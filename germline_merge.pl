#!/usr/bin/perl
# Merge a blood variant report with a tumor variant report to determine if 
# variants are germline or somatic.
#
# TODO:
#    - Add some output formatting options, including TSV / PP options.
#    - Add filtering to pull up update only by position, gene, ID, etc. Can 
#      leverage vcfExtractor filters to limit data if need be.
#    - Add dbSNP annotations?
#    - Add SNP only filtering, which can be turned into an ID verify kind of 
#      thing if we really want. 
#
# 12/7/2017 - D Sims
###############################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Data::Dump;
use Sort::Versions;

my $scriptname = basename($0);
my $version = "v0.4.030118";
my $vaf_diff = 20;

my $description = <<"EOT";
Input a blood VCF and a tumor VCF, and get a comparison table of VAFs between the two.
May be helpful for determining which variants are somatic and which are germline.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options]
    -b, --blood     Blood VCF to be processed.
    -t, --tumor     Tumor VCF to be processed.

    -d, --diff      VAF difference (as an INT) threshold for annotation (default: $vaf_diff%).
    -o, --output    Send output to custom file.  Default is STDOUT.
    -v, --version   Version information
    -h, --help      Print this help information
EOT

my $help;
my $ver_info;
my $outfile;
my $blood_vcf;
my $tumor_vcf;

GetOptions( 
    "blood|b=s"     => \$blood_vcf,
    "tumor|t=s"     => \$tumor_vcf,
    "diff|d=i"      => \$vaf_diff,
    "output|o=s"    => \$outfile,
    "version|v"     => \$ver_info,
    "help|h"        => \$help )
or die $usage;

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help if $help;
version if $ver_info;

if (! ($blood_vcf and $tumor_vcf)) {
    print "ERROR: you must input both a blood VCF and tumor VCF file to be ",
        "processed!\n";
    print "$usage\n";
    exit 1;
}

# Make sure enough args passed to script
#if ( scalar( @ARGV ) < 1 ) {
	#print "ERROR: Not enough arguments passed to script!\n\n";
	#print "$usage\n";
	#exit 1;
#}

# Verify that we have VCF Extractor in our path. should be a non-issue if one
# gets whole repo with this script, but just in case.
unless ( qx(which vcfExtractor.pl) ) {
    print "ERROR: You must have vcfExtractor.pl installed and in your \$PATH. This ",
        "should already be a part of the biofx repo (https://github.com/drmrgd/biofx-utils) ",
        "from which this script comes.\n";
    exit 1;
}

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
    print "Writing output to $outfile...\n";
	open( $out_fh, ">", $outfile ) 
        || die "Can't open the output file '$outfile' for writing: $!";
} else {
	$out_fh = \*STDOUT;
}
#########--------------------- END ARG Parsing ------------------------#########

my $blood_data = parse_vcf($blood_vcf);
my $tumor_data = parse_vcf($tumor_vcf);
my $married_data = marry_results($blood_data, $tumor_data);

# Add blood only, tumor only, etc annotations to data.
find_uniq($married_data);

print_results($married_data, $out_fh);

sub find_uniq {
    # Read the married results in, and determine if each variant is germline (i.e.
    # the same call and roughly the same VAF (+/- 5% VAF), or if only in tumor. 
    # Also flag variants that are in blood only, which is probably indicative of 
    # artifact.
    my $data = shift;
    my %results;

    for my $var (keys %$data) {
        my ($bvaf, $tvaf) = @{$$data{$var}}[3,5];
        if ($bvaf eq 'ND') {
            splice(@{$$data{$var}}, 7, 0, 'Tumor Only');
        }
        elsif ($tvaf eq 'ND') {
            splice(@{$$data{$var}}, 7, 0, 'Blood Only');
        }
        elsif (abs($bvaf - $tvaf) > $vaf_diff) {
            splice(@{$$data{$var}}, 7, 0, 'Diff');
        } else {
            splice(@{$$data{$var}}, 7, 0, '-');
        }
    }
    return;
}

sub parse_vcf {
    my $vcf = shift;
    my %data;
    my @wanted_fields = qw(position ref alt vaf cov refcov altcov varid gene transcript
        cds aa location function gene_class variant_class);
    # TODO: Need to verify we have an IR annotated VCF or this is not going to work.
    open(my $vcf_data, "-|", "vcfExtractor.pl -Nna $vcf")
       or die "ERROR: Can't parse VCF file.\n";

   while(<$vcf_data>) {
       chomp;
       next unless /^chr/;
       my @fields = split;
       my $varid = join(':', @fields[0..2]);
       @{$data{$varid}}{@wanted_fields} = @fields;
   }
   return \%data;
}

sub marry_results {
    my ($blood_data, $tumor_data) = @_;
    my %married_data;
    my @wanted = qw(position ref alt vaf cov transcript gene cds aa varid location 
        function);
        
    # Work with data in context of tumor, and then got back to verify we didn't
    # miss any blood results.  Presumably tumor will be more.
    for my $pos (keys %$tumor_data) {
        my @outdata = @{$$tumor_data{$pos}}{@wanted};
        if ($blood_data->{$pos}) {
            splice(@outdata, 3, 0, @{$$blood_data{$pos}}{'vaf'}, @{$$blood_data{$pos}}{'cov'});
        } else {
            splice(@outdata, 3, 0, 'ND', 'ND');
        }
        $married_data{$pos} = \@outdata;
    }

    # Check to see if we need to add any blood variants that didn't get picked
    # up in tumor for some reason.
    my @missing = grep{ ! $married_data{$_} } keys %$blood_data;
    if (@missing) {
        for my $pos (@missing) {
            my @outdata = @{$$blood_data{$pos}}{@wanted};
            splice(@outdata, 5, 0, 'ND', 'ND');
            $married_data{$pos} = \@outdata;
        }
    }
    return \%married_data;
}

sub print_results {
    my ($data, $outfh) = @_;
    my @header = qw(Position REF ALT Blood_VAF Blood_Cov Tumor_VAF Tumor_Cov Annot Transcript Gene 
        CDS AA VarID Location Function);
    select $outfh;
    print join(',', @header), "\n";
    for my $var (sort {versioncmp($a, $b)} keys %$data) {
        print join(',', @{$$data{$var}}), "\n";
    }
}
