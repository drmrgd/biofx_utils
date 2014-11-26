#!/usr/bin/perl
# Read in all the variant call files, filter and truncate the data, and associate with the run ID.  Set
# up to be able to find the common variants in all the runs.  Add column to indicate if it's one of the 
# variants located within the NIST CEPH standard VCF file.
# 
# TODO:
#   - Clean up output file generation.  Don't make an output file unless it's needed!
#   - Better error checking for input files. erroneous output files from frist todo list item are 
#     breaking the program when globbed in
#
# Created: 2/27/2013 - Dave Sims
#
##########################################################################################################
use warnings;
use strict;
use autodie;

use List::Util qw{ min max sum };
use Getopt::Long qw{ :config auto_abbrev bundling no_ignore_case };
use Sort::Versions;
use File::Basename;
use Data::Dump;

my $scriptname = basename($0);
my $version = "v4.4.0_112614";
my $description = <<EOT;
Program to read in all of the variant call table files from an Ion Torrent run on CEPH, and report out
the ID and number of times each variant is seen.  This is used to track the number of variants reported 
in that sample in order to identify variants that should be excluded from the reportable range of the 
assay.
EOT

my $usage = <<EOT;
USAGE: $scriptname [options] <input_files>
    -C, --Classic    Use TVCv3.2.1 formatted data files.  Will be deprecated.
    -V, --VCF        Input is a VCF file instead of TVC Tabular File
    -c, --coverage   Coverage cutoff (default is 450)
    -p, --preview    Write output only to STDOUT to preview the results.
    -o, --output     Write output to custom file (DEFAULT: "CEPH_###_Run_Variant_Tally.tsv")
    -v, --version    Print version information
    -h, --help       Print help information
EOT

my $verinfo;
my $help;
my $covfilter = 450;
my $preview;
my $output;
my $classic;
my $vcf_input;

# Set up some commandline opts
GetOptions( "Classic|C"       => \$classic,
            "VCF|V"           => \$vcf_input,
            "coverage|c=i"    => \$covfilter,
            "preview|p"       => \$preview,
            "output|o=s"      => \$output,
            "version|v"       => \$verinfo,
            "help|h"          => \$help,
    ) or die "\n$usage";

sub help {
	printf "%s - %s\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version {
    printf "%s - %s\n", $scriptname, $version;
    exit;
}

version if $verinfo;
help if $help;

my @filesList;
my $total_runs;
if ( ! @ARGV ) {
    print "ERROR: no files loaded for analysis\n";
    exit 1;
} else {
    @filesList = @ARGV;
    $total_runs = scalar @filesList;
}

# Set up Tally output file.
#TODO: Clean this up.  Don't generate output files unless necessary.
my ($out_fh, $outfile);
($output) ? ($outfile = $output) : ($outfile = "CEPH_".$total_runs."_Run_Variant_Tally.tsv");
if ( $preview ) {
    $out_fh = \*STDOUT;
} else {
    open( $out_fh, ">", $outfile ) || die "Can't open the output file '$outfile' for writing: $!";
    print "Writing data to '$outfile'\n";
}

######--------------------------------- END COMMANDLINE ARGS ----------------------------------------######
my ( %all_variants, %var_freq, %var_cov );

# Set up fields to use to capture data from the tables.  Set up here so that we can mod it easy later when the 
# tables change.
my %fields = (
    "v3.2.1"  => { "varid" => [qw(0 1 6 7)],
                   "data"  => [qw(0 1 3 6 7 8 10)],
                 },
    "v4.0.2"  => { "varid" => [qw(0 1 15 16)],
                   "data"  => [qw(0 1 13 15 16 6 18)],
                 },
    # Add in VCF entry; use data retrieved from vcfExtractor
    "vcf"    => { "varid"  => [qw()],
                  "data"   => [qw()],
              },
);

my $ts_version;
if ($classic) {
    $ts_version = "v3.2.1";
}
elsif ($vcf_input) {
    $ts_version = "VCF";
}
else {
    $ts_version = "v4.0.2";
}

if ( $vcf_input ) {
    #print "Processing data as VCF file...\n";
    proc_vcf( \@filesList, \$ts_version );
} else {
    proc_datatable( \@filesList, \%fields, \$ts_version );
}

# Get some field width data
my ( $rwidth, $awidth ) = field_width( \%all_variants );

# Get statistics about each variant and generate a formated hash table to print out the results with 
my %results;
foreach my $variant ( keys %all_variants) {
	my $count = @{$all_variants{$variant}};
    my ( $chr, $pos, $ref, $alt ) = split( /:/, $variant );
	
	# Get min, max, median coverage and frequency info	
	my ( $minCov, $maxCov, $meanCov ) = stats( \@{$var_cov{$variant}} );
	my ( $minFreq, $maxFreq, $meanFreq ) = stats( \@{$var_freq{$variant}} );
    my $detection_freq = sprintf( "%0.2f", ($count/$total_runs)*100);

    $results{$variant} = [$chr, $pos, $ref, $alt, "$count/$total_runs", $minCov, $meanCov, $maxCov, $minFreq, $meanFreq, $maxFreq];

}

# Print out the collected summary data 
my @header_cols = qw{ Chr Position Ref Var Count MinCov MedCov MaxCov MinVAF MedVAF MaxVAF };
my $header = sprintf( "%-8s %-12s %-${rwidth}s %-${awidth}s %-10s %-7s %-7s %-7s %9s %9s %9s", @header_cols ); 
my $title = "Frequency of detected variants with at least $covfilter reads in $total_runs CEPH runs";

print $out_fh "$title\n\n";
print $out_fh "$header\n";

my $format = "%-8s %-12d %-${rwidth}s %-${awidth}s %-10s %-7d %-7d %-7d %8.2f%% %8.2f%% %8.2f%%\n";
for my $variant ( sort { 
        versioncmp( $results{$b}->[4], $results{$a}->[4] ) || 
        versioncmp( $results{$a}->[0], $results{$b}->[0] ) || 
        versioncmp( $results{$a}->[1], $results{$b}->[1] ) 
    } keys %results ) {
        printf $out_fh $format, @{$results{$variant}};
}

sub proc_datatable {
    my $files = shift;
    my $data_fields = shift;
    my $version = shift;

    # Get fields to use for spliting the tables
    print "Processing '$$version' data...\n";
    my @varid_index = @{$$data_fields{$$version}{'varid'}};
    my @field_index = @{$$data_fields{$$version}{'data'}};

    for my $file ( @$files ) {
        open( my $in_fh, "<", $file );
        my $header = <$in_fh>;
        if ( $header =~ /^#+.*VCF/ ) {
            print "ERROR: file '$file' appears to be a VCF file.  You should use the -V option to process.\n";
            exit 1;
        }
        while (<$in_fh>) {
            next if ( /Chrom/ || /Absent/ || /No Call/ );
            my @fields = split;
            if ( $fields[$field_index[6]] > $covfilter ) {
                my $varid = join( ':', @fields[@varid_index] );
                push( @{$all_variants{$varid}}, [@fields[@field_index]] );
                push( @{$var_freq{$varid}}, $fields[$field_index[5]] );
                push( @{$var_cov{$varid}}, $fields[$field_index[6]] );
            }
        }
        close $in_fh;
    }
    return;
}

sub proc_vcf {
    # Process VCF file input.  vcfExtractor program required
    my $files = shift;
    my $version = shift;

    print "Processing '$$version' data...\n";

    # Check to see that we have vcfExtractor in our $PATH
    my $vcfextractor_path;
    unless ( ($vcfextractor_path = qx(which vcfExtractor)) || ($vcfextractor_path = qx(which vcfExtractor.pl)) ) {
        print "ERROR: 'vcfExtractor' required when inputting VCF files, but program not found.  Check to see that this is in your path\n";
        exit 1;
    }
    chomp $vcfextractor_path;

    #if ( ! `which vcfExtractor` && ! `which vcfExtractor.pl` ) {
        #print "ERROR: 'vcfExtractor' required when inputting VCF files, but program not found.  Check to see that this is in your path\n";
        #exit 1;
    #}

    for my $file ( @$files ) {
         my $cmd = qq($vcfextractor_path --noref --NOCALL $file 2>/dev/null);
         open( my $data_fh, "-|", $cmd ) || die "Can't open the stream: $!";
         while (<$data_fh>) {
             next until ( /^chr/ );
             my @fields = split;
             if ( $fields[6] > $covfilter ) {
                 my $varid = join( ':', @fields[0,1,2] );
                 push( @{$all_variants{$varid}}, [$fields[0],'---',@fields[1,2,5,6]] );
                 push( @{$var_freq{$varid}}, $fields[5] );
                 push( @{$var_cov{$varid}}, $fields[6] );
             }
         }
         close $data_fh;
     }
     return;
}

sub stats {
    # Get some stats on variant metrics
	my $input = shift;
	my $min = min( @$input );
	my $max = max( @$input );
	my $sum = sum( @$input );
	my $mean = $sum/@$input;
	return ( $min, $max, $mean );
}

sub field_width {
    #set dynamic field with for formatting later.
    my $hash_ref = shift;

    my $rwidth = 0;
    my $awidth = 0;

    for my $var ( keys %$hash_ref ) {
        my ($chr, $start, $ref, $alt) = split( /:/, $var );
        $rwidth = length( $ref ) + 3  if ( length( $ref ) > $rwidth ); 
        $awidth = length( $alt ) + 3  if ( length( $alt ) > $awidth );
    }
    return( $rwidth, $awidth );
}
