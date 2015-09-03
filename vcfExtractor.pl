#!/usr/bin/perl
# Script to pull out column information from a VCF file.  Can also grab variant information
# based on a position lookup using the '-p' option.  
#
# Need to make sure to install the latest version of VCF Tools to avoid generic Perl error 
# message in output.  Can build from source, or I built a .deb file to installed it on 
# Debian systems.
# 
# D Sims 2/21/2013
#################################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use List::Util qw( sum min max );
use Sort::Versions;
use JSON::XS; 
use Data::Dump;
use File::Basename;
use Term::ANSIColor;

use constant 'DEBUG' => 0;
print colored("*" x 50, 'bold yellow on_black'), "\n";
print colored("\t\tDEVELOPMENT VERSION", 'bold yellow on_black'), "\n";
print colored("*" x 50, 'bold yellow on_black'), "\n\n";

my $scriptname = basename($0);
my $version = "v4.5.4_090315-dev";
my $description = <<"EOT";
Program to extract fields from an Ion Torrent VCF file generated by TVCv4.0+.  By default the program 
will extract the following fields:

     CHROM:POS REF ALT Filter Filter_Reason VAF TotCov RefCov AltCov COSID

This can only be modified currently by the hardcoded variable '\$vcfFormat'.

This version of the program also supports extracting only variants that match a position query based
on using the following string: chr#:position. Additionally, Hotspot annotated variants (i.e. those
variants that have a COSMIC, OM, or other annotation in the TVC output), can be searched using the
Hotspot ID (e.g. COSM476).  

Multiple positions can be searched by listed each separated by a space and wrapping the whole query in
quotes:

        vcfExtractor -p "chr17:29553485 chr17:29652976" <vcf_file>

For batch processing, a lookup file with the positions (one on each line in the same format as
above) can be passed with the '-f' option to the script:

        vcfExtractor -l lookup_file <vcf_file>
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] [-f {1,2,3}] <input_vcf_file>

    Program Options
    -o, --output    Send output to custom file.  Default is STDOUT.
    -a, --annot     Add IR and Oncomine OVAT annotation information to output if available.
    -V, --Verbose   Print additional information during processing.
    -v, --version   Version information
    -h, --help      Print this help information

    Filter and Output Options
    -p, --pos       Output only variants at this position.  Format is "chr<x>:######" 
    -c, --cosid     Look for variant with matching COSMIC ID (or other Hotspot ID)
    -l, --lookup    Read a list of variants from a file to query the VCF. 
    -f, --fuzzy     Less precise (fuzzy) position match. Strip off n digits from the position string.
                    MUST be used with a query option (e.g. -p, -c, -l), and can not trim more than 3 
                    digits from string.
    -n, --noref     Output reference calls.  Ref calls filtered out by default
    -N, --NOCALL    Remove 'NOCALL' entries from output
    -O, --OVAT      Only report Oncomine Annotated Variants.
    -H, --HS        Print out only variants that have a Hotspot ID.
EOT

my $help;
my $ver_info;
my $outfile;
my $positions;
my $lookup;
my $fuzzy;
my $noref;
my $nocall;
my $hsids;
my $annots;
my $ovat_filter;
my $verbose;
my $hotspot;

GetOptions( "output|o=s"    => \$outfile,
            "annot|a"       => \$annots,
            "OVAT|O"        => \$ovat_filter,
            "cosid|c=s"     => \$hsids,
            "NOCALL|N"      => \$nocall,
            "pos|p=s"       => \$positions,
            "lookup|l=s"    => \$lookup,
            "fuzzy|f=i"     => \$fuzzy,
            "noref|n"       => \$noref,
            "version|v"     => \$ver_info,
            "Verbose|V"     => \$verbose,
            "HS|H"          => sub{ print colored( "Filtering by hotspots is not yet implemented. skipping...\n\n", "bold yellow on_black") ; return },
            #"HS|H"          => \$hotspot,
            "help|h"        => \$help )
        or die "\n$usage";

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

# Set up some colored output flags and warn / error variables
my $warn  = colored( "WARN:", 'bold yellow on_black');
my $err   = colored( "ERROR:", 'bold red on_black');
my $info  = colored( "INFO:", 'bold cyan on_black');
my $on    = colored( "On", 'bold green on_black');
my $off   = colored( "Off", 'bold red on_black');
my @warnings;

# Check for vcftools; we can't run without it...for now.
if ( ! qx(which vcftools) ) {
    print "$err Required package vcftools is not installed on this system.  Install vcftools ('vcftools.sourceforge.net') and try again.\n";
    exit 1;
}

# Double check that fuzzy option is combined intelligently with a position lookup.
if ( $fuzzy ) {
    if ( $fuzzy > 3 ) {
        print "\n$err Can not trim more than 3 digits from the query string.\n\n";
        print $usage;
        exit 1;
    }
    elsif ( $lookup ) {
        print "\n$warn fuzzy lookup in batch mode may produce a lot of results! Continue? ";
        chomp( my $response = <STDIN> );
        exit if ( $response =~ /[(n|no)]/i );
        print "\n";
    }
    elsif ( ! $positions && ! $hsids ) {
        print "$err must include position or hotspot ID query with the '-f' option\n\n";
        print $usage;
        exit 1;
    }
}

# Throw a warning if using the ovat filter without asking for OVAT annotations.
if ( $ovat_filter && ! $annots ) {
    print "$info Requested Oncomine annotation filter without adding the OVAT annotations. Auto adding the OVAT annotations!\n" if $verbose;
    $annots=1;
}

# Make sure enough args passed to script
if ( scalar( @ARGV ) < 1 ) {
    print "$err No VCF file passed to script!\n\n";
    print $usage;
    exit 1;
}

# Parse the lookup file and add variants to the postions list if processing batch-wise
my @valid_hs_ids = qw( BT COSM OM OMINDEL MCH );
if ($lookup) {
    my $query_list = batch_lookup(\$lookup) if $lookup;
    if ( grep { $$query_list =~ /$_\d+/ } @valid_hs_ids ) {
        $hsids = $$query_list;
    } 
    elsif ( $$query_list =~ /chr[0-9YX]+/ ) {
        $positions = $$query_list;
    }
    else {
        print "$err Issue with lookup file.  Check and retry\n"; 
        exit 1;
    }
}

# Double check the query (position / cosid) format is correct
my (@coords, @cosids);
if ( $positions ) {
    print "Checking position id lookup...\n" if $verbose;
    @coords = split( /\s+/, $positions );
    for my $coord ( @coords ) {
        if ( $coord !~ /\Achr[0-9YX]+:\d+$/i ) {
            print "$err '$coord' not valid. Please use the following format for position queries: 'chr#:position'\n";
            exit 1;
        }
    }
} 
elsif ( $hsids ) {
    print "Checking variant id lookup...\n" if $verbose;
    @cosids = split( /\s+/, $hsids );
    for my $varid ( @cosids ) { 
        if ( ! grep { $varid =~ /$_\d+/ } @valid_hs_ids ) {
            print "$err '$varid' is not a valid hotspot query term! Valid lookups are:\n";
            print "\t${_}###\n" for @valid_hs_ids;
            exit;
        }
    }
}

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
	open( $out_fh, ">", $outfile ) || die "Can't open the output file '$outfile' for writing: $!";
} else {
	$out_fh = \*STDOUT;
}
#########------------------------------ END ARG Parsing ---------------------------------#########
my $inputVCF = shift;

# Check VCF file and options to make sure they're valid
open ( my $vcf_fh, "<", $inputVCF );
my @header = grep { /^#/ } <$vcf_fh>;
if ( $header[0] !~ /VCFv4/ ) {
    print "$err '$inputVCF' does not appear to be a valid VCF file or does not have a header.\n\n";
    print "$usage\n";
    exit 1;
}

# Crude check for TVC3.2 or TVC4.0+ VCF file.  Still need to refine this
if ( grep { /^##INFO.*Bayesian_Score/ } @header ) {
    print "$warn '$inputVCF' appears to be from TVCv3.2, and the 'tvc32' option was not selected.  The file may not be processed correctly.\n";
    die "Pre TVCv4.0 VCF file detected.  These files are not longer supported by this utility\n";
}

# Trigger OVAT annot capture if available
my $ovat_annot;
( grep { /OncomineVariantAnnotation/ } @header ) ? ($ovat_annot = 1) : ($ovat_annot = 0);
die "$err OVAT output selected, but no OVAT annotations found in VCF file!\n" if ( $annots && $ovat_annot == 0 ); 
close $vcf_fh;

# Get the data from VCF Tools
my $vcfFormat;
if ($annots) {
    $vcfFormat = "'%CHROM:%POS\t%REF\t%ALT\t%FILTER\t%INFO/FR\t%INFO/OID\t%INFO/OPOS\t%INFO/OREF\t%INFO/OALT\t%INFO/OMAPALT\t%INFO/FUNC\t[%GTR\t%AF\t%FRO\t%RO\t%FAO\t%AO\t%DP]\n'";
} else {
    $vcfFormat = "'%CHROM:%POS\t%REF\t%ALT\t%FILTER\t%INFO/FR\t%INFO/OID\t%INFO/OPOS\t%INFO/OREF\t%INFO/OALT\t%INFO/OMAPALT\t---\t[%GTR\t%AF\t%FRO\t%RO\t%FAO\t%AO\t%DP]\n'";
}

my @extracted_data = qx/ vcf-query $inputVCF -f $vcfFormat /;

# Read in the VCF file data and create a hash
my %vcf_data = parse_data( \@extracted_data );

my $ovat_stat;
($ovat_filter) ? ($ovat_stat = "$on!") : ($ovat_stat = "$off.");
print "$info OVAT filter status: $ovat_stat\n" if $verbose;

# Filter and format extracted data or just format and print it out.
if ( @cosids ) {
    filter_data(\%vcf_data, \@cosids);
}
elsif ( @coords ) {
    filter_data(\%vcf_data, \@coords); 
} 
elsif ( $ovat_filter ) {
    ovat_filter(\%vcf_data);
}
else {
    format_output(\%vcf_data);
}

sub parse_data {
    # Extract the VCF information and create a hash of the data.  
    my $data = shift;
    my %parsed_data;

    for ( @$data ) {
        my ( $pos, $ref, $alt, $filter, $reason, $oid, $opos, $oref, $oalt, $omapalt, $func, $gtr, $af, $fro, $ro, $fao, $ao, $dp ) = split( /\t/ );

        # IR generates CNV and Fusion entries that are not compatible.  
        # Skip for now; implement a sub for each later.
        next if ( $alt =~ /[.><\]\d+]/ ); 

        # Clean up filter reason string
        $reason =~ s/^\.,//;

        # Filter out vars we don't want to print out later anyway.
        next if $reason eq "NODATA";  # Don't print out NODATA...nothing to learn there.
        $filter = "NOCALL" if ( $gtr =~ m|\./\.| );
        next if ( $nocall && $filter eq "NOCALL" );
        next if ( $noref && $gtr eq '0/0' );

        # Create some arrays to hold the variant data in case we have MNP calls here
        my @alt_array     = split( /,/, $alt );
        my @oid_array     = split( /,/, $oid );
        my @opos_array    = split( /,/, $opos );
        my @oref_array    = split( /,/, $oref );
        my @oalt_array    = split( /,/, $oalt );
        my @omapalt_array = split( /,/, $omapalt );
        my @fao_array     = split( /,/, $fao );
        my @ao_array      = split( /,/, $ao );

        my @indices;
        for my $alt_index ( 0..$#alt_array ) {
            my $alt_var = $alt_array[$alt_index];

            # Get the normalizedRef, normalizedAlt, and normalizedPos values from the REF and ALT fields so that
            # we can map the FUNC block.
            my @coords = split(/:/, $pos);
            my %norm_data = normalize_variant(\$ref, \$alt_var, $coords[1]);

            my @array_pos = grep { $omapalt_array[$_] eq $alt_var } 0..$#omapalt_array;
            for my $index ( @array_pos ) {
                (my $parsed_pos = $pos) =~ s/(chr\d+:).*/$1$opos_array[$index]/; 
                
                # Stupid bug with de novo and hotspot merge that can create two duplicate entries for the same
                # variant but one with and one without a HS (also different VAF, coverage,etc). Try this to 
                # capture only HS entry.
                #my $var_id = join( ":", $parsed_pos, $oref_array[$index], $oalt_array[$index] );
                my $var_id = join( ":", $parsed_pos, $oref_array[$index], $oalt_array[$index], $filter );
                my $cosid = $oid_array[$index];
                if ( $cosid ne '.' && exists $parsed_data{$var_id} ) {
                   delete $parsed_data{$var_id}; 
                }

                # Grab the OVAT annotation information from the FUNC block if possible.
                my ($ovat_gc, $ovat_vc, $gene_name, $hgvs);
                if ( $func eq '.' ) {
                    push( @warnings, "$warn could not find FUNC entry for '$pos'\n") if $annots;
                    $ovat_vc = $ovat_gc = $gene_name = "NULL";
                } else {
                    ($ovat_gc, $ovat_vc, $gene_name, $hgvs) = get_ovat_annot(\$func, \%norm_data) unless $func eq '---'; 
                }

                # Start the var string.
                push( @{$parsed_data{$var_id}},
                    $parsed_pos,
                    $norm_data{'normalizedRef'},
                    $norm_data{'normalizedAlt'},
                    $filter, 
                    $reason, 
                    $gtr );

                my ($vaf, $tot_coverage);
                # Check to see if call is result of long indel assembler and handle appropriately. 
                if ( $fao_array[$alt_index] eq '.' ) {
                    $tot_coverage = $ao_array[$alt_index] + $ro;
                    $vaf = vaf_calc( \$filter, \$dp, \$ro, \$ao_array[$alt_index] );
                    push(@{$parsed_data{$var_id}}, $vaf, $tot_coverage, $ro, $ao_array[$alt_index], $cosid);
                } else {
                    my @cleaned_fao_array = grep { $_ ne '.' } @fao_array;
                    $tot_coverage = sum( @cleaned_fao_array ) + $fro;
                    $vaf = vaf_calc( \$filter, \$tot_coverage, \$fro, \$fao_array[$alt_index] );
                    push( @{$parsed_data{$var_id}}, $vaf, $tot_coverage, $fro, $fao_array[$alt_index], $cosid );
                }

                if ($noref) {
                    if ( ${$parsed_data{$var_id}}[6] ne '.' ) {
                        if ( ${$parsed_data{$var_id}}[6] < 1 ) {
                                delete $parsed_data{$var_id};
                                next;
                        }
                    }
                }
                push( @{$parsed_data{$var_id}}, $gene_name, $hgvs, $ovat_gc, $ovat_vc ) if $annots;
            }
        }
    }
    #dd \%parsed_data;
    #exit;
    return %parsed_data;
}

sub get_ovat_annot {
    # If this is IR VCF, add in the OVAT annotation. 
    no warnings;
    my $func = shift; 
    my $norm_data = shift; 
    $$func =~ tr/'/"/;

    my %data = (
        oncomineGeneClass     => '---',
        oncomineVariantClass  => '---',
        gene                  => '---',
        transcript            => '---',
        protein               => '---',
        coding                => '---',
        function              => '---',
        normalizedRef         => '---',
        normalizedAlt         => '---',
        location              => '---',
    );

    my $json_annot = JSON::XS->new->decode($$func);

    # Create some default mappings since only 1 block per entry (unless MNP) and gene, location, etc is the same.
    $data{'transcript'} = $$json_annot[0]->{'transcript'};
    $data{'gene'}       = $$json_annot[0]->{'gene'};
    $data{'location'}   = $$json_annot[0]->{'location'};

    my $match;
    for my $func_block ( @$json_annot ) {
        # if there is a normalizedRef entry, then let's map the func block appropriately
        if ($$func_block{'normalizedRef'}) {
            if ($$func_block{'normalizedRef'} eq $$norm_data{'normalizedRef'} && $$func_block{'normalizedAlt'} eq $$norm_data{'normalizedAlt'}) {
                %data = %$func_block;
                $match = 1;
                last;
            } 
        } 
    }

    my $gene_class    = $data{'oncomineGeneClass'}    // '---';
    my $variant_class = $data{'oncomineVariantClass'} // '---';
    my $gene_name     = $data{'gene'}                 // '---';
    my $protein       = $data{'protein'}              // '---';
    my $hgvs          = $data{'coding'}               // '---';
    my $function      = $data{'function'}             // '---';
    my $ref           = $data{'normalizedRef'}        // '---';
    my $alt           = $data{'normalizedAlt'}        // '---';
    my $exon;
    ($data{'location'} eq 'exonic') ? ($exon = $data{'exon'}) : ($exon = 'intronic');

    if (DEBUG) {
        print "======================  DEBUG  =======================\n\n";
        print "gc       => $gene_class\n";
        print "vc       => $variant_class\n";
        print "gene     => $gene_name\n";
        print "ref      => $ref\n";
        print "alt      => $alt\n";
        print "AA       => $protein\n";
        print "HGVS     => $hgvs\n";
        print "function => $function\n";
        print "exon     => $exon\n";
        print "======================================================\n\n";
    }
    return ($gene_class, $variant_class, $gene_name, $hgvs);
}

sub normalize_variant {
    # Borrowed from ThermoFisher's vcf.py script to convert IR VCFs. Trim from both ends until only unique
    # sequence left
    my $ref = shift;
    my $alt = shift;
    my $pos = shift;
    my ($norm_ref, $norm_alt);

    my ($rev_ref, $rev_alt, $position_delta) = rev_and_trim($ref, $alt);
    ($norm_ref, $norm_alt, $position_delta) = rev_and_trim(\$rev_ref, \$rev_alt);

    my $adj_position = $position_delta + $pos;
    return ( 'normalizedRef' => $norm_ref, 'normalizedAlt' => $norm_alt, 'normalizedPos' => $adj_position );
}

sub rev_and_trim {
    # Borrowed from ThermoFisher's vcf.py script to convert IR VCFs
    my $ref = shift;
    my $alt = shift;

    my $position_delta = 0;
    my @rev_ref = split(//, reverse($$ref));
    my @rev_alt = split(//, reverse($$alt));

    while (@rev_ref > 1 && @rev_alt > 1 && $rev_ref[0] eq $rev_alt[0]) {
        shift @rev_ref;
        shift @rev_alt;
        $position_delta++;
    }
    return (join('',@rev_ref), join('', @rev_alt), $position_delta);
}

sub vaf_calc {
    # Determine the VAF
    my $filter = shift;
    my $tcov = shift;
    my $rcov = shift;
    my $acov = shift;

    my $vaf;

    if ( $$filter eq "NOCALL" ) { 
        $vaf = '.';
    }
    elsif( $$filter eq "NODATA" || $$tcov == 0) {
        $vaf = 0;
    }
    else {
        $vaf = sprintf( "%.2f", 100*($$acov / $$tcov) );
    }
    return $vaf;
}

sub filter_data {
    # Filtered out extracted dataset.
    my $data = shift;
    my $filter = shift;

    my %filtered_data;
    my @fuzzy_pos;
    my %counter;

    if ( $fuzzy ) {
        my $re = qr/(.*).{$fuzzy}/;
        @fuzzy_pos = map { /$re/ } @$filter;
        for my $query ( @fuzzy_pos ) {
            for ( sort keys %$data ) {
                next if ( $ovat_filter && $$data{$_}[11] eq '---' );
                if ( $$data{$_}[0] =~ /$query.{$fuzzy}/ ) {
                    push( @{$filtered_data{$query}},  [@{$$data{$_}}] );
                    $counter{$query} = 1;
                }
            }
        }
    } 
    else {
        for my $variant ( keys %$data ) {
            if ($hsids) {
                if ( my ($query) = grep { ($_) eq $$data{$variant}[10] } @$filter ) {
                    @{$filtered_data{$variant}} = @{$$data{$variant}};
                    $counter{$query} = 1;
                }
            } else {
                if ( my ($query) = grep { ($_) eq $$data{$variant}[0] } @$filter ) {
                    @{$filtered_data{$variant}} = @{$$data{$variant}};
                    $counter{$query} = 1;
                }
            }
        }
    }

    # Send results to output report formatter
    ($ovat_filter) ? ovat_filter(\%filtered_data) : format_output(\%filtered_data);

    my $term;
    ($hsids) ? ($term = "with Hotspot ID:") : ($term = "at position:");

    if ( $fuzzy ) {
        for my $query ( @fuzzy_pos ) {
            my $string = $query . ( '*' x $fuzzy );
            printf $out_fh "\n>>> No variant found $term %s <<<\n", $string if ( ! exists $counter{$query} );
        }
    } else {
        for my $query ( @$filter ) {
            print $out_fh "\n>>> No variant found $term $query <<<\n" if ( ! exists $counter{$query} );
        } 
    }
}

sub get_hotspots {
    my $data = shift;

    print "Filtering by Hotspot ID is not yet implemented!  Skipping this....\n";
    return;

    print "$info Filtering out non-hotspots\n" if $verbose;
    for my $variant (keys %$data ) {
        delete $$data{$variant} if $$data{$variant}->[9] eq '.';
    }

    format_output($data);
    print {$out_fh} "\n>>> No Hotspots detected! <<<\n" unless $data;
}

sub ovat_filter {
    # Filter out calls that are not oncomine reportable variants
    my $data = shift;

    # Add in the OVAT filter here
    if ( $ovat_filter ) {
        print "$info Running ovat filter\n" if $verbose;
        for my $variant ( keys %$data ) {
            delete $$data{$variant} if $$data{$variant}->[12] =~ /^(-|NULL)/;
        }
    }

    format_output($data);
    print {$out_fh} "\n>>> No Oncomine Annotated Variants Found! <<<\n" unless %$data;
}

sub format_output {
    # XXX
    # Format and print out the results
    my $data = shift;
    my ($w1, $w2, $w3, $w4) ;
    # Set up the output header
    my ($format, @header);
    if ($annots) {
        ($w1, $w2, $w3, $w4) = field_width( $data );
        $format = "%-19s %-${w1}s %-${w2}s %-10s %-${w3}s %-10s %-10s %-10s %-10s %-${w4}s %-14s %-12s %-21s %s\n";
        @header = qw( CHROM:POS REF ALT Filter Filter_Comment VAF TotCov RefCov AltCov COSID Gene HGVS oncomineGeneClass oncomineVariantClass );
    } else {
        ($w1, $w2, $w3) = field_width( $data );
        $format = "%-19s %-${w1}s %-${w2}s %-10s %-${w3}s %-10s %-10s %-10s %-10s %s\n";
        @header = qw( CHROM:POS REF ALT Filter Filter_Commment VAF TotCov RefCov AltCov COSID );
    }
    printf {$out_fh} $format, @header;

    # Need to parse the data stucture differently if fuzzy lookup
    if ( $fuzzy ) {
        for my $variant ( sort { versioncmp( $a, $b ) }  keys %$data ) {
            for my $common_var ( @{$$data{$variant}} ) {
                printf $out_fh $format, @$common_var[0..4,6..14];
            }
        }
    } else {
        for my $variant ( sort { versioncmp( $a, $b ) } keys %$data ) {
            printf {$out_fh} $format, @{$$data{$variant}}[0..4,6..14];
        }
    }

    # Wrap up
    if ( @warnings && $verbose ) {
        print "\n";
        print $_ for @warnings;
    }
}

sub field_width {
    # Get the longest field width for formatting later.
    my $data_ref = shift;
    my $ref_width = 0;
    my $var_width = 0;
    my $filter_width= 0;
    my $hgvs_width = 0;

    if ( $fuzzy ) {
        for my $variant ( keys %$data_ref ) {
            for ( @{$$data_ref{$variant}} ) {
                my $ref_len = length( $$_[1] );
                my $alt_len = length( $$_[2] );
                my $filter_len = length( $$_[4] );
                my $hgvs_len = length($$_[10]);
                $ref_width = $ref_len if ( $ref_len > $ref_width );
                $var_width = $alt_len if ( $alt_len > $var_width );
                $filter_width = $filter_len if ( $filter_len > $filter_width );
                $hgvs_width = $hgvs_len if ( $hgvs_len > $hgvs_width );
            }
        }
    } else {
        for my $variant ( keys %$data_ref ) {
            my $ref_len = length( $$data_ref{$variant}[1] );
            my $alt_len = length( $$data_ref{$variant}[2] );
            my $filter_len = length( $$data_ref{$variant}[4] );
            my $hgvs_len = length( $$data_ref{$variant}[10] );
            $ref_width = $ref_len if ( $ref_len > $ref_width );
            $var_width = $alt_len if ( $alt_len > $var_width );
            $filter_width = $filter_len if ( $filter_len > $filter_width );
            $hgvs_width = $hgvs_len if ( $hgvs_len > $hgvs_width );
        }
    }

    ( $filter_width > 13 ) ? ($filter_width += 4) : ($filter_width = 17);
    return ( $ref_width + 4, $var_width + 4, $filter_width, $hgvs_width + 2);
}

sub batch_lookup {
    # Process a lookup file, and load up @query_list
    my $file = shift;
    my @query_list;

    open( my $fh, "<", $$file ) or die "Can't open the lookup file: $!";
    chomp( @query_list = <$fh> );
    close $fh;

    my $query_string = join( ' ', @query_list );
    return \$query_string;
}
