#!/usr/bin/perl
# Script to parse Ion Torrent VCF files and output columnar information.  Can also filter based
# on several different criteria.  Designed to work with TVC v4.0+ VCF files.
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
my $scriptname = basename($0);
my $version = "v7.4.111617";

print colored("*" x 75, 'bold yellow on_black'), "\n";
print colored("\tDEVELOPMENT VERSION ($version) OF VCF EXTRACTOR", 'bold yellow on_black'), "\n";
print colored("*" x 75, 'bold yellow on_black'), "\n\n"; 

my $description = <<"EOT";
Parse and filter an Ion Torrent VCF file.  By default, this program will output 
a simple table in the following format:

     CHROM:POS REF ALT Filter Filter_Reason VAF TotCov RefCov AltCov COSID

However, using the '-a' option, we can add Ion Reporter (IR) annotations to the 
output, assuming the data was run through IR.  

In addition to simple output, we can also filter the data based on teh following criteria:
    - Non-reference calls can be omitted with '-n' option.
    - Only calls with Hotspot IDs can be output with '-H' option.
    - Only calls with an OVAT annotation can be output with '-O' option..
    - NOCALL variants can removed with the '-N' option.
    - Calls matching a specific gene or genes can be acquired with the '-g' option.
    - Calls matching a specific Hotspot ID can be acquired with the '-c' option.
    
This program can also output variants that match a position query based on using
the following string: 
    chr#:position. 
If a position is not quite known or we are dealing with a 0 vs 1 based position 
rule, we can perform a fuzzy lookup by using the '-f' option and the number of 
right-most digits in the position that we are unsure of (e.g. -f1 => 1234*, 
-f2 => 123**).

We can also use batch files to lookup multiple positions or hotspots using 
the '-l' option.  

        vcfExtractor -l lookup_file <vcf_file>
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] [-f {1,2,3}] <input_vcf_file>

    Program Options
    -a, --annot     Add IR and Oncomine OVAT annotation information to output if available.
    -V, --Verbose   Print additional information during processing.
    -c, --cfdna     Data is from a cfDNA run, and some metrics and thresholds will be different.
    -o, --output    Send output to custom file.  Default is STDOUT.
    -v, --version   Version information
    -h, --help      Print this help information

    Filter and Output Options
    -p, --pos       Output only variants at this position.  Format is "chr<x>:######" 
    -i, --id        Look for variant with matching variant ID (COSMIC ID or other Hotspot ID)
    -g, --gene      Filter variant calls by gene id. Can input a single value or comma 
                    separated list of gene ids to query. Can only be used with 
                    the '--annot' option as the annotations have to come from IR.
    -l, --lookup    Read a list of variants from a file to query the VCF. 
    -f, --fuzzy     Less precise (fuzzy) position match. Strip off n digits from the position string.
                    MUST be used with a query option (e.g. -p, -c, -l), and can not trim more than 3 
                    digits from string.
    -n, --noref     Output reference calls.  Ref calls filtered out by default
    -N, --NOCALL    Remove 'NOCALL' entries from output
    -O, --OVAT      Only report Oncomine Annotated Variants.
    -H, --HS        Print out only variants that have a Hotspot ID (NOT YET IMPLEMENTED).
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
my $gene;
my $hotspots;
my $cfdna;
my $debug_pos;  # undocumented.

GetOptions( "output|o=s"    => \$outfile,
            "cfdna|c"       => \$cfdna,
            "annot|a"       => \$annots, 
            "OVAT|O"        => \$ovat_filter,
            "id|i=s"        => \$hsids,
            "NOCALL|N"      => \$nocall,
            "pos|p=s"       => \$positions,
            "lookup|l=s"    => \$lookup,
            "fuzzy|f=i"     => \$fuzzy,
            "noref|n"       => \$noref,
            "gene|g=s"      => \$gene,
            "debug_pos=s"   => \$debug_pos, #undocumented.
            "version|v"     => \$ver_info,
            "Verbose|V"     => \$verbose,
            "HS|H"          => \$hotspots,
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
my @warnings;

# Check for vcftools; we can't run without it...for now.
if ( ! qx(which vcftools) ) {
    print "$err Required package 'vcftools' is not installed on this system. ",
        "Install vcftools ('vcftools.sourceforge.net') and try again.\n";
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
if ( $ovat_filter ) {
    if ($cfdna) {
        die "$err Can no combine the cfDNA option with the OVAT option at this time ",
            "as it would seem that OVAT is not yet supported by \nthe pipeline. This ",
            "might be added as a feature later on.\n";
    }
    elsif ( ! $annots ) {
    print "$info Requested Oncomine annotation filter without adding the OVAT ",
        "annotations. Auto adding the OVAT annotations!\n" if $verbose;
    $annots=1;
    }
}

# Implementing a debug position method to help with development.  This is an
# undocumented method that will allow for one to input a position and output
# the parsed hash of data and rest of method on only that position alone.
if ($debug_pos) {
    print '-'x25 . '  DEBUG  ' . '-'x25, "\n";
    print "Outputting data for position $debug_pos only.\n";
    print '-'x59, "\n";
}

# Make sure enough args passed to script
if ( scalar( @ARGV ) < 1 ) {
    print "$err No VCF file passed to script!\n\n";
    print $usage;
    exit 1;
}

# Parse the lookup file and add variants to the postions list if processing batch-wise
# TODO: this is going to break with cfDNA panel since TF decided to use AA changes
# as the ID.
my @valid_hs_ids = qw( BT COSM OM OMINDEL MCH PM_COSM PM_B PM_D PM_MCH PM_E CV );
if ($lookup) {
    my $query_list = batch_lookup(\$lookup);
    if ( grep { $$query_list =~ /$_\d+/ } @valid_hs_ids ) {
        $hsids = $$query_list;
    } 
    elsif ( $$query_list =~ /chr[0-9YX]+/ ) {
        $positions = $$query_list;
    }
    else {
        print "$err Issue with lookup file. Check and retry\n"; 
        exit 1;
    }
}

# Double check the query (position / cosid) format is correct
my (@coords, @cosids,@filter_list);
if ( $positions ) {
    print "Checking position id lookup...\n" if $verbose;
    @coords = split( /\s+/, $positions );
    for my $coord ( @coords ) {
        if ( $coord !~ /\Achr[0-9YX]+:\d+$/i ) {
            print "$err '$coord' not valid. Please use the following format for ",
                "position queries: 'chr#:position'\n";
            exit 1;
        }
    }
    push(@filter_list, 'position');
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
    push(@filter_list, 'hsid');
}

# Add gene list for query if we have one.
if ($gene) {
    die "$err Can not use the '--gene' option without the '--annot' option!\n" unless $annots;
    push(@filter_list, 'gene');
}
my @query_genes = map{ uc $_ } split(/,/,$gene) if $gene;

# Setting up hash of filters.  Right now, only accept one type as combinations 
# are probably redundant. We might find an  excuse to this later, so, so keep 
# the data struct.  Pass the array to the filter function later.
my %vcf_filters = (
    'gene'     => \@query_genes,
    'hsid'     => \@cosids,
    'position' => \@coords,
);
print "Applied filters: \n" and dd \%vcf_filters if DEBUG or $verbose;

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
	open( $out_fh, ">", $outfile ) 
    || die "Can't open the output file '$outfile' for writing: $!";
} else {
	$out_fh = \*STDOUT;
}

#########------------------------------ END Arg Parsing and validation ---------------------------------#########
my $inputVCF = shift;

# Check VCF file and options to make sure they're valid
open ( my $vcf_fh, "<", $inputVCF );
my @header = grep { /^#/ } <$vcf_fh>;
die "$err '$inputVCF' does not appear to be a valid VCF file or does not have a 
    header.\n" unless @header;
close $vcf_fh;

# Crude check for TVC3.2 or TVC4.0+ VCF file.  Still need to refine this
if ( grep { /^##INFO.*Bayesian_Score/ } @header ) {
    print "$warn '$inputVCF' appears to be from TVCv3.2, and the 'tvc32' option 
        was not selected.  The file may not be processed correctly.\n";
    die "Pre TVCv4.0 VCF file detected. These files are not longer supported 
        by this utility\n";
}

# Trigger IR / OVAT annot capture if available
my ($ir_annot,$ovat_annot);
( grep { /OncomineVariantAnnotation/ } @header ) ? ($ovat_annot = 1) : ($ovat_annot = 0);
( grep { /IonReporterExportVersion/ } @header ) ? ($ir_annot = 1) : ($ir_annot = 0);

if ( $annots && $ir_annot == 0 ) {
    die "$err IR output selected, but VCF does not appear to have been run through IR!\n";
}

# Figure out if these are cfDNA files to help with downstream options.
if (grep { /^##INFO=<ID=MMDP/ } @header) {
        if (! $cfdna) {
            print "$err You appear to be running a VCF derived from the taqseq ",
                "cfDNA panel, but have not set the --cfdna option.\n";
            exit 1
        }
        elsif ($ovat_filter) {
            print "$err OVAT filtered output selected, but VCF appears to have been run ",
            "through the TagSeq cfDNA pipeline, and and there is no OVAT annotation\n",
            "available at this time!\n";
            exit 1;
        }
}

# Get the data from VCF Tools
my @wanted_fields = qw(%CHROM:%POS %REF %ALT %FILTER %INFO/FR %INFO/OID %INFO/OPOS 
    %INFO/OREF %INFO/OALT %INFO/OMAPALT --- [%GTR %AF %FRO %RO %FAO %AO %DP]);

if ($annots) {
    $wanted_fields[10] = '%INFO/FUNC';
}

if ($cfdna) {
    $wanted_fields[12] = '%MAF';
    $wanted_fields[13] = '%MRO';
    $wanted_fields[15] = '%MAO';
    $wanted_fields[17] = '%MDP]';
}

my $vcf_format = join('\t', @wanted_fields);
my @extracted_data = qx( vcf-query $inputVCF -f "$vcf_format\n" );

# Read in the VCF file data and create a hash
my %vcf_data = parse_data( \@extracted_data );

# Filter parsed data.
my $filtered_vcf_data = filter_data(\%vcf_data, \%vcf_filters);

# if fuzzy results, make them look like regular
$filtered_vcf_data = flatten_fuzzy_results($filtered_vcf_data) if $fuzzy; 

# Finally print it all out.
format_output($filtered_vcf_data, \%vcf_filters);

# Wrap up
if ( @warnings && $verbose ) {
    print "\n";
    print $_ for @warnings;
}

sub parse_data {
    # Extract the VCF information and create a hash of the data.  
    # TODO: Want to add amplicon coverage data and LOD% data for cfDNA panel. Can
    # add values to end of @$data array which will be undef if not available?
    my $data = shift;
    my %parsed_data;

    for ( @$data ) {
        my ( $pos, $ref, $alt, $filter, $reason, $oid, $opos, $oref, $oalt, 
            $omapalt, $func, $gtr, $af, $fro, $ro, $fao, $ao, $dp ) = split( /\t/ );

        # Limit processing to just one position and output more metrics so that we
        # can figure out what's going on.
        if ($debug_pos) {
            next unless $pos eq $debug_pos;
            print_debug_output([split(/\t/)]);
        }
        #__exit__('319','Post debugger output');
        
        # IR generates CNV and Fusion entries that are not compatible.  
        # TODO: for cfDNA can we integrate these calls now?  
        next if ( $alt =~ /[.><\]\d+]/ ); 

        # Clean up filter reason string
        $reason =~ s/^\.,//;

        # Filter out vars we don't want to print out later anyway.
        next if $reason eq "NODATA";
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

            # Get the normalizedRef, normalizedAlt, and normalizedPos values from 
            # the REF and ALT fields so that we can map the FUNC block.
            my @coords = split(/:/, $pos);
            my %norm_data = normalize_variant(\$ref, \$alt_var, $coords[1]);

            my @array_pos = grep { $omapalt_array[$_] eq $alt_var } 0..$#omapalt_array;
            for my $index ( @array_pos ) {
                (my $parsed_pos = $pos) =~ s/(chr\d+:).*/$1$norm_data{'normalizedPos'}/; 
                
                my $var_id = join( ":", $parsed_pos, $oref_array[$index], $oalt_array[$index] );
                my $cosid = $oid_array[$index];
                # Stupid bug with de novo and hotspot merge that can create two 
                # duplicate entries for the same variant but one with and one 
                # without a HS (also different VAF, coverage,etc). Try this to 
                # capture only HS entry.
                if ( $cosid ne '.' && exists $parsed_data{$var_id} ) {
                   delete $parsed_data{$var_id}; 
                }

                # Grab the OVAT annotation information from the FUNC block if possible.
                my ($ovat_gc, $ovat_vc, $gene_name, $transcript, $hgvs, $protein, $function, $exon);
                if ( $func eq '.' ) {
                    push( @warnings, "$warn could not find FUNC entry for '$pos'\n") if $annots;
                    $ovat_vc = $ovat_gc = $gene_name = "NULL";
                } else {
                    ($ovat_gc, $ovat_vc, $gene_name, $transcript, $hgvs, $protein, $function, 
                        $exon) = get_ovat_annot(\$func, \%norm_data) unless $func eq '---'; 
                }

                # Start the var string.
                push( @{$parsed_data{$var_id}},
                    $parsed_pos,
                    $norm_data{'normalizedRef'},
                    $norm_data{'normalizedAlt'},
                    $filter, 
                    $reason, 
                    $gtr );

                # Check to see if call is result of long indel assembler and handle appropriately. 
                my ($vaf, $tot_coverage);
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
                
                # Filter out reference calls if we have turned on the noref filter. Have to leave the NOCALL 
                # calls if we have left those in, and have to deal with sub 1% VAFs for cf DNA assay.
                my $calc_vaf = ${$parsed_data{$var_id}}[6];
                if ( $calc_vaf ne '.' ) {
                    if ($calc_vaf == 0) {
                        delete $parsed_data{$var_id} and next if $noref;
                    }
                    elsif ($calc_vaf < 1) {
                        delete $parsed_data{$var_id} and next if ! $cfdna; 
                    }
                }

                # Now handle in two steps.  Add IR annots if there, and then if wanted ovat annots, add them too.
                push(@{$parsed_data{$var_id}}, $gene_name, $transcript, $hgvs, $protein, $exon, $function) if $annots;
                push(@{$parsed_data{$var_id}}, $ovat_gc, $ovat_vc) if $annots and $ovat_annot;
            }
        }
    }
    #dd \%parsed_data;
    #exit;
    return %parsed_data;
}

sub get_ovat_annot {
    # If this is IR VCF, add in the OVAT annotation. 
    #no warnings;
    my ($func, $norm_data) = @_;
    my %data;
    my @wanted_elems = qw(oncomineGeneClass oncomineVariantClass gene transcript 
        protein coding function normalizedRef normalizedAlt location exon);

    $$func =~ tr/'/"/;
    my $json_annot = JSON::XS->new->decode($$func);

    my $match;
    for my $func_block ( @$json_annot ) {
        # if there is a normalizedRef entry, then let's map the func block 
        # appropriately....as long as we're not running the cfDNA panel, which
        # has new "decoy" alleles that I think screw this all up.  
        if (! $cfdna && $$func_block{'normalizedRef'}) {
            if ($$func_block{'normalizedRef'} eq $$norm_data{'normalizedRef'} && $$func_block{'normalizedAlt'} eq $$norm_data{'normalizedAlt'}) {
                %data = %$func_block;
                $match = 1;
                last;
            } 
        } 
        else {
            @data{@wanted_elems} = @{$func_block}{@wanted_elems};
        }
    }

    my $gene_class    = $data{'oncomineGeneClass'}    // '---';
    my $variant_class = $data{'oncomineVariantClass'} // '---';
    my $gene_name     = $data{'gene'}                 // '---';
    my $protein       = $data{'protein'}              // '---';
    my $hgvs          = $data{'coding'}               // '---';
    my $transcript    = $data{'transcript'}           // '---';
    my $function      = $data{'function'}             // '---';
    my $ref           = $data{'normalizedRef'}        // '---';
    my $alt           = $data{'normalizedAlt'}        // '---';
    my $location;
    if ($data{'location'} eq 'exonic') {
        $location = "Exon$data{'exon'}";
    } else {
        $location = $data{'location'};
    }
    $location //= '---';

    # Sometimes, for reasons I'm not quite sure of, there can be an array for the functional annotation.  I think it's 
    # safe to take the most severe of the list and to use.  
    ($function) = grep {/(missense|nonsense)/} @$function if ref $function eq 'ARRAY';

    if (DEBUG) {
        print "======================  DEBUG  =======================\n\n";
        print "gc       => $gene_class\n";
        print "vc       => $variant_class\n";
        print "gene     => $gene_name\n";
        print "ref      => $ref\n";
        print "alt      => $alt\n";
        print "AA       => $protein\n";
        print "tscript  => $transcript\n";
        print "HGVS     => $hgvs\n";
        print "function => $function\n";
        print "location => $location\n";
        print "======================================================\n\n";
    }
    return ($gene_class, $variant_class, $gene_name, $transcript, $hgvs, $protein, $function, $location );
}

sub normalize_variant {
    # Borrowed from ThermoFisher's vcf.py script to convert IR VCFs. Trim from both ends until only unique
    # sequence left
    my ($ref,$alt,$pos) = @_;
    my ($norm_ref, $norm_alt);

    my ($rev_ref, $rev_alt, $position_delta) = rev_and_trim($ref, $alt);
    ($norm_ref, $norm_alt, $position_delta) = rev_and_trim(\$rev_ref, \$rev_alt);

    my $adj_position = $position_delta + $pos;
    return ( 'normalizedRef' => $norm_ref, 'normalizedAlt' => $norm_alt, 
        'normalizedPos' => $adj_position );
}

sub rev_and_trim {
    # Borrowed from ThermoFisher's vcf.py script to convert IR VCFs
    my ($ref, $alt) = @_;
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
    my ($filter, $tcov, $rcov, $acov) = @_;
    my $vaf;

    if ( $$filter eq "NOCALL" ) { 
        $vaf = '.';
    }
    elsif( $$filter eq "NODATA" || $$tcov == 0) {
        $vaf = 0;
    }
    else {
        if ($cfdna) { 
            $vaf = sprintf( "%.4f", 100*($$acov / $$tcov) );
        } else {
            $vaf = sprintf( "%.2f", 100*($$acov / $$tcov) );
        }
    }
    return $vaf;
}

sub filter_data {
    # Filter extracted VCF data and return a hash of filtered data.
    my ($data, $filter) = @_;
    my %filtered_data;
    my @fuzzy_pos;
    my %counter;

    my $on  = colored( "On", 'bold green on_black');
    my $off = colored( "Off", 'bold red on_black');

    # First run OVAT filter; no need to push big list of variants through other filters.
    if ($verbose) {
        print "$info OVAT filter status: ";
        ($ovat_filter) ? print "$on!\n" : print "$off.\n";
        print "$info Hotspot ID filter status: ";
        ($hotspots) ? print "$on!\n" : print "$off.\n";
        print "$info NOCALLs output to results: ";
        ($nocall) ? print "$off!\n" : print "$on.\n";
        print "$info Reference calls output to results: ";
        ($noref) ? print "$off!\n" : print "$on.\n";
    }
    $data = ovat_filter($data) if $ovat_filter;
    $data = hs_filtered($data) if $hotspots;

    # Determine filter to run, and if none, just return the full set of data.
    my @selected_filters = grep { scalar @{$$filter{$_}} > 0 } keys %$filter;
    if (@selected_filters > 1) {
        print "ERROR: Using more than one type of filter is redundant and not 
            accepted at this time. (Filters chosen: ";
        print join(',', @selected_filters), " )\n";
        exit 1;
    }
    return $data unless @selected_filters;

    # Now run full filter tree.
    if ( $fuzzy ) {
        my $re = qr/(.*).{$fuzzy}/;
        @fuzzy_pos = map { /$re/ } @{$$filter{position}};
        for my $query (@fuzzy_pos) {
            for (sort keys %$data) {
                push(@{$filtered_data{$query}}, [@{$$data{$_}}]) if ($$data{$_}[0] =~ /$query.{$fuzzy}$/);
            }
        }
    } 
    elsif ($selected_filters[0] eq 'gene') {
        print "$info Running the gene filter...\n" if $verbose;
        for my $variant (keys %$data) {
            if ( grep { $$data{$variant}[11] eq $_ } @{$$filter{$selected_filters[0]}}) {
                @{$filtered_data{$variant}} = @{$$data{$variant}};
            }
        }
    }
    elsif ($selected_filters[0] eq 'hsid') {
        print "$info Running the HSID filter...\n" if $verbose;
        for my $variant (keys %$data) {
            if ( grep { $$data{$variant}[10] eq $_ } @{$$filter{$selected_filters[0]}}) {
                @{$filtered_data{$variant}} = @{$$data{$variant}};
            }
        }
    }
    elsif ($selected_filters[0] eq 'position') {
        print "$info Running the Position filter...\n" if $verbose;
        for my $variant (keys %$data) {
            if (grep { $$data{$variant}[0] eq $_ } @{$$filter{$selected_filters[0]}}) {
                @{$filtered_data{$variant}} = @{$$data{$variant}};
            }
        }
    }
    return \%filtered_data;
}

sub ovat_filter {
    # Filter out calls that are not oncomine reportable variants
    my $data = shift;
    print "$info Running ovat filter\n" if $verbose;
    for my $variant ( keys %$data ) {
        delete $$data{$variant} if $$data{$variant}->[18] eq '---';
    }
    return $data;
}

sub hs_filtered {
    # Filter out calls that are not oncomine reportable variants
    my $data = shift;
    print "$info Running Hotspots filter\n" if $verbose;
    for my $variant ( keys %$data ) {
        delete $$data{$variant} if $$data{$variant}->[10] eq '.';
    }
    return $data;
}

sub flatten_fuzzy_results {
    # Make the fuzzy results look more like non-fuzzy match results to make reporting easier.
    my $data = shift;
    my %results;
    for my $fuzzy_position (keys %$data) {
        my $counter = 0;
        for my $var (@{$$data{$fuzzy_position}}) {
            $counter++;
            $results{"$fuzzy_position:$counter"} = $var;
        }
    }
    return \%results;
}

sub format_output {
    # Format and print out the results
    # w1 => REF(index:1), w2 => ALT(index:2), w3 => Filter comment(index:4), w4 => CDS(index:13), w5 => AA(index:14)
    # TODO: Fix the header for cfDNA to output the molecular tag headers since we're using that data.
    my ($data,$filter_list) = @_;
    my $ref_width = 8;
    my $alt_width = 8;
    my $filter_width = 17;
    my $cds_width = 7;
    my $aa_width = 7;

    select $out_fh;

    # Set up the output header and the correct format string to use.
    my ($format, @header);
    if ($nocall) {
        ($ref_width, $alt_width) = field_width($data,[1,2]) if %$data;
        $format = "%-17s %-${ref_width}s %-${alt_width}s %-8s %-8s %-8s %-8s %-12s\n";
        @header = qw( CHROM:POS REF ALT VAF TotCov RefCov AltCov COSID );
    } else {
        ($ref_width,$alt_width,$filter_width) = field_width($data,[1,2,4]) if %$data;
        $filter_width = 17 if $filter_width < 17;
        $format = "%-17s %-${ref_width}s %-${alt_width}s %-8s %-${filter_width}s %-8s %-8s %-8s %-10s %-12s\n";
        @header = qw( CHROM:POS REF ALT Filter Filter_Reason VAF TotCov RefCov AltCov COSID );
    }
    if ($annots) {
        ($cds_width,$aa_width) = field_width($data,[13,14]) if %$data;
        push(@header, qw(Gene Transcript CDS AA Location Function));
        $format =~ s/\n$/ %-14s %-15s %-${cds_width}s %-${aa_width}s %-12s %-22s\n/;
        if ($ovat_annot) {
            push(@header,qw(oncomineGeneClass oncomineVariantClass));
            $format =~ s/\n$/ %-21s %-21s \n/;
        }
    }
    printf $format, @header;

    # Handle null result reporting depending on the filter used.
    if (! %$data) {
        print "\n>>> No Oncomine Annotated Variants Found! <<<\n" and exit if $ovat_filter;
        print "\n>>> No Variants Found for Gene(s): " . join(', ', @{$$filter_list{gene}}), "! <<<\n" and exit if @{$$filter_list{gene}};
        print "\n>>> No Variants Found for Hotspot ID(s): " . join(', ', @{$$filter_list{hsid}}), "! <<<\n" and exit if @{$$filter_list{hsid}};
        if (@{$$filter_list{position}}) {
            my @positions;
            for my $query (@{$$filter_list{position}}) {
                ($fuzzy) ? push(@positions, (substr($query, 0, -$fuzzy) . '*'x$fuzzy)) : push(@positions, $query);
            }
            print "\n>>> No variant found at position(s): ", join(', ', @positions), "! <<<\n" and exit;
        } 
    } else {
        # TODO: Ouput formatting is not rigtht for cfDNA data.  I think there is some missing formatting.
        my @output_data;
        for my $variant ( sort { versioncmp( $a, $b ) } keys %$data ) {
            ($nocall) ? (@output_data = @{$$data{$variant}}[0,1,2,6..18]) : (@output_data = @{$$data{$variant}}[0..4,6..18]);
            @output_data[9..13] = map { $_ //= 'NULL' } @output_data[9..13];
            printf $format, @output_data;
        }
    }
}

sub field_width {
    # Load in a hash of data and an array of indices for which we want field width info, and
    # output an array of field widths to use in the format string.
    my ($data,$indices) = @_;
    my @return_widths;
    for my $pos (@$indices) {
        my @elems = map { ${$$data{$_}}[$pos] } keys %$data;
        push(@return_widths, get_longest(\@elems)+2);
    }
    return @return_widths;
}

sub get_longest {
    my $array = shift;
    my @lens = map { length($_) } @$array;
    my @sorted_lens = sort { versioncmp( $b, $a) } @lens;
    return $sorted_lens[0];
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

sub print_debug_output {
    # DEBUG: Can add position to filter and output a hash of data to help.
    my $data = shift;
    my @fields = qw(pos ref alt filter reason oid opos oref oalt omapalt func 
        gtr af fro ro fao ao dp);
    my %foo;

    @foo{@fields} = map{chomp; $_} @$data;
    print '='x25, "  DEBUG  ", "="x25, "\n";
    dd \%foo;
    print '='x59, "\n";
}

sub __exit__ {
    my ($line, $msg) = @_;
    print "\n\n";
    print colored("Got exit message at line $line with message: $msg", 'bold white on_green');
    print "\n";
    exit;
}
