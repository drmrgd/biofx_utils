#!/usr/bin/perl
# Script to parse Ion Torrent VCF files and output columnar information.  Can 
# also filter based on several different criteria.  Designed to work with TVC 
# v4.0+ VCF files.
# 
# D Sims 2/21/2013
################################################################################
use warnings;
use strict;
use version;
use autodie;
use feature 'state';

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use List::Util qw( sum min max );
use Cwd qw(abs_path);
use Sort::Versions;
use JSON::XS; 
use Data::Dump;
use File::Basename;
use Term::ANSIColor;
use Time::Piece;

use constant 'DEBUG' => 0;
my $scriptname = basename($0);
my $version = "v7.28.032218";

print colored("*" x 75, 'bold yellow on_black'), "\n";
print colored("\tDEVELOPMENT VERSION ($version) OF VCF EXTRACTOR", 
    'bold yellow on_black'), "\n";
print colored("*" x 75, 'bold yellow on_black'), "\n\n"; 

my $description = <<"EOT";
Parse an Ion Torrent or Ion Reporter VCF file. By default, this program will
output a simple table in the following format:

     CHROM:POS REF ALT Filter Filter_Reason VAF TotCov RefCov AltCov COSID

However, using the '-a' option, we can add Ion Reporter (IR) annotations, 
including OVAT annotation, to the output, assuming the data was run through IR.

In addition to simple output, we can also filter the data based on the 
following criteria:
    - Non-reference calls can be omitted with '-n' option.
    - NOCALL variants can removed with the '-N' option.
    - Only calls with Hotspot IDs can be output with '-i' option.
    - Only calls with an OVAT annotation can be output with '-O' option..
    - Calls matching a specific gene or genes can be acquired with the '-g' 
      option.
    - Calls matching a specific Hotspot ID can be acquired with the '-c' 
      option.
    
This program has also been updated to include parsing of TaqSeq based cfDNA 
assay results, and in addition to the LOD metric, outputs amplicon coverage as
TotalCov, with the rest of the coverage being derived from Molecular Family 
data.

This program can also output variants that match a position query based on 
using the following string: 
    chr#:position. 
If a position is not quite known or we are dealing with a 0 vs 1 based 
position rule, we can perform a fuzzy lookup by using the '-f' option and the 
number of right-most digits in the position that we are unsure of (e.g. 
-f1 => 1234*, -f2 => 123**).

We can also use batch files to lookup multiple positions or hotspots using 
the '-l' option.  

        vcfExtractor -l lookup_file <vcf_file>

In addition to the Perl modules, Sort::Versions, JSON::XS, and Term::ANSIcolor,
which may not be standard in your distribution, the program will require the
vcftools ('vcftools.sourceforge.net') package to be installed and in your 
\$PATH.

Finally, output can be either in the form of an easy to read table (default), or
as a new VCF file, containing the data that was left after filtering. This can be
useful for cases where one still needs a VCF format, but does not want all of the 
extra stuff that is associated with a TVC VCF file.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <input_vcf_file>

    Program Options
    -a, --annot     Add IR and Oncomine OVAT annotation information to output 
                    if available.
    -q, --quiet     Print additional information during processing.
    -c, --cfdna     Data is from a cfDNA run, and some metrics and thresholds 
                    will be different.
    -o, --output    Send output to custom file.  Default is STDOUT.
    -v, --version   Version information
    -h, --help      Print this help information

    Filter and Output Options
    -p, --pos       Output only variants at this position.  Format is 
                    "chr<x>:######" 
    -i, --id        Look for variant with matching variant ID (COSMIC ID or 
                    other Hotspot ID)
    -g, --gene      Filter variant calls by gene id. Can input a single value 
                    or comma 
                    separated list of gene ids to query. Can only be used with 
                    the '--annot' option as the annotations have to come from 
                    IR.
    -l, --lookup    Read a list of search terms from a file to query the VCF.
                    Note that you will not have to use the '-p', '-g', or '-i'
                    option (for example) when using this option as the type 
                    will automatically be detected from the file. Also note 
                    that just like the other options, mixed query types is not
                    supported and will result in an error.
    -f, --fuzzy     Less precise (fuzzy) position match. Strip off n digits 
                    from the position string. MUST be used with a query option
                    (e.g. -p, -c, -l), and can not trim more than 3 digits from
                    string.
    -n, --noref     Output reference calls.  Ref calls filtered out by default
    -N, --NOCALL    Remove 'NOCALL' entries from output
    -O, --OVAT      Only report Oncomine Annotated Variants.
    -H, --HS        Print out only variants that have a Hotspot ID (NOT YET 
                    IMPLEMENTED).
    -V, --VCF       Output data as a new, simple VCF file, rather than a table.
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
my $quiet;
my $gene;
my $hotspots;
my $cfdna;
my $debug_pos;  # undocumented.
my $vcf_output;

GetOptions( 
    "output|o=s"    => \$outfile,
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
    "quiet|q"       => \$quiet,
    "HS|H"          => \$hotspots,
    "VCF|V"         => \$vcf_output,
    "help|h"        => \$help )
or die "\n$usage";

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version_info {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help if $help;
version_info if $ver_info;

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

# Double check that fuzzy option is combined intelligently with a position 
# lookup.
if ( $fuzzy) {
    if ( $fuzzy > 3 ) {
        print "\n$err Can not trim more than 3 digits from the query string.\n";
        exit 1;
    }
    elsif ( $lookup ) {
        print "\n$warn fuzzy lookup in batch mode may produce a lot of ",
            "results! Continue? ";
        chomp( my $response = <STDIN> );
        exit if ( $response =~ /[(n|no)]/i );
        print "\n";
    }
    elsif ( ! $positions && ! $hsids ) {
        print "$err must include position or hotspot ID query with the '-f' ",
            "option\n\n";
        print $usage;
        exit 1;
    }
}

# Throw a warning if using the ovat filter without asking for OVAT annotations.
if ( $ovat_filter ) {
    if ($cfdna) {
        die "$err Can no combine the cfDNA option with the OVAT option at ",
            "this time as it would seem that OVAT is not yet supported by\n",
            "the pipeline. This might be added as a feature later on.\n";
    }
    elsif ( ! $annots ) {
    print "$info Requested Oncomine annotation filter without adding the OVAT ",
        "annotations. Auto adding the OVAT annotations!\n" if $quiet;
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

# Make sure VCF has been passed to script.
if ( scalar( @ARGV ) < 1 ) {
    print "$err No VCF file passed to script!\n\n";
    print $usage;
    exit 1;
}

# Parse the lookup file and add variants to the postions list if processing 
# batch-wise
my @valid_hs_ids = qw( BT COSM OM OMINDEL MCH PM_COSM PM_B PM_D PM_MCH PM_E 
    CV MAN );
if ($lookup) {
    if ($positions or $hsids) {
        print "$err You can not use individual filters with the lookup file ",
            "option.\n";
        exit 1;
    }

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
my (@coords, @cosids, @filter_list);
if ( $positions ) {
    print "Checking position id lookup...\n" if $quiet;
    @coords = split( /\s+/, $positions );
    for my $coord ( @coords ) {
        if ( $coord !~ /\Achr[0-9YX]+:\d+$/i ) {
            print "$err '$coord' not valid. Please use the following format ",
                "for position queries: 'chr#:position'\n";
            exit 1;
        }
    }
    push(@filter_list, 'position');
} 
elsif ( $hsids ) {
    if ($cfdna) {
        print "You can not use HS lookup with the cfDNA option at this time ",
            "since the variant ID is an amino acid sequence and difficult\n",
            "to work with.\n";
        exit 1;
    }
    print "Checking variant id lookup...\n" if $quiet;
    @cosids = split( /\s+/, $hsids );
    for my $varid ( @cosids ) { 
        if ( ! grep { $varid =~ /$_\d+/ } @valid_hs_ids ) {
            print "$err '$varid' is not a valid hotspot query term! Valid ",
                "lookups are:\n";
            print "\t${_}###\n" for @valid_hs_ids;
            exit;
        }
    }
    push(@filter_list, 'hsid');
}

# Add gene list for query if we have one.
if ($gene) {
    die "$err Can not use the '--gene' option without the '--annot' ",
        "option!\n" unless $annots;
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
print "Applied filters: \n" and dd \%vcf_filters if DEBUG or $quiet;

# Set up output type
my $output_type;
($vcf_output) ? ($output_type = 'vcf') : ($output_type = 'table');

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
	open( $out_fh, ">", $outfile ) 
        || die "Can't open the output file '$outfile' for writing: $!";
    print "Writing output to '$outfile.'\n";
} else {
	$out_fh = \*STDOUT;
}

#########--------------- END Arg Parsing and validation ---------------#########
my $input_vcf = shift;

# Check VCF file and options to make sure they're valid
open ( my $vcf_fh, "<", $input_vcf);
my @header = grep { /^#/ } <$vcf_fh>;
die "$err '$input_vcf' does not appear to be a valid VCF file or does not ",
    "have a header.\n" unless @header;
close $vcf_fh;

# Check for TVC4.4+ VCF file. We can no longer process older files.
my ($vcf_version) = map { /^##source="tvc (.*?) \(.*$/ } @header;
$vcf_version =~ s/-/./; 
if (version->parse($vcf_version) < version->parse('4.4.2')) {
    die "$err Pre TVCv4.4 VCF file detected.  These are no longer supported by",
        " this utility and must be manually parsed.\n";
}

# Trigger IR / OVAT annot capture if available
my ($ir_annot,$ovat_annot);
( grep { /OncomineVariantAnnotation/ } @header ) 
    ? ($ovat_annot = 1) 
    : ($ovat_annot = 0);

( grep { /IonReporterExportVersion/ } @header ) 
    ? ($ir_annot = 1) 
    : ($ir_annot = 0);

if ( $annots && $ir_annot == 0 ) {
    die "$err IR output selected, but VCF does not appear to have been run ",
    "through IR!\n";
}

# Figure out if these are cfDNA files to help with downstream options.
if (grep { /^##INFO=<ID=LOD/ } @header) {
        if (! $cfdna) {
            print "$err You appear to be running a VCF derived from a taqseq ",
                "cfDNA panel, but have not set the --cfdna option.\n";
            exit 1
        }
        elsif ($ovat_filter) {
            print "$err OVAT filtered output selected, but VCF appears to have",
            " been run through the TagSeq cfDNA pipeline, and and there is no ",
            "OVAT annotation\navailable at this time!\n";
            exit 1;
        }
} else {
    if ($cfdna) {
        print "$err cfDNA option (--cfdna) selected, but VCF does not appear ",
        "to be from a cfDNA run!\n";
        exit 1; 
    }
}

# Get the data from VCF Tools
my @wanted_fields = qw(%CHROM:%POS %REF %ALT %FILTER %INFO/FR %INFO/OID 
    %INFO/OPOS %INFO/OREF %INFO/OALT %INFO/OMAPALT --- --- [%GTR %AF %FRO %RO 
    %FAO %AO %DP]);

if ($annots) {
    $wanted_fields[10] = '%INFO/FUNC';
}

if ($cfdna) {
    $wanted_fields[11] = '%INFO/LOD';
    $wanted_fields[13] = '%MAF';
    $wanted_fields[14] = '%MRO';
    $wanted_fields[16] = '%MAO';
}

my $vcf_format = join('\t', @wanted_fields);
my @extracted_data = qx( vcf-query $input_vcf -f "$vcf_format\n" );

# Read in the VCF file data and create a hash
my %vcf_data = parse_data( \@extracted_data );

# Filter parsed data.
my $filtered_vcf_data = filter_data(\%vcf_data, \%vcf_filters);

# Finally print it all out.
# XXX
format_output($filtered_vcf_data, \%vcf_filters, \@header, $output_type);

# Wrap up
if ( @warnings && $quiet) {
    print "\n";
    print $_ for @warnings;
}

sub parse_data {
    # Extract the VCF information and create a hash of the data.  
    my $data = shift;
    my %parsed_data;

    for ( @$data ) {
        chomp;
        my ( $pos, $ref, $alt, $filter, $reason, $oid, $opos, $oref, $oalt, 
            $omapalt, $func, $lod, $gtr, $af, $fro, $ro, $fao, $ao, 
            $dp ) = split( /\t/ );

        # Limit processing to just one position and output more metrics so
        # that we can figure out what's going on.
        if ($debug_pos) {
            next unless $pos eq $debug_pos;
            print_debug_output([split(/\t/)]);
        }

        # IR generates CNV and Fusion entries that are not compatible.  
        next if ( $alt =~ /[.><\]\d+]/ ); 

        # Clean up filter reason string
        $reason =~ s/^\.,//;

        # Filter out vars we don't want to print out later anyway.
        next if $reason eq "NODATA";
        
        # Sometimes not getting correct 'NOCALL' flag; manually set if need to.
        $filter = "NOCALL" if (($gtr =~ m|\./\.|) or ($reason eq 'REJECTION'));
        next if ( $nocall && $filter eq "NOCALL" );
        next if ( $noref && $gtr eq '0/0' );

        # Create some arrays to hold the variant data in case we have MNP calls
        my @alt_array     = split( /,/, $alt );
        my @oid_array     = split( /,/, $oid );
        my @opos_array    = split( /,/, $opos );
        my @oref_array    = split( /,/, $oref );
        my @oalt_array    = split( /,/, $oalt );
        my @omapalt_array = split( /,/, $omapalt );
        my @fao_array     = split( /,/, $fao );
        my @ao_array      = split( /,/, $ao );
        my @lod_array     = split( /,/, $lod );

        for my $alt_index ( 0..$#alt_array ) {

            # NOTE: New vars for flagging.
            my $caller = 'tvc'; # set default as TVC and change later.

            my $alt_var = $alt_array[$alt_index];

            # Get the normalizedRef, normalizedAlt, and normalizedPos values
            # from the REF and ALT fields so that we can map the FUNC block.
            my @coords = split(/:/, $pos);
            my %norm_data = normalize_variant(\$ref, \$alt_var, $coords[1]);

            my @array_pos = grep { $omapalt_array[$_] eq $alt_var } 0..$#omapalt_array;
            for my $index ( @array_pos ) {
                my @var_data; # Temp holder array for collected variant data.
                (my $parsed_pos = $pos) =~ s/(chr\d+:).*/$1$norm_data{'normalizedPos'}/; 
                
                my $var_id = join( ":", $parsed_pos, $oref_array[$index], 
                    $oalt_array[$index] );
                my $cosid = $oid_array[$index];

                # Stupid bug with de novo and hotspot merge that can create two 
                # duplicate entries for the same variant but one with and one 
                # without a HS (also different VAF, coverage,etc). Try this to 
                # capture only HS entry.
                if ($parsed_data{$var_id}) {
                    # if data already there, check if it's from LIA and replace
                    if ($parsed_data{$var_id}->[-1] eq 'lia') {
                        delete $parsed_data{$var_id};
                    }
                    # Pedmatch has duplicate hotspots, both of which are 
                    # populating the output. If we already have a variant 
                    # entry and it's the PM version replace with COSMIC 
                    # version.
                    elsif ($parsed_data{$var_id}->[8] =~ /^PM/) {
                        delete $parsed_data{$var_id};
                    }

                    # If TVC Duplicate Hotspot bug, 
                    delete $parsed_data{$var_id} if ( $cosid eq '.');
                } 

                # Start the var string.
                push( @var_data,
                    $parsed_pos,
                    $norm_data{'normalizedRef'},
                    $norm_data{'normalizedAlt'},
                    $filter, 
                    $reason
                );

                # Check to see if call is result of long indel assembler and 
                # handle appropriately. 
                my ($vaf, $tot_coverage);
                if ( $fao_array[$alt_index] eq '.' ) {
                    # Result is from LIA
                    $caller = 'lia';
                    $tot_coverage = $ao_array[$alt_index] + $ro;
                    $vaf = vaf_calc( \$filter, \$dp, \$ro, 
                        \$ao_array[$alt_index] );
                    
                    push(@var_data, $vaf, $lod_array[$alt_index],
                        $tot_coverage, $ro, $ao_array[$alt_index], $cosid);
                } else {
                    my @cleaned_fao_array = grep { $_ ne '.' } @fao_array;
                    $tot_coverage = sum( @cleaned_fao_array ) + $fro;
                    $vaf = vaf_calc( \$filter, \$tot_coverage, \$fro, 
                        \$fao_array[$alt_index] );
                    
                    push(@var_data, $vaf, $lod_array[$alt_index],
                        $tot_coverage, $fro, $fao_array[$alt_index], $cosid );
                }
                
                # Filter out reference calls if we have turned on the noref 
                # filter. Have to leave the NOCALL calls if we have left those 
                # in, and have to deal with sub 1% VAFs for cfDNA assay.
                my $calc_vaf = $var_data[5];
                if ( $calc_vaf ne '.' ) {
                    if ($calc_vaf == 0) {
                        next if $noref;
                    }
                    elsif ($calc_vaf < 1) {
                        next if ! $cfdna; 
                    }
                }

                # Make some data changes for output if we have cfDNA and not 
                # conventional panel / assay.
                if ($cfdna) {
                    $var_data[7] = $dp;
                }
                else {
                    # We don't have cfDNA, so remove the LOD field from output.
                    splice(@var_data, 6, 1);
                }
                
                # Grab the OVAT annotation information from the FUNC block if 
                # possible.
                my ($ovat_gc, $ovat_vc, $gene_name, $transcript, $hgvs, 
                    $protein, $function, $exon);
                if ( $func eq '.' ) {
                    push( @warnings, "$warn could not find FUNC entry for '$pos'\n") if $annots;
                    $ovat_vc = $ovat_gc = $gene_name = "NULL";
                } 
                else {
                    ($ovat_gc, $ovat_vc, $gene_name, $transcript, $hgvs,
                        $protein, $function, $exon) = 
                        get_ovat_annot(\$func, \%norm_data) unless $func eq '---'; 
                }

                # Now handle in two steps.  Add IR annots if there, and then 
                # if wanted ovat annots, add them too.
                
                push(@var_data, $gene_name, $transcript, $hgvs, $protein, $exon,
                    $function) if $annots;
                push(@var_data, $ovat_vc) if $annots and $ovat_annot;
                push(@var_data, $caller);

                # Load data into hash.
                $parsed_data{$var_id} = \@var_data;
            }
        }
    }
    
    #dd \%parsed_data;
    #__exit__(__LINE__, "Finsished parsing data and generating variant hash.");
    return %parsed_data;
}

sub get_ovat_annot {
    # If this is IR VCF, add in the OVAT annotation. 
    my ($func, $norm_data) = @_;
    my %data;
    my @wanted_elems = qw(oncomineGeneClass oncomineVariantClass gene 
        transcript protein coding function normalizedRef normalizedAlt 
        location exon);

    $$func =~ tr/'/"/;
    my $json_annot = JSON::XS->new->decode($$func);
    
    for my $func_block ( @$json_annot ) {
        # In some cases, overlapping transcripts / genes can induce multiple 
        # FUNC block entries, which makes things a mess. Skip over any FUNC
        # block entry that does not have the correct transcript ID based on our
        # canonical transcript list.  Since some primary transcripts have 
        # changed over time, will get an array of ids to map.
        my $can_tran = get_can_tran($func_block->{'gene'});

        # No transcript ID in FUNC block on occassion for some reason  
        $func_block->{'transcript'} //= $can_tran->[0];
        next if $can_tran eq 'None' 
            or ! grep { $func_block->{'transcript'} eq $_ } @$can_tran;

        # If there is "normalized" data, then we got a positive variant call; 
        # map the appropriate elems.  If not, then likely ref call, and map 
        # what you can, fill in the rest later.
        if ($$func_block{'normalizedRef'}) {
            if ($$func_block{'normalizedRef'} eq $$norm_data{'normalizedRef'} 
                && $$func_block{'normalizedAlt'} eq $$norm_data{'normalizedAlt'}
                && $$func_block{'normalizedPos'} eq $$norm_data{'normalizedPos'}) {
                %data = %$func_block;
                last;
            } else {
                @data{qw(gene transcript location exon)} = 
                    @{$func_block}{qw(gene transcript location exon)};
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
    $data{'location'} //= '---';
    
    if ($data{'location'} eq 'exonic') {
        $location = "Exon$data{'exon'}";
    } else {
        $location = $data{'location'};
    }

    # Sometimes, for reasons I'm not quite sure of, there can be an array for 
    # the functional annotation in each func block entry.  I think it's safe to
    # take the most severe of the list and to use.  

    # TODO: Refine this.  I can't figure out when and in what context we get 
    # this multiple function annotations bit. So, for now, let's just print 
    # them all out and see what the trend looks like, and the figure it out 
    # from there.  
    ($function) = join('|', @$function) if ref $function eq 'ARRAY';

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
    return ($gene_class, $variant_class, $gene_name, $transcript, $hgvs, 
        $protein, $function, $location );
}

sub normalize_variant {
    # Borrowed from ThermoFisher's vcf.py script to convert IR VCFs. Trim from 
    # both ends until only unique sequence left
    my ($ref,$alt,$pos) = @_;
    my ($norm_ref, $norm_alt);

    my ($rev_ref, $rev_alt, $position_delta) = rev_and_trim($ref, $alt);
    ($norm_ref, $norm_alt, $position_delta) = rev_and_trim(\$rev_ref, 
        \$rev_alt);

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

    my $on  = colored( "On", 'bold green on_black');
    my $off = colored( "Off", 'bold red on_black');

    if ($quiet) {
        print "$info OVAT filter status: ";
        ($ovat_filter) ? print "$on!\n" : print "$off.\n";
        print "$info Hotspot ID filter status: ";
        ($hotspots) ? print "$on!\n" : print "$off.\n";
        print "$info NOCALLs output to results: ";
        ($nocall) ? print "$off!\n" : print "$on.\n";
        print "$info Reference calls output to results: ";
        ($noref) ? print "$off!\n" : print "$on.\n";
    }
    
    # First run OVAT filter; no need to push big list of variants through
    # other filters.
    $data = ovat_filter($data) if $ovat_filter;
    $data = hs_filtered($data) if $hotspots;

    # Determine filter to run, and if none, just return the full set of data.
    my @selected_filters = grep { scalar @{$$filter{$_}} > 0 } keys %$filter;
    if (@selected_filters > 1) {
        print "ERROR: Using more than one type of filter is redundant and not ",
            "accepted at this time. Filters chosen: ";
        print join(', ', @selected_filters), "\n";
        exit 1;
    }
    return $data unless @selected_filters;

    # If we're running a fuzzy position lookup, need to configure things a bit 
    # first
    my @fuzzy_pos;
    if ( $fuzzy ) {
        my $re = qr/(.*).{$fuzzy}/;
        @fuzzy_pos = map { /$re/ } @{$$filter{position}};
        @{$$filter{'fuzzy'}} = @fuzzy_pos;
        $selected_filters[0] = 'fuzzy';
    } 

    # Now run filters
    run_filter($selected_filters[0], \@{$$filter{$selected_filters[0]}}, 
        $fuzzy, $data, \%filtered_data);

    return \%filtered_data;
}

sub run_filter {
    # Run the established filter on the data and return the result.
    my ($term, $filter_vals, $fuzzy, $data, $results) = @_;
    my %mapped_terms = (
        'position' => 0,
        'gene'     => 11,
        'hsid'     => 10,
    );

    my $index = $mapped_terms{$term};
    # Things shift around if using cfDNA...next iteration store data in a hash!
    $index-- unless ($cfdna or (grep { $term eq $_ } qw(position fuzzy)));

    ($term eq 'hsid') ? ($term = uc($term)) : ($term = ucfirst($term));
    print "$info Running the $term filter...\n";

    for my $variant (keys %$data) {
        # if we're running a fuzzy lookup, we need a regex match and have to 
        # do it a little differently
        if ($term eq 'Fuzzy') {
            if ( grep { $$data{$variant}[0] =~/$_.{$fuzzy}$/ } @$filter_vals ) {
                @{$$results{$variant}} = @{$$data{$variant}};
            }
        } else {
            if ( grep { $$data{$variant}[$index] eq $_ } @$filter_vals ) {
                @{$$results{$variant}} = @{$$data{$variant}};
            }
        }
    }
    return $results;
}

sub ovat_filter {
    # Filter out calls that are not oncomine reportable variants
    my $data = shift;
    print "$info Running ovat filter\n" if $quiet;
    for my $variant ( keys %$data ) {
        delete $$data{$variant} if $$data{$variant}->[18] eq '---';
    }
    return $data;
}

sub hs_filtered {
    # Filter out calls that are not oncomine reportable variants
    my $data = shift;
    print "$info Running Hotspots filter\n" if $quiet;
    for my $variant ( keys %$data ) {
        delete $$data{$variant} if $$data{$variant}->[10] eq '.';
    }
    return $data;
}

sub format_output {
    # Format and print out the results
    # Fields in output data:
    #     0. chr:pos             10. Gene
    #     1. REF                 11. Transcript
    #     2. ALT                 12. CDS
    #     3. Filter (NOCALL)     13. AA
    #     4. Filter Reason       14. Location
    #     5. VAF                 15. Function
    #     6. TotCov              16. OVAT VC
    #     7. RefCov              17. tvc/lia flag
    #     8. AltCov
    #     9. VarID

    my ($data, $filter_list, $header, $output_type) = @_;

    # DEBUG
    # Dump first entry for help in figuring out array indices
    #for (keys %$data) {
        #my $i = 0;
        #for my $v (@{$$data{$_}}) {
            #print "$i: $v\n";
            #$i++;
        #}
        #exit;
    #}
    
    if ($output_type eq 'vcf') {
        make_vcf($header, $data);
    }
   
    # Default starting values.
    my $ref_width    = 5;
    my $alt_width    = 5;
    my $varid_width  = 10;
    my $filter_width = 17;
    my $cds_width    = 7;
    my $aa_width     = 7;
    my $func_width   = 9;

    if (%$data) {
        my ($calc_ref_width, $calc_alt_width, $calc_varid_width) = 
            field_width($data, [1,2,9]);

        # Have to pre-declare and set to 0, else we will get warning when no 
        # opt
        my ($calc_filter_width, $calc_cds_width, $calc_aa_width, 
            $calc_func_width) = (0)x4;
        ($calc_filter_width) = field_width($data, [4]) unless $nocall;
        if ($annots) {
            # Need to figure out the index positions of cds, aa, and func,
            # depending on whether we have a cfDNA assay or not.
            my @i;
            ($cfdna) ? (@i = [13,14,16]) : (@i = [12,13,15]);
            ($calc_cds_width, $calc_aa_width, $calc_func_width) = 
                field_width($data, @i);
        }

        # Use calculated values unless defaults are bigger.
        $ref_width = $calc_ref_width unless $ref_width > $calc_ref_width;
        $alt_width = $calc_alt_width unless $alt_width > $calc_alt_width;
        $varid_width = $calc_varid_width unless $varid_width > $calc_varid_width;
        $filter_width = $calc_filter_width unless $filter_width > $calc_filter_width;
        $cds_width = $calc_cds_width unless $cds_width > $calc_cds_width;
        $aa_width = $calc_aa_width unless $aa_width > $calc_aa_width;
        $func_width = $calc_func_width unless $func_width > $calc_func_width;
    }
    
    # Easier to store all formatter elements in a hash for string construction?
    my %string_formatter = (
        'CHROM:POS'             => '%-17s',
        'REF'                   => "%-${ref_width}s",
        'ALT'                   => "%-${alt_width}s",
        'VAF'                   => "%-9s",
        'TotCov'                => "%-8s",
        'RefCov'                => "%-8s",
        'AltCov'                => '%-8s',
        'VarID'                 => "%-${varid_width}s",
        'Filter'                => '%-8s',
        'Filter_Reason'         => "%-${filter_width}s",
        'Gene'                  => '%-13s',
        'Transcript'            => '%-15s',
        'CDS'                   => "%-${cds_width}s",
        'AA'                    => "%-${aa_width}s",
        'Location'              => '%-13s',
        'Function'              => "%-${func_width}s",
        #'oncomineGeneClass'     => '%-21s',
        'oncomineVariantClass'  => '%s', # Since last field don't set a width.
        'LOD'                   => '%-7s',
    );

    # Set up the output header and the correct format string to use.
    my @header = qw( CHROM:POS REF ALT VAF TotCov RefCov AltCov VarID );
    splice(@header, 3, 0, ('Filter', 'Filter_Reason')) unless ($nocall);

    # Add in the LOD field if running cfDNA.  But, VAF position can change
    # depending on whether outputting filter or not.  So, find that first.
    my ($vaf_index) = grep { $header[$_] eq 'VAF' } 1..$#header;
    splice(@header, $vaf_index+1, 0, 'LOD') if ($cfdna);

    if ($annots) {
        push(@header, qw(Gene Transcript CDS AA Location Function));
        # Add OVAT annots and expand function column width if we have OVAT 
        # annots.
        push(@header, 'oncomineVariantClass') if $ovat_annot;
    }

    select $out_fh;
    my $format_string = join(' ', @string_formatter{@header}) . "\n";
    printf $format_string, @header;

    if (%$data) {
        my @output_data;
        for my $variant ( sort { versioncmp( $a, $b ) } keys %$data ) {
            # if not outputting nocall, remove fields 3 and 4; always remove 
            # genotype field.
            ($nocall) 
                ? (@output_data = @{$$data{$variant}}[0,1,2,5..17]) 
                : (@output_data = @{$$data{$variant}});

            # Fill in undef slots with NULL
            @output_data[9..13] = map { $_ //= 'NULL' } @output_data[9..13];
            printf $format_string, @output_data;
        }
    } else {
        # Handle null result reporting depending on the filter used.
        if ($ovat_filter) {
            print "\n>>> No Oncomine Annotated Variants Found! <<<\n";
            exit;
        }
        if (@{$$filter_list{'gene'}}) {
            print "\n>>> No Variants Found for Gene(s): " . 
                join(', ', @{$$filter_list{gene}}), "! <<<\n";
            exit;
        }
        if (@{$$filter_list{'hsid'}}) {
            print "\n>>> No Variants Found for Hotspot ID(s): " . 
                join(', ', @{$$filter_list{hsid}}), "! <<<\n";
            exit;
        }
        if (@{$$filter_list{position}}) {
            my @positions;
            for my $query (@{$$filter_list{position}}) {
                ($fuzzy) 
                    ? push(@positions, 
                        (substr($query, 0, -$fuzzy) . '*'x$fuzzy)) 
                    : push(@positions, $query);
            }
            print "\n>>> No variant found at position(s): ", 
                join(', ', @positions), "! <<<\n" and exit;
        } 
    }
}

sub field_width {
    # Load in a hash of data and an array of indices for which we want field 
    # width info, and output an array of field widths to use in the format 
    # string.
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
    my @sorted_lens = sort { versioncmp($b, $a) } @lens;
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
        lod gtr af fro ro fao ao dp);
    my %foo;

    @foo{@fields} = map{chomp; $_} @$data;

    print '='x25, "  DEBUG  ", "="x25, "\n";
    dd \%foo;
    print '='x59, "\n";
}

sub __load_tscript_file {
    my $tscript_file = shift;
    my %ids;

    open (my $fh, "<", $tscript_file);
    while (<$fh>) {
        chomp;
        next if /^#/;
        my @elems = split(',');
        @{$ids{$elems[0]}} = @elems[1..$#elems];
    }
    return %ids;
}

sub get_can_tran {
    # IR Annotation can get strange with overlapping genes. Figure out the 
    # canonical transcript for filtering, based on transcript list in resources
    # dir in the current package.
    my $gene = shift;
    my $tscript_file = abs_path(dirname($0)) . '/resources/refseq.txt';
    die "Error: Can not find the canonical transcript file!\n" 
        unless -f $tscript_file;

    # Use 'state' variable feature to prevent us reinitializing the hash every
    # single time we want to do a look up, while still allowing us to scope this
    # only when really needed (i.e. not if we don't have an IR file).
    state %tscript_ids;
    %tscript_ids = __load_tscript_file($tscript_file) if ! %tscript_ids;

    if ($tscript_ids{$gene}) {
        return $tscript_ids{$gene};
    } else {
        return 'None';
    }
}

sub __make_header {
    my $header = shift;
    my @captured;
    no warnings;
    
    my @wanted = qw(##fileformat ##fileDate ##source ##reference ##sampleGender
        ##sampleDisease ##contig ##INFO=<ID=AF ##INFO=<ID=AO ##INFO=<ID=DP 
        ##INFO=<ID=RO ##FORMAT=<ID=AF ##FORMAT=<ID=AO ##FORMAT=<ID=DP 
        ##FORMAT=<ID=RO ##FORMAT=<ID=GT );

    for my $line (@$header) {
        push(@captured, $line) if grep { $line =~ /$_/ } @wanted;
    }
    push(@captured, $header->[-1]);

    # Update the date to today's
    my $today = localtime->strftime('%Y-%m-%d');
    ($captured[1] =  $captured[1]) =~ s/\d{8}/$today/;
    return \@captured;
}

sub __make_vcf_line {
    my $data = shift;

    my ($chr, $pos) = split(/:/, $data->[0]);
    my $vaf = $data->[5];
    my $gt;
    if ($vaf eq '.') {
        $gt = './.';
    }
    elsif ($vaf > 50) {
        $gt = '1/1';
    }
    else {
        $gt = '0/1';
    }

    my $info = sprintf('AF=%s;DP=%s;RO=%s;AO=%s', @$data[5..8]);
    my $sample_data = join(':', $gt, @$data[5..8]);

    my $vcf_line = join("\t", $chr, $pos, $data->[9], $data->[1], $data->[2], '.',
        $data->[3], $info, 'GT:AF:DP:RO:AO', $sample_data) . "\n"; 
    return $vcf_line;
}

sub make_vcf {
    # Output data in VCF format instead of a pretty table.
    my ($header, $data) = @_;
    my $final_header = __make_header($header);

    my @vcf_lines;
    for (sort { versioncmp( $a, $b ) } keys %$data) {
        push(@vcf_lines, __make_vcf_line($$data{$_}));
    }

    print {$out_fh} $_ for @$final_header;
    print {$out_fh} $_ for @vcf_lines;
    exit;
}

sub __exit__ {
    my ($line, $msg) = @_;
    $msg //= '';
    print "\n\n";
    print colored("Got exit message at line $line with message: $msg", 
        'bold white on_green');
    print "\n";
    exit;
}
