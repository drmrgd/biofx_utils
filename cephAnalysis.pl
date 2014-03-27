#!/usr/bin/perl
# Read in all the variant call files, filter and truncate the data, and associate with the run ID.  Set
# up to be able to find the common variants in all the runs.  Add column to indicate if it's one of the 
# variants located within the NIST CEPH standard VCF file.
#
# Created: 2/27/2013 - Dave Sims
#
##########################################################################################################

use warnings;
use strict;
use autodie;

use List::Util qw{ min max sum };
use Getopt::Long qw{ :config no_ignore_case };
use File::Basename;
use Data::Dump;

#( my $scriptname = $0 ) =~ s/^(.*\/)+//;
my $scriptname = basename($0);
my $version = "v4.0.032714";
my $description = <<EOT;
Program to read in all of the variant call table files from an Ion Torrent run on CEPH, and report out
the ID and number of times each variant is seen.  This is used to track the number of variants reported 
in that sample in order to identify variants that should be excluded from the reportable range of the 
assay.
EOT

my $usage = <<EOT;
USAGE: $scriptname [options] <input_files>
    -C, --Classic    Use TVCv3.2.1 formatted data files.  Will be deprecated.
    -c, --coverage   Coverage cutoff (default is 450)
    -p, --preview    Write output only to STDOUT to preview the results.
    -o, --output     Write output to file (default is STDOUT)
    -v, --version    Print version information
    -h, --help       Print help information
EOT

my $verinfo;
my $help;
my $covfilter = 450;
my $preview;
my $output;
my $classic;

# Set up some commandline opts
GetOptions( "Classic"       => \$classic,
            "coverage=i"    => \$covfilter,
            "preview"       => \$preview,
            "output=s"      => \$output,
            "version"       => \$verinfo,
            "help"          => \$help,
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
my $totRuns;
if ( ! @ARGV ) {
    print "ERROR: no files loaded for analysis\n";
    exit 1;
} else {
    @filesList = @ARGV;
    $totRuns = scalar @filesList;
}

# Set up Tally output file.
my $outfile = "CEPH_".$totRuns."_Run_Variant_Tally.tsv";
my $out_fh;
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
    "v3.2.1"  => { "varid" => [qw(2 0 1 6 7)],
                   "data"  => [qw(0 1 2 3 6 7 8 10)],
                 },
    "v4.0.2"  => { "varid" => [qw(12 0 1 15 16)],
                   "data"  => [qw(0 1 12 13 15 16 6 18)],
                 }
);

proc_datatable( \@filesList, \%fields );

sub proc_datatable {
    my $files = shift;
    my $data_fields = shift;

    # Get fields to use for spliting the tables
    my $version = $classic ? 'v3.2.1' : 'v4.0.2';
    print "Processing '$version' data...\n";
    my @varid_index = @{$$data_fields{$version}{'varid'}};
    my @field_index = @{$$data_fields{$version}{'data'}};

    for my $file ( @$files ) {
        open( my $in_fh, "<", $file );
        while (<$in_fh>) {
            next if ( /Chrom/ || /Absent/ || /No Call/ );
            my @fields = split;
            if ( $fields[$field_index[7]] > $covfilter ) {
                my $varid = join( ':', @fields[@varid_index] );
                push( @{$all_variants{$varid}}, [@fields[@field_index]] );
                #push( @{$all_variants{$varid}->{@fields[@field_index[0..5]]}}, [@fields[@field_index[6,7]]] );
                #push( @{$all_variants{$varid}}, [$fields[$field_index[6]], $fields[$field_index[7]]] ); 
                push( @{$var_freq{$varid}}, $fields[$field_index[6]] );
                push( @{$var_cov{$varid}}, $fields[$field_index[7]] );
            }
        }
        close $in_fh;
    }
    return;
}

# Get some field width data
my ( $rwidth, $awidth ) = field_width( \%all_variants );
#print "rwidth: $rwidth\nawidth: $awidth\n";

#dd \%all_variants;
#exit;

# Get statistics about each variant and generate a formated hash table to print out the results with 
my %results;
foreach my $variant ( keys %all_variants) {
	my $count = @{$all_variants{$variant}};
    my ( $gene, $chr, $pos, $ref, $alt ) = split( /:/, $variant );
	
	# Get min, max, median coverage and frequency info	
	my ( $minCov, $maxCov, $meanCov ) = stats( \@{$var_cov{$variant}} );
	my ( $minFreq, $maxFreq, $meanFreq ) = stats( \@{$var_freq{$variant}} );
    my $detection_freq = sprintf( "%0.2f", ($count/$totRuns)*100);

    my $format = "%-8s %-8s %-12d %-${rwidth}s %-${awidth}s %-10s %-7d %-7d %-7d %8.2f%% %8.2f%% %8.2f%%";
    my $varLine = sprintf( $format, $gene, $chr, $pos, $ref, $alt, "$count/$totRuns", $minCov, $meanCov, $maxCov, $minFreq, $meanFreq, $maxFreq );

    $results{$variant} = [$varLine, $detection_freq];
}

# Print out the collected summary data 
my @header_cols = qw{ Gene Chr Position Ref Var Count MinCov MedCov MaxCov MinVAF MedVAF MaxVAF };
my $header = sprintf( "%-8s %-8s %-12s %-${rwidth}s %-${awidth}s %-10s %-7s %-7s %-7s %9s %9s %9s", @header_cols ); 
my $title = "Frequency of detected variants with at least $covfilter reads in $totRuns CEPH runs";

print $out_fh "$title\n\n";
print $out_fh "$header\n";

for my $variant ( sort { $results{$b}[1] <=> $results{$a}[1] } keys %results ) {
    print $out_fh $results{$variant}[0], "\n";
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
        my ($ref, $alt) = $var =~ /.*?(\w+):(\w+)$/;
        #print "ref: $ref\nalt: $alt\n";
        
        $rwidth = length( $ref ) + 3  if ( length( $ref ) > $rwidth ); 
        $awidth = length( $alt ) + 3  if ( length( $alt ) > $awidth );
    }
    return( $rwidth, $awidth );
}
