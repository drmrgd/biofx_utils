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
use List::Util qw{ min max sum };
use Data::Dumper;

( my $scriptname = $0 ) =~ s/^(.*\/)+//;
my $version = "v3.0.0";
my $description = <<EOT;
Program to read in all of the variant call table files from an Ion Torrent run on CEPH, and report out
the ID and number of times each variant is seen.  This is used to track the number of variants reported 
in that sample in order to identify variants that should be excluded from the reportable range of the 
assay.
EOT

my $usage = <<EOT;
USAGE: $scriptname [options] <input_files>
    -d, --dev        Development Version that works with format of TVCv4.0
    -c, --coverage   Coverage cutoff (default is 450)
    -p, --preview    Write output only to STDOUT to preview the results.
    -o, --output     Write output to file (default is STDOUT)
    -v, --version    Print version information
    -h, --help       Print help information
EOT

sub help {
	printf "%s - %s\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

my $verinfo;
my $help;
my $covfilter = 450;
my $preview;
my $output;
my $dev_version;

# Set up some commandline opts
while ( scalar( @ARGV ) > 0 ) {
	last if ( $ARGV[0] !~ /^-/ );
	my $opt = shift;
    if    ( $opt eq '-d' || $opt eq '--dev' )      { $dev_version = 1; }
	elsif ( $opt eq '-o' || $opt eq '--output' )   { $output = shift; }
	elsif ( $opt eq '-p' || $opt eq '--preview' )  { $preview = 1; }
	elsif ( $opt eq '-c' || $opt eq '--coverage' ) { $covfilter = shift; }
	elsif ( $opt eq '-v' || $opt eq '--version' )  { $verinfo = 1; }
	elsif ( $opt eq '-h' || $opt eq '--help' )     { $help = 1; }
	else {
		print "$scriptname: Invalid Option '$opt'\n";
		print $usage;
		exit 1;
	}
}

if ( $verinfo ) {
	print "$scriptname\t$version\n";
	exit;
}

help if $help;

if ( $covfilter !~ /\d+/ ) {
	print "ERROR: coverage cutoff must be numeric!\n";
	exit 1;
}

my @filesList;
die "No files loaded for analysis!\n" if ( ! defined @ARGV );

my $totRuns;
foreach my $file ( @ARGV ) {
	push( @filesList, $file);
	$totRuns = scalar( @filesList );
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
my ( %varCount, %trimHash );
my ( %varFreq, %varCov );

# Get the variant call data from each input file and start collected in the data in a hash
foreach my $varFile ( @filesList ) {
	( my $filename = $varFile ) =~ s/\.txt//;

	open( VARFILE, "<", $varFile ) || die "Error opening file: $varFile. $!";
	
	while (<VARFILE>) {
        next if ( /Chrom/ || /Absent/ || /No Call/ );
        my @line = split( /\t/, $_ );
        if ( ! defined $dev_version ) {
            #next if (/Chrom/); 
            #my @line = split;
            #my @line = split( /\t/, $_ );
            
            if ( $line[10] > $covfilter ) {
                my $varID = join( ":", @line[2,0,1,6,7] ); #Gene:Chr:Pos:REF:ALT
                $varCount{$varID}++;
                push( @{$trimHash{$varID}}, @line[0,1,2,3,6,7,8,10] );
                push( @{$varFreq{$varID}}, $line[8] );
                push( @{$varCov{$varID}}, $line[10] );
            }
        } else {
            # Set up correct fields for TVCv4.0
            if ( $line[18] > $covfilter ) {
                my $varID = join( ":", @line[12,0,14,15,16] );
                $varCount{$varID}++;
                push( @{$trimHash{$varID}}, @line[0,14,12,13,15,16,6,18] );
                push( @{$varFreq{$varID}}, $line[6] );
                push( @{$varCov{$varID}}, $line[18] );
                #print "$varID\n";
                #next;
            }
        }
	}
}
close( VARFILE );

# Get some field width data
my ( $rwidth, $awidth ) = field_width( \%trimHash );

# Get statistics about each variant and generate a formated hash table to print out the results with 
my %results;
foreach my $variant ( sort keys %varCount ) {
	my $count = $varCount{$variant};
    my ( $gene, $chr, $pos, $ref, $alt ) = split( /:/, $variant );
	
	# Get min, max, median coverage and frequency info	
	my ( $minCov, $maxCov, $meanCov ) = stats( \@{$varCov{$variant}} );
	my ( $minFreq, $maxFreq, $meanFreq ) = stats( \@{$varFreq{$variant}} );
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
        my ($ref,$alt) = $var =~ /.*:(\w+):(\w+)$/;
        
        $rwidth = length( $ref ) + 3  if ( length( $ref ) > $rwidth ); 
        $awidth = length( $alt ) + 3  if ( length( $alt ) > $awidth );
    }
    return( $rwidth, $awidth );
}
