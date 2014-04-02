#!/usr/bin/perl
# Read in either a file containing a human genome position or a position from the CLI, and get the sequence
# from the UCSC DAS server.  The input can either be in the following forms:
#
#     chrx:start,stop
#     chrx:start-stop
#     chrx	start	stop
#     x	start	stop
#     x:start,stop
#
# 6/15/2013 - D Sims
##############################################################################################################

use warnings;
use strict;
use Data::Dumper;
use LWP::Simple;
use XML::Twig;
use Sort::Versions;

( my $scriptname = $0 ) =~ s/^(.*\/)+//;
my $version = "v0.7.040214"; # add version sort to output
my $description = <<"EOT";
Program to retrieve sequence from the UCSC DAS server.  Enter sequence coordinates in the form of 'chr:start-stop',
and the output will be sequence from hg19, padded by 10 bp.  Extra padding can be added with the '-p' option.  
The input is flexible and can accomodate the following formats:

    chrx:start,stop
    chrx:start-stop
    chrx:start..stop (direct cp / paste from COSMIC)
    x start   stop
    x:start,stop

The sequence can be fed from a file list or can be entered one by one on the CLI.
EOT

my $usage = <<"EOT";
USAGE: $0 [options] <chr:start-stop>
    -p, --pad       Pad the output sequence (Default is 10bp).
    -o, --output    Send output to custom file.  Default is STDOUT.
    -v, --version   Version information
    -h, --help      Print this help information
EOT

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version_info {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

my $help = 0;
my $verInfo = 0;
my $outfile;
my $padding = 10;

while ( scalar( @ARGV ) > 0 ) {
	last if ( $ARGV[0] !~ /^-/ );
	my $opt = shift;
	if    ( $opt eq '-v' || $opt eq '--version' )  { $verInfo = 1; }
    elsif ( $opt eq '-p' || $opt eq '--pad' )      { $padding = shift; }
	elsif ( $opt eq '-h' || $opt eq '--help' )     { $help = 1; }
	elsif ( $opt eq '-o' || $opt eq '--output' )   { $outfile = shift; }
	else {
		print "Invalid option: '$opt'\n";
		print "$usage\n";
		exit 1;
	}
}

help if $help;
version_info if $verInfo;

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
	open( $out_fh, ">", $outfile ) || die "Can't open the output file '$outfile' for writing: $!";
} else {
	$out_fh = \*STDOUT;
}

#########------------------------------ END ARG Parsing ---------------------------------#########

my @queries;
my %lookup;

# Read in position file and create a query line. Position has to be in format chr*\tstart\tstop
while (<>) {
	chomp;
	my $queryline;
	
	# Get query components to set up the query string.
	my ( $chr, $start, $end )  =  $_ =~ /^(?:chr)?(\d+).(\d+)(?:\D(?:[\-,\. ])?(\d+))?$/;
    $end = $start if ( ! $end );

    # Add some padding to the start and end positions
    $start -= $padding;
    $end += $padding;

	# Make sure that the sequence of start and end is right before loading.  Gets a silent error otherwise!
	( $end > $start ) ? $queryline = "chr$chr:$end,$start" : $queryline = "chr$chr:$start,$end";
	push( @queries, $queryline );
}

my $URL = "http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=";

# Query the UCSC DAS server and extract sequence from the resulting XML file.
my %result;
foreach ( @queries ) {
	my $query = $URL . $_;
	my $twig = XML::Twig->new();
	$twig->parse( LWP::Simple::get( $query ) );
	for my $seq ( $twig->findnodes( '//DNA' ) ) {
		( my $returnSeq = $seq->text ) =~ s/\R/\n/g;
	    $result{$_} = $returnSeq;
	}
}

# Print out the formatted results
for ( sort { versioncmp( $a, $b ) } keys %result ) {

    ( my $seq_id = $_) =~ s/,/-/g;
    print $seq_id."\n".uc($result{$_})."\n";

    # TODO: Configure fancy MOI output
    # Fancy formatting with a box around the MOI...needs work!
	#( my $anchorStr = $result{$_} ) =~ s/(.{10})(.)(.{10})/$1 \[$2\] $3/;
	#printf "%-30s%s\n", $_, $anchorStr;
    
}
