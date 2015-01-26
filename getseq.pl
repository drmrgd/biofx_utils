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
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Data::Dump;
use LWP::Simple;
use XML::Twig;
use Sort::Versions;

my $scriptname = basename($0);
my $version = "v1.0.0_012615";
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
USAGE: $scriptname [options] <chr:start-stop>
    -p, --pad       Pad the output sequence (Default is 10bp).
    -b, --batch     Load up a batch file of positions to search.
    -o, --output    Send output to custom file.  Default is STDOUT.
    -v, --version   Version information
    -h, --help      Print this help information
EOT

my $help;
my $ver_info;
my $outfile;
my $padding = 10;
my $batch_file;

GetOptions( "padding|p=i"   => \$padding,
            "batch|b=s"     => \$batch_file,
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

# Make sure enough args passed to script
if ( scalar( @ARGV ) < 1 && ! $batch_file ) {
    print "ERROR: need to enter at least one coordinate to search!\n\n";
    print "$usage\n";
    exit 1;
}

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
	open( $out_fh, ">", $outfile ) || die "Can't open the output file '$outfile' for writing: $!";
} else {
	$out_fh = \*STDOUT;
}

#########------------------------------ END ARG Parsing ---------------------------------#########
my @queries;

if ($batch_file) {
    my @batch_list = proc_batch( $batch_file );
    for my $coord (@batch_list) {
        my $formatted_query = gen_queries($coord);
        push(@queries, $formatted_query);
    }
} else {
    my $input_query = shift;
    my $formatted_query = gen_queries($input_query);
    push(@queries, $formatted_query);
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
    print {$out_fh} ">$seq_id" . uc($result{$_}) . "\n";
}

sub gen_queries {
    my $input_string = shift;

    my ($chr, $start, $end) = $input_string =~ /^(?:chr)?(\d+).(\d+)(?:\D(?:[\-,\. ])?(\d+))?$/;
    $end = $start unless $end;
    $start -= $padding;
    $end += $padding;

    my $formatted_query;
    # TODO:  Bug with the string if someone has put the end before the start.
    ($end > $start) ? $formatted_query = "chr$chr:$end,$start" : $formatted_query = "chr$chr:$start,$end";

    return $formatted_query;
}

sub proc_batch {
    # read in a batch file if many positons needed
    my $input_file = shift;

    open(my $fh, "<", $input_file);
    my @query_list = map { chomp($_); $_ } <$fh>;
    close $fh;

    return @query_list;
}
