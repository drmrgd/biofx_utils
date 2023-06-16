#!/usr/bin/env perl
# Long overdue update to sequence retrieval script. This one is going to use all
# local resources rather than the ucsc method from the original.  While this
# requires having some local files like the reference sequence, which is a bit
# on the large side, as well as samtools installed, this way is better and
# faster than the original.  
#
# TODO: Right now we only have a preliminary draft version without many
# features. Work on adding a few comfort features to make this a little better:
#   1. Batch queries from a file.
#   2. Batch queries as a comma separated list.
#   3. Better output formatting and reporting. 
#   4. Ability to output just the retrieved sequences for programmatic use.
#
# 4/28/2022 - D Sims
################################################################################
use strict;
use warnings;
use autodie;

use Data::Dump;
use File::Basename;
use Getopt::Long;
use Sort::Versions;

my $version = "3.0.061623";

my $reference = "$ENV{HOME}/Dropbox/reference/hg19/hg19.fasta";

my $scriptname = basename($0);
my $description = <<"EOT";
Input a position and retrieve sequence from a reference file. Query should be 
formatted as chromosome:start-end.  Note that the addition of 'chr' to the 
beginning of the chromosome string is optional.  Also, if one just wants a 
single base in the return, one can omit the end position in the query.

EOT

my $usage = <<"EOT";
$scriptname [options] <position>

Options:
    -r, --reference    Reference sequence to use. Default: $reference.
    -v, --version      Print the version information and exit.
    -h, --help         Print this help text and exit.

EOT

my $help;
my $verinfo;

GetOptions(
    'r|reference=s'     => \$reference,
    'v|version'         => \$verinfo,
    'h|help'            => \$help
) or die "$usage";

sub version {print "$scriptname - v$version\n";}
sub help { version(); print "$description\n$usage\n"; }

version() and exit if $verinfo;
help() and exit if $help;

my $input_query = shift;
my $formatted_query = parse_query($input_query);

sub parse_query {
    my $query = shift;

    my ($chr, $start, $end) = $query =~ /^((?:chr)?[12]?[0-9XY]):(\d+)(?:-(\d+))?/;
    unless($chr and $start) {
        die "ERROR: Input query `$query` does not appear to be properly " .
            "formatted! Can not process!\n";
    }
    $end //= $start;
    my $formatted_query = sprintf("%s:%s-%s", 
        ($chr =~ /^chr/) ? $chr : "chr$chr",
        $start,
        $end
    );
    return $formatted_query;
}

open(my $stream, "-|", "samtools faidx $reference $formatted_query");
# my $ret_seq = (map{chomp; $_} <$stream>)[1];
my @ret_seq = map{ chomp; $_ } <$stream>;
my $outseq = join('', @ret_seq[1..$#ret_seq]);
print "[ $formatted_query ]: $outseq\n";
