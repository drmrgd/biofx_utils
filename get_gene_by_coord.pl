#!/usr/bin/env perl
# Simple script to read in a chr:pos coord and output the gene at that 
# position.
#
# 12/14/2018 - D Sims
################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Data::Dump;
use Sort::Versions;
use List::Util qw(min max);
use Parallel::ForkManager;
use MCE::Shared;

my $gene_file = dirname($0) . "/resources/gene_coordinates.txt";

my $scriptname = basename($0);
my $version = "v0.8.121718";
my $description = <<"EOT";
Input GRCh37 (hg19) chromosome and position in the format of: 

    chr1:1234567
    
and get back a Hugo Gene name. Can also input a comma separated list of coords,
or a batchfile of coords to lookup.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <coords>
    -b, --batchfile  File of coordinates, one per line, to retrieve.
    -o, --output     Send output to custom file.  Default is STDOUT.
    -v, --version    Version information
    -h, --help       Print this help information
EOT

my $help;
my $ver_info;
my $outfile;
my $batchfile;

GetOptions( 
    "batchfile|b=s" => \$batchfile,
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
if ( scalar( @ARGV ) < 1 ) {
    unless ($batchfile) {
        print "ERROR: Not enough arguments passed to script!\n\n";
        print "$usage\n";
        exit 1;
    }
}

# Build up a lookup table of gene and positions that will be very fast to query
my $gene_dict = read_gene_file($gene_file);

# Create a list of coordinates to look up.
my @coords;
if ($batchfile) {
    open(my $bfh, "<", $batchfile);
    @coords = map{ chomp; $_ } <$bfh>;
    close $bfh;
} else {
    @coords = split(/,/, shift);
}

# Write output to either indicated file or STDOUT
my $out_fh;
if ( $outfile ) {
	open( $out_fh, ">", $outfile );
} else {
	$out_fh = \*STDOUT;
}
################------ END Arg Parsing and Script Setup ------#################
my $pm = new Parallel::ForkManager(48);

# Since we can not make a hash of these data since there won't be unique keys,
# we will use the MCE::Shared module to create memory links to store these data
# in an array.  
my $data = MCE::Shared->array;

for my $coord (@coords) {
    $pm->start and next;
    my $gene = map_gene(\$coord);
    $data->push("$coord|$gene");
    $pm->finish;
}
$pm->wait_all_children;

print_results($data, $out_fh);


sub print_results {
    # NOTE: See the MCE::Shared::Array docs for details how this works. Basically
    # we're storing the results array in a sort of hash where the index is a key
    # and the return result is the value.  There is an iterator method that will
    # then allow use to traverse this structure using the index as the "key" and
    # the results as the "value"
    my ($results, $outfh) = @_;
    my @missing;

    my $iter = $results->iterator;
    while (my ($i, $v) = $iter->()) {
        my ($coord, $gene) = split(/\|/, $v);
        if ($gene eq 'None') {
            push(@missing, $coord);
        } else {
            print {$outfh} "$coord,$gene\n";
        }
    }

    if (@missing) {
        print "WARNING: The following coordinates did not map to a gene: \n";
        print "\t$_\n" for @missing;
    }
}

sub map_gene {
    my $coords = shift;
    my $gene;

    my ($chr, $pos) = split(/:/, $$coords);
    if (exists $gene_dict->{$chr}) {
        for my $range (keys %{$gene_dict->{$chr}}) {
            my ($x, $y) = split(/-/, $range);
            if ($x < $pos and $y > $pos) {
                return get_ent($gene_dict->{$chr}{$range}, $pos);
            }
        }
    }
}

sub get_ent {
    my ($entries, $pos) = @_;
    for my $entry (@$entries) {
        my ($chr, $start, $end, $gene, $strand) = split(/;/, $entry);

        # Pad the start and stop by 200 bp to account for upstream and 
        # downstream regulatory variants that might pop up.
        if (($start-150) < $pos && ($end+150) > $pos) {
            return $gene;
        }
    }
    # Didn't find a result. Probably a coord typo
    return "None";
}

sub read_gene_file {
    my $gene_file = shift;
    my @gene_info;
    my %positions;

    my %data;
    open(my $fh, "<", $gene_file);
    while (<$fh>) {
        chomp;
        my ($chr, $start, $end, $gene, $strand) = split(/\t/);

        # Filter out some things we don't want for most purposes.
        # 1. MiRs
        next if $gene =~ /^MIR/;

        # for now, keep all data in the file.  However, there are many chrom.
        # variants (e.g. 'chr17_ctg5_hap1', 'chr17_gl000205_random', etc.) that
        # we probably don't need.  Also quite a few genes we probably won't ever
        # need.
        push(@{$positions{$chr}}, $start, $end);
        push(@gene_info, [$chr, $start, $end, $gene, $strand]);
    }
    close $fh;
    return build_lookup(\%positions, \@gene_info);
}

sub build_lookup {
    my ($positions, $genes) = @_;
    my %lookup_table;
    my $bin_width = 50000000;

    # Build some bins into which we'll put the genes of interest. Use the max
    # and min positions, along with the $bin_width var, to determine how many
    # bins we need.  Might need to adjust a bit.
    for my $chr ( sort { versioncmp($a, $b) } keys %$positions) {
        my $min = min(@{$positions->{$chr}});
        my $max = max(@{$positions->{$chr}});
        my $floor = int($min/$bin_width) * $bin_width;
        while ($floor < $max) {
            my $range = sprintf( "%s-%s", $floor+1, $floor += $bin_width);
            $lookup_table{$chr}->{$range} = [];
        }
    }

    # Now that we have our bins, let's put the genes into place.
    for my $gene (@$genes) {
        if ($lookup_table{$gene->[0]}) {
            insert_gene($lookup_table{$gene->[0]}, $gene);
        }
    }
    return \%lookup_table;
}

sub insert_gene {
    my ($regions, $gene) = @_;
    for my $range (keys %$regions) {
        my ($start, $end) = split(/-/, $range);
        if ($gene->[1] > $start and $gene->[2] < $end) {
            push(@{$regions->{$range}}, join(';', @$gene)) and return;
        }
    }
}
