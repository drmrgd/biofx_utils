#!/usr/bin/perl
# Read in a FASTQ file and print out a read length histogram table based on the bin size specified on 
# the commandline.  The default size is 5.  Right now can only read in 1 FASTQ file at a time.  More to
# come later on when I can figure out how to store them.
#
# 5/6/13
use warnings;
use strict;
use autodie;

use File::Basename;
use Getopt::Long qw(:config bundling auto_abbrev no_ignore_case);
use Data::Dump;
use Cwd 'abs_path';
use List::Util qw(sum max min);
use POSIX qw(ceil floor);
use Sort::Versions;
use Statistics::R;

my $scriptname = basename($0);
my $version = "v2.0.0_090116";
my $description = <<"EOT";
Program to read BAM file and create a read length histogram using the Statistics::R module. Can load a 
FASTQ file directly so that we don't need to convert first.  
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <sequence_file>
    -f, --fastq      Sample is a FASTQ file and no conversion needs to be done.
    -s, --sample     Sample name (DEFAULT: FASTQ filename)
    -b, --bin        Bin size (DEFAULT: 10)
    -z, --zoom       Zoom in on bins by only outputting data for which there are bins.  Default is range from 0-210.
    -o, --output     Output file (DEFAULT: rlh_data.txt)
    -v, --version    Version Information
    -h, --help       Display the help information
EOT

my $help;
my $ver_info;
my $outfile = "rlh_data.txt";
my $binwidth = 10;
my $sample;
my $fastq_file;
my $zoom;

GetOptions( "fastq|f"     => \$fastq_file,
            "sample|s=s"  => \$sample,
            "bin|b=i"     => \$binwidth,
            "zoom|z"      => \$zoom,
            "output|o=s"  => \$outfile,
            "version|v"   => \$ver_info,
            "help|h"      => \$help )
        or die $usage;

sub help {
    printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
    exit;
}

sub version {
    printf "%s - %s\n", $scriptname, $version;
}

help if $help;
version if $ver_info;


###------------------------------ END COMMANDLINE ARGS --------------------------###
my $seqfile = shift;
unless ($seqfile) {
    print "ERROR: no BAM or FASTQ files loaded!\n";
    print "$usage\n";
    exit 1;
}

# Check the file type is correct and covert the BAM file to a FASTQ if need be.
my $fastq;
if ($fastq_file) {
    file_check(\$seqfile, 'fastq');
    $fastq = $seqfile;
} else {
    file_check(\$seqfile, 'bam');
    $fastq = convert_bam($seqfile);
}

my $sample_name;
($sample) ? ($sample_name = $sample) : (($sample_name = $fastq) =~ s/\.fastq//);

my (%readlength, $total_reads);
open( my $input_fh, "<", $fastq) || die "FASTQ file '$fastq' can not be opened: $!\n";
while (<$input_fh>) {
    chomp;
    if ( $_ =~ /^@.*:\d+:\d+$/ ) {
        $total_reads++;
        my $seq = <$input_fh>;
        my $len = length($seq);
        $readlength{$len}++;
    }
}
close( $input_fh );

# Call rlh sub and print out results
my %rlhData;
($zoom) 
    ? (%rlhData = read_length_histogram( \%readlength, $binwidth, $outfile))
    : (%rlhData = static_bin_histogram(\%readlength, $binwidth, $outfile));

# Draw histo with R
plot_histogram($sample_name,$total_reads);

sub file_check {
    my ($file,$type) = @_;
    my ($extension) = $$file =~ /\.(\w+)$/;
    if ($type eq 'fastq') {
        die "ERROR: '$$file' does not appear to be a FASTQ file! Did you accidentally use the '-f' option?\n" unless $extension eq 'fastq';
    }
    elsif ($type eq 'bam') {
        die "ERROR: '$$file' does not appear to be a BAM file! Did you forget to add the '-f' option?\n" unless $extension eq 'bam';
    }
    return;
}

sub static_bin_histogram {
    my ($rl_data, $binwidth, $outfile) = @_;
    my @bins = map{ $_."-".($_+$binwidth) } (15..210);
    my %result = map{ $_ => 0 } @bins;

    while (my ($length, $count) = each %$rl_data) {
        for my $bin (@bins) {
            my ($start, $end) = split(/-/, $bin);
            if ($length > $start && $length < $end) {
                $result{$bin} += $count;
            }
        }
    }
    open(my $out_fh, ">", $outfile);
    print {$out_fh} "Length\tCount\n";
    for my $lengths( sort{versioncmp( $a, $b )} keys %result) {
        print {$out_fh} "$lengths\t$result{$lengths}\n";
    }
    return %result;
}

sub read_length_histogram {
    # Read in %readlength data and output RLH table
    my ($rlh,$binwidth,$outfile) = @_;
	my ( %result, %returnHash,$maxLength,$minLength );

    print "Calculating read length counts based on a bin size of '$binwidth'...";
    open (my $out_fh, ">", $outfile);

	while ( my ( $length, $count ) = each %$rlh) {
		push( my @lengths, $length );
		$maxLength = max( @lengths );
		$minLength = min( @lengths );
		$result{ceil( ( $length + 1 ) / $binwidth ) - 1} += $count;
	}
	
	print {$out_fh} "Length\tCount\n";
	for my $lengths ( sort { $a <=> $b } keys %result ) {
		my $bin = $lengths * $binwidth;
		my $binstring = $bin. "-" . ($bin + ( $binwidth - 1 ));
		print {$out_fh} $binstring, "\t$result{$lengths}\n";
		$returnHash{$binstring} = $result{$lengths};
	}
    print "Done!\n";
	return %returnHash;
}

sub plot_histogram {
    #my $sample_name = shift;
    my ($sample_name, $total_reads) = @_;
    my $data_file = abs_path($outfile);
    print "Generating histogram plot...";

    my $R = Statistics::R->new();
    $R->run( q/library(ggplot2)/ );
    $R->run( q/library(scales)/ );
    $R->run( qq/data<-data.frame(read.table("$data_file",sep="\t",header=TRUE))/ );
    $R->run( q/data$ordered_bins<-factor(data$Length,as.character(data$Length))/ );

	my $rlh_plot = <<EOF;
    p <- ggplot(data, aes(x=ordered_bins, y=Count)) +
         geom_bar(stat="identity") + 
         labs(title="$sample_name Read Length Histogram (Total Reads: $total_reads)", x="Bins", y="Count") +
         theme_bw() + 
         theme(axis.text.x = element_text(angle=45, hjust=1)) + 
         scale_y_continuous(labels=comma)
EOF
    $R->run(qq/$rlh_plot/);
    $R->run(qq/ggsave(file="${sample_name}_rlh.pdf",width=10, height=5, dpi=100, plot = p)/);
    print "Done!\nHistogram '${sample_name}_rlh.pdf' generated for this file.\n";
}

sub convert_bam {
    my $bam_file = shift;
    (my $fastq_file = $bam_file) =~ s/\.bam$/\.fastq/;
    print "Coverting '$bam_file' to FASTQ file: $fastq_file...";
    die "ERROR: Samtools not found in \$PATH!" unless `which samtools`;
    `samtools bam2fq $bam_file 2> /dev/null > $fastq_file`;
    print "Done!\n";
    return $fastq_file;
}
