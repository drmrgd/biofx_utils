#!/usr/bin/perl
# Based on advice from Jingwei, this will fix an RNA BAM file that is missing the correct components to import into IR
# when we save the wrong version from the backend.  
use strict;
use warnings;
use autodie;

sub usage {
    print "USAGE: $0 <rna_bam.bam>\n";
    exit;
}

my $bam_file = shift;
$bam_file or ( print "ERROR: You must load an RNA BAM file!\n" and usage() );
usage() if $bam_file eq '-h';

my @missing_header_elem = qw(FO:TACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCATCGATCGATGTACAGCTACGTACGTCTGAGCA KS:TCAGTCTGGCAACGGCGAT PG:tmap);

open(my $new_header, ">", 'tmp.header');
open( my $bam_header, "-|", "samtools view -H $bam_file");
while (<$bam_header>) {
    my @elems = split;
    if ($elems[0] eq '@RG') {
        print {$new_header} join("\t",@elems, @missing_header_elem), "\n";
    } else {
        print {$new_header} $_;
    }
}
close $new_header;

(my $new_bamfile = $bam_file) =~ s/\.bam/_reheader.bam/;
print "Fixing the RNA BAM header...";
system("samtools reheader tmp.header $bam_file > $new_bamfile") and die $!;
print "Done!\n";

unlink('tmp.header');
