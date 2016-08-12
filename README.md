Miscellaneous Bioinformatics Scripts
--
* <b>cephAnalysis.pl</b>
  Script to read in tabular variant call information for multiple CEPH runs and generate some metrics about the number of variants seen
  and their frequency

* <b>getseq.pl</b>
  Script to query the UCSC DAS server and pull out sequence information given a position string.  Can use either
  a single string entry:

           <<<chrx:123435
  
  or can use a file with positions in the same format and get a batch output.  Script only beta version and needs 
  some more tweaks when I have time.

* <b>vcfExtractor.pl</b>
  Script to parse Ion Torrent specific VCF files and pull out variant data.  This works with TS v4.2 and v5.0
  files.  Some parts are compatible with Ion Reporter generated VCF files.

* <b>genemed_compare.pl</b>
  Starting with either the XLS or CSV derived from the XLS report, compare GeneMed output from different versions
  of the Torrent Variant Caller (TVC) and output a Venn Diagram and Tables of the results

* <b>readlength_histogram.pl</b>
  Read in an Ion Torrent BAM file and generate a readlength histogram plot from the sample.  This script will require
  the `Statistics::R` perl module, as well as, the `ggplot2` library in R.  You can find [Statistics::R here](http://search.cpan.org/~fangly/Statistics-R/lib/Statistics/R.pm).
  You can find the  most excellent [ggplot2 package here](http://ggplot2.org).
