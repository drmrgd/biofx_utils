Miscellaneous Bioinformatics Scripts
--
This is a set of tools I developed to do some basic, but helpful things.

* <b>getseq.pl</b>
  Script to query the UCSC DAS server and pull out sequence information given a position string.  Can use either
  a single string entry:

           <<<chrx:123435
  
  or can use a file with positions in the same format and get a batch output.  Script only beta version and needs 
  some more tweaks when I have time.

* <b>vcfExtractor.pl</b>
  Script to parse Ion Torrent specific VCF files and pull out variant data.  This works with TS v4.2 and v5.0
  files.  Some parts are compatible with Ion Reporter generated VCF files.

* <b>readlength_histogram.pl</b>
  Read in an Ion Torrent BAM file and generate a readlength histogram plot from the sample.  This script will require
  the `Statistics::R` perl module, as well as, the `ggplot2` library in R.  You can find [Statistics::R here](http://search.cpan.org/~fangly/Statistics-R/lib/Statistics/R.pm).
  You can find the  most excellent [ggplot2 package here](http://ggplot2.org).

* <b>map_ref.py</b>
    Map coordinates of two reference assemblies (e.g. hg18 and hg19) together in order.  This utility requires Konstantin's excellent python pyliftover library ([github](https://github.com/konstantint/pyliftover)), which leverages the UCSC liftOver utility for mapping reference assemblies. 

* <b>germline_merge.pl</b>
    Combine OCA Ion Reporter blood and tumor VCFs to generate a tumor / normal comparison file.  Quite academic for now, but may have more features as the need grows.

* <b>get_clinvar_variant_data.py</b>
    Input a file, comma separated list, or a single ClinVar ID and get a table of variant information derived from ClinVar using the eutils API functionality of NCBI. No filtering possible for now, but will be added later.  Requires: Python3
