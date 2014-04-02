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
