#!/usr/bin/env Rscript
# Quickie script to get some gene level information based on some inputs.
# Can get more Filters info with a `listFilters(ensembl)` call, or more 
# attributes with a `listAttributes(ensembl)` call.  If you don't want to use
# the `hsapiens_gene_ensembl` dataset, can find others with 
# `listDatasets(useMart("ensembl"))` call.
# 
# 7/25/2019 - D sims
################################################################################
library(argparse)
library(biomaRt)

version <- '2.0.012721'

valid_terms <- c('refseq_mrna','ensemble_transcript_id', 'refseq_mrna', 
          'entrezgene_id', 'hgnc_symbol', 'ensembl_gene_id')

parser <- ArgumentParser(
    description=paste0('Input a refseq ID and output some other ids and terms ',
        'relevant to the gene of interest. One can also use the "--Field" ',
        'argument in the event a refseq ID is not available, to search using ',
        'that term instead (i.e. ensembl_transcript_id). Available terms are ',
        '"ensembl_transcript_id", "refseq_mrna", "entrezgene_id", or ',
        '"hgnc_symbol".')
)
parser$add_argument(
    "refseq", 
    metavar='<STR> term to lookup (e.g. RefSeq ID).',
    nargs='?',
    help="RefSeq ID, or comma separated list of IDs to search"
)
parser$add_argument(
    "-f", "--file", 
    metavar='<FILE> File of terms to search.',
    help="Batch file of RefSeq IDs, one per line, to search"
)
parser$add_argument(
    "-F", "--Field", 
    metavar='<STR> Field',
    choices=valid_terms,
    help="Field to use for searching. Default is refseq_mrna."
)
parser$add_argument(
    "-w", "--wanted",
    metavar='<STR> Wanted Terms',
    help=paste0('In the event that only a subset of the terms are desired in the ',
        'output (like in the event we just want to map an ensemble_gene_id ',
        'to a hgnc_symbol), limit the output to these terms. If more than one ',
        'term is desired, enter each separated by a comma.')
)
parser$add_argument(
    "-o", "--outfile",
    metavar='<FILE> Output File',
    help='Write results to file instead of stdout.'
)
parser$add_argument(
    "-v", "--version", 
    action="version", 
    version=version,
    help="Print the version information and exit."
)
args <- parser$parse_args()

if (is.null(args$Field)) {
    search_term <- 'refseq_mrna'
} else {
    search_term <- args$Field
}

if (is.null(args$refseq) && is.null(args$file)) {
    stop("You must input either a RefSeq ID or a batch file of RefSeq IDs")
}

if (!is.null(args$refseq)) {
    refseq_queries <- strsplit(args$refseq, ',')
} else {
    # We are processing a batchfile
    refseq_queries <- scan(args$file, character(), quote="")
}

ensembl <- useMart("ensembl", host="grch37.ensembl.org", 
    dataset="hsapiens_gene_ensembl")
wanted <- c("ensembl_transcript_id", "refseq_mrna", "entrezgene_id", 
            "ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position",
            "end_position", "strand", "ensembl_exon_id")

if (! is.null(args$wanted)) {
    want_list <- unlist(strsplit(args$wanted, ','))
    for (elem in want_list) {
        if (! elem %in% valid_terms) {
            message(sprintf("ERROR: term '%s' not valid. Choose from:", elem))
            message(sprintf("\t> %s\n", valid_terms))
            stop()
        }
    }
    wanted <- want_list
}

results <- getBM(attributes=wanted, filters=search_term, values=refseq_queries,
    mart=ensembl)

if (! is.null(args$outfile)) {
    print("Writing results to 'results.txt'")
    write.csv(results, file="results.txt", quote=FALSE, row.names=FALSE)
} else {
    (results)
}
