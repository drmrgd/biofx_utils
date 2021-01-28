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

version <- '2.5.012821'

valid_terms <- c('refseq_mrna','ensembl_transcript_id', 'refseq_mrna', 
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
    metavar='<STR> Search Term',
    nargs='?',
    help="RefSeq ID, or comma separated list of IDs to search"
)
parser$add_argument(
    '-r', '--ref_version',
    metavar='<STR> Ref Version',
    default='grch37',
    choices=c("grch37", "grch38"),
    help=paste0("Version of human refernece to use. Input either 'grch37' ",
                "or 'grch38'. Default is 'grch37'")
)
parser$add_argument(
    "-f", "--file", 
    metavar='<FILE> Batchfile',
    help="Batch file of RefSeq IDs, one per line, to search"
)
parser$add_argument(
    "-F", "--Field", 
    default='refseq_mrna',
    metavar='<STR> Field',
    help=paste0("Field to use for searching. Default is refseq_mrna. Input '?' ",
                "to get a list of all valid terms.")
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


validate_terms <- function(term) {
    # Check to see that the query and wanted terms are valid, and if not, of if
    # we input '?', return a list of valid terms.
    # if (term != '?' && ! term %in% valid_terms) {
    if (term == '?') {
        message("Valid terms:")
        message(sprintf("\t> %s\n", valid_terms))
        quit()
    } 
    else if (! term %in% valid_terms) {
        message(sprintf("ERROR: term '%s' not valid.", term))
        validate_terms('?')
    } else {
        invisible(return)
    }
}


# Validate the query field, or output a list of valid terms if we're not sure.
validate_terms(args$Field)

wanted <- c("ensembl_transcript_id", "refseq_mrna", "entrezgene_id", 
            "ensembl_gene_id", "hgnc_symbol", "chromosome_name",
            "start_position", "end_position", "strand", "ensembl_exon_id")

# Set up and validated the desired output terms.
if (! is.null(args$wanted)) {
    want_list <- unlist(strsplit(args$wanted, ','))
    for (elem in want_list) {
        validate_terms(elem)
    }
    wanted <- want_list
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

# Create ensembl connection.
if (args$ref_version == 'grch38') {
    message("[INFO]: Using GRCh38 as the reference.\n")
    ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
} else {
    message("[INFO]: Using GRCh37 as the reference.\n")
    ensembl <- useMart("ensembl", host="grch37.ensembl.org", 
        dataset="hsapiens_gene_ensembl")
}

results <- getBM(attributes=wanted, filters=args$Field, values=refseq_queries,
    mart=ensembl)

# Print the output to either stdout or a file.
if (! is.null(args$outfile)) {
    print(paste0("Writing results to '", args$outfile, "'."))
    write.csv(results, file=args$outfile, quote=FALSE, row.names=FALSE)
} else {
    (results)
}
