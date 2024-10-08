#' Function to use BLAST to align ASV sequences across multiple datasets.
#'
#' Note that currently no attempt is made to resolve dependencies for this function so you will need to
#' have rBLAST from github (\code{devtools::install_github("mhahsler/rBLAST")}) and Biostrings/seqinr
#' already installed. I'll get to at some point.
#'
#' @param seqs1 A named list of ASV sequences. Generally the names should be ASV numbers.
#' @param seqs2 An optional second named list of ASV sequences. If db is left NULL then this is used
#' to make a blast database and seqs1 will be blasted against it before it is removed.
#' @param cutoff A percentage identity required to return a match. Defaults to 97.5
#' @param db An optional database to use. Note if this is supplied then seqs2 is ignored.
#' @param maxMatches An optional number of maximum matches to return per sequence in seqs1.
#' Defaults to Inf, which will return all matches above the cutoff.
#' @param cores Number of cores optionally to run in parallel. Defaults to 1 if "mc.cores" is not set.
#' @keywords BLAST
#' @import parallel
#' @return A dataframe of blast results
#'
#' @examples
#'
#' if ("rBLAST" %in% installed.packages() && rBLAST::has_blast()) {
#'   seqs1 <- list(
#'     ASV9 = paste0(
#'       "GTGCCAGCAGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGT",
#'       "AAAGGGTGCGTAGGCGGGTCTTTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCT",
#'       "TTGATACTGAAGATCTTGAGTTCGGGAGAGGTGAGTGGAACTGCGAGTGTAGAGGTGAAATTCGTA",
#'       "GATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGATACTGACGCTGAGGCACGAAAGCGT",
#'       "GGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAATGCCAGCCGTTAGTGGGTT",
#'       "TACTCACTAGTGGCGCAGCTAACGCTTTAAGCATTCCGCCTGGGGAGTACGGTCGCAAGATTAAAACTCAAATGAATTGACGG"
#'     ),
#'     ASV10 = paste0(
#'       "GTGCCAGCCGCCGCGGTAATACGAAGGGGGCTAGCGTTGCTCGGAATCACTGGGCGTAAAGGGTGCGTAGGCGGGTCT",
#'       "TTAAGTCAGGGGTGAAATCCTGGAGCTCAACTCCAGAACTGCCTTTGATACTGAAGATCTTGAGTTCGGGAGAGGTGAGT",
#'       "GGAACTGCGAGTGTAGAGGTGAAATTCGTAGATATTCGCAAGAACACCAGTGGCGAAGGCGGCTCACTGGCCCGATACT",
#'       "GACGCTGAGGCACGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAATGCCAGCC",
#'       "GTTAGTGGGTTTACTCACTAGTGGCGCAGCTAACGCTTTAAGCATTCCGCCTGGGGAGTACGGTCGCAAGATTAAAACT",
#'       "CAAATGAATTGACGG"
#'     )
#'   )
#'   seqs2 <- list(
#'     ASV1 = paste0(
#'       "GTGCCAGCAGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAAGT",
#'       "TGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATTTCCT",
#'       "GTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGAGGTGC",
#'       "GAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCAACTAGCCGTTGGGAGCCTTG",
#'       "AGCTCTTAGTGGCGCAGCTAACGCATTAAGTTGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGG"
#'     ),
#'     ASV2 = paste0(
#'       "GTGCCAGCCGCCGCGGTAATACAGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCGCGTAGGTGGTTTGTTAA",
#'       "GTTGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCAAAACTGACAAGCTAGAGTATGGTAGAGGGTGGTGGAATT",
#'       "TCCTGTGTAGCGGTGAAATGCGTAGATATAGGAAGGAACACCAGTGGCGAAGGCGACCACCTGGACTGATACTGACACTGA",
#'       "GGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCAACTAGCCGTTGGGAG",
#'       "CCTTGAGCTCTTAGTGGCGCAGCTAACGCATTAAGTTGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAA",
#'       "TTGACGG"
#'     )
#'   )
#'   blastASVs(seqs1, seqs2, cutoff = 0)
#' }
#'
#' @export


blastASVs <- function(seqs1, seqs2, cutoff = 97.5, db = NULL, maxMatches = Inf,
                      cores = getOption("mc.cores", 1)) {
  delete_db <- FALSE
  if (is.null(db)) {
    tmp_file <- tempdir()
    timestamp <- round(as.numeric(Sys.time()))
    db <- paste0(tmp_file, "/db_", timestamp, "/db")
    seqinr::write.fasta(seqs2, names(seqs2),
      file.out = paste0(tmp_file, "/", timestamp, ".fasta")
    )
    rBLAST::makeblastdb(
      file = paste0(tmp_file, "/", timestamp, ".fasta"),
      db_name = db, dbtype = "nucl"
    )
    delete_db <- TRUE
  }
  db_obj <- rBLAST::blast(db)
  res <- do.call(rbind, parallel::mclapply(seq_along(seqs1), function(i) {
    sq <- seqs1[[i]]
    xsq <- Biostrings::DNAStringSet(sq)
    res <- predict(db_obj, xsq, BLAST_args = paste0("-perc_identity ", cutoff))
    colnames(res)[1:2] <- c("seq1_id", "db_id")
    if (nrow(res) >= 1) {
      res$seq1_id <- names(seqs1)[i]
      res$match_rank <- nrow(res) + 1 - rank(res$pident, ties.method = "average")
    } else {
      res <- structure(
        list(
          seq1_id = names(seqs1)[i], db_id = "NONE", pident = 0,
          length = 0L, mismatch = 0L, gapopen = 0L, qstart = 1L,
          qend = 0L, sstart = 0L, send = 0L, evalue = 0, bitscore = 0L,
          match_rank = 1
        ),
        row.names = 1L, class = "data.frame"
      )
    }
    res <- res[as.numeric(rownames(res)) <= maxMatches, ]
    return(res)
  }, mc.cores = cores))
  if (delete_db) {
    unlink(paste0(tmp_file, "/", timestamp, ".fasta"))
    unlink(paste0(tmp_file, "/db_", timestamp), recursive = TRUE)
  }
  return(res)
}
