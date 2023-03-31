#' getCDSfromGTF
#'
#' Get the full spliced coding sequences for each transcript
#'
#' @param gtf the path to the gtf file
#' @param seqs the sequences that are associated with the gtf file
#' @param attrsep in the attributes column of gtf, how are attributes separated? Default "; ".
#'
#' @export
#'
getCDSfromGTF <- function (gtf, seqs, attrsep = "; ")
{
  # check gtf file, import or use imported data frame
  cat("Getting Gtf file \n")
  gtfClasses = c("character", "character", "character", "integer",
                 "integer", "character", "character", "character", "character")
  if (is.character(gtf)) {
    gtf_dat = read.table(gtf, sep = "\t", as.is = TRUE, quote = "",
                         header = FALSE, comment.char = "#", nrows = -1, colClasses = gtfClasses)
  }
  else if (is.data.frame(gtf)) {
    stopifnot(ncol(gtf) == 9)
    if (!all(unlist(lapply(gtf, class)) == gtfClasses)) {
      stop("one or more columns of gtf have the wrong class")
    }
    gtf_dat <- gtf
    rm(gtf)
  }
  else {
    stop("gtf must be a file path or a data frame")
  }
  colnames(gtf_dat) = c("seqname", "source", "feature", "start",
                        "end", "score", "strand", "frame", "attributes")

  # get Coding sequences
  gtf_dat = gtf_dat[gtf_dat[, 3] == "CDS", ]

  # does any of the CDS gtf entries have non NA start and end
  stopifnot(!any(is.na(gtf_dat$start)), !any(is.na(gtf_dat$end)))

  # check fasta file given the unique sequences in gtf file
  # check if all chromosomes are matching
  cat("Checking Fasta file \n")
  chrs = unique(gtf_dat$seqname)
  if (is.character(seqs)) {
    fafiles = list.files(seqs)
    lookingFor = paste0(chrs, ".fa")
  }
  else {
    fafiles = names(seqs)
    lookingFor = chrs
  }
  if (!(all(lookingFor %in% fafiles))) {
    stop("all chromosomes in gtf must have corresponding sequences in seqs")
  }

  # get merge cds per each transcript
  cat("Getting CDS sequences \n")
  seqlist <- lapply(chrs, function(chr) {

    cat("Processing sequence:", chr, "\n")

    # get the current chromosome from the gtf and fasta files
    dftmp = gtf_dat[gtf_dat[, 1] == chr, ]
    fullseq = seqs[which(names(seqs) == chr)]

    # get cds reads from the start and end given in gtf file
    these_seqs = subseq(rep(fullseq, times = nrow(dftmp)), start = dftmp$start, end = dftmp$end)

    # get transcript names
    transcripts = getAttributeField(dftmp$attributes, field = "transcript_id", attrsep = attrsep)
    if (substr(transcripts[1], 1, 1) == "\"") {
      x <- transcripts
      transcripts <- substr(x, 2, nchar(x) - 1)
    }

    # aggregate by transcripts
    unique_transcripts <- unique(transcripts)
    names(these_seqs) <- transcripts
    new_seqs <- sapply(unique_transcripts, function(trans){
      temp <- these_seqs[transcripts == trans]
      DNAStringSet(unlist(temp))
    })
    names(new_seqs) <- NULL
    new_seqs <- do.call(c, new_seqs)

    # reverse complement sequences, if negative strand
    transcripts_strand <- data.frame(transcripts = transcripts, strand = dftmp$strand)
    transcripts_strand <- transcripts_strand[!duplicated(transcripts_strand),]
    revstrand = which(transcripts_strand$strand == "-")
    new_seqs[revstrand] = reverseComplement(new_seqs[revstrand])

    # return sequences
    names(new_seqs) <- transcripts_strand$transcripts
    new_seqs
  })

  # return sequences
  do.call(c, seqlist)
}

#' scanSeqWindow
#'
#' Scan a set of sequences using sliding windows and nucleotide pattern
#'
#' @param sequences a DNAStringSet object
#' @param window_size the size/length of the sliding windows
#' @param pattern the nucleotide pattern
#'
#' @import data.table
#' @import Biostrings
#'
#' @export
#'
scanSeqWindow <- function (sequences, window_size = 20, pattern)
{
  if(is.null(pattern) | is.null(sequences))
    stop("You have to provide both the sequences and the searched pattern")

  if(class(sequences) != "DNAStringSet")
    stop("You have to provide both the sequences and the searched pattern")

  if(nchar(pattern) > window_size)
    stop("The pattern cannot be longer than the window size =", window_size)

  # search for the pattern and count
  seq_scan <- lapply(seq_along(1:length(sequences)), function(i){
    cat("Sequence:", names(sequences)[i], "\n")
    seq_matched <- matchPattern(pattern, sequences[[i]])
    matched <- rep(0,length(sequences[[i]]))
    matched[start(seq_matched)] <- 1
    count_pattern <- frollsum(matched, window_size, align = "center")
    start_count_pattern <- 1:length(count_pattern)
    data.frame(chr = names(sequences)[i], start = start_count_pattern, end = start_count_pattern, score = count_pattern/window_size)
  })

  # return
  do.call(rbind,seq_scan)
}
