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
    these_seqs <- subseq(rep(fullseq, times = nrow(dftmp)), start = dftmp$start, end = dftmp$end)

    # get transcript names
    transcripts = getAttributeField(dftmp$attributes, field = "transcript_id", attrsep = attrsep)
    if (substr(transcripts[1], 1, 1) == "\"") {
      x <- transcripts
      transcripts <- substr(x, 2, nchar(x) - 1)
    }

    # aggregate sequences and locations by transcripts
    unique_transcripts <- unique(transcripts)
    names(these_seqs) <- transcripts
    new_seqs <- sapply(unique_transcripts, function(trans){
      temp <- these_seqs[transcripts == trans]
      DNAStringSet(unlist(temp))
    })
    names(new_seqs) <- NULL
    new_seqs <- do.call(c, new_seqs)
    locations <- sapply(unique_transcripts, function(trans){
      dftmp_trans <- dftmp[transcripts == trans,]
      dftmp_trans_seq <- mapply(function(x,y) return(x:y), dftmp_trans$start, dftmp_trans$end)
      as.vector(unlist(dftmp_trans_seq))
    })

    # reverse complement sequences, if negative strand
    transcripts_strand <- data.frame(transcripts = transcripts, strand = dftmp$strand)
    transcripts_strand <- transcripts_strand[!duplicated(transcripts_strand),]
    revstrand <- which(transcripts_strand$strand == "-")
    new_seqs[revstrand] <- reverseComplement(new_seqs[revstrand])
    locations[revstrand] <- lapply(locations[revstrand], rev)
    locations <- lapply(locations, function(loc) {
      data.frame(locations = loc, entry = rep(chr, length(loc)))
    })

    # return sequences
    names(new_seqs) <- transcripts_strand$transcripts
    return(list(strings = new_seqs, locations = locations))
  })

  # get all strings
  strings <- do.call(c, lapply(seqlist, function(x) x$strings))

  # get and combine all locations with fasta entries
  locations <- do.call(c, lapply(seqlist, function(x) x$locations))

  # return sequences
  return(list(strings = strings, locations = locations))
}

#' scanCDSSeqWindow
#'
#' Scan a set of sequences using sliding windows and nucleotide pattern
#'
#' @param sequences a DNAStringSet object
#' @param
#' @param annotation the annotation file
#' @param window_size the size/length of the sliding windows
#' @param pattern the nucleotide pattern
#'
#' @import data.table
#' @import Biostrings
#'
#' @export
#'
scanCDSSeqWindow <- function (sequences, locations, annotation, window_size = 20, pattern)
{
  if(is.null(pattern) | is.null(sequences))
    stop("You have to provide both the sequences and the searched pattern")

  if(class(sequences) != "DNAStringSet")
    stop("You have to provide both the sequences and the searched pattern")

  if(nchar(pattern) > window_size)
    stop("The pattern cannot be longer than the window size =", window_size)

  # search for the pattern and count
  cat("Scanning sequences \n")
  break_seq <- floor(length(sequences)/50)
  seq_scan <- lapply(seq_along(1:length(sequences)), function(i){
    # cat("ID:", i, "Sequence:", names(sequences)[i], "\n")
    if(!(i %% break_seq)){
      prog <- floor(i/break_seq)
      cat('\r', "Working  [", strrep("#",prog),
          strrep(" ", 50-prog), "] ", 2*prog, "%", sep="")
    }
    seq <- sequences[[i]]
    loc <- locations[[i]]
    seq_matched <- matchPattern(pattern, seq)
    matched <- rep(0,length(seq))
    matched[start(seq_matched)] <- 1
    count_pattern <- frollsum(matched, window_size, align = "center")
    start_count_pattern <- loc$locations
    window_size_adj <- floor(window_size/nchar(pattern))
    data.frame(chr = loc$entry, start = start_count_pattern, end = start_count_pattern, score = count_pattern/window_size_adj)
  })
  cat(" \n\r")

  # combine
  cat("Combine Scans \n")
  seq_scan <- do.call(rbind,seq_scan)
  seq_scan <- na.omit(seq_scan)

  # sort for each chromosome
  cat("Sorting fasta entries \n")
  seq_scan_sorted <- NULL
  all_chr <- unique(seq_scan$chr)
  cat("Total number of fasta entries: ", length(all_chr), "\n")
  break_seq <- floor(length(all_chr)/10)
  for(i in 1:length(all_chr)){
    if(!(i %% break_seq)){
      prog <- floor(i/break_seq)
      cat('\r', "Working  [", strrep("#",prog),
          strrep(" ", 10-prog), "] ", 2*prog, "%", sep="")
    }
    temp <- seq_scan[seq_scan$chr == all_chr[i],]
    seq_scan_sorted <- rbind(seq_scan_sorted,
                             temp[order(seq_scan$start, decreasing = FALSE),])
  }

  # return
  seq_scan_sorted
}
