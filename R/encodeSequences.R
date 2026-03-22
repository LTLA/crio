#' Encode and decode sequences
#'
#' Encode sequences to a 2-bit integer representation, and decode the integers back to sequences.
#'
#' @param sequences Character vector of nucleotide sequences.
#' @param encoded Integer vector of 2-bit encodings of the sequences. 
#' @param lengths Integer vector containing the lengths of sequences.
#' This should be the same length as \code{encoded}, otherwise it is recycled to the length of \code{encoded}.
#'
#' @return For \code{encodeSequences}, an integer vector containing the 2-bit integer encoding for each sequence.
#' 
#' For \code{decodeSequences}, a character vector containing the decoded sequences.
#'
#' @author Aaron Lun
#'
#' @examples
#' enc <- encodeSequences(c("AAAA", "CCCC", "TTTT", "GGGG"))
#' enc
#' decodeSequences(enc, lengths=4)
#'
#' @seealso
#' \code{\link{readMolecules}}, which uses a 2-bit encoding for the UMI sequences.
#'
#' @export 
encodeSequences <- function(sequences) {
    encode_sequences(sequences)
}

#' @export
#' @rdname encodeSequences
decodeSequences <- function(encoded, lengths) {
    decode_sequences(encoded, rep(lengths, length.out=length(encoded)))
}
