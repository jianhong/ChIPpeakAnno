#' translate pattern from IUPAC Extended Genetic Alphabet to regular expression
#' 
#' translate pattern containing the IUPAC nucleotide ambiguity codes to regular
#' expression.  For example,Y->[C|T], R-> [A|G], S-> [G|C], W-> [A|T], K->
#' [T|U|G], M-> [A|C], B-> [C|G|T], D-> [A|G|T], H-> [A|C|T], V-> [A|C|G] and
#' N-> [A|C|T|G].
#' 
#' 
#' @param pattern a character vector with the IUPAC nucleotide ambiguity codes
#' @return a character vector with the pattern represented as regular
#' expression
#' @author Lihua Julie Zhu
#' @seealso countPatternInSeqs, summarizePatternInPeaks
#' @keywords misc
#' @export
#' @examples
#' 
#'     pattern1 = "AACCNWMK"
#'     translatePattern(pattern1)
#' 
translatePattern <- function(pattern)
{
    pattern = toupper(pattern)
    pattern = gsub("Y","[C|T]", pattern)
    pattern = gsub("R", "[A|G]", pattern)
    pattern = gsub("S", "[G|C]", pattern)
    pattern = gsub("W", "[A|T]", pattern)
    pattern = gsub("K", "[T|U|G]", pattern)
    pattern = gsub("M", "[A|C]", pattern)
    pattern = gsub("B", "[C|G|T]", pattern)
    pattern = gsub("D", "[A|G|T]", pattern)
    pattern = gsub("H", "[A|C|T]", pattern)
    pattern = gsub("V", "[A|C|G]", pattern)
    pattern = gsub("N", "[A|C|T|G]", pattern)	
    pattern
}
