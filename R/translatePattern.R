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