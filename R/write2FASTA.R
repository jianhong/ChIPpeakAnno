write2FASTA <-
function(mySeq, file="", width=80)
{
	descriptions=rownames(mySeq)
	sequences=mySeq$sequence
	names(sequences) = descriptions
	writeXStringSet(as(sequences,"XStringSet"), file, width=width)
}
