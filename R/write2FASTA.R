write2FASTA <-
function(mySeq, file="", width=80)
{
	descriptions=rownames(mySeq)
	sequences=mySeq$sequence
	names(sequences) = descriptions
	#numF = length(descriptions)
	
	#ff = lapply(seq_len(numF),function(i)
	{
	#	desc <- descriptions[i]
        #seq <- sequences[i]
        #list(desc = desc, seq = seq)
    #})
	write.XStringSet(as(sequences,"XStringSet"), file, width=width)
}
