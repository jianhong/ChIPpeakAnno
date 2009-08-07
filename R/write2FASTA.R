write2FASTA <-
function(mySeq, file="", width=80)
{
	descriptions=rownames(mySeq)
	sequences=mySeq$sequence
	numF = length(descriptions)
	ff = lapply(seq_len(numF),function(i)
	{
		desc <- descriptions[i]
        seq <- sequences[i]
        list(desc = desc, seq = seq)
    })
	writeFASTA(ff, file=file, width=width)
}