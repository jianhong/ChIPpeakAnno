getUniqueGOidCount <-
function(goList)
{
	x = goList[order(goList)]
	duplicated.go<-duplicated(x)
	unique.go <- unique(x)
	go.count<- numeric()

	count <- 1
	j <- 1

	for (i in 2:length(duplicated.go))
	{
		
		if (duplicated.go[i] == FALSE)
		{
			go.count[j] <- count # previous GO
			j <- j + 1
			count <- 1
		}
		else
		{
			count <- count + 1
		}
	}
	go.count[j] <- count #last GO
	list(GOterm=unique.go, GOcount=go.count)
}

