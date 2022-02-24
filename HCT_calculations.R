#Set the folder containing the results as working directory
setwd("~/")

#Enter the path to the generated file "fimo.tsv"
result <- "~/"

#enter the path of the folder conatining scanned results
list <- list.files("",pattern = ".tsv")

for(i in 1:length(list))
{
  #Store the path of the working directory as add
  scan <- read.table(file = list[i],header = TRUE,stringsAsFactors = FALSE,fill = TRUE)
  scan <- scan[order(scan$sequence_name),]
  tally <- as.data.frame(matrix(0,ncol = length(unique(scan$motif_alt_id)),nrow = length(unique(scan$motif_id))))
  row.names(tally) <- unique(scan$motif_id)
  colnames(tally) <- unique(scan$motif_alt_id)
  for(j in 1:ncol(tally))
  {
    print(paste(i,j,sep="_"))
    sub <- scan[scan$motif_alt_id %in% colnames(tally)[j],]
    tf <- unique(sub$motif_id)
    for(k in 1:length(tf))
    {
      p <- which(row.names(tally) == tf[k])
      sub_tf <- sub[sub$motif_id %in% tf[k],]
      scan_pos <- sub_tf[sub_tf$stop == "+",]
      scan_neg <- sub_tf[sub_tf$stop == "-",]
      s1 <- 1
      pos <- 1
      while(pos <= nrow(scan_pos))
      {
        a <- scan_pos$start[pos]+5
        a <- which(scan_pos$sequence_name > a)
        if(length(a) > 0)
        {
          pos <- a[1]
          s1 <- s1+1
        }else{
          pos <- nrow(scan_pos)+1
        }
      }
      neg <- 1
      s2 <- 1
      while(neg <= nrow(scan_neg))
      {
        a <- scan_neg$start[neg]+5
        a <- which(scan_neg$sequence_name > a)
        if(length(a) > 0)
        {
          neg <- a[1]
          s2 <- s2+1
        }else{
          neg <- nrow(scan_neg)+1
        }
      }
      tally[p,j] <- pmax(s1,s2)
    }
  }
  #Store the path of the folder that will contain results as result
  table <- as.data.frame(as.table(as.matrix(tally)))
  table <- table[table$Freq > 0,]
  colnames(table) <- c("TFs","Genes","Sites")
  result <- paste(result,list[i],sep = "")
  write.table(table, file = result,row.names = FALSE,quote = FALSE)
}

    
          
        
        
        
        
        
        
        


