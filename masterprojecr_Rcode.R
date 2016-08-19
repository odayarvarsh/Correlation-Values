###Reading the two table datas
groupa=read.table("GroupA_BRCA_miR.txt",header=T,row.names=1,stringsAsFactors=F)
trimmed_data=read.table("Trimmed_BRCA_mRNA.txt",header=T,row.names=1,stringsAsFactors=F)
dim(groupa)
dim(trimmed_data)

#Trim the data to 25 columns as opposed to 49 
common_groupa=find.common.rna(groupa,trimmed_data)
common_groupamir=common_groupa[1:dim(groupa)[1],]
common_groupamr=common_groupa[-(1:dim(groupa)[1]),]

######find.cross.cor####

#generalized function to find the correlation between miRNA's of each rows 
#use k+1 function in for loop to make sure that same row is not being correlated with itself
#otherwise, all correlation values will be 1, and that cannot happen as that is perefect correlation-this would mean that all miRNA's in each row are the same 


groupavec=vector()

for (i in 1:nrow(common_groupamir))
{
  for (k in i:nrow(common_groupamir))
  {
    groupavec=c(groupavec,cor(as.numeric(common_groupamir[i,]), as.numeric(common_groupamir[k+1,])))
  } 
}
find.cross.cor=function(common_groupamir){groupavec=vector()
rownames_i=vector()
rownames_k=vector()
for (i in 1:(nrow(common_groupamir)-1))
{
  for (k in (i+1):nrow(common_groupamir))
  {
    groupavec=c(groupavec,cor(as.numeric(common_groupamir[i,]), as.numeric(common_groupamir[k,])))
    rownames_i=c(rownames_i, (row.names(common_groupamir[i,])))
    rownames_k=c(rownames_k, (row.names(common_groupamir[k,])))
  } 
}
crosscor=rbind(rbind(groupavec,rownames_i),rownames_k)
index=order(crosscor[1,],decreasing=TRUE)
vecc=crosscor[,index]
top_20=vecc[,1:20]
return(t(top_20))
}

############

crosscor=find.cross.cor(common_groupamir)
