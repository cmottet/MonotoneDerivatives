source("init.R")
input.file1  = file.path(input_dir,"data/Hill_Horror_dist.Rdata")
input.file2  = file.path(output_dir,"data/TransformedHorrorDist_TailProb.Rdata")
load(input.file1);load(input.file2)

NorderM = length(orderM)
NorderL = length(orderL)
NmaxParam = NorderL + NorderM
Data = NULL

for (j in 1:length(p))
  Data = rbind(Data, transform(optim[[j]],
                               p_index = j,
                               RelErr  = (primal_uB - (1-p[j]))/(1-p[j]),
                               nparam  = as.factor(rowSums(optim[[j]][,1:NmaxParam])),
                               emptyM  = rowSums(optim[[j]][,NorderL:NmaxParam])==0 
  ))

ddply(Data,.(p_index), subset, min(RelErr) == RelErr, select = c(1:8,13,21))
ddply(Data,.(p_index,nparam), subset, min(RelErr) == RelErr, select = c(1:8,13,21))

FindJandL_byp = ddply(Data,.(p_index), function(Data) with(Data,which(min(RelErr) == RelErr)))
FindJandL_bypN = ddply(Data,.(p_index,nparam), function(Data) with(Data,which(min(RelErr) == RelErr)))

Data[FindJandL_bypN$V1,]
FindJandL = ddply(Data,.(p_index), function(Data) with(Data,which(min(RelErr) == RelErr)))


