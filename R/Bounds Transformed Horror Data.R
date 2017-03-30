remove(list = ls())


library(latex2exp)
library(ggplot2)
library(magrittr)
load("data/Hill_Horror_dist.Rdata")
load("data/MinMaxTransformedHorrorDist_TailProb.Rdata")


tmp <- which(apply(optim[,13:17],1,function(x)sum(!is.na(x))>0 ))
length(which(apply(optim[,13:17],1,function(x)sum(!is.na(x))>0 )))
Data <- optim[tmp,]

J1 <- m
J2 <- d
NJ1 <- length(J1)+1 # To account for 0
NJ2 <- length(J2)-1
N <- NJ1 + NJ2
ParamNames <- names(Data)[1:N]


Data <- transform(Data,
                 nparam = as.factor(rowSums(Data[,1:N])),
                 J1empty = rowSums(Data[,(NJ2+1):N])==0,
                 J2empty = rowSums(Data[,1:NJ1])==0,
                 RelErr = with(Data,apply(direction*(cbind(lB,uB) - (1-P))/(1-P),1 ,min)),
                 D  = as.factor(D),
                 Est = rowSums(Data[,1:NJ2]) >= 1,
                 param.inc = apply(Data[,1:N], 1,function(mat){paste(subset(ParamNames,unlist(mat)),collapse = ",")})
)

maxErrRel <- expand.grid(D = 0:4, P = P) %>% transform(maxErrRel = P/(1-P))

ggplot(data = subset(Data,direction ==1), aes(x = D, y = RelErr,text = param.inc,size = J1empty)) +
  geom_jitter(width = 0.3,color = alpha("black",1/5)) + 
  geom_hline(data = maxErrRel, aes(yintercept = maxErrRel))+
  facet_wrap(~P,ncol = 4,nrow = 4,
             labeller = as_labeller(c("0.9" = "p = 90",
                                      "0.99" = " p = 99",
                                      "0.999" = " p = 99.9", 
                                      "0.9999" = " p = 99.99")))   + 
  coord_trans(y = "log10", limy = c(0.24,1.5e4)) + 
  scale_y_continuous(breaks =   c(0.25,1,2.5,10,25,1e2,250,1e3,2500,1e4), 
                     labels = c(0.25,1,2.5,10,25,1e2,250,1e3,2500,1e4)) + 
  scale_x_discrete(breaks =   seq(0,5,by = 2)) +
  labs(x = TeX('$D$'), y = TeX("Relative error"), size = TeX('Subset $\\mathit{J_1}$')) + 
  scale_size_discrete(labels = c("Not empty","Empty"))

ggsave("pics/Bounds_Transformed_Horror_Data.png",
       width =10,
       height = 5) # Save plot in png

