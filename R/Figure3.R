library(latex2exp)
library(ggplot2)
library(magrittr)

load("data/log_Hill_data_optimparam.Rdata")
load("data/MaxLogdHorrorDist_TailProb.Rdata")

optim <- optim[-which(apply(optim, 1, function(x)all(is.na(x)))),]
J1 <- 0:3
J2 <- 3:1
NJ1 <- length(J1)+1 # To account for 0
NJ2 <- length(J2)-1
N <- NJ1 + NJ2
ParamNames <- names(optim)[1:N]

Data <- transform(optim,
                  nparam = as.factor(rowSums(optim[,1:N])),
                  J1empty = rowSums(optim[,(NJ2+1):N])==0,
                  J2empty = rowSums(optim[,1:NJ1])==0,
                  RelErr = with(optim,apply((cbind(lB) - (1-P))/(1-P),1 ,min)),
                  D  = as.factor(D),
                  Est = rowSums(optim[,1:NJ2]) >= 1,
                  param.inc = apply(optim[,1:N], 1,function(mat){paste(subset(ParamNames,unlist(mat)),collapse = ",")})
)

###
### Plot Figure 3
###

convexData <- subset(Data, (D==2) & (!d3) & (d2) & d1 & m0 & (!m1) & (!m2) & (!m3) )
maxErrRel <- expand.grid(D = 0:5, P = P) %>% transform(maxErrRel = P/(1-P))

plot <- ggplot(data = Data, aes(x = D, y = RelErr,text = param.inc,size = J1empty)) +
  geom_jitter(width = 0.3,color = alpha("black",1/5)) +
  geom_hline(data = maxErrRel, aes(yintercept = maxErrRel))+
#  geom_point(data = convexData, aes(x = D, y = RelErr, color="Convex Approach") ) +
  facet_wrap(~P,ncol = 4,nrow = 4,
             labeller = as_labeller(c("0.9" = "p = 90",
                                      "0.99" = " p = 99",
                                      "0.999" = " p = 99.9",
                                      "0.9999" = " p = 99.99")))   +
  coord_trans(y = "log10", limy = c(0.01,1.5e4)) +
  scale_y_continuous(breaks =   c(0.25,1,2.5,10,25,1e2,250,1e3,2500,1e4),
                     labels = c(0.25,1,2.5,10,25,1e2,250,1e3,2500,1e4)) +
  scale_x_discrete(breaks =   seq(0,5,by = 2)) +
  labs(x = TeX('$D$'), y = TeX("Relative error"), size = TeX('Subset $\\mathit{J_1}$')) +
  labs(x = TeX('$D$'), y = TeX("Relative error"), size = "", color="") +
  scale_size_discrete(labels = c("At least one moment constraint","No moment constraint"))

plot

ggsave(plot, file = "pics/Figure3.pdf", width = 10,height = 5,dpi=300)

###
### Get best subsets J1* and J2*
###
library(plyr)
Table <- ddply(Data, .(P,D), .fun = function(x){
  index_min <- which.min(x$RelErr)
  row_min <- x[index_min, ]

  J2_star <- J2[unlist(row_min[1:length(J2)])] %>% toString()
  J1_star <- J1[unlist(row_min[length(J2) + 1:length(J1)])] %>% toString()

  if (length(J2_star) ==0) J2_star = ' '
  if (length(J1_star) ==0) J2_star = ' '

  Z_star <- format(min(row_min$uB,row_min$lB), scientific=TRUE, digits =3)
  data.frame(J1 = paste0("{",J1_star,"}"), J2 = paste0("{",J2_star,"}"), Z_star, row_min$RelErr)
})
Table$P <- Table$P*100
names(Table) <- c("P", "D", "J1", "J2","Optimal Objective Value", "Relative Error")

# Create the LaTex Tables
library(xtable)
library(dplyr)
display <- c(rep("d",4), "s","E")
for (p in P*100){
  file <- paste0("tables/Transformed_Hill_Horror_TailProb_best_",p,".tex")
  table <- subset(Table, p == P)%>% select(-P)
  align  <- c(rep("|c",ncol(table)+1),"|")

  xtable(table, align = align, digits = 3) %>% print(include.rownames=FALSE, type = "latex",file = file)

  # We need to replace the tables by subtables (our way around uses unix commands)
  command_1 <- paste0("sed -i '' -- 's/begin{table}\\[ht\\]/begin{subtable}{\\\\textwidth}/g' ",file)
  system(command_1)

  command_2 <- paste0("sed -i '' -- 's/end{table}/end{subtable}/g' ",file)
  system(command_2)

  # We need to add the subcaption to the subtables (our way around uses unix commands)
  command_3 <- paste0("sed -i '' -- 's/D \\& J1 \\& J2/\\$D\\$ \\& \\$\\\\mathcal J\\_1\\^\\*\\$ \\& \\$\\\\mathcal J\\_2\\^\\*\\$/g' ",file)
  system(command_3)

  # Make the right column names
  command_4 <- paste0("sed -i '' -- 's/end{tabular}/end{tabular}\\\\subcaption{$p = ",p," $}/g' ",file)
  system(command_4)

}
