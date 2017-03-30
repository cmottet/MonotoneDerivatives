# install.packages("filehash")
# install.packages("tikzDevice")

# library(tikzDevice)
# library(systemfit)
# options(tikzMetricPackages = c("\\usepackage[utf8]{inputenc}", 
#                                "\\usepackage[T1]{fontenc}",
#                                "\\usetikzlibrary{calc}",
#                                "\\usepackage{amssymb}"))
# tikz("text.tex",width = 10,height = 5,standAlone = TRUE,
#      packages = c("\\usepackage{tikz}",
#                   "\\usepackage[active,tightpage,psfixbb]{preview}",
#                   "\\PreviewEnvironment{pgfpicture}",
#                   "\\setlength\\PreviewBorder{0pt}",
#                   "\\usepackage{amssymb}"))
# par(mar = c(3,5,3,5))

source("init.R")
library(latex2exp)

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
                               emptyM  = as.factor(rowSums(optim[[j]][,(NorderL):NmaxParam])==0) 
  ))

ggplot(data = Data, aes(x = nparam, y = RelErr, size = emptyM)) +
       geom_jitter(width = 0.3,color = alpha("black",1/3)) + 
       facet_wrap(~p_index,ncol = 4,nrow = 4,
                  labeller = as_labeller(c("1" = "p = 90","2" = " p = 99","3" = " p = 99.9", "4" = " p = 99.99"))) + 
       coord_trans(y = "log10", limy = c(0.24,1.5e4)) + 
       scale_y_continuous(breaks =   c(0.25,1,2.5,10,25,1e2,250,1e3,2500,1e4), 
                          labels = c(0.25,1,2.5,10,25,1e2,250,1e3,2500,1e4)) + 
       scale_x_discrete(breaks =   seq(0,NmaxParam,by = 2)) + 
       scale_size_discrete(labels = c("Not empty","Empty"))+
       labs(x = TeX('$N$'),
            y = TeX("Relative error"),
            size = TeX('Subset $\\mathit{J}$'))

output.file = file.path(out_pic,"Paper Plots/Bounds_Transformed_Horror_Data.png")
ggsave(output.file,width =10,height = 5) # Save plot in png
