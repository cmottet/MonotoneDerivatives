remove(list = ls())

library(DistributionPty)
library(ggplot2)
load("data/Hill_Horror_dist.Rdata")


Tsample <- log(sample)
p <- c(0.9,0.99,0.999,0.9999)
q <-  qlhorror(p)
x <- seq(min(Tsample),8.01,by = 0.01)

Data  <- data.frame(X = Tsample)
Data2 <- data.frame(X = x, Y = dlhorror(x))

marks <- data.frame(X = q)

ggplot(Data,aes(x= X))  + 
  geom_histogram(aes(y = ..density..),colour = "black",fill = "white",binwidth = 0.1) +
  geom_segment(data = marks,aes(x = X,y = 0,xend = X, yend = 0.4),linetype = "longdash")+
  labs(list(x = "Log(X) -  Sample Values", y = "", title = "")) +
  geom_text(data = marks,aes(x = X, y =0.3, 
                             label = c("90th percentile","99th percentile","99.90th percentile","99.99th percentile"),
                             angle = 90,vjust = 1.5),cex = 3)+
  ylim(0,0.6) + 
  geom_line(data = Data2,aes(x = X, y = Y,label = "True"),colour = "black",size = 1.2)+
  theme(axis.title.x = element_text(size = rel(0.8)))

ggsave("pics/Transformed_Histogram_Horror_Plot.png",width =5,height = 5) # Save plot in png

