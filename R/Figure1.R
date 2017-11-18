remove(list = ls())

library(ggplot2)
library(DistributionPty)
load("data/log_Hill_Horror_dist.Rdata")

p <- c(0.9,0.99,0.999,0.9999)
q <-  qlhorror(p)
x <- seq(qlhorror(0),7.5,by = 0.1)

Data <- data.frame(X = sample)
Data2 <- data.frame(X = x, Y = dlhorror(x))

marks <- data.frame(X = q)

plot <- ggplot(Data,aes(x= X))  +
  geom_histogram(aes(y = ..density..),colour = "black",fill = "white",binwidth = 0.5) +
  geom_segment(data = marks,aes(x = X,y = 0,xend = X, yend = 0.3),linetype = "longdash")+
  labs(list(x = "Sample values", y = "", title = "")) +
  geom_text(data = marks,aes(x = X, y =0.2,
                             label = c("90th percentile","99th percentile",
                                       "99.9th percentile","99.99th percentile"),
                             angle = 90,vjust = 1.5),cex = 3)+
  ylim(0,0.34) +
  theme(axis.title.x = element_text(size = rel(0.8))) +
  geom_line(data = Data2,aes(x = X, y = Y,label = "True"),colour = "black",size = 1.1)+
  scale_x_continuous(breaks = q,labels = signif(q,2))

plot 

ggsave(plot, file = "pics/Figure1.pdf", width = 5,height = 5,dpi=300)
