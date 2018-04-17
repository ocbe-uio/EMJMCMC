library(ggplot2)
library(MASS)

S1 <- matrix(c(1,0,0,1),nrow=2,byrow=TRUE)
mu1 <- c(4,3)
S2 <- matrix(c(1,0,0,1),nrow=2,byrow=TRUE)
mu2 <- c(-1,-5)
S3 <- matrix(c(1,0,0,1),nrow=2,byrow=TRUE)
mu3 <- c(3,-2)
S4 <- matrix(c(1,0,0,1),nrow=2,byrow=TRUE)
mu4 <- c(-1,0)
S5 <- matrix(c(1,0,0,1),nrow=2,byrow=TRUE)
mu5 <- c(0,7)


n <- 1000
p1 <- 0.20
n1 <- rbinom(1,size=n,prob=p1)  ## how many from first distribution?
n2 <- rbinom(1,size=n,prob=p1)
n3 <- rbinom(1,size=n,prob=p1)
n4 <- rbinom(1,size=n,prob=p1)
n5 <- n - n1 - n2 - n3 - n4
val1 <- mvrnorm(n1,mu=mu1,Sigma=S1)
val2 <- mvrnorm(n2,mu=mu2,Sigma=S2)
val3 <- mvrnorm(n3,mu=mu3,Sigma=S3)
val4 <- mvrnorm(n4,mu=mu4,Sigma=S4)
val5 <- mvrnorm(n4,mu=mu5,Sigma=S5)
allval <- rbind(val1,val2,val3,val4,val5)      ## combine
allval <- allval[sample(n,n),]  ## scramble order

df = data.frame(allval); colnames(df) = c("x","y")

#write.csv(df,"x.csv")

commonTheme = list(labs(color="Density",fill="Density",
                        x="",
                        y=""),
                   theme_bw(),
                   theme(legend.position=c(0,1),
                         legend.justification=c(0,1)))

ggplot(data = df, aes(x=x, y=y) ) +
  stat_density2d(aes(fill=..level..,alpha=..level..),geom='polygon',colour='black')+
  scale_fill_continuous(low="green",high="red") +
  guides(alpha="none") +
  commonTheme  + 
  #annotate("rect", xmin = 3.80, xmax = 4.15, ymin = 2.75, ymax = 3.25,
  #                                                      colour = "gray", fill = "white")+ 
  annotate("text", x=3.95, y=3, label= paste(expression(chi[k])),parse = TRUE,cex = 10)  +
  #annotate("rect", xmin = 3.35, xmax = 3.70, ymin = 2.05, ymax = 2.55,
  #         colour = "gray", fill = "white")+ 
  annotate("text", x=3.20, y=4.30, label= paste(expression(gamma)),parse = TRUE,cex = 10)  +
  #annotate("rect", xmin = 5.50, xmax = 5.85, ymin = 1.75, ymax = 2.25,
  #         colour = "gray", fill = "white")+ 
  annotate("text", x=5.70, y=2, label= paste(expression(chi[0])),parse = TRUE,cex = 10)  +
  #annotate("rect", xmin = -1.20, xmax = -0.83, ymin = -0.25, ymax = 0.25,
  #         colour = "gray", fill = "white")+ 
  annotate("text", x=-1, y=0, label= paste(expression({{chi[k]}*{textstyle("*")}})),parse = TRUE,cex = 10)+
  #annotate("rect", xmin = -0.60, xmax = -0.25, ymin = -0.35, ymax = 0.15,
  #         colour = "gray", fill = "white")+
  annotate("text", x=0.00, y=-0.1, label= paste(expression(gamma^{""[""[""[""[textstyle("*")]]]]})),parse = TRUE,cex = 10)+
  #annotate("rect", xmin = -1.20, xmax = -0.83, ymin = 2.25, ymax = 2.75,
  #       colour = "gray", fill = "white")+ 
  annotate("text", x=-1, y=2.5, label= paste(expression({{chi[0]}*{textstyle("*")}})),parse = TRUE,cex = 10)+
  geom_segment(aes(x=0.05, y=-0.05, xend=5.50, yend=2.00), 
               arrow = arrow(), colour = "red")+
  geom_segment(aes(x=3.15, y=4.30, xend=-0.83, yend=2.40), 
               arrow = arrow(), colour = "red")+
  geom_segment(aes(x=5.50, y=2.15, xend=4.18, yend=2.90), 
               arrow = arrow(), colour = "blue")+
  geom_segment(aes(x=3.75, y=3.00, xend=3.30, yend=4.30), 
               arrow = arrow(), colour = "green")+
  geom_segment(aes(x=-1.00, y=1.95, xend=-1.00, yend=0.25), 
               arrow = arrow(), colour = "blue")+
  geom_segment(aes(x=-0.80, y=0.00, xend=-0.20, yend=-0.25), 
               arrow = arrow(), colour = "green")+theme(axis.text.x=element_blank(), axis.text.y=element_blank())

