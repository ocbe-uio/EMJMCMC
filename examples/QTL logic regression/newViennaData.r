data = read.csv("C:/Users/Aliaksandr/Downloads/Genotypingdata-EUmlauflogicregression.csv",colClasses = rep('factor', 19))

dat = as.data.frame(model.matrix( ~ 0+., data[,-1]))[,-2]
names(dat)
names(data)

write.csv(x = dat,file = "C:/Users/Aliaksandr/Downloads/qtlalzg.csv")