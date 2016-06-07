#tmp=read.table('data/data_exemple1.csv',sep=',',header=T) #from Epic Italy
#attach(tmp)
#model.fit=glm(Outcome~Mvalue,family=binomial(link='logit'))
#cof=summary(model.fit)$coefficients[,1]
#x=rnorm(100,mean=mean(Mvalue),sd=sd(Mvalue))
#cof[2]=-0.5 #modified to have a weaker pvalue
#y=rbinom(100,1,1/(1+exp(-cof[1]-cof[2]*x)))
#summary(glm(y~x, family=binomial(link='logit')))

y=c(rep(0,50),rep(1,50))
x=c(rnorm(50,2,0.7),rnorm(50,1.5,0.7))

example_methy=data.frame(Outcome=y,Mvalue=x)
write.table(example_methy,file='data/ex_methy.csv',sep=',',col.names=T,row.names=F,quote=F)

