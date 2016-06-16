#############
# Example methylation
#############
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




#############
# Example independent SNPs
#############

require(waffect)

n=500
geno=matrix(NA,nrow=n,ncol=10)
maf=c(.2,.1,.3,.1,.1,.2,.4,.3,.2,.4)
for(s in 1:10) geno[,s]=sample(0:2,size=n,replace=TRUE,prob=c((1-maf[s])^2,2*maf[s]*(1-maf[s]),maf[s]^2))
geno=as.data.frame(geno)
names(geno)=paste('snp',1:10,sep='')
corrplot(corr=cor(geno),method='square',type='upper')
attach(geno)

f0=0.1; OR=1.6
pi=rep(f0,n)
pi=1/(1+exp(-log(f0) -log(OR)*snp1 - log(OR)*snp8))
pheno=waffect(prob=pi,count=250,label=c(1,0))

save(geno,pheno,file='data/ex_snps.Rda')

#Checking
load('data/ex_snps.Rda')
one_test=function(x) fisher.test(x,y=pheno)$p
(pvalues=apply(geno,2,one_test))
plot(-log10(pvalues),type='h')
abline(h=-log10(0.05/10),lty=2)


#############
# Example dependent SNPs
#############
require(corrplot)
n=500

maf=0.2; snp1=sample(0:2,size=n,replace=TRUE,prob=c((1-maf)^2,2*maf*(1-maf),maf^2))
snp2=snp1; snp2[sample(x=1:n,150,replace=FALSE)]=1
snp3=snp2; snp3[sample(x=1:n,100,replace=FALSE)]=0
maf=0.1; snp4=sample(0:2,size=n,replace=TRUE,prob=c((1-maf)^2,2*maf*(1-maf),maf^2))
snp5=snp4; snp5[sample(x=1:n,150,replace=FALSE)]=1

#maf=0.3; snp6=sample(0:2,size=n,replace=TRUE,prob=c((1-maf)^2,2*maf*(1-maf),maf^2))
snp6=snp4; snp6[sample(x=1:n,150,replace=FALSE)]=1
snp7=snp6; snp7[sample(x=1:n,150,replace=FALSE)]=1
snp8=snp7; snp8[sample(x=1:n,100,replace=FALSE)]=0
maf=0.1; snp9=sample(0:2,size=n,replace=TRUE,prob=c((1-maf)^2,2*maf*(1-maf),maf^2))
snp10=snp9; snp10[sample(x=1:n,150,replace=FALSE)]=1

geno_d=2-data.frame(SNP1=snp1,SNP2=snp2,SNP3=snp3,SNP4=snp4,SNP5=snp5,SNP6=snp6,SNP7=snp7,SNP8=snp8,SNP9=snp9,SNP10=snp10)
corrplot(corr=cor(geno_d),method='square',type='upper')

save(geno_d,pheno,file='data/ex_snps_d.Rda')


#################
# Exercise, 1000 genomes data
#################
#rm(list = ls())
require(waffect)

system('plink --file data/chrX_629 --recodeA --out data/geno')
geno=read.table('data/geno.raw',header=TRUE)
geno=geno[,7:1006]

#f0_s=1.2; OR1_s=1.5; OR2_s=1.7 #2, 1.5
f0_s=0.4; OR1_s=1.6
#pi_s=1/(1+exp(-log(f0_s)-log(OR1_s)*geno[,150]-log(OR2_s)*geno[,950]))
pi_s=1/(1+exp(-log(f0_s)-log(OR1_s)*geno[,150]))
pheno=waffect(prob=pi_s,label=c(1,0),count=314)

#save(geno,pheno, pi_s, f0_s, OR1_s, OR2_s,file='data/exercise_geno6.Rda')
save(geno,pheno, pi_s, f0_s, OR1_s,file='data/exercise_geno.Rda')

###################################################################################

int=c




#DA QUI IN POI PULIRE
##################
# Exercise, 1000 genomes data, environmental exposure
##################

#########
#Power study case-contol
#########
load('data/exercise_geno.Rda')
cores <- detectCores()
N=500; p=ncol(geno); n=nrow(geno)

e=sample(c(0,1),prob = c(0.5,0.5),size=n,replace=TRUE)
maf=0.2; geno[,150]=sample(0:1,size=n,replace=TRUE,prob=c((1-maf)^2,2*maf*(1-maf)+(maf)^2))
pi_e=1/(1+exp(-log(0.05/(1-0.05))-log(1.5)*geno[,150]-log(2)*e-log(1.5)*geno[,150]*e))

require(waffect)
one.sim_caco=function(i){
  perm=waffect(pi_e,label=c(1,0),count=315)
  out=summary(glm(perm~geno[,150]*e,family=binomial(link='logit')))$coef[4,4]
  return(out)
}
pvalue_150_caco=mclapply(1:N, one.sim_caco, mc.cores=cores)
pvalue_150_caco=as.numeric(pvalue_150_caco)
mean(pvalue_150_caco<=0.05)

#########
#Power study case only
#########
one.sim_ca=function(i){
  perm=waffect(pi_e,label=c(1,0),count=315)
  new.data=geno[which(perm==1),]
  out=summary(glm(e~geno[,150],family=binomial(link='logit')))$coef[2,4]
  return(out)
}
pvalue_150_ca=mclapply(1:N, one.sim_ca, mc.cores=cores)
pvalue_150_ca=as.numeric(pvalue_150_ca)
mean(pvalue_150_ca<=0.05)
