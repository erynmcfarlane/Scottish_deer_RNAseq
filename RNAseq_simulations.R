### want to simulate data RNAseq and genotype (or Q score?) data
## want to try:
##1) smurf (sparse multitype regularized feature)
###2) corncob loop?
####3) hierarchical model in Stan, as described in Harrison et al. 2020 MER

### simulate the data
### following https://aosmith.rbind.io/2018/01/09/simulate-simulate-part1/

library(purrr) # v. 0.3.4
library(broom) # v. 0.5.6
library(dplyr) # v. 1.0.0
library(ggplot2) # v. 3.3.1
library(dirmult)

set.seed(42)

##let's simulate 30 individuals, 1000 genes and two tissue types
### this would give a dataframe of 30 rows, and 2000 columns (gene 1- heart, gene 1 - muscle... gene 1000 - heart, gene 1000 - muscle)
n=50 ### number of individuals
num_genes<-1000 ### number of genes
genes<-as.factor(rep(c(1:num_genes), 2))
tissues<-as.factor(c(rep("heart", n), rep("muscle", n)))
Q<-rep(rbeta(n, 0.25,0.25))
### alternatively to Q, I can simulate diagnostic markers
### NB, the Q and the genotypes aren't consistent, I've simulated them separately

genotypes<-matrix(rep(sample(0:2, 1000, replace=TRUE, prob=c(rdirichlet(n=1,alpha=c(0.5,0.3,0.2)))), n), nrow=n, ncol=num_genes)

### stack the genotypes just because each individual is sampled twice, once for each tissue.
genotypes_twice<-rbind(genotypes, genotypes)
genotypes_twice<-as.data.frame(genotypes_twice)

age<-round(rnorm(n, mean=2.5), 1)
### expression is count data, but can only sum to 30 000 000 per individual 
###matrix, length = n, width = 2* genes

deername<-as.factor(rep(1:n,2))
alpha<-rep(1, num_genes)
## I've simulated using a dirichlet, but then I've multipled by 30... so these numbers can be interpreted in millions. 
heart_expression<-(rdirichlet(n=n, alpha=c(alpha)))*30000000
muscle_expression<-(rdirichlet(n=n, alpha=c(alpha)))*30000000

### this is my data frame, where each column is a gene, and each individual is represented in two rows, once for heart and one for muscle
expression<-rbind(heart_expression, muscle_expression)

expression.df<-as.data.frame(expression)
phenotypes<-cbind(deername,Q, tissues, age, genotypes_twice) 
names(phenotypes)

data<-cbind(phenotypes, expression)

library(smurf)
formu<-list()
#### what's my response variable going to be?

### right now this is just for expression of the first gene.
### need to write a loop to do expression of each of the genes.
for(i in 1:3){
 formu<-expression[,i]~p(deername, pen="lasso")+
   p(Q, pen="lasso")+
     p(tissues, pen="flasso")+
    p(age, pen="none")
  expression.fit<-glmsmurf(formula=formu, family=gaussian(), data=phenotypes, pen.weights="glm.stand", lambda = "is.aic")
}
summary(expression.fit)  
plot(expression.fit)



###let's try to do a multivariate lasso?
library(glmnet)

fit<-glmnet(phenotypes, expression.df, family="mgaussian")
plot(fit, xvar="lambda", label=TRUE, type.coef="2norm")


### I have no idea how to interpret this, but I did get it to work!
### needs to be said, m-gaussion has to be the wrong distribution for my data. 

### let's look at corncob?
library(corncob)
library(magrittr)
library(phyloseq)


### not going to include the genotype data in the phenotypes - this isn't a sparse model approach. 
phyloseq::otu_table(expression, taxa_are_rows=FALSE)->expression.phyloseq
phyloseq::sample_data(phenotypes[,1:4])->phenotypes.phyloseq
phyloseq(expression.phyloseq, phenotypes.phyloseq)->data.phyloseq


da_analysis<-differentialTest(formula = ~ Q+age+tissue+deername, 
                              phi.formula = ~ Q, 
                              formula_null= ~1,
                              test="LTR", 
                              boot=FALSE, 
                              data=data.phyloseq,
                              fdr_cutoff=0.05)

#### this isn't working even with the example data and I don't know why
soil <- soil_phylum_small
corncob<-bbdml(formula=sp1~1, 
               phi.formula=~1, 
               data=soil)
