##
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")

## Loading data
myY<- read.table("pheno_for_gapit.txt", header = TRUE)

myG<- read.delim("final_AD_gwas_snp.hmp.txt", header = FALSE)

myCV <- read.table("AD_CV_structure.txt", header= TRUE) ### don't forget to delete one Q from CV file before running

## GLM model
myGAPIT<- GAPIT(
  Y=myY[,1:7],
  G=myG,
  PCA.total =6,
  model ="GLM",
  Multiple_analysis = TRUE
)

## Blink model
myGAPIT<- GAPIT(
  Y=myY[,1:8],
  G=myG,
  CV=myCV,
  model ="Blink",
  Multiple_analysis = TRUE,
cutOff=0.05,
Geno.View.output =FALSE,
Random.model = TRUE,
N.sig=75,
PCA.View.output = FALSE,
Phenotype.View= FALSE,
PCA.3d = FALSE
)

## CMLM
myGAPIT<- GAPIT(
  Y=myY[,1:8],
  G=myG,
  group.from = 1,
  group.to = 317,
  group.by = 10,
  model ="CMLM",
  Multiple_analysis = TRUE,
  cutOff=0.05,
Geno.View.output =FALSE,
Random.model = TRUE,
N.sig=100,
PCA.View.output = FALSE,
Phenotype.View= FALSE
)

### ECMLM model
ECMLM <- GAPIT(
Y=myY,
G=myG,
CV=myCV,
kinship.cluster=c("average", "complete", "ward"),
kinship.group=c("Mean", "Max"),
group.from=1,
group.to=223,
group.by=20
)
#### GAPIT all parameters
  Y = NULL, #phenotype
  G = NULL, #hapmap genotype
  GD = NULL, #numeric genotype
  GM = NULL, #genotype map information
  KI = NULL, #kinship
  Z = NULL, #Z matrix for MLM, cMLM, encMLM
  CV = NULL, #corvariance matrix
  Aver.Dis=1000,
  buspred = FALSE, #Prediction option for after GWAS MABLUP
  bin.from = 10000, #SUPER 
  bin.to = 10000, #SUPER
  bin.by = 10000, #SUPER
  cutOff = 0.05, #threshold for significant
  CV.Extragenetic = 0, # the top number of no-inheritance columns in CV
  effectunit = 1, #Simulation phenotype
  file.output = TRUE, #output option
  FDRcut = FALSE, # filter pseudo QTN based on cutOff in blink
  group.from = 1000000,#MLM
  group.to = 1000000,#MLM
  group.by = 50,#MLM
  Geno.View.output = TRUE,#genotype analysis option
  h2 = NULL, #simulation phenotype heritability
  inclosure.from = 10, #SUPER
  inclosure.to = 10, #SUPER
  inclosure.by = 10, #SUPER
  Inter.Plot = FALSE, #Interactive plot option
  Inter.type = c("m","q"), #Interactive plot type for Manhattan and QQ plots
  kinship.cluster = "average", #cMLM
  kinship.group = 'Mean',#cMLM
  kinship.algorithm = "Zhang",#cMLM
  lmpred = FALSE, #option for linear model prediction or MABLUP prediction, that could be set as multiple parameters
  model = "MLM",# model or method in GWAS or GS
  maxOut = 100, # power for top number of markers in the GWAS
  memo = NULL, #label for marking
  Model.selection = FALSE,# optimum number of CV and PCAs
  Multi_iter = FALSE, #Multiple step for FarmCPU and BLink
  Major.allele.zero = FALSE, #convert hapmap file to numeric file, set major marker as 0
  Multiple_analysis = TRUE, #option for multiple Manhattan and QQ plots
  num_regwas = 10,# the max number of  markers in second GWAS
  NQTN = NULL, #Simulation phenotype, number of QTN
  N.sig=NULL, #Random.model, Number of significant markers..N.sig=10 for getting PVE for snp 
  NJtree.group = NULL, #NJtree set number of cluster group
  NJtree.type = c("fan","unrooted"),#NJtree type
  output.numerical = FALSE,# option for output numeric files
  output.hapmap = FALSE, # option for output hapmap files
  QTN.position = NULL, #Simulation phenotype, QTN position in the order of map file
  QTNDist = "normal",
  Random.model = TRUE, #Random.model to calculate PVE
  sangwich.top = NULL, #SUPER
  sangwich.bottom = NULL,#SUPER
  SNP.P3D = TRUE,
  SNP.effect = "Add",
  SNP.impute = "Middle",
  SNP.test = TRUE,
  SNP.MAF = 0,
  SNP.FDR = 1,  
  PCA.total = 0, # PCA number
  PCA.col = NULL, #indicater colors for individuals in PCA plot
  PCA.3d = FALSE, #3D PCA plot option
  PCA.View.output = TRUE, #option for PCA plot
  Phenotype.View= TRUE, # option for phenotype view plot