setwd("C:/Users/j_sar/Desktop/guppies/Brendan_guppy_vcf_new/guppy_vcf")
vcf<- read.table("Chr_12.noheader.vcf",header=TRUE,sep="\t",stringsAsFactors = FALSE)

#Columns corresponding to each species/sex
gambusia.DNA<-c("G_holbrooki_2F","G_holbrooki_3F","G_holbrooki_4F","G_holbrooki_6M","G_holbrooki_7M","G_holbrooki_8M")
latipinna.M.DNA<-c("P_latipinna_3M","P_latipinna_4M","P_latipinna_7M")
latipinna.F.DNA<-c("P_latipinna_1F","P_latipinna_2F","P_latipinna_6F")
picta.F.DNA<-c("P_picta_1F","P_picta_2F","P_picta_4F")
picta.M.DNA<-c("P_picta_3M","P_picta_5M","P_picta_6M")
wingei.M.DNA<-c("P_wingei_7M","P_wingei_5M","P_wingei_6M")
wingei.F.DNA<-c("P_wingei_13F","P_wingei_14F","P_wingei_15F")

# Convert from vcf format to 0,1,2 genotype format (where "1" is heterozygous)
procfunc = function(x) {
  t1 = sub(":.*","",x)
  if(t1 == "0/0") return(0)
  else if(t1 == "0/1" | t1 == "1/0") return(1)
  else if(t1 == "1/1") return(2)
  else return(-1)
}

pos = vcf$POS 
geno = data.frame(apply(vcf[,10:ncol(vcf)],1:2,procfunc),stringsAsFactors = FALSE)
names(geno) = names(vcf[,10:ncol(vcf)])
geno$POS<-vcf$POS
write.table(geno,file="guppies.geno",append=FALSE,quote=FALSE,sep="\t",row.names=FALSE) #Save genotype format file



geno<- read.table("guppies.geno",header=TRUE,sep="\t",stringsAsFactors = FALSE)

# Remove any SNPs that are missing data for all samples from one species/sex
geno.nodels<-geno[apply(geno[,c(picta.M.DNA)],1,sum)!=-3 & apply(geno[,c(picta.F.DNA)],1,sum)!=-3 
                  & apply(geno[,c(wingei.M.DNA)],1,sum)!=-3 & apply(geno[,c(wingei.F.DNA)],1,sum)!=-3 
                  & apply(geno[,c(latipinna.M.DNA,latipinna.F.DNA)],1,sum)!=-6& apply(geno[,c(gambusia.DNA)],1,sum)!=-6,]



#Retain only the important columns
geno.pw<-geno.nodels[,c("POS",picta.M.DNA,picta.F.DNA,wingei.M.DNA,wingei.F.DNA, latipinna.M.DNA,latipinna.F.DNA,gambusia.DNA)]

#Need at least two picta males & two picta females to identify sex-linked SNPs 
geno.pw<-geno.pw[apply(geno.pw[,c(picta.M.DNA)],1,function(x) length(which(x=="-1")))<2,]
geno.pw<-geno.pw[apply(geno.pw[,c(picta.F.DNA)],1,function(x) length(which(x=="-1")))<2,]
geno.pw<-geno.pw[apply(geno.pw[,c(wingei.M.DNA)],1,function(x) length(which(x=="-1")))<2,]
geno.pw<-geno.pw[apply(geno.pw[,c(wingei.F.DNA)],1,function(x) length(which(x=="-1")))<2,]

#For each SNP, calculate frequencies of possible genotypes for each species/sex
geno.pw$picta.males.hetfreq<-apply(geno.pw[,c(picta.M.DNA)],1,function(x) length(which(x=="1")))/apply(geno.pw[,c(picta.M.DNA)],1,function(x) length(which(x!="-1")))
geno.pw$picta.males.00freq<-apply(geno.pw[,c(picta.M.DNA)],1,function(x) length(which(x=="0")))/apply(geno.pw[,c(picta.M.DNA)],1,function(x) length(which(x!="-1")))
geno.pw$picta.males.11freq<-apply(geno.pw[,c(picta.M.DNA)],1,function(x) length(which(x=="2")))/apply(geno.pw[,c(picta.M.DNA)],1,function(x) length(which(x!="-1")))
geno.pw$picta.females.hetfreq<-apply(geno.pw[,c(picta.F.DNA)],1,function(x) length(which(x=="1")))/apply(geno.pw[,c(picta.F.DNA)],1,function(x) length(which(x!="-1")))
geno.pw$picta.females.00freq<-apply(geno.pw[,c(picta.F.DNA)],1,function(x) length(which(x=="0")))/apply(geno.pw[,c(picta.F.DNA)],1,function(x) length(which(x!="-1")))
geno.pw$picta.females.11freq<-apply(geno.pw[,c(picta.F.DNA)],1,function(x) length(which(x=="2")))/apply(geno.pw[,c(picta.F.DNA)],1,function(x) length(which(x!="-1")))
geno.pw$wingei.males.hetfreq<-apply(geno.pw[,c(wingei.M.DNA)],1,function(x) length(which(x=="1")))/apply(geno.pw[,c(wingei.M.DNA)],1,function(x) length(which(x!="-1")))
geno.pw$wingei.males.00freq<-apply(geno.pw[,c(wingei.M.DNA)],1,function(x) length(which(x=="0")))/apply(geno.pw[,c(wingei.M.DNA)],1,function(x) length(which(x!="-1")))
geno.pw$wingei.males.11freq<-apply(geno.pw[,c(wingei.M.DNA)],1,function(x) length(which(x=="2")))/apply(geno.pw[,c(wingei.M.DNA)],1,function(x) length(which(x!="-1")))
geno.pw$wingei.females.hetfreq<-apply(geno.pw[,c(wingei.F.DNA)],1,function(x) length(which(x=="1")))/apply(geno.pw[,c(wingei.F.DNA)],1,function(x) length(which(x!="-1")))
geno.pw$wingei.females.00freq<-apply(geno.pw[,c(wingei.F.DNA)],1,function(x) length(which(x=="0")))/apply(geno.pw[,c(wingei.F.DNA)],1,function(x) length(which(x!="-1")))
geno.pw$wingei.females.11freq<-apply(geno.pw[,c(wingei.F.DNA)],1,function(x) length(which(x=="2")))/apply(geno.pw[,c(wingei.F.DNA)],1,function(x) length(which(x!="-1")))
geno.pw$outgroups.hetfreq<-apply(geno.pw[,c(latipinna.M.DNA,latipinna.F.DNA,gambusia.DNA)],1,function(x) length(which(x=="1")))/apply(geno.pw[,c(latipinna.M.DNA,latipinna.F.DNA,gambusia.DNA)],1,function(x) length(which(x!="-1")))
geno.pw$outgroups.00freq<-apply(geno.pw[,c(latipinna.M.DNA,latipinna.F.DNA,gambusia.DNA)],1,function(x) length(which(x=="0")))/apply(geno.pw[,c(latipinna.M.DNA,latipinna.F.DNA,gambusia.DNA)],1,function(x) length(which(x!="-1")))
geno.pw$outgroups.11freq<-apply(geno.pw[,c(latipinna.M.DNA,latipinna.F.DNA,gambusia.DNA)],1,function(x) length(which(x=="2")))/apply(geno.pw[,c(latipinna.M.DNA,latipinna.F.DNA,gambusia.DNA)],1,function(x) length(which(x!="-1")))

#All picta females must be homozygous for same allele
geno.filter<-geno.pw[geno.pw$picta.females.00freq==1 | geno.pw$picta.females.11freq==1, ]

#All wingei females must be homozygous for same allele
geno.filter<-geno.filter[geno.filter$wingei.females.00freq==1 | geno.filter$wingei.females.11freq==1, ]

#All picta males must be heterozygous or homozygous for same allele as females
geno.filter<-geno.filter[geno.filter$picta.males.hetfreq==1 | (geno.filter$picta.males.hetfreq==0 & geno.filter$picta.males.00freq==geno.filter$picta.females.00freq), ]

#All wingei males must be heterozygous or homozygous for same allele as females
geno.filter<-geno.filter[geno.filter$wingei.males.hetfreq==1 | (geno.filter$wingei.males.hetfreq==0 & geno.filter$wingei.males.00freq==geno.filter$wingei.females.00freq), ]

#Outgroups must be fixed for same allele
geno.filter<-geno.filter[geno.filter$outgroups.00freq==1 | geno.filter$outgroups.11freq==1, ]

#identify derived allele based on outgroup
geno.filter$derivedAllele<-ifelse(geno.filter$outgroups.00freq==1,1,0)

#Identify which chromosomes have derived alleles
geno.filter$picta.X.derived<-ifelse(geno.filter$picta.females.00freq==geno.filter$outgroups.00freq,0,1)
geno.filter$wingei.X.derived<-ifelse(geno.filter$wingei.females.00freq==geno.filter$outgroups.00freq,0,1)
geno.filter$picta.Y.derived<-ifelse(geno.filter$picta.males.hetfreq==1,
                                    ifelse(geno.filter$picta.X.derived==0,1,0),
                                    ifelse(geno.filter$picta.males.00freq==geno.filter$outgroups.00freq,0,1))
geno.filter$wingei.Y.derived<-ifelse(geno.filter$wingei.males.hetfreq==1,
                                    ifelse(geno.filter$wingei.X.derived==0,1,0),
                                    ifelse(geno.filter$wingei.males.00freq==geno.filter$outgroups.00freq,0,1))

# Either 2 or 3 of Xw, Yw, Xp, and Yp must have derived allele to be topologically informative
geno.filter$NumDerived<-geno.filter$picta.X.derived+geno.filter$wingei.X.derived+geno.filter$picta.Y.derived+geno.filter$wingei.Y.derived
geno.filter<-geno.filter[geno.filter$NumDerived>1 & geno.filter$NumDerived<4,]
             
nrow(geno.filter)
                       
 
# Identify topology associated with each SNP
geno.filter$derivedAlleCat<-0
# Cat 1 = (0,1,0,1) -> homologous SDR
geno.filter$derivedAlleCat<-ifelse(geno.filter$wingei.X.derived==0 & geno.filter$wingei.Y.derived==1 & geno.filter$picta.X.derived == 0 & geno.filter$picta.Y.derived ==1, 1, geno.filter$derivedAlleCat)
# Cat 2 = (1,0,1,0) -> homologous SDR
geno.filter$derivedAlleCat<-ifelse(geno.filter$wingei.X.derived==1 & geno.filter$wingei.Y.derived==0 & geno.filter$picta.X.derived == 1 & geno.filter$picta.Y.derived ==0, 2, geno.filter$derivedAlleCat)
# Cat 3 = (0,0,1,1) -> independent evolution / hemizygosity
geno.filter$derivedAlleCat<-ifelse(geno.filter$wingei.X.derived==0 & geno.filter$wingei.Y.derived==0 & geno.filter$picta.X.derived == 1 & geno.filter$picta.Y.derived ==1, 3, geno.filter$derivedAlleCat)
# Cat 4 = (1,1,0,0) -> independent evolution / hemizygosity
geno.filter$derivedAlleCat<-ifelse(geno.filter$wingei.X.derived==1 & geno.filter$wingei.Y.derived==1 & geno.filter$picta.X.derived == 0 & geno.filter$picta.Y.derived ==0, 4, geno.filter$derivedAlleCat)
# Cat 5 = (1,1,1,0) -> turnover Xw->Yw, hemizygous wingei
geno.filter$derivedAlleCat<-ifelse(geno.filter$wingei.X.derived==1 & geno.filter$wingei.Y.derived==1 & geno.filter$picta.X.derived == 1 & geno.filter$picta.Y.derived ==0, 5, geno.filter$derivedAlleCat)
# Cat 6 = (1,0,1,1) -> turnover Xp->Yp, hemizygous picta
geno.filter$derivedAlleCat<-ifelse(geno.filter$wingei.X.derived==1 & geno.filter$wingei.Y.derived==0 & geno.filter$picta.X.derived == 1 & geno.filter$picta.Y.derived ==1, 6, geno.filter$derivedAlleCat)
# Cat 7 = (1,1,0,1) -> turnover Yw->Xw, hemizygous wingei
geno.filter$derivedAlleCat<-ifelse(geno.filter$wingei.X.derived==1 & geno.filter$wingei.Y.derived==1 & geno.filter$picta.X.derived == 0 & geno.filter$picta.Y.derived ==1, 7, geno.filter$derivedAlleCat)
# Cat 8 = (0,1,1,1) -> turnover Yp->Xp, hemizygous picta
geno.filter$derivedAlleCat<-ifelse(geno.filter$wingei.X.derived==0 & geno.filter$wingei.Y.derived==1 & geno.filter$picta.X.derived == 1 & geno.filter$picta.Y.derived ==1, 8, geno.filter$derivedAlleCat)
# Cat 9 = (1,0,0,1) -> nonsense
geno.filter$derivedAlleCat<-ifelse(geno.filter$wingei.X.derived==1 & geno.filter$wingei.Y.derived==0 & geno.filter$picta.X.derived == 0 & geno.filter$picta.Y.derived ==1, 9, geno.filter$derivedAlleCat)
# Cat 10 = (0,1,1,0) -> nonsense
geno.filter$derivedAlleCat<-ifelse(geno.filter$wingei.X.derived==0 & geno.filter$wingei.Y.derived==1 & geno.filter$picta.X.derived == 1 & geno.filter$picta.Y.derived ==0, 10, geno.filter$derivedAlleCat)

#Counts of all possible topologies
aggregate(geno.filter$derivedAlleCat,by=list(geno.filter$derivedAlleCat),length)




