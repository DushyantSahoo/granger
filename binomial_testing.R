#! /usr/bin/Rscript

args = commandArgs()
print(args[1])
# Help
if (identical(args[1],"--help")){
  print("example: binomial_testing list_of_nii_files.txt mask.nii.gz top_voxels 2 2 outfile ")
  print("example: binomial_testing arg1 arg2 arg3 arg4 arg5 arg6")
  print("arg1: list of nii files stored in text")
  print("arg2: mask file (avg152T1_brain.nii.gz is given along the source code)")
  print("arg3: top voxels")
  print("arg4: 1 if one sided or 2 if two sided")
  print("arg5: 1 if symmetric or 2 if non-symmetric (if two sided testing)")
  print("arg6: output folder location")
  print("Output-")
  print("1) log of q_value (nii file)")
  print("2) log of p_value (nii file)")
}

stop("Stop")
# script for binomial testing
library(AnalyzeFMRI)
library(fastICA)
library(oro.nifti)
library(R.matlab)
library(base)
library(fmsb)

# Read data table containing file names
mydata = read.table(args[1])


BB<-readNIfTI(args[2], reorient = FALSE)
B<-BB@.Data
temp2 <- 1
for ( i in 1:91){
  for (j in 1:109){
    for (k in 1:91){
      if (B[i,j,k]>0 ){
        temp2 <- temp2+1
      }
    }
  }
}

# parameter intialization
example <- readNIfTI(mydata[1],reorient = FALSE)
example_data <- example@.Data
time <- dim(example_data)[4]
bigdata <- array(0,dim=c(dim(example_data)[1],dim(example_data)[2],dim(example_data)[3],time))
bigdata1 <- array(0,dim=c(dim(example_data)[1],dim(example_data)[2],dim(example_data)[3],time))

# Function for binomial testing
binomial <- function(Time,top,input1,input2,Data,siz){
  for (timer in 1:Time){
    
    # one sided
    if (input1 == 1){
      simp1 <- array(0 , dim=c(length(Data),temp2))
      for (sublength in 1:length(Data)){
        
        # Read data file
        AA<-readNIfTI(Data[sublength], reorient = FALSE)
        A1<-AA@.Data
        
        A <- A1[,,,timer]
        
        xlength <- AA@pixdim[2]
        ylength <- AA@pixdim[3]
        zlength <- AA@pixdim[4]
        value <- dim(A)
        a<-value[1]
        b<-value[2]
        c<-value[3]
        
        # transform 3d data to 1d data
        temp<- 1
        simp <- array(0 , dim=c(siz))
        for ( i in 1:a){
          for (j in 1:b){
            for (k in 1:c){
              if (B[i,j,k]>0 ){
                simp[temp] <- A[i,j,k]
                temp <- temp+1
              }
            }
          }
        }
        xc <- 1
        thres <- quantile(simp,1-(top/100))
        for ( ii in 1:(siz)){
          if (simp[ii]>thres )
          {
            simp1[sublength,ii] <- 1
            xc = xc + 1
          }
          else {
            simp1[sublength,ii] <- 0
          }
        }
      }
      
      sum1 <- colSums(simp1)
      for ( ii in 1:(siz)){
        sum1[ii] <- binom.test(sum1[ii], length(sub), p = top/100 , alternative = "greater" , conf.level = 1-(top/100))$p.value
      }
      final1 <- p.adjust(sum1, method = "fdr", n = length(sum1))
      temp1 <- 1
      new2 <- array(0 , dim=c(a,b,c))
      new4 <- array(0 , dim=c(a,b,c))
      for ( i in 1:a){
        for (j in 1:b){
          for (k in 1:c){
            if (B[i,j,k]>0 ){
              new2[i,j,k] <- -log10(final1[temp1] )
              new4[i,j,k] <- -log10(sum1[temp1] )
              temp1 <- temp1+1
            }
          }
        }
      }
    } 
    # 2 sided symmteric
    else if (input1 == 1 && input2 == 1){
      simp1 <- array(0 , dim=c(length(Data),temp2))
      simp2 <- array(0 , dim=c(length(Data),temp2))
      for (sublength in 1:length(Data)){
        
        # Read data file
        AA<-readNIfTI(Data[sublength], reorient = FALSE)
        A1<-AA@.Data
        
        A <- A1[,,,timer]
        
        xlength <- AA@pixdim[2]
        ylength <- AA@pixdim[3]
        zlength <- AA@pixdim[4]
        value <- dim(A)
        a<-value[1]
        b<-value[2]
        c<-value[3]
        
        # transform 3d data to 1d data
        temp<- 1
        simp <- array(0 , dim=c(siz))
        for ( i in 1:a){
          for (j in 1:b){
            for (k in 1:c){
              if (B[i,j,k]>0 ){
                simp[temp] <- A[i,j,k]
                temp <- temp+1
              }
            }
          }
        }
        xc <- 1
        thres <- quantile(simp,1-(top/200))
        thres1 <- quantile(simp,(top/200))
        for ( ii in 1:(siz)){
          if (simp[ii]>thres )
          {
            simp1[sublength,ii] <- 1
            xc = xc + 1
          }
          else {
            simp1[sublength,ii] <- 0
          }
        }
        xc <- 1
        for ( ii in 1:(siz)){
          if (simp[ii]<(thres1))
          {
            simp2[sublength,ii] <- 1
            xc = xc + 1
          }
          else {
            simp2[sublength,ii] <- 0
          }
        }
      }
      
      sum1 <- colSums(simp1)
      sum2 <- colSums(simp2)
      for ( ii in 1:(siz)){
        sum2[ii] <- binom.test(sum2[ii], length(sub), p = top/100 , alternative = "greater" , conf.level = 1-(top/100))$p.value
        
        sum1[ii] <- binom.test(sum1[ii], length(sub), p = top/100 , alternative = "greater" , conf.level = 1-(top/100))$p.value
      }
      final1 <- p.adjust(sum1, method = "fdr", n = length(sum1))
      final2 <- p.adjust(sum2, method = "fdr", n = length(sum2))
      temp1 <- 1
      new1 <- array(0 , dim=c(a,b,c))
      new2 <- array(0 , dim=c(a,b,c))
      new3 <- array(0 , dim=c(a,b,c))
      new4 <- array(0 , dim=c(a,b,c))
      for ( i in 1:a){
        for (j in 1:b){
          for (k in 1:c){
            if (B[i,j,k]>0 ){
              new1[i,j,k] <- log10(final2[temp1])
              new2[i,j,k] <- -log10(final1[temp1] )
              new3[i,j,k] <- log10(sum2[temp1])
              new4[i,j,k] <- -log10(sum1[temp1] )
              temp1 <- temp1+1
            }
          }
        }
      }
      new2 = new2 + new1
      new4 = new4 + new3
    }
    # 2 sided non-symmteric
    else {
      simp1 <- array(0 , dim=c(length(Data),temp2))
      simp2 <- array(0 , dim=c(length(Data),temp2))
      for (sublength in 1:length(Data)){
        
        # Read data file
        AA<-readNIfTI(Data[sublength], reorient = FALSE)
        A1<-AA@.Data
        
        A <- A1[,,,timer]
        
        xlength <- AA@pixdim[2]
        ylength <- AA@pixdim[3]
        zlength <- AA@pixdim[4]
        value <- dim(A)
        a<-value[1]
        b<-value[2]
        c<-value[3]
        
        # transform 3d data to 1d data
        temp<- 1
        simp <- array(0 , dim=c(siz))
        for ( i in 1:a){
          for (j in 1:b){
            for (k in 1:c){
              if (B[i,j,k]>0 ){
                simp[temp] <- A[i,j,k]
                temp <- temp+1
              }
            }
          }
        }
        xc <- 1
        new1 <- abs(simp)
        thres <- quantile(new1,1-(top/100))
        for ( ii in 1:(siz)){
          if (simp[ii]>thres )
          {
            simp1[sublength,ii] <- 1
            xc = xc + 1
          }
          else {
            simp1[sublength,ii] <- 0
          }
        }
        xc <- 1
        for ( ii in 1:(siz)){
          if (simp[ii]<(-thres))
          {
            simp2[sublength,ii] <- 1
            xc = xc + 1
          }
          else {
            simp2[sublength,ii] <- 0
          }
        }
      }
      
      sum1 <- colSums(simp1)
      sum2 <- colSums(simp2)
      for ( ii in 1:(siz)){
        sum2[ii] <- binom.test(sum2[ii], length(sub), p = top/100 , alternative = "greater" , conf.level = 1-(top/100))$p.value
        
        sum1[ii] <- binom.test(sum1[ii], length(sub), p = top/100 , alternative = "greater" , conf.level = 1-(top/100))$p.value
      }
      final1 <- p.adjust(sum1, method = "fdr", n = length(sum1))
      final2 <- p.adjust(sum2, method = "fdr", n = length(sum2))
      temp1 <- 1
      new1 <- array(0 , dim=c(a,b,c))
      new2 <- array(0 , dim=c(a,b,c))
      new3 <- array(0 , dim=c(a,b,c))
      new4 <- array(0 , dim=c(a,b,c))
      for ( i in 1:a){
        for (j in 1:b){
          for (k in 1:c){
            if (B[i,j,k]>0 ){
              new1[i,j,k] <- log10(final2[temp1])
              new2[i,j,k] <- -log10(final1[temp1] )
              new3[i,j,k] <- log10(sum2[temp1])
              new4[i,j,k] <- -log10(sum1[temp1] )
              temp1 <- temp1+1
            }
          }
        }
      }
      new2 = new2 + new1
      new4 = new4 + new3
    }
    
    
    bigdata[,,,timer]<-new2
    bigdata1[,,,timer]<-new4
  }
  return(list(bigdata,bigdata1))
}
if (args[4] == 1){
  f.write.analyze(mat = binomial2(time,args[3],args[4],args[5],mydata,temp2)[1],file = paste(args[6],logq_value,sep=""),size = "float",pixdim=c(xlength,ylength,zlength))
  f.write.analyze(mat = binomial2(time,args[3],args[4],args[5],mydata,temp2)[2],file = paste(args[6],logp_value,sep=""),size = "float",pixdim=c(xlength,ylength,zlength))
} else if (args[4] == 2 && args[5] == 1){
  f.write.analyze(mat = binomial1(time,args[3],args[4],args[5],mydata,temp2)[1],file = paste(args[6],logq_value,sep=""),size = "float",pixdim=c(xlength,ylength,zlength))
  f.write.analyze(mat = binomial1(time,args[3],args[4],args[5],mydata,temp2)[2],file = paste(args[6],logp_value,sep=""),size = "float",pixdim=c(xlength,ylength,zlength))
} else if (args[4] == 2 && args[5] == 2){
  f.write.analyze(mat = binomial(time,args[3],args[4],args[5],mydata,temp2)[1],file = paste(args[6],logq_value,sep=""),size = "float",pixdim=c(xlength,ylength,zlength))
  f.write.analyze(mat = binomial(time,args[3],args[4],args[5],mydata,temp2)[2],file = paste(args[6],logp_value,sep=""),size = "float",pixdim=c(xlength,ylength,zlength))
} else {
  print("Incorrect input")
}
