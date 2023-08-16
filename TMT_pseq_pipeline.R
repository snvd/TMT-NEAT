#Pipeline for SLN, zero imputation, and IRS on TMT data followed by DE analysis using PoissonSeq
#Last updated: Jan 06, 2023 by CMS

TMT_pseq_pipeline <- function(workdir, datafile, metadatafile, exp, REGEX, SLN, PTM, DE, stat, qval,compsfile){

#make sure the exp name is syntactically valid if not using REGEXP
if(REGEX == "no") {
  exp <- make.names(exp)
}

print(paste0("Working directory: ", workdir))
setwd(workdir)
message("Loading data files...")

#load the data file
#make the peptide ids the row names
data <- read.delim(datafile,row.names = "id", stringsAsFactors = FALSE)

colnames(data)[1] <- "Proteins"

#read in the metadata
metadataorg <- read.table(metadatafile,sep="\t",header=TRUE,stringsAsFactors=FALSE)

#remove blanks
noblanks <- metadataorg[!grepl('blank',metadataorg$name,ignore.case=TRUE),]
metadata <- noblanks
plex <- length(metadata$sample[metadata$run==1]) #number of TMT channels used
runs <- max(metadata$run) #number of LC-MS/MS runs for this experiment
reps <- metadata$rep #number of replicates per condition
numrefs <- sum(grepl("Ref",metadata$name,ignore.case=TRUE))/runs

message("Beginning data processing...")
if (PTM == "P"){
  #get rid of peptides with no phospho sites
  data = data[!(data$Number.of.Phospho..STY.==""),]
  #get the raw intensities for differential abundance
  intensities = data[,grepl('Reporter.intensity',colnames(data))]
  rawintensities = intensities[!grepl('corrected',colnames(intensities))]
  rawintensities = rawintensities[!grepl('count',colnames(rawintensities))]
  
  if(exp=="None"){
    intensitiesK=rawintensities
  }else{
    intensitiesK = rawintensities[,grepl(exp,colnames(rawintensities))] 
  }
  #remove blanks
  metadata3 = rep(metadataorg$name,each=3)
  intensitiesK = intensitiesK[,!grepl('blank',metadata3,ignore.case=TRUE)]

  
  message("Separating multiplicities...")
  #first we need to make a new table where we 1) take only the columns we need and 2) separate peptides by multiplicity
  #this assumes columns are always named the same thing
  currentrow=1
  for (i in 1:length(row.names(data))){
    myprotein = data[i,]
    #check number of phospho sites
    #if there is a semicolon, that means we have multiple sites and need to split into multiple rows
    if (length(grep(";",myprotein$Number.of.Phospho..STY.))>0){
      #get the info for this protein
      mydata = data.frame(myprotein$Proteins,myprotein$Positions.within.proteins,myprotein$Protein,myprotein$Fasta.headers,
                          myprotein$Localization.prob,myprotein$Number.of.Phospho..STY.,
                          myprotein$Sequence.window, myprotein$Modification.window, myprotein$Peptide.window.coverage,
                          myprotein$Phospho..STY..Probabilities, row.names(myprotein), stringsAsFactors=FALSE)
      colnames(mydata)=c("Proteins","Positions.within.proteins","Protein","Fasta.headers",
                         "Localization.prob",
                         "Phospho.Site","Sequence.window",'Modification.window',
                         'Peptide.window.coverage','Phospho.site.probs','original.id')
      #split the phospho sites by semincolon
      sites = strsplit(myprotein$Number.of.Phospho..STY.,";")
      mysites = sites[[1]]
      for (j in 1:length(mysites)){
        #if the site # is >3, we need to make it 3
        if (mysites[j]>3){
          currentsite = 3
        } else{
          currentsite=mysites[j]
        }
        myintensities = intensitiesK[i,grepl(paste('_',currentsite,sep=""),colnames(intensitiesK))]
        colnames(myintensities)=paste(rep("R",length(reps)),reps,sep="_")
        if (i==1){
          newdata = data.frame(mydata,myintensities,stringsAsFactors=FALSE)
          row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",mysites[j],sep="")
          currentrow = currentrow+1
        } else{
          newdata[currentrow,] =  data.frame(mydata,myintensities)
          row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",mysites[j],sep="")
          currentrow = currentrow+1
        }
      }
    }else{
      #get the info for this protein
      mydata = data.frame(myprotein$Proteins,myprotein$Positions.within.proteins,myprotein$Protein,myprotein$Fasta.headers,
                          myprotein$Localization.prob,myprotein$Number.of.Phospho..STY.,
                          myprotein$Sequence.window, myprotein$Modification.window, myprotein$Peptide.window.coverage,
                          myprotein$Phospho..STY..Probabilities, row.names(myprotein), stringsAsFactors=FALSE)
      colnames(mydata)=c("Proteins","Positions.within.proteins","Protein","Fasta.headers",
                         "Localization.prob",
                         "Phospho.Site","Sequence.window",'Modification.window',
                         'Peptide.window.coverage','Phospho.site.probs','original.id')
      #get the correct intensity values for this protein
      if (myprotein$Number.of.Phospho..STY.>3){
        currentsite = 3
      } else{
        currentsite=myprotein$Number.of.Phospho..STY.
      }
      myintensities = intensitiesK[i,grepl(paste('_',currentsite,sep=""),colnames(intensitiesK))]
      colnames(myintensities)=paste(rep("R",length(reps)),reps,sep="_")
      if (i==1){
        newdata = data.frame(mydata,myintensities,stringsAsFactors=FALSE)
        row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",myprotein$Number.of.Phospho..STY.,sep="")
        currentrow = currentrow+1
      } else{
        newdata[currentrow,] =  data.frame(mydata,myintensities)
        row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",myprotein$Number.of.Phospho..STY.,sep="")
        currentrow = currentrow+1
      }
    }
  }
}else if(PTM=="U"){
  #get rid of peptides with no ubiquitin sites
  data <- data[!(data$Number.of.GlyGly..K.==""),]
  #get the raw intensities for differential abundance
  intensities <- data[,grepl('Reporter.intensity',colnames(data))]
  rawintensities <- intensities[!grepl('corrected',colnames(intensities))]
  rawintensities <- rawintensities[!grepl('count',colnames(rawintensities))]
  
  if(exp=="None"){
    intensitiesK <- rawintensities
  }else{
    intensitiesK <- rawintensities[,grepl(exp,colnames(rawintensities))] 
  }
  #remove blanks
  metadata3 <- rep(metadataorg$name,each=3)
  intensitiesK <- intensitiesK[,!grepl('blank',metadata3,ignore.case=TRUE)]
  
  message("Separating multiplicities...")
  #first we need to make a new table where we 1) take only the columns we need and 2) separate peptides by multiplicity
  #this assumes columns are always named the same thing
  currentrow <- 1
  for (i in 1:length(row.names(data))){
    myprotein <- data[i,]
    #check number of Ubiquitin sites
    #if there is a semicolon, that means we have multiple sites and need to split into multiple rows
    if (length(grep(";",myprotein$Number.of.GlyGly..K.))>0){
      #get the info for this protein
      mydata = data.frame(myprotein$Proteins,myprotein$Positions.within.proteins,myprotein$Protein,myprotein$Fasta.headers,
                          myprotein$Localization.prob,myprotein$Number.of.GlyGly..K.,
                          myprotein$Sequence.window, myprotein$Modification.window, myprotein$Peptide.window.coverage,
                          myprotein$GlyGly..K..Probabilities, row.names(myprotein), stringsAsFactors=FALSE)
      colnames(mydata)=c("Proteins","Positions.within.proteins","Protein","Fasta.headers",
                         "Localization.prob",
                         "GlyGly.Site","Sequence.window",'Modification.window',
                         'Peptide.window.coverage','GlyGly.site.probs','original.id')
      #split the Ubiquitin sites by semincolon
      sites = strsplit(myprotein$Number.of.GlyGly..K.,";")
      mysites = sites[[1]]
      for (j in 1:length(mysites)){
        #if the site # is >3, we need to make it 3
        if (mysites[j]>3){
          currentsite = 3
        } else{
          currentsite=mysites[j]
        }
        myintensities = intensitiesK[i,grepl(paste('_',currentsite,sep=""),colnames(intensitiesK))]
        colnames(myintensities)=paste(rep("R",length(reps)),reps,sep="_")
        if (i==1){
          newdata = data.frame(mydata,myintensities,stringsAsFactors=FALSE)
          row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",mysites[j],sep="")
          currentrow = currentrow+1
        } else{
          newdata[currentrow,] =  data.frame(mydata,myintensities)
          row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",mysites[j],sep="")
          currentrow = currentrow+1
        }
      }
    }else{
      #get the info for this protein
      mydata <- data.frame(myprotein$Proteins,myprotein$Positions.within.proteins,myprotein$Protein,myprotein$Fasta.headers,
                          myprotein$Localization.prob,myprotein$Number.of.GlyGly..K.,
                          myprotein$Sequence.window, myprotein$Modification.window, myprotein$Peptide.window.coverage,
                          myprotein$GlyGly..K..Probabilities, row.names(myprotein), stringsAsFactors=FALSE)
      colnames(mydata) <- c("Proteins","Positions.within.proteins","Protein","Fasta.headers",
                         "Localization.prob",
                         "GlyGly.Site","Sequence.window",'Modification.window',
                         'Peptide.window.coverage','GlyGly.site.probs','original.id')
      #get the correct intensity values for this protein
      if (myprotein$Number.of.GlyGly..K.>3){
        currentsite <- 3
      } else{
        currentsite <- as.numeric(myprotein$Number.of.GlyGly..K.)
      }
      myintensities <- intensitiesK[i,grepl(paste('_',currentsite,sep=""),colnames(intensitiesK))]
      colnames(myintensities) <- paste(rep("R",length(reps)),reps,sep="_")
      if (i==1){
        newdata = data.frame(mydata,myintensities,stringsAsFactors=FALSE)
        row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",myprotein$Number.of.GlyGly..K.,sep="")
        currentrow = currentrow+1
      } else{
        newdata[currentrow,] =  data.frame(mydata,myintensities)
        row.names(newdata)[currentrow]=paste(row.names(newdata)[currentrow],".",myprotein$Number.of.GlyGly..K.,sep="")
        currentrow = currentrow+1
      }
    }
  }
}else{
  #get the raw intensities for differential abundance
  intensities = data[,grepl('Reporter.intensity',colnames(data))]
  rawintensities = intensities[!grepl('corrected',colnames(intensities))]
  rawintensities = rawintensities[!grepl('count',colnames(rawintensities))]
  intensitiesK = rawintensities[,grepl(exp,colnames(rawintensities))]
  if(exp=="None"){
    intensitiesK=rawintensities
  }else{
    intensitiesK = rawintensities[,grepl(exp,colnames(rawintensities))] 
  }
  #remove blanks
  metadata3 = metadataorg$name
  intensitiesK = intensitiesK[,!grepl('blank',metadata3,ignore.case=TRUE)]

  #get info for all proteins
  mydata = data.frame(data$Proteins,data$Majority.protein.IDs,data$Peptide.counts..all,
                      data$Peptide.counts..unique.,data$Fasta.headers,data$Number.of.proteins,data$Peptides,
                      data$Unique.peptides, data$Sequence.coverage....,data$Unique.sequence.coverage....,
                      data$Sequence.lengths, data$Peptide.IDs, data$Evidence.IDs,
                      row.names(data), stringsAsFactors=FALSE)
  colnames(mydata)=c("Proteins","Majority.protein.IDs","Peptide.counts.all","Peptide.counts.unique",
                     "Fasta.headers",
                     "Number.of.proteins","Peptides","Unique.peptides",'Sequence.coverage',
                     'Unique.sequence.coverage','Sequence.lengths','Peptide.IDs','Evidence.IDs','original.id')
  newdata = data.frame(mydata,intensitiesK)
}

#replace column names with sample names from the metadata
colnames(newdata)[(dim(mydata)[2]+1):dim(newdata)[2]] <- paste0(metadata$name,"_",metadata$rep)

#get rid of all contaminants
if (PTM=="P" | PTM=="U"){
  newdata <-  newdata[!grepl("CON",newdata$Protein,ignore.case=TRUE),]
  newdata <-  newdata[!grepl("REV",newdata$Protein, ignore.case=TRUE),]
}else{
  newdata <-  newdata[!grepl("CON",newdata$Majority.protein.IDs,ignore.case=TRUE),]
  newdata <-  newdata[!grepl("REV",newdata$Majority.protein.IDs, ignore.case=TRUE),]
}

#finally add log2 transformed values because apparently people like those
intensities <- newdata[,(dim(mydata)[2]+1):dim(newdata)[2]]
newdatalog <- data.frame(newdata,log2(intensities+1),stringsAsFactors=FALSE)
colnames(newdatalog) <- c(colnames(newdata),paste("log2(",metadata$name,"_",metadata$rep,")",sep=""))

#save table
newdatalog |>
  rownames_to_column(var = "UID") |>
  write.csv(file = 'prenormalized_data.csv', row.names = FALSE)

message("Removing proteins in less than 50% of the runs...")

#we need to get rid of proteins that are not expressed in at least half of the runs
#to do this we can just pull the reference lanes from all the runs and weed out those with 50% missing (0) intensity value
#if there are no references, we need to impute them using the protein mean intensity

if (numrefs==0){
  refs = rowMeans(intensities)
  } else {
  refs = newdata[,grepl('Ref',colnames(newdata),ignore.case=TRUE)]
  }

if (numrefs>1) {
  for (i in 1:numrefs) {
    if (i==1) {
      refsums = data.frame(rowSums(refs[,(numrefs*(i-1)+1):(numrefs*i)]))
      } else {
      refsums = data.frame(refsums,rowSums(refs[,(numrefs*(i-1)+1):(numrefs*i)]))
      }
    }
  #count how many times the protein is zero in the reference lane
  zerosums = rowSums(refsums[,]==0)
  #remove the proteins that are zero in >=50-1% of the runs
  nozeros = newdata[zerosums <= (floor(dim(refsums)[2]/2-1)),]
  } else if (numrefs==1) {
  #remove the proteins that are zero in all samples
  nozeros <- newdata[rowSums(refs==0)==0,]
  } else {
  nozeros <- newdata[refs!=0,]
  }
#save
nozeros %>%
  rownames_to_column(var = "UID") %>%
  write.csv(file = 'prenormalized_data_in_at_least_half_of_runs.csv', row.names = FALSE)
finalintensities = nozeros[,(dim(mydata)[2]+1):dim(newdata)[2]]

#plot boxplot before normalization
#QC plots
colors = c("blue","green","orange","yellow","red","purple","white","blue","green","orange","yellow","red","purple","white")
png(filename='boxplot_log2.png',width=5000,height=2000,res=300)
invisible(b <- boxplot(log2(finalintensities+1),col=colors[unlist(lapply(1:runs,function(x) rep(x,plex)))],ylab="log2(Intensity)",cex.axis=0.75,las=2))
print(b)
dev.off()

if(SLN=="Yes"){
  message("Sample loading normalization...")
  #peform sample loading normalization (SLN)
  if (runs>1){
    for (i in 1:runs){
      #get intensity values for this run
      myints = finalintensities[,(1+(i-1)*plex):(i*plex)]
      #column normalize
      sums = colSums(myints)
      meansums = mean(sums)
      normfactor = sums/meansums
      myintensitiesnorm = myints
      for (j in 1:dim(myints)[2]){
        myintensitiesnorm[,j] = ceiling(myints[,j]/normfactor[j])
      }
      #make a table of normalized intensities
      if (i==1){
        normintensities = myintensitiesnorm
      }else{
        normintensities = data.frame(normintensities,myintensitiesnorm)
      }
    }
    write.csv(normintensities,"sample_loading_normalization.csv")
  }else{
    myints <- finalintensities
    #column normalize
    sums <- colSums(myints)
    meansums <- mean(sums)
    normfactor <- sums/meansums
    myintensitiesnorm <- myints
    for (j in 1:dim(myints)[2]){
      myintensitiesnorm[,j] <-  ceiling(myints[,j]/normfactor[j])
    }
    normintensities <- myintensitiesnorm
    myintensitiesnorm %>%
      rownames_to_column(var = "UID") %>%
      write.csv(file = "sample_loading_normalization.csv", row.names = FALSE)
  }
}else{
  normintensities <- finalintensities
}

normintensitiesimpall <- normintensities

#perform IRS
if (runs>1){
  message("Internal reference normalization...")
  refs = as.matrix(normintensitiesimpall[,grepl('Ref',colnames(normintensitiesimpall),ignore.case=TRUE)])
  if (numrefs>1){
    for (i in 1:runs){
      if (i==1){
        refscomb = data.frame(rowMeans(refs[,(numrefs*(i-1)+1):(numrefs*i)]))
      }else{
        refscomb = data.frame(refscomb,rowMeans(refs[,(numrefs*(i-1)+1):(numrefs*i)]))
      }
    }
  }else{
    refscomb=refs
  }
  refslog = log(refscomb+1)
  zerosums = rowSums(refscomb[,]==0)
  irsaverage <- rowSums(refslog)/(dim(refscomb)[2]-zerosums)
  irsaverage = exp(irsaverage)
  normfactorref = irsaverage/refscomb
  for (i in 1:dim(normfactorref)[2]){
    normfactorref[normfactorref[,i]=="Inf",i]=NA
  }
  for (i in 1:dim(refscomb)[2]){
    #get intensity values for this run
    myintensities = normintensitiesimpall[,colnames(normintensitiesimpall)%in%paste(metadata$name[metadata$run==i],"_",metadata$rep[metadata$run==i],sep="")]
    mynormfactorref = normfactorref[,i]
    #row normalize each protein against the reference
    for (j in 1:dim(myintensities)[2]){
      myintensities[,j] = ceiling(myintensities[,j]*mynormfactorref)
    }
    #make a table of normalized intensities
    if (i==1){
      normintensitiesIRS = myintensities
    }else{
      normintensitiesIRS = data.frame(normintensitiesIRS,myintensities)
    }
  }
  write.csv(normintensitiesIRS,'IRS_normalized_values.csv')
  finalimpintensitiesIRS=normintensitiesIRS
  
  metadata %>%
    filter(!grepl(pattern = "Ref", x = .[]$name), .preserve = TRUE) -> metadata
}else{
  finalimpintensitiesIRS <- normintensitiesimpall
}

#QC plots
# Normalized intensity box plot
colors = c("blue","green","orange","yellow","red","purple","white","blue","green","orange","yellow","red","purple","white")
png(filename='boxplot_log2_norm.png',width=5000,height=2000,res=300)
invisible(b <- boxplot(log2(finalimpintensitiesIRS+1),col=colors[unlist(lapply(1:runs,function(x) rep(x,plex)))],ylab="log2(Intensity)",cex.axis=0.75,las=2))
print(b)
dev.off()

#pca
if (numrefs>0){
  pcaresults = prcomp(t(na.omit(finalimpintensitiesIRS[,!grepl("Ref",colnames(finalimpintensitiesIRS),ignore.case=TRUE)])))
  #plot PCA
  g <- ggbiplot(pcaresults, groups = metadata$name[!grepl("Ref",metadata$name,ignore.case=TRUE)],
                labels=colnames(finalimpintensitiesIRS)[grep("Ref",colnames(finalimpintensitiesIRS),invert=TRUE,ignore.case=TRUE)],
                var.axes=FALSE,labels.size=3,ellipse=TRUE)
  if (length(unique(metadata$name))>2){
    g <- g+scale_color_brewer(palette="Set1")
  }
  g<- g+ theme(text = element_text(size=14))
  png(filename='PCA.png',width=5000,height=2000,res=300)
  print(g)
  dev.off()
}else{
  pcaresults = prcomp(t(finalimpintensitiesIRS))
  #plot PCA
  g <- ggbiplot(pcaresults, groups = metadata$name,
                labels=colnames(finalimpintensitiesIRS),
                var.axes=FALSE,labels.size=3,ellipse=TRUE)
 # if (length(unique(metadata$name))>2){
#    g <- g+scale_color_brewer(palette="Paired")
#  }
  g<- g+ theme(text = element_text(size=14))
  png(filename='PCA.png',width=4000,height=2000,res=300)
  print(g)
  dev.off()
}

if(DE=="Yes"){
  if (length(compsfile)==0) {
    print("Comparing everything versus everything (comps.xlsx file not provided)")
    #make list of all pairwise comparisons
    mysamples = unique(metadata$name)
    comps = data.frame()
     currentrow = 1
     for (i in 1:(length(mysamples)-1)){
       for (j in (i+1):length(mysamples)){
         comps[currentrow,1] <- mysamples[i]
         comps[currentrow,2] <- mysamples[j]
#         comps[j,1] = paste0(mysamples[i],"_vs_",mysamples[j])
         currentrow = currentrow+1
       }
     }
  } else {
    print("Comparing sample pairs depicted in 'comps.xlsx' file")
    #read in list of comparisons
    comps <- read.xlsx(compsfile)
  }

  #perform PoissonSeq
  message("Differential expression analysis...")
  pseqdata <- finalimpintensitiesIRS
  newwb <- createWorkbook()
  newwb2 <- createWorkbook()
  for (i in 1:dim(comps)[1]){
    #get the intensities for this comparison
    pseqdata |>
      select(matches(comps[i,1]),matches(comps[i,2])) |>
      na.omit() -> pdata

    #make indicator variable y
    y <- c(rep(1,ncol(pdata|>select(matches(comps[i,1])))),
             rep(2,ncol(pdata|>select(matches(comps[i,2])))))
    
    #perform PSeq
    pseq<-PS.Main(dat=list(n=pdata,y=y,type="twoclass",pair=FALSE,gname=row.names(pdata)),para=list(ct.sum=0,ct.mean=0))
    
    # Let's create results table. Here, fold-change is calculated from previously normalized intensities.
    pseq |>
      select(2,4,5) |>
      rename(UID = gname) |>
      inner_join(nozeros |>
                   select(Proteins:original.id) |>
                   rownames_to_column(var="UID"), by = "UID") |>
      inner_join(pdata |>
                   mutate(UID = rownames(pdata)), by = "UID") %>%
      mutate(log2FC = log2(rowMeans(across(starts_with(comps[i,2])))/rowMeans(across(starts_with(comps[i,1])))),
             GeneID = gsub(pattern = ";(.+)", replacement = "", x = .[]$Proteins, perl = T), .after = fdr) -> myresults
    
    #save
    mycomp = paste0(comps[i,2],"_vs_",comps[i,1])
    if (nchar(mycomp)>31){
      mysheet = paste0(abbreviate(comps[i,2],minlength=15),"_vs_",abbreviate(comps[i,1],minlength=15))
    }else{
      mysheet=paste0(comps[i,2],"_vs_",comps[i,1])
    }
    addWorksheet(wb = newwb2, sheetName = mysheet, gridLines = TRUE)
    writeDataTable(wb=newwb2, sheet=mysheet,x=myresults,tableStyle="none",
                   rowNames=FALSE,withFilter=FALSE,
                   bandedRows=FALSE,bandedCols=FALSE)
    
    #table formatting for MA and volcano plots
    myresults %>%
      mutate(mean_mock = rowMeans(select(pdata,matches(comps[i,1]))),
             mean_treatment = rowMeans(select(pdata, matches(comps[i,2])))) -> myData
    myData$DE <- "NS"
    myData$DE[myData$pval<=qval & myData$log2FC < 0] <- "Down"
    myData$DE[myData$pval<=qval & myData$log2FC > 0] <- "Up"
    myData$DE <- factor(myData$DE, levels = c("Up", "Down", "NS"))
    
    #make MA plot
    if (stat=="q"){
      MA <- ggplot(myData,
                   aes(x = log2(mean_mock*mean_treatment)*0.5,
                       y = log2FC))+
        geom_point(color=alpha('black', 0.3), shape=21, size=2, aes(fill=factor(DE))) +
        scale_fill_manual(values=alpha(c('#FD6467','#56B4E9','#FDFD96'),0.8)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = alpha("black", 0.7), linewidth = 1)+
        ggtitle(paste0(comps[i,2],"/",comps[i,1]," (",sum(pseq$fdr<qval)," DE elements)"))+
        xlab(label = bquote("A (Average"~log[2]~"Intensity)"))+
        ylab(label = bquote("M ("*log[2]~"Fold-change)"))+
        labs(fill = element_blank())+
        theme_classic()+
        theme(axis.title.x = element_text(size =20),
              axis.title.y = element_text(size =20),
              legend.text = element_text(size=15),
              axis.text.x= element_text(size=15),
              axis.text.y= element_text(size=15))
      
      png(filename = paste0(comps[i,2],"_vs_",comps[i,1],"_MA_plot_",qval,".png"), width = 2500, height = 1500, res = 300)
      plot(MA)
      dev.off()
      
    } else {
      
      MA <- ggplot(myData,
                   aes(x = log2(mean_mock*mean_treatment)*0.5,
                       y = log2FC))+
        geom_point(color=alpha('black', 0.3), shape=21, size=2, aes(fill=factor(DE))) +
        scale_fill_manual(values=alpha(c('#FD6467','#56B4E9','#FDFD96'),0.8)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = alpha("black", 0.7), linewidth = 1)+
        ggtitle(paste0(comps[i,2],"/",comps[i,1]," (",sum(pseq$pval<qval)," DE elements)"))+
        xlab(label = bquote("A (Average"~log[2]~"Intensity)"))+
        ylab(label = bquote("M ("*log[2]~"Fold-change)"))+
        labs(fill = element_blank())+
        theme_classic()+
        theme(axis.title.x = element_text(size =20),
              axis.title.y = element_text(size =20),
              legend.text = element_text(size=15),
              axis.text.x= element_text(size=15),
              axis.text.y= element_text(size=15))
      
      png(filename = paste0(comps[i,2],"_vs_",comps[i,1],"_MA_plot_",qval,".png"), width = 2500, height = 1500, res = 300)
      plot(MA)
      dev.off()
    }
    #make volcano plot
    signum = sum(pseq$pval<qval)
    if (stat=="q"){
      volcanoPlot <- ggplot(myData, aes(x = log2FC, y = -log10(fdr)))+
        geom_point(color='black', shape=21, size=2, alpha=0.7, aes(fill=factor(DE))) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1, alpha = 0.7)+
        scale_fill_manual(values=c('#FD6467','#56B4E9','#FDFD96' ))+
        ggtitle(paste0(comps[i,2],"/",comps[i,1]," (",sum(pseq$fdr<qval)," DE elements)"))+
        xlab(label = bquote(log[2]~"Fold-change"))+
        ylab(label = bquote(-log[10]~"qvalue"))+
        labs(fill = element_blank())+
        ylim(0,3)+
        xlim(-3,3)+
        theme_classic()+
        theme(legend.position = "top")
      
      png(filename=paste0(comps[i,2],"_vs_",comps[i,1],"_volcano_plot_",qval,".png"),width=2500,height=2000,res=300)
      print(volcanoPlot)
      dev.off()
      
    }else{
      
      volcanoPlot <- ggplot(myData, aes(x = log2FC, y = -log10(pval)))+
        geom_point(color='black', shape=21, size=2, alpha=0.7, aes(fill=factor(DE))) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1, alpha = 0.7)+
        scale_fill_manual(values=c('#FD6467','#56B4E9','#FDFD96' ))+
        ggtitle(paste0(comps[i,2],"/",comps[i,1]," (",sum(pseq$pval<qval)," DE elements)"))+
        xlab(label = bquote(log[2]~"Fold-change"))+
        ylab(label = bquote(-log[10]~"pvalue"))+
        labs(fill = element_blank())+
        ylim(0,3)+
        xlim(-3,3)+
        theme_classic()+
        theme(legend.position = "top")
  
      png(filename=paste0(comps[i,2],"_vs_",comps[i,1],"_volcano_plot_",qval,".png"),width=2500,height=2000,res=300)
      plot(volcanoPlot)
      dev.off()
    }
    
    #make pvalue and qvalue histogram

      png(filename=paste0(comps[i,2],"_vs_",comps[i,1],"_qval_hist.png"),width=2000,height=2000,res=300)
      h <- hist(x=pseq$fdr,breaks=100)
      plot(h,main="q-value distribution histogram", xlab="q-value")
      dev.off()

      png(filename=paste0(comps[i,2],"_vs_",comps[i,1],"_pval_hist.png"),width=2000,height=2000,res=300)
      h <- hist(x=pseq$pval,breaks=100)
      plot(h,main="p-value distribution histogram",xlab="p-value")
      dev.off()

    #get differentially expressed genes and save
    if (stat=="q"){
      myresults %>%
        filter(fdr<qval, .preserve = TRUE) -> mypros
    }else{
      myresults %>%
        filter(pval<qval, .preserve = TRUE) -> mypros
    }

    #save
    addWorksheet(wb = newwb, sheetName = mysheet, gridLines = TRUE)
    writeDataTable(wb=newwb, sheet=mysheet,x=mypros,tableStyle="none",
                   rowNames=FALSE,withFilter=FALSE,
                   bandedRows=FALSE,bandedCols=FALSE)
  }
  
  #write workbook
  if (stat=="q"){
    saveWorkbook(newwb, paste("Pseq_all_comps_q",qval,".xlsx",sep=""),overwrite=TRUE)
  }else{
    saveWorkbook(newwb, paste("Pseq_all_comps_p",qval,".xlsx",sep=""),overwrite=TRUE)
  }
  saveWorkbook(newwb2, "Pseq_all_comps.xlsx",overwrite=TRUE)
}else{
  #save expression values with annotation information
  finalimpintensitiesIRS = finalimpintensitiesIRS[order(row.names(finalimpintensitiesIRS)),]
  nozeros = nozeros[order(row.names(nozeros)),]
  myresults = data.frame(nozeros[row.names(nozeros)%in%row.names(finalimpintensitiesIRS),1:dim(mydata)[2]],finalimpintensitiesIRS)
  write.csv(myresults,'Normalized_values.csv')
}

message("Finished! Please close TMT-NEAT window")
}
