
##############################
#### Analysis in MSstats
##############################

##############################
## Load MSstats package
##############################
library(MSstats)

##############################
## Read Skyline report
##############################
raw <- read.csv(file="Choi2017_DDA_Skyline_input.csv")
head(raw)

annotation <- read.csv('Choi2017_DDA_Skyline_annotation.csv')
annotation

##############################
## Make MSstats required format
##############################
quant <- SkylinetoMSstatsFormat(raw, 
                                annotation = annotation,
                                removeProtein_with1Feature = TRUE)

head(quant)


##############################
## dataProcess
## including Normalization, decide censored cutoff, protein-level summarization
##############################

processed.quant <- dataProcess(quant,
                               normalization = 'equalizeMedians',
                               summaryMethod="TMP",
                               cutoffCensored="minFeature",
                               censoredInt="0",
                               MBimpute=TRUE,
                               maxQuantileforCensored=0.999)

save(processed.quant, file='processed.quant.rda')

##############################
## Data visualization
##############################

dataProcessPlots(processed.quant, type="QCplot", 
                 ylimDown=0, 
                 which.Protein = 'allonly',
                 width=7, height=7,  
                 address="Choi2017_DDA_Skyline_")

dataProcessPlots(processed.quant, type="Profileplot", 
                 ylimDown=0, 
                 originalPlot = TRUE,
                 summaryPlot = TRUE,
                 width=7, height=7,  
                 address="Choi2017_DDA_Skyline_")

dataProcessPlots(processed.quant, type="Conditionplot", 
                 ylimDown=0, 
                 width=7, height=7,  
                 address="Choi2017_DDA_Skyline_")


##############################
## Model-based comparison + adjust p-value
##############################

comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)
comparison <- rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison)<-c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")


test.MSstats <- groupComparison(contrast.matrix=comparison, data=processed.quant)
test.MSstats <- test.MSstats$ComparisonResult

##############################
## save the result
##############################

save(test.MSstats, file='test.MSstats.rda')
write.csv(test.MSstats, file='Choi2017_DDA_Skyline_testResult_byMSstats.csv')


##############################
## Visualization of result
##############################
groupComparisonPlots(data=test.MSstats, type="VolcanoPlot",
                     width=6, height=6,
                     address="Choi2017_DDA_Skyline_")

groupComparisonPlots(data=test.MSstats, type="ComparisonPlot",
                     width=6, height=6,
                     address="Choi2017_DDA_Skyline_")

