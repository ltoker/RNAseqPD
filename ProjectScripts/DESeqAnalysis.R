source("SetUp.R")
packageF("DESeq2")

name = "Parkome"
load(paste0(GeneralResultsPath, "Parkome.RData"))

MetaCovar = "sex + age_years + rin + PMI_hours + Batch"
Model = as.formula(paste0("~Profile +", MetaCovar))

Aspects = c("B", "M", "C")

source("ProjectScripts/ProjectFunctions.R")

BeforeMGPcorrection <- lapply(studyFinal, function(Cohort){
  Meta <- Cohort$Metadata 
  Counts = Cohort$countMatrix[,-1]
  DESeqOut <- DESeq2RUN(data =  Counts, Meta = Meta, model = Model)
  DESeqResults <- GetDESeq2Results(DESeqOut, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)
  DESeqResultsDF <- GetOneSidedPval(ResultsObj = DESeqResults)
  EnrichListPDdown <- gsr(scores = DESeqResultsDF, scoreColumn = "DownPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = Aspects)
  EnrichListPDup <- gsr(scores = DESeqResultsDF, scoreColumn = "UpPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = Aspects)
  list(DESeqOut = DESeqOut,
       DESeqResultsDF = DESeqResultsDF,
       EnrichListPDdown = EnrichListPDdown,
       EnrichListPDup = EnrichListPDup)
})


#Repeat but correcting just for Oligo MGP
OligModel <- as.formula(paste0("~Profile + ", MetaCovar, " + Oligo_Genes"))
MicrogliaModel <- as.formula(paste0("~Profile + ", MetaCovar, " + Microglia_Genes"))
OligMicrogliaModel <- as.formula(paste0("~Profile + ", MetaCovar, " + Microglia_Genes + Oligo_Genes"))

AfterOligocorrection <- lapply(studyFinal, function(Cohort){
  Meta <- Cohort$Metadata 
  Counts = Cohort$countMatrix[,-1]
  DESeqOut <- DESeq2RUN(data =  Counts, Meta = Meta, model = OligModel)
  DESeqResults <- GetDESeq2Results(DESeqOut, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)
  
  DESeqResultsDF <- GetOneSidedPval(ResultsObj = DESeqResults)
  EnrichListPDdown <- gsr(scores = DESeqResultsDF, scoreColumn = "DownPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = Aspects)
  EnrichListPDup <- gsr(scores = DESeqResultsDF, scoreColumn = "UpPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = Aspects)
  list(DESeqOut = DESeqOut,
       DESeqResultsDF = DESeqResultsDF,
       EnrichListPDdown = EnrichListPDdown,
       EnrichListPDup = EnrichListPDup)
})

AfterMicrogliacorrection <- lapply(studyFinal, function(Cohort){
  Meta <- Cohort$Metadata 
  Counts = Cohort$countMatrix[,-1]
  DESeqOut <- DESeq2RUN(data =  Counts, Meta = Meta, model = MicrogliaModel)
  DESeqResults <- GetDESeq2Results(DESeqOut, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)
  
  DESeqResultsDF <- GetOneSidedPval(ResultsObj = DESeqResults)
  EnrichListPDdown <- gsr(scores = DESeqResultsDF, scoreColumn = "DownPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = Aspects)
  EnrichListPDup <- gsr(scores = DESeqResultsDF, scoreColumn = "UpPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = Aspects)
  list(DESeqOut = DESeqOut,
       DESeqResultsDF = DESeqResultsDF,
       EnrichListPDdown = EnrichListPDdown,
       EnrichListPDup = EnrichListPDup)
})

AfterOligMicrogliacorrection <- lapply(studyFinal, function(Cohort){
  Meta <- Cohort$Metadata 
  Counts = Cohort$countMatrix[,-1]
  DESeqOut <- DESeq2RUN(data =  Counts, Meta = Meta, model = OligMicrogliaModel)
  DESeqResults <- GetDESeq2Results(DESeqOut, coef = "Profile_PD_vs_Cont", alpha = 0.05, indepFilter = TRUE)
  
  DESeqResultsDF <- GetOneSidedPval(ResultsObj = DESeqResults)
  EnrichListPDdown <- gsr(scores = DESeqResultsDF, scoreColumn = "DownPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = Aspects)
  EnrichListPDup <- gsr(scores = DESeqResultsDF, scoreColumn = "UpPval", bigIsBetter = F, logTrans = T, annotation = GenericHumanAnno, aspects = Aspects)
  list(DESeqOut = DESeqOut,
       DESeqResultsDF = DESeqResultsDF,
       EnrichListPDdown = EnrichListPDdown,
       EnrichListPDup = EnrichListPDup)
})
#Compare results from both cohorts before Olig adjustment
CompareResultsAll(data1 = BeforeMGPcorrection$NBB$DESeqResultsDF, data2 = BeforeMGPcorrection$Norway$DESeqResultsDF, name1 = "NBB",
               name2 = "Norway", colorCol = "Cohort",
               Title = "Before Oligo correction")


#Compare results from both cohorts after Olig adjustment
CompareResultsAll(data1 = AfterOligMicrogliacorrection$NBB$DESeqResultsDF, data2 = AfterOligocorrection$Norway$DESeqResultsDF, name1 = "NBB",
               name2 = "Norway", colorCol = "Cohort",
               Title = "After Microglia/Oligo correction")



rm(datas, cpmCountFiltered, cpmCountFiltered, ensembl, estimates, ExpDataAll, ExpDataCPM)


########Get adjusted log2 counts for genes #####################
AdjCovarPV <- data.frame(Cov = c("sex_M_vs_F", 
                                 "Batch_2_vs_1", "Batch_3_vs_1", "Batch_4_vs_1",
                                 "age_years","PMI_hours", "rin",
                                 "Oligo_Genes"),
                         adjType = c(rep("base", 4), rep("mean", 4)))

AdjCovarNBB <- data.frame(Cov = c("sex_M_vs_F", "Batch_1_vs_0",
                                  "Batch_2_vs_0", "Batch_3_vs_0", "Batch_4_vs_0",
                                  "age_years","PMI_hours", "rin",
                                  "Oligo_Genes", "Microglia_Genes"),
                          adjType = c(rep("base", 5), rep("mean", 5)))


PlotAdjCounts <- function(Gene){
  EnsemblID = BeforeMGPcorrection$Norway$DESeqResultsDF %>% filter(GeneSymbol == Gene) %>% .$EnsemblID
  PV <- plotCounts(BeforeMGPcorrection$Norway$DESeqOut, gene = EnsemblID, intgroup = "Profile", returnData = T) %>% data.frame
  PV %<>% mutate(Adj_MGP = GetAdjCountDESeq(dds = AfterOligocorrection$Norway$DESeqOut, Gene = EnsemblID, adjCov = AdjCovarPV),
                       Adj = GetAdjCountDESeq(dds = BeforeMGPcorrection$Norway$DESeqOut, Gene = EnsemblID, adjCov = AdjCovarPV %>% filter(!Cov %in%  c("Oligo_Genes"))),
                       count = log2(count),
                       CommonName = rownames(attr(AfterOligocorrection$Norway$DESeqOut, "modelMatrix")),
                       Cohort = "PW")
  PV$RNA2 <- studyFinal$Norway$Metadata$RNA2[match(PV$CommonName, studyFinal$Norway$Metadata$CommonName)]
  
  
  
  NBB <- plotCounts(BeforeMGPcorrection$NBB$DESeqOut, gene = EnsemblID, intgroup = "Profile", returnData = T) %>% data.frame
  NBB %<>% mutate(Adj_MGP = GetAdjCountDESeq(dds = AfterOligMicrogliacorrection$NBB$DESeqOut, Gene = EnsemblID, adjCov = AdjCovarNBB),
                        Adj = GetAdjCountDESeq(dds = BeforeMGPcorrection$NBB$DESeqOut, Gene = EnsemblID, adjCov = AdjCovarNBB %>% filter(!Cov %in%  c("Oligo_Genes", "Microglia_Genes"))),
                        count = log2(count),
                        CommonName = rownames(attr(AfterOligMicrogliacorrection$NBB$DESeqOut, "modelMatrix")),
                        Cohort = "NBB")
  NBB$RNA2 <- studyFinal$NBB$Metadata$RNA2[match(NBB$CommonName, studyFinal$NBB$Metadata$CommonName)]
  

  
  
  Gene_All <- rbind(PV, NBB) %>% data.frame() %>% gather(key = "CountType", value = "GeneCount", Adj_MGP, Adj)
  
  Gene_All$condition <- relevel(Gene_All$Profile,ref = "Cont")
  
  levels(Gene_All$condition) <- c("Cont", "PD")
  
  Gene_All$CountType <- factor(Gene_All$CountType, levels = c("Adj", "Adj_MGP"))
  levels(Gene_All$CountType) <- c("NotMGPAdjusted", "MGPAdjusted")
  Gene_All$Cohort <- factor(Gene_All$Cohort, levels = c("PW", "NBB", "Dimitriu")) 
  
  Plot <- ggplot(Gene_All, aes(condition, GeneCount, color = condition)) +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_blank()) +
    labs(x = "", y = paste0("log2(",Gene, ")")) +
    geom_boxplot(outlier.shape = NA, aes(fill = condition), alpha = 0.5) +
    geom_jitter(width = 0.2) +
    scale_color_manual(values =  c("dodgerblue4", "chocolate1"), name = "Group") +
    scale_fill_manual(values =  c("dodgerblue4", "chocolate1"), name = "Group") +
    facet_grid(CountType~Cohort, scales = "free_y")
  ggsave(paste0(Gene,"geneLevels.pdf"), plot = Plot, device = "pdf", width =  7, height =5, dpi = 300, useDingbats=F)

  Plot <- ggplot(Gene_All %>% filter(CountType == "MGPAdjusted"), aes(condition, GeneCount, color = condition)) +
    theme_bw(base_size = 14) +
    theme(panel.grid = element_blank()) +
    labs(x = "", y = paste0("log2(",Gene, ")")) +
    geom_boxplot(outlier.shape = NA, aes(fill = condition), alpha = 0.5) +
    geom_jitter(width = 0.2) +
    scale_color_manual(values =  c("dodgerblue4", "chocolate1"), name = "Group") +
    scale_fill_manual(values =  c("dodgerblue4", "chocolate1"), name = "Group") +
    facet_wrap(~Cohort, scales = "free")
  ggsave(paste0(Gene,"geneLevelsMGPonly.pdf"), plot = Plot, device = "pdf", width =  7, height =3, dpi = 300, useDingbats=F)
  print(Plot)
}

PlotAdjCounts("PTPRH")
PlotAdjCounts("JUP")
PlotAdjCounts("DLG2")
PlotAdjCounts("KCNH1")
PlotAdjCounts("MAP4K4")
PlotAdjCounts("FRMD5")
PlotAdjCounts("CNTNAP2")

#Get the corrected counts for all genes
AdjustedPV <- sapply(AfterOligocorrection$Norway$DESeqResultsDF$GeneSymbol, function(Gene){
  EnsemblID = AfterOligocorrection$Norway$DESeqResultsDF %>% filter(GeneSymbol == Gene) %>% .$EnsemblID
  GetAdjCountDESeq(dds = AfterOligocorrection$Norway$DESeqOut, Gene = EnsemblID, adjCov = AdjCovarPV)
}) %>% t %>% data.frame

sampleNamesPV <- rownames(attr(AfterOligocorrection$Norway$DESeqOut, "modelMatrix"))
names(AdjustedPV) <- studyFinal$Norway$Metadata$RNA2[match(sampleNamesPV, studyFinal$Norway$Metadata$CommonName)]
AdjustedPV %<>% mutate(GeneSymbol = rownames(.))


AdjustedNBB <- sapply(AfterOligMicrogliacorrection$NBB$DESeqResultsDF$GeneSymbol, function(Gene){
  EnsemblID = AfterOligMicrogliacorrection$NBB$DESeqResultsDF %>% filter(GeneSymbol == Gene) %>% .$EnsemblID
  GetAdjCountDESeq(dds = AfterOligMicrogliacorrection$NBB$DESeqOut, Gene = EnsemblID, adjCov = AdjCovarNBB)
}) %>% t %>% data.frame

sampleNamesNBB <- rownames(attr(AfterOligMicrogliacorrection$NBB$DESeqOut, "modelMatrix"))
names(AdjustedNBB) <- studyFinal$NBB$Metadata$RNA2[match(sampleNamesNBB, studyFinal$NBB$Metadata$CommonName)]
AdjustedNBB %<>% mutate(GeneSymbol = rownames(.))

saveRDS(list(AdjustedPV = AdjustedPV,
             AdjustedNBB = AdjustedNBB), file = "AdjestedCounts.Rds")


save.image(file = paste0(GeneralResultsPath, "DESeqAnalysisParkome.Rdata"))
