# Script to run interaction importance analysis on HCC
library(tidyverse)
setwd('/work/schnablelab/jdavis132/HIPS')
# Modifies dataframe genotypePairs with columns genotype1 and genotype2 at a minimum that identify all pairwise comparisons between genotypes
# Remove dashes from genotype names before using so we can split the comparisons in the tukey step
getSignificantCrossovers <- function(data, pheno, environments)
{
  phenotype <- paste0(pheno, '.sp')
  phenotypeMean <- paste0(pheno, 'Mean')
  phenotypeRank <- paste0(pheno, 'Rank')
  phenotypeAdjustedP <- paste0(pheno, 'AdjP')
  phenotypeSigDiff <- paste0(pheno, 'SigDiff')
  phenotypeRankChange <- paste0(pheno, 'RC')
  phenotypeScore <- paste0(pheno, 'Score')
  phenotypeComparedEnvs <- paste0(pheno, 'ComparedEnvs')
  
  for(env in environments)
  {
    envSuffix <- paste0('.E', env)
    
    environmentData <- data %>%
      filter(environmentCode==env) %>%
      select(genotype, all_of(phenotype))
    environmentData <- environmentData[complete.cases(environmentData), ]
    
    if(length(environmentData[[phenotype]]) < 1|length(unique(environmentData$genotype)) < 1){next}
    
    anova <- aov(as.formula(paste(phenotype, ' ~ genotype')), data = environmentData)
    
    tukey <- TukeyHSD(anova)$genotype %>%
      as_tibble(rownames = 'genotypes') %>%
      rowwise() %>%
      mutate(genotype1 = str_split_i(genotypes, '-', 1),
             genotype2 = str_split_i(genotypes, '-', 2)) %>%
      rename('{phenotypeAdjustedP}' := `p adj`) %>%
      mutate('{phenotypeSigDiff}' := .data[[phenotypeAdjustedP]] < 0.05) %>%
      select(c(genotypes, genotype1, genotype2,  all_of(c(phenotypeAdjustedP, phenotypeSigDiff))))
    
    environmentDataSummary <- data %>%
      filter(environmentCode==env) %>%
      group_by(genotype) %>%
      summarise('{phenotypeMean}' := mean(.data[[phenotype]], na.rm = TRUE)) %>%
      mutate('{phenotypeRank}' := dense_rank(desc(.data[[phenotypeMean]]))) %>%
      select(c(all_of(phenotypeRank), genotype))
    
    envG1Suffix <- paste0(envSuffix, '.G1')
    envG2Suffix <- paste0(envSuffix, '.G2')
    
    genotypePairs <- left_join(genotypePairs, tukey, join_by(genotype1==genotype1, genotype2==genotype2), keep = FALSE, suffix = c('', '')) %>%
      left_join(tukey, join_by(genotype1==genotype2, genotype2==genotype1), keep = FALSE, suffix = c('', '.t2')) %>%
      rowwise() %>%
      mutate('{phenotypeAdjustedP}' := case_when(!is.na(.data[[phenotypeAdjustedP]]) ~ .data[[phenotypeAdjustedP]],
                                                 !is.na(.data[[paste0(phenotypeAdjustedP, '.t2')]]) ~ .data[[paste0(phenotypeAdjustedP, '.t2')]]),
             '{phenotypeSigDiff}' := case_when(!is.na(.data[[phenotypeSigDiff]]) ~ .data[[phenotypeSigDiff]],
                                               !is.na(.data[[paste0(phenotypeSigDiff, '.t2')]]) ~ .data[[paste0(phenotypeSigDiff, '.t2')]])) %>%
      distinct(genotype1, genotype2, .keep_all = TRUE) %>%
      select(!ends_with('.t2')) %>%
      rename('{phenotypeAdjustedP}{envSuffix}' := .data[[phenotypeAdjustedP]],
             '{phenotypeSigDiff}{envSuffix}' := .data[[phenotypeSigDiff]]) %>%
      full_join(environmentDataSummary, join_by(genotype1==genotype), keep = FALSE, suffix = c('', ''), relationship = 'many-to-one') %>%
      rename('{phenotypeRank}{envG1Suffix}' := .data[[phenotypeRank]]) %>%
      full_join(environmentDataSummary, join_by(genotype2==genotype), keep = FALSE, suffix = c('', ''), relationship = 'many-to-one') %>%
      rename('{phenotypeRank}{envG2Suffix}' := .data[[phenotypeRank]]) %>%
      filter(!is.na(genotype1) & !is.na(genotype2))
  }
  
  cols <- colnames(genotypePairs)
  for(i in 1:(totalEnvironments - 1))
  {
    envI <- environments[i]
    if(is.na(envI)){next} 
    envISuffix <- paste0('.E', envI)
    envIG1Rank <- paste0(phenotypeRank, envISuffix, '.G1')
    envIG2Rank <- paste0(phenotypeRank, envISuffix, '.G2')
    envISigDiff <- paste0(phenotypeSigDiff, envISuffix)
    
    if(length(setdiff(c(envISigDiff, envIG1Rank, envIG2Rank), cols)) > 0){next}
    
    for(j in (i + 1):totalEnvironments)
    {
      envJ <- environments[j]
      if(is.na(envJ)){next} 
      envJSuffix <- paste0('.E', envJ)
      envJG1Rank <- paste0(phenotypeRank, envJSuffix, '.G1')
      envJG2Rank <- paste0(phenotypeRank, envJSuffix, '.G2')
      envJSigDiff <- paste0(phenotypeSigDiff, envJSuffix)
      
      if(length(setdiff(c(envJSigDiff, envJG1Rank, envJG2Rank), cols)) > 0){next}
      
      envPairSuffix <- paste0('.E', envI, '-', envJ)
      envPairRankChange <- paste0(phenotypeRankChange, envPairSuffix)
      envPairScore <- paste0(phenotypeScore, envPairSuffix)
      
      genotypePairs <- genotypePairs %>%
        rowwise() %>%
        mutate('{envPairRankChange}' := ((.data[[envIG1Rank]] - .data[[envIG2Rank]])/(.data[[envJG1Rank]] - .data[[envJG2Rank]])) < 0) %>%
        mutate('{envPairScore}' := case_when(!.data[[envPairRankChange]] ~ 0, 
                                             .data[[envPairRankChange]] ~ .data[[envISigDiff]] + .data[[envJSigDiff]]))
    }
  }
  
  genotypePairs <- genotypePairs %>%
    rowwise() %>%
    mutate('{phenotypeScore}' := rowSums(across(contains(phenotypeScore)), na.rm = TRUE), 
           '{pheno}ComparedEnvs' := rowSums(!is.na(across(contains(phenotypeAdjustedP))), na.rm = TRUE)) %>%
    mutate('{phenotypeScore}Normalized' := .data[[phenotypeScore]]/((.data[[phenotypeComparedEnvs]]*(.data[[phenotypeComparedEnvs]] - 1))))
  return(genotypePairs)
}

hybrids <- hybrids <- read.csv('HYBRIDS_2022_2023_SPATIALLYCORRECTED.csv') %>%
  filter(location!='') %>% 
  mutate(nitrogenTreatment = factor(nitrogenTreatment, levels = c('Low', 'Medium', 'High'))) %>%
  rowwise() %>%
  mutate(across(where(is.numeric), ~case_when(.==-Inf ~ NA, .default = .)))

envsPerHybrid <- tibble(hybrid = unique(hybrids$genotype), numEnvs = NULL)
for(i in 1:length(unique(envsPerHybrid$hybrid)))
{
  hybridData <- filter(hybrids, genotype==envsPerHybrid$hybrid[i])
  envsPerHybrid$numEnvs[i] <- length(unique(hybridData$environment))
}

singleEnvHybrids <- envsPerHybrid$hybrid[envsPerHybrid$numEnvs<4]
hybrids <- filter(hybrids, !(genotype %in% singleEnvHybrids)) %>%
  rowwise() %>%
  mutate(genotype = str_remove_all(genotype, '-'))

phenotypes <- c("combineTestWeight", "combineMoisture", "flagLeafHeight", "earHeight", "yieldPerAcre", 
                'GDDToAnthesis', 'GDDToSilk', 'anthesisSilkingIntervalGDD', 'kernelRowNumber', 'earWidth',
                'earLength', 'shelledCobWidth', 'shelledCobMass', 'kernelMassPerEar', 'kernelsPerEar', 'hundredKernelMass',
                'earFillLength', 'kernelsPerRow')
#phenotypes <- c('yieldPerAcre')
# How often are interactions between a pair of hybrids crossover interactions AND represent significant differences in the phenotype?
genotypePairs <- tibble(genotype1 = NULL, genotype2 = NULL)
hybridEnvs <- hybrids %>%
  group_by(environment) %>%
  summarise(environmentCode = cur_group_id())
hybrids <- full_join(hybrids, hybridEnvs, join_by(environment), keep = FALSE, suffix = c('', ''))

hybridsSigCrossovers <- hybrids %>%
  rowwise() %>%
  mutate(genotype = str_replace_all(genotype, '-', ' ')) %>%
  filter(!is.na(genotype))

allGenotypes <- unique(hybridsSigCrossovers$genotype)
totalGenotypes <- length(allGenotypes)

for(i in 1:(totalGenotypes - 1))
{
  df <- tibble(genotype1 = allGenotypes[i], genotype2 = allGenotypes[(i + 1):totalGenotypes])
  genotypePairs <- bind_rows(genotypePairs, df)
}

environments <- unique(hybrids$environmentCode)
totalEnvironments <- length(environments)

for(i in phenotypes)
{
  genotypePairs <- full_join(genotypePairs,
                             getSignificantCrossovers(hybridsSigCrossovers, i, environments),
                             join_by(genotype1, genotype2),
                             keep = FALSE,
                             suffix = c('', ''))
}

genotypePairs <- genotypePairs %>%
  rowwise() %>%
  mutate(genotype1 = case_when(genotype1=='SYNGENTA NK06593120EZ1' ~ 'SYNGENTA NK0659-3120-EZ1',
                               genotype1=='SYNGENTA NK07603111' ~ 'SYNGENTA NK0760-3111', 
                               .default = genotype1),
         genotype2 = case_when(genotype2=='SYNGENTA NK06593120EZ1' ~ 'SYNGENTA NK0659-3120-EZ1',
                               genotype2=='SYNGENTA NK07603111' ~ 'SYNGENTA NK0760-3111', 
                               .default = genotype2))

write.csv(genotypePairs, 'significantCrossovers.csv')
