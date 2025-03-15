library(stringr)
library(dplyr)
library(tidyr)
library(tibble)
library(viridis)
# Modifies dataframe genotypePairs with columns genotype1 and genotype2 at a minimum that identify all pairwise comparisons between genotypes
# Remove dashes from genotype names before using so we can split the comparisons in the tukey step
getSignificantCrossovers <- function(data, pheno, environment, genotypePairs)
{
  phenotype <- paste0(pheno)
  phenotypeMean <- paste0(pheno, 'Mean')
  phenotypeRank <- paste0(pheno, 'Rank')
  phenotypeAdjustedP <- paste0(pheno, 'AdjP')
  phenotypeSigDiff <- paste0(pheno, 'SigDiff')
  phenotypeRankChange <- paste0(pheno, 'RC')
  phenotypeScore <- paste0(pheno, 'Score')
  phenotypeComparedEnvs <- paste0(pheno, 'ComparedEnvs')
  envs <- unique(data[[environment]])
  totalEnvironments <- length(envs)
  for(env in envs)
  {
    envSuffix <- paste0('.E', env)

    environmentData <- data %>%
      filter(.data[[environment]]==env) %>%
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
      filter(.data[[environment]]==env) %>%
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

  # cols <- colnames(genotypePairs)
  for(i in 1:(totalEnvironments - 1))
  {
    envI <- envs[i]
    if(is.na(envI)){next}
    envISuffix <- paste0('.E', envI)
    envIG1Rank <- paste0(phenotypeRank, envISuffix, '.G1')
    envIG2Rank <- paste0(phenotypeRank, envISuffix, '.G2')
    envISigDiff <- paste0(phenotypeSigDiff, envISuffix)

    # if(length(setdiff(c(envISigDiff, envIG1Rank, envIG2Rank), cols)) > 0){next}

    for(j in (i + 1):totalEnvironments)
    {
      envJ <- envs[j]
      if(is.na(envJ)){next}
      envJSuffix <- paste0('.E', envJ)
      envJG1Rank <- paste0(phenotypeRank, envJSuffix, '.G1')
      envJG2Rank <- paste0(phenotypeRank, envJSuffix, '.G2')
      envJSigDiff <- paste0(phenotypeSigDiff, envJSuffix)

      # if(length(setdiff(c(envJSigDiff, envJG1Rank, envJG2Rank), cols)) > 0){next}

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


plotInteractionImportanceGrid <- function(significantInteractionsData = sigCrossovers, 
                                          performanceData = hybrids, trait, traitLabel, 
                                          legendTitle = str_wrap('Normalized Interaction Importance Score', 25),
                                          legendPosition = 'right', legendTextAngle = 0, 
                                          legendTextHJust = 1, xAxisLabelAngle = 0)
{
  phenotype <- trait
  phenotypeSpatial <- paste0(phenotype)
  phenotypeLabel <- traitLabel
  phenotypeScoreNormalized <- paste0(phenotype, 'ScoreNormalized')
  
  df1 <- significantInteractionsData %>%
    select(c(genotype1, genotype2, all_of(phenotypeScoreNormalized)))
  df2 <- df1 %>%
    rowwise() %>%
    mutate(genotype1A = genotype2, 
           genotype2 = genotype1) %>%
    select(c(genotype1A, genotype2, all_of(phenotypeScoreNormalized))) %>%
    rename(genotype1 = genotype1A)
  
  df <-bind_rows(df1, df2) %>% 
    filter(!is.na(.data[[phenotypeScoreNormalized]]))
  
  blues <- lm(as.formula(paste0(phenotypeSpatial, ' ~ environment + genotype')), data = performanceData) 
  blues <- blues$coefficients
  blues <- as_tibble(blues, rownames = 'genotype') %>%
    filter(str_detect(genotype, 'genotype')) %>%
    rename(blue = value) %>%
    mutate(rank = dense_rank(blue)) %>%
    rowwise() %>% 
    mutate(genotype = str_remove(genotype, 'genotype')) %>%
    select(genotype, rank)
  
  df <- left_join(df, blues, join_by(genotype1==genotype), keep = FALSE, suffix = c('', ''), relationship = 'many-to-one') %>%
    rename(rankG1 = rank) %>%
    left_join(blues, join_by(genotype2==genotype), keep = FALSE, suffix = c('', '')) %>%
    rename(rankG2 = rank) %>%
    select(c(genotype1, genotype2, rankG1, rankG2, all_of(phenotypeScoreNormalized)))
  
  heatmap <- ggplot(df, aes(rankG1, rankG2, fill = .data[[phenotypeScoreNormalized]])) + 
    geom_tile() + 
    scale_fill_viridis(direction = -1) +
    labs(x = 'Hybrid Rank', y = 'Hybrid Rank', fill = legendTitle, title = phenotypeLabel) + 
    theme_minimal() +
    theme(text = element_text(color = 'black', size = 9),
          axis.text.x = element_text(color = 'black', size = 9, angle = xAxisLabelAngle),
          axis.text = element_text(color = 'black', size = 9),
          legend.text = element_text(color = 'black', size = 9, angle = legendTextAngle, hjust = legendTextHJust, vjust = 1),
          axis.line = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_blank(), 
          legend.position = legendPosition,
          legend.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 9))
  return(heatmap)
}

