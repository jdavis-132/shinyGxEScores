library(shiny)
library(DT)
library(readr)
source('R/Functions.R')

shinyServer
(
  function(input, output, session)
  {
    output$introNote <- renderText({ "Determine which empirical GxE interactiopsn are likely to impact selection decisions depending on the environment(s) used for selection. If you use this tool, please cite Davis et al. (2025): https://www.biorxiv.org/content/early/2025/01/24/2025.01.21.634104. \n \n Performance and genotype-by-environment interactions of one hybrid pair across three environment pairs. Box plots show population-level yield in the given environment ($n =$ 163 -- 175 plots). Boxes indicate the range from the 25\textsuperscript{th} -- 75\textsuperscript{th} percentile of values. Black lines within the boxplots indicate the median value. Whiskers indicate the most extreme values within 1.5 times the interquartile range and black points indicate the values of data points outside that range. The colored lines connect the mean yields of each hybrid between environments. Asterisks denote a significant difference in mean yields between the two hybrids in the environment at a significance level of alpha = 0.05$ according to Tukey’s HSD. The top panel shows a pair of environments where the hybrid pair did not change rank order and received an interaction importance score of 0. In the middle panel, the hybrids changed rank order in this pair of environments and had a significant difference in mean yield in one of the two environments, and they received an interaction importance score of 1. The bottom panel shows an environment pair where the hybrids changed rank order and had significant differences in mean yield in both environments, and they received an interaction importance score of 2. \textbf{C)} Incidence matrix indicating the frequency with which two hybrids exhibited an interaction for yield between two environments that represents a potentially important change in the selection decision between environments. For each hybrid pair, interaction scores were summed across all environment pairs and divided by the total score possible for the hybrid pair based on the number of environments in which both hybrids were present. Hybrids are ranked in order of ascending yield BLUP values fitting the environment as a fixed effect." })
    output$interactionImportanceScoreExample <- renderImage(
      {
        filename <- normalizePath(file.path('./images/example.png'))
        list(src = filename)
      }, 
      deleteFile = FALSE)
    # Display running message
    jobStatus <- reactiveValues(running=FALSE) 
    output$runningMessage <- renderText({
        if (jobStatus$running) {
          "Calculating... "
        } else {
          ""
        }
    })
    
    observeEvent(input$go, {
        jobStatus$running <- TRUE
        updateActionButton(session, "go", label = "Running")
        
        # Prep data
        data <- input$df$datapath %>%
          read_csv() %>%
          select(all_of(c(input$envCol, input$traitCol, input$genotypeCol))) %>% 
          rename(genotype = .data[[input$genotypeCol]],
                 environment = .data[[input$envCol]]) %>%
          rowwise() %>%
          mutate(environment = str_remove_all(environment, ' ')) %>% 
          mutate(genotype = str_replace_all(genotype, '-', ' ')) %>%
          filter(!is.na(genotype) & !is.na(environment))
        genotypes <- unique(data$genotype)
        num_genotypes <- length(genotypes) - 1
        print(num_genotypes)
        genotypePairs <- tibble()
        for(i in 1:num_genotypes)
        {
          df <- tibble(genotype1 = genotypes[i], genotype2 = genotypes[(i + 1):num_genotypes])
          genotypePairs <- bind_rows(genotypePairs, df)
        }
        
        data_summary <- data %>% 
          group_by(environment, genotype) %>% 
          summarise(traitMean = mean(.data[[input$traitCol]], na.rm = TRUE)) 
        allEnvironments <- unique(data_summary$environment)
        # Run computations
        interactionScores <- getSignificantCrossovers(data, input$traitCol, 'environment', genotypePairs)
        print('scores done')
        importantInteractionScores <- interactionScores %>% 
          select(genotype1, genotype2, contains('Score')) %>% 
          select(!contains('Normalized')) %>% 
          select(genotype1,  genotype2, contains('.E')) %>%
          pivot_longer(contains('Score'), names_to = 'interaction', values_to = 'score', names_prefix = paste0(input$traitCol, 'Score.E')) %>% 
          filter(score >= 1) %>% 
          rowwise() %>%
          mutate(environment1 = str_remove(interaction, paste0(input$traitCol, 'Score.E')) %>%
                   str_split_i('-', 1),
                 environment2 = str_remove(interaction, paste0(input$traitCol, 'Score.E')) %>%
                   str_split_i('-', 2)) %>%
          select(!interaction) %>%
          left_join(data_summary, join_by(genotype1==y$genotype, environment1==y$environment)) %>%
          rename(meanG1E1 = traitMean) %>%
          left_join(data_summary, join_by(genotype2==y$genotype, environment1==y$environment)) %>%
          rename(meanG2E1 = traitMean) %>%
          left_join(data_summary, join_by(genotype1==y$genotype, environment2==y$environment)) %>%
          rename(meanG1E2 = traitMean) %>%
          left_join(data_summary, join_by(genotype2==y$genotype, environment2==y$environment)) %>%
          rename(meanG2E2 = traitMean)
        # render outputs
        grid <- plotInteractionImportanceGrid(significantInteractionsData = interactionScores, performanceData = data, trait = input$traitCol, 
                                              traitLabel = input$traitCol)
        output$interactionImportanceGrid <- renderPlot(grid)
        output$importantInteractions <- renderDT(importantInteractionScores)
        jobStatus$running <- FALSE
        updateActionButton(session, "go", label = "GO")
      }
    )
  }
)
