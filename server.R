library(shiny)
library(DT)
library(readr)
source('R/Functions.R')

shinyServer
(
  function(input, output, session)
  {
    output$introNote <- renderText({ "Determine which empirical GxE interactiopsn are likely to impact selection decisions depending on the environment(s) used for selection. If you use this tool, please cite Davis et al. (2025): https://www.biorxiv.org/content/early/2025/01/24/2025.01.21.634104. Genotypes are arranged in order of ascending BLUP values, i.e., the genotype with the highest BLUP has the highest rank value." })
    # Display running message
    jobStatus <- reactiveValues(running=FALSE) 
    output$runningMessage <- renderText({
      if (jobStatus$running) {
        return("Calculating... ")
      } else {
        return(NULL)
      }
    })
    
    observeEvent(input$go, {
        jobStatus$running = TRUE
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
          summarise(traitMean = mean(.data[[input$traitCol]])) 
        allEnvironments <- unique(data_summary$environment)
        # Run computations
        interactionScores <- getSignificantCrossovers(data, input$traitCol, 'environment', genotypePairs)
        print('scores done')
        importantInteractionScores <- interactionScores %>% 
          select(genotype1, genotype2, contains('Score')) %>% 
          select(!contains('Normalized')) %>% 
          pivot_longer(contains('Score'), names_to = 'interaction', values_to = 'score') %>% 
          filter(score >= 1) #%>% 
          rowwise() %>% 
          mutate(environment1 = str_remove(interaction, paste0(input$traitCol, 'Score.E')) %>% 
                   str_split_i('-.E', 1), 
                 environment2 = str_remove(interaction, paste0(input$traitCol, 'Score.E')) %>% 
                   str_split_i('-.E', 2)) %>% 
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
        output$interactionImportanceScoreExample <- renderImage(
          {
            filename <- normalizePath(file.path('./images/example.png'))
            list(src = filename)
          }, 
          deleteFile = FALSE)
        jobStatus$running = FALSE
        updateActionButton(session, "go", label = "GO")
      }
    )
  }
)
