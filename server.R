library(shiny)
library(DT)
source('R/Functions.R')

shinyServer
(
  function(input, output, session)
  {
    output$introNote <- renderText({ "Determine which empirical GxE interactiopsn are likely to impact selection decisions depending on the environment(s) used for selection. If you use this tool, please cite Davis et al. (2025): https://www.biorxiv.org/content/early/2025/01/24/2025.01.21.634104." })
    
    # Display running message
    jobStatus <- reactiveValues(running=FALSE) 
    output$runningMessage <- renderUI({
      if (jobStatus$running) {
        return(tags$div(
          class = "alert alert-info",
          "Calculating... "
        ))
      } else {
        return(NULL)
      }
    })
    
    observeEvent
    (
      input$go,
      {
        jobStatus$running = TRUE
        updateActionButton(session, "go", label = "Running")
        
        # Prep data
        data <- read.csv(input$df) %>%
          select(all_of(c(input$envCol, input$traitCol, input$genotypeCol))) %>% 
          filter(!is.na(.data[[input$envCol]] & !is.n(.data[[input$genotypeCol]])))
        genotypes <- unique(data[[input$genotypeCol]])
        num_genotypes <- length(genotypes)
        genotypePairs <- tibble()
        for(i in 1:(num_genotypes - 1))
        {
          df <- tibble(genotype1 = allGenotypes[i], genotype2 = allGenotypes[(i + 1):num_genotypes])
          genotypePairs <- bind_rows(genotypePairs, df)
        }
        environments <- unique(data[[input$envCol]])
        
        data_summary <- data %>% 
          group_by(all_of(c(input$envCol, input$genotypeCol))) %>% 
          summarise(traitMean = mean(.data[[input$traitCol]]))
        
        # Run computations
        interactionScores <- getSignificantCrossovers(data, input$traitCol, environments)
        importantInteractionScores <- interactionScores %>% 
          select(genotype1, genotype2, contains('Score')) %>% 
          select(!contains('Normalized')) %>% 
          pivot_longer(contains('Score'), names_to = 'interaction', values_to = 'score') %>% 
          filter(score > 0) %>% 
          rowwise() %>% 
          mutate(environment1 = str_remove(interaction, paste0(input$traitCol, 'Score.E')) %>% 
                   str_split_i('-.E', 1), 
                 environment2 = str_remove(interaction, paste0(input$traitCol, 'Score.E')) %>% 
                   str_split_i('-.E', 2)) %>% 
          select(!interaction) %>% 
          left_join(data_summary, join_by(genotype1==y[[input$genotypeCol]], environment1==y[[input$envCol]])) %>% 
          rename(meanG1E1 = traitMean) %>% 
          left_join(data_summary, join_by(genotype2==y[[input$genotypeCol]], environment1==y[[input$envCol]])) %>% 
          rename(meanG2E1 = traitMean) %>% 
          left_join(data_summary, join_by(genotype1==y[[input$genotypeCol]], environment2==y[[input$envCol]])) %>% 
          rename(meanG1E2 = traitMean) %>% 
          left_join(data_summary, join_by(genotype2==y[[input$genotypeCol]], environment2==y[[input$envCol]])) %>% 
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
