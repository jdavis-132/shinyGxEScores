library(shiny)

shinyUI
(
  fluidPage
  (
    titlePanel('Welcome to ShinyGxEScores!'),
    textOutput("introNote"),
    sidebarPanel
    (
      fileInput("df", "Choose csv file with multi-environment data", accept = ".csv"),
      textInput("envCol", "Name of environment identifier column", "environment", placeholder = "environment"),
      textInput("traitCol", "Name of trait column to analyze GxE interactions for", "yieldPerAcre", placeholder = "yieldPerAcre"),
      textInput("genotypeCol", "Name of genotype identifier column", "genotype", placeholder = "genotype"),
      actionButton("go", "GO"), 
      width = 2
    ),
    
    mainPanel
    (
      # Display running message
      uiOutput("runningMessage"),
      imageOutput("interactionImportanceScoreExample"),
      plotOutput("interactionImportanceGrid"), 
      dataTableOutput("importantInteractions")
    )
  )
)
