# shinyGxEScores
RShiny app to score importance of GxE to selection decisions

## To Run Locally: 
Open R and run the following code. Example data is available in data/.
```
if (!require('shiny'))
{
  install.packages('shiny')
}
library('shiny')
runGitHub('shinyGxEScores', 'jdavis-132', launch.browser=TRUE)
```
