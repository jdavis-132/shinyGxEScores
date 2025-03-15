# shinyGxEScores
RShiny app to score importance of GxE to selection decisions
If you use this tool in your work, please cite Davis et al. (2025): https://www.biorxiv.org/content/early/2025/01/24/2025.01.21.634104. 

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
## To see it hosted on AWS: 
Go to http://13.58.27.199:3838/shinyGxEScores/. 
