
library(knitr)

knitr::purl(input = "01-data.Rmd",
            output = "./R/01-data.R",
            documentation = 0)

knitr::purl(input = "03-cointegrated-var.Rmd",
            output = "./R/03-cointegrated-var.R",
            documentation = 0)
