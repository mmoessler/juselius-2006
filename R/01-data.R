
# Load data
library(readxl)
data <- read_excel("./CointegratedVARModelHandbook/book.xls")
data <- data[,colnames(data)[-1]]



# Number of observations
TT <- nrow(data)

# Define time variable
data$time <- seq(as.Date("1973-01-01"), as.Date("2003-01-01"), by="quarter")



# Create "transitory intervention" dummy at 1975:4
data$Dt754 <- rep(0,TT)
data$Dt754[which(data$time == as.Date("1975-10-01"))] <- 1
data$Dt754[which(data$time == as.Date("1976-01-01"))] <- -0.5
data$Dt754[which(data$time == as.Date("1976-04-01"))] <- -0.5



# Create "permanent intervention" dummy at 1976:4
data$Dp764 <- rep(0,TT)
data$Dp764[which(data$time == as.Date("1976-10-01"))] <- 1



# Create "mean shift" dummy at 1983:1
data$Ds831 <- rep(0,TT)
data$Ds831[which(data$time %in% seq(as.Date("1983-01-01"), as.Date("2003-01-01"), by="quarter"))] <- 1



# Create difference of "mean shift" dummy, i.e., "permanent intervention" dummy, at 1983:1
data$Dp831 <- rep(0,TT)
data$Dp831[which(data$time == as.Date("1983-01-01"))] <- 1



# Define lag of "permanent intervention" dummy, at 1983:1
data$Dp832 <- rep(0,TT)
data$Dp832[which(data$time == as.Date("1983-04-01"))] <- 1

