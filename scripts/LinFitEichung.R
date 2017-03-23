library("Hmisc")
source("functions.R")

x <- c(13.37, 17.44, 22.10, 32.06, 44.23)
y1 <- c(178.96, 237.79, 301.68, 418.56, 588.46)

erry1 <- c(11, 13, 16, 19, 24)


drawCIlim(x, y1, erry1,  xlab = "Energy/keV", ylab = "channel")

grid()

text(x=40, y=100, "y=(13.1+-0.3)*channel/keV*x+(7+-6)*channel", cex=0.8)
fm1 <- lm(y1 ~ x, weights=1/erry1^2)
summary(fm1)
abline(fm1, col = "red")


#*******************************************************************************
