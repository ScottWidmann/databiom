# My first script

1+1


cells <- c(1, 2, 3, 4)
rnames <- c("R1", "R2")
cnames <- c("C1", "C2")
mymatrix <- matrix(cells, nrow=2, ncol=2, byrow=TRUE, dimnames=list(rnames, cnames))


wave <- c(700, 520, 475, 200)
light <- c("Red", "Green", "Blue", "Ultraviolet")
vis <- c(TRUE, TRUE, TRUE, FALSE)
mydf <- data.frame(wave, light, vis)
names(mydf) <- c("Wavelength", "Color", "Visible")


write.table(mymatrix, "mymatrix.txt", sep = "\t", quote=FALSE)
write.csv(mydf, "mydf.csv", row.names = FALSE, quote = FALSE)



mycsv <- read.csv("processed.cleveland.data.csv")
head(mycsv)

library(tidyverse)
mycsv <- read_csv("processed.cleveland.data.csv")
head(mycsv)

#Summary stats
summary(mycsv)

#Using the pipe operator ( %>% )
mycsv %>% summary()

#Drop missing values and get summary stats
mycsv %>% drop_na() %>% summary()

#Drop NAs and create a new variable
mycsv_noNAs <- mycsv %>% drop_na()
mycsv_noNAs
mycsv 

#Install DESeq2 Package via BiocManager
BiocManager::install("DESeq2")

library(DESeq2)
