#!/usr/bin/env Rscript

library(ggplot2)

args<-commandArgs(TRUE)

fn=args[1]

cov <- read.delim(fn, header=FALSE)
names(cov) <- c("CHR", "POSITION", "COVERAGE")

png(filename="coverage.pdf")
ggplot(cov, aes(POSITION, COVERAGE)) + geom_point() + geom_hline(yintercept = 30, color = "green") +
ggtitle(paste("Depth of Coverage"))
dev.off()
