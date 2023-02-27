#
# Get mapping stats from parsed flagstat results and plot them
#
# TODO: Make this a side-by-side barplot split by 10 reference and put maybe (?)4 plots on one page
#

library(stringr)
library(gtools)

header <- c("reference", "index", "mapped")

## statsBefore<-read.table("beforeFilter.flagstat.txt",sep="\t", stringsAsFactors = F)
#statsAfter <- read.table("afterFilter.flagstat.txt", sep = "\t", stringsAsFactors = F)
statsAfter <- read.table(snakemake@input[[1]], sep = "\t", stringsAsFactors = F)

colnames(statsAfter) <- header

# statsBefore$mapped<-as.numeric(str_split_fixed(as.character(statsBefore$mapped), " ", n=2)[,1]) # Keep only first column from mapped
statsAfter$mapped <- as.numeric(str_split_fixed(as.character(statsAfter$mapped), " ", n = 2)[, 1])

statsAfter <- statsAfter[mixedorder(statsAfter$index, decreasing = TRUE), ]

pdf(snakemake@output[[1]], height = 10) # Rename by removing everything after last "_"
  for (reference in unique(statsAfter$reference)) {
      barplot(statsAfter[statsAfter$reference == reference, 3],
        names.arg = statsAfter[statsAfter$reference == reference, 2],
        cex.names = 0.7, las = 2, xlim = c(0, max(statsAfter[statsAfter$reference == reference, 3]) * 1.15),
        horiz = T, col = rainbow(n = ncol(statsAfter[statsAfter$reference == reference, ])), space = 0.5, beside = F,
        xlab = "Number of mappings after filtering", main = paste0("Mapping to ", gsub("(.*)_.*", "\\1", reference))
      )
  }
dev.off()