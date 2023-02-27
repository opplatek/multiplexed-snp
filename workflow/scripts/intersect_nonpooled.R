# script to calculate overlaps for SNP analysis
# this script works with extract_snps.sh bash script

# library(xlsx)
library(gtools) # Using mixedsort http://www.inside-r.org/packages/cran/gtools/docs/mixedsort
###
# setwd("/home/jan/Projects/vets/2015/results/kaja/bwamem_PVRL4e2_10062015")
###
#files <- list.files(pattern = "*.extract.snp") # store names of all input files
files <- snakemake@input[[1]]
ofile <- snakemake@output[[1]]
###
# R doesn't have function fo cbind data.frames with different length so we need to create a function
cbind.fill <- function(...) { # create a function to cbind different length data.frames
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function(x) {
    rbind(x, matrix(, n - nrow(x), ncol(x)))
  }))
}
###
# Remove empty files from the main list of input files. Most likely no variants identified
emptyFiles <- NULL
for (i in 1:length(files)) # loop for calculate the overlap for all input files
{
  # print(files[i])
  if (file.info(files[i])[[1]] == 0) {
    emptyFiles <- rbind(emptyFiles, files[i])
  }
}

files <- files[!(files %in% emptyFiles)]

for (i in 1:length(files)) # loop for calculate the overlap for all input files
{
  print(files[i])
  # Read table and if empty create dummy matrix
  if (file.info(files[i])[[1]] == 0) {
    overlap <- as.data.frame(matrix(ncol = 4, nrow = 1))
  } else {
    overlap <- read.table(file = files[i], sep = "\t", fill = T, header = T, stringsAsFactors = F)
  }

  ############################ New overlap##################################
  # calculate overlaps between all columns
  tmp_overlap <- rbind(as.matrix(overlap[, 1]), as.matrix(overlap[, 2]), as.matrix(overlap[, 3])) # Merge all columns to one big column
  tmp_overlap <- as.matrix(tmp_overlap[tmp_overlap != "-", ]) # Remove all "-" if there are any
  tmp_overlap <- as.matrix(tmp_overlap[tmp_overlap != "", ]) # Remove all empty if there are any
  tmp_overlap <- as.matrix(tmp_overlap[!is.na(tmp_overlap)]) # remove NA values
  tmp_overlap <- trimws(tmp_overlap, which = c("both")) # Remove whitespaces at the beginning and end if present http://stackoverflow.com/questions/2261079/how-to-trim-leading-and-trailing-whitespace-in-r
  # tmp_overlap<-as.matrix(sort(tmp_overlap[,1], na.last = NA)) # Sort the big column by numbers and exclude NA values
  tmp_overlap <- as.matrix(mixedsort(tmp_overlap, na.last = NA)) # Sort the big column by numbers and exclude NA values
  # To use sort correctly the table must be factor-less, therefore matrix or similar
  tmp_overlap2 <- as.data.frame(unique(tmp_overlap)) # Extract only unique positions
  # tmp_overlap2<-as.data.frame(mixedsort(tmp_overlap2[,1], na.last = NA)) # Sort the big column by numbers and exclude NA values
  tmp_overlap <- as.data.frame(tmp_overlap) # Transform list do data.frame
  # tmp_overlap2<-as.data.frame(tmp_overlap2) # Transform list do data.frame

  for (tmp_pos in 1:nrow(tmp_overlap2)) { # For all unique positions
    #  tmp_overlap2[tmp_pos,"sum"]<-sum(tmp_overlap==tmp_overlap2[tmp_pos,1]) # Sum their occurences and put them to new columns "sum
    tmp_overlap2[tmp_pos, "sum"] <- sum(unlist(tmp_overlap[, 1]) %in% unlist(tmp_overlap2[tmp_pos, 1])) # Sum their occurences and put them to new columns "sum"
  }

  # All occurences
  all <- tmp_overlap2[tmp_overlap2[, 2] == 3, 1]
  all <- as.matrix(all)
  all <- as.data.frame(mixedsort(all[, 1])) # sort the unique values (SNPs)

  # Three occurencces
  # justThree<-tmp_overlap2[tmp_overlap2[,2]>=3,1]
  # calculate overlaps between at least two columns
  justTwo <- tmp_overlap2[tmp_overlap2[, 2] >= 2, 1]
  justTwo <- as.matrix(justTwo)
  justTwo <- as.data.frame(mixedsort(justTwo[, 1])) # sort the unique values (SNPs)

  ############################ New overlap##################################

  ############################ Old overlap##################################
  # #calculate overlaps between all columns
  # i12<-t(intersect(overlap[,1],overlap[,2])) #calculate overlaps between columns
  # i13<-t(intersect(overlap[,1],overlap[,3]))
  # i23<-t(intersect(overlap[,2],overlap[,3]))
  # all<-as.data.frame(intersect(i12,i23)) #calculate overlap of all columns
  # all <- as.data.frame(all[!is.na(all)]) #remove NA values
  #
  # #calculate overlaps between at least two columns
  # justTwo<-as.matrix(cbind(i12, i13)) #merge two overlaps
  # justTwo<-as.matrix(t(cbind(justTwo, i23))) #merge merged overlaps with last overlap
  # justTwo <- justTwo[!is.na(justTwo)] #remove NA values
  # justTwo<-as.data.frame(unique(justTwo)) #extract only unique values (SNPs)
  # justTwo<-as.data.frame(sort(justTwo[,1])) #sort the unique values (SNPs)
  ############################ Old overlap##################################

  # merge all results to one final table
  final_table <- as.data.frame(cbind.fill(overlap, all, justTwo)) # merge overlaps with our function
  colnames(final_table)[length(final_table) - 1] <- "overlap_all" # rename columns
  colnames(final_table)[length(final_table)] <- "overlap_two" # rename columns

  # There might be problems if there is only one SNP in the results - create table with one row which
  # is later during conversion as.data.frame turned into one column data.frame instead of one row
  # data.frame. To avoid this we add one additional dummy column which we will delete later
  if (nrow(final_table) == 1) { # Check if there is only one row in the table and if True add one row
    final_table[2, ] <- NA # Add dummy row
    DELETE <- T # Save information that we should delete it afterwards
  } else {
    DELETE <- F # If there are more than one row keep information DELETE=F
  }

  # final adjustments to the final table
  final_table <- sapply(final_table, as.character) # store data.frame as character to replace NA values
  final_table[is.na(final_table)] <- "" # replace NA values with nothing
  final_table[final_table == ""] <- "- -" # replace empty cells with "- -" - will be suitable for later formating
  final_table <- as.data.frame(final_table)
  final_table <- sapply(final_table, as.factor) # store data.frame back to factor

  if (DELETE == T) { # If we had only one SNP - one row in original table we added dummy row and we have to
    # delete it
    final_table <- as.data.frame(final_table[-2, ]) # remove added dummy row
    final_table <- t(final_table) # and transpose the table to keep it in columns instead of rows
    rownames(final_table)[1] <- 1 # Rename rowname because otherwise it would keep final_table[-2, ] in the
    # name and it would case troubles in calculation
  }

  # if final table is empty create dummy matrix
  if (ncol(final_table) < 4) {
    final_table <- as.data.frame(matrix(ncol = 4, nrow = 1))
  }

  tmp <- final_table[, 4]
  tmp[tmp == "- -"] <- "-"
  final_table[, 4] <- tmp
#  write.table(final_table, file = files[i], sep = "\t", col.names = NA) # write the final table
  write.table(final_table, file = ofile, sep = "\t", col.names = NA, quote = F) # write the final table
}

# Write table for empty input files
for (emptyFile in emptyFiles)
{
  sink(file = emptyFile)
  sink(file = ofile)
    print("No results, empty input files. Probably no variants identified")
  sink()
}
