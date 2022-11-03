library(dplyr)
library(ggplot2)
library(ggalluvial)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	stop("Please, provide a path to the pipeline output folder")
}


input_path <- paste(args[1], 'final_analysis', 'cladeAssignment.tsv', sep='/')
output_path <- paste(args[1], 'final_analysis', 'cladeAssignments.pdf', sep='/')


t1 <- read.table(input_path, sep = "\t", header = TRUE)
t2 <- names(sort(table(t1[, "UShER"]), decreasing = TRUE))
t1[, "UShER"] <- factor(t1[, "UShER"], levels = t2, labels = t2)
t2 <- names(sort(table(t1[, "Nextclade"]), decreasing = TRUE))
t1[, "Nextclade"] <- factor(t1[, "Nextclade"], levels = t2, labels = t2)

cladeAssignment <- t1[order(as.numeric(t1[, "UShER"]), as.numeric(t1[, "Nextclade"])), ]
rm(t1, t2)

cladeFlow <- as.data.frame(cladeAssignment %>% group_by(UShER, Nextclade) %>% summarise(n = n()))

pdf(output_path, height = 11, width = 7)
ggplot(as.data.frame(cladeFlow),
	aes(y = n, axis1 = UShER, axis2 = Nextclade, fill = UShER)) +
	    geom_alluvium() +
	    geom_stratum() +
	    geom_text(stat = "stratum", aes(label = after_stat(stratum)), min.y = 50) +
	    scale_x_discrete(limits = c("UShER", "Nextclade"))
dev.off()
