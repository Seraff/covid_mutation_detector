library(dplyr)
library(ggplot2)
library(ggalluvial)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	stop("Please, provide a path to the pipeline output folder")
}


input_path <- paste(args[1], 'final_analysis', 'lineageAssignmentChanges.tsv', sep='/')
output_path <- paste(args[1], 'final_analysis', 'lineageAssignments.pdf', sep='/')


t1 <- read.table(input_path, sep = "\t", header = TRUE)
t2 <- names(sort(table(t1[, "Nextclade_pango"]), decreasing = TRUE))
t1[, "Nextclade_pango"] <- factor(t1[, "Nextclade_pango"], levels = t2, labels = t2)
t2 <- names(sort(table(t1[, "UShER"]), decreasing = TRUE))
t1[, "UShER"] <- factor(t1[, "UShER"], levels = t2, labels = t2)

lineageAssignment <- t1[order(as.numeric(t1[, "Nextclade_pango"]), as.numeric(t1[, "UShER"])), ]
rm(t1, t2)

lineageFlow <- as.data.frame(lineageAssignment %>% group_by(Nextclade_pango, UShER) %>% summarise(n = n()))

pdf(output_path, height = 11, width = 7)
ggplot(as.data.frame(lineageFlow),
	aes(y = n, axis1 = Nextclade_pango, axis2 = UShER, fill = Nextclade_pango)) +
	    geom_alluvium() +
	    geom_stratum() +
	    geom_text(stat = "stratum", aes(label = after_stat(stratum)), min.y = 50) +
	    scale_x_discrete(limits = c("Nextclade_pango", "UShER"))
dev.off()
