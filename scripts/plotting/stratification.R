## Author: Michal Kolář

sink("stratification.Rout", split = TRUE)

rm(list = ls())

cmdArgs <- commandArgs()
self_rel_path <- unlist(strsplit(cmdArgs[grep('^--file', cmdArgs)][1], split='='))[2]
self_abs_path <- file.path(getwd(), self_rel_path)
self_dir_path <- dirname(self_abs_path)
root_path <- dirname(dirname(dirname(self_abs_path)))

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	stop("Please, provide a path to the pipeline output folder")
}

paths <- c()
paths["output_dir"] <- file.path(args[1], 'final_analysis')
paths["output_pdf_dir"] <- file.path(paths["output_dir"], "pdf")
paths["output_png_dir"] <- file.path(paths["output_dir"], "png")
paths["output_tables_dir"] <- file.path(paths["output_dir"], "tables")

paths["metadata_tsv"] <- file.path(args[1], 'data', 'metadata.tsv')
paths["pangolin_csv"] <- file.path(args[1], 'pangolin', 'lineage_report.csv')
paths["report_json"] <- file.path(args[1], 'report.json')

dir.create(paths["output_dir"], showWarnings=FALSE)
dir.create(paths["output_pdf_dir"], showWarnings=FALSE)
dir.create(paths["output_png_dir"], showWarnings=FALSE)
dir.create(paths["output_tables_dir"], showWarnings=FALSE)

library(multcomp)
library(jsonlite)
library(lattice)
library(BiocParallel)

##----------------##
##-- parameters --##
##----------------##

bpparam <- MulticoreParam(8)

##-- make an independent sample group from NA samples<
if (impute <-  FALSE) {
    unknown <- "_unknown_"
} else {
    unknown <- NA
}

##-- minimal lineage and mutation frequency in the dataset
minFreqLineage  <- 2
minFreqMutation <- 3

##-- regions west to east
region <- na.omit(c("Karlovarsky", "Plzensky", "Ustecky", "Stredocesky", "Jihocesky", "Praha", "Liberecky", "Vysocina",
                    "Kralovehradecky","Pardubicky", "Jihomoravsky", "Olomoucky", "Moravskoslezsky", "Zlinsky", unknown))
names(region) <- region


##---------------##
##-- functions --##
##---------------##

logisticRegression <- function(d, explanatory, data) {
    ## https://www.datacamp.com/community/tutorials/logistic-regression-R
    ## https://stats.stackexchange.com/questions/60352/comparing-levels-of-factors-after-a-glm-in-r
    ## http://r.qcbs.ca/workshop06/book-en/binomial-glms.html

    glm.model <- as.formula(paste(d, "~ 0 +", paste(explanatory, collapse = "+")))
    glm.fit   <- glm(glm.model, data = data, family = binomial(link = "logit"))

    glm.contrast <- summary(glht(glm.fit, mcp(Week = "Tukey",
                                              Region = "Tukey",
                                              Sex = "Tukey",
                                              Age = "Tukey",
                                              Vaccination = "Tukey",
                                              Select = "Tukey")))

    glm.pvalues <- glm.contrast$test$pvalues
    names(glm.pvalues) <- names(glm.contrast$test$coefficients)

    return(glm.pvalues)
}

plotProfile <- function(p, data) {

    ex <- names(p)
    reC<- grep("^C", p[[1]], value = TRUE)
    names(reC) <- reC
    reL<- grep("^L", p[[1]], value = TRUE)
    names(reL) <- reL

    if (length(reC)) {
        pData <- sapply(reC, function(r) sapply(split(data[, r],  data[, ex]), mean))
        pData <- data.frame(pData, factor(rownames(pData), levels = rownames(pData)))
        colnames(pData)[ncol(pData)] <- ex

        ylim <- c(0, min(1, 1.05 * max(pData[reC])))

        pdf(paste0(paths["output_pdf_dir"], "/", ex, "_nextclade.pdf"), width = 7, height = 4)
        for (r in seq(along = reC))
            print(dotplot(as.formula(paste(reC[r], "~", ex)), data = pData,
                          ylab = gsub("C_", "Nextclade ", reC[r]),
                          pch = 17, col = r, type = "b", ylim = ylim, scales=list(x=list(rot=45))))
        dev.off()
    }

    if (length(reL)) {
        pData <- sapply(reL, function(r) sapply(split(data[, r],  data[, ex]), mean))
        pData <- data.frame(pData, factor(rownames(pData), levels = rownames(pData)))
        colnames(pData)[ncol(pData)] <- ex

        ylim <- c(0, min(1, 1.05 * max(pData[reL])))

        pdf(paste0(paths["output_pdf_dir"], "/", ex, "_pangolin.pdf"), width = 7, height = 4)
        for (r in seq(along = reL))
            print(dotplot(as.formula(paste(reL[r], "~", ex)), data = pData,
                          ylab = gsub("L_", "Pangolin ", reL[r]),
                          pch = 17, col = r, type = "b", ylim = ylim, scales=list(x=list(rot=45))))
        dev.off()
    }

    return(invisible(TRUE))
}

plotProfileMutation <- function(p, data) {

    ex <- names(p)
    re <- p[[1]]
    names(re) <- re

    if (length(re)) {
        pData <- sapply(re, function(r) sapply(split(data[, r],  data[, ex]), mean, na.rm = TRUE))
        pData <- data.frame(pData, factor(rownames(pData), levels = rownames(pData)))
        colnames(pData)[ncol(pData)] <- ex

        ylim <- c(0, min(1, 1.05 * max(pData[re])))

        pdf(paste0(paths["output_pdf_dir"], "/", ex, "_mutation.pdf"), width = 7, height = 4)
        for (r in seq(along = re))
            print(dotplot(as.formula(paste(re[r], "~", ex)), data = pData,
                          ylab = gsub("stop", "*", gsub("del","-", gsub("_", ":", re[r]))),
                          pch = 17, col = r, type = "b", ylim = ylim, scales=list(x=list(rot=45))))
        dev.off()
    }

    return(invisible(TRUE))
}


##------------##
##-- script --##
##------------##

##-- check if there is anything to stratify on
metadata  <- read.csv(paths["metadata_tsv"], sep = "\t")
covCat <- c("gender", "district", "region", "week_number", "Vaccinated", "Vakcina", "N_doses", "Reinfection", "H_18", "H_50", "Imported")
covCon <- c("age", "D_last_dose")

for (c1 in covCat) {
    print(c1)
    print(table(metadata[, c1], useNA = "always"))
}

for (c1 in covCon) {
    print(c1)
    metadata[, c1] <- as.numeric(metadata[, c1])
    print(stem(metadata[, c1]))
}

##-- use only informative covariates
covCat <- c("gender", "region", "week_number", "Vaccinated", "Vakcina", "N_doses", "Reinfection", "Imported")
covCon <- c("age", "D_last_dose")


##-- vaccination status
## valid vaccinations

validVaccine <- rep(NA, times = nrow(metadata))

## Comirnaty, SPIKEVAX, VAXZEVRIA
## later than 14 days after the second dose and before 9 months (270 days) after second dose or indefinite after booster dose
selVaccine <- metadata[, "Vakcina"] %in% c("Comirnaty", "Comirnaty 5-11", "SPIKEVAX", "VAXZEVRIA")

sel2ndDose <- metadata[, "N_doses"] == 2
selTime    <- (metadata[, "D_last_dose"] >= 14) & (metadata[, "D_last_dose"] <= 270)
selValid2ndDose <- selVaccine & sel2ndDose & selTime

sel3rdDose <- metadata[, "N_doses"] == 3
selValid3rdDose <- selVaccine & sel3rdDose

validVaccine[selValid2ndDose] <- TRUE
validVaccine[selValid3rdDose] <- TRUE

## COVID-19 Vaccine Janssen
## later than 14 days after the first dose and before 2 months (60 days) after the first dose
## before 9 months (270 days) after second dose
## or indefinite after booster dose
selVaccine <- metadata[, "Vakcina"] %in% c("COVID-19 Vaccine Janssen")

sel1stDose <- metadata[, "N_doses"] == 1
selTime    <- (metadata[, "D_last_dose"] >= 14) & (metadata[, "D_last_dose"] <= 60)
selValid1stDose <- selVaccine & sel1stDose & selTime

sel2ndDose <- metadata[, "N_doses"] == 2
selTime    <- (metadata[, "D_last_dose"] >= 0) & (metadata[, "D_last_dose"] <= 270)
selValid2ndDose <- selVaccine & sel2ndDose & selTime

sel3rdDose <- metadata[, "N_doses"] == 3
selValid3rdDose <- selVaccine & sel3rdDose

validVaccine[selValid1stDose] <- TRUE
validVaccine[selValid2ndDose] <- TRUE
validVaccine[selValid3rdDose] <- TRUE

## invalid vaccinations (all other vaccinations)
vaccineInfo <- metadata[, "Vaccinated"] %in% c("Yes", "No")

validVaccine[vaccineInfo & is.na(validVaccine)] <- FALSE

table(validVaccine, useNA = "always")

metadata <- cbind(metadata, Valid_Vaccination = validVaccine)


##-- coarse grain age
age <- cut(metadata[, "age"], breaks = c(-1, 15, 30, 40, 50, 65, max(metadata[, "age"], na.rm = TRUE) + 1),
           labels = c("15-", "(15,30)", "(30,40)", "(40,50)", "(50,65)", "65+"))

table(age, useNA = "always")

metadata <- cbind(metadata, Age = age)


##-- sample select variable
select <- rep(NA, times = nrow(metadata))

selReinfectionInfo <- metadata[, "Reinfection"] %in% c("Yes", "No")
selImportedInfo    <- metadata[, "Imported"]    %in% c("Yes", "No")
selH18Info         <- metadata[, "H_18"]        %in% c("Yes", "No")
selH50Info         <- metadata[, "H_50"]        %in% c("Yes", "No")

selAnyAvailable <- (selReinfectionInfo | selImportedInfo | selH18Info | selH50Info)
select[selAnyAvailable] <- "Basal"

selReinfection <- metadata[, "Reinfection"] %in% c("Yes")
select[selReinfection] <- "Reinfection"

selImported <- metadata[, "Imported"] %in% c("Yes")
select[selImported] <- "Imported"

selH18 <- metadata[, "H_18"] %in% c("Yes")
select[selH18] <- "HospUnder18"

selH50 <- metadata[, "H_50"] %in% c("Yes")
select[selH50] <- "HospUnder50"

table(select, useNA = "always")

metadata <- cbind(metadata, Select = select)


##-- explanatory variables
group <- metadata[, c("fasta_id", "week_number", "region", "gender", "Age", "Valid_Vaccination", "Select")]

for (v1 in sort(colnames(group)[-1]))
    for (v2 in sort(colnames(group)[-1]))
        if (v1 < v2) {
            print(v1)
            print(v2)
            print(table(group[, v1], group[, v2], useNA = "always"))
        }

##-- age
if (impute) {
    warning(paste("Replacing",
                  sum(is.na(group[, "Age"]))/length(group[, "Age"]),
                  "missing age values by (40,50). The mean age of other samples is",
                  mean(metadata[, "age"], na.rm = TRUE)))

    group[is.na(group[, "Age"]), "Age"] <- "(40,50)"
}

table(group[, "Age"], useNA = "always")

##-- gender: only male and female allowed
table(group[, "gender"], useNA = "always")

group[, "gender"] <- factor(group[, "gender"], levels = c("Female", "Male"))
table(group[, "gender"], useNA = "always")

##-- region: only actual regions allowed
table(group[, "region"], useNA = "always")

group[group[, "region"] == "", "region"] <- unknown

group[, "region"] <- factor(group[, "region"], levels = region)
table(group[, "region"], useNA = "always")

##-- week_number: only actual weeks kept
table(group[, "week_number"], useNA = "always")

group[, "week_number"] <- factor(group[, "week_number"])

##-- Valid_Vaccination
table(group[, "Valid_Vaccination"], useNA = "always")

vaccinated <- c("Unvaccinated", "Vaccinated")[group[, "Valid_Vaccination"] + 1]
vaccinated[is.na(vaccinated)] <- unknown

group[, "Valid_Vaccination"] <- factor(vaccinated)

table(group[, "Valid_Vaccination"], useNA = "always")

##-- Select
table(group[, "Select"], useNA = "always")

group[is.na(group[, "Select"]), "Select"] <- unknown

group[, "Select"] <- factor(group[, "Select"])

table(group[, "Select"], useNA = "always")

##-- contingency tables
for (v1 in sort(colnames(group)[-1]))
    for (v2 in sort(colnames(group)[-1]))
        if (v1 < v2) {
            print(v1)
            print(v2)
            print(table(group[, v1], group[, v2], useNA = "always"))
        }

colnames(group) <- c("Isolate", "Week", "Region", "Sex", "Age", "Vaccination", "Select")
group <- unique(group)

summary(group)


##-------------------------------------------------##
##-- logistic redression for lineages and clades --##
##-------------------------------------------------##
nextclade <- read.delim(paths["metadata_tsv"], sep = "\t")[, 1:2]
colnames(nextclade) <- c("Isolate", "Clade")

pangolin  <- read.delim(paths["pangolin_csv"], sep = ",")[, 1:2]
colnames(pangolin) <- c("Isolate", "Lineage")

t1 <- match(group[, "Isolate"], nextclade[, "Isolate"])
if (!all(group[, "Isolate"] == nextclade[t1, "Isolate"]))
    stop("Nextclade mismatch, ...")
nextclade <- nextclade[t1, ]
nextclade <- model.matrix(~0+Clade, data = nextclade)
colnames(nextclade) <- substr(gsub("Clade", "C_", colnames(nextclade)), 1, 5)

t1 <- match(group[, "Isolate"], pangolin[, "Isolate"])
if (!all(group[, "Isolate"] == pangolin[t1, "Isolate"]))
    stop("Pangolin mismatch, ...")
pangolin <- pangolin[t1, ]
pangolin <- model.matrix(~0+Lineage, data = pangolin)
colnames(pangolin) <- gsub("Lineage", "L_", colnames(pangolin))

data <- cbind(group, nextclade, pangolin)

explanatory <- colnames(group)[-1]

response <- c(colnames(nextclade), colnames(pangolin))
response <- names(which(colSums(data[, response]) >= minFreqLineage))
names(response) <- response

data <- unique(data)
data <- data[, c(explanatory, response)]
summary(data)

glm.result <- bplapply(response, logisticRegression, explanatory = explanatory, data = data, BPPARAM = bpparam)
glm.result <- do.call(cbind, glm.result)

##-- pseudo Bnferroni (combined with Tukey, ...)
glm.signif <- which(glm.result < 0.05 / length(response), arr.ind = TRUE)
glm.signif <- data.frame(Response    = colnames(glm.result)[glm.signif[, 2]],
                         Contrast    = rownames(glm.result)[glm.signif[, 1]],
                         Explanatory = gsub(":.*", "", rownames(glm.result)[glm.signif[, 1]]),
                         PValue = glm.result[glm.signif])

##-- plot the significant results
glm.plot <- lapply(split(glm.signif[, "Response"], glm.signif[, "Explanatory"]), function(x) sort(unique(x)))
lapply(seq(along = glm.plot), function(x) plotProfile(glm.plot[x], data = data))

write.table(glm.signif, file = file.path(paths["output_tables_dir"], "significantChangesInVariants.tsv"), sep = "\t", row.names = FALSE)


##---------------------------------------##
##-- logistic regression for mutations --##
##---------------------------------------##
mutation <- fromJSON(file(paths["report_json"]))
t1 <- sapply(mutation, function(x) x$cnt) >= minFreqMutation
mutation <- lapply(mutation[t1], function(x) x$seq_ids)
t2 <- unique(unlist(mutation))
t3 <- matrix(0, nrow = length(t2), ncol = length(mutation))
rownames(t3) <- t2
colnames(t3) <- names(mutation)

for (m in seq(along = mutation)) {

    mCol <- match(names(mutation)[m], colnames(t3))
    mRow <- match(mutation[[m]],      rownames(t3))
    t3[mRow, mCol] <- 1
}

colnames(t3) <- gsub("\\*", "stop", gsub("-", "del", gsub(":", "_", colnames(t3))))

t4 <- match(group[, "Isolate"], rownames(t3))

data <- data.frame(group, t3[t4, ])

##-- get rid of "all NA" rows
sel <- apply(data, 1, function(x) !all(is.na(x[colnames(t3)])))
data <- data[sel, ]

explanatory <- colnames(group)[-1]

response <- colnames(t3)
response <- names(which(colSums(data[, response], na.rm = TRUE) >= minFreqMutation))
names(response) <- response

data <- data[, c(explanatory, response)]
summary(data)

glm.result <- bplapply(response, logisticRegression, explanatory = explanatory, data = data, BPPARAM = bpparam)
glm.result <- do.call(cbind, glm.result)

##-- pseudo Bnferroni (combined with Tukey, ...)
glm.signif <- which(glm.result < 0.05 / length(response), arr.ind = TRUE)
glm.signif <- data.frame(Response    = colnames(glm.result)[glm.signif[, 2]],
                         Contrast    = rownames(glm.result)[glm.signif[, 1]],
                         Explanatory = gsub(":.*", "", rownames(glm.result)[glm.signif[, 1]]),
                         PValue = glm.result[glm.signif])

##-- plot the significant results
glm.plot <- lapply(split(glm.signif[, "Response"], glm.signif[, "Explanatory"]), function(x) sort(unique(x)))
lapply(seq(along = glm.plot), function(x) plotProfileMutation(glm.plot[x], data = data))

write.table(glm.signif, file = file.path(paths["output_tables_dir"], "significantChangesInMutations.tsv"), sep = "\t", row.names = FALSE)


##------------------##
##-- session info --##
##------------------##
rm(age, c1, m, mCol, mRow, region, sel, sel1stDose, sel2ndDose, sel3rdDose, selAnyAvailable, selH18, selH18Info, selH50, selH50Info,
   selImported, selImportedInfo, selReinfection, selReinfectionInfo, selTime, selVaccine, selValid1stDose, selValid2ndDose,
   selValid3rdDose, select, t1, t2, t3, t4, v1, v2, vaccinated, vaccineInfo, validVaccine)

save.image("stratification.RData")
sessionInfo()
warnings()
traceback()
sink()
