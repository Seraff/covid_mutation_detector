sink("mutationTable.Rout", split = TRUE)

dir.create("./Pdf")
dir.create("./Png")
dir.create("./Tables")

library(gdata)
library(seqinr)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(jsonlite)
library(stringr)

##-- parameters
genomeLength <- 29903
region <- c("Ceska republika", "Jihocesky", "Jihomoravsky", "Karlovarsky", "Kralovehradecky", "Liberecky",
            "Moravskoslezsky", "Olomoucky", "Pardubicky", "Plzensky", "Praha", "Stredocesky", "Ustecky",
            "Vysocina", "Zlinsky", "_unknown_")
names(region) <- region

minFreqLineage  <- 2
minProbLineage  <- 0.01
minFreqMutation <- 3
minProbMutation <- 0.10

heatmapSplit <- c(0, 1/9, 1/3, 1)

##-- functions
replaceDeletion <- function(x, from, to) {

    if (from %in% x[, "mutation"]) {
        x1 <- matrix(rep(x[x[, "mutation"] == from, ], length(to)), byrow = TRUE, ncol = ncol(x))
        x1[, 1] <- to
        colnames(x1) <- colnames(x)
        x <- x[x[, "mutation"] != from, ]
        x <- rbind(x, x1)
    }

    return(x)
}

createTableSimple <- function(r, f) {

    if (r == "Ceska republika") {
        id <- metadata[, "fasta_id"]
    } else if (r == "_unknown_") {
        id <- metadata[metadata[, "region"] == "", "fasta_id"]
    } else {
        id <- metadata[metadata[, "region"] == r, "fasta_id"]
    }

    if (length(id) < 1)
        return(data.frame(Factor = character(0), Freq = numeric(0), Prob = numeric(0), Count = numeric(0)))

    dat <- metadata[id, f]
    tab <- data.frame(table(unlist(dat), useNA = "ifany"))
    colnames(tab)[1] <- "Factor"

    tab <- cbind(tab, Prob = tab[, "Freq"] / length(id), Count = length(id))
    tab <- tab[order(tab[, "Freq"], decreasing = TRUE), ]
    return(tab)
}

createTableList <- function(r, f) {

    if (r == "Ceska republika") {
        id <- metadata[, "fasta_id"]
    } else if (r == "_unknown_") {
        id <- metadata[metadata[, "region"] == "", "fasta_id"]
    } else {
        id <- metadata[metadata[, "region"] == r, "fasta_id"]
    }

    if (length(id) < 1)
        return(data.frame(Factor = character(0), Freq = numeric(0), Prob = numeric(0), Count = numeric(0)))

    dat <- strsplit(metadata[id, f], ",")
    tab <- data.frame(table(unlist(dat), useNA = "ifany"))
    colnames(tab)[1] <- "Factor"

    tab <- cbind(tab, Prob = tab[, "Freq"] / length(id), Count = length(id))
    tab <- tab[order(tab[, "Freq"], decreasing = TRUE), ]
    return(tab)
}

selectLevelsToShow <- function(factorTable, minFreq, minProb) {

    factorSel <- sort(unique(unlist(lapply(factorTable, function(x) {
        as.character(x[( (x[, "Freq"] >= minFreq) & (x[, "Prob"] >= minProb) ), 1])
    }))))
    t1 <- factorTable[["Ceska republika"]]
    t1 <- t1[t1[, "Factor"] %in% factorSel, ]
    factorSel <- as.character(t1[order(t1[, "Freq"], decreasing = TRUE), "Factor"])

    return(factorSel)
}

summarizeLevelsToShow <- function(r, factorTable, factorSel) {

    tab <- factorTable[[r]]

    s1 <- sum(tab[, "Freq"])
    s2 <- sum(tab[, "Prob"])

    tab <- tab[tab[, "Factor"] %in% factorSel, ]
    tab[, "Factor"] <- factor(tab[, "Factor"], levels = c(factorSel, "Other"))
    tab <- rbind(tab,
                 data.frame(Factor = "Other", Freq = s1 - sum(tab[, "Freq"]), Prob = s2 - sum(tab[, "Prob"]),
                            Count = tab[1, "Count"]))

    tab <- cbind(tab, Region = r)
    return(tab)
}

createTableSimpleWeek <- function(w, f) {

    if (w == "Vsechny tydny") {
        id <- metadata[, "fasta_id"]
    } else if (w == "_unknown_") {
        id <- metadata[metadata[, "week"] == "tyden_NA", "fasta_id"]
    } else {
        id <- metadata[metadata[, "week"] == w, "fasta_id"]
    }

    dat <- metadata[id, f]
    tab <- data.frame(table(unlist(dat), useNA = "ifany"))
    colnames(tab)[1] <- "Factor"

    tab <- cbind(tab, Prob = tab[, "Freq"] / length(id), Count = length(id))
    tab <- tab[order(tab[, "Freq"], decreasing = TRUE), ]
    return(tab)
}

createTableListWeek <- function(w, f) {

    if (w == "Vsechny tydny") {
        id <- metadata[, "fasta_id"]
    } else if (w == "_unknown_") {
        id <- metadata[metadata[, "week"] == "tyden_NA", "fasta_id"]
    } else {
        id <- metadata[metadata[, "week"] == w, "fasta_id"]
    }

    dat <- strsplit(metadata[id, f], ",")
    tab <- data.frame(table(unlist(dat), useNA = "ifany"))
    colnames(tab)[1] <- "Factor"

    tab <- cbind(tab, Prob = tab[, "Freq"] / length(id), Count = length(id))
    tab <- tab[order(tab[, "Freq"], decreasing = TRUE), ]
    return(tab)
}

selectWeekLevelsToShow <- function(factorTable, minFreq, minProb) {

    factorSel <- sort(unique(unlist(lapply(factorTable, function(x) {
        as.character(x[( (x[, "Freq"] >= minFreq) & (x[, "Prob"] >= minProb) ), 1])
    }))))
    t1 <- factorTable[["Vsechny tydny"]]
    t1 <- t1[t1[, "Factor"] %in% factorSel, ]
    factorSel <- as.character(t1[order(t1[, "Freq"], decreasing = TRUE), "Factor"])

    return(factorSel)
}

summarizeWeekLevelsToShow <- function(w, factorTable, factorSel) {

    tab <- factorTable[[w]]

    s1 <- sum(tab[, "Freq"])
    s2 <- sum(tab[, "Prob"])

    tab <- tab[tab[, "Factor"] %in% factorSel, ]
    tab[, "Factor"] <- factor(tab[, "Factor"], levels = c(factorSel, "Other"))
    tab <- rbind(tab,
                 data.frame(Factor = "Other", Freq = s1 - sum(tab[, "Freq"]), Prob = s2 - sum(tab[, "Prob"]),
                            Count = tab[1, "Count"]))

    tab <- cbind(tab, Week = w)
    return(tab)
}


##-- variant definition from outbreak.info

## list of detected variants, VOC, VOI and VUM (tazbyk-joCpy5-jycqet)
detected <- system("cut -f2 -d',' datafreeze-*_12_weeks_good.pangolin-*-UShER-*.csv | sort -u | grep -v lineage | grep -v None", intern = TRUE)

detected <- setdiff(detected,
                    c("BA.4.1.5", "BA.4.1.8",
                      "BA.2.75.2",
                      "BA.5.1.10", "BA.5.1.5", "BA.5.1.9", "BA.5.10", "BA.5.8", "BA.5.9",
                      "BF.10", "BF.11", "BF.12", "BF.14", "BF.15","BF.6","BF.7", "BF.8",
                      "BK.1",
                      "XAJ", "XAN",
                      "Unassigned"))

## https://outbreak.info/situation-reports
VOC <- c("Delta", "Omicron", "Alpha", "Beta", "Gamma", "BA.1", "BA.2", "BA.4", "BA.5")
VOI <- c("Lambda", "Mu", "BA.3")##, "BA.2.75")
VUM <- c("C.1.2", "Eta", "Iota", "Kappa", "B.1.617.3", "Theta", "B.1.1.318-related", "C.36.3-related")

MUC   <- c()
MUI   <- c()

variant <- sort(unique(c(detected, VOC, VOI, VUM)))

if (! (dir.exists("VariantSignatures/") & all(variant %in% gsub(".tsv$", "", dir("VariantSignatures/")))) ) {
    missingVariant <- setdiff(variant, gsub(".tsv$", "", dir("VariantSignatures/")))
    ## use web scraper to download characteristic mutations for each missing variant
    system(paste("export PYTHONPATH=/Users/kolarmi/opt/miniconda3/lib/python3.9/site-packages; python3 scrape.py", paste(missingVariant, collapse = " ")))
}

t2 <- lapply(variant, function(v) {
    v1 <- read.delim(paste0("./VariantSignatures/", v, ".tsv"), sep = "\t", header = FALSE)

    if (as.numeric(gsub(",", "", strsplit(v1[1, 3], "/")[[1]][2])) < 500)
        warning(paste("Variant", v, "has less than 500 isolates.\n"))

    v2 <- data.frame(mutation = v1[, 1], lineage = v, prevalence = as.numeric(v1[, 2]))
})

t2 <- do.call(rbind, t2)
t2[, "mutation"] <- toupper(t2[, "mutation"])

t3 <- reshape(t2, idvar = "mutation", timevar = "lineage", direction = "wide")
t3 <- as.matrix(t3)

## match mutation format
t3 <- replaceDeletion(x = t3, from = "ORF1A:DEL1-1",       to = "ORF1A:M1-")
t3 <- replaceDeletion(x = t3, from = "ORF1A:DEL17-17",     to = "ORF1A:S17-")
t3 <- replaceDeletion(x = t3, from = "ORF1A:DEL82-86",   to = c("ORF1a:G82-", "ORF1A:H83-", "ORF1A:V84-", "ORF1A:M85-", "ORF1A:V86-"))
t3 <- replaceDeletion(x = t3, from = "ORF1A:DEL141-143",   to = c("ORF1a:K141-", "ORF1A:S142-", "ORF1A:F143-"))
t3 <- replaceDeletion(x = t3, from = "ORF1A:DEL1698-1698", to = "ORF1A:A1698-")
t3 <- replaceDeletion(x = t3, from = "ORF1A:DEL2084-2084", to = "ORF1a:L2084-")
t3 <- replaceDeletion(x = t3, from = "ORF1A:DEL3674-3676", to = c("ORF1a:L3674-", "ORF1A:S3675-", "ORF1A:G3676-"))
t3 <- replaceDeletion(x = t3, from = "ORF1A:DEL3675-3677", to = c("ORF1A:S3675-", "ORF1A:G3676-", "ORF1A:F3677-"))
t3 <- replaceDeletion(x = t3, from = "ORF1A:DEL3676-3678", to = c("ORF1A:G3676-", "ORF1A:F3677-", "ORF1A:K3678-"))

t3 <- replaceDeletion(x = t3, from = "ORF1B:DEL1129-1129", to = "ORF1B:G1129-")
t3 <- replaceDeletion(x = t3, from = "ORF1B:DEL1353-1353", to = "ORF1B:K1353-")
t3 <- replaceDeletion(x = t3, from = "ORF1B:DEL1545-1545", to = "ORF1B:T1545-")

t3 <- replaceDeletion(x = t3, from = "ORF3A:DEL96-96",     to = "ORF3A:L96-")
t3 <- replaceDeletion(x = t3, from = "ORF3A:DEL255-256",   to = c("ORF3A:V255-", "ORF3A:V256-"))
t3 <- replaceDeletion(x = t3, from = "ORF3A:DEL256-256",   to = "ORF3a:V256-")
t3 <- replaceDeletion(x = t3, from = "ORF3A:DEL256-257",   to = c("ORF3A:V256-", "ORF3a:N257-"))
t3 <- replaceDeletion(x = t3, from = "ORF3A:DEL257-257",   to = "ORF3a:N257-")

t3 <- replaceDeletion(x = t3, from = "ORF6:DEL2-3",        to = c("ORF6:F2-", "ORF6:H3-"))
t3 <- replaceDeletion(x = t3, from = "ORF6:DEL3-3",        to = "ORF6:H3-")

t3 <- replaceDeletion(x = t3, from = "ORF7A:DEL72-77",     to = c("ORF7A:K72-", "ORF7A:H73-", "ORF7A:V74-", "ORF7A:Y75-", "ORF7A:H76-", "ORF7A:L77-"))
t3 <- replaceDeletion(x = t3, from = "ORF7A:DEL101-104",   to = c("ORF7A:F101-", "ORF7A:L102-", "ORF7A:I103-", "ORF7A:V104-"))

t3 <- replaceDeletion(x = t3, from = "ORF8:DEL67-68",      to = c("ORF8:S67-", "ORF8:K68-"))
t3 <- replaceDeletion(x = t3, from = "ORF8:DEL87-87",      to = c("ORF8:T87-"))
t3 <- replaceDeletion(x = t3, from = "ORF8:DEL119-120",    to = c("ORF8:D119-", "ORF8:F120-"))
t3 <- replaceDeletion(x = t3, from = "ORF8:DEL120-121",    to = c("ORF8:F120-", "ORF8:I121-"))
t3 <- replaceDeletion(x = t3, from = "ORF8:DEL121-121",    to = "ORF8:I121-")

t3 <- replaceDeletion(x = t3, from = "N:DEL3-3",           to = "N:D3-")
t3 <- replaceDeletion(x = t3, from = "N:DEL31-33",         to = c("N:E31-", "N:R32-", "N:S33-"))
t3 <- replaceDeletion(x = t3, from = "N:DEL204-204",       to = "N:G204-")
t3 <- replaceDeletion(x = t3, from = "N:DEL209-209",       to = "N:R209-")
t3 <- replaceDeletion(x = t3, from = "N:DEL215-215",       to = "N:G215-")

t3 <- replaceDeletion(x = t3, from = "S:DEL25-27",         to = c("S:P25-","S:P26-","S:A27-"))
t3 <- replaceDeletion(x = t3, from = "S:DEL69-70",         to = c("S:H69-", "S:V70-"))
t3 <- replaceDeletion(x = t3, from = "S:DEL136-144",       to = c("S:C136-", "S:N137-", "S:D138-", "S:P139-", "S:F140-", "S:L141-", "S:G142-", "S:V143-", "S:Y144-"))
t3 <- replaceDeletion(x = t3, from = "S:DEL140-143",       to = c("S:F140-", "S:L141-", "S:G142-", "S:V143-"))
t3 <- replaceDeletion(x = t3, from = "S:DEL141-143",       to = c("S:L141-", "S:G142-", "S:V143-"))
t3 <- replaceDeletion(x = t3, from = "S:DEL143-145",       to = c("S:V143-", "S:Y144-", "S:Y145-"))
t3 <- replaceDeletion(x = t3, from = "S:DEL144-144",       to = "S:Y144-")
t3 <- replaceDeletion(x = t3, from = "S:DEL144-145",       to = c("S:Y144-", "S:Y145-"))
t3 <- replaceDeletion(x = t3, from = "S:DEL157-158",       to = c("S:F157-", "S:R158-"))
t3 <- replaceDeletion(x = t3, from = "S:DEL157-159",       to = c("S:F157-", "S:R158-", "S:V159-"))
t3 <- replaceDeletion(x = t3, from = "S:DEL212-212",       to = "S:L212-")
t3 <- replaceDeletion(x = t3, from = "S:DEL241-243",       to = c("S:L241-", "S:L242-", "S:A243-"))
t3 <- replaceDeletion(x = t3, from = "S:DEL243-244",       to = c("S:A243-", "S:L244-"))
t3 <- replaceDeletion(x = t3, from = "S:DEL247-253",       to = c("S:S247-", "S:Y248-", "S:L249-", "S:T250-", "S:P251-", "S:G252-", "S:D253-"))
t3 <- replaceDeletion(x = t3, from = "S:DEL681-681",       to = "S:P681-")
t3 <- replaceDeletion(x = t3, from = "S:DEL1072-1072",     to = "S:E1072-")

t3[, 1] <- toupper(t3[, 1])
t3[, 1] <- gsub("A:", "a:", t3[, 1])
t3[, 1] <- gsub("B:", "b:", t3[, 1])

##-- identify any remaining deletions and insertions
grep("del", t3[, 1], ignore.case = TRUE, value = TRUE)
grep("ins", t3[, 1], ignore.case = TRUE, value = TRUE)

rownames(t3) <- NULL
colnames(t3) <- gsub("prevalence.", "", colnames(t3))

t3[is.na(t3)] <- 0

##-- overlapping deletions may create duplicated rows. Sum them.
t4 <- t3[duplicated(t3[, 1]), 1]

for (m in t4) {
    m1 <- t3[t3[, 1] == m, -1]
    m2 <- apply(m1, 2, function(x) sum(as.numeric(x)))
    m3 <- matrix(as.character(m2), byrow = TRUE, nrow = nrow(m1), ncol = ncol(m1))

    t3[t3[, 1] == m, -1] <- m3
    rm(m1, m2, m3)
}

variant <- unique(t3)
variant <- variant[, c("mutation", VOC, VOI, VUM, setdiff(detected, c(VOC, VOI, VUM)))]

## add MUI, MUC, MU.CZ
if (length(MUC) | length(MUI)) {

    variant <- cbind(variant, MUC = "0", MUI = "0")

    for (m in setdiff(c(MUC, MUI), variant[, 1])) {
        variant <- rbind(variant, "0")
        variant[nrow(variant), 1] <- m
    }
    variant[variant[, 1] %in% MUC ,   "MUC"]   <- "1"
    variant[variant[, 1] %in% MUI ,   "MUI"]   <- "1"
}

##-- read the data in
metadata <- read.csv("datafreeze-2022-09-23_12_weeks_good.tsv", sep = "\t")

##-- order by fasta_id
metadata <- metadata[order(metadata$fasta_id), ]

##-- keep only metadata with fasta file associated
metadata <- metadata[!is.na(metadata$fasta_id), ]

##-- mutations from nextclade output
nextclade <- read.delim("datafreeze-2022-09-23_12_weeks_good.nextclade.2.6.0.tsv", sep = "\t")

selFasta <- (nextclade$seqName %in% metadata$fasta_id)
table(selFasta)

nextclade <- nextclade[selFasta, ]

selOverallStatus <- ( (nextclade$qc.overallStatus == "good") | (nextclade$qc.overallStatus == "mediocre") )
table(selOverallStatus)

selTotalMissing <- (nextclade$totalMissing < 0.05 * genomeLength)
table(selTotalMissing)

sel <- selOverallStatus & selTotalMissing
table(sel)
##--> keep all sequences to match others

mutation <- nextclade[, c("seqName", "substitutions", "deletions", "insertions", "aaSubstitutions", "aaDeletions")]

mutation <- cbind(id = sapply(strsplit(mutation$seqName, " "), "[[", 1), mutation)
mutation <- mutation[!duplicated(mutation$id), ]
rownames(mutation) <- mutation$id

##-- pangolin lineage from pangolin UShER output
pangolin <- read.delim("datafreeze-2022-09-23_12_weeks_good.pangolin-4.1.2-UShER-.csv", sep = ",")[, 1:3]

selFasta <- (pangolin$taxon %in% metadata$fasta_id)
table(selFasta)

selNoConflict <- ( (pangolin$conflict < 0.01) | is.na(pangolin$conflict) )
table(selNoConflict)

## keep also the conflicting sequences, they differ in Delta sublineage assignment, only
##pangolin <- pangolin[selNoConflict, ]

colnames(pangolin) <- c("id", "lineage", "conflict")
pangolin <- pangolin[!duplicated(pangolin$id), ]
rownames(pangolin) <- pangolin$id

sel <- sort(intersect(mutation$id, pangolin$id))

mutation <- mutation[sel, ]
pangolin <- pangolin[sel, ]

mutation <- cbind(mutation, lineage = pangolin[, "lineage"])

rm(nextclade, pangolin, sel, selNoConflict, selOverallStatus, selTotalMissing)

##-- metadata
selNextclade <- metadata$fasta_id %in% rownames(mutation)
table(selNextclade)

metadata <- metadata[selNextclade, ]

metadata <- metadata[!duplicated(metadata$fasta_id), ]
rownames(metadata) <- metadata$fasta_id

##-- common stuff
sel <- sort(intersect(rownames(metadata), rownames(mutation)))

metadata <- metadata[sel, ]
mutation <- mutation[sel, ]

metadata <- cbind(metadata, mutation[, -1])
rm(mutation, sel, selNextclade, selFasta)

##-- tables
lineageTable <-  lapply(region, createTableSimple, f = "lineage")

##-- lineage plot
lineageSel     <- selectLevelsToShow(lineageTable, minFreq = minFreqLineage, minProb = minProbLineage)
lineageSummary <- lapply(region, summarizeLevelsToShow, factorTable = lineageTable, factorSel = lineageSel)
lineageSummary <- do.call(rbind, lineageSummary)
colnames(lineageSummary)[1] <- "Lineage"

t1 <- read.xls("palettes.xls", sheet = "palette variant anotation")
colVariant <- t1[, "Colour"]

t2 <- setdiff(lineageSel, t1[, "Variant"])
colVariant <- c(colVariant, brewer.pal(length(t2), "Greys"))

names(colVariant) <- c(t1[, "Variant"], t2)
colVariant <- colVariant[levels(lineageSummary$Lineage)]
rm(t1, t2)

pdf("./Pdf/plotLineage.pdf", width = 8.25, height = 5.7)
ggplot(data=lineageSummary, aes(x = Region, y = Prob, fill = Lineage)) +
    geom_col(colour="white", size = 0.1) +
    geom_text(aes(x = Region, y = 0.025, hjust = 1, label = Count)) +
    coord_flip() +
    scale_fill_manual(values=colVariant)
dev.off()

##-- summary tables
lineageSummary <- split(lineageSummary, lineageSummary$Region)
lineageSel     <- c(lineageSel, "Other")
lineageWrite   <- data.frame(Lineage = lineageSel,
                             lapply(lineageSummary, function(x) {
                                 rownames(x) <- x[, "Lineage"];
                             x1 <- x[match(lineageSel, rownames(x)), "Prob"]
                                 x1[is.na(x1)] <- 0
                                 return(signif(x1, 3))
                              }),
                              Total = sapply(lineageSummary["Ceska republika"], function(x) {
                                  rownames(x) <- x[, "Lineage"];
                                  x1 <- x[match(lineageSel, rownames(x)), "Freq"]
                                  x1[is.na(x1)] <- 0
                                  return(x1)
                              }), check.names = FALSE)
colnames(lineageWrite)[ncol(lineageWrite)] <- "Total"
write.table(lineageWrite, file = "./Tables/tableLineage.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
cat(c("\"Count\"", sapply(lineageSummary, function(x) x[1, "Count"])[colnames(lineageWrite)[-c(1, ncol(lineageWrite))]], "\n"),
    sep = "\t", append = TRUE, file = "./Tables/tableLineage.tsv")

##-- all lineages in the Czech republic
lineageAll        <- selectLevelsToShow(lineageTable, minFreq = 0, minProb = 0)
lineageAllSummary <- lapply(region, summarizeLevelsToShow, factorTable = lineageTable, factorSel = lineageAll)

lineageAll      <- c(lineageAll, "Other")
lineageAllWrite <- data.frame(Lineage = lineageAll,
                              lapply(lineageAllSummary, function(x) {
                                  rownames(x) <- x[, "Factor"];
                                  x1 <- x[match(lineageAll, rownames(x)), "Prob"]
                                  x1[is.na(x1)] <- 0
                                  return(signif(x1, 3))
                              }),
                              Total = sapply(lineageAllSummary["Ceska republika"], function(x) {
                                  rownames(x) <- x[, "Factor"];
                                  x1 <- x[match(lineageAll, rownames(x)), "Freq"]
                                  x1[is.na(x1)] <- 0
                                  return(x1)
                              }), check.names = FALSE)
colnames(lineageAllWrite)[ncol(lineageAllWrite)] <- "Total"
lineageAllWrite[nrow(lineageAllWrite), "Total"] <- ""
write.table(lineageAllWrite, file = "./Tables/tableAllLineage.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
cat(c("\"Count\"", sapply(lineageAllSummary, function(x) x[1, "Count"])[colnames(lineageAllWrite)[-c(1, ncol(lineageAllWrite))]], "\n"),
    sep = "\t", append = TRUE, file = "./Tables/tableAllLineage.tsv")


##-- time course of the mutations
## add year_week to metadata
week <- paste0("tyden_", strftime(as.Date(metadata$collection_date, format = "%Y-%m-%d"), format = "%G_%V"))
metadata <- cbind(metadata, week)

week <- sort(unique(week))
if ("tyden_NA" %in% week)
    week <- c("tyden_NA", setdiff(week, "tyden_NA"))
week <- c("Vsechny tydny", week)
week <- gsub("tyden_NA", "_unknown_", week)
names(week) <- week

lineageTableWeek <-  lapply(week, createTableSimpleWeek, f = "lineage")

##-- lineage plot
lineageSel     <- selectWeekLevelsToShow(lineageTableWeek, minFreq = minFreqLineage, minProb = minProbLineage)
lineageSummary <- lapply(week, summarizeWeekLevelsToShow, factorTable = lineageTableWeek, factorSel = lineageSel)
lineageSummary <- do.call(rbind, lineageSummary)
colnames(lineageSummary)[1] <- "Lineage"

t1 <- read.xls("palettes.xls", sheet = "palette variant anotation")
colVariant <- t1[, "Colour"]

t2 <- setdiff(lineageSel, t1[, "Variant"])
colVariant <- c(colVariant, brewer.pal(length(t2), "Greys"))

names(colVariant) <- c(t1[, "Variant"], t2)
colVariant <- colVariant[levels(lineageSummary$Lineage)]
rm(t1, t2)

pdf("./Pdf/plotLineageWeek.pdf", width = 8.25, height = 5.7)
ggplot(data=lineageSummary, aes(x = Week, y = Prob, fill = Lineage)) +
    geom_col(colour="white", size = 0.1) +
    geom_text(aes(x = Week, y = 0.025, hjust = 1, label = Count)) +
    coord_flip() +
    scale_fill_manual(values=colVariant)
dev.off()

##-- tables
lineageSummary <- split(lineageSummary, lineageSummary$Week)
lineageSel     <- c(lineageSel, "Other")
lineageWrite   <- data.frame(Lineage = lineageSel,
                             lapply(lineageSummary, function(x) {
                                 rownames(x) <- x[, "Lineage"];
                                 x1 <- x[match(lineageSel, rownames(x)), "Prob"]
                                 x1[is.na(x1)] <- 0
                                 return(signif(x1, 3))
                             }),
                             Total = sapply(lineageSummary["Vsechny tydny"], function(x) {
                                 rownames(x) <- x[, "Lineage"];
                                 x1 <- x[match(lineageSel, rownames(x)), "Freq"]
                                 x1[is.na(x1)] <- 0
                                 return(x1)
                             }), check.names = FALSE)
colnames(lineageWrite)[ncol(lineageWrite)] <- "Total"
##lineageWrite[nrow(lineageWrite), "Total"] <- ""
write.table(lineageWrite, file = "./Tables/tableLineageWeek.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
cat(c("\"Count\"", sapply(lineageSummary, function(x) x[1, "Count"])[colnames(lineageWrite)[-c(1, ncol(lineageWrite))]], "\n"),
    sep = "\t", append = TRUE, file = "./Tables/tableLineageWeek.tsv")

##-- all lineages in the Czech republic
lineageAll        <- selectWeekLevelsToShow(lineageTableWeek, minFreq = 0, minProb = 0)
lineageAllSummary <- lapply(week, summarizeWeekLevelsToShow, factorTable = lineageTableWeek, factorSel = lineageAll)

lineageAll      <- c(lineageAll, "Other")
lineageAllWrite <- data.frame(Lineage = lineageAll,
                              lapply(lineageAllSummary, function(x) {
                                  rownames(x) <- x[, "Factor"];
                                  x1 <- x[match(lineageAll, rownames(x)), "Prob"]
                                  x1[is.na(x1)] <- 0
                                  return(signif(x1, 3))
                              }),
                              Total = sapply(lineageAllSummary["Vsechny tydny"], function(x) {
                                  rownames(x) <- x[, "Factor"];
                                  x1 <- x[match(lineageAll, rownames(x)), "Freq"]
                                  x1[is.na(x1)] <- 0
                                  return(x1)
                              }), check.names = FALSE)
colnames(lineageAllWrite)[ncol(lineageAllWrite)] <- "Total"
##lineageAllWrite[nrow(lineageAllWrite), "Total"] <- ""
write.table(lineageAllWrite, file = "./Tables/tableAllLineageWeek.tsv", sep = "\t", col.names = TRUE, row.names = FALSE)
cat(c("\"Count\"", sapply(lineageAllSummary, function(x) x[1, "Count"])[colnames(lineageAllWrite)[-c(1, ncol(lineageAllWrite))]], "\n"),
    sep = "\t", append = TRUE, file = "./Tables/tableAllLineageWeek.tsv")

##-- create output from Martin's data files in subdir mkolisko
##-- by region
aaMutationSummaryWide  <- as.matrix(read.csv("mkolisko/table_regions.reduced.csv", row.names = 1, sep = ","))

##-- remove Orf9b
aaMutationSummaryWide <- aaMutationSummaryWide[!grepl("ORF9b:", rownames(aaMutationSummaryWide)), ]

##-- I do not need pocet column
aaMutationSummaryWide  <- aaMutationSummaryWide[, -ncol(aaMutationSummaryWide)]
aaMutationSummaryCount <- as.numeric(sapply(strsplit(colnames(aaMutationSummaryWide), "_"), "[[", 2))
colnames(aaMutationSummaryWide) <- names(aaMutationSummaryCount) <-
    gsub("X", "Neurceno", sapply(strsplit(colnames(aaMutationSummaryWide), "_"), "[[", 1))

##-- add the Czech Republic
aaMutationSummaryWide <- cbind(aaMutationSummaryWide,
                               `Ceska republika` = colSums(t(aaMutationSummaryWide) * aaMutationSummaryCount) / sum(aaMutationSummaryCount))
aaMutationSummaryCount <- c(aaMutationSummaryCount, `Ceska republika` = sum(aaMutationSummaryCount))

##-- "_unknown_" is "Neurceno" now
region <- gsub("_unknown_", "Neurceno", region)
names(region) <- region

##-- if some region has no data, add NA column
t1 <- setdiff(region, colnames(aaMutationSummaryWide))

for (t2 in t1) {
    aaMutationSummaryWide <- cbind(aaMutationSummaryWide, NA)
    aaMutationSummaryCount <- c(aaMutationSummaryCount, 0)
    colnames(aaMutationSummaryWide)[ncol(aaMutationSummaryWide)] <- t2
    names(aaMutationSummaryCount)[length(aaMutationSummaryCount)] <- t2
}
rm(t1, t2)

##-- sort
aaMutationSummaryWide  <- aaMutationSummaryWide[, region]
aaMutationSummaryCount <- aaMutationSummaryCount[region]

otherMutations <- aaMutationSummaryWide["other", ]

aaMutationSummaryCount <- aaMutationSummaryCount[setdiff(names(aaMutationSummaryCount), "other")]
aaMutationSummaryWide  <- aaMutationSummaryWide[setdiff(rownames(aaMutationSummaryWide), "other"), ]

protein <- factor(gsub(":.*", "", rownames(aaMutationSummaryWide)),
                  levels = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF10", "ORF14"))
positionInProtein <- as.numeric(sapply(strsplit(rownames(aaMutationSummaryWide), ":"), function(x) substr(x[2], 2, nchar(x[2]) - 1)))
aaMutationSummaryOrder <- order(protein, positionInProtein)
names(protein) <- names(positionInProtein) <- names(aaMutationSummaryOrder) <- rownames(aaMutationSummaryWide)

t1 <- read.xls("palettes.xls", sheet = "palette protein anotation")
colProtein <- t1[, "Colour"]
t2 <- setdiff(levels(protein), t1[, "Protein"])
colProtein <- c(colProtein, brewer.pal(length(t2), "Greys"))
names(colProtein) <- c(t1[, "Protein"], t2)
colProtein <- colProtein[levels(protein)]
rm(t1, t2)

colFun  <- colorRamp2(heatmapSplit, c("white", "cornflowerblue", "yellow", "red"))
colFun2 <- colorRamp2(c(0, 1/3, 2/3, 1) * max(otherMutations, na.rm = TRUE), c("#EDF8FB", "#B3CDE3", "#8C96C6", "#88419D"))
##sub("(\\d)[^0-9]+$", "\\1", df1$v1)
variantHeatmap <- variant[match(rownames(aaMutationSummaryWide), variant[, "mutation"]),
                          c(intersect(colnames(variant), c(lineageSel)))]
variantHeatmap[is.na(variantHeatmap)] <- 0
rownames(variantHeatmap) <- rownames(aaMutationSummaryWide)
variantHeatmap <- as.matrix(variantHeatmap)
storage.mode(variantHeatmap) <- "numeric"

##-- sort in a more meaningfull manner
variantHeatmap <- variantHeatmap[, c("BA.2","BA.2.9", "BA.2.12.1",
                                     "BA.4", "BA.4.1", "BA.4.6",
                                     "BA.5", "BA.5.1", "BA.5.1.2", "BA.5.1.3", "BA.5.2", "BA.5.2.1", "BF.1", "BF.5", "BA.5.2.3",
                                     "BA.5.3.1", "BE.1", "BE.1.1")]

pdf("./Pdf/heatmapAAMutation.pdf", height = 25, width = 11.7)
Heatmap(aaMutationSummaryWide, name = "Frekvence", rect_gp = gpar(col = "white", lwd = 2), na_col = "white", col = colFun,
        cluster_columns = FALSE,
        row_order = aaMutationSummaryOrder,
        row_split = protein,
        left_annotation  = rowAnnotation(Protein = protein, col = list(Protein = colProtein), show_legend = FALSE),
        top_annotation = columnAnnotation(`Dalsi mutace` = otherMutations, col = list(`Dalsi mutace` = colFun2), na_col = "white"),
        bottom_annotation = columnAnnotation(`Pocet vzorku` = anno_text(aaMutationSummaryCount, just = "right", location = 1))) +
    Heatmap(variantHeatmap, name = "VOC/VOI/VUM", rect_gp = gpar(col = "white", lwd = 2), na_col = "white",
            col = colorRamp2(c(0, 0.5, 1), c("#EEEEEE", "#EEEEEE", "#222222")),
            cluster_columns = FALSE, width = 12)
dev.off()

##-- spike glycosylation
glyco <- c("FUL","BMA","XYL","SIA","NDG","NAG","MAN","FUC","GAL")

t2 <- lapply(glyco, function(v) {
    v1 <- read.delim(paste0("./VariantGlycoSig/", v, ".tsv"), sep = "\t", header = FALSE)
    v2 <- data.frame(position = v1[, 1], glycosylation = v, occurence = as.numeric(v1[, 2]))
})

t2 <- do.call(rbind, t2)
t2[, "position"] <- toupper(t2[, "position"])

t2$position <- str_replace_all(t2$position,
                               c(
                                   "ALA"="A", "ARG"="R", "ASN"="N", "ASP"="D",
                                   "CYS"="C", "GLU"="E", "GLN"="Q", "GLY"="G",
                                   "HIS"="H", "ILE"="I", "LEU"="L", "LYS"="K",
                                   "MET"="M", "PHE"="F", "PRO"="P", "SER"="S",
                                   "THR"="T", "TRP"="W", "TYR"="Y", "VAL"="V"
                                   ))

t3 <- reshape(t2, idvar = "position", timevar = "glycosylation", direction = "wide")
t3 <- as.matrix(t3)

rownames(t3) <- NULL
colnames(t3) <- gsub("occurence.", "", colnames(t3))
t3[is.na(t3)] <- 0

variant_glyco <- unique(t3)
variant_glyco <- variant_glyco[, c("position",glyco)]

variantHeatmap_glyco <- variant_glyco[match(sub("(\\d)[^0-9]+$", "\\1",rownames(aaMutationSummaryWide)), variant_glyco[, "position"]),
                                      c(intersect(colnames(variant_glyco), c(glyco)))]
variantHeatmap_glyco[is.na(variantHeatmap_glyco)] <- 0
rownames(variantHeatmap_glyco) <- rownames(aaMutationSummaryWide)
variantHeatmap_glyco <- as.matrix(variantHeatmap_glyco)
storage.mode(variantHeatmap_glyco) <- "numeric"

variantHeatmap_glyco <- cbind(variantHeatmap_glyco,Glykosylace=
                                                     c(ifelse(rowSums(variantHeatmap_glyco)>1,1,0)))
variantHeatmap_glyco <- as.matrix(variantHeatmap_glyco[,ncol(variantHeatmap_glyco)])

##-- sort in a more meaningfull manner
##variantHeatmap_glyco <- variantHeatmap_glyco[, c("FUL","GAL","NDG","XYL","BMA","SIA","MAN","FUC","NAG")]

sel <- grep("S:", rownames(aaMutationSummaryWide), value = TRUE)

pdf("./Pdf/heatmapAAMutation_Spike_glyco.pdf", height = 8.4, width = 11.7)
Heatmap(aaMutationSummaryWide[sel, ], name = "Frekvence", rect_gp = gpar(col = "white", lwd = 2), na_col = "white", col = colFun,
        cluster_columns = FALSE,
        row_order = order(as.numeric(gsub("S:.(\\d+).", "\\1",  sel))),
        row_split = protein[sel],
        left_annotation  = rowAnnotation(Protein = protein[sel], col = list(Protein = colProtein), show_legend = FALSE),
        bottom_annotation = columnAnnotation(`Pocet vzorku` = anno_text(aaMutationSummaryCount, just = "right", location = 1))) +
    Heatmap(variantHeatmap[sel, ], name = "VOC/VOI/VUM", rect_gp = gpar(col = "white", lwd = 2), na_col = "white",
            col = colorRamp2(c(0, 0.5, 1), c("#EEEEEE", "#EEEEEE", "#222222")),
            cluster_columns = FALSE, width = 12) +
    Heatmap(variantHeatmap_glyco[sel, ], name = "Glykosylace", rect_gp = gpar(col = "white", lwd = 2), na_col = "white",
            col = colorRamp2(c(0, 0.5, 1), c("#EEEEEE", "#EEEEEE", "#660066")),show_heatmap_legend = FALSE,
            cluster_columns = FALSE, width = 1)
dev.off()

##-- by week
aaMutationSummaryWide  <- as.matrix(read.csv("mkolisko/table_weeks.reduced.csv", row.names = 1, sep = ","))

##-- remove Orf9b
aaMutationSummaryWide <- aaMutationSummaryWide[!grepl("ORF9b:", rownames(aaMutationSummaryWide)), ]

##-- I do not need pocet column
aaMutationSummaryWide  <- aaMutationSummaryWide[, -ncol(aaMutationSummaryWide)]

aaMutationSummaryCount <- as.numeric(sapply(strsplit(colnames(aaMutationSummaryWide), "_"), "[[", 2))
colnames(aaMutationSummaryWide) <- names(aaMutationSummaryCount) <-
    gsub("NA", "_unknown_", sapply(strsplit(colnames(aaMutationSummaryWide), "_"), "[[", 1))
colnames(aaMutationSummaryWide) <- names(aaMutationSummaryCount) <-
    gsub("X(\\d{2})\\.(\\d{4})", "tyden_\\2_\\1", colnames(aaMutationSummaryWide))

##-- add the All weeks
aaMutationSummaryWide <- cbind(aaMutationSummaryWide,
                               `Vsechny tydny` = colSums(t(aaMutationSummaryWide) * aaMutationSummaryCount) / sum(aaMutationSummaryCount))
aaMutationSummaryCount <- c(aaMutationSummaryCount, `Vsechny tydny` = sum(aaMutationSummaryCount))

##-- sort alphabetically
aaMutationSummaryWide  <- aaMutationSummaryWide[, week]
aaMutationSummaryCount <- aaMutationSummaryCount[week]

otherMutations <- aaMutationSummaryWide["other", ]

aaMutationSummaryCount <- aaMutationSummaryCount[setdiff(names(aaMutationSummaryCount), "other")]
aaMutationSummaryWide  <- aaMutationSummaryWide[setdiff(rownames(aaMutationSummaryWide), "other"), ]

protein <- factor(gsub(":.*", "", rownames(aaMutationSummaryWide)),
                  levels = c("ORF1a", "ORF1b", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF9b", "ORF10", "ORF14"))
positionInProtein <- as.numeric(sapply(strsplit(rownames(aaMutationSummaryWide), ":"), function(x) substr(x[2], 2, nchar(x[2]) - 1)))
aaMutationSummaryOrder <- order(protein, positionInProtein)
names(protein) <- names(positionInProtein) <- names(aaMutationSummaryOrder) <- rownames(aaMutationSummaryWide)

t1 <- read.xls("palettes.xls", sheet = "palette protein anotation")
colProtein <- t1[, "Colour"]
t2 <- setdiff(levels(protein), t1[, "Protein"])
colProtein <- c(colProtein, brewer.pal(length(t2), "Greys"))
names(colProtein) <- c(t1[, "Protein"], t2)
colProtein <- colProtein[levels(protein)]
rm(t1, t2)

colFun  <- colorRamp2(heatmapSplit, c("white", "cornflowerblue", "yellow", "red"))
colFun2 <- colorRamp2(c(0, 1/3, 2/3, 1) * max(otherMutations), c("#EDF8FB", "#B3CDE3", "#8C96C6", "#88419D"))

variantHeatmap <- variant[match(rownames(aaMutationSummaryWide), variant[, "mutation"]),
                          c(intersect(colnames(variant), c(lineageSel)))]
variantHeatmap[is.na(variantHeatmap)] <- 0
rownames(variantHeatmap) <- rownames(aaMutationSummaryWide)
variantHeatmap <- as.matrix(variantHeatmap)
storage.mode(variantHeatmap) <- "numeric"

##-- sort in a more meaningfull manner
variantHeatmap <- variantHeatmap[, c("BA.2","BA.2.9", "BA.2.12.1",
                                     "BA.4", "BA.4.1", "BA.4.6",
                                     "BA.5", "BA.5.1", "BA.5.1.2", "BA.5.1.3", "BA.5.2", "BA.5.2.1", "BF.1", "BF.5", "BA.5.2.3",
                                     "BA.5.3.1", "BE.1", "BE.1.1")]


pdf("./Pdf/heatmapAAMutationWeek.pdf", height = 21, width = 11.7)
Heatmap(aaMutationSummaryWide, name = "Frekvence", rect_gp = gpar(col = "white", lwd = 2), na_col = "white", col = colFun,
        cluster_columns = FALSE,
        row_order = aaMutationSummaryOrder,
        row_split = protein,
        left_annotation  = rowAnnotation(Protein = protein, col = list(Protein = colProtein), show_legend = FALSE),
        top_annotation = columnAnnotation(`Dalsi mutace` = otherMutations, col = list(`Dalsi mutace` = colFun2)),
        bottom_annotation = columnAnnotation(`Pocet vzorku` = anno_text(aaMutationSummaryCount, just = "right", location = 1))) +
    Heatmap(variantHeatmap, name = "VOC/VOI/VUM", rect_gp = gpar(col = "white", lwd = 2), na_col = "white",
            col = colorRamp2(c(0, 0.5, 1), c("#EEEEEE", "#EEEEEE", "#222222")),
            cluster_columns = FALSE, width = 12)
dev.off()

##-- Spike glycosylation week
variantHeatmap_glyco <- variant_glyco[match(sub("(\\d)[^0-9]+$", "\\1",rownames(aaMutationSummaryWide)), variant_glyco[, "position"]),
                                      c(intersect(colnames(variant_glyco), c(glyco)))]
variantHeatmap_glyco[is.na(variantHeatmap_glyco)] <- 0
rownames(variantHeatmap_glyco) <- rownames(aaMutationSummaryWide)
variantHeatmap_glyco <- as.matrix(variantHeatmap_glyco)
storage.mode(variantHeatmap_glyco) <- "numeric"

##-- sort in a more meaningfull manner
##variantHeatmap_glyco <- variantHeatmap_glyco[, c("FUL","GAL","NDG","XYL","BMA","SIA","MAN","FUC","NAG")]
variantHeatmap_glyco <- cbind(variantHeatmap_glyco,Glykosylace=
                              c(ifelse(rowSums(variantHeatmap_glyco)>1,1,0)))
variantHeatmap_glyco <- as.matrix(variantHeatmap_glyco[,ncol(variantHeatmap_glyco)])

sel <- grep("S:", rownames(aaMutationSummaryWide), value = TRUE)

pdf("./Pdf/heatmapAAMutationWeek_Spike_glyco_week.pdf", height = 8.4, width = 11.7)
Heatmap(aaMutationSummaryWide[sel, ], name = "Frekvence", rect_gp = gpar(col = "white", lwd = 2), na_col = "white", col = colFun,
        cluster_columns = FALSE,
        row_order = order(as.numeric(gsub("S:.(\\d+).", "\\1",  sel))),
        row_split = protein[sel],
        left_annotation  = rowAnnotation(Protein = protein[sel], col = list(Protein = colProtein), show_legend = FALSE),
        bottom_annotation = columnAnnotation(`Pocet vzorku` = anno_text(aaMutationSummaryCount, just = "right", location = 1))) +
    Heatmap(variantHeatmap[sel, ], name = "VOC/VOI/VUM", rect_gp = gpar(col = "white", lwd = 2), na_col = "white",
            col = colorRamp2(c(0, 0.5, 1), c("#EEEEEE", "#EEEEEE", "#222222")),
            cluster_columns = FALSE, width = 12) +
    Heatmap(variantHeatmap_glyco[sel, ], name = "Glykosylace", rect_gp = gpar(col = "white", lwd = 2), na_col = "white",
            col = colorRamp2(c(0, 0.5, 1), c("#EEEEEE", "#EEEEEE", "#660066")), show_heatmap_legend = FALSE,
            cluster_columns = FALSE, width = 1)
dev.off()

##-- Spike
sel <- grep("S:", rownames(aaMutationSummaryWide), value = TRUE)

pdf("./Pdf/heatmapAAMutationWeek_Spike.pdf", height = 8.4, width = 11.7)
Heatmap(aaMutationSummaryWide[sel, ], name = "Frekvence", rect_gp = gpar(col = "white", lwd = 2), na_col = "white", col = colFun,
        cluster_columns = FALSE,
        row_order = order(as.numeric(gsub("S:.(\\d+).", "\\1",  sel))),
        row_split = protein[sel],
        left_annotation  = rowAnnotation(Protein = protein[sel], col = list(Protein = colProtein), show_legend = FALSE),
        bottom_annotation = columnAnnotation(`Pocet vzorku` = anno_text(aaMutationSummaryCount, just = "right", location = 1))) +
    Heatmap(variantHeatmap[sel, ], name = "VOC/VOI/VUM", rect_gp = gpar(col = "white", lwd = 2), na_col = "white",
            col = colorRamp2(c(0, 0.5, 1), c("#EEEEEE", "#EEEEEE", "#222222")),
            cluster_columns = FALSE, width = 12)
dev.off()

##-- session info
save.image("mutationTable.RData")
sessionInfo()
warnings()
traceback()
sink()
