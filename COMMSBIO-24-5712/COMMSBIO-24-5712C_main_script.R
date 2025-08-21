### Single-cell DNA methylation analysis tool Amethyst resolves distinct non-CG methylation patterns in human astrocytes and oligodendrocytes ###
# Author: Lauren Rylaarsdam, PhD 
# Date: July 2023 - August 2025

# ***Please cite our manuscript if you use Amethyst or any of the following code.***
# This script encompasses analysis performed for manuscript COMMSBIO-24-5712C.
# Accompanying data to this script can be found under: NCBI GEO GSE303678.

# This analysis was done over the course of multiple years, rounds of revisions, and updates (v0.0.0.9000 - v1.0.2).
# Because of this, there may be minor inconsistencies with the tutorials published for the current version.
# Please see Github for the most up-to-date workflows and function parameters. 

# Load libraries
library(amethyst)
library(data.table)
library(ggplot2)
library(ggrepel)
library(tibble)
library(tidyr)
library(plyr)
library(dplyr)
library(future)
library(furrr)
library(purrr)
library(stringr)
library(pheatmap)
library(cowplot)
library(webr)
library(rrvgo)

###############################################################################
### BRAIN DATA PROCESSING ###

# PART ONE: DONE ON SERVER
# load object and perform basic filtering 
brain <- createObject()
# add data from cellInfo file
brain <- addCellInfo(brain, file = "/secret/path/sciMETv2_3842F_reAnalysis/sciMETv2_3842F.all.cellInfo.txt")
# add batch annotation information
brain <- addAnnot(brain, "/secret/path/sciMETv2_3842F_reAnalysis/sciMETv2_3842F.all.batch.annot", name = "batch")
# filter out cells with high/low coverage or mch_pct beyond what is biologically reasonable
brain@metadata <- brain@metadata |> dplyr::filter(cov > 6000000 & cov < 100000000 & mch_pct < 12) |> dplyr::filter(!(batch %in% c("SLN_1", "SLN_2")))
brain@h5paths <- data.frame(barcode = rownames(brain@metadata), 
                            path = rep("/secret/path/sciMETv2_3842F_reAnalysis/sciMETv2_3842F_merged.h5", length(rownames(brain@metadata))))

# index chromosomes and calculate methylation levels over fixed genomic windows
brain@index[["chr_cg"]] <- indexChr(brain, type = "CG", threads = 20)
brain@index[["chr_ch"]] <- indexChr(brain, type = "CH", threads = 20)
brain@genomeMatrices[["cg_100k_score"]] <- makeWindows(brain, stepsize = 100000, type = "CG", metric = "score", threads = 20, index = "chr_cg", nmin = 2, species = "human")
brain@genomeMatrices[["ch_100k_pct"]] <- makeWindows(brain, stepsize = 100000, type = "CH", metric = "percent", threads = 20, index = "chr_ch", nmin = 5, species = "human")

brain@genomeMatrices[["ch_100k_pct"]] <- brain@genomeMatrices[["ch_100k_pct"]][, colnames(brain@genomeMatrices[["ch_100k_pct"]]) %in% rownames(brain@metadata)]
dimEstimate(brain, genomeMatrices = "ch_100k_pct", dims = 40, threshold = 0.989) # 23
brain@genomeMatrices[["cg_100k_score"]] <- brain@genomeMatrices[["cg_100k_score"]][, colnames(brain@genomeMatrices[["cg_100k_score"]]) %in% rownames(brain@metadata)]
dimEstimate(brain, genomeMatrices = "cg_100k_score", dims = 40, threshold = 0.982) # 26
brain@reductions[["irlba"]] <- runIrlba(brain, genomeMatrices = c("ch_100k_pct", "cg_100k_score"), dims = c(23, 26), replaceNA = c(0, 0))
brain@reductions[["irlba"]] <- regressCovBias(brain, reduction = "irlba")
brain <- runCluster(brain, k = 25, method = "louvain", reduction = "irlba", colname = "cluster_id")
brain@reductions[["umap_100k_cg_ch"]] <- runUmap(brain, neighbors = 25, dist = 0.1, method = "euclidean", reduction = "irlba") 

# clustering w/ cg and ch 100kb windows only for figure 2
brain@reductions[["irlba_cg_100k"]] <- brain@reductions[["irlba"]][, c(24:49)]
brain@reductions[["irlba_ch_100k"]] <- brain@reductions[["irlba"]][, c(1:23)]
set.seed(4321)
brain@reductions[["umap_irlba_cg_100k"]] <- runUmap(brain, neighbors = 25, dist = 0.1, method = "euclidean", reduction = "irlba_cg_100k") 
set.seed(1234)
brain@reductions[["umap_irlba_ch_100k"]] <- runUmap(brain, neighbors = 25, dist = 0.1, method = "euclidean", reduction = "irlba_ch_100k")

# calculate extended mCH 
protein_coding <- unique(brain@ref |> dplyr::filter(type == "gene" & gene_type == "protein_coding") |> dplyr::pull(gene_name))
brain@genomeMatrices[["gene_ch_extended"]] <- makeWindows(brain, type = "CH", genes = protein_coding, metric = "percent", index = "chr_ch", threads = 20) 
saveRDS(brain@genomeMatrices[["gene_ch_extended"]], "brain_mCH_pct_20k_genes.RData")

# make cluster tracks
brain_cluster500bpwindows <- calcSmoothedWindows(brain, type = "CG", threads = 30, step = 500, species = "human", 
                                           index = "chr_cg", groupBy = "cluster_id", returnSumMatrix = T, returnPctMatrix = T)
brain@genomeMatrices[["cg_cluster_tracks"]] <- brain_cluster500bpwindows[["pct_matrix"]] 

# calculate DMRs from cluster tracks
brain_cg_dmrs <- testDMR(sumMatrix = brain_cluster500bpwindows[["sum_matrix"]], eachVsAll = TRUE, nminTotal = 10, nminGroup = 10)
brain_cg_filtered_dmrs <- filterDMR(brain_cg_dmrs, method = "bonferroni", filter = TRUE, pThreshold = 0.01, logThreshold = 2)
brain_cg_collapsed_dmrs <- collapseDMR(brain, brain_cg_filtered_dmrs, maxDist = 2000, minLength = 2000, reduce = T, annotate = T) 
brain_cg_collapsed_dmrs$type <- dplyr::recode(brain_cg_collapsed_dmrs$member_id,
                                          "1" = "Oligo",
                                          "2" = "Oligo",
                                          "3" = "Oligo",
                                          "4" = "Astro",
                                          "5" = "OPC",
                                          "6" = "Inh_CGE",
                                          "7" = "Exc_L2-4_RORB",
                                          "8" = "Exc_L4-6_LRRK1",
                                          "9" = "Micro",
                                          "10" = "Exc_L6_TLE4",
                                          "11" = "Exc_L1-3_CUX2",
                                          "12" = "Exc_L5-6_PDZRN4", 
                                          "13" = "Inh_MGE",
                                          "14" = "Exc_L4-5_FOXP2") 
# brain@results[["brain_cg_collapsed_dmrs"]] <- brain_cg_collapsed_dmrs

brain_top_cg_dmrs <- brain_cg_collapsed_dmrs |> 
  dplyr::group_by(member_id, direction) |> 
  dplyr::arrange(dmr_padj, .by_group = TRUE) |> dplyr::mutate(rank_padj = 1:n()) |>
  dplyr::arrange(desc(abs(dmr_logFC)), .by_group = TRUE) |> dplyr::mutate(rank_logFC = 1:n()) |>
  rowwise() |> dplyr::mutate(total_rank = sum(rank_padj, rank_logFC)) |> 
  group_by(member_id, direction) |> slice_min(n = 5, order_by = total_rank) |>
  dplyr::mutate(location = paste0(chr, "_", (dmr_start - 1000), "_", (dmr_end - 1000))) |> dplyr::arrange(direction)

# classification with genehancer (brain)
genehancer <- data.table(read.delim("~/Desktop/mount/hg38_genehancer.class.bed", header=FALSE))
setnames(genehancer, c("chr", "start", "end", "class"))
setkey(genehancer, chr, start, end)
brain_dmrs <- copy(brain_cg_collapsed_dmrs)
setkey(brain_dmrs, chr, dmr_start, dmr_end)
brain_overlap_result <- foverlaps(brain_dmrs, genehancer, type="any", nomatch=NA)
# Remove rows that are duplicated by test, direction, dmr_start, dmr_end, and class
brain_overlap_result <- unique(brain_overlap_result, by = c("member_id", "direction", "dmr_start", "dmr_end", "class"))
brain_overlap_result <- brain_overlap_result |> dplyr::group_by(type, direction, class) |> dplyr::summarise(n = n(), dmr_sum = sum(end - start))
brain_overlap_result <- brain_overlap_result |> dplyr::mutate(other = ifelse(is.na(class), "YES", "NO"))
ggplot(brain_overlap_result, aes(y = type, x = n, fill = class)) + geom_col() + facet_grid(vars(direction), scales = "free_y") + scale_fill_manual(values = c("#eb4034", "#52a1a3", "#f6cb52")) + theme_classic()
brain_overlap_result |> dplyr::group_by(class) |>  dplyr::summarise(n = sum(n)) |>dplyr::mutate(percent = 100 * n / sum(n)) # 83.1% are not NA

# topGO of DMR results (cg brain - example: hyomethylated genes in Exc L2-4 RORB)
background <- unique(brain@ref |> dplyr::filter(type == "gene" & seqid != "chrM" & seqid != "chrY") |> dplyr::pull(gene_name)) 
brain_resultElim <- as.list(c("Astro", "Oligo", "OPC", "Micro", "Inh_CGE", "Inh_MGE", "Exc_L1-3_CUX2", "Exc_L2-4_RORB", "Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L5-6_PDZRN4", "Exc_L6_TLE4"))
direction <- as.list(c("hypo", "hyper"))

for (i in brain_resultElim) {
  for (j in direction) {
    tryCatch({
      query <- unlist(strsplit(brain_cg_collapsed_dmrs$gene_names[brain_cg_collapsed_dmrs$type == i & brain_cg_collapsed_dmrs$direction == j], ", "))
      GOdata <- new("topGOdata", 
                    description = "GO Enrichment Analysis", 
                    ontology = "BP", 
                    allGenes = setNames(factor(as.integer(background %in% query), levels = c(0, 1)), background),
                    geneSel = function(x) x == 1, 
                    nodeSize = 10, 
                    annot = annFUN.org, 
                    mapping = "org.Hs.eg.db", 
                    ID = "symbol")
      fisher_res <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
      brain_resultElim[[i]][[j]] <- GenTable(GOdata, Fisher = fisher_res, topNodes = 500, numChar = 60)
      brain_resultElim[[i]][[j]] <- brain_resultElim[[i]][[j]] |> dplyr::filter(Fisher < 0.01 & Significant > 5) |> dplyr::mutate(fold_change = Significant/Expected, Fisher = as.numeric(Fisher))
      brain_resultElim[[i]][[j]] <- brain_resultElim[[i]][[j]] |> dplyr::filter(fold_change > 2)
      brain_resultElim[[i]][[j]] <- janitor::clean_names(brain_resultElim[[i]][[j]])
      brain_resultElim[[i]][[j]]$type <- as.character(i)
      brain_resultElim[[i]][[j]]$direction <- as.character(j)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}

brain_resultElim[1:12] <- NULL
for (i in c(1:12)) {
  brain_resultElim[[i]] <- do.call(rbind, brain_resultElim[[i]])
}
brain_resultElim <- do.call(rbind, brain_resultElim)
brain_resultElim <- brain_resultElim |> group_by(type, direction) |> dplyr::mutate(rank_pval = rank(fisher, ties.method = "min"), 
                                                                                   rank_logfc = rank(-fold_change, ties.method = "min"),
                                                                                   rank_total = rank_pval + rank_logfc) 
ggplot(brain_resultElim |> dplyr::filter(annotated < 150 & expected > 0.3), aes(x = fold_change, y = -log10(fisher), color = type)) + 
  geom_point() + scale_color_manual(values = brain_type_colors) + theme_classic() + ggrepel::geom_text_repel(aes(label = term)) 

# pie plot
brain_simMatrix <- calculateSimMatrix(brain_resultElim$go_id, orgdb="org.Hs.eg.db", ont="BP", method="Rel")
brain_scores <- setNames(-log10(brain_resultElim$fisher), brain_resultElim$go_id)
brain_reducedTerms <- reduceSimMatrix(brain_simMatrix, brain_scores, threshold=0.9, orgdb="org.Hs.eg.db")
treemapPlot(brain_reducedTerms)
PieDonut(brain_reducedTerms, aes(parentTerm, term, count=score), showRatioDonut = F, showPieName = F, r0 = 0.45, r1 = 1.0)

###############################################################################
#### BRAIN LOCAL ####
# basic QC figures
dimFeature(brain, colorBy = cluster_id, pointSize = 0.3, colors = makePalette(option = 19, n = 14)) +
  geom_text_repel(aes(label = cluster_id), color = "black", data =brain@metadata |> dplyr::group_by(cluster_id) |> dplyr::summarise(umap_x = mean(umap_x), umap_y = mean(umap_y)))

# correlate to ref
brain@genomeMatrices[["gene_ch"]] <- brain@genomeMatrices[["gene_ch"]][, colnames(brain@genomeMatrices[["gene_ch"]]) %in% rownames(brain@metadata)]
ref <- readRDS("~/Library/CloudStorage/OneDrive-OregonHealth&ScienceUniversity/amethyst/atlas data/Luo 2022/PMC9004682_5972genes_BA10_ref.RData")
mmc9 <- read_excel("~/Downloads/mmc9.xlsx", skip = 3) # PMC9004682

brain@genomeMatrices[["gene_ch_cluster"]] <- aggregateMatrix(brain, "gene_ch", "cluster_id")
both <- merge(ref, brain@genomeMatrices[["gene_ch_cluster"]], by = 0) |> tibble::column_to_rownames(var = "Row.names")
both <- both[rownames(both) %in% mmc9$`Gene Name`, ]
cor <- cor(both)
cor <- cor[c(1:17), c(18:31)]
pheatmap(cor)

# markers for cell type annotation
heatMap(brain, genes = c("MOG", "MBP", "MAG"), nrow = 3, track = "cg_cluster_tracks") # Oligodendrocytes
heatMap(brain, genes = c("PDGFRA", "OLIG1"), nrow = 2, matrix = "cg_cluster_tracks") # Oligodendrocyte progenitors
heatMap(brain, genes = c("GFAP", "SLC1A2", "AQP4", "S100B", "ALDH1L1", "AGT"), nrow = 3, matrix = "cg_cluster_tracks") # Astrocytes
heatMap(brain, genes = c("C1QA", "TMEM119", "CXCR1", "CSF1R"), nrow = 4, matrix = "cg_cluster_tracks") # Microglia
heatMap(brain, genes = c("GAD1", "SLC32A1", "PVALB", "DLX1", "DLX2", "LHX6"), nrow = 3, matrix = "cg_cluster_tracks") # Gabaergic neurons
heatMap(brain, genes = c("SLC17A7", "GRIN2B", "BCL11B", "SATB2", "TBR1", "NRGN"), nrow = 3, matrix = "cg_cluster_tracks") # Glutamatergic neurons

# assign
annotation <- tidyr::pivot_longer(as.data.frame(cor) %>% rownames_to_column("group2"), cols = c(2:(length(colnames(cor))+1)), names_to = "group1", values_to = "cor") %>% 
  dplyr::group_by(group1) %>% dplyr::arrange(desc(cor), .by_group = TRUE) %>% top_n(n = 1)

# call cell type
brain@metadata[["type"]] <- dplyr::recode(brain@metadata[["cluster_id"]],
                                          "1" = "Oligo",
                                          "2" = "Oligo",
                                          "3" = "Oligo",
                                          "4" = "Astro",
                                          "5" = "OPC",
                                          "6" = "Inh_CGE",
                                          "7" = "Exc_L2-4_RORB",
                                          "8" = "Exc_L4-6_LRRK1",
                                          "9" = "Micro",
                                          "10" = "Exc_L6_TLE4",
                                          "11" = "Exc_L1-3_CUX2",
                                          "12" = "Exc_L5-6_PDZRN4", 
                                          "13" = "Inh_MGE",
                                          "14" = "Exc_L4-5_FOXP2") 
brain@metadata$type <- factor(brain@metadata$type, levels = c("Astro", "Oligo", "OPC", "Micro", "Inh_CGE", "Inh_MGE", "Exc_L1-3_CUX2", "Exc_L2-4_RORB", "Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L5-6_PDZRN4", "Exc_L6_TLE4"))
brain_type_colors <- c("Astro" =  "#C41858",
                 "Micro" =  "#F79C46",
                 "Oligo" =   "#E95E4E",
                 "OPC" = "#fc9286",
                 "Inh_CGE" = "#E5C05E",
                 "Inh_MGE" = "#B2C897",
                 "Exc_L1-3_CUX2" = "#aed6d6",
                 "Exc_L2-4_RORB" =  "#5CBCC4",
                 "Exc_L4-5_FOXP2" =   "#176e6b",
                 "Exc_L4-6_LRRK1" =    "#1E2C34",  
                 "Exc_L5-6_PDZRN4" =   "#8C1E4C",
                 "Exc_L6_TLE4" =  "#d66997")

# create aggregated cluster tracks - need to make function to do this automatically
brain@tracks[["cg_type_tracks"]] <- copy(brain@tracks[["cg_cluster_tracks"]])
setnames(brain@tracks[["cg_type_tracks"]], c("chr", "start", "end", "Oligo_1", "Oligo_2", "Oligo_3", "Astro", "OPC", "Inh_CGE", "Exc_L2-4_RORB", "Exc_L4-6_LRRK1", "Micro", "Exc_L6_TLE4", "Exc_L1-3_CUX2", "Exc_L5-6_PDZRN4", "Inh_MGE", "Exc_L4-5_FOXP2"))
brain@tracks[["cg_type_tracks"]] <- brain@tracks[["cg_type_tracks"]][, Oligo := round(rowMeans(.SD, na.rm = TRUE), 2), .SDcols = c("Oligo_1", "Oligo_2", "Oligo_3")]
brain@tracks[["cg_type_tracks"]] <- brain@tracks[["cg_type_tracks"]][, c("Oligo_1", "Oligo_2", "Oligo_3") := NULL]

### mCH DMR analysis ###
# calculate mCH DMRs from cluster tracks
ch_sum_matrix <- calcSmoothedWindows(brain, type = "CH", threads = 30, step = 500, species = "human", index = "chr_ch", groupBy = "type", returnSumMatrix = T, returnPctMatrix = T)
brain@tracks[["ch_type_tracks"]] <- ch_sum_matrix[["pctm_matrix"]]
brain_ch_dmrs <- testDMR(sumMatrix = ch_sum_matrix[["sum_matrix"]], eachVsAll = TRUE, nminTotal = 20, nminGroup = 10)
brain_ch_filtered_dmrs <- filterDMR(brain_ch_dmrs, method = "bonferroni", filter = TRUE, pThreshold = 0.01, logThreshold = 1.5)
brain_ch_collapsed_dmrs <- collapseDMR(brain, brain_ch_filtered_dmrs, maxDist = 2000, minLength = 2500, reduce = T, annotate = T) 
brain_ch_collapsed_dmrs$location <- paste0(brain_ch_collapsed_dmrs$chr, "_", brain_ch_collapsed_dmrs$dmr_start, "_", brain_ch_collapsed_dmrs$dmr_end)
brain_ch_collapsed_dmrs <- brain_ch_collapsed_dmrs |> dplyr::group_by(member_id, direction) |> 
  dplyr::arrange(dmr_padj, .by_group = TRUE) |> dplyr::mutate(rank_padj = 1:n()) |>
  dplyr::arrange(desc(abs(dmr_logFC)), .by_group = TRUE) |> dplyr::mutate(rank_logFC = 1:n()) |>
  rowwise() |> dplyr::mutate(total_rank = sum(rank_padj, rank_logFC)) |> 
  dplyr::mutate(location = paste0(chr, "_", (dmr_start), "_", (dmr_end)))
# brain@results[["brain_ch_collapsed_dmrs"]] <- brain_ch_collapsed_dmrs

brain_ch_dmr_counts <- full_join(brain_ch_collapsed_dmrs |> dplyr::group_by(member_id, direction) |> dplyr::summarise(n_dmrs = n()) |> dplyr::rename("type" = "member_id"),
                                 brain@metadata |> dplyr::group_by(type) |> dplyr::summarise(mch_pct = mean(mch_pct)),
                                 by = "type")
brain_ch_dmr_counts$type <- factor(brain_ch_dmr_counts$type, levels = c("Astro", "Oligo", "OPC", "Micro", "Inh_CGE", "Inh_MGE",  "Exc_L1-3_CUX2", "Exc_L2-4_RORB", "Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L5-6_PDZRN4", "Exc_L6_TLE4"))
ggplot(brain_ch_dmr_counts, aes(x = type, y = n_dmrs, fill = type)) +  geom_col() + scale_y_log10() +
  scale_fill_manual(values = brain_type_colors) + theme_classic() + facet_wrap(vars(direction))

brain_ch_collapsed_dmrs <- brain_ch_collapsed_dmrs |> dplyr::mutate(class = ifelse(member_id %in% c("Astro", "Oligo", "OPC"), "non-neuron", "neuron"))
table(brain_ch_collapsed_dmrs$class[brain_ch_collapsed_dmrs$direction == "hyper"]) # neuron: 255782 non-neuron: 239

# make tracks
ccre <- as.data.table(rtracklayer::import("~/Downloads/encodeCcreCombined.bb"))
setnames(ccre, "seqnames", "chr")
tracks <- split(ccre, ccre$ucscLabel)
ccre_track <- ccre[,.(chr, start, end)]

# plot top dmrs
toplot <- brain_ch_collapsed_dmrs |> dplyr::filter(direction == "hyper" & dmr_length > 5000) |> dplyr::group_by(member_id) |> top_n(-1, wt = dmr_padj) |> pull(location)
heatMap(brain, regions = toplot, trackOverhang = 20000, track = "ch_type_tracks", nrow = 5)

dotM(brain, c("SATB1", "SATB2","GAD1", "GAD2", "CUX2", "RORB", "FOXP2", "LRRK1", "PDZRN4", "TLE4"), groupBy = "type", sizeMatrix = "gene_ch", colorMatrix = "gene_ch") + scale_size(range = c(1, 12))
dotM(brain, c("CTTNBP2", "APP", "GABBR2", "KIF5C", "CELF4"), groupBy = "type", sizeMatrix = "gene_ch", colorMatrix = "gene_ch") + scale_size(range = c(1, 10))
dimM(brain, genes = c("CTTNBP2", "APP", "GABBR2", "KIF5C"), matrix = "gene_ch", blend = F, nrow = 1, squish = 6, colors = c("#dbdbdb", "#cccccc", "#265A61", "#0C1E20"), reduction = "umap_100k_cg_ch")
heatMap(brain, genes = c("QKI", "TBR1", "SLC17A6", "MBP"), nrow = 2, track = "cg_type_tracks")
heatMap(brain, track = "cg_cluster_tracks", regions = brain_top_dmrs$location[brain_top_dmrs$test == 4], nrow = 2, legend = F)
histograM(brain, genes = c("QKI", "TBR1", "SLC17A6"), track = "cg_type_tracks")

###############################################################################
### PBMC DATA PROCESSING ###

# PART ONE: DONE ON SERVER
# load object and perform basic filtering 
pbmc <- createObject()
# pbmc@metadata <- readRDS("~/Desktop/mount/pbmc_metadata_tmp.RData")
# add data from cellInfo file
pbmc <- addCellInfo(pbmc, file = "/secret/path/sciMET_PBMC2_reanalysis/PBMC2.cellInfo.txt")
# add batch annotation information
pbmc <- addAnnot(pbmc, "/secret/path/sciMET_PBMC2_reanalysis/PBMC2_all.paper_types.annot", name = "type")
# filter out cells with high/low coverage or mch_pct beyond what is biologically reasonable
pbmc@metadata <- pbmc@metadata |> dplyr::filter(cov > 10000000 & cov < 40000000) 
pbmc@h5paths <- data.frame(barcode = rownames(pbmc@metadata), 
                            path = rep("/secret/path/sciMET_PBMC2_reanalysis/PBMC2.h5", length(rownames(pbmc@metadata))))

# index chromosomes and calculate methylation levels over fixed genomic windows
pbmc@index[["chr_cg"]] <- indexChr(pbmc, type = "CG", threads = 20)
pbmc@genomeMatrices[["cg_100k_score"]] <- makeWindows(pbmc, stepsize = 100000, type = "CG", metric = "score", threads = 30, index = "chr_cg", nmin = 2)

# CLUSTERING TEST - CG SCORE 100KB WINDOWS
dimEstimate(pbmc, genomeMatrices = "cg_100k_score", dims = c(40), threshold = 0.99) # 24
pbmc@reductions[["irlba_100k_score"]] <- runIrlba(pbmc, genomeMatrices = c("cg_100k_score"), dims = c(24), replaceNA = c(0))
pbmc@reductions[["irlba_100k_score"]] <- regressCovBias(pbmc, reduction = "irlba_100k_score", method = "lm")
set.seed(1234)
pbmc <- runCluster(pbmc, k = 30, method = "louvain", reduction = "irlba_100k_score", colname = "100k_cluster_id")
set.seed(1234)
pbmc@reductions[["umap_irlba_100k_score"]] <- runUmap(pbmc, neighbors = 30, dist = 0.1, method = "euclidean", reduction = "irlba_100k_score") 

# calculate methylation levels over twist regions
pbmc@genomeMatrices[["twist_min1kb_score"]] <- makeWindows(pbmc, bed = "/secret/path/sciMET_CAP/merged_probes_Twist_Methylome_V1_hg38_noalt.min1kbp.bed", type = "CG", metric = "score", threads = 30, index = "chr_cg", nmin = 2)
pbmc@genomeMatrices[["twist_min1kb_score"]] <- pbmc@genomeMatrices[["twist_min1kb_score"]][!sapply(rownames(pbmc@genomeMatrices[["twist_min1kb_score"]]), function(name) length(strsplit(name, "_")[[1]]) > 3 || grepl("chrY", name)), ]
pbmc@genomeMatrices[["twist_min1kb_score"]] <- pbmc@genomeMatrices[["twist_min1kb_score"]][, colnames(pbmc@genomeMatrices[["twist_min1kb_score"]]) %in% rownames(pbmc@metadata)] 
pbmc@genomeMatrices[["twist_min1kb_score"]] <- pbmc@genomeMatrices[["twist_min1kb_score"]][rowSums(!is.na(pbmc@genomeMatrices[["twist_min1kb_score"]])) >= 200, ]

# estimate number of necessary dimensions and perform reduction
dimEstimate(pbmc, genomeMatrices = c("twist_min1kb_score"), dims = c(40), threshold = 0.99) # 36
pbmc@reductions[["irlba_twist_score"]] <- runIrlba(pbmc, genomeMatrices = c("twist_min1kb_score"), dims = c(36), replaceNA = c(0))
pbmc@reductions[["irlba_twist_score"]] <- regressCovBias(pbmc, reduction = "irlba_twist_score", method = "lm")
set.seed(1234)
pbmc <- runCluster(pbmc, k = 30, method = "louvain", reduction = "irlba_twist_score", colname = "twist_cluster_id")
set.seed(1234)
pbmc@reductions[["umap_irlba_twist_score"]] <- runUmap(pbmc, neighbors = 30, dist = 0.1, method = "euclidean", reduction = "irlba_twist_score") 

# pbmc - cluster by dmrs
download.file("https://raw.githubusercontent.com/lrylaarsdam/amethyst/main/vignettes/pbmc_vignette/pbmc_lowconfidence_dmrs.bed", "./github_pbmc_lowconf_dmr.bed") 
pbmc@genomeMatrices[["dmr_score"]] <- makeWindows(pbmc, "CG", bed = "github_pbmc_lowconf_dmr.bed", metric = "score", threads = 20)
pbmc@genomeMatrices[["dmr_score"]] <- pbmc@genomeMatrices[["dmr_score"]][rowSums(!is.na(pbmc@genomeMatrices[["dmr_score"]])) >= 300, ]
pbmc@reductions[["irlba_dmr_score"]] <- runIrlba(pbmc, genomeMatrices = c("dmr_score"), dims = c(18), replaceNA = c(0))
pbmc@reductions[["irlba_dmr_score"]] <- regressCovBias(pbmc, reduction = "irlba_dmr_score")
set.seed(12345)
pbmc <- runCluster(pbmc, k = 30, reduction = "irlba_dmr_score", method = "louvain", colname = "dmr_cluster_id")
set.seed(12345)
pbmc@reductions[["umap_irlba_dmr_score"]] <- runUmap(pbmc, neighbors = 30, dist = 0.1, method = "euclidean", reduction = "irlba_dmr_score") 

# calculate promoter mCG levels
pbmc@ref <- makeRef("hg38") # Add reference annotation file
protein_coding <- unique(pbmc@ref |> dplyr::filter(type == "gene" & gene_type == "protein_coding" & seqid != "chrM") |> dplyr::pull(gene_name))
pbmc@genomeMatrices[["cg_promoter_pct"]] <- makeWindows(pbmc, genes = protein_coding, promoter = TRUE,  type = "CG", metric = "percent", threads = 20, index = "chr_cg", nmin = 2) 
pbmc@genomeMatrices[["cg_promoter_pct"]] <- pbmc@genomeMatrices[["cg_promoter_pct"]][rowSums(!is.na(pbmc@genomeMatrices[["cg_promoter_pct"]])) >= 100, ]

# calculate high conf DMR for annotation
download.file("https://raw.githubusercontent.com/lrylaarsdam/amethyst/main/vignettes/pbmc_vignette/pbmc_highconfidence_dmrs.bed", "./github_pbmc_highconf_dmr.bed") 
pbmc@genomeMatrices[["high_dmr_pct"]] <- makeWindows(pbmc, "CG", bed = "github_pbmc_highconf_dmr.bed", metric = "percent", threads = 30)
pbmc@genomeMatrices[["high_dmr_pct_cluster"]] <- aggregateMatrix(pbmc, matrix = "high_dmr_pct", groupBy = "dmr_cluster_id")

# assess locally; rename
download.file("https://raw.githubusercontent.com/lrylaarsdam/amethyst/main/vignettes/pbmc_vignette/pbmc_ref.RData", "~/Downloads/pbmc_ref.RData") 
pbmc_ref <- readRDS("~/Downloads/pbmc_ref.RData")
pbmc_cor <- cor(merge(pbmc_ref, pbmc@genomeMatrices[["high_dmr_pct_cluster"]], by = 0) |> tibble::column_to_rownames(var = "Row.names"), use = "pairwise.complete.obs")
pbmc_cor <- pbmc_cor[c(1:ncol(pbmc_ref)), c((ncol(pbmc_ref) + 1)):ncol(pbmc_cor)]
pheatmap(pbmc_cor)
pbmc_annotation <- tidyr::pivot_longer(as.data.frame(pbmc_cor) %>% rownames_to_column("group2"), cols = c(2:(length(colnames(pbmc_cor))+1)), names_to = "group1", values_to = "cor") %>% dplyr::group_by(group1) %>% dplyr::arrange(desc(cor), .by_group = TRUE) %>% top_n(n = 1)

pbmc@metadata[["type"]] <- dplyr::recode(pbmc@metadata[["dmr_cluster_id"]],
                                         "1" = "T_CD8+", 
                                         "2" = "B", 
                                         "3" = "T_CD4+",
                                         "4" = "Mono_2",
                                         "5" = "Mono_1", 
                                         "6" = "NK", 
                                         "7" = "T_CD4+", 
                                         "8" = "T_CD4+",
                                         "9" = "Mono_1") 

pbmc_type_colors <- c("B" =  "#EE7868",
                      "NK" =  "#f5c66e",
                      "Mono_1" = "#c43b72",
                      "Mono_2" =   "#55075D",
                      "T_CD4+" =   "#B5DCA5",
                      "T_CD8+" =  "#03514F")

# aggregate matrix by type
pbmc@genomeMatrices[["cg_promoter_pct_type"]] <- aggregateMatrix(pbmc, matrix = "cg_promoter_pct", groupBy = "type")
pbmc@genomeMatrices[["cg_promoter_z_type"]] <- as.data.frame(scale(pbmc@genomeMatrices[["cg_promoter_pct_type"]], center = TRUE, scale = TRUE))
pbmc_promoters <- inner_join(tidyr::pivot_longer(pbmc@genomeMatrices[["cg_promoter_pct_type"]] |> rownames_to_column(var = "gene"), cols = -c("gene"), values_to = "pct", names_to = "type"),
                             tidyr::pivot_longer(pbmc@genomeMatrices[["cg_promoter_z_type"]] |> rownames_to_column(var = "gene"), cols = -c("gene"), values_to = "z", names_to = "type"),
                             by = c("gene", "type"))

# calculate type DMRs
pbmc_type_500bpwindows <- calcSmoothedWindows(pbmc, type = "CG", threads = 20, step = 500, genome = "hg38", index = "chr_cg", groupBy = "type", returnSumMatrix = T, returnPctMatrix = T)
pbmc_dmrs <- testDMR(sumMatrix = pbmc_type_500bpwindows[["sum_matrix"]], eachVsAll = TRUE, nminTotal = 10, nminGroup = 10)
pbmc_filtered_dmrs <- filterDMR(pbmc_dmrs, method = "bonferroni", filter = TRUE, pThreshold = 0.01, logThreshold = 2)
pbmc_collapsed_dmrs <- collapseDMR(pbmc, pbmc_filtered_dmrs, maxDist = 2000, minLength = 1000, reduce = T, annotate = T)
pbmc_collapsed_dmrs <- pbmc_collapsed_dmrs |> dplyr::group_by(member_id, direction) |> 
  dplyr::arrange(dmr_padj, .by_group = TRUE) |> dplyr::mutate(rank_padj = 1:n()) |>
  dplyr::arrange(desc(abs(dmr_logFC)), .by_group = TRUE) |> dplyr::mutate(rank_logFC = 1:n()) |>
  rowwise() |> dplyr::mutate(total_rank = sum(rank_padj, rank_logFC), location = paste0(chr, "_", dmr_start, "_", dmr_end))

pbmc@results[["dmrs_cg_type_collapsed"]] <- pbmc_collapsed_dmrs
pbmc@tracks[["cg_type_tracks"]] <- pbmc_type_500bpwindows[["pct_matrix"]]

pbmc_top_dmrs <- pbmc_collapsed_dmrs |> 
  dplyr::group_by(member_id, direction) |> 
  dplyr::arrange(dmr_padj, .by_group = TRUE) |> dplyr::mutate(rank_padj = 1:n()) |>
  dplyr::arrange(desc(abs(dmr_logFC)), .by_group = TRUE) |> dplyr::mutate(rank_logFC = 1:n()) |>
  rowwise() |> dplyr::mutate(total_rank = sum(rank_padj, rank_logFC)) |> 
  group_by(member_id, direction) |> slice_min(n = 1, order_by = total_rank) |>
  dplyr::mutate(location = paste0(chr, "_", (dmr_start - 1000), "_", (dmr_end - 1000))) |> dplyr::arrange(direction)

# classification with genehancer (pbmc)
pbmc_dmrs <- copy(data.table(pbmc_collapsed_dmrs))
setkey(pbmc_dmrs, chr, dmr_start, dmr_end)
pbmc_overlap_result <- foverlaps(pbmc_dmrs, genehancer, type="any", nomatch=NA)
# Remove rows that are duplicated by test, direction, dmr_start, dmr_end, and class
pbmc_overlap_result <- unique(pbmc_overlap_result, by = c("member_id", "direction", "dmr_start", "dmr_end", "class"))
pbmc_overlap_result <- pbmc_overlap_result |> dplyr::group_by(member_id, direction, class) |> dplyr::summarise(n = n(), dmr_sum = sum(end - start))
pbmc_overlap_result <- pbmc_overlap_result |> dplyr::mutate(other = ifelse(is.na(class), "YES", "NO")) 
ggplot(pbmc_overlap_result, aes(y = member_id, x = n, fill = class)) + geom_col() + facet_grid(vars(direction), scales = "free_y") + scale_fill_manual(values = c("#eb4034", "#52a1a3", "#f6cb52")) + theme_classic()
pbmc_overlap_result |> dplyr::group_by(class) |>  dplyr::summarise(n = sum(n)) |>dplyr::mutate(percent = 100 * n / sum(n)) # 91.2% are not NA

# GO of top dmrs (pbmc)
# topGO of results

pbmc_resultElim <- as.list(c("T_CD4+", "T_CD8+", "NK", "Mono_1", "Mono_2", "B"))
direction <- as.list(c("hypo", "hyper"))
for (i in pbmc_resultElim) {
  for (j in direction) {
    tryCatch({
      query <- unlist(strsplit(pbmc_collapsed_dmrs$gene_names[pbmc_collapsed_dmrs$member_id == i & pbmc_collapsed_dmrs$direction == j], ", "))
      GOdata <- new("topGOdata", 
                    description = "GO Enrichment Analysis", 
                    ontology = "BP", 
                    allGenes = setNames(factor(as.integer(background %in% query), levels = c(0, 1)), background),
                    geneSel = function(x) x == 1, 
                    nodeSize = 10, 
                    annot = annFUN.org, 
                    mapping = "org.Hs.eg.db", 
                    ID = "symbol")
      fisher_res <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
      pbmc_resultElim[[i]][[j]] <- GenTable(GOdata, Fisher = fisher_res, topNodes = 500, numChar = 60)
      pbmc_resultElim[[i]][[j]] <- pbmc_resultElim[[i]][[j]] |> dplyr::filter(Fisher < 0.01 & Significant > 5) |> dplyr::mutate(fold_change = Significant/Expected, Fisher = as.numeric(Fisher))
      pbmc_resultElim[[i]][[j]] <- pbmc_resultElim[[i]][[j]] |> dplyr::filter(fold_change > 2)
      pbmc_resultElim[[i]][[j]] <- janitor::clean_names(pbmc_resultElim[[i]][[j]])
      pbmc_resultElim[[i]][[j]]$type <- as.character(i)
      pbmc_resultElim[[i]][[j]]$direction <- as.character(j)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
}
pbmc_resultElim[1:6] <- NULL
for (i in c(1:6)) {
  pbmc_resultElim[[i]] <- do.call(rbind, pbmc_resultElim[[i]])
}
pbmc_resultElim <- do.call(rbind, pbmc_resultElim)
pbmc_resultElim <- pbmc_resultElim |> group_by(type, direction) |> dplyr::mutate(rank_pval = rank(fisher, ties.method = "min"), 
                                                                                   rank_logfc = rank(-fold_change, ties.method = "min"),
                                                                                   rank_total = rank_pval + rank_logfc) 

###############################################################################
### BENCHMARKING ###

# for amethyst, running separate scripts to get /usr/bin/time metrics. The scripts are pasted below.
library(amethyst)
library(data.table)
library(dplyr)
library(tibble)
library(tidyr)
library(plyr)
library(future)
library(furrr)
library(purrr)

# sciMETv2_3842F_benchmark_CG_index.R
brain <- readRDS("/secret/path/amethyst/manuscript_analysis_v3/benchmark/amethyst/sciMETv2_3842F_index_metadata_paths.RData")
t1_amethyst <- Sys.time() # "2025-03-28 11:10:41 PDT"
brain@index[["chr_cg"]] <- indexChr(obj = brain, type = "CG", threads = 10) # not including chrY
t2_amethyst <- Sys.time() # "2025-03-28 11:13:16 PDT"
save.image("/secret/path/amethyst/manuscript_analysis_v3/benchmark/amethyst/sciMETv2_3842F_benchmark_CG_index.RData")

# sciMETv2_3842F_benchmark_CH_index.R
brain <- readRDS("/secret/path/amethyst/manuscript_analysis_v3/benchmark/amethyst/sciMETv2_3842F_index_metadata_paths.RData")
t1_amethyst <- Sys.time() # "2025-03-28 11:11:39 PDT"
brain@index[["chr_ch"]] <- indexChr(obj = brain, type = "CH", threads = 10) # not including chrY
t2_amethyst <- Sys.time() # "2025-03-28 11:56:33 PDT"
save.image("/secret/path/amethyst/manuscript_analysis_v3/benchmark/amethyst/sciMETv2_3842F_benchmark_CH_index.RData")

# sciMETv2_3842F_benchmark_CG_100k.R
brain <- readRDS("/secret/path/amethyst/manuscript_analysis_v3/benchmark/amethyst/sciMETv2_3842F_index_metadata_paths.RData")
t1_amethyst <- Sys.time() # "2025-03-28 10:37:23 PDT"
brain@genomeMatrices[["cg_100k_pct"]] <- makeWindows(obj = brain, index = "chr_cg", type = "CG", stepsize = 100000, metric = "percent", threads = 10, nmin = 2)
t2_amethyst <- Sys.time() # "2025-03-28 10:47:29 PDT"
save.image("/secret/path/amethyst/manuscript_analysis_v3/benchmark/amethyst/sciMETv2_3842F_benchmark_CG_100k.RData")

# sciMETv2_3842F_benchmark_CH_100k.R
brain <- readRDS("/secret/path/amethyst/manuscript_analysis_v3/benchmark/amethyst/sciMETv2_3842F_index_metadata_paths.RData")
t1_amethyst <- Sys.time() #"2025-03-28 10:37:29 PDT"
brain@genomeMatrices[["ch_100k_pct"]] <- makeWindows(obj = brain, index = "chr_ch", type = "CH", stepsize = 100000, metric = "percent", threads = 10, nmin = 5)
t2_amethyst <- Sys.time() # "2025-03-28 13:05:05 PDT"
save.image("/secret/path/amethyst/manuscript_analysis_v3/benchmark/amethyst/sciMETv2_3842F_benchmark_CH_100k.RData")

# sciMETv2_3842F_benchmark_CG_DMR.R
brain <- readRDS("/secret/path/amethyst/manuscript_analysis_v3/benchmark/amethyst/sciMETv2_3842F_index_metadata_paths.RData")
brain@metadata$major_type <- dplyr::case_when(
  grepl("Exc", brain@metadata$type) ~ "exc",
  grepl("Inh", brain@metadata$type) ~ "inh",
  TRUE ~ NA)
neurons <- subsetObject(brain, cells = rownames(brain@metadata |> dplyr::filter(!is.na(major_type))))
comparisons <- data.frame(
  stringsAsFactors = FALSE,
  name = c("exc_vs_inh"),
  A = c("exc"),
  B = c("inh")
)

t1_amethyst <- Sys.time() # "2025-03-31 17:05:58 PDT"
neuron_500bp_windows <- calcSmoothedWindows(obj = neurons, type = "CG", threads = 10, step = 500, smooth = 3, genome = "hg38", index = "chr_cg", futureType = "multicore", groupBy = "major_type")
t2_amethyst <- Sys.time() # "2025-03-31 17:54:46 PDT" 48.80223 mins
neuron_dmrs <- testDMR(neuron_500bp_windows[["sum_matrix"]], comparisons = comparisons, nminTotal = 10, nminGroup = 10) 
t3_amethyst <- Sys.time() # "2025-03-31 17:57:10 PDT" 
neuron_dmrs_high <- filterDMR(neuron_dmrs, method = "bonferroni", filter = TRUE, pThreshold = 0.01, logThreshold = 1.5)
neuron_dmrs_med <- filterDMR(neuron_dmrs, method = "BH", filter = TRUE, pThreshold = 0.01, logThreshold = 1.25)
t4_amethyst <- Sys.time() # "2025-03-31 17:57:15 PDT"
neurons@ref <- makeRef("hg38")
t5_amethyst <- Sys.time() # "2025-03-31 17:57:49 PDT"
neuron_collapsed_dmrs_high <- collapseDMR(neurons, neuron_dmrs_high, maxDist = 1000, minLength = 2000, reduce = T, annotate = T)
neuron_collapsed_dmrs_med <- collapseDMR(neurons, neuron_dmrs_med, maxDist = 1000, minLength = 1000, reduce = T, annotate = T)
t6_amethyst <- Sys.time() # "2025-03-31 17:57:49 PDT" # 3.05552 mins

save.image("/secret/path/amethyst/manuscript_analysis_v3/benchmark/amethyst/sciMETv2_3842F_benchmark_CG_DMR.RData")

### comparison of multiple dimensionality reduction and clustering methods - done on cluster for accurate time comparison ###
# sciMETv2_3842F_benchmark_CG_irlba.R
brain <- readRDS("sciMETv2_3842F_cg_100k.RData") 
Sys.time() # "2025-04-03 12:44:04 PDT"
brain@reductions[["irlba_cg_100k_score_benchmark"]] <- runIrlba(brain, genomeMatrices = c("cg_100k_score"), dims = c(26), replaceNA = c(0))
brain@reductions[["umap_irlba_cg_100k_score_benchmark"]] <- runUmap(brain, neighbors = 25, dist = 0.1, method = "euclidean", reduction = "irlba_cg_100k_score_benchmark") 
brain <- runCluster(brain, k = 25, reduction = "irlba_cg_100k_score_benchmark", method = "louvain", colname = "irlba_cg_louvain") 
Sys.time() # "2025-04-03 12:44:23 PDT"

# pca chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://cran.r-project.org/web/packages/LearnPCA/vignettes/Vig_03_Step_By_Step_PCA.pdf
# sciMETv2_3842F_benchmark_CG_pca.R
brain_cg_100k_score_pca <- data.frame(t(brain@genomeMatrices[["cg_100k_score"]])) # already centered and scaled
brain_cg_100k_score_pca[is.na(brain_cg_100k_score_pca)] <- 0

t1_amethyst <- Sys.time() #"2025-04-16 12:38:16 PDT"
brain_cg_100k_score_pca <- stats::prcomp(brain_cg_100k_score_pca, center = F, scale. = F)
brain@reductions[["pca_cg_100k_score_benchmark"]] <- brain_cg_100k_score_pca$x[, 1:26]
t2_amethyst <- Sys.time() # "2025-04-16 12:41:00 PDT"
brain@reductions[["umap_pca_cg_100k_score_benchmark"]] <- runUmap(brain, neighbors = 25, dist = 0.1, method = "euclidean", reduction = "pca_cg_100k_score_benchmark") 
brain <- runCluster(brain, k = 25, reduction = "pca_cg_100k_score_benchmark", method = "louvain", colname = "pca_cg_louvain") 

# testing updated clustering
brain <- runCluster(brain, k = 25, reduction = "irlba", method = "louvain", colname = "louvain") 
brain <- runCluster(brain, k = 25, reduction = "irlba", method = "leiden", colname = "leiden")
mclust::adjustedRandIndex(brain@metadata$louvain, brain@metadata$leiden)

p1 <- dimFeature(brain, colorBy = louvain, reduction = "umap", colors = makePalette(14, 14)) + ggtitle("louvain cluster assignment")
p2 <- dimFeature(brain, colorBy = leiden, reduction = "umap", colors = makePalette(14, 11)) + ggtitle("leiden cluster assignment")
plot_grid(p1, p2)

###############################################################################
### MOFA+ ### https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/getting_started_R.html

# cd /secret/path/amethyst/manuscript_analysis_v3/benchmark/MOFA
# load("/secret/path/amethyst/manuscript_analysis_v1-2/manuscript_workspace_240712.RData")

library("MOFA2")

data <- brain@genomeMatrices[c("cg_100k_score")]
data <- lapply(names(data), function(name) {
  x <- data[[name]] |> 
    tibble::rownames_to_column(var = "feature") |> 
    tidyr::pivot_longer(cols = -"feature", names_to = "sample", values_to = "value") |>
    dplyr::filter(!is.na(value)) |>
    dplyr::mutate(view = name) 
  return(x)
})
data <- data.table::rbindlist(data)
MOFAobject <- MOFA2::create_mofa(data)

# plot_data_overview(MOFAobject) # already filtered
# saveRDS(MOFAobject, "sciMETv2_3842F_mofa_object.RData")
# the following is in rscript sciMETv2_3842F_run_mofa.R; only MOFAobject is loaded

# set parameters
data_opts <- get_default_data_options(MOFAobject) # using defaults - data is already scaled
model_opts <- get_default_model_options(MOFAobject)
train_opts <- get_default_training_options(MOFAobject)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path("/secret/path/amethyst/manuscript_analysis_v3/benchmark/MOFA/brain_cg100kscore_mofa_model_allfeatures.hdf5")

t1_mofa_allfeatures <- Sys.time() # "2025-03-27 10:51:04 PDT"
set.seed(1234)
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)
t2_mofa_allfeatures <- Sys.time() # "2025-03-27 10:54:21 PDT"
total_time <- t2_mofa_allfeatures - t1_mofa_allfeatures # Time difference of 3.290394 mins

# save.image("sciMETv2_3842F_run_mofa_result.RData")
# downstream analysis: https://www.bioconductor.org/packages/release/bioc/vignettes/MOFA2/inst/doc/downstream_analysis.html

model_allfeatures <- load_model("~/Desktop/mount/MOFA/brain_cg100kscore_mofa_model_allfeatures.hdf5")
samples_metadata(model_allfeatures) <- inner_join(samples_metadata(model_allfeatures), 
                                      brain@metadata |> tibble::rownames_to_column(var = "sample"), by = "sample")
plot_data_scatter(model_allfeatures, view = "cg_100k_score", factor = 1, features = 5, add_lm = TRUE, color_by = "type")

set.seed(1234)
t1 <- Sys.time()
model_allfeatures <- run_umap(model_allfeatures, n_neighbors = 25, min_dist = 0.1, metric = "euclidean")
t2 <- Sys.time() # Time difference of 4.33674 secs
plot_dimred(model_allfeatures, method = "UMAP", color_by = "type", dot_size = 1) + ggtitle("MOFA+; all features")

#####################
### MOFA using variable features as suggested ###
# identify variable features
brain@genomeMatrices[["cg_100k_score_type"]] <- aggregateMatrix(brain, matrix = "cg_100k_score", groupBy = "type")
keep <- rowSums(!is.na(brain@genomeMatrices[["cg_100k_score"]])) >= 1300 # select for windows where large majority of cells have values
var_100k_cg_windows <- brain@genomeMatrices[["cg_100k_score_type"]][keep, ] |> 
  tibble::rownames_to_column(var = "window") |> 
  tidyr::pivot_longer(cols = -c("window"), values_to = "score_m", names_to = "type") |>
  dplyr::filter(!is.na(score_m)) |>
  dplyr::group_by(window) |>
  dplyr::filter(n() == 12) |> 
  dplyr::summarise(cg_mean = mean(score_m, na.rm = T), cg_sd = sd(score_m, na.rm = T)) |> 
  dplyr::arrange(desc(cg_sd))
var_100k_cg_windows <- var_100kwindows$window[1:3000]

data <- brain@genomeMatrices[c("cg_100k_score")]
data[["cg_100k_score"]] <- data[["cg_100k_score"]][var_100k_cg_windows, ]

data <- lapply(names(data), function(name) {
  x <- data[[name]] |> 
    tibble::rownames_to_column(var = "feature") |> 
    tidyr::pivot_longer(cols = -"feature", names_to = "sample", values_to = "value") |>
    dplyr::filter(!is.na(value)) |>
    dplyr::mutate(view = name) 
  return(x)
})
data <- data.table::rbindlist(data)
MOFAobject <- MOFA2::create_mofa(data)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path("/secret/path/amethyst/manuscript_analysis_v3/benchmark/MOFA/brain_cg100kscore_mofa_model_varfeatures.hdf5")

t3_mofa_varfeatures <- Sys.time() # "2025-03-24 10:09:59 PDT"
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)
t4_mofa_varfeatures <- Sys.time() # "2025-03-24 10:10:19 PDT"
total_time <- t4_mofa_varfeatures - t3_mofa_varfeatures # Time difference of 19.88988 secs

model_varfeatures <- load_model("~/Desktop/mount/MOFA/brain_cg100kscore_mofa_model_varfeatures.hdf5")
samples_metadata(model_varfeatures) <- inner_join(samples_metadata(model_varfeatures), 
                                                  brain@metadata |> tibble::rownames_to_column(var = "sample"), by = "sample")
plot_data_scatter(model_varfeatures, view = "cg_100k_score", factor = 1, features = 5, add_lm = TRUE, color_by = "type")

set.seed(4321)
model_varfeatures <- run_umap(model_varfeatures, n_neighbors = 25, min_dist = 0.1, metric = "euclidean")
plot_dimred(model_varfeatures, method = "UMAP", color_by = "type", dot_size = 1) + ggtitle("MOFA+; 3k var 100kb CG windows")

brain@reductions[["umap_mofa_100k_vmrs"]] <- data.frame(row.names = model_varfeatures@samples_metadata[["sample"]], 
                                                        dim_x = model_varfeatures@dim_red[["UMAP"]][["UMAP1"]],
                                                        dim_y = model_varfeatures@dim_red[["UMAP"]][["UMAP2"]])

#####################
### MOFA using DMRs for pbmc ###
pbmc_data <- pbmc@genomeMatrices[c("dmr_score")]
pbmc_data <- lapply(names(pbmc_data), function(name) {
  x <- pbmc_data[[name]] |> 
    tibble::rownames_to_column(var = "feature") |> 
    tidyr::pivot_longer(cols = -"feature", names_to = "sample", values_to = "value") |>
    dplyr::filter(!is.na(value)) |>
    dplyr::mutate(view = name) 
  return(x)
})
pbmc_data <- data.table::rbindlist(pbmc_data)
pbmc_MOFAobject <- MOFA2::create_mofa(pbmc_data)

# set parameters
data_opts <- get_default_data_options(pbmc_MOFAobject) # using defaults - data is already scaled
model_opts <- get_default_model_options(pbmc_MOFAobject)
train_opts <- get_default_training_options(pbmc_MOFAobject)

pbmc_MOFAobject <- prepare_mofa(
  object = pbmc_MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

outfile = file.path("/secret/path/amethyst/manuscript_analysis_v3/benchmark/MOFA/pbmc_cg_dmr_score_mofa_model_allfeatures.hdf5")
set.seed(1234)
t1_mofa <- Sys.time() # "2025-04-08 15:38:24 PDT"
pbmc_MOFAobject.trained <- run_mofa(pbmc_MOFAobject, outfile, use_basilisk = TRUE)
t2_mofa <- Sys.time() # "2025-04-08 15:57:38 PDT"

pbmc_model_varfeatures <- run_umap(pbmc_MOFAobject.trained, n_neighbors = 30, min_dist = 0.1, metric = "euclidean")
pbmc@reductions[["umap_mofa_vmrs"]] <- data.frame(row.names = pbmc_model_varfeatures@samples_metadata[["sample"]], 
                                                  dim_x = pbmc_model_varfeatures@dim_red[["UMAP"]][["UMAP1"]],
                                                  dim_y = pbmc_model_varfeatures@dim_red[["UMAP"]][["UMAP2"]])

###############################################################################
### scAI ### https://github.com/sqjin/scAI
# https://htmlpreview.github.io/?https://github.com/sqjin/scAI/blob/master/examples/walkthrough_mESC_dataset.html

# Requires paired transcriptomic and epigenomic data in same cell

###############################################################################
### methscan ### https://anders-biostat.github.io/MethSCAn/tutorial.html
# /usr/bin/time -v apptainer exec methscan.sif methscan prepare sciMETv2_3842F_files/*.cov sciMETv2_3842F_compact &> methscan_prepare_stats.log &
# /usr/bin/time -v apptainer exec methscan.sif methscan smooth sciMETv2_3842F_compact &> methscan_smooth_stats.log & # already filtered
# /usr/bin/time -v apptainer exec methscan.sif methscan scan --threads 10 sciMETv2_3842F_compact sciMETv2_3842F_VMRs.bed &> methscan_smooth_stats.log &
# /usr/bin/time -v apptainer exec methscan.sif methscan matrix --threads 10 sciMETv2_3842F_VMRs.bed sciMETv2_3842F_compact sciMETv2_3842F_VMR_matrix &> methscan_matrix_stats.log &

### brain ###  
# cg_sciMETv2_3842F_cluster.R
meth_mtx <- read.csv("/secret/path/amethyst/manuscript_analysis_v3/benchmark/methscan/cg_sciMETv2_3842F_VMR_matrix/mean_shrunken_residuals.csv.gz", row.names=1) |> as.matrix()

prcomp_iterative <- function(x, n=10, n_iter=50, min_gain=0.001, ...) {
  mse <- rep(NA, n_iter)
  na_loc <- is.na(x)
  x[na_loc] = 0  # zero is our first guess
  
  for (i in 1:n_iter) {
    prev_imp <- x[na_loc]  # what we imputed in the previous round
    # PCA on the imputed matrix
    pr <- prcomp_irlba(x, center = F, scale. = F, n = n, ...)
    # impute missing values with PCA
    new_imp <- (pr$x %*% t(pr$rotation))[na_loc]
    x[na_loc] <- new_imp
    # compare our new imputed values to the ones from the previous round
    mse[i] = mean((prev_imp - new_imp) ^ 2)
    # if the values didn't change a lot, terminate the iteration
    gain <- mse[i] / max(mse, na.rm = T)
    if (gain < min_gain) {
      message(paste(c("\n\nTerminated after ", i, " iterations.")))
      break
    }
  }
  pr$mse_iter <- mse[1:i]
  pr
}

t1_methscan <- Sys.time() # "2025-04-17 07:58:11 PDT"
pca <- meth_mtx |> scale(center = T, scale = F) |> prcomp_iterative(n = 26)  # same as CG 100k windows
t2_methscan <- Sys.time() # "2025-04-17 08:02:37 PDT"
pca_tbl <- data.frame(pca$x) 
rownames(pca_tbl) <- gsub(pattern = ".methscan.CG", replacement = "", x = rownames(meth_mtx))
t3_methscan <- Sys.time() # "2025-04-17 08:02:37 PDT"

methscan_obj <- createObject()
t4_methscan <- Sys.time()
umap_dims <- as.data.frame(umap::umap(pca_tbl, method = "naive", dims = 2, n_components = 2, n_neighbors = 25, min_dist = .1, metric = "euclidean")$layout)
t5_methscan <- Sys.time() # Time difference of 4.164483 secs

brain@reductions[["umap_methscan_vmrs"]] <- umap_dims
colnames(brain@reductions[["umap_methscan_vmrs"]]) <- c("dim_x", "dim_y")

### pbmc ###  
pbmc_meth_mtx <- read.csv("/secret/path/amethyst/manuscript_analysis_v3/benchmark/methscan/pbmc_files_compact_VMR_matrix/mean_shrunken_residuals.csv.gz", row.names=1) |> as.matrix()
pbmc_pca <- pbmc_meth_mtx |> scale(center = T, scale = F) |> prcomp_iterative(n = 18)  # same as CG DMRs
pbmc_pca_tbl <- data.frame(pbmc_pca$x) 
rownames(pbmc_pca_tbl) <- gsub(pattern = ".methscan.CG", replacement = "", x = rownames(pbmc_meth_mtx))
pbmc_methscan_obj <- createObject()
pbmc_umap_dims <- as.data.frame(umap::umap(pbmc_pca_tbl, method = "naive", dims = 2, n_components = 2, n_neighbors = 30, min_dist = .1, metric = "euclidean")$layout)
pbmc@reductions[["umap_methscan_vmrs"]] <- pbmc_umap_dims
colnames(pbmc@reductions[["umap_methscan_vmrs"]]) <- c("dim_x", "dim_y")


###############################################################################
### scMET ### https://rpubs.com/cakapourani/scmet-analysis
# cd "/secret/path/amethyst/manuscript_analysis_v3/benchmark/scmet"
library("scMET")

### CG 100kb windows ###
## make Y data matrix
future::plan(future::multicore, workers = 10)
barcodes <- as.list(rownames(brain@h5paths))
paths <- as.list(brain@h5paths$path)
chr_groups <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
index <- "chr_cg"
type <- "CG"
stepsize <- 100000

by_chr <- list()
for (chr in chr_groups) {
  
  # get index positions
  sites <- brain@index[[index]][[chr]] # get chr index for h5 file
  
  # add up sum c and sum t in member cells
  windows <- furrr::future_pmap(.l = list(paths, barcodes), .f = function(path, barcode) {
    tryCatch({
      barcode_name <- sub("\\..*$", "", barcode)
      h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", barcode),
                                                 start = sites$start[sites$cell_id == barcode],
                                                 count = sites$count[sites$cell_id == barcode])) # read in 1 chr at a time
      h5 <- h5[pos %% stepsize == 0, pos := pos + 1] # otherwise sites exactly divisible by stepsize will be their own window
      h5 <- h5[, Feature := paste0(chr, "_", plyr::round_any(pos, stepsize, floor), "_", plyr::round_any(pos, stepsize, ceiling))]
      meth_window <- h5[, .(Cell = barcode_name, total_reads = sum(c + t, na.rm = TRUE), met_reads = sum(c, na.rm = TRUE)), by = Feature]
    }, error = function(e) {
      cat("Error processing data for barcode", barcode, ":", conditionMessage(e), "\n")
      return(NULL)  # Return NULL or any other value indicating failure
    })
  }, .progress = TRUE)
  
  # Reduce merge the data.tables per cell in chunks (way faster for large datasets)
  by_chr[[chr]] <- do.call(rbind, windows)
  cat("\nCompleted ", chr,"\n")
}

cg_windows <- rbindlist(by_chr)
cg_windows <- cg_windows[, if(.N >= 1000) .SD, by = Feature] # at least 1000/1346 of cells have the feature

### make covariate matrix
cg_sites <- read.csv("hs38d1_noalt.CG.bed", header = F, sep = '\t')
colnames(cg_sites) <- c("chr", "pos", "seq", "strand")
cg_sites <- data.table(cg_sites)
cg_sites <- cg_sites[pos %% stepsize == 0, pos := pos + 1] # otherwise sites exactly divisible by stepsize will be their own window
cg_sites <- cg_sites[, Feature := paste0(chr, "_", plyr::round_any(pos, stepsize, floor), "_", plyr::round_any(pos, stepsize, ceiling))]
cg_density <- cg_sites[, .(cg_density = sum(pos != 0)/100000), by = Feature]
X <- model.matrix(~ cg_density, data = cg_density |> tibble::column_to_rownames(var = "Feature"))
X[, "cg_density"] <- scale(X[, "cg_density"])
X <- X[unique(cg_windows$Feature), ]

# save.image("sciMETv2_3842F_benchmark_scmet_preprocessing.RData") # after removing everything besides cg_windows and X
# running script below w/ /usr/time/bin -v sciMETv2_3842F_benchmark_scmet_fit_obj.R

t1_scmet <- Sys.time() # "2025-03-21 16:34:49 PDT"
fit_obj_cg <- scmet(Y = cg_windows, X = X, L = 4, iter = 1000, seed = 12, n_cores = 10) #
t2_scmet <- Sys.time() # "2025-03-22 01:46:38 PDT"; Time difference of 9.196847 hours
# saveRDS(fit_obj_cg , "scmet_brain_fit_obj_cg.RData")
t3_scmet <- Sys.time() # "2025-03-22 19:49:17 PDT"
fit_obj_cg <- scmet_hvf(scmet_obj = fit_obj_cg, delta_e = 0.75, evidence_thresh = 0.8, efdr = 0.1)
t4_scmet <- Sys.time() # "2025-03-22 19:49:20 PDT"

# save.image("/secret/path/amethyst/manuscript_analysis_v3/benchmark/scmet/sciMETv2_3842F_benchmark_scmet_result.RData")
# end of sciMETv2_3842F_benchmark_scmet_fit_obj.R script

## visualize output locally
scmet_plot_efdr_efnr_grid(obj = fit_obj_cg, task = "hvf")
gg1 <- scmet_plot_mean_var(obj = fit_obj_cg, y = "gamma", task = NULL, show_fit = TRUE)
gg2 <- scmet_plot_mean_var(obj = fit_obj_cg, y = "epsilon", task = NULL, show_fit = TRUE)
cowplot::plot_grid(gg1, gg2, ncol = 2)

scmet <- createObject()
scmet@metadata <- brain@metadata
scmet@genomeMatrices[["cg_100k_score_scmet_hvf"]] <- brain@genomeMatrices[["cg_100k_score"]][rownames(fit_obj_cg$hvf$summary[1:3000, ]), ] # top 3k features for comparison purposes; only 9 identified as HVF
t1 <- Sys.time()
scmet@reductions[["irlba"]] <- runIrlba(scmet, genomeMatrices = c("cg_100k_score_scmet_hvf"), dims = c(26), replaceNA = c(0))
t2 <- Sys.time() # 0.04040568 mins
scmet <- runUmap(scmet, neighbors = 25, dist = 0.1, method = "euclidean", reduction = "irlba") 
t3 <- Sys.time() # 0.1538821 mins
scmet <- runCluster(scmet, k = 25, reduction = "irlba")
dimFeature(scmet, colorBy = type)

###############################################################################
### BPRMeth ### https://www.bioconductor.org/packages/release/bioc/vignettes/BPRMeth/inst/doc/BPRMeth_vignette.html
### Melissa ### https://www.bioconductor.org/packages/devel/bioc/vignettes/Melissa/inst/doc/run_melissa.html

# same authors; starting structure is the same. Seems to me that melissa adds imputation to what BPRMeth already does.

library('BPRMeth')
library('Melissa')

barcodes <- as.list(rownames(brain@h5paths))
paths <- as.list(brain@h5paths$path)
chr_groups <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
index <- "chr_cg"
type <- "CG"
stepsize <- 100000

future::plan(future::multicore, workers = 5)

by_chr <- list()
for (chr in chr_groups) {
  
  # get index positions
  sites <- brain@index[[index]][[chr]] # get chr index for h5 file
  
  # add up sum c and sum t in member cells
  windows <- furrr::future_pmap(.l = list(paths, barcodes), .f = function(path, barcode) {
    tryCatch({
      barcode_name <- sub("\\..*$", "", barcode)
      h5 <- data.table::data.table(rhdf5::h5read(path, name = paste0(type, "/", barcode),
                                                 start = sites$start[sites$cell_id == barcode],
                                                 count = sites$count[sites$cell_id == barcode])) # read in 1 chr at a time
      h5 <- h5[pos %% stepsize == 0, pos := pos + 1] # otherwise sites exactly divisible by stepsize will be their own window
      h5 <- h5[, feature := paste0(chr, "_", plyr::round_any(pos, stepsize, floor), "_", plyr::round_any(pos, stepsize, ceiling))]
      h5 <- h5[, cell := barcode_name]
    }, error = function(e) {
      cat("Error processing data for barcode", barcode, ":", conditionMessage(e), "\n")
      return(NULL)  # Return NULL or any other value indicating failure
    })
  }, .progress = TRUE)
  
  # Reduce merge the data.tables per cell in chunks (way faster for large datasets)
  by_chr[[chr]] <- do.call(rbind, windows)
  cat("\nCompleted ", chr,"\n")
}
melissa_data <- rbindlist(by_chr)

# Split by 'cell', then within each, split by 'feature'
melissa_list <- split(melissa_data, by = "cell", flatten = FALSE)
melissa_list <- lapply(melissa_list, function(dt) split(dt, by = "feature", flatten = FALSE))

# Assuming `Y` is a data.table
# Convert to a list format for BPRMeth
melissa_data <- split(melissa_data, melissa_data$cell) #(see scmet data preparation for cg_windows construction)
melissa_data <- lapply(melissa_data, function(x) {x = x |> tibble::column_to_rownames(var = "feature") |> dplyr::select(total_reads, met_reads)})
melissa_data <- lapply(melissa_data, as.matrix)

t1_bprmeth <- Sys.time()
set.seed(1234)
bpr_clusters <- cluster_profiles_vb(X = bpr_data, K = 12, model = "binomial", alpha_0 = .5, beta_0 = .1) # Unable to converge using a variety of parameters
t2_bprmeth <- Sys.time()

###############################################################################
### scMelody ### https://github.com/TQBio/scMelody/tree/main

# not able to implement "package" due to poor documentation and errors in source code.

scmelody <- brain@genomeMatrices[["cg_100k_score"]] |> 
  tibble::rownames_to_column(var = "window") |> 
  tidyr::pivot_longer(cols = -window, names_to = "cell", values_to = "met") |>
  dplyr::filter(!is.na(met))
scmelody <- split(scmelody, scmelody$cell)
scmelody <- lapply(scmelody, function(df) {
  df %>%
    tidyr::separate(window, into = c("chr", "start", "end"), sep = "_", convert = TRUE) %>%
    mutate(location = paste(start, end, sep = "_")) %>%
    select(chr, location, met)
})


cl <- makeCluster(8) # the number of CPU cores for computing
clusterExport(cl,"mat_F",envir = environment())
t1_scmelody <- Sys.time()
results <- parLapply(cl,1:length(scmelody),get_res,scmelody) #data_list:single-cell methylation profiles organized as above
t2_scmelody <- Sys.time()
stopCluster(cl)
res_cor <- do.call(cbind,results)

###############################################################################
### Epiclomal ### https://github.com/molonc/Epiclomal
# remotes::install_github("molonc/Epiclomal/REpiclomal")

# only produced 13 VMRs after days

###############################################################################
### Benchmarking results - time/memory stats ###
benchmarking$step_major_class <- factor(benchmarking$step_major_class, levels = c("preparation", "generate_features", "clustering", "dmr"))
benchmarking$package <- factor(benchmarking$package, 
                               levels = rev(c("amethyst", "amethyst_facet", "allcools", "episcanpy", "methylstar", "methscan", "scmet", "melissa", "bpr_meth", "epiclomal", "deepcpg", "scmeformer", "maple", "mofa2", "seurat_v5", "couplecoc", "liger", "scai", "scmelody")))
ggplot(benchmarking |> dplyr::filter(modality != "CH" & !(step_sub_class == "pca" & package == "amethyst")), aes(x =  log(1 + wall_clock_time_mins), y = package, fill = -sub_step_numeric)) + geom_col() + 
  facet_wrap(vars(step_major_class), nrow = 1, scales = "free_y")

ggplot(benchmarking |> dplyr::filter(modality != "CH" & !(step_sub_class == "pca" & package == "amethyst")), aes(x = sub_step_numeric, y = package, size =  wall_clock_time_mins, color = max_res_set_size_kb)) + geom_point() + 
  facet_wrap(vars(step_major_class), nrow = 1, scales = "free_y") + scale_size(breaks = c(0.1, 1, 10, 100, 1000), range = c(1, 10)) + scale_color_viridis_c()

benchmarking_summary <- benchmarking |> dplyr::filter(!(step_sub_class == "pca" & package == "amethyst")) |> dplyr::group_by(package, step_major_class, modality) |> dplyr::summarise(time_mins = sum(wall_clock_time_mins), set_size_kb = mean(max_res_set_size_kb, na.rm = T))
ggplot(benchmarking_summary |> dplyr::filter(modality != "CH"), aes(x = step_major_class, y = package, color = set_size_kb, size = time_mins)) + 
  geom_point() + scale_color_viridis() + scale_size(breaks = c(0.1, 1, 10, 100, 1000), range = c(1, 10))

ggplot(benchmarking |> dplyr::distinct(package, max_met_test), aes(x = 1, y = package, size = max_met_test)) + 
  geom_point() + scale_color_viridis() + scale_size(breaks = c(1, 10, 100, 1000, 10000, 100000, 500000), range = c(.1, 8))

###############################################################################
### Benchmarking results - DMR analysis exc vs inh neurons ###

amethyst_dmrs_high_confidence <- readRDS("~/Desktop/mount/amethyst/sciMETv2_3842F_cg_exc_vs_inh_dmrs_high_confidence.RData")
amethyst_dmrs_high_confidence <- amethyst_dmrs_high_confidence |> dplyr::mutate(location = paste0(chr, "_", dmr_start, "_", dmr_end)) |>
  dplyr::arrange(dmr_padj) |> dplyr::mutate(p_rank = 1:n()) |> dplyr::arrange(desc(abs(dmr_logFC))) |> dplyr::mutate(logFC_rank = 1:n()) |>
  dplyr::rowwise() |> dplyr::mutate(total_rank = sum(p_rank, logFC_rank))

amethyst_dmrs <- readRDS("~/Desktop/mount/amethyst/sciMETv2_3842F_cg_exc_vs_inh_dmrs_medium_confidence.RData")
amethyst_dmrs$location <- paste0(amethyst_dmrs$chr, "_", amethyst_dmrs$dmr_start, "_", amethyst_dmrs$dmr_end)
amethyst_dmrs <- amethyst_dmrs |> dplyr::filter(chr != "chrY")

methscan_dmrs <- read.csv2("~/Desktop/mount/methscan/cg_sciMETv2_3842F_exc_vs_inh_DMRs_1kb.bed", sep = '\t', header = F)
colnames(methscan_dmrs) <- c("chr", "dmr_start", "dmr_end", "t_statistic", "n_sites", "n_cells_group1", 
                             "n_cells_group2", "meth_frac_group1", "meth_frac_group2", "low_group_label", "p", "adjusted_p")
methscan_dmrs$length <- methscan_dmrs$dmr_end - methscan_dmrs$dmr_start
methscan_dmrs$location <- paste0("chr", methscan_dmrs$chr, "_", methscan_dmrs$dmr_start, "_", methscan_dmrs$dmr_end)
methscan_dmrs$chr <- paste0("chr", methscan_dmrs$chr)
methscan_dmrs$chr[methscan_dmrs$chr == "chr0"] <- "chrX"
methscan_dmrs$bonferroni_p <- p.adjust(method = "bonferroni", methscan_dmrs$p)
methscan_dmrs <- methscan_dmrs |> dplyr::filter(n_sites > 10 & n_cells_group1 > 31 & n_cells_group2 > 18)

allcools_dmrs <- read.csv2("~/Desktop/mount/allcools/dmr/allcools_sciMETv2_3842F_exc_vs_inh_dmrs.csv", sep = '\t', header = T)
allcools_dmrs$dmr_padj <- p.adjust(allcools_dmrs$dmr_pval, method = "BH") # bonferroni not possible
allcools_dmrs$length <- allcools_dmrs$dmr_end - allcools_dmrs$dmr_start
allcools_dmrs$location <- paste0(allcools_dmrs$dmr_chrom, "_", allcools_dmrs$dmr_start, "_", allcools_dmrs$dmr_end)
allcools_dmrs <- allcools_dmrs |> dplyr::filter(length < 1000000)
allcools_dmrs <- allcools_dmrs |> dplyr::rename(chr = dmr_chrom)

library(GenomicRanges)
library(ggvenn)
library(ggridges)

# Convert to GRanges
gr_list <- list(
  #amethyst_high = GRanges(seqnames = amethyst_dmrs_high_confidence$chr, ranges = IRanges(start = amethyst_dmrs_high_confidence$dmr_start, end = amethyst_dmrs_high_confidence$dmr_end)),
  amethyst_med = GRanges(seqnames = amethyst_dmrs$chr, ranges = IRanges(start = amethyst_dmrs$dmr_start, end = amethyst_dmrs$dmr_end)),
  methscan = GRanges(seqnames = methscan_dmrs$chr, ranges = IRanges(start = methscan_dmrs$dmr_start, end = methscan_dmrs$dmr_end)),
  allcools = GRanges(seqnames = allcools_dmrs$chr, ranges = IRanges(start = allcools_dmrs$dmr_start, end = allcools_dmrs$dmr_end))
)

# 1. Create a union of all GRanges
all_regions <- GenomicRanges::reduce(c(gr_list$amethyst_med, gr_list$methscan, gr_list$allcools))

### counting overlaps
# 2. For each region in the union, check overlaps with each set (as TRUE/FALSE - treating each DMR as 1 locus)
venn_df_dmr <- tibble(
  region = paste0(seqnames(all_regions), ":", start(all_regions), "-", end(all_regions)),
  Amethyst = countOverlaps(all_regions, gr_list$amethyst_med) > 0,
  MethScan = countOverlaps(all_regions, gr_list$methscan) > 0,
  AllCools = countOverlaps(all_regions, gr_list$allcools) > 0
)

### counting overlaps by bp number instead (based on reviewer comment)
all_regions_expanded <- unlist(tile(all_regions, width = 1)) # expand so each bp is treated separately
cg_sites <- read.csv("OneDrive - Oregon Health & Science University/amethyst/manuscript v4/hs38d1_noalt.CG.bed", sep = '\t', header = F, col.names = c("chr", "pos", "seq", "strand"))
cg_sites = GRanges(seqnames = cg_sites$chr,
                   ranges = IRanges(start = cg_sites$pos, end = cg_sites$pos),
                   strand = cg_sites$strand)
all_regions_expanded_filtered <- subsetByOverlaps(all_regions_expanded, cg_sites) # filter by real CG sites

venn_df_dms <- tibble(
  region = paste0(seqnames(all_regions_expanded_filtered), ":", start(all_regions_expanded_filtered), "-", end(all_regions_expanded_filtered)),
  Amethyst = countOverlaps(all_regions_expanded_filtered, gr_list$amethyst_med) > 0,
  MethScan = countOverlaps(all_regions_expanded_filtered, gr_list$methscan) > 0,
  AllCools = countOverlaps(all_regions_expanded_filtered, gr_list$allcools) > 0
)

# 3. Plot with ggvenn
ggvenn(venn_df_dmr, c("Amethyst", "MethScan", "AllCools"), fill_color = c("#487ca3","#b80245","#7558a1"))
ggvenn(venn_df_dms, c("Amethyst", "MethScan", "AllCools"), fill_color = c("#487ca3","#b80245","#7558a1"))

### GO results of the different DMR methods
allcools_dmrs_annotated <- merge(
  setDT(allcools_dmrs),
  foverlaps(
    setkey(setDT(allcools_dmrs), chr, dmr_start, dmr_end),
    setkey(setDT(rename(
      distinct(dplyr::select(dplyr::filter(obj@ref, type == "gene"), seqid, start, end, gene_name)),
      chr = "seqid", gene = "gene_name"
    )), chr, start, end),
    by.x = c("chr", "dmr_start", "dmr_end"),
    by.y = c("chr", "start", "end"),
    type = "any"
  )[, .(gene_names = paste(unique(gene), collapse = ", ")), by = .(chr, dmr_start, dmr_end)],
  by = c("chr", "dmr_start", "dmr_end"),
  all.x = TRUE
)

methscan_dmrs_annotated <- merge(
  setDT(methscan_dmrs),
  foverlaps(
    setkey(setDT(methscan_dmrs), chr, dmr_start, dmr_end),
    setkey(setDT(rename(
      distinct(dplyr::select(dplyr::filter(obj@ref, type == "gene"), seqid, start, end, gene_name)),
      chr = "seqid", gene = "gene_name"
    )), chr, start, end),
    by.x = c("chr", "dmr_start", "dmr_end"),
    by.y = c("chr", "start", "end"),
    type = "any"
  )[, .(gene_names = paste(unique(gene), collapse = ", ")), by = .(chr, dmr_start, dmr_end)],
  by = c("chr", "dmr_start", "dmr_end"),
  all.x = TRUE
)

hypo_exc_results <- list(
  amethyst_query = unlist(strsplit(amethyst_dmrs$gene_names[amethyst_dmrs$direction == "hypo"], ", ")),
  methscan_query = unlist(strsplit(methscan_dmrs_annotated$gene_names[methscan_dmrs_annotated$low_group_label == "exc"], ", ")),
  allcools_query = unlist(strsplit(allcools_dmrs_annotated$gene_names[allcools_dmrs_annotated$dmr_state_sciMETv2_3842F_excitatory == 1], ", "))
)

background <- unique(brain@ref |> dplyr::filter(type == "gene" & seqid != "chrM" & seqid != "chrY") |> dplyr::pull(gene_name)) 
package_test_resultElim <- as.list(c("amethyst_query", "methscan_query", "allcools_query"))

for (i in package_test_resultElim) {
  tryCatch({
    query <- hypo_exc_results[[i]]
    GOdata <- new("topGOdata", 
                  description = "GO Enrichment Analysis", 
                  ontology = "BP", 
                  allGenes = setNames(factor(as.integer(background %in% query), levels = c(0, 1)), background),
                  geneSel = function(x) x == 1, 
                  nodeSize = 10, 
                  annot = annFUN.org, 
                  mapping = "org.Hs.eg.db", 
                  ID = "symbol")
    fisher_res <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
    package_test_resultElim[[i]] <- GenTable(GOdata, Fisher = fisher_res, topNodes = 500, numChar = 60)
    package_test_resultElim[[i]] <- package_test_resultElim[[i]] |> dplyr::filter(Fisher < 0.01 & Significant > 5) |> dplyr::mutate(fold_change = Significant/Expected, Fisher = as.numeric(Fisher))
    package_test_resultElim[[i]] <- package_test_resultElim[[i]] |> dplyr::filter(fold_change > 2)
    package_test_resultElim[[i]] <- janitor::clean_names(package_test_resultElim[[i]])
    package_test_resultElim[[i]]$method <- as.character(i)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

package_test_resultElim[1:3] <- NULL
package_test_resultElim <- do.call(rbind, package_test_resultElim)
package_test_resultElim <- package_test_resultElim |> dplyr::group_by(method) |> dplyr::mutate(rank_pval = rank(fisher, ties.method = "min"), 
                                                                    rank_logfc = rank(-fold_change, ties.method = "min"),
                                                                    rank_total = rank_pval + rank_logfc) 
ggplot(package_test_resultElim  |> dplyr::filter(annotated < 150), aes(x = fold_change, y = -log10(fisher), color = method)) + 
  geom_point() + theme_classic() + ggrepel::geom_text_repel(aes(label = term)) + facet_wrap(vars(method), scales = "free_x")

### Plot all DMRs
all_regions_df <- rbindlist(list((amethyst_dmrs_high_confidence |> dplyr::select(chr, dmr_start, dmr_end) |> dplyr::mutate(package = "amethyst_high", ymin = 2, ymax = 3)),
                                 (amethyst_dmrs |> dplyr::select(chr, dmr_start, dmr_end) |> dplyr::mutate(package = "amethyst_med", ymin = 2, ymax = 3)),
                                 (methscan_dmrs |> dplyr::select(chr, dmr_start, dmr_end) |> dplyr::mutate(package = "methscan", ymin = 1, ymax = 2)),
                                 (allcools_dmrs |> dplyr::select(chr, dmr_start, dmr_end) |> dplyr::mutate(package = "allcools", ymin = 0, ymax = 1))))
all_regions_df$chr <- factor(all_regions_df$chr, levels = c(paste0("chr", 1:22), "chrX"))
dmr_tracks <- list(amethyst_dmrs |> dplyr::select(chr, dmr_start, dmr_end) |> dplyr::rename(start = dmr_start, end = dmr_end),
                   methscan_dmrs |> dplyr::select(chr, dmr_start, dmr_end) |> dplyr::rename(start = dmr_start, end = dmr_end),
                   allcools_dmrs |> dplyr::select(chr, dmr_start, dmr_end) |> dplyr::rename(start = dmr_start, end = dmr_end))
names(dmr_tracks ) <- c("amethyst", "methscan", "allcools")
all_regions_df$length <- all_regions_df$dmr_end - all_regions_df$dmr_start
all_regions_df$package <- factor(all_regions_df$package, levels = rev(c("amethyst_high", "amethyst_med", "methscan", "allcools")))

ggplot(all_regions_df |> dplyr::filter(chr == "chr11" & dmr_start > 0 & dmr_end < 2500000)) + 
  geom_rect(aes(xmin = dmr_start, xmax = dmr_end, ymin = ymin, ymax = ymax, fill = package)) + 
  facet_grid(rows = vars(chr)) + theme_classic() + scale_fill_manual(values = rev(c("#487ca3","#b80245","#7558a1")))

ggplot(all_regions_df |> dplyr::filter(package != "amethyst_high"), aes(x = length, y = package, fill = package)) + theme_classic() + geom_density_ridges(alpha = 0.8) + scale_x_continuous(trans = 'log10') + scale_fill_manual(values = c("#487ca3","#b80245","#7558a1"))

brain@tracks[["cg_type_neurons"]] <- brain@tracks[["cg_type_tracks"]][, c(1:3, 6:8, 10:14)]
heatMap(brain,  track = "cg_type_neurons", regions = c(getGeneCoords(brain@ref, "LINGO1"), getGeneCoords(brain@ref, "SLC6A1")), 
        order = c("Inh_CGE", "Inh_MGE", "Exc_L1-3_CUX2", "Exc_L2-4_RORB", "Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L5-6_PDZRN4", "Exc_L6_TLE4"),
        extraTracks = dmr_tracks, extraTrackColors = rev(c("#487ca3","#b80245","#7558a1")), remove = "ENS", trackOverhang = 50000,
        colors = c("#005eff", "#9daec9", "#cccccc", "#dbdbdb"))

# classification with genehancer (brain)
genehancer <- data.table(read.delim("~/Desktop/mount/hg38_genehancer.class.bed", header=FALSE))
setnames(genehancer, c("chr", "start", "end", "class"))
setkey(genehancer, chr, start, end)
dmrs <- copy(all_regions_df)
setkey(dmrs, chr, dmr_start, dmr_end)
overlap_result <- foverlaps(dmrs, genehancer, type="any", nomatch=NA)
# Remove rows that are duplicated by test, direction, dmr_start, dmr_end, and class
overlap_result <- unique(overlap_result, by = c("chr", "dmr_start", "dmr_end", "class", "package"))
overlap_result <- overlap_result |> dplyr::group_by(package, class) |> dplyr::summarise(n = n())
ggplot(overlap_result, aes(y = n, x = package, fill = class)) + geom_col(position = "fill") + scale_fill_manual(values = c("#eb4034", "#52a1a3", "#f6cb52")) + theme_classic()

###############################################################################
### Glia hyper-mCH analysis ###
### step one - calculate DMGs
genematrix <- as.matrix(brain@genomeMatrices[["gene_ch_extended"]])
set.seed(1111)
# To prevent over-representation bias, average to 33 cells
downsampled_expr_list <- lapply(split(rownames(brain@metadata), brain@metadata$type), function(samples) {
  if (length(samples) <= 33) { # 33 is the smallest number of groups
    return(genematrix[, samples, drop = FALSE])
  }
  group_indices <- split(samples, rep(1:33, length.out = length(samples)))
  as.data.frame(sapply(group_indices, function(idx) rowMeans(genematrix[, idx, drop = FALSE], na.rm = TRUE)))
})
genematrix <- do.call(cbind, downsampled_expr_list)

membership <- data.frame(
  column = colnames(genematrix),
  type = sub("\\..*", "", colnames(genematrix))
)
results <- furrr::future_map(.x = rownames(genematrix), .f = function(gene) {
  gene_results <- list()  # Initialize outside the loop
  for (id in unique(membership$type)) {
    members <- membership$column[membership$type == id]
    nonmembers <- membership$column[membership$type != id]
    
    member_values <- as.numeric(genematrix[gene, members])
    nonmember_values <- as.numeric(genematrix[gene, nonmembers])
    
    tryCatch(
      {
        gene_results[[id]] <- data.frame(
          "p.val" = stats::wilcox.test(
            x = member_values,
            y = nonmember_values)$p.value,
          "gene" = gene,
          "type" = id,
          mean_1 = mean(member_values, na.rm = TRUE),
          mean_2 = mean(nonmember_values, na.rm = TRUE),
          n_members = sum(!is.na(member_values)),
          n_nonmembers = sum(!is.na(nonmember_values))
        ) |> dplyr::mutate(
          logFC = log2(mean_1 / mean_2),
          direction = ifelse(mean_1 > mean_2, "hyper", "hypo")
        )
      },
      error = function(e) {
        cat("Error processing gene:", gene, "and cluster:", id, "\n")
        gene_results[[id]] <- NA
      }
    )
    cat(paste0("\nFinished ", id))
  }
  gene_results
}, .progress = TRUE)

brain_ch_dmgs <- do.call(rbind, lapply(results, function(x) do.call(rbind, x)))
brain_ch_dmgs <- brain_ch_dmgs |> dplyr::group_by(type) |> dplyr::mutate(p.adj = stats::p.adjust(p.val, method = "bonferroni")) |> dplyr::select(p.val, p.adj, everything())
# saveRDS(brain_ch_dmgs, "brain_mCH_test_20k_genes.RData")

gene_lengths <- brain@ref |> dplyr::filter(type == "gene") |> dplyr::mutate(gene_length = (end - start))
brain_ch_dmgs$gene_length <- gene_lengths$gene_length[match(brain_ch_dmgs$gene, gene_lengths$gene_name)]
brain_ch_dmgs <- brain_ch_dmgs |> dplyr::filter(p.adj < 0.05 & abs(logFC) > 1 & n_members > 20 & n_nonmembers > 200 & gene_length > 3000) |> dplyr::group_by(gene, direction) |> dplyr::mutate(n_hits = n())
brain_ch_dmgs$direction[brain_ch_dmgs$direction == "hypermethylated"] <- "hyper"
brain_ch_dmgs$direction[brain_ch_dmgs$direction == "hypomethylated"] <- "hypo"
brain_ch_dmgs <- brain_ch_dmgs |> dplyr::filter(!(direction == "hyper" & mean_1 < 2))
#brain@results[["brain_ch_dmgs"]] <- brain_ch_dmgs
glia_hyper_ch_dmgs <- brain_ch_dmgs |> dplyr::filter(direction == "hyper" & type %in% c("Oligo", "Astro", "OPC", "Micro"))
glia_hyper_ch_dmrs <- brain_ch_collapsed_dmrs |> dplyr::filter(direction == "hyper" & type %in% c("Micro", "OPC", "Astro", "Oligo"))
brain_hyper_ch_dmgs <-  brain_ch_dmgs |> dplyr::filter(direction == "hyper")
brain_hyper_ch_dmgs <- brain_hyper_ch_dmgs|> dplyr::mutate(class = ifelse(type %in% c("Astro", "Oligo", "OPC"), "non-neuron", "neuron"))

brain_ch_dm_counts <- full_join(
  full_join(brain_ch_collapsed_dmrs  |> dplyr::group_by(type, direction) |> dplyr::summarise(n_dmrs = n()),
            brain_ch_dmgs |> dplyr::group_by(type, direction) |> dplyr::summarise(n_dmgs = n(), pct_unique = (sum(n_hits == 1)*100 / n_dmgs)),
            by = c("type", "direction")),
  brain@metadata |> dplyr::group_by(type) |> dplyr::summarise(mch_pct = mean(mch_pct)), by = "type")

brain_ch_dm_counts$n_dmgs[is.na(brain_ch_dm_counts$n_dmgs)] <- 0

# bar chart showing DMR or DMG counts
ggplot(brain_ch_dm_counts, aes(x = type, y = n_dmrs, fill = type)) +  geom_col() + scale_y_log10() +
  scale_fill_manual(values = brain_type_colors) + theme_classic() + facet_wrap(vars(direction))

### finding glia hyper-mCH genes that were identified in Lister et al. 2013 ###
glia_hyper_mch_lister <- janitor::clean_names(read_excel("~/Library/CloudStorage/OneDrive-OregonHealth&ScienceUniversity/amethyst/atlas data/Lister 2013/Table S3 - hyper mCH genes in glia.xlsx", sheet = "Table S3", skip = 2)[c(1:6)])
glia_hyper_mch_lister$gene <- toupper(glia_hyper_mch_lister$gene)
validated <- intersect(glia_hyper_ch_dmgs$gene, glia_hyper_mch_lister$gene)
heatMap(brain, track = "ch_type_tracks", genes = validated, nrow = 4)

gene_avg <- data.frame(
  gene = "gene_average",
  cell_id = colnames(brain@genomeMatrices[["gene_ch_extended"]]),
  pct = colMeans(brain@genomeMatrices[["gene_ch_extended"]], na.rm = TRUE),
  norm = colMeans(brain@genomeMatrices[["gene_ch_norm"]], na.rm = TRUE)
)
glia_hyper_mch_validated <- rbind(inner_join(tidyr::pivot_longer(brain@genomeMatrices[["gene_ch_extended"]][validated, ] |> rownames_to_column(var = "gene"), cols = -c("gene"), values_to = "pct", names_to = "cell_id"),
                                             tidyr::pivot_longer(brain@genomeMatrices[["gene_ch_norm"]][validated, ] |> rownames_to_column(var = "gene"), cols = -c("gene"), values_to = "norm", names_to = "cell_id"), by = c("gene", "cell_id")),gene_avg)
glia_hyper_mch_validated <- inner_join(glia_hyper_mch_validated, brain@metadata |> rownames_to_column(var = "cell_id") |> dplyr::select(cell_id, type), by = c("cell_id")) |> dplyr::group_by(type, gene) |> dplyr::summarise(pct = mean(pct, na.rm = T), norm = mean(norm, na.rm =T))

glia_hyper_venn_df <- tibble(
  gene = paste0(seqnames(all_regions), ":", start(all_regions), "-", end(all_regions)),
  Amethyst = countOverlaps(all_regions, gr_list$amethyst_med) > 0,
  MethScan = countOverlaps(all_regions, gr_list$methscan) > 0,
  AllCools = countOverlaps(all_regions, gr_list$allcools) > 0
)
ggvenn(list('this_study' = c(glia_hyper_ch_dmgs$gene), 
            'lister' = c(glia_hyper_mch_lister$gene)), c('this_study', 'lister'), fill_color = c("#487ca3","#b80245"))

### correlating hyper-mCH genes to RNA expression using Luo 2022 atlas ###
# relation to RNA data https://www.cell.com/cell-genomics/fulltext/S2666-979X(22)00027-1?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2666979X22000271%3Fshowall%3Dtrue
luo2022rna <- read.csv2("~/Library/CloudStorage/OneDrive-OregonHealth&ScienceUniversity/amethyst/atlas data/Luo 2022/Luo2022_34734gene_RNAexpr_by_majortype.txt", header = T, sep = '\t')                                                                    
luo2022rna <- `rownames<-`(as.data.frame(lapply(luo2022rna, as.numeric)), rownames(luo2022rna))
colnames(luo2022rna) <- gsub(pattern = "\\.", replacement = "-", x = colnames(luo2022rna))
luo2022rna <- luo2022rna[rowSums(luo2022rna[, sapply(luo2022rna, is.numeric)] == 0) < 10, ]
# there is still variation in total expression. Normalized cell count matrix, but re-normalizing group count matrix
luo2022rna <- sweep(luo2022rna, 2, colSums(luo2022rna), FUN = "/") * 10000

# look at DMGs
rna <- left_join(luo2022rna |> rownames_to_column(var = "gene"),  brain_hyper_ch_dmgs |> dplyr::ungroup() |> dplyr::select(type, gene), by = "gene") |> dplyr::rename(hyper_mch_in_type = type)
rna$hyper_mch_in_type[is.na(rna$hyper_mch_in_type)] <- "none"

# look at DMRs
# brain_ch_dmrs_genes <- brain_ch_collapsed_dmrs |> dplyr::filter(direction == "hyper") |> mutate(gene = str_split(gene_names, ",\\s*")) |> unnest(gene) |> mutate(gene = str_trim(gene)) |> dplyr::distinct(type, gene)
# rna <- left_join(luo2022rna |> rownames_to_column(var = "gene"), brain_ch_dmrs_genes |> dplyr::ungroup() |> dplyr::select(type, gene), by = "gene") |> dplyr::select(-'Micro-Endo_TYROBP') |> dplyr::rename(hyper_mch_in_type = type)
# rna$hyper_mch_in_type[is.na(rna$hyper_mch_in_type)] <- "none"

# normalize each gene to relative expression levels
# rna[, sapply(rna, is.numeric)] <- sweep(rna[, sapply(rna, is.numeric)], 1, rowMeans(rna[, sapply(rna, is.numeric)], na.rm = TRUE), FUN = "/")
# that looks weird - alternatively, make sure gene sets are comparable. (they are not - just noting that mean exp is higher between gene sets)
# ggplot(rna |> tidyr::pivot_longer(cols = -c("gene", "hyper_mch_in_type"), names_to = "type", values_to = "rna"), aes(x = hyper_mch_in_type, y = rna)) + geom_jitter(size = 0.1) + geom_violin()
# rna |> tidyr::pivot_longer(cols = -c("gene", "hyper_mch_in_type"), names_to = "type", values_to = "rna") |> dplyr::group_by(hyper_mch_in_type) |> dplyr::summarize(mean_exp = mean(rna, na.rm = T))

rna_summary <- rna |> dplyr::group_by(hyper_mch_in_type) |> dplyr::select(-gene) |> dplyr::summarise_all(mean, na.rm = T) |> column_to_rownames(var = "hyper_mch_in_type")
pheatmap(rna_summary, 
         color = colorRampPalette(rev(c("#D7191C", "#fd9061","#FDAE61", "#FFFFBF", "#ABD9E9", "#2C7BB6")))(100),
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation")

###############################################################################
### Glia hyper-mCH analysis continued - how does mCG compare in these genes? ###
### repeat the process above for mCG

brain@genomeMatrices[["gene_cg_pct"]] <- makeWindows(brain, type = "CG", genes = protein_coding, metric = "percent", index = "chr_cg", threads = 30) 
genematrix <- as.matrix(brain@genomeMatrices[["gene_cg_pct"]])
set.seed(1111)
# To prevent over-representation bias, average to 33 cells
downsampled_expr_list <- lapply(split(rownames(brain@metadata), brain@metadata$type), function(samples) {
  if (length(samples) <= 33) { # 33 is the smallest number of groups
    return(genematrix[, samples, drop = FALSE])
  }
  group_indices <- split(samples, rep(1:33, length.out = length(samples)))
  as.data.frame(sapply(group_indices, function(idx) rowMeans(genematrix[, idx, drop = FALSE], na.rm = TRUE)))
})
genematrix <- do.call(cbind, downsampled_expr_list)

membership <- data.frame(
  column = colnames(genematrix),
  type = sub("\\..*", "", colnames(genematrix))
)
results <- furrr::future_map(.x = rownames(genematrix), .f = function(gene) {
  gene_results <- list()  # Initialize outside the loop
  for (id in unique(membership$type)) {
    members <- membership$column[membership$type == id]
    nonmembers <- membership$column[membership$type != id]
    
    member_values <- as.numeric(genematrix[gene, members])
    nonmember_values <- as.numeric(genematrix[gene, nonmembers])
    
    tryCatch(
      {
        gene_results[[id]] <- data.frame(
          "p.val" = stats::wilcox.test(
            x = member_values,
            y = nonmember_values)$p.value,
          "gene" = gene,
          "type" = id,
          mean_1 = mean(member_values, na.rm = TRUE),
          mean_2 = mean(nonmember_values, na.rm = TRUE),
          n_members = sum(!is.na(member_values)),
          n_nonmembers = sum(!is.na(nonmember_values))
        ) |> dplyr::mutate(
          logFC = log2(mean_1 / mean_2),
          direction = ifelse(mean_1 > mean_2, "hyper", "hypo")
        )
      },
      error = function(e) {
        cat("Error processing gene:", gene, "and cluster:", id, "\n")
        gene_results[[id]] <- NA
      }
    )
    cat(paste0("\nFinished ", id))
  }
  gene_results
}, .progress = TRUE)

brain_cg_dmgs <- do.call(rbind, lapply(results, function(x) do.call(rbind, x)))
brain_cg_dmgs <- brain_cg_dmgs |> dplyr::group_by(type) |> dplyr::mutate(p.adj = stats::p.adjust(p.val, method = "bonferroni")) |> dplyr::select(p.val, p.adj, everything())
brain_cg_dmgs$gene_length <- gene_lengths$gene_length[match(brain_cg_dmgs$gene, gene_lengths$gene_name)]
# saveRDS(brain_cg_dmgs, "/secret/path/amethyst/manuscript_analysis_v3/brain_mCG_test_20k_genes.RData")
brain_cg_dmgs <- brain_cg_dmgs |> dplyr::filter(n_members > 20 & n_nonmembers > 200 & gene_length > 3000) |> dplyr::group_by(gene, direction) |> dplyr::mutate(n_hits = n())

brain_dmgs <- inner_join(brain_cg_dmgs, brain_ch_dmgs, by = c("gene", "type"), suffix = c(".cg", ".ch"))
brain_dmgs <- left_join(brain_dmgs, brain_hyper_ch_dmgs |> dplyr::select(gene, type) |> dplyr::rename(hyper_mch_in_type = type), by = "gene") |> dplyr::filter(type == hyper_mch_in_type)
brain_dmgs$type <- factor(brain_dmgs$type, levels = rev(brain_order))

brain_type_summary_metrics <- brain@metadata |> dplyr::group_by(type) |> dplyr::summarise(mean_mch_pct = mean(mch_pct), mean_mcg_pct = mean(mcg_pct))
brain_type_summary_metrics$type <- factor(brain_type_summary_metrics$type, levels = rev(brain_order))

# plot results - let's be honest chatGPT helped with the lm_labels function
library(broom)
lm_labels <- brain_dmgs |> dplyr::group_by(type) |>
  do({
    fit <- lm(mean_1.cg ~ mean_1.ch, data = .)
    tidy_fit <- tidy(fit)
    glance_fit <- glance(fit)
    
    intercept <- tidy_fit$estimate[1]
    slope <- tidy_fit$estimate[2]
    r2 <- glance_fit$r.squared
    p <- tidy_fit$p.value[2]
    
    tibble(
      label = paste0(
        "y = ", round(slope, 2), "x + ", round(intercept, 2),
        "\nR² = ", round(r2, 2),
        "\np = ", signif(p, 2)
      ),
      x = 10, 
      y = 0.95
    )
  })

ggplot() + 
  geom_point(data = brain_dmgs, aes(x = mean_1.ch, y = mean_1.cg, color = type), size = 0.5, alpha = 0.5) + 
  geom_smooth(data = brain_dmgs, aes(group = type, x = mean_1.ch, y = mean_1.cg, color = type), method = "lm", se = FALSE) +
  geom_vline(data = brain_type_summary_metrics |> dplyr::filter(!(type %in% c("Micro", "OPC"))), aes(xintercept = mean_mch_pct), linetype = "dashed", color = "gray40") +
  geom_hline(data = brain_type_summary_metrics |> dplyr::filter(!(type %in% c("Micro", "OPC"))), aes(yintercept = mean_mcg_pct), linetype = "dashed", color = "gray40") +
  geom_text(data = lm_labels, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 3, hjust = 0) +
  scale_color_manual(values = brain_type_colors) + theme_classic() + facet_wrap(vars(type), nrow = 5) 


###############################################################################
### ATLAS DATA ANALYSIS ###
# initial processing was done for sciMETv3 paper https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00355-0
# ultima + illumina data
atlas <- createObject() # combined CG 100k score and CG 100k pct matrices from Illumina NIH data and Ultima NIH data that were already calculated
atlas@metadata <- readRDS("~/Desktop/mount/atlas_illuminaplusultima_metadata.RData")

atlas@reductions[["irlba"]] <- runIrlba(atlas, genomeMatrices = c("cg_100k_score", "ch_100k_pct"), dims = c(10, 12))
atlas@reductions[["irlba_regressed"]] <- regressCovBias(atlas)
atlas <- runCluster(atlas, k = 50)
atlas@reductions[["umap"]] <- runUmap(atlas, neighbors = 30, dist = 0.01, method = "euclidean")

protein_coding <- unique(atlas@ref |> dplyr::filter(type == "gene" & gene_type == "protein_coding") |> dplyr::pull(gene_name))
atlas@genomeMatrices[["gene_ch"]] <- makeWindows(atlas, type = "CH", genes = protein_coding, metric = "percent", index = "chr_ch", threads = 10, save = T) 
# saved as ("/secret/path/amethyst/manuscript_analysis_v3/atlas_illuminaplusultima_gene_mch/atlas_CH_percent_20013genes_nmin2.RData")
atlas@genomeMatrices[["gene_ch_norm_cttnpb2_brox"]] <- sweep(atlas@genomeMatrices[["gene_ch"]][c("BROX", "CTTNBP2"), ], 2, colMeans(atlas@genomeMatrices[["gene_ch"]], na.rm = T), FUN="/")

atlas@metadata[["type"]] <- dplyr::recode(atlas@metadata[["cluster_id"]],
                                          "1" = "Exc_L4-6_LRRK1", 
                                          "2" = "Oligo",
                                          "3" = "Exc_L2-4_RORB_1",
                                          "4" = "Exc_L2-4_RORB_3",
                                          "5" = "Exc_L5-6_PDZRN4", 
                                          "6" = "Exc_L2-4_RORB_1",
                                          "7" = "Inh_CGE_VIP",
                                          "8" = "Exc_L6_TLE4",
                                          "9" = "Inh_CGE_VIP", 
                                          "10" = "Astro_2",
                                          "11" = "Inh_MGE_CALB1",
                                          "12" = "TBD_1", 
                                          "13" = "Oligo",
                                          "14" = "Astro_1",
                                          "15" = "Micro_1",
                                          "16" = "Exc_L4-5_FOXP2",
                                          "17" = "Exc_L2-4_RORB_1",
                                          "18" = "Exc_L2-4_RORB_2",
                                          "19" = "Inh_MGE_PVALB",
                                          "20" = "TBD_2",
                                          "21" = "Micro_2",
                                          "22" = "Astro_2",
                                          "23" = "Exc_L4-5_FOXP2",
                                          "24" = "OPC",
                                          "25" = "Exc_L2-4_RORB_3",
                                          "26" = "TBD_3",
                                          "27" = "Inh_MGE_UNC5B",
                                          "28" = "Exc_L6_TSHZ2",
                                          "29" = "Oligo",
                                          "30" = "Exc_L2-4_RORB_1",
                                          "31" = "Exc_L2-4_RORB_1")

# remove TBD2 and TBD3 - they are low cov and have weird distributions
ggplot(atlas@metadata, aes(y = type, x = cov)) + geom_jitter(size = 0.1) + geom_violin()
atlas@metadata <- atlas@metadata |> dplyr::filter(!(type %in% c("TBD_2", "TBD_3")))
atlas_order <- c("Astro_1", "Astro_2", "Oligo", "OPC", "Micro_1", "Micro_2", "TBD_1", 
                 "Inh_CGE_VIP", "Inh_MGE_PVALB", "Inh_MGE_UNC5B", "Inh_MGE_CALB1", 
                 "Exc_L2-4_RORB_1", "Exc_L2-4_RORB_2", "Exc_L2-4_RORB_3","Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L5-6_PDZRN4","Exc_L6_TSHZ2", "Exc_L6_TLE4")
atlas@metadata$type <- factor(atlas@metadata$type, levels = atlas_order)
ecker_type_colors <- c("Exc_L2-4_RORB_1" = "#e86c5f",
                       "Exc_L2-4_RORB_2" = "#6b021c",
                       "Exc_L2-4_RORB_3" = "#db5375",
                       "Exc_L4-5_FOXP2" = "#eec170",
                       "Exc_L4-6_LRRK1" = "#f58549",
                       "Exc_L6_TSHZ2" ="#ccd0b5",
                       "Exc_L6_TLE4" = "#ddc993",
                       "Exc_L5-6_PDZRN4" =  "#f2a65a",
                       "Inh_MGE_CALB1" = "#351842",
                       "Inh_CGE_VIP" =   "#bfadc7",
                       "Inh_MGE_PVALB" =  "#9db6d1",
                       "Inh_MGE_UNC5B" = "#a17f9c",
                       "TBD_1" = "#bbd7d7",
                       "TBD_2" = "#74adad",
                       "TBD_3" = "#2b5c6b",
                       "Astro_1" = "#C5D86D",
                       "Astro_2" = "#8EC379",
                       "Micro_1" =  "#57AE84",
                       "Micro_2" = "#29978B",
                       "OPC" = "#287976",                
                       "Oligo" = "#265A61")

###############################################################################
#### Atlas-scale dataset downsample tests #### 
# individual scripts ran separately for accurate memory metrics; see /Users/rylaarsd/Library/CloudStorage/OneDrive-OregonHealth&ScienceUniversity/amethyst/manuscript v1-2/downsample_test_using100kb/downsample_test.R
downsample <- read_excel("Library/CloudStorage/OneDrive-OregonHealth&ScienceUniversity/amethyst/manuscript v1-2/downsample_test_using100kb/downsample_test_outputs_using100kb.xlsx")
downsample <- janitor::clean_names(downsample) |> dplyr::filter(!is.na(context))

# Define the conversion function
convert_to_seconds <- function(time_str) {
  if (is.na(time_str) || time_str == "") return(NA)
  
  num <- as.numeric(sub(" .*", "", time_str))
  unit <- sub(".* ", "", time_str)
  
  if (unit == "secs") {
    return(num)
  } else if (unit == "mins") {
    return(num * 60)
  } else if (unit == "hours") {
    return(num * 3600)
  } else if (unit == "days") {
    return(num * 86400)
  } else {
    stop("Unknown time unit")
  }
}

# Apply the conversion function to columns 15-19
downsample <- downsample %>%
  mutate(across(15:20, ~ map_dbl(.x, convert_to_seconds)))
downsample <- downsample |> pivot_longer(cols = c(15:20), names_to = "param", values_to = "seconds") |> dplyr::select(-c(3:6, 8:14)) |> 
  dplyr::mutate(max_res_set_size_gb = maximum_resident_set_size_kbytes/1000000, cells = as.numeric(cells), norm_seconds = (seconds * 1000000000 / (mean_cov * cells))) |> dplyr::filter(!is.na(seconds))
downsample$param <- factor(downsample$param, levels = rev(c("t_indexing", "t_100kb_windows", "t_clustering", "t_calc_smoothed_windows", "t_dmr", "t_1kgenes")))

ggplot(downsample, aes(x = log10(cells), y = log10(seconds), color = param)) + geom_point(size = 2) + geom_line(linewidth = 1.5, aes(linetype = context)) +
  theme_classic() + scale_color_manual(values = rev(c("#264653", "#298B83", "#9CB478", "#edbe66", "#f28e4b", "#db4f40"))) 

downsample_summary <- downsample |> dplyr::filter(param %in% c("t_indexing", "t_100kb_windows")) |> dplyr::group_by(cells, context) |> dplyr::summarise(mem = mean(max_res_set_size_gb), cov = mean(mean_cov), norm_seconds = sum(norm_seconds))
ggplot(downsample_summary, aes(x = log10(cells), y = mem, color = context)) + geom_point(aes(size = norm_seconds)) + geom_line(linewidth = 1) + theme_classic() + scale_size(breaks = c(307, 1000, 2418), range = c(1, 7)) + 
  scale_color_manual(values = c("#373f51", "#f54e42")) + scale_y_log10()

###############################################################################
### Glia mCH / Trinucleotide analysis w/ atlas data ###
# find context-specific methylation values - did on server
nmin <- 4
future::plan(future::multicore, workers = 8)
barcodes <- as.list(rownames(atlas@h5paths))
paths <- as.list(atlas@h5paths$path)
chr_groups <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

by_chr <- list()
for (chr in chr_groups) {
  
  key <- read.table(paste0("/secret/path/amethyst/bedfiles/hg38_CHN_allgenes_", chr, ".bed"), header = F, sep = '\t')
  colnames(key) <- c("chr", "pos", "context", "sense", "gene", "gene_sense")
  key <- key |> dplyr::filter(gene %in% brain_hyper_ch_dmgs$gene) # did all of chr1 and chrX genes before switching to this
  setDT(key)
  
  # get index positions
  sites <- atlas@index[["chr_ch"]][[chr]] # get chr index for h5 file
  key_tmp <- copy(key)
  key_tmp <- key_tmp[, c("chr", "sense", "gene_sense") := NULL]
  setkey(key_tmp, pos)
  
  # add up sum c and sum t in member cells
  windows <- furrr::future_pmap(.l = list(paths, barcodes), .f = function(path, barcode) {
    tryCatch({
      barcode_name <- sub("\\..*$", "", barcode)
      result <- data.table::data.table(rhdf5::h5read(path, name = paste0("CH/", barcode),
                                                     start = sites$start[sites$cell_id == barcode],
                                                     count = sites$count[sites$cell_id == barcode])) # read in 1 chr at a time
      setkey(result, pos)
      result <- result[key_tmp, on = .(pos = pos), nomatch = 0L, .(chr, pos = x.pos, c, t, gene, context)]
      result <- result[!is.na(context)]
      result <- result[, .(cell_id =
                             round(((sum(c != 0) * 100 )/ (sum(c != 0) + sum(t != 0))), 2),
                           n = sum(c + t, na.rm = TRUE)),
                       by = .(gene, context)]
      result <- result[n >= nmin, .(cell_id, gene, context)]
      setnames(result, "cell_id", barcode_name)
    }, error = function(e) {
      cat("Error processing data for barcode", barcode, ":", conditionMessage(e), "\n")
      return(NULL)  # Return NULL or any other value indicating failure
    })
  }, .progress = TRUE)
  
  # Reduce merge the data.tables per cell in chunks (way faster for large datasets)
  windows <- split(windows, ceiling(seq_along(windows)/1000))
  windows_merged <- Reduce(function(x, y) merge(x, y, by = c("gene", "context"), all = TRUE, sort = FALSE),
                           furrr::future_map(windows, ~ Reduce(function(x, y) merge(x, y, by = c("gene", "context"), all = TRUE, sort = FALSE), .x), .progress = TRUE))
  
  saveRDS(windows_merged, paste0(chr, "_ultimaillumina_CHN.RData"))
  
  rm(key_tmp, windows, windows_merged)
  gc() 
  
  cat(paste0("\nCompleted ", unique(by_chr[[chr]][["gene"]]), " in ", chr))
}
windows_merged <- do.call(rbind, by_chr)
future::plan(future::sequential)

membership <- merge(atlas@metadata |> dplyr::select(type),
                    data.frame(mch_gene_avg = colMeans(atlas@genomeMatrices[["gene_ch"]], na.rm = T)), 
                               by = 0) |> tibble::column_to_rownames(var = "Row.names")
groups <- as.list(unique(membership$type))

chn <- list.files(pattern = "ultimaillumina_CHN.RData")
chn <- lapply(chn, function(x) {
  aggregated <- list()
  y <- readRDS(x)
  y <- as.data.frame(y)
  for (i in groups) {
    members <- membership |> dplyr::filter(type == i)
    mch_gene_avg <- mean(members$mch_gene_avg)
    members <- rownames(members)
    members <- data.frame(gene = y$gene, 
                          context = y$context, 
                          mch = rowMeans(y[, members], na.rm = TRUE), 
                          n_obs = rowSums(!is.na(y[, members])))
    members <- members |> dplyr::mutate(norm_mch = mch / mch_gene_avg, type = paste0(i))
    aggregated[[i]] <- members
  }
  return(aggregated)
})

chn <- do.call(rbind, unlist(chn, recursive = FALSE))
chn$class <- sapply(strsplit(as.character(chn$type), split = "[-_]"), function(x) x[1])
chn <- chn |> dplyr::filter(!is.na(mch) & class != "TBD" & mch < 70 & !(context %in% c("CTY", "CYA")))
# saveRDS(chn, "/secret/path/amethyst/manuscript_analysis_v3/atlas_illuminaplusultima_gene_mch/atlas_trinucleotide_data_5987hypermCHgenes_plus_chr1chrXall.RData")

##### LOCAL #####
chn$gene_length <- gene_lengths$gene_length[match(chn$gene, gene_lengths$gene_name)]
chn <- chn |> dplyr::filter(gene_length > 3000 & n_obs > 20)
chn <- left_join(chn, brain_hyper_ch_dmgs |> dplyr::select(gene, type) |> dplyr::rename(hyper_mch_in_type = type), by = "gene")
chn$hyper_mch_in_type[is.na(chn$hyper_mch_in_type)] <- "none"
chn_summary <- chn |> dplyr::group_by(context, type, hyper_mch_in_type) |> dplyr::summarise(mean_mch = mean(mch), mean_norm_mch = mean(norm_mch))
chn_summary$type <- factor(chn_summary$type, levels = atlas_order)
chn_summary$hyper_mch_in_type <- factor(chn_summary$hyper_mch_in_type, levels = c(brain_order, "none"))
chn_summary$context <- factor(chn_summary$context, levels = c("CAA", "CCA", "CTA", "CAC", "CCC", "CTC", "CAG", "CCG", "CTG", "CAT", "CCT", "CTT"))

ggplot(chn_summary, aes(y = type, x = hyper_mch_in_type, size = mean_mch, color = mean_norm_mch)) + geom_point() +  
  theme_classic()  + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(vars(context), nrow = 4) + scale_size(range = c(0, 3)) + scale_color_viridis_c(limits = c(0, 6), oob = scales::squish)

ggplot(chn_summary, aes(y = type, x = hyper_mch_in_type, fill = mean_norm_mch)) + geom_tile() +  
  theme_classic()  + theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  facet_wrap(vars(context), nrow = 4) + scale_fill_gradientn(colors = c("#440154", "#3b528b", "#21918c", "#5ec962", "#93D34E", "#C8DD39", "#fde725", "#fff6a8"))

# context prevalence
# on server: cut -f3 hg38_CHN_allgenes.bed | sort | uniq -c | sort -nr # allgenes = 20,013 protein coding genes
prevalence <- data.frame(
  stringsAsFactors = FALSE,
                         context = c("CTG","CAG","CTT","CAA","CCA","CAT",
                                     "CCT","CTC","CAC","CCC","CTA","CCG"),
                           count = c(118700539L,118700539L,117552696L,111787887L,
                                     107708924L,107512216L,103589601L,
                                     98900294L,88434749L,76256002L,75544877L,
                                     16240996L)
              )
ggplot(prevalence, aes(x = context, y = count, fill = context)) + geom_col() + scale_fill_manual(values = rev(makePalette(12, 12))) + theme_classic()  + theme(axis.text.x = element_text(angle = 90, hjust=1))

# brox and cttnbp2 only
ggplot(chn |> dplyr::filter(gene %in% c("QKI", "SATB2", "SATB1", "BROX", "PCDHGC3", "CTTNBP2")), aes(y = type, x = context, fill = norm_mch)) + geom_tile() +  
  theme_classic()  + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  facet_wrap(vars(gene), nrow = 2) + scale_size(range = c(0, 3)) + scale_fill_viridis_c(limits = c(0, 8), oob = scales::squish)

# cac and ccg focus
cacccg <- chn |> dplyr::filter(context %in% c("CAC", "CCG")) |> dplyr::filter(norm_mch < 10)
cacccg$type <- factor(cacccg$type, levels = atlas_order)
ggplot(cacccg, aes(color = type, fill = type, x = type, y = log(norm_mch + 1))) + facet_wrap(vars(context)) +  
  geom_jitter(size = 0.1, alpha = 0.3) + geom_boxplot(outlier.size = 0, outlier.alpha = 0, alpha = 0.3) +
  theme_classic() + scale_fill_manual(values = ecker_type_colors) + scale_color_manual(values = ecker_type_colors)  + theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(legend.position = "none")
ggplot(cacccg, aes(color = type, fill = type, x = type, y = log(mch + 1))) + facet_wrap(vars(context), nrow = 1) +  
  geom_boxplot(outlier.size = 0, outlier.alpha = 0, alpha = 0.7) +
  theme_classic() + scale_fill_manual(values = ecker_type_colors) + scale_color_manual(values = ecker_type_colors)  + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(legend.position = "none")
ggplot(cacccg, aes(color = type, fill = type, x = type, y = log(norm_mch + 1))) + facet_wrap(vars(context), nrow = 1) +  
  geom_boxplot(outlier.size = 0, outlier.alpha = 0, alpha = 0.7) +
  theme_classic() + scale_fill_manual(values = ecker_type_colors) + scale_color_manual(values = ecker_type_colors)  + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  theme(legend.position = "none")

# top hits per trinucleotide context
context_incidence <- read.table("~/Downloads/hg38_CHN_allgenes_context_incidence.bed") # stored in /secret/path/amethyst/bedfiles/hg38_CHN_allgenes_context_incidence.bed
colnames(context_incidence) <- c("gene", "context", "sites")
chn <- left_join(chn, context_incidence, by = c("gene", "context"))
chn$di_context <- substr(chn$context, 1, 2)

top_aggregated <- chn |> dplyr::filter(gene %in% brain_hyper_ch_dmgs$gene & !(type %in% c("OPC", "Micro_1", "Micro_2")) & !(hyper_mch_in_type == "none") & sites > 50) 
top_aggregated <- top_aggregated |> dplyr::ungroup() |> dplyr::group_by(context, type) |> top_n(n = 1, wt = mch) 
top_aggregated$type <- factor(top_aggregated$type, levels = atlas_order)

ggplot(top_aggregated, aes(y = type, x = mch, color = context, shape = context)) + geom_point(width = 0.4) + geom_text_repel(aes(label = gene)) + theme_classic() +
  scale_color_manual(values = rev(makePalette(12, 12))) + scale_shape_manual(values = c(rep(c(17, 15, 16, 18), 3)))

ggplot(top_aggregated, aes(x = norm_mch, y = mch, color = type, size = norm_mch, shape = context)) + geom_jitter(width = 0.1) + geom_text_repel(aes(label = gene)) + theme_classic() +
  scale_color_manual(values = ecker_type_colors) + scale_shape_manual(values = c(rep(c(17, 15, 16, 18), 3))) + facet_wrap(vars(di_context))

# which trinucleotide context most strongly relates to RNA
# pretty much all except CCG. maybe suggesting it is least biologically relevant, or more noisy
common_genes <- Reduce(intersect, list(chn$gene, rownames(luo2022rna), brain_hyper_ch_dmgs$gene))
luo2022rna_filtered <- luo2022rna[common_genes, ]

chn_filtered <- chn %>% dplyr::filter(gene %in% common_genes)
chn_filtered <- split(chn_filtered, chn_filtered$context)
chn_filtered <- lapply(chn_filtered, function(x) {
  x <- as.data.frame(x)
  x <- x %>%
    dplyr::group_by(gene, type) %>%
    dplyr::summarise(mch = mean(mch, na.rm = TRUE), .groups = "drop") %>%  # Or median, or first, etc.
    tidyr::pivot_wider(
      id_cols = gene,
      names_from = type,
      values_from = mch,
      names_repair = "minimal"
    )
  x <- x |> tibble::column_to_rownames(var = "gene")
  x <- x[common_genes, ]
  return(x)
})

cor_results <- lapply(chn_filtered, function(mch_mat) {
  sapply(colnames(luo2022rna_filtered), function(rna_col) {
    sapply(colnames(mch_mat), function(mch_col) {
      cor(luo2022rna_filtered[, rna_col],
          mch_mat[, mch_col],
          use = "pairwise.complete.obs",
          method = "pearson")
    })
  })
})

cor_df <- imap_dfr(cor_results, function(mat, context) {
  as.data.frame(mat) %>%
    rownames_to_column("mch_cluster") %>%
    pivot_longer(-mch_cluster, names_to = "rna_cluster", values_to = "correlation") %>%
    mutate(context = context)
})

cor_df$mch_cluster <- factor(cor_df$mch_cluster, levels = atlas_order)
cor_df$context <- factor(cor_df$context, levels = c("CAA", "CCA", "CTA", "CAC", "CCC", "CTC", "CAG", "CCG", "CTG", "CAT", "CCT", "CTT"))
ggplot(cor_df, aes(x = rna_cluster, y = mch_cluster, fill = correlation)) + 
  geom_tile() +facet_wrap(vars(context), nrow = 4) + scale_fill_viridis_c() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

mch_to_rna_key <- data.frame(
  stringsAsFactors = FALSE,
  methylation_atlas_type = c("Exc_L4-6_LRRK1",
                             "Oligo","Exc_L2-4_RORB_1",
                             "Exc_L2-4_RORB_3","Exc_L5-6_PDZRN4",
                             "Inh_CGE_VIP","Exc_L6_TLE4",
                             "Astro_2","Inh_MGE_CALB1",
                             "TBD_1","Astro_1","Micro_1",
                             "Exc_L4-5_FOXP2","Exc_L2-4_RORB_2",
                             "Inh_MGE_PVALB","Micro_2","OPC",
                             "Inh_MGE_UNC5B","Exc_L6_TSHZ2"),
  rna_atlas_type = c("L4-6_LRRK1",
                     "Oligo_MBP","L2-4_RORB",
                     "L2-4_RORB","L5-6_PDZRN4","CGE_VIP",
                     "L6_TLE4","Astro_FGF3R",
                     "MGE_CALB1",NA,"Astro_FGF3R",NA,
                     "L4-5_FOXP2","L2-4_RORB","MGE_PVALB",
                     NA,NA,"MGE_UNC5B","L6_TSHZ2"),
  brain_type = c("Exc_L4-6_LRRK1",
                 "Oligo","Exc_L2-4_RORB",
                 "Exc_L2-4_RORB","Exc_L5-6_PDZRN4",
                 "Inh_CGE","Exc_L6_TLE4",
                 "Astro","Inh_MGE",
                 NA,"Astro","Micro",
                 "Exc_L4-5_FOXP2","Exc_L2-4_RORB",
                 "Inh_MGE","Micro","OPC",
                 "Inh_MGE","Exc_L6_TSHZ2"))

expression_results <- imap_dfr(chn_filtered, ~as.data.frame(.x) %>%
                                 rownames_to_column("gene") %>%
                                 pivot_longer(-gene, names_to = "methylation_atlas_type", values_to = "mCHH") %>%
                                 mutate(context = .y)) %>%
  left_join(mch_to_rna_key, by = "methylation_atlas_type") %>%
  inner_join(luo2022rna_filtered %>%
               rownames_to_column("gene") %>%
               pivot_longer(-gene, names_to = "rna_atlas_type", values_to = "RNA"),
             by = c("gene", "rna_atlas_type"))
ggplot(expression_results |> dplyr::filter(methylation_atlas_type %in% c("Astro_1")), 
       aes(x = log(RNA), y = mCHH, color = context)) + geom_smooth()

# trinucleotide % in hyper-mCH genes per type
chn$corresponding_brain_type <- mch_to_rna_key$brain_type[match(chn$type, mch_to_rna_key$methylation_atlas_type)]
chn_hyper_mch_only <- chn |> dplyr::filter(corresponding_brain_type == hyper_mch_in_type)
chn_hyper_mch_only <- chn_hyper_mch_only |> dplyr::group_by(type, context) |> dplyr::summarise(mean_mch = mean(mch))
plot_chn_hyper_mch_only <- ggplot(chn_hyper_mch_only, aes(y = type, x = mean_mch, fill = context)) + geom_col(position = "fill") + scale_fill_manual(values = rev(makePalette(12, 12)))
chn_not_hyper_mch <- chn |> dplyr::filter(hyper_mch_in_type == "none")
chn_not_hyper_mch <- chn_not_hyper_mch |> dplyr::group_by(type, context) |> dplyr::summarise(mean_mch = mean(mch))
plot_chn_not_hyper_mch <- ggplot(chn_not_hyper_mch , aes(y = type, x = mean_mch, fill = context)) + geom_col(position = "fill") + scale_fill_manual(values = rev(makePalette(12, 12)))
plot_grid(plot_chn_hyper_mch_only, plot_chn_not_hyper_mch, nrow = 2)

###############################################################################
### Glia mCH analysis between individuals
atlas@genomeMatrices[["gene_ch"]] <- readRDS("/secret/path/amethyst/manuscript_analysis_v3/atlas_illuminaplusultima_gene_mch/atlas_CH_percent_20013genes_nmin2.RData")
atlas@metadata$type_sample <- paste0(atlas@metadata$type, "_", atlas@metadata$sample)

individuals <- c("BA46_6596", "BA46_6927", "BA46_6996", "BA46_6998")
cell_types <- unique(atlas@metadata$type)

# Create all combinations of 1 individual vs 3 others
comparison_list <- list()

for (cell in cell_types) {
  comb <- combn(individuals, 4, simplify = FALSE)
  for (group in comb) {
    combos <- combn(group, 1, simplify = FALSE)
    for (a in combos) {
      b <- setdiff(group, a)
      comparison_list[[length(comparison_list) + 1]] <- data.table(
        name = paste0(cell, "_", a),
        A = paste0(cell, "_", a),
        B = paste(paste0(cell, "_", b), collapse = ",")
      )
    }
  }
}

comparison_list <- data.frame(rbindlist(comparison_list))

genematrix <- as.matrix(atlas@genomeMatrices[["gene_ch"]])
genematrix <- genematrix[sum(!is.na(genematrix)) > 1000, ]
membership <- atlas@metadata |> dplyr::select("type_sample")
gene_results <- list()

# variant of findClusterMarkers function
library(furrr)
future::plan(future::multicore, workers = 10)

results <- future_map(1:nrow(comparison_list), function(i) {
  
  # set up - extract corresponding barcode values to group i & group name
  group_a_types <- unlist(strsplit(comparison_list[i, "A"], ','))
  group_b_types <- unlist(strsplit(comparison_list[i, "B"], ','))
  name <- comparison_list[i, "name"]
  group_a <- rownames(membership[membership$type_sample %in% group_a_types, , drop = FALSE])
  group_b <- rownames(membership[membership$type_sample %in% group_b_types, , drop = FALSE])
  if (length(group_a) < 50 || length(group_b) < 50) {
    return(NULL)
  }
  group_a_vals <- genematrix[, group_a, drop = FALSE]
  group_b_vals <- genematrix[, group_b, drop = FALSE]
  
  # For current group, apply test to each gene
  gene_results <- purrr::map_dfr(rownames(genematrix), function(gene) {
    tryCatch({
      x <- group_a_vals[gene, ]
      y <- group_b_vals[gene, ]
      
      if (sum(!is.na(x)) < 10 || sum(!is.na(y)) < 10) return(NULL)
      
      mean_1 <- mean(x, na.rm = TRUE)
      mean_2 <- mean(y, na.rm = TRUE)
      
      data.frame(
        p.val = wilcox.test(x, y)$p.value,
        gene = gene,
        comparison = name,
        mean_1 = mean_1,
        mean_2 = mean_2,
        logFC = log2(mean_1 / mean_2),
        direction = ifelse(mean_1 > mean_2, "hyper", "hypo"),
        n_group_a = sum(!is.na(x)),
        n_group_b = sum(!is.na(y))
      )
    }, error = function(e) NULL)
  })
  
  return(gene_results)
}, .progress = TRUE)

atlas_mch_gene_individual_results <- dplyr::bind_rows(results)
atlas_mch_gene_individual_results$gene_length <- gene_lengths$gene_length[match(atlas_mch_gene_individual_results$gene, gene_lengths$gene_name)]
atlas_mch_gene_individual_results$location <- atlas@ref$location[match(atlas_mch_gene_individual_results$gene, atlas@ref$gene_name)]

future::plan(future::sequential)
# saveRDS(atlas_mch_gene_individual_results, "/secret/path/amethyst/manuscript_analysis_v4/atlas_20k_gene_ch_pct_typexindividual_results.RData")

# post-process filtering
atlas_mch_gene_individual_results <- atlas_mch_gene_individual_results |> dplyr::group_by(comparison) |> dplyr::mutate(p.adj = stats::p.adjust(p.val, method = "bonferroni")) 
atlas_mch_gene_individual_results <- atlas_mch_gene_individual_results |> dplyr::mutate(significant = ifelse((p.adj < 0.05 & abs(logFC) > 1 & n_group_a > 20 & n_group_b > 20 & gene_length > 3000), "YES", "NO"))
atlas_mch_gene_individual_results$significant[atlas_mch_gene_individual_results$direction == "hyper" & atlas_mch_gene_individual_results$mean_1 < 2]  <- "NO"
atlas_mch_gene_individual_results <- atlas_mch_gene_individual_results |> dplyr::mutate(individual = str_extract(comparison, "BA46_\\d+"),type = str_remove(comparison, "_BA46_\\d+")) |> dplyr::filter(type != "TBD_1") 
atlas_mch_gene_individual_results <- tidyr::separate(data = atlas_mch_gene_individual_results, col = "location", into = c("chr", "start", "end"), sep = "_", remove = FALSE) |> dplyr::filter(significant == "YES")

atlas_mch_gene_individual_results$type <- factor(atlas_mch_gene_individual_results$type, levels = c("Exc_L2-4_RORB_1", "Exc_L2-4_RORB_2", "Exc_L2-4_RORB_3", "Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L5-6_PDZRN4", "Exc_L6_TLE4", "Exc_L6_TSHZ2", "Inh_MGE_CALB1", "Inh_MGE_PVALB", "Inh_MGE_UNC5B", "Inh_CGE_VIP", "Astro_1", "Astro_2","Micro_1","Micro_2",  "OPC", "Oligo"))
atlas_mch_gene_individual_results$chr <- factor(atlas_mch_gene_individual_results$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))
atlas_mch_gene_individual_results$sex <- factor(atlas_mch_gene_individual_results$sex, levels = c("M", "F"))

key <- data.frame(
  stringsAsFactors = FALSE,
        individual = c("BA46_6596", "BA46_6927", "BA46_6996", "BA46_6998"),
               age = c(38L, 48L, 27L, 40L),
               sex = c("F", "M", "F", "M")
)
atlas_mch_gene_individual_results <- left_join(atlas_mch_gene_individual_results, key, by = "individual")
atlas_mch_gene_individual_results <- atlas_mch_gene_individual_results |> dplyr::filter(chr != "chrY") |> dplyr::group_by(gene) |> dplyr::mutate(n_hits = n(), start = as.numeric(start), end = as.numeric(end))

# saveRDS(atlas_mch_gene_individual_results, "/secret/path/amethyst/manuscript_analysis_v4/atlas_20k_gene_ch_pct_typexindividual_results_filtered.RData")
# writexl::write_xlsx(atlas_mch_gene_individual_results, "~/OneDrive - Oregon Health & Science University/amethyst/manuscript v4/Table S7 - atlas CH DMGs between individuals.xlsx")
# atlas@results[["atlas_mch_gene_individual_results"]] <- atlas_mch_gene_individual_results

# local - plot # DMGs
ggplot(atlas_mch_gene_individual_results |> dplyr::group_by(gene, direction, individual) |>
         dplyr::summarise(n_hits = n(), mean_logFC = mean(logFC), mean_padj = pmin(-log10(mean(p.adj)), 220)), 
       aes(x = mean_logFC, y = n_hits, color = individual, size = mean_padj)) + 
  geom_point() + scale_color_manual(values = c("#143642","#EC9928","#128B8C","#A82023")) +
  theme_classic() + geom_text_repel(aes(label = gene)) 

ggplot(atlas_mch_gene_individual_results |> dplyr::mutate(p.adj = ifelse(direction == "hypo", log10(p.adj), -log10(p.adj))), 
       aes(x = p.adj, y = chr, color = type, shape = individual)) + geom_jitter() + geom_text_repel(aes(label = gene)) + 
  theme_classic() + scale_color_manual(values = ecker_type_colors) + facet_wrap(vars(sex)) + scale_shape_manual(values = c(17, 15, 16, 18))

# local - plot disribution across genome
chromosome_sizes <- data.frame(chr = c("chr1", "chr2", "chr3", "chr4", "chr5", 
                                              "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", 
                                              "chr12", "chr13", "chr14", "chr15", "chr16", 
                                              "chr17", "chr18", "chr19", "chr20", "chr21", 
                                              "chr22", "chrX", "chrY"), 
                               size = c(248956422, 242193529, 198295559, 190214555, 
                                        181538259, 170805979, 159345973, 145138636, 
                                        138394717, 133797422, 135086622, 133275309, 
                                        114364328, 107043718, 101991189, 90338345, 
                                        83257441, 80373285, 58617616, 64444167, 46709983, 
                                        50818468, 156040895, 57227415))
chromosome_sizes$chr <- factor(chromosome_sizes$chr, levels = c(paste0("chr", 1:22), "chrX", "chrY"))

ggplot() +
  geom_rect(data = chromosome_sizes, aes(xmin = 0, xmax = size, ymin = 0, ymax = 1), fill = "grey90", inherit.aes = FALSE) +
  geom_rect(data = atlas_mch_gene_individual_results |> dplyr::distinct(chr, gene, n_hits, start, end), 
            aes(xmin = start-100000, xmax = end+100000, ymin = 0, ymax = 1, color = n_hits)) + scale_color_viridis_b() +
  facet_grid(rows = vars(chr), space = "free_x", scales = "free_x") + scale_shape_manual(values = c(17, 15, 16, 18)) + theme_classic()

# local - plot USP9X, KDM6A
atlas@tracks[["type_individual_ch_pct"]] <- readRDS("~/Downloads/individualxtype_smoothed_ch_pct.RData")
atlas@tracks[["type_individual_ch_pct"]] <- atlas@tracks[["type_individual_ch_pct"]][, .SD, .SDcols = !grepl("^TBD_", colnames(atlas@tracks[["type_individual_ch_pct"]]))]

order1 <- paste0(rep(atlas_order, each = 4), "_", individuals)
order2 <- rev(paste0(rev(atlas_order), "_", rep(individuals, each = 18)))
heatMap(obj = atlas, 
        track = "type_individual_ch_pct", 
        order = order2,
        genes = "USP9X", 
        trackOverhang = 200000, colorMax = 15)

###############################################################################
### Glia mCG analysis between individuals - to check if USP9X comes up
# repeat above analysis with CG to see if USP9X comes up as significant
# cd secret/path/amethyst/manuscript_analysis_v3/atlas_illuminaplusultima_gene_mcg
atlas@genomeMatrices[["gene_cg"]] <- makeWindows(atlas, type = "CG", genes = protein_coding, metric = "percent", index = "chr_cg", threads = 30, save = T) 

# Create all combinations of 1 individual vs 3 others
comparison_list <- list()

for (cell in cell_types) {
  comb <- combn(individuals, 4, simplify = FALSE)
  for (group in comb) {
    combos <- combn(group, 1, simplify = FALSE)
    for (a in combos) {
      b <- setdiff(group, a)
      comparison_list[[length(comparison_list) + 1]] <- data.table(
        name = paste0(cell, "_", a),
        A = paste0(cell, "_", a),
        B = paste(paste0(cell, "_", b), collapse = ",")
      )
    }
  }
}

comparison_list <- data.frame(rbindlist(comparison_list))

genematrix <- as.matrix(atlas@genomeMatrices[["gene_cg"]])
genematrix <- genematrix[sum(!is.na(genematrix)) > 1000, ]
membership <- atlas@metadata |> dplyr::select("type_sample")
gene_results <- list()

# variant of findClusterMarkers function
library(furrr)
results <- future_map(1:nrow(comparison_list), function(i) {
  
  # set up - extract corresponding barcode values to group i & group name
  group_a_types <- unlist(strsplit(comparison_list[i, "A"], ','))
  group_b_types <- unlist(strsplit(comparison_list[i, "B"], ','))
  name <- comparison_list[i, "name"]
  group_a <- rownames(membership[membership$type_sample %in% group_a_types, , drop = FALSE])
  group_b <- rownames(membership[membership$type_sample %in% group_b_types, , drop = FALSE])
  if (length(group_a) < 50 || length(group_b) < 50) {
    return(NULL)
  }
  group_a_vals <- genematrix[, group_a, drop = FALSE]
  group_b_vals <- genematrix[, group_b, drop = FALSE]
  
  # For current group, apply test to each gene
  gene_results <- purrr::map_dfr(rownames(genematrix), function(gene) {
    tryCatch({
      x <- group_a_vals[gene, ]
      y <- group_b_vals[gene, ]
      
      if (sum(!is.na(x)) < 10 || sum(!is.na(y)) < 10) return(NULL)
      
      mean_1 <- mean(x, na.rm = TRUE)
      mean_2 <- mean(y, na.rm = TRUE)
      
      data.frame(
        p.val = wilcox.test(x, y)$p.value,
        gene = gene,
        comparison = name,
        mean_1 = mean_1,
        mean_2 = mean_2,
        logFC = log2(mean_1 / mean_2),
        direction = ifelse(mean_1 > mean_2, "hyper", "hypo"),
        n_group_a = sum(!is.na(x)),
        n_group_b = sum(!is.na(y))
      )
    }, error = function(e) NULL)
  })
  
  return(gene_results)
}, .progress = TRUE)

atlas_mcg_gene_individual_results <- dplyr::bind_rows(results)
atlas_mcg_gene_individual_results$gene_length <- gene_lengths$gene_length[match(atlas_mcg_gene_individual_results$gene, gene_lengths$gene_name)]
atlas_mcg_gene_individual_results$location <- atlas@ref$location[match(atlas_mcg_gene_individual_results$gene, atlas@ref$gene_name)]
# saveRDS(atlas_mcg_gene_individual_results, "/secret/path/amethyst/manuscript_analysis_v4/atlas_20k_gene_cg_pct_typexindividual_results.RData")

# post-process filtering
atlas_mcg_gene_individual_results <- atlas_mcg_gene_individual_results |> dplyr::group_by(comparison) |> dplyr::mutate(p.adj = stats::p.adjust(p.val, method = "bonferroni")) 
atlas_mcg_gene_individual_results <- atlas_mcg_gene_individual_results |> dplyr::mutate(significant = ifelse((p.adj < 0.05 & abs(logFC) > 1 & n_group_a > 20 & n_group_b > 20 & gene_length > 3000), "YES", "NO"))
atlas_mcg_gene_individual_results <- atlas_mcg_gene_individual_results |> dplyr::mutate(individual = str_extract(comparison, "BA46_\\d+"),type = str_remove(comparison, "_BA46_\\d+")) |> dplyr::filter(type != "TBD_1") 
atlas_mcg_gene_individual_results <- tidyr::separate(data = atlas_mcg_gene_individual_results, col = "location", into = c("chr", "start", "end"), sep = "_", remove = FALSE) |> dplyr::filter(significant == "YES")
# saveRDS(atlas_mcg_gene_individual_results, "/secret/path/amethyst/manuscript_analysis_v4/atlas_20k_gene_cg_pct_typexindividual_results_filtered.RData")
# note: 114 results but USP9X not one of them

# check corresponding USP9X CG patterns - local
atlas@tracks[["type_individual_cg_pct"]] <- readRDS("~/Downloads/individualxtype_smoothed_cg_pct.RData")
atlas@tracks[["type_individual_cg_pct"]] <- atlas@tracks[["type_individual_cg_pct"]][, .SD, .SDcols = !grepl("^TBD_", colnames(atlas@tracks[["type_individual_cg_pct"]]))]

order1 <- paste0(rep(atlas_order, each = 4), "_", individuals)
order2 <- rev(paste0(rev(atlas_order), "_", rep(individuals, each = 18)))
heatMap(obj = atlas, 
        track = "type_individual_cg_pct", 
        order = order2,
        genes = "USP9X", 
        trackOverhang = 200000,
        colors = c("#005eff", "#9daec9", "#cccccc", "#dbdbdb", "grey90"))

# check how many inter-individual hyper-mCH chrX genes have been previously reported to escape x-inactivation?
# https://pubmed.ncbi.nlm.nih.gov/23828890/ supplementary table 4
confirmed_x_inactivation <- read_excel("~/Downloads/1237905s_tables_new.xlsx", sheet = 4, skip = 2) |> dplyr::select(Gene, Note) |> janitor::clean_names()
atlas_mch_individual_results_unique_chrx <- atlas_mch_gene_individual_results |> dplyr::filter(chr == "chrX") |> dplyr::distinct(gene)
atlas_mch_individual_results_unique_chrx <- left_join(atlas_mch_individual_results_unique_chrx, confirmed_x_inactivation, by = "gene")
atlas_mch_individual_results_unique_chrx |> dplyr::group_by(note) |> dplyr::summarise(n = n())

# Confirmed escapee                                  20
# Escapee genes with female hypo-mCG and hypo-mCH     3
# Inactivated genes                                   3
# Novel predicted escapee                             1
# NA                                                 30

27*100/57 # 47.36842

###############################################################################
### CODE TO GENERATE FIGURE 1-5; S1-S4 ###

###############################################################################
### FIGURE 1 - Amethyst overview + benchmarking
# Figure 1a: conceptual illustration

# Fig 1b: max test dataset size
ggplot(benchmarking |> dplyr::distinct(package, max_met_test), aes(x = 1, y = package, size = max_met_test)) + 
  geom_point() + scale_color_viridis() + scale_size(breaks = c(50, 500, 5000, 50000, 500000), range = c(.1, 10)) + theme_classic()

# Fig 1c: time comparison
ggplot(benchmarking |> dplyr::filter(modality != "CH" & !(step_sub_class == "pca" & package == "amethyst")), 
       aes(x =  wall_clock_time_mins, y = package, fill = max_res_set_size_kb)) + geom_col() + theme_classic() + xlim(0,3000) +
  facet_wrap(vars(step_major_class), nrow = 1) + scale_fill_gradientn(colors = c("#d7cfe3", "#beaed6", "#937db7", "#695091", "#452d6b"), breaks = c(0, 4000000, 8000000, 12000000, 16000000))

###############################################################################
### FIGURE 2 - Dataset introduction / cell type annotation ###

brain_order <- c("Astro", "Oligo", "OPC", "Micro", "Inh_CGE", "Inh_MGE", "Exc_L1-3_CUX2", "Exc_L2-4_RORB", "Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L5-6_PDZRN4", "Exc_L6_TLE4")
pbmc_order <- c("T_CD4+", "T_CD8+", "NK", "Mono_1", "Mono_2", "B")

# Fig 2a - pbmc clustering
p1 <- dimFeature(pbmc, colorBy = type, colors = pbmc_type_colors, pointSize = 0.1, reduction = "umap_irlba_100k_score")
p2 <- dimFeature(pbmc, colorBy = type, colors = pbmc_type_colors, pointSize = 0.1, reduction = "umap_methscan_vmrs")
p3 <- dimFeature(pbmc, colorBy = type, colors = pbmc_type_colors, pointSize = 0.1, reduction = "umap_mofa_vmrs")
p4 <- dimFeature(pbmc, colorBy = type, colors = pbmc_type_colors, pointSize = 0.1, reduction = "umap_irlba_dmr_score")
plot_grid(p1, p2, p3, p4, nrow = 1)

# Fig 2b - brain clustering
p1 <- dimFeature(brain, colorBy = type, colors = brain_type_colors, pointSize = 0.1, reduction = "umap_irlba_cg_100k") 
#p2 <- dimFeature(brain, colorBy = type, colors = brain_type_colors, pointSize = 0.1, reduction = "umap_irlba_ch_100k") 
p2 <- dimFeature(brain, colorBy = type, colors = brain_type_colors, pointSize = 0.1, reduction = "umap_methscan_vmrs") 
p3 <- dimFeature(brain, colorBy = type, colors = brain_type_colors, pointSize = 0.1, reduction = "umap_mofa_100k_vmrs") 
p4 <- dimFeature(brain, colorBy = type, colors = brain_type_colors, pointSize = 0.1, reduction = "umap_100k_cg_ch") 
plot_grid(p1, p2, p3, p4, nrow = 1)

# silhouette score analysis - checked that rownames of matrix and metadata are in same order
mean(cluster::silhouette(dist = dist(brain@reductions[["umap_irlba_cg_100k"]]), x = as.integer(brain@metadata$type))[, 'sil_width']) # 0.5695381
mean(cluster::silhouette(dist = dist(brain@reductions[["umap_irlba_ch_100k"]]), x = as.integer(brain@metadata$type))[, 'sil_width']) # 0.4974
mean(cluster::silhouette(dist = dist(brain@reductions[["umap_100k_cg_ch"]]), x = as.integer(brain@metadata$cluster_id))[, 'sil_width']) # 0.6694814

# Fig 2c - pbmc marker genes
heatMap(pbmc, genes = c("CD3G", "CD8A", "KIR2DL4", "SPI1", "MPO", "CD19"), nrow = 6, order = rev(pbmc_order),
        track = "cg_type_tracks", trackOverhang = 2000, arrowOverhang = 1000 , colors = c("#005eff", "#9daec9", "#cccccc", "#dbdbdb"))
histograM(pbmc, genes = "CD2", order = pbmc_order, baseline = 0,
        track = "cg_type_tracks", trackOverhang = 2000, arrowOverhang = 1000 , colors = c("#005eff", "#9daec9", "#cccccc", "#dbdbdb"), trackScale = 10, arrowScale = 0.05)

# Fig 2d - brain marker genes (cg)
heatMap(brain, genes = c("SLC17A7", "DLX1", "C1QA", "MAG", "GFAP"), nrow = 5, track = "cg_type_tracks", order = brain_order,
        trackOverhang = 2000, arrowOverhang = 1000, colors = c("#005eff", "#9daec9", "#cccccc", "#dbdbdb"))
histograM(brain, genes = "TBR1", track = "cg_type_tracks", trackOverhang = 2000, arrowOverhang = 1000, baseline = 0, order = rev(brain_order),
          colors = c("#005eff", "#9daec9", "#cccccc", "#dbdbdb"), arrowScale = 0.025)

# Fig 2e - brain marker genes (ch)
histograM(brain, genes = c("GAD1", "SATB2"), track = "ch_type_tracks", trackOverhang = 2000, arrowOverhang = 1000, baseline = 0, order = rev(brain_order),
          colors = brain_type_colors, arrowScale = 0.025, nrow = 1) # changed function code to fill by group

# Fig 2f - pbmc marker genes. since writing this code, dotM has been modified to accept two matrices.
genes <- c("CD2", "CD3D", "CD3G", "GZMK", "KIR2DL4", "MPO", "SPI1", "MPEG1", "CD19", "CD79A")
ggplot(pbmc_promoters |> dplyr::filter(gene %in% genes)) + 
  geom_point(aes(x = factor(type, levels = pbmc_order), y = factor(gene, levels = rev(genes)), size = pct, color = z)) + scale_size_continuous(range = c(1, 12)) + theme_classic() +
  scale_color_gradientn(colors = c("#005eff", "grey60", "#ffa600")) 

# Fig 2g - brain marker genes (mCH). since writing this code, dotM has been modified to accept two matrices.
genes <- c("SATB1", "SATB2" ,"GAD1", "GAD2", "CUX2", "RORB", "FOXP2", "LRRK1", "PDZRN4", "TLE4")
brain@genomeMatrices[["gene_ch_norm"]] <- sweep(brain@genomeMatrices[["gene_ch_extended"]], 2, colMeans(brain@genomeMatrices[["gene_ch_extended"]], na.rm = T), FUN="/") # actually global %mCH in figure; changed to this for glia analysis
brain_mch <- inner_join(tidyr::pivot_longer(aggregateMatrix(brain, matrix = "gene_ch", groupBy = "type") |> rownames_to_column(var = "gene"), cols = -c("gene"), values_to = "pct", names_to = "type"),
                        tidyr::pivot_longer(aggregateMatrix(brain, matrix = "gene_ch_norm", groupBy = "type") |> rownames_to_column(var = "gene"), cols = -c("gene"), values_to = "norm", names_to = "type"),
                        by = c("gene", "type"))
ggplot(brain_mch |> dplyr::filter(gene %in% genes)) + 
  geom_point(aes(y = factor(type, levels = brain_order), x = factor(gene, levels = genes), size = pct, color = norm)) + scale_size_continuous(range = c(1, 12)) + theme_classic() +
  scale_color_gradientn(colors = c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725")) + scale_size(range = c(1, 12))

# Fig 2h - brain mch across SATB2 and GAD1
dimM(brain, genes = c("SATB2", "GAD1"), matrix = "gene_ch_imputed", blend = T, reduction = "umap_100k_cg_ch", pointSize = 0.4)

# statistics
test <- merge(brain@metadata, as.data.frame(t(brain@genomeMatrices[["gene_ch_extended"]][c("SATB2", "GAD1"), ])), by = 0)
shapiro.test(test$GAD1) # W = 0.77476, p-value < 2.2e-16
wilcox.test(test$GAD1[test$major_type == "exc"], test$GAD1[test$major_type != "exc"]) # W = 54280, p-value < 2.2e-16
test |> dplyr::mutate(group = case_when(major_type == "exc" ~ "exc", TRUE ~ "non-exc")) |> group_by(group) |> dplyr::summarise(
  mean_GAD1 = mean(GAD1, na.rm = TRUE), n = n(),
  se_GAD1 = sd(GAD1, na.rm = TRUE) / sqrt(sum(!is.na(GAD1))))

shapiro.test(test$SATB2) # W = 0.71413, p-value < 2.2e-16
wilcox.test(test$SATB2[test$major_type == "inh"], test$SATB2[test$major_type != "inh"]) # W = 57048, p-value < 2.2e-16
test |> dplyr::mutate(group = case_when(major_type == "inh" ~ "inh", TRUE ~ "non-inh")) |> group_by(group) |> dplyr::summarise(
  mean_SATB2 = mean(SATB2, na.rm = TRUE), n = n(),
  se_SATB2 = sd(SATB2, na.rm = TRUE) / sqrt(sum(!is.na(SATB2))))

# Fig 2i - gene body %mCH correlation to Ecker data
#### re-naming cell types according to ecker reference Luo 2022
ref <- readRDS("~/Library/CloudStorage/OneDrive-OregonHealth&ScienceUniversity/amethyst/atlas data/Luo 2022/PMC9004682_5972genes_BA10_ref.RData")
type_ch <- aggregateMatrix(brain, "gene_ch", "type")
both <- merge(ref, type_ch, by = 0) |> tibble::column_to_rownames(var = "Row.names")
both <- both[rownames(both) %in% mmc9$`Gene Name`, ] # DEGs from BA10 cell types
cor <- cor(both)
cor <- cor[c(1:17), c(18:29)]
pheatmap(cor)

###############################################################################
#### FIGURE 3 - DMR analysis ###

# Fig 3a - pbmc top mCG hit (heatMap)
heatMap(pbmc, regions = c("chr7_130938000_130962500", "chr19_54905000_54909000", "chr14_105853000_105865500"), trackOverhang = 2000, extraTracks = list(ccre_track), 
        track = "cg_type_tracks", remove = "MIR|ENS", order = rev(pbmc_order), colors = c("#005eff", "#9daec9", "#cccccc", "#dbdbdb"), arrowOverhang = 500)
pbmc_dmr_hypo_tracks <- split(pbmc_collapsed_dmrs |> dplyr::filter(direction == "hypo") |> ungroup() |> dplyr::select(chr, dmr_start, dmr_end, type) |> dplyr::rename(start = dmr_start, end = dmr_end) |> data.table(), pbmc_collapsed_dmrs[pbmc_collapsed_dmrs$direction == "hypo", ]$type)
heatMap(pbmc, regions = c("chr7_130938000_130962500", "chr19_54905000_54909000", "chr14_105853000_105865500"), trackOverhang = 2000, extraTracks = pbmc_dmr_hypo_tracks[pbmc_order], 
        track = "cg_type_tracks", remove = "MIR|ENS", order = rev(pbmc_order), colors = c("#005eff", "#9daec9", "#cccccc", "#dbdbdb"), arrowOverhang = 500)

# Fig 3b - brain top mCG hit (heatMap)
heatMap(brain, regions = c("chr9_79573500_79586000", "chr14_78169500_78174500", "chr22_37979000_37989500"), order = brain_order, trackOverhang = 2000, nrow = 3,
        track = "cg_type_tracks", remove = "MIR|ENS", colors = c("#005eff", "#9daec9", "#cccccc", "#dbdbdb"), arrowOverhang = 500, extraTracks = list(ccre_track))
brain_cg_dmr_hypo_tracks <- split(brain_cg_collapsed_dmrs |> dplyr::filter(direction == "hypo") |> ungroup() |> dplyr::select(chr, dmr_start, dmr_end, type) |> dplyr::rename(start = dmr_start, end = dmr_end) |> data.table(), brain_cg_collapsed_dmrs[brain_cg_collapsed_dmrs$direction == "hypo", ]$type)
heatMap(brain, regions = c("chr9_79573500_79586000", "chr14_78169500_78174500", "chr22_37979000_37989500"), order = brain_order, trackOverhang = 2000, nrow = 3,
        track = "cg_type_tracks", remove = "MIR|ENS", colors = c("#005eff", "#9daec9", "#cccccc", "#dbdbdb"), arrowOverhang = 500, extraTracks = brain_cg_dmr_hypo_tracks[rev(brain_order)])

# Fig 3c - brain top mCH hit (heatMap)
brain_ch_dmr_hyper_tracks <- split(brain_ch_collapsed_dmrs[brain_ch_collapsed_dmrs$direction == "hyper", ] |> ungroup() |> dplyr::select(chr, dmr_start, dmr_end, type) |> dplyr::rename(start = dmr_start, end = dmr_end) |> data.table(), brain_ch_collapsed_dmrs[brain_ch_collapsed_dmrs$direction == "hyper", "type"])
heatMap(brain, track = "ch_type_tracks", regions = c(getGeneCoords(brain@ref, "GRIK3"), "chr19_13064000_13023500", getGeneCoords(brain@ref, "SATB2")), order = brain_order, trackOverhang = 100000, 
        extraTracks = brain_ch_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"), remove = "ENS|MIR|RPL|RNU")

# Fig 3d - brain top mCH hit (histograM)
histograM(brain, track = "ch_type_tracks", regions = c("chr7_101916500_102082500"), order = rev(brain_order), trackOverhang = 100000,  remove = "ENS|MIR|RPL|RNU", baseline = "mean")
ggplot() +
  geom_rect(data = brain_ch_collapsed_dmrs |> dplyr::filter(chr == "chr7", dmr_start > (101916500 - 100000), dmr_end < (102082500 + 100000)),
            aes(xmin = dmr_start, xmax = dmr_end, ymin = 0, ymax = 1, fill = direction)) +
  geom_rect(data = brain@ref[seqid == "chr7" & start > (101916500 - 100000) & gene_name == "CUX1" & end < (102082500 + 100000) & type == "exon", ], 
            aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = type)) +
  facet_wrap(~ type, ncol = 1, scales = "free_y") + scale_fill_manual(values = c("black", "#ffa600", "#0073bf")) + xlim((101916500 - 100000), (102082500 + 100000)) + theme_classic()

histograM(brain, track = "ch_type_tracks", regions = c("chr15_77635500_77818000"), order = rev(brain_order), trackOverhang = 100000,  remove = "ENS|MIR|RPL|RNU", baseline = "mean", arrowScale = 0.005)
ggplot() + 
  geom_rect(data = brain_ch_collapsed_dmrs |> dplyr::filter(chr == "chr15", dmr_start > (77635500 - 100000), dmr_end < (77818000 + 100000)),
            aes(xmin = dmr_start, xmax = dmr_end, ymin = 0, ymax = 1, fill = direction)) + 
  geom_rect(data = brain@ref[seqid == "chr15" & start > (77635500 - 100000) & gene_name == "LINGO1" & end < (77818000 + 100000) & type == "exon", ], 
            aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = type)) +
  scale_fill_manual(values = c("black", "#ffa600", "#0073bf")) +
  facet_wrap(~ factor(type, levels = rev(brain_order)), ncol = 1, scales = "free_y") + theme_classic()

# Fig 3e - pbmc DMR genehancer enrichment
ggplot(pbmc_overlap_result |> dplyr::mutate(n = ifelse(direction == "hypo", (n * -1), n)), 
       aes(y = factor(member_id, levels = rev(pbmc_order)), x = n, fill = class)) + geom_col() + 
  scale_fill_manual(values = c("#eb4034", "#52a1a3", "#f6cb52")) + theme_classic() + scale_x_continuous(breaks = c(-1500, -1000, -500, 0, 500, 1000))

# Fig 3f - brain DMR genehancer enrichment
ggplot(brain_overlap_result |> dplyr::mutate(n = ifelse(direction == "hypo", (n * -1), n)),
       aes(y = factor(type, levels = brain_order), x = n, fill = class)) + geom_col() + 
  scale_fill_manual(values = c("#eb4034", "#52a1a3", "#f6cb52")) + theme_classic() + scale_x_continuous(breaks = c(-3000, -2500, -2000, -1500, -1000, -500, 0, 500, 1000))

# Fig 3g - pbmc hypo mCG GO terms
ggplot(pbmc_resultElim[sample(x = 1:nrow(x = pbmc_resultElim)), ] |> dplyr::filter(annotated < 150 & expected > 0.3 & direction == "hypo"), aes(x = -log10(fisher), y = log10(fold_change), color = type, size = 1 / rank_total)) + 
  geom_point() + scale_color_manual(values = pbmc_type_colors) + theme_classic() + ggrepel::geom_text_repel(aes(label = term))

# Fig 3h - brain hypo mCG GO terms
ggplot(brain_resultElim[sample(x = 1:nrow(x = brain_resultElim)), ] |> dplyr::filter(annotated < 150 & expected > 0.3 & direction == "hypo"), aes(x = -log10(fisher), y = log10(fold_change), color = type, size = 1 / rank_total)) + 
  geom_point() + scale_color_manual(values = brain_type_colors) + theme_classic() + ggrepel::geom_text_repel(aes(label = term)) 

# Fig 3i - colored manually

# Fig 3j - (benchmark) venn of benchmark DMR results
ggvenn(venn_df_dms, c("Amethyst", "MethScan", "AllCools"), fill_color = c("#487ca3","#b80245","#7558a1"))

# Fig 3k - (benchmark) ridgeplot of DMR sizes
ggplot(all_regions_df |> dplyr::filter(package != "amethyst_high"), aes(x = length, y = package, fill = package)) + theme_classic() + geom_density_ridges(alpha = 0.8) + scale_x_continuous(trans = 'log10') + scale_fill_manual(values = c("#487ca3","#b80245","#7558a1"))

# Fig 3l - (benchmark) genehancer comparison for DMRs
dmrs <- copy(all_regions_df)
setkey(dmrs, chr, dmr_start, dmr_end)
overlap_result <- foverlaps(dmrs, genehancer, type="any", nomatch=NA)
overlap_result <- unique(overlap_result, by = c("chr", "dmr_start", "dmr_end", "class", "package"))
overlap_result <- overlap_result |> dplyr::group_by(package, class) |> dplyr::summarise(n = n())
ggplot(overlap_result, aes(y = n, x = package, fill = class)) + geom_col(position = "fill") + scale_fill_manual(values = c("#eb4034", "#52a1a3", "#f6cb52")) + theme_classic()

# Fig 3m - (benchmark) heatMap example of DMRs
brain@tracks[["cg_type_neurons"]] <- brain@tracks[["cg_type_tracks"]][, c(1:3, 6:8, 10:14)]
heatMap(brain,  track = "cg_type_neurons", regions = c(getGeneCoords(brain@ref, "LINGO1"), getGeneCoords(brain@ref, "SLC6A1")), 
        order = c("Inh_CGE", "Inh_MGE", "Exc_L1-3_CUX2", "Exc_L2-4_RORB", "Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L5-6_PDZRN4", "Exc_L6_TLE4"),
        extraTracks = dmr_tracks, extraTrackColors = rev(c("#487ca3","#b80245","#7558a1")), remove = "ENS", trackOverhang = 50000,
        colors = c("#005eff", "#9daec9", "#cccccc", "#dbdbdb"))

# Fig 3n - (benchmark) zoomed-in heatMap examples
heatMap(brain,  track = "cg_highres_lingo1_subset", regions = getGeneCoords(brain@ref, "LINGO1-AS1"), trackOverhang = 2500, colors = c("#0073bf", "grey60", "#eb4034"),
        extraTracks = dmr_tracks, extraTrackColors = rev(c("#487ca3","#b80245","#7558a1")), 
        order = c("Inh_CGE", "Inh_MGE", "Exc_L1-3_CUX2", "Exc_L2-4_RORB", "Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L5-6_PDZRN4", "Exc_L6_TLE4"))
brain@tracks[["cg_highres_slc6a1_subset"]] <- readRDS("~/Desktop/mount/sciMETv2_3842F_cg_highres_tracks_slc6a1.RData")
heatMap(brain,  track = "cg_highres_slc6a1_subset", regions = c("chr3_10949000_10950000"), trackOverhang = 5000, colors = c("#0073bf", "grey60", "#eb4034"),
        extraTracks = dmr_tracks, extraTrackColors = rev(c("#487ca3","#b80245","#7558a1")), 
        order = c("Inh_CGE", "Inh_MGE", "Exc_L1-3_CUX2", "Exc_L2-4_RORB", "Exc_L4-5_FOXP2", "Exc_L4-6_LRRK1", "Exc_L5-6_PDZRN4", "Exc_L6_TLE4"))

###############################################################################
### FIGURE 4 - Glia hyper-mCH analysis ###

# Fig 4a - total numbers
ggplot(brain_ch_dm_counts, aes(x = log10(n_dmrs), y = log10(n_dmgs), color = type, size = pct_unique)) + geom_point() + 
  facet_wrap(vars(direction)) + scale_color_manual(values = brain_type_colors) + theme_classic() + xlim(0, 5.2) + ylim(0, 4.2)

# Fig 4b - overlap w/ Lister
ggvenn(list('this_study' = c(glia_hyper_ch_dmgs$gene), 
            'lister' = c(glia_hyper_mch_lister$gene)), c('this_study', 'lister'), fill_color = c("#7558a1","#f5985b"))

# Fig 4b - dot plot of genes overlapping w/ Lister
# plotting both in custom dot plot
glia_hyper_mch_validated$gene <- factor(glia_hyper_mch_validated$gene, levels = c("gene_average", "PCDHGC3", "PCDHGC4", "PCDHGC5", "TBR1", setdiff(validated, c("PCDHGC3", "PCDHGC4", "PCDHGC5", "TBR1"))))
ggplot(glia_hyper_mch_validated, aes(x = gene, y = type)) + geom_point(aes(size = pct, color = norm)) +
  theme_classic() + scale_size(range = c(.1, 10)) + scale_color_viridis_c() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Fig 4c - top astro results (heatMap)
heatMap(brain, regions = getGeneCoords(brain@ref, gene = "BMI1"), order = brain_order, trackOverhang = 10000, track = "ch_type_tracks", extraTracks = brain_ch_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))
heatMap(brain, regions = getGeneCoords(brain@ref, gene = "CTTNBP2"), order = brain_order, trackOverhang = 5000, track = "ch_type_tracks", extraTracks = brain_ch_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))
heatMap(brain, regions = getGeneCoords(brain@ref, gene = "MAP3K3"), order = brain_order, trackOverhang = 5000, track = "ch_type_tracks", extraTracks = brain_ch_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))
heatMap(brain, regions = getGeneCoords(brain@ref, gene = "POU4F1"), order = brain_order, trackOverhang = 12500, track = "ch_type_tracks", extraTracks = brain_ch_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))

# Fig 4d - top astro results (dimM)
# plot collective increase over astro hyper-mCH genes
astro_genemodule <- tidyr::pivot_longer(brain@genomeMatrices[["gene_ch_norm"]][glia_hyper_ch_dmgs$gene[glia_hyper_ch_dmgs$type == "Astro"], ] |> rownames_to_column(var = "gene"), cols = -c("gene"), values_to = "norm", names_to = "cell_id") |> 
  dplyr::group_by(cell_id) |> dplyr::summarise(astro_genemodule = mean(norm, na.rm = TRUE)) |> column_to_rownames(var = "cell_id")
brain <- addMetadata(brain, astro_genemodule)
dimM(brain, genes = c("BMI1"), matrix = "gene_ch_extended", pointSize = 0.3, blend = F, squish = 8, colors = c("#dbdbdb", "#cccccc", "#265A61", "#0C1E20"), reduction = "umap_100k_cg_ch")
dimFeature(brain, colorBy = astro_genemodule, pointSize = 0.3, reduction = "umap_100k_cg_ch") + scale_color_gradientn(name = "astro", colors = c("black", "turquoise", "gold", "red"), limits = c(0, 4)) + ggtitle("astro module")

# Fig 4e - top oligo results (heatMap)
heatMap(brain, regions = c("chr5_141496500_141517500"), order = brain_order, trackOverhang = 15000, track = "ch_type_tracks", extraTracks = brain_ch_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))
heatMap(brain, regions = c("chr1_222710000_222742000"), order = brain_order, trackOverhang = 0000, track = "ch_type_tracks", extraTracks = brain_ch_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))
heatMap(brain, regions = c("chr14_28779000_28783500"), order = brain_order, trackOverhang = 5000, track = "ch_type_tracks", extraTracks = brain_ch_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))
heatMap(brain, regions = c("chr5_88815500_88818000"), order = brain_order, trackOverhang = 10000, track = "ch_type_tracks", extraTracks = brain_ch_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))

# Fig 4f - top oligo results (dimM)
# plot collective increase over oligo hyper-mCH genes
oligo_genemodule <- tidyr::pivot_longer(brain@genomeMatrices[["gene_ch_norm"]][glia_hyper_ch_dmgs$gene[glia_hyper_ch_dmgs$type == "Oligo"], ] |> rownames_to_column(var = "gene"), cols = -c("gene"), values_to = "norm", names_to = "cell_id") |> 
  dplyr::group_by(cell_id) |> dplyr::summarise(oligo_genemodule = mean(norm, na.rm = TRUE)) |> column_to_rownames(var = "cell_id")
brain <- addMetadata(brain, oligo_genemodule)
dimM(brain, genes = c("BROX"), matrix = "gene_ch_extended", pointSize = 0.3, blend = F, squish = 7, colors = c("#dbdbdb", "#cccccc", "#265A61", "#0C1E20"), reduction = "umap_100k_cg_ch")
dimFeature(brain, colorBy = oligo_genemodule, pointSize = 0.3, reduction = "umap_100k_cg_ch") + scale_color_gradientn(name = "oligo", colors = c("black", "turquoise", "gold", "red"), limits = c(0, 4)) + ggtitle("oligo module")

# Fig 4g - correlation to RNA
pheatmap(rna_summary, 
         color = colorRampPalette(rev(c("#D7191C", "#fd9061", "#FDAE61", "#FFFFBF", "#ABD9E9", "#2C7BB6")))(100),
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation")

###############################################################################
### FIGURE 5 - Atlas data analysis ###

# Fig 5a - atlas data w/ cell type key
dimFeature(atlas, colorBy = type, colors = ecker_type_colors, pointSize = 0.1,reduction = "umap") + theme(legend.position = "none")
ggplot(atlas@metadata, aes(x = log(cov))) + geom_histogram(bins = 40) + xlim(13, 17) + ylim(0, 10000)

# Fig 5b - atlas data %mCH
dimFeature(atlas, colorBy = log(mch_pct), pointSize = 0.1) + theme(legend.position = "none") + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"), limits = c(0, 2.5), oob = scales::squish)

# Fig 5c - atlas data cell type composition
sampleComp(atlas, groupBy = "sample", colorBy = "type", colors = ecker_type_colors) + facet_wrap(vars(seq), nrow = 1)
# iLISI score to test for batch effects - https://github.com/immunogenomics/LISI
lisi_scores <- lisi::compute_lisi(atlas@reductions[["umap"]], atlas@metadata, label_colnames = c('sample', 'type', 'seq'), perplexity = 50) # 50 also used in runCluster
summary(lisi_scores$sample) #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.  1.005   2.175   2.724   2.684   3.202   3.999 
summary(lisi_scores$type) #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 1.000   1.000   1.011   1.150   1.101   6.039 
summary(lisi_scores$seq) # summary(lisi_scores$seq) Min. 1st Qu.  Median    Mean 3rd Qu.    Max.  1.000   1.495   1.727   1.672   1.904   2.000

# Fig 5d - atlas data by sample
dimFeature(atlas, colorBy = type, colors = ecker_type_colors, pointSize = 0.1) + theme(legend.position = "none") + facet_wrap(vars(sample), nrow = 1)

# Fig 5e - time/memory stats
ggplot(downsample, aes(x = log10(cells), y = log10(seconds), color = param)) + geom_point(size = 2) + geom_line(linewidth = 1.5, aes(linetype = context)) + theme_classic() + scale_color_manual(values = rev(c("#264653", "#298B83", "#9CB478", "#edbe66", "#f28e4b", "#db4f40"))) 
ggplot(downsample_summary, aes(x = log10(cells), y = mem, color = context)) + geom_point(aes(size = norm_seconds)) + geom_line(linewidth = 1) + theme_classic() + scale_size(breaks = c(307, 1000, 2418), range = c(1, 7)) + scale_color_manual(values = c("#373f51", "#f54e42")) + scale_y_log10()

# Fig 5f - dimM Brox vs CTTBP2
genem <- as.data.frame(t(atlas@genomeMatrices[["gene_ch_norm_cttnpb2_brox"]])) |> filter(!is.na(BROX) & !is.na(CTTNBP2)) |>
  mutate(BROX = pmin(BROX, 4) / 4, 
         CTTNBP2 = pmin(CTTNBP2, 4) / 4) # use pmin to scale to a max for extreme outliers
plot <- merge(genem, atlas@metadata, by = 0) |> dplyr::mutate(mix = grDevices::rgb(red = BROX, green = CTTNBP2, blue = 0, maxColorValue = 1))
p1 <- ggplot2::ggplot(plot, ggplot2::aes(x = .data[["umap_x"]], y = .data[["umap_y"]], color = mix)) + ggplot2::geom_point(size = .1) + 
  ggplot2::theme_classic() + ggplot2::guides(colour = guide_legend(override.aes = list(size = 3))) + ggplot2::scale_color_identity()  + noAxes()

legend <- expand.grid(red = seq(0, 1, by = 0.1), green = seq(0, 1, by = 0.1))
legend <- within(legend, mix <- grDevices::rgb(green = green, red = red, blue = 0, maxColorValue = 1))
p2 <- ggplot2::ggplot(legend, ggplot2::aes(x=red, y=green)) +
  ggplot2::geom_tile(aes(fill=mix), color="white") +
  ggplot2::scale_fill_identity() +
  ggplot2::theme_classic() +
  ggplot2::labs(x = "BROX", y = "CTTNBP2", title = "Legend")
p2 <- p2 + coord_fixed(ratio = 1/1)
cowplot::plot_grid(p1, p2, rel_widths = c(3,1))

# statistics 05/01/25
test <- merge(as.data.frame(t(atlas@genomeMatrices[["gene_ch_norm_cttnpb2_brox"]])), atlas@metadata, by = 0)
test <- test |> dplyr::mutate(astro = ifelse(type %in% c("Astro_1", "Astro_2"), "YES", "NO"), oligo = ifelse(type == "Oligo", "YES", "NO"))
shapiro.test(sample(test$CTTNBP2, size = 5000)) # W = 0.70252, p-value < 2.2e-16
wilcox.test(test$CTTNBP2[test$astro == "YES"], test$CTTNBP2[test$astro == "NO"]) # W = 2054098960, p-value < 2.2e-16
test |> dplyr::group_by(astro) |> dplyr::summarise(mean = mean(CTTNBP2, na.rm = TRUE), se = sd(CTTNBP2, na.rm = TRUE) / sqrt(sum(!is.na(CTTNBP2))))

shapiro.test(sample(test$BROX, size = 5000)) # W = 0.50974, p-value < 2.2e-16
wilcox.test(test$BROX[test$oligo == "YES"], test$BROX[test$oligo == "NO"]) # W = 482049861, p-value < 2.2e-16
test |> dplyr::group_by(oligo) |> dplyr::summarise(mean = mean(BROX, na.rm = TRUE), se = sd(BROX, na.rm = TRUE) / sqrt(sum(!is.na(BROX))))

# Fig 5g - bar chart of mean %mCHH in each trinucleotide context for hyper and non hyper-mCH genes
plot_grid(plot_chn_hyper_mch_only, plot_chn_not_hyper_mch, nrow = 2)

# Fig 5h - bar chart of n DMGs between individuals identified per chromosome
ggplot(atlas_mch_gene_individual_results |> dplyr::group_by(chr, sex, individual, direction) |> dplyr::summarise(n = n()) |> dplyr::mutate(n = ifelse(direction == "hypo", (-1*n), n)), 
       aes(x = sex, y = n, fill = individual)) + geom_col() + facet_grid(cols = vars(chr), rows = vars(direction)) + theme_classic() + scale_fill_manual(values = c("#143642","#EC9928","#128B8C","#A82023"))

# Fig 5i - top mCH DMGs per individual per type
ggplot(atlas_mch_gene_individual_results |> dplyr::mutate(p.adj = -log10(p.adj)), aes(x = logFC, y = p.adj, color = individual)) + 
  geom_point() + scale_color_manual(values = c("#143642","#EC9928","#128B8C","#A82023")) + theme_classic() + 
  new_scale_color() + geom_text_repel(aes(label = gene, color = type)) + scale_color_manual(values = ecker_type_colors) 

# Fig 5j - mCH type x individual heatMap of USP9X
order2 <- rev(paste0(rev(atlas_order), "_", rep(individuals, each = 18)))
heatMap(obj = atlas, 
        track = "type_individual_ch_pct", 
        order = order2,
        genes = "USP9X", 
        trackOverhang = 200000, colorMax = 15)

###############################################################################
### FIGURE S1 - Doublet analysis ###

set.seed(111)
# Remove doublets
doublets <- pbmc
dbobj <- makeDoubletObject(doublets, simFraction=0.25, threads = 10, genomeMatrices=c("cg_100k_score"))
dbobj@reductions[["irlba"]] <- runIrlba(dbobj, genomeMatrices=c('cg_100k_score'), dims=c(10))
result <- buildDoubletModel(dbobj, reduction = "irlba", method="rf")
dbobj <- predictDoubletScores(dbobj, reduction = "irlba", model = result$model)
dbobj <- runCluster(dbobj, k = 30, reduction = "irlba")
dbobj <- runUmap(dbobj, neighbors = 30, dist = 0.1, method = "euclidean", reduction = "irlba")
doublets <- addDoubletScores(doublets, dbobj)
doublets@metadata <- doublets@metadata |> dplyr::mutate(doublet = ifelse(doublet_score.y > 0.5, "YES", "NO"))
# save.image("/secret/path/amethyst/doublet_test/doublet_workspace.RData")

dimFeature(doublets, colorBy = cluster_id, colors = makePalette(19, 10))
dimFeature(dbobj, colorBy = cluster_id, colors =  makePalette(19, 15))
dimFeature(dbobj, colorBy = doublet_score) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red"))
dimFeature(dbobj, colorBy = doublet_score) + scale_color_gradientn(colors = c("black", "turquoise", "gold", "red")) + facet_wrap(vars(doublet_info))
dimFeature(doublets, colorBy = doublet, colors = c("grey80", "red"))
ggplot(dbobj@metadata, aes(x = reorder(doublet_info, doublet_score), y = doublet_score)) + geom_boxplot(outlier.size = 0, outlier.alpha = 0) + geom_jitter(size = 0.3) + theme_classic()

shapiro.test(dbobj@metadata$doublet_score) # W = 0.66428, p-value < 0.00000000000000022
wilcox.test(dbobj@metadata$doublet_score ~ dbobj@metadata$doublet_info) # W = 2427350, p-value < 0.00000000000000022

dbobj@metadata |> dplyr::group_by(doublet_info) |>  dplyr::summarise(mean = mean(doublet_score, na.rm = TRUE),se = sd(doublet_score, na.rm = TRUE) / sqrt(n()))
dbobj@metadata |> dplyr::mutate(pass = ifelse(doublet_score > 0.5, "FAIL", "PASS")) |> dplyr::group_by(doublet_info, pass) |> dplyr::summarise(n())
chisq.test(data.table(c(695, 90), c(14, 3124))) # X-squared = 3284.8, df = 1, p-value < 0.00000000000000022

###############################################################################
### FIGURE S2 - GAD1 and SATB2 results, not imputed ###

dimM(brain, genes = c("SATB2", "GAD1"), matrix = "gene_ch_imputed", blend = T, reduction = "umap_100k_cg_ch", pointSize = 0.4)
dimM(brain, genes = c("SATB2", "GAD1"), matrix = "gene_ch", blend = F, reduction = "umap_100k_cg_ch", pointSize = 0.4, squish = 10, colors = c("#dbdbdb", "#cccccc", "#265A61", "#0C1E20"))

# statistics
test <- merge(brain@metadata, as.data.frame(t(brain@genomeMatrices[["gene_ch_extended"]][c("SATB2", "GAD1"), ])), by = 0)
shapiro.test(test$GAD1) # W = 0.77476, p-value < 2.2e-16
wilcox.test(test$GAD1[test$major_type == "exc"], test$GAD1[test$major_type != "exc"]) # W = 54280, p-value < 2.2e-16
test |> dplyr::mutate(group = case_when(major_type == "exc" ~ "exc", TRUE ~ "non-exc")) |> group_by(group) |> dplyr::summarise(
  mean_GAD1 = mean(GAD1, na.rm = TRUE), n = n(),
  se_GAD1 = sd(GAD1, na.rm = TRUE) / sqrt(sum(!is.na(GAD1))))

shapiro.test(test$SATB2) # W = 0.71413, p-value < 2.2e-16
wilcox.test(test$SATB2[test$major_type == "inh"], test$SATB2[test$major_type != "inh"]) # W = 57048, p-value < 2.2e-16
test |> dplyr::mutate(group = case_when(major_type == "inh" ~ "inh", TRUE ~ "non-inh")) |> group_by(group) |> dplyr::summarise(
  mean_SATB2 = mean(SATB2, na.rm = TRUE), n = n(),
  se_SATB2 = sd(SATB2, na.rm = TRUE) / sqrt(sum(!is.na(SATB2))))

###############################################################################
### FIGURE S3 - mCG and mCH relationship ###
brain_cg_dmr_hyper_tracks <- split(brain_cg_collapsed_dmrs[brain_cg_collapsed_dmrs$direction == "hyper", ] |> ungroup() |> dplyr::select(chr, dmr_start, dmr_end, type) |> dplyr::rename(start = dmr_start, end = dmr_end) |> data.table(), brain_cg_collapsed_dmrs[brain_cg_collapsed_dmrs$direction == "hyper", "type"])
brain_cg_dmr_hypo_tracks <- split(brain_cg_collapsed_dmrs[brain_cg_collapsed_dmrs$direction == "hypo", ] |> ungroup() |> dplyr::select(chr, dmr_start, dmr_end, type) |> dplyr::rename(start = dmr_start, end = dmr_end) |> data.table(), brain_cg_collapsed_dmrs[brain_cg_collapsed_dmrs$direction == "hypo", "type"])

# Fig S3a - top astro results (heatMap)
heatMap(brain, regions = getGeneCoords(brain@ref, gene = "BMI1"), order = brain_order, trackOverhang = 10000, track = "cg_type_tracks", extraTracks = brain_cg_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))
heatMap(brain, regions = getGeneCoords(brain@ref, gene = "CTTNBP2"), order = brain_order, trackOverhang = 5000, track = "cg_type_tracks", extraTracks = brain_cg_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))
heatMap(brain, regions = getGeneCoords(brain@ref, gene = "MAP3K3"), order = brain_order, trackOverhang = 5000, track = "cg_type_tracks", extraTracks = brain_cg_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))
heatMap(brain, regions = getGeneCoords(brain@ref, gene = "POU4F1"), order = brain_order, trackOverhang = 12500, track = "cg_type_tracks", extraTracks = brain_cg_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))

# Fig S3b - top oligo results (heatMap)
heatMap(brain, regions = c("chr5_141496500_141517500"), order = brain_order, trackOverhang = 15000, track = "cg_type_tracks", extraTracks = brain_cg_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))
heatMap(brain, regions = c("chr1_222710000_222742000"), order = brain_order, trackOverhang = 000, track = "cg_type_tracks", extraTracks = brain_cg_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))
heatMap(brain, regions = c("chr14_28779000_28783500"), order = brain_order, trackOverhang = 5000, track = "cg_type_tracks", extraTracks = brain_cg_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))
heatMap(brain, regions = c("chr5_88815500_88818000"), order = brain_order, trackOverhang = 10000, track = "cg_type_tracks", extraTracks = brain_cg_dmr_hyper_tracks[setdiff(rev(brain_order), "Micro")], extraTrackColors = setdiff(rev(brain_type_colors), "#F79C46"))

# Fig S3c - relationship of mCG and mCH in hyper-mCH genes for each group
ggplot() + 
  geom_point(data = brain_dmgs, aes(x = mean_1.ch, y = mean_1.cg, color = type), size = 0.5) + 
  geom_smooth(data = brain_dmgs, aes(group = type, x = mean_1.ch, y = mean_1.cg, color = type), method = "lm", se = FALSE) +
  geom_vline(data = brain_type_summary_metrics |> dplyr::filter(!(type %in% c("Micro", "OPC"))), aes(xintercept = mean_mch_pct), linetype = "dashed", color = "gray40") +
  geom_hline(data = brain_type_summary_metrics |> dplyr::filter(!(type %in% c("Micro", "OPC"))), aes(yintercept = mean_mcg_pct), linetype = "dashed", color = "gray40") +
  geom_text(data = lm_labels, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 3, hjust = 0) +
  scale_color_manual(values = brain_type_colors) + theme_classic() + facet_wrap(vars(type), nrow = 5) 

###############################################################################
### FIGURE S4 - Atlas-scale brain dataset; mCG patterns over USP9X by cell type x individual ###
order2 <- rev(paste0(rev(atlas_order), "_", rep(individuals, each = 18)))
heatMap(obj = atlas, 
        track = "type_individual_cg_pct", 
        order = order2,
        genes = "USP9X", 
        trackOverhang = 200000,
        colors = c("#005eff", "#9daec9", "#cccccc", "#dbdbdb", "grey90"))
heatMap(obj = atlas, 
        track = "type_individual_cg_pct", 
        order = order2,
        genes = "KDM6A", 
        trackOverhang = 200000,
        colors = c("#005eff", "#9daec9", "#cccccc", "#dbdbdb", "grey90"))


