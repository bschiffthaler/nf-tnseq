suppressPackageStartupMessages({
  library(here) # Path management
  library(readr) # Better file reading
  library(Rsubread) # Alignment feature counting
  library(Tnseq) # Statistical testing
  library(stringr) # Matching sample table to files
  library(stringdist) # Matching sample table to files
  library(pander) # Render non-interactive tables
  library(dplyr) # Data wrangling
  library(tidyr) # Data wrangling
  library(ggplot2) # Plotting
  library(plotly) # Plotting
  library(ggrepel) # Plotting
  library(writexl) # Writing xlsx documents
  library(pheatmap) # Heatmaps
  library(reactable) # Interactive tables
  library(pheatmap)
  library(Gviz)
  library(corrplot)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  set_here(args[1])
}

my_theme <-theme_bw() + 
  theme(text = element_text(size = 18))

#' ## Read input data
#' We read the available data. Here's what's what:
#' 
#' * ta_sites - The locations of TA sites within the selected features (typically genes)
#' * sample_info - The sample metadata. conditions, replicates etc.
#' * bamfiles - Sequence alignments to the selected genome
#' * gois - Gene of interest that will be highlighted in plots (Optional)
#' 
ta_sites <- here("genome","genes_tasites.saf")
ta_df <- read_tsv(ta_sites,
                  col_names = c("GeneID", "Chr", "Start", "End", "Strand"))

sample_info <- read_csv(here("datafiles","sampledata.csv"), col_types = 'ccffffc')

bamfiles <- dir(here("analysis","03-align"),
                pattern = "*\\.bam$",
                full.names = TRUE)

if (file.exists(here("datafiles", "genes_of_interest.txt"))) {
  gois <- read_lines(here("datafiles", "genes_of_interest.txt"))
} else {
  gois <- character()
}

#' Check for order disparity between mapping files and the sample information data
#' and reorder if needed. **Manually verify that the order is correct**
closest <- function(a, b) {
  which.min(stringdist(a, b))
}
ind <- sapply(bamfiles, closest, sample_info$Name)
sample_info <- sample_info[ind, ]
pander(data.frame(A = sample_info$Name, B = basename(bamfiles)))


#' ## Prepare output directory structure
suppressWarnings({
  dir.create(here("analysis/","04-statistical-analysis"), recursive = TRUE)
  dir.create(here("analysis/","04-statistical-analysis", "plots"), recursive = TRUE)
  dir.create(here("analysis/","04-statistical-analysis", "tables"), recursive = TRUE)
})

#' ## Count number of reads in features
#' Here we count how many reads overlap each feature (i.e. a TA site). The default
#' is to allow a read to be counted for multiple TA sites. If you set
#' `allowMultiOverlap = FALSE` then only reads which are overlapping only a single
#' TA site would be counted. Change this only if you know what you are doing. In
#' the following example Read 3 would never be counted since it overlaps two
#' TA sites and Read 2 could not be counted since it overlaps a single TA site in 
#' two separate genes (on either strand).
#' 
#'     Genome    -----TATA------
#'     Gene5'3'  <<<<<<<<<<<<<<<
#'     Gene3'5'  >>>>>>>
#'     Read 1    =====
#'     Read 2       ====
#'     Read 3        =====
#'     Read 4           =====
#'     
sr_counts <- featureCounts(bamfiles,
                           annot.ext = ta_df,
                           useMetaFeatures = FALSE,
                           allowMultiOverlap = TRUE)
sr_counts$countsNames <- dimnames(sr_counts$counts)
sr_counts$counts <- unname(sr_counts$counts)

#' We sum up over technical replicates even if there are none (or only one)
sr_counts$tr_counts <- t(apply(sr_counts$counts, 1, tapply,
                                      sample_info$TechnicalReplicate, sum))
sample_info_tr <- do.call(rbind,
        lapply(split(sample_info, sample_info$TechnicalReplicate), function(f) {
          x <- head(f, 1)
          x$Name <- paste(x$Condition, x$BiologicalReplicate, x$Pool, sep="-")
          x
        })
)
sr_counts$tr_counts <- sr_counts$tr_counts[, sample_info_tr$TechnicalReplicate]


#' ## Check percent disruption
x <- as_tibble(sr_counts$counts)
colnames(x) <- sample_info$Name
x$Gene <- sr_counts$annotation$GeneID
x %>% 
  pivot_longer(-Gene) %>%
  left_join(sample_info, by = c("name" = "Name")) -> x
x %>% 
  group_by(Gene, Condition) %>%
  summarise(PercentDisrupted = sum(value > 0) / length(value) * 100) %>%
  ggplot(aes(x = PercentDisrupted, fill = Condition)) +
  geom_histogram(bins = 50) +
  facet_wrap(~Condition) +
  my_theme -> pl
ggsave(here("analysis", "04-statistical-analysis/", "plots", "percent_disruption.pdf"),
       width = 12, height = 10, plot = pl)
ggplotly(pl)

#' ## Perform statistical testing if we have more than one condition
#' 
#' See https://dx.doi.org/10.1186%2Fs12859-017-1745-2 for details.
#' 
tnseq <- TnseqDiff(countData = sr_counts$tr_counts,
                   geneID = ta_df$GeneID,
                   location = ta_df$Start,
                   pool = sample_info_tr$Pool,
                   condition = sample_info_tr$Condition)

#' Adjust P values for multiple testing. We use the Benjamini Hochberg procedure.
#' See https://en.wikipedia.org/wiki/False_discovery_rate for details.
tnseq$resTable$padj <- p.adjust(tnseq$resTable$pvalue, method = 'BH')

#' Order by most significant genes. NOTE that the pvalue can never be actually 0,
#' but sometimes it is lower than the rounding error of the computer and gets 
#' represented as a 0. In your results, these should be described as $p<<10^{-150}$
#' to avoid confusion.
tnseq$resTable <- tnseq$resTable[order(tnseq$resTable$padj),]

#' Write results tables as XLSX and CSV
write_xlsx(tnseq$resTable,
           here("analysis", "04-statistical-analysis/", "tables",
                "results.xlsx"))
write_csv(tnseq$resTable,
          here("analysis", "04-statistical-analysis/", "tables",
               "results.csv"))
write_xlsx(tnseq$est.insertion,
           here("analysis", "04-statistical-analysis/", "tables",
                "est_insertion.xlsx"))
write_csv(tnseq$est.insertion,
          here("analysis", "04-statistical-analysis/", "tables",
               "est_insertion.csv"))

reactable::reactable(tnseq$resTable, filterable = TRUE, compact = TRUE,
                     striped = TRUE)

#' ## Plots
#' 
#' ### Sample PCA
#' A principal component analysis of samples based on read counts from TA insertions
#' using `limma::voom` for data transformation. The transformation is not optimal 
#' for this kind of data and should be taken with a grain of salt.
gene_sums <- apply(sr_counts$tr_counts, 2, tapply, sr_counts$annotation$GeneID, sum)
voom_data <- limma::voom(gene_sums)$E
pca <- prcomp(t(voom_data))
tibble(PC1 = pca$x[, 1], PC2 = pca$x[, 2]) %>%
  cbind(sample_info_tr) %>%
  ggplot(aes(x = PC1, y = PC2, col = Condition,
             pch = BiologicalReplicate)) +
  geom_point(size = 4) + 
  my_theme -> pl
ggsave(here("analysis", "04-statistical-analysis/", "plots", "sample_pca.pdf"),
       width = 12, height = 10, plot = pl)
ggplotly(pl)

#' ### Technical replicate correlation plot
gene_sums <- apply(sr_counts$counts, 2, tapply, sr_counts$annotation$GeneID, sum)
voom_data <- limma::voom(gene_sums)$E
colnames(voom_data) <- sample_info$Name
cors <- cor(voom_data)
corrplot(cors, type = "lower", is.corr = FALSE)

#' ### Volcano plot
#' Top 20 genes are labeled in the print version of this plot

fc_name <- sym(colnames(tnseq$resTable)[7])
as_tibble(tnseq$resTable) %>%
  ggplot(aes(x = !!fc_name,
             y = -log10(padj + .Machine$double.eps),
             label = ID)) +
  geom_label_repel(data = head(tnseq$resTable, 20), direction = 'both') +
  geom_point() +
  my_theme -> pl
ggsave(here("analysis", "04-statistical-analysis/", "plots", "volcano_plot.pdf"),
       width = 12, height = 10, plot = pl)
suppressWarnings(ggplotly(pl))

get_goidx <- function(rnames, gois) {
  if (length(gois) > 1) {
    gois_idx <- sapply(str_to_lower(rnames), str_detect, str_to_lower(gois))
    if (any(colSums(gois_idx) > 1)) {
      warning("Ambiguous gene names in gene of interest! Please provide unique names")
    }
    return (apply(gois_idx, 1, which))
  } else {
    return (which(str_detect(str_to_lower(rnames), str_to_lower(gois))))
  }
}

#' ### Volcano plot with GsOI highlighted
if (length(gois) > 0) {
  gois_idx <- get_goidx(tnseq$resTable$ID, gois)
  as_tibble(tnseq$resTable) -> j
  j$GeneOfInterest <- FALSE
  j$GeneOfInterest[gois_idx] <- TRUE
  j %>%
    ggplot(aes(x = !!fc_name,
               y = -log10(padj + .Machine$double.eps),
               label = ID,
               col = GeneOfInterest)) +
    geom_text_repel(data = j[gois_idx, ], direction = 'both') +
    geom_point() +
    my_theme -> pl
  ggsave(here("analysis", "04-statistical-analysis/", "plots", "volcano_gois_plot.pdf"),
         width = 12, height = 10, plot = pl)
  suppressWarnings(ggplotly(pl))
}

#' ### Volcano plot colored by chromosome
as_tibble(tnseq$resTable) %>%
  left_join(ta_df, by = c("ID" = "GeneID")) %>%
  ggplot(aes(x = !!fc_name,
             y = -log10(padj + .Machine$double.eps),
             label = ID,
             col = Chr)) +
  geom_point() +
  my_theme -> pl
ggsave(here("analysis", "04-statistical-analysis/", "plots", "volcano_chr_plot.pdf"),
       width = 12, height = 10, plot = pl)
suppressWarnings(ggplotly(pl))

#' ### Distribution of adjusted P-values
ggplot(tnseq$resTable, aes(x = padj)) +
  geom_histogram(bins = 100) +
  my_theme -> pl
ggsave(here("analysis", "04-statistical-analysis/", "plots", "distribution_of_padj.pdf"),
       width = 12, height = 10, plot = pl)
suppressWarnings(ggplotly(pl))

#' ### Count data heatmap
gene_sums <- apply(sr_counts$tr_counts, 2, tapply, sr_counts$annotation$GeneID, sum)
voom_data <- limma::voom(gene_sums)$E
top50_idx <- head(tnseq$resTable, 50)[, "ID"]
ann_col <- as.data.frame(sample_info_tr[, c("Condition", "BiologicalReplicate")])
rownames(ann_col) <- sample_info_tr$Name
colnames(voom_data) <- rownames(ann_col)
pdf(here("analysis", "04-statistical-analysis/", "plots", "count_heatmap.pdf"),
    width = 12, height = 10)
pheatmap(voom_data[top50_idx, ], annotation_col = ann_col)
dev.off()

#' ### Count data heatmap for genes of interest
if (length(gois) > 0) {
  gois_idx <- get_goidx(rownames(voom_data), str_to_lower(gois))
  pdf(here("analysis", "04-statistical-analysis/", "plots", "count_heatmap_gois.pdf"),
      width = 12, height = 10)
  pheatmap(voom_data[gois_idx, ,drop=FALSE], annotation_col = ann_col,
           cluster_rows = ifelse(length(gois) > 1, TRUE, FALSE))
  dev.off()
}

#' ### Coverage region plots for genes of interest
#' 
#' Setup for Gviz
options(ucscChromosomeNames=FALSE) 
gff <- rtracklayer::readGFFAsGRanges(here("genome", "genome.gff3"))
gff$feature <- gff$type
gff$transcript <- gff$Parent
gff$gene <- gff$Parent
gff$exon <- gff$ID
gff$symbol <- gff$Name

dna <- Biostrings::readDNAStringSet(here("genome","genome.fa"))
names(dna) <- sapply(str_split(names(dna), ' '), '[[', 1)
sinfo <- Seqinfo(names(dna), seqlengths = Biostrings::width(dna))

if (length(gois) > 0) {
  gois_idx <- get_goidx(str_to_lower(gff$ID), str_to_lower(gois))
  for (f in gois_idx) {
    r <- gff[f]
    chr <- chr <- as.character(unique(seqnames(r)))
    s <- start(r)
    e <- end(r)
    ta_sub_df <- filter(ta_df, GeneID == gff[f]$ID)
    ta_gr <- tibble(start = ta_sub_df$Start, end = ta_sub_df$End, chromosome = ta_sub_df$Chr,
                    strand = ta_sub_df$Strand)
    gtrack <- GenomeAxisTrack(name = as.character(seqnames(gff[f])))
    atrack <- GeneRegionTrack(gff[gff$feature=="CDS"], name = "Gene/CDS")
    altrack <- lapply(seq_along(bamfiles), function(f){
      DataTrack(bamfiles[f], name = sample_info$Name[f], stacking = "hide",
                type="heatmap", ylim = c(0, 100))
    })
    tatrack <- DetailsAnnotationTrack(ta_gr, name = "TA Sites", stacking = "full")
    plotTracks(c(gtrack, tatrack, atrack, altrack), from = s, to = e, chromosome = chr,
               transcriptAnnotation = "gene", shape = "arrow", extend.right = 0.1,
               extend.left = 0.1, background.title = "black")
    pdf(here("analysis", "04-statistical-analysis/", "plots",
             paste0("coverage_",gff[f]$ID,".pdf")), width = 12, height = 10)
    plotTracks(c(gtrack, tatrack, atrack, altrack), from = s, to = e, chromosome = chr,
               transcriptAnnotation = "gene", shape = "arrow", extend.right = 0.1,
               extend.left = 0.1, background.title = "black")
    dev.off()
  }
}
