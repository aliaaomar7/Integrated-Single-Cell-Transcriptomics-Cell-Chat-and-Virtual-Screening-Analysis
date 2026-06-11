dir.create("C:/Users/aliaa/Desktop/IPF1")

download.file(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891627/suppl/GSM3891627_IPF1_All_matrix.mtx.gz",
  "C:/Users/aliaa/Desktop/IPF1/matrix.mtx.gz"
)


options(timeout = 300)
download.file(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891627/suppl/GSM3891627_IPF1_All_matrix.mtx.gz",
  "C:/Users/aliaa/Desktop/IPF1/matrix.mtx.gz",
  mode = "wb"
)

download.file(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891627/suppl/GSM3891627_IPF1_All_barcodes.tsv.gz",
  "C:/Users/aliaa/Desktop/IPF1/barcodes.tsv.gz",
  mode = "wb"
)

download.file(
  "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891627/suppl/GSM3891627_IPF1_All_genes.tsv.gz",
  "C:/Users/aliaa/Desktop/IPF1/genes.tsv.gz",
  mode = "wb"
)

file.rename(
  "C:/Users/aliaa/Desktop/IPF1/genes.tsv.gz",
  "C:/Users/aliaa/Desktop/IPF1/features.tsv.gz"
)

IPF_I <- Read10X(data.dir = "C:/Users/aliaa/Desktop/IPF1")

# IPF2
dir.create("C:/Users/aliaa/Desktop/IPF2")

download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891629/suppl/GSM3891629_IPF2_All_barcodes.tsv.gz",
              "C:/Users/aliaa/Desktop/IPF2/barcodes.tsv.gz", mode = "wb")

download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891629/suppl/GSM3891629_IPF2_All_genes.tsv.gz",
              "C:/Users/aliaa/Desktop/IPF2/features.tsv.gz", mode = "wb")

download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891629/suppl/GSM3891629_IPF2_All_matrix.mtx.gz",
              "C:/Users/aliaa/Desktop/IPF2/matrix.mtx.gz", mode = "wb")
# IPF3
dir.create("C:/Users/aliaa/Desktop/IPF3")

download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891631/suppl/GSM3891631_IPF3_All_barcodes.tsv.gz",
              "C:/Users/aliaa/Desktop/IPF3/barcodes.tsv.gz", mode = "wb")

download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891631/suppl/GSM3891631_IPF3_All_features.tsv.gz",
              "C:/Users/aliaa/Desktop/IPF3/features.tsv.gz", mode = "wb")

download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891631/suppl/GSM3891631_IPF3_All_matrix.mtx.gz",
              "C:/Users/aliaa/Desktop/IPF3/matrix.mtx.gz", mode = "wb")

#NML1

dir.create("C:/Users/aliaa/Desktop/NML1")

download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891621/suppl/GSM3891621_NML1_All_barcodes.tsv.gz",
              "C:/Users/aliaa/Desktop/NML1/barcodes.tsv.gz", mode = "wb")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891621/suppl/GSM3891621_NML1_All_genes.tsv.gz",
              "C:/Users/aliaa/Desktop/NML1/features.tsv.gz", mode = "wb")

download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891621/suppl/GSM3891621_NML1_All_matrix.mtx.gz",
              "C:/Users/aliaa/Desktop/NML1/matrix.mtx.gz", mode = "wb")
# NML2
dir.create("C:/Users/aliaa/Desktop/NML2")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891623/suppl/GSM3891623_NML2_All_barcodes.tsv.gz",
              "C:/Users/aliaa/Desktop/NML2/barcodes.tsv.gz", mode = "wb")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891623/suppl/GSM3891623_NML2_All_genes.tsv.gz",
              "C:/Users/aliaa/Desktop/NML2/features.tsv.gz", mode = "wb")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891623/suppl/GSM3891623_NML2_All_matrix.mtx.gz",
              "C:/Users/aliaa/Desktop/NML2/matrix.mtx.gz", mode = "wb")

# NML3
dir.create("C:/Users/aliaa/Desktop/NML3")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891625/suppl/GSM3891625_NML3_All_barcodes.tsv.gz",
              "C:/Users/aliaa/Desktop/NML3/barcodes.tsv.gz", mode = "wb")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891625/suppl/GSM3891625_NML3_All_features.tsv.gz",
              "C:/Users/aliaa/Desktop/NML3/features.tsv.gz", mode = "wb")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3891nnn/GSM3891625/suppl/GSM3891625_NML3_All_matrix.mtx.gz",
              "C:/Users/aliaa/Desktop/NML3/matrix.mtx.gz", mode = "wb")
#load

library(Seurat)
library(tidyverse)
set.seed(42)

base_dir <- "C:/Users/aliaa/Desktop"

# Load
CON_I   <- Read10X(data.dir = file.path(base_dir, "NML1"))
CON_II  <- Read10X(data.dir = file.path(base_dir, "NML2"))
CON_III <- Read10X(data.dir = file.path(base_dir, "NML3"))
IPF_I   <- Read10X(data.dir = file.path(base_dir, "IPF1"))
IPF_II  <- Read10X(data.dir = file.path(base_dir, "IPF2"))
IPF_III <- Read10X(data.dir = file.path(base_dir, "IPF3"))


# Create objects
CON_I   <- CreateSeuratObject(counts = CON_I,   project = "CON", min.cells = 3, min.features = 200)
CON_II  <- CreateSeuratObject(counts = CON_II,  project = "CON", min.cells = 3, min.features = 200)
CON_III <- CreateSeuratObject(counts = CON_III, project = "CON", min.cells = 3, min.features = 200)
IPF_I   <- CreateSeuratObject(counts = IPF_I,   project = "IPF", min.cells = 3, min.features = 200)
IPF_II  <- CreateSeuratObject(counts = IPF_II,  project = "IPF", min.cells = 3, min.features = 200)
IPF_III <- CreateSeuratObject(counts = IPF_III, project = "IPF", min.cells = 3, min.features = 200)


# Mitochondrial %
CON_I   <- PercentageFeatureSet(CON_I,   pattern = "^MT-", col.name = "percent.mt")
CON_II  <- PercentageFeatureSet(CON_II,  pattern = "^MT-", col.name = "percent.mt")
CON_III <- PercentageFeatureSet(CON_III, pattern = "^MT-", col.name = "percent.mt")
IPF_I   <- PercentageFeatureSet(IPF_I,   pattern = "^MT-", col.name = "percent.mt")
IPF_II  <- PercentageFeatureSet(IPF_II,  pattern = "^MT-", col.name = "percent.mt")
IPF_III <- PercentageFeatureSet(IPF_III, pattern = "^MT-", col.name = "percent.mt")

# Filter
CON_I   <- subset(CON_I,   subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
CON_II  <- subset(CON_II,  subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
CON_III <- subset(CON_III, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
IPF_I   <- subset(IPF_I,   subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & nCount_RNA < 40000 & percent.mt < 10)
IPF_II  <- subset(IPF_II,  subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)
IPF_III <- subset(IPF_III, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 10)


# Subsample to 3000 cells per donor
CON_I   <- subset(CON_I,   cells = sample(Cells(CON_I),   min(3000, ncol(CON_I))))
CON_II  <- subset(CON_II,  cells = sample(Cells(CON_II),  min(3000, ncol(CON_II))))
CON_III <- subset(CON_III, cells = sample(Cells(CON_III), min(3000, ncol(CON_III))))
IPF_I   <- subset(IPF_I,   cells = sample(Cells(IPF_I),   min(3000, ncol(IPF_I))))
IPF_II  <- subset(IPF_II,  cells = sample(Cells(IPF_II),  min(3000, ncol(IPF_II))))
IPF_III <- subset(IPF_III, cells = sample(Cells(IPF_III), min(3000, ncol(IPF_III))))

# Normalize + HVG
CONIPF.list <- list(
  CON_I = CON_I, CON_II = CON_II, CON_III = CON_III,
  IPF_I = IPF_I, IPF_II = IPF_II, IPF_III = IPF_III
)
rm(CON_I, CON_II, CON_III, IPF_I, IPF_II, IPF_III)
gc()

CONIPF.list <- lapply(X = CONIPF.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Integration
features        <- SelectIntegrationFeatures(object.list = CONIPF.list)
CONIPF.anchors  <- FindIntegrationAnchors(object.list = CONIPF.list, anchor.features = features)
CONIPF.combined <- IntegrateData(anchorset = CONIPF.anchors)

CONIPF.combined <- IntegrateData(anchorset = CONIPF.anchors)

saveRDS(CONIPF.combined, file = "C:/Users/aliaa/Desktop/CONIPF_integrated.rds")

#scaling

DefaultAssay(CONIPF.combined) <- "integrated"

CONIPF.combined <- ScaleData(CONIPF.combined, verbose = FALSE)
CONIPF.combined <- RunPCA(CONIPF.combined, npcs = 50, verbose = FALSE)
ElbowPlot(CONIPF.combined, ndims = 50)

#dimensional reduction
CONIPF.combined <- RunUMAP(CONIPF.combined, reduction = "pca", dims = 1:20)
CONIPF.combined <- FindNeighbors(CONIPF.combined, reduction = "pca", dims = 1:20)
CONIPF.combined <- FindClusters(CONIPF.combined, resolution = 0.1)
DimPlot(CONIPF.combined, reduction = "umap", label = TRUE)

DimPlot(CONIPF.combined, reduction = "umap", group.by = "orig.ident")
DimPlot(CONIPF.combined, reduction = "umap", split.by = "orig.ident")
#annotation
DefaultAssay(CONIPF.combined) <- "RNA"
FeaturePlot(CONIPF.combined, features = c("EPCAM", "SFTPC", "AGER", "SCGB1A1"),
            cols = c("lightgrey", "blue"))
FeaturePlot(CONIPF.combined, features = c("CLDN5", "CCL21", "PECAM1", "EMCN"),
            label = TRUE, cols = c("lightgrey", "blue"))

FeaturePlot(CONIPF.combined, features = c("COL1A2", "LUM", "PDGFRA", "PDGFRB"),
            label = TRUE, cols = c("lightgrey", "blue"))
FeaturePlot(CONIPF.combined, features = c("PTPRC", "CD52", "AIF1", "TRBC2"),
            label = TRUE, cols = c("lightgrey", "blue"))

CONIPF.combined <- RenameIdents(CONIPF.combined,
                                `0` = "Macrophages",
                                `1` = "Macrophages",
                                `2` = "Endothelial",
                                `3` = "Epithelial",
                                `4` = "T/NK Cells",
                                `5` = "Fibroblasts",
                                `6` = "Fibroblasts",
                                `7` = "Immune",
                                `8` = "Endothelial",
                                `9` = "T Cells",
                                `10` = "Unknown",
                                `11` = "Unknown",
                                `12` = "Unknown",
                                `13` = "Unknown"
)
DimPlot(CONIPF.combined, reduction = "umap", label = TRUE)
DimPlot(CONIPF.combined, reduction = "umap", split.by = "orig.ident", label = TRUE)

saveRDS(CONIPF.combined, file = "C:/Users/aliaa/Desktop/CONIPF_annotated.rds")

#biomarkers

DefaultAssay(CONIPF.combined) <- "RNA"


CONIPF.combined <- JoinLayers(CONIPF.combined)
markers <- FindAllMarkers(CONIPF.combined, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) %>% print(n = 40)

CONIPF.combined <- RenameIdents(CONIPF.combined,
                                "Macrophages" = "Monocyte-derived Macrophages",
                                "Endothelial" = "Endothelial",
                                "Epithelial"  = "Epithelial",
                                "T/NK Cells"  = "NK Cells",
                                "Fibroblasts" = "Fibroblasts",
                                "Immune"      = "Dendritic Cells",
                                "T Cells"     = "Plasma B Cells",
                                "Unknown"     = "Mast Cells"
)
DimPlot(CONIPF.combined, reduction = "umap", label = TRUE, repel = TRUE)

saveRDS(CONIPF.combined, file = "C:/Users/aliaa/Desktop/CONIPF_annotated.rds")

#proportion comparison plot — this shows how cell type composition shifts between Normal and IPF
library(ggplot2)
library(dplyr)

# Get cell type identities into metadata
CONIPF.combined$cell_type <- Idents(CONIPF.combined)

# Calculate proportions
props <- CONIPF.combined@meta.data %>%
  group_by(orig.ident, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(prop = n / sum(n) * 100)

# Custom colors
cell_colors <- c(
  "Monocyte-derived Macrophages" = "#E63946",
  "Endothelial"                  = "#F4A261",
  "Epithelial"                   = "#2A9D8F",
  "NK Cells"                     = "#457B9D",
  "Fibroblasts"                  = "#1D3557",
  "Dendritic Cells"              = "#A8DADC",
  "Plasma B Cells"               = "#9B5DE5",
  "Mast Cells"                   = "#F15BB5"
)
# Rename conditions for cleaner labels
props$Condition <- ifelse(props$orig.ident == "CON", "Normal", "IPF")

ggplot(props, aes(x = Condition, y = prop, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.6, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = cell_colors) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101)) +
  labs(
    title = "Cell Type Composition in Normal vs IPF Lung",
    subtitle = "GSE132771 | 3 donors per condition | CCA integration",
    x = NULL,
    y = "Percentage of Cells (%)",
    fill = "Cell Type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey50", size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.text = element_text(color = "black", size = 13),
    axis.line = element_line(color = "black"),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave("C:/Users/aliaa/Desktop/cell_proportion_CON_vs_IPF.pdf", 
       width = 8, height = 6, dpi = 300)

#dotplot

DefaultAssay(CONIPF.combined) <- "RNA"

DotPlot(CONIPF.combined, 
        features = c("S100A8", "FCN1",
                     "EMCN", "PECAM1",       
                     "EPCAM", "SFTPC",        
                     "KLRD1", "PRF1",         
                     "DCN", "COL1A2",         
                     "CD1C", "FCER1A",        
                     "IGHG1", "IGKC",         
                     "TPSB2"),               
        cols = c("#f0f0f0", "#E63946"),
        dot.scale = 10,
        scale = TRUE) +
  RotatedAxis() +
  geom_point(aes(size = pct.exp), shape = 21, colour = "black", stroke = 0.3) +
  scale_size(range = c(1, 10), name = "% Expressed") +
  labs(title = "Canonical Marker Gene Expression Across Cell Types",
       subtitle = "Dot size = % cells expressing | Color = scaled average expression",
       x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey50", size = 9),
    axis.text.x   = element_text(size = 10, face = "italic", angle = 45, hjust = 1),
    axis.text.y   = element_text(size = 10),
    legend.title  = element_text(size = 9, face = "bold"),
    panel.grid.major = element_line(color = "grey92")
  )


#featureplot
p1 <- FeaturePlot(CONIPF.combined, 
                  features = "COL1A2",
                  split.by = "orig.ident",
                  min.cutoff = "q10",
                  cols = c("lightgrey", "#E63946"),
                  pt.size = 0.3) &
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_blank())
p2 <- FeaturePlot(CONIPF.combined, 
                  features = "DCN",
                  split.by = "orig.ident",
                  min.cutoff = "q10",
                  cols = c("lightgrey", "#1D3557"),
                  pt.size = 0.3) &
  theme(plot.title = element_text(face = "bold"),
        axis.title = element_blank())

p1 / p2 + 
  plot_annotation(
    title = "Fibroblast Marker Expression: Normal vs IPF",
    subtitle = "COL1A2 and DCN — min.cutoff q10 to enhance contrast",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "grey50", size = 10)
    )
  )

#vlnplot

VlnPlot(CONIPF.combined,
        features = c("COL1A2"),
        idents = "Fibroblasts",
        split.by = "orig.ident",
        split.plot = TRUE,
        cols = c("#457B9D", "#E63946"),
        pt.size = 0) +
  labs(title = "COL1A2 Expression in Fibroblasts: Normal vs IPF")

#visual

DefaultAssay(CONIPF.combined) <- "RNA"

p <- FeaturePlot(CONIPF.combined,
                 features = c("PTPRC", "COL1A2", "CLDN5", "EPCAM"),
                 ncol = 2,
                 min.cutoff = "q10",
                 cols = c("lightgrey", "#E63946"),
                 pt.size = 0.3,
                 label = TRUE,
                 label.size = 3,
                 repel = TRUE) &
  theme(
    plot.title   = element_text(face = "italic", size = 11, color = "#1D3557"),
    axis.title   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA)
  )

p + plot_annotation(
  title    = "Lineage Marker Expression Across Cell Types",
  subtitle = "PTPRC = immune | COL1A2 = mesenchymal | CLDN5 = endothelial | EPCAM = epithelial",
  theme = theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey50", size = 9)
  )
)

ggsave("C:/Users/aliaa/Desktop/featureplot_markers.pdf", 
       width = 10, height = 9, dpi = 300)


library(ggrepel)
DefaultAssay(CONIPF.combined) <- "RNA"


Idents(CONIPF.combined) <- "cell_type"

fib_de <- FindMarkers(CONIPF.combined,
                      ident.1 = "IPF",
                      ident.2 = "CON",
                      group.by = "orig.ident",
                      subset.ident = "Fibroblasts",
                      min.pct = 0.25,
                      logfc.threshold = 0.25)


fib_de$gene <- rownames(fib_de)
fib_de$significance <- "Not Significant"

fib_de$significance[fib_de$avg_log2FC > 0.5 & fib_de$p_val_adj < 0.05] <- "Up in IPF"
fib_de$significance[fib_de$avg_log2FC < -0.5 & fib_de$p_val_adj < 0.05] <- "Up in Normal"

top_up   <- fib_de %>% filter(significance == "Up in IPF") %>% top_n(10, avg_log2FC)

top_down <- fib_de %>% filter(significance == "Up in Normal") %>% top_n(10, -avg_log2FC)
top_genes <- rbind(top_up, top_down)

ggplot(fib_de, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_point(data = top_genes, size = 2.5, alpha = 1) +
  geom_text_repel(data = top_genes,
                  aes(label = gene),
                  size = 3,
                  fontface = "italic",
                  max.overlaps = 20,
                  box.padding = 0.5) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
  scale_color_manual(values = c(
    "Up in IPF"       = "#E63946",
    "Up in Normal"    = "#457B9D",
    "Not Significant" = "grey70"
  )) +
  scale_x_continuous(limits = c(-3, 3)) +
  labs(
    title    = "Differential Gene Expression in Fibroblasts: IPF vs Normal",
    subtitle = "Genes with |log2FC| > 0.5 and adj.p < 0.05 are highlighted",
    x        = "Average log2 Fold Change (IPF vs Normal)",
    y        = "-log10 Adjusted P-value",
    color    = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey50", size = 9),
    legend.position = "top",
    axis.text     = element_text(color = "black")
  )

fib_de %>% 
  filter(significance == "Up in IPF") %>% 
  arrange(desc(avg_log2FC)) %>% 
  head(15) %>% 
  select(gene, avg_log2FC, p_val_adj)

p <- FeaturePlot(CONIPF.combined,
                 features = "CTHRC1",
                 split.by = "orig.ident",
                 min.cutoff = "q10",
                 cols = c("lightgrey", "#E63946"),
                 pt.size = 0.5) &
  theme(plot.title = element_text(face = "bold.italic"),
        axis.title = element_blank())
p + plot_annotation(
  title = "CTHRC1 Expression: Normal vs IPF",
  subtitle = "Pathological myofibroblast marker — exclusively elevated in IPF fibroblasts",
  theme = theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(color = "grey50", size = 9)
  )
)

ggsave("C:/Users/aliaa/Desktop/CTHRC1_featureplot.pdf",
       width = 10, height = 5, dpi = 300)

install.packages("enrichR")
library(enrichR)

# Set database
dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human")


ipf_genes <- fib_de %>%
  filter(significance == "Up in IPF") %>%
  arrange(desc(avg_log2FC)) %>%
  pull(gene)
# Run enrichment
enrich_results <- enrichr(ipf_genes, dbs)

# Plot GO
plotEnrich(enrich_results[["GO_Biological_Process_2023"]], 
           showTerms = 15,
           numChar = 50,
           title = "GO Biological Process — IPF Fibroblast Upregulated Genes") +
  theme(plot.title = element_text(face = "bold", size = 12),
        axis.text.y = element_text(size = 9))

ggsave("C:/Users/aliaa/Desktop/GO_enrichment.pdf", width = 11, height = 7, dpi = 300)

ipf_genes_clean <- ipf_genes[!grepl("^IGH|^IGK|^IGL|^IGHA|^IGHG", ipf_genes)]
enrich_results_clean <- enrichr(ipf_genes_clean, dbs)
go_results <- enrich_results_clean[["GO_Biological_Process_2023"]]

go_df <- go_results %>%
  head(12) %>%
  mutate(Term = gsub(" \\(GO:.*\\)", "", Term),
         Term = str_to_title(Term),
         GeneCount = as.numeric(sapply(strsplit(Overlap, "/"), `[`, 1)))

ggplot(go_df, aes(x = reorder(Term, GeneCount),
                  y = GeneCount,
                  fill = P.value)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_gradient(low = "#E63946", high = "#A8DADC",
                      name = "P value") +
  labs(title    = "GO Biological Process Enrichment",
       subtitle = "Upregulated genes in IPF vs Normal Fibroblasts",
       x = NULL,
       y = "Number of Genes") +
  theme_classic(base_size = 12) +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "grey50", size = 9),
        axis.text.y   = element_text(size = 10))

file.exists("C:/Users/aliaa/Desktop/CONIPF_integrated.rds")

file.exists("C:/Users/aliaa/Desktop/CONIPF_annotated.rds")