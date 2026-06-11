##  ─────────────────────────────────
##  1. run_cellchat_inference() contained hardcoded "cellchat_CON" references
##     instead of the generic argument "cc" — function now fully generic.
##  2. identifyOverExpressedInteractions() was missing from the function body.
##  3. pak::pak() call was placed mid-analysis — moved to Section 0.
##  4. orig.ident subsetting: confirmed CON = "CON", IPF = "IPF" from
##     CreateSeuratObject(project = ...) calls in IPF.R.
##  5. getMaxWeight() for shared pathways was passed a character vector —
##     corrected to loop per pathway individually.
##  6. rankNet(return.data = TRUE) — not available in all CellChat versions;
##     added safe fallback via tryCatch.
##  7. group.cellType lineage block now fully populated with your cell types.
##
################################################################################


## ══════════════════════════════════════════════════════════════════════════════
## 0.  ENVIRONMENT SETUP
## ══════════════════════════════════════════════════════════════════════════════

## -- Install (run once, then comment out) -------------------------------------
# install.packages("pak")
# pak::pak("jinworks/CellChat")          # CellChat v2
# pak::pak("immunogenomics/presto")      # fast Wilcoxon (required by CellChat)
# pak::pak("NMF")
# pak::pak("ggalluvial")

## -- Load libraries -----------------------------------------------------------
library(CellChat)
library(Seurat)
library(ggplot2)
library(patchwork)
library(NMF)             # communication pattern analysis
library(ggalluvial)      # river (alluvial) plots
library(circlize)        # chord diagrams
library(ComplexHeatmap)  # comparative heatmaps
library(dplyr)

options(stringsAsFactors = FALSE)
set.seed(42)

## -- Output directory ---------------------------------------------------------
out_dir <- "C:/Users/aliaa/Desktop/CellChat_Results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("Output directory:", out_dir, "\n")


## ══════════════════════════════════════════════════════════════════════════════
## 1.  LOAD & INSPECT SEURAT OBJECT
## ══════════════════════════════════════════════════════════════════════════════

seurat_obj <- readRDS("C:/Users/aliaa/Desktop/CONIPF_annotated.rds")

## Confirm metadata structure
cat("=== Metadata columns ===\n");      print(colnames(seurat_obj@meta.data))
cat("\n=== Active identity table ===\n"); print(table(Idents(seurat_obj)))
cat("\n=== orig.ident (samples) ===\n"); print(table(seurat_obj@meta.data$orig.ident))
cat("\n=== Cell types × condition ===\n")
print(table(Idents(seurat_obj), seurat_obj@meta.data$orig.ident))

# ── Step 1: What do the active Idents actually look like?
print(head(Idents(seurat_obj)))       # are these names or numbers?
print(levels(Idents(seurat_obj)))     # all unique values


# ── Step 2: Does the UMAP reduction exist?
print(names(seurat_obj@reductions))   # must contain "umap"

# ── Step 3: What metadata columns exist?
print(colnames(seurat_obj@meta.data))

# ── Step 4: What does the cell_type column contain?
print(table(seurat_obj@meta.data$cell_type))

# Step 1: confirm Idents are the cell type names (already confirmed above)
# Add them permanently into metadata so group.by works reliably
seurat_obj$cell_type <- as.character(Idents(seurat_obj))

# Verify it now exists
print(table(seurat_obj$cell_type))  # should show 8 cell types with counts

# Step 2: absolute minimum plot — no colors, no options
# If THIS is blank, there is a deeper Seurat v5 layer issue
DimPlot(seurat_obj)
my_cols <- c(
  "Monocyte-derived Macrophages" = "#E63946",
  "Endothelial"                  = "#F4A261",
  "Epithelial"                   = "#2A9D8F",
  "NK Cells"                     = "#457B9D",
  "Fibroblasts"                  = "#1D3557",
  "Dendritic Cells"              = "#A8DADC",
  "Plasma B Cells"               = "#9B5DE5",
  "Mast Cells"                   = "#F15BB5"
)

p_umap <- DimPlot(seurat_obj,
                  group.by = "cell_type",   # now exists in metadata
                  label    = TRUE,
                  repel    = TRUE) +
  scale_color_manual(values = my_cols) +    # safer than cols= for named vectors
  ggtitle("Cell types — full dataset (CON + IPF)") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

print(p_umap)

# ── Save cell_type to metadata permanently (from last fix)
seurat_obj$cell_type <- as.character(Idents(seurat_obj))

# ── Split by condition
seurat_CON <- subset(seurat_obj, subset = orig.ident == "CON")
seurat_IPF <- subset(seurat_obj, subset = orig.ident == "IPF")

cat("CON cells:", ncol(seurat_CON), "\n")
cat("IPF cells:", ncol(seurat_IPF), "\n")

# Sanity check — confirm cell types survived the split
cat("\nCON cell types:\n"); print(table(Idents(seurat_CON)))
cat("\nIPF cell types:\n"); print(table(Idents(seurat_IPF)))

# ── Create CellChat objects
cellchat_CON <- createCellChat(
  object   = seurat_CON,
  group.by = "cell_type",
  assay    = "RNA"
)

cellchat_IPF <- createCellChat(
  object   = seurat_IPF,
  group.by = "cell_type",
  assay    = "RNA"
)

print(cellchat_CON)
print(cellchat_IPF)

# ── Assign full human LR database to both
CellChatDB.use        <- CellChatDB.human
cellchat_CON@DB <- CellChatDB.use
cellchat_IPF@DB <- CellChatDB.use

cat("Database assigned.\n")
cat("Categories:", unique(CellChatDB.human$interaction$annotation), "\n")


library(future)

run_cellchat_inference <- function(cc, label) {
  
  cat("\n===", label, "- starting inference ===\n")
  
  cc <- subsetData(cc)
  cat("  subsetData done\n")
  
  plan("multisession", workers = 4)
  cc <- identifyOverExpressedGenes(cc)
  plan("sequential")
  cat("  identifyOverExpressedGenes done\n")
  
  cc <- identifyOverExpressedInteractions(cc)
  cat("  identifyOverExpressedInteractions done\n")
  
  # projectData removed in CellChat v2 — use raw.use = TRUE instead
  cc <- computeCommunProb(cc, type = "triMean", raw.use = TRUE)
  cat("  computeCommunProb done\n")
  
  cc <- filterCommunication(cc, min.cells = 10)
  cat("  filterCommunication done\n")
  
  cc <- computeCommunProbPathway(cc)
  cat("  computeCommunProbPathway done\n")
  
  cc <- aggregateNet(cc)
  cat("  aggregateNet done\n")
  
  cat("\n  Active pathways:", length(cc@netP$pathways), "\n")
  cat("  Interaction matrix:", paste(dim(cc@net$count), collapse = " x "), "\n\n")
  
  return(cc)
}

seurat_obj <- readRDS("C:/Users/aliaa/Desktop/CONIPF_annotated.rds")
seurat_obj$cell_type <- as.character(Idents(seurat_obj))

# Split
seurat_CON <- subset(seurat_obj, subset = orig.ident == "CON")
seurat_IPF <- subset(seurat_obj, subset = orig.ident == "IPF")

# Recreate CellChat objects
cellchat_CON <- createCellChat(object = seurat_CON, group.by = "cell_type", assay = "RNA")
cellchat_IPF <- createCellChat(object = seurat_IPF, group.by = "cell_type", assay = "RNA")

# Assign database
cellchat_CON@DB <- CellChatDB.human
cellchat_IPF@DB <- CellChatDB.human

cat("Ready.\n")

# Now run inference
cellchat_CON <- run_cellchat_inference(cellchat_CON, "CON")
cellchat_IPF <- run_cellchat_inference(cellchat_IPF, "IPF")

saveRDS(cellchat_CON, file.path(out_dir, "cellchat_CON_inferred.rds"))
saveRDS(cellchat_IPF, file.path(out_dir, "cellchat_IPF_inferred.rds"))
cat("Both saved.\n")
cellchat_CON <- netAnalysis_computeCentrality(cellchat_CON, slot.name = "netP")
cellchat_IPF <- netAnalysis_computeCentrality(cellchat_IPF, slot.name = "netP")
cat("Centrality done.\n")


groupSize_CON <- as.numeric(table(cellchat_CON@idents))
groupSize_IPF <- as.numeric(table(cellchat_IPF@idents))

pdf(file.path(out_dir, "03_aggregate_network_CON.pdf"), width = 11, height = 5.5)
par(mfrow = c(1, 2), xpd = TRUE)

netVisual_circle(cellchat_CON@net$count,
                 vertex.weight = groupSize_CON,
                 weight.scale  = TRUE,
                 label.edge    = FALSE,
                 title.name    = "CON — Number of interactions")
netVisual_circle(cellchat_CON@net$weight,
                 vertex.weight = groupSize_CON,
                 weight.scale  = TRUE,
                 label.edge    = FALSE,
                 title.name    = "CON — Interaction strength")
dev.off()

pdf(file.path(out_dir, "03_aggregate_network_IPF.pdf"), width = 11, height = 5.5)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat_IPF@net$count,
                 vertex.weight = groupSize_IPF,
                 weight.scale  = TRUE,
                 label.edge    = FALSE,
                 title.name    = "IPF — Number of interactions")
netVisual_circle(cellchat_IPF@net$weight,
                 vertex.weight = groupSize_IPF,
                 weight.scale  = TRUE,
                 label.edge    = FALSE,
                 title.name    = "IPF — Interaction strength")
dev.off()
cat("Circle plots saved.\n")

par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat_CON@net$count,
                 vertex.weight = groupSize_CON,
                 weight.scale  = TRUE,
                 label.edge    = FALSE,
                 title.name    = "CON - Number of interactions")
netVisual_circle(cellchat_CON@net$weight,
                 vertex.weight = groupSize_CON,
                 weight.scale  = TRUE,
                 label.edge    = FALSE,
                 title.name    = "CON - Interaction strength")

# IPF
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat_IPF@net$count,
                 vertex.weight = groupSize_IPF,
                 weight.scale  = TRUE,
                 label.edge    = FALSE,
                 title.name    = "IPF - Number of interactions")
netVisual_circle(cellchat_IPF@net$weight,
                 vertex.weight = groupSize_IPF,
                 weight.scale  = TRUE,
                 label.edge    = FALSE,
                 title.name    = "IPF - Interaction strength")

pathways_CON <- cellchat_CON@netP$pathways
pathways_IPF <- cellchat_IPF@netP$pathways
shared     <- intersect(pathways_CON, pathways_IPF)
gained_IPF <- setdiff(pathways_IPF, pathways_CON)
lost_IPF   <- setdiff(pathways_CON, pathways_IPF)

cat("CON pathways  :", length(pathways_CON), "\n")
cat("IPF pathways  :", length(pathways_IPF), "\n")
cat("Shared        :", length(shared), "\n")
cat("Gained in IPF :", gained_IPF, "\n")
cat("Lost in IPF   :", lost_IPF, "\n")

# Merge for comparative analysis
object.list     <- list(CON = cellchat_CON, IPF = cellchat_IPF)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))
cat("Merged.\n")

# Overall interaction comparison
p_count  <- compareInteractions(cellchat_merged, show.legend = FALSE,
                                group = c(1,2), measure = "count")
p_weight <- compareInteractions(cellchat_merged, show.legend = FALSE,
                                group = c(1,2), measure = "weight")

p_count + p_weight

# Information flow — ranks all pathways by strength CON vs IPF
rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)

rankNet(cellchat_merged, mode = "comparison", stacked = FALSE, do.stat = TRUE)

# Differential circle plots (red = gained in IPF, blue = lost)
par(mfrow = c(1,2), xpd = TRUE)
netVisual_diffInteraction(cellchat_merged, comparison = c(1,2),
                          measure = "count", weight.scale = TRUE,
                          title.name = "Differential interactions (count)")
netVisual_diffInteraction(cellchat_merged, comparison = c(1,2),
                          measure = "weight", weight.scale = TRUE,
                          title.name = "Differential interactions (strength)")

# Signaling role scatter — who are the dominant senders/receivers?
par(mfrow = c(1,2))
netAnalysis_signalingRole_scatter(cellchat_CON, title = "CON — Signaling Roles")
netAnalysis_signalingRole_scatter(cellchat_IPF, title = "IPF — Signaling Roles")

# Extract the bubble plot data directly
bubble_data <- netVisual_bubble(cellchat_merged,
                                sources.use = c("Monocyte-derived Macrophages",
                                                "Fibroblasts"),
                                targets.use = c("Fibroblasts", "Epithelial",
                                                "Endothelial"),
                                comparison = c(1,2),
                                thresh = 0.01,
                                remove.isolate = TRUE,
                                max.dataset = 2,
                                min.dataset = 1,
                                return.data = TRUE)

# Check what's inside
df <- bubble_data$communication
head(df)
nrow(df)

library(ggplot2)
library(dplyr)


df_plot <- df %>%
  mutate(
    # Short readable label: Ligand → Receptor
    lr_label = paste0(ligand, " → ", receptor),
    # Clean source-target label
    axis_label = paste0(source, "\n→ ", target),
    # -log10 pval for size
    sig_size = ifelse(pval == 1, 1, ifelse(pval == 2, 2, 3)),
    # Dataset as factor
    dataset = factor(dataset, levels = c("CON", "IPF"))
  )

top_lr <- df_plot %>%
  group_by(lr_label) %>%
  summarise(mean_prob = mean(prob.original)) %>%
  arrange(desc(mean_prob)) %>%
  slice_head(n = 25) %>%
  pull(lr_label)

df_top <- df_plot %>% filter(lr_label %in% top_lr)

# Order LR pairs by mean probability
lr_order <- df_top %>%
  group_by(lr_label) %>%
  summarise(mean_prob = mean(prob.original)) %>%
  arrange(mean_prob) %>%
  pull(lr_label)

df_top$lr_label <- factor(df_top$lr_label, levels = lr_order)
df_top$group.names <- factor(df_top$group.names,
                             levels = sort(unique(df_top$group.names)))

p <- ggplot(df_top, aes(x = group.names, y = lr_label)) +
  
  # Bubbles
  geom_point(aes(size  = prob.original,
                 color = prob.original,
                 shape = dataset),
             alpha = 0.85) +
  
  # Scales
  scale_color_gradientn(
    colors = c("#313695","#4575b4","#74add1","#abd9e9",
               "#fee090","#fdae61","#f46d43","#d73027","#a50026"),
    name = "Communication\nProbability"
  ) +
  scale_size_continuous(range = c(1.5, 8), name = "Probability") +
  scale_shape_manual(values = c("CON" = 16, "IPF" = 17),
                     name = "Condition") +
  
  # Facet by condition
  facet_grid(. ~ dataset, scales = "free_x", space = "free") +
  
  # Labels
  labs(
    title    = "Top Ligand-Receptor Interactions: CON vs IPF",
    subtitle = "Fibroblasts & Macrophages → Target Cells | p < 0.01",
    x = "Cell-Cell Interaction",
    y = "Ligand → Receptor"
  ) +
  
  # Theme
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.x   = element_text(angle = 40, hjust = 1, size = 9, color = "black"),
    axis.text.y   = element_text(size = 9.5, color = "black"),
    axis.title    = element_text(face = "bold", size = 11),
    strip.text    = element_text(face = "bold", size = 12,
                                 color = "white"),
    strip.background = element_rect(fill = c("#2166ac"), color = NA),
    legend.title  = element_text(face = "bold", size = 10),
    legend.text   = element_text(size = 9),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.4),
    panel.grid.major.x = element_line(color = "grey93", linewidth = 0.4),
    panel.spacing = unit(0.8, "lines")
  ) +
  guides(
    color = guide_colorbar(barwidth = 1, barheight = 6),
    size  = "none"
  )

# ── Save ──────────────────────────────────────────────────────────────────────
ggsave(file.path(out_dir, "04_bubble_publication.pdf"),
       p, width = 14, height = 10, dpi = 300)
ggsave(file.path(out_dir, "04_bubble_publication.png"),
       p, width = 14, height = 10, dpi = 300)
print(p)

# Extract outgoing signaling data for both conditions
library(reshape2)
library(ggplot2)
library(dplyr)

# Get the signaling role matrices
mat_CON <- cellchat_CON@netP$prob
mat_IPF <- cellchat_IPF@netP$prob

# Compute outgoing strength per cell type per pathway
outgoing_CON <- apply(mat_CON, 3, function(x) rowSums(x))
outgoing_IPF <- apply(mat_IPF, 3, function(x) rowSums(x))

# Normalize
outgoing_CON_norm <- t(apply(outgoing_CON, 1, function(x) x / max(x, na.rm=TRUE)))
outgoing_IPF_norm <- t(apply(outgoing_IPF, 1, function(x) x / max(x, na.rm=TRUE)))

# Convert to long format
df_CON <- melt(outgoing_CON_norm) %>%
  rename(CellType = Var1, Pathway = Var2, Strength = value) %>%
  mutate(Condition = "CON")

df_IPF <- melt(outgoing_IPF_norm) %>%
  rename(CellType = Var1, Pathway = Var2, Strength = value) %>%
  mutate(Condition = "IPF")

# Combine
df_all <- rbind(df_CON, df_IPF)

# Focus on top pathways by mean strength in IPF
top_pathways <- df_IPF %>%
  group_by(Pathway) %>%
  summarise(mean_str = mean(Strength, na.rm=TRUE)) %>%
  arrange(desc(mean_str)) %>%
  slice_head(n = 20) %>%
  pull(Pathway)

df_plot <- df_all %>%
  filter(Pathway %in% top_pathways) %>%
  mutate(
    Condition = factor(Condition, levels = c("CON", "IPF")),
    Pathway = factor(Pathway, levels = top_pathways),
    CellType = gsub("Monocyte-derived Macrophages", "Macrophages", CellType)
  )

# ── Build publication heatmap ─────────────────────────────────────────────────
p_heat <- ggplot(df_plot, aes(x = CellType, y = Pathway, fill = Strength)) +
  geom_tile(color = "white", linewidth = 0.4) +
  facet_grid(. ~ Condition) +
  scale_fill_gradientn(
    colors = c("white", "#fee8c8", "#fdbb84", "#e34a33", "#99000d"),
    na.value = "white",
    name = "Relative\nStrength",
    limits = c(0, 1)
  ) +
  labs(
    title    = "Outgoing Signaling Strength: CON vs IPF",
    subtitle = "Top 20 pathways | Normalized per pathway",
    x = NULL, y = "Signaling Pathway"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.x   = element_text(angle = 40, hjust = 1, size = 10,
                                 color = "black", face = "bold"),
    axis.text.y   = element_text(size = 10, color = "black"),
    axis.title.y  = element_text(face = "bold", size = 11),
    strip.text    = element_text(face = "bold", size = 13, color = "white"),
    strip.background = element_rect(fill = "#2c3e50", color = NA),
    legend.title  = element_text(face = "bold", size = 10),
    legend.text   = element_text(size = 9),
    panel.spacing = unit(1, "lines"),
    panel.border  = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )

ggsave(file.path(out_dir, "05_outgoing_heatmap_publication.pdf"),
       p_heat, width = 12, height = 8, dpi = 300)
ggsave(file.path(out_dir, "05_outgoing_heatmap_publication.png"),
       p_heat, width = 12, height = 8, dpi = 300)
print(p_heat)


# ── Differential heatmap: IPF minus CON ──────────────────────────────────────
# Get shared pathways only
shared_paths <- intersect(colnames(outgoing_CON_norm), 
                          colnames(outgoing_IPF_norm))

# Align cell types
shared_cells <- intersect(rownames(outgoing_CON_norm), 
                          rownames(outgoing_IPF_norm))

diff_mat <- outgoing_IPF_norm[shared_cells, shared_paths] - 
  outgoing_CON_norm[shared_cells, shared_paths]

# Focus on most changed pathways
top_diff <- names(sort(apply(abs(diff_mat), 2, max), decreasing = TRUE))[1:25]

df_diff <- melt(diff_mat[, top_diff]) %>%
  rename(CellType = Var1, Pathway = Var2, Delta = value) %>%
  mutate(
    CellType = gsub("Monocyte-derived Macrophages", "Macrophages", 
                    as.character(CellType)),
    Direction = ifelse(Delta > 0, "Increased in IPF", "Decreased in IPF")
  )

# Order pathways by mean delta
path_order <- df_diff %>%
  group_by(Pathway) %>%
  summarise(mean_delta = mean(Delta)) %>%
  arrange(mean_delta) %>%
  pull(Pathway)

df_diff$Pathway  <- factor(df_diff$Pathway, levels = path_order)

p_diff <- ggplot(df_diff, aes(x = CellType, y = Pathway, fill = Delta)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(data = df_diff %>% filter(abs(Delta) > 0.3),
            aes(label = round(Delta, 2)),
            size = 2.8, color = "white", fontface = "bold") +
  scale_fill_gradientn(
    colors  = c("#2166ac","#4393c3","#92c5de","#f7f7f7",
                "#f4a582","#d6604d","#b2182b"),
    limits  = c(-1, 1),
    na.value = "grey90",
    name    = "Δ Strength\n(IPF − CON)"
  ) +
  labs(
    title    = "Differential Outgoing Signaling: IPF vs Control",
    subtitle = "Red = increased in IPF | Blue = decreased in IPF",
    x = NULL, y = "Signaling Pathway"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.x   = element_text(angle = 40, hjust = 1, size = 10,
                                 color = "black", face = "bold"),
    axis.text.y   = element_text(size = 10, color = "black"),
    axis.title.y  = element_text(face = "bold", size = 11),
    legend.title  = element_text(face = "bold", size = 10),
    panel.border  = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
  )

ggsave(file.path(out_dir, "06_differential_outgoing_heatmap.pdf"),
       p_diff, width = 11, height = 9, dpi = 300)
ggsave(file.path(out_dir, "06_differential_outgoing_heatmap.png"),
       p_diff, width = 11, height = 9, dpi = 300)
print(p_diff)


# ── Signaling role scatter: who sends vs receives ─────────────────────────────
gg1 <- netAnalysis_signalingRole_scatter(cellchat_CON) +
  ggtitle("Control") +
  theme_classic(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title  = element_text(face = "bold", size = 11),
    axis.text   = element_text(color = "black", size = 10),
    legend.position = "bottom"
  )

gg2 <- netAnalysis_signalingRole_scatter(cellchat_IPF) +
  ggtitle("IPF") +
  theme_classic(base_size = 12) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title  = element_text(face = "bold", size = 11),
    axis.text   = element_text(color = "black", size = 10),
    legend.position = "bottom"
  )

p_scatter <- (gg1 | gg2) +
  plot_annotation(
    title = "Signaling Roles of Cell Types: CON vs IPF",
    subtitle = "X-axis: outgoing strength | Y-axis: incoming strength",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
    )
  )

ggsave(file.path(out_dir, "07_signaling_roles_scatter.pdf"),
       p_scatter, width = 13, height = 6, dpi = 300)
ggsave(file.path(out_dir, "07_signaling_roles_scatter.png"),
       p_scatter, width = 13, height = 6, dpi = 300)
print(p_scatter)


# Use circle layout instead of chord - more stable
key_pathways <- c("COLLAGEN", "MIF", "CXCL", "CCL")

for(pw in key_pathways) {
  pdf(file.path(out_dir, paste0("08_circle_", pw, ".pdf")), 
      width = 12, height = 6)
  par(mfrow = c(1,2), xpd = TRUE, mar = c(3,3,5,3))
  
  if(pw %in% pathways_CON) {
    netVisual_aggregate(cellchat_CON, signaling = pw, layout = "circle")
    title(paste("CON:", pw), line = 0.5, cex.main = 1.4, font.main = 2)
  } else {
    plot.new()
    text(0.5, 0.5, paste("CON:\nNo significant\n", pw, "signaling"),
         cex = 1.3, font = 2, col = "grey50")
  }
  
  if(pw %in% pathways_IPF) {
    netVisual_aggregate(cellchat_IPF, signaling = pw, layout = "circle")
    title(paste("IPF:", pw), line = 0.5, cex.main = 1.4, font.main = 2)
  } else {
    plot.new()
    text(0.5, 0.5, paste("IPF:\nNo significant\n", pw, "signaling"),
         cex = 1.3, font = 2, col = "grey50")
  }
  
  dev.off()
  cat("Saved:", pw, "\n")
}

# View COLLAGEN in RStudio
par(mfrow = c(1,2), xpd = TRUE, mar = c(3,3,5,3))
netVisual_aggregate(cellchat_CON, signaling = "COLLAGEN", layout = "circle")
title("CON: COLLAGEN Signaling", line = 0.5, cex.main = 1.4, font.main = 2)
netVisual_aggregate(cellchat_IPF, signaling = "COLLAGEN", layout = "circle")
title("IPF: COLLAGEN Signaling", line = 0.5, cex.main = 1.4, font.main = 2)



# ── Extract pathway-level communication data ──────────────────────────────────
library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)

# Helper: extract edge list for a pathway from a cellchat object
get_pathway_edges <- function(cc, pathway) {
  if (!pathway %in% cc@netP$pathways) return(NULL)
  prob <- cc@netP$prob[,,pathway]
  df <- melt(prob) %>%
    rename(Source = Var1, Target = Var2, Probability = value) %>%
    filter(Probability > 0) %>%
    mutate(
      Source = gsub("Monocyte-derived Macrophages", "Macrophages", as.character(Source)),
      Target = gsub("Monocyte-derived Macrophages", "Macrophages", as.character(Target))
    )
  return(df)
}

# ── Get data for all pathways ─────────────────────────────────────────────────
mif_con  <- get_pathway_edges(cellchat_CON, "MIF")  %>% mutate(Condition = "CON", Pathway = "MIF")
mif_ipf  <- get_pathway_edges(cellchat_IPF, "MIF")  %>% mutate(Condition = "IPF", Pathway = "MIF")
cxcl_ipf <- get_pathway_edges(cellchat_IPF, "CXCL") %>% mutate(Condition = "IPF", Pathway = "CXCL")
spp1_ipf <- get_pathway_edges(cellchat_IPF, "SPP1") %>% mutate(Condition = "IPF", Pathway = "SPP1")

# ── Figure A: MIF CON vs IPF heatmap ─────────────────────────────────────────
mif_all <- rbind(mif_con, mif_ipf) %>%
  mutate(Condition = factor(Condition, levels = c("CON", "IPF")))

# Get all cell types for consistent axis
all_cells <- union(unique(mif_all$Source), unique(mif_all$Target))

# Complete grid
mif_complete <- mif_all %>%
  tidyr::complete(Source, Target, Condition, fill = list(Probability = 0)) %>%
  filter(Source != Target)

p_mif <- ggplot(mif_complete, 
                aes(x = Target, y = Source, fill = Probability)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(data = mif_complete %>% filter(Probability > 0.005),
            aes(label = sprintf("%.3f", Probability)),
            size = 2.5, color = "white", fontface = "bold") +
  facet_wrap(~ Condition) +
  scale_fill_gradientn(
    colors = c("white","#fee8c8","#fdbb84","#e34a33","#99000d"),
    name = "Communication\nProbability",
    limits = c(0, max(mif_complete$Probability))
  ) +
  labs(
    title    = "MIF Signaling Pathway: CON vs IPF",
    subtitle = "Macrophage inhibitory factor — inflammatory & fibrotic regulator",
    x = "Target Cell Type", y = "Source Cell Type"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40",
                                 face = "italic"),
    axis.text.x   = element_text(angle = 40, hjust = 1, size = 9,
                                 color = "black", face = "bold"),
    axis.text.y   = element_text(size = 9, color = "black", face = "bold"),
    axis.title    = element_text(face = "bold", size = 10),
    strip.text    = element_text(face = "bold", size = 12, color = "white"),
    strip.background = element_rect(fill = "#1d3557", color = NA),
    legend.title  = element_text(face = "bold", size = 9),
    panel.border  = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
    panel.spacing = unit(1, "lines")
  )

# ── Figure B: CXCL IPF heatmap ───────────────────────────────────────────────
cxcl_complete <- cxcl_ipf %>%
  tidyr::complete(Source, Target, fill = list(Probability = 0)) %>%
  filter(Source != Target)

p_cxcl <- ggplot(cxcl_complete,
                 aes(x = Target, y = Source, fill = Probability)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(data = cxcl_complete %>% filter(Probability > 0.005),
            aes(label = sprintf("%.3f", Probability)),
            size = 2.8, color = "white", fontface = "bold") +
  scale_fill_gradientn(
    colors = c("white","#d0e8f5","#74b9e7","#1a6faf","#08306b"),
    name = "Communication\nProbability"
  ) +
  labs(
    title    = "CXCL Signaling Pathway: IPF Only",
    subtitle = "Chemokine signaling — absent in control, gained in IPF",
    x = "Target Cell Type", y = "Source Cell Type"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40",
                                 face = "italic"),
    axis.text.x   = element_text(angle = 40, hjust = 1, size = 9,
                                 color = "black", face = "bold"),
    axis.text.y   = element_text(size = 9, color = "black", face = "bold"),
    axis.title    = element_text(face = "bold", size = 10),
    strip.text    = element_text(face = "bold", size = 12, color = "white"),
    strip.background = element_rect(fill = "#457b9d", color = NA),
    legend.title  = element_text(face = "bold", size = 9),
    panel.border  = element_rect(color = "grey70", fill = NA, linewidth = 0.5)
  )

# ── Figure C: SPP1 IPF heatmap ───────────────────────────────────────────────
spp1_complete <- spp1_ipf %>%
  tidyr::complete(Source, Target, fill = list(Probability = 0)) %>%
  filter(Source != Target)

p_spp1 <- ggplot(spp1_complete,
                 aes(x = Target, y = Source, fill = Probability)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(data = spp1_complete %>% filter(Probability > 0.001),
            aes(label = sprintf("%.3f", Probability)),
            size = 2.8, color = "white", fontface = "bold") +
  scale_fill_gradientn(
    colors = c("white","#e8d5f5","#c39bd3","#7d3c98","#4a235a"),
    name = "Communication\nProbability"
  ) +
  labs(
    title    = "SPP1 Signaling Pathway: IPF Only",
    subtitle = "Osteopontin — macrophage-driven fibrotic signal, absent in control",
    x = "Target Cell Type", y = "Source Cell Type"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40",
                                 face = "italic"),
    axis.text.x   = element_text(angle = 40, hjust = 1, size = 9,
                                 color = "black", face = "bold"),
    axis.text.y   = element_text(size = 9, color = "black", face = "bold"),
    axis.title    = element_text(face = "bold", size = 10),
    strip.text    = element_text(face = "bold", size = 12, color = "white"),
    strip.background = element_rect(fill = "#6b2fa0", color = NA),
    legend.title  = element_text(face = "bold", size = 9),
    panel.border  = element_rect(color = "grey70", fill = NA, linewidth = 0.5)
  )

# ── Save individually ─────────────────────────────────────────────────────────
ggsave(file.path(out_dir, "08_MIF_heatmap.pdf"),
       p_mif, width = 11, height = 5, dpi = 300)
ggsave(file.path(out_dir, "08_CXCL_heatmap.pdf"),
       p_cxcl, width = 7, height = 5, dpi = 300)
ggsave(file.path(out_dir, "08_SPP1_heatmap.pdf"),
       p_spp1, width = 7, height = 5, dpi = 300)

# ── Save as combined panel ────────────────────────────────────────────────────
p_combined <- (p_mif) / (p_cxcl | p_spp1) +
  plot_annotation(
    title = "Key Signaling Pathways in IPF vs Control",
    subtitle = "MIF (shared, rewired) | CXCL & SPP1 (IPF-specific, gained)",
    tag_levels = "A",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
    )
  )

ggsave(file.path(out_dir, "08_pathway_panel.pdf"),
       p_combined, width = 13, height = 11, dpi = 300)
ggsave(file.path(out_dir, "08_pathway_panel.png"),
       p_combined, width = 13, height = 11, dpi = 300)

print(p_mif)
print(p_cxcl)
print(p_spp1)




# ── NMF Communication Patterns ───────────────────────────────────────────────
# Step 1: Select optimal number of patterns (k) for CON
par(mfrow = c(1,1))
selectK(cellchat_CON, pattern = "outgoing")
pdf(file.path(out_dir, "selectK_CON_outgoing.pdf"), width=6, height=4)
selectK(cellchat_CON, pattern = "outgoing")
dev.off()

pdf(file.path(out_dir, "selectK_IPF_outgoing.pdf"), width=6, height=4)
selectK(cellchat_IPF, pattern = "outgoing")
dev.off()

pdf(file.path(out_dir, "selectK_CON_incoming.pdf"), width=6, height=4)
selectK(cellchat_CON, pattern = "incoming")
dev.off()

pdf(file.path(out_dir, "selectK_IPF_incoming.pdf"), width=6, height=4)
selectK(cellchat_IPF, pattern = "incoming")
dev.off()

selectK(cellchat_IPF, pattern = "outgoing")

# ── Run NMF patterns ──────────────────────────────────────────────────────────
set.seed(42)
cellchat_CON <- identifyCommunicationPatterns(
  cellchat_CON,
  pattern = "outgoing",
  k = 4,
  width = 8, height = 6
)

set.seed(42)
cellchat_IPF <- identifyCommunicationPatterns(
  cellchat_IPF,
  pattern = "outgoing",
  k = 5,
  width = 8, height = 6
)

# ── Also run incoming patterns ────────────────────────────────────────────────
set.seed(42)
cellchat_CON <- identifyCommunicationPatterns(
  cellchat_CON,
  pattern = "incoming",
  k = 4,
  width = 8, height = 6
)

set.seed(42)
cellchat_IPF <- identifyCommunicationPatterns(
  cellchat_IPF,
  pattern = "incoming",
  k = 5,
  width = 8, height = 6
)

cat("NMF done for both conditions.\n")


# Check actual structure
str(cellchat_CON@netP$pattern$outgoing, max.level = 3)



# ── Correct extraction and publication figures ─────────────────────────────────

build_pattern_plots <- function(cc, condition, k_out, k_in) {
  
  # ── Outgoing ────────────────────────────────────────────────────────────────
  cell_out <- cc@netP$pattern$outgoing$pattern$cell %>%
    mutate(CellGroup = gsub("Monocyte-derived Macrophages", 
                            "Macrophages", as.character(CellGroup)))
  
  sig_out <- cc@netP$pattern$outgoing$pattern$signaling
  
  # Normalize contribution per pattern
  cell_out <- cell_out %>%
    group_by(Pattern) %>%
    mutate(Score = Contribution / max(Contribution)) %>%
    ungroup()
  
  sig_out <- sig_out %>%
    group_by(Pattern) %>%
    mutate(Score = Contribution / max(Contribution)) %>%
    ungroup()
  
  # Top 12 pathways per pattern
  top_sig <- sig_out %>%
    group_by(Signaling) %>%
    summarise(max_score = max(Score)) %>%
    arrange(desc(max_score)) %>%
    slice_head(n = 20) %>%
    pull(Signaling)
  
  sig_top <- sig_out %>% filter(Signaling %in% top_sig)
  
  # Cell heatmap
  p_cell_out <- ggplot(cell_out, 
                       aes(x = Pattern, y = CellGroup, fill = Score)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(data = cell_out %>% filter(Score > 0.25),
              aes(label = sprintf("%.2f", Score)),
              size = 3.2, color = "white", fontface = "bold") +
    scale_fill_gradientn(
      colors = c("white","#c6dbef","#6baed6","#2171b5","#084594"),
      name = "Relative\nScore", limits = c(0,1)
    ) +
    labs(
      title = paste0(condition, " | Outgoing Patterns (k=", k_out, ")\nCell Type Contributions"),
      x = NULL, y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.text.x  = element_text(face = "bold", size = 11, color = "black"),
      axis.text.y  = element_text(size = 10, color = "black"),
      legend.title = element_text(face = "bold", size = 9),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5)
    )
  
  # Pathway heatmap
  p_sig_out <- ggplot(sig_top,
                      aes(x = Pattern, y = Signaling, fill = Score)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(data = sig_top %>% filter(Score > 0.6),
              aes(label = sprintf("%.2f", Score)),
              size = 2.8, color = "white", fontface = "bold") +
    scale_fill_gradientn(
      colors = c("white","#fee8c8","#fdbb84","#e34a33","#99000d"),
      name = "Relative\nScore", limits = c(0,1)
    ) +
    labs(
      title = paste0(condition, " | Outgoing Patterns (k=", k_out, ")\nSignaling Pathway Contributions"),
      x = NULL, y = "Signaling Pathway"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.text.x  = element_text(face = "bold", size = 11, color = "black"),
      axis.text.y  = element_text(size = 9, color = "black"),
      legend.title = element_text(face = "bold", size = 9),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5)
    )
  
  # ── Incoming ────────────────────────────────────────────────────────────────
  cell_in <- cc@netP$pattern$incoming$pattern$cell %>%
    mutate(CellGroup = gsub("Monocyte-derived Macrophages",
                            "Macrophages", as.character(CellGroup))) %>%
    group_by(Pattern) %>%
    mutate(Score = Contribution / max(Contribution)) %>%
    ungroup()
  
  sig_in <- cc@netP$pattern$incoming$pattern$signaling %>%
    group_by(Pattern) %>%
    mutate(Score = Contribution / max(Contribution)) %>%
    ungroup()
  
  top_sig_in <- sig_in %>%
    group_by(Signaling) %>%
    summarise(max_score = max(Score)) %>%
    arrange(desc(max_score)) %>%
    slice_head(n = 20) %>%
    pull(Signaling)
  
  sig_in_top <- sig_in %>% filter(Signaling %in% top_sig_in)
  
  p_cell_in <- ggplot(cell_in,
                      aes(x = Pattern, y = CellGroup, fill = Score)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(data = cell_in %>% filter(Score > 0.25),
              aes(label = sprintf("%.2f", Score)),
              size = 3.2, color = "white", fontface = "bold") +
    scale_fill_gradientn(
      colors = c("white","#c7e9c0","#74c476","#238b45","#00441b"),
      name = "Relative\nScore", limits = c(0,1)
    ) +
    labs(
      title = paste0(condition, " | Incoming Patterns (k=", k_in, ")\nCell Type Contributions"),
      x = NULL, y = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.text.x  = element_text(face = "bold", size = 11, color = "black"),
      axis.text.y  = element_text(size = 10, color = "black"),
      legend.title = element_text(face = "bold", size = 9),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5)
    )
  
  p_sig_in <- ggplot(sig_in_top,
                     aes(x = Pattern, y = Signaling, fill = Score)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(data = sig_in_top %>% filter(Score > 0.6),
              aes(label = sprintf("%.2f", Score)),
              size = 2.8, color = "white", fontface = "bold") +
    scale_fill_gradientn(
      colors = c("white","#dadaeb","#9e9ac8","#6a51a3","#3f007d"),
      name = "Relative\nScore", limits = c(0,1)
    ) +
    labs(
      title = paste0(condition, " | Incoming Patterns (k=", k_in, ")\nSignaling Pathway Contributions"),
      x = NULL, y = "Signaling Pathway"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
      axis.text.x  = element_text(face = "bold", size = 11, color = "black"),
      axis.text.y  = element_text(size = 9, color = "black"),
      legend.title = element_text(face = "bold", size = 9),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5)
    )
  
  return(list(cell_out  = p_cell_out,
              sig_out   = p_sig_out,
              cell_in   = p_cell_in,
              sig_in    = p_sig_in))
}

# ── Generate plots ────────────────────────────────────────────────────────────
plots_CON <- build_pattern_plots(cellchat_CON, "CON", k_out=4, k_in=4)
plots_IPF <- build_pattern_plots(cellchat_IPF, "IPF", k_out=5, k_in=5)

# ── Save outgoing panel ───────────────────────────────────────────────────────
p_out <- (plots_CON$cell_out | plots_CON$sig_out) /
  (plots_IPF$cell_out | plots_IPF$sig_out) +
  plot_annotation(
    title    = "Outgoing Communication Patterns: CON vs IPF",
    subtitle = "NMF decomposition | Blue = cell contributions | Red = pathway contributions",
    tag_levels = "A",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
    )
  )

ggsave(file.path(out_dir, "10_NMF_outgoing.pdf"),
       p_out, width = 16, height = 14, dpi = 300)
ggsave(file.path(out_dir, "10_NMF_outgoing.png"),
       p_out, width = 16, height = 14, dpi = 300)

# ── Save incoming panel ───────────────────────────────────────────────────────
p_in <- (plots_CON$cell_in | plots_CON$sig_in) /
  (plots_IPF$cell_in | plots_IPF$sig_in) +
  plot_annotation(
    title    = "Incoming Communication Patterns: CON vs IPF",
    subtitle = "NMF decomposition | Green = cell contributions | Purple = pathway contributions",
    tag_levels = "A",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
    )
  )

ggsave(file.path(out_dir, "10_NMF_incoming.pdf"),
       p_in, width = 16, height = 14, dpi = 300)
ggsave(file.path(out_dir, "10_NMF_incoming.png"),
       p_in, width = 16, height = 14, dpi = 300)

print(p_out)
print(p_in)

print(plots_CON$cell_out)

print(plots_CON$sig_out)

print(plots_IPF$cell_out)

print(plots_IPF$sig_out)

print(plots_CON$cell_in)
print(plots_CON$sig_in)
print(plots_IPF$cell_in)
print(plots_IPF$sig_in)

install.packages("cowplot")


cellchat_CON <- readRDS("C:/Users/aliaa/Desktop/CellChat_Results/cellchat_CON_inferred.rds")
cellchat_IPF <- readRDS("C:/Users/aliaa/Desktop/CellChat_Results/cellchat_IPF_inferred.rds")

groupSize_CON <- as.numeric(table(cellchat_CON@idents))
groupSize_IPF <- as.numeric(table(cellchat_IPF@idents))

my_cols <- c(
  "Monocyte-derived Macrophages" = "#E63946",
  "Endothelial"                  = "#F4A261",
  "Epithelial"                   = "#2A9D8F",
  "NK Cells"                     = "#457B9D",
  "Fibroblasts"                  = "#1D3557",
  "Dendritic Cells"              = "#E9C46A",
  "Plasma B Cells"               = "#9B5DE5",
  "Mast Cells"                   = "#F15BB5"
)

out_dir <- "C:/Users/aliaa/Desktop/CellChat_Results"
library(CellChat)
# Shared max for comparable scaling
max_count  <- max(max(cellchat_CON@net$count),  max(cellchat_IPF@net$count))
max_weight <- max(max(cellchat_CON@net$weight), max(cellchat_IPF@net$weight))

# Panel A — CON count
png(file.path(out_dir, "tmp_A.png"), width=700, height=650, res=100)
par(mar=c(2,2,4,2))
netVisual_circle(cellchat_CON@net$count, vertex.weight=groupSize_CON,
                 weight.scale=TRUE, edge.weight.max=max_count,
                 edge.width.max=12, label.edge=FALSE, color.use=my_cols,
                 title.name="A   CON — Number of Interactions",
                 vertex.label.cex=0.85, arrow.size=0.12)
dev.off()

# Panel B — IPF count
png(file.path(out_dir, "tmp_B.png"), width=700, height=650, res=100)
par(mar=c(2,2,4,2))
netVisual_circle(cellchat_IPF@net$count, vertex.weight=groupSize_IPF,
                 weight.scale=TRUE, edge.weight.max=max_count,
                 edge.width.max=12, label.edge=FALSE, color.use=my_cols,
                 title.name="B   IPF — Number of Interactions",
                 vertex.label.cex=0.85, arrow.size=0.12)
dev.off()

# Panel C — CON strength
png(file.path(out_dir, "tmp_C.png"), width=700, height=650, res=100)
par(mar=c(2,2,4,2))
netVisual_circle(cellchat_CON@net$weight, vertex.weight=groupSize_CON,
                 weight.scale=TRUE, edge.weight.max=max_weight,
                 edge.width.max=12, label.edge=FALSE, color.use=my_cols,
                 title.name="C   CON — Interaction Strength",
                 vertex.label.cex=0.85, arrow.size=0.12)
dev.off()

# Panel D — IPF strength
png(file.path(out_dir, "tmp_D.png"), width=700, height=650, res=100)
par(mar=c(2,2,4,2))
netVisual_circle(cellchat_IPF@net$weight, vertex.weight=groupSize_IPF,
                 weight.scale=TRUE, edge.weight.max=max_weight,
                 edge.width.max=12, label.edge=FALSE, color.use=my_cols,
                 title.name="D   IPF — Interaction Strength",
                 vertex.label.cex=0.85, arrow.size=0.12)
dev.off()

# Combine into one publication figure
library(cowplot)
library(ggplot2)

pA <- ggdraw() + draw_image(file.path(out_dir, "tmp_A.png"))
pB <- ggdraw() + draw_image(file.path(out_dir, "tmp_B.png"))
pC <- ggdraw() + draw_image(file.path(out_dir, "tmp_C.png"))
pD <- ggdraw() + draw_image(file.path(out_dir, "tmp_D.png"))

final <- plot_grid(pA, pB, pC, pD, nrow=2, ncol=2)

ggsave(file.path(out_dir, "03_aggregate_network_CON_IPF.pdf"),
       final, width=14, height=13, dpi=300)

cat("Circle plot saved.\n")

file.exists(file.path(out_dir, "tmp_A.png"))
file.exists(file.path(out_dir, "tmp_B.png"))
file.exists(file.path(out_dir, "tmp_C.png"))
file.exists(file.path(out_dir, "tmp_D.png"))

install.packages("magick")
library(magick)
library(cowplot)

pA <- ggdraw() + draw_image(file.path(out_dir, "tmp_A.png"))
pB <- ggdraw() + draw_image(file.path(out_dir, "tmp_B.png"))
pC <- ggdraw() + draw_image(file.path(out_dir, "tmp_C.png"))
pD <- ggdraw() + draw_image(file.path(out_dir, "tmp_D.png"))

final <- plot_grid(pA, pB, pC, pD, nrow=2, ncol=2)

ggsave(file.path(out_dir, "03_aggregate_network_CON_IPF.pdf"),
       final, width=14, height=13, dpi=300)

cat("Done\n")

library(magick)
img_C <- image_read(file.path(out_dir, "tmp_C.png"))
print(img_C)

library(magick)
library(cowplot)

pA <- ggdraw() + draw_image(file.path(out_dir, "tmp_A.png"))
pB <- ggdraw() + draw_image(file.path(out_dir, "tmp_B.png"))
pC <- ggdraw() + draw_image(file.path(out_dir, "tmp_C.png"))
pD <- ggdraw() + draw_image(file.path(out_dir, "tmp_D.png"))

final <- plot_grid(pA, pB, pC, pD, nrow=2, ncol=2, scale=0.95)

ggsave(file.path(out_dir, "03_aggregate_network_CON_IPF.pdf"),
       final, width=16, height=16, dpi=300)

cat("Done\n")


library(cowplot)
library(magick)

pA <- ggdraw() + draw_image(file.path(out_dir, "tmp_A.png")) +
  draw_label("A", x=0.02, y=0.98, hjust=0, vjust=1, fontface="bold", size=20) +
  draw_label("CON — Number of Interactions", x=0.5, y=0.98, hjust=0.5, vjust=1, fontface="bold", size=14, color="grey20")

pB <- ggdraw() + draw_image(file.path(out_dir, "tmp_B.png")) +
  draw_label("B", x=0.02, y=0.98, hjust=0, vjust=1, fontface="bold", size=20) +
  draw_label("IPF — Number of Interactions", x=0.5, y=0.98, hjust=0.5, vjust=1, fontface="bold", size=14, color="grey20")

pC <- ggdraw() + draw_image(file.path(out_dir, "tmp_C.png")) +
  draw_label("C", x=0.02, y=0.98, hjust=0, vjust=1, fontface="bold", size=20) +
  draw_label("CON — Interaction Strength", x=0.5, y=0.98, hjust=0.5, vjust=1, fontface="bold", size=14, color="grey20")

pD <- ggdraw() + draw_image(file.path(out_dir, "tmp_D.png")) +
  draw_label("D", x=0.02, y=0.98, hjust=0, vjust=1, fontface="bold", size=20) +
  draw_label("IPF — Interaction Strength", x=0.5, y=0.98, hjust=0.5, vjust=1, fontface="bold", size=14, color="grey20")

final <- plot_grid(pA, pB, pC, pD, nrow=2, ncol=2, scale=0.95) +
  theme(plot.background = element_rect(fill="white", color=NA))

ggsave(file.path(out_dir, "03_aggregate_network_CON_IPF.pdf"),
       final, width=16, height=16, dpi=300)
ggsave(file.path(out_dir, "03_aggregate_network_CON_IPF.png"),
       final, width=16, height=16, dpi=300)
cat("Done\n")




library(CellChat)
library(dplyr)

cellchat_IPF <- readRDS("C:/Users/aliaa/Desktop/CellChat_Results/cellchat_IPF_inferred.rds")
cellchat_CON <- readRDS("C:/Users/aliaa/Desktop/CellChat_Results/cellchat_CON_inferred.rds")

# Extract all significant LR pairs in IPF
lr_ipf <- subsetCommunication(cellchat_IPF)

# Filter for your key pathways
key <- lr_ipf %>%
  filter(pathway_name %in% c("SPP1", "CXCL", "MIF", "COLLAGEN", "CCL")) %>%
  arrange(desc(prob)) %>%
  select(source, target, ligand, receptor, pathway_name, prob, pval)
View(key)


# Summarize receptors cleanly
receptor_summary <- key %>%
  mutate(receptor_clean = case_when(
    grepl("CXCR4", receptor) ~ "CXCR4",
    grepl("CD44",  receptor) ~ "CD44",
    grepl("ACKR1", receptor) ~ "ACKR1",
    grepl("ITGB1", receptor) ~ "ITGB1",
    TRUE ~ receptor
  )) %>%
  group_by(receptor_clean, pathway_name) %>%
  summarise(
    n_interactions = n(),
    mean_prob      = round(mean(prob), 4),
    max_prob       = round(max(prob), 4),
    .groups        = "drop"
  ) %>%
  arrange(desc(max_prob))

print(as.data.frame(receptor_summary))

# SPP1 macrophage → fibroblast specifically
spp1_mf <- key %>%
  filter(pathway_name == "SPP1",
         source == "Monocyte-derived Macrophages",
         target == "Fibroblasts")
print(as.data.frame(spp1_mf))


# Run this to get your network pharmacology target confirmation
# This is what formally justifies CD44 and CXCR4 to a reviewer

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(STRINGdb)
library(igraph)

targets <- c("CD44", "CXCR4", "MIF", "SPP1", 
             "CXCL12", "CD74", "ITGB1", "CCL2")




library(igraph)

ppi <- read.table("C:/Users/aliaa/Desktop/string_interactions.tsv",
                  header = TRUE, sep = "\t")

# Check columns
cat("Columns:", colnames(ppi), "\n")

# Build network — columns 1 and 2 are the interacting proteins
g <- graph_from_data_frame(ppi[, 1:2], directed = FALSE)

# Three hub metrics
deg <- degree(g)
bet <- betweenness(g, normalized = TRUE)
clo <- closeness(g, normalized = TRUE)

# Combined ranking table
hub_df <- data.frame(
  gene        = names(deg),
  degree      = as.numeric(deg),
  betweenness = round(as.numeric(bet), 3),
  closeness   = round(as.numeric(clo), 3)
)
hub_df$hub_score <- scale(hub_df$degree) + 
  scale(hub_df$betweenness) + 
  scale(hub_df$closeness)
hub_df <- hub_df[order(-hub_df$hub_score), ]
rownames(hub_df) <- NULL

cat("\n=== Final hub ranking ===\n")
print(hub_df)

write.csv(hub_df, "C:/Users/aliaa/Desktop/hub_genes_ranked.csv",
          row.names = FALSE)
cat("\nSaved.\n")
