################-- Libraries --################
library(dplyr)
library(Seurat) #v4
library(ggplot2)
library(stringr)
library(org.Hs.eg.db)
library(SingleCellExperiment)
library(data.table)
library(VennDiagram)
library("ggpubr")
library(corrplot)
library(svglite)
library(ComplexHeatmap)

################ --Folders and functions --################
generate_folder <- function(foldername) {
    workDir <- getwd()
    subDir <- foldername
    results_path <- file.path(workDir, subDir)
    if (file.exists(subDir)) {
    } else {
        dir.create(results_path)
    }
    return(results_path)
}

SAMPLE_PATH_SQV <- "./2.CellRanger_GRCh38/"
SAMPLE_PATH_CR <- "./2.CellRanger_GRCh38_SARS_CoV/"

results_path <- "4.SC_Analysis/"
generate_folder(results_path)

theme_Publication <- function(base_size = 14, base_family = "arial") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
            plot.title = element_text(
                face = "bold",
                size = rel(1.2), hjust = 0.5
            ),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold", size = rel(1)),
            axis.title.y = element_text(angle = 90, vjust = 2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour = "black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour = "#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size = unit(0.4, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(size = 10, face = "bold"),
            plot.margin = unit(c(10, 5, 5, 5), "mm"),
            strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
            strip.text = element_text(face = "bold")
        ))
}

scale_fill_Publication <- function(...) {
    library(scales)
    discrete_scale("fill", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ...)
}

scale_colour_Publication <- function(...) {
    library(scales)
    discrete_scale("colour", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ...)
}

run_normalization <- function(immune.combined, SD = TRUE, SCTRANSFORM = FALSE) {
    print("STATUS: normalizing data")
    all.genes <- rownames(immune.combined)
    if (SD == TRUE) {
        immune.combined <- ScaleData(immune.combined,
            vars.to.regress = "percent.mt",
            features = all.genes, verbose = FALSE
        )
    } else if (SCTRANSFORM == TRUE) {
        immune.combined <- SCTransform(immune.combined,
            vars.to.regress = "percent.mt",
            verbose = FALSE
        )
    }
    return(immune.combined)
}

feature_reduction <- function(immune.combined, result_folder) {
    print("STATUS: performing PCA and UMAP")

    results_path <- generate_folder(result_folder)
    # unlink(paste0(result_folder, "/*"))

    immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
    immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)

    DimPlot(immune.combined,
        reduction = "umap",
        group.by = "orig.ident",
        label = FALSE
    )
    ggsave(paste0(result_folder, "UMAP_samples.png"), height = 5, width = 8, dpi = 500)

    DimPlot(immune.combined,
        reduction = "umap",
        group.by = "ZIKA",
        cols = c("#ddd2d2", "#2d0446"),
        label = FALSE
    )
    ggsave(paste0(result_folder, "UMAP_ZIKA.png"), height = 5, width = 8, dpi = 500)


    return(immune.combined)
}

initial_SC_Seurat_CR <- function(data_dir, sampleID, cells, mock = FALSE) {
    pbmc.data <- Read10X(data.dir = data_dir)

    pbmc <- CreateSeuratObject(
        counts = pbmc.data,
        project = sampleID,
        min.cells = 0,
        min.features = 0
    )
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave(paste0(results_folder, sampleID, "_VnPlotMt.png"), dpi = 300)

    pbmc <- subset(pbmc, cells = cells)
    VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave(paste0(results_folder, sampleID, "_VnPlotMt_Filtered.png"), dpi = 300)

    message("STATUS: Number of cells cell ranger ", dim(pbmc@meta.data)[1])
    # print ("hiv positive")
    z <- pbmc["sars-cov2-wh1", ]
    x1 <- z@meta.data[z@meta.data$nFeature_RNA != 0, ]
    COVIDdata <- c()
    for (cell in colnames(pbmc)) {
        if (cell %in% rownames(x1)) {
            COVIDdata <- c(COVIDdata, "COVID")
        } else {
            COVIDdata <- c(COVIDdata, "noCOVID")
        }
    }
    pbmc$COVID <- COVIDdata

    
    pbmc <- NormalizeData(pbmc,
        normalization.method = "LogNormalize",
        scale.factor = 10000
    )

    print ("end")
    # Find features (genes) that are highly variable from cell-to-cell
    pbmc <- FindVariableFeatures(pbmc,
        selection.method = "vst",
        nfeatures = 2000
    )

    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(pbmc), 10)

    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(pbmc)
    LabelPoints(plot = plot1, points = top10, repel = TRUE)
    ggsave(paste0(results_folder, sampleID, "_HighlyVariableGenes.png"), dpi = 300)

    pbmc <- ScaleData(pbmc)
    pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
    pbmc <- RunUMAP(pbmc, dims = 1:10)

    VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
    ggsave(paste0(results_folder, sampleID, "_PCA_Loadings.png"), dpi = 300)

    DimPlot(pbmc, reduction = "pca")
    ggsave(paste0(results_folder, sampleID, "_PCA.png"), dpi = 300)

    DimPlot(pbmc, reduction = "umap")
    ggsave(paste0(results_folder, sampleID, "_UMAP.png"), dpi = 300)

    FeaturePlot(pbmc, reduction = "umap", features = c("sars-cov2-wh1"))
    ggsave(paste0(results_folder, sampleID, "_COVIDgenes.png"), width = 6, height = 4, dpi = 300)

    VlnPlot(pbmc, features = c("sars-cov2-wh1"))
    ggsave(paste0(results_folder, sampleID, "_COVIDgenes_VnPlot.png"), dpi = 300)

    z <- pbmc["sars-cov2-wh1", ]
    x1 <- z@meta.data[z@meta.data$nFeature_RNA != 0, ]

    if (length(rownames(x1)) > 0) {
        print(paste0("STATUS: Generating sars-cov2-wh1 figure for ", sampleID))
        pbmccovid <- subset(pbmc, cells = rownames(x1))

        print(pbmccovid)
        VlnPlot(pbmccovid, features = c("sars-cov2-wh1"))
        ggsave(paste0(results_folder, sampleID, "_COVID_VnPlotMt_onlyCOVID:cells.png"), dpi = 300)
    }

    return(pbmc)
}

initial_SC_Seurat_SQV <- function(data_dir, sampleID, cells, mock = FALSE) {
    pbmc.data <- Read10X(data.dir = data_dir)


    pbmc <- CreateSeuratObject(
        counts = pbmc.data,
        project = sampleID,
        min.cells = 0,
        min.features = 0
    )
    pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
    VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave(paste0(results_folder, sampleID, "_VnPlotMt.png"), dpi = 300)

    # pbmc <- subset(pbmc, cells = cells)
    VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    ggsave(paste0(results_folder, sampleID, "_VnPlotMt_Filtered.png"), dpi = 300)

    message("STATUS: Number of cells SQV ", dim(pbmc@meta.data)[1])

    z <- pbmc["sars-CoV2-WuhanHu1", ]
    x1 <- z@meta.data[z@meta.data$nFeature_RNA != 0, ]
    print(head(z))
    COVIDdata <- c()
    for (cell in colnames(pbmc)) {
        if (cell %in% rownames(x1)) {
            COVIDdata <- c(COVIDdata, "COVID")
        } else {
            COVIDdata <- c(COVIDdata, "noCOVID")
        }
    }
    pbmc$COVID <- COVIDdata
    pbmc <- NormalizeData(pbmc,
        normalization.method = "LogNormalize",
        scale.factor = 10000
    )


    # Find features (genes) that are highly variable from cell-to-cell
    pbmc <- FindVariableFeatures(pbmc,
        selection.method = "vst",
        nfeatures = 2000
    )

    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(pbmc), 10)

    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(pbmc)
    LabelPoints(plot = plot1, points = top10, repel = TRUE)
    ggsave(paste0(results_folder, sampleID, "_HighlyVariableGenes.png"), dpi = 300)

    pbmc <- ScaleData(pbmc)
    pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc))
    pbmc <- RunUMAP(pbmc, dims = 1:10)

    VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
    ggsave(paste0(results_folder, sampleID, "_PCA_Loadings.png"), dpi = 300)

    DimPlot(pbmc, reduction = "pca")
    ggsave(paste0(results_folder, sampleID, "_PCA.png"), dpi = 300)

    DimPlot(pbmc, reduction = "umap")
    ggsave(paste0(results_folder, sampleID, "_UMAP.png"), dpi = 300)

    FeaturePlot(pbmc, reduction = "umap", features = c("sars-CoV2-WuhanHu1"))
    ggsave(paste0(results_folder, sampleID, "_COVIDgenes.png"), width = 6, height = 4, dpi = 300)

    VlnPlot(pbmc, features = c("sars-CoV2-WuhanHu1"))
    ggsave(paste0(results_folder, sampleID, "_COVIDgenes_VnPlot.png"), dpi = 300)


    z <- pbmc["sars-CoV2-WuhanHu1", ]
    x1 <- z@meta.data[z@meta.data$nFeature_RNA != 0, ]

    if (length(rownames(x1)) > 0) {
        print(paste0("STATUS: Generating COVID figure for ", sampleID))
        pbmccovid <- subset(pbmc, cells = rownames(x1))

        print(pbmccovid)
        VlnPlot(pbmccovid, features = c("sars-CoV2-WuhanHu1"))
        ggsave(paste0(results_folder, sampleID, "_COVID_VnPlotMt_onlyCOVIDcells.png"), dpi = 300)
    }

    return(pbmc)
}

scatterplots4remaining <- function(overlap, SRR25cr, SRR25sqv, figname) {
    SRR25crhivov <- subset(SRR25cr, cells = overlap)
    zcrhivov <- SRR25crhivov["sars-cov2-wh1", ]
    gfpcrhivov <- zcrhivov@meta.data[zcrhivov@meta.data$nFeature_RNA != 0, ]

    SRR25sqvhivov <- subset(SRR25sqv, cells = overlap)
    zsqvhivov <- SRR25sqvhivov["sars-CoV2-WuhanHu1", ]
    gfpsqvhivov <- zsqvhivov@meta.data[zsqvhivov@meta.data$nFeature_RNA != 0, ]

    df <- cbind(gfpsqvhivov[, "nCount_RNA"], gfpcrhivov[, "nCount_RNA"])
    rownames(df) <- rownames(gfpsqvhivov)
    colnames(df) <- c("sqv", "cr")
    df <- as.data.frame(df)
    ggscatter(df,
        x = "sqv", y = "cr", size = 5, shape=21,
        add = "reg.line", conf.int = TRUE, font.label = c(16, "plain"),
        cor.coef = TRUE, cor.method = "kendall",
        xlab = "SVQ UMI counts", ylab = "Cellranger UMI Counts",
        cor.coeff.args = list(method = "kendall", size=6, label.x = 10, label.sep = "\n"),
        ggthemes=theme_minimal()) + font("xy.text", size = 14) + 
        font("xlab", size = 16) + font("ylab", size = 16) +
        scale_x_continuous(breaks = get_breaks(n = 10)) +
        scale_y_continuous(breaks = get_breaks(n = 10)) +
        rotate_x_text(angle = 90, hjust = 1, vjust = 0.5)
    ggsave(paste0(figname, "_scatter.svg"), dpi = 300, height = 4, width = 4, units = "in")
    ggsave(paste0(figname, "_scatter.png"), dpi = 300, height = 4, width = 4, units = "in")
    ggsave(paste0(figname, "_scatter.pdf"), dpi = 300, height = 4, width = 4, units = "in")

    message("STATUS: mean sqv remaining cells ", mean(gfpsqvhivov[, "nCount_RNA"]))
    message("STATUS: mean cr remaining cells ", mean(gfpcrhivov[, "nCount_RNA"]))
    # dfmelt1 <- melt(df)
    # dfmelt2 <- rbind(dfmelt, dfmelt1)
    # p <- ggplot(dfmelt2, aes(x = value, color = variable)) +
    #     geom_density()
    # ggsave(paste0(figname, "_density.png"), dpi = 300, height = 4, width = 4, units = "in")
}

upset_plots <- function(remaining, results_folder, at=c(0, -1000, -2000, -3000, -4000), 
    labels=c(0, 1000, 2000, 3000, 4000)) {
    # m <- list_to_matrix(remaining)
    m <- make_comb_mat(remaining)
    cs <- comb_size(m)
    ss <- set_size(m)

    ht <- UpSet(m,
        set_order = c("CellRanger", "SVQ bowtie2", "SVQ bbmap"),
        lwd = 2, top_annotation = HeatmapAnnotation(
            "Intersection of cells\ndetected with\nSars-CoV-2" = anno_barplot(cs,
                ylim = c(0, max(cs) * 1.1),
                border = FALSE,
                gp = gpar(fill = "black", fontsize = 6),
                height = unit(4, "cm"),
                axis_param = list(gp = gpar(fontsize = 12))
            ),
            annotation_name_side = "left",
            annotation_name_rot = 90
        ),
        right_annotation = rowAnnotation(
            "Sars-CoV-2+ cells per\nalignment method" = anno_barplot(-ss,
                baseline = 0,
                axis_param = list(
                    at = at,
                    labels = labels,
                    labels_rot = 90, gp = gpar(fontsize = 12, fill = "black")
                ),
                border = FALSE,
                gp = gpar(fill = "black"),
                width = unit(2, "cm")
            ), set_name = anno_text(set_name(m),
                location = 0.05,
                just = "left", gp = gpar(fontsize = 14),
                width = max_text_width(set_name(m)) + unit(4, "mm")
            )
        ),
        left_annotation = NULL, show_row_names = FALSE
    )

    svglite(file.path(results_folder, "Upsetplot.svg"), width = 5.5, height = 3.4)
    ht <- draw(ht)
    od <- column_order(ht)
    decorate_annotation("Intersection of cells\ndetected with\nSars-CoV-2", {
        grid.text(cs[od],
            x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
            default.units = "native", just = c("left", "bottom"),
            gp = gpar(fontsize = 14, col = "#404040"), rot = 0
        )
    })
    dev.off()

    pdf(file.path(results_folder, "Upsetplot.pdf"), width = 5.5, height = 3.4)
    ht <- draw(ht)
    od <- column_order(ht)
    decorate_annotation("Intersection of cells\ndetected with\nSars-CoV-2", {
        grid.text(cs[od],
            x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
            default.units = "native", just = c("left", "bottom"),
            gp = gpar(fontsize = 14, col = "#404040"), rot = 0
        )
    })
    dev.off()

    png(file.path(results_folder, "Upsetplot.png"), width = 5.5, height = 3.4, units = "in", res = 500)
    ht <- draw(ht)
    od <- column_order(ht)
    decorate_annotation("Intersection of cells\ndetected with\nSars-CoV-2", {
        grid.text(cs[od],
            x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
            default.units = "native", just = c("left", "bottom"),
            gp = gpar(fontsize = 14, col = "#404040"), rot = 0
        )
    })
    dev.off()

}
################ -- Processing 1dpi_CoV2_HHT --###############

# # -- Pull out intersecting cells to include in analysis (want to include only overlapping cells)
pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_CR, "1dpi_CoV2_HHT","outs", "filtered_feature_bc_matrix"))

pbmccr <- CreateSeuratObject(
    counts = pbmc.data,
    project = "1dpi_CoV2_HHT",
    min.cells = 0,
    min.features = 0
)

pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_SQV, "1dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix_bowtie2_global"))

pbmcsqvbow2global <- CreateSeuratObject(
    counts = pbmc.data,
    project = "1dpi_CoV2_HHT",
    min.cells = 0,
    min.features = 0
)

pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_SQV, "1dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix_bbmap_global"))

pbmcsqvbbmapglobal <- CreateSeuratObject(
    counts = pbmc.data,
    project = "1dpi_CoV2_HHT",
    min.cells = 0,
    min.features = 0
)

message("STATUS: total cells in cellranger analysis ", length(rownames(pbmccr@meta.data)))
message("STATUS: total cells in SQV bowtie2 global analysis ", length(rownames(pbmcsqvbow2global@meta.data)))
message("STATUS: total cells in SQV bbmap global analysis ", length(rownames(pbmcsqvbbmapglobal@meta.data)))

cells <- intersect(intersect(rownames(pbmccr@meta.data), rownames(pbmcsqvbow2global@meta.data)),
        rownames(pbmcsqvbbmapglobal@meta.data))

#Number cells can generally very based on read counts per cell this is just a way to standardize cells in the analysis 

message("STATUS: overlapping cells in all 3 analyses ", length(cells))

results_folder_base <- "4.SC_Analysis_CR/"
generate_folder(results_folder_base)
results_folder <- file.path(results_folder_base, "1dpi_CoV2_HHT/")
generate_folder(results_folder)

cr <- initial_SC_Seurat_CR(file.path(SAMPLE_PATH_CR, "1dpi_CoV2_HHT","outs", "filtered_feature_bc_matrix"),
    sampleID = "1dpi_CoV2_HHT", cells = cells
)

results_folder_base <- "4.SC_Analysis/"
generate_folder(results_folder_base)

results_folder <- file.path(results_folder_base, "1dpi_CoV2_HHT", "bowtie2_global/")
generate_folder(results_folder)
sqvbow2global <- initial_SC_Seurat_SQV(file.path(SAMPLE_PATH_SQV, "1dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix_bowtie2_global"),
    sampleID = "1dpi_CoV2_HHT", cells = cells
)

results_folder <- file.path(results_folder_base, "1dpi_CoV2_HHT", "bbmap_global/")
generate_folder(results_folder)
sqvbbmapglobal <- initial_SC_Seurat_SQV(file.path(SAMPLE_PATH_SQV, "1dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix_bbmap_global"),
    sampleID = "1dpi_CoV2_HHT", cells = cells
)

zcr <- cr["sars-cov2-wh1", ]
crcells <- zcr@meta.data[zcr@meta.data$nFeature_RNA != 0, ]

zsqvhivbow2global <- sqvbow2global["sars-CoV2-WuhanHu1", ]
sqvbow2globalcells <- zsqvhivbow2global@meta.data[zsqvhivbow2global@meta.data$nFeature_RNA != 0, ]
zsqvhivbbmapglobal <- sqvbbmapglobal["sars-CoV2-WuhanHu1", ]
sqvbbmapglobalcells <- zsqvhivbbmapglobal@meta.data[zsqvhivbbmapglobal@meta.data$nFeature_RNA != 0, ]

vd <- venn.diagram(
    x = list(rownames(crcells), rownames(sqvbow2globalcells)),
    category.names = c("cellranger", "sqv bowtie2 global"), sub.cex = 20,
    filename = file.path(results_folder_base, "1dpi_CoV2_HHT","4.VD_1dpi_CoV2_HHT_COVIDpositivebow2global.png"),
    output = TRUE, imagetype = "png"
)
vd <- venn.diagram(
    x = list(rownames(crcells), rownames(sqvbbmapglobalcells)),
    category.names = c("cellranger", "sqv bbmap global"), sub.cex = 20,
    filename = file.path(results_folder_base, "1dpi_CoV2_HHT","4.VD_1dpi_CoV2_HHT_COVIDpositivebbmapglobal.png"),
    output = TRUE, imagetype = "png"
)

overlapbow2global <- intersect(rownames(sqvbow2globalcells), rownames(crcells))
overlapbbmapglobal <- intersect(rownames(sqvbbmapglobalcells), rownames(crcells))

scatterplots4remaining(overlapbow2global, cr, sqvbow2global, file.path(results_folder_base, "1dpi_CoV2_HHT", "4.Correlationbow2global1dpi_CoV2_HHT"))
scatterplots4remaining(overlapbbmapglobal, cr, sqvbbmapglobal, file.path(results_folder_base, "1dpi_CoV2_HHT", "4.Correlationbbbmapglobal1dpi_CoV2_HHT"))

#-- Upset Plot 

# Remaining upset plot 
remaining <- list(
    "CellRanger" = rownames(crcells),
    "SVQ bowtie2" = rownames(sqvbow2globalcells),
    "SVQ bbmap" = rownames(sqvbbmapglobalcells)
)
upset_plots(remaining, file.path(results_folder_base, "1dpi_CoV2_HHT"))


################ -- Processing 2dpi_CoV2_HHT -- ################

# # -- Pull out intersecting cells to include in analysis (want to include only overlapping cells)
pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_CR, "2dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix"))

pbmccr <- CreateSeuratObject(
    counts = pbmc.data,
    project = "2dpi_CoV2_HHT",
    min.cells = 0,
    min.features = 0
)

pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_SQV, "2dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix_bowtie2_global"))

pbmcsqvbow2global <- CreateSeuratObject(
    counts = pbmc.data,
    project = "2dpi_CoV2_HHT",
    min.cells = 0,
    min.features = 0
)

pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_SQV, "2dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix_bbmap_global"))

pbmcsqvbbmapglobal <- CreateSeuratObject(
    counts = pbmc.data,
    project = "2dpi_CoV2_HHT",
    min.cells = 0,
    min.features = 0
)

message("STATUS: total cells in cellranger analysis ", length(rownames(pbmccr@meta.data)))
message("STATUS: total cells in SQV bowtie2 global analysis ", length(rownames(pbmcsqvbow2global@meta.data)))
message("STATUS: total cells in SQV bbmap global analysis ", length(rownames(pbmcsqvbbmapglobal@meta.data)))

cells <- intersect(
        intersect(rownames(pbmccr@meta.data), rownames(pbmcsqvbow2global@meta.data)
    ), rownames(pbmcsqvbbmapglobal@meta.data)
)

# Number cells can generally very based on read counts per cell this is just a way to standardize cells in the analysis

message("STATUS: overlapping cells in all 3 analyses ", length(cells))

results_folder_base <- "4.SC_Analysis_CR/"
generate_folder(results_folder_base)
results_folder <- file.path(results_folder_base, "2dpi_CoV2_HHT/")
generate_folder(results_folder)

cr <- initial_SC_Seurat_CR(file.path(SAMPLE_PATH_CR, "2dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix"),
    sampleID = "2dpi_CoV2_HHT", cells = cells
)

results_folder_base <- "4.SC_Analysis/"
generate_folder(results_folder_base)

results_folder <- file.path(results_folder_base, "2dpi_CoV2_HHT", "bowtie2_global/")
generate_folder(results_folder)
sqvbow2global <- initial_SC_Seurat_SQV(file.path(SAMPLE_PATH_SQV, "2dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix_bowtie2_global"),
    sampleID = "2dpi_CoV2_HHT", cells = cells
)

results_folder <- file.path(results_folder_base, "2dpi_CoV2_HHT", "bbmap_global/")
generate_folder(results_folder)
sqvbbmapglobal <- initial_SC_Seurat_SQV(file.path(SAMPLE_PATH_SQV, "2dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix_bbmap_global"),
    sampleID = "2dpi_CoV2_HHT", cells = cells
)

zcr <- cr["sars-cov2-wh1", ]
crcells <- zcr@meta.data[zcr@meta.data$nFeature_RNA != 0, ]

zsqvhivbow2global <- sqvbow2global["sars-CoV2-WuhanHu1", ]
sqvbow2globalcells <- zsqvhivbow2global@meta.data[zsqvhivbow2global@meta.data$nFeature_RNA != 0, ]
zsqvhivbbmapglobal <- sqvbbmapglobal["sars-CoV2-WuhanHu1", ]
sqvbbmapglobalcells <- zsqvhivbbmapglobal@meta.data[zsqvhivbbmapglobal@meta.data$nFeature_RNA != 0, ]

vd <- venn.diagram(
    x = list(rownames(crcells), rownames(sqvbow2globalcells)),
    category.names = c("cellranger", "sqv bowtie2 global"), sub.cex = 20,
    filename = file.path(results_folder_base, "2dpi_CoV2_HHT", "4.VD_2dpi_CoV2_HHT_COVIDpositivebow2global.png"),
    output = TRUE, imagetype = "png"
)
vd <- venn.diagram(
    x = list(rownames(crcells), rownames(sqvbbmapglobalcells)),
    category.names = c("cellranger", "sqv bbmap global"), sub.cex = 20,
    filename = file.path(results_folder_base, "2dpi_CoV2_HHT", "4.VD_2dpi_CoV2_HHT_COVIDpositivebbmapglobal.png"),
    output = TRUE, imagetype = "png"
)

overlapbow2global <- intersect(rownames(sqvbow2globalcells), rownames(crcells))
overlapbbmapglobal <- intersect(rownames(sqvbbmapglobalcells), rownames(crcells))

scatterplots4remaining(overlapbow2global, cr, sqvbow2global, file.path(results_folder_base, "2dpi_CoV2_HHT", "4.Correlationbow2global2dpi_CoV2_HHT"))
scatterplots4remaining(overlapbbmapglobal, cr, sqvbbmapglobal, file.path(results_folder_base, "2dpi_CoV2_HHT", "4.Correlationbbbmapglobal2dpi_CoV2_HHT"))

#--Upset Plot

library(ComplexHeatmap)

# Remaining upset plot
remaining <- list(
    "CellRanger" = rownames(crcells),
    "SVQ bowtie2" = rownames(sqvbow2globalcells),
    "SVQ bbmap" = rownames(sqvbbmapglobalcells)
)
upset_plots(remaining, file.path(results_folder_base, "2dpi_CoV2_HHT"), 
    at=c(0, -5000, -10000, -15000), labels= c(0, 5000, 10000, 15000))


################ -- Processing 3dpi_CoV2_HHT --################

# # -- Pull out intersecting cells to include in analysis (want to include only overlapping cells)
pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_CR, "3dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix"))

pbmccr <- CreateSeuratObject(
    counts = pbmc.data,
    project = "3dpi_CoV2_HHT",
    min.cells = 0,
    min.features = 0
)

pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_SQV, "3dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix_bowtie2_global"))

pbmcsqvbow2global <- CreateSeuratObject(
    counts = pbmc.data,
    project = "3dpi_CoV2_HHT",
    min.cells = 0,
    min.features = 0
)

pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_SQV, "3dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix_bbmap_global"))

pbmcsqvbbmapglobal <- CreateSeuratObject(
    counts = pbmc.data,
    project = "3dpi_CoV2_HHT",
    min.cells = 0,
    min.features = 0
)

message("STATUS: total cells in cellranger analysis ", length(rownames(pbmccr@meta.data)))
message("STATUS: total cells in SQV bowtie2 global analysis ", length(rownames(pbmcsqvbow2global@meta.data)))
message("STATUS: total cells in SQV bbmap global analysis ", length(rownames(pbmcsqvbbmapglobal@meta.data)))

cells <- intersect(
        intersect(rownames(pbmccr@meta.data), rownames(pbmcsqvbow2global@meta.data)),
        rownames(pbmcsqvbbmapglobal@meta.data)
)

# Number cells can generally very based on read counts per cell this is just a way to standardize cells in the analysis

message("STATUS: overlapping cells in all 3 analyses ", length(cells))

results_folder_base <- "4.SC_Analysis_CR/"
generate_folder(results_folder_base)
results_folder <- file.path(results_folder_base, "3dpi_CoV2_HHT/")
generate_folder(results_folder)

cr <- initial_SC_Seurat_CR(file.path(SAMPLE_PATH_CR, "3dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix"),
    sampleID = "3dpi_CoV2_HHT", cells = cells
)

results_folder_base <- "4.SC_Analysis/"
generate_folder(results_folder_base)

results_folder <- file.path(results_folder_base, "3dpi_CoV2_HHT", "bowtie2_global/")
generate_folder(results_folder)
sqvbow2global <- initial_SC_Seurat_SQV(file.path(SAMPLE_PATH_SQV, "3dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix_bowtie2_global"),
    sampleID = "3dpi_CoV2_HHT", cells = cells
)

results_folder <- file.path(results_folder_base, "3dpi_CoV2_HHT", "bbmap_global/")
generate_folder(results_folder)
sqvbbmapglobal <- initial_SC_Seurat_SQV(file.path(SAMPLE_PATH_SQV, "3dpi_CoV2_HHT", "outs", "filtered_feature_bc_matrix_bbmap_global"),
    sampleID = "3dpi_CoV2_HHT", cells = cells
)

zcr <- cr["sars-cov2-wh1", ]
crcells <- zcr@meta.data[zcr@meta.data$nFeature_RNA != 0, ]

zsqvhivbow2global <- sqvbow2global["sars-CoV2-WuhanHu1", ]
sqvbow2globalcells <- zsqvhivbow2global@meta.data[zsqvhivbow2global@meta.data$nFeature_RNA != 0, ]
zsqvhivbbmapglobal <- sqvbbmapglobal["sars-CoV2-WuhanHu1", ]
sqvbbmapglobalcells <- zsqvhivbbmapglobal@meta.data[zsqvhivbbmapglobal@meta.data$nFeature_RNA != 0, ]

vd <- venn.diagram(
    x = list(rownames(crcells), rownames(sqvbow2globalcells)),
    category.names = c("cellranger", "sqv bowtie2 global"), sub.cex = 20,
    filename = file.path(results_folder_base, "3dpi_CoV2_HHT", "4.VD_3dpi_CoV2_HHT_COVIDpositivebow2global.png"),
    output = TRUE, imagetype = "png"
)
vd <- venn.diagram(
    x = list(rownames(crcells), rownames(sqvbbmapglobalcells)),
    category.names = c("cellranger", "sqv bbmap global"), sub.cex = 20,
    filename = file.path(results_folder_base, "3dpi_CoV2_HHT", "4.VD_3dpi_CoV2_HHT_COVIDpositivebbmapglobal.png"),
    output = TRUE, imagetype = "png"
)

overlapbow2global <- intersect(rownames(sqvbow2globalcells), rownames(crcells))
overlapbbmapglobal <- intersect(rownames(sqvbbmapglobalcells), rownames(crcells))

scatterplots4remaining(overlapbow2global, cr, sqvbow2global, file.path(results_folder_base, "3dpi_CoV2_HHT", "4.Correlationbow2global3dpi_CoV2_HHT"))
scatterplots4remaining(overlapbbmapglobal, cr, sqvbbmapglobal, file.path(results_folder_base, "3dpi_CoV2_HHT", "4.Correlationbbbmapglobal3dpi_CoV2_HHT"))

#--Upset Plot
library(ComplexHeatmap)

# Remaining upset plot
remaining <- list(
    "CellRanger" = rownames(crcells),
    "SVQ bowtie2" = rownames(sqvbow2globalcells),
    "SVQ bbmap" = rownames(sqvbbmapglobalcells)
)

upset_plots(remaining, file.path(results_folder_base, "3dpi_CoV2_HHT"),
    at = c(0, -10000, -20000, -30000), labels = c(0, 10000, 20000, 30000)
)
