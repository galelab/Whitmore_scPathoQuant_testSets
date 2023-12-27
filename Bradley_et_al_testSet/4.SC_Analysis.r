################-- Libraries--################
library(dplyr)
library(Seurat) #v4
library(ggplot2)
library(stringr)
library(org.Hs.eg.db)
library(SingleCellExperiment)
library(data.table)
library(RColorBrewer)
library(VennDiagram)
library("ggpubr")
library(svglite)
library(corrplot)

################-- Folders and functions--###############

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

SAMPLE_PATH_SQV <- "./2.CellRanger_hg38/"
SAMPLE_PATH_CR <- "./2.CellRanger_hg38_HIV_Redo/"

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
    z <- pbmc["hiv-positive", ]
    print("HIV")
    print(head(z))
    x1 <- z@meta.data[z@meta.data$nFeature_RNA != 0, ]
    HIVdata <- c()
    for (cell in colnames(pbmc)) {
        if (cell %in% rownames(x1)) {
            HIVdata <- c(HIVdata, "HIV")
        } else {
            HIVdata <- c(HIVdata, "noHIV")
        }
    }
    pbmc$HIV <- HIVdata

    print ("GFP")
    z <- pbmc["EGFP-positive", ]
    print(head(z))
    gfp_x1 <- z@meta.data[z@meta.data$nFeature_RNA != 0, ]
    message("STATUS: Number of cells with EGFP cell ranger ", dim(gfp_x1)[1])
    print(length(rownames(gfp_x1)))
    print(length(rownames(gfp_x1)))

    if (length(rownames(gfp_x1)) > 0) {
        vd <- venn.diagram(
            x = list(rownames(gfp_x1), rownames(x1)),
            category.names = c("GFP", "HIV"), sub.cex = 24,
            filename = paste0(results_folder, "VD_GFP_HIV.png"),
            output = TRUE, imagetype = "png"
        )
    }

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

    FeaturePlot(pbmc, reduction = "umap", features = c("hiv-positive"))
    ggsave(paste0(results_folder, sampleID, "_HIVgenes.png"), width = 6, height = 4, dpi = 300)

    FeaturePlot(pbmc, reduction = "umap", features = c("EGFP-positive"))
    ggsave(paste0(results_folder, sampleID, "_GFP.png"), width=6, height=4, dpi = 300)

    VlnPlot(pbmc, features = c("hiv-positive"))
    ggsave(paste0(results_folder, sampleID, "_HIVgenes_VnPlot.png"), dpi = 300)

    z <- pbmc["hiv-positive", ]
    x1 <- z@meta.data[z@meta.data$nFeature_RNA != 0, ]

    if (length(rownames(x1)) > 0) {
        print(paste0("STATUS: Generating HIV figure for ", sampleID))
        pbmchiv <- subset(pbmc, cells = rownames(x1))

        print(pbmchiv)
        VlnPlot(pbmchiv, features = c("hiv-positive"))
        ggsave(paste0(results_folder, sampleID, "_hiv_VnPlotMt_onlyhivcells.png"), dpi = 300)
    }

    return(pbmc)
}

initial_SC_Seurat_SQV <- function(data_dir, sampleID, cells, mock = FALSE) {
    pbmc.data <- Read10X(data.dir = data_dir)


    pbmc <- CreateSeuratObject(
        counts = pbmc.data,
        project = "HIV",
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

    z <- pbmc["HIV-1-pNL4-3", ]
    x1 <- z@meta.data[z@meta.data$nFeature_RNA != 0, ]
    print("STATUS: hiv")
    HIVdata <- c()
    for (cell in colnames(pbmc)) {
        if (cell %in% rownames(x1)) {
            HIVdata <- c(HIVdata, "HIV")
        } else {
            HIVdata <- c(HIVdata, "noHIV")
        }
    }
    pbmc$HIV <- HIVdata
    z <- pbmc["EGFP-positive", ]
    print("STATUS: EGP")
    gfp_x1 <- z@meta.data[z@meta.data$nFeature_RNA != 0, ]
    message("STATUS: Number of cells with EGFP SQV ", dim(gfp_x1)[1])
    if (length(rownames(gfp_x1)) > 0) {
        vd <- venn.diagram(
            x = list(rownames(gfp_x1), rownames(x1)),
            category.names = c("GFP", "HIV"), sub.cex = 20,
            filename = paste0(results_folder, "VD_GFP_HIV.png"),
            output = TRUE, imagetype = "png"
        )
    }
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

    FeaturePlot(pbmc, reduction = "umap", features = c("HIV-1-pNL4-3"))
    ggsave(paste0(results_folder, sampleID, "_HIVgenes.png"), width = 6, height = 4, dpi = 300)

    FeaturePlot(pbmc, reduction = "umap", features = c("EGFP-positive"))
    ggsave(paste0(results_folder, sampleID, "_GFP.png"), width = 6, height = 4, dpi = 300)


    VlnPlot(pbmc, features = c("HIV-1-pNL4-3"))
    ggsave(paste0(results_folder, sampleID, "_HIVgenes_VnPlot.png"), dpi = 300)


    z <- pbmc["HIV-1-pNL4-3", ]
    x1 <- z@meta.data[z@meta.data$nFeature_RNA != 0, ]

    if (length(rownames(x1)) > 0) {
        print(paste0("STATUS: Generating HIV figure for ", sampleID))
        pbmchiv <- subset(pbmc, cells = rownames(x1))

        print(pbmchiv)
        VlnPlot(pbmchiv, features = c("HIV-1-pNL4-3"))
        ggsave(paste0(results_folder, sampleID, "_hiv_VnPlotMt_onlyhivcells.png"), dpi = 300)
    }

    return(pbmc)
}

scatterplots4GFP <- function(gfpcr, gfpsqv, SRR25cr, SRR25sqv, figname) {
    SRR25crgfp <- subset(SRR25cr, cells = rownames(gfpcr))
    zcrhiv <- SRR25crgfp["hiv-positive", ]
    gfpcrhiv <- zcrhiv@meta.data[zcrhiv@meta.data$nFeature_RNA != 0, ]

    SRR25sqvgfp <- subset(SRR25sqv, cells = rownames(gfpsqv))
    zsqvhiv <- SRR25sqvgfp["HIV-1-pNL4-3", ]
    gfpsqvhiv <- zsqvhiv@meta.data[zsqvhiv@meta.data$nFeature_RNA != 0, ]

    df <- cbind(gfpsqvhiv[, "nCount_RNA"], gfpcrhiv[, "nCount_RNA"])
    rownames(df) <- rownames(gfpsqvhiv)
    colnames(df) <- c("sqv", "cr")
    df <- as.data.frame(df)
    ggscatter(df,
        x = "sqv", y = "cr", size = 5, shape=21,
        add = "reg.line", conf.int = TRUE,
        cor.coef = TRUE, cor.method = "kendall",
        xlab = "SVQ UMI counts", ylab = "Cellranger UMI Counts",  cor.coeff.args = list(method = "kendall", size=6, label.x = 10, label.sep = "\n"),
        ggthemes=theme_minimal()) + font("xy.text", size = 14) + 
        font("xlab", size = 16) + font("ylab", size = 16) +
        scale_x_continuous(breaks = get_breaks(n = 10)) +
        scale_y_continuous(breaks = get_breaks(n = 10)) +
        rotate_x_text(angle = 90, hjust = 1, vjust = 0.5)
    ggsave(paste0(figname, ".svg"), dpi = 300, height = 4, width = 4, units = "in")
    ggsave(paste0(figname, ".png"), dpi = 300, height = 4, width = 4, units = "in")
    ggsave(paste0(figname, ".pdf"), dpi = 300, height = 4, width = 4, units = "in")

    message("STATUS: mean sqv ", length(rownames(gfpsqvhiv)), " cells ", mean(gfpsqvhiv[, "nCount_RNA"]))
    message("STATUS: mean cr ", length(rownames(gfpcrhiv)), " cells ", mean(gfpcrhiv[, "nCount_RNA"]))
    dfmelt <- melt(df)
    dfmelt$variable <- paste(dfmelt$variable, "GFP", sep = "_")
    # p <- ggplot(dfmelt2, aes(x = value, color = variable)) +
    #     geom_density()
    # ggsave(paste0(figname, "_density.svg"), dpi = 300, height = 4, width = 4, units = "in")
}

scatterplots4remaining <- function(overlap, SRR25cr, SRR25sqv, figname) {
    SRR25crhivov <- subset(SRR25cr, cells = overlap)
    zcrhivov <- SRR25crhivov["hiv-positive", ]
    gfpcrhivov <- zcrhivov@meta.data[zcrhivov@meta.data$nFeature_RNA != 0, ]

    SRR25sqvhivov <- subset(SRR25sqv, cells = overlap)
    zsqvhivov <- SRR25sqvhivov["HIV-1-pNL4-3", ]
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

########### --Processing HIV replicate SRR6825024 --###########

#-- Pull out intersecting cells to include in analysis (want to include only overlapping cells)
pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_CR, "SRR6825024",
    "outs", "filtered_feature_bc_matrix"))

pbmccr <- CreateSeuratObject(
    counts = pbmc.data,
    project = "HIV",
    min.cells = 0,
    min.features = 0
)

pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_SQV, "SRR6825024",
    "outs", "filtered_feature_bc_matrix_bowtie2_global"))

pbmcsqvbow2global <- CreateSeuratObject(
    counts = pbmc.data,
    project = "HIV",
    min.cells = 0,
    min.features = 0
)

pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_SQV, "SRR6825024",
    "outs", "filtered_feature_bc_matrix_bbmap_global"))

pbmcsqvbbmapglobal <- CreateSeuratObject(
    counts = pbmc.data,
    project = "HIV",
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
results_folder <- file.path(results_folder_base, "SRR24/")
generate_folder(results_folder)

SRR25cr <- initial_SC_Seurat_CR(file.path(SAMPLE_PATH_CR, "SRR6825024", "outs", "filtered_feature_bc_matrix"),
    sampleID = "HIV_REP_2", cells = cells
)

results_folder_base <- file.path("4.SC_Analysis", "SRR24/")
generate_folder(results_folder_base)

results_folder <- file.path(results_folder_base, "bowtie2_global/")
generate_folder(results_folder)
SRR25sqvbow2global <- initial_SC_Seurat_SQV(file.path(SAMPLE_PATH_SQV, "SRR6825024", "outs", "filtered_feature_bc_matrix_scPathoQuant_bowtie2"),
    sampleID = "HIV_REP_2", cells = cells
)

results_folder <- file.path(results_folder_base, "bbmap_global/")
generate_folder(results_folder)
SRR25sqvbbmapglobal <- initial_SC_Seurat_SQV(file.path(SAMPLE_PATH_SQV, "SRR6825024", "outs", "filtered_feature_bc_matrix_scPathoQuant_bbmap"),
    sampleID = "HIV_REP_2", cells = cells
)

zcr <- SRR25cr["EGFP-positive", ]
gfpcr <- zcr@meta.data[zcr@meta.data$nFeature_RNA != 0, ]

zsqvbow2global <- SRR25sqvbow2global["EGFP-positive", ]
gfpsqvbow2global <- zsqvbow2global@meta.data[zsqvbow2global@meta.data$nFeature_RNA != 0, ]
zsqvbbmapglobal <- SRR25sqvbbmapglobal["EGFP-positive", ]
gfpsqvbbmapglobal <- zsqvbbmapglobal@meta.data[zsqvbbmapglobal@meta.data$nFeature_RNA != 0, ]

message(
    "STATUS: Number of cells that are GFP and HIV positive based on CellRanger and sqv bowtie2 global ",
    length(intersect(rownames(gfpsqvbow2global), rownames(gfpcr))), " sqv=", length(rownames(gfpsqvbow2global)), " Cellranger=", length(rownames(gfpcr))
)
message(
    "STATUS: Number of cells that are GFP and HIV positive based on CellRanger and sqv bbmap global ",
    length(intersect(rownames(gfpsqvbbmapglobal), rownames(gfpcr))), " sqv=", length(rownames(gfpsqvbbmapglobal)), " Cellranger=", length(rownames(gfpcr))
)

scatterplots4GFP(gfpcr, gfpsqvbow2global[rownames(gfpcr), ], SRR25cr, SRR25sqvbow2global, file.path(results_folder_base, "4.CorrelationGFPHIVbow2global_SRR24test"))
scatterplots4GFP(gfpcr, gfpsqvbbmapglobal[rownames(gfpcr), ], SRR25cr, SRR25sqvbbmapglobal, file.path(results_folder_base, "4.CorrelationGFPHIVbbmapglobal_SRR24test"))

message(all.equal(rownames(gfpcr), rownames(gfpsqvbbmapglobal[rownames(gfpcr), ])))
message(all.equal(rownames(gfpsqvbow2global[rownames(gfpcr), ]), rownames(gfpsqvbbmapglobal[rownames(gfpcr), ])))

gfptable <- cbind(gfpcr[, "nCount_RNA"], gfpsqvbow2global[rownames(gfpcr), "nCount_RNA"], gfpsqvbbmapglobal[rownames(gfpcr), "nCount_RNA"])

rownames(gfptable) <- rownames(gfpcr)
colnames(gfptable) <- c("CellRanger", "SQV bowtie2", "SQV bbmap")
testRes <- cor.mtest(gfptable, method = "kendall", conf.level = 0.95)
M <- cor(gfptable, method = "kendall")

purporange <- colorRampPalette(brewer.pal(8, "PuOr"))(8)

png(file.path(results_folder_base, "4.CorrplotEGFPSR24.png"), width = 5.4, height = 4.5, units="in", res=300)
corrplot(M,
    type = "lower", order = "original", tl.col = "black",
    cl.ratio = 0.4, cl.cex = 1, col = purporange, # col = COL2("PuOr", 10),
    addCoef.col = "white", cl.pos = "r", number.cex=2,
    pch.cex = 1, pch.col = "white", diag = FALSE
)
dev.off()

svglite(file.path(results_folder_base, "4.CorrplotEGFPSR24.svg"), width = 5.4, height = 4.5)
corrplot(M,
    type = "lower", order = "original", tl.col = "black",
    cl.ratio = 0.4, cl.cex = 1, col = purporange, # col = COL2("PuOr", 10),
    addCoef.col = "white", cl.pos = "r", number.cex=2,
    pch.cex = 1, pch.col = "white", diag = FALSE
)
dev.off()
zcrhiv <- SRR25cr["hiv-positive", ]
gfpcrhivremain <- zcrhiv@meta.data[zcrhiv@meta.data$nFeature_RNA != 0, ]

zsqvhivbow2global <- SRR25sqvbow2global["HIV-1-pNL4-3", ]
gfpsqvhivremainzsqvhivbow2global <- zsqvhivbow2global@meta.data[zsqvhivbow2global@meta.data$nFeature_RNA != 0, ]
zsqvhivbbmapglobal <- SRR25sqvbbmapglobal["HIV-1-pNL4-3", ]
gfpsqvhivremainzsqvhivbbmapglobal <- zsqvhivbbmapglobal@meta.data[zsqvhivbbmapglobal@meta.data$nFeature_RNA != 0, ]

crcells <- setdiff(rownames(gfpcrhivremain), rownames(gfpcr))
sqvcellsbow2global <- setdiff(rownames(gfpsqvhivremainzsqvhivbow2global), rownames(gfpsqvbow2global))
sqvcellsbbmapglobal <- setdiff(rownames(gfpsqvhivremainzsqvhivbbmapglobal), rownames(gfpsqvbbmapglobal))

vd <- venn.diagram(
    x = list(crcells, sqvcellsbow2global),
    category.names = c("cellranger", "sqv bowtie2 global"), sub.cex = 20,
    filename = file.path(results_folder_base, "4.VD_SRR25_HIVpositivebow2global.png"),
    output = TRUE, imagetype = "png"
)
vd <- venn.diagram(
    x = list(crcells, sqvcellsbbmapglobal),
    category.names = c("cellranger", "sqv bbmap global"), sub.cex = 20,
    filename = file.path(results_folder_base, "4.VD_SRR25_HIVpositivebbmapglobal.png"),
    output = TRUE, imagetype = "png"
)

overlapbow2global <- intersect(sqvcellsbow2global, crcells)
overlapbbmapglobal <- intersect(sqvcellsbbmapglobal, crcells)

scatterplots4remaining(overlapbow2global, SRR25cr, 
    SRR25sqvbow2global,
    file.path(results_folder_base, "4.Correlationbow2globalRemainSRR24"))
scatterplots4remaining(overlapbbmapglobal, SRR25cr,
    SRR25sqvbbmapglobal,
    file.path(results_folder_base, "4.CorrelationbbbmapglobalRemainSRR24"))

# -Upset Plot
library(ComplexHeatmap)

# Remaining upset plot
remaining <- list(
    "CellRanger" = rownames(gfpcrhivremain),
    "SVQ bowtie2" = rownames(gfpsqvhivremainzsqvhivbow2global),
    "SVQ bbmap" = rownames(gfpsqvhivremainzsqvhivbbmapglobal)
)

# m <- list_to_matrix(remaining)
m <- make_comb_mat(remaining)
cs <- comb_size(m)
ss <- set_size(m)

ht <- UpSet(m,
    set_order = c("CellRanger", "SVQ bowtie2", "SVQ bbmap"),
    lwd = 2, 
    top_annotation = HeatmapAnnotation(
        "Intersection of cells\ndetected with HIV" = anno_barplot(cs,
            ylim = c(0, max(cs) * 1.1),
            border = FALSE,
            gp = gpar(fill = "black"),
            axis = TRUE, 
            height = unit(4, "cm"),
            axis_param = list(gp = gpar(fontsize = 12))
        ), 
        annotation_name_side = "left",
        annotation_name_rot = 90
    ),
    right_annotation = rowAnnotation(
        "HIV + cells per\nalignment method" = anno_barplot(-ss,
            baseline = 0,
            axis_param = list(
                at = c(0, -2000, -4000, -6000, -8000),
                labels = c(0, 2000, 4000, 6000, 8000),
                labels_rot =90, gp = gpar(fontsize = 12,fill = "black")
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
# png(file.path(results_folder_base, "Upsetplotv5.png"), width = 5.5, height = 3.4, units = "in", res=500)
# svglite(file.path(results_folder_base, "Upsetplotv2.svg"), width = 5.5, height = 3.4)
pdf(file.path(results_folder_base, "Upsetplot.pdf"), width = 5.5, height = 3.4)
ht <- draw(ht)
od <- column_order(ht)
decorate_annotation("Intersection of cells\ndetected with HIV", {
    grid.text(cs[od],
        x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
        default.units = "native", just = c("left", "bottom"),
        gp = gpar(fontsize = 14, col = "#404040"), rot = 0
    )
})
dev.off()

# png(file.path(results_folder_base, "Upsetplotv5.png"), width = 5.5, height = 3.4, units = "in", res=500)
svglite(file.path(results_folder_base, "Upsetplot.svg"), width = 5.5, height = 3.4)
# pdf(file.path(results_folder_base, "Upsetplot.pdf"), width = 5.5, height = 3.4)
ht <- draw(ht)
od <- column_order(ht)
decorate_annotation("Intersection of cells\ndetected with HIV", {
    grid.text(cs[od],
        x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
        default.units = "native", just = c("left", "bottom"),
        gp = gpar(fontsize = 14, col = "#404040"), rot = 0
    )
})
dev.off()

png(file.path(results_folder_base, "Upsetplot.png"), width = 5.5, height = 3.4, units = "in", res=500)
# svglite(file.path(results_folder_base, "Upsetplotv2.svg"), width = 5.5, height = 3.4)
# pdf(file.path(results_folder_base, "Upsetplot.pdf"), width = 5.5, height = 3.4)
ht <- draw(ht)
od <- column_order(ht)
decorate_annotation("Intersection of cells\ndetected with HIV", {
    grid.text(cs[od],
        x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
        default.units = "native", just = c("left", "bottom"),
        gp = gpar(fontsize = 14, col = "#404040"), rot = 0
    )
})
dev.off()
bc <- read.csv(
    "/share/lwhitmo/projects/scViralQuantTesting/CD4_HIV_DukeData/scViralQuantResuls_bbmap_global/SRR6825024/_tmp/barcode_umi_read_table.csv",
    header = TRUE
)
bowtie_bbmap <- intersect(rownames(gfpsqvhivremainzsqvhivbbmapglobal), rownames(gfpsqvhivremainzsqvhivbow2global))
bowtie_bbmap <- setdiff(bowtie_bbmap, rownames(gfpcrhivremain))
bctmp <- bc[bc$cell_barcode %in% bowtie_bbmap, ]
write.table(
    bctmp, file.path(results_folder_base, "bowtiebbmap_overlaping_positive.csv"),
    quote = FALSE, col.names = FALSE, row.names = FALSE
)

bbmap_cr <- intersect(
    rownames(gfpsqvhivremainzsqvhivbbmapglobal), rownames(gfpcrhivremain)
)

bbmap_cr <- setdiff(bbmap_cr, rownames(gfpsqvhivremainzsqvhivbow2global))

bctmp <- bc[bc$cell_barcode %in% bbmap_cr, ]

write.table(
    bctmp, file.path(results_folder_base, "bbmap_cr_overlaping_positive.csv"),
    quote = FALSE, col.names = FALSE, row.names = FALSE
)

cronly <- setdiff(
    rownames(gfpcrhivremain),
    c(rownames(gfpsqvhivremainzsqvhivbow2global), rownames(gfpsqvhivremainzsqvhivbbmapglobal))
)
bctmp <- bc[bc$cell_barcode %in% cronly, ]

write.table(
    bctmp, file.path(results_folder_base, "cellranger_overlaping_positive.csv"),
    quote = FALSE, col.names = FALSE, row.names = FALSE
)

bbmaponly <- setdiff(
    rownames(gfpsqvhivremainzsqvhivbbmapglobal),
    c(rownames(gfpsqvhivremainzsqvhivbow2global), rownames(gfpcrhivremain))
)
bctmp <- bc[bc$cell_barcode %in% bbmaponly, ]

write.table(
    bctmp, file.path(results_folder_base, "bbmap_overlaping_positive.csv"),
    quote = FALSE, col.names = FALSE, row.names = FALSE
)
########### --Processing HIV replicate SRR6825025 --###########
## Pull out intersecting cells to include in analysis (want to include only overlapping cells) 
pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_CR, "SRR6825025","outs", "filtered_feature_bc_matrix"))

pbmccr <- CreateSeuratObject(
    counts = pbmc.data,
    project = "HIV",
    min.cells = 0,
    min.features = 0
)

pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_SQV, "SRR6825025", "outs", "filtered_feature_bc_matrix_scPathoQuant_bowtie2"))

pbmcsqvbow2global <- CreateSeuratObject(
    counts = pbmc.data,
    project = "HIV",
    min.cells = 0,
    min.features = 0
)

pbmc.data <- Read10X(data.dir = file.path(SAMPLE_PATH_SQV, "SRR6825025", "outs", "filtered_feature_bc_matrix_scPathoQuant_bbmap"))

pbmcsqvbbmapglobal <- CreateSeuratObject(
    counts = pbmc.data,
    project = "HIV",
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
results_folder <- file.path(results_folder_base, "SRR25/")
generate_folder(results_folder)

SRR25cr <- initial_SC_Seurat_CR(file.path(SAMPLE_PATH_CR, "SRR6825025","outs", "filtered_feature_bc_matrix"),
    sampleID = "HIV_REP_2", cells = cells
)

results_folder_base <- file.path("4.SC_Analysis", "SRR25/")
generate_folder(results_folder_base)

results_folder <- file.path(results_folder_base, "bowtie2_global/")
generate_folder(results_folder)
SRR25sqvbow2global <- initial_SC_Seurat_SQV(file.path(SAMPLE_PATH_SQV, "SRR6825025", "outs", "filtered_feature_bc_matrix_scPathoQuant_bowtie2"),
    sampleID = "HIV_REP_2", cells = cells
)

results_folder <- file.path(results_folder_base, "bbmap_global/")
generate_folder(results_folder)
SRR25sqvbbmapglobal <- initial_SC_Seurat_SQV(file.path(SAMPLE_PATH_SQV, "SRR6825025", "outs", "filtered_feature_bc_matrix_scPathoQuant_bbmap"),
    sampleID = "HIV_REP_2", cells = cells
)

zcr <- SRR25cr["EGFP-positive", ]
gfpcr <- zcr@meta.data[zcr@meta.data$nFeature_RNA != 0, ]

zsqvbow2global <- SRR25sqvbow2global["EGFP-positive", ]
gfpsqvbow2global <- zsqvbow2global@meta.data[zsqvbow2global@meta.data$nFeature_RNA != 0, ]
zsqvbbmapglobal <- SRR25sqvbbmapglobal["EGFP-positive", ]
gfpsqvbbmapglobal <- zsqvbbmapglobal@meta.data[zsqvbbmapglobal@meta.data$nFeature_RNA != 0, ]

message(
    "STATUS: Number of cells that are GFP and HIV positive based on CellRanger and sqv bowtie2 global ",
    length(intersect(rownames(gfpsqvbow2global), rownames(gfpcr))), " sqv=", length(rownames(gfpsqvbow2global)), " Cellranger=", length(rownames(gfpcr))
)
message(
    "STATUS: Number of cells that are GFP and HIV positive based on CellRanger and sqv bbmap global ",
    length(intersect(rownames(gfpsqvbbmapglobal), rownames(gfpcr))), " sqv=", length(rownames(gfpsqvbbmapglobal)), " Cellranger=", length(rownames(gfpcr))
)

scatterplots4GFP(gfpcr, gfpsqvbow2global, SRR25cr, SRR25sqvbow2global, file.path(results_folder_base, "4.CorrelationGFPHIVbow2global_SRR25"))
scatterplots4GFP(gfpcr, gfpsqvbbmapglobal, SRR25cr, SRR25sqvbbmapglobal, file.path(results_folder_base, "4.CorrelationGFPHIVbbmapglobal_SRR25"))

message(all.equal(rownames(gfpcr), rownames(gfpsqvbbmapglobal)))
message(all.equal(rownames(gfpsqvbow2global), rownames(gfpsqvbbmapglobal)))

gfptable <- cbind(gfpcr[, "nCount_RNA"], gfpsqvbow2global[, "nCount_RNA"], gfpsqvbbmapglobal[, "nCount_RNA"])

rownames(gfptable) <-rownames(gfpsqvbbmapglobal)
colnames(gfptable) <- c("CellRanger","SQV bowtie2", "SQV bbmap")
testRes <- cor.mtest(gfptable,method="kendall", conf.level = 0.95)
M <- cor(gfptable,method="kendall")
purporange <- colorRampPalette(brewer.pal(8, "PuOr"))(8)

svglite(file.path(results_folder_base, "4.CorrplotEGFP25.svg"), width = 5.4, height = 4.5)
corrplot(M,
    type = "lower", order = "original", tl.col = "black",
    cl.ratio = 0.4, cl.cex = 1, col = purporange, # col = COL2("PuOr", 10),
    addCoef.col = "white", cl.pos = "r", number.cex=2,
    pch.cex = 1, pch.col = "white", diag = FALSE
)
dev.off()

png(file.path(results_folder_base, "4.CorrplotEGFP25.png"), width = 5.4, height = 4.5, units="in", res=300)
corrplot(M,
    type = "lower", order = "original", tl.col = "black",
    cl.ratio = 0.4, cl.cex = 1, col = purporange, # col = COL2("PuOr", 10),
    addCoef.col = "white", cl.pos = "r",  number.cex=2,
    pch.cex = 1, pch.col = "white", diag = FALSE
)
dev.off()
zcrhiv <- SRR25cr["hiv-positive", ]
gfpcrhivremain <- zcrhiv@meta.data[zcrhiv@meta.data$nFeature_RNA != 0, ]

zsqvhivbow2global <- SRR25sqvbow2global["HIV-1-pNL4-3", ]
gfpsqvhivremainzsqvhivbow2global <- zsqvhivbow2global@meta.data[zsqvhivbow2global@meta.data$nFeature_RNA != 0, ]
zsqvhivbbmapglobal <- SRR25sqvbbmapglobal["HIV-1-pNL4-3", ]
gfpsqvhivremainzsqvhivbbmapglobal <- zsqvhivbbmapglobal@meta.data[zsqvhivbbmapglobal@meta.data$nFeature_RNA != 0, ]

crcells <- setdiff(rownames(gfpcrhivremain), rownames(gfpcr))
sqvcellsbow2global <- setdiff(rownames(gfpsqvhivremainzsqvhivbow2global), rownames(gfpsqvbow2global))
sqvcellsbbmapglobal <- setdiff(rownames(gfpsqvhivremainzsqvhivbbmapglobal), rownames(gfpsqvbbmapglobal))

vd <- venn.diagram(
    x = list(crcells, sqvcellsbow2global),
    category.names = c("cellranger", "sqv bowtie2 global"), sub.cex = 20,
    filename = file.path(results_folder_base, "4.VD_SRR25_HIVpositivebow2global.png"),
    output = TRUE, imagetype = "png"
)
vd <- venn.diagram(
    x = list(crcells, sqvcellsbbmapglobal),
    category.names = c("cellranger", "sqv bbmap global"), sub.cex = 20,
    filename = file.path(results_folder_base, "4.VD_SRR25_HIVpositivebbmapglobal.png"),
    output = TRUE, imagetype = "png"
)

overlapbow2global <- intersect(sqvcellsbow2global, crcells)
overlapbbmapglobal <- intersect(sqvcellsbbmapglobal, crcells)

scatterplots4remaining(overlapbow2global, SRR25cr, SRR25sqvbow2global, file.path(results_folder_base, "4.Correlationbow2globalRemainSRR25"))
scatterplots4remaining(overlapbbmapglobal, SRR25cr, SRR25sqvbbmapglobal, file.path(results_folder_base, "4.CorrelationbbbmapglobalRemainSRR25"))

# --Upset Plot
library(ComplexHeatmap)

# Remaining upset plot 
remaining <- list(
    "CellRanger" = rownames(gfpcrhivremain),
    "SVQ bowtie2" = rownames(gfpsqvhivremainzsqvhivbow2global),
    "SVQ bbmap" = rownames(gfpsqvhivremainzsqvhivbbmapglobal)
)

# m <- list_to_matrix(remaining)
m <- make_comb_mat(remaining)
cs <- comb_size(m)
ss <- set_size(m)

ht <- UpSet(m,
 set_order = c("CellRanger", "SVQ bowtie2", "SVQ bbmap"),
 lwd=2, top_annotation = HeatmapAnnotation(
        "Intersection of cells\ndetected with HIV" = anno_barplot(cs, 
            ylim = c(0, max(cs)*1.1),
            border = FALSE, 
            gp = gpar(fill = "black", fontsize=14), 
            height = unit(4, "cm"),
            axis_param = list(gp = gpar(fontsize = 12))
        ),
        annotation_name_side = "left", 
        annotation_name_rot = 90),
        right_annotation = rowAnnotation(
        "HIV + cells per\nalignment method" = anno_barplot(-ss, 
            baseline = 0,
            axis_param = list(
                at = c(0, -200, -400, -600),
                labels = c(0, 200, 400, 600),
                labels_rot = 90, gp = gpar(fontsize=12)),
            border = FALSE, 
            gp = gpar(fontsize = 14,fill = "black"),
            width = unit(2, "cm")
         ), set_name = anno_text(set_name(m), 
             location = 0.05, 
             just = "left", gp = gpar(fontsize=14),
             width = max_text_width(set_name(m)) + unit(4, "mm"))
            ),
    left_annotation = NULL, show_row_names = FALSE
)

# png(file.path(results_folder_base, "Upsetplot.png"), width = 5.5, height = 3.4, units = "in", res=300)
svglite(file.path(results_folder_base, "Upsetplot.svg"), width = 5.5, height = 3.5)
# pdf(file.path(results_folder_base, "Upsetplot.pdf"), width = 5.5, height = 3.4)
ht <- draw(ht)
od <- column_order(ht)
decorate_annotation("Intersection of cells\ndetected with HIV", {
    grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
        default.units = "native", just = c("left", "bottom"), 
        gp = gpar(fontsize = 14, col = "#404040"), rot = 0)})
dev.off()

png(file.path(results_folder_base, "Upsetplot.png"), width = 5.5, height = 3.4, units = "in", res=300)
# svglite(file.path(results_folder_base, "Upsetplot.svg"), width = 5.5, height = 3.5)
# pdf(file.path(results_folder_base, "Upsetplot.pdf"), width = 5.5, height = 3.4)
ht <- draw(ht)
od <- column_order(ht)
decorate_annotation("Intersection of cells\ndetected with HIV", {
    grid.text(cs[od],
        x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
        default.units = "native", just = c("left", "bottom"),
        gp = gpar(fontsize = 14, col = "#404040"), rot = 0
    )
})
dev.off()

# png(file.path(results_folder_base, "Upsetplot.png"), width = 5.5, height = 3.4, units = "in", res = 300)
# svglite(file.path(results_folder_base, "Upsetplot.svg"), width = 5.5, height = 3.5)
pdf(file.path(results_folder_base, "Upsetplot.pdf"), width = 5.5, height = 3.4)
ht <- draw(ht)
od <- column_order(ht)
decorate_annotation("Intersection of cells\ndetected with HIV", {
    grid.text(cs[od],
        x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"),
        default.units = "native", just = c("left", "bottom"),
        gp = gpar(fontsize = 14, col = "#404040"), rot = 0
    )
})
dev.off()

bc <- read.csv(
    "/share/lwhitmo/projects/scViralQuantTesting/CD4_HIV_DukeData/scViralQuantResuls_bbmap_global/SRR6825025/_tmp/barcode_umi_read_table.csv",
    header = TRUE
)
bowtie_bbmap <- intersect(rownames(gfpsqvhivremainzsqvhivbbmapglobal), rownames(gfpsqvhivremainzsqvhivbow2global))
bowtie_bbmap <- setdiff(bowtie_bbmap, rownames(gfpcrhivremain))
bctmp <- bc[bc$cell_barcode %in% bowtie_bbmap, ]
write.table(
    bctmp, file.path(results_folder_base, "bowtiebbmap_overlaping_positive.csv"),
    quote = FALSE, col.names = FALSE, row.names = FALSE
)

bbmap_cr <- intersect(
    rownames(gfpsqvhivremainzsqvhivbbmapglobal), rownames(gfpcrhivremain)
)

bbmap_cr <- setdiff(bbmap_cr, rownames(gfpsqvhivremainzsqvhivbow2global))

bctmp <- bc[bc$cell_barcode %in% bbmap_cr, ]

write.table(
    bctmp, file.path(results_folder_base, "bbmap_cr_overlaping_positive.csv"),
    quote = FALSE, col.names = FALSE, row.names = FALSE
)

cronly <- setdiff(
    rownames(gfpcrhivremain),
    c(rownames(gfpsqvhivremainzsqvhivbow2global), rownames(gfpsqvhivremainzsqvhivbbmapglobal))
)
bctmp <- bc[bc$cell_barcode %in% cronly, ]

write.table(
    bctmp, file.path(results_folder_base, "cellranger_overlaping_positive.csv"),
    quote = FALSE, col.names = FALSE, row.names = FALSE
)

bbmaponly <- setdiff(
    rownames(gfpsqvhivremainzsqvhivbbmapglobal),
    c(rownames(gfpsqvhivremainzsqvhivbow2global), rownames(gfpcrhivremain))
)
bctmp <- bc[bc$cell_barcode %in% bbmaponly, ]

write.table(
    bctmp, file.path(results_folder_base, "bbmap_overlaping_positive.csv"),
    quote = FALSE, col.names = FALSE, row.names = FALSE
)
