# overall_water_climate_analysis.R performs rank-based correlation on the overall data for Le et al. (2023)
# Author(s): Edward Le
# Statistical/Code Validation By: Joseph Janssen

library(tidyverse)
library(reshape2)
library(ggcorrplot)
library(testit)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(hash)

plotTitlesHash <- hash()
plotTitlesHash[["mean_sp"]] <- "SP"
plotTitlesHash[["mean_ai"]] <- "AI"
plotTitlesHash[["mean_si"]] <- "SI"
plotTitlesHash[["bfi"]] <- "BFI"
plotTitlesHash[["q5Frac"]] <- expression("Normalized" ~ Q[5])
plotTitlesHash[["low_fdc"]] <- "Low-FDC"
plotTitlesHash[["mean_low_flow_dur"]] <- "Low-flow duration"
plotTitlesHash[["q95Frac"]] <- expression("Normalized" ~ Q[95])
plotTitlesHash[["high_fdc"]] <- "High-FDC"
plotTitlesHash[["mean_high_flow_dur"]] <- "High-flow duration"

masterDir <- "../"
dataDir <- paste(masterDir, "Data/", sep = "")
plotDir <- paste(masterDir, "Plots/", sep = "")
spearmanPlotDir <- paste(plotDir, "Spearman Correlational Analysis/", sep = "")

dir.create(plotDir, showWarnings = FALSE)
dir.create(spearmanPlotDir, showWarnings = FALSE)

# Data name (presumes data folder already exists)
combinedDataName <- paste(dataDir, "overall_water_climate_data.csv", sep = "")


# Read and seperate out data
combinedData <- read_csv(combinedDataName)
combinedData <- subset(combinedData, select = -c(catchment, longitude, latitude))

# Common constants
tag <- ""
sigLevel <- 0.05
hiddenColumns <- c("mean_ai", "mean_si")
focalVariableColumn <- "mean_sp"

aridInPhaseData <- combinedData[which(combinedData$mean_ai >= quantile(combinedData$mean_ai, 0.75) & combinedData$mean_si >= quantile(combinedData$mean_si, 0.75)), ]
aridOutOfPhaseData <- combinedData[which(combinedData$mean_ai >= quantile(combinedData$mean_ai, 0.75) & combinedData$mean_si <= quantile(combinedData$mean_si, 0.25)), ]
wetOutOfPhaseData <- combinedData[which(combinedData$mean_ai <= quantile(combinedData$mean_ai, 0.25) & combinedData$mean_si <= quantile(combinedData$mean_si, 0.25)), ]

# Input: colNames (vector of column names), columnLabelsHash (hashmap of labels)
# Output: Vector of labels
get_labels <- function(colNames, columnLabelsHash) {
    output <- c()
    for (currName in colNames) {
        output <- c(output, columnLabelsHash[[currName]])
    }
    return(output)
}

# Input: colNames (list of dataframes), columnLabelsHash (list of labels)
# Output: List of labels with sample sizes
append_sample_sizes_to_column_labels <- function(dataList, labelsList) {
    assert(length(dataList) == length(labelsList))
    appendedLabelsList <- list()
    for (i in seq_len(length(dataList))) {
        appendedLabelsList <- append(appendedLabelsList, paste(labelsList[[i]], "\n(n = ", nrow(dataList[[i]]), ")", sep = ""))
    }
    return(appendedLabelsList)
}

# Input: data (dataframe), correlations (dataframe), focalVariableColumn (string)
# Output: Dataframe of scaled values
scale_data <- function(data, correlations, focalVariableColumn) {
    focalVariableStdev <- sd(data[[focalVariableColumn]])
    scaledData <- data.frame()
    colsToLoop <- colnames(data)[colnames(data) != focalVariableColumn]
    for (currCol in colsToLoop) {
        currStdev <- sd(data[[currCol]])
        assert(is.numeric(correlations$value[correlations$Var1 == currCol]))
        correlations$value[correlations$Var1 == currCol] <- correlations$value[correlations$Var1 == currCol] * (currStdev / focalVariableStdev)
        scaledData <- bind_rows(scaledData, correlations[correlations$Var1 == currCol, ])
    }
    return(scaledData)
}

# Input: dataList (list of dataframes), titleList (list of titles), variablesToUse (vector of variables),
# focalVariableColumn (string), scaleData (boolean)
# Output: dataframe of correlations
get_correl <- function(dataList, titleList, variablesToUse, focalVariableColumn, scaleData) {
    correlData <- data.frame()
    for (i in seq_len(length(dataList))) {
        currData <- dataList[[i]][variablesToUse]
        currCorr <- cor(currData, method = "spearman")
        currCorr <- round(currCorr, 2)
        currCorr <- melt(currCorr)
        currCorr <- currCorr[currCorr$Var2 == focalVariableColumn, ]
        currCorr <- currCorr[!currCorr$Var1 == focalVariableColumn, ]
        currCorr$Var2 <- titleList[[i]]
        if (scaleData) {
            currCorr <- scale_data(currData, currCorr, focalVariableColumn)
        }
        correlData <- bind_rows(correlData, currCorr)
    }
    return(correlData)
}

# Input: dataList (list of dataframes), titleList (list of titles), variablesToUse (vector of variables),
# focalVariableColumn (string)
# Output: Dataframe with insignificant correlations
get_insig_correls <- function(dataList, titleList, variablesToUse, focalVariableColumn) {
    pValueDataList <- data.frame()
    for (i in seq_len(length(dataList))) {
        currData <- dataList[[i]][variablesToUse]
        currPValues <- round(cor_pmat(currData, method = "spearman"), 2)
        currPValues <- melt(currPValues)
        currPValues <- currPValues[currPValues$Var2 == focalVariableColumn, ]
        currPValues <- currPValues[!currPValues$Var1 == focalVariableColumn, ]
        currPValues$Var2 <- titleList[[i]]
        currPValues$value[currPValues$value >= sigLevel] <- TRUE
        currPValues$value[currPValues$value < sigLevel] <- NA
        currPValues <- drop_na(currPValues)
        pValueDataList <- bind_rows(pValueDataList, currPValues)
    }
    return(pValueDataList)
}

# Input: dataList (list of dataframes), titleList (list of titles), variablesToUse (vector of variables),
# focalVariableColumn (string), columnLabelsHash (hashmap of labels), plotTitle (string), scaleData (boolean)
# Output: Plot of Spearman correlations
plot_correl_rank <- function(dataList, titleList, variablesToUse, focalVariableColumn, columnLabelsHash, plotTitle, scaleData) {
    correlData <- get_correl(dataList, titleList, variablesToUse, focalVariableColumn, scaleData)

    yAxisOrder <- rev(variablesToUse[variablesToUse != focalVariableColumn])
    yAxisLabels <- get_labels(yAxisOrder, columnLabelsHash)

    # Plot data
    plot <- NULL
    if (scaleData) {
        fillLabel <- "Scaled Spearman\nCorrelation"
        plot <- ggplot(data = correlData, aes(Var2, Var1, fill = value >= 0)) +
            geom_tile(color = "white") +
            scale_fill_manual(
                values = c("pink", "lightblue"),
                labels = c("Negative Scaled\nCorrelation", "Positive Scaled\nCorrelation"),
                name = fillLabel
            ) +
            geom_text(aes(Var2, Var1, label = round(value, 2)), color = "black", size = 4)
    } else {
        pValueData <- get_insig_correls(dataList, titleList, variablesToUse, focalVariableColumn)
        fillLabel <- expression("Spearman Correlation (ρ)")
        plot <- ggplot(data = correlData, aes(Var2, Var1, fill = value)) +
            geom_tile(color = "white") +
            scale_fill_gradient2(high = "blue", mid = "white", low = "red", midpoint = 0, limits = c(-1, 1), name = fillLabel) +
            geom_text(aes(Var2, Var1, label = round(value, 2)), color = "black", size = 4) +
            geom_point(data = pValueData, shape = 4, aes(Var2, Var1), size = 20)
    }
    plot <- plot +
        scale_y_discrete(limits = yAxisOrder, labels = yAxisLabels) +
        ggtitle(plotTitle) +
        theme_classic() +
        theme(
            axis.title = element_blank(),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = element_text(
                angle = 0,
                hjust = 0.5
            ),
            legend.position = "bottom",
            legend.title.align = 0.5,
            legend.text.align = 0.5
        )

    return(plot)
}

# Input: dataList (list of dataframes), titleList (list of titles), variablesToUse (vector of variables),
# focalVariableColumn (string), columnLabelsHash (hashmap of labels), filename (string)
# Output: Saves plot of combined Spearman correlations
create_correl_plot_panels <- function(dataList, titleList, variablesToUse, focalVariableColumn, columnLabelsHash, filename) {
    combinedPlot <- plot_grid(
        plot_correl_rank(dataList, titleList, variablesToUse, focalVariableColumn, columnLabelsHash, "(A)", FALSE),
        # Remove all data for panel (b)
        plot_correl_rank(dataList[2:length(dataList)], titleList[2:length(dataList)], variablesToUse, focalVariableColumn, columnLabelsHash, "(B)", TRUE)
    )
    ggsave(filename, combinedPlot, height = 8, width = 14)
}


dataList <- list(combinedData, aridInPhaseData, aridOutOfPhaseData, wetOutOfPhaseData)
labelsList <- list("All\nClimates", "Arid\nIn-Phase", "Arid\nOut-of-Phase", "Wet\nOut-of-Phase")
labelsList <- append_sample_sizes_to_column_labels(dataList, labelsList)

# Create Low Flow Plots
variablesToUse <- c("mean_sp", "bfi", "q5Frac", "low_fdc", "mean_low_flow_dur")
lowFlowFileName <- paste(spearmanPlotDir, "fig02_low_flow_spearman.jpeg")
create_correl_plot_panels(dataList, labelsList, variablesToUse, focalVariableColumn, plotTitlesHash, lowFlowFileName)

# Create High Flow Plots
variablesToUse <- c("mean_sp", "q95Frac", "high_fdc", "mean_high_flow_dur")
highFlowFileName <- paste(spearmanPlotDir, "fig03_high_flow_spearman.jpeg")
create_correl_plot_panels(dataList, labelsList, variablesToUse, focalVariableColumn, plotTitlesHash, highFlowFileName)