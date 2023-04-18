# annual_water_climate_analysis.R analyzes annual data for Le et al. (2023)
# Author(s): Edward Le, John Hammond
# Statistical/Code Validation By: Joseph Janssen

library(tidyverse)
library(ggplot2)
library(rnaturalearthhires)
library(rgeos)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(testit)
library(hash)
library(RColorBrewer)
library(stringr)

# Global variables
spRichnessThresh <- 0.5
smallestNumberOfYearsInACatchment <- 99999

plotTitlesHash <- hash()
plotTitlesHash[["sp"]] <- "SP"
plotTitlesHash[["ai"]] <- "AI"
plotTitlesHash[["ln(ai)"]] <- "Ln(AI)"
plotTitlesHash[["si"]] <- "SI"
plotTitlesHash[["bfi"]] <- "BFI"
plotTitlesHash[["q5Frac"]] <- expression("Normalized" ~ Q[5], "")
plotTitlesHash[["low_fdc"]] <- "Low-FDC"
plotTitlesHash[["mean_low_flow_dur"]] <- "Low-flow duration (days)"
plotTitlesHash[["q95Frac"]] <- expression("Normalized" ~ Q[95], "")
plotTitlesHash[["high_fdc"]] <- "High-FDC"
plotTitlesHash[["mean_high_flow_dur"]] <- "High-flow duration (days)"

stdPlotTitlesHash <- hash()
stdPlotTitlesHash[["sp"]] <- "Std. SP"
stdPlotTitlesHash[["ai"]] <- "Std. AI"
stdPlotTitlesHash[["ln(ai)"]] <- "Std. Ln(AI)"
stdPlotTitlesHash[["si"]] <- "Std. SI"
stdPlotTitlesHash[["bfi"]] <- "Std. BFI"
stdPlotTitlesHash[["q5Frac"]] <- expression("Std. normalized" ~ Q[5], "")
stdPlotTitlesHash[["low_fdc"]] <- "Std. Low-FDC"
stdPlotTitlesHash[["mean_low_flow_dur"]] <- "Std. low-flow duration"
stdPlotTitlesHash[["q95Frac"]] <- expression("Std. normalized" ~ Q[95], "")
stdPlotTitlesHash[["high_fdc"]] <- "Std. High-FDC"
stdPlotTitlesHash[["mean_high_flow_dur"]] <- "Std. high-flow duration"

geoPlotTitlesHash <- hash()
geoPlotTitlesHash[["sp"]] <- "SP"
geoPlotTitlesHash[["bfi"]] <- "(A) BFI"
geoPlotTitlesHash[["q5Frac"]] <- expression("(B) Normalized" ~ Q[5], "")
geoPlotTitlesHash[["low_fdc"]] <- "(C) Low-FDC"
geoPlotTitlesHash[["mean_low_flow_dur"]] <- "(D) Low-flow duration"
geoPlotTitlesHash[["q95Frac"]] <- expression("(E) Normalized" ~ Q[95], "")
geoPlotTitlesHash[["high_fdc"]] <- "(F) High-FDC"
geoPlotTitlesHash[["mean_high_flow_dur"]] <- "(G) High-flow duration"

# Counters
catchmentsLt15Years <- c()
catchmentsSnowRich <- c()
catchmentsSnowPoor <- c()

# Initialize dataframes
allYearsNoStandardizationDf <- data.frame()
allYearsSnowRichNoStandardizationDf <- data.frame()
allYearsSnowPoorNoStandardizationDf <- data.frame()
allYearsSnowRichStandardizedDf <- data.frame()
allYearsSnowPoorStandardizedDf <- data.frame()

predictors <- c("sp", "ai", "si")
predictorsWithLnAI <- c("sp", "ai", "ln(ai)", "si")
spPredictorVector <- c("sp")
responses <- c("bfi", "q5Frac", "low_fdc", "mean_low_flow_dur", "q95Frac", "high_fdc", "mean_high_flow_dur")

aiLabels <- c("Wetter", "Drier")
siLabels <- c("More Out-of-Phase", "More In-Phase")

# Column Names
spColName <- "sp"
aridityColName <- "ai"
seasonalityColName <- "si"

# Directory management
# Setwd (if necessary; depends on IDE)
# setwd(<fill_in>)
masterDir <- "../"
dataDir <- paste(masterDir, "Data/", sep = "")
annualDataDir <- paste(dataDir, "Annual Data/", sep = "")
statisticsDir <- paste(masterDir, "Statistics/", sep = "")

plotDir <- paste(masterDir, "Plots/", sep = "")
correlationalPlotDir <- paste(plotDir, "Pearson Correlational Analysis/", sep = "")

dir.create(plotDir, showWarnings = FALSE)
dir.create(correlationalPlotDir, showWarnings = FALSE)

overallCorrelationRowLabels <- c(
    "(A)", "", "", "",
    "(B)", "", "", "",
    "(C)", "", "", "",
    "(D)", "", "", "",
    "(E)", "", "", "",
    "(F)", "", "", "",
    "(G)", "", "", ""
)
temporalCorrelationRowLabelsAi <- c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)", "(G)")
temporalCorrelationRowLabelsSi <- c("", "", "", "", "", "", "")
annualDataPattern <- "_annual_water_climate_data.csv"
catchmentNumberPattern <- "([0-9]+)_"

overallDataFilename <- "overall_water_climate_data.csv"

# Input:  df (dataframe), predictor (string), response, methodString (one of "pearson", "kendall", "spearman"), titlesHash (hashmap of titles)
# splitColumn (which column to split), splitColumnLabels (label for that column), rowLabel (label per row)
# Output: Correlation plot (combined annual data)
correl_single_predictor_single_response <- function(df, predictor, response, methodString, titlesHash, splitCol = NULL, splitColLabels = NULL, rowLabel = NULL) {
    # For handling correlation coefficient of Ln(AI) plot
    if (predictor == "ln(ai)" && is.null(splitCol)) {
        df[predictor] <- log(df$ai)
        currPlot <- ggscatter(df,
            x = predictor, y = response,
            add = "reg.line",
            cor.coef = TRUE, cor.method = methodString,
            cor.coeff.args = list(
                mapping = aes(label = ..r.label..), digits = 2,
                label.x.npc = "left", label.y.npc = "top"
            ),
            alpha = 0.005
        ) +
            xlab(titlesHash[[predictor]]) +
            ylab(titlesHash[[response]])
    }
    # Plotting standard correlational analysis
    else if (is.null(splitCol)) {
        currPlot <- ggscatter(df,
            x = predictor, y = response,
            add = "reg.line",
            cor.coef = TRUE, cor.method = methodString,
            cor.coeff.args = list(
                mapping = aes(label = ..r.label..), digits = 2,
                label.x.npc = "left", label.y.npc = "top"
            ),
            alpha = 0.005
        ) +
            xlab(titlesHash[[predictor]]) +
            ylab(titlesHash[[response]])
    }
    # Plotting standardized correlational analysis; seperated by group
    else {
        # SplitCol splits the selected column based on some critical point
        criticalColPoint <- 0
        splitColCat <- paste(splitCol, "Cat", sep = "")
        df[splitColCat] <- ""
        df[df[splitCol] < criticalColPoint, splitColCat] <- splitColLabels[1]
        df[df[splitCol] >= criticalColPoint, splitColCat] <- splitColLabels[2]
        df[splitColCat] <- factor(df[[splitColCat]])

        currPlot <- ggscatter(df,
            x = predictor, y = response,
            color = splitColCat,
            palette = c("red", "blue"),
            add = "reg.line",
            cor.coef = TRUE, cor.method = methodString,
            cor.coeff.args = list(
                mapping = aes_string(
                    color = splitColCat,
                    label = quote(..r.label..)
                ),
                show.legend = FALSE,
                digits = 2,
                label.x.npc = "left",
                label.y.npc = "top"
            ),
            alpha = 0.1
        ) +
            xlab(titlesHash[[predictor]]) +
            ylab(titlesHash[[response]])
        currPlot <- ggpar(currPlot, legend.title = titlesHash[[splitCol]])
    }
    if (!is.null(rowLabel)) {
        currPlot <- currPlot + ggtitle(rowLabel)
    }
    return(currPlot)
}

# Input:  df (dataframe), predictor (string), response (string), methodString (string), titlesHash (hashmap of titles)
# Output: Geographical catchment-based correlation plot
correl_single_predictor_single_response_geo_by_catchment <- function(df, predictor, response, methodString, titlesHash) {
    northAmerica <- ne_states(c("united states of america", "canada"), returnclass = "sf")
    plotTitle <- titlesHash[[response]]
    dataColumnsNeeded <- c("catchment", predictor, response)
    dfCatchmentGeoms <- distinct(df[c("catchment", "latitude", "longitude")])
    dfForGeoPlot <- df[dataColumnsNeeded]
    dfForGeoPlot <- group_by(dfForGeoPlot, catchment)
    dfForGeoPlot <- summarize(dfForGeoPlot, cor = cor(.data[[predictor]], .data[[response]],
        method = methodString
    ))

    # Validate joins worked properly by not dropping catchments
    assert(nrow(dfForGeoPlot) == nrow(dfCatchmentGeoms))
    dfForGeoPlot <- inner_join(dfCatchmentGeoms, dfForGeoPlot, by = "catchment")
    assert(nrow(dfCatchmentGeoms) == nrow(dfCatchmentGeoms))

    methodSymbol <- "()"
    if (methodString == "pearson") {
        methodSymbol <- "(R)"
    }

    geoCorrelPlot <- ggplot(data = northAmerica) +
        geom_sf() +
        coord_sf(xlim = c(-145, -50), ylim = c(25, 72), expand = FALSE) +
        geom_point(data = dfForGeoPlot, aes(x = longitude, y = latitude, fill = cor), shape = 21) +
        scale_size_manual(values = c(0.75)) +
        scale_fill_gradient2(high = "blue", mid = "white", low = "red", midpoint = 0, limits = c(-1, 1)) +
        labs(fill = paste("Correlation with", titlesHash[[predictor]], methodSymbol), title = plotTitle) +
        theme(
            plot.title = element_text(hjust = 0.5),
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank()
        )
    return(geoCorrelPlot)
}

# Input: df (dataframe), predictorsVector (vector of strings), responsesVector (vector of strings),
# methodString (one of "pearson", "kendall", "spearman"), titlesHash (hashmap of titles),
# splitColumn (which column to split), splitColumnLabels (label for that column), rowLabels (label per row)
# Output: List of list of plots
correl_predictors_responses_plot_list <- function(df, predictorsVector, responsesVector, methodString, titlesHash, splitCol = NULL, splitColLabels = NULL, rowLabels = NULL) {
    allResponseCorrelPlotsVector <- list()
    plotNumber <- 1
    for (currResponse in responsesVector) {
        correlPlotVectorSingleResponse <- list()
        for (currPredictor in predictorsVector) {
            currRowLabel <- getRowLabel(rowLabels, plotNumber)
            vectorElemName <- paste(eval(currPredictor), eval(currResponse), sep = "-")
            correlPlotVectorSingleResponse[[vectorElemName]] <- correl_single_predictor_single_response(df, currPredictor, currResponse, methodString, titlesHash, splitCol, splitColLabels, currRowLabel)
            plotNumber <- plotNumber + 1
        }
        allResponseCorrelPlotsVector[[eval(currResponse)]] <- correlPlotVectorSingleResponse
    }
    return(allResponseCorrelPlotsVector)
}

# Input:  df (dataframe), predictorsVector (vector of strings), responsesVector (vector of strings),
# titlesHash (hashmap of titles), methodString (one of "pearson", "kendall", "spearman")
# Output: List of list of plots (geographical/catchment-based correlations)
correl_predictors_responses_geo_plot_list <- function(df, predictorsVector, responsesVector, methodString, titlesHash) {
    allResponseCorrelPlotsVector <- list()
    for (currResponse in responsesVector) {
        correlPlotVectorSingleResponse <- list()
        for (currPredictor in predictorsVector) {
            vectorElemName <- paste(eval(currPredictor), eval(currResponse), sep = "-")
            correlPlotVectorSingleResponse[[vectorElemName]] <- correl_single_predictor_single_response_geo_by_catchment(df, currPredictor, currResponse, methodString, titlesHash)
        }
        allResponseCorrelPlotsVector[[eval(currResponse)]] <- correlPlotVectorSingleResponse
    }
    return(allResponseCorrelPlotsVector)
}

# Input: plotList (list of plots), ncol (integer), nrow (integer),
# legendPosition (See ggplot2 documentation for legend positions)
# Output: An arranged plot
arrangePlots <- function(plotList, ncol, nrow, legendPosition = "bottom") {
    if (is.null(ncol) && is.null(nrow)) {
        arrangedPlot <- ggpubr::ggarrange(
            plotlist = plotList,
            common.legend = TRUE,
            legend = legendPosition
        )
    } else if (!is.null(ncol) && is.null(nrow)) {
        arrangedPlot <- ggpubr::ggarrange(
            ncol = ncol,
            plotlist = plotList,
            common.legend = TRUE,
            legend = legendPosition
        )
    } else if (is.null(ncol) && !is.null(nrow)) {
        arrangedPlot <- ggpubr::ggarrange(
            nrow = nrow,
            plotlist = plotList,
            common.legend = TRUE,
            legend = legendPosition
        )
    } else {
        arrangedPlot <- ggpubr::ggarrange(
            ncol = ncol,
            nrow = nrow,
            plotlist = plotList,
            common.legend = TRUE,
            legend = legendPosition
        )
    }
    return(arrangedPlot)
}

# Input: ggplot, plot title (string)
# Output: ggplot decorated with plot title
addPlotTitle <- function(arrangedPlot, plotTitle) {
    plotTitle <- ggdraw() +
        draw_label(plotTitle)
    arrangedPlot <- plot_grid(plotTitle, arrangedPlot, nrow = 2, rel_heights = c(0.05, 1))
    return(arrangedPlot)
}

# Input: vector of plot titles (or null), positional argument for plot
# Output: plot title (null or string)
getRowLabel <- function(rowLabels, position) {
    if (is.null(rowLabels)) {
        return(NULL)
    }
    return(rowLabels[position])
}

# Input: inputPlotList (List of List of plots), predictorsVector, responsesVector, plotTitle (string), filename (string), width and height of plot (in)
# ncol (integer for number of rows; default set to length of predictors), nrow (integer for number of rows; default set to length of responses),
# legendPosition (See ggplot2 documentation for legend positions)
# Output: Arranged plot (for further arrangement)
plot_list_of_plots <- function(inputPlotList, predictorsVector, responsesVector, plotTitle, filename, width, height, ncol = length(predictorsVector),
                               nrow = length(responsesVector), legendPosition = "bottom") {
    plotList <- list()
    for (response in responsesVector) {
        plotList <- append(plotList, inputPlotList[[response]])
    }
    arrangedPlot <- arrangePlots(plotList, ncol, nrow, legendPosition)
    if (!is.null(plotTitle)) {
        arrangedPlot <- addPlotTitle(arrangedPlot, plotTitle)
    }
    arrangedPlot <- arrangedPlot + theme_light() + theme(
        text = element_text(color = "black"),
        axis.title = element_text(family = "sans", ),
        legend.text = element_text(family = "sans", ),
        legend.title = element_text(family = "sans", ),
        strip.text = element_text(family = "sans", )
    )
    ggsave(filename, arrangedPlot, width = width, height = height, units = "in", dpi = 320)
    return(arrangedPlot)
}

# Get annual datafiles
catchmentFilenames <- list.files(path = annualDataDir, pattern = annualDataPattern)

# Read overall data
overallData <- read_csv(paste(dataDir, overallDataFilename, sep = ""))

# Populate dataframes
for (currCatchmentFilename in catchmentFilenames) {
    # Read data
    currCatchmentFilepath <- paste(annualDataDir, currCatchmentFilename, sep = "")
    currDf <- read_csv(currCatchmentFilepath)
    currDf <- drop_na(currDf)

    # Get catchment number
    currCatchmentString <- str_extract(currCatchmentFilename, catchmentNumberPattern)
    currCatchmentString <- str_replace(currCatchmentString, "_", "")
    currCatchment <- strtoi(currCatchmentString, 10)

    # Warn against low-year data
    if (nrow(currDf) < 15) {
        print(paste(currCatchmentFilename, "does not have >= 15 years of data"))
        catchmentsLt15Years <- c(catchmentsLt15Years, currCatchmentFilename)
        smallestNumberOfYearsInACatchment <- min(smallestNumberOfYearsInACatchment, nrow(currDf))
    }
    if (!(currCatchment %in% overallData$catchment)) {
        stop("Catchment does not belong in overall data")
    }
    assert(unique(currDf$catchment) == currCatchment)

    # Bind unstandardized data
    allYearsNoStandardizationDf <- bind_rows(allYearsNoStandardizationDf, currDf)

    # Center data
    standardizedCurrDf <- data.frame(currDf)
    standardizedCurrDf[predictors] <- scale(standardizedCurrDf[predictors])
    standardizedCurrDf[responses] <- scale(standardizedCurrDf[responses])

    # Seperate data by overall snow richness
    currCatchmentAllYearsMeanSP <- overallData[overallData$catchment == currCatchment, ]$mean_sp
    if (currCatchmentAllYearsMeanSP >= spRichnessThresh) {
        allYearsSnowRichNoStandardizationDf <- bind_rows(allYearsSnowRichNoStandardizationDf, currDf)
        allYearsSnowRichStandardizedDf <- bind_rows(allYearsSnowRichStandardizedDf, standardizedCurrDf)
        catchmentsSnowRich <- c(catchmentsSnowRich, currCatchment)
    } else {
        allYearsSnowPoorNoStandardizationDf <- bind_rows(allYearsSnowPoorNoStandardizationDf, currDf)
        allYearsSnowPoorStandardizedDf <- bind_rows(allYearsSnowPoorStandardizedDf, standardizedCurrDf)
        catchmentsSnowPoor <- c(catchmentsSnowPoor, currCatchment)
    }
}

# All Data Correlations
allSitesNoStandardizationPlots <- correl_predictors_responses_plot_list(allYearsNoStandardizationDf, predictorsWithLnAI, responses, "pearson", plotTitlesHash, rowLabels = overallCorrelationRowLabels)
allSitesNoStandardizationPlotFilename <- paste(correlationalPlotDir, "figS1_unstandardized_correlations_pearson.jpeg", sep = "")
plot_list_of_plots(allSitesNoStandardizationPlots, predictorsWithLnAI, responses, NULL, allSitesNoStandardizationPlotFilename, 12, 15)

# Geographical Correlations
# Snow Rich
snowRichStandardizedGeoPlots <- correl_predictors_responses_geo_plot_list(allYearsSnowRichNoStandardizationDf, spPredictorVector, responses, "pearson", geoPlotTitlesHash)
snowRichStandardizedGeoPlotsFilename <- paste(correlationalPlotDir, "fig04_snow_rich_correlations_geo_pearson.jpeg", sep = "")
snowRichStandardizedGeoPlotsPanel <- plot_list_of_plots(snowRichStandardizedGeoPlots, spPredictorVector, responses, NULL, snowRichStandardizedGeoPlotsFilename, 12, 15, NULL, NULL, "top")

# # Snow Poor (not reported)
# snowPoorStandardizedGeoPlots <- correl_predictors_responses_geo_plot_list(allYearsSnowPoorNoStandardizationDf, spPredictorVector, responses, "pearson", geoPlotTitlesHash)
# snowPoorStandardizedGeoPlotsFilename <- paste(correlationalPlotDir, "snow_poor_correlations_geo_pearson.jpeg", sep = "")
# snowPoorStandardizedGeoPlotsPanel <- plot_list_of_plots(snowPoorStandardizedGeoPlots, spPredictorVector, responses, NULL, snowPoorStandardizedGeoPlotsFilename, 12, 15, NULL, NULL, "top")

# Snow Rich
# AI
snowRichStandardizedPlotsAi <- correl_predictors_responses_plot_list(allYearsSnowRichStandardizedDf, spPredictorVector, responses, "pearson", stdPlotTitlesHash, "ai", aiLabels, temporalCorrelationRowLabelsAi)
snowRichStandardizedPlotsAiTitle <- "Drier (std. AI ≥ 0)  vs. wetter (std. AI < 0) years"
snowRichStandardizedPlotsAiFilename <- paste(correlationalPlotDir, "fig05_snow_rich_correlations_pearson_ai.jpeg", sep = "")
snowRichStandardizedPlotsAiPanel <- plot_list_of_plots(snowRichStandardizedPlotsAi, spPredictorVector, responses, snowRichStandardizedPlotsAiTitle, snowRichStandardizedPlotsAiFilename, 12, 15)
# SI
snowRichStandardizedPlotsSi <- correl_predictors_responses_plot_list(allYearsSnowRichStandardizedDf, spPredictorVector, responses, "pearson", stdPlotTitlesHash, "si", siLabels, temporalCorrelationRowLabelsSi)
snowRichStandardizedPlotsSiTitle <- "In-phase (std. SI ≥ 0) vs. out-of-phase (std. SI < 0) years"
snowRichStandardizedPlotsSiFilename <- paste(correlationalPlotDir, "fig05_snow_rich_correlations_pearson_si.jpeg", sep = "")
snowRichStandardizedPlotsSiPanel <- plot_list_of_plots(snowRichStandardizedPlotsSi, spPredictorVector, responses, snowRichStandardizedPlotsSiTitle, snowRichStandardizedPlotsSiFilename, 12, 15)

snowRichCorrelationsAiSiFilename <- paste(correlationalPlotDir, "fig05_snow_rich_correlations_pearson.jpeg", sep = "")
snowRichCorrelationsAiSi <- ggpubr::ggarrange(snowRichStandardizedPlotsAiPanel, snowRichStandardizedPlotsSiPanel)
ggsave(snowRichCorrelationsAiSiFilename,
    snowRichCorrelationsAiSi,
    height = 20,
    width = 15,
    units = "in",
    dpi = 320
)

# # Snow Poor (not reported)
# # AI
# snowPoorStandardizedPlotsAi <- correl_predictors_responses_plot_list(allYearsSnowPoorStandardizedDf, spPredictorVector, responses, "pearson", stdPlotTitlesHash, "ai", aiLabels, temporalCorrelationRowLabelsAi)
# snowPoorStandardizedPlotsAiTitle <- "Drier (std. AI ≥ 0)  vs. wetter (std. AI < 0) years"
# snowPoorStandardizedPlotsAiFilename <- paste(correlationalPlotDir, "snow_poor_correlations_pearson_ai.jpeg", sep = "")
# snowPoorStandardizedPlotsAiPanel <- plot_list_of_plots(snowPoorStandardizedPlotsAi, spPredictorVector, responses, snowPoorStandardizedPlotsAiTitle, snowPoorStandardizedPlotsAiFilename, 12, 15)
# # SI
# snowPoorStandardizedPlotsSi <- correl_predictors_responses_plot_list(allYearsSnowPoorStandardizedDf, spPredictorVector, responses, "pearson", stdPlotTitlesHash, "si", siLabels, temporalCorrelationRowLabelsSi)
# snowPoorStandardizedPlotsSiTitle <- "In-phase (std. SI ≥ 0) vs. out-of-phase (std. SI < 0) years"
# snowPoorStandardizedPlotsSiFilename <- paste(correlationalPlotDir, "snow_poor_correlations_pearson_si.jpeg", sep = "")
# snowPoorStandardizedPlotsSiPanel <- plot_list_of_plots(snowPoorStandardizedPlotsSi, spPredictorVector, responses, snowPoorStandardizedPlotsSiTitle, snowPoorStandardizedPlotsSiFilename, 12, 15)

# snowPoorStandardizedPlotsAiSiFilename <- paste(correlationalPlotDir, "snow_poor_correlations_pearson.jpeg", sep = "")
# snowPoorStandardizedPlotsAiSi <- ggpubr::ggarrange(snowPoorStandardizedPlotsAiPanel, snowPoorStandardizedPlotsSiPanel)
# ggsave(snowPoorStandardizedPlotsAiSiFilename,
#     snowPoorStandardizedPlotsAiSi,
#     height = 20,
#     width = 15,
#     units = "in",
#     dpi = 320
# )