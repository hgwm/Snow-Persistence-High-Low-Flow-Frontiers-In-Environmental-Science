# plotter_sites_signatures_climatic_indices.R plots overall signature/climatic data for all available catchments in the Le et al. (2023)
# Input:    - Overall signature/climatic data (assumed to come pre-loaded with latitudes and longitudes)
# Output:   - Plots for all signatures found within the trend and signature/climatic data
# Author(s): Edward Le, Joe Janssen
# Note: "cat" refers to the "category" of data (e.g., signatures or climatic indices)

library(ggplot2)
library(sf)
library(stringr)
library(rnaturalearth)
library(rnaturalearthdata)
library(readr)
library(sp)
library(maps)
library(ggspatial)
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(gsubfn)
library(hash)
library(ggeasy)
library(patchwork)
library(ggsn)
library(BAMMtools)
library(ggnewscale)
library(scales)
library(viridis)
library(ggnewscale)
library(cowplot)
library(egg)
library(plyr)

filenameTag <- ""

# Set working directory
setwd("./")

# Create plot folders
plotDir <- "../Plots/Geo/"
dir.create("../Plots", showWarnings = FALSE)
dir.create("../Plots/Geo", showWarnings = FALSE)

# Read attribute data
overallAttributes <- read_csv("../Data/overall_water_climate_data.csv")
# Drop na
overallAttributes <- drop_na(overallAttributes)

northAmerica <- ne_states(c("united states of america", "canada"), returnclass = "sf")

panelOrderHash <- hash()
panelOrderHash[["mean_sp"]] <- "a"
panelOrderHash[["mean_ai"]] <- "b"
panelOrderHash[["mean_si"]] <- "c"
panelOrderHash[["bfi"]] <- "d"
panelOrderHash[["q5Frac"]] <- "e"
panelOrderHash[["low_fdc"]] <- "f"
panelOrderHash[["mean_low_flow_dur"]] <- "g"
panelOrderHash[["q95Frac"]] <- "h"
panelOrderHash[["high_fdc"]] <- "i"
panelOrderHash[["mean_high_flow_dur"]] <- "j"

# Titles List:
# Unfortunately, expression() is buggy with string/expression concatenation for Q5/Q95 subscripting; duplicate
# data on ordering here
plotTitlesHash <- hash()
plotTitlesHash[["mean_sp"]] <- "(A) Snow persistence"
plotTitlesHash[["mean_ai"]] <- "(B) Aridity index"
plotTitlesHash[["mean_si"]] <- "(C) Seasonality index"
plotTitlesHash[["bfi"]] <- "(D) Baseflow index"
plotTitlesHash[["q5Frac"]] <- expression("(E) Normalized" ~ Q[5], "")
plotTitlesHash[["low_fdc"]] <- "(F) Low flow duration curve slope\n(0.05-0.30)"
plotTitlesHash[["mean_low_flow_dur"]] <- "(G) Mean low-flow duration (days)"
plotTitlesHash[["q95Frac"]] <- expression("(H) Normalized" ~ Q[95], "")
plotTitlesHash[["high_fdc"]] <- "(I) High flow duration curve slope\n(0.70-0.95)"
plotTitlesHash[["mean_high_flow_dur"]] <- "(J) Mean high-flow duration (days)"

# Figure prefix ordering
attributeGeoPlotFigures <- "fig01"

untransformedVariables <- c("mean_sp", "mean_si", "bfi", "q5Frac", "mean_low_flow_dur", "mean_high_flow_dur")

# Helper functions

# Create attribute plot filenames
# Dependencies: panelOrderHash encoded with panel orders
# Input: destPlotDir (string), sigName (string), prefixString (string), suffixString (string)
# Ouput: Formatted filename string
create_sig_filenames_string <- function(destPlotDir, sigName, prefixString, suffixString) {
  paste(destPlotDir, prefixString, panelOrderHash[[sigName]], "_", sigName, suffixString, sep = "")
}

# Find minimum, midpoint, and maximum number from a distribution
# Will omit min and 25% to middle if they are the same (when rounded to two decimal places untransformed or log transform)
# Input: nums (Vector of numbers)
# Output: Vector of min, middle, and max values
get_min_mids_max <- function(nums) {
  min <- min(nums)
  quarterToMiddle <- min(nums) + ((max(nums) - min(nums)) / 4)
  middle <- min(nums) + ((max(nums) - min(nums)) / 2)
  threeQuartersToMiddle <- max(nums) - ((max(nums) - min(nums)) / 4)
  max <- max(nums)
  minMidsMax <- c(
    # Min
    min,
    # 25% to middle
    quarterToMiddle,
    # Middle
    middle,
    # 75%
    threeQuartersToMiddle,
    # Max
    max
  )
  if (round(min, 2) == round(quarterToMiddle, 2) | (round(exp(min), 2) == round(exp(quarterToMiddle), 2))) {
    minMidsMax <- c(
      # Min
      min,
      # Middle
      middle,
      # 75%
      threeQuartersToMiddle,
      # Max
      max
    )
  }
  return(minMidsMax)
}

# Formats min, mid, max with two decmial places
# Input: nums (vector of numbers), decPlaces (integer)
# Output: formatted vector of strings
format_list_decimal_places <- function(nums, decPlaces = 2) {
  nums <- lapply(nums, round, decPlaces)
  nums <- lapply(nums, format, nsmall = decPlaces)
  return(nums)
}

# Input: spatialData (ne_states data), catDataFrame (datatrame), currCat (string),
# minMidsMax (vector of numbers), minMidsMaxUntransformedLabels (vector of strings),
# plotTitlesHash (hasmap of titles)
# Output: Geospatial plot of data
longTermPlotter <- function(spatialData, catDataFrame, currCat, minMidsMax, minMidsMaxUntransformedLabels, plotTitlesHash) {
  longTermPlot <- ggplot(data = spatialData) +
    geom_sf() +
    theme(
      legend.text = element_text(size = 11),
      legend.key.height = unit(0.5, "cm"),
      legend.spacing.y = unit(0.25, "cm")
    ) +
    coord_sf(xlim = c(-145, -50), ylim = c(25, 72), expand = FALSE) +
    geom_point(
      data = catDataFrame,
      na.rm = TRUE,
      shape = 21,
      aes_string(
        x = "longitude",
        y = "latitude",
        fill = currCat
      ),
      inherit.aes = FALSE
    ) +
    scale_fill_gradient2(
      low = "#d7191c",
      mid = "#ffffbf",
      high = "#2c7bb6",
      midpoint = as.numeric(minMidsMax[ceiling(length(minMidsMax) / 2)]),
      breaks = minMidsMax,
      labels = minMidsMaxUntransformedLabels,
    ) +
    labs(
      title = plotTitlesHash[[currCat]],
      fill = NULL
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, ),
      plot.margin = grid::unit(c(0.01, 0.01, 0.01, 0.01), "null"),
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )
  return(longTermPlot)
}

# Plot Characteristics
catPlotList <- list()

# First col of attribute data is catchment number, 2nd col is longitude, 3rd is latitude
categories <- colnames(overallAttributes)[4:length(colnames(overallAttributes))]
print(categories)
catLen <- length(categories)
for (i in 1:catLen) {
  currCat <- categories[i]
  if (currCat %in% keys(plotTitlesHash)) {
    catDataFrame <- data.frame(overallAttributes[c("catchment", "latitude", "longitude", currCat)])

    # Min/Midpoints/Max breaks; midpoint defined as middle between min and max
    minMidsMax <- get_min_mids_max(catDataFrame[, currCat])
    # Round and then pad 0s, turn into strings
    minMidsMaxUntransformedLabels <- format_list_decimal_places(minMidsMax)

    # Special treatment for low-flow durations and high-flow durations due to skew
    if (grepl("dur", currCat, fixed = TRUE)) {
      # Round to nearest 10
      midpoint <- round_any(minMidsMax[ceiling(length(minMidsMax) / 2)], 10)

      catDataFrame[catDataFrame[, currCat] > midpoint, currCat] <- midpoint
      minMidsMax <- get_min_mids_max(catDataFrame[, currCat])
      minMidsMaxUntransformedLabels <- format_list_decimal_places(minMidsMax)
      minMidsMaxUntransformedLabels[length(minMidsMaxUntransformedLabels)] <- paste(minMidsMaxUntransformedLabels[length(minMidsMaxUntransformedLabels)], "+", sep = "")
    }

    # Apply log transformation for visual clarity
    if (!(currCat %in% untransformedVariables)) {
      catDataFrame[, currCat] <- log(catDataFrame[, currCat])
      # Min/Midpoints/Max breaks; midpoint defined as middle between min and max
      minMidsMax <- get_min_mids_max(catDataFrame[, currCat])
      # Untransform min, mid, max for labels
      minMidsMaxUntransformedLabels <- exp(minMidsMax)
      # Round and then pad 0s, turn into strings
      minMidsMaxUntransformedLabels <- format_list_decimal_places(minMidsMaxUntransformedLabels)
    }
    # Long-Term Value Plot
    longTermPlot <- longTermPlotter(northAmerica, catDataFrame, currCat, minMidsMax, minMidsMaxUntransformedLabels, plotTitlesHash)
    longTermPlotFileName <- create_sig_filenames_string(
      destPlotDir = plotDir,
      sigName = currCat,
      prefixString = attributeGeoPlotFigures,
      suffixString = paste(filenameTag, ".jpeg", sep = "")
    )
    ggsave(
      filename = longTermPlotFileName,
      plot = longTermPlot,
      dpi = 320
    )

    # Convert letter order to number
    catPlotOrder <- which(letters == panelOrderHash[[currCat]])
    catPlotList[[catPlotOrder]] <- longTermPlot
  }
}

# Plot all panels as one
allCatPlotsName <- paste(plotDir, attributeGeoPlotFigures, filenameTag, ".jpeg", sep = "")
allCatPlots <- grid.arrange(grobs = catPlotList, ncol = 2)
ggsave(allCatPlotsName,
  allCatPlots,
  height = 11,
  width = 8.5,
  dpi = 320
)