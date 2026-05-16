## -----------------------------------------------------------------------------
## Shared utilities for all plot_*.R scripts.
##
## Every plot_*.R script should start with:
##     source("plotting_utils.R")
##
## This file provides:
##   - Library imports
##   - load_data(fig)  : reads CSV files from results/fig/
##   - init_settings(...): defines global color/shape/alpha aesthetic mappings
##                         used across all plotting functions
## -----------------------------------------------------------------------------

options(width = 300)

library(tidyverse)
library(latex2exp)
library(RColorBrewer)
library(kableExtra)

## Read all CSV result files for a given fig number, concatenated.
load_data <- function(fig) {
    idir <- sprintf("../results/%s", fig)
    ifile.list <- list.files(idir)
    results <- do.call("rbind", lapply(ifile.list, function(ifile) {
        df <- read_delim(sprintf("%s/%s", idir, ifile), delim=",", col_types=cols(), guess_max=2)
    }))
    return(results)
}

## Define the global aesthetic mappings used by all plotting functions.
## Writes to the global environment (via <<-) so that downstream functions
## can use method.values, method.labels, color.scale, etc. directly.
##
## Args:
##   idx.exclude     : optional integer indices to drop from method.values
##                     (position-based; applied after present_methods filter)
##   names_ACODE     : if TRUE, prefix method labels with "ACODE ("
##   present_methods : optional character vector; restricts method.values to
##                     methods actually present in the data
init_settings <- function(idx.exclude=NULL, names_ACODE=FALSE, present_methods=NULL) {
  cbPalette <<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#6e57d2", "red", "#DEB887")

  method.values <<- c(
    "lb_simes", "lb_storey_simes", "lb_fisher",
    "lb_wmw", "lb_higher_k2",
    "lb_lmp", "lb_auto", "disc_bh"
  )

  if (names_ACODE) {
    method.labels <<- c(
      "ACODE (Simes)", "ACODE (Storey-Simes)", "ACODE (Fisher)",
      "ACODE (WMW)", "ACODE (WMW, k=3)",
      "ACODE (LMP,G-hat)", "ACODE (Adaptive)",
      "BH (identification)"
    )
  } else {
    method.labels <<- c(
      "Simes", "Storey-Simes", "Fisher",
      "WMW", "LMPR (Lehmann, k=3)",
      "LMPR (G-hat)", "Adaptive",
      "BH (identification)"
    )
  }

  classifier.values <<- c("occ-auto", "bc-auto", "auto")
  classifier.labels <<- c("One-Class", "Binary", "Automatic")

  data.values <<- c("pendigits", "creditcard", "cover", "shuttle", "mammography", "aloi")
  data.labels <<- c("Pendigits", "Creditcard", "Covertype", "Shuttle", "Mammography", "ALOI")

  ## Named aesthetic mappings, keyed by method code.
  color.by.method <- c(
    lb_simes        = cbPalette[1],
    lb_storey_simes = cbPalette[1],
    lb_fisher       = cbPalette[4],
    lb_wmw          = cbPalette[3],
    lb_higher_k2    = cbPalette[6],
    lb_lmp          = cbPalette[7],
    lb_auto         = cbPalette[8],
    disc_bh         = cbPalette[9]
  )

  shape.by.method <- c(
    lb_simes        = 2,
    lb_storey_simes = 6,
    lb_fisher       = 3,
    lb_wmw          = 1,
    lb_higher_k2    = 0,
    lb_lmp          = 9,
    lb_auto         = 8,
    disc_bh         = 5
  )

  alpha.by.method <- c(
    lb_simes        = 0.75,
    lb_storey_simes = 0.75,
    lb_fisher       = 0.75,
    lb_wmw          = 0.75,
    lb_higher_k2    = 0.75,
    lb_lmp          = 0.75,
    lb_auto         = 1,
    disc_bh         = 0.75
  )

  ## Optionally restrict to methods actually present in the data
  if (!is.null(present_methods)) {
    keep <- method.values %in% present_methods
    method.values <<- method.values[keep]
    method.labels <<- method.labels[keep]
  }

  ## Apply position-based exclusion (after present_methods filter)
  if (!is.null(idx.exclude) && length(idx.exclude) > 0) {
    method.values <<- method.values[-idx.exclude]
    method.labels <<- method.labels[-idx.exclude]
  }

  ## Build position-indexed scales (legacy use)
  color.scale <<- unname(color.by.method[method.values])
  shape.scale <<- unname(shape.by.method[method.values])
  alpha.scale <<- unname(alpha.by.method[method.values])

  ## Defensive check: catch typos / unmapped methods early
  if (any(is.na(color.scale))) stop("Missing color mapping for: ", paste(method.values[is.na(color.scale)], collapse=", "))
  if (any(is.na(shape.scale))) stop("Missing shape mapping for: ", paste(method.values[is.na(shape.scale)], collapse=", "))
  if (any(is.na(alpha.scale))) stop("Missing alpha mapping for: ", paste(method.values[is.na(alpha.scale)], collapse=", "))

  ## Label-indexed scales: robust to dropped/unused factor levels in plots
  ## that match by Method label rather than by position.
  color.by.label <<- setNames(unname(color.by.method[method.values]), method.labels)
  shape.by.label <<- setNames(unname(shape.by.method[method.values]), method.labels)
  alpha.by.label <<- setNames(unname(alpha.by.method[method.values]), method.labels)
}
