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

## Define the aesthetic mappings used by power-style plots (Fisher / WMW /
## Shirashi variants). Writes to the global environment via <<- so plotting
## functions can use method.values, method.labels, colors, shapes directly.
init_settings_power <- function() {
    cbPalette <<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                    "#0072B2", "#6e57d2", "red", "#DEB887")

    method.values <<- c("Fisher", "WMW",
                        "Shirashi_oracle",
                        "Shirashi_ghat_betamix",
                        "Shirashi_ghat_betamix_inc")
    method.labels <<- c("Fisher", "WMW",
                        "LMPR (oracle)",
                        "LMPR (G-hat)",
                        "LMPR (G-hat, monotone)")

    ## Named aesthetic mappings, keyed by method label.
    ## (Names avoid base::colors / grDevices::colors binding collision.)
    color.by.label <<- c("Fisher"                 = cbPalette[4],
                         "WMW"                    = cbPalette[3],
                         "LMPR (oracle)"          = cbPalette[6],
                         "LMPR (G-hat)"           = cbPalette[7],
                         "LMPR (G-hat, monotone)" = cbPalette[7])

    shape.by.label <<- c("Fisher"                 = 3,
                         "WMW"                    = 1,
                         "LMPR (oracle)"          = 8,
                         "LMPR (G-hat)"           = 15,
                         "LMPR (G-hat, monotone)" = 16)
}

## Variant for the lower-bound plots (figA12 family). The non-monotone
## "Shirashi_ghat_betamix" run is the only Shirashi variant present, and is
## labelled as the monotone version (matches the original plot_lb_1).
init_settings_lb <- function() {
    cbPalette <<- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                    "#0072B2", "#6e57d2", "red", "#DEB887")

    method.values <<- c("Fisher", "WMW",
                        "Shirashi_oracle",
                        "Shirashi_ghat_betamix")
    method.labels <<- c("Fisher", "WMW",
                        "LMPR (oracle, monotone)",
                        "LMPR (G-hat, monotone)")

    color.by.label <<- c("Fisher"                 = cbPalette[4],
                         "WMW"                    = cbPalette[3],
                         "LMPR (oracle, monotone)"= cbPalette[6],
                         "LMPR (G-hat, monotone)" = cbPalette[7])

    shape.by.label <<- c("Fisher"                 = 3,
                         "WMW"                    = 1,
                         "LMPR (oracle, monotone)"= 8,
                         "LMPR (G-hat, monotone)" = 16)
}
