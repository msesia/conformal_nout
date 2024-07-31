## Load necessary libraries
library(tidyverse)    ## For data manipulation and visualization
library(progress)     ## For displaying progress bars
library(nout)         ## Custom library (assuming this contains necessary functions)

## Source utility functions for data generation and experiments
source("../R/utils_data.R")
source("../R/utils_experiments.R")
source("../R/utils_g.R")

###########################
## Experiment parameters ##
###########################

## Flag to determine if input should be parsed from command line
parse_input <- TRUE

if(parse_input) {
    ## Reading command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    # Checking if the correct number of arguments is provided
    if (length(args) < 6) {
        stop("Insufficient arguments provided. Expected 6 arguments.")
    }
    ## Assigning command line arguments to variables
    setup <- as.integer(args[1])
    n_cal <- as.integer(args[2])
    n_test <- as.integer(args[3])
    seed <- as.integer(args[4])
    alternative <- args[5]
    prop_out <- as.numeric(args[6])
} else {
    ## Use default values
    setup <- 2
    n_cal <- 500
    n_test <- 200
    seed <- 1
    alternative <- "lehmann_k2"
    prop_out <- 0.25
}

## Print the values to verify they are correctly assigned
cat("n_cal:", n_cal, "\n")
cat("n_test:", n_test, "\n")
cat("seed:", seed, "\n")
cat("alternative:", alternative, "\n")
cat("prop_out:", prop_out, "\n")

## Generate a unique and interpretable file name based on the input parameters
output_file <- paste0("results/", "setup", setup, "/",
  "n_cal_", n_cal, "_",
  "n_test_", n_test, "_",
  "seed_", seed, "_",
  "alt_", alternative, "_",
  "prop_out_", prop_out, ".txt"
)

## Print the output file name to verify
cat("Output file name:", output_file, "\n")

## Number of repetitions for each experimental setting
n_exp <- 10
    
## Make tibble with experiment meta-data
header <- tibble(n_cal=n_cal, n_test=n_test, alternative=alternative, prop_out=prop_out)

##########################
## Experiment functions ##
##########################

## Function to run a single experiment
## Args:
##   seed: Random seed for reproducibility
## Returns:
##   A tibble containing the results of the experiment
run_experiment <- function(i) {
    random_state = seed*1000 + i
    set.seed(random_state)  ## Set seed for reproducibility

    ## Generate calibration and test data with specified parameters
    data <- generate_cal_test_scores(n_cal = n_cal, n_test = n_test, prop_out = prop_out, alternative = alternative)

    ## Calculate true number of outliers
    n.out <- sum(data$outlier.test)

    ## Apply global testing methods to the generated data
    res <- run_outlier_enumeration(data, alternative = alternative) |> select(Method, Lower)

    ## Combine the results with experiment metadata
    results <- tibble(Seed = random_state) |> cbind(header) |> cbind(res)

    return(results)
}

## Function to run multiple experiments and gather results
## Args:
##   n_exp: Number of repetitions for each experimental setting
## Returns:
##   A tibble containing the combined results of all experiments
run_multiple_experiments <- function(n_exp) {
    results_df <- data.frame()  # Initialize an empty data frame to store cumulative results

    # Print a progress bar header
    cat("Running experiments\n")
    pb <- txtProgressBar(min = 0, max = n_exp, style = 3)  # Initialize progress bar

    # Loop over each repetition
    for (i in 1:n_exp) {
        result <- run_experiment(i)  # Run experiment and get the result

        # Convert the result to a data frame (if it's not already)
        result_df <- as.data.frame(result)

        # Add the result to the cumulative data frame
        results_df <- rbind(results_df, result_df)

        # Write the cumulative results to the CSV file
        write.csv(results_df, output_file, row.names = FALSE)

        setTxtProgressBar(pb, i)  # Update progress bar
    }

    close(pb)  # Close the progress bar

    return(results_df)  # Return the cumulative results data frame
}


#####################
## Run experiments ##
#####################

## Run the experiments with specified parameters
results <- run_multiple_experiments(n_exp)


if(FALSE) {
##################
## Plot results ##
##################

    ## Define the significance level
    alpha <- 0.1

    ## Calculate power for different methods and proportions of outliers
    power_results <- results %>%
        group_by(Method, Prop_Out) %>%
        summarize(
            Power = mean(p.value < alpha),  ## Calculate power as the proportion of p-values below 0.05
            SE = sqrt((Power * (1 - Power)) / n())  ## Calculate standard error of the power estimate
        )

    ## Plot the power for different methods and proportions of outliers
    power_results %>%
        ggplot(aes(x = Prop_Out, y = Power, color = Method, shape = Method)) +  ## Set aesthetics
        geom_line() +  ## Add lines
        geom_point() +  ## Add points
        geom_errorbar(aes(ymin = Power - 2 * SE, ymax = Power + 2 * SE), width = 0.02) +  ## Add error bars
        geom_hline(yintercept = alpha, linetype = 2) +  ## Add horizontal line at significance level
        theme_minimal(base_size = 15) +  ## Set minimal theme for the plot
        labs(
            title = "Power of Different Methods for Various Proportions of Outliers",
            subtitle = sprintf("N-cal: %d, N-test: %d, Alternative distribution: %s", n_cal, n_test, alternative),
            x = "Proportion of Outliers",
            y = "Power",
            color = "Method"
        ) +
        ylim(0, 1)  ## Set y-axis limits

}
