## -----------------------------------------------------------------------------
## Required results : tabA2
## Produces         : Table A2
## -----------------------------------------------------------------------------


rm(list=ls())

results_file <- "../results/tabA2/simes_perm_results.csv"

dir.create("tables", showWarnings = FALSE, recursive = TRUE)


out_df <- read.csv(results_file)
MS              <- out_df$m
SIZE_SIMES      <- out_df$size_simes
SIZE_SIMES_PERM <- out_df$size_simes_perm
CRIT_VAL        <- out_df$crit_val

## -----------------------------------------------------------------------------
## Table
## -----------------------------------------------------------------------------
tab <- round(rbind(MS, SIZE_SIMES, SIZE_SIMES_PERM, CRIT_VAL), 3)
print(tab)

## LaTeX output
fmt <- function(x) sprintf("%.3f", x)
ncol_tab <- length(MS)
col_spec <- paste0("r|", paste(rep("r", ncol_tab), collapse=""))

latex_lines <- c(
    "\\begin{tabular}{", col_spec, "}",
    "\\toprule",
    paste0("$m$ & ", paste(MS, collapse=" & "), " \\\\"),
    "\\midrule",
    paste0("$\\mathbb{P}(\\phi_S^{\\mathrm{Simes}}=1 \\mid H_S)$ & ",
           paste(fmt(SIZE_SIMES), collapse=" & "), " \\\\"),
    paste0("$\\mathbb{P}(\\phi_S^{\\mathrm{SimesPerm}}=1 \\mid H_S)$ & ",
           paste(fmt(SIZE_SIMES_PERM), collapse=" & "), " \\\\"),
    paste0("$\\alpha^{\\mathrm{Perm}}$ & ",
           paste(fmt(CRIT_VAL), collapse=" & "), " \\\\"),
    "\\hline",
    "\\end{tabular}"
)
latex_str <- paste(c(paste0(latex_lines[1], latex_lines[2], latex_lines[3]),
                     latex_lines[4:length(latex_lines)]),
                   collapse="\n")
cat(latex_str, "\n")

## Also write LaTeX to file
writeLines(latex_str, "tables/tabA2_simes_perm.tex")
