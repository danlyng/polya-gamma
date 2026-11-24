# ---------------------------------------------------------------------
# Compare Rust PG(1,z) sampler with BayesLogit reference
# ---------------------------------------------------------------------
library(readr)      # fast file I/O
library(BayesLogit) # reference implementation
library(ggplot2)    # plots
library(dplyr)      # simple data wrangling
library(e1071)      # skewness / kurtosis  (base R works as well)

set.seed(1)         # reproducibility
N      <- 1e5       # sample size produced by both samplers
z_grid <- c(0.5, 1, 2, 3.2, 5)   # <- change or reduce as you wish
outdir <- "~/projects.nosync/polya-gamma"  # where Rust files live

# ------------- helper functions --------------------------------------

mean_theory <- function(z) tanh(z/2) / (2*z)

one_z_check <- function(z) {
  file <- file.path(outdir, sprintf("out_%g.txt", z))
  rust  <- read_csv(file, col_names = FALSE, show_col_types = FALSE)[[1]]
  stopifnot(length(rust) == N)

  ref   <- BayesLogit::rpg(num = N, h = 1, z = z)

  # basic numbers -----------------------------------------------------
  df_sum <- tibble(
    which = c("Rust", "Reference"),
    mean  = c(mean(rust), mean(ref)),
    var   = c(var(rust),  var(ref)),
    sd    = c(sd(rust),   sd(ref)),
    skew  = c(skewness(rust), skewness(ref)),
    kurt  = c(kurtosis(rust), kurtosis(ref))
  )
  cat("\n------------------------------------------------------------\n")
  cat(sprintf("PG(1, %.3g)  —  sample size = %d\n", z, N))
  print(df_sum, digits = 4)

  # theoretical mean
  cat(sprintf("Theoretical mean  = %.6f\n", mean_theory(z)))

  # formal tests ------------------------------------------------------
  cat("\nTwo–sample t-test:\n")
  print(t.test(rust, ref))

  cat("\nKolmogorov–Smirnov test:\n")
  print(ks.test(rust, ref))

  # plots -------------------------------------------------------------
  # QQ-plot
  qqplot(rust, ref, pch = 20,
         main = sprintf("QQ plot  PG(1, %.2f): Rust vs BayesLogit", z),
         xlab = "Rust quantiles", ylab = "Reference quantiles")
  abline(0, 1, col = "red")

  # Overlaid kernel densities
  d <- bind_rows(
    data.frame(value = rust, which = "Rust"),
    data.frame(value = ref,  which = "Reference")
  )
  plt <- ggplot(d, aes(x = value, color = which)) +
      geom_density(adjust = 1.2, size = 1) +
      theme_minimal() +
      labs(title = sprintf("PG(1, %.2f) – density comparison", z))
}

# ------------- run the checks ----------------------------------------
for (z in z_grid) one_z_check(z)