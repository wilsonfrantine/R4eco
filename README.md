# Ecological Analysis Made Easy

R4eco is a small R package designed to simplify ecological analyses by providing a range of functions for diverse tasks. While the package is continuously expanding, it currently includes functions for Renkonen's similarity index (S) calculation, bootstrap implementation, and convenient plotting of UPGMA dendrograms.

## Installation

To install R4eco, use the following code:

```{r
# Install R4eco
remotes::install_github("wilsonfrantine/R4eco")

# Load the package
library(R4eco)

```

## # Example Data and Model

```{r}
d <- data.frame(
  Type = rep(c("Forest", "Regeneration", "Restoration"), each = 12),
  Landscape = rep(paste0("L", 1:12), times = 3),
  Mean_NDVI_SD_500 = rnorm(36, mean = 0.2, sd = 0.02),
  hill_q0 = sample(5:14, 36, replace = TRUE),
  Abundance = sample(5:50, 36, replace = TRUE)
)

modelX <- lme4::lmer(formula = hill_q0 ~ Mean_NDVI_SD_500 * Type + (1|Landscape), data = d)

```
## lmerPredictionPlot Examples

```{r}
# Plot linear prediction
lmerPredictionPlot(model = modelX, type = "linear.prediction")

# Customizing plot with color and fill scales
lmerPredictionPlot(model = modelX) +
  scico::scale_color_scico_d(palette = "batlow", begin = 0.1, end = 0.7) +
  scico::scale_fill_scico_d(palette = "batlow", begin = 0.1, end = 0.7)

# Adding facet_wrap to the plot
lmerPredictionPlot(model = modelX) +
  scico::scale_color_scico_d(palette = "batlow", begin = 0.1, end = 0.7) +
  scico::scale_fill_scico_d(palette = "batlow", begin = 0.1, end = 0.7) +
  ggplot2::facet_wrap(~Type)
```

## Renkonen Similarity Calculation

To begin, you can generate sample data using the simulate_data() function: `simulate_data()`.

```{r}
# Generate sample data
simulate_data()

```

### Calculate Renkonen similarity with or without bootstrapping:

```{r}
# Calculate Renkonen similarity
df <- simulate_data(5, 9)
S <- renkonen(df)
head(S)

# Calculate Renkonen similarity with bootstrap
df <- simulate_data()
S <- renkonen(df, boot = 1000)
head(S)
```

## UPGMA Dendrogram Plotting

Visualize Renkonen similarity as a UPGMA dendrogram:

```{r}
# Plot UPGMA dendrogram
df <- simulate_data()
S <- renkonen(df, boot = 1000)
trees <- upgma(S)
plot_upgma(trees)
```

## Consensus Dendrogram
Compute a consensus dendrogram for saving as Newick format or using in other software:

```{r}
# Compute consensus dendrogram
df <- simulate_data()
S <- renkonen(df, boot = 1000)
trees <- upgma(S)
tree_consensus <- tree_consensus(trees, type = "mcc")
```

## Contact
For assistance or inquiries, please reach out to:
Wilson Frantine-Silva
wilsonfrantine@gmail.com

Enjoy using R4eco to enhance your ecological analyses!
