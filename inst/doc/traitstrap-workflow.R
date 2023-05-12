## ---- setup, include = FALSE--------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(traitstrap)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)

# wes anderson colours
col_palettes <- list(
  GrandBudapest1 = c("#FD6467", "#5B1A18", "#D67236")
)
theme_set(theme_minimal(base_size = 12))

## ----hex, echo=FALSE, out.width='15%', out.extra='style="float:right; padding:8px"'----
knitr::include_graphics("../man/figures/Traitstrap_hex.png")

## ----true-dist, echo=FALSE, out.width='90%'-----------------------------------
knitr::include_graphics("true-dist.png")

## ----comm-boot, echo=FALSE, out.width='90%'-----------------------------------
knitr::include_graphics("comm-boot.png")

## ----data-prank, echo=FALSE, eval=TRUE----------------------------------------
community <- community |>
  mutate(Taxon = case_when(
    Taxon == "equisetum scirpoides" ~ "enquistetum scirpoides",
    Taxon == "micranthes hieracifolia" ~ "maitneranthes hieracifolia",
    Taxon == "bistorta vivipara" ~ "bistorta vigdis",
    Taxon == "stellaria humifusa" ~ "stelfordaria humifusa",
    Taxon == "oxyria digyna" ~ "oxyria tanyna",
    Taxon == "silene acaulis" ~ "silene acaudis",
    TRUE ~ Taxon
  ))

trait <- trait |>
  mutate(Taxon = case_when(
    Taxon == "equisetum scirpoides" ~ "enquistetum scirpoides",
    Taxon == "micranthes hieracifolia" ~ "maitneranthes hieracifolia",
    Taxon == "bistorta vivipara" ~ "bistorta vigdis",
    Taxon == "stellaria humifusa" ~ "stelfordaria humifusa",
    Taxon == "oxyria digyna" ~ "oxyria tanyna",
    Taxon == "silene acaulis" ~ "silene acaudis",
    TRUE ~ Taxon
  ))

# shorten trait names
trait <- trait |>
  mutate(Trait = recode(Trait, "Leaf_Thickness_Ave_mm" = "Thickness_mm"))

## ----comm-data, echo=FALSE, eval=TRUE-----------------------------------------
community

## ----trait-data, echo=FALSE, eval=TRUE----------------------------------------
trait

## ----trait-fill, echo=TRUE, eval=TRUE-----------------------------------------
trait_filling <- trait_fill(
  # input data (mandatory)
  comm = community,
  traits = trait,

  # specifies columns in your data (mandatory)
  abundance_col = "Cover",
  taxon_col = "Taxon",
  trait_col = "Trait",
  value_col = "Value",

  # specifies sampling hierarchy
  scale_hierarchy = c("Site", "PlotID"),

  # min number of samples
  min_n_in_sample = 9
)

trait_filling

## ----trait-fill2, echo=TRUE, eval=FALSE---------------------------------------
#  trait_filling2 <- trait_fill(
#    comm = community,
#    traits = trait,
#    abundance_col = "Cover",
#  
#    # defining taxonomic hierarchy
#    taxon_col = c("Taxon", "Genus"),
#    trait_col = "Trait",
#    value_col = "Value",
#    scale_hierarchy = c("Site", "PlotID"),
#    min_n_in_sample = 3,
#  
#    # specifying experimental design
#    treatment_col = "Treatment",
#    treatment_level = "Site",
#  )

## ----non-parap-boot, echo=TRUE, eval=TRUE-------------------------------------
# run nonparametric bootstrapping
np_bootstrapped_moments <- trait_np_bootstrap(
  trait_filling, 
  nrep = 50
)

np_bootstrapped_moments

## ----summarize, echo=TRUE, eval=TRUE------------------------------------------
# summarizes bootstrapping output
sum_boot_moment <- trait_summarise_boot_moments(
  np_bootstrapped_moments
)

sum_boot_moment

## ----fit-dist, echo=TRUE, eval=TRUE-------------------------------------------
# fit distributions
fitted_distributions <- trait_fit_distributions(
  filled_traits = trait_filling, 
  distribution_type = "lognormal"
)

fitted_distributions

## ----fit-dist2, echo=TRUE, eval=FALSE-----------------------------------------
#  # fit several types of distributions
#  fitted_distributions <- trait_fit_distributions(
#    filled_traits = trait_filling,
#    distribution_type = c(Plant_Height_cm = "normal",
#                             Wet_Mass_g = "lognormal")
#  )
#  
#  fitted_distributions

## ----para-boot, echo=TRUE, eval=TRUE------------------------------------------
# run parametric bootstrapping
p_bootstrapped_moments <- trait_parametric_bootstrap(
  fitted_distributions = fitted_distributions,
  nrep = 50
)

p_bootstrapped_moments

## ----raw-dist-np, echo=TRUE, eval=TRUE----------------------------------------
# run nonparametric bootstrapping
raw_dist_np <- trait_np_bootstrap(
  filled_traits = trait_filling,
  raw = TRUE
)

raw_dist_np

## ----plot-raw-dist-np, echo=TRUE, eval=TRUE, fig.width = 6--------------------
ggplot(raw_dist_np, aes(x = log(Value), fill = Site)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = col_palettes$GrandBudapest1) +
  labs(x = "log(trait value)") +
  facet_wrap(facets = vars(Trait), scales = "free")

## ----multivariate-prep, warning=FALSE-----------------------------------------
library(FD)

multivariate_traits <- trait_fill(
  comm = community |>
    # to make the example faster, we'll only use a subset of plots
    filter(
      Site == "1",
      PlotID %in% c("A", "B")
    ),
  traits = trait,
  scale_hierarchy = c("Site", "PlotID"),
  taxon_col = "Taxon",
  value_col = "Value",
  trait_col = "Trait",
  abundance_col = "Cover",
  complete_only = TRUE,
  leaf_id = "ID"
)

## ----multivariate-bootstrap, warning=FALSE------------------------------------
boot_multi <- trait_multivariate_bootstrap(
  filled_traits = multivariate_traits,
  nrep = 5, # number of reps is set low for demo purposes
  sample_size = 200,
  id = "ID",
  fun = function(x) {
    dbFD(
      x = x,
      calc.FRic = FALSE,
      calc.FDiv = FALSE,
      calc.CWM = FALSE,
      stand.x = FALSE,
      scale.RaoQ = FALSE
    )
  }
)

# data wrangling
raoq_est <- boot_multi |>
  mutate(result = map(result, as.data.frame)) |>
  unnest(result)

## ----multivariate-plot, warning=FALSE, fig.width = 6--------------------------
ggplot(data = raoq_est, mapping = aes(x = RaoQ, fill = PlotID)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = col_palettes$GrandBudapest1) +
  xlim(c(0, 6))

## ----multivariate-tpd, echo=TRUE, warning=FALSE, message=FALSE, fig.width = 6----
library(TPD)

boot_tpd <- trait_multivariate_bootstrap(multivariate_traits,
  id = "ID",
  nrep = 5, # Note that the number of reps is set low for demo purposes
  fun = function(x) {
    TPDs(
      species = rep(1, 200),
      traits = x
    ) |>
      REND(TPDs = _)
  }
)

# wrangling data
tpd <- boot_tpd |>
  mutate(result = map(result, as.data.frame)) |>
  unnest(result) |>
  pivot_longer(
    cols = species.FRichness:species.FDivergence,
    values_to = "value",
    names_to = "metric"
  ) |>
  mutate(metric = str_remove(metric, "species\\."))


# and make plot
ggplot(data = tpd) +
  geom_violin(
    mapping = aes(y = value, x = PlotID, fill = PlotID),
    alpha = 0.5,
    draw_quantiles = c(0.025, 0.975)
  ) +
  scale_fill_manual(values = col_palettes$GrandBudapest1) +
  facet_wrap(facets = vars(metric), scales = "free")

## ----coverage-plot, echo=TRUE, eval=TRUE, fig.width = 6-----------------------
# show coverage plot
autoplot(trait_filling) +
  scale_fill_manual(values = col_palettes$GrandBudapest1) +
  theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5))

## ----missing-traits, echo=TRUE, eval=TRUE-------------------------------------
# list missing traits
trait_missing(
  filled_trait = trait_filling,
  comm = community
)

