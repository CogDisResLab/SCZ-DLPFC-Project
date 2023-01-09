# Look at Kinase module gene expressions

library(tidyverse)

kinase_expression <- read_csv("raw/Kinase-Expression.csv") |>
  drop_na(Name) |>
  select(
    Name,
    HGNC_Symbol,
    Lewis_2017_L5,
    Lewis_2015_L5,
    Deep_Neurons,
    Superficial_Neurons,
    Superficial_Deep_Neurons,
    MtSinaiDLPFC,
    Iwamoto_BA46_SCZ
  )

erk_genes <-
  c("MEKK2",
    "MEKK3",
    "TPL2",
    "ERK5",
    "RAF1",
    "RAFA",
    "RAFB",
    "MEK1",
    "MEK2",
    "ERK1",
    "ERK2")
jnk_genes <-
  c(
    "MEKK1",
    "MEKK2",
    "MEKK3",
    "DLK",
    "MLK2",
    "TPL2",
    "ASK1",
    "TAK1",
    "TAO1",
    "TAO2",
    "MEK4",
    "MEK7",
    "JNK1",
    "JNK2",
    "JNK3"
  )
p38_genes <-
  c(
    "MEKK1",
    "MEKK2",
    "MEKK3",
    "DLK",
    "MLK2",
    "TPL2",
    "ASK1",
    "TAK1",
    "TAO1",
    "TAO2",
    "MEK3",
    "MEK6",
    "P38A",
    "P38B",
    "P38G",
    "P38D"
  )

p38_module <- kinase_expression |>
  filter(Name %in% p38_genes) |>
  mutate(Name = str_glue("{Name} ({HGNC_Symbol})")) |>
  select(-HGNC_Symbol) |>
  pivot_longer(where(is.numeric), names_to = "Dataset", values_to = "LogFC") |>
  unique() |>
  write_csv("data/P38-Module-Gene-Expression.csv")

erk_module <- kinase_expression |>
  filter(Name %in% erk_genes) |>
  mutate(Name = str_glue("{Name} ({HGNC_Symbol})")) |>
  select(-HGNC_Symbol) |>
  pivot_longer(where(is.numeric), names_to = "Dataset", values_to = "LogFC") |>
  unique() |>
  write_csv("data/ERK-Module-Gene-Expression.csv")

jnk_module <- kinase_expression |>
  filter(Name %in% jnk_genes) |>
  mutate(Name = str_glue("{Name} ({HGNC_Symbol})")) |>
  select(-HGNC_Symbol) |>
  pivot_longer(where(is.numeric), names_to = "Dataset", values_to = "LogFC") |>
  unique() |>
  write_csv("data/JNK-Module-Gene-Expression.csv")

make_heatmap <- function(data, kinase) {
  g <-
    ggplot(data, aes(
      y = Name,
      x = Dataset,
      fill = LogFC,
      label = LogFC
    ))

  p <- g + geom_tile(width = 0.9,
                     height = 0.9,
                     color = "grey80") +
    geom_text() +
    scale_fill_gradient2(
      low = "darkred",
      high = "darkgreen",
      breaks = seq(-0.3, 0.6, 0.3),
      limits = c(-0.3, 0.6),
      na.value = "grey60"
    ) +
    scale_x_discrete(
      breaks = c(
        "Lewis_2015_L5",
        "Lewis_2017_L5",
        "Superficial_Neurons",
        "Superficial_Deep_Neurons",
        "Deep_Neurons",
        "MtSinaiDLPFC",
        "Iwamoto_BA46_SCZ"
      ),
      limits = c(
        "Lewis_2015_L5",
        "Lewis_2017_L5",
        "Superficial_Neurons",
        "Superficial_Deep_Neurons",
        "Deep_Neurons",
        "MtSinaiDLPFC",
        "Iwamoto_BA46_SCZ"
      ),
      labels = c(
        "Lewis (2015)",
        "Lewis (2017)",
        "Superficial Neurons",
        "Superficial Deep Neurons",
        "Deep Neurons",
        "Mt. Sinai DLPFC",
        "Iwamoto DLPFC"
      )
    ) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    xlab("Datasets") + ylab("Kinase (HGNC Symbol)") +
    ggtitle(str_glue("Gene expression for {kinase} module"))

  p
}

modules <- list(P38 = p38_module,
                JNK = jnk_module,
                ERK = erk_module)

figures <- map2(modules, names(modules),
                ~ make_heatmap(.x, .y)) |>
  {
    \(l) map2(l,
              names(l),
              ~ ggsave(
                file.path("figures",
                          str_glue("{.y}-Expression-Heatmap.png")),
                plot = .x,
                width = 12,
                height = 8,
                bg = "white"

              ))
  }()

p38_matrix <- p38_module |>
  pivot_wider(names_from = "Dataset",
              values_from = "LogFC") |>
  column_to_rownames("Name") |>
  as.matrix()

p38_cor <- p38_matrix |>
  cor(use = "complete.obs", method = "pearson")

jnk_matrix <- jnk_module |>
  pivot_wider(names_from = "Dataset",
              values_from = "LogFC") |>
  column_to_rownames("Name") |>
  as.matrix()

jnk_cor <- jnk_matrix |>
  cor(use = "complete.obs", method = "pearson")

erk_matrix <- erk_module |>
  pivot_wider(names_from = "Dataset",
              values_from = "LogFC") |>
  column_to_rownames("Name") |>
  as.matrix()

erk_cor <- erk_matrix |>
  cor(use = "complete.obs", method = "pearson")

