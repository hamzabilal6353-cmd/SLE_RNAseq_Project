# Load plotting packages
library(enrichplot)
library(ggplot2)

# Load saved pathway results
go_upregulated <- readRDS(
  "results/tables/GO_BP_upregulated_object.rds"
)

go_downregulated <- readRDS(
  "results/tables/GO_BP_downregulated_object.rds"
)

# Function for wrapping long pathway names
wrap_pathway_names <- function(x, width = 32) {
  
  vapply(
    x,
    function(label) {
      paste(
        strwrap(label, width = width),
        collapse = "\n"
      )
    },
    character(1)
  )
}

# Clean upregulated pathway plot
go_up_clean <- enrichplot::dotplot(
  go_upregulated,
  showCategory = 10,
  font.size = 11
) +
  scale_y_discrete(
    labels = wrap_pathway_names
  ) +
  labs(
    title = "GO Biological Processes Enriched in SLE",
    subtitle = "Upregulated genes: SLE versus healthy controls",
    x = "Gene ratio",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.y = element_text(
      size = 10,
      lineheight = 0.9
    ),
    plot.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      15, 25, 15, 15
    )
  )

# Display and save upregulated plot
print(go_up_clean)

ggsave(
  "results/figures/GO_BP_upregulated_clean.png",
  plot = go_up_clean,
  width = 12,
  height = 9,
  dpi = 300,
  bg = "white"
)

# Clean downregulated pathway plot
go_down_clean <- enrichplot::dotplot(
  go_downregulated,
  showCategory = 10,
  font.size = 11
) +
  scale_y_discrete(
    labels = wrap_pathway_names
  ) +
  labs(
    title = "GO Biological Processes Reduced in SLE",
    subtitle = "Downregulated genes: SLE versus healthy controls",
    x = "Gene ratio",
    y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.y = element_text(
      size = 10,
      lineheight = 0.9
    ),
    plot.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      15, 25, 15, 15
    )
  )

# Display and save downregulated plot
print(go_down_clean)

ggsave(
  "results/figures/GO_BP_downregulated_clean.png",
  plot = go_down_clean,
  width = 12,
  height = 9,
  dpi = 300,
  bg = "white"
)

cat(
  "Clean upregulated plot saved:",
  file.exists(
    "results/figures/GO_BP_upregulated_clean.png"
  ),
  "\n"
)

cat(
  "Clean downregulated plot saved:",
  file.exists(
    "results/figures/GO_BP_downregulated_clean.png"
  ),
  "\n"
)

print("Clean pathway plots completed successfully")