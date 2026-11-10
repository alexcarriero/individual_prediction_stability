library(arrow)
library(dplyr)
library(ggplot2)

evaluation <- open_dataset(
  "results/evaluation_v2/dgm=linear/sample_size=458"
) |>
  filter(
    method %in% c("independent", "bootstrap", "bootstrap_635"),
    model %in% c("glmnet", "random_forest")
  ) |>
  select(method, model, iter, id, pi_width) |>
  collect()

reference <- evaluation |>
  filter(method == "independent") |>
  select(model, iter, id, reference = pi_width)

plot_data <- evaluation |>
  filter(method != "independent") |>
  left_join(reference, by = c("model", "iter", "id")) |>
  mutate(
    model = factor(
      model,
      c("glmnet", "random_forest"),
      c("Elastic Net", "Random Forest")
    ),
    method = factor(
      method,
      c("bootstrap", "bootstrap_635"),
      c("Non-parametric\nbootstrap", "Out-of-bag\nbootstrap")
    )
  )

plot_limits <- range(c(plot_data$reference, plot_data$pi_width), na.rm = TRUE)

figure <- ggplot(plot_data, aes(reference, pi_width)) +
  geom_bin_2d(bins = 100) +
  geom_abline(slope = 1, intercept = 0, colour = grey(0.5), linewidth = 0.3) +
  facet_grid(rows = vars(method), cols = vars(model)) +
  coord_equal(xlim = plot_limits, ylim = plot_limits) +
  scale_fill_viridis_c(trans = "log10") +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  labs(
    title = "95% prediction range width: Bootstrap vs independent sampling",
    x = "Independent sampling estimate",
    y = "Bootstrap estimate",
    fill = "Count"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = grey(0.95), colour = NA),
    legend.position = "inside",
    legend.justification = c(1, 0.01),
    legend.background = element_rect(fill = "white", color = "transparent"),
    plot.title = element_text(
      size = 16,
      face = "bold",
      margin = margin(b = 10, t = 10)
    ),
    plot.title.position = "plot",
    legend.key.size = unit(0.4, "cm"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )
figure

ggsave(
  "results/social_media_figure.png",
  figure,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)
