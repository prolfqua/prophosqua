library(tidyverse)
data <-read_tsv("sim2_results.tsv")

# Add a unique identifier for each row

# Convert data to long format for plotting
# here instead of  2 cols w TP and FP  we add a column for Type and take counts appart
data_long <- data %>%
  gather(key = "Type", value = "Count", TP, FP)

data_long <- data_long |>
  group_by(method, Type) |>
  mutate(ID = row_number()) |> ungroup()

# Plot
pdf("StackedBars_comparison_sim2.pdf",3,3)
ggplot(data_long, aes(x = factor(ID), y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Group Size", y = "Count", fill = "Type") +
  scale_x_discrete(labels = data$GrpSize) +
  theme_minimal() +
  facet_wrap(~method, ncol = 1)
dev.off()
