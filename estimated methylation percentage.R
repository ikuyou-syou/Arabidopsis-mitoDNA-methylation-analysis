library(ggplot2)
library(dplyr)
# Define the data
data <- premitomethyl

# Filter data
filtered_data <- data %>% 
  dplyr::filter(readssum >= 10)

# Group the data by and assign a group ID
filtered_data <- filtered_data %>% 
  mutate(group = (end - 1) %/% 1000 + 1) 

# Calculate the methylation level for each group 
smoothed_data <- filtered_data %>%
  group_by(group, C_type) %>%
  summarise(
    methylated_reads = sum(methylated_reads),
    unmethylated_reads = sum(unmethylated_reads),
    methylation_level = (methylated_reads / (methylated_reads + unmethylated_reads)) * 100,
    end = mean(end)  
  ) %>%
  arrange(end)

# Plot the data
al2 <- ggplot(smoothed_data, aes(x = end, y = methylation_level, color = C_type)) +
  geom_line(size = 0.75) +
  scale_color_manual(values = c("CG" = "darkorchid4", "CHG" = "cyan4", "CHH" = "yellow3")) +
  labs(x = "Position", y = "Methylation(%)") + ylim(0, 100) +
  theme_bw() +
  theme(panel.border = element_rect(fill = NA, color = "black", size = 1, linetype = "solid"))

# Display the plot
print(al2)


