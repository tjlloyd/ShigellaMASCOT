######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library(coda)  # needed to calculate ESS values
library("methods")
require(ggtree)
library(ggpubr)
library(treeio)

# clear workspace
rm(list = ls())

# Set the directory to the directory of the file
setwd('/Users/yourdirectoryfilepathhere')

# Define colors for locations
color_palette <- c(
  "location" = "black",
  "anotherlocation" = 'orange',
  "morelocations" = "skyblue",
  "evenmorelocations" = "#CC79A7",
  "lastlocation" = 'darkgreen'
)



# Read the tree data
tree = read.beast("yourtreefilehere.tree")

# Define locations each location should be listed in quotes
loc = c("locationhere", "anotherlocation", "morelocations", "evenmorelocations", "lastlocation")
tree@data$loc = vector("character", length(tree@data$locationhere))
tree@data$site = vector("character", length(tree@data$locationhere))
tree@data$entropy = vector("numeric", length(tree@data$locationhere))

# read in metadata.csv
metadata = read.table("metadata.csv", header = TRUE, sep = ",")
# remove duplicates PNUSAE rows if needed
#metadata = metadata[!duplicated(metadata$PNUSAE),]

# Calculate the most probable location and entropy for each node
for (ii in seq(1, length(tree@data$Locationhere))) {
  prob_vals = numeric(length(loc))
  entr = 0
  for (j in seq(1, length(loc))) {
    prob_vals[j] = as.numeric(tree@data[[loc[j]]][ii])
    if (prob_vals[j] != 0) {
      entr = entr + prob_vals[j] * log(prob_vals[j])
    }
  }
  tree@data$loc[ii] = loc[which.max(prob_vals)]
  tree@data$entropy[ii] = -entr
}
for (ii in seq(1, length(tree@phylo$edge[,1]))) {
  jj = as.numeric(tree@data$node[ii])
  if (jj <= length(tree@phylo$tip.label)){
    tmp = strsplit(tree@phylo$tip.label[jj], "\\|")[[1]]
    # find metadata$PNUSAE == tmp[[1]], if it's there, assgnsite as metadata$site
    if (length(which(metadata$PNUSAE == tmp[[1]])) > 0){
      if (!grepl("Unknown", metadata$site[which(metadata$PNUSAE == tmp[[1]])])){
        print(tmp)
        tree@data$site[ii] = metadata$site[which(metadata$PNUSAE == tmp[[1]])]  
        print(tree@data$site[ii])
      }
    }
  }
}


# Function to convert year fractions to date format
year_fraction_to_date <- function(year_fraction) {
  year <- as.integer(year_fraction)
  month <- round((year_fraction - year) * 12 + 1)
  as.Date(paste(year, month, "01", sep = "-"), format = "%Y-%m-%d")
}

# Plot the tree with correct coloring and x-axis labels in months
p = ggtree(tree, aes(color = loc, fill = loc, alpha = -entropy),size = 0.8, mrsd = "2024-03-01") + 
  geom_tippoint(aes(fill = loc), size = 1.75) + 
  theme_tree() +
  theme_minimal() +
  coord_flip() +
  scale_x_reverse(breaks = seq(2022, 2024, by = 1/3),
                  labels = function(x) format(year_fraction_to_date(x), "%b %Y")
  ) +
  scale_y_discrete(breaks = NULL) +
  scale_fill_manual(name = "Location", values = color_palette, labels = c("locationhere", "anotherlocation", "morelocations", "evenmorelocations", "lastlocation")) +
  scale_color_manual(name = "Location" , values = color_palette, labels = c("locationhere", "anotherlocation", "morelocations", "evenmorelocations", "lastlocation")) + 
  scale_shape_manual() +
  scale_size_manual() +
  scale_alpha(name = "Uncertainty") +
  theme(legend.position = c(0.8, 0.7), legend.box = "horizontal") 

# Display the plot
plot(p)

#export the plot
ggsave(filename = "yourfinalplotfilename.pdf" , plot = p , height = 6 , width = 8)
