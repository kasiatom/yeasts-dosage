require(readr)
require(dplyr)
require(ggplot2)
require(ggdendro)
require(ape)

set.seed(0)
# setwd("$HOME/distances") ## put here name of the working directory with all distance metrics files

sdiff <- read_delim("1011.sdiff.summary",
  delim = "\t", escape_double = FALSE,
  trim_ws = TRUE)
n_sdiff <-colnames(sdiff) 
n_sdiff[1] <-"IID1"
colnames(sdiff) <- n_sdiff

king <- read_delim("1011-king.kin0", delim = "\t",
  escape_double = FALSE, trim_ws = TRUE)
n_king <-colnames(king) 
n_king[2] <-"IID1"
n_king[1] <-"IID2"
colnames(king) <- n_king
 
## plink 1.9 1-ibs
ibs_distance <- read_delim("1011-ibs-distance.txt",
  delim = "\t", escape_double = FALSE,
  col_names = FALSE, trim_ws = TRUE)

colnames(ibs_distance) <- c("IID1", "IID2", "DISTANCE")

ibs_distance <- ibs_distance %>%
  filter(DISTANCE > 0)

quan <- quantile(ibs_distance$DISTANCE, probs=seq(0,1, 0.05))
cutoff <- quan[["10%"]]

mean_ibs <- ibs_distance

## plot before
plot <- mean_ibs %>%
  ggplot(aes(x=DISTANCE)) +
  geom_histogram(bins=50, fill = "orange", color = "grey") + 
  theme_bw() +
  ggtitle("PAIRWISE IBS DIASTANCEs - ALL STRAINS")
ggsave("distances-all.png", device = "png")
###

i <- 301
j <- 300
while (i >= 100 & j > 0){
  mean_ibs <-mean_ibs %>%
    group_by(IID1) %>%
    filter(DISTANCE != 0) %>%
    mutate("MEAN" = mean(DISTANCE)) %>%
    mutate("SD" = sd(DISTANCE)) %>%
    mutate("MIN" = min(DISTANCE)) %>%
    mutate("MAX" = max(DISTANCE)) %>%
    mutate("CLOSE" = if_else(DISTANCE <= cutoff, 1, 0)) %>%
    mutate(CLOSE = sum(CLOSE)) %>%
    ungroup()
  

  md <-min(mean_ibs$DISTANCE)
  ids <- mean_ibs %>% 
    filter(DISTANCE == md) %>%
    arrange(IID1)%>%
    as.data.frame()

  ids <- ids[1:2,] %>%
    arrange(CLOSE, rev(MEAN)) %>%
    select(IID1) 
  
  id <- ids[2,1]

  mean_ibs <- mean_ibs %>%       
    filter(IID1 != id & IID2 != id)
 
  j <- max(mean_ibs$CLOSE) 
  i <- length(unique(c(mean_ibs$IID1, mean_ibs$IID2)))
} 
###
## plot after
plot1 <- mean_ibs %>%
  ggplot(aes(x=DISTANCE)) +
  geom_histogram(bins=50, fill = "orange", color = "grey") + 
  theme_bw() +
  ggtitle("PAIRWISE IBS DIASTANCES - SELECTED STRAINS")

ggsave("distances-selected.png", device = "png")

ids <- unique(c(mean_ibs$IID1, mean_ibs$IID2))
write(ids, "selected-ids.txt")

## other metrics
merged <- sdiff %>%
  left_join(., king, by = c("IID1", "IID2"))

merged <- merged %>%
  left_join(., ibs_distance, by = c("IID1", "IID2"))


plot <-merged %>%
  ggplot(aes(x=KINSHIP)) +
  geom_histogram(bins=50, fill = "orange", color = "grey") + 
  theme_bw()
ggsave("kinship.png", device = "png")

plot <- merged %>%
  ggplot(aes(x=DIFF_CT/OBS_CT)) +
  geom_histogram(bins=50, fill = "orange", color = "grey") + 
  theme_bw()
ggsave("sdiff.png", device = "png")

plot <-merged %>%
  ggplot(aes(x=DIFF_CT/OBS_CT, y = DISTANCE)) +
  geom_point(shape = 21, fill = "orange", color = "grey") +
  theme_bw()
ggsave("sdiff-ibsdist.png", device = "png")


plot <-merged %>%
  ggplot(aes(x=DISTANCE, y = KINSHIP)) +
  geom_point(shape = 21, fill = "orange", color = "grey") +
  theme_bw()
ggsave("ibsdist-kinship.png", device = "png")


plot <- merged %>%
  ggplot(aes(x=DISTANCE, y = HETHET)) +
  geom_point(shape = 21, fill = "orange", color = "grey") +
  theme_bw()
ggsave("ibsdist-hethet.png", device = "png")

plot <- merged %>%
  ggplot(aes(x=DISTANCE, y = IBS0)) +
  geom_point(shape = 21, fill = "orange", color = "grey") +
  theme_bw()

ggsave("ibsdist-ibs0.png", device = "png")

## cluster
mibs <- read_delim("1011-dist.mdist", delim = "\t",
                   escape_double = FALSE, col_names = FALSE,
                   trim_ws = TRUE)
labels <- read_delim("1011-dist.mdist.id",
  delim = "\t", escape_double = FALSE,
  col_names = FALSE, trim_ws = TRUE)

colnames(mibs) <- labels$X1
rownames(mibs) <- labels$X1


ward <- hclust(as.dist(mibs), method = "ward.D2")
ward_d <- as.dendrogram(ward)
complete <- hclust(as.dist(mibs), method = "complete")
complete_d <- as.dendrogram(complete)
upgma <- hclust(as.dist(mibs), method = "average")
upgma_d <- as.dendrogram(complete)


#ggdendrogram(hc, size = 0.2, color = "orange") + 
#  theme(axis.text=element_text(size=4),
#        axis.title=element_text(size=4)) 

pdf(file = "ward-dendogram.pdf", width = 8.27, height = 11.69)
par("mar")
par(mar=c(2, 0, 0, 0))

nodePar <- list(lab.cex = 0.3, pch = c(NA, 19), 
                cex = 0.2, col = "orange", cex.axis = 0.4, cex.lab = 0.3)
plot(ward_d, nodePar = nodePar, horiz = TRUE, edgePar = list(col = "orange", lwd = 0.5))
dev.off()

pdf(file = "complete-dendogram.pdf", width = 8.27, height = 11.69)
par("mar")
par(mar=c(2, 0, 0, 0))

nodePar <- list(lab.cex = 0.3, pch = c(NA, 19), 
                cex = 0.2, col = "orange", cex.axis = 0.4, cex.lab = 0.3)
plot(complete_d, nodePar = nodePar, horiz = TRUE, edgePar = list(col = "orange", lwd = 0.5))
dev.off()

pdf(file = "upgma-dendogram.pdf", width = 8.27, height = 11.69)
par("mar")
par(mar=c(2, 0, 0, 0))

nodePar <- list(lab.cex = 0.3, pch = c(NA, 19), 
                cex = 0.2, col = "orange", cex.axis = 0.4, cex.lab = 0.3)
plot(upgma_d, nodePar = nodePar, horiz = TRUE, edgePar = list(col = "orange", lwd = 0.5))
dev.off()

groups <- cutree(upgma, h = 0.02)

my_strains <-c()
for (i in seq(1, 98, 1)) {
  s <- sample(names(groups[groups == i]), size =1)
  my_strains <-c(my_strains, s)
}

write(my_strains, "selected-ids-dendogram.txt")


## njt


nj_tree <- nj(as.dist(mibs))
pdf(file = "nj.pdf", width = 8.27, height = 11.69)  
par("mar")
par(mar=c(2, 0, 0, 0))
tips <- nj_tree$tip.label
tips <- as.data.frame(tips)
tips <- tips %>% mutate("color" = if_else(tips %in% my_strains, "red", "blue"))
label_colours <- tips$color

nj_plot <- plot(nj_tree, cex = 0.3, edge.color = "grey", node.color = "grey", tip.color = label_colours, type = "fan")
dev.off()


subtrees <- subtrees(nj_tree)
