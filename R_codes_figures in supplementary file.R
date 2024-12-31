
## packages
library(ggplot2)
library(ggrepel)

#### Figure S2: Distribution of contents in pH, soil organic matter, nitrogen, aluminium, sodium, potassium, phosphorus, and density of particles by soil depths (0-10, 10-20, and 20-30 cm) in two transects in an Eastern Amazonian forest.

solo<- read.table("supplementary_figure_complete_soil_samples_dataset.csv", header=TRUE, sep="")

solo <- solo[which (solo$Prof %in% c("0_10", "10_20", "20_30")),]

par(mfrow=c(4,2))
par(mar = c(2,2,2,2))
dotchart(solo$pH, main="pH", groups = as.factor (solo$Prof))
dotchart(solo$MO, main="Organic Matter", groups = as.factor (solo$Prof))
dotchart(solo$N, main="Nitrogen", groups = as.factor (solo$Prof))
dotchart(solo$Al, main="Aluminium", groups = as.factor (solo$Prof))
dotchart(solo$Na_Cl, main="Sodium", groups = as.factor (solo$Prof))
dotchart(solo$K, main="Potassium", groups = as.factor (solo$Prof))
dotchart(solo$P, main="Phosphorus", groups = as.factor (solo$Prof))
dotchart(solo$dens_part, main="Density of particles", groups = as.factor (solo$Prof))

dev2bitmap (file = "dotchart_solo_prof", type = "png16m",
            height = 7, width = 7, res = 300,
            units = "in",
            method = "postscript")


#### Figure S3: Principal Component Analysis of morphological traits for 48 below-ground ant species recorded in two transects in an Eastern Amazonian forest

traits.std <- read.table("dataset_standardized_ant_traits.txt", header = TRUE, sep="\t")
groups <- read.table("supplementary_figure_species_classified_in_groups.txt", header = TRUE, sep="\t")
groups$group <- as.factor(groups$group)

pca<- vegan::rda(traits.std[, 2:11], center = TRUE)
summary(pca)
ylabel<- "PCA 2 (48.4%)"
xlabel<- "PCA 1 (21.8%)"

with(groups,levels (group))
colvec <- c("gray","black", "red") #colocar na ordem 

### PCA plot

BD_unit <- data.frame(pca$CA$u)[, 1:2]
BD_unit$ssp <- row.names(pca.scores)
BD_unit$ssp == groups$species
BD_unit$code <-  groups$ID_Rony
BD_unit$group <- as.factor(groups$group)

BD_arrow <- data.frame(pca$CA$v[,1:2])
BD_arrow$arrow <- gsub(".std","",row.names(BD_arrow))


cor_vec <- "red"
cor_ssp <- 1
adj2 <- 0.65


ggplot(data=BD_unit, aes(x = PC1, y = PC2, colour=group))+
  geom_hline(yintercept = 0, colour="gray", linetype="dashed")+ geom_vline(xintercept = 0, colour="gray", linetype="dashed")+
  geom_segment(data=BD_arrow, aes(x = 0, xend = PC1*adj2, y = 0, yend = PC2*adj2), arrow = arrow(length = unit(0.25, "cm")), colour = "darksalmon", size = .8)+
  geom_point(size=3, alpha=.6, shape = 16)+
  geom_text(data=BD_arrow, aes(x = PC1*adj2, y = PC2*adj2, label = arrow), size = 4, hjust=-.1, vjust=-.5, colour=cor_vec)+ #Segmento
  geom_text_repel(aes(label = code), hjust = .1, vjust=.6, segment.size = 0.02, segment.color = "gray90", colour=cor_ssp, size = 2.5, box.padding = -0.01, alpha=0.8) + #Especies
  scale_x_continuous(sec.axis = sec_axis(~./adj2))+
  scale_y_continuous(sec.axis = sec_axis(~./adj2))+
  scale_color_manual(values = colvec)+
  labs(x=xlabel, y=ylabel)+
  theme_classic(base_size = 13)+
  theme(legend.position = "none",
        legend.background = element_rect(color = "black", size = .5, linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.1),
        axis.text.x.bottom = element_text(color = cor_ssp),
        axis.text.y.left = element_text(color = cor_ssp),
        axis.text.x.top = element_text(color = cor_vec),
        axis.text.y.right = element_text(color = cor_vec))

#ggsave("Figure_S3.tiff", width = 9, height = 7, dpi = 300, compression = "lzw")


# Additional data visualization: plotting species colored by groups (ground-dwelling, leaf-litter habitat or subterranean species)

p2_pca <- ggplot(data=BD_unit, aes(x = PC1, y = PC2, colour=group))+
  geom_hline(yintercept = 0, colour="gray", linetype="dashed")+ geom_vline(xintercept = 0, colour="gray", linetype="dashed")+
  geom_point(size=3, shape = 16)+
  geom_text_repel(aes(label = code), hjust = .1, vjust=.6, segment.size = 0.02, segment.color = "gray90", colour=cor_ssp, size = 2.5, box.padding = -0.01, alpha=0.8) + #Especies
  scale_color_manual(values = colvec)+
  labs(x=xlabel, y=ylabel)+
  theme_classic(base_size = 13)+
  theme(legend.position = "none",
        legend.background = element_rect(color = "black", size = .5, linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.1),
        axis.text.x.bottom = element_text(color = cor_ssp),
        axis.text.y.left = element_text(color = cor_ssp))
p2_pca
#ggsave("PCA_2.tiff", width = 9, height = 7, dpi = 300, compression = "lzw")


# Additional data visualization: biplot (plotting species colored by groups: ground-dwelling, leaf-litter habitat or subterranean species)

p3_pca <- ggplot(data=BD_unit, aes(x = PC1, y = PC2, colour=group))+
  geom_hline(yintercept = 0, colour="gray", linetype="dashed")+ geom_vline(xintercept = 0, colour="gray", linetype="dashed")+
  geom_point(size=3, shape = 16)+
  geom_segment(data=BD_arrow, aes(x = 0, xend = PC1*adj2, y = 0, yend = PC2*adj2), arrow = arrow(length = unit(0.25, "cm")), colour = "red", size = .8)+
  geom_text(data=BD_arrow, aes(x = PC1*adj2, y = PC2*adj2, label = arrow), size = 4, hjust=-.1, vjust=-.5, colour=cor_vec)+ #Segment
  scale_x_continuous(sec.axis = sec_axis(~./adj2))+
  scale_y_continuous(sec.axis = sec_axis(~./adj2))+
  scale_color_manual(values = colvec)+
  labs(x=xlabel, y=ylabel)+
  theme_classic(base_size = 13)+
  theme(legend.position = "none",
        legend.background = element_rect(color = "black", size = .5, linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.1),
        axis.text.x.bottom = element_text(color = cor_ssp),
        axis.text.y.left = element_text(color = cor_ssp),
        axis.text.x.top = element_text(color = cor_vec),
        axis.text.y.right = element_text(color = cor_vec))
p3_pca

library(cowplot)
plot_grid(p2_pca+labs(x=NULL), p3_pca, ncol = 1)
#ggsave("PCA_3.tiff", width = 8, height = 11, dpi = 300, compression = "lzw")

### PCA biplot (species points not colored)

S3 <- ggplot(data = BD_unit, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, colour = "gray", linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "gray", linetype = "dashed") +
  geom_point(data = BD_unit, aes(x = PC1, y = PC2), size = 3, alpha = 0.2, shape = 16) + 
  geom_text_repel(data = BD_unit, aes(x=PC1, y=PC2, label= traits.std$species_code), fontface="italic", colour=1,
                  segment.color = '#cccccc', size = 3, alpha=.7, max.overlaps = 15) +
  geom_text(data = BD_arrow, aes(x = PC1, y = PC2, label = rownames(BD_arrow)), col = 'red') +
  geom_segment(data = BD_arrow, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'darkred') +
  labs(x = xlabel, y = ylabel) +
  theme_bw(base_size = 14) + 
  theme (panel.grid = element_blank())

S3



#### Figure S4: correlations between soil variables

#Continuous figure function:
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point(aes(shape=Area), size=2, alpha=0.6) +
    geom_smooth(method=lm, formula = y ~ x, color = 1, se = F)
  p
}

#Diagonal figure function:
my_fn2 <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_histogram(bins = 10, colour="gray20")
  p
}

#Correlation figure function:
cor_fun <- function(data, mapping, method="pearson", ndp=2, sz=7, stars=FALSE, ...){
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  corr <- cor.test(x, y, method=method)
  est <- corr$estimate
  lb.size <- sz
  
  if(stars){
    stars <- c("***", "**", "*", "")[findInterval(corr$p.value, c(0, 0.001, 0.01, 0.05, 1))]
    lbl <- paste0(round(est, ndp), stars)
  }else{
    lbl <- round(est, ndp)
  }
  ggplot(data=data, mapping=mapping) + 
    annotate("text", x=mean(x, na.rm=TRUE), y=mean(y, na.rm=TRUE), label=lbl, size=lb.size,...)+
    theme(panel.grid = element_blank())
}


## - Plot the Figure ----
MyVar<- c("pH", "C", "MO", "N", "C_N", "Al", "Na_Cl", "K", "P", "dens_part")
var_pairs <- data.frame(Area=df$Area, df[,MyVar])
colnames(var_pairs) <- c("Area", "pH", "C", "OM", "N", "C:N", "Al", "Na", "K", "P", "Dp")

library(GGally)
library(ggplot2)

p1 <- ggpairs(var_pairs, columns = 2:11,
              lower = list(continuous = my_fn),
              diag = list(continuous = my_fn2), 
              upper = list(continuous = cor_fun))

cor_m <- p1+theme_minimal()+
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(),
        axis.line.x.bottom = element_line(color = "black", size = .8),
        axis.text.x.bottom = element_text(color = "black"),
        axis.line.y.left = element_line(color = "black", size = .8),
        axis.text.y.left = element_text(color = "black"),
        panel.border = element_rect(linetype = "solid", colour = "gray50", fill = NA),
        strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12),
        strip.background = element_rect(colour="black", fill="white"))
cor_m

#ggsave("Fig_S1.png", width = 15, height = 10, dpi = 400)
#ggsave("Fig_S1.tiff", width = 15, height = 10, dpi = 400, compression = 'lzw')


#### Figure S5: Distribution of contents in sand, silt, and clay in two studied transects

soil<- read.table("supplementary_figures_soil_granulometry.txt", header = TRUE, sep="")
soil$fTransect<- as.factor (soil$Transect)
par(mfrow=c(2,2))
dotchart(soil$Areia, groups= soil$fTransect, main="Sand content", ylab="")
dotchart(soil$Argila, groups= soil$fTransect, main="Clay content", xlab="")
dotchart(soil$Silte, groups= soil$fTransect, main="Silt content", xlab="")


#### Figure S6: Relationships between observed species richness of below-ground ants and soil variables 

library(ggpubr)
MyNames <- c("P", "Na", "Dp", "pH", "Al", "N", "OM", "K")
MyVar <- c("P", "Na_Cl", "dens_part", "pH", "Al", "N", "MO", "K")
df2 <- df[, MyVar]
df2$Sobs <- df$Sobs
df2$Area <- df$Area


gg_var <- function(coNum){
  p1 <- ggplot(df2, aes(x = df2[,coNum], y = Sobs)) +
    stat_smooth(formula = y ~ x, method = "lm", fullrange=T, alpha = 0.2, size = 1.5, color=1)+
    geom_point(aes(shape=Area), size=3, alpha=0.7)+
    stat_regline_equation(
      aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), formula = y ~ x) +
    labs(x = MyNames[coNum], y = "Sobs")+
    theme_classic(base_size = 14)+
    theme(legend.position = "none",
          plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line (size=.4),
          axis.line.x.bottom = element_line(color = "black", size = .6),
          axis.text.x.bottom = element_text(color = "black"),
          axis.line.y.left = element_line(color = "black", size = .6),
          axis.text.y.left = element_text(color = "black"))
  return(p1)
}

library(cowplot)
plot_grid(gg_var(1), gg_var(2)+labs(y = NULL),
          gg_var(3), gg_var(4)+labs(y = NULL),
          gg_var(5), gg_var(6)+labs(y = NULL),
          gg_var(7), gg_var(8)+labs(y = NULL),
          ncol = 2, align = "hv")

#ggsave("Fig_S2.2.png", width = 12, height = 14, dpi = 400)
#ggsave("Fig_S2.2.pdf", width = 12, height = 14)
#ggsave("Fig_S2.2.tiff", width = 12, height = 14, dpi = 400, compression = 'lzw')


