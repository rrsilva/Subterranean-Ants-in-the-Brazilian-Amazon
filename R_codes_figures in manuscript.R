
#############################################################################################################
#                                                                                                           #
#                                    Figures in the main text                                               #
#                                                                                                           #
#############################################################################################################

### Plotting below figures depends on objects created on codes for statistical analyses.

#### Fig. 2. The relative importance of each variable for best fitting to the relationship between species richness and soil properties.

# Saving the results to export
mod_func <- function(mod_ac, name_v){
  models_tot <- mod_ac$msTable
  tab_f2 <- data.frame(models_tot)
  tab_f2 <- tab_f2[tab_f2$delta<=2,]
  tab_f2$adjR <- selec[1:nrow(tab_f2),'adjR^2']
  tab_f2$varY <- rep(name_v, nrow(tab_f2))
  var_names <- rownames(data.frame(mod_ac$sw))
  var_names <- var_names[order(var_names)]
  tab_f2$varX <- row.names(tab_f2)
  tab_f2$names_varX <- rep(paste(paste(1:length(var_names), var_names, sep = "="), collapse = " "), nrow(tab_f2))
  return(tab_f2)
}

imp_func <- function(selecao, name_v){
  imp2 <- data.frame(sw(selecao))
  imp2$var <- row.names(imp2)
  imp2 <- data.frame(imp2, rep(name_v, nrow(imp2)))
  colnames(imp2) <- c("importance", "varX", "varY")
  return(imp2)
}


name_var <- "Sobs"
tab_mod_temp <- mod_func(mod_ac=mod_acaic, name_v=name_var)
tab_mod_fin <- tab_mod_temp

tab_impor_temp <- imp_func(selecao=selec, name_v=name_var)
tab_impor_fin <- tab_impor_temp


import_var_a <- data.frame(sw(selec))
row.names(import_var_a)
import_var_a$nomes_var_a <- c("P", "Na", "Dp", "pH", "Al", "N", "OM", "K")

import_var_a$import_a <- round(import_var_a$sw.selec., 3)
import_var_a2 <- import_var_a[order(import_var_a$import_a, decreasing = T),]

### Figure 2 (a)
library(ggplot2)
p1 <- ggplot(data=import_var_a2, aes(x=nomes_var_a, y=import_a)) +
  geom_bar(position="dodge", stat="identity", colour="black") + 
  geom_hline(aes(yintercept=0.8),color="red", linetype="dashed", size=1)+
  coord_flip()+
  labs(x = "Variables", y = "Importance", tag = "A")+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1))+
  scale_x_discrete(limits=rev(c("P", "Na", "Dp", "pH", "Al", "N", "OM", "K")), labels=rev(c("P", "Na", "Dp", "pH", "Al", "N", "OM", "K")))+
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

p1


#interaction model (see Fig. 2 (b) ====

MyVar <- c("Na_Cl", "P", "Sobs")
x_y_ric <- df[, MyVar]
MyNames <- c("x", "y", "z")
colnames(x_y_ric) <- MyNames

h <- c(MASS::bandwidth.nrd(x_y_ric$x), MASS::bandwidth.nrd(x_y_ric$y))
dens <- MASS::kde2d(x_y_ric$x, x_y_ric$y, h = h, n = 100, lims = c(range(x_y_ric$x), range(x_y_ric$y)))
df1 <- expand.grid(x = dens$x, y = dens$y)
df1$z <- as.vector(dens$z)
n_bins <- 8


p4 <- ggplot(data = df1, aes(x = x, y = y, z = z)) +
  geom_contour_filled(alpha=0.5, bins = n_bins)+
  geom_contour(color = "black", size = 0.1, bins = n_bins)+
  geom_point(data = x_y_ric, aes(x = x, y = y, size=z)) +
  scale_fill_brewer("Sobs_dens", palette = 'YlOrRd')+
  scale_size_continuous("Sobs", breaks = c(1,3,7,16))+
  xlab("Na") +
  ylab("P")+
  theme_classic(base_size = 14)+
  theme(legend.position = "right",
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
p4

library(visreg)
fit <- lm(Sobs ~ Na_Cl * P, data=df)
val_adj <- 1.0000001
val_adj <- 1.0051
bd2 <- data.frame(x1=c(0.05,2.2,2.735,2.71,2.71,2.71),
                  y1=c(1.2,1.2,3.6,7,10.2,13.6),
                  riq=c(0,4,8,12,16,20)) 
#https://kwstat.github.io/pals/reference/brewer.html
#saimos do default de palette que era do tipo Diverging ("bwr") para uma sequencial ("YlOrRd") https://matplotlib.org/stable/tutorials/colors/colormaps.html
p4.1 <- visreg2d(fit, x="Na_Cl", y="P", plot.type="gg") +
  geom_contour(aes(z=z), color="gray40", binwidth = 4.1, alpha=.6)+
  scale_fill_distiller("Sobs", palette= "YlOrRd", direction = 1, breaks = c(0,4,8,12,16,20)) +
  geom_point(data=df, aes(x=Na_Cl, y=P, shape=Area), alpha=0.8, size = 3)+
  geom_text(data=bd2, aes(x=x1, y=y1, label = riq), colour="gray40", size = 6) + #Species name
  labs(x= "Na", y= "P", tag = "B")+
  lims(x=c(-((max(df$Na_Cl)*val_adj)-max(df$Na_Cl)), (max(df$Na_Cl)*val_adj)),
       y=c(-((max(df$P)*val_adj)-max(df$P)),(max(df$P)*val_adj)))+
  theme_classic(base_size = 14)+
  theme(legend.position = "right",
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line (size=.4),
        axis.line.x.bottom = element_line(color = "black", size = .6),
        axis.text.x.bottom = element_text(color = "black"),
        axis.line.y.left = element_line(color = "black", size = .6),
        axis.text.y.left = element_text(color = "black"))+
  guides(shape="none")

p4.1

library(gridExtra)
library(cowplot)
plot_grid(arrangeGrob(p1+
                        theme(legend.position = "none"),
                      p4.1,
                      ncol=2, widths = c(1,3)))

#ggsave("Fig_3.1.png", width = 11, height = 6, dpi = 400)
#ggsave("Fig_3.1.pdf", width = 11, height = 6)
#ggsave("Fig_3.1.tiff", width = 11, height = 6, dpi = 400, compression = 'lzw')


#### Figure 3: RDA

# Morphological composition

MyVar<- c("WL", "HL.std", "HW.std", "SL.std", "ML.std", "EL.std", "IO.std", "DEM.std", "FL.std", "TL.std")
morf <- vegan::decostand(df[,MyVar], method="hellinger")
var <- vegan::decostand(df[,c("Na_Cl", "Al", "pH", "dens_part", "P", "K", "MO", "N")], method = "standardize")

RDA1 <- vegan::rda(morf, var)

variaveis <- t(vegan::intersetcor(RDA1)[,1:3])
compo1 <- as.data.frame(summary(RDA1)$cont)
eigenv1 <- t(as.vector(compo1[1,1:3]))
cont1 <- t(as.vector(compo1[2,1:3]*100))
acum1 <- t(as.vector(compo1[3,1:3]*100))
teste <- as.data.frame(anova(RDA1))
p_valor <- as.data.frame(c(teste[1,4],NA,NA))

(tabela1 <- t(cbind(variaveis, eigenv1, cont1, acum1, p_valor)))
#write.table(tabela,"RDA_formigas_tabela.txt") #export results

plot(RDA1, choices = c(1, 2), display = c("sp", "lc", "cn"))

##Figure ggplot2
value1 <- summary(RDA1)

#Units:
rda_units <- data.frame(value1$sites)[,1:2]
rda_units$area <- df$Area

#Species:
rda_ssp <- data.frame(value1$species)[,1:2]
rda_ssp$species_names <-gsub(".std","",rownames(rda_ssp))

#Arrows:
rda_arrows <- data.frame(value1$biplot)[,1:2]
paste(rownames(rda_arrows), collapse = "','")
rda_arrows$vari <- c('Na','Al','pH','Dp','P','K','OM','N')

#Adjust for arrows
mult <- 0.4 #Fator de ajuste do x2 e y2. Mudar at? ter uma rela??o boa entre os dois

cor_ssp <- "gray50" #Color species
cor_name_seta <- "red" #Color names variables


library(ggrepel)
rda1 <- ggplot(rda_units, aes(x = RDA1, y = RDA2))+
  geom_hline(yintercept = 0, colour="gray", linetype="dashed")+ geom_vline(xintercept = 0, colour="gray", linetype="dashed")+
  geom_segment(data = rda_arrows, aes(x = 0, xend = RDA1*mult, y = 0, yend = RDA2*mult), arrow = arrow(length = unit(0.25, "cm")), colour = "red", size = .8, alpha=0.2)+
  geom_point(aes(shape=area), size = 3, alpha=0.5)+
  geom_point(data = rda_ssp, aes(x = RDA1, y = RDA2), size = 1.5, shape=3, colour = 1, alpha=.6)+ #Species point
  geom_text_repel(data = rda_ssp, aes(x = RDA1, y = RDA2, label = species_names), segment.color = "gray", colour=1, size = 3) + #Species name
  geom_text_repel(data = rda_arrows, aes(x = RDA1*mult, y = RDA2*mult, label = vari), direction = "x", segment.color = 1, colour="red", size = 5) +
  labs(x = "RDA 1 (13.22%)", y = "RDA 2 (9.84%)", tag = "B")+
  scale_x_continuous(sec.axis = sec_axis(~./mult))+
  scale_y_continuous(sec.axis = sec_axis(~./mult))+
  theme_classic(base_size = 12)+
  theme(legend.position = c(0.12, 0.92),
        legend.background = element_rect(color = "black", size = .5, linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(),
        axis.line.x.bottom = element_line(color = "black", size = .8),
        axis.text.x.bottom = element_text(color = "black"),
        axis.line.y.left = element_line(color = "black", size = .8),
        axis.text.y.left = element_text(color = "black"),
        axis.line.x.top = element_line(color = cor_name_seta, size = .8),
        axis.text.x.top = element_text(color = cor_name_seta),
        axis.ticks.x.top = element_line(color = cor_name_seta),
        axis.line.y.right = element_line(color = cor_name_seta, size = .8),
        axis.text.y.right = element_text(color = cor_name_seta),
        axis.ticks.y.right = element_line(color = cor_name_seta))
rda1


# community composition

ssp1 <- vegan::decostand(my.sample, method="hellinger")

RDA2 <- vegan::rda(ssp1, var)

variaveis <- t(vegan::intersetcor(RDA2)[,1:3])
compo1 <- as.data.frame(summary(RDA2)$cont)
eigenv1 <- t(as.vector(compo1[1,1:3]))
cont1 <- t(as.vector(compo1[2,1:3]*100))
acum1 <- t(as.vector(compo1[3,1:3]*100))
teste <- as.data.frame(anova(RDA2))
p_valor <- as.data.frame(c(teste[1,4],NA,NA))

(tabela2 <- t(cbind(variaveis, eigenv1, cont1, acum1, p_valor)))
#write.table(rbind(tabela1, tabela2),"res_RDA.txt") # export results

plot(RDA2, choices = c(1, 2), display = c("sp", "lc", "cn"))

##Figure ggplot2
value2 <- summary(RDA2)

#Units:
rda_units1 <- data.frame(value2$sites)[,1:2]
rda_units1$area <- df$Area

#Species:
rda_ssp1 <- data.frame(value2$species)[,1:2]
#rda_ssp1$species_names <-gsub(".std","",rownames(rda_ssp))
row.names(rda_ssp1)
rda_ssp1$species_names <- c("Br.M01", "Ca.bre", "Ca.glo", "Ca.inc", "Ca.M04", "Ca.M05", "Cr.obs", "Cr.sot", "Cy.ham", "Cy.lae", "Cy.pel", "Di.sex", "Do.imi", "Fu.elo", "Gn.ple", "Hy.M04", "Hy.M05", "Hy.M06", "Hy.M07", "Mo.flo", "My.for", "Ne.spl", "Ny.ful", "Oc.neo", "Oc.bal", "Oc.ihe", "Pa.har", "Ph.fra", "Ph.D01", "Ph.F01", "Ph.F02", "Ph.F03", "Ph.sco", "Ph.M26", "Ra.arh", "Ra.M01", "Rh.M01", "Ro.M01", "So.M20", "So.M15", "So.M19", "So.M21", "So.M22", "St.car", "Ty.mei", "Ty.pus", "Ty.rog")


#Arrows:
rda_arrows1 <- data.frame(value2$biplot)[,1:2]
paste(rownames(rda_arrows), collapse = "','")
rda_arrows1$vari <- c('Na','Al','pH','Dp','P','K','OM','N')

#Adjust for arrows
mult1 <- 1.1 #Fator de ajuste do x2 e y2. Mudar at? ter uma rela??o boa entre os dois
mult2 <- 1.5


rda2 <- ggplot(rda_units1, aes(x = RDA1, y = RDA2))+
  geom_hline(yintercept = 0, colour="gray", linetype="dashed")+ geom_vline(xintercept = 0, colour="gray", linetype="dashed")+
  geom_segment(data = rda_arrows1[-3,], aes(x = 0, xend = RDA1*mult1, y = 0, yend = RDA2*mult1), arrow = arrow(length = unit(0.25, "cm")), colour = "red", size = .8, alpha=0.2)+ #N sig
  geom_segment(data = rda_arrows1[3,], aes(x = 0, xend = RDA1*mult1, y = 0, yend = RDA2*mult1), arrow = arrow(length = unit(0.25, "cm")), colour = "red", size = .8)+ #Sig
  geom_point(aes(shape=area), size = 3, alpha=0.5)+
  #geom_point(data = rda_ssp1, aes(x = RDA1*mult1, y = RDA2*mult1), size = 1.5, shape=16, colour = cor_name_seta, alpha=.4)+ #Species point
  geom_point(data = rda_ssp1, aes(x = RDA1*mult2, y = RDA2*mult2), size = 1.5, shape=3, colour = cor_ssp, alpha=.6)+ #Species point
  geom_text_repel(data = rda_ssp1, aes(x = RDA1*mult2, y = RDA2*mult2, label = species_names), segment.color = "gray", colour=1, size = 1.8, fontface="italic") +
  geom_text_repel(data = rda_arrows1, aes(x = RDA1*mult1, y = RDA2*mult1, label = vari), direction = "x", segment.color = 1, colour="red", size = 5) +
  labs(x = "RDA 1 (7.65%)", y = "RDA 2 (6.15%)")+
  scale_x_continuous(sec.axis = sec_axis(~./mult1))+
  scale_y_continuous(sec.axis = sec_axis(~./mult1))+
  theme_classic(base_size = 12)+
  theme(legend.position = c(0.12, 0.92),
        legend.background = element_rect(color = "black", size = .5, linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(),
        axis.line.x.bottom = element_line(color = "black", size = .8),
        axis.text.x.bottom = element_text(color = "black"),
        axis.line.y.left = element_line(color = "black", size = .8),
        axis.text.y.left = element_text(color = "black"),
        axis.line.x.top = element_line(color = cor_name_seta, size = .8),
        axis.text.x.top = element_text(color = cor_name_seta),
        axis.ticks.x.top = element_line(color = cor_name_seta),
        axis.line.y.right = element_line(color = cor_name_seta, size = .8),
        axis.text.y.right = element_text(color = cor_name_seta),
        axis.ticks.y.right = element_line(color = cor_name_seta))
rda2

mult2 <- 4
rda2 <- ggplot(rda_units1, aes(x = RDA1, y = RDA2))+
  geom_hline(yintercept = 0, colour="gray", linetype="dashed")+ geom_vline(xintercept = 0, colour="gray", linetype="dashed")+
  geom_segment(data = rda_arrows1[-3,], aes(x = 0, xend = RDA1, y = 0, yend = RDA2), arrow = arrow(length = unit(0.25, "cm")), colour = 1, size = .8, alpha=0.2)+ #N sig
  geom_segment(data = rda_arrows1[3,], aes(x = 0, xend = RDA1, y = 0, yend = RDA2), arrow = arrow(length = unit(0.25, "cm")), colour = 1, size = .8)+ #Sig
  geom_point(aes(shape=area), size = 3, alpha=0.5)+
  #geom_point(data = rda_ssp1, aes(x = RDA1*mult1, y = RDA2*mult1), size = 1.5, shape=16, colour = cor_name_seta, alpha=.4)+ #Species point
  geom_point(data = rda_ssp1, aes(x = RDA1*mult2, y = RDA2*mult2), size = 1, shape=3, colour = "red", alpha=.6)+ #Species point
  geom_text_repel(data = rda_ssp1, aes(x = RDA1*mult2, y = RDA2*mult2, label = species_names), segment.color = "red", min.segment.length = 10, colour="red", size = 2, fontface="italic") +
  geom_text_repel(data = rda_arrows1, aes(x = RDA1, y = RDA2, label = vari), direction = "x", segment.color = 1, colour=1, size = 5) +
  labs(x = "RDA 1 (7.65%)", y = "RDA 2 (6.15%)", tag = "A")+
  scale_x_continuous(sec.axis = sec_axis(~./mult2))+
  scale_y_continuous(sec.axis = sec_axis(~./mult2))+
  theme_classic(base_size = 12)+
  theme(legend.position = c(0.12, 0.92),
        legend.background = element_rect(color = "black", size = .5, linetype = "solid"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(),
        axis.line.x.bottom = element_line(color = "black", size = .8),
        axis.text.x.bottom = element_text(color = "black"),
        axis.line.y.left = element_line(color = "black", size = .8),
        axis.text.y.left = element_text(color = "black"),
        axis.line.x.top = element_line(color = cor_name_seta, size = .8),
        axis.text.x.top = element_text(color = cor_name_seta),
        axis.ticks.x.top = element_line(color = cor_name_seta),
        axis.line.y.right = element_line(color = cor_name_seta, size = .8),
        axis.text.y.right = element_text(color = cor_name_seta),
        axis.ticks.y.right = element_line(color = cor_name_seta))
rda2



library(gridExtra)
library(cowplot)
plot_grid(arrangeGrob(rda2+
                        theme(legend.position = "none"),
                      rda1+
                        theme(legend.position = "none"),
                      ncol=1))

#ggsave("Fig_4.1.pdf", width = 8, height = 12)
#ggsave("Fig_4.1.png", width = 8, height = 12, dpi = 500)
#ggsave("Fig_4.1.tiff", width = 8, height = 12, dpi = 500, compression = 'lzw')


