#
dir ()

#### Statiscal analyses

# Datasets for analyses

traits.std <- read.table("dataset_standardized_ant_traits.txt", header = TRUE, sep="\t")
my.sample <- read.table("dataset_ant_community.txt", header = TRUE, sep="\t")
df <- read.table("dataset_variables_IN_soil.txt", header = TRUE, sep="\t") # contains also standardized soil variables and ant traits used in models (i.e., PC1-PC2 axes)
head(df)


#### soils properties between the two transects ####

kruskal.test(df$Sobs ~df$Area) # yes for richness

kruskal.test(df$Areia ~df$Area) # yes for granulometric
kruskal.test(df$Silte ~df$Area)
kruskal.test(df$Argila ~df$Area)

kruskal.test(df$pH ~df$Area) 
kruskal.test(df$C ~df$Area)
kruskal.test(df$MO ~df$Area)
kruskal.test(df$N ~df$Area) # yes for nitrogen
kruskal.test(df$C_N ~df$Area) # yes for carbon:nitrogen ratio
kruskal.test(df$Al ~df$Area) # yes for aluminium
kruskal.test(df$Na_Cl ~df$Area) # yes for sodium
kruskal.test(df$K ~df$Area) 
kruskal.test(df$P ~df$Area) # yes for phosphorus
kruskal.test(df$P ~df$dens_part) 

#### Relationships between ant species richness, taxonomic composition and soil properties

library(MASS)
library(fmsb)
MyVar<- c("pH", "C", "MO", "N", "C_N", "Al", "Na_Cl", "K", "P", "dens_part")
source("vif_fun.r")

vif_func(in_frame=df[,MyVar], thresh=5, trace=T) 
# selected for statistical models: pH + N + MO + Al + Na_Cl + K + P + dens_part
# removed: C and C_N, 



### Statistical models using linear mixed models (LMMs) ====
library(nlme)
m1 <- lme(Sobs ~ pH.std + N.std + MO.std + Al.std + Na_Cl.std + K.std + P.std + dens_part.std,
          random=~1|Area, method = "REML", data=df)

library(MuMIn)
selec <- dredge(m1, rank="AICc", extra = "adjR^2")
(mod_acaic <- summary(model.avg(selec))) # AIC result
sw(selec)
r.squaredGLMM(m1)


#### Relationships between ant morphology and soil properties

m2 <- lme(PC1 ~ pH.std + MO.std + N.std + Al.std + Na_Cl.std + K.std + P.std + dens_part.std, random=~1|Area, method = "REML", data=df)

selec <- dredge(m2, rank="AICc", extra = "adjR^2")
(mod_acaic <- summary(model.avg(selec))) # PC1 - morphology
sw(selec)
r.squaredGLMM(m2)

name_var <- "PC1"
tab_mod_temp <- mod_func(mod_ac=mod_acaic, name_v=name_var)
tab_mod_fin <- data.frame(rbind(tab_mod_fin, tab_mod_temp))
tab_impor_temp <- imp_func(selecao=selec, name_v=name_var)
tab_impor_fin <- data.frame(rbind(tab_impor_fin, tab_impor_temp))


m3 <- lme(PC2 ~ pH.std + MO.std + N.std + Al.std + Na_Cl.std + K.std + P.std + dens_part.std, random=~1|Area, method = "REML", data=df)

selec <- dredge(m3, rank="AICc", extra = "adjR^2")
(mod_acaic <- summary(model.avg(selec))) # PC2 - morphology
sw(selec)
r.squaredGLMM(m3)

name_var <- "PC2"
tab_mod_temp <- mod_func(mod_ac=mod_acaic, name_v=name_var)
tab_mod_fin <- data.frame(rbind(tab_mod_fin, tab_mod_temp))
tab_impor_temp <- imp_func(selecao=selec, name_v=name_var)
tab_impor_fin <- data.frame(rbind(tab_impor_fin, tab_impor_temp))


#write.table(tab_impor_fin, "tabela.txt",sep="\t", row.names = FALSE) # export results; change r object for each model


#### Morphological composition ====
MyVar<- c("WL", "HL.std", "HW.std", "SL.std", "ML.std", "EL.std", "IO.std", "DEM.std", "FL.std", "TL.std")

#Ordered by R2
library(vegan)
adonis2(df[,MyVar] ~ df$Al.std, permutations=9999, method = "euclidian", strata = df$Area) #0.12229
adonis2(df[,MyVar] ~ df$Na_Cl.std, permutations=9999, method = "euclidian", strata = df$Area) #0.07787
adonis2(df[,MyVar] ~ df$dens_part.std, permutations=9999, method = "euclidian", strata = df$Area) #0.04106
adonis2(df[,MyVar] ~ df$MO.std, permutations=9999, method = "euclidian", strata = df$Area) #0.03669
adonis2(df[,MyVar] ~ df$P.std, permutations=9999, method = "euclidian", strata = df$Area) #0.03602
adonis2(df[,MyVar] ~ df$pH.std, permutations=9999, method = "euclidian", strata = df$Area) #0.01852
adonis2(df[,MyVar] ~ df$N.std, permutations=9999, method = "euclidian", strata = df$Area) #0.01575
adonis2(df[,MyVar] ~ df$K.std, permutations=9999, method = "euclidian", strata = df$Area) #0.01491


result1 <- adonis2(df[,MyVar] ~ (df$Al.std + df$Na_Cl.std + df$dens_part.std + df$MO.std + df$P.std + df$pH.std + df$N.std + df$K.std), permutations=9999, method = "euclidian", strata = df$Area)
result1
round(data.frame(result1$aov.tab),2)


#### Taxonomic Composition ====
#Ordered by R2
adonis2(my.sample ~ df$Na_Cl.std, permutations=9999, method = "jaccard", strata = df$Area) #0.06776
adonis2(my.sample ~ df$Al.std, permutations=9999, method = "jaccard", strata = df$Area) #0.0631
adonis2(my.sample ~ df$pH.std, permutations=9999, method = "jaccard", strata = df$Area) #0.05602
adonis2(my.sample ~ df$dens_part.std, permutations=9999, method = "jaccard", strata = df$Area) #0.05085
adonis2(my.sample ~ df$P.std, permutations=9999, method = "jaccard", strata = df$Area) #0.03785
adonis2(my.sample ~ df$K.std, permutations=9999, method = "jaccard", strata = df$Area) #0.03248
adonis2(my.sample ~ df$MO.std, permutations=9999, method = "jaccard", strata = df$Area) #0.02869
adonis2(my.sample ~ df$N.std, permutations=9999, method = "jaccard", strata = df$Area) #0.02211


result2 <- adonis2(my.sample ~ (df$Na_Cl.std + df$Al.std + df$pH.std + df$dens_part.std + df$P.std + df$K.std + df$MO.std + df$N.std),
                 permutations=9999, method = "jaccard", strata = df$Area) #Import?ncia biol?gica e feito em blocos
result2

res_adonis <- data.frame(var=c(rownames(result2$aov.tab),rownames(result1$aov.tab)), R2ord=c(0.06776,0.0631,0.05602,0.05085,0.03785,0.03248,0.02869,0.02211,NA,NA,
                                                                                             0.12229, 0.07787, 0.04106, 0.03669, 0.03283, 0.02632, 0.01943,0.01535,NA,NA), rbind(result2$aov.tab, result1$aov.tab))

#write.table(res_adonis, "res_model_PERMANOVA.txt",sep="\t", row.names = FALSE) # export results
rownames(result2$aov.tab)


#### modelos para cada trait ####
#WL
m <- lme(WL ~ pH.std + N.std + MO.std + Al.std + Na_Cl.std + K.std + P.std + dens_part.std, random=~1|Area, method = "REML", data=df)

selec <- dredge(m, rank="AICc", extra = "adjR^2")
(mod_acaic <- summary(model.avg(selec)))
sw(selec)
r.squaredGLMM(m)

name_var <- "WL"
tab_mod_temp <- mod_func(mod_ac=mod_acaic, name_v=name_var)
tab_mod_fin <- data.frame(rbind(tab_mod_fin, tab_mod_temp))
tab_impor_temp <- imp_func(selecao=selec, name_v=name_var)
tab_impor_fin <- data.frame(rbind(tab_impor_fin, tab_impor_temp))


#HL.std
m <- lme(HL.std ~ pH.std + N.std + MO.std + Al.std + Na_Cl.std + K.std + P.std + dens_part.std,
         random=~1|Area, method = "REML", data=df)
r.squaredGLMM(m)
selec <- dredge(m, rank="AICc", extra = "adjR^2")
(mod_acaic <- summary(model.avg(selec)))
sw(selec)

name_var <- "HL.std"
tab_mod_temp <- mod_func(mod_ac=mod_acaic, name_v=name_var)
tab_mod_fin <- data.frame(rbind(tab_mod_fin, tab_mod_temp))
tab_impor_temp <- imp_func(selecao=selec, name_v=name_var)
tab_impor_fin <- data.frame(rbind(tab_impor_fin, tab_impor_temp))


#ML.std
m <- lme(ML.std ~ pH.std + N.std + MO.std + Al.std + Na_Cl.std + K.std + P.std + dens_part.std,
         random=~1|Area, method = "REML", data=df)
r.squaredGLMM(m)
selec <- dredge(m, rank="AICc", extra = "adjR^2")
(mod_acaic <- summary(model.avg(selec)))
importance(selec)

name_var <- "ML.std"
tab_mod_temp <- mod_func(mod_ac=mod_acaic, name_v=name_var)
tab_mod_fin <- data.frame(rbind(tab_mod_fin, tab_mod_temp))
tab_impor_temp <- imp_func(selecao=selec, name_v=name_var)
tab_impor_fin <- data.frame(rbind(tab_impor_fin, tab_impor_temp))


#SL.std
m <- lme(SL.std ~ pH.std + N.std + MO.std + Al.std + Na_Cl.std + K.std + P.std + dens_part.std,
         random=~1|Area, method = "REML", data=df)
r.squaredGLMM(m)
selec <- dredge(m, rank="AICc", extra = "adjR^2")
(mod_acaic <- summary(model.avg(selec)))
sw(selec)

name_var <- "SL.std"
tab_mod_temp <- mod_func(mod_ac=mod_acaic, name_v=name_var)
tab_mod_fin <- data.frame(rbind(tab_mod_fin, tab_mod_temp))
tab_impor_temp <- imp_func(selecao=selec, name_v=name_var)
tab_impor_fin <- data.frame(rbind(tab_impor_fin, tab_impor_temp))


#EL.std
m <- lme(EL.std ~ pH.std + N.std + MO.std + Al.std + Na_Cl.std + K.std + P.std + dens_part.std,
         random=~1|Area, method = "REML", data=df)
r.squaredGLMM(m)
selec <- dredge(m, rank="AICc", extra = "adjR^2")
(mod_acaic <- summary(model.avg(selec)))
sw(selec)

name_var <- "EL.std"
tab_mod_temp <- mod_func(mod_ac=mod_acaic, name_v=name_var)
tab_mod_fin <- data.frame(rbind(tab_mod_fin, tab_mod_temp))
tab_impor_temp <- imp_func(selecao=selec, name_v=name_var)
tab_impor_fin <- data.frame(rbind(tab_impor_fin, tab_impor_temp))


#FL.std
m <- lme(FL.std ~ pH.std + N.std + MO.std + Al.std + Na_Cl.std + K.std + P.std + dens_part.std,
         random=~1|Area, method = "REML", data=df)
r.squaredGLMM(m)
selec <- dredge(m, rank="AICc", extra = "adjR^2")
(mod_acaic <- summary(model.avg(selec)))
sw(selec)

name_var <- "FL.std"
tab_mod_temp <- mod_func(mod_ac=mod_acaic, name_v=name_var)
tab_mod_fin <- data.frame(rbind(tab_mod_fin, tab_mod_temp))
tab_impor_temp <- imp_func(selecao=selec, name_v=name_var)
tab_impor_fin <- data.frame(rbind(tab_impor_fin, tab_impor_temp))


#TL.std
m <- lme(TL.std ~ pH.std + N.std + MO.std + Al.std + Na_Cl.std + K.std + P.std + dens_part.std,
         random=~1|Area, method = "REML", data=df)
r.squaredGLMM(m)
selec <- dredge(m, rank="AICc", extra = "adjR^2")
(mod_acaic <- summary(model.avg(selec)))
sw(selec)

name_var <- "TL.std"
tab_mod_temp <- mod_func(mod_ac=mod_acaic, name_v=name_var)
tab_mod_fin <- data.frame(rbind(tab_mod_fin, tab_mod_temp))
tab_impor_temp <- imp_func(selecao=selec, name_v=name_var)
tab_impor_fin <- data.frame(rbind(tab_impor_fin, tab_impor_temp))


#IO.std
m <- lme(IO.std ~ pH.std + N.std + MO.std + Al.std + Na_Cl.std + K.std + P.std + dens_part.std,
         random=~1|Area, method = "REML", data=df)
r.squaredGLMM(m)
selec <- dredge(m, rank="AICc", extra = "adjR^2")
(mod_acaic <- summary(model.avg(selec)))
sw(selec)

name_var <- "IO.std"
tab_mod_temp <- mod_func(mod_ac=mod_acaic, name_v=name_var)
tab_mod_fin <- data.frame(rbind(tab_mod_fin, tab_mod_temp))
tab_impor_temp <- imp_func(selecao=selec, name_v=name_var)
tab_impor_fin <- data.frame(rbind(tab_impor_fin, tab_impor_temp))


#HW.std
m <- lme(HW.std ~ pH.std + N.std + MO.std + Al.std + Na_Cl.std + K.std + P.std + dens_part.std,
         random=~1|Area, method = "REML", data=df)
r.squaredGLMM(m)
selec <- dredge(m, rank="AICc", extra = "adjR^2")
(mod_acaic <- summary(model.avg(selec)))
sw(selec)

name_var <- "HW.std"
tab_mod_temp <- mod_func(mod_ac=mod_acaic, name_v=name_var)
tab_mod_fin <- data.frame(rbind(tab_mod_fin, tab_mod_temp))
tab_impor_temp <- imp_func(selecao=selec, name_v=name_var)
tab_impor_fin <- data.frame(rbind(tab_impor_fin, tab_impor_temp))


#DEM.std
m <- lme(DEM.std ~ pH.std + N.std + MO.std + Al.std + Na_Cl.std + K.std + P.std + dens_part.std,
         random=~1|Area, method = "REML", data=df)
r.squaredGLMM(m)
selec <- dredge(m, rank="AICc", extra = "adjR^2")
(mod_acaic <- summary(model.avg(selec)))
sw(selec)

name_var <- "DEM.std"
tab_mod_temp <- mod_func(mod_ac=mod_acaic, name_v=name_var)
tab_mod_fin <- data.frame(rbind(tab_mod_fin, tab_mod_temp))
tab_impor_temp <- imp_func(selecao=selec, name_v=name_var)
tab_impor_fin <- data.frame(rbind(tab_impor_fin, tab_impor_temp))


#write.table(tab_mod_fin, "res_models_selection.txt",sep="\t", row.names = FALSE) # para pegar resultados
#write.table(tab_impor_fin, "res_models_importance.txt",sep="\t", row.names = FALSE) # para pegar resultados






