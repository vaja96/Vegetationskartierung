#########
# Vegetationskartierung der Drachenwiese, vor Ort und mit multispektralen Orthophotos
# Wie gut lassen sich folgende Vegetationsmerkmale der Drachenwiese anhand raeumlich hoch
# aufgeloester Multispektralbilder kartieren?
# Pflanzengesellschaften als floristische Gradienten und in Klassen
#########

########
# Floristische Gradienten: Ordination + Random Forest Regression
# Ergaenzung zu den vorgegebenen Skripten: 
# Datenvorbereitung, Klassifikation, Regression 
# von Michael Ewald fuer den Kurs Vegetationskartierung und Vegetationsaufnahme fuer das Sommersemester 2026 am KIT

# Floristische Gradienten werden mittels Ordination (NMDS) ergaenzt, genauso wie in Regression mit EIVE-N
# Dann als kontinuierliche Zielvariable per Random Forest modelliert und in die Flaeche uebertragen
# Extraktion der Spektralwerte, isopam + RF/SVM für diskrete Klassen, RF-Regression für EIVE-N als Beispiel
#######

#######
# Autorin: Vanessa Jasmin Herr
#######
# Start des Codes
# libraries laden


library(vegan)         # Ordination (NMDS)
library(dplyr)          # Datenoperationen
library(terra)           # Rasterdatenverarbeitung
library(randomForest)    # Random Forest Modelle
library(VSURF)            # Variablenselektion
library(ggplot2)           # Abbildungen
library(pdp)                # Partial Dependence Plots
library(tidyr)             # Visualisierung

# Artentabelle einlesen wie in Klassifikation
sp <- read.csv2("Arten_final26.csv", row.names = 1)
sp <- sp[, -which(colnames(sp) == "p19")]   # p19: keine Koordinaten
sp <- t(sp)                                  # Plots in Zeilen, Arten in Spalten
sp <- as.data.frame(sp)

# Da der Plot 34 nur eine Art besitzt und diese die Bray-Curtis Distanz zu stark verzerrt, wird diese Art ausgeschlossen (siehe Klassifikation und Regression)
sp <- sp[rownames(sp) != "p34", ]


#######
# NMDS Ordination
# reduziert hochdimensionale Artzusammensetzung auf wenige floristische Gradienten (Achsen) die die Unterschiede zwischen den Aufnahmen moeglichst gut abbilden
#######

set.seed(84) # Zahl gewaehlt wie im Klassifikationsskript
nmds <- metaMDS(sp, distance = "bray", k = 2, trymax = 100)
nmds

# Stress-Wert pruefen (Faustregel: < 0.2 akzeptabel, < 0.1 gut)
nmds$stress
stressplot(nmds) # 0.22 > 0.2

# knapp ueber der Grenze von 0.2 - die floristischen Gradienten sind also mit einer gewissen Unsicherheit behaftet
# Probieren mit k = 3

set.seed(84) # Zahl gewaehlt wie im Klassifikationsskript
nmds <- metaMDS(sp, distance = "bray", k = 3, trymax = 100)
nmds

# Stress-Wert pruefen (Faustregel: < 0.2 akzeptabel, < 0.1 gut)
nmds$stress
stressplot(nmds) # 0.16 < 0.2

#######
# Stress-Scree-Plot: Stress in Abhaengigkeit von der Dimensionszahl k 
# Faustregel: wo die Kurve deutlich abflacht (Ellbogen) bringt eine weitere Dimension nur noch marginale Verbesserung
#######

set.seed(84)
k_werte <- 1:6
stress_werte <- numeric(length(k_werte))
 
for (i in seq_along(k_werte)) {
  m <- metaMDS(sp, distance = "bray", k = k_werte[i],
               trymax = 100, trace = FALSE)
  stress_werte[i] <- m$stress
}
 
plot(k_werte, stress_werte, type = "b", pch = 19,
     xlab = "Anzahl Dimensionen (k)", ylab = "Stress",
     main = "")
abline(h = 0.2, col = "red", lty = 2)   # kritische Schwelle einzeichnen
text(k_werte, stress_werte, labels = round(stress_werte, 3), pos = 3)

# Varianzanteil je Achse in der 3D Loesung
# metaMDS führt eine PC Rotation durch: die Site Scores sind beteits nach abnehmender Varianz sortiert (NMDS1 > NMDS2 > NMDS3)
# Quantifikation direkt über die Varianz der Achsen-Scores

set.seed(84)
nmds3 <- metaMDS(sp, distance = "bray", k = 3, trymax = 100)
 
sc3 <- scores(nmds3, display = "sites")
var_je_achse <- apply(sc3, 2, var)
anteil_je_achse <- round(var_je_achse / sum(var_je_achse) * 100, 1)
 
anteil_je_achse


# Procrustes Vergleich, ob die dritte Achse verzerrt
set.seed(84)
nmds2 <- metaMDS(sp, distance = "bray", k = 2, trymax = 100)
 
proc <- procrustes(nmds2, nmds3, choices = c(1, 2))
summary(proc)
plot(proc, main = "Procrustes: k=2 vs. erste 2 Achsen von k=3")
 
# statistischer Test: ist die Uebereinstimmung besser als Zufall?
set.seed(84)
prot <- protest(nmds2, nmds3, choices = c(1, 2))
prot

#######
# Es wir mit k =3 weiter gemacht
# k = 4 wird nicht in Betracht gezogen, weil die Stichprobengroesse zu klein wird um stabile Loesungen zu finden
# zusaetzlich wird der Aufwand immer groesser und die die oekologische Interpretierbarkeit laesst nach
#######

sp <- read.csv2("Arten_final26.csv", row.names = 1)
sp <- sp[, -which(colnames(sp) == "p19")]   # p19: keine Koordinaten
sp <- t(sp)                                  # Plots in Zeilen, Arten in Spalten
sp <- as.data.frame(sp)
 
# Plot mit nur einer Art (p34) verzerrt die Bray-Curtis-Distanz stark
# und wurde auch in Klassifikation / Regression ausgeschlossen
sp <- sp[rownames(sp) != "p34", ]
 

# NMDS Ordination mit k = 3
set.seed(84)
nmds <- metaMDS(sp, distance = "bray", k = 3, trymax = 100)
nmds
 
# Stress-Wert pruefen (Faustregel: < 0.2 akzeptabel, < 0.1 gut)
nmds$stress
stressplot(nmds)

# Artenscores nachtraeglich ergaenzen (falls "Species: scores missing")
sp_bereinigt <- sp[, colSums(sp) > 0]
sp_scores <- wascores(scores(nmds, display = "sites"), sp_bereinigt)
nmds$species <- sp_scores

# scores anschauen für die oekologische Interpretation
# NMDS1
species_scores <- as.data.frame(scores(nmds, display = "species"))

species_scores[order(species_scores$NMDS1), ]

#NMDS2
species_scores[order(species_scores$NMDS2), ]

#NMDS3
species_scores[order(species_scores$NMDS3), ]

# Ordinationsdiagramme zur Kontrolle (paarweise Achsenkombinationen)
plot(nmds, choices = c(1, 2), type = "t", main = "NMDS Achse 1 vs. 2")
plot(nmds, choices = c(1, 3), type = "t", main = "NMDS Achse 1 vs. 3")
plot(nmds, choices = c(2, 3), type = "t", main = "NMDS Achse 2 vs. 3")
 

# Achsenwerte (Scores) je Plot extrahieren
scores_df <- as.data.frame(scores(nmds, display = "sites"))
scores_df$plot <- rownames(scores_df)

####### 
# Mit Spektraldaten verknuepfen wie bei Regression (eive, uavdata)
#######
uavdata <- read.csv2("extracted_values_Drachenwiese26.csv")
 
scores_df <- arrange(scores_df, plot)
uavdata   <- arrange(uavdata, plot)

# Nur Plots behalten, die in beiden Datensaetzen vorhanden sind (Plot 19 sollte rausfallen, weil für den keine GPS Koordinaten da sind)
common_plots <- intersect(scores_df$plot, uavdata$plot)
scores_df <- scores_df[scores_df$plot %in% common_plots, ]
uavdata   <- uavdata[uavdata$plot %in% common_plots, ]
 
identical(scores_df$plot, uavdata$plot)   # sollte TRUE sein
## --> insgesamt also nur 46 Plotaufnahmen

#########
# Umweltvariablen (Spektralwerte) auf die Ordination projizieren
# Ergaenzend zu den VSURF-Ergebnissen (diese betrachten nur eine Achse) um eine Gesamtschau aller 12 Variablen zeitgleich

# Reihenfolge der Plots exakt an die NMDS Site Scores angleichen
uav_env <- uavdata[match(rownames(scores(nmds, display = "sites")), uavdata$plot), 1:12]

fit <- envfit(nmds, uav_env, permutations = 999)
fit   # r2 und Pr(>r) je Variable pruefen - nur signifikante Vektoren sind interpretierbar

plot(nmds, display = "sites", type = "t")
plot(fit, col = "red", p.max = 0.05)   # nur signifikante Vektoren mit p < 0.05 anzeigen

### Plot der fuer die Auswertung genommen wird
plot(nmds, display = "sites", type = "p")   # nur Punkte statt Plotnamen
plot(fit, col = "red", p.max = 0.05)

# 10 Arten mit hoechsten und niedrigsten Werten jeder Achse ausgeben
species_scores <- as.data.frame(scores(nmds, display="species"))

species_scores$Art <- rownames(species_scores)

species_scores %>%
  arrange(NMDS1) %>%
  head(10)

species_scores %>%
  arrange(desc(NMDS1)) %>%
  head(10)


#######
# Trainings- und Testdatensatz analog zu Regression
# Ein split, welcher gemeinsam für alle drei NMDS Achsen verwendet wird
#######

set.seed(84) # wie bei Regression
n <- nrow(uavdata)
test_id <- sample(seq_len(n), size = round(0.2 * n))
 
uavdata_test  <- uavdata[test_id, ]
uavdata_train <- uavdata[-test_id, ]
scores_test   <- scores_df[test_id, ]
scores_train  <- scores_df[-test_id, ]

#######
# Random Forest Regresseion für alle drei NMDS Achsen
# Ablauf genau wie für EIVE-N in Regression
# wiederverwendbare Funktion, damit nicht 3 mal kopiert werden muss
#######
rf_gradient_workflow <- function(zielvariable_name,
                                 uavdata_train,
                                 uavdata_test,
                                 scores_train,
                                 scores_test) {
  
  # Zielvariable
  y_train <- scores_train[[zielvariable_name]]
  y_test  <- scores_test[[zielvariable_name]]

  #############
  # Variablenselektion auf Trainingsdaten
  #############
  vsurf_res <- VSURF(
    x = uavdata_train[, 1:12],
    y = y_train
  )
  
  sel_vars <- colnames(uavdata_train[, 1:12])[
    vsurf_res$varselect.pred
  ]
  
  print(sel_vars)
  
  ############
  # Random Forest mit selektierten Variablen auf Trainingsdaten
  ############
  rf_sel <- randomForest(
    x = uavdata_train[, sel_vars, drop = FALSE],
    y = y_train,
    importance = TRUE,
    keep.forest = TRUE
  )
  
  print(rf_sel)
  
  ###########
  # Variablenwichtigkeit 
  ###########

  importances <- round(importance(rf_sel), 2)
  print(
    importances[
      order(importances[, "%IncMSE"], decreasing = TRUE),
      "%IncMSE"
    ]
  )
  ##########
  # Vorhersage Trainingsdaten
  ##########

  trainpred <- predict(
    rf_sel,
    newdata = uavdata_train[, sel_vars, drop = FALSE]
  )
  
  ###########
  # Vorhersage Testdaten
  ##########

  testpred <- predict(
    rf_sel,
    newdata = uavdata_test[, sel_vars, drop = FALSE]
  )
  
  ###########
  # Guetemasse Training
  ##########

  rsqrt_train <- cor(y_train, trainpred)^2
  rmse_train  <- sqrt(mean((y_train - trainpred)^2))
  nrmse_train <- rmse_train /
    (max(y_train) - min(y_train)) * 100
  
  ##########
  # Guetemasse Test
  ##########
  
  rsqrt_test <- cor(y_test, testpred)^2
  rmse_test  <- sqrt(mean((y_test - testpred)^2))
  nrmse_test <- rmse_test /
    (max(y_test) - min(y_test)) * 100
  
  guete <- data.frame(
    zielvariable = zielvariable_name,
    sel_vars = paste(sel_vars, collapse = ", "),
    r2_train = round(rsqrt_train, 3),
    rmse_train = round(rmse_train, 3),
    nrmse_train = round(nrmse_train, 1),
    r2_test = round(rsqrt_test, 3),
    rmse_test = round(rmse_test, 3),
    nrmse_test = round(nrmse_test, 1)
  )
  
  list(
    modell = rf_sel,
    sel_vars = sel_vars,
    guete = guete
  )
}
 

# Fuer jede der drei Achsen einmal durchlaufen lassen
ergebnis_nmds1 <- rf_gradient_workflow("NMDS1", uavdata_train, uavdata_test,
                                        scores_train, scores_test)
ergebnis_nmds2 <- rf_gradient_workflow("NMDS2", uavdata_train, uavdata_test,
                                        scores_train, scores_test)
ergebnis_nmds3 <- rf_gradient_workflow("NMDS3", uavdata_train, uavdata_test,
                                        scores_train, scores_test)
 
# Guetemasse aller drei Achsen gemeinsam vergleichen
guete_gesamt <- rbind(ergebnis_nmds1$guete,
                       ergebnis_nmds2$guete,
                       ergebnis_nmds3$guete)
guete_gesamt


# Visualisierung 
# R^2 Training vs. Test nebeneinander
guete_r2_long <- guete_gesamt %>%
  select(zielvariable, r2_train, r2_test) %>%
  pivot_longer(cols = c(r2_train, r2_test),
               names_to = "datensatz", values_to = "r2") %>%
  mutate(datensatz = recode(datensatz,
                             r2_train = "Training",
                             r2_test = "Test"))
ggplot(guete_r2_long, aes(x = zielvariable, y = r2, fill = datensatz)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  labs(title = "Modellguete je floristischem Gradienten",
       subtitle = "Grosse Luecke zwischen Training und Test = Hinweis auf instabiles Modell",
       x = "NMDS-Achse", y = expression(R^2), fill = "Datensatz") +
  theme_bw()


# Plot fuer Auswertung
ggplot(guete_r2_long, aes(x = zielvariable, y = r2, fill = datensatz)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  labs(title = "",
       x = "NMDS-Achse", y = expression(R^2), fill = "Datensatz") +
  theme_bw()
 
ggplot(guete_r2_long, aes(x = zielvariable, y = r2, fill = datensatz)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey40") +
  labs(
    x = "NMDS-Achse",
    y = expression(R^2),
    fill = "Datensatz"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

# nRMSE Training vs. Test nebeneinander (normierter Fehler, % vom Wertebereich)
guete_nrmse_long <- guete_gesamt %>%
  select(zielvariable, nrmse_train, nrmse_test) %>%
  pivot_longer(cols = c(nrmse_train, nrmse_test),
               names_to = "datensatz", values_to = "nrmse") %>%
  mutate(datensatz = recode(datensatz,
                             nrmse_train = "Training",
                             nrmse_test = "Test"))
 
ggplot(guete_nrmse_long, aes(x = zielvariable, y = nrmse, fill = datensatz)) +
  geom_col(position = "dodge") +
  labs(title = "Normierter Vorhersagefehler (nRMSE) je floristischem Gradienten",
       x = "NMDS-Achse", y = "nRMSE (%)", fill = "Datensatz") +
  theme_bw()
 
####### 
# Karten erstellen wie in Regression
#  Verwendet das bereits auf 1x1 m aggregierte Rasterbild mit 12 Baendern
#    (r/g/b/RE/NIR/NDVI je mean + sd), z.B. aus Klassifikation:
#    "Drachenwiese26_aggregated_merged.tif"
# Hinzu kommt noch die Baum- und Schattenmaskierung aus mask_trees_shadows.R, diese nutzt ein DEM und dann werden Werte unter den Baeumen rausgestrichen, weil die Schatten die Ergebnisse verfaelschen
#######
dws.agg <- rast("dw26_sued_aggregated.tif")
dwn.agg <- rast("dw26_nord_aggregated.tif")
dw.agg <- merge(dws.agg, dwn.agg)

writeRaster(dw.agg, filename = "Drachenwiese26_aggregated_merged.tif",
            overwrite = TRUE)

#ms <- rast("Drachenwiese26_aggregated_merged.tif")

# Baummaske 
DEM_nord <- rast("drachenwiese_nord25_DEM.tif")
refn <- rast("dw26_nord_aggregated.tif")
DEM_nord <- resample(DEM_nord, refn, method = "max")
treemask.nord <- ifel(DEM_nord > 162.25, 1, NA)

DEM_sued <- rast("drachenwiese_sued25_DEM.tif")
refs <- rast("dw26_sued_aggregated.tif")
DEM_sued <- resample(DEM_sued, refs, method = "max")
treemask.sued <- ifel(DEM_sued > 163.5, 1, NA)

overallext <- ext(xmin(refs), xmax(refn), ymin(refs), ymax(refn))
treemask.sued <- extend(treemask.sued, overallext)
treemask.nord <- resample(treemask.nord, treemask.sued)
tm <- merge(treemask.sued, treemask.nord)
writeRaster(tm, filename = "Drachenwiese26_treemask.tif", overwrite = TRUE)

# Schattenmaske 
ms_roh <- rast("Drachenwiese26_aggregated_merged.tif")
RGB <- subset(ms_roh, c(1, 2, 3))
RGBsum <- sum(RGB)
shadowmask <- ifel(RGBsum < 42000, 1, NA)
writeRaster(shadowmask, filename = "Drachenwiese26_shadowmask.tif", overwrite = TRUE)

# Beide Masken anwenden 
tm <- rast("Drachenwiese26_treemask.tif")
sm <- rast("Drachenwiese26_shadowmask.tif")
ms <- mask(ms_roh, tm, maskvalue = 1)
ms <- mask(ms, sm, maskvalue = 1)
plotRGB(ms, stretch = "lin", main = "Maskiertes Orthofoto (ohne Baeume/Schatten)")
writeRaster(ms, filename = "Drachenwiese26_aggregated_masked.tif", overwrite = TRUE)


pmap_nmds1 <- predict(ms, ergebnis_nmds1$modell, na.rm = TRUE)
pmap_nmds2 <- predict(ms, ergebnis_nmds2$modell, na.rm = TRUE)
pmap_nmds3 <- predict(ms, ergebnis_nmds3$modell, na.rm = TRUE)
 
par(mfrow = c(1, 3))
plot(pmap_nmds1, main = "NMDS1")
plot(pmap_nmds2, main = "NMDS2")
plot(pmap_nmds3, main = "NMDS3")
par(mfrow = c(1, 1))
 
writeRaster(pmap_nmds1, filename = "NMDS1_Drachenwiese_2026.tif", overwrite = TRUE)
writeRaster(pmap_nmds2, filename = "NMDS2_Drachenwiese_2026.tif", overwrite = TRUE)
writeRaster(pmap_nmds3, filename = "NMDS3_Drachenwiese_2026.tif", overwrite = TRUE)
 

########### 
# Korrelation der NMDS Achsen mit den eiv e Zeigerwerten
# pruefen ob sich die statistischen Gradienten (NMDS) mit bekannten oekologischen Gradienten (Feuchte, Naehrstoffe, Reatkionszahl/pH decken)
##########

eive <- read.csv("EIVE_plot_means.csv")
eive <- eive[eive$plot %in% scores_df$plot, ]
eive <- arrange(eive, plot)

identical(eive$plot, scores_df$plot)   # sollte TRUE sein

cor_tab <- data.frame(
  Achse = rep(c("NMDS1", "NMDS2", "NMDS3"), each = 3),
  EIVE  = rep(c("Feuchte (M)", "Naehrstoffe (N)", "Reaktion/pH (R)"), times = 3),
  r = c(
    cor(scores_df$NMDS1, eive$cwm_M), cor(scores_df$NMDS1, eive$cwm_N), cor(scores_df$NMDS1, eive$cwm_R),
    cor(scores_df$NMDS2, eive$cwm_M), cor(scores_df$NMDS2, eive$cwm_N), cor(scores_df$NMDS2, eive$cwm_R),
    cor(scores_df$NMDS3, eive$cwm_M), cor(scores_df$NMDS3, eive$cwm_N), cor(scores_df$NMDS3, eive$cwm_R)
  )
)
cor_tab$r <- round(cor_tab$r, 2)
cor_tab

##########
# Artenscores für alle drei Achsen vervollstaendigen

species_scores <- as.data.frame(scores(nmds, display = "species"))
species_scores$Art <- rownames(species_scores)

# NMDS1
species_scores %>% arrange(NMDS1) %>% select(Art, NMDS1) %>% head(10)
species_scores %>% arrange(desc(NMDS1)) %>% select(Art, NMDS1) %>% head(10)

# NMDS2
species_scores %>% arrange(NMDS2) %>% select(Art, NMDS2) %>% head(10)
species_scores %>% arrange(desc(NMDS2)) %>% select(Art, NMDS2) %>% head(10)

# NMDS3
species_scores %>% arrange(NMDS3) %>% select(Art, NMDS3) %>% head(10)
species_scores %>% arrange(desc(NMDS3)) %>% select(Art, NMDS3) %>% head(10)


##########################################################################################
# Plots fuer Bericht
y_train_1 <- scores_train$NMDS1
y_test_1  <- scores_test$NMDS1

y_train_2 <- scores_train$NMDS2
y_test_2  <- scores_test$NMDS2

y_train_3 <- scores_train$NMDS3
y_test_3  <- scores_test$NMDS3

library(ggplot2)
library(dplyr)

# Vorhersagen für Trainings- und Testdaten
pred_train_1 <- predict(ergebnis_nmds1$modell,
                        newdata = uavdata_train[, ergebnis_nmds1$sel_vars, drop = FALSE])
pred_test_1 <- predict(ergebnis_nmds1$modell,
                       newdata = uavdata_test[, ergebnis_nmds1$sel_vars, drop = FALSE])

pred_train_2 <- predict(ergebnis_nmds2$modell,
                        newdata = uavdata_train[, ergebnis_nmds2$sel_vars, drop = FALSE])
pred_test_2 <- predict(ergebnis_nmds2$modell,
                       newdata = uavdata_test[, ergebnis_nmds2$sel_vars, drop = FALSE])

pred_train_3 <- predict(ergebnis_nmds3$modell,
                        newdata = uavdata_train[, ergebnis_nmds3$sel_vars, drop = FALSE])
pred_test_3 <- predict(ergebnis_nmds3$modell,
                       newdata = uavdata_test[, ergebnis_nmds3$sel_vars, drop = FALSE])


# Daten zusammenführen
plot_data <- bind_rows(
  data.frame(Achse = "NMDS1",
             Beobachtung = y_train_1,
             Vorhersage = pred_train_1,
             Datensatz = "Training"),
  data.frame(Achse = "NMDS1",
             Beobachtung = y_test_1,
             Vorhersage = pred_test_1,
             Datensatz = "Test"),
  
  data.frame(Achse = "NMDS2",
             Beobachtung = y_train_2,
             Vorhersage = pred_train_2,
             Datensatz = "Training"),
  data.frame(Achse = "NMDS2",
             Beobachtung = y_test_2,
             Vorhersage = pred_test_2,
             Datensatz = "Test"),
  
  data.frame(Achse = "NMDS3",
             Beobachtung = y_train_3,
             Vorhersage = pred_train_3,
             Datensatz = "Training"),
  data.frame(Achse = "NMDS3",
             Beobachtung = y_test_3,
             Vorhersage = pred_test_3,
             Datensatz = "Test")
)


# Plot
ggplot(plot_data, aes(x = Beobachtung, y = Vorhersage,
                      shape = Datensatz)) +
  geom_point(size = 2.5) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~ Achse, nrow = 1) +
  labs(
    x = "Beobachteter NMDS-Score",
    y = "Vorhergesagter NMDS-Score",
    shape = "Datensatz"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

##########
# NMDS + EIV fuer grafische Darstellung
##########

# NMDS-Scores
sites <- as.data.frame(
  scores(nmds, display = "sites")
)

sites$plot <- rownames(sites)

# EIV-Daten auf gleiche Plots beschraenken
eive_plot <- eive %>%
  select(plot, cwm_M, cwm_N, cwm_R)

# Zusammenfuehren
nmds_eiv <- sites %>%
  inner_join(eive_plot, by = "plot")

###########
# Korrelation EIV <-> NMDS
###########

cor_eiv <- data.frame(
  Variable = c(
    "Feuchte (M)",
    "Nährstoffe (N)",
    "Reaktion/pH (R)"
  ),
  
  M = c(
    cor(nmds_eiv$NMDS1, nmds_eiv$cwm_M),
    cor(nmds_eiv$NMDS2, nmds_eiv$cwm_M),
    cor(nmds_eiv$NMDS3, nmds_eiv$cwm_M)
  ),
  
  N = c(
    cor(nmds_eiv$NMDS1, nmds_eiv$cwm_N),
    cor(nmds_eiv$NMDS2, nmds_eiv$cwm_N),
    cor(nmds_eiv$NMDS3, nmds_eiv$cwm_N)
  ),
  
  R = c(
    cor(nmds_eiv$NMDS1, nmds_eiv$cwm_R),
    cor(nmds_eiv$NMDS2, nmds_eiv$cwm_R),
    cor(nmds_eiv$NMDS3, nmds_eiv$cwm_R)
  )
)

rownames(cor_eiv) <- c("NMDS1", "NMDS2", "NMDS3")

cor_eiv

##############
# Funktion fuer EIV-Pfeile
##############

get_eiv_vectors <- function(x, y, data) {
  
  vars <- c(
    "cwm_M",
    "cwm_N",
    "cwm_R"
  )
  
  labels <- c(
    "Feuchte (M)",
    "Nährstoffe (N)",
    "Reaktion/pH (R)"
  )
  
  arrows <- data.frame(
    Variable = labels,
    x = sapply(vars, function(v) {
      cor(data[[x]], data[[v]], use = "complete.obs")
    }),
    y = sapply(vars, function(v) {
      cor(data[[y]], data[[v]], use = "complete.obs")
    })
  )
  
  arrows
}


##########
# NMDS-Plot mit EIV-Pfeilen
##########

plot_nmds_eiv <- function(x, y) {
  
  arrows_eiv <- get_eiv_vectors(
    x = x,
    y = y,
    data = nmds_eiv
  )
  
  # Pfeile etwas verlaengern
  arrow_mult <- 1.5
  
  arrows_eiv$xend <- arrows_eiv$x * arrow_mult
  arrows_eiv$yend <- arrows_eiv$y * arrow_mult
  
  ggplot(
    nmds_eiv,
    aes(
      x = .data[[x]],
      y = .data[[y]]
    )
  ) +
    
    geom_point(size = 2.5) +
    
    geom_segment(
      data = arrows_eiv,
      aes(
        x = 0,
        y = 0,
        xend = xend,
        yend = yend
      ),
      arrow = arrow(length = unit(0.25, "cm")),
      inherit.aes = FALSE
    ) +
    
    geom_text(
      data = arrows_eiv,
      aes(
        x = xend,
        y = yend,
        label = Variable
      ),
      inherit.aes = FALSE,
      vjust = -0.5
    ) +
    
    labs(
      x = x,
      y = y
    ) +
    
    theme_bw()
}


p12_eiv <- plot_nmds_eiv("NMDS1", "NMDS2")
p13_eiv <- plot_nmds_eiv("NMDS1", "NMDS3")
p23_eiv <- plot_nmds_eiv("NMDS2", "NMDS3")

p12_eiv
p13_eiv
p23_eiv
