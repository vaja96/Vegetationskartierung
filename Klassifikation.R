#########
# Vegetationskartierung der Drachenwiese, vor Ort und mit multispektralen Orthophotos
# Wie gut lassen sich folgende Vegetationsmerkmale der Drachenwiese anhand räumlich hoch
# aufgeloester Multispektralbilder kartieren?
# Pflanzengesellschaften als floristische Gradienten und in Klassen
#########



#######
# Hier jetzt die Klassen (diskrete Zielvariable), nach dem Dokument Klassifikation
# ergaenzend ist noch die Baum und Schattenmaskierung angewandt damit die Karte nicht für Flaechen mit Baumkronen oder Schatten eine Klasse vorhersagt
# Die gleiche Datenbasis wie beim Ordinationsskript
# Klassen- und Gradientenansatz sollen vergleichbar sein (gleiche Plots, gleicher Seed)
# von Michael Ewald fuer den Kurs Vegetationskartierung und Vegetationsaufnahme fuer das Sommersemester 2026 am KIT
#######

# libraries
library(vegan)         # Ordination (NMDS)
library(dplyr)          # Datenoperationen
library(terra)           # Rasterdatenverarbeitung
library(randomForest)    # Random Forest Modelle
library(VSURF)            # Variablenselektion
library(ggplot2)           # Abbildungen
library(pdp)                # Partial Dependence Plots
library(isopam)
library(caret)
library(lattice)


############
# Was passiert in diesem Skript
# 1. Artentabelle einlesen und transponieren
# 2. Isopam-Klassifikation (bestimmt Klassenzahl selbst)
# 3. Verknuepfung mit Spektraldaten, p34 raus
# 4. Train Test Split (Seed84 wie Ordinationsskript)
# 5. Random Forest Klassifikation
# 6. SVM Klassifikation (Vergleichsmodell)
# 7. Validierung beider Modelle über Konfusionsmatrix
# 8. Baum- und Schattenmaskierung
# 9. Kartierung mit dem maskierten Raster
# 10. Guetemasse fuer den Vergleich mit dem Gradienten Ansatz

# Artentabelle einlesen
sp <- read.csv2("Arten_final26.csv", row.names = 1)
sp <- sp[, -which(colnames(sp) == "p19")]   # keine Koordinaten für p19
sp <- t(sp)
sp <- as.data.frame(sp)


#################
#  Isopam Klassifikation der Vegetationsaufnahmen
# hierarchisches Partionierungsverfahren, welches speziell für Vegetationsaufnahmen entwickelt wurde (Schmidtlein et al.2010)
# isopam bestimmt die Klassenzahl selbst (Datengetrieben), sie wird dementsprechen nicht von uns vorgegeben
#################
ip <- isopam(sp, wordy = FALSE)
 
# "Synthetische" Vegetationstabelle zur inhaltlichen Charakterisierung der Klassen
# Stetigkeiten in %
it <- isotab(ip)
it$tab
 
# Klassenzuordnung extrahieren (Level 2 = feinere Aufteilung, hier 3 Klassen)
classes <- as.factor(ip$flat$lev.2)
classes <- classes[order(names(classes))]
table(classes) #Klassengroesse pruefen

# Isopam hat hier einen hierarchischen Clusterbaum mit 2 Ebenen erzeugt
# Klasse 1 teilt sich in zwei Unterklassen auf (1.1 und 1.2)
# Klasse 2 ist separat
# Hierarchei sichtbar machen: welche Level-2-Klasse gehoert zu welcher Level-1-Gruppe?
hierarchie <- data.frame(
  plot = names(ip$flat$lev.1),
  level1 = ip$flat$lev.1,
  level2 = ip$flat$lev.2
)
table(hierarchie$level1, hierarchie$level2)

###
# Charakterarten je Klasse bestimmen
# it$tb enthaelt die Stetigkeit jeder Art in jeder Klasse mit Signifikanzsternchen (bspw. 80%**)
# Sternchen werden jetzt entfert, Werte numerisch umgewandelt und je Klasse die Arten mit der größten Differenz zwischen der eigenen Stetigkeit und hoechsten Stetigkeiten in den anderen Klassen ausgegeben
# das sind dann genau die Arten, die eine Klasse am staerksten von den anderen unterscheiden --> Charakterarten
### 
tab_numeric <- as.data.frame(lapply(as.data.frame(it$tab), function(x) {
  as.numeric(gsub("[^0-9.]", "", as.character(x)))
}))
rownames(tab_numeric) <- rownames(it$tab)
colnames(tab_numeric) <- colnames(it$tab)
 
top_n <- 10
class_namen <- colnames(tab_numeric)
 
charakterarten <- lapply(class_namen, function(cl) {
  andere <- setdiff(class_namen, cl)
  df <- data.frame(
    Art = rownames(tab_numeric),
    Stetigkeit_eigene_Klasse = tab_numeric[[cl]],
    Stetigkeit_max_andere = apply(tab_numeric[, andere, drop = FALSE], 1, max)
  )
  df$Differenz <- df$Stetigkeit_eigene_Klasse - df$Stetigkeit_max_andere
  df <- df[order(-df$Differenz), ]
  head(df, top_n)
})
names(charakterarten) <- paste0("Klasse_", class_namen)

# Ausgabe je Klasse
charakterarten[["Klasse_1.1"]]   
charakterarten[["Klasse_1.2"]]   
charakterarten[["Klasse_2.0"]]     
 
# Zum direkten Vergleich alle drei Klassen untereinander ausgeben
for (i in seq_along(charakterarten)) {
  cat("\n---", names(charakterarten)[i], "---\n")
  print(charakterarten[[i]])
}

#####################
# Spektraldaten laden und mit den Klassen verknuepfen
#####################
 
uavdata <- read.csv2("extracted_values_Drachenwiese26.csv")
uavdata <- arrange(uavdata, plot)
 
# p34 (nur eine Art, keine sinnvolle Klassenzuordnung) ausschliessen
if ("p34" %in% uavdata$plot) {
  uavdata <- uavdata[-which(uavdata$plot == "p34"), ]
}
 
# Sicherstellen, dass beide Datensaetze exakt dieselbe Plot-Reihenfolge haben
identical(uavdata$plot, names(classes))
# Falls FALSE: Klassen auf die in uavdata vorhandenen Plots einschraenken
classes <- classes[names(classes) %in% uavdata$plot]
uavdata <- uavdata[uavdata$plot %in% names(classes), ]
uavdata <- arrange(uavdata, plot)
classes <- classes[order(names(classes))]
 
identical(uavdata$plot, names(classes))   # ist jetzt TRUE
 
#################
# Trainings und Testdatensatz
# selber seed (84) wie in Ordinationsskript --> identischer split
# Klassen und Gradientenguete so direkt vergleichbar
################
set.seed(84)
n <- nrow(uavdata)
test_id <- sample(seq_len(n), size = round(0.2 * n))
 
uavdata_test  <- uavdata[test_id, ]
uavdata_train <- uavdata[-test_id, ]
classes_test  <- classes[test_id]
classes_train <- classes[-test_id]

####################
# Random Forest Klassifikation
# 10 wiederholungen der Kreuzvalidierung
# --> Abschaetzung der Modellguete waehrend des Trainings (model tuning)
#####################
 
ctrl <- trainControl(method = 'repeatedcv', number = 10, repeats = 10)
 
forest <- train(y = classes_train, x = uavdata_train[, 1:12],
                 method = 'rf', trControl = ctrl)
forest
 
# Variablenwichtigkeit (Mean Decrease Gini), welche Spektralvariablen tragen am meisten zur Klassentrennung bei
forest$finalModel$importance
 
#############
# Support vector machine (polynom kernel) 
# als 2. Modell neben Random Forest --> Einfluss Methodenwahl
###########

svm <- train(y = classes_train,
             x = uavdata_train[, 1:12],
             method = 'svmPoly',
             trControl = ctrl,
             preProcess = c("center", "scale"),
             tuneLength = 3)
svm

# Variablenwichtigkeit über ROC Kurve je Klasse
varImp(svm)

###############
# Unabhaengige Validierung mit dem Testdatensatz für beide Modelle
# Konfusionsmatrix um vorhergesagte und beobachtete Klassen gegenüber zu stellen (Accuracy, Kappa, Sensitivity, Specificity je Klasse)
###############
 
svm_pred <- predict(svm, uavdata_test[, 1:12])
conf_mat_svm <- confusionMatrix(svm_pred, classes_test)
print(conf_mat_svm)
 
# Zum Vergleich auch fuer Random Forest
forest_pred <- predict(forest, uavdata_test[, 1:12])
conf_mat_rf <- confusionMatrix(forest_pred, classes_test)
print(conf_mat_rf)

 
# Kurzer Vergleich der beiden Modelle
data.frame(
  Modell = c("Random Forest", "SVM"),
  Accuracy = c(conf_mat_rf$overall["Accuracy"], conf_mat_svm$overall["Accuracy"]),
  Kappa = c(conf_mat_rf$overall["Kappa"], conf_mat_svm$overall["Kappa"])
)
 
#####
# Baum- und Schattenmaskierung 
#####
 
DEM_nord <- rast("drachenwiese_nord25_DEM.tif")
refn <- rast("dw26_nord_aggregated.tif")
DEM_nord <- resample(DEM_nord, refn, method = "max")
treemask.nord <- ifel(DEM_nord > 162.25, 1, NA)
 
DEM_sued <- rast("drachenwiese_sued25_DEM.tif")
refs <- rast("dw26_sued_aggregated.tif")
DEM_sued <- resample(DEM_sued, refs, method = "max")
treemask.sued <- ifel(DEM_sued > 163.5, 1, NA)
 
# Beide Teilmasken auf eine gemeinsame Ausdehnung bringen und zusammenfuegen
overallext <- ext(xmin(refs), xmax(refn), ymin(refs), ymax(refn))
treemask.sued <- extend(treemask.sued, overallext)
treemask.nord <- resample(treemask.nord, treemask.sued)
tm <- merge(treemask.sued, treemask.nord)
 
plot(tm, main = "Baummaske (1 = Baum/Strauch)")
writeRaster(tm, filename = "Drachenwiese26_treemask.tif", overwrite = TRUE)

# Schattenmaske 
ms_roh <- rast("Drachenwiese26_aggregated_merged.tif")
RGB <- subset(ms_roh, c(1, 2, 3))   # r_mean, g_mean, b_mean
RGBsum <- sum(RGB)
plot(RGBsum, main = "Summe R+G+B (zur Schattenerkennung)")
 
shadowmask <- ifel(RGBsum < 42000, 1, NA)
plot(shadowmask, main = "Schattenmaske (1 = Schatten)")
writeRaster(shadowmask, filename = "Drachenwiese26_shadowmask.tif", overwrite = TRUE)
 
# beide Masken auf Orthofoto anwenden
tm <- rast("Drachenwiese26_treemask.tif")
sm <- rast("Drachenwiese26_shadowmask.tif")
 
ms <- mask(ms_roh, tm, maskvalue = 1)   # Baumpixel -> NA
ms <- mask(ms, sm, maskvalue = 1)        # Schattenpixel -> NA
 
plotRGB(ms, stretch = "lin",
        main = "Maskiertes Orthofoto (ohne Baeume/Schatten)")
writeRaster(ms, filename = "Drachenwiese26_aggregated_masked.tif", overwrite = TRUE)
 
#############
# Kartierung
# Klassifikation in die Flaeche übertragen, verwendet das maskierte Raster aus dem Schritt davor
#############

# Random Forest Karte
pred_rf <- terra::predict(ms, model = forest, na.rm = TRUE)
coltab(pred_rf) <- data.frame(value = 1:3,
                               color = c("#FDE725", "#21908C", "#440154"))
plot(pred_rf, type = "classes", legend = "bottomright",
     main = "Klassenkarte - Random Forest (maskiert)")
 
# SVM Karte
pred_svm <- terra::predict(ms, model = svm, na.rm = TRUE)
coltab(pred_svm) <- data.frame(value = 1:3,
                                color = c("#FDE725", "#21908C", "#440154"))
plot(pred_svm, type = "classes", legend = "bottomright",
     main = "Klassenkarte - SVM (maskiert)")
 
writeRaster(pred_rf,  filename = "Klassen_RF_Drachenwiese_2026_masked.tif",  overwrite = TRUE)
writeRaster(pred_svm, filename = "Klassen_SVM_Drachenwiese_2026_masked.tif", overwrite = TRUE)


# Karte fuer Auswertung
# Random Forest Karte
pred_rf <- terra::predict(ms, model = forest, na.rm = TRUE)

# Klassen benennen
levels(pred_rf) <- data.frame(
  value = 1:3,
  class = c("1.1", "1.2", "2.0")
)

coltab(pred_rf) <- data.frame(
  value = 1:3,
  color = c("#FDE725", "#21908C", "#440154")
)


# SVM Karte
pred_svm <- terra::predict(ms, model = svm, na.rm = TRUE)

# Klassen benennen
levels(pred_svm) <- data.frame(
  value = 1:3,
  class = c("1.1", "1.2", "2.0")
)

coltab(pred_svm) <- data.frame(
  value = 1:3,
  color = c("#FDE725", "#21908C", "#440154")
)


# Beide nebeneinander plotten
par(mfrow = c(1, 2))

plot(
  pred_rf,
  type = "classes",
  legend = "bottomright",
  main = "Random Forest"
)

plot(
  pred_svm,
  type = "classes",
  legend = "bottomright",
  main = "SVM"
)

par(mfrow = c(1, 1))


# Speichern
writeRaster(
  pred_rf,
  filename = "Klassen_RF_Drachenwiese_2026_masked.tif",
  overwrite = TRUE
)

writeRaster(
  pred_svm,
  filename = "Klassen_SVM_Drachenwiese_2026_masked.tif",
  overwrite = TRUE
)

########### 
# Guetemasse, Accuracy und Kappa aus den beiden Konfusionsmatritzen, den R2 und RMSE Werten der Floristischen Gradienten (Ordinationsskript) gegenueberstellen um die eigentliche Fragestellung: Wie gut lassen sich 
# Pflanzengesellschaften als floristische Gradienten und in Klassen kartieren?
###########


conf_mat_rf$overall[c("Accuracy", "Kappa")]
conf_mat_svm$overall[c("Accuracy", "Kappa")]
