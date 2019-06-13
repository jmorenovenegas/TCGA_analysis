source("modulesR/0_loadLibraries.R")

loadpkg('dplyr')
load('Data/sample_data.RData')

#Cargamos las muestras clinicas pertenecientes a pacientes triple negativo y luminal procedentes de TCGA.
tnbc_samples <- sample_data %>% dplyr::filter(`ER Status` == "Negative" & `PR Status` == "Negative" & `HER2 Final Status` == "Negative" & `PAM50 mRNA` != "Luminal A")
luminal_samples <- sample_data %>% dplyr::filter(`PAM50 mRNA` == "Luminal A")
mydata <- rbind(tnbc_samples, luminal_samples)
mydata <- mydata %>% select(`Complete TCGA ID`, `OS event`, `OS Time`)
mydata <- mydata %>% mutate(`OS event`= as.numeric(`OS event`), `OS Time`= as.numeric(`OS Time`)) 

# Cargamos el objeto voom obtenido en el analisis de expresion diferencial y cargar la matriz de expresion diferencial. El objeto voom presenta una matriz numerica con los valores de expresion genetica nomalizados en escala logCPM (log2 counts per million).
loadpkg("limma")
load('Data/voom_object.RData')

# Cargar los genes significativos. El archivo csv procede del analisis de expresion diferencial realizado previamente.
sigGenes <- as.character(read.csv('Results/TNBCvsLumA_dea.csv'
                                  , header = TRUE, sep = ',', quote = '"')$Gene)

# Extraemos la matriz de datos de expresion diferencial del objeto voom y la filtramos para quedarnos con los genes significativos.
expr <- v$E[sigGenes,]
expr <- t(expr)
ID <- row.names(expr)
expr <- as.data.frame(expr, row.names = 1:length(row.names(expr)))

#Hacemos un inner join de la matriz de expresion diferencial y las muestras de pacientes.
expr <- mutate(expr, `Complete TCGA ID` = ID)
mydata <- inner_join(expr, mydata)

loadpkg("gplots")
loadpkg("factoextra")
loadpkg("pROC")

# Visualizamos la matriz de expresion mediante un heatmap.

pdf(file = "Results/TNBCvsLumA_heatmap.pdf", width = 9, height = 4.5)
heatmap.2(scale(expr[,-127]), col = redgreen(75), trace = 'none', density.info = 'none')
dev.off()

set.seed(123)
k.clust <- kmeans(scale(expr[,-127]), 2)

# Visualizar resultado del clustering

pdf(file = "Results/TNBCvsLumA_clustering.pdf", width = 9, height = 4.5)
fviz_cluster(k.clust, expr[,-127],geom = 'point', main = "kmeans clustering")
dev.off()

mydata <- mutate(mydata, predicted_label = as.factor(k.clust$cluster))
data.table::setnames(mydata, "OS Time", "time")
data.table::setnames(mydata, "OS event", "status")

loadpkg("survival")
loadpkg("survminer")

# Visualizar curvas de supervivencia de los clusters

set.seed(123)
pdf(file = "Results/predicted_groups_survival_curves.pdf", width = 9, height = 4.5)
ggsurvplot(survfit(Surv(time ,status) ~ predicted_label, data = dplyr::filter(mydata, time <= 1200))
           , conf.int = FALSE, censor = TRUE, main = "Survival"
           , pval = TRUE, xlab="Time in years"
           , legend.title = 'Kaplan-Meier curve by predicted phenotype'
           , legend.labs = c('Predicted TNBC', 'Predicted LUM')
           , palette = "simpsons"
           , xscale = 365.25)
dev.off()