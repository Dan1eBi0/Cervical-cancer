###############################################################################
###  Cervical cancer behavior risk - downloaded data set from                 #
###  https://archive.ics.uci.edu/ml/datasets/Cervical+Cancer+Behavior+Risk    #
###  Author: Daniela Reis                                                     #
###  Data: 10/12/2021                                                         #                                                                   #
###############################################################################


pacotes <- c("readr","plotly","tidyverse", "knitr","kableExtra","car","rgl","gridExtra",
             "PerformanceAnalytics","reshape2","rayshader","psych","pracma",
             "polynom","rqPen","ggrepel","factoextra","cluster", "fpc", "pls")

if(sum(as.numeric(!pacotes %in% installed.packages())) != 0){
  instalador <- pacotes[!pacotes %in% installed.packages()]
  for(i in 1:length(instalador)) {
    install.packages(instalador, dependencies = T)
    break()}
  sapply(pacotes, require, character = T) 
} else {
  sapply(pacotes, require, character = T) 
}


cancer_table <- read.csv(file = "sobar-72.csv",
                         header = TRUE, sep = ",", quote = "\"",
                         dec = ".", fill = TRUE, comment.char = "") #Reading .csv file

glimpse(cancer_table) #Some visualizations


################################### PCA ########################################

rho_ca <- cor(cancer_table[,1:19]) #Creating correlation matrix

chart.Correlation(cancer_table[,1:19])#Viewing correlations

#Building a heat map
rho_ca %>% 
  melt() %>% 
  ggplot() +
  geom_tile(aes(x = Var1, y = Var2, fill = value)) +
  geom_text(aes(x = Var1, y = Var2, label = round(x = value, digits = 3)),
            size = 4) +
  labs(x = NULL,
       y = NULL,
       fill = "Correlações") +
  scale_fill_gradient2(low = "dodgerblue4", 
                       mid = "white", 
                       high = "violetred4",
                       midpoint = 0) +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 0))


KMO(r = rho_ca) #KMO (Kaiser-Meyer-Olkin) Statistics


#Bartlett sphericity test (below). If the p-value is less than 0.05, there are correlations that
#are statistically different from zero, when comparing the correlation matrix with an 
#identity matrix of the same dimension.
cortest.bartlett(R = rho_ca, n = 72) 

#pchisq(q = 828.6449, df = 171, lower.tail = FALSE) 


eigenvalues_rho <- eigen(rho_ca)#Obtainging the eigeinvalues

eigenvalues_rho$values



polinomio_caracteristico <- charpoly(rho_ca) #characteristic polynomial
polinomio_caracteristico 

polinomio_caracteristico_invertido <- as.polynomial(rev(polinomio_caracteristico))

solve(polinomio_caracteristico_invertido) # Rho's eigenvalues


var_compartilhada <- (eigenvalues_rho$values/sum(eigenvalues_rho$values))#Shared variance
var_compartilhada

var_cumulativa <- cumsum(var_compartilhada)
var_cumulativa

principais_componentes <- 1:sum(eigenvalues_rho$values)#Possible factors
principais_componentes

#Joining!
data.frame(principais_componentes = paste0("PC", principais_componentes),
           autovalor = eigenvalues_rho$values,
           var_compartilhada = var_compartilhada,
           var_cumulativa = var_cumulativa) -> relatorio_eigen


relatorio_eigen %>% #First results
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = T, 
                font_size = 12)


relatorio_eigen %>% 
  ggplot(aes(x = principais_componentes, 
             y = var_compartilhada,
             group = 1,
             label = paste0(round(var_compartilhada * 100,
                                  digits = 2), "%"))) +
  geom_col(fill = "dodgerblue4", color = "black") +
  geom_line(color = "darkgoldenrod3",
            size = 1.2) +
  geom_point(size = 2) +
  geom_text(size = 3, vjust = 2, color = "white") +
  labs(x = "Componentes Principais",
       y = "Variância Compartilhada") +
  theme_bw()



eigenvalues_rho$vectors #Eigenvectors


data.frame(eigenvalues_rho$vectors) %>% 
  rename(PC1 = X1, PC2 = X2, PC3 = X3, PC4 = X4, PC5 = X5, PC6 = X6, PC7 = X7, PC8 = X8, PC9 = X9, PC10 = X10, PC11 = X11, PC12 = X12, PC13 = X13, PC14 = X14, PC15 = X15, PC16 = X16, PC17 = X17, PC18 = X18, PC19 = X19) %>% 
  mutate(var = names(cancer_table[1:19])) %>% 
  melt(id.vars = "var") %>% 
  mutate(var = factor(var)) %>% 
  ggplot(aes(x = var, y = value, fill = var)) +
  geom_bar(stat = "identity", color = "black") +
  facet_wrap(~variable) +
  labs(x = NULL, y = NULL, fill = "Legenda:") +
  scale_fill_viridis_d() +
  theme_bw()


L2 <- diag(eigenvalues_rho$values)
L2

prova_01 <- t(eigenvalues_rho$vectors) %*% rho_ca %*% eigenvalues_rho$vectors
round(x = prova_01,
      digits = 14)


scores_fatoriais <- t(eigenvalues_rho$vectors)/sqrt(eigenvalues_rho$values)
scores_fatoriais

cancer_table_std <- cancer_table %>% #Standardizing the database
  scale() %>% 
  data.frame()

fatores <- list()

for(i in 1:nrow(scores_fatoriais)){
  fatores[[i]] <- rowSums(x = sweep(x = cancer_table_std, 
                                    MARGIN = 2, 
                                    STATS = scores_fatoriais[i,], 
                                    FUN = `*`))
}


#warnings()

fatores_df <- data.frame((sapply(X = fatores, FUN = c)))
fatores_df

fatores_df %>%
  rename(F1 = X1,
         F2 = X2,
         F3 = X3,
         F4 = X4,
         F5 = X5,
         F6 = X6,
         F7 = X7,
         F8 = X8,
         F9 = X9,
         F10 = X10,
         F11 = X11,
         F12 = X12,
         F13 = X13,
         F14 = X14,
         F15 = X15,
         F16 = X16,
         F17 = X17,
         F18 = X18,
         F19 = X19) %>% 
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = T, 
                font_size = 12)


round(x = cor(fatores_df), 
      digits = 14)


cancer_table_final <-  cbind(cancer_table,
                             fatores_df) %>% 
  rename(F1 = X1,
         F2 = X2,
         F3 = X3,
         F4 = X4,
         F5 = X5,
         F6 = X6,
         F7 = X7,
         F8 = X8,
         F9 = X9,
         F10 = X10,
         F11 = X11,
         F12 = X12,
         F13 = X13,
         F14 = X14,
         F15 = X15,
         F16 = X16,
         F17 = X17,
         F18 = X18,
         F19 = X19) 


correlacoes_entre_fatores <- cor(cancer_table_final[,1:39])#Correlations of values

correlacoes_entre_fatores %>% 
  melt() %>% 
  filter(Var1 %in% c("F1","F2","F3","F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12", "F13", "F14", "F15", "F16", "F17", "F18", "F19") &
           Var2 %in% c("behavior_sexualRisk","behavior_eating",
                       "behavior_personalHygine", "intention_aggregation",
                       "intention_commitment", "attitude_consistency",
                       "attitude_spontaneity", "norm_significantPerson",
                       "norm_fulfillment", "perception_vulnerability",
                       "perception_severity", "motivation_strength",
                       "motivation_willingness", "socialSupport_emotionality", 
                       "socialSupport_appreciation", "socialSupport_instrumental",
                       "empowerment_knowledge", "empowerment_abilities", 
                       "empowerment_desires")) %>% 
  ggplot() +
  geom_tile(aes(x = Var1, y = Var2, fill = value)) +
  geom_text(aes(x = Var1, y = Var2, label = round(x = value, digits = 3)),
            size = 4) +
  labs(x = NULL,
       y = NULL,
       fill = "Correlações") +
  scale_fill_gradient2(low = "dodgerblue4", 
                       mid = "white", 
                       high = "brown4",
                       midpoint = 0) +
  theme(panel.background = element_rect("white"),
        panel.grid = element_line("grey95"),
        panel.border = element_rect(NA),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 0))


correlacoes_entre_fatores %>% 
  melt() %>% 
  filter(Var1 %in% c("F1","F2","F3","F4", "F5", "F6", "F7", "F8", "F9", "F10", "F11", "F12", "F13", "F14", "F15", "F16", "F17", "F18", "F19") &
           Var2 %in% c("behavior_sexualRisk","behavior_eating",
                       "behavior_personalHygine", "intention_aggregation",
                       "intention_commitment", "attitude_consistency",
                       "attitude_spontaneity", "norm_significantPerson",
                       "norm_fulfillment", "perception_vulnerability",
                       "perception_severity", "motivation_strength",
                       "motivation_willingness", "socialSupport_emotionality", 
                       "socialSupport_appreciation", "socialSupport_instrumental",
                       "empowerment_knowledge", "empowerment_abilities", 
                       "empowerment_desires")) %>% 
  dcast(Var1 ~ Var2) %>% 
  column_to_rownames("Var1") -> correlacoes_entre_fatores_df


correlacoes_entre_fatores_df %>% 
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = T, 
                font_size = 6)


correlacoes_entre_fatores_df %>% 
  mutate(eigenvalues = behavior_sexualRisk^2 + behavior_eating^2 +
           behavior_personalHygine^2 +  intention_aggregation^2 +
           intention_commitment^2 + attitude_consistency^2 +
           attitude_spontaneity^2 + norm_significantPerson^2 +
           norm_fulfillment^2 + perception_vulnerability^2 +
           perception_severity^2 + motivation_strength^2 +
           motivation_willingness^2 + socialSupport_emotionality^2 + 
           socialSupport_appreciation^2 + socialSupport_instrumental^2 +
           empowerment_knowledge^2 + empowerment_abilities^2 + 
           empowerment_desires^2) %>% 
  kable() %>%
  kable_styling(bootstrap_options = "striped", 
                full_width = T, 
                font_size = 6)

correlacoes_entre_fatores_df %>% 
  mutate(eigenvalues = behavior_sexualRisk^2 + behavior_eating^2 +
           behavior_personalHygine^2 +  intention_aggregation^2 +
           intention_commitment^2 + attitude_consistency^2 +
           attitude_spontaneity^2 + norm_significantPerson^2 +
           norm_fulfillment^2 + perception_vulnerability^2 +
           perception_severity^2 + motivation_strength^2 +
           motivation_willingness^2 + socialSupport_emotionality^2 + 
           socialSupport_appreciation^2 + socialSupport_instrumental^2 +
           empowerment_knowledge^2 + empowerment_abilities^2 + 
           empowerment_desires^2) %>% 
  filter(eigenvalues > 1) %>% 
  select(-eigenvalues) %>% 
  t() %>% 
  data.frame() %>% 
  square() %>% 
  mutate(comunalidades = rowSums(.))


correlacoes_entre_fatores_df %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column("variáveis") %>% 
  ggplot(aes(x = F1, y = F2, label = variáveis)) +
  geom_point(color = "dodgerblue4",
             size = 2) +
  geom_text_repel() +
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "red") +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "red") +
  expand_limits(x= c(-1.25, 0.25), y=c(-0.25, 1)) +
  labs(x = paste("Dimensão 1", paste0("(",round(var_compartilhada[1] * 100, 
                                                digits = 2),"%)")),
       y = paste("Dimensão 2", paste0("(",round(var_compartilhada[2] * 100, 
                                                digits = 2),"%)"))) 
theme_bw()

########################### REGRESSÃO LOGÍSTICA ################################