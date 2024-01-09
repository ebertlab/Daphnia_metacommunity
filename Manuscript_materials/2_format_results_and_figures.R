# Long-term species coexistence in a metacommunity is compatible with 
# niche differences and strong competition at a local scale
# Study by: Maxime Dubart, P.David, F.Ben-Ami, C.R. Haag, P.D. Fields, V.I. Pajunen & D.Ebert
# title2: "Summary model 15b"
# author of script: "Maxime Dubart"
# date: "25 01 2018"
# updated "04 01 2024" by dieter ebert


library(rjags)
library(ggplot2)
library(cowplot)

trunc <- function(x, prec = 2) formatC(x, digits = prec, format = "f")

### Read results ----------
# # 1000 saved iterations (dev : 21020), bin 2000, ite 4000, thin 2 (bin in included in niter)
# js = as.mcmc.list(list(as.mcmc(readRDS("./results/jsample.full.allsp.model_15_34522179.RDS")$BUGSoutput$sims.matrix),
#                         as.mcmc(readRDS("./results/jsample.full.allsp.model_15_94462948.RDS")$BUGSoutput$sims.matrix)))
# 
# # 2500 saved (dev = 21110), bin = 5000, it? = 5000, thin = 2, (bin included in niter)
# js = as.mcmc.list(list(as.mcmc(readRDS("./results/jsample.full.allsp.model_15_93694702.RDS")$BUGSoutput$sims.matrix),
#                         as.mcmc(readRDS("./results/jsample.full.allsp.model_15_80227941.RDS")$BUGSoutput$sims.matrix)))

# 5000 saved iterations (dev : 7660.35556), bin 10000, ite 10000, thin 2
js = as.mcmc.list(list(as.mcmc(readRDS("./results/jsample.full.allsp.model_15_90837942.RDS")$BUGSoutput$sims.matrix),
                       as.mcmc(readRDS("./results/jsample.full.allsp.model_15_81433744.RDS")$BUGSoutput$sims.matrix),
                       as.mcmc(readRDS("./results/jsample.full.allsp.model_15_81704949.RDS")$BUGSoutput$sims.matrix),
                       as.mcmc(readRDS("./results/jsample.full.allsp.model_15_26143440.RDS")$BUGSoutput$sims.matrix),
                       as.mcmc(readRDS("./results/jsample.full.allsp.model_15_80784603.RDS")$BUGSoutput$sims.matrix),
                       as.mcmc(readRDS("./results/jsample.full.allsp.model_15_92051553.RDS")$BUGSoutput$sims.matrix)))


# ggmcmc(ggs(jsL), file = "ggmcmc_longispina.pdf")
jsL = data.frame()
jsM = data.frame()
jsP = data.frame()

for(i in 1:length(js)){
  jsL = rbind(jsL, js[[i]][,c(1,4,10,5,11,16,19,25,31,37,20,26,32,38,43)])
  jsM = rbind(jsM, js[[i]][,c(2,6,12,7,13,17,21,27,33,39,22,28,34,40,43)])
  jsP = rbind(jsP, js[[i]][,c(3,8,14,9,15,18,23,29,35,41,24,30,36,42,43)])
}


jsLL = jsL[,-15]
colnames(jsLL) = c("alpha","p1_e","p1_c","p2_e","p2_c","d","beta1_PC1_e","beta1_PC2_e","beta1_PC1_c","beta1_PC2_c",
                   "beta2_PC1_e","beta2_PC2_e","beta2_PC1_c","beta2_PC2_c")
# write.table(jsLL, "./simulated/data/jsL", row.names = F, col.names = T)

jsMM = jsM[,-15]
colnames(jsMM) = c("alpha","p1_e","p1_c","p2_e","p2_c","d","beta1_PC1_e","beta1_PC2_e","beta1_PC1_c","beta1_PC2_c",
                   "beta2_PC1_e","beta2_PC2_e","beta2_PC1_c","beta2_PC2_c")
# write.table(jsMM, "./simulated/data/jsM", row.names = F, col.names = T)

jsPP = jsP[,-15]
colnames(jsPP) = c("alpha","p1_e","p1_c","p2_e","p2_c","d","beta1_PC1_e","beta1_PC2_e","beta1_PC1_c","beta1_PC2_c",
                   "beta2_PC1_e","beta2_PC2_e","beta2_PC1_c","beta2_PC2_c")
# write.table(jsPP, "./simulated/data/jsP", row.names = F, col.names = T)

### Does e-c cross zero in winter ?
quantile(exp(jsMM$p2_e)-exp(jsMM$p2_c),c(0.025,0.5,0.975))
quantile(exp(jsLL$p2_e)-exp(jsLL$p2_c),c(0.025,0.5,0.975))
quantile(exp(jsPP$p2_e)-exp(jsPP$p2_c),c(0.025,0.5,0.975))

## Transform to readable format :

jsM[,2:5] <- exp(jsM[,2:5])
jsM <- t(apply(jsM, 2, quantile, c(0.5,0.025,0.975)))

jsL[,2:5] <- exp(jsL[,2:5])
jsL <- t(apply(jsL, 2, quantile, c(0.5,0.025,0.975)))

jsP[,2:5] <- exp(jsP[,2:5])
jsP <- t(apply(jsP, 2, quantile, c(0.5,0.025,0.975)))

siz <- 15
df <- as.data.frame(rbind(cbind(rep('D. longispina',siz), 
                                c('both','sum.','sum.','win.','win.','both','sum','sum','sum','sum','win.','win.','win.','win.','dev'), 
                                c('alpha','e','c','e','c','d','PC1_es','PC2_es','PC1_cs','PC2_cs','PC1_ew','PC2_ew','PC1_cw','PC2_cw','dev') ,
                                (jsL) ),
                          cbind(rep('D. magna',siz), 
                                c('both','sum.','sum.','win.','win.','both','sum','sum','sum','sum','win.','win.','win.','win.','dev'), 
                                c('alpha','e','c','e','c','d','PC1_es','PC2_es','PC1_cs','PC2_cs','PC1_ew','PC2_ew','PC1_cw','PC2_cw','dev') ,
                                (jsM) ) ,
                          cbind(rep('D. pulex',siz), 
                                c('both','sum.','sum.','win.','win.','both','sum','sum','sum','sum','win.','win.','win.','win.','dev'), 
                                c('alpha','e','c','e','c','d','PC1_es','PC2_es','PC1_cs','PC2_cs','PC1_ew','PC2_ew','PC1_cw','PC2_cw','dev') ,
                                (jsP) )))
rownames(df) <- NULL
colnames(df)[1:3] <- c("species","period","parameter")
df[,4] <- as.numeric(as.character(df[,4]))
df[,5] <- as.numeric(as.character(df[,5]))
df[,6] <- as.numeric(as.character(df[,6]))

df2 <- df
df <- df2
str(df)
df$species = factor(df$species)
df$species = factor(df$species,levels(df$species[c(2,1,3)]))


###.   ----

## Figure 1 ------
 #### Plot colonization distances  ------
dec_exp = ggplot(data = data.frame(x = seq(0,200)), aes(x = x))+
  stat_function(fun = function(x) exp(- 0.0348 * x), linetype = "dashed", linewidth = 0.6)+
  geom_vline(xintercept = log(2)/0.0348, linetype = "dashed", linewidth = 0.6)+
  stat_function(fun = function(x) exp(- 0.0425 * x), linetype = "dotted", linewidth = 0.6)+
  geom_vline(xintercept = log(2)/0.0425, linetype = "dotted", linewidth = 0.6)+
  stat_function(fun = function(x) exp(- 0.0432 * x), linewidth = 0.6)+
  geom_vline(xintercept = log(2)/0.0432, linewidth = 0.6)+
  geom_hline(yintercept = 0.5, color = "grey", linetype = "dashed")+
  labs(x = "Distance (m)", y = "Contribution")

dec_exp

#### Plot detection probability colonization distances  ------
pd_plot <- plot_grid(
  ggplot(df[df$parameter %in% c('d'),], aes(`50%`, x=species))+
    geom_errorbar(aes(ymax=`97.5%`, ymin = `2.5%`), width = 0.25,position = position_dodge(0.4))+
    geom_jitter(position = position_dodge(0.4), size = 3)+
    labs(y = "Detection probability")+scale_y_continuous(limits = c(0,1))+
    theme(legend.position = "none", 
          axis.text.x = element_text(face = "italic")) + #, axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+
    labs(x = ""),
  ggplot(df[df$parameter %in% c('alpha'),], aes(log(2)/ `50%`, x=species))+
    geom_errorbar(aes(ymax=log(2)/ `97.5%`, ymin = log(2)/ `2.5%`), width = 0.25,position = position_dodge(0.4))+
    geom_jitter(position = position_dodge(0.4), size = 3)+
    labs(y = "Median colonization distance (m)")+
    scale_y_continuous(limits = c(0,50))+
    theme(legend.position = "none", 
          axis.text.x = element_text(face = "italic")) + #, axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))+
    labs(x = "")+
    annotation_custom(ggplotGrob(dec_exp), xmin = 0.8, xmax = 3.5, ymin = 27, ymax = 50),
  labels = c("A","B")
)
pd_plot


df_rates = df[df$parameter %in% c('e','c'),]
df_rates$pam = paste(df_rates$period, df_rates$parameter,sep="")

df_rate_scaled_by_month = df_rates
df_rate_scaled_by_month[df_rate_scaled_by_month$period == "win.",4:6] = df_rate_scaled_by_month[df_rate_scaled_by_month$period == "win.",4:6] / 10
df_rate_scaled_by_month[df_rate_scaled_by_month$period == "sum.",4:6] = df_rate_scaled_by_month[df_rate_scaled_by_month$period == "sum.",4:6] / 2

scaled_rates = ggplot(df_rate_scaled_by_month, aes(`50%`, x=species, shape = pam ))+
  geom_errorbar(aes(ymax=`97.5%`, ymin = `2.5%`), width = 0.25,position = position_dodge(0.4))+
  geom_jitter(position = position_dodge(0.4), size = 2)+
  scale_shape_manual(values = c(19,17,1,2) ,labels = c("Colonization rate (Summer)",
                                                       "Extinction rate (Summer)",
                                                       "Colonization rate (Winter)",
                                                       "Extinction rate (Winter)"))+
  labs(y = "rates (e,c)", x = "", shape = "")+
  scale_y_continuous(limits = c(0,0.3))+
  theme(legend.position = c(0.895,0.85), axis.text.x = element_text(face = "italic"))

#scaled_rates
#ggsave(scaled_rates, filename = "./figures/base_demo_rescaled.jpeg", height = 4, width = 7.5, dpi = 600)


#### Plot grid with detection prob., col. distance and extinction/colonization rates  ------

base_param = plot_grid(
  pd_plot,
  ggplot(df_rates, aes(`50%`, x=species, shape = pam ))+
    geom_errorbar(aes(ymax=`97.5%`, ymin = `2.5%`), width = 0.25,position = position_dodge(0.4))+
    geom_jitter(position = position_dodge(0.4), size = 3)+
    scale_shape_manual(values = c(19,17,1,2) ,labels = c("Colonization rate (Summer)",
                                                         "Extinction rate (Summer)",
                                                         "Colonization rate (Winter)",
                                                         "Extinction rate (Winter)"))+
    labs(y = "Extinction / Colonization rate", x = "", shape = "")+
    scale_y_continuous(limits = c(0,0.5))+
    theme(legend.position = c(0.85,0.85), 
          axis.text.x = element_text(face = "italic")), ncol = 1, labels = c("","C"))
base_param

#### save Figure 1 (filename = "base_demo") ----
ggsave(base_param, filename = "./figures/base_demo.png", height = 7.5, width = 7.5,  dpi = 600)
ggsave(base_param, filename = "./figures/base_demo.pdf", height = 7.5, width = 7.5,  dpi = 600)


###.   ----

## Figure 2: environmental interactions  ------

cbbPalette <- c("gray0", "gray44", "white", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cov_pc1 <- ggplot(df[df$parameter %in% c('PC1_es','PC1_cs','PC1_ew','PC1_cw'),], aes(`50%`, x= parameter, fill=species))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_bar(stat='identity', position = position_dodge())+
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),position = position_dodge(.9), width = 0.4)+
  scale_fill_manual(values=cbbPalette, name = "Species")+
  theme(legend.position = c(0.1,0.83), legend.text = element_text(face="italic"))+
  labs(y = "Median effect (95 % CI)", x = "")+
  scale_y_continuous(limits=c(-2.1,2))+
  scale_x_discrete(labels = c("Summer col.","Winter col.","Summer ext.","Winter ext."))

cov_pc2 <- ggplot(df[df$parameter %in% c('PC2_es','PC2_cs','PC2_ew','PC2_cw'),], aes(`50%`, x= parameter, fill=species))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_bar(stat='identity', position = position_dodge())+
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),position = position_dodge(.9), width = 0.4)+
  scale_fill_manual(values=cbbPalette)+
  theme(legend.position = "none")+
  labs(y = "Median effect (95 % CI)", x = "")+
  scale_y_continuous(limits=c(-2.1,2))+
  scale_x_discrete(labels = c("Summer col.","Winter col.","Summer ext.","Winter ext."))

cov = plot_grid(cov_pc1, cov_pc2, ncol = 1, labels = c("A","B"))
cov

#### save Figure 2. ------
ggsave(filename = "./figures/env_interactions_cov_15b.jpeg", cov, width = 7, height = 10, dpi = 600)
ggsave(filename = "./figures/env_interactions_cov_15bPDF.pdf", cov, width = 7, height = 10, dpi = 600)
ggsave(filename = "./figures/env_interactions_cov_15bPNG.png", cov, width = 7, height = 10, dpi = 600)

####.    -----
## Interaction matrix. ----

tempoL = jsLL
tempoL$beta1_PC1_e_scaled = ((tempoL$p1_e - tempoL$beta1_PC1_e) - (tempoL$p1_e + tempoL$beta1_PC1_e)) / tempoL$p1_e
tempoL$beta1_PC2_e_scaled = ((tempoL$p1_e - tempoL$beta1_PC2_e) - (tempoL$p1_e + tempoL$beta1_PC2_e)) / tempoL$p1_e
tempoL$beta1_PC1_c_scaled = ((tempoL$p1_c - tempoL$beta1_PC1_c) - (tempoL$p1_c + tempoL$beta1_PC1_c)) / tempoL$p1_c
tempoL$beta1_PC2_c_scaled = ((tempoL$p1_c - tempoL$beta1_PC2_c) - (tempoL$p1_c + tempoL$beta1_PC2_c)) / tempoL$p1_c

tempoL$beta2_PC1_e_scaled = ((tempoL$p2_e - tempoL$beta2_PC1_e) - (tempoL$p2_e + tempoL$beta2_PC1_e)) / tempoL$p2_e
tempoL$beta2_PC2_e_scaled = ((tempoL$p2_e - tempoL$beta2_PC2_e) - (tempoL$p2_e + tempoL$beta2_PC2_e)) / tempoL$p2_e
tempoL$beta2_PC1_c_scaled = ((tempoL$p2_c - tempoL$beta2_PC1_c) - (tempoL$p2_c + tempoL$beta2_PC1_c)) / tempoL$p2_c
tempoL$beta2_PC2_c_scaled = ((tempoL$p2_c - tempoL$beta2_PC2_c) - (tempoL$p2_c + tempoL$beta2_PC2_c)) / tempoL$p2_c

Lscaled = t(as.data.frame((apply(tempoL, 2, quantile, c(0.025,0.5,0.975))))[,15:22])

tempoL = jsMM
tempoL$beta1_PC1_e_scaled = ((tempoL$p1_e - tempoL$beta1_PC1_e) - (tempoL$p1_e + tempoL$beta1_PC1_e)) / tempoL$p1_e
tempoL$beta1_PC2_e_scaled = ((tempoL$p1_e - tempoL$beta1_PC2_e) - (tempoL$p1_e + tempoL$beta1_PC2_e)) / tempoL$p1_e
tempoL$beta1_PC1_c_scaled = ((tempoL$p1_c - tempoL$beta1_PC1_c) - (tempoL$p1_c + tempoL$beta1_PC1_c)) / tempoL$p1_c
tempoL$beta1_PC2_c_scaled = ((tempoL$p1_c - tempoL$beta1_PC2_c) - (tempoL$p1_c + tempoL$beta1_PC2_c)) / tempoL$p1_c

tempoL$beta2_PC1_e_scaled = ((tempoL$p2_e - tempoL$beta2_PC1_e) - (tempoL$p2_e + tempoL$beta2_PC1_e)) / tempoL$p2_e
tempoL$beta2_PC2_e_scaled = ((tempoL$p2_e - tempoL$beta2_PC2_e) - (tempoL$p2_e + tempoL$beta2_PC2_e)) / tempoL$p2_e
tempoL$beta2_PC1_c_scaled = ((tempoL$p2_c - tempoL$beta2_PC1_c) - (tempoL$p2_c + tempoL$beta2_PC1_c)) / tempoL$p2_c
tempoL$beta2_PC2_c_scaled = ((tempoL$p2_c - tempoL$beta2_PC2_c) - (tempoL$p2_c + tempoL$beta2_PC2_c)) / tempoL$p2_c

Mscaled = t(as.data.frame((apply(tempoL, 2, quantile, c(0.025,0.5,0.975))))[,15:22])

tempoL = jsPP
tempoL$beta1_PC1_e_scaled = ((tempoL$p1_e - tempoL$beta1_PC1_e) - (tempoL$p1_e + tempoL$beta1_PC1_e)) / tempoL$p1_e
tempoL$beta1_PC2_e_scaled = ((tempoL$p1_e - tempoL$beta1_PC2_e) - (tempoL$p1_e + tempoL$beta1_PC2_e)) / tempoL$p1_e
tempoL$beta1_PC1_c_scaled = ((tempoL$p1_c - tempoL$beta1_PC1_c) - (tempoL$p1_c + tempoL$beta1_PC1_c)) / tempoL$p1_c
tempoL$beta1_PC2_c_scaled = ((tempoL$p1_c - tempoL$beta1_PC2_c) - (tempoL$p1_c + tempoL$beta1_PC2_c)) / tempoL$p1_c

tempoL$beta2_PC1_e_scaled = ((tempoL$p2_e - tempoL$beta2_PC1_e) - (tempoL$p2_e + tempoL$beta2_PC1_e)) / tempoL$p2_e
tempoL$beta2_PC2_e_scaled = ((tempoL$p2_e - tempoL$beta2_PC2_e) - (tempoL$p2_e + tempoL$beta2_PC2_e)) / tempoL$p2_e
tempoL$beta2_PC1_c_scaled = ((tempoL$p2_c - tempoL$beta2_PC1_c) - (tempoL$p2_c + tempoL$beta2_PC1_c)) / tempoL$p2_c
tempoL$beta2_PC2_c_scaled = ((tempoL$p2_c - tempoL$beta2_PC2_c) - (tempoL$p2_c + tempoL$beta2_PC2_c)) / tempoL$p2_c

Pscaled = t(as.data.frame((apply(tempoL, 2, quantile, c(0.025,0.5,0.975))))[,15:22])

siz <- 8
dfscaled <- as.data.frame(rbind(cbind(rep('D. longispina',siz), 
                                      c('sum','sum','sum','sum','win.','win.','win.','win.'), 
                                      c('PC1_es','PC2_es','PC1_cs','PC2_cs','PC1_ew','PC2_ew','PC1_cw','PC2_cw') ,
                                      (Lscaled) ),
                                cbind(rep('D. magna',siz), 
                                      c('sum','sum','sum','sum','win.','win.','win.','win.'), 
                                      c('PC1_es','PC2_es','PC1_cs','PC2_cs','PC1_ew','PC2_ew','PC1_cw','PC2_cw') ,
                                      (Mscaled) ) ,
                                cbind(rep('D. pulex',siz), 
                                      c('sum','sum','sum','sum','win.','win.','win.','win.'), 
                                      c('PC1_es','PC2_es','PC1_cs','PC2_cs','PC1_ew','PC2_ew','PC1_cw','PC2_cw') ,
                                      (Pscaled) )))

rownames(dfscaled) <- NULL
colnames(dfscaled)[1:3] <- c("species","period","parameter")
dfscaled[,4] <- as.numeric(as.character(dfscaled[,4]))
dfscaled[,5] <- as.numeric(as.character(dfscaled[,5]))
dfscaled[,6] <- as.numeric(as.character(dfscaled[,6]))

#dfscaled$species = factor(dfscaled$species, levels = levels(dfscaled$species)[c(2,1,3)])
dfscaled$species = factor(dfscaled$species)

#cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("gray0", "gray44", "white", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cov_pc1_scaled <- ggplot(dfscaled[dfscaled$parameter %in% c('PC1_es','PC1_cs','PC1_ew','PC1_cw'),], aes(`50%`, x= parameter, fill=species))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_bar(stat='identity', position = position_dodge())+
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),position = position_dodge(.9), width = 0.4)+
  scale_fill_manual(values=cbbPalette, name = "Species")+
  theme(legend.position = c(0.05,0.8), legend.text = element_text(face="italic"))+
  labs(y = "Median effect (95 % CI)", x = "")+
  scale_y_continuous(limits=c(-2.2,1.7))+
  scale_x_discrete(labels = c("Summer col.","Winter col.","Summer ext.","Winter ext."))

cov_pc2_scaled <- ggplot(dfscaled[dfscaled$parameter %in% c('PC2_es','PC2_cs','PC2_ew','PC2_cw'),], aes(`50%`, x= parameter, fill=species))+
  geom_hline(yintercept = 0, linetype="dashed")+
  geom_bar(stat='identity', position = position_dodge())+
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`),position = position_dodge(.9), width = 0.4)+
  scale_fill_manual(values=cbbPalette)+
  theme(legend.position = "none", panel.background = element_rect("white", "black"),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "lightgray"))+
#  theme_light(legend.position = "none")+
  labs(y = "Median effect (95 % CI)", x = "")+
  scale_y_continuous(limits=c(-2.2,1.7))+
  scale_x_discrete(labels = c("Summer col.","Winter col.","Summer ext.","Winter ext."))

cov_scaled = plot_grid(cov_pc1_scaled, cov_pc2_scaled, ncol = 1, labels = c("A","B"))
#cov_scaled
#ggsave(filename = "./figures/environmental_images/cov_15b_scaled.jpeg", cov_scaled, width = 7, height = 10, dpi = 600)
#ggsave(filename = "./figures/environmental_images/cov_15b_scaled_pdf.pdf", cov_scaled, width = 7, height = 10, dpi = 600)



## Figure 4: Interaction matrix. ----

gamma = data.frame()
gamma_scaled = data.frame()

for(i in 1:length(js)){
  gamma = rbind(gamma, js[[i]][,44:79])
  gamma_scaled = rbind(gamma_scaled, js[[i]][,c(4:15,44:79)])
}
# gamma = rbind(js[[1]][, 44:79], js[[2]][, 44:79], js[[3]][, 44:79], js[[4]][, 44:79])
# # gamma = rbind(js[[1]][, 44:79], js[[2]][, 44:79], js[[3]][, 44:79])


gamma_for_simu = gamma[,c(-1,-2,-9,-10,-17,-18,-19,-20,-27,-28,-35,-36)]
# write.table(gamma_for_simu, "./simulated/data/gamma", row.names = F, col.names = T)

gamma = apply(gamma, 2, quantile, c(0.025,0.5,0.975))

## Or with transformed param : :
# gamma = apply(gamma_scaled, 2, quantile, c(0.025,0.5,0.975))

tmpStructure = structure(rep(0, 3*3*3), .Dim = c(3,3,3))

gammas = list(winter = list(e = tmpStructure, c = tmpStructure), summer = list(e = tmpStructure, c = tmpStructure))

# Subset things
gammas$summer$c[,,1] = matrix(gamma[,grep("gamma_c[1", colnames(gamma), fixed = T)][2,], byrow = F, ncol = 3)
gammas$summer$c[,,2] = matrix(gamma[,grep("gamma_c[1", colnames(gamma), fixed = T)][1,], byrow = F, ncol = 3)
gammas$summer$c[,,3] = matrix(gamma[,grep("gamma_c[1", colnames(gamma), fixed = T)][3,], byrow = F, ncol = 3)

gammas$winter$c[,,1] = matrix(gamma[,grep("gamma_c[2", colnames(gamma), fixed = T)][2,], byrow = F, ncol = 3)
gammas$winter$c[,,2] = matrix(gamma[,grep("gamma_c[2", colnames(gamma), fixed = T)][1,], byrow = F, ncol = 3)
gammas$winter$c[,,3] = matrix(gamma[,grep("gamma_c[2", colnames(gamma), fixed = T)][3,], byrow = F, ncol = 3)

gammas$summer$e[,,1] = matrix(gamma[,grep("gamma_e[1", colnames(gamma), fixed = T)][2,], byrow = F, ncol = 3)
gammas$summer$e[,,2] = matrix(gamma[,grep("gamma_e[1", colnames(gamma), fixed = T)][1,], byrow = F, ncol = 3)
gammas$summer$e[,,3] = matrix(gamma[,grep("gamma_e[1", colnames(gamma), fixed = T)][3,], byrow = F, ncol = 3)

gammas$winter$e[,,1] = matrix(gamma[,grep("gamma_e[2", colnames(gamma), fixed = T)][2,], byrow = F, ncol = 3)
gammas$winter$e[,,2] = matrix(gamma[,grep("gamma_e[2", colnames(gamma), fixed = T)][1,], byrow = F, ncol = 3)
gammas$winter$e[,,3] = matrix(gamma[,grep("gamma_e[2", colnames(gamma), fixed = T)][3,], byrow = F, ncol = 3)

library(reshape2)
### when showing exp() (true values)
interactions_on_c_summer = melt(exp(gammas$summer$c))
interactions_on_e_summer = melt(exp(gammas$summer$e))
interactions_on_c_winter = melt(exp(gammas$winter$c))
interactions_on_e_winter = melt(exp(gammas$winter$e))


sp = c("D. longispina","D. magna","D. pulex")

for(i in 1:length(sp)){
  interactions_on_c_summer$Var1[interactions_on_c_summer$Var1 == i] = sp[i]
  interactions_on_c_summer$Var2[interactions_on_c_summer$Var2 == i] = sp[i]
  
  interactions_on_e_summer$Var1[interactions_on_e_summer$Var1 == i] = sp[i]
  interactions_on_e_summer$Var2[interactions_on_e_summer$Var2 == i] = sp[i]
  
  interactions_on_c_winter$Var1[interactions_on_c_winter$Var1 == i] = sp[i]
  interactions_on_c_winter$Var2[interactions_on_c_winter$Var2 == i] = sp[i]
  
  interactions_on_e_winter$Var1[interactions_on_e_winter$Var1 == i] = sp[i]
  interactions_on_e_winter$Var2[interactions_on_e_winter$Var2 == i] = sp[i]
}

sets = unique(interactions_on_c_summer[,c(1,2)])
interactions_on_c_summer$labels = NA
interactions_on_e_summer$labels = NA
interactions_on_c_winter$labels = NA
interactions_on_e_winter$labels = NA

for(i in 1:dim(sets)[1]){
  interactions_on_c_summer[interactions_on_c_summer$Var1 == sets[i,'Var1'] & interactions_on_c_summer$Var2 ==  sets[i,'Var2'] & interactions_on_c_summer$Var3 == 1, 'labels'] = paste(trunc(interactions_on_c_summer[interactions_on_c_summer$Var1 == sets[i,'Var1'] & interactions_on_c_summer$Var2 ==  sets[i,'Var2'] & interactions_on_c_summer$Var3 == 1, 'value'], prec = 2),
                                                                                                                                                                                      '\n (',
                                                                                                                                                                                      trunc(interactions_on_c_summer[interactions_on_c_summer$Var1 == sets[i,'Var1'] & interactions_on_c_summer$Var2 ==  sets[i,'Var2'] & interactions_on_c_summer$Var3 == 2, 'value'], prec = 2),
                                                                                                                                                                                      ' ; ',
                                                                                                                                                                                      trunc(interactions_on_c_summer[interactions_on_c_summer$Var1 == sets[i,'Var1'] & interactions_on_c_summer$Var2 ==  sets[i,'Var2'] & interactions_on_c_summer$Var3 == 3, 'value'], prec = 2),
                                                                                                                                                                                      ")", sep='')
  
  interactions_on_e_summer[interactions_on_e_summer$Var1 == sets[i,'Var1'] & interactions_on_e_summer$Var2 ==  sets[i,'Var2'] & interactions_on_e_summer$Var3 == 1, 'labels'] = paste(trunc(interactions_on_e_summer[interactions_on_e_summer$Var1 == sets[i,'Var1'] & interactions_on_e_summer$Var2 ==  sets[i,'Var2'] & interactions_on_e_summer$Var3 == 1, 'value'], prec = 2),
                                                                                                                                                                                      '\n (',
                                                                                                                                                                                      trunc(interactions_on_e_summer[interactions_on_e_summer$Var1 == sets[i,'Var1'] & interactions_on_e_summer$Var2 ==  sets[i,'Var2'] & interactions_on_e_summer$Var3 == 2, 'value'], prec = 2),
                                                                                                                                                                                      ' ; ',
                                                                                                                                                                                      trunc(interactions_on_e_summer[interactions_on_e_summer$Var1 == sets[i,'Var1'] & interactions_on_e_summer$Var2 ==  sets[i,'Var2'] & interactions_on_e_summer$Var3 == 3, 'value'], prec = 2),
                                                                                                                                                                                      ")", sep='')
  
  interactions_on_c_winter[interactions_on_c_winter$Var1 == sets[i,'Var1'] & interactions_on_c_winter$Var2 ==  sets[i,'Var2'] & interactions_on_c_winter$Var3 == 1, 'labels'] = paste(trunc(interactions_on_c_winter[interactions_on_c_winter$Var1 == sets[i,'Var1'] & interactions_on_c_winter$Var2 ==  sets[i,'Var2'] & interactions_on_c_winter$Var3 == 1, 'value'], prec = 2),
                                                                                                                                                                                      '\n (',
                                                                                                                                                                                      trunc(interactions_on_c_winter[interactions_on_c_winter$Var1 == sets[i,'Var1'] & interactions_on_c_winter$Var2 ==  sets[i,'Var2'] & interactions_on_c_winter$Var3 == 2, 'value'], prec = 2),
                                                                                                                                                                                      ' ; ',
                                                                                                                                                                                      trunc(interactions_on_c_winter[interactions_on_c_winter$Var1 == sets[i,'Var1'] & interactions_on_c_winter$Var2 ==  sets[i,'Var2'] & interactions_on_c_winter$Var3 == 3, 'value'], prec = 2),
                                                                                                                                                                                      ")", sep='')
  
  interactions_on_e_winter[interactions_on_e_winter$Var1 == sets[i,'Var1'] & interactions_on_e_winter$Var2 ==  sets[i,'Var2'] & interactions_on_e_winter$Var3 == 1, 'labels'] =   paste(trunc(interactions_on_e_winter[interactions_on_e_winter$Var1 == sets[i,'Var1'] & interactions_on_e_winter$Var2 ==  sets[i,'Var2'] & interactions_on_e_winter$Var3 == 1, 'value'], prec = 2),
                                                                                                                                                                                        '\n (',
                                                                                                                                                                                        trunc(interactions_on_e_winter[interactions_on_e_winter$Var1 == sets[i,'Var1'] & interactions_on_e_winter$Var2 ==  sets[i,'Var2'] & interactions_on_e_winter$Var3 == 2, 'value'], prec = 2),
                                                                                                                                                                                        ' ; ',
                                                                                                                                                                                        trunc(interactions_on_e_winter[interactions_on_e_winter$Var1 == sets[i,'Var1'] & interactions_on_e_winter$Var2 ==  sets[i,'Var2'] & interactions_on_e_winter$Var3 == 3, 'value'], prec = 2),
                                                                                                                                                                                        ")", sep='')
}

interactions_on_c_summer$value = ifelse(interactions_on_c_summer$value > 1, 1, -1)
interactions_on_c_summer = interactions_on_c_summer[interactions_on_c_summer$Var3 == 1 & interactions_on_c_summer$Var1 != interactions_on_c_summer$Var2,]
interactions_on_c_summer$value = c(0,0,1,0,0,0)


interaction_cs_plot = ggplot(interactions_on_c_summer, aes(Var2, Var1, fill = as.factor(value)))+
  geom_tile(color = "white") +
  theme_minimal()+ # minimal theme
  scale_fill_manual(values = c("#e55539","grey","#79a4e5"), limits = c("-1","0","1"))+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20))+
  scale_y_discrete(labels = rev(c("M","L","P")),limits = rev(c("D. magna","D. longispina","D. pulex")))+
  scale_x_discrete(labels = c("M","L","P"),limits = (c("D. magna","D. longispina","D. pulex")), position = "top")+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = labels), color = "black", size = 6)+
  labs(y="on ...",x="Effect of ...")

interactions_on_e_summer$value = ifelse(interactions_on_e_summer$value > 1, -1, 1)
interactions_on_e_summer = interactions_on_e_summer[interactions_on_e_summer$Var3 == 1 & interactions_on_e_summer$Var1 != interactions_on_e_summer$Var2,]
interactions_on_e_summer$value = c(0,-1,0,0,-1,1)

interaction_es_plot = ggplot(interactions_on_e_summer, aes(Var2, Var1, fill = as.factor(value)))+
  geom_tile(color = "white")+
  scale_fill_manual(values = c("#e55539","grey","#79a4e5"), limits = c("-1","0","1"))+
  theme_minimal()+ # minimal theme
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20))+
  coord_fixed()+
  scale_y_discrete(labels = rev(c("M","L","P")),limits = rev(c("D. magna","D. longispina","D. pulex")))+
  scale_x_discrete(labels = c("M","L","P"),limits = (c("D. magna","D. longispina","D. pulex")), position = "top")+
  geom_text(aes(Var2, Var1, label = labels), color = "black", size = 6)+
  labs(y="on ...",x="Effect of ...")

interactions_on_c_winter$value = ifelse(interactions_on_c_winter$value > 1, 1, -1)
interactions_on_c_winter = interactions_on_c_winter[interactions_on_c_winter$Var3 == 1 & interactions_on_c_winter$Var1 != interactions_on_c_winter$Var2,]
interactions_on_c_winter$value = c(1,1,1,0,1,0)

interaction_cw_plot = ggplot(interactions_on_c_winter, aes(Var2, Var1, fill = as.factor(value)))+
  geom_tile(color = "white")+
  scale_fill_manual(values = c("#e55539","grey","#79a4e5"), limits = c("-1","0","1"))+
  theme_minimal()+ # minimal theme
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20))+
  coord_fixed()+
  scale_y_discrete(labels = rev(c("M","L","P")),limits = rev(c("D. magna","D. longispina","D. pulex")))+
  scale_x_discrete(labels = c("M","L","P"),limits = (c("D. magna","D. longispina","D. pulex")), position = "top")+
  geom_text(aes(Var2, Var1, label = labels), color = "black", size = 6)+
  labs(y="on ...",x="Effect of ...")

interactions_on_e_winter$value = ifelse(interactions_on_e_winter$value > 1, -1, 1)
interactions_on_e_winter = interactions_on_e_winter[interactions_on_e_winter$Var3 == 1 & interactions_on_e_winter$Var1 != interactions_on_e_winter$Var2,]
interactions_on_e_winter$value = c(-1,0,-1,0,0,-1)

interaction_ew_plot = ggplot(interactions_on_e_winter, aes(Var2, Var1, fill = as.factor(value)))+
  geom_tile(color = "white")+
  scale_fill_manual(values = c("#e55539","grey","#79a4e5"), limits = c("-1","0","1"))+
  theme_minimal()+ # minimal theme
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),
        axis.text = element_text(size = 18), 
        axis.title = element_text(size = 20))+
  coord_fixed()+
  scale_y_discrete(labels = rev(c("M","L","P")),limits = rev(c("D. magna","D. longispina","D. pulex")))+
  scale_x_discrete(labels = c("M","L","P"),limits = (c("D. magna","D. longispina","D. pulex")), position = "top")+
  geom_text(aes(Var2, Var1, label = labels), color = "black", size = 6)+
  labs(y="on ...",x="Effect of ...")



pg_interactions = plot_grid(interaction_cs_plot,interaction_es_plot,interaction_cw_plot,interaction_ew_plot,
                            labels =  c("A", "B","C", "D"),
                            label_size = 20, vjust = 0.95, scale = 0.9)
pg_interactions
 
#### save Figure 4: interactions (filename = interactions) ----
ggsave(pg_interactions, filename = "./figures/interactions.jpeg", width = 14, height = 14, dpi = 600)
ggsave(pg_interactions, filename = "./figures/interactions_pdf.pdf", width = 14, height = 14, dpi = 600)
ggsave(pg_interactions, filename = "./figures/interactions_png.png", width = 14, height = 14, dpi = 600)


###.    -----

### Figure 3   -----

#### Figure 3: Extinctions and colonization  -----

upperLimit_col = .4

colpalette = c("#a1d99b","#00441b")

plot_prob_PCA_space = function(dt, tit, lim, l = colpalette[1], h = colpalette[2]){
  ggplot(dt, aes(V2,V3, color= 1-exp(-1*(V1))))+geom_point()+
    labs(title=tit)+
    scale_color_continuous(low = l, high = h, limits = lim)+
    # scale_color_gradient()+
    labs(color = "Prob.")+
    theme()+labs(y="PC2",x="PC1")+
    scale_x_reverse()+
    scale_y_reverse()
}


envt = as.matrix(readRDS("data/environment_PCA.RDS")[,1:2])

m = as.data.frame(cbind((exp(log(jsM[2,1]) + envt %*% jsM[c(7,8),1])), envt[,1],envt[,2]))
m_g1<-plot_prob_PCA_space(m, "Spring extinction", c(0,1))
e_sum = m

m = as.data.frame(cbind((exp(log(jsM[3,1]) + envt %*% jsM[c(9,10),1])), envt[,1],envt[,2]))
m_g2<-plot_prob_PCA_space(m, "Spring colonization", c(0,upperLimit_col), colpalette[2], colpalette[1])
c_sum = m

m = as.data.frame(cbind((exp(log(jsM[4,1]) + envt %*% jsM[c(11,12),1])), envt[,1],envt[,2]))
m_g3<-plot_prob_PCA_space(m, "Winter extinction", c(0,1))
e_win = m

m = as.data.frame(cbind((exp(log(jsM[5,1]) + envt %*% jsM[c(13,14),1])), envt[,1],envt[,2]))
m_g4<-plot_prob_PCA_space(m, "Winter colonization", c(0,upperLimit_col), colpalette[2], colpalette[1])
c_win = m

# plot_grid(g1,g2,g3,g4)

al = 1 - (e_win+e_sum)/(c_win+c_sum)
al[,2:3] <- e_win[,2:3]
al[al[,1] < 0,1] <- 0
alM = al
# al[al[,1] > 0,1] <- 1
gMa<-ggplot(al, aes(V2,V3, color= 1-(V1)))+geom_point(aes(alpha = 0.05+V1*0.4))+
  labs(title="Expected Occ. (D.magna)", y="PC2",x="PC1")+
  scale_color_continuous(limits=c(0,1))+
  scale_alpha_continuous(range=c(0.2,1))+
  theme(legend.position = "none")+
  scale_y_reverse()+
  scale_x_reverse()


alwin = 1-(e_win/c_win)
alwin[,2:3] = e_win[,2:3]
alwin[alwin[,1] < 0, 1] = 0
gMaWin = ggplot(alwin, aes(V2,V3, color= 1-(V1)))+geom_point(aes(alpha = 0.05+V1*0.4))+
  labs(title="Expected Occ. (D.magna)",y="PC2",x="PC1")+
  scale_color_continuous(limits=c(0,1))+
  scale_alpha_continuous(range=c(0.2,1))+
  theme(legend.position = "none")+
  scale_y_reverse()+
  scale_x_reverse()

alsum = 1-(e_sum/c_sum)
alsum[,2:3] = e_sum[,2:3]
alsum[alsum[,1] < 0, 1] = 0
gMaSum = ggplot(alsum, aes(V2,V3, color= 1-(V1)))+geom_point(aes(alpha = 0.05+V1*0.4))+
  labs(title="Expected Occ. (D.magna)", y="PC2",x="PC1")+
  scale_alpha_continuous(range=c(0.2,1))+
  scale_color_continuous(limits=c(0,1))+
  theme(legend.position = "none")+
  scale_y_reverse()+
  scale_x_reverse()

gMaPeriods = plot_grid(gMaWin, gMaSum, ncol = 2)



m = as.data.frame(cbind((exp(log(jsP[2,1]) + envt %*% jsP[c(7,8),1])), envt[,1],envt[,2]))
p_g1<-plot_prob_PCA_space(m, "Spring extinction", c(0,1))
e_sum = m

m = as.data.frame(cbind((exp(log(jsP[3,1]) + envt %*% jsP[c(9,10),1])), envt[,1],envt[,2]))
p_g2<-plot_prob_PCA_space(m, "Spring colonization", c(0,upperLimit_col),colpalette[2], colpalette[1])
c_sum = m

m = as.data.frame(cbind((exp(log(jsP[4,1]) + envt %*% jsP[c(11,12),1])), envt[,1],envt[,2]))
p_g3<-plot_prob_PCA_space(m, "Winter extinction", c(0,1))
e_win = m

m = as.data.frame(cbind((exp(log(jsP[5,1]) + envt %*% jsP[c(13,14),1])), envt[,1],envt[,2]))
p_g4<-plot_prob_PCA_space(m, "Winter colonization", c(0,upperLimit_col), colpalette[2], colpalette[1])
c_win = m

# plot_grid(g1,g2,g3,g4)

al = 1 - (e_win+e_sum)/(c_win+c_sum)
al[,2:3] <- e_win[,2:3]
al[al[,1] < 0,1] <- 0
alP = al
# al[al[,1] > 0,1] <- 1
gPu<-ggplot(al, aes(V2,V3, color=  1 - (V1)))+geom_point(aes(alpha = 0.05+V1*0.4))+
  labs(title="Expected Occ. (D.pulex)", y="PC2",x="PC1")+
  scale_alpha_continuous(range=c(0.2,1))+
  scale_color_continuous(limits=c(0,1))+
  theme(legend.position = "none")+
  scale_y_reverse()+
  scale_x_reverse()

alwin = 1-(e_win/c_win)
alwin[,2:3] = e_win[,2:3]
alwin[alwin[,1] < 0, 1] = 0
gPuWin = ggplot(alwin, aes(V2,V3, color= 1-(V1)))+geom_point(aes(alpha = 0.05+V1*0.4))+labs(title="Expected Occ. (D.pulex)")+
  labs(y="PC2",x="PC1")+scale_color_continuous(limits=c(0,1))+
  scale_alpha_continuous(range=c(0.3,1))+
  theme(legend.position = "none")+
  scale_y_reverse()+
  scale_x_reverse()

alsum = 1-(e_sum/c_sum)
alsum[,2:3] = e_sum[,2:3]
alsum[alsum[,1] < 0, 1] = 0
gPuSum = ggplot(alsum, aes(V2,V3, color= 1-(V1)))+geom_point(aes(alpha = 0.05+V1*0.4))+labs(title="Expected Occ. (D.pulex)")+
  labs(y="PC2",x="PC1")+scale_color_continuous(limits=c(0,1))+
  scale_alpha_continuous(range=c(0.3,1))+
  theme(legend.position = "none")+
  scale_y_reverse()+
  scale_x_reverse()

gPuPeriods = plot_grid(gPuWin, gPuSum, ncol = 2)



m = as.data.frame(cbind((exp(log(jsL[2,1]) + envt %*% jsL[c(7,8),1])), envt[,1],envt[,2]))
l_g1<-plot_prob_PCA_space(m, "Spring extinction", c(0,1))
e_sum = m

m = as.data.frame(cbind((exp(log(jsL[3,1]) + envt %*% jsL[c(9,10),1])), envt[,1],envt[,2]))
l_g2<-plot_prob_PCA_space(m, "Spring colonization", c(0,upperLimit_col),colpalette[2], colpalette[1])
c_sum = m

m = as.data.frame(cbind((exp(log(jsL[4,1]) + envt %*% jsL[c(11,12),1])), envt[,1],envt[,2]))
l_g3<-plot_prob_PCA_space(m, "Winter extinction", c(0,1))
e_win = m

m = as.data.frame(cbind((exp(log(jsL[5,1]) + envt %*% jsL[c(13,14),1])), envt[,1],envt[,2]))
l_g4<-plot_prob_PCA_space(m, "Winter colonization", c(0,upperLimit_col), colpalette[2], colpalette[1])
c_win = m

# plot_grid(g1,g2,g3,g4)

al = 1 - (e_win+e_sum)/(c_win+c_sum)
al[,2:3] <- e_win[,2:3]
al[al[,1] < 0,1] <- 0
alL = al
# al[al[,1] > 0,1] <- 1
gLo<-ggplot(al, aes(V2,V3, color= 1-(V1)))+geom_point(aes(alpha = 0.05+V1*0.4))+
  labs(title="Expected Occ. (D.longi)", y="PC2",x="PC1")+
  scale_color_continuous(limits=c(0,1))+
  scale_alpha_continuous(range=c(0.2,1))+
  scale_x_reverse()+
  scale_y_reverse()
lege <- get_legend(gLo)


alwin = 1-(e_win/c_win)
alwin[,2:3] = e_win[,2:3]
alwin[alwin[,1] < 0, 1] = 0
gLoWin = ggplot(alwin, aes(V2,V3, color= 1-(V1)))+geom_point(aes(alpha = 0.05+V1*0.4))+labs(title="Expected Occ. (D.longi)")+
  labs(y="PC2",x="PC1")+scale_color_continuous(limits=c(0,1))+
  scale_alpha_continuous(range=c(0.3,1))+
  theme(legend.position = "none")+
  scale_x_reverse()+
  scale_y_reverse()

alsum = 1-(e_sum/c_sum)
alsum[,2:3] = e_sum[,2:3]
alsum[alsum[,1] < 0, 1] = 0
gLoSum = ggplot(alsum, aes(V2,V3, color= 1-(V1)))+
  geom_point(aes(alpha = 0.05+V1*0.4))+labs(title="Expected Occ. (D.longi)")+
  labs(y="PC2",x="PC1")+
  scale_alpha_continuous(range=c(0.3,1))+
  scale_color_continuous(limits=c(0,1))+theme(legend.position = "none")+
  scale_x_reverse()+
  scale_y_reverse()

gLoPeriods = plot_grid(gLoWin, gLoSum, ncol = 2)

# Plot Ext

legend_ext = get_legend(m_g3)
ext_grid = plot_grid(m_g1+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title="",x = "", y = ""), 
                     m_g3+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title="",x = "", y = ""),
                     l_g1+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title="",x = "", y = ""),
                     l_g3+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title="",x = "", y = ""),
                     p_g1+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title="",x = "", y = ""),
                     p_g3+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title="",x = "", y = ""),
                     ncol = 2)#,
# labels = c("M","M","P","P","L","L") )



legend_col = get_legend(m_g2)
col_grid = plot_grid(m_g2+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title="",x = "", y = ""), 
                     m_g4+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title="",x = "", y = ""),
                     l_g2+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title="",x = "", y = ""),
                     l_g4+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title="",x = "", y = ""),
                     p_g2+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title="",x = "", y = ""),
                     p_g4+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title="",x = "", y = ""),
                     ncol = 2)#,

plot_grid(ext_grid, col_grid, ncol = 2)



#### Figure 3: Expectations & Observations -----

expec = plot_grid(gMa+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title = "", y = "", x = ""),
                  gLo+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title = "",y = "", x = ""),
                  gPu+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title = "",y = "", x = ""),
                  ncol = 1)


# Gather observations
envt = as.matrix(readRDS("data/environment_PCA.RDS")[,1:3])
magna = readRDS("data/occupancies_magna_82-17_core.RDS")
longi = readRDS("data/occupancies_longispina_82-17_core.RDS")
pulex = readRDS("data/occupancies_pulex_82-17_core.RDS")

magna = cbind(magna[[1]], magna[[2]])
longi = cbind(longi[[1]], longi[[2]])
pulex = cbind(pulex[[1]], pulex[[2]])

magna_1 = apply(magna, 1, mean, na.rm=T)
longi_1 = apply(longi, 1, mean, na.rm=T)
pulex_1 = apply(pulex, 1, mean, na.rm=T)

magna_envt = as.data.frame(cbind(magna_1, envt[,1], envt[,2]))
longi_envt = as.data.frame(cbind(longi_1, envt[,1], envt[,2]))
pulex_envt = as.data.frame(cbind(pulex_1, envt[,1], envt[,2]))

colnames(magna_envt)[1] = colnames(longi_envt)[1] = colnames(pulex_envt)[1] <- c("Obs") 

Observations_magna = ggplot(magna_envt, aes(V2,V3, color= 1-(Obs)))+
  geom_point(aes(alpha = 0.05+Obs*0.4))+labs(title="Expected Occ. (D.longi)")+
  scale_alpha_continuous(range=c(0.3,1))+
  scale_color_continuous(limits=c(0,1))+theme(legend.position = "none")+
  scale_x_reverse()+
  scale_y_reverse()

Observations_pulex = ggplot(pulex_envt, aes(V2,V3, color= 1-(Obs)))+
  geom_point(aes(alpha = 0.05+Obs*0.4))+labs(title="Expected Occ. (D.longi)")+
  scale_alpha_continuous(range=c(0.3,1))+
  scale_color_continuous(limits=c(0,1))+theme(legend.position = "none")+
  scale_x_reverse()+
  scale_y_reverse()

Observations_longi = ggplot(longi_envt, aes(V2,V3, color= 1-(Obs)))+
  geom_point(aes(alpha = 0.05+Obs*0.4))+labs(title="Expected Occ. (D.longi)")+
  scale_alpha_continuous(range=c(0.3,1))+
  scale_color_continuous(limits=c(0,1))+theme(legend.position = "none")+
  scale_x_reverse()+
  scale_y_reverse()

Observ = plot_grid(Observations_magna+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title = "", y = "", x = ""),
                   Observations_longi+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title = "",y = "", x = ""),
                   Observations_pulex+theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank())+labs(title = "",y = "", x = ""),
                   ncol = 1)

ext_col_expt = plot_grid(ext_grid, col_grid, expec, Observ, ncol = 4, rel_widths = c(2,2,1,1))
ext_col_expt

#### save Figure 3 (complete)   ----
ggsave(ext_col_expt, filename = "./figures/col_ext_expect.jpeg", dpi = 600, width = 14.4)
ggsave(ext_col_expt, filename = "./figures/col_ext_expect.pdf", dpi = 600, width = 14.4)

# plot_grid(gMaPeriods, gLoPeriods, gPuPeriods, ncol = 1)
#### not used  


###.    ----
### Expected & observed occupancies by season (2 PC's) ------

magna = readRDS("data/occupancies_magna_82-17_core.RDS")
magna_1 = apply(magna[[1]], 1, mean, na.rm=T)
magna_2 = apply(magna[[2]], 1, mean, na.rm=T)
magna_envt = as.data.frame(cbind(magna_1, magna_2, envt[,1], envt[,2]))
colnames(magna_envt)[1:2] <- c("ObsSumm","ObsWin")

gMaObsSumm<-ggplot(magna_envt, aes(V3,V4, color= (ObsSumm)))+
  geom_point(size = 1, alpha = 0.1+(magna_envt$ObsSumm/max(magna_envt$ObsSumm)*0.9))+
  labs(title="Obs Occ. (Spring)", y="PC2",x="PC1")+
  scale_color_continuous(limits=c(0,1))+
  theme(legend.position = "none")+
  scale_y_reverse()+
  scale_x_reverse()

gMaObsWin<-ggplot(magna_envt, aes(V3,V4, color= (ObsWin)))+
  geom_point(size = 1, alpha = 0.1+(magna_envt$ObsWin/max(magna_envt$ObsWin)*0.9))+
  labs(title="Obs Occ. (Winter)",y="PC2",x="PC1")+
  scale_color_continuous(limits=c(0,1))+
  theme(legend.position = "none")+
  scale_y_reverse()+
  scale_x_reverse()
# 

pulex = readRDS("data/occupancies_pulex_82-17_core.RDS")
pulex_1 = apply(pulex[[1]], 1, mean, na.rm=T)
pulex_2 = apply(pulex[[2]], 1, mean, na.rm=T)
pulex_envt = as.data.frame(cbind(pulex_1, pulex_2, envt[,1], envt[,2]))
colnames(pulex_envt)[1:2] <- c("ObsSumm","ObsWin") 

gPuObsSumm<-ggplot(pulex_envt, aes(V3,V4, color= (ObsSumm)))+
  geom_point(size = 1, alpha = 0.1+(pulex_envt$ObsSumm/max(pulex_envt$ObsSumm)*0.9))+
  labs(title="Obs Occ. (Spring)",y="PC2",x="PC1")+
  scale_color_continuous(limits=c(0,1))+
  theme(legend.position = "none")+
  scale_y_reverse()+
  scale_x_reverse()

gPuObsWin<-ggplot(pulex_envt, aes(V3,V4, color= (ObsWin)))+
  geom_point(size = 1, alpha = 0.1+(pulex_envt$ObsWin/max(pulex_envt$ObsWin)*0.9))+
  labs(title="Obs Occ. (Winter)",y="PC2",x="PC1")+
  scale_color_continuous(limits=c(0,1))+
  theme(legend.position = "none")+
  scale_y_reverse()+
  scale_x_reverse()

longi = readRDS("data/occupancies_longispina_82-17_core.RDS")
longi_1 = apply(longi[[1]], 1, mean, na.rm=T)
longi_2 = apply(longi[[2]], 1, mean, na.rm=T)
longi_envt = as.data.frame(cbind(longi_1, longi_2, envt[,1], envt[,2]))
colnames(longi_envt)[1:2] <- c("ObsSumm","ObsWin") 

gLoObsSumm<-ggplot(longi_envt, aes(V3,V4, color= (ObsSumm)))+
  geom_point(size = 1, alpha = 0.1+(longi_envt$ObsSumm/max(longi_envt$ObsSumm)*0.9))+
  labs(title="Obs Occ. (Spring)",y="PC2",x="PC1")+
  scale_color_continuous(limits=c(0,1))+
  theme(legend.position = "none")+
  scale_y_reverse()+
  scale_x_reverse()

gLoObsWin<-ggplot(longi_envt, aes(V3,V4, color= (ObsWin)))+
  geom_point(size = 1, alpha = 0.1+(longi_envt$ObsWin/max(longi_envt$ObsWin)*0.9))+
  labs(title="Obs Occ. (Winter)",y="PC2",x="PC1")+
  scale_color_continuous(limits=c(0,1))+
  theme(legend.position = "none")+
  scale_y_reverse()+
  scale_x_reverse()

plot_grid(gMa, gMaObsSumm,gMaObsWin,
          gLo+theme(legend.position = "none"), gLoObsSumm,gLoObsWin,
          gPu, gPuObsSumm, gPuObsWin,
          ncol = 3)






