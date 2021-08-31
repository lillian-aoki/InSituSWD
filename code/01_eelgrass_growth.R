# Code files for In Situ Seagrass Wasting Disease manuscript
# 01_specific_productivity_model

# Last updated 2021-08-30 by Lillian Aoki

# This script imports, models, and visualizes data from the eelgrass wasting disease in situ marking experiment,
# conducted in 2019 by OJG and LRA

# Outputs include Fig 1, 2, 3, S1, S2 in the manuscript

# Load libraries ####
library(nlme)
library(lme4)
library(glmmTMB)
library(tidyverse)
library(DHARMa)
library(sjPlot)
library(ggeffects)
library(quantreg)
library(broom)
library(car)
library(here)

# Read in data ####
blades <- read.csv("data/BladesSpProdAnalysis_NSFEeLISA_correct_split.csv")
blades <- na.omit(blades)
# use only leaf 1-5 (few shoots and more than 5 leaves and they don't grow much at all)
blades <- subset(blades,Leaf<6)
blades$LeafRank <- as.factor(blades$Leaf)
blades$Prevalence <- as.factor(blades$Prevalence)
# need a shootID 
blades$ShootId <- paste(blades$Month, blades$Transect, blades$Shoot, sep="_")
# If we want to include Blade Area as a continuous predictor, need to center and scale it
blades$sBladeArea <- scale(blades$BladeArea,center = TRUE,scale=TRUE)
# Exclude blades with 0 growth that are the oldest blade
try <- blades[blades$Leaf==blades$TotalLeaves & blades$SpProdMm2==0,]
blades_2 <- anti_join(blades,try)
# Drop additional zeros (42 total)
blades_drop <- blades_2[blades_2$SpProdMm2>0,]

# Build prevalence model ####
# Run a model with Prevalence, LeafRank, Month, and BladeArea (scaled) as fixed effects
# ShootId as a random effect
glmm_final <- glmmTMB(SpProdMm2~Prevalence+LeafRank+Month+sBladeArea+(1|ShootId),
                      family = Gamma(link="log"),
                      data=blades_drop,
                      dispformula = ~LeafRank)
# Simulate residuals and check validate model
E.sim <- simulateResiduals(glmm_final)
plot(E.sim)
# Residuals look pretty good - no dispersion issues! Still have outliers but they are real data, so don't drop.
# The K-S deviation (non-normality) is also not that surprising with such a large sample size. 
# Plot against model covariates
plot(E.sim$scaledResiduals~blades_drop$sBladeArea)
plot(E.sim$scaledResiduals~blades_drop$Prevalence)
plot(E.sim$scaledResiduals~blades_drop$LeafRank)
plot(E.sim$scaledResiduals~as.factor(blades_drop$Month))
# Resids against Leaf Rank still has a slight pattern, but it's really hard to eliminate that (partly due to the unbalanced samples)
# This is an acceptable model.
summary(glmm_final)
drop1(glmm_final, test = "Chisq")
# Prevalence and blade rank are significant, interaction is not.

# Visualize the prevalence model ####
# Note, using sjPlot() package, to get the model terms to plot in the order we want, it's easiest to refit the model
glmm_plot <- glmmTMB(SpProdMm2~LeafRank+Prevalence+Month+sBladeArea+(1|ShootId),
                     family = Gamma(link="log"),
                     data=blades_drop,
                     dispformula = ~LeafRank)
term_names <- c("Disease \npresence","Leaf rank \n(2nd youngest)","Leaf rank \n(3rd youngest)",
                "Leaf rank \n(4th youngest)", "Leaf rank \n(5th youngest)")
p1 <- plot_model(glmm_plot,
                 terms = c("Prevalence [1]","LeafRank [2,3,4,5]"),
                 order.terms = c(5,1,2,3,4),
                 type="std",
                 axis.labels = rev(term_names),
                 title="",
                 show.p = TRUE,
                 show.values = TRUE,
                 value.offset = 0.2,
                 digits = 3,
                 colors="black"
)
p1_e <- p1 + scale_y_log10(limits=c(0.001,1))+
  theme_bw()
p1_e
ggsave(here("figures","FigureS1.jpg"),width=6,height=4)

pred <- ggpredict(glmm_final,terms=c("Prevalence","LeafRank"))
pred$Prev <- pred$x
ggplot(pred,aes(x=group,y=predicted))+
  geom_jitter(data=blades_drop,aes(x=LeafRank,y=SpProdMm2,color=Prevalence),#width=0.3,
              alpha=0.25,position=position_dodge(width = 0.5),show.legend = FALSE,size=0.75)+
  geom_errorbar(aes(ymin=conf.low,ymax=conf.high,group=x),width=0.15,size=0.75, show.legend = FALSE,
                position=position_dodge(width = 0.5),color="grey20")+
  geom_point(shape=21,position=position_dodge(width = 0.5),size=1.5,aes(fill=x),color="grey20",stroke=0.75)+
  scale_y_log10(limits=c(0.0004,6))+
  scale_fill_manual(values=c("#98CAE1","#FDB366"),labels=c("Healthy","Diseased"))+
  scale_x_discrete(labels=c("1st\n(youngest)","2nd","3rd","4th","5th\n(oldest)"))+
  scale_color_manual(values=c("#98CAE1","#FDB366"))+
  #scale_color_manual(values=c("grey20","grey20"))+
  #scale_color_manual(values=c("seagreen","wheat4"),labels=c("Healthy","Diseased"))+
  #labs(fill="")+
  xlab("Blade Rank")+
  ylab("Specific productivty \n(% new leaf area per day)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(size=9),
        axis.title = element_text(size=10),
        #   legend.key = element_rect(size=4),
        legend.key.size = unit(0.25,"cm"),
        legend.text = element_text(size=8),
        legend.title=element_blank(),   
        #     legend.title = element_text(size=9),
        legend.position = c(0.85,0.85),
        legend.spacing = unit(2,"pt"),
        #legend.background = element_rect(color="black"),
        legend.margin = margin(t=1,r=1,b=1,l=1,unit="pt"),
        plot.margin = margin(t=1,r=1,b=1,l=1,unit="pt"))
ggsave(here("figures","Figure1.jpg"),width=3.3,height=2)

# Build severity model ####
# use the same data set as before, but limit to leaf 3 to capture range of severities
leaf3 <- subset(blades,Leaf==3)
# Set the quantiles of interest
qs3 <- c(0.10,0.25,0.50,0.75,0.9)
# run the quantile regression
qr1 <- rq(Mm2LD~Severity,data=leaf3,tau=qs3)
# Look at model output
summary(qr1)

# Visualize the severity model ####
# For first plot, need to create a new data frame with slopes (coefficient estimates) and intercepts from the model
coef1 <- qr1$coefficients
rq_coef1 <- data.frame(intercept = coef1[1, ], coef = coef1[2, ],
                       tau_flag =colnames(coef1))
rq_coef1$xend <- -rq_coef1$intercept/rq_coef1$coef
rq_coef1$yend <- rq_coef1$intercept+rq_coef1$coef*rq_coef1$xend


SevQ <- ggplot(leaf3,aes(x=Severity,y=Mm2LD))+
  geom_point(size=1.75,color="grey40")+
  geom_segment(data=rq_coef1,aes(x=0,y=intercept, xend=xend, yend=yend,color=tau_flag),size=1)+
  scale_color_manual(values=c("#4A7BB7","#98CAE1","#EAECCC","#FDB366","#DD3D2D"),
                     labels=c("10th","25th","50th","75th","90th"))+
  #  scale_color_viridis_d(labels=c("10th","25th","50th","75th","90th"),direction = -1)+
  xlab("Disease severity (% diseased blade area)")+
  ylab(expression("Blade growth (mm"^"2"~"d"^"-1"~")"))+
  labs(col="Quantile")+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),
        panel.grid = element_blank(),
        legend.position = c(0.85,0.7),
        legend.title = element_text(size=9),
        legend.key.size = unit(0.25,"cm"),
        strip.text=element_text(size=9),
        axis.text=element_text(size=9),
        axis.title=element_text(size=10),
        legend.margin = margin(t=1,r=1,b=1,l=1,unit="pt"),
        plot.margin = margin(t=1,r=1,b=1,l=1,unit="pt"))
SevQ
ggsave(here("figures","Figure2.jpg"),width=3.3,height=2)

# Second plot, re-calculated the regressions for additional quantiles (0.02 intervals)
# This is just to get a smoother line (i.e. 50 points instead of 5 points) 
# But we can calculate however we want
qs4 <- 1:49/50
qr2 <- rq(Mm2LD~Severity,data=leaf3,tau=qs4)
# Bootstrapped the confidence intervals for the coefficient estimates
dat <- broom::tidy(qr2,se.type="boot",conf.int=TRUE,conf.level=0.95) %>%
  filter(!grepl("(Intercept)", term))

ggplot(dat,aes(x=tau,y=estimate))+
  #geom_point(color="grey40", size = 3)+ 
  geom_line(color="grey40", size = 1)+ 
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=0.25, fill="grey50")+
  geom_hline(yintercept = 0,linetype="dashed")+
  xlab("Quantile")+
  ylab("Severity coefficient estimate")+
  theme_bw()+
  theme(panel.grid = element_blank())
ggsave(here("figures","FigureS2.jpg"), width=6, height=4)

## Compare lesion and leaf growth rates ####
lesion <- read.csv("data/LesionAreaChangeNSF.csv")
month_order <- c("June","July")
lesion_leaf <- left_join(lesion,blades,by=c("ShootId","Month","Transect","Shoot","Subshoot","LeafNum"="Leaf"))
lesion_leaf$lesion_mm2_d <- lesion_leaf$DelLesionArea/lesion_leaf$Days
lesion_leaf$Month <- ordered(lesion_leaf$Month,levels=month_order)

t.test(lesion_leaf$lesion_mm2_d,lesion_leaf$Mm2LD,paired=TRUE)

# Visualize lesion and leaf rates ####
p1 <- ggplot(lesion_leaf,aes(x=lesion_mm2_d,y=Mm2LD))+geom_point(size=1.75, color="grey40")+
  geom_abline(slope=1,intercept=0,linetype="dashed")+
  #scale_color_manual(values = c("steelblue","darkgoldenrod"))+
  xlab(expression("Lesion growth (mm"^"2"~"d"^"-1"~")"))+
  ylab(expression("Blade growth (mm"^"2"~"d"^"-1"~")"))+
  theme_bw()+
  theme(axis.text = element_text(size=9),
        axis.title = element_text(size=10),
        panel.grid = element_blank(),
        plot.margin = margin(t=1,r=1,b=1,l=1,unit="pt"))
p1
ggsave(here("figures","Figure3.jpg"),width=3.3, height=2)
