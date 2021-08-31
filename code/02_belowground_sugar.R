# Code files for In Situ Seagrass Wasting Disease manuscript
# 02_belowground_sugar

# Last updated 2021-08-30 by Lillian Aoki

# This script compares disease and belowground sugar content from the in situ marking experiment, 
# conducted in 2019 by Olivia Graham, Lillian Aoki, Tiffany Stephens, Joshua Stokes

# Outputs include Fig 4 in the manuscript

# Load libraries ####
library(tidyverse)
library(here)

# Read in data ####
df.all <- read.csv("data/BelowgroundDiseaseSugar.csv")

df.density <- read.csv("data/SeagrassDensity.csv")
df.density.summ <- df.density%>%
  group_by(site, transect_no)%>%
  summarize(shoots=mean(shoots.m2))

#df.all <- left_join(df.sugar, df.disease, by = c("site", "plant_id"))
df.all$site <- str_trim(df.all$site)
df.canopy <- read.csv("data/CanopyHeight.csv")
df.final <- left_join(df.all,df.canopy,by=c("site","transect_no","plant_id"))
df.final <- left_join(df.final,df.density.summ,by=c("site","transect_no"))

# Build belowground model ####
M1 <- lm(sugar_mgg~severity+shoots+blade_length_cm+site+severity:site,
         data=df.final)
E.M1 <- resid(M1,type="response")
F.M1 <- fitted(M1)
plot(E.M1~F.M1)
hist(E.M1,breaks = 10)
qqnorm(E.M1)
qqline(E.M1)
## residuals look good
# plot against model covariates as well
plot(E.M1~df.final$shoots)
plot(E.M1~df.final$severity)
plot(E.M1~df.final$blade_length_cm)
plot(E.M1~as.factor(df.final$site))
# still look good
# need to check for multicollinearity
car::vif(M1)
# No issues with VIF
# So, acceptable to include all the terms in the model.
# Model terms - look at the significance of individual model terms
drop1(M1,test="Chisq")
# AIC doesn't improve by dropping the 3-way interaction or site and the model changes significantly
# Therefore we keep all terms in the model.
# Assess significance of terms
summary(M1)

# Visualize sugar model ####
newdat <- with(df.final,seq(from=min(severity), to = max(severity), length.out = 100))
df.new <- data.frame(severity=newdat,blade_length_cm=median(df.final$blade_length_cm),
                     shoots=median(df.final$shoots),site=rep(unique(df.final$site),each=100))
pred <- predict(M1,newdata=df.new,se.fit=TRUE,level=0.05)
df.pred <- data.frame(predicted=pred$fit, predicted.se=pred$se.fit)
df.pred <- cbind(df.pred,df.new)

ggplot(df.pred,aes(x=severity,y=predicted,color=site))+
  #  geom_ribbon(aes(ymax=predicted+predicted.se*1.96,ymin=predicted-predicted.se*1.96,fill=site,color=site),alpha=0.45)+
  geom_line(aes(color=site),size=1,show.legend = FALSE)+
  geom_point(data=df.final,aes(x=severity,y=sugar_mgg,shape=site), size=3)+
  #facet_wrap(~site,scales = "free_x")+
  xlab("Disease severity (% diseased blade area)")+
  ylab("Sugar content \n(mg per g rhizome DW)")+
  #scale_color_viridis_d()+
  labs(color="")+
  #  scale_color_manual(values=c("#4A7BB7","#EAECCC","#DD3D2D"))+
  scale_color_manual(values=c("#98CAE1","#EAECCC","#FDB366"),name="",
                     labels=c("4th of July Beach","False Bay","Indian Cove"))+
  scale_shape_discrete(name="",labels=c("4th of July Beach","False Bay","Indian Cove"))+
  theme_bw()+
  theme(strip.background=element_blank(),
        panel.grid = element_blank(),
        axis.title = element_text(size=10),
        axis.text=element_text(size=9),
        legend.text=element_text(size=9),
        legend.title=element_blank(),
        legend.position = c(0.75,0.85),
        legend.key.size = unit(0.25,"cm"),
        strip.text=element_text(size=9),
        #legend.background = element_rect(color="grey50"),
        legend.margin = margin(t=0,r=1,b=2,l=1,unit="pt"),
        plot.margin = margin(t=1,r=1,b=1,l=1,unit="pt"))
#ggsave(here("figures/high_res","Figure3.tiff"),width = 3.3,height=2)
ggsave(here("figures","Figure4.jpg"),width = 3.3,height=2)
