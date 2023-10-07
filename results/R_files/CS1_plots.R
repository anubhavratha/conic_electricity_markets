library(data.table)
library(ggplot2)
library(MASS)
library(RColorBrewer)
library(dbplyr)
library(tikzDevice)

BoxColors = c("#076256","#fc452d","#90a4ae")
savepdf <- function(file, width=17.55, height=12.08)
{
  fname <- paste("../plots/rplots/",file,".pdf",sep="")
  pdf(fname, width=width/2.54, height=height/2.54,
      pointsize=10)
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
}


#Getting cost 
cost_data = fread("../cost_vs_wcap_100scen.csv")
tikzDevice::tikz("../plots/rplots/cost_vs_wcap_100scens.tex", width=6, height=3)
#savepdf("cost_plot")
p_cost <- ggplot(data=cost_data, aes(x = as.factor(ceiling(resShare)), y = costTot, fill=as.factor(model))) + 
  geom_bar(stat = 'identity', position = position_dodge())+
  facet_wrap(~congestion)+
  scale_fill_manual(values=BoxColors, labels=c("Mcc", "R1", "R2"))+
  xlab("RES Share (pc of Peak Load)")+
  ylab("In-sample Expected Cost (kUSD)")+
  guides(fill=guide_legend(ncol=3))+
  labs(fill="Category")+
  theme(
    legend.position="top",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color="grey80",linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    strip.background = element_rect(color="black",fill=NA, linetype="solid")
  )
p_cost
#print(p_cost)
dev.off()

savepdf("cost_compare_stochastic")
cost_data_subset = subset(cost_data,model!="R1")
p_cost_sto <- ggplot(data=cost_data_subset, aes(x = resShare, y = costTot, fill=as.factor(model))) + 
  geom_bar(stat = 'identity', position = position_dodge())+
  facet_wrap(~congestion)+
  scale_fill_manual(values=BoxColors, labels=c("Mcc", "R2"))+
  xlab("RES Share (% of Peak Load)")+
  ylab("In-sample Expected Cost (kUSD)")+
  guides(fill=guide_legend(ncol=2))+
  #  labs(fill="Category")+
  theme(
    legend.position="top",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color="grey80",linetype = "dashed"),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    strip.background = element_rect(color="black",fill=NA, linetype="solid")
  )
print(p_cost_sto)
dev.off()

### OUT OF SAMPLE COST COMPARISON With real-time penalty = 10%
#Getting results table
rawdata  = fread("../new_oos_cost_comparison_500scens.csv")
tikzDevice::tikz("../plots/rplots/oos_cost_comparison.tex", width=6, height=3)
#Boxplot for cost distribution
plot_costcomp <- ggplot(data=rawdata, aes(x=factor(networkconfig), y=oos_cost, color=factor(model)))+
  geom_boxplot(outlier.shape = 1, outlier.size=3, lwd=1)+
  #geom_point(position=position_jitterdodge(), alpha=0.3, size=0.7)+
  theme_bw(base_size=16)+
 # facet_wrap(~networkconfig)+
  labs(color="Model")+
  scale_color_manual(values=BoxColors)+
  xlab("Network Configuration")+
  ylab("Out-of-sample expected cost")+
#  ylim(200000,500000)+
  theme(
#    legend.position="top",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major.y = element_line(color="grey80",linetype = "dashed"),
    panel.grid.major.x = element_line(color="grey80",linetype = "dashed"),
    panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    strip.background = element_rect(color="black",fill=NA, linetype="solid")
  )
plot_costcomp

dev.off()


### OUT OF SAMPLE COST COMPARISON With Different Values of real-time penalties
#Getting results table
rawdata  = fread("../oos_costdata_with_rtpenvals.csv")
tikzDevice::tikz("../plots/rplots/oos_cost_comparison.tex", width=6, height=3)
#Boxplot for cost distribution
plot_costcomp <- ggplot(data=rawdata, aes(x=factor(rt_pen), y=oos_cost, color=factor(model)))+
  geom_boxplot(outlier.shape = NA )+
  geom_point(position=position_jitterdodge(), alpha=0.4, size=0.5)+
  theme_bw(base_size=16)+
  facet_wrap(~networkconfig)
labs(color="Model")+ 
  xlab("Network Configuration")+
  ylab("Out-of-sample expected cost")+
  scale_y_continuous(labels = scales::scientific)+
  ylim(150000,500000)
plot_costcomp

dev.off()


### Prices plot on the network
#Getting node and coordinates data - also price values
fdata=fread("out_data_24node.csv")

# #Getting the nodal points on the coordinate axis
ggplot(data=fdata, aes(x=x,y=y))+
  geom_point(color="BLUE")

# #Plot for node sizes depending on variance + transparent graph
ggplot(data=fdata)+
  geom_point(mapping=aes(x=x,y=y,size=el_price), colour="BLUE")+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA))


#Inputing kernel data for the KDE
price_kerdata = fread("price_kernel_plot.csv")

### --- PLOTS FOR THE UNPENALIZED CASE --- ###
#Exploring data values for sanity check
#ggplot(data=unpen_data, aes(x=fgval))+
# geom_histogram()+
# scale_x_log10()


savepdf <- function(file, width=18, height=27)
{
  fname <- paste("Rplots/",file,".pdf",sep="")
  pdf(fname, width=width/2.54, height=height/2.54,
      pointsize=10)
  par(mgp=c(2.2,0.45,0), tcl=-0.4, mar=c(3.3,3.6,1.1,1.1))
}

#compute normalized legend entries
max_kde = round(max(price_kerdata$fgval), digits=2)
min_kde = round(min(price_kerdata$fgval), digits=2)
max_price = round(max(fdata$el_price), digits=2)
min_price = round(min(fdata$el_price), digits=2)
norm_fact = (max_price)/max_kde

savepdf("price_kernel_plot")
p1 <- ggplot(data=price_kerdata)+
  geom_contour_filled(aes(x=xval,y=yval, z=fgval), show.legend = FALSE, breaks=pretty(0:max_kde+1, n=6))+
  scale_fill_brewer(palette='Oranges', guide = guide_colorsteps(
    frame.colour = "black",
    ticks.colour = "black", # you can also remove the ticks with NA
    barwidth=1,
    barheight=8,
    show.limits = FALSE),
    name = "EL Price",
    labels = round(norm_fact*pretty(0:max_kde+1, n=6)[2:5]))+
    theme_void()+
    theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.margin=grid::unit(c(0,0,0,0), "mm"),
    legend.position = c(0.3,0.25)
  )+
  xlim(c(0,18))+
  ylim(c(3,30))
print(p1)
dev.off()


tikzDevice::tikz("Rplots/legend_price_ker.tex", width=6, height=3)
#savepdf("legend_price_kernel_plot")
p2 <- ggplot(data=price_kerdata)+
  geom_contour_filled(aes(x=xval,y=yval, z=fgval), show.legend = TRUE, breaks=pretty(0:max_kde+1, n=6))+
  scale_fill_brewer(palette='Oranges', guide = guide_colorsteps(
    frame.colour = "black",
    ticks.colour = "black", # you can also remove the ticks with NA
    barwidth=1,
    barheight=8,
    show.limits = FALSE),
    name = "EL Price",
    labels = round(norm_fact*pretty(0:max_kde+1, n=6)[2:5]))+
  theme_void()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    plot.margin=grid::unit(c(0,0,0,0), "mm"),
    legend.position = "right"
  )+
  xlim(c(0,18))+
  ylim(c(3,30))
p2
dev.off()


