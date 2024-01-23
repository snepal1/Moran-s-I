
rm(list=ls())

#### Moran's I ####

library(rgdal)
library(spdep)
library(dplyr)

# correlogram function
icorrelogram <- function(locations, z, binsize, maxdist){
  distbin <- seq(0,maxdist,by=binsize)
  Nbin <- length(distbin)-1
  moran.results <- data.frame("dist"= rep(NA, Nbin), "Morans.i"=NA,"null.lower"=NA, "null.upper"=NA)
  
  for (i in 1:Nbin){
    d.start<-distbin[i] 
    d.end<-distbin[i+1]
    neigh <- dnearneigh(x=locations, d1=d.start, d2=d.end, longlat=F)
    wts <- nb2listw(neighbours=neigh, style='B', zero.policy=T)
    mor.i <- moran.mc(x=z, listw=wts, nsim=1000, alternative="greater", zero.policy=T, na.action = na.omit)  # note alternative is for P-value, so only 'significant if positive autocorrelation
    
    moran.results[i, "probability"]<-mor.i$p.value
    moran.results[i, "dist"]<-(d.end+d.start)/2 
    moran.results[i, "Morans.i"]<-mor.i$statistic 								                #observed moran's i
    moran.results[i, "null.lower"]<-quantile(mor.i$res, probs = 0.025,na.rm = T)  #95% null envelope	
    moran.results[i, "null.upper"]<-quantile(mor.i$res, probs = 0.975,na.rm = T)#95% null envelope
  }
  return(moran.results)
}



# read in data
#filename_df <- readOGR("D:/data_groundtruth", layer = 'biomass_both_point')

filename_df <- readOGR("S:\\Biomass\\Biomass_data", layer = 'Both_biomass')

filename<-subset(filename_df,UNIT==39)

fig_dataRA<-list(cbind(filename@data, filename@coords)%>% 
                   dplyr::select(coords.x1,coords.x2,Biomass_20,Biomass_19))

coords <- as.matrix(data.frame('x' = fig_dataRA[[1]]$coords.x1, 'y' = fig_dataRA[[1]]$coords.x2))

filename1<-subset(filename_df,UNIT==40)

fig_dataRB<-list(cbind(filename1@data, filename1@coords)%>% 
                   dplyr::select(coords.x1,coords.x2,Biomass_20,Biomass_19))

coords1 <- as.matrix(data.frame('x' = fig_dataRB[[1]]$coords.x1, 'y' = fig_dataRB[[1]]$coords.x2))


filename2<-subset(filename_df,UNIT==43)

fig_dataRC<-list(cbind(filename2@data, filename2@coords)%>% 
                   dplyr::select(coords.x1,coords.x2,Biomass_20,Biomass_19))

coords2 <- as.matrix(data.frame('x' = fig_dataRC[[1]]$coords.x1, 'y' = fig_dataRC[[1]]$coords.x2))


filename3<-subset(filename_df,UNIT==44)

fig_dataRD<-list(cbind(filename3@data, filename3@coords)%>% 
                   dplyr::select(coords.x1,coords.x2,Biomass_20,Biomass_19))

coords3 <- as.matrix(data.frame('x' = fig_dataRD[[1]]$coords.x1, 'y' = fig_dataRD[[1]]$coords.x2))





filename4<-subset(filename_df,UNIT==45)

fig_dataRE<-list(cbind(filename4@data, filename4@coords)%>% 
                   dplyr::select(coords.x1,coords.x2,Biomass_20,Biomass_19))

coords4 <- as.matrix(data.frame('x' = fig_dataRE[[1]]$coords.x1, 'y' = fig_dataRE[[1]]$coords.x2))




# calculate moran's I and simulated random confidence interval for all response variables
################## RNA_A
final_out <- data.frame()
for(i in 1:length(fig_dataRA)){
  year <- names(fig_dataRA[i])
  temp_out <- fig_dataRA[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords, data, binsize = 50.1, maxdist = 1000)
    corr$year = year
    corr$variable = variable
    final_out <- rbind(final_out, corr)
  }
}


final_out1 <- data.frame()
for(i in 1:length(fig_dataRA)){
  temp_out <- fig_dataRA[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out1 <- rbind(final_out1, corr)
  }
}



fig_dataRNA_A<-rbind(final_out,final_out1)
#plot correlograms

############################## RNAB#####################


final_out2 <- data.frame()
for(i in 1:length(fig_dataRB)){
  year <- names(fig_dataRB[i])
  temp_out <- fig_dataRB[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords1, data, binsize = 50.1, maxdist = 1000)
    corr$year = year
    corr$variable = variable
    final_out2 <- rbind(final_out2, corr)
  }
}


final_out3 <- data.frame()
for(i in 1:length(fig_dataRB)){
  year <- names(fig_dataRB[i])
  temp_out <- fig_dataRB[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords1, data, binsize = 50.1, maxdist = 1000)
    corr$year = year
    corr$variable = variable
    final_out3 <- rbind(final_out3, corr)
  }
}



fig_dataRNA_B<-rbind(final_out2,final_out3)
#plot correlograms
################################################# RNA-C ###################################

############################## #####################

final_out4 <- data.frame()
for(i in 1:length(fig_dataRC)){
  year <- names(fig_dataRC[i])
  temp_out <- fig_dataRC[[i]] %>% dplyr:: select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords2, data, binsize = 50.1, maxdist = 1000)
    corr$year = year
    corr$variable = variable
    final_out4 <- rbind(final_out4, corr)
  }
}


final_out5 <- data.frame()
for(i in 1:length(fig_dataRC)){
  year <- names(fig_dataRC[i])
  temp_out <- fig_dataRC[[i]] %>% dplyr:: select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords2, data, binsize = 50.1, maxdist = 1000)
    corr$year = year
    corr$variable = variable
    final_out5 <- rbind(final_out5, corr)
  }
}



fig_dataRNA_c<-rbind(final_out4,final_out5)
#plot correlograms

############################################## RNA-D #############################################################
final_out6 <- data.frame()
for(i in 1:length(fig_dataRD)){
  year <- names(fig_dataRD[i])
  temp_out <- fig_dataRD[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords3, data, binsize = 50.1, maxdist =1000)
    corr$year = year
    corr$variable = variable
    final_out6 <- rbind(final_out6, corr)
  }
}


final_out7 <- data.frame()
for(i in 1:length(fig_dataRD)){
  year <- names(fig_dataRD[i])
  temp_out <- fig_dataRD[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords3, data, binsize = 50.1, maxdist = 1000)
    corr$year = year
    corr$variable = variable
    final_out7 <- rbind(final_out7, corr)
  }
}

fig_dataRNA_d<-rbind(final_out6,final_out7)


final_out8 <- data.frame()
for(i in 1:length(fig_dataRE)){
  year <- names(fig_dataRE[i])
  temp_out <- fig_dataRE[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords4, data, binsize = 50.1, maxdist =1000)
    corr$year = year
    corr$variable = variable
    final_out8 <- rbind(final_out8, corr)
  }
}


final_out9 <- data.frame()
for(i in 1:length(fig_dataRE)){
  year <- names(fig_dataRE[i])
  temp_out <- fig_dataRE[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords4, data, binsize = 50.1, maxdist = 1000)
    corr$year = year
    corr$variable = variable
    final_out9 <- rbind(final_out9, corr)
  }
}

fig_dataRNA_E<-rbind(final_out8,final_out9)


#########################
library(ggplot2)

#plot correlograms


(figA <- ggplot(fig_dataRNA_A) +
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_point(aes(y = Morans.i, x = dist, color = variable), size = 2) +
    geom_line(aes(y = Morans.i, x = dist, color = variable), size = 1) +
    geom_line(aes(y = null.lower, x = dist, color = variable), size = 1, linetype = 'dotted') +
    geom_line(aes(y = null.upper, x = dist, color = variable), size = 1, linetype = 'dotted') +
    scale_color_brewer(palette = 'Set1') + ggtitle("A) UNIT-39")+
    theme(plot.title =element_text(hjust =0))+
    theme_classic(base_line_size = 1, base_rect_size = 1) + 
    theme(legend.position = 'NONE', 
          legend.title = element_blank(), 
          legend.text = element_text(size = 20),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 14, angle = 90),
          strip.text = element_text(size = 24),
          panel.spacing = unit(0.3, units = 'in')) +
    xlab("") +
    scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 200), limits = c(0, 1000)) +
    ylab("") +
    scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.40), minor_breaks = NULL, limits = c(-0.5, 0.8) ) )


(figB <- ggplot(fig_dataRNA_B) +
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_point(aes(y = Morans.i, x = dist, color = variable), size = 2) +
    geom_line(aes(y = Morans.i, x = dist, color = variable), size = 1) +
    geom_line(aes(y = null.lower, x = dist, color = variable), size = 1, linetype = 'dotted') +
    geom_line(aes(y = null.upper, x = dist, color = variable), size = 1, linetype = 'dotted') +
    scale_color_brewer(palette = 'Set1')+ ggtitle("B) UNIT-40")+
    theme(plot.title =element_text(hjust =0))+
    theme_classic(base_line_size = 1, base_rect_size = 1) + 
    theme(legend.position = 'NONE', 
          legend.title = element_blank(), 
          legend.text = element_text(size = 20),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 14, angle = 90),
          strip.text = element_text(size = 24),
          panel.spacing = unit(0.3, units = 'in')) +
    xlab("") +
    scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 200), limits = c(0, 1000)) +
    ylab("") +
    scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.40), minor_breaks = NULL, limits = c(-0.5, 0.8) ) )




(figC <- ggplot(fig_dataRNA_c) +
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_point(aes(y = Morans.i, x = dist, color = variable), size = 2) +
    geom_line(aes(y = Morans.i, x = dist, color = variable), size = 1) +
    geom_line(aes(y = null.lower, x = dist, color = variable), size = 1, linetype = 'dotted') +
    geom_line(aes(y = null.upper, x = dist, color = variable), size = 1, linetype = 'dotted') +
    ggtitle("C) UNIT-43")+
    theme(plot.title =element_text(size=30))+
    scale_color_brewer(palette = 'Set1') +
    theme_classic(base_line_size = 1, base_rect_size = 1) + 
    theme(legend.position = 'NONE', 
          legend.title = element_blank(), 
          legend.text = element_text(size = 20),
          axis.ticks.length = unit(0.25, "cm"),
          axis.ticks.x=element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 14, angle = 90),
          strip.text = element_text(size = 24),
          panel.spacing = unit(0.3, units = 'in'))+
    xlab("Lag-distance (m)") +
    scale_x_continuous(breaks = seq(from = 0, to = 1100, by = 200), limits = c(0, 1100) ) +
    ylab("Moran's I") +
    scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.40), minor_breaks = NULL, limits = c(-0.6, 0.9) ) )



(figD <- ggplot(fig_dataRNA_d) +
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_point(aes(y = Morans.i, x = dist, color = variable), size = 2) +
    geom_line(aes(y = Morans.i, x = dist, color = variable), size = 1) +
    geom_line(aes(y = null.lower, x = dist, color = variable), size = 1, linetype = 'dotted') +
    geom_line(aes(y = null.upper, x = dist, color = variable), size = 1, linetype = 'dotted') +
    ggtitle("D) UNIT-44")+
    theme(plot.title =element_text(hjust =0))+
    scale_color_brewer(palette = 'Set1') +
    theme_classic(base_line_size = 1, base_rect_size = 1) + 
    theme(legend.position = 'NONE', 
          legend.title = element_blank(), 
          legend.text = element_text(size = 20),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_text(size=14),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_text(size = 24),
          panel.spacing = unit(0.3, units = 'in')) +
    xlab("") +
    scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 200), limits = c(0, 1000)) +
    ylab("Moran's I") +
    scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.40), minor_breaks = NULL, limits = c(-0.5, 0.8) ) )

(figE <- ggplot(fig_dataRNA_E) +
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_point(aes(y = Morans.i, x = dist, color = variable), size = 2) +
    geom_line(aes(y = Morans.i, x = dist, color = variable), size = 1) +
    geom_line(aes(y = null.lower, x = dist, color = variable), size = 1, linetype = 'dotted') +
    geom_line(aes(y = null.upper, x = dist, color = variable), size = 1, linetype = 'dotted') +
    ggtitle("E) UNIT-45")+
    theme(plot.title =element_text(hjust =0))+
    scale_color_brewer(palette = 'Set1') +
    theme_classic(base_line_size = 1, base_rect_size = 1) + 
    theme(legend.position = 'NONE', 
          legend.title = element_blank(), 
          legend.text = element_text(size = 20),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_text(size=14),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_text(size = 24),
          panel.spacing = unit(0.3, units = 'in')) +
    xlab("") +
    scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 200), limits = c(0, 1000)) +
    ylab("Moran's I") +
    scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.40), minor_breaks = NULL, limits = c(-0.5, 0.8) ) )

#######################################################################################################################

# read in data
#filename_df <- readOGR("D:/data_groundtruth", layer = 'biomass_both_point')

filename_df <- readOGR("S:\\Biomass\\Biomass_data", layer = 'Both_biomass')

filename<-subset(filename_df,UNIT==38)

fig_dataRA<-list(cbind(filename@data, filename@coords)%>% 
                   dplyr::select(coords.x1,coords.x2,Biomass_20,Biomass_19))

coords <- as.matrix(data.frame('x' = fig_dataRA[[1]]$coords.x1, 'y' = fig_dataRA[[1]]$coords.x2))

filename1<-subset(filename_df,UNIT==41)

fig_dataRB<-list(cbind(filename1@data, filename1@coords)%>% 
                   dplyr::select(coords.x1,coords.x2,Biomass_20,Biomass_19))

coords1 <- as.matrix(data.frame('x' = fig_dataRB[[1]]$coords.x1, 'y' = fig_dataRB[[1]]$coords.x2))


filename2<-subset(filename_df,UNIT==42)

fig_dataRC<-list(cbind(filename2@data, filename2@coords)%>% 
                   dplyr::select(coords.x1,coords.x2,Biomass_20,Biomass_19))

coords2 <- as.matrix(data.frame('x' = fig_dataRC[[1]]$coords.x1, 'y' = fig_dataRC[[1]]$coords.x2))


filename3<-subset(filename_df,UNIT==47)

fig_dataRD<-list(cbind(filename3@data, filename3@coords)%>% 
                   dplyr::select(coords.x1,coords.x2,Biomass_20,Biomass_19))

coords3 <- as.matrix(data.frame('x' = fig_dataRD[[1]]$coords.x1, 'y' = fig_dataRD[[1]]$coords.x2))





filename4<-subset(filename_df,UNIT==48)

fig_dataRE<-list(cbind(filename4@data, filename4@coords)%>% 
                   dplyr::select(coords.x1,coords.x2,Biomass_20,Biomass_19))

coords4 <- as.matrix(data.frame('x' = fig_dataRE[[1]]$coords.x1, 'y' = fig_dataRE[[1]]$coords.x2))




# calculate moran's I and simulated random confidence interval for all response variables
################## RNA_A

final_out <- data.frame()
for(i in 1:length(fig_dataRA)){
  year <- names(fig_dataRA[i])
  temp_out <- fig_dataRA[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords, data, binsize = 50.1, maxdist = 1000)
    corr$year = year
    corr$variable = variable
    final_out <- rbind(final_out, corr)
  }
}


final_out1 <- data.frame()
for(i in 1:length(fig_dataRA)){
  temp_out <- fig_dataRA[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out1 <- rbind(final_out1, corr)
  }
}



fig_dataRNA_A<-rbind(final_out,final_out1)
#plot correlograms

############################## RNAB#####################


final_out2 <- data.frame()
for(i in 1:length(fig_dataRB)){
  year <- names(fig_dataRB[i])
  temp_out <- fig_dataRB[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords1, data, binsize = 50.1, maxdist = 1000)
    corr$year = year
    corr$variable = variable
    final_out2 <- rbind(final_out2, corr)
  }
}


final_out3 <- data.frame()
for(i in 1:length(fig_dataRB)){
  year <- names(fig_dataRB[i])
  temp_out <- fig_dataRB[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords1, data, binsize = 50.1, maxdist = 1000)
    corr$year = year
    corr$variable = variable
    final_out3 <- rbind(final_out3, corr)
  }
}



fig_dataRNA_B<-rbind(final_out2,final_out3)
#plot correlograms
################################################# RNA-C ###################################

############################## #####################

final_out4 <- data.frame()
for(i in 1:length(fig_dataRC)){
  year <- names(fig_dataRC[i])
  temp_out <- fig_dataRC[[i]] %>% dplyr:: select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords2, data, binsize = 50.1, maxdist = 1000)
    corr$year = year
    corr$variable = variable
    final_out4 <- rbind(final_out4, corr)
  }
}


final_out5 <- data.frame()
for(i in 1:length(fig_dataRC)){
  year <- names(fig_dataRC[i])
  temp_out <- fig_dataRC[[i]] %>% dplyr:: select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords2, data, binsize = 50.1, maxdist = 1000)
    corr$year = year
    corr$variable = variable
    final_out5 <- rbind(final_out5, corr)
  }
}



fig_dataRNA_c<-rbind(final_out4,final_out5)
#plot correlograms

############################################## RNA-D #############################################################
final_out6 <- data.frame()
for(i in 1:length(fig_dataRD)){
  year <- names(fig_dataRD[i])
  temp_out <- fig_dataRD[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords3, data, binsize = 50.1, maxdist =1000)
    corr$year = year
    corr$variable = variable
    final_out6 <- rbind(final_out6, corr)
  }
}


final_out7 <- data.frame()
for(i in 1:length(fig_dataRD)){
  year <- names(fig_dataRD[i])
  temp_out <- fig_dataRD[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords3, data, binsize = 50.1, maxdist = 1000)
    corr$year = year
    corr$variable = variable
    final_out7 <- rbind(final_out7, corr)
  }
}

fig_dataRNA_d<-rbind(final_out6,final_out7)


final_out8 <- data.frame()
for(i in 1:length(fig_dataRE)){
  year <- names(fig_dataRE[i])
  temp_out <- fig_dataRE[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords4, data, binsize = 50.1, maxdist =1000)
    corr$year = year
    corr$variable = variable
    final_out8 <- rbind(final_out8, corr)
  }
}


final_out9 <- data.frame()
for(i in 1:length(fig_dataRE)){
  year <- names(fig_dataRE[i])
  temp_out <- fig_dataRE[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords4, data, binsize = 50.1, maxdist = 1000)
    corr$year = year
    corr$variable = variable
    final_out9 <- rbind(final_out9, corr)
  }
}

fig_dataRNA_E<-rbind(final_out8,final_out9)


#########################
library(ggplot2)

#plot correlograms


(figf <- ggplot(fig_dataRNA_A) +
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_point(aes(y = Morans.i, x = dist, color = variable), size = 2) +
    geom_line(aes(y = Morans.i, x = dist, color = variable), size = 1) +
    geom_line(aes(y = null.lower, x = dist, color = variable), size = 1, linetype = 'dotted') +
    geom_line(aes(y = null.upper, x = dist, color = variable), size = 1, linetype = 'dotted') +
    scale_color_brewer(palette = 'Set1') + ggtitle("F) UNIT-38")+
    theme(plot.title =element_text(hjust =0))+
    theme_classic(base_line_size = 1, base_rect_size = 1) + 
    theme(legend.position = c(0.8,1), 
          legend.title = element_blank(), 
          legend.text = element_text(size = 20),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_text(size = 24),
          panel.spacing = unit(0.3, units = 'in')) +
    xlab("") +
    scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 200), limits = c(0, 1000)) +
    ylab("") +
    scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.40), minor_breaks = NULL, limits = c(-0.5, 0.8) ) )




(figg <- ggplot(fig_dataRNA_B) +
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_point(aes(y = Morans.i, x = dist, color = variable), size = 2) +
    geom_line(aes(y = Morans.i, x = dist, color = variable), size = 1) +
    geom_line(aes(y = null.lower, x = dist, color = variable), size = 1, linetype = 'dotted') +
    geom_line(aes(y = null.upper, x = dist, color = variable), size = 1, linetype = 'dotted') +
    ggtitle("G) UNIT-41")+
    theme(plot.title =element_text(hjust =0))+
    scale_color_brewer(palette = 'Set1') +
    theme_classic(base_line_size = 1, base_rect_size = 1) + 
    theme(legend.position = "NONE", 
          legend.title = element_blank(), 
          legend.text = element_text(size = 20),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_text(size = 24),
          panel.spacing = unit(0.3, units = 'in')) +
    xlab("") +
    scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 200), limits = c(0, 1000)) +
    ylab("") +
    scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.40), minor_breaks = NULL, limits = c(-0.5, 0.8) ) )




(figh <- ggplot(fig_dataRNA_c) +
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_point(aes(y = Morans.i, x = dist, color = variable), size = 2) +
    geom_line(aes(y = Morans.i, x = dist, color = variable), size = 1) +
    geom_line(aes(y = null.lower, x = dist, color = variable), size = 1, linetype = 'dotted') +
    geom_line(aes(y = null.upper, x = dist, color = variable), size = 1, linetype = 'dotted') +
    ggtitle("H) UNIT-42")+
    theme(plot.title =element_text(size=30))+
    scale_color_brewer(palette = 'Set1') +
    theme_classic(base_line_size = 1, base_rect_size = 1) + 
    theme(legend.position = 'NONE', 
          legend.title = element_blank(), 
          legend.text = element_text(size = 20),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_text(size = 24),
          panel.spacing = unit(0.3, units = 'in')) +
    xlab("") +
    scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 200), limits = c(0, 1000)) +
    ylab("") +
    scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.40), minor_breaks = NULL, limits = c(-0.5, 0.8) ) )


(figi <- ggplot(fig_dataRNA_d) +
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_point(aes(y = Morans.i, x = dist, color = variable), size = 2) +
    geom_line(aes(y = Morans.i, x = dist, color = variable), size = 1) +
    geom_line(aes(y = null.lower, x = dist, color = variable), size = 1, linetype = 'dotted') +
    geom_line(aes(y = null.upper, x = dist, color = variable), size = 1, linetype = 'dotted') +
    ggtitle("I) UNIT-47")+
    theme(plot.title =element_text(hjust =0))+
    scale_color_brewer(palette = 'Set1') +
    theme_classic(base_line_size = 1, base_rect_size = 1) + 
    theme(legend.position = 'NONE', 
          legend.title = element_blank(), 
          legend.text = element_text(size = 20),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_text(size = 24),
          panel.spacing = unit(0.3, units = 'in')) +
    xlab("") +
    scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 200), limits = c(0, 1000)) +
    ylab("") +
    scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.40), minor_breaks = NULL, limits = c(-0.5, 0.8) ) )


(figj <- ggplot(fig_dataRNA_E) +
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_point(aes(y = Morans.i, x = dist, color = variable), size = 2) +
    geom_line(aes(y = Morans.i, x = dist, color = variable), size = 1) +
    geom_line(aes(y = null.lower, x = dist, color = variable), size = 1, linetype = 'dotted') +
    geom_line(aes(y = null.upper, x = dist, color = variable), size = 1, linetype = 'dotted') +
    ggtitle("J) UNIT-48")+
    theme(plot.title =element_text(hjust =0))+
    scale_color_brewer(palette = 'Set1') +
    theme_classic(base_line_size = 1, base_rect_size = 1) + 
    theme(legend.position = 'NONE', 
          legend.title = element_blank(), 
          legend.text = element_text(size = 20),
          axis.ticks.length = unit(0.25, "cm"),
          axis.text.x = element_text(size=14),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          strip.text = element_text(size = 24),
          panel.spacing = unit(0.3, units = 'in')) +
    xlab("") +
    scale_x_continuous(breaks = seq(from = 0, to = 1000, by = 200), limits = c(0, 1000)) +
    ylab("Moran's I") +
    scale_y_continuous(breaks = seq(from = -1, to = 1, by = 0.40), minor_breaks = NULL, limits = c(-0.5, 0.8) ) )

library(cowplot)
library(grid)
library(gridExtra)

grids_bs <- plot_grid(figA,figf,figB,figg,figC,figh,figD,figi,figE,figj,
                      ncol = 2, align = "v")

y.grob <- textGrob("Moran's I", 
                   gp=gpar( col="black", fontsize=30), rot=90)

x.grob <- textGrob("Lag-distance (m)", 
                   gp=gpar( col="black", fontsize=30))

grid.arrange(arrangeGrob(grids_bs, left = y.grob, bottom = x.grob))





