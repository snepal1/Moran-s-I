
rm(list=ls())

#### Load the required libraries all at once ################
#define vector of packages to load ########################
some_packages <- c('ggplot2', 'dplyr', 'rgdal','spdep','mapview',
                   'grid','gridExtra')

#load all packages at once
lapply(some_packages, library, character.only=TRUE)

#library(rgdal)
#library(spdep)
#library(dplyr)
#library(mapview)

# Read in the data as the point shape file ####################################

filename_df <- readOGR("D:\\Moran's\\Moran-s-I", layer = 'Both_biomass')

#mapview(filename_df)

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


#mapview(filename_df)
########### Subset the data to required units #################################################
filename<-subset(filename_df,UNIT==39)
filename1<-subset(filename_df,UNIT==40)
filename2<-subset(filename_df,UNIT==43)
filename3<-subset(filename_df,UNIT==44)
filename4<-subset(filename_df,UNIT==45)

############# write a function to extract the spatial and non-spatial data that will go
#in Moran's I calculation function'
process_data <- function(data) {
  # Combine data and coordinates
  fig_dataRA<-list(cbind(data@data, data@coords)%>% 
                     dplyr::select(coords.x1,coords.x2,Biomass_20,Biomass_19))
  
  # Extract coordinates and create a matrix
  coords <- as.matrix(data.frame('x' = fig_dataRA[[1]]$coords.x1, 'y' = fig_dataRA[[1]]$coords.x2))
  
  # Return both fig_dataRA and coords in a list
  return(list(fig_dataRA = fig_dataRA, coords = coords))
}

# Return the co-ordinate and the data for unit 39:
result <- process_data(filename)
fig_dataRA <- result$fig_dataRA
coords <- result$coords
##### Return the co-ordinate and the data for unit 40:
result <- process_data(filename1)
fig_dataRB <- result$fig_dataRA
coords1 <- result$coords

##### Return the co-ordinate and the data for unit 43:
result <- process_data(filename2)
fig_dataRC <- result$fig_dataRA
coords2 <- result$coords
##### Return the co-ordinate and the data for unit 44:
result <- process_data(filename3)
fig_dataRD <- result$fig_dataRA
coords3 <- result$coords
##### Return the co-ordinate and the data for unit 45:
result <- process_data(filename4)
fig_dataRE <- result$fig_dataRA
coords4 <- result$coords


##################### Check the neighborhood object plot #############################

iweights <- function(locations){
  # Create an empty plot to initialize the plotting area
  plot(1, type = "n", xlab = "UTMX", ylab = "UTMY", xlim = c(min(locations[, 1]), max(locations[, 1])),
       ylim = c(min(locations[, 2]), max(locations[, 2])))
  
  # Add the spatial weights using the dnearneigh function
  neigh <- dnearneigh(x = locations, d1 = 0, d2 = 1000, longlat = FALSE)
  P <- plot(neigh, coords = locations, add = TRUE, col = "red")
  
  return(P)
}

iweights(coords2)
############# calculate moran's I and simulated 95 %  confidence envelope ##################################
################## 

final_out <- data.frame()
for(i in 1:length(fig_dataRA)){
  temp_out <- fig_dataRA[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords, data, binsize = 50.1, maxdist = 1000)
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


############### Put it back as a data frame #####################
fig_dataRNA_A<-rbind(final_out,final_out1)
#str(fig_dataRNA_A)

############################## RNAB #####################


final_out2 <- data.frame()
for(i in 1:length(fig_dataRB)){
  temp_out <- fig_dataRB[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords1, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out2 <- rbind(final_out2, corr)
  }
}


final_out3 <- data.frame()
for(i in 1:length(fig_dataRB)){
  temp_out <- fig_dataRB[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords1, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out3 <- rbind(final_out3, corr)
  }
}

############### Put it back as a data frame #####################
fig_dataRNA_B<-rbind(final_out2,final_out3)


################################################# RNA-C ###################################

############################## #############################################################################

final_out4 <- data.frame()
for(i in 1:length(fig_dataRC)){
  temp_out <- fig_dataRC[[i]] %>% dplyr:: select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords2, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out4 <- rbind(final_out4, corr)
  }
}


final_out5 <- data.frame()
for(i in 1:length(fig_dataRC)){
  temp_out <- fig_dataRC[[i]] %>% dplyr:: select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords2, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out5 <- rbind(final_out5, corr)
  }
}

############### Put it back as a data frame #####################
fig_dataRNA_c<-rbind(final_out4,final_out5)


############################################## RNA-D #############################################################
final_out6 <- data.frame()
for(i in 1:length(fig_dataRD)){
  temp_out <- fig_dataRD[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords3, data, binsize = 50.1, maxdist =1000)
    corr$variable = variable
    final_out6 <- rbind(final_out6, corr)
  }
}


final_out7 <- data.frame()
for(i in 1:length(fig_dataRD)){
  temp_out <- fig_dataRD[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords3, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out7 <- rbind(final_out7, corr)
  }
}
############### Put it back as a data frame #####################
fig_dataRNA_d<-rbind(final_out6,final_out7)

############################################## RNA-E #############################################################
final_out8 <- data.frame()
for(i in 1:length(fig_dataRE)){
  temp_out <- fig_dataRE[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords4, data, binsize = 50.1, maxdist =1000)
    corr$variable = variable
    final_out8 <- rbind(final_out8, corr)
  }
}


final_out9 <- data.frame()
for(i in 1:length(fig_dataRE)){
  temp_out <- fig_dataRE[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords4, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out9 <- rbind(final_out9, corr)
  }
}

############### Put it back as a data frame #####################
fig_dataRNA_E<-rbind(final_out8,final_out9)


############################ Creating the combined plots usign ggplot2 #########################################
#plot correlograms
library(ggplot2)

(figA <- ggplot(fig_dataRNA_A)
  + theme_classic() 
  + geom_hline(yintercept = 0)+
    geom_point(aes(y = Morans.i, x = dist, color = variable), size = 2)+
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

################################## Let us do it for the HID units#####################################################################################

############################# Subset the data to the required unit only ###########################################

filename5<-subset(filename_df,UNIT==38)
filename6<-subset(filename_df,UNIT==41)
filename7<-subset(filename_df,UNIT==42)
filename8<-subset(filename_df,UNIT==47)
filename9<-subset(filename_df,UNIT==48)

################################## Use the function created initially to process the data#########################
# Return the co-ordinate and the data for unit 38:
result <- process_data(filename5)
fig_dataRA5 <- result$fig_dataRA
coords5 <- result$coords
##### Return the co-ordinate and the data for unit 41:
result <- process_data(filename6)
fig_dataRB6 <- result$fig_dataRA
coords6 <- result$coords

##### Return the co-ordinate and the data for unit 42:
result <- process_data(filename7)
fig_dataRC7 <- result$fig_dataRA
coords7 <- result$coords
##### Return the co-ordinate and the data for unit 47:
result <- process_data(filename8)
fig_dataRD8 <- result$fig_dataRA
coords8 <- result$coords
##### Return the co-ordinate and the data for unit 48:
result <- process_data(filename9)
fig_dataRE9 <- result$fig_dataRA
coords9 <- result$coords


# calculate moran's I and simulated random confidence interval for all response variables
################## RNA_A

final_out5 <- data.frame()
for(i in 1:length(fig_dataRA5)){
  temp_out <- fig_dataRA5[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords5, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out5 <- rbind(final_out5, corr)
  }
}


final_out5A <- data.frame()
for(i in 1:length(fig_dataRA5)){
  temp_out <- fig_dataRA5[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords5, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out5A <- rbind(final_out5A, corr)
  }
}



fig_dataRNA_A<-rbind(final_out5,final_out5A)


############################## RNAB##################################################################


final_out6 <- data.frame()
for(i in 1:length(fig_dataRB6)){
  temp_out <- fig_dataRB6[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords6, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out6 <- rbind(final_out6, corr)
  }
}


final_out6A <- data.frame()
for(i in 1:length(fig_dataRB6)){
  temp_out <- fig_dataRB6[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords6, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out6A <- rbind(final_out6A, corr)
  }
}



fig_dataRNA_B<-rbind(final_out6,final_out6A)
#plot correlograms
################################################# RNA-C ###################################

############################## #####################

final_out7 <- data.frame()
for(i in 1:length(fig_dataRC7)){
  temp_out <- fig_dataRC7[[i]] %>% dplyr:: select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords7, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out7 <- rbind(final_out7, corr)
  }
}


final_out7A <- data.frame()
for(i in 1:length(fig_dataRC7)){
  temp_out <- fig_dataRC7[[i]] %>% dplyr:: select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords7, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out7A <- rbind(final_out7A, corr)
  }
}



fig_dataRNA_c<-rbind(final_out7,final_out7A)
#plot correlograms

############################################## RNA-D #############################################################
final_out8 <- data.frame()
for(i in 1:length(fig_dataRD8)){
  temp_out <- fig_dataRD8[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords8, data, binsize = 50.1, maxdist =1000)
    corr$variable = variable
    final_out8 <- rbind(final_out8, corr)
  }
}


final_out8A <- data.frame()
for(i in 1:length(fig_dataRD8)){
  temp_out <- fig_dataRD8[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords8, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out8A <- rbind(final_out8A, corr)
  }
}

fig_dataRNA_d<-rbind(final_out8,final_out8A)


final_out9 <- data.frame()
for(i in 1:length(fig_dataRE9)){
  temp_out <- fig_dataRE9[[i]] %>% dplyr::select(Biomass_20)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords9, data, binsize = 50.1, maxdist =1000)
    corr$variable = variable
    final_out9 <- rbind(final_out9, corr)
  }
}


final_out9A <- data.frame()
for(i in 1:length(fig_dataRE9)){
  temp_out <- fig_dataRE9[[i]] %>% dplyr::select(Biomass_19)
  for(j in 1:ncol(temp_out)){
    data <- temp_out[,j]
    variable <- colnames(temp_out[j])
    corr <- icorrelogram(locations = coords9, data, binsize = 50.1, maxdist = 1000)
    corr$variable = variable
    final_out9A <- rbind(final_out9A, corr)
  }
}

fig_dataRNA_E<-rbind(final_out9,final_out9A)


######################### Plot the correlograms ##############



(figf <- ggplot(fig_dataRNA_A) +
    theme_classic() +
    geom_hline(yintercept = 0) +
    geom_point(aes(y = Morans.i, x = dist, color = variable), size = 2) +
    geom_line(aes(y = Morans.i, x = dist, color = variable), size = 1) +
    geom_line(aes(y = null.lower, x = dist, color = variable), size = 1, linetype = 'dotted') +
    geom_line(aes(y = null.upper, x = dist, color = variable), size = 1, linetype = 'dotted') +
    scale_color_brewer(palette = 'Set1') + 
    scale_color_manual(name = "Legend Title", 
                      values = c("Biomass_19" = "red", "Biomass_20" = "#0072B2"), 
                      labels = c("1934", "2016")) +
   ggtitle("F) UNIT-38")+
    theme(plot.title =element_text(hjust =0))+
    theme_classic(base_line_size = 1, base_rect_size = 1) + 
    theme(legend.position = c(0.8,1), 
          legend.title = element_blank(), 
          legend.text=element_text(size=20),
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

#library(cowplot)
#library(grid)
#library(gridExtra)

grids_bs <- plot_grid(figA,figf,figB,figg,figC,figh,figD,figi,figE,figj,
                      ncol = 2, align = "v")

y.grob <- textGrob("Moran's I", 
                   gp=gpar( col="black", fontsize=30), rot=90)

x.grob <- textGrob("Lag-distance (m)", 
                   gp=gpar( col="black", fontsize=30))

grid.arrange(arrangeGrob(grids_bs, left = y.grob, bottom = x.grob))

################################ End of the run ###########################################################################


