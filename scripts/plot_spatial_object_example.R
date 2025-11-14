# example analysis for Wen Rui plotting a spatial file (geometry plus raster)
# two plots: contineuous axis and categorial exis

# load packages
require(terra);require(sf);require(data.table);require(ggplot2)
require(scales)

# read in raster map (this is the p_z file send by Siattwinda)
r1 <- rast('D:/ESA/02 phd projects/05 siatwinda/02 data/phosphate.tif')

# read in the polygon map of zambia (the gzambia file send by Siattwinda)
s1 <- st_read('D:/ESA/02 phd projects/05 siatwinda/02 data/zambia.gpkg')

# disaggregate the raster (to simplify plotting, not needed, but it makes testing faster)
r2 <- aggregate(r1, fact = 25, fun = 'mean')

# to do (difficult) analyses for all rasters it is pratical to convert to data.table first
d1 <- as.data.frame(r2,xy=TRUE)
d1 <- as.data.table(d1)

# convert unit from mg/100kg to mg/kg
d1[,p2 := phosphate * 0.01]

# here a fake addition to set a P dose depending on Pmehlig
d1[p2 <= 7,pdose := 85]
d1[p2 > 7 & p2 <= 10,pdose := 55]
d1[p2 > 10,pdose := 5]

# some computers print the degrees sign wrongly in ggplot, then this is a work around
xlabs <- seq(22,34,2)
ylabs <- seq(-8,-18,-2)

# here a single plot using a continuous legend
p1 <- ggplot(data = s1) + geom_sf(color = "black", fill = "gray92") +
  geom_tile(data = d1,aes(x=x,y=y,fill= pdose)) +
  scale_fill_viridis_c()+ theme_bw() + #theme_void() +
  theme(legend.position = c(0.25,0.85),legend.direction="horizontal") +
  scale_x_continuous(breaks = xlabs, labels = paste0(xlabs,'°E')) +
  scale_y_continuous(breaks = ylabs, labels = paste0(abs(ylabs),'°S')) +
  xlab("") + ylab("") +
  ggtitle("Recommende P dose", subtitle = "scenario 1") +
  coord_sf(crs = 4326) +
  theme(legend.text=element_text(size=8),legend.title = element_text(size=10))+
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))

# here a single plot using a discete legend
pbreak <- c(0,20,60,1000)
plabel <- c('<20','20-60','>60')
p2 <- ggplot(data = s1) + geom_sf(color = "black", fill = "gray92") +
  geom_tile(data = d1,aes(x=x,y=y,fill= cut(pdose, pbreak,labels = plabel))) +
  scale_fill_viridis_d()+ theme_bw() + #theme_void() +
  theme(legend.position = c(0.25,0.85),legend.direction="horizontal") +
  scale_x_continuous(breaks = xlabs, labels = paste0(xlabs,'°E')) +
  scale_y_continuous(breaks = ylabs, labels = paste0(abs(ylabs),'°S')) +
  labs(fill = 'P dose (kg /ha)')+
  xlab("") + ylab("") +
  ggtitle("Recommende P dose", subtitle = "scenario 1") +
  coord_sf(crs = 4326)+
  theme(legend.text=element_text(size=8),legend.title = element_text(size=10))+
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5))

# you can combine both plots using patchwork
require(patchwork)

# combine both figures
p3 <- p1 + p2
ggsave(plot = p3, filename = 'D:/ESA/02 phd projects/05 siatwinda/02 data/test_pdose.png',width = 25,height = 14,units='cm')
