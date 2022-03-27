# Installing and loading libraries
library(tidyverse)
library(sf)
# install.packages("rnaturalearth")
library(rnaturalearth)
# install.packages("rnaturalearthdata")
library(rnaturalearthdata)
# devtools::install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)
library(rgeos)
library(fs)
# library(devtools)
# install_github("jgx65/hierfstat")
library(hierfstat)
#install.packages("viridis")
library(viridis)
library(ggmap)

# Setting the ggplot2 theme to bw
theme_set(theme_bw())

# Loading in the data to display
#world <- ne_coastline(scale="large",returnclass="sp")
world <- ne_countries(scale="large",returnclass="sp")

# Because of the Eurocentric nature of the R community
# most plots cannot span "0" in longitude, despite this
# projection being relevant to half of the globe.
# Fortunately, a solution is available at:
#https://rpubs.com/valentin/pacific-centered-map-voronoi-tessellation

# We need to cut out the bit of the world we are not interested in
box_cut <- bbox2SP(n=90,s=-90,w=-70,e=35,proj4string = world@proj4string)
world_crop <- gDifference(world, box_cut)

# Now we will further pair to just lookng at the longitudes we are interested in
sf_use_s2(FALSE) # After updating sf to > 1.0 need this to make it work
final_crop <- world_crop %>% st_as_sf() %>% 
  st_shift_longitude() %>% 
  st_crop(c(xmin = st_bbox(.)[["xmin"]],
            xmax = st_bbox(.)[["xmax"]],
            ymin = -65,
            ymax = 65))

ggplot(data = final_crop) +
  geom_sf() +
  coord_sf(expand = FALSE)

# Creating haplotype "rings" following https://www.r-graph-gallery.com/128-ring-or-donut-plot.html
# Reading in data:
data <- read_tsv("/Users/aleal62p/Dropbox (Otago University)/spermwhale_taranaki/haplotype_data_21Nov2021.txt")

# An example of what data looks like. Haplotype is in left most
# column, followed by count of each haplotype in each population
# in the subsequent columns. Count rather than freq data is 
# necessary for downstream FST calculations
#Haplotype `Crozet Is.` `Aldabra Is` Seychelles Mauritius Kerguelen
#<chr>            <dbl>        <dbl>      <dbl>     <dbl>     <dbl>
#  1 A                    1            4          2         7         0
#2 B                    1            0          0         0         0
#3 C                    0            5          5        27         0
#4 D                    0            0          0         0         0
#5 E                    0            0          0         0         0
#6 F                    0            0          0         0         0
#7 H                    0            0          0         0         0
#8 J                    0            0          6         0         0
#9 K                    0            0          0         0         0
#10 L                    0            0          0         0         0

# Defining colour palette

haplotype_colours <- c("#00B0F2",
                       "#FF0001",
                       "#92D24F",
                       "#FFC000",
                       "#00B050",
                       "#FFFF02",
                       "#0072C1",
                       "#6E309E",
                       "#00215D",
                       "#FDFFBD",
                       "#CE00CE",
                       "#C00000",
                       "#984907",
                       "#01CB9A",
                       "#2F33FE",
                       "#CDCC00",
                       "#FF0069",
                       "#EBF0DF",
                       "#01FFFD",
                       "#C5D9F2",
                       "#C967FE",
                       "#604A7D",
                       "#FEE9D9",
                       "#E7E0EE",
                       "#7E7E7E",
                       "#0C0C0C",
                       "#4B4429",
                       "#32869D",
                       "#02CBFF",
                       "#DBEDF4",
                       "#C05C0A",
                       "#CDBFDA",
                       "#F4DBDB",
                       "#F1F1F1",
                       "#F802FA",
                       "#DFD9C2",
                       "#E4B8B8",
                       "#A7A7A7",
                       "#C5E0B4",
                       "#FFFFFF",
                       "#FBD5B5",
                       "#385723")

names(haplotype_colours) <- c("A",
                       "B",
                       "C",
                       "J",
                       "E",
                       "N",
                       "KK",
                       "D",
                       "O",
                       "I",
                       "H",
                       "SEAUS7",
                       "SW.M2",
                       "M",
                       "S",
                       "F",
                       "G",
                       "Q",
                       "T",
                       "Z",
                       "EE",
                       "GG",
                       "K",
                       "P",
                       "AA",
                       "CC",
                       "FF",
                       "HH",
                       "II",
                       "NN",
                       "SEAUS8",
                       "L",
                       "R",
                       "U",
                       "BB",
                       "DD",
                       "LL",
                       "OO",
                       "SW.K1",
                       "New",
                       "JJ",
                       "SEAUS1")

# Stepping through each location
for (i in 2:dim(data)[2]) {
  # Grabbing name
  region_name <- colnames(data)[i]
  # Pulling out the data and haplotype labels
  data_subset <- data[-which(data[,i]==0),c(1,i)]
  # Renaming to be consistent with example
  names(data_subset) <- c("category","count")
  # Calculating proportions
  data_subset$proportion <- data_subset$count / sum(data_subset$count)
  # Cumulative proportions
  data_subset$ymax <- cumsum(data_subset$proportion)
  # Bottom of each triangle
  data_subset$ymin <- c(0, head(data_subset$ymax, n=-1))
  # label position
  data_subset$labelPosition <- (data_subset$ymax + data_subset$ymin) / 2
  
  # subset colours for plotting
  colour_subset <- haplotype_colours[which(names(haplotype_colours) %in% data_subset$category)]
  
  # subset label categories if they are too narrow
  category_subset <- data_subset$category
  category_subset[data_subset$proportion<0.1] <- ""
  
  # Plotting
  ggplot(data_subset, aes(ymax=ymax, ymin=ymin, xmax=1, xmin=0, fill=category)) +
    geom_rect(color="black") +
    scale_fill_manual(values=colour_subset) +
    geom_text(x=0.5, aes(y=labelPosition, label=category_subset), size=12) +
    coord_polar(theta="y") +
    xlim(c(-1, 4)) +
    theme_void() +
    theme(legend.position = "none")
  
  ggsave(paste("/Users/aleal62p/Dropbox (Otago University)/spermwhale_taranaki/",gsub(" ","_",path_sanitize(region_name,"")),".pdf",sep=""))

}

# Pairwise FST values
# Need to convert data into hierfstat compatible format
pop_key <- rbind(seq(0,length(names(data))-1),names(data))
hap_key <- rbind(seq(1,dim(data)[1]),as.matrix(data[1:dim(data)[1],1])[,1])
hierfstat_dat <- NULL
for (i in 2:dim(data)[2]) {
  tempdata <- cbind(hap_key[1,],data[,i])
  tempdata <- tempdata[which(tempdata[,2]!=0),]
  for (j in 1:dim(tempdata)[1]) {
    tempheirfstat <- cbind((i-1),tempdata[j,1])
    for (k in 1:tempdata[j,2]) {
      hierfstat_dat <- rbind(hierfstat_dat,tempheirfstat)
    }
  }
}

# Need to convert into numerical dataframe
hierfstat_dat <- as.data.frame(hierfstat_dat)
names(hierfstat_dat) <- c("Pop","loc-1")
hierfstat_dat[,1] <- as.numeric(hierfstat_dat[,1])
hierfstat_dat[,2] <- as.numeric(hierfstat_dat[,2])

WCfst <- pairwise.WCfst(dat = hierfstat_dat,diploid=FALSE)

# Converting to long form for heatmap plotting
WCfst_long <- NULL
for (i in 1:(dim(WCfst)[2]-1)) {
  for (j in (i+1):dim(WCfst)[2]) {
    temp_long <- cbind(pop_key[2,(i+1)],pop_key[2,(j+1)],WCfst[i,j])
    WCfst_long <- rbind(WCfst_long,temp_long)
  }
}

# Exact test for population differentiation
pvalue <- NULL
for (i in 2:(dim(data)[2]-1)) {
  for (j in (i+1):dim(data)[2]) {
    temppvalue <- fisher.test(cbind(data[,i],data[,j]),simulate.p.value=TRUE)$p.value
    pvalue <- c(pvalue,temppvalue)
  }
}

# Binding them together
WCfst_long <- cbind(WCfst_long,pvalue)
WCfst_long <- as.data.frame(WCfst_long)
WCfst_long[,3] <- as.numeric(WCfst_long[,3])
WCfst_long[,4] <- as.numeric(WCfst_long[,4])

# "zeroing" out FST values if negative, or p-values not significant
WCfst_heatmap <- WCfst_long %>% mutate(V3=ifelse(V3<0,0,V3))
WCfst_heatmap$V1 <- factor(WCfst_heatmap$V1,levels=rev(c(pop_key[-1,-1])))
WCfst_heatmap$V2 <- factor(WCfst_heatmap$V2,levels=c(pop_key[-1,-1]))

# Plotting the heatmap for regional populations only
WCfst_regional <- WCfst_heatmap %>% filter(V1!="TOTAL") %>% filter(V1!="Other Pacific Ocean") %>% filter(V1!="New Zealand total") %>% filter(V1!="Australia total") %>% filter(V1!="SE Australia") %>% filter(V1!="SW Australia") %>% filter(V1!="Other Indian Ocean") %>% filter(V1!="Southern Hemisphere")%>% filter(V2!="TOTAL") %>% filter(V2!="Other Pacific Ocean") %>% filter(V2!="New Zealand total") %>% filter(V2!="Australia total") %>% filter(V2!="SE Australia") %>% filter(V2!="SW Australia") %>% filter(V2!="Other Indian Ocean") %>% filter(V2!="Southern Hemisphere")

# Obtaining colour scheme for the labels of the different regions
label_colours_x <- c(rep("dark green",7),rep("dark orange",6),rep("red",3),rep("blue",14))
label_colours_y <- rev(c(rep("dark green",8),rep("dark orange",6),rep("red",3),rep("blue",13)))

ggplot(WCfst_regional,aes(V2,V1,fill=V3)) + geom_tile() + scale_fill_viridis(discrete=FALSE) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(fill = expression(F[ST])) +
  theme(axis.text.x = element_text(colour = label_colours_x)) +
  theme(axis.text.y = element_text(colour = label_colours_y))

ggsave("/Users/aleal62p/Dropbox (Otago University)/spermwhale_taranaki/regional_fsts.pdf",width=175,height=150,units="mm")

# Plotting the heatmap for 2018 + summary pops
WCfst_sum <- WCfst_heatmap %>% filter(V1=="Other Pacific Ocean" | V1=="NZ females" | V1=="NZ males (inc. 2018)" | V1=="2018 stranding" | V1=="SE Australia" | V1=="SW Australia" | V1=="Other Indian Ocean")  %>% filter(V2=="Other Pacific Ocean" | V2=="NZ females" | V2=="NZ males (inc. 2018)" | V2=="2018 stranding" | V2=="SE Australia" | V2=="SW Australia" | V2=="Other Indian Ocean") 

ggplot(WCfst_sum,aes(V2,V1,fill=V3)) + geom_tile() + scale_fill_viridis(discrete=FALSE) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(fill = expression(F[ST]))

ggsave("/Users/aleal62p/Dropbox (Otago University)/spermwhale_taranaki/sum_fsts.pdf",width=175,height=150,units="mm")

# Outputting 2018 stranding's statistics
WCfst_long %>% filter(V1=="2018 stranding" | V2=="2018 stranding")

# Outputting NZ male's statistics
WCfst_long %>% filter(V1=="NZ males (inc. 2018)" | V2=="NZ males (inc. 2018)") %>% mutate_if(is.numeric,~round(.,digits=7))

# Outputting NZ female's statistics
WCfst_long %>% filter(V1=="NZ females" | V2=="NZ females")


# Plotting finer scale New Zealand
# Taranaki male lat/long details for unique haps
taranaki <- rbind(c("T7052","-39.5797","174.1152","O"),
                  c("T7053","-39.5797","174.1152","C"),
                  c("T7054","-39.5797","174.1152","New"),
                  c("T7055","-39.5797","174.1152","Z"),
                  c("T7059","-39.5797","174.1152","M"),
                  c("T7061","-39.5797","174.1152","B"),
                  c("T7064","-39.5797","174.1152","A"))
                  
taranaki <- as.data.frame(taranaki)                  
names(taranaki) <- c("Sample_ID","Latitude","Longitude","Haplotype")
taranaki[,2] <- as.numeric(taranaki[,2])
taranaki[,3] <- as.numeric(taranaki[,3])

# Obtaining map
nzmap <- map_data("nz")

# Loading in data from Alexander et al. (2016) with lat/longs for NZ
alexanderspermwhales <- read_tsv("/Users/aleal62p/Dropbox (Otago University)/spermwhale_taranaki/Sperm_whale_DNA_profiles_AA_16Aug2014.txt")

nzonly <- alexanderspermwhales %>% filter(Region=="NZ") 
nzonly$Latitude <- as.numeric(nzonly$Latitude)
nzonly$Longitude <- as.numeric(nzonly$Longitude)

nzonly_NA <- nzonly %>% filter(is.na(Latitude) | is.na(Longitude))
nzonly <- nzonly %>% filter(!is.na(Latitude) & !is.na(Longitude))
nzonly_Chathams <- nzonly %>% filter(Longitude<0)
nzonly_noChathams <- nzonly %>% filter(Longitude>0)

nz_haplotype_colours <- haplotype_colours[names(haplotype_colours) %in% unique(c(taranaki$Haplotype, nzonly$Haplotype, nzonly_NA$Haplotype))]

ggplot(nzmap) +
  geom_polygon(mapping=aes(x=long,y=lat,group=group),color="black",fill="light grey") +
  coord_quickmap() +
  theme_bw(base_size = 15) +
  geom_point(taranaki,mapping=aes(x=Longitude,y=Latitude,fill=Haplotype),shape=21,size=5) +
  geom_jitter(nzonly_noChathams,mapping=aes(x=Longitude,y=Latitude,fill=Haplotype),shape=21,size=5,width=0.25) +
  scale_fill_manual(name="Haplotype",values=nz_haplotype_colours) +
  theme(legend.position = "left") +
  labs(y= "Latitude", x = "Longitude")

# Copy the image directly into powerpoing  
