library(ggplot2)
library(senamhiR)
############################################Interpolacion spline pp_Moquegua
###################################################################
###################################################################
###################################################################
path <- "C:/Users/ASUS/Desktop/Moquegua"
path
setwd("C:/Users/ASUS/Desktop/Moquegua")
datos <- list.files(path=path,pattern = ".*\\.txt$", full.names = F)
datos_full <- list.files(path=path,pattern = ".*\\.txt$", full.names = T)

code <- substr(datos,5,10)

all_data <- station_search()

#all_data[all_data$StationID== "000235",]


library(dplyr)

metadata <- NULL
for (i in 1:length(code)){
  metadata <- rbind(metadata,all_data[all_data$StationID == code[i],])
}
metadata

ggplot(data=metadata,aes(x=Longitude,y=Latitude,label=Station))+
  geom_point()+geom_label()

all_data$Latitude



sum_na <- function(x) {if (all(is.na(x)) == TRUE | sum(is.na(x))*100/length(x) >20) {NA} else {sum(x, na.rm = TRUE)}}

mean_na <- function(x) {if (all(is.na(x)) == TRUE | sum(is.na(x))*100/length(x) >20) {NA} else {mean(x, na.rm = TRUE)}}

i=1

final_table <- NULL

for (i in 1:length(datos_full)){
  data <- read.table(datos_full[i],na.strings = c(-99.9,-999,-888,888,999,"NA"))
  names(data) <- c("Año","Mes","Dia","Pp","Tmax","Tmin")  ###Cambiar nombre de columnas
  data$fecha <- as.Date(paste0(data$Año,"-",data$Mes,"-",data$Dia))  ###Crear fechas
  data$Tmax[data$Tmin > data$Tmax] <- NA     ###Pequeño control de calidad
  data$Tmin[data$Tmin > data$Tmax] <- NA
  data$Pp[data$Pp < 0] <- NA
  data$Tmed <- (data$Tmax + data$Tmin)/2    
  data$Hr <- ((6.11*(10**((7.5*data$Tmin)/(237.3+data$Tmin)))) / (6.11*(10**((7.5*data$Tmax)/(237.3+data$Tmax)))))*100 ###Hr
  
  
  data <- data %>% filter(fecha >= as.Date("2003-01-01") & fecha <= as.Date("2012-12-31"))%>%  group_by(Año,Mes) %>% summarise(Pp= sum_na(Pp),Tmed= mean_na(Tmed),Tmax=mean_na(Tmax),Tmin= mean_na(Tmin),Hr = mean_na(Hr))
  
  Fecha_month <- as.Date(paste0(data$Año,"-",data$Mes,"-",01))
  
  data <- data.frame(Fecha = Fecha_month,Latitud = metadata$Latitude[i],Longitud = metadata$Longitude[i],Altitud = metadata$Altitude[i],Estacion = metadata$Station[i],data[,-c(1,2)])
  
  final_table <- rbind(final_table,data)
  
}


View(metadata)
View(final_table)

#6120/(51*12)

#observaciones_regular <- kriging::kriging(x = observaciones$lon, 
#                                          y = observaciones$lat,
#                                         response = observaciones$tmax_media,polygons = p)

#library(terra)
#rast("https://data.chc.ucsb.edu/products/CHIRTSmonthly/CHIRTSmax.CDR/CHIRTSmax.1983.01.tif")



interpolado_spline <- NULL
interpolado_kriging <- NULL
i=1
for (i in 1:length(unique(final_table$Fecha))){
  #####Filtrando fecha 
  final_table2 <- (final_table[final_table$Fecha == unique(final_table$Fecha)[i],])
  
  #####Eliminando NAs
  final_table2 <- final_table2[is.na(final_table2$Hr) ==F,]
  
  ######Creando rejilla a interpolar
  y = seq(-19.25,-14.25,0.5)
  x = seq(-72.75,-68.25,0.5)
  length(y)
  length(x)
  
  #unique(mon_data$Latitud)[order(unique(mon_data$Latitud))]
  #unique(mon_data$Longitud)
  grid1 <- expand.grid(x,y)
  
  coord<- cbind(final_table2$Longitud,final_table2$Latitud)
  Pp<- final_table2$Hr
  
  
  #############SPLINE
  
  library(fields)
  fit<- Tps(coord,Pp)
  library(gstat)
  pred <- predict(fit, grid1)
  
  spline <- data.frame(grid1,pred,Fecha= unique(final_table$Fecha)[i])
  names(spline) <- c("Longitud","Latitud","Pp","Fecha")
  
  
  spline$Pp[spline$Pp < 0] <- 0
  spline$Pp[spline$Pp > 100] <- 100
  
  interpolado_spline <- rbind(interpolado_spline,spline)
  
  #################################Kriging
  
  
  grid <- expand.grid(x,y)
  names(grid) <- c("Longitud","Latitud")
  
  library(sf)
  
  ####Puntos iniciales en fomato sf
  sf<- st_as_sf(final_table2, coords = c("Longitud", "Latitud"), remove = FALSE)
  #plot(sf)
  
  ####Puntos finales en formato sf
  sf2<- st_as_sf(grid, coords = c("Longitud", "Latitud"), remove = T)
  #plot(sf2)
  
  sf <- as(sf, "Spatial")
  #plot(sf)
  ####Convirtiendo puntos finales a grillas
  sf2 <- as(sf2, "Spatial")
  #plot(sf2)
  
  v_emp_ok = variogram(Hr ~ 1, sf)
  
  library(automap)
  v_mod_ok = autofitVariogram(Hr ~ 1, sf)
  
  
  #plot(v_emp_ok)
  #plot(v_mod_ok)
  
  
  
  g = gstat(formula = Hr ~ 1, model = v_mod_ok$var_model, data = sf)
  z = predict(g,sf2)
  
  
  kriging <- data.frame(Longitud = sf2@coords[,1],Latitud = sf2@coords[,2],Pp = z$var1.pred,Fecha= unique(final_table$Fecha)[i])
  
  
  kriging$Pp[kriging$Pp < 0] <- 0
  kriging$Pp[kriging$Pp > 100] <- 100
  
  interpolado_kriging <- rbind(interpolado_kriging,kriging)
  
  
}

head(interpolado_kriging)
head(interpolado_spline)
head(interpolado)

interpolado <- data.frame(interpolado_spline,Kriging = interpolado_kriging$Pp)
names(interpolado)[3] <- "Spline"
#interpolado <- interpolado

2022-1984+1

library(terra)
###Leyendo datos mensuales chirpsC:\Users\ASUS\Desktop\Tesissss\obs_data
data <- rast("C:/Users/ASUS/Downloads/POWER_Regional_Monthly_1984_2022.nc")
data
data[data <=-99] <- NA

###Leyendo shapefile departamental
shp <- terra::vect("C:/Users/ASUS/Downloads/DEPARTAMENTOS_inei_geogpsperu_suyopomalia.shp")
###Filtrando shapefile cajamarca
#shp_caj <- shp[shp$NOMBDEP == "CAJAMARCA"] 
shp_caj <- shp[shp$NOMBDEP == "MOQUEGUA"] 
###Estableciendo mismo sistema de coordenadas
terra::crs(data) <- terra::crs(shp_caj)

###Recortando datos chirps para cajamarca
data_crop <- terra::crop(data,shp_caj,snap="out",mask=T)

#plot(data_crop[[which(fecha == as.Date("2011-11-01"))]])
# Generar la secuencia del 1 al 130
numeros <- seq(1, 507)


# Filtrar los múltiplos de 13
numeros_sin_multiplos_13 <- numeros[which(numeros %% 13 != 0)]

# Mostrar la secuencia resultante
numeros_sin_multiplos_13

data_crop <- data_crop[[numeros_sin_multiplos_13]]

###Extrayendo fechas
fecha <- seq(as.Date("1984-01-01"),as.Date("2022-12-31"),by="month")


###Cambiando nombre de capas
names(data_crop) <- fecha


###Filtrando periodo 1991 al 2020
data_crop <- data_crop[[which(fecha >= as.Date("2003-01-01") & fecha <= as.Date("2012-12-31"))]]

###Grafica basica de algun mes
plot(data_crop[[119]])

###Paquete para manipulacion de datos grillados en formato dataframe
library(tidyterra)

###Conviertiendo datos grillados a dataframe
df_data <- as_tibble(data_crop,xy=T)

####dataframe a formato largo
df2_data <-df_data %>% 
  pivot_longer(
    cols = !c(x,y), 
    names_to = "Fecha", 
    values_to = "Pp",
    values_drop_na = TRUE)

df2_data$Fecha <- as.Date(df2_data$Fecha)

library(lubridate) ###Paquete para manipulacion de fechas
library(dplyr)  ###Pquete para maenjo de datos

###Calculando climatología mensual
mon_data <- df2_data %>% group_by(x,y,year(Fecha),month(Fecha),month(Fecha,label=T,abbr=F)) %>%
  summarise(Pp= mean(Pp))
names(mon_data) <- c("Longitud","Latitud","Año","Mes","Mes_lab","Chirps")
mon_data$Fecha <- as.Date(paste0(mon_data$Año,"-",mon_data$Mes,"-",01))

#####Ajustando coordenadas
mon_data$Longitud <- round(mon_data$Longitud,3)
mon_data$Latitud <- round(mon_data$Latitud,3)

#unique(mon_data$Longitud)
#unique(mon_data$Latitud)
interpolado$Longitud <- round(interpolado$Longitud,3)
interpolado$Latitud <- round(interpolado$Latitud,3)

library(dplyr)
#min(table_f$Pp)
table_f<- mon_data %>% 
  inner_join(interpolado, by = c("Fecha","Longitud","Latitud")) 
View(table_f)
table_f
library(Metrics)
stat_table1 <- na.omit(table_f) %>% group_by(Longitud,Latitud,Mes,Mes_lab) %>% 
  summarise(R = cor(Chirps,Spline,method = "spearman",
                    use = "pairwise.complete.obs"),
            pval = cor.test(Chirps,Spline,method = "spearman",use = "pairwise.complete.obs",exact=F)$p.value,
            RMSE = rmse(Spline, Chirps),
            MAE = mae(Spline, Chirps),
            MSE = mse(Spline, Chirps),
            MAD = mdae(Spline, Chirps),
            BIAS = Fgmutils::bias(Spline, Chirps))

stat_table2 <- table_f %>% group_by(Longitud,Latitud,Mes,Mes_lab) %>% 
  summarise(R = cor(Chirps,Kriging,method = "spearman",
                    use = "pairwise.complete.obs"),
            pval = cor.test(Chirps,Kriging,method = "spearman",use = "pairwise.complete.obs",exact=F)$p.value,
            RMSE = rmse(Kriging, Chirps),
            MAE = mae(Kriging, Chirps),
            MSE = mse(Kriging, Chirps),
            MAD = mdae(Kriging, Chirps),
            BIAS = Fgmutils::bias(Kriging, Chirps))


stat_table1$Interpolacion = "Spline"
stat_table2$Interpolacion = "Kriging"


shp2 <- shp[shp$NOMBDEP != "MOQUEGUA", ]

unit <- rbind(stat_table1,stat_table2)


setwd("C:/Users/ASUS/Desktop/Moquegua")
#names(metadata)
library(ggspatial)
library(paletteer)

min(unit$RMSE)
max(unit$RMSE)


g1 <- ggplot(data=shp_caj)+
  geom_raster(data=unit,aes(x=Longitud,y=Latitud,fill=RMSE))+
  geom_sf(data=shp2,color="black",size=1,fill="white")+
  geom_sf(fill="transparent",size=1,col="black")+facet_wrap(~Interpolacion*Mes_lab,nrow=4,dir="v")  +
  guides(fill=guide_colorsteps()) +
  scale_fill_gradientn(colours=(paletteer_c("grDevices::Geyser", 30)),breaks=seq(3,30,1)) +theme_bw()+coord_sf(xlim=c(-72,-69.5), ylim=c(-18,-15.5))+theme(plot.title = element_text(hjust = 0.5,face = "bold",size=14))+
  theme_light()+
  theme(panel.grid.major = element_line(color="black", linetype="dashed", size=0.25),legend.key.width  = unit(4,"cm"),
        panel.ontop = TRUE,
        panel.background = element_rect(color=NA, fill=NA),  
        legend.key = element_rect(color="black"),  
        plot.title = element_text(hjust=0.5, face="bold"),
        strip.background = element_rect(fill="#FFFACD",color=NA),
        strip.text = element_text(color="black", face="bold",size=15),
        legend.direction = "horizontal",legend.position = "bottom"
  ) +
  ggtitle("RMSE Kriging y Spline - Nasa Power : Humedad Relativa")+xlab("")+ylab("")+scale_x_continuous(breaks = seq(-72, -69.5, by = 1), labels = scales::number_format(accuracy = 0.1))+
  annotation_scale(location = "bl", width_hint = 0.28) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),height = unit(0.8, "cm"),
                         width = unit(0.8, "cm"),
                         style = north_arrow_fancy_orienteering)+geom_point(data = metadata,aes(y=Latitude,x=Longitude), size = 1.5,shape=12,col="#5D478B")


ggsave(g1,filename=paste0("map_RMSE_spline__kriging_hr.png"), width = 11, height =10,dpi=700)


min(unit$MAE)
max(unit$MAE)


g1 <- ggplot(data=shp_caj)+
  geom_raster(data=unit,aes(x=Longitud,y=Latitud,fill=MAE))+
  geom_sf(data=shp2,color="black",size=1,fill="white")+
  geom_sf(fill="transparent",size=1,col="black")+facet_wrap(~Interpolacion*Mes_lab,nrow=4,dir="v")  +
  guides(fill=guide_colorsteps()) +
  scale_fill_gradientn(colours=(paletteer_c("ggthemes::Sunset-Sunrise Diverging", 30)),breaks=seq(2,29,1)) +theme_bw()+coord_sf(xlim=c(-72,-69.5), ylim=c(-18,-15.5))+theme(plot.title = element_text(hjust = 0.5,face = "bold",size=14))+
  theme_light()+
  theme(panel.grid.major = element_line(color="black", linetype="dashed", size=0.25),legend.key.width  = unit(4,"cm"),
        panel.ontop = TRUE,
        panel.background = element_rect(color=NA, fill=NA),  
        legend.key = element_rect(color="black"),  
        plot.title = element_text(hjust=0.5, face="bold"),
        strip.background = element_rect(fill="#FFFACD",color=NA),
        strip.text = element_text(color="black", face="bold",size=15),
        legend.direction = "horizontal",legend.position = "bottom"
  ) +
  ggtitle("MAE Kriging y Spline - Nasa Power : Humedad Relativa")+xlab("")+ylab("")+scale_x_continuous(breaks = seq(-72, -69.5, by = 1), labels = scales::number_format(accuracy = 0.1))+
  annotation_scale(location = "bl", width_hint = 0.28) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),height = unit(0.8, "cm"),
                         width = unit(0.8, "cm"),
                         style = north_arrow_fancy_orienteering)+geom_point(data = metadata,aes(y=Latitude,x=Longitude), size = 1.5,shape=12,col="#5D478B")


ggsave(g1,filename=paste0("map_MAE_spline__kriging_hr.png"), width = 11, height =10,dpi=700)


min(unit$MAD)
max(unit$MAD)


g1 <- ggplot(data=shp_caj)+
  geom_raster(data=unit,aes(x=Longitud,y=Latitud,fill=MAD))+
  geom_sf(data=shp2,color="black",size=1,fill="white")+
  geom_sf(fill="transparent",size=1,col="black")+facet_wrap(~Interpolacion*Mes_lab,nrow=4,dir="v")  +
  guides(fill=guide_colorsteps()) +
  scale_fill_gradientn(colours=rev(
    paletteer_c("grDevices::GnBu", 30)),breaks=seq(2,29,1)) +theme_bw()+coord_sf(xlim=c(-72,-69.5), ylim=c(-18,-15.5))+theme(plot.title = element_text(hjust = 0.5,face = "bold",size=14))+
  theme_light()+
  theme(panel.grid.major = element_line(color="black", linetype="dashed", size=0.25),legend.key.width  = unit(4,"cm"),
        panel.ontop = TRUE,
        panel.background = element_rect(color=NA, fill=NA),  
        legend.key = element_rect(color="black"),  
        plot.title = element_text(hjust=0.5, face="bold"),
        strip.background = element_rect(fill="#FFFACD",color=NA),
        strip.text = element_text(color="black", face="bold",size=15),
        legend.direction = "horizontal",legend.position = "bottom"
  ) +
  ggtitle("MAD Kriging y Spline - Nasa Power : Humedad Relativa")+xlab("")+ylab("")+scale_x_continuous(breaks = seq(-72, -69.5, by = 1), labels = scales::number_format(accuracy = 0.1))+
  annotation_scale(location = "bl", width_hint = 0.28) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),height = unit(0.8, "cm"),
                         width = unit(0.8, "cm"),
                         style = north_arrow_fancy_orienteering)+geom_point(data = metadata,aes(y=Latitude,x=Longitude), size = 1.5,shape=12,col="#5D478B")


ggsave(g1,filename=paste0("map_MAd_spline__kriging_hr.png"), width = 11, height =10,dpi=700)


#########################################33333


min(unit$MSE)
max(unit$MSE)


g1 <- ggplot(data=shp_caj)+
  geom_raster(data=unit,aes(x=Longitud,y=Latitud,fill=MSE))+
  geom_sf(data=shp2,color="black",size=1,fill="white")+
  geom_sf(fill="transparent",size=1,col="black")+facet_wrap(~Interpolacion*Mes_lab,nrow=4,dir="v")  +
  guides(fill=guide_colorsteps()) +
  scale_fill_gradientn(colours=rev(
    paletteer_c("grDevices::Earth", 30)),breaks=seq(10,900,50)) +theme_bw()+coord_sf(xlim=c(-72,-69.5), ylim=c(-18,-16))+theme(plot.title = element_text(hjust = 0.5,face = "bold",size=14))+
  theme_light()+
  theme(panel.grid.major = element_line(color="black", linetype="dashed", size=0.25),legend.key.width  = unit(4.5,"cm"),
        panel.ontop = TRUE,
        panel.background = element_rect(color=NA, fill=NA),  
        legend.key = element_rect(color="black"),  
        plot.title = element_text(hjust=0.5, face="bold"),
        strip.background = element_rect(fill="#FFFACD",color=NA),
        strip.text = element_text(color="black", face="bold",size=15),
        legend.direction = "horizontal",legend.position = "bottom"
  ) +
  ggtitle("MSE Kriging y Spline - Nasa Power : Humedad Relativa")+xlab("")+ylab("")+scale_x_continuous(breaks = seq(-72, -69.5, by = 1), labels = scales::number_format(accuracy = 0.1))+
  annotation_scale(location = "bl", width_hint = 0.28) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),height = unit(0.8, "cm"),
                         width = unit(0.8, "cm"),
                         style = north_arrow_fancy_orienteering)+geom_point(data = metadata,aes(y=Latitude,x=Longitude), size = 1.5,shape=12,col="#5D478B")


ggsave(g1,filename=paste0("map_MSE_spline__kriging_hr.png"), width = 10.5, height =10,dpi=700)
