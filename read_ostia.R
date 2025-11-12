# source("read_ostia.R")

xx      <- read.csv('DATA/database_ofrd-mopa.csv', head=TRUE)
mhw.date <- strptime(xx$date, format="%Y-%m-%d", tz='GMT')

 plot(mhw.date,xx$ostia_cdr_NWS_sst,ty='l')
lines(mhw.date,xx$ostia_nrt_NWS_sst,col=2)

nws <- xx$ostia_cdr_NWS_sst
nws[is.na(nws)] <- xx$ostia_nrt_NWS_sst[is.na(nws)]
 plot(mhw.date,nws,ty='l')

nws_a <- xx$ostia_cdr_NWS_anomaly
nws_a[is.na(nws_a)] <- xx$ostia_nrt_NWS_anomaly[is.na(nws_a)]
 plot(mhw.date,nws_a,ty='l')

###  extract ssts   ##########################################
regions <- c("ostia_cdr_Global_.90S_to_90N._sst",
             "ostia_cdr_Global_.60S_to_60N._sst",
             "ostia_cdr_North_Atlantic_.CR._sst",
             "ostia_cdr_North_Atlantic_.MO._sst",
             "ostia_cdr_UKV_sst", 
             "ostia_cdr_NWS_sst",      
             "ostia_cdr_North_Sea_sst",
             "ostia_cdr_N_Hem_.20N_to_60N._sst",
             "ostia_cdr_Tropics_.20S_to_20N._sst",
             "ostia_cdr_S_Hem_.60S_to_20S._sst")
regions.short <- gsub("ostia_cdr_","",regions)
regions.short <- gsub("_sst","",regions.short)   

l.mhw <- list()
for(i in seq_along(regions)) {

    l.mhw[[i]]       <- xx[,regions[i]]
    i.na             <- is.na(l.mhw[[i]])
    l.mhw[[i]][i.na] <- xx[i.na, gsub("cdr","nrt",regions[i])]
}
names(l.mhw) <- regions.short
l.mhw$date <- mhw.date

plot(mhw.date,l.mhw[['NWS']],ty='l')

###  add anomalies  ##########################################
regions <- c("ostia_cdr_Global_.90S_to_90N._anomaly",
             "ostia_cdr_Global_.60S_to_60N._anomaly",
             "ostia_cdr_North_Atlantic_.CR._anomaly",
             "ostia_cdr_North_Atlantic_.MO._anomaly",
             "ostia_cdr_UKV_anomaly", 
             "ostia_cdr_NWS_anomaly",      
             "ostia_cdr_North_Sea_anomaly",
             "ostia_cdr_N_Hem_.20N_to_60N._anomaly",
             "ostia_cdr_Tropics_.20S_to_20N._anomaly",
             "ostia_cdr_S_Hem_.60S_to_20S._anomaly")
regions.short <- gsub("ostia_cdr_","",regions)
regions.short <- gsub("_anomaly","",regions.short)   

l.mhw_an <- list()
for(i in seq_along(regions)) {

    l.mhw_an[[i]]       <- xx[,regions[i]]
    i.na             <- is.na(l.mhw_an[[i]])
    l.mhw_an[[i]][i.na] <- xx[i.na, gsub("cdr","nrt",regions[i])]
}
names(l.mhw_an) <- regions.short
l.mhw_an$date <- mhw.date

plot(mhw.date,l.mhw_an[['NWS']],ty='l')

###  check climatology ##########################################
plot(mhw.date,l.mhw[['NWS']]-l.mhw_an[['NWS']],ty='l')


save(file="DATA/ostia_cdr_nrt_regions.RData", l.mhw, l.mhw_an)

#
