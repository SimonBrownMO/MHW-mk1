# source("digest_mass_files.R") 

library(ncdf4)
library(PCICt)
source("/home/users/simon.brown/extremes/quantile_modelling/NSQM/libs/general_stats.R")

st.out <- "../DATA/mass_dump/../1980-2100-Region_mean_timeseries_Daily_r001i1p00000.RData"

flist_m <- "19800101_Region_mean_timeseries_Daily.nc"
flist_m <- list.files("../DATA/mass_dump/", pattern="*Region_mean_timeseries_Daily.nc", full.names=TRUE) 

m.sst    <- NULL #  $ reg_sst_ave    : chr "Regional mean sea surface temperature"
m.mld    <- NULL #  $ reg_mld_ave    : chr "Regional mean mixed layer depth"
m.time   <- NULL #  $ time_centered  : chr "time_centered"

ireg  <- 11 # required region and mask id
imask <- 3  # or 6?
st1   <- "O9_ap977_ar095_au084_r001i1p00000"

for(i in seq_along(flist_m)) {
    nc1      <- nc_open(flist_m[i])
    vlist    <- nc1$var
    if(i==1) {
        lname    <- lapply(vlist, '[[', "longname")
        r1       <- ncvar_get(nc1, 'reg_id')
        m1       <- ncvar_get(nc1, 'mask_id')
    }
    x1       <- ncvar_get(nc1, 'reg_sst_ave')
    d1       <- ncvar_get(nc1, 'reg_mld_ave')
    # t1       <- ncvar_get(nc1, 'time_centered')
    m.sst    <- cbind(m.sst, x1) 
    m.mld    <- cbind(m.mld, d1) 

    time1    <- vlist$time_centered$dim[[1]]
    origin   <- tail(unlist(strsplit(time1$units,'seconds since ')),1)
    m.time   <- c(m.time,time1$vals)  # time1$calendar  "360_day"

    nc_close(nc1)
    if(i %% 100 == 0)cat("Done with file ",i," of ",length(flist_m),"\n")
}

m.time2 <- as.PCICt( m.time, cal=time1$calendar ,origin=origin)
m.years <- as.integer(format(m.time2,"%Y"))
m.doy   <- as.integer(strftime(m.time2, format = "%j"))
m.date  <- m.years + (m.doy-0.5)/(unlist(lapply(m.years,days_in_year))) 

i0 <- which(r1[,1]==ireg & m1[,1]==imask) 
pcht <- paste(st1,"Region",ireg,"mask",imask )
plot(as.POSIXct(m.time2), m.sst[i0,], pch=46, main=pcht)
plot(as.POSIXct(m.time2), m.sst[i0,], ty='l', main=pcht)

m.sst <- list(sst=m.sst, mld=m.mld, time=m.date, years=m.years, date=m.time2, doy=m.doy)

save(file=st.out, m.sst)    
