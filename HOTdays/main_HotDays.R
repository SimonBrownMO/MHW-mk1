# source("main_HotDays.R")

# sbatch --time=210 --mem=30000 -o $PWD/main.R.rout --wrap="Rscript $PWD/main.R"
# sbatch --qos=long --time=4000 --mem=30000 -o $PWD/main.R.rout --wrap="Rscript $PWD/main.R"
# 2023-06-02 start to doJointMakeStationary.R took 2.4h

cat("###############################################################################",cr)
sysdate0   <- Sys.time()
stsysdate0 <- print(sysdate0)
cat("starting main_HotDays.R",cr)
cat("###############################################################################",cr,cr)

cargs <- commandArgs()

TESTINBATCH <- grep('-f',cargs)
if(!any(TESTINBATCH)) {
    cat("NOT IN a BATCH environment",cr)
    ### variables that would be passed in by script
    # MSref.file <- "/data/users/lowen/extremes/heatwaves/HadUKGrid/dur-clim/CPM5km/v_1/UK/15/MakeStationary/MSref/ukgd_cpm85_5k_x147y44.MSref.2024-10-30-122429.RData"
    MSref.file <- "/data/users/laura.owen/extremes/heatwaves/HadUKGrid/dur-clim/CPM5km/v2/UK/01/MakeStationary/MSref/ukgd_cpm85_5k_x147y46.MSref.2025-04-25.RData"
} else {
    cat("In a BATCH environment",cr)
    cargs      <- commandArgs(trailingOnly = TRUE)
    MSref.file <- cargs[[1]]
}

cat(cr,"MSref.file",cr)
cat(MSref.file,cr,cr)
st.pwd    <- system("pwd", intern=TRUE)

### setup variables #####################################
    source(paste(st.pwd,"/../setup_all.R",sep=''))
    # MSref.file loaded in setup_HotDays.R
    source(paste(st.pwd,"/setup_HotDays.R",sep=''))
    cat("HotDays main: Completed setup_HotDays.R",cr,"#########################",cr,cr)
    cat("HotDays main: st_Term",HDconfig$files$st_Term,cr)
    # reload_HD_ref(HDconfig$files$st_HDconfig)

# load('/data/users/hadsx/extremes/heatwaves/HadUKGrid/joint/v_j2c/CPM/HotDays/vHD2c/ukgd_cpm85_5k_x146y44.HDref.2024-09-12-164818.RData',verb=TRUE)
# load(HDconfig$files$msfiles$st_msgpd,verb=TRUE)
# load(HDconfig$files$msfiles$st_msdata01,verb=TRUE)

 if (!file.exists(HDconfig$files$st_Term)) {
    ### extract events #####################################
        source(paste(st.pwd,"/doExtractEvents.R",sep=''))
        iobs <- which(data01$isobs==1)
        if(DOEVENTSPLOT) fn_plotEvents(events01$obs, data01[iobs,], nplot=3, savefile=sub('Events/ukgd','Events/plots/ukgd',sub('.RData','.pdf',st_events)) )
        cat("HotDays main: Completed doExtractEvents.R",cr,"#########################",cr,cr)

    ### doFitInitEvent.R ##################################### MEDIUM-SLOW
        source(paste(st.pwd,"/doFitInitEvent.R",sep=''))
        cat("HotDays main: Completed doFitInitEvent.R",cr,"#########################",cr,cr)

    ### doFitInitValue.R ##################################### FAST
        source(paste(st.pwd,"/doFitInitValue.R",sep=''))
        cat("HotDays main: Completed doFitInitValue.R",cr,"#########################",cr,cr)

    ### doFitHt.R ##################################### FAST
        source(paste(st.pwd,"/doFitHt.R",sep=''))
        cat("HotDays main: Completed doFitHt.R",cr,"#########################",cr,cr)

    ### doFitTerm.R ##################################### FAST
        source(paste(st.pwd,"/doFitTerm.R",sep=''))
        # save_metadata(fin=paste(st.pwd,"/"),fout=stmetadata)       ### need to add appending
        cat("HotDays main: Completed doFitTerm.R",cr,"#########################",cr,cr)


        cat("st_Term",HDconfig$files$st_Term,cr,cr)

    cat("###############################################################################",cr)
    sysdate    <- Sys.time()
    stsysdate  <- print(sysdate)

    cat("Time taken",cr)
    print(sysdate-sysdate0)
    cat("All Done main_HotDays.R",cr)
    cat("###############################################################################",cr,cr)

 } else {

    cat("File exists",cr,HDconfig$files$st_Term,cr)

 }
#
