gfr_anomaly <- function(gfr){
    # Generate decimal dates
    gfr$gfrdate <- with(gfr, ymd_to_decimal(paste0(year,'-',month,'-',day)) )

    # **********************************************************************************************************************
    # ** GFR ANOMOLIES, ADDED JULY 2021																					**
    # ** GFR, GFRdate, and RRT date anomolies: The following hardcode changes to U25eGFR value/method and or GFR dates	**
    # ** are the result of an extensive review of participant longitudinal GFR profiles and emails to the CCCs.			**
    # ** More information regarding these anomolies is located at ~cpierce\KIDMAC\GFR\anomolies\							**
    # **																													**
    # ** 20AUG2021, CBP: Code updated to include specific problem in the condition of the "if..do" statement 	
    # **
    # ** 21OCT2021, Code transported to R (questions & bug reports to sunjae@jhu.edu)
    # *********************************************************************************************************************;

    #  hard core corrections for analysis as of 2021.07.23;

    # patterns of U25eGFR in transplanted
    gfr[gfr$kid==102008 & gfr$visit==101 & gfr$krtstatus==0, ]$gfrdate <- 2015.87
    gfr[gfr$kid==102008 & gfr$visit==101 & gfr$krtstatus==0, ]$krtstatus <- 21

    gfr[gfr$kid==105009 & gfr$visit==70 & gfr$cyc_ifcc==0.70, ]$u25egfr <- 51.8
    gfr[gfr$kid==105009 & gfr$visit==70 & gfr$cyc_ifcc==0.70, ]$u25egfrmeth <- 2
    gfr[gfr$kid==105009 & gfr$visit==70 & gfr$cyc_ifcc==0.70, ]$cyc_ifcc <- NA

    gfr <- gfr[-which(gfr$kid==105014 & gfr$visit==130), ]

    gfr[gfr$kid==106005 & gfr$visit==30 & gfr$scr==2.90, ]$u25egfr <- NA
    gfr[gfr$kid==106005 & gfr$visit==30 & gfr$scr==2.90, ]$u25egfrmeth <- NA
    gfr[gfr$kid==106005 & gfr$visit==30 & gfr$scr==2.90, ]$e2012gfr <- NA

    gfr[gfr$kid==115011 & gfr$visit==21 & gfr$krtstatus==21, ]$gfrdate <- 2009.68
    gfr[gfr$kid==115011 & gfr$visit==21 & gfr$krtstatus==21, ]$krtstatus <- 0

    gfr <- gfr[-which(gfr$kid==151003 & gfr$visit==71), ]

    gfr[gfr$kid==168008 & gfr$visit==50 & gfr$cyc_ifcc==0.20, ]$u25egfr <- 33.6
    gfr[gfr$kid==168008 & gfr$visit==50 & gfr$cyc_ifcc==0.20, ]$u25egfrmeth <- 2
    gfr[gfr$kid==168008 & gfr$visit==50 & gfr$cyc_ifcc==0.20, ]$cyc_ifcc <- NA

    gfr[gfr$kid==210011 & gfr$visit==41  & gfr$krtstatus==21, ]$gfrdate <- 2014.24
    gfr[gfr$kid==210011 & gfr$visit==41  & gfr$krtstatus==21, ]$krtstatus <- 0


    # *from looking at tails of U25eGFR using all child-visits by krtstatus;
    gfr[gfr$kid==112007 & gfr$visit==10  & gfr$cyc_ifcc==0.32, ]$u25egfr <- 53.4
    gfr[gfr$kid==112007 & gfr$visit==10  & gfr$cyc_ifcc==0.32, ]$u25egfrmeth <- 2
    gfr[gfr$kid==112007 & gfr$visit==10  & gfr$cyc_ifcc==0.32, ]$cyc_ifcc <- NA


    # extreme values of scr and cystatin;
    gfr[gfr$kid==174013 & gfr$visit==10  & gfr$scr==0.20, ]$u25egfr <- 52.45
    gfr[gfr$kid==174013 & gfr$visit==10  & gfr$scr==0.20, ]$u25egfrmeth <- 3
    gfr[gfr$kid==174013 & gfr$visit==10  & gfr$scr==0.20, ]$scr <- NA

    gfr[gfr$kid==175011 & gfr$visit==40  & gfr$cyc_ifcc==3.03, ]$u25egfr <- 75.06
    gfr[gfr$kid==175011 & gfr$visit==40  & gfr$cyc_ifcc==3.03, ]$u25egfrmeth <- 2
    gfr[gfr$kid==175011 & gfr$visit==40  & gfr$cyc_ifcc==3.03, ]$cyc_ifcc <- NA


    # EXTREME VALUES OF studentized residuals of U25eGFR (within children based on random regression with logUPCR modifying level and slope);
    gfr[gfr$kid==128004 & gfr$visit==30  & gfr$scr==0.70, ]$u25egfr <- 38.5
    gfr[gfr$kid==128004 & gfr$visit==30  & gfr$scr==0.70, ]$u25egfrmeth <- 3
    gfr[gfr$kid==128004 & gfr$visit==30  & gfr$scr==0.70, ]$scr <- NA

    gfr[gfr$kid==150004 & gfr$visit==70  & gfr$cyc_ifcc==0.64, ]$u25egfr <- 37.8
    gfr[gfr$kid==150004 & gfr$visit==70  & gfr$cyc_ifcc==0.64, ]$u25egfrmeth <- 2
    gfr[gfr$kid==150004 & gfr$visit==70  & gfr$cyc_ifcc==0.64, ]$cyc_ifcc <- NA

    gfr[gfr$kid==161008 & gfr$visit==60  & gfr$cyc_ifcc==0.37, ]$u25egfr <- 33.5
    gfr[gfr$kid==161008 & gfr$visit==60  & gfr$cyc_ifcc==0.37, ]$u25egfrmeth <- 2
    gfr[gfr$kid==161008 & gfr$visit==60  & gfr$cyc_ifcc==0.37, ]$cyc_ifcc <- NA

    # ** END ANOMOLIES *****************************************************************************************************;

    return(gfr)
}
