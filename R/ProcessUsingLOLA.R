ProcessUsingLOLA <- function() {
  #dbPath = system.file("extdata", "hg19", package="LOLA")

  dbPath = "/media/H_driver/2016/Yang/Results/data/groups/lab_bock/shared/resources/regions/LOLACore/mm10"

  regionDB = loadRegionDB(dbLocation=dbPath)

  #read triple overlap peaks
  regionSetA = readBed("/media/H_driver/2016/Yang/Results/WT_triple_overlap.bed")

  #read each peak
  regionSet.WT.ASXL1=readBed("/media/H_driver/2016/Yang/BedFromPeng/10_WT_BM_ASXL1_peaks_from_PeakCall4Yang.bed")
  regionSet.WT.SMC1=readBed("/media/H_driver/2016/Yang/BedFromPeng/11_2470IUPUI_WT_BM_SMC1_peaks_from_PeakCall4Yang.bed")
  regionSet.WT.Rad21=readBed("/media/H_driver/2016/Yang/BedFromPeng/13_2470IUPUI_WT_BM_Rad21_peaks_from_PeakCall4Yang.bed")


  db.841<-regionDB$regionGRL

  userSets = GRangesList(regionSet.WT.ASXL1,regionSet.WT.SMC1,regionSet.WT.Rad21)

  db.841.2<-c(db.841,userSets)

  userUniverse<-buildRestrictedUniverse(userSets)

  userUniverse.844<-buildRestrictedUniverse(db.841.2)

  re.SMC1A.RAD21.only<-makeVennDiagram(re.out[c(2,4)],NameOfPeaks=c("SMC1A", "RAD21"),totalTest=35000)

  ol.SMC1A.RAD21.only <- findOverlapsOfPeaks(re.out[c(2,4)])

  ol.SMC1A.RAD21.only.2<-ol.SMC1A.RAD21.only$peaklist$`11_2470IUPUI_WT_BM_SMC1_peaks.bed///13_2470IUPUI_WT_BM_Rad21_peaks.bed`

  checkUniverseAppropriateness(ol.SMC1A.RAD21.only.2,userUniverse.844, cores = 1,fast = FALSE)

  locResults2 = runLOLA(userSets,userUniverse,regionDB, cores=1)
  locResults = runLOLA(regionSetA,userUniverse,regionDB, cores=1)

  #SMC1A and RAD21 overlap



  locResults.SMC1A.RAD21 = runLOLA(ol.SMC1A.RAD21.only.2,userUniverse,regionDB, cores=4)

  head(locResults.SMC1A.RAD21)

  pValue<-as.data.frame(10^(-locResults.SMC1A.RAD21$pValueLog))

  temp<-as.data.frame(locResults.SMC1A.RAD21)

  locResults.SMC1A.RAD21.with.p<-cbind(temp[,1:3],pValue,temp[,4:23])

  colnames(locResults.SMC1A.RAD21.with.p)[4]="pValue"

  #It takes extensive computation resource, be carefully to run
  locResults.SMC1A.RAD21.based.large.universe = runLOLA(ol.SMC1A.RAD21.only.2,userUniverse.844,regionDB, cores=4)

  write.table(locResults.SMC1A.RAD21.with.p,
              file=paste0(out.dir.name,"SMC1A.RAD21.with.Other.Antibody.csv"),sep=",",
              quote = FALSE,row.names = FALSE,col.names = TRUE)

  write.table(locResults,file=paste0(out.dir.name,"Other.csv"),sep=",",quote = FALSE,row.names = FALSE,col.names = TRUE)

  listRegionSets(regionDB)
  extractEnrichmentOverlaps(locResults,regionSetA, regionDB)

  loadCaches("/media/H_driver/2016/Yang/Results/LOLACoreCaches_latest.tgz")

  save.image(file =paste0(out.dir.name,"re_search.RData"))

  # regionSetB = readBed("lola_vignette_data/setB_100.bed")
  # regionSetC = readBed("lola_vignette_data/setC_100.bed")
  # activeDHS = readBed("lola_vignette_data/activeDHS_universe.bed")
  # data("sample_universe", package="LOLA")
  #
  # dbPath = system.file("extdata", "hg19", package="LOLA")
  # regionDB = loadRegionDB(dbLocation=dbPath)
  # data("sample_universe", package="LOLA")
  # data("sample_input", package="LOLA")
  #
  # getRegionSet(regionDB, collections="ucsc_example", filenames="vistaEnhancers.bed")
  # getRegionSet(dbPath, collections="ucsc_example", filenames="vistaEnhancers.bed")
  #
  # res = runLOLA(userSets, userUniverse, regionDB, cores=1)
  # locResult = res[2,]
}
