ReadIntercatch <- function(file){
  IC <- read.table(file ,sep=",", col.names=as.character(1:33), fill=T)
  HI <- subset(IC,X1=='HI')[,1:12]
  names(HI) <- c("RecordType", "Country", "Year", "SeasonType", "Season", "Fleet",
                 "AreaType", "FishingArea", "DepthRange", "UnitEffort", "Effort",
                 "AreaQualifier")
  SI <- subset(IC,X1=='SI')[,1:24]
  names(SI) <- c("RecordType", "Country", "Year", "SeasonType", "Season", "Fleet",
                 "AreaType", "FishingArea", "DepthRange", "Species", "Stock", "CatchCategory",
                 "ReportingCategory", "DataToFrom", "Usage", "SamplesOrigin", "Qualityflag",
                 "UnitCaton", "CATON", "OffLandings", "VarCaton", "InfoFleet",
                 "InfoStockCoordinator", "InfoGeneral")
  if(sum(!SI$UnitCaton %in% c("t","kg"))>0)
    stop("Invalid UnitCaton: only 't' or 'kg' allowed")
  x <- ifelse(SI$UnitCaton == "kg", 0.001, 1)
  SI$CATON <- SI$CATON * x
  SI$UnitCaton <- "t"
  SD <- subset(IC,X1=='SD')[,1:33]
  names(SD) <- c("RecordType", "Country", "Year", "SeasonType", "Season", "Fleet",
                 "AreaType", "FishingArea", "DepthRange", "Species", "Stock", "CatchCategory",
                 "ReportingCategory", "Sex", "CANUMtype", "AgeLength", "PlusGroup",
                 "SampledCatch", "NumSamplesLngt", "NumLngtMeas", "NumSamplesAge", "NumAgeMeas",
                 "unitMeanWeight", "unitCANUM", "UnitAgeOrLength", "UnitMeanLength", "Maturity",
                 "NumberCaught", "MeanWeight", "MeanLength", "VarNumLanded", "VarWeightLanded",
                 "VarLngthLanded")
  if(sum(!SD$unitMeanWeight %in% c("g","kg"))>0)
    stop("Invalid unitMeanWeight: only 'g' or 'kg' allowed")
  x <- ifelse(SD$unitMeanWeight == "g", 0.001, 1)
  SD$MeanWeight <- SD$MeanWeight * x
  SD$unitMeanWeight <- 'kg'
  if(sum(!SD$unitCANUM %in% c("k","m","n"))>0)
    stop("Invalid unitCANUM: only 'k', 'm' or 'n' allowed")
  x <- ifelse(SD$unitCANUM == "m", 0.001, ifelse(SD$unitCANUM == "n", 1000, 1))
  SD$NumberCaught <- SD$NumberCaught * x
  SD$unitCANUM <- 'k'
  if(sum(!SD$UnitAgeOrLength %in% c("cm","mm"))>0)
    stop("Invalid UnitAgeOrLength: only 'cm' or 'mm' allowed for this datacall")
  x <- ifelse(SD$UnitAgeOrLength == "mm", 0.1, 1)
  SD$AgeLength <- as.numeric(as.character(SD$AgeLength)) * x
  SD$unitMeanWeight <- 'cm'
  return(list(HI,SI,SD))
}

