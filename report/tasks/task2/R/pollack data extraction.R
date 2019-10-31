library(RODBC)

Q1 <- 
"
select 'Survey' collate Latin1_General_CI_AS as dataType
   ,a.fldCruiseName as dataSource
   ,datepart(yy,fldDateTimeHaul) as [year]
   ,datepart(ww,fldDateTimeHaul) as [week]
   ,a.fldFishSex collate Latin1_General_CI_AS as sex
   ,b.fldICESRectangle collate Latin1_General_CI_AS as area
   ,a.fldFishLength as lenCls
   ,a.fldResult1 collate Latin1_General_CI_AS as age
   ,a.fldFishWholeWeight as indWt
   ,convert(char,a.fldFishMaturity) collate Latin1_General_CI_AS as matStage
from fss_survey..tblDataBiologicalSamples a
   join fss_survey..tblDataStationLogs b
      on a.fldCruiseName = b.fldCruiseName
      and a.fldCruiseStationNumber = b.fldCruiseStationNumber
where a.fldMainSpeciesCode = 'POL'

union all 

select 'Landings' as dataType
   ,'Lan' + convert(varchar,datepart(yy,b.SampleDate)) as dataSource
   ,datepart(yy,b.SampleDate) as [year]
   ,datepart(ww,b.SampleDate) as [week]
   ,case when g.Sex in ('M','F') then upper(g.sex) end as sex
   ,c.areas as area
	,cast(g.FishLength*10 as int) as lenCls
	,g.Age as age
	,case when e.PresentationTypeID='Gutted' then g.FishWeight*f.CONVERSION_FACTOR else g.FishWeight end as indWt
	,case when g.MaturityCode = 'none' then NULL else convert(char,g.MaturityCode) end as matStage
from vmfsssql02.stockman2015.dbo.SamplingTrip a
   join vmfsssql02.stockman2015.dbo.SampleHeader b
      on a.TripID = b.TripID
   left join IntercatchAreasLut c
      on b.ICESDivision = c.ICESDIV
   left join vmfsssql02.stockman2015.dbo.MD_Port d
      on a.LandingPlaceID = d.port_id
   join vmfsssql02.stockman2015.dbo.SampleDetail e
      on b.SampleHeaderID = e.SampleHeaderID
   join vmfsssql02.stockman2015.dbo.SPECIESLUT f
      on e.SpeciesID = f.SPECIESLUT_ID
   join vmfsssql02.stockman2015.dbo.SampleAgeData g
      on e.SampleDetailID = g.SampleDetailID
where f.FISHLATINNAME = 'Pollachius pollachius'

union all 
select 'Discards' as dataType
   ,'Dis' + convert(varchar,datepart(yy, a.[Departure Date]))
   ,datepart(yy, a.[Departure Date]) as [year]
   ,datepart(ww, a.[Departure Date]) as [week]
   ,case when f.Sex = 'F' then 'F'
      when f.sex = 'M' then 'M' end as sex
   ,h.ices_code as area
   ,convert(int,f.[Length]) * 10 lenCls
  ,f.[Age] as age
   ,case when convert(int,f.[Weight]) between 1 and 99999 then convert(int,f.[Weight]) end as indWt
   ,convert(char,maturity) as MatStage
from Discards..CRUISe a
   left join Discards..FSSPort b
      on a.ArrPort_id = b.port_id
   left join Discards..FSSVessel c
      on a.fss_vessel_id = c.vessel_id
   join Discards..HAUl d
      on a.[Cruise Code] = d.[Cruise Code]
   join Discards..[SAMPLE-HEADEr] e
      on d.[Cruise Code] = e.[Cruise code]
      and d.Haul = e.[Haul code]
   join Discards..[SAMPLe] f
      on e.sample_header_id = f.sample_header_id
   join Discards..SPECIES g
      on f.Species = g.Species
   left join Discards..FSSICES h
      on d.ICES_id = h.ices_id
where a.[cruise code] not like 'eco%'
   and f.[Weight] > 0 
   and g.[Scientific Name] = 'Pollachius pollachius'
"

Q2 <- "
 select datepart(yy,b.SampleDate) as [year]
   ,b.SampleHeaderID
   ,case when substring(b.GearCode,1,3) in ('GN ','GND','GNS','GTR') then 'GNS'
      when substring(b.GearCode,1,3) in ('OTB','PTB','TRR','TWR') then 'OTB'
      when substring(b.GearCode,1,3) in ('PS ','SB ','SDN','SSC') then 'SSC'
      when substring(b.GearCode,1,3) in ('TBB') then 'TBB'
      else 'OTH' end as gear
	,cast(g.LengthClass*10 as int) as lenCls
	,cast(sum(g.NumberAtLength) as int) as lenNum
from vmfsssql02.stockman2015.dbo.SamplingTrip a
   join vmfsssql02.stockman2015.dbo.SampleHeader b
      on a.TripID = b.TripID
   left join vmfsssql02.stockman2015.dbo.MD_Port d
      on a.LandingPlaceID = d.port_id
   join vmfsssql02.stockman2015.dbo.SampleDetail e
      on b.SampleHeaderID = e.SampleHeaderID
   join vmfsssql02.stockman2015.dbo.SPECIESLUT f
      on e.SpeciesID = f.SPECIESLUT_ID
   join vmfsssql02.stockman2015.dbo.SampleLengthData g
      on e.SampleDetailID = g.SampleDetailID
where b.CATCHTYPEID = 1
   and g.NumberAtLength > 0 
   and UseAsMeasuredOnlyAsWell = 1
   and f.FISHLATINNAME = 'Pollachius pollachius'
group by datepart(yy,b.SampleDate)
   ,b.SampleHeaderID
   ,case when substring(b.GearCode,1,3) in ('GN ','GND','GNS','GTR') then 'GNS'
      when substring(b.GearCode,1,3) in ('OTB','PTB','TRR','TWR') then 'OTB'
      when substring(b.GearCode,1,3) in ('PS ','SB ','SDN','SSC') then 'SSC'
      when substring(b.GearCode,1,3) in ('TBB') then 'TBB'
      else 'OTH' end
   ,cast(g.LengthClass*10 as int)

"

channel <- odbcDriverConnect("Driver=SQL Server; Server=VMFSSDEV02; 
                             Database=Stockman2015_views")
 bio <- sqlQuery(channel,Q1)
 len <- sqlQuery(channel,Q2)
close(channel)


source('\\\\GalwayFS03\\Fishdata\\Data for ICESWG\\2018\\Intercatch functions\\ReadIntercatch.R')

ic2016 <- ReadIntercatch(
  '\\\\GalwayFS03\\Fishdata\\Data for ICESWG\\2017\\WGCSE\\pol.27.67\\intercatch_pol.27.67.csv')
ic2017 <- ReadIntercatch(
  '\\\\GalwayFS03\\Fishdata\\Data for ICESWG\\2018\\WGCSE\\pol.27.67\\intercatch_pol.27.67_landings.csv')

lan <- rbind(ic2016[[3]],ic2017[[3]])

save(bio,len,lan,file='pollack.RData')

     