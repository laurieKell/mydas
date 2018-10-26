## Datasets    

The main dataset on commercial landings are from the [STECF](https://stecf.jrc.ec.europa.eu/dd/effort/graphs-quarter), as new data become available the database will be updated.

## Summary

The data are in a postgres database and can be summarised and analysed using R, see 

[pdf](https://github.com/laurieKell/mydas/blob/master/tasks/task1/R/stockprioritisation.pdf), 
[source](https://github.com/laurieKell/mydas/blob/master/tasks/task1/R/stockprioritisation.Rmd)

## PostgreSQL access

[database description](https://github.com/tunafish72/mydas/blob/master/Databasedescripton.pdf),
host = 'postgresql-seascope.csffkpr9jjjn.eu-west-2.rds.amazonaws.com'
dbname='mydasDev'
port = 5432
user = 'MydasApplication'
password = 'gmit2017!'
## Shiny application

[shiny application](http://35.177.86.42:3838/mydas/), 
user: mydas, pwd:gmit1
