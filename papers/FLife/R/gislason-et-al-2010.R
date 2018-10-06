#' Data from [Gislason et al 2010](http://onlinelibrary.wiley.com/doi/10.1111/j.1467-2979.2009.00350.x/full)
#'
#' Gislason, H., Daan, N., Rice, J. C., & Pope, J. G. (2010). 
#' Size, growth, temperature and the natural mortality of marine fish. 
#' Fish and Fisheries, 11(2), 149-158.)

GislasonEtAl2010Data <- object('GislasonEtAl2010Data')

GislasonEtAl2010Data$read <- function(directory='.'){
  read.table(file.path(directory,"gislason-et-al-2010.tsv"),header=T,sep='\t')
}

