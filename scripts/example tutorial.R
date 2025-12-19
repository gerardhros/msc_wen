# example script for learning R and doing analysis
# Hiiiiiiii
# add a second line example

  # install library data.table
  require(data.table)
  
# objects in R

  # most simple one is a vector
  d = 5
  
  # vector with a length
  d = c(5,10,6,18)
  
  # data.frame (is a table)
  d = data.frame(c1 = LETTERS[1:5],c2 = 1:5,c3 = c(T,T,T,F,F),c4 = c(1.5,2.5,3.2,4,3))
  
  # data types: character, number, integer, boolean
  
  # matrix
  d = as.matrix(d)
  
  # data.table
  d = data.table(c1 = LETTERS[1:5],c2 = 1:5,c3 = c(T,T,T,F,F),c4 = c(1.5,2.5,3.2,4,3))
  
  # the data.frame calculation
  d$c5 <- d$c4 + 5 * 2 + d$c2 / 2 + sqrt(d$c4)
  
  # the same calculation in data.table
  d[, c5 := c4 + 5 * 2 + c2/2 + sqrt(c4)]
  
  # for data.table you need to know a few assumptions
  d[c2 < 3 & c1 =='A'] # filtering
  d[,c5 := c4 * 25]
  d[c2 < 3 & c1 =='A', c6 := c4* sqrt(c4)] # filtering
  
  # doing calculations on groups of data
  d = data.table(c1 = LETTERS[c(1:10,1:10)])
  d[,c2 := rnorm(.N,mean=25,sd=12)]
  d[,mean(c2)]
  d[,list(mean = mean(c2),sd = sd(c2)),by=c1]
  d$c1[c(1,5,8)] <- NA_real_
  
  d = data.table(c1 = LETTERS[c(1:10,1:10)])
  d[,c2 := rnorm(.N,mean=25,sd=12)]
  
  e = data.table(c1n = LETTERS[1:8],c8 = 5, c9 = 'hello',c10= T)
  
  f <- merge(d,e,by.x = 'c1',by.y = 'c1n',all.x=TRUE)
  
  # here is the cheatsheet for data.table
  # https://raw.githubusercontent.com/rstudio/cheatsheets/master/datatable.pdf
  
  # for plotting using ggplot2
  # for nice plotting you can make use of patchwork
  
  # spatial analysis and mapping
  # use sf package for geometries (Shape files) and terra for raster files
  
  
  