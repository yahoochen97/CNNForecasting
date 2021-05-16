# define horizons in the loyo process
horizons = c('0',
               '7',
               '14',
               '21',
               '28',
               '42',
               '56')

# number of candidate points for loyo
search_sizes = c(100,20)

# all years in loyo process
cv_years = c(1992,1994,1996,1998,2000,2002,2004,2006,2008,2010,2012,2014,2016)

# model type
TYPES = c('GP','LM')

for (i in 2:length(TYPES)) {
  MAX_AVE_NLZ = c()
  HORIZONS = c()
  OPT_IDX = c()
  TYPE = TYPES[i]
  search_size = search_sizes[i]
  for (a in 1:length(horizons)) {
    averaged_LL = rep(0, search_size)
    for (b in 1:search_size){
      year_nlZs = c()
      for (cv_year in cv_years) {
        input_file = paste('nlZs/', TYPE, '_' , cv_year, 'day', horizons[a], '_', b,'.csv',sep='')
        nlZ <- read.csv(input_file)
        nlZ = nlZ$x
        year_nlZs = c(year_nlZs, nlZ)
      }
      averaged_LL[b] = -mean(year_nlZs)
    }
    # obtain the maximal index of the averaged ll array
    MAX_AVE_NLZ = c(MAX_AVE_NLZ, max(averaged_LL))
    HORIZONS = c(HORIZONS, strtoi(horizons[a]))
    OPT_IDX = c(OPT_IDX, which(averaged_LL==max(averaged_LL)))
  }
  # write results to csv
  output_file = paste('results/', TYPE, '_opthyp.csv',sep='')
  result <- data.frame(HORIZONS, MAX_AVE_NLZ, OPT_IDX)
  names(result) <- tolower(names(result))
  write.csv(result,output_file)
}
