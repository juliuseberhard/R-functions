################################################################################
# Some Functions
# Author: Julius Eberhard
# Last Edit: 2016-11-21
################################################################################

# Package Check and Installation -----------------------------------------------
# reasonable for multiple packages only

CheckPack <- function(packages  # vector of packages
                      ) {

  for (i in 1:length(packages))
    if (!packages[i] %in% installed.packages()[, 1])
      install.packages(packages[[i]])

}

# Leap Year Check --------------------------------------------------------------
# (Author: Forester, quantitative-ecology.blogspot.de)

IsLeapyear <- function(year
                       ) {
  
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))

}
