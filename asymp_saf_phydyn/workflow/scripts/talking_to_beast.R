##-----------------------------------------------------
## Util functions for beast
##------------------------------------------------------

library(readr)
library(lubridate)

to_date <- function(x, mrs) {
  round_date(date_decimal(decimal_date(ymd(mrs)) - x), unit = "day")
}

to_num <- function(x, mrs) {
  decimal_date(ymd(mrs)) - decimal_date(ymd(x))
}

read_trace <- function(traceFile, burninFrac){
  df <- read_table(traceFile, comment = "#") %>%
    mutate(file = traceFile)
  
  if (burninFrac > 0) {
    n <- dim(df)[1]
    df <- df[-(1:ceiling(burninFrac * n)), ]
  }
  
  return(df)
}

  