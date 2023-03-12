## reloads packages and source file
reload_source <- function(){
  if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
  if (!require('readr')) install.packages('readr'); library('readr')
  if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
  if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
  if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
  if (!require('lubridate')) install.packages('lubridate'); library('lubridate')
  if (!require('DescTools')) install.packages('DescTools'); library('DescTools') # for Gtest
  if (!require('gtools')) install.packages('gtools'); library('gtools') 
  if (!require('survival')) install.packages('survival'); library('survival') 
  if (!require('kableExtra')) install.packages('kableExtra'); library('kableExtra') 
  if (!require('utils')) install.packages('utils'); library('utils') 
  if (!require('mgcv')) install.packages('mgcv'); library('mgcv') 
  if (!require('lubridate')) install.packages('lubridate'); library('lubridate') 
  if (!require('msm')) install.packages('msm'); library('msm') 
  if (!require('ipw')) install.packages('ipw'); library('ipw') 
  if (!require('patchwork')) install.packages('patchwork'); library('patchwork') 
  if (!require('sas7bdat')) install.packages('sas7bdat'); library('sas7bdat') 
  if (!require('haven')) install.packages('haven'); library('haven') 
  if (!require('janitor')) install.packages('haven'); library('janitor') 
  if (!require('rlist')) install.packages('rlist'); library('rlist') 
  source(here::here("source/utils.R"))
}

##' 
##' Utility function for making a crosstab that shows 
##' both proportion in each category and the 
##' parenthetical percentage.
##' 
##' @param data
##' @param row_column
##' @param col_column
##' @param grp_column
##' 
##' @return a dataframe with the appropriate data
##' 
summary_tbl <- function(data, row_column, col_column, grp_column){
  
  require(janitor)
  require(rlist)
  tbl <- data%>%tabyl(!!ensym(row_column),!!ensym(col_column),!!ensym(grp_column)) %>% 
    #group_by(!!ensym(row_column),
    #                     !!ensym(col_column))%>%
    #summarize(num=n()) %>% 
    #pivot_wider(names_from=!!ensym(col_column), 
    #            values_from =num)%>%
    replace(., is.na(.), 0) %>% 
    adorn_totals()%>%
    #adorn_totals("col")%>%
    adorn_percentages("col")%>%
    adorn_pct_formatting(digits=1)%>%
    adorn_ns(position="front") 
  
  tbl <- tbl[c("H1","H3","B","NA_")] %>% list.cbind()
  
  return(tbl)
}


summary_tbl2 <- function(data, row_column, col_column){
  
  require(janitor)
  require(rlist)
  tbl <- data%>%tabyl(!!ensym(row_column),!!ensym(col_column)) %>% 
    #group_by(!!ensym(row_column),
    #                     !!ensym(col_column))%>%
    #summarize(num=n()) %>% 
    #pivot_wider(names_from=!!ensym(col_column), 
    #            values_from =num)%>%
    replace(., is.na(.), 0) %>% 
    adorn_totals()%>%
    #adorn_totals("col")%>%
    adorn_percentages("col")%>%
    adorn_pct_formatting(digits=1)%>%
    adorn_ns(position="front") 
  return(tbl)
}

