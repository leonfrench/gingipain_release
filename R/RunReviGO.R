library(here)
library(readr)
library(magrittr)
library(dplyr)
library(rvest)
library(httr)

#call with string with GO ID and uncorrected p-value. One diretion at a time.
#These are expected to be the first two columns of the input table.  
#some documentation at http://revigo.irb.hr/invokeRevigoAndFillFields.html
#GO:0009268 1e-14
#GO:0010447 1e-14
run_revigo <- function(goID_and_pvalue_table) {
  Sys.sleep(.5) #delay half a second to prevent overloading revigo
  if(ncol(goID_and_pvalue_table) != 2 ) stop('Input data frame is not two columns wide, should be go ID and p-value only (eg. GO:0009268 1e-14)')
  url <- "http://revigo.irb.hr/revigo.jsp"
  fd <- list(
    submit ="submitToRevigo",
    goList  = format_tsv(goID_and_pvalue_table, col_names = F),
    #goList  = "GO:0006614	2.139547041149198e-33
    #GO:0006613	6.532814810756933e-33",
    #goList  = "GO:0006613 1e-14
    #GO:0006614 1e-14",
    cutoff = "0.70",
    isPValue = "yes",
    whatIsBetter = "higher",
    goSizes = "0", #0 representing "whole UniProt (default)", alternatives are "Homo sapiens"= 9606 or another NCBI species code
    measure = "SIMREL" 
  )
  h1 <- handle('') #reset the handle
  resp <- POST(url, body=fd, encode="form",  handle=h1) #use verbose() for checking
  
  #for debugging
  #html(content(resp, as = "text"))
  #write(content(resp, as = "text"), "/Users/lfrench/Downloads/z.txt")
  #write(content(resp, as = "text"), "/Users/lfrench/Downloads/z.html")
  
  #clean up of the table
  extracted_data_frame <- html_table(content(resp))[[1]]
  extracted_data_frame <- extracted_data_frame[-1,]
  colnames(extracted_data_frame) <- extracted_data_frame[1,]
  extracted_data_frame <- extracted_data_frame[-1,]
  extracted_data_frame %<>% as_tibble()
  extracted_data_frame %<>% mutate(frequency = gsub(" %", "", frequency))
  extracted_data_frame %<>% type.convert()
  extracted_data_frame
}
