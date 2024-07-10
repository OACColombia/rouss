require(tidyverse)
#' Convert time (Hour Minute Second) to decimal number.
#'
#' @description
#' This function helps to convert time information (Hour Minute Second) to decimal number. It was extracted from "Best Practices for Using eBird Data" by Matthew Strimas-Mackey, Wesley M. Hochachka, Viviana Ruiz-Gutierrez, Orin J. Robinson, Eliot T. Miller, Tom Auer, Steve Kelling, Daniel Fink, & Alison Johnston
#'
#' @param x A character vector of hour minute second
#'
#' @return A numeric vector of decimal time
#' @export
#'
#' @examples
#' x <- "7 6 5"
#' time_to_decimal(x)
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}
