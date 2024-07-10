require(tidyverse)
#' Convert time (Hour Minute Second) to decimal number.
#'
#' @description
#' This function helps to convert time information (Hour Minute Second) to decimal number.
#' It was extracted from "Best Practices for Using eBird Data" by Matthew Strimas-Mackey, Wesley M. Hochachka, Viviana Ruiz-Gutierrez, Orin J. Robinson, Eliot T. Miller, Tom Auer, Steve Kelling, Daniel Fink, & Alison Johnston.
#' This function requires the package `lubridate`, included in `tidyverse`.
#'
#' @param x A character vector of hour minute second
#'
#' @return A numeric vector of decimal time
#' @export
#'
#' @examples
#' x <- "5:30:00 AM"
#' time_to_decimal(x)
time_to_decimal <- function(x) {
  x <- lubridate::hms(x, quiet = TRUE)
  lubridate::hour(x) + lubridate::minute(x) / 60 + lubridate::second(x) / 3600
}
