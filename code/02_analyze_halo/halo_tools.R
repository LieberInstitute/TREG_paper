
#' Filter Halo Data By Region
#'
#' @param halo_data halo data from one Sample
#' @param bad_region_df Problem regions defined by X_min and X_max, Y_min and Y_max 
#'
#' @return halo_data annotaed with region_filter and problem column
#' @export
#'
#' @examples
#' halo_test_filter <- filter_halo_regions(halo_test, bad_regions)
#' halo_test_filter %>% count(region_filter, problem)
filter_halo_regions <- function(halo_data, bad_region_df){
  
  halo_data <- halo_data %>%
    mutate(problem = NA, region_filter = FALSE)
  
  if(nrow(bad_region_df) == 0) return(halo_data)
  
  for(row in 1:nrow(bad_region_df)){
    halo_data <- .filter_halo_region(halo_data, bad_region_df[row,])
  }
  
  halo_data <- halo_data %>% select(-filterX, -filterY)
  return(halo_data)
}


.filter_halo_region <- function(halo_data, bad_region){
  # message(bad_region$problem)
  h2 <- halo_data %>%
    mutate(filterX = halo_data$XMax <= bad_region$X_max & halo_data$XMax >= bad_region$X_min,
           filterY = halo_data$YMax <= bad_region$Y_max & halo_data$YMax >= bad_region$Y_min,
           region_filter = region_filter|(filterX & filterY),
           problem = ifelse((filterX & filterY), bad_region$problem, problem))
  return(h2)
}

