library(GEOquery)
library(glue)
library(readr)

cfg <- list( 
    "track_dir" = "./test/data/track_1"
)
load_geo <- function(){
    fpath <- file.path(cfg$track_dir,"geo.txt")
    geo_files <- read_lines(fpath)
    geo_info <- map(geo_files, ~ getGEO(.x, destdir=cfg$track_dir))
    geo_file <- map(geo_files, ~ getGEOfile(.x, destdir=cfg$track_dir ))
    geo_parse <- map(
        geo_files, 
        ~ parseGEO(file.path(cfg$track_dir,glue("{.x}.soft")), destdir=cfg$track_dir )
    )
    geo_download <- map(
        geo_files, 
        ~ getGEOSuppFiles(.x, baseDir = cfg$track_dir)
    )
}