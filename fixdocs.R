library(xml2)
library(tidyverse)
library(pkgdown)
build_site(document = T)
fs <- list.files(path = "docs/articles/simstandard_tutorial_files/figure-html",
                 pattern = "*.png")
fs <- paste0("docs/articles/simstandard_tutorial_files/figure-html/", fs)

file.rename(fs, str_replace_all(fs, "png", "svg"))


tutorial <-
  xml2::read_html("docs/articles/simstandard_tutorial.html")
img <- xml_find_all(tutorial, ".//img")
src <- xml_attr(img, "src")
newsrc <- str_replace_all(src, "png", "svg")
for (i in seq_len(newsrc)) {
  xml_set_attr(img[i], "src", newsrc[i])
}


write_html(tutorial, "docs/articles/simstandard_tutorial.html")
