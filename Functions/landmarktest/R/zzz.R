#' @importFrom RCurl getURL
#' 
#' 
.onAttach <- function(libname = find.package("landmarktest"), pkgname = "landmarktest") {
    sanitising_fun <- RCurl::getURL("https://raw.githubusercontent.com/TGuillerme/dispRity/master/R/sanitizing.R", ssl.verifypeer = FALSE)
    eval(parse(text = sanitising_fun))
}