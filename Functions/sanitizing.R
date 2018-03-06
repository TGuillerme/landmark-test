#####
## TEMPORARY: get sanitizing functions
#####
script <- RCurl::getURL("https://raw.githubusercontent.com/TGuillerme/dispRity/master/R/sanitizing.R", ssl.verifypeer = FALSE)
eval(parse(text = script))