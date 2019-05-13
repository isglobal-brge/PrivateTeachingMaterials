#
# META-ANALYSIS
#

install.packages(c("meta", "rmeta", "metafor", "RCurl",
                   "altmeta", "RSImed"))

library(RCurl)

# P-curve
script <- getURL("https://raw.githubusercontent.com/nicebread/p-checker/master/p-curve.R", ssl.verifypeer = FALSE)
eval(parse(text = script))