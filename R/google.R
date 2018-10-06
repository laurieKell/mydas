library(FLCore)

library(googledrive)
fls=drive_find()
drive_download("ray.RData",path="ray.RData")
load("ray.RData")


download.file("https://drive.google.com/open?id=1ASdwzY5TPwIpsbWcbYM00ZE_JjeGh7Dg", "turbot.RData")


