library(shinylive)
shinylive::export("app", "site")
httpuv::runStaticServer("site", port = 8008)
