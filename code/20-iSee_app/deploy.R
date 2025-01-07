library("rsconnect")

source("token.R")

options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appFiles = c("app.R", "sce.Rdata", "initial.R"),
    appName = "snRNAseq_NAc",
    account = "libd",
    server = "shinyapps.io"
)
