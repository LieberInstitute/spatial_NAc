library("rsconnect")
library("here")
library("withr")

## Or you can go to your shinyapps.io account and copy this
## Here we do this to keep our information hidden.
load(here("code", "06_deploy_app", ".deploy_info.Rdata"), verbose = TRUE)
rsconnect::setAccountInfo(
    name = deploy_info$name,
    token = deploy_info$token,
    secret = deploy_info$secret
)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Deploy the app, that is, upload it to shinyapps.io
rsconnect::deployApp(
    appDir = here("code", "06_deploy_app"),
    appFiles = c(
        "app.R",
        "spe_shiny.rds",
        with_dir(here("code", "06_deploy_app"), dir("www", full.names = TRUE))
    ),
    appName = "Spatial_NAc",
    account = "libd",
    server = "shinyapps.io"
)
