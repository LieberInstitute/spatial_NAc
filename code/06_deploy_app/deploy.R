library("rsconnect")
library("here")
library("withr")

## Or you can go to your shinyapps.io account and copy this
## Here we do this to keep our information hidden.
# load(here("code", "06_deploy_app", ".deploy_info.Rdata"), verbose = TRUE)
# rsconnect::setAccountInfo(
#     name = deploy_info$name,
#     token = deploy_info$token,
#     secret = deploy_info$secret
# )

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Boost the maximum size to 5 GB
## https://support.posit.co/hc/en-us/articles/219449487--How-much-data-can-I-upload-to-shinyapps-io#:~:text=For%20the%20latter%20plans%2C%20note,option%20first%3A%20options(rsconnect.
# Error in `enforceBundleLimits()`:
# ! `appDir` (/Users/leocollado/Dropbox/Code/spatial_NAc/code/06_deploy_app) is too large to be deployed.
# ✖ The maximum size is 3145728000 bytes.
# ✖ This directory is at least 4688682982 bytes.
# ℹ Remove some files or adjust the rsconnect.max.bundle.size option.
# Run `rlang::last_trace()` to see where the error occurred.
options(rsconnect.max.bundle.size = 5 * 1024^3)

## Deploy the app, that is, upload it to shinyapps.io
rsconnect::deployApp(
    appDir = here("code", "06_deploy_app"),
    appFiles = c(
        "app.R",
        with_dir(here("code", "06_deploy_app"), dir("spe_shiny", full.names = TRUE)),
        with_dir(here("code", "06_deploy_app"), dir("www", full.names = TRUE))
    ),
    appName = "Spatial_NAc",
    account = "libd",
    server = "shinyapps.io"
)
