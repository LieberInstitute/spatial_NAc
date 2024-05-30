library(getopt)
library(sessioninfo)

# Import command-line parameters
spec <- matrix(
    c(
        c("final_step", "k"),
        c("f", "k"),
        rep("1", 2),
        c("character", "numeric"),
        c(
            "Using array coordinates after either 'imagej' or 'samui",
            "Number of clusters"
        )
    ),
    ncol = 5
)
opt <- getopt(spec)

print("Using the following parameters:")
print(opt)

session_info()

## This script was made using slurmjobs version 1.2.2
## available from http://research.libd.org/slurmjobs/
