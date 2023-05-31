#   First, had to go into 'nextflow.config' and set singularity to true instead
#   of docker. First run of the pipeline needs a huge amount of memory (e.g.
#   ~60G) to build the singularity image from the docker version

module load nextflow/22.10.7
module load singularity/3.2.1

WORK_DIR=$(git rev-parse --show-toplevel)/processed-data/03_webatlas_tests/work

nextflow run webatlas-pipeline/main.nf \
    -params-file ./config.yaml \
    -entry Full_pipeline \
    -w $WORK_DIR
