#   Modelled after https://lcolladotor.github.io/bioc_team_ds/organizing-your-work.html#qsva-example

#!/bin/bash
#$ -cwd
#$ -l rnet,mem_free=2G,h_vmem=2G,h_fsize=400G
#$ -o transfer.log
#$ -e transfer.log
#$ -t 1-8
#$ -tc 8

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module list

## Locate file
SAMPLE=$(awk "NR==${SGE_TASK_ID}" samples.txt)

## Locate directories
ORIGINALDIR=/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/01_spaceranger_reorg/
NEWDIR=${ORIGINALDIR}/${SAMPLE}-1
CPDIR=${ORIGINALDIR}/${SAMPLE}/${SAMPLE}
#DESTDIR=/dcs05/lieber/marmaypag/dcl01Legacy_LIBD001/keri/



# ${SAMPLE} /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/01_spaceranger_reorg

#mkdir ${NEWDIR}

## Copy from dcl01 to dcs05
#rsync -rltgvh --chown=:lieber_marmaypag ${CPDIR}/ ${NEWDIR}/

mv ${CPDIR}/ ${NEWDIR}/

## Label as trash the files that were moved
#mv ${ORIGINALDIR} ${TRASHDIR}/

## Create link
#ln -s ${DESTDIR} ${ORIGINALDIR}


## Check settings
#nfs4_getfacl ${DESTDIR}

echo "**** Job ends ****"
date
