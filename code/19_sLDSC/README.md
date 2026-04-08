The scripts in this directory can be used to performed sLDSC on variant annotations corresponding to specific gene sets. 

The R scripts below can be used to process the gene sets that we would like to derive variant annotations corresponding to. 
- 01-aggregate_*.R
    Where * refers to the gene set that is being used. For discrete cell clusters and SpDs, this corresponds to the aggregated log-normalized gene expression, for continuous spatial factors or NMF, we consider the correlation between the factor and spatial/ snRNA-seq log-normalized gene expression respectively. 
- 02-specificity_score_*.R
    Where * refers to the gene set that is being used. For discrete cell clusters and SpDs, for each gene normalize the expression across all reported cell types or SpDs to ensure we are selecting genes that are not just highly expressed in a particular cell type/SpD, but specific to it. Then pick the top x% of genes. For the NMFs/ MCPs, pick the top x% of genes ranked by correlation. 
- 03-bed_*.R
    Where * refers to the gene set that is being used. For each category, Cell type/ SpD/ NMF/ MCP create a bed file that includes the selected genes, and choose a window around the gene body. Here we picked +/- 100 Kb

Once we have the gene sets, we can run sLDSC in the following steps.
1. Create variant annotations: Example: python /dcs04/lieber/shared/statsgen/LDSC/base/scripts/make_annot.py --bed-file /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC/snRNA_seq/input_files/bedfiles/Astrocyte_A.bed --bimfile /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.1.bim --annot-file /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC/snRNA_seq/input_files/Astrocyte_A/chr.1.annot.gz

To run this step successfully, you will need
- Bed files corresponding to the selected genes and their respective windows within which you would like to annotate variants
- Bim files corresponding to 1000G_EUR_Phase3_plink, which will contain variant positions and IDs. 

2. Compute LD scores: Example: python /dcs04/lieber/shared/statsgen/LDSC/base/scripts/ldsc.py --l2 --thin-annot --ld-wind-cm 1 --bfile /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_EUR_Phase3_plink/1000G.EUR.QC.1 --anno /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC/snRNA_seq/input_files/Inh_D/chr.1.annot.gz --out /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC/snRNA_seq/input_files/Inh_D/chr.1 --print-snps /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/hapmap3_snps/hm.1.snp

To run this step succesfully, we need, 
- The annotations generated in the previous step
- List of SNPs that we will include, in this case hapmap3

3. Run sLDSC: Example: python /dcs04/lieber/shared/statsgen/LDSC/base/scripts/ldsc.py --h2 /dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain/adhd.gz --w-ld-chr /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --ref-ld-chr /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC/snRNA_seq/input_files/Inh_D/chr.,/dcs04/lieber/shared/statsgen/LDSC/base/baseline2/baselineLD. --overlap-annot --frqfile-chr /dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_Phase3_frq/1000G.EUR.QC. --out /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC/snRNA_seq/sLDSC_coefficients/Inh_D/adhd.out --print-coefficients

To run this step, we need, 
- The GWAS summary stats
- Regression weights downloaded from sLDSC
- Baseline functional annotations that we want to include
- Output location for coefficients

Post processing: 04-make_plots_*.R Read in the results from sLDSC including the partitioned hertiability estimate and the standardized effect size for the annotation. Then we examine the coefficient z_score, compute a corresponding p-value for the coefficient, and adjust for multiple hypothesis correction. 
