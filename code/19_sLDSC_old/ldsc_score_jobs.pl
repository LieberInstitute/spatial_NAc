use warnings;

# change path below to your correct path
$ldsc="/dcs04/lieber/shared/statsgen/LDSC/base/scripts/ldsc.py";
$referenceDir="/dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_EUR_Phase3_plink/";
$hapmap3="/dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/hapmap3_snps/";

open(OUT, ">ldsc_score_jobs_10x_visium.txt");

opendir($dh, "/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC_old/10x_visium/input_files/");
@outdir = grep { !/^\./ && !/bedfiles/} readdir($dh);
closedir $dh;

foreach $out (@outdir){
  foreach $i(1..22){
	$bfile=$referenceDir."1000G.EUR.QC.".$i;
	$anno="/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC_old/10x_visium/input_files/$out/chr.".$i.".annot.gz";
	$hapmapSNP=$hapmap3."hm.".$i.".snp";
	print OUT "python $ldsc --l2 --thin-annot --ld-wind-cm 1 ";
	print OUT "--bfile $bfile ";
	print OUT "--anno $anno ";
	print OUT "--out /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC_old/10x_visium/input_files/$out/chr.$i ";
	print OUT "--print-snps $hapmapSNP\n";
  }
}

close(OUT);