use warnings;
use File::Basename;

# change path below to your correct path
$make_annot="/dcs04/lieber/shared/statsgen/LDSC/base/scripts/make_annot.py";
$referenceDir="/dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_EUR_Phase3_plink/";
$inputDir="/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC_old/snRNA_seq/input_files";
open(OUT, ">ldsc_anno_jobs_snRNA_seq.txt");
foreach $bedfile (<$inputDir/bedfiles/*.bed>){
	$out = basename($bedfile);
	$out =~ s/\.bed//;
	$outdir = $inputDir."/".$out;
	mkdir $outdir;
	foreach $i(1..22){
		$bimfile=$referenceDir."1000G.EUR.QC.".$i.".bim";
		print OUT "python $make_annot ";
		print OUT "--bed-file $bedfile ";
		print OUT "--bimfile $bimfile ";
		print OUT "--annot-file ${outdir}/chr.${i}.annot.gz\n";
	}
}
close(OUT);