use warnings;
use Cwd;

# change path below to your correct path
$ldsc = "/dcs04/lieber/shared/statsgen/LDSC/base/scripts/ldsc.py";
$baseline2 = "/dcs04/lieber/shared/statsgen/LDSC/base/baseline2/baselineLD.";
$weights = "/dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.";
$freq = "/dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_Phase3_frq/1000G.EUR.QC.";

open(OUT, ">ldsc_h2_jobs_snRNA_seq.txt");

opendir($dh, "/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC_old/snRNA_seq/input_files/");
@outdirs = grep {!/\.csv$/i && !/\.tsv$/i && !/^\./ && !/bedfiles/} readdir($dh);
closedir $dh;

# change path below to your correct path
chdir "/dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain";
$cwd = cwd();
foreach $outdir (@outdirs){
	foreach $gwas(<*.gz>){
		$out = $gwas;
		$out =~ s/.gz$/.out/;
		print OUT "python $ldsc --h2 ${cwd}/$gwas ";
		print OUT "--w-ld-chr $weights ";
		print OUT "--ref-ld-chr /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC_old/snRNA_seq/input_files/${outdir}/chr.,$baseline2 ";
		print OUT "--overlap-annot ";
		print OUT "--frqfile-chr $freq ";
		print OUT "--out /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC_old/snRNA_seq/sLDSC_coefficients/${outdir}/$out ";
		print OUT "--print-coefficients\n";	
	}
}
close(OUT);
