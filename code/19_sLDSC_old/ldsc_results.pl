use warnings;
use Cwd;

$cwd = cwd();

# trait meta
open(IN, "/dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain/trait_names_keys");
while(<IN>){
	chomp;
	@tokens=split(/\t/,$_);
	$filename=$tokens[0];
	$filename =~ s/\.gz/\.out\.results/;
	#$type{$filename}=$tokens[1];
	$trait{$filename}=$tokens[2];
	print "$filename\t$trait{$filename}\n";
}
close(IN);

opendir($dh, "/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC_old/10x_visium/sLDSC_coefficients/");
@outdirs = grep { !/^\./} readdir($dh);
closedir $dh;

$count=0;

open(OUT, ">/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC_old/10x_visium/ldsc_results.txt");
print OUT "cell\ttrait\tProp._SNPs	Prop._h2	Prop._h2_std_error	Enrichment	Enrichment_std_error	Enrichment_p	Coefficient	Coefficient_std_error	Coefficient_z-score\n";
foreach $dir (@outdirs){
	$count++;
	print "$count\n";
	$group=$dir;
	chdir "/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/19_sLDSC_old/10x_visium/sLDSC_coefficients/${dir}";
	foreach $file(<*results>){
		if(not defined $trait{$file}){
			print "$file\n";
		}
		print OUT "$group\t$trait{$file}\t";
		open(IN,$file);
		<IN>;
		$line=<IN>;
		chomp $line;
		@tokens=split(/\t/,$line);
		shift @tokens;
		print OUT join("\t", @tokens), "\n";
		close(IN);
	}
	chdir "..";
}
close(OUT);