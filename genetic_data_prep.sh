cd /Users/yu-hang/Library/CloudStorage/OneDrive-Personal/R_workstation/thyroid cancer/GWASTutorial/results
../tools/plink2 --bfile ../01_Dataset/1KG.EAS.ap.split.rare002.common015.missing --missing --out ../results/plink_results
head plink_results.smiss
head plink_results.vmiss
../tools/plink2 --bfile ../01_Dataset/1KG.EAS.ap.split.rare002.common015.missing --freq --out ../results/plink_results
../tools/plink2 --bfile 1KG.EAS.auto.snp.norm.nodup.split.rare002.common015.missing --hardy --out ../results/plink_results
../tools/plink2 --bfile 1KG.EAS.auto.snp.norm.nodup.split.rare002.common015.missing --maf 0.01 --geno 0.02 --mind 0.02 --hwe 1e-6 --indep-pairwise 50 5 0.2 --out ../results/plink_results
../tools/plink2 --bfile 1KG.EAS.auto.snp.norm.nodup.split.rare002.common015.missing --extract plink_results.prune.in --het --out ../results/plink_results
../tools/plink2 --bfile 1KG.EAS.auto.snp.norm.nodup.split.rare002.common015.missing --extract plink_results.prune.in --genome --out ../results/plink_results
../tools/plink2 --bfile 1KG.EAS.auto.snp.norm.nodup.split.rare002.common015.missing --chr 22 --r2-phased --out ../results/plink_results
../tools/plink2 --bfile 1KG.EAS.auto.snp.norm.nodup.split.rare002.common015.missing --maf 0.01 --geno 0.02 --mind 0.02 --hwe 1e-6 --indep-pairwise 50 5 0.2 --keep-allele-order --make-bed --out ../results/plink_results/sample_clean