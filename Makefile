data/raw_dds/20240808_4h-abx_pathogen.Rds: src/loadData.R data/raw_counts/20240808_4h-abx_pathogen.tsv
	Rscript src/loadData.R -i data/raw_counts/20240808_4h-abx_pathogen.tsv -l data/gene_info/H37Rv_gene-lengths.csv --pathogen -s Drug Dose Replicate -d Drug Dose

data/raw_dds/20240816_24h-axenic_pathogen.Rds: src/loadData.R data/raw_counts/20240816_24h-axenic_pathogen.tsv
	Rscript src/loadData.R -i data/raw_counts/20240816_24h-axenic_pathogen.tsv -l data/gene_info/H37Rv_gene-lengths.csv --pathogen -s Drug Replicate -d Drug

data/raw_dds/20240502_pel-timecourse-6donor_pathogen.Rds: src/loadData.R data/raw_counts/20240502_pel-timecourse-6donor_pathogen.tsv
	Rscript src/loadData.R -i data/raw_counts/20240502_pel-timecourse-6donor_pathogen.tsv -l data/gene_info/H37Rv_gene-lengths.csv --pathogen -s Drug Day Donor Replicate -d Drug Day

data/raw_dds/20240502_pel-timecourse-6donor_host.Rds: src/loadData.R data/raw_counts/20240502_pel-timecourse-6donor_host.Rds
	Rscript src/loadData.R -i data/raw_counts/20240502_pel-timecourse-6donor_host.Rds -s Drug Donor Replicate -d Drug -m Day d1

raw_dds: data/raw_dds/20240808_4h-abx_pathogen.Rds data/raw_dds/20240816_24h-axenic_pathogen.Rds data/raw_dds/20240502_pel-timecourse-6donor_pathogen.Rds data/raw_dds/20240502_pel-timecourse-6donor_host.Rds


data/clean_dds/%.Rds: src/runQC.R data/raw_dds/%.Rds
	Rscript src/runQC.R -i $(word 2, $^)
clean_dds: data/clean_dds/20240808_4h-abx_pathogen.Rds data/clean_dds/20240816_24h-axenic_pathogen.Rds data/clean_dds/20240502_pel-timecourse-6donor_pathogen.Rds data/clean_dds/20240502_pel-timecourse-6donor_host.Rds

data/DE_results/%.Rds: src/runDE.R data/clean_dds/%.Rds data/comparisons/comparisons_%.txt
	Rscript src/runDE.R -i $(word 2, $^) -c $(word 3, $^)
DE_dds: data/DE_results/20240808_4h-abx_pathogen.Rds data/DE_results/20240816_24h-axenic_pathogen.Rds data/DE_results/20240502_pel-timecourse-6donor_pathogen.Rds data/DE_results/20240502_pel-timecourse-6donor_host.Rds

data/corr_results/20240502_pel-timecourse-6donor_corr.csv: src/plotCorr.R data/DE_results/20240502_pel-timecourse-6donor_pathogen.Rds data/DE_results/20240502_pel-timecourse-6donor_host.Rds 
	Rscript src/plotCorr.R -i $(word 2, $^) $(word 3, $^) -m Drug Day Donor -c gef_d1_vs_DMSO_d1 pel_d1_vs_DMSO_d1 sara_d1_vs_DMSO_d1 -d gef_vs_DMSO pel_vs_DMSO sara_vs_DMSO
corr: data/corr_results/20240502_pel-timecourse-6donor_corr.csv
