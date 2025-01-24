data/raw_dds/20240808_4h-abx_pathogen.Rds: src/loadData.R data/raw_counts/20240808_4h-abx_pathogen.tsv
	Rscript src/loadData.R -i data/raw_counts/20240808_4h-abx_pathogen.tsv -l data/gene_info/H37Rv_gene-lengths.csv --pathogen -s Drug Dose Replicate -d Drug Dose

data/raw_dds/20240816_24h-axenic_pathogen.Rds: src/loadData.R data/raw_counts/20240816_24h-axenic_pathogen.tsv
	Rscript src/loadData.R -i data/raw_counts/20240816_24h-axenic_pathogen.tsv -l data/gene_info/H37Rv_gene-lengths.csv --pathogen -s Drug Replicate -d Drug

data/raw_dds/20240502_pel-timecourse-6donor_pathogen.Rds: src/loadData.R data/raw_counts/20240502_pel-timecourse-6donor_pathogen.tsv
	Rscript src/loadData.R -i data/raw_counts/20240502_pel-timecourse-6donor_pathogen.tsv -l data/gene_info/H37Rv_gene-lengths.csv --pathogen -s Drug Day Donor Replicate -d Drug Day

data/raw_dds/20240502_pel-timecourse-6donor_host.Rds: src/loadData.R data/raw_counts/20240502_pel-timecourse-6donor_host.Rds
	Rscript src/loadData.R -i data/raw_counts/20240502_pel-timecourse-6donor_host.Rds -s Drug Donor Replicate -d Drug

raw_dds: data/raw_dds/20240808_4h-abx_pathogen.Rds data/raw_dds/20240816_24h-axenic_pathogen.Rds data/raw_dds/20240502_pel-timecourse-6donor_pathogen.Rds data/raw_dds/20240502_pel-timecourse-6donor_host.Rds


data/QC_dds/%.Rds: src/runQC.R data/raw_dds/20240808_4h-abx_pathogen.Rds
	Rscript src/runQC.R -i data/raw_dds/20240808_4h-abx_pathogen.Rds

data/clean_dds/%.Rds: src/runDE.R data/raw_dds/*.Rds data/*_comparisons.txt
	Rscript src/runDE.R -i data/raw_dds/%.Rds -c data/%_comparisons.txt
