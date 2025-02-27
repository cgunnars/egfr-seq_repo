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

data/corr_results/20240502_pel-timecourse-6donor_corr.csv: src/plotCorr.R data/clean_dds/20240502_pel-timecourse-6donor_pathogen.Rds data/clean_dds/20240502_pel-timecourse-6donor_host.Rds 
	Rscript src/plotCorr.R -i $(word 2, $^) $(word 3, $^) -m Drug Day Donor Replicate -c gef_d1_vs_DMSO_d1 pel_d1_vs_DMSO_d1 sara_d1_vs_DMSO_d1 -d gef_vs_DMSO pel_vs_DMSO sara_vs_DMSO
corr: data/corr_results/20240502_pel-timecourse-6donor_corr.csv

fig/heatmap/20240502_pel-timecourse-6donor_pel-early.pdf: src/plotHeatmap.R data/DE_results/20240502_pel-timecourse-6donor_pathogen.Rds
	Rscript src/plotHeatmap.R -i $(word 2, $^) -c pel_d1_vs_DMSO_d1 -s pel_d1 DMSO_d1 phago_4h -o pel-early

fig/heatmap/20240502_pel-timecourse-6donor_gef-early.pdf: src/plotHeatmap.R data/DE_results/20240502_pel-timecourse-6donor_pathogen.Rds
	Rscript src/plotHeatmap.R -i $(word 2, $^) -c gef_d1_vs_DMSO_d1 -s gef_d1 DMSO_d1 phago_4h -o gef-early

fig/heatmap/20240502_pel-timecourse-6donor_sara-early.pdf: src/plotHeatmap.R data/DE_results/20240502_pel-timecourse-6donor_pathogen.Rds
	Rscript src/plotHeatmap.R -i $(word 2, $^) -c sara_d1_vs_DMSO_d1 -s sara_d1 DMSO_d1 phago_4h -o sara-early

fig/heatmap/20240502_pel-timecourse-6donor_all-early.pdf: src/plotHeatmap.R data/DE_results/20240502_pel-timecourse-6donor_pathogen.Rds
	Rscript src/plotHeatmap.R -i $(word 2, $^) -c pel_d1_vs_DMSO_d1 -s pel_d1 DMSO_d1 sara_d1 gef_d1 -o all-early

fig/heatmap/20240502_pel-timecourse-6donor_EGFR-early.pdf: src/plotHeatmap.R data/DE_results/20240502_pel-timecourse-6donor_pathogen.Rds
	Rscript src/plotHeatmap.R -i $(word 2, $^) -c pel_d1_vs_DMSO_d1 gef_d1_vs_DMSO_d1 -s pel_d1 DMSO_d1 gef_d1 -o EGFR-early

fig/heatmap/20240502_pel-timecourse-6donor_timecourse.pdf: src/plotHeatmap.R data/DE_results/20240502_pel-timecourse-6donor_pathogen.Rds
	Rscript src/plotHeatmap.R -i $(word 2, $^) -c DMSO_d1_vs_phago_4h DMSO_d3_vs_DMSO_d1 DMSO_d2_vs_DMSO_d1 -s DMSO_d3 DMSO_d2 DMSO_d1 phago_4h -o timecourse

fig/heatmap/20240502_pel-timecourse-6donor_timecourse-pel.pdf: src/plotHeatmap.R data/DE_results/20240502_pel-timecourse-6donor_pathogen.Rds
	Rscript src/plotHeatmap.R -i $(word 2, $^) -c DMSO_d1_vs_phago_4h pel_d1_vs_phago_4h pel_d3_vs_DMSO_d3 pel_d2_vs_DMSO_d2 pel_d1_vs_DMSO_d1 -s DMSO_d3 DMSO_d2 DMSO_d1 pel_d3 pel_d2 pel_d1 phago_4h -o timecourse-pel

fig/heatmap/20240502_pel-timecourse-6donor_host_all-early.pdf: src/plotHeatmap.R data/DE_results/20240502_pel-timecourse-6donor_host.Rds
	Rscript src/plotHeatmap.R -i $(word 2, $^) -c pel_vs_DMSO gef_vs_DMSO sara_vs_DMSO -s pel sara gef DMSO -o host_all-early

fig/heatmap/20240808_4h-abx_pathogen_lap.pdf: src/plotHeatmap.R data/DE_results/20240808_4h-abx_pathogen.Rds
	Rscript src/plotHeatmap.R -i $(word 2, $^) -c lap_5_vs_gef_5 lap_5_vs_pel_5 lap_25_vs_gef_25 lap_25_vs_pel_25 lap_25_vs_lap_5 -s lap_5 lap_25 gef_5 gef_25 pel_5 pel_25 -o lap
fig/heatmap/20240808_4h-abx_pathogen_var.pdf: src/plotHeatmap.R data/DE_results/20240808_4h-abx_pathogen.Rds
	Rscript src/plotHeatmap.R -i $(word 2, $^) -c var_5_vs_gef_5 var_5_vs_pel_5 var_25_vs_gef_25 var_25_vs_pel_25 var_25_vs_var_5 -s var_5 var_25 gef_5 gef_25 pel_5 pel_25 -o var
fig/heatmap/20240808_4h-abx_pathogen_dose.pdf: src/plotHeatmap.R data/DE_results/20240808_4h-abx_pathogen.Rds
	Rscript src/plotHeatmap.R -i $(word 2, $^) -c var_25_vs_var_5 gef_25_vs_gef_5 pel_25_vs_pel_5 -s var_5 var_25 gef_5 gef_25 pel_25 pel_5 lap_5 lap_25 -o dose


fig/heatmap/joint_axenic-6donor_pel.pdf: src/plotJointHeatmap.R data/DE_results/20240502_pel-timecourse-6donor_pathogen.Rds data/DE_results/20240816_24h-axenic_pathogen_pel_vs_DMSO_full.csv
	Rscript src/plotJointHeatmap.R -i $(word 2, $^) -c pel_d1_vs_DMSO_d1 -s pel_d1 DMSO_d1 -d $(word 3, $^) -o joint_axenic-6donor_pel

heatmaps: fig/heatmap/20240502_pel-timecourse-6donor_pel-early.pdf fig/heatmap/20240502_pel-timecourse-6donor_all-early.pdf fig/heatmap/20240502_pel-timecourse-6donor_host_all-early.pdf fig/heatmap/20240502_pel-timecourse-6donor_EGFR-early.pdf fig/heatmap/20240502_pel-timecourse-6donor_gef-early.pdf fig/heatmap/20240502_pel-timecourse-6donor_timecourse.pdf fig/heatmap/20240502_pel-timecourse-6donor_timecourse-pel.pdf
abx_heatmaps: fig/heatmap/20240808_4h-abx_pathogen_lap.pdf fig/heatmap/20240808_4h-abx_pathogen_var.pdf fig/heatmap/20240808_4h-abx_pathogen_dose.pdf

joint_heatmaps: fig/heatmap/joint_axenic-6donor_pel.pdf


fig/biplot/joint_axenic-6donor_pel.pdf: src/plotbiFC.R data/DE_results/20240502_pel-timecourse-6donor_pathogen_pel_d1_vs_DMSO_d1_full.csv data/DE_results/20240816_24h-axenic_pathogen_pel_vs_DMSO_full.csv 
	Rscript src/plotbiFC.R -i $(word 2, $^) -d $(word 3, $^) -o joint_axenic-6donor_pel
fig/biplot/joint_abx-6donor_pel.pdf: src/plotbiFC.R data/DE_results/20240502_pel-timecourse-6donor_pathogen_pel_d1_vs_DMSO_d1_full.csv data/DE_results/20240808_4h-abx_pathogen_pel_25_vs_pel_5_full.csv
	Rscript src/plotbiFC.R -i $(word 2, $^) -d $(word 3, $^) -o joint_abx-6donor_pel
fig/biplot/joint_abx-6donor_gef.pdf: src/plotbiFC.R data/DE_results/20240502_pel-timecourse-6donor_pathogen_gef_d1_vs_DMSO_d1_full.csv data/DE_results/20240808_4h-abx_pathogen_gef_25_vs_gef_5_full.csv
	Rscript src/plotbiFC.R -i $(word 2, $^) -d $(word 3, $^) -o joint_abx-6donor_gef
fig/biplot/joint_axenic-6donor_gef.pdf: src/plotbiFC.R data/DE_results/20240502_pel-timecourse-6donor_pathogen_gef_d1_vs_DMSO_d1_full.csv data/DE_results/20240816_24h-axenic_pathogen_pel_vs_DMSO_full.csv
	Rscript src/plotbiFC.R -i $(word 2, $^) -d $(word 3, $^) -o joint_axenic-6donor_gef
fig/biplot/joint_intracellular_gef-pel.pdf: src/plotbiFC.R data/DE_results/20240502_pel-timecourse-6donor_pathogen_pel_d1_vs_DMSO_d1_full.csv data/DE_results/20240502_pel-timecourse-6donor_pathogen_gef_d1_vs_DMSO_d1_full.csv
	Rscript src/plotbiFC.R -i $(word 2, $^) -d $(word 3, $^) -o joint_intracellular_gef-pel
biplots: fig/biplot/joint_axenic-6donor_pel.pdf fig/biplot/joint_abx-6donor_pel.pdf fig/biplot/joint_abx-6donor_gef.pdf fig/biplot/joint_axenic-6donor_gef.pdf fig/biplot/joint_intracellular_gef-pel.pdf
