# Define high-level variables -- name of experiment, should match name of raw .tsv
axenic_stem=20240816_24h-axenic_pathogen
abx_stem=20240808_4h-abx_pathogen
intra6h_stem=20240502_pel-timecourse-6donor_host
intra6p_stem=20240502_pel-timecourse-6donor_pathogen
intrap_stem=20240405_pel-timecourse_pathogen
# Named conditions used for testing
axenic_conds := gef pel 
intra_conds := gef_d1 pel_d1 sara_d1
intrahost_conds := gef pel sara
intraaxenic_conds := gef_d1 pel_d1
# name of design --> composite design term Drug_Day or Drug_Dose 
axenic_design := Drug
abx_design := Drug Dose
intra6p_design := Drug Day
intra6h_design := Drug
intrap_design := Drug Day

all_stems := $(axenic_stem) $(abx_stem) $(intra6h_stem) $(intra6p_stem) $(intrap_stem)

figs: fig1 fig2 fig3 fig4 fig5
sfigs: sfig_dose-upset sfig_regulators sfig_phago sfig_degmethod relative_heatmap_host 
tables: DE_tables combined_host_tables combined_pathogen_tables 

### PREPROCESSING into a standard format
data/raw_dds/$(abx_stem).Rds: src/loadData.R data/raw_counts/$(abx_stem).tsv
	Rscript $< -i $(word 2, $^) -l data/gene_info/H37Rv_gene-lengths.csv --pathogen -s $(abx_design) Replicate -d $(abx_design) 

data/raw_dds/$(axenic_stem).Rds: src/loadData.R data/raw_counts/$(axenic_stem).tsv
	Rscript $< -i $(word 2, $^) -l data/gene_info/H37Rv_gene-lengths.csv --pathogen -s $(axenic_design) Replicate -d $(axenic_design)

data/raw_dds/$(intra6p_stem).Rds: src/loadData.R data/raw_counts/$(intra6p_stem).tsv
	Rscript $< -i $(word 2, $^) -l data/gene_info/H37Rv_gene-lengths.csv --pathogen -s $(intra6p_design) Donor Replicate -d $(intra6p_design)

data/raw_dds/$(intra6h_stem).Rds: src/loadData.R data/raw_counts/$(intra6h_stem).Rds
	Rscript $< -i $(word 2, $^) -s $(intra6h_design) Donor Replicate -d $(intra6h_design) -m Day d1

data/raw_dds/$(intrap_stem).Rds: src/loadData.R data/raw_counts/$(intrap_stem).tsv
	Rscript $< -i $(word 2, $^) -l data/gene_info/H37Rv_gene-lengths.csv --pathogen -s $(intrap_design) Well -d $(intrap_design)	
raw := $(foreach n, $(all_stems), $(addprefix data/raw_dds/, $(addprefix $n, .Rds)))
raw_dds: $(raw)

### RUN QC
data/clean_dds/$(abx_stem).Rds: src/runQC.R data/raw_dds/$(abx_stem).Rds
	Rscript $< -i $(word 2, $^) -d gef_25_S15 pel_5_S13 pel_25_S12
data/clean_dds/%.Rds: src/runQC.R data/raw_dds/%.Rds
	Rscript $< -i $(word 2, $^)
### also generates QC plots and plots of expression similarity
fig/QC/%_post.pdf: data/clean_dds/%.Rds
	@if test -f $@; then :; else\ #remake using clean_dds as dependency
		rm -f $<; \
		make $<; \
	fi
fig/QC/%_pre.pdf: data/clean_dds/%.Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/QC/%_post-heatmap.pdf: data/clean_dds/%.Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/QC/%_pre-heatmap.pdf: data/clean_dds/%.Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
clean_dds := $(foreach n, $(all_stems), $(addprefix data/clean_dds/, $(addprefix $(n), .Rds)))

data/clean_dds/$(intrap_stem)_df.csv: src/makeGeneDf.R data/clean_dds/$(intrap_stem).Rds
	Rscript $< -i $(intrap_stem) -s $(intrap_design) Donor
data/clean_dds/$(intra6p_stem)_df.csv: src/makeGeneDf.R data/clean_dds/$(intra6p_stem).Rds
	Rscript $< -i $(intra6p_stem) -s $(intra6p_design) Donor
data/clean_dds/$(intra6h_stem)_df.csv: src/makeGeneDf.R data/clean_dds/$(intra6h_stem).Rds
	Rscript $< -i $(intra6h_stem) -s $(intra6h_design) Donor
data/clean_dds/$(axenic_stem)_df.csv: src/makeGeneDf.R data/clean_dds/$(axenic_stem).Rds
	Rscript $< -i $(axenic_stem) -s $(axenic_design) Replicate 
data/clean_dds/$(abx_stem)_df.csv: src/makeGeneDf.R data/clean_dds/$(abx_stem).Rds
	Rscript $< -i $(abx_stem) -s $(abx_design) Replicate 
clean_df := $(clean_dds:%.Rds=%_df.csv)

clean: $(clean_df) $(clean_dds) 


### RUN DE 
data/DE_results/%.Rds: src/runDE.R data/clean_dds/%.Rds data/comparisons/comparisons_%.txt
	Rscript $< -i $(word 2, $^) -c $(word 3, $^)
## also generates volcano plot and DE results table
fig/DE_results/$(abx_stem)_volcano_%.pdf: data/DE_results/$(abx_stem).Rds #must exist in comparisons.txt to rebuild correctly
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/DE_results/$(axenic_stem)_volcano_%.pdf: data/DE_results/$(axenic_stem).Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/DE_results/$(intra6p_stem)_volcano_%.pdf: data/DE_results/$(intra6p_stem).Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/DE_results/$(intrap_stem)_volcano_%.pdf: data/DE_results/$(intrap_stem).Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/DE_results/$(intra6h_stem)_volcano_%.pdf: data/DE_results/$(intra6h_stem).Rds	
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(abx_stem)_%_full.csv: data/DE_results/$(abx_stem).Rds 
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(axenic_stem)_%_full.csv: data/DE_results/$(axenic_stem).Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(intra6p_stem)_%_full.csv: data/DE_results/$(intra6p_stem).Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(intrap_stem)_%_full.csv: data/DE_results/$(intrap_stem).Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(intra6h_stem)_%_full.csv: data/DE_results/$(intra6h_stem).Rds	
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
## use list of comparisons to generate a list of results tables that should exist
axenic_de := $(foreach n, $(shell cat data/comparisons/comparisons_$(axenic_stem).txt), $(addprefix data/DE_results/$(axenic_stem)_, $(addsuffix _full.csv, $n)))
abx_de := $(foreach n, $(shell cat data/comparisons/comparisons_$(abx_stem).txt), $(addprefix data/DE_results/$(abx_stem)_, $(addsuffix _full.csv, $n)))
6p_de := $(foreach n, $(shell cat data/comparisons/comparisons_$(intra6p_stem).txt), $(addprefix data/DE_results/$(intra6p_stem)_, $(addsuffix _full.csv, $n)))
6h_de := $(foreach n, $(shell cat data/comparisons/comparisons_$(intra6h_stem).txt), $(addprefix data/DE_results/$(intra6h_stem)_, $(addsuffix _full.csv, $n)))
p_de := $(foreach n, $(shell cat data/comparisons/comparisons_$(intrap_stem).txt), $(addprefix data/DE_results/$(intrap_stem)_, $(addsuffix _full.csv, $n)))
DE_tables: $(axenic_de) $(abx_de) $(6p_de) $(6h_de) $(p_de)

## DE objects
DE := $(foreach n, $(all_stems), $(addprefix data/DE_results/, $(addprefix $n, .Rds)))
DE_dds: $(DE)


## FIGURE 1 CHEMICAL INFO
fig/chem_info/tanimoto.svg: src/getSMILES.py data/lux_data/compiled-clean-data.xlsx
	python $<
fig/chem_info/EGFR-spec-heatmap.svg: fig/chem_info/tanimoto.svg 
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/chem_info/EGFR-kd-ctrl.svg: src/getSMILES.py data/lux_data/compiled-clean-data.xlsx
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/chem_info/scaffold_grid.pdf: src/getSMILES.py 
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi

fig1: fig/chem_info/tanimoto.svg fig/chem_info/EGFR-spec-heatmap.svg fig/chem_info/EGFR-kd-ctrl.svg fig/chem_info/scaffold_grid.pdf

## FIGURE 2 ANTIBIOTIC EFFECTS
fig2: fig/abx-dose/pca_all_1_2.pdf fig/abx-dose/var-heatmap.pdf fig/abx-dose/lap-heatmap.pdf  
fig3: fig/abx-dose/pca_nolap_1_2.pdf fig/abx-dose/pca_nolap_1_3.pdf fig/abx-dose/pca_nolap_2_3.pdf fig/abx-dose/non-lap_loadings.pdf  fig/abx-dose/pel-dose-heatmap.pdf fig/abx-dose/gef-dose-heatmap.pdf  
sfig_dose-upset: fig/abx-dose/lap-var_upset.pdf fig/abx-dose/lap_upset.pdf fig/abx-dose/var_upset.pdf

## abx-dose.R generates figures 2 and 3
## generates upset plot, curated heatmaps, pca with and without lapatinib
fig/abx-dose/lap-var_upset.pdf: src/abx-dose.R $(axenic_de) 
	Rscript $<
fig/abx-dose/var_upset.pdf: fig/abx-dose/lap-var_upset.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/abx-dose/lap_upset.pdf: fig/abx-dose/lap-var_upset.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/abx-dose/dose_upset.pdf: fig/abx-dose/lap-var_upset.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/abx-dose/var-heatmap.pdf: fig/abx-dose/lap-var_upset.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/abx-dose/lap-heatmap.pdf: fig/abx-dose/lap-var_upset.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/abx-dose/pca_nolap_1_2.pdf: fig/abx-dose/lap-var_upset.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/abx-dose/pca_nolap_1_3.pdf: fig/abx-dose/lap-var_upset.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/abx-dose/pca_nolap_2_3.pdf: fig/abx-dose/lap-var_upset.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/abx-dose/pca_all_1_2.pdf: fig/abx-dose/lap-var_upset.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi

## FIGURE 4 intraCELLULAR 
fig4: $(relative_heatmaps_pathogen) iMod_enrich_unique fig/growth/RNAseq_lux.pdf fig/combined_bar/$(intra6p_stem)_fill.pdf fig/combined_bar/$(intra6p_stem)_n.pdf fig/combined_bar/$(intra6p_stem)_vennlikely.pdf

fig/growth/RNAseq_lux.pdf: src/plotDay1.R #data/lux_data/pel-clean_data.xlsx
	Rscript $< 

## Make comparisons of "DE / not DE" for each drug relative to the other two
combined_p_conds := pel_d1_gef_d1_sara_d1 gef_d1_pel_d1_sara_d1 sara_d1_pel_d1_gef_d1
combined_h_conds := pel_gef_sara gef_pel_sara sara_pel_gef
relative_heatmap_host: $(addprefix fig/relative_heatmap/relative_heatmap_, $(addsuffix .pdf, $(addprefix $(intra6h_stem)_,  $(combined_h_conds))))
relative_heatmap_pathogen: $(addprefix fig/relative_heatmap/relative_heatmap_, $(addsuffix .pdf, $(addprefix $(intra6p_stem)_,  $(combined_p_conds)))) 

fig/relative_heatmap/relative_heatmap_$(intra6p_stem)_pel_d1_gef_d1_sara_d1.pdf: src/relative_heatmap.R data/DE_results/$(intra6p_stem).Rds
	Rscript $< -i $(intra6p_stem) -r DMSO_d1 -c pel_d1 gef_d1 sara_d1 -g Drug_Day
fig/relative_heatmap/relative_heatmap_$(intra6p_stem)_gef_d1_pel_d1_sara_d1.pdf: src/relative_heatmap.R data/DE_results/$(intra6p_stem).Rds
	Rscript $< -i $(intra6p_stem) -r DMSO_d1 -c gef_d1 pel_d1 sara_d1 -g Drug_Day
fig/relative_heatmap/relative_heatmap_$(intra6p_stem)_sara_d1_pel_d1_gef_d1.pdf: src/relative_heatmap.R data/DE_results/$(intra6p_stem).Rds
	Rscript $< -i $(intra6p_stem) -r DMSO_d1 -c sara_d1 pel_d1 gef_d1 -g Drug_Day
data/DE_results/combined/$(intra6p_stem)_%.csv: fig/relative_heatmap/relative_heatmap_$(intra6p_stem)_%.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
relative_heatmaps: $(relative_heatmap_host) $(relative_heatmap_pathogen)
combined_pathogen_tables: $(addprefix data/DE_results/combined/$(intra6p_stem)_, $(addsuffix .csv, $(combined_p_conds)))
combined_host_tables: $(addprefix data/DE_results/combined/$(intra6h_stem)_, $(addsuffix .csv, $(combined_h_conds))) 


fig/combined_bar/$(intra6p_stem)_fill.pdf: src/combined_bar_venn.R combined_pathogen_tables
	Rscript $< -e $(intra6p_stem) -c pel_d1 gef_d1 sara_d1 -r DMSO_d1 -g Drug_Day
fig/combined_bar/$(intra6p_stem)_n.pdf: fig/combined_bar/$(intra6p_stem)_fill.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/combined_bar/$(intra6p_stem)_vennlikely.pdf: fig/combined_bar/$(intra6p_stem)_fill.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(intra6p_stem)_%_d1_unique.txt: fig/combined_bar/$(intra6p_stem).pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(intra6p_stem)_likely_shared.txt: fig/combined_bar/$(intra6p_stem).pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
combined_bar: fig/combined_bar/$(intra6p_stem)_fill.pdf fig/combined_bar/$(intra6h_stem)_fill.pdf

fig/relative_heatmap/all_drugs_iModulon.pdf: src/geneListToGSEA.R data/DE_results/$(intra6p_stem)_likely_shared.txt 
	Rscript $< -c all -m drugs -n 320 -e $(intra6p_stem) 
data/enrich/all_drugs_%.csv: fig/relative_heatmap/all_drugs_iModulon.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
iMod_enrich_unique: fig/relative_heatmap/all_drugs_iModulon.pdf

### SUPPLEMENTARY FIGURE POST PHAGOCYTOSIS 
phago_plots := $(addprefix fig/time-dependent/phago_, $(addsuffix _hm.pdf, $(addprefix $(intra6p_stem)_, $(intra_conds))))
sfig_phago: $(phago_plots) fig/time-dependent/phago_$(intrap_stem)_pel_d1_hm.pdf

fig/time-dependent/phago_$(intra6p_stem)_%_hm.pdf: src/plotPhago.R data/DE_results/$(intra6p_stem).Rds	
	Rscript $< -i $(intra6p_stem) -c $*
fig/time-dependent/phago_$(intrap_stem)_%_hm.pdf: src/plotPhago.R data/DE_results/$(intrap_stem).Rds
	Rscript $< -i $(intrap_stem) -c pel_d1 

### SUPPLEMENTARY FIGURE SIGMA AND TCS
fig/regulators/sigma_$(intra6p_stem)_%_d1.pdf: src/plotRegulators.R data/DE_results/$(intra6p_stem)_%_d1_vs_DMSO_d1_full.csv data/clean_dds/$(intra6p_stem)_df.csv
	Rscript $< -i $(intra6p_stem) -c $* -r DMSO
fig/regulators/sigma_%_pel_d1.pdf: src/plotRegulators.R data/DE_results/%_pel_d1_vs_DMSO_d1_full.csv data/clean_dds/%_df.csv
	Rscript $< -i $* -c pel -r DMSO
regulators_plot := $(addprefix fig/regulators/sigma_$(intra6p_stem)_, $(addsuffix .pdf, $(intra_conds))) fig/regulators/sigma_$(intrap_stem)_pel_d1.pdf
sfig_regulators: $(regulators_plot)

### SUPPLEMENTARY FIGURE GSEA 
axenic_gsea := $(addprefix fig/gsea/gsea_, $(addsuffix _$(axenic_stem)_axenic.pdf, $(axenic_conds)))
intra_gsea := $(addprefix fig/gsea/gsea_, $(addsuffix _$(intra6p_stem)_intra.pdf, $(intra_conds)))
fig/gsea/gsea_%_$(axenic_stem)_axenic.pdf: src/plotGSEA.R data/DE_results/$(axenic_stem).Rds
	Rscript $< -i data/DE_results/$(axenic_stem)_$*_vs_DMSO_full.csv -g marR sigE kstR -o $*_$(axenic_stem)_axenic
fig/gsea/gsea_%_$(intra6p_stem)_intra.pdf: src/plotGSEA.R data/DE_results/$(intra6p_stem).Rds
	Rscript $< -i data/DE_results/$(intra6p_stem)_$*_vs_DMSO_d1_full.csv -g marR sigE kstR -o $*_$(intra6p_stem)_intra
fig/gsea/gsea_pel_d1_$(intrap_stem)_intra.pdf: src/plotGSEA.R data/DE_results/$(intrap_stem).Rds
	Rscript $< -i data/DE_results/$(intrap_stem)_pel_d1_vs_DMSO_d1_full.csv -g marR sigE kstR -o pel_d1_$(intrap_stem)_intra
data/enrich/gsea_%.csv: fig/gsea/gsea_%.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
gsea: $(axenic_gsea) $(intra_gsea) fig/gsea/gsea_pel_d1_$(intrap_stem)_intra.pdf

### SUPPLEMENTARY FIGURE HOST INFORMATION
fig/relative_heatmap/relative_heatmap_$(intra6h_stem)_pel_gef_sara.pdf: src/relative_heatmap.R data/DE_results/$(intra6h_stem).Rds 
	Rscript $< -i $(intra6h_stem) -r DMSO -c pel gef sara -g Drug
data/DE_results/combined/$(intra6h_stem)_%.csv: fig/relative_heatmap/relative_heatmap_$(intra6h_stem)_%.pdf 
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/relative_heatmap/relative_heatmap_$(intra6h_stem)_gef_pel_sara.pdf: src/relative_heatmap.R data/DE_results/$(intra6h_stem).Rds 
	Rscript $< -i $(intra6h_stem) -r DMSO -c gef pel sara -g Drug
fig/relative_heatmap/relative_heatmap_$(intra6h_stem)_sara_pel_gef.pdf: src/relative_heatmap.R data/DE_results/$(intra6h_stem).Rds 
	Rscript $< -i $(intra6h_stem) -r DMSO -c sara pel gef -g Drug

fig/combined_bar/$(intra6h_stem)_fill.pdf: src/combined_bar_venn.R combined_host_tables
	Rscript $< -e $(intra6h_stem) -c pel gef sara -r DMSO -g Drug

fig/combined_bar/$(intra6h_stem)_n.pdf: fig/combined_bar/$(intra6h_stem)_fill.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/combined_bar/$(intra6h_stem)_vennlikely.pdf: fig/combined_bar/$(intra6h_stem)_fill.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(intra6h_stem)_%_unique.txt: fig/combined_bar/$(intra6h_stem).pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(intra6h_stem)_likely_shared.txt: fig/combined_bar/$(intra6h_stem).pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi

## FIGURE 5 AXENIC EFFECTS
fig5: fig/DE_results/$(axenic_stem)_volcano_pel_vs_DMSO.pdf fig/DE_results/$(axenic_stem)_volcano_gef_vs_DMSO.pdf $(axenic_joint_heatmaps) iMod_enrich_axenic	

# AXENIC EFFECT CATEGORIZATION PLOT + HEATMAP
fig/axenic_heatmap/axenic_heatmap_%_d1_$(axenic_stem).pdf: src/axenic_heatmap.R data/DE_results/$(axenic_stem).Rds data/DE_results/$(intra6p_stem).Rds
	Rscript $< -i $(intra6p_stem) -a $(axenic_stem) -c $*_d1 -d $* -v DMSO_d1 -w DMSO -g Drug_Day -j Drug
fig/axenic_heatmap/axenic_heatmap_%_d1_$(abx_stem).pdf: src/axenic_heatmap.R data/DE_results/$(abx_stem).Rds data/DE_results/$(intra6p_stem).Rds
	Rscript $< -i $(intra6p_stem) -a $(abx_stem) -c $*_d1 -d $*_25 -v DMSO_d1 -w $*_5 -g Drug_Day -j Drug_Dose
fig/axenic_heatmap/axenic_heatmap_pel_d1_$(intrap_stem).pdf: src/axenic_heatmap.R data/DE_results/$(intrap_stem).Rds data/DE_results/$(intra6p_stem).Rds
	Rscript $< -i $(intra6p_stem) -a $(intrap_stem) -c pel_d1 -d pel_d1 -v DMSO_d1 -w DMSO_d1 -g Drug_Day -j Drug_Day
data/DE_results/combined/combined_intraaxenic_%.csv: fig/axenic_heatmap/axenic_heatmap_%.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/axenic_heatmap/venn_biplot_%.pdf: fig/axenic_heatmap/axenic_heatmap_%.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi

axenic_joint_heatmaps := $(addprefix fig/axenic_heatmap/axenic_heatmap_, $(addsuffix _d1_$(axenic_stem).pdf, $(axenic_conds)))
abx_joint_heatmaps := $(addprefix fig/axenic_heatmap/axenic_heatmap_, $(addsuffix _d1_$(abx_stem).pdf, $(axenic_conds)))
joint_heatmaps: $(axenic_joint_heatmaps) $(abx_joint_heatmaps) fig/axenic_heatmap/axenic_heatmap_pel_d1_$(intrap_stem).pdf

combined_intraaxenic := $(addprefix data/DE_results/combined/combined_intraaxenic_, $(addsuffix _$(axenic_stem).csv, $(intraaxenic_conds)))
combined_intraabx := $(addprefix data/DE_results/combined/combined_intraaxenic_, $(addsuffix _$(abx_stem).csv, $(intraaxenic_conds)))
combined_intraaxenic_tables: $(combined_intraaxenic) $(combined_intraabx)

# TEST ENRICHMENT WHILE EXCLUDING AXENIC EFFECTS
fig/relative_heatmap/%_intraaxenic_allcomps_iModulon.pdf: src/geneListToGSEA.R data/DE_results/$(intra6p_stem)_likely_shared.txt data/DE_results/combined/combined_intraaxenic_%_$(axenic_stem).csv
	Rscript $< -c $* -m intraaxenic -e $(intra6p_stem) $(axenic_stem)
data/enrich/pel_d1_intraaxenic_%.csv: fig/relative_heatmap/pel_d1_intraaxenic_allcomps_iModulon.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/enrich/gef_d1_intraaxenic_%.csv: fig/relative_heatmap/gef_d1_intraaxenic_allcomps_iModulon.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi

fig/relative_heatmap/%_single_iModulon.pdf: src/geneListToGSEA.R
	Rscript $< -c $* -m single -e $(intrap_stem)
data/enrich/%_single_pel_d1.csv: fig/relative_heatmap/%_single_iModulon.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
iMod_enrich_axenic: fig/relative_heatmap/pel_d1_intraaxenic_allcomps_iModulon.pdf fig/relative_heatmap/gef_d1_intraaxenic_allcomps_iModulon.pdf
	

## SUPPLEMENTARY FIGURE -- WALD TEST SCHEME
sfig_degmethod_intra: fig/deg-method/sfig6a_pel_intra.pdf fig/deg-method/sfig6b_pel_intra.pdf fig/deg-method/sfig6c_pel_intra.pdf
sfig_degmethod_axenic: fig/deg-method/sfig6a_pel_axenic.pdf
sfig_degmethod: sfig_degmethod_intra sfig_degmethod_axenic

fig/deg-method/sfig6a_%.pdf: src/testDEMethods.R data/DE_results/$(intra6p_stem).Rds
	Rscript $< -i data/DE_results/$(intra6p_stem).Rds -c Drug_Day_pel_d1_vs_DMSO_d1 -o $*
fig/deg-method/sfig6b_%.pdf: fig/deg-method/sfig6a_%.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/deg-method/sfig6c_%.pdf: fig/deg-method/sfig6a_%.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
