AXENIC_STEM=20240816_24h-axenic_pathogen
ABX_STEM=20240808_4h-abx_pathogen
INTRA6H_STEM=20240502_pel-timecourse-6donor_host
INTRA6P_STEM=20240502_pel-timecourse-6donor_pathogen
INTRAP_STEM=20240405_pel-timecourse_pathogen

AXENIC_CONDS := gef pel
INTRA_CONDS := gef_d1 pel_d1 sara_d1
INTRAHOST_CONDS := gef pel sara
INTRAAXENIC_CONDS := gef_d1 pel_d1

ALL_STEMS := $(AXENIC_STEM) $(ABX_STEM) $(INTRA6H_STEM) $(INTRA6P_STEM) $(INTRAP_STEM)

figs: fig1 fig2 fig3 fig4 fig5
sfigs: sfig_dose-upset sfig_regulators sfig_phago sfig_degmethod relative_heatmap_host 
tables: DE_tables combined_host_tables combined_pathogen_tables 

### PREPROCESSING into a standard format
data/raw_dds/$(ABX_STEM).Rds: src/loadData.R data/raw_counts/$(ABX_STEM).tsv
	Rscript $< -i $(word 2, $^) -l data/gene_info/H37Rv_gene-lengths.csv --pathogen -s Drug Dose Replicate -d Drug Dose

data/raw_dds/$(AXENIC_STEM).Rds: src/loadData.R data/raw_counts/$(AXENIC_STEM).tsv
	Rscript $< -i $(word 2, $^) -l data/gene_info/H37Rv_gene-lengths.csv --pathogen -s Drug Replicate -d Drug

data/raw_dds/$(INTRA6P_STEM).Rds: src/loadData.R data/raw_counts/$(INTRA6P_STEM).tsv
	Rscript $< -i $(word 2, $^) -l data/gene_info/H37Rv_gene-lengths.csv --pathogen -s Drug Day Donor Replicate -d Drug Day

data/raw_dds/$(INTRA6H_STEM).Rds: src/loadData.R data/raw_counts/$(INTRA6H_STEM).Rds
	Rscript $< -i $(word 2, $^) -s Drug Donor Replicate -d Drug -m Day d1

data/raw_dds/$(INTRAP_STEM).Rds: src/loadData.R data/raw_counts/$(INTRAP_STEM).tsv
	Rscript $< -i $(word 2, $^) -l data/gene_info/H37Rv_gene-lengths.csv --pathogen -s Drug Day Well -d Drug Day	
raw := $(foreach n, $(ALL_STEMS), $(addprefix data/raw_dds/, $(addprefix $n, .Rds)))
raw_dds: $(raw)

### RUN QC
data/clean_dds/$(ABX_STEM).Rds: src/runQC.R data/raw_dds/$(ABX_STEM).Rds
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
clean := $(foreach n, $(ALL_STEMS), $(addprefix data/clean_dds/, $(addprefix $(n), .Rds)))
clean_dds: $(clean) 

### RUN DE 
data/DE_results/%.Rds: src/runDE.R data/clean_dds/%.Rds data/comparisons/comparisons_%.txt
	Rscript $< -i $(word 2, $^) -c $(word 3, $^)
## also generates volcano plot and DE results table
fig/DE_results/$(ABX_STEM)_volcano_%.pdf: data/DE_results/$(ABX_STEM).Rds #must exist in comparisons.txt to rebuild correctly
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/DE_results/$(AXENIC_STEM)_volcano_%.pdf: data/DE_results/$(AXENIC_STEM).Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/DE_results/$(INTRA6P_STEM)_volcano_%.pdf: data/DE_results/$(INTRA6P_STEM).Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/DE_results/$(INTRAP_STEM)_volcano_%.pdf: data/DE_results/$(INTRAP_STEM).Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/DE_results/$(INTRA6H_STEM)_volcano_%.pdf: data/DE_results/$(INTRA6H_STEM).Rds	
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(ABX_STEM)_%_full.csv: data/DE_results/$(ABX_STEM).Rds 
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(AXENIC_STEM)_%_full.csv: data/DE_results/$(AXENIC_STEM).Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(INTRA6P_STEM)_%_full.csv: data/DE_results/$(INTRA6P_STEM).Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(INTRAP_STEM)_%_full.csv: data/DE_results/$(INTRAP_STEM).Rds
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(INTRA6H_STEM)_%_full.csv: data/DE_results/$(INTRA6H_STEM).Rds	
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
## use list of comparisons to generate a list of results tables that should exist
axenic_de := $(foreach n, $(shell cat data/comparisons/comparisons_$(AXENIC_STEM).txt), $(addprefix data/DE_results/$(AXENIC_STEM)_, $(addsuffix _full.csv, $n)))
abx_de := $(foreach n, $(shell cat data/comparisons/comparisons_$(ABX_STEM).txt), $(addprefix data/DE_results/$(ABX_STEM)_, $(addsuffix _full.csv, $n)))
6p_de := $(foreach n, $(shell cat data/comparisons/comparisons_$(INTRA6P_STEM).txt), $(addprefix data/DE_results/$(INTRA6P_STEM)_, $(addsuffix _full.csv, $n)))
6h_de := $(foreach n, $(shell cat data/comparisons/comparisons_$(INTRA6H_STEM).txt), $(addprefix data/DE_results/$(INTRA6H_STEM)_, $(addsuffix _full.csv, $n)))
p_de := $(foreach n, $(shell cat data/comparisons/comparisons_$(INTRAP_STEM).txt), $(addprefix data/DE_results/$(INTRAP_STEM)_, $(addsuffix _full.csv, $n)))
DE_tables: $(axenic_de) $(abx_de) $(6p_de) $(6h_de) $(p_de)

## DE objects
DE := $(foreach n, $(ALL_STEMS), $(addprefix data/DE_results/, $(addprefix $n, .Rds)))
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

## FIGURE 4 INTRACELLULAR 
fig4: $(relative_heatmaps_pathogen) iMod_enrich_unique fig/growth/RNAseq_lux.pdf fig/combined_bar/$(INTRA6P_STEM)_fill.pdf fig/combined_bar/$(INTRA6P_STEM)_n.pdf fig/combined_bar/$(INTRA6P_STEM)_vennlikely.pdf

fig/growth/RNAseq_lux.pdf: src/plotDay1.R #data/lux_data/pel-clean_data.xlsx
	Rscript $< 

## Make comparisons of "DE / not DE" for each drug relative to the other two
combined_p_conds := pel_d1_gef_d1_sara_d1 gef_d1_pel_d1_sara_d1 sara_d1_pel_d1_gef_d1
combined_h_conds := pel_gef_sara gef_pel_sara sara_pel_gef
relative_heatmap_host: $(addprefix fig/relative_heatmap/relative_heatmap_, $(addsuffix .pdf, $(addprefix $(INTRA6H_STEM)_,  $(combined_h_conds))))
relative_heatmap_pathogen: $(addprefix fig/relative_heatmap/relative_heatmap_, $(addsuffix .pdf, $(addprefix $(INTRA6P_STEM)_,  $(combined_p_conds)))) 

fig/relative_heatmap/relative_heatmap_$(INTRA6P_STEM)_pel_d1_gef_d1_sara_d1.pdf: src/relative_heatmap.R data/DE_results/$(INTRA6P_STEM).Rds
	Rscript $< -i $(INTRA6P_STEM) -r DMSO_d1 -c pel_d1 gef_d1 sara_d1 -g Drug_Day
fig/relative_heatmap/relative_heatmap_$(INTRA6P_STEM)_gef_d1_pel_d1_sara_d1.pdf: src/relative_heatmap.R data/DE_results/$(INTRA6P_STEM).Rds
	Rscript $< -i $(INTRA6P_STEM) -r DMSO_d1 -c gef_d1 pel_d1 sara_d1 -g Drug_Day
fig/relative_heatmap/relative_heatmap_$(INTRA6P_STEM)_sara_d1_pel_d1_gef_d1.pdf: src/relative_heatmap.R data/DE_results/$(INTRA6P_STEM).Rds
	Rscript $< -i $(INTRA6P_STEM) -r DMSO_d1 -c sara_d1 pel_d1 gef_d1 -g Drug_Day
data/DE_results/combined/$(INTRA6P_STEM)_%.csv: fig/relative_heatmap/relative_heatmap_$(INTRA6P_STEM)_%.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
relative_heatmaps: $(relative_heatmap_host) $(relative_heatmap_pathogen)
combined_pathogen_tables: $(addprefix data/DE_results/combined/$(INTRA6P_STEM)_, $(addsuffix .csv, $(combined_p_conds)))
combined_host_tables: $(addprefix data/DE_results/combined/$(INTRA6H_STEM)_, $(addsuffix .csv, $(combined_h_conds))) 


fig/combined_bar/$(INTRA6P_STEM)_fill.pdf: src/combined_bar_venn.R combined_pathogen_tables
	Rscript $< -e $(INTRA6P_STEM) -c pel_d1 gef_d1 sara_d1 -r DMSO_d1 -g Drug_Day
fig/combined_bar/$(INTRA6P_STEM)_n.pdf: fig/combined_bar/$(INTRA6P_STEM)_fill.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/combined_bar/$(INTRA6P_STEM)_vennlikely.pdf: fig/combined_bar/$(INTRA6P_STEM)_fill.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(INTRA6P_STEM)_%_d1_unique.txt: fig/combined_bar/$(INTRA6P_STEM).pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(INTRA6P_STEM)_likely_shared.txt: fig/combined_bar/$(INTRA6P_STEM).pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi

fig/combined_bar/$(INTRA6H_STEM)_fill.pdf: src/combined_bar_venn.R combined_host_tables
	Rscript $< -e $(INTRA6H_STEM) -c pel gef sara -r DMSO -g Drug

fig/combined_bar/$(INTRA6H_STEM)_n.pdf: fig/combined_bar/$(INTRA6H_STEM)_fill.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/combined_bar/$(INTRA6H_STEM)_vennlikely.pdf: fig/combined_bar/$(INTRA6H_STEM)_fill.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(INTRA6H_STEM)_%_unique.txt: fig/combined_bar/$(INTRA6H_STEM).pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
data/DE_results/$(INTRA6H_STEM)_likely_shared.txt: fig/combined_bar/$(INTRA6H_STEM).pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
combined_bar: fig/combined_bar/$(INTRA6P_STEM)_fill.pdf fig/combined_bar/$(INTRA6H_STEM)_fill.pdf

### SUPPLEMENTARY FIGURE POST PHAGOCYTOSIS 
phago_plots := $(addprefix fig/time-dependent/phago_, $(addsuffix _hm.pdf, $(addprefix $(INTRA6P_STEM)_, $(INTRA_CONDS))))
sfig_phago: $(phago_plots) fig/time-dependent/phago_$(INTRAP_STEM)_pel_d1_hm.pdf

fig/time-dependent/phago_$(INTRA6P_STEM)_%_hm.pdf: src/plotPhago.R data/DE_results/$(INTRA6P_STEM).Rds	
	Rscript $< -i $(INTRA6P_STEM) -c $*
fig/time-dependent/phago_$(INTRAP_STEM)_%_hm.pdf: src/plotPhago.R data/DE_results/$(INTRAP_STEM).Rds
	Rscript $< -i $(INTRAP_STEM) -c pel_d1 

### SUPPLEMENTARY FIGURE SIGMA AND TCS
fig/regulators/sigma_%_d1.pdf: src/plotRegulators.R
	Rscript $< -i 20240502_pel-timecourse-6donor_pathogen -c $* -r DMSO
regulators_plot := $(addprefix fig/regulators/sigma_, $(addsuffix .pdf, $(INTRA_CONDS)))
sfig_regulators: $(regulators_plot)

### SUPPLEMENTARY FIGURE GSEA 
axenic_gsea := $(addprefix fig/gsea/gsea_, $(addsuffix _$(AXENIC_STEM)_axenic.pdf, $(AXENIC_CONDS)))
intra_gsea := $(addprefix fig/gsea/gsea_, $(addsuffix _$(INTRA6P_STEM)_intra.pdf, $(INTRA_CONDS)))
fig/gsea/gsea_%_$(AXENIC_STEM)_axenic.pdf: src/plotGSEA.R data/DE_results/$(AXENIC_STEM).Rds
	Rscript $< -i data/DE_results/$(AXENIC_STEM)_$*_vs_DMSO_full.csv -g marR sigE kstR -o $*_$(AXENIC_STEM)_axenic
fig/gsea/gsea_%_$(INTRA6P_STEM)_intra.pdf: src/plotGSEA.R data/DE_results/$(INTRA6P_STEM).Rds
	Rscript $< -i data/DE_results/$(INTRA6P_STEM)_$*_vs_DMSO_d1_full.csv -g marR sigE kstR -o $*_$(INTRA6P_STEM)_intra
fig/gsea/gsea_pel_d1_$(INTRAP_STEM)_intra.pdf: src/plotGSEA.R data/DE_results/$(INTRAP_STEM).Rds
	Rscript $< -i data/DE_results/$(INTRAP_STEM)_pel_d1_vs_DMSO_d1_full.csv -g marR sigE kstR -o pel_d1_$(INTRAP_STEM)_intra
data/enrich/gsea_%.csv: fig/gsea/gsea_%.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
gsea: $(axenic_gsea) $(intra_gsea) fig/gsea/gsea_pel_d1_$(INTRAP_STEM)_intra.pdf

### SUPPLEMENTARY FIGURE HOST INFORMATION
fig/relative_heatmap/relative_heatmap_$(INTRA6H_STEM)_pel_gef_sara.pdf: src/relative_heatmap.R data/DE_results/$(INTRA6H_STEM).Rds 
	Rscript $< -i $(INTRA6H_STEM) -r DMSO -c pel gef sara -g Drug
data/DE_results/combined/$(INTRA6H_STEM)_%.csv: fig/relative_heatmap/relative_heatmap_$(INTRA6H_STEM)_%.pdf 
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/relative_heatmap/relative_heatmap_$(INTRA6H_STEM)_gef_pel_sara.pdf: src/relative_heatmap.R data/DE_results/$(INTRA6H_STEM).Rds 
	Rscript $< -i $(INTRA6H_STEM) -r DMSO -c gef pel sara -g Drug
fig/relative_heatmap/relative_heatmap_$(INTRA6H_STEM)_sara_pel_gef.pdf: src/relative_heatmap.R data/DE_results/$(INTRA6H_STEM).Rds 
	Rscript $< -i $(INTRA6H_STEM) -r DMSO -c sara pel gef -g Drug


## FIGURE 5 AXENIC EFFECTS
fig5: fig/DE_results/$(AXENIC_STEM)_volcano_pel_vs_DMSO.pdf fig/DE_results/$(AXENIC_STEM)_volcano_gef_vs_DMSO.pdf $(axenic_joint_heatmaps) iMod_enrich_axenic	

# produces plot + data table containing categorizations
axenic_joint_heatmaps := $(addprefix fig/axenic_heatmap/axenic_heatmap_, $(addsuffix _d1_$(AXENIC_STEM).pdf, $(AXENIC_CONDS)))
abx_joint_heatmaps := $(addprefix fig/axenic_heatmap/axenic_heatmap_, $(addsuffix _d1_$(ABX_STEM).pdf, $(AXENIC_CONDS)))
joint_heatmaps: $(axenic_joint_heatmaps) $(abx_joint_heatmaps) fig/axenic_heatmap/axenic_heatmap_pel_d1_$(INTRAP_STEM).pdf
fig/axenic_heatmap/axenic_heatmap_%_d1_$(AXENIC_STEM).pdf: src/axenic_heatmap.R data/DE_results/$(AXENIC_STEM).Rds
	Rscript $< -i $(INTRA6P_STEM) -a $(AXENIC_STEM) -c $*_d1 -d $* -v DMSO_d1 -w DMSO -g Drug_Day -j Drug
data/DE_results/combined/combined_intraaxenic_%.csv: fig/axenic_heatmap/axenic_heatmap_%.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
fig/axenic_heatmap/axenic_heatmap_%_d1_$(ABX_STEM).pdf: src/axenic_heatmap.R data/DE_results/$(ABX_STEM).Rds
	Rscript $< -i $(INTRA6P_STEM) -a $(ABX_STEM) -c $*_d1 -d $*_25 -v DMSO_d1 -w $*_5 -g Drug_Day -j Drug_Dose
fig/axenic_heatmap/axenic_heatmap_pel_d1_$(INTRAP_STEM).pdf: src/axenic_heatmap.R data/DE_results/$(INTRAP_STEM).Rds
	Rscript $< -i $(INTRA6P_STEM) -a $(INTRAP_STEM) -c pel_d1 -d pel_d1 -v DMSO_d1 -w DMSO_d1 -g Drug_Day -j Drug_Day

combined_intraaxenic := $(addprefix data/DE_results/combined/combined_intraaxenic_, $(addsuffix _$(AXENIC_STEM).csv, $(INTRAAXENIC_CONDS)))
combined_intraabx := $(addprefix data/DE_results/combined/combined_intraaxenic_, $(addsuffix _$(ABX_STEM).csv, $(INTRAAXENIC_CONDS)))
combined_intraaxenic_tables: $(combined_intraaxenic) $(combined_intraabx)


fig/relative_heatmap/%_intraaxenic_allcomps_iModulon.pdf: src/geneListToGSEA.R data/DE_results/$(INTRA6P_STEM)_likely_shared.txt data/DE_results/combined/combined_intraaxenic_%_$(AXENIC_STEM).csv
	Rscript $< -c $* -m intraaxenic -e $(INTRA6P_STEM) $(AXENIC_STEM)
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
	Rscript $< -c $* -m single -e $(INTRAP_STEM)
data/enrich/%_single_pel_d1.csv: fig/relative_heatmap/%_single_iModulon.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi
iMod_enrich_axenic: fig/relative_heatmap/pel_d1_intraaxenic_allcomps_iModulon.pdf fig/relative_heatmap/gef_d1_intraaxenic_allcomps_iModulon.pdf
	
fig/relative_heatmap/all_drugs_iModulon.pdf: src/geneListToGSEA.R data/DE_results/$(INTRA6P_STEM)_likely_shared.txt 
	Rscript $< -c all -m drugs -n 320 -e $(INTRA6P_STEM) 
data/enrich/all_drugs_%.csv: fig/relative_heatmap/all_drugs_iModulon.pdf
	@if test -f $@; then :; else\
		rm -f $<; \
		make $<; \
	fi

iMod_enrich_unique: fig/relative_heatmap/all_drugs_iModulon.pdf

## SUPPLEMENTARY FIGURE -- WALD TEST SCHEME
sfig_degmethod_intra: fig/deg-method/sfig6a_pel_intra.pdf fig/deg-method/sfig6b_pel_intra.pdf fig/deg-method/sfig6c_pel_intra.pdf
sfig_degmethod_axenic: fig/deg-method/sfig6a_pel_axenic.pdf
sfig_degmethod: sfig_degmethod_intra sfig_degmethod_axenic

fig/deg-method/sfig6a_%.pdf: src/testDEMethods.R data/DE_results/$(INTRA6P_STEM).Rds
	Rscript $< -i data/DE_results/$(INTRA6P_STEM).Rds -c Drug_Day_pel_d1_vs_DMSO_d1 -o $*
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
