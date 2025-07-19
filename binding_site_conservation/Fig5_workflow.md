# Cross-species RNA-binding Protein Binding Site Conservation Analysis Workflow

## Cross-species homologous transcript screening

```bash
nohup bash script/mashmap.sh & # Use mashmap for rapid cross-species homologous transcript search

nohup bash script/uniq.sh & # Remove duplicates from mashmap alignment results, count conserved transcript pairs between two species, use human as reference, count transcript ids
```

## HAR entry filtering

```bash
nohup bash script/grep.sh & # After alignment, filter cross-species RBP binding sites based on core transcripts

nohup bash script/filter.sh & # Filter entries, keep only cross-species conserved RBP entries, core to core.filted

merge_bed.sh # Merge bed files

grep_bed.sh # Filter bed files at transcript level
```

## Homologous sequence alignment

```bash
cat /data1/MuSIC2_new_model/MotifLearnModel/11species_RNA_raw_fasta/*/!(*cut200nt*).rna > refs.rna 

Rscript script/extract.R -i core.filted.tsv -f refs.rna -o seqs_test/

bash script/align.sh 
```

## Statistical analysis of binding site conservation

```bash
bash script/count.sh

nohup bash script/count_box.sh &

nohup bash script/count_heatmap.sh &

Rscript script/plot_side.R -H count/count_heatmap.txt -B count/count_box.txt -c 2 -o figs/count_side_v2.pdf 

Rscript script/plot_side_cluster.R -H count/count_heatmap.txt -B count/count_box.txt -c 2 -o figs/count_side_v3.pdf # Clustering analysis
```

## Binding site examples

```bash
nohup bash script/example_filter.sh &

nohup bash script/focus.sh &

nohup bash script/update_bed_focus.sh &

nohup bash script/igv_ref.sh &

nohup bash script/igv_peak.sh &
```

## RBP family conservation

```bash
#Rscript script/plot_box.R -i count/count_box.txt -c 2 -o figs/count_box_col2.pdf

Rscript script/plot_box_family.R -i count/count_box.txt -c 2 -p '^YTHD' -o figs/count_box_famaliy_YTHD.pdf

Rscript script/plot_box_family.R -i count/count_box.txt -c 2 -p '^CPSF' -o figs/count_box_famaliy_CPSF.pdf
```

## Species cumulative conservation

```bash
nohup bash script/focus_CPSF.sh & ## har to bed

nohup bash script/update_bed_focus_CPSF.sh & ## bed cross-species align site update

nohup bash script/count_box_focus_CPSF.sh &

Rscript script/plot_box_family.R -i update_bed_focus/all/count_box_famaliy_CPSF.txt -c 2 -p '^CPSF' -o figs/count_box_famaliy_CPSF_all.pdf

nohup bash script/count_box_focus_CPSF_extend80.sh &

Rscript script/plot_box_family.R -i update_bed_focus/all/count_box_focus_CPSF_extend80.txt -c 2 -p '^CPSF' -o figs/count_box_famaliy_CPSF_all_extend80.pdf

nohup bash /home/bioinfo/02_project/05_RBP/05_homo/script/count_box_focus_CPSF_extend80_species_1.sh &

nohup bash /home/bioinfo/02_project/05_RBP/05_homo/script/count_box_focus_CPSF_extend80_species_2.sh &

nohup bash /home/bioinfo/02_project/05_RBP/05_homo/script/count_box_focus_CPSF_extend80_species_3.sh &

nohup bash /home/bioinfo/02_project/05_RBP/05_homo/script/count_box_focus_CPSF_extend80_species_4.sh &

nohup bash /home/bioinfo/02_project/05_RBP/05_homo/script/count_box_focus_CPSF_extend80_species_5.sh &

nohup bash /home/bioinfo/02_project/05_RBP/05_homo/script/count_box_focus_CPSF_extend80_species_6.sh &

nohup bash /home/bioinfo/02_project/05_RBP/05_homo/script/count_box_focus_CPSF_extend80_species_7.sh &

nohup bash /home/bioinfo/02_project/05_RBP/05_homo/script/count_box_focus_CPSF_extend80_species_8.sh &

nohup bash /home/bioinfo/02_project/05_RBP/05_homo/script/count_box_focus_CPSF_extend80_species_9.sh &

nohup bash /home/bioinfo/02_project/05_RBP/05_homo/script/count_box_focus_CPSF_extend80_species_10.sh &

Rscript script/plot_forest_family_ggplot.R -o figs/count_forest_famaliy_CPSF_all_extend80_ggplot.pdf -c 2 -i update_bed_focus/all/count_box_famaliy_CPSF_extend80_species_1.txt,update_bed_focus/all/count_box_famaliy_CPSF_extend80_species_2.txt,update_bed_focus/all/count_box_famaliy_CPSF_extend80_species_3.txt,update_bed_focus/all/count_box_famaliy_CPSF_extend80_species_4.txt,update_bed_focus/all/count_box_famaliy_CPSF_extend80_species_5.txt,update_bed_focus/all/count_box_famaliy_CPSF_extend80_species_6.txt,update_bed_focus/all/count_box_famaliy_CPSF_extend80_species_7.txt,update_bed_focus/all/count_box_famaliy_CPSF_extend80_species_8.txt,update_bed_focus/all/count_box_famaliy_CPSF_extend80_species_9.txt,update_bed_focus/all/count_box_famaliy_CPSF_extend80_species_10.txt
```
