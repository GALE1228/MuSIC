


## 跨物种同源转录本筛选
nohup bash script/mashmap.sh & #使用mashmap进行跨物种同源转录本快速搜索
nohup bash script/uniq.sh & # mashmap比对结果去重复，统计两物种保守的转录本对，以人为基准，统计转录本id

## HAR条目筛选
nohup bash script/grep.sh & # align 之后，根据core转录本筛选跨物种的RBP结合位点
nohup bash script/filter.sh & # 过滤条目，只保留跨物种保守的RBP条目，core to core.filted
merge_bed.sh # 合并bed文件
grep_bed.sh # 筛选bed文件为各转录本水平

## 同源序列比对
cat /data1/MuSIC2_new_model/MotifLearnModel/11species_RNA_raw_fasta/*/!(*cut200nt*).rna > refs.rna 
Rscript script/extract.R -i core.filted.tsv -f refs.rna -o seqs_test/
bash script/align.sh 

## 统计结合位点保守性
bash script/count.sh
nohup bash script/count_box.sh &
nohup bash script/count_heatmap.sh &
Rscript script/plot_side.R -H count/count_heatmap.txt -B count/count_box.txt -c 2 -o figs/count_side_v2.pdf 
Rscript script/plot_side_cluster.R -H count/count_heatmap.txt -B count/count_box.txt -c 2 -o figs/count_side_v3.pdf # 聚类分析

## 结合位点示例
nohup bash script/example_filter.sh &
nohup bash script/focus.sh &
nohup bash script/update_bed_focus.sh &
nohup bash script/igv_ref.sh &
nohup bash script/igv_peak.sh &

## RBP家族保守性
#Rscript script/plot_box.R -i count/count_box.txt -c 2 -o figs/count_box_col2.pdf
Rscript script/plot_box_family.R -i count/count_box.txt -c 2 -p '^YTHD' -o figs/count_box_famaliy_YTHD.pdf
Rscript script/plot_box_family.R -i count/count_box.txt -c 2 -p '^CPSF' -o figs/count_box_famaliy_CPSF.pdf

## 物种累加保守性
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
