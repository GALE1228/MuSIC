

# train cross-species RBP model for your species
python main.py --train --rbp_name TARDBP_HUMAN --file_path data/186rbp_dataset --gpuid 0 --cross --species_name mouse --smooth_rate 0.64

# predict variants using trained model
python main.py --infer --rbp_name TARDBP_HUMAN --infer_fasta_path mouse_variants.fa --gpuid 0 --cross --species_name mouse --smooth_rate 0.64

# predict wt using trained model
python main.py --infer --rbp_name TARDBP_HUMAN --infer_fasta_path mouse_wt.fa --gpuid 0 --cross --species_name mouse --smooth_rate 0.64

python SNV_impact_prediction/merge_snv_impact.py 
    \ --variants TARDBP_HUMAN_cross_mouse_music_mouse_variants.inference
    \ --wt TARDBP_HUMAN_cross_mouse_music_mouse_wt.inference
    \ --output SNV_impact_prediction