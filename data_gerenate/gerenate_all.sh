for p in $(cat data/cross_species_rbp.list); do
    echo $p
    python data_gerenate/RNAfold_annotation_gerenate_h5.py --data_path "data/train_validate_data" --rbp_name $p
done

