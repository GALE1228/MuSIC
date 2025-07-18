import h5py

h5_file = "data/predict_data/test.h5"
# Open .h5 file
with h5py.File(h5_file, "r") as h5f:
    rna_names = h5f["rna_names"][:]  # Load RNA names
    one_hot_matrices = h5f["one_hot_matrices"][:]  # Load one-hot matrices

    # Print data
    print("RNA Names:", [name.decode('utf-8') for name in rna_names])  # Convert to string
    print("One-hot Matrices Shape:", one_hot_matrices.shape)