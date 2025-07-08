import h5py

h5_file = "data/predict_data/test.h5"
# 打开 .h5 文件
with h5py.File(h5_file, "r") as h5f:
    rna_names = h5f["rna_names"][:]  # 加载 RNA 名称
    one_hot_matrices = h5f["one_hot_matrices"][:]  # 加载 One-hot 矩阵

    # 打印数据
    print("RNA Names:", [name.decode('utf-8') for name in rna_names])  # 转换为字符串
    print("One-hot Matrices Shape:", one_hot_matrices.shape)