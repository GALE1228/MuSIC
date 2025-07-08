import torch
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, TensorDataset
import pandas as pd
import os
from .one_hot_encode_decode import *
import h5py


# 检查one-hot编码是否正确
def load_data_test(data_file, matrices_to_concatenate, max_length=200):
    df = pd.read_csv(data_file, sep="\t", header=None)
    row = df.iloc[0]
    rna_name = row.iloc[0]
    
    # 根据需要拼接的矩阵构建One-hot编码
    one_hot_matrices = []
    
    # 根据命令行输入选择需要拼接的矩阵
    if 'seq_4' in matrices_to_concatenate:
        seq_4 = row.iloc[1]
        seq_4_oh = convert_one_hot_seq_4(seq_4 ,max_length)
        seq_4_deoh = decode_seq_4(seq_4_oh)
        print(seq_4_deoh == seq_4)
        one_hot_matrices.append(seq_4_oh)

    if 'str_2' in matrices_to_concatenate:
        str_2 = row.iloc[2]
        str_2_oh = convert_one_hot_str_2(str_2 ,max_length)
        str_2_deoh = decode_str_2(str_2_oh)
        print(str_2_deoh == str_2)
        one_hot_matrices.append(str_2_oh)

    if 'str_4' in matrices_to_concatenate:
        str_4 = row.iloc[3]
        str_4_oh = convert_one_hot_str_4(str_4 ,max_length)
        str_4_deoh = decode_str_4(str_4_oh)
        print(str_4_deoh == str_4)
        one_hot_matrices.append(str_4_oh)

    if 'str_7' in matrices_to_concatenate:
        str_7 = row.iloc[4]
        str_7_oh = convert_one_hot_str_7(str_7 ,max_length)
        str_7_deoh = decode_str_7(str_7_oh)
        print(str_7_deoh == str_7)
        one_hot_matrices.append(str_7_oh)

    if 'seq_str_8' in matrices_to_concatenate:
        seq_str_8 = row.iloc[5]
        seq_str_8_oh = convert_one_hot_seq_str_8(seq_str_8 ,max_length)
        seq_str_8_deoh = decode_seq_str_8(seq_str_8_oh)
        print(seq_str_8_deoh == seq_str_8)
        one_hot_matrices.append(seq_str_8_oh)

    if 'seq_str_16' in matrices_to_concatenate:
        seq_str_16 = row.iloc[6]
        seq_str_16_oh = convert_one_hot_seq_str_16(seq_str_16 ,max_length)
        seq_str_16_deoh = decode_seq_str_16(seq_str_16_oh)
        print(seq_str_16_deoh == seq_str_16)
        one_hot_matrices.append(seq_str_16_oh)

    if 'seq_str_28' in matrices_to_concatenate:
        seq_str_28 = row.iloc[7]
        seq_str_28_oh = convert_one_hot_seq_str_28(seq_str_28 ,max_length)
        seq_str_28_deoh = decode_seq_str_28(seq_str_28_oh)
        print(seq_str_28_deoh == seq_str_28)
        one_hot_matrices.append(seq_str_28_oh)

    # 拼接选择的One-hot矩阵
    combined_oh = combine_one_hot_matrix(one_hot_matrices, axis=0)
    # print(one_hot_matrices)
    
    # 输出每个矩阵的形状以及解码验证
    print(f"0\t{rna_name}\t" + "\t".join([str(mat.shape) for mat in one_hot_matrices]))
    
    
    return combined_oh

def load_data_save_h5(data_file, output_file, max_length):
    """
    该函数读取数据文件，将其中的 `seq_4` 和 `str_2` 特征进行 one-hot 编码，
    然后将它们合并并保存为 HDF5 格式。

    参数:
        data_file (str): 输入的数据文件路径，包含 `seq_4` 和 `str_2` 列。
        output_file (str): 输出文件路径，将保存为 HDF5 格式。
        max_length (int): 用于 one-hot 编码时序列的最大长度。
    """
    
    dirname = os.path.dirname(output_file)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    
    # 读取数据
    df = pd.read_csv(data_file, sep="\t", header=None)
    
    all_combined_matrices = []  # 用于保存每一行的合并后的矩阵
    all_rna_names = []  # 保存每一行的 RNA 名称
    
    for index, row in df.iterrows():
        rna_name = row.iloc[0]
        one_hot_matrices = []

        seq_4 = row.iloc[1]
        seq_4_oh = convert_one_hot_seq_4(seq_4, max_length)
        one_hot_matrices.append(seq_4_oh)

        str_2 = row.iloc[2]
        str_2_oh = convert_one_hot_str_2(str_2, max_length)
        one_hot_matrices.append(str_2_oh)

        combined_oh = combine_one_hot_matrix(one_hot_matrices, axis=0)

        all_rna_names.append(rna_name)
        all_combined_matrices.append(combined_oh)

    print("combined feature matrix shape after one-hot encode is ", np.array(all_combined_matrices).shape)
    GREEN = "\033[32m"
    RESET = "\033[0m"

    # 保存到 HDF5 文件
    with h5py.File(output_file, "w") as h5f:
        h5f.create_dataset("rna_names", data=np.array(all_rna_names, dtype='S'))
        h5f.create_dataset("one_hot_matrices", data=np.array(all_combined_matrices, dtype=np.float32))
    print(f"Data saved to {GREEN}{output_file}{RESET}")
