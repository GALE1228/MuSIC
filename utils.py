import torch
from sklearn.model_selection import train_test_split
from torch.utils.data import Dataset, DataLoader, TensorDataset
import pandas as pd
import os
import h5py
import numpy as np
from data_gerenate.RNAfold_annotation_gerenate_h5 import *

def init_fasta_headers(fasta_path, replacement='|'):
    import re

    lines = []
    with open(fasta_path, 'r') as infile:
        for line in infile:
            if line.startswith('>'):
                clean_line = re.sub(r'[\s\t]+', replacement, line.strip())
                lines.append(clean_line + '\n')
            else:
                lines.append(line)

    with open(fasta_path, 'w') as outfile:
        outfile.writelines(lines)

    return fasta_path

def load_h5_file(h5_file):
    with h5py.File(h5_file, "r") as h5f:
        rna_names = h5f["rna_names"][:]  # 加载 RNA 名称
        one_hot_matrices = h5f["one_hot_matrices"][:]  # 加载 One-hot 矩阵

    return rna_names,one_hot_matrices

def split_dataset(data ,labels ,test_size):
    X_train, X_test, y_train, y_test = train_test_split(data, labels, test_size=test_size, random_state=42 ,stratify=labels)
    return X_train, X_test, y_train, y_test

def create_dataloader(X, y, y_smooth, batch_size):
    tensor_x = torch.Tensor(X)
    tensor_y = torch.Tensor(y)
    tensor_y_smooth = torch.Tensor(y_smooth)
    dataset = TensorDataset(tensor_x, tensor_y, tensor_y_smooth)
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)
    
    return dataloader

def create_infer_dataloader(X, y, rna_names, batch_size=32):
    tensor_x = torch.Tensor(X)  # 转换为Tensor
    tensor_y = torch.Tensor(y)  # 转换为Tensor
    dataset = TensorDataset(tensor_x, tensor_y)
    
    # 创建数据加载器
    dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=False)

    # 返回数据加载器和 rna_names
    return dataloader, rna_names

def load_model(model, model_path, device):
    model.load_state_dict(torch.load(model_path))
    model.to(device)
    model.eval()
    return model

def make_directory(path, foldername, verbose=1):
    """make a directory"""

    if not os.path.isdir(path):
        os.mkdir(path)
        print("making directory: " + path)

    outdir = os.path.join(path, foldername)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
        print("making directory: " + outdir)
    return outdir

def smooth_onehot_label(one_hot_label, smooth_rate):

    total_count = one_hot_label.size(0)
    class_1_distribution = (one_hot_label == 1).sum().item() / total_count
    smoothed_labels = torch.zeros_like(one_hot_label)
    for i in range(one_hot_label.size(0)):
        if one_hot_label[i] == 1:
            smoothed_labels[i] = smooth_rate * one_hot_label[i] + (1 - smooth_rate) * class_1_distribution
        else:
            smoothed_labels[i] = smooth_rate * one_hot_label[i] + (1 - smooth_rate) * (1 - class_1_distribution)

    smoothed_labels = smoothed_labels.cpu().numpy()
    # print("smoothed_labels:", smoothed_labels[1:5])
    return smoothed_labels

def train_dataset(file_path, rbp_name, batch_size, smooth_rate):
    # H5 文件路径
    train_positive_h5_file = f"{file_path}/{rbp_name}/positive_data/train.h5"
    train_negative_h5_file = f"{file_path}/{rbp_name}/negative_data/train.h5"
    
    # 检查是否存在 H5 文件，如果不存在则调用 process_rnafold_data
    if not os.path.exists(train_positive_h5_file) or not os.path.exists(train_negative_h5_file):
        print(f"{train_positive_h5_file} or {train_negative_h5_file} not found, generating H5 files.")
        process_train_rnafold_data(file_path, rbp_name)

    # 加载正负数据集
    train_positive_rna_names, train_positive_dataset = load_h5_file(train_positive_h5_file)
    train_negative_rna_names, train_negative_dataset = load_h5_file(train_negative_h5_file)

    # 合并数据集
    train_data = np.concatenate((train_positive_dataset, train_negative_dataset), axis=0)

    # 创建训练集和测试集的标签
    train_labels = np.concatenate((np.ones(len(train_positive_dataset)), np.zeros(len(train_negative_dataset))), axis=0)

    print(train_labels.shape)
    smoothed_label = smooth_onehot_label(torch.tensor(train_labels) , smooth_rate)
    print("smoothed_label shape:", smoothed_label.shape)

    # 划分训练集和测试集
    X_train, X_test, y_train, y_test, y_train_smooth, y_test_smooth = train_test_split(
    train_data, train_labels, smoothed_label, test_size=0.2, random_state=42
    )

    train_loader = create_dataloader(X_train, y_train, y_train_smooth, batch_size)
    test_loader = create_dataloader(X_test, y_test, y_test_smooth, batch_size)


    return train_loader, test_loader

def validation_dataset(file_path , rbp_name , batch_size , smooth_rate):
    
    test_positive_h5_file = f"{file_path}/{rbp_name}/positive_data/test.h5"
    test_negative_h5_file = f"{file_path}/{rbp_name}/negative_data/test.h5"

    if not os.path.exists(test_positive_h5_file) or not os.path.exists(test_negative_h5_file):
        print(f"{test_positive_h5_file} or {test_negative_h5_file} not found, generating H5 files.")
        process_validation_rnafold_data(file_path, rbp_name)

    test_positive_rna_names , test_positive_dataset = load_h5_file(test_positive_h5_file)
    test_negative_rna_names , test_negative_dataset = load_h5_file(test_negative_h5_file)

    test_data = np.concatenate((test_positive_dataset, test_negative_dataset), axis=0)
    test_labels = np.concatenate((np.ones(len(test_positive_dataset)), np.zeros(len(test_negative_dataset))), axis=0)

    print(test_labels.shape)
    smoothed_label = smooth_onehot_label(torch.tensor(test_labels) , smooth_rate)
    print("smoothed_label shape:", smoothed_label.shape)

    data_loader = create_dataloader(test_data, test_labels, smoothed_label , batch_size)

    return data_loader

def inference_dataset(fasta_path , batch_size):
    print(os.path.basename(fasta_path))
    fasta_path = init_fasta_headers(fasta_path)
    inference_h5_file = os.path.splitext(fasta_path)[0] + ".h5"
    print(inference_h5_file)

    if not os.path.exists(inference_h5_file):
        print(f"{inference_h5_file} not found, generating H5 files.")
        process_rnafold_infer_data(fasta_path)

    inference_rna_names , inference_dataset = load_h5_file(inference_h5_file)

    test_data = inference_dataset
    test_labels = np.ones(len(inference_dataset))

    # 创建数据加载器
    data_loader, rna_names_all = create_infer_dataloader(test_data, test_labels, inference_rna_names, batch_size)

    return data_loader, rna_names_all

def save_validations(out_dir, filename, dataname, predictions, label, met):
    evals_dir = make_directory(out_dir, f"out/evals")
    metrics_path = os.path.join(evals_dir, filename+'.metrics')
    probs_path = os.path.join(evals_dir, filename+'.probs')
    with open(metrics_path,"w") as f:
        print("{:s}\t{:.3f}\t{:.3f}\t{:.3f}\t{:d}\t{:d}\t{:d}\t{:d}".format(
            dataname,
            met.acc,
            met.auc,
            met.prc,
            met.tp,
            met.tn,
            met.fp,
            met.fn,
        ), file=f)
    with open(probs_path,"w") as f:
        for i in range(len(predictions)):
            print("{:.3f}\t{}".format(predictions[i], label[i]), file=f)
    print("Evaluation file:", metrics_path)
    print("Prediction file:", probs_path)

def save_infers(out_dir, filename, rna_names_all, y_all, p_all):
    # 创建目录
    evals_dir = make_directory(out_dir, "out/infer")
    probs_path = os.path.join(evals_dir, filename + '.inference')
    
    # 打开文件进行写入
    with open(probs_path, "w") as f:
        for i in range(len(y_all)):
            rna_name = rna_names_all[i]
            
            # 检查是否是字节串，如果是，则解码
            if isinstance(rna_name, bytes):
                rna_name = rna_name.decode('utf-8')
            
            y_line = "{:f}".format(y_all[i])
            p_line = "{:f}".format(p_all[i][0])
            f.write(f"{rna_name}\t{y_line}\t{p_line}\n")

    print(f"Prediction file saved to: {probs_path}")

