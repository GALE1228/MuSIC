import argparse
from .annotation_tools import *
from .data_utils import load_data_save_h5
import os

def process_train_rnafold_data(data_path, rbp_name):
    """
    处理 RNAfold 数据，并执行折叠和注释。
    
    参数:
        data_path (str): 数据所在路径
        rbp_name (str): RBP 名称，如 "AGO2_MiClip"
    """
    # 数据类型（positive 和 negative）和任务类型（train 和 test）
    dts = ["positive", "negative"]

    # 遍历所有数据类型和任务类型的组合
    for dt in dts:
        # 执行 RNAfold
        rnafold_result = run_rnafold(data_path, dt, rbp_name, tt = "train")
        # 生成注释文件的路径
        seq_str_anno_file = f"{data_path}/{rbp_name}/{dt}_data/train_annotation.tsv"
        # 进行 RNAfold 结果的处理和注释
        process_rnafold_and_annotate(rnafold_result, seq_str_anno_file)
        output_file = f"{data_path}/{rbp_name}/{dt}_data/train.h5"
        load_data_save_h5(seq_str_anno_file, output_file, max_length=200)

    print("RNAfold processing, annotation complete and generate h5 file for MuSIC.")

def process_validation_rnafold_data(data_path, rbp_name):
    """
    处理 RNAfold 数据，并执行折叠和注释。
    
    参数:
        data_path (str): 数据所在路径
        rbp_name (str): RBP 名称，如 "AGO2_MiClip"
    """
    # 数据类型（positive 和 negative）和任务类型（train 和 test）
    dts = ["positive", "negative"]

    # 遍历所有数据类型和任务类型的组合
    for dt in dts:
        # 执行 RNAfold
        rnafold_result = run_rnafold(data_path, dt, rbp_name ,tt = "test")
        # 生成注释文件的路径
        seq_str_anno_file = f"{data_path}/{rbp_name}/{dt}_data/test_annotation.tsv"

        print("Begin to process RNAfold result and generate annotation file.")

        process_rnafold_and_annotate(rnafold_result, seq_str_anno_file)
        output_file = f"{data_path}/{rbp_name}/{dt}_data/test.h5"
        load_data_save_h5(seq_str_anno_file, output_file, max_length=200)

    print("RNAfold processing, annotation complete and generate h5 file for MuSIC.")


def process_rnafold_infer_data(fasta_filepath):
    """
    处理没有label的数据集，进行inference 或者 计算 HAR 等等
    """
    # 执行 RNAfold
    rnafold_result = run_infer_rnafold(fasta_filepath)
    # 生成注释文件的路径
    seq_str_anno_file = os.path.splitext(fasta_filepath)[0] + "_annotation.tsv"

    print("Begin to process RNAfold result and generate annotation file.")

    process_rnafold_and_annotate(rnafold_result, seq_str_anno_file)
    output_file = os.path.splitext(fasta_filepath)[0] + ".h5"
    load_data_save_h5(seq_str_anno_file, output_file, max_length=200)

    print("RNAfold processing, annotation complete and generate h5 file for MuSIC.")
