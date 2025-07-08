import argparse
from annotation_tools import *
from data_utils import load_data_save_h5
import os


def process_rnafold_infer_data(fasta_filepath):
    """
    处理没有label的数据集，进行inference 或者 计算 HAR 等等
    """
    # 执行 RNAfold
    rnafold_result = "/data3/MuSIC2_new_model/MotifLearnModel/MuSIC/data/11species_all_trans_predict_data/HUMAN/HUMAN_cut200nt_all_fold.result"
    # 生成注释文件的路径
    seq_str_anno_file = os.path.splitext(fasta_filepath)[0] + "_annotation.tsv"

    print("Begin to process RNAfold result and generate annotation file.")

    process_rnafold_and_annotate(rnafold_result, seq_str_anno_file)
    output_file = os.path.splitext(fasta_filepath)[0] + ".h5"
    load_data_save_h5(seq_str_anno_file, output_file, max_length=200)

    print("RNAfold processing, annotation complete and generate h5 file for MuSIC.")


process_rnafold_infer_data("/data3/MuSIC2_new_model/MotifLearnModel/MuSIC/data/11species_all_trans_predict_data/HUMAN/HUMAN_cut200nt_all.fa")