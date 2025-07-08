import subprocess
import tempfile
import argparse
import os
import shutil
from tqdm import tqdm  # 导入tqdm库来显示进度条
import threading
import time

import glob


def monitor_folding_progress(output_file, total_sequences, pbar):
    """
    监控 RNAfold 计算进度，不再依赖 .ps 文件，而是从 output_file 统计 '>' 头的个数。

    参数:
        output_file (str): RNAfold 结果文件路径。
        total_sequences (int): 输入 RNA 序列的总数。
        pbar (tqdm): 进度条对象。
    """
    folded_rna_count = 0
    while folded_rna_count < total_sequences:
        # 统计 output_file 中 ">" 的个数
        if os.path.exists(output_file):
            with open(output_file, "r") as f:
                folded_rna_count = sum(1 for line in f if line.startswith(">"))

        pbar.update(folded_rna_count - pbar.n)  # 更新进度条
        time.sleep(1)  # 每隔1秒检查一次


def run_rnafold(data_path, data_type, rbp_name, tt):
    """
    处理 RNAfold 计算，类似于 Bash 脚本中的 RNAfold 调用，输出结果为 RNAfold 结果文本。

    参数:
        data_path (str): 数据存储的路径
        data_type (str): 数据类型（如：训练数据，测试数据等）
        rbp_name (str): RBP 名称
        tt (str): 输入文件名（不包括扩展名）
    """
    # RNAfold 命令路径
    RNAfold = "RNAfold"

    # 输入文件路径
    input_file_path = os.path.join(data_path, rbp_name, f"{data_type}_data", f"{tt}.fa")

    # 颜色代码
    GREEN = "\033[32m"
    RESET = "\033[0m"
    print(f"Start using RNAfold folding {GREEN}{input_file_path}{RESET}")

    # 创建目标文件夹
    target_directory = os.path.join(data_path, "RNAfold_results", rbp_name)
    os.makedirs(target_directory, exist_ok=True)

    # 统计输入文件中序列的数量
    with open(input_file_path, "r") as infile:
        total_sequences = sum(1 for line in infile if line.startswith(">"))

    # 结果文件路径
    output_file = os.path.join(target_directory, f"{tt}_fold.result")

    print(f"Total sequences: {total_sequences}")
    print(f"Folding sequences...")

    # 使用 tqdm 显示进度条
    with tqdm(total=total_sequences, desc="Folding RNA", unit="seq") as pbar:
        # 启动监控线程，实时追踪 output_file 中的 ">" 行数
        monitor_thread = threading.Thread(target=monitor_folding_progress, args=(output_file, total_sequences, pbar))
        monitor_thread.daemon = True
        monitor_thread.start()

        # 定义定期删除 .ps 文件的函数
        def delete_ps_files_periodically():
            while True:
                for file_path in glob.iglob(os.path.join(target_directory, "*.ps")):
                    try:
                        os.remove(file_path)
                    except Exception:
                        # 不输出任何信息，直接忽略删除失败的情况
                        pass
                # 每隔 2 秒执行一次删除操作
                import time
                time.sleep(2)

        # 启动定期删除 .ps 文件的线程
        ps_deletion_thread = threading.Thread(target=delete_ps_files_periodically)
        ps_deletion_thread.daemon = True
        ps_deletion_thread.start()

        # 执行 RNAfold 计算
        with open(input_file_path, "r") as infile, open(output_file, "w") as outfile:
            subprocess.run(
                [RNAfold, "-p", "--noPS"],
                stdin=infile,
                stdout=outfile,
                cwd=target_directory
            )

        monitor_thread.join()  # 等待监控线程结束
        # 让主线程等待一段时间，确保删除线程有机会完成最后一次删除
        import time
        time.sleep(2)

    # 输出统计信息
    print(f"Folding complete: {total_sequences}/{total_sequences} sequences folded.")
    print(f"Folding complete. Result saved to: {output_file}")

    return output_file

def run_infer_rnafold(fasta_filepath):
    """
    处理 RNAfold 计算，类似于 Bash 脚本中的 RNAfold 调用，输出结果为 RNAfold 结果文本。

    参数:
        fasta_filepath: fasta 文件的路径
    """
    RNAfold = "/data3/software/viennarna/viennarna-install-2.4.14/bin/RNAfold"
    GREEN = "\033[32m"
    RESET = "\033[0m"

    print(f"Start using RNAfold folding {GREEN}{fasta_filepath}{RESET}")

    target_directory = os.path.dirname(fasta_filepath)
    os.makedirs(target_directory, exist_ok=True)

    # 生成 RNAfold 结果文件路径
    output_file = os.path.splitext(fasta_filepath)[0] + "_fold.result"

    # 统计输入文件中序列的数量
    with open(fasta_filepath, "r") as infile:
        total_sequences = sum(1 for line in infile if line.startswith(">"))

    print(f"Total sequences to fold: {total_sequences}")

    # 使用 tqdm 显示进度条
    with tqdm(total=total_sequences, desc="Folding RNA", unit="seq") as pbar:
        # 启动监控线程，实时追踪 output_file 中的 ">" 行数
        monitor_thread = threading.Thread(target=monitor_folding_progress, args=(output_file, total_sequences, pbar))
        monitor_thread.daemon = True
        monitor_thread.start()

        # 启动一个新线程来定期删除 .ps 文件
        def delete_ps_files_periodically():
            while True:
                for file_path in glob.iglob(os.path.join(target_directory, "*.ps")):
                    try:
                        os.remove(file_path)
                    except Exception:
                        # 不输出任何信息，直接忽略删除失败的情况
                        pass
                # 每隔 2 秒执行一次删除操作
                import time
                time.sleep(2)

        ps_deletion_thread = threading.Thread(target=delete_ps_files_periodically)
        ps_deletion_thread.daemon = True
        ps_deletion_thread.start()

        # 执行 RNAfold 计算
        with open(fasta_filepath, "r") as infile, open(output_file, "w") as outfile:
            subprocess.run(
                [RNAfold, "-p", "--noPS"],
                stdin=infile,
                stdout=outfile,
                cwd=target_directory
            )

        monitor_thread.join()  # 等待监控线程结束
        # 虽然线程是守护线程，但这里可以让主线程等待一小段时间确保删除线程有机会执行最后一次删除
        import time
        time.sleep(2)

    print(f"Folding complete: {total_sequences}/{total_sequences} sequences folded.")
    print(f"Folding complete. Result saved to: {output_file}")

    return output_file

    
def annotate_str(mfe_structure_line):
    """
    注释单个 MFE 结构行，保存到临时文件后调用 C 程序进行处理，并提取注释结果。

    参数:
        mfe_structure_line (str): 点括号格式的 MFE 结构。

    返回:
        str: 注释后的字符串（C 程序输出）。
    """
    try:
        # 使用 tempfile 创建临时文件
        with tempfile.NamedTemporaryFile(delete=True, mode="w") as structure_file, \
             tempfile.NamedTemporaryFile(delete=True, mode="r") as annotation_file:

            # 将 MFE 结构保存到临时文件
            structure_file.write(mfe_structure_line + "\n")
            structure_file.flush()  # 确保写入文件

            # 调用 C 程序对临时文件进行注释
            c_program = "./data_gerenate/fold_str_annotion/parse_secondary_structure_v2"  # C 程序的路径
            subprocess.run(
                [c_program, structure_file.name, annotation_file.name],
                check=True)
            # 从注释文件中读取结果
            annotated_structure = annotation_file.read().strip()

        return annotated_structure

    except subprocess.CalledProcessError as e:
        print(f"Error running C program: {e}")
        return None

    except Exception as e:
        print(f"Unexpected error: {e}")
        return None

def process_rnafold_and_annotate(input_file, output_file):
    """
    从 RNAfold 输出文件中提取第 1 行（标题行）、第 2 行（RNA 序列行）和第 3 行（MFE 结构行），
    对结构进行注释，并将结构转换为 2-letter 格式，然后将结果保存为 TSV 格式。

    参数:
        input_file (str): 输入 RNAfold 输出文件路径。
        output_file (str): 输出处理后结果文件路径（应以 .tsv 结尾）。
    """
    # 定义7-letter到4-letter和4-letter到2-letter结构的转换字典
    struct_7_convert_4 = {'B': 'M', 'E': 'U', 'H': 'L', 'L': 'P', 'M': 'M', 'R': 'P', 'T': 'M'}
    struct_4_convert_2 = {'P': 'P', 'L': 'U', 'U': 'U', 'M': 'U'}

    # 检查输出目录是否存在，不存在则创建
    dirname = os.path.dirname(output_file)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    # 打开输入文件并处理每个 RNA 记录块
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        lines = infile.readlines()

        # 每个 RNA 记录块的行数为 6 行
        for i in range(0, len(lines), 6):
            # 检查当前块是否有足够的行
            if i + 2 < len(lines):
                header_line = lines[i].strip()  # 第 1 行
                sequence_line = lines[i + 1].strip()  # 第 2 行
                mfe_structure_line = lines[i + 2].strip().split()[0]  # 第 3 行，去掉能量值

                # 对结构进行注释（调用注释函数）
                annotated_structure = annotate_str(mfe_structure_line)

                # 转换结构：7-letter -> 4-letter -> 2-letter
                struct_4 = []
                struct_2 = []

                for i in range(len(annotated_structure)):
                    # Convert the 7-letter structure to 4-letter and then to 2-letter
                    struct_4.append(struct_7_convert_4[annotated_structure[i]])
                    struct_2.append(struct_4_convert_2[struct_7_convert_4[annotated_structure[i]]])

                # Join the results into strings
                struct_4 = ''.join(struct_4)
                struct_2 = ''.join(struct_2)

                # Write the final results: header, sequence, original structure, annotated structure
                outfile.write(f"{header_line}\t{sequence_line}\t{struct_2}\n")

    RED = "\033[31m"
    RESET = "\033[0m"

    print(f"Sequence and Structure Annotation Files {RED}{output_file}{RESET} has been saved.")
    
    os.remove(input_file)
