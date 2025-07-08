import numpy as np

"""
# 使用示例
rna_sequences = ["CAGGGUCAGGAUCCGCGUGCUUGCUAGUCCACCUUACCGCCGGAUGGGCAGCCGCUGCUGAGCGCCUGCUAGGUGGUCCGACGCGCAGCUCCGAGCUGUCCUGGACUGGCGCUUUUUAGGAGCGAAGGGGAACCCCGCAGGAGACCAGGGCCCUGAACUCAGGGGCUUCGCCACUGAUUGUCCAAACGCAAUUCUUGUACGAGUCUGCGGCCAACCGAGAAUUGUGGCUGGACAUCUGUGGCUGAGCUCCGGGCGCAACAGGGGCGGGGGCCCCAGGGACAGGGCUCAGCGCGGGCGAGACCUCUCGGGGCGGCGGCUGGCUGGCCCCAGCGCGAGGAUCGCGGUCCCGGCCCGCGCGCACAGAGACGCCGGUCCCUGCCACCAGCGCCGCCAUCCGCAU",
                 "AAGCGUUCAAGCUCAACACCCACUACCUAAAAAAUCCCAAACAUAUAACUGAACUCCUCACACCCAAUUGGACCAAUCUAUCACCCUAUAGAAGAACUAAUGUUAGUAUAAGUAACAUGAAAACAUUCUCCUCCGCAUAAGCCUGCGUCAGAUUAAAACACUGAACUGACAAUUAACAGCCCAAUAUCUACAAUCAACCAACAAGUCAUUAUUACCCUCACUGUCAACCCAACACAGGCAUGCUCAUAAGGAAAGGUUAAAAAAAGUAAAAGGAACUCGGCAAAUCUUACCCCGCCUGUUUACCAAAAACAUCACCUCUAGCAUCACCAGUAUUAGAGGCACCGCCUGCCCAGUGACACAUGUUUAACGGCCGCGGUACCCUAACCGUGCAAAGGUAGCA"]
max_length = 400

# 1. 将RNA序列转换为one-hot编码
one_hot_encoded = convert_one_hot(rna_sequences, max_length=max_length)
print("One-hot编码结果：")
print(one_hot_encoded)
print(one_hot_encoded)

# 2. 解码回RNA序列
decoded_sequences = [decodeRNA(seq) for seq in one_hot_encoded]
print("\n解码后的RNA序列：")
print(decoded_sequences)

"""
# seq_4
def convert_one_hot_seq_4(sequence, max_length):
    # 假设sequence是一个单一的RNA序列字符串
    seq = sequence.upper()  # 转换为大写
    seq_length = len(seq)
    
    if seq_length > max_length:
        seq = seq[:max_length]
        seq_length = max_length
    
    # 创建一个4行seq_length列的零矩阵
    one_hot = np.zeros((4, seq_length))  # 4个通道，分别对应 A, C, G, U

    # 为每个核苷酸分配one-hot编码
    index = [j for j in range(seq_length) if seq[j] == 'A']
    one_hot[0, index] = 1
    index = [j for j in range(seq_length) if seq[j] == 'C']
    one_hot[1, index] = 1
    index = [j for j in range(seq_length) if seq[j] == 'G']
    one_hot[2, index] = 1
    index = [j for j in range(seq_length) if seq[j] == 'U']  # RNA中使用U代替T
    one_hot[3, index] = 1

    # 如果给定了max_length，进行零填充
    if max_length:
        offset1 = int((max_length - seq_length) / 2)
        offset2 = max_length - seq_length - offset1
        if offset1:
            one_hot = np.hstack([np.zeros((4, offset1)), one_hot])
        if offset2:
            one_hot = np.hstack([one_hot, np.zeros((4, offset2))])

    return one_hot  # 返回一个单一的(4, max_length)矩阵
# seq_4
def decode_seq_4(m):
    na = ["A", "C", "G", "U"]  # RNA中的核苷酸
    seq = ""
    for i in range(m.shape[1]):  # 遍历每列
        # 查找该列中值为1的行索引（即对应核苷酸的位置）
        var = np.where(m[:, i] == 1)[0]
        if len(var) == 1:  # 确保找到唯一的匹配
            seq = seq + na[var[0]]  # 将对应的核苷酸添加到序列中
    return seq

"""

# 使用示例
sequences = ["UUUUUPPPPPUUUUUUUUPPUPPPPPPPUUPPPPUPPPPPPUUPPUPPPPUUUUUPPPPPPPPPPUUUPPPPPUUUUUUPPPPPUUUPPPPUUUUUPPPPPPUUUPPPPUPPUUUUPPPPPPUPPPPUPPPPUPPPPPPUUUUUUUPUUUUUUUUPUUUUUUPPPPPPUPPPPPPPPPPUUPPPPPPPUUUUPPPPPPPUUPPPPPPUUUUUPPPPPPPPUUUUUPPPPPPPPPPPPPPPUPPPPPPPUPPPPUPPUUUUUUUUUUUUPPUPPPPPPPPPPPPPPUUUUUUUPPPPPPPUUUUUUUUUUUUUPPPPPPPPPPUUUUUUUUPPPPPUUUUUUUUUUPPPPPUUUUUUUPPPPPPUUUPPPPPUUUUUUUUUUPPPPPPPPPPPPPPPPPPU",
             "UUPPPPPPPPPPPPPUUUUUUUPPPPUPPPUUPPPPPPUPPPUUPPPPPPPPPPPUUUUUPPPPPPPUPPPUUUPPPUPPPPUPPPPPPPPPPPUUPPPPUPPUUPPPPPUUUPPPPPPPPPPUUUUPPPPPPPPPPPUUPPPPPPUUUUPPPPPPUUUUUUUUPPPPPPPPPPPPPPPPPPPPPPUUUPPPPPPPPPPUUUUUUPPPPPPUUUUPPPPPPUUUUUUUUUUUUUUUUUUUUUUUUUUUUUPPPPPPPPPUPPPPPPPUUUUPPPPPPPPPPPPUUUPPPPPPPPUPPPUPPPUUUPPPPUUUUUPPPPPPPPUPPPPPPPPPPUPPPUUUPPPUUPPPPPUUUUUUUUUUPPPPPUUPPPUPPPUPPPPPPPPPPUPPPPPPPUUPPPUU"
             ]
max_length = 400  # 假设最大长度为6

# 1. 将序列转换为one-hot编码
one_hot_encoded = convert_up_one_hot(sequences, max_length=max_length)
print("One-hot编码结果：")
print(one_hot_encoded)

# 2. 解码回原始序列
decoded_sequences = [decode_up(seq) for seq in one_hot_encoded]
print("\n解码后的序列：")
print(decoded_sequences)

"""
# str_2
def convert_one_hot_str_2(sequence, max_length):
    # 假设sequence是一个单一的字符串
    seq = sequence.upper()  # 转换为大写
    seq_length = len(seq)

    if seq_length > max_length:
        seq = seq[:max_length]
        seq_length = max_length

    # 创建一个2行seq_length列的零矩阵
    one_hot = np.zeros((2, seq_length))  # 2个通道，分别对应 U 和 P

    # 为每个字符分配one-hot编码
    index_u = [j for j in range(seq_length) if seq[j] == 'U']
    one_hot[0, index_u] = 1  # 'U' 对应 [1, 0]
    
    index_p = [j for j in range(seq_length) if seq[j] == 'P']
    one_hot[1, index_p] = 1  # 'P' 对应 [0, 1]

    # 如果给定了max_length，进行零填充
    if max_length:
        offset1 = int((max_length - seq_length) / 2)
        offset2 = max_length - seq_length - offset1
        if offset1:
            one_hot = np.hstack([np.zeros((2, offset1)), one_hot])  # 左边填充
        if offset2:
            one_hot = np.hstack([one_hot, np.zeros((2, offset2))])  # 右边填充

    return one_hot  # 返回一个单一的(2, max_length)矩阵
# str_2
def decode_str_2(m):
    seq = ""
    for i in range(m.shape[1]):  # 遍历每列
        # 找到值为1的通道
        if m[0, i] == 1:
            seq += 'U'
        elif m[1, i] == 1:
            seq += 'P'
    return seq


"""
# 输入序列
sequences = ["PLUMPULMMULP" , "PPP"]
max_length = 14  # 假设最大长度为6

# 1. 将序列转换为one-hot编码
one_hot_encoded = convert_plum_one_hot(sequences, max_length=max_length)
print("One-hot编码结果：")
print(one_hot_encoded)

# 2. 解码回原始序列
decoded_sequences = [decode_plum(seq) for seq in one_hot_encoded]
print("\n解码后的序列：")
print(decoded_sequences)
"""
# str_4
def convert_one_hot_str_4(sequence, max_length):
    # 字母到one-hot编码的映射
    letter_to_index = {'P': 0, 'L': 1, 'U': 2, 'M': 3}
    
    # 假设sequence是一个单一的字符串
    seq = sequence.upper()  # 转换为大写
    seq_length = len(seq)

    if seq_length > max_length:
        seq = seq[:max_length]
        seq_length = max_length
    
    one_hot = np.zeros((4, seq_length))  # 4个通道，分别对应 P, L, U, M
    
    # 为每个字母分配one-hot编码
    for i, letter in enumerate(seq):
        if letter in letter_to_index:
            one_hot[letter_to_index[letter], i] = 1  # 选择对应的通道并赋值为1

    # 如果给定了max_length，进行零填充
    if max_length:
        offset1 = int((max_length - seq_length) / 2)
        offset2 = max_length - seq_length - offset1
        if offset1:
            one_hot = np.hstack([np.zeros((4, offset1)), one_hot])  # 左边填充
        if offset2:
            one_hot = np.hstack([one_hot, np.zeros((4, offset2))])  # 右边填充

    return one_hot  # 返回一个单一的(4, max_length)矩阵
# str_4
def decode_str_4(m):
    index_to_letter = {0: 'P', 1: 'L', 2: 'U', 3: 'M'}
    seq = ""
    for i in range(m.shape[1]):  # 遍历每列
        # 找到值为1的通道
        var = np.where(m[:, i] == 1)[0]
        if len(var) == 1:  # 确保找到唯一的匹配
            seq += index_to_letter[var[0]]  # 将对应的字母添加到序列中
    return seq

"""
# 输入序列
sequences = ["BHTLM", "RHTLM", "BEEE", "RRR"]
max_length = 6  # 假设最大长度为6

# 1. 将序列转换为one-hot编码
one_hot_encoded = convert_behmlrt_one_hot(sequences, max_length=max_length)
print("One-hot编码结果：")
print(one_hot_encoded)

# 2. 解码回原始序列
decoded_sequences = [decode_behmlrt(seq) for seq in one_hot_encoded]
print("\n解码后的序列：")
print(decoded_sequences)
"""
# str_7
def convert_one_hot_str_7(sequence, max_length):
    # 字母到one-hot编码的映射
    letter_to_index = {'B': 0, 'E': 1, 'H': 2, 'L': 3, 'M': 4, 'R': 5, 'T': 6}
    
    # 假设sequence是一个单一的字符串
    seq = sequence.upper()  # 转换为大写
    seq_length = len(seq)
    if seq_length > max_length:
        seq = seq[:max_length]
        seq_length = max_length
    # 创建一个7行seq_length列的零矩阵
    one_hot = np.zeros((7, seq_length))  # 7个通道，分别对应 B, E, H, L, M, R, T
    
    # 为每个字母分配one-hot编码
    for i, letter in enumerate(seq):
        if letter in letter_to_index:
            one_hot[letter_to_index[letter], i] = 1  # 选择对应的通道并赋值为1

    # 如果给定了max_length，进行零填充
    if max_length:
        offset1 = int((max_length - seq_length) / 2)
        offset2 = max_length - seq_length - offset1
        if offset1:
            one_hot = np.hstack([np.zeros((7, offset1)), one_hot])  # 左边填充
        if offset2:
            one_hot = np.hstack([one_hot, np.zeros((7, offset2))])  # 右边填充

    return one_hot  # 返回一个单一的(7, max_length)矩阵
# str_7
def decode_str_7(m):
    index_to_letter = {0: 'B', 1: 'E', 2: 'H', 3: 'L', 4: 'M', 5: 'R', 6: 'T'}
    seq = ""
    for i in range(m.shape[1]):  # 遍历每列
        # 找到值为1的通道
        var = np.where(m[:, i] == 1)[0]
        if len(var) == 1:  # 确保找到唯一的匹配
            seq += index_to_letter[var[0]]  # 将对应的字母添加到序列中
    return seq

# seq_str_8
def convert_one_hot_seq_str_8(sequence, max_length):
    # 字母到one-hot编码的映射
    letter_to_index = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7}
    
    # 假设sequence是一个单一的字符串
    seq = sequence.upper()  # 转换为大写
    seq_length = len(seq)
    if seq_length > max_length:
        seq = seq[:max_length]
        seq_length = max_length
    # 创建一个8行seq_length列的零矩阵
    one_hot = np.zeros((8, seq_length))  # 8个通道，分别对应 A, B, C, D, E, F, G, H
    
    # 为每个字母分配one-hot编码
    for i, letter in enumerate(seq):
        if letter in letter_to_index:
            one_hot[letter_to_index[letter], i] = 1  # 选择对应的通道并赋值为1

    # 如果给定了max_length，进行零填充
    if max_length:
        offset1 = int((max_length - seq_length) / 2)
        offset2 = max_length - seq_length - offset1
        if offset1:
            one_hot = np.hstack([np.zeros((8, offset1)), one_hot])  # 左边填充
        if offset2:
            one_hot = np.hstack([one_hot, np.zeros((8, offset2))])  # 右边填充

    return one_hot  # 返回一个单一的(8, max_length)矩阵
# seq_str_8
def decode_seq_str_8(m):
    index_to_letter = {0: 'A', 1: 'B', 2: 'C', 3: 'D', 4: 'E', 5: 'F', 6: 'G', 7: 'H'}
    seq = ""
    for i in range(m.shape[1]):  # 遍历每列
        # 找到值为1的通道
        var = np.where(m[:, i] == 1)[0]
        if len(var) == 1:  # 确保找到唯一的匹配
            seq += index_to_letter[var[0]]  # 将对应的字母添加到序列中
    return seq

"""
# 输入序列
sequences = ["ABCD", "EFGH", "IJKL", "MNOP"]
max_length = 8  # 假设最大长度为8

# 1. 将序列转换为one-hot编码
one_hot_encoded = convert_behmlrt_16_one_hot(sequences, max_length=max_length)
print("One-hot编码结果：")
print(one_hot_encoded)

# 2. 解码回原始序列
decoded_sequences = [decode_behmlrt_16(seq) for seq in one_hot_encoded]
print("\n解码后的序列：")
print(decoded_sequences)
"""
# seq_str_16
def convert_one_hot_seq_str_16(sequence, max_length):
    # 字母到one-hot编码的映射
    letter_to_index = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 
                       'I': 8, 'J': 9, 'K': 10, 'L': 11, 'M': 12, 'N': 13, 'O': 14, 'P': 15}
    
    # 处理单个序列
    seq = sequence.upper()  # 转换为大写
    seq_length = len(seq)

    if seq_length > max_length:
        seq = seq[:max_length]
        seq_length = max_length
    
    one_hot = np.zeros((16, seq_length))  # 创建一个 16 通道的矩阵
    
    # 为每个字母分配one-hot编码
    for i, letter in enumerate(seq):
        if letter in letter_to_index:
            one_hot[letter_to_index[letter], i] = 1  # 选择对应的通道并赋值为1

    # 如果给定了max_length，进行零填充
    if max_length:
        offset1 = int((max_length - seq_length) / 2)  # 左侧填充
        offset2 = max_length - seq_length - offset1  # 右侧填充
        if offset1:
            one_hot = np.hstack([np.zeros((16, offset1)), one_hot])  # 左边填充
        if offset2:
            one_hot = np.hstack([one_hot, np.zeros((16, offset2))])  # 右边填充

    return one_hot  # 返回一个(16, max_length)的矩阵
# seq_str_16
def decode_seq_str_16(m):
    index_to_letter = {0: 'A', 1: 'B', 2: 'C', 3: 'D', 4: 'E', 5: 'F', 6: 'G', 7: 'H',
                       8: 'I', 9: 'J', 10: 'K', 11: 'L', 12: 'M', 13: 'N', 14: 'O', 15: 'P'}
    seq = ""
    for i in range(m.shape[1]):  # 遍历每列
        # 找到值为1的通道
        var = np.where(m[:, i] == 1)[0]
        if len(var) == 1:  # 确保找到唯一的匹配
            seq += index_to_letter[var[0]]  # 将对应的字母添加到序列中
    return seq

# seq_str_28
def convert_one_hot_seq_str_28(sequence, max_length):
    one_hot_seq = []
    
    # 28字母到one-hot编码的映射，包括大小写字母
    letter_to_index = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7, 
                       'I': 8, 'J': 9, 'K': 10, 'L': 11, 'M': 12, 'N': 13, 'O': 14, 'P': 15,
                       'Q': 16, 'R': 17, 'S': 18, 'T': 19, 'U': 20, 'V': 21, 'W': 22, 'X': 23,
                       'Y': 24, 'Z': 25, 'a': 26, 'b': 27}
    
    # 转换为大写字母并过滤掉非字母字符
    seq = ''.join([char for char in sequence if char.isalpha()])
    seq_length = len(seq)

    if seq_length > max_length:
        seq = seq[:max_length]
        seq_length = max_length
    
    one_hot = np.zeros((28, seq_length))  # 创建28个通道的矩阵
    
    # 为每个字母分配one-hot编码
    for i, letter in enumerate(seq):
        if letter in letter_to_index:
            one_hot[letter_to_index[letter], i] = 1  # 将对应位置设为1

    # 如果给定了max_length，进行零填充
    if max_length:
        offset1 = int((max_length - seq_length) / 2)
        offset2 = max_length - seq_length - offset1
        if offset1:
            one_hot = np.hstack([np.zeros((28, offset1)), one_hot])  # 左边填充
        if offset2:
            one_hot = np.hstack([one_hot, np.zeros((28, offset2))])  # 右边填充

    # 返回一个numpy数组，包含所有序列的one-hot编码
    return one_hot


def decode_seq_str_28(m):
    # 反向映射字典
    index_to_letter = {0: 'A', 1: 'B', 2: 'C', 3: 'D', 4: 'E', 5: 'F', 6: 'G', 7: 'H', 
                       8: 'I', 9: 'J', 10: 'K', 11: 'L', 12: 'M', 13: 'N', 14: 'O', 15: 'P',
                       16: 'Q', 17: 'R', 18: 'S', 19: 'T', 20: 'U', 21: 'V', 22: 'W', 23: 'X',
                       24: 'Y', 25: 'Z', 26: 'a', 27: 'b'}
    
    seq = ""
    for i in range(m.shape[1]):  # 遍历每列
        # 查找值为1的通道
        var = np.argmax(m[:, i])  # 找到最大值的索引，避免np.where的性能问题
        seq += index_to_letter[var]  # 将对应的字母添加到序列中

    return seq


def combine_one_hot_matrix(matrices, axis=0):
    combined_matrix = np.concatenate(matrices, axis=axis)
    # print(f"Combined shape: {combined_matrix.shape}")
    return combined_matrix