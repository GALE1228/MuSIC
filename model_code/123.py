from PIL import Image
import numpy as np
import os

# 假设你已经加载了 chars
chars = np.load("/data3/MuSIC2_new_model/MotifLearnModel/MuSIC/model_code/acgu.npz", allow_pickle=True)['data']

# 输出文件夹
os.makedirs("acgu_chars_output", exist_ok=True)

# 逐张保存图片
for i, img_array in enumerate(chars):
    img = Image.fromarray(img_array.astype(np.uint8))
    img.save(f"acgu_chars_output/char_{i}.png")