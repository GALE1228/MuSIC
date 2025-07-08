from utils import *
import numpy as np
import argparse
import torch
import torch.nn as nn
from model_code.model import CNN
from train_code.train_loop import train ,validate ,inference ,compute_saliency ,compute_saliency_img ,compute_high_attention_region ,process_shap_labels,save_shap_fig
from model_code.GradualWarmupScheduler import GradualWarmupScheduler
import logging
import os
import shap
import matplotlib.pyplot as plt


def main(args):

    file_path = args.file_path
    out_dir = args.out_dir
    rbp = args.rbp_name
    batch_size = args.batch_size
    gpuid = args.gpuid
    nepochs = args.num_epochs
    weight_decay = args.weight_decay
    pos_weight = args.pos_weight
    early_stopping = args.early_stopping
    exp_name = args.exp_name
    learn_rate = args.learn_rate
    
    species_name = args.species_name


    device = torch.device(f"cuda:{gpuid}" if torch.cuda.is_available() else "cpu")
    print("device gpu ID is", gpuid)

    if args.cross:
        identity = f"{rbp}_{exp_name}_cross_{species_name}"
        smooth_rate = args.smooth_rate
    else:
        identity = f"{rbp}_{exp_name}_within"
        smooth_rate = 1

    print("rbp_clip information is", identity)

    best_model_path = f"{out_dir}/out/model/{identity}_best.pth"
    os.makedirs(os.path.dirname(best_model_path), exist_ok=True)
    

    if args.train:

        """
        测试平滑label的训练代码

        """
        train_loader , test_loader = train_dataset(file_path , rbp , batch_size ,smooth_rate)

        model = CNN().to(device)

        optimizer = torch.optim.Adam(model.parameters(), lr=learn_rate, betas=(0.9, 0.999), weight_decay=weight_decay)
        scheduler = GradualWarmupScheduler(
            optimizer, multiplier=8, total_epoch=float(nepochs), after_scheduler=None)
        criterion = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(pos_weight))

        best_auc = 0
        best_acc = 0
        best_epoch = 0

        log_dir = f"{out_dir}/out/logs"
        print("log dir is ", log_dir)
        os.makedirs(log_dir, exist_ok=True)
        log_file = os.path.join(log_dir, f"{identity}.txt")
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file , mode='w'),
                logging.StreamHandler()
            ]
        )

        # 定义颜色的 ANSI 转义码
        COLOR_GREEN = '\033[92m'
        COLOR_RED = '\033[91m'
        COLOR_RESET = '\033[0m'

        for epoch in range(1, nepochs + 1):
            t_met = train(model, device, train_loader, criterion, optimizer, batch_size=64 , smooth_rate=smooth_rate)
            v_met, _, _ = validate(args, model, device, test_loader, criterion , smooth_rate=smooth_rate)
            scheduler.step()
            lr = scheduler.get_lr()[0]
            
            color_best = 'green'  # 默认颜色为绿色
            
            # 如果当前的 AUC 比最好的 AUC 高，更新最佳模型
            if best_auc < v_met.auc:
                best_auc = v_met.auc
                best_acc = v_met.acc
                best_epoch = epoch
                color_best = 'red'  # 如果是最好的 AUC，使用红色
                # 确保传递正确的文件路径
                torch.save(model.state_dict(), best_model_path)

            # 检查是否满足早停条件
            if epoch - best_epoch > early_stopping:
                print(f"Early stop at {epoch}, {exp_name}")
                break

            # 打印日志信息，并添加颜色
            if color_best == 'red':
                color = COLOR_RED
            else:
                color = COLOR_GREEN

            logging.info(f'{color}Train Epoch: {epoch} avg.loss: {t_met.other[0]:.4f} '
                        f'Acc: {t_met.acc:.2f}, AUC: {t_met.auc:.4f}, lr: {lr:.6f}{COLOR_RESET}')
            
            logging.info(f'{color}Test Epoch: {epoch} avg.loss: {v_met.other[0]:.4f} '
                        f'Acc: {v_met.acc:.2f}, AUC: {v_met.auc:.4f} '
                        f'({best_auc:.4f} best){COLOR_RESET}')
            
        # 打印最终结果
        logging.info("%s auc: %.4f acc: %.4f", "TEST", best_auc, best_acc)

        filename = best_model_path.format("best")
        print("Loading model: {}".format(filename))
        model.load_state_dict(torch.load(filename))

    if args.validate:
        data_loader = validation_dataset(file_path , rbp , batch_size ,smooth_rate)
        print("Test  set:", len(data_loader.dataset))

        best_model = CNN().to(device)
        best_model = load_model(best_model, best_model_path, device)
        print("load best model path is ", best_model_path)

        criterion = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(pos_weight))
        met, y_all, p_all = validate(args, best_model, device, data_loader, criterion , smooth_rate)

        p_name = identity

        # out_dir = f"./{exp_name}"

        print("> eval {} auc: {:.4f} acc: {:.4f}".format(p_name, met.auc, met.acc))

        save_validations(out_dir, identity, p_name, p_all, y_all, met)

    if args.gerenate_h5:
        fasta_path = args.infer_fasta_path
        print("Inference fasta file path :", fasta_path)
        data_loader, rna_names_all = inference_dataset(fasta_path , batch_size)
        print("Test  set:", len(data_loader.dataset))

    if args.infer:
        fasta_path = args.infer_fasta_path
        print("Inference fasta file path :", fasta_path)
        data_loader, rna_names_all = inference_dataset(fasta_path , batch_size)
        print("Test  set:", len(data_loader.dataset))

        best_model = CNN().to(device)
        best_model = load_model(best_model, best_model_path, device)
        print("load best model path is ", best_model_path)

        criterion = nn.BCEWithLogitsLoss(pos_weight=torch.tensor(pos_weight))
        p_all, y_all, rna_names_out = inference(args, best_model, device, data_loader , rna_names_all)

        identity = identity+"_"+ os.path.basename(fasta_path).replace(".fa","") 
        # out_dir = f"./{exp_name}"

        print(p_all.shape)
        save_infers(out_dir, identity, rna_names_out ,y_all, p_all)

    if args.saliency:

        fasta_path = args.infer_fasta_path
        print("Inference fasta file path :", fasta_path)
        data_loader, rna_names_all = inference_dataset(fasta_path , batch_size)
        print("Test  set:", len(data_loader.dataset))

        best_model = CNN().to(device)
        best_model = load_model(best_model, best_model_path, device)
        print("load best model path is ", best_model_path)


        identity = identity+"_"+ os.path.basename(fasta_path).replace(".fa","") 
        # out_dir = f"./{exp_name}"

        compute_saliency(args, out_dir, best_model, device, data_loader, identity)

    if args.saliency_img:
        fasta_path = args.infer_fasta_path
        print("Inference fasta file path :", fasta_path)
        data_loader, rna_names_all = inference_dataset(fasta_path , batch_size)
        print("Test  set:", len(data_loader.dataset))

        best_model = CNN().to(device)
        best_model = load_model(best_model, best_model_path, device)
        print("load best model path is ", best_model_path)


        identity = identity+"_"+ os.path.basename(fasta_path).replace(".fa","") 
        # out_dir = f"./{exp_name}"

        compute_saliency_img(args, out_dir, best_model, device, data_loader, identity, rna_names_all)
    
    if args.har:
        fasta_path = args.infer_fasta_path
        print("Inference fasta file path :", fasta_path)
        data_loader, rna_names_all = inference_dataset(fasta_path , batch_size)
        print("Test  set:", len(data_loader.dataset))

        best_model = CNN().to(device)
        best_model = load_model(best_model, best_model_path, device)
        print("load best model path is ", best_model_path)


        identity = identity+"_"+ os.path.basename(fasta_path).replace(".fa","") 
        # out_dir = f"./{exp_name}"

        compute_high_attention_region(args, out_dir, best_model, device, data_loader, identity)



    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process and concatenate one-hot encoded matrices.")
    
    parser.add_argument('--file_path', type=str)
    parser.add_argument('--infer_fasta_path', type=str)
    parser.add_argument('--rbp_name', type=str)
    parser.add_argument('--out_dir', type=str, default="music", help="Output directory for the results.")

    parser.add_argument('--gerenate_h5', action='store_true', help='RNAfold fasta file and generate h5 file')
    parser.add_argument('--train', action='store_true', help='train mode')
    parser.add_argument('--validate', action='store_true', help='validate mode')
    parser.add_argument('--infer', action='store_true', help='infer mode')
    parser.add_argument('--saliency', action='store_true', help='compute saliency mode')
    parser.add_argument('--saliency_img', action='store_true', help='compute saliency and plot image mode')
    parser.add_argument('--har', action='store_true', help='compute highest attention region')
    parser.add_argument('--cross', action='store_true', help='cross species application')


    parser.add_argument('--gpuid', type=int, default=0, help="Specify which GPU to use, e.g., 0 for the first GPU.")
    parser.add_argument('--smooth_rate', type=float, default=1, help="Rate for smoothing labels.")
    parser.add_argument('--batch_size', type=int, default=64, help="Batch size for training or testing.")
    parser.add_argument('--num_epochs', type=int, default=200, help="Number of epochs for training.")
    parser.add_argument('--weight_decay', type=float, default=1e-6, help="Weight decay (L2 regularization) factor.")
    parser.add_argument('--pos_weight', type=float, default=2, help="Weight factor for the positive class in imbalanced datasets.")
    parser.add_argument('--learn_rate', type=float, default=0.001, help="Learning rate for training.")
    parser.add_argument('--early_stopping', type=int, default=20, help="Number of epochs with no improvement before early stopping.")
    parser.add_argument('--exp_name', type=str, default="music", help="Name for the experiment, used for saving and logging.")
    parser.add_argument('--species_name', type=str, default="human", help="Name for the cross species, used for saving and logging.")
    
    args = parser.parse_args()
    main(args)


