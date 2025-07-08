import torch
import torch.optim as optim
import torch.nn as nn
from .metrics_utils import MLMetrics
import numpy as np
import argparse, os, copy
from tqdm import tqdm
import matplotlib.pyplot as plt



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

def mat2str(m):
    string=""
    if len(m.shape)==1:
        for j in range(m.shape[0]):
            string+= "%.3f," % m[j]
    else:
        for i in range(m.shape[0]):
            for j in range(m.shape[1]):
                string+= "%.3f," % m[i,j]
    return string

def train(model, device, train_loader, criterion, optimizer, batch_size, smooth_rate):
    model.train()
    met = MLMetrics(objective='binary')
    previous_loss_cls = 0.0
    previous_loss_smooth = 0.0
    for batch_idx, (x0, y0, y_s0) in enumerate(train_loader):
        x, y, y_s = x0.float().to(device), y0.to(device).float(), y_s0.to(device).float()

        if y0.sum() == 0 or y0.sum() == batch_size:
            continue

        optimizer.zero_grad()
        output = model(x)
        output = output.squeeze(1)

        loss_cls = criterion(output, y)
        loss_smooth = criterion(output, y_s)

        loss_cls.backward(retain_graph=True)
        loss_smooth.backward()

        gradients_loss_cls = [p.grad for p in model.parameters()]
        gradients_loss_smooth = [p.grad for p in model.parameters()]

        delta_smooth = loss_smooth.item() / previous_loss_smooth if previous_loss_smooth > 0 else 1
        delta_cls = loss_cls.item() / previous_loss_cls if previous_loss_cls > 0 else 1

        weight_cross_entropy = 1 / delta_cls
        weight_smooth = 1 / delta_smooth

        weight_cross_entropy_norm = weight_cross_entropy / (weight_cross_entropy + weight_smooth)
        weight_smooth_norm = weight_smooth / (weight_cross_entropy + weight_smooth)

        final_gradients = []
        for grad_lc, grad_smooth in zip(gradients_loss_cls, gradients_loss_smooth):
            final_gradients.append(weight_cross_entropy_norm * grad_lc + weight_smooth_norm * grad_smooth)

        prob = torch.sigmoid(output)

        y_np = y.to(device='cpu', dtype=torch.long).detach().numpy()
        p_np = prob.to(device='cpu').detach().numpy()

        loss_total = smooth_rate * loss_cls + (1 - smooth_rate) * loss_smooth

        met.update(y_np, p_np, [loss_total.item()])
        
        # loss_total.backward(retain_graph=True)
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5)

        for param, final_grad in zip(model.parameters(), final_gradients):
            param.grad = final_grad

        optimizer.step()

        previous_loss_cls = loss_cls.item()
        previous_loss_smooth = loss_smooth.item()

    return met

def train_rbp_smi(model, device, train_loader, criterion, optimizer ,batch_size, rbp_smi):
    model.train()
    met = MLMetrics(objective='binary')
    for batch_idx, (x0, y0) in enumerate(train_loader):
        x, y = x0.float().to(device), y0.to(device).float()
        if y0.sum() ==0 or y0.sum() ==batch_size:
            continue
        optimizer.zero_grad()
        output = model(x)
        output = output.squeeze(1)

        output_max = torch.max(output)
        output_min = torch.min(output)
    
        # print(f"train Output max: {output_max.item()}, train Output min: {output_min.item()}")
        loss = rbp_smi * criterion(output, y)
        prob = torch.sigmoid(output)

        y_np = y.to(device='cpu', dtype=torch.long).detach().numpy()
        p_np = prob.to(device='cpu').detach().numpy()
        met.update(y_np, p_np,[loss.item()])
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5)
        optimizer.step()

    return met

def validate(args, model, device, test_loader, criterion ,smooth_rate):
    model.eval()
    y_all = []
    p_all = []
    l_all = []
    with torch.no_grad():
        for batch_idx, (x0, y0, y_s0) in enumerate(test_loader):
            x, y, y_s = x0.float().to(device), y0.to(device).float(), y_s0.to(device).float()
            output  = model(x)
            output = output.squeeze(1)

            loss_cls = criterion(output, y)
            loss_smooth = criterion(output, y_s)
            prob = torch.sigmoid(output)
            loss_total = smooth_rate * loss_cls + (1 - smooth_rate) * loss_smooth
            y_np = y.to(device='cpu', dtype=torch.long).numpy()
            p_np = prob.to(device='cpu').numpy()
            l_np = loss_total.item()

            y_all.append(y_np)
            p_all.append(p_np)
            l_all.append(l_np)

    y_all = np.concatenate(y_all)
    p_all = np.concatenate(p_all)
    l_all = np.array(l_all)
    
    met = MLMetrics(objective='binary')
    met.update(y_all, p_all,[l_all.mean()])
    

    
    return met, y_all, p_all

def inference(args, model, device, test_loader, rna_names_all):
    model.eval()
    p_all = []
    y_all = []
    rna_names_out = []  # 保存 RNA 名称
    with torch.no_grad():
        for batch_idx, (x0, y0) in enumerate(test_loader):
            x, y = x0.float().to(device), y0.to(device).float()
            output = model(x)
            prob = torch.sigmoid(output)

            p_np = prob.to(device='cpu').numpy()
            p_all.append(p_np)

            y_np = y.to(device='cpu').numpy()
            y_all.append(y_np)

            # 获取 RNA 名称并追加
            rna_names_out.extend(rna_names_all[batch_idx * test_loader.batch_size : (batch_idx + 1) * test_loader.batch_size])

    p_all = np.concatenate(p_all)
    y_all = np.concatenate(y_all)
    return p_all, y_all, rna_names_out

def compute_saliency(args, out_dir, model, device, test_loader, identity):
    from model_code.smoothgrad import GuidedBackpropSmoothGrad
    model.eval()

    saliency_dir = make_directory(out_dir, "out/saliency")
    saliency_path = os.path.join(saliency_dir, identity+'.sal')

    # sgrad = SmoothGrad(model, device=device)
    sgrad = GuidedBackpropSmoothGrad(model, device=device)
    sal = ""
    for batch_idx, (x0, y0) in enumerate(test_loader):
        X, Y = x0.float().to(device), y0.to(device).float()
        output = model(X)
        prob = torch.sigmoid(output)
        p_np = prob.to(device='cpu').detach().numpy().squeeze(-1)
        guided_saliency = sgrad.get_batch_gradients(X, Y)

        # print(f"Shape of guided_saliency: {guided_saliency.shape}") (N, 6, 200)
        N, _, NS = guided_saliency.shape # (N, 6, 200)
        
        for i in range(N):
            inr = batch_idx*args.batch_size + i
            str_sal = mat2str(np.squeeze(guided_saliency[i]))
            sal += "{}\t{:.6f}\t{}\n".format(inr, p_np[i], str_sal)
            
    f = open(saliency_path,"w")
    f.write(sal)
    f.close()
    print(saliency_path)

def compute_saliency_img(args, out_dir, model, device, test_loader, identity , rna_names_all):
    from model_code.smoothgrad import GuidedBackpropSmoothGrad
    from model_code import visualize


    def saliency_img(X, mul_saliency, outdir="results"):
        """generate saliency image

        Args:
            X ([np.ndarray]): raw input(L x 6)
            mul_saliency ([np.ndarray]): [description]
            outdir (str, optional): [description]. Defaults to "results".
        """
        X =X.T
        mul_saliency =mul_saliency.T

        s_str = mul_saliency[:, 4:]
        s_str = (s_str - s_str.min()) / (s_str.max() - s_str.min())  # Normalize structure saliency

        # Replace the original structure saliency with the normalized one
        mul_saliency[:, 4:] = s_str

        # Visualize the saliency image
        visualize.plot_saliency(
            X.T, 
            mul_saliency.T, 
            nt_width=100, 
            norm_factor=5, 
            str_null=None,  # No longer needed
            outdir=outdir
        )


    prefix_n = len(str(len(test_loader.dataset)))
    make_directory(out_dir, "out/saliency_imgs/")
    imgs_dir = make_directory(out_dir, "out/saliency_imgs/"+identity)
    imgs_path = imgs_dir + '/{:0' + str(prefix_n) + 'd}_{}_ {:.3f}.pdf'
    saliency_path = os.path.join(imgs_dir, 'all.saliency')

    # sgrad = SmoothGrad(model, device=device)
    sgrad = GuidedBackpropSmoothGrad(model, device=device, magnitude=1)
    for batch_idx, (x0, y0) in enumerate(test_loader):
        X, Y = x0.float().to(device), y0.to(device).float()
        output = model(X)
        prob = torch.sigmoid(output)
        p_np = prob.to(device='cpu').detach().numpy().flatten()
        guided_saliency  = sgrad.get_batch_gradients(X, Y)
        mul_saliency = copy.deepcopy(guided_saliency)
        # print("guided_saliency.shape",guided_saliency.shape)
        mul_saliency[:,:,:4] = guided_saliency[:,:,:4] * X[:,:,:4]
        N, _, NS = guided_saliency.shape # (N, 6, 200)
        sal = ""
        for i in tqdm(range(N)):
            inr = batch_idx*args.batch_size + i
            rna_name = rna_names_all[inr]
            str_sal = mat2str(np.squeeze(guided_saliency[i]))
            sal += "{}\t{:.6f}\t{}\n".format(inr, p_np[i], str_sal)
            img_path = imgs_path.format(inr, rna_name, p_np[i])
            # import pdb; pdb.set_trace()
            saliency_img(
                X[i].to(device='cpu').detach().numpy(),  # 直接使用 1D 索引
                mul_saliency[i].to(device='cpu').numpy(), 
                outdir=img_path)
    if not os.path.exists(saliency_path):     
        f = open(saliency_path,"w")
        f.write(sal)
        f.close()
        print(saliency_path)


def compute_high_attention_region(args, out_dir, model, device, test_loader, identity):
    from model_code.smoothgrad import GuidedBackpropSmoothGrad
    model.eval()
    har_dir = make_directory(out_dir, "out/har")
    har_path = os.path.join(har_dir, identity+'.har')

    L = 20
    har = ""
    # sgrad = SmoothGrad(model, device=device)
    sgrad = GuidedBackpropSmoothGrad(model, device=device)
    for batch_idx, (x0, y0) in enumerate(test_loader):
        X, Y = x0.float().to(device), y0.to(device).float()
        output = model(X)
        prob = torch.sigmoid(output)
        p_np = prob.to(device='cpu').detach().numpy().squeeze()
        guided_saliency  = sgrad.get_batch_gradients(X, Y)

        attention_region = guided_saliency.sum(dim=1, keepdim=True)[:, 0, :].to(device='cpu').numpy() # (N, 1, 200)
        N,NS = attention_region.shape # (N, 101)
        for i in range(N):
            inr = batch_idx*args.batch_size + i
            iar = attention_region[i]
            ar_score = np.array([ iar[j:j+L].sum() for j in range(NS-L+1)])
            # import pdb; pdb.set_trace()
            highest_ind = np.argmax(iar)
            har += "{}\t{:.6f}\t{}\t{}\n".format(inr, p_np[i], highest_ind, highest_ind+L)

    f = open(har_path,"w")
    f.write(har)
    f.close()
    print(har_path)


def process_shap_labels(rna_names_all , num ,shap_values): 
    rna_names_all = [name.decode().lstrip('> ') if isinstance(name, bytes) else name.lstrip('> ') for name in rna_names_all]
    labels = rna_names_all[:num]
    labels= np.tile(labels, (len(shap_values), 1)).T
    return labels

def save_shap_fig(out_dir, identity , dpi):
    shap_dir = make_directory(out_dir, "out/shap")
    shap_path = os.path.join(shap_dir, identity+'.pdf')
    plt.savefig(shap_path, bbox_inches='tight', dpi=dpi)
    plt.close()  # 关闭当前的绘图窗口
    print(shap_path," saved")
