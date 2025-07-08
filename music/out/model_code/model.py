import torch
import torch.nn as nn
import torch.nn.functional as F
from .resnet import ResidualBlock1D, ResidualBlock2D
from .se import SEBlock
from .Conv2d import Conv2d

class CNN(nn.Module):
    def __init__(self):
        super(CNN, self).__init__()
        h_k = 0
        self.n_features = 6
        h_p, h_k = 1, 3
        base_channel = 16
        self.conv    = Conv2d(in_channels = 1, out_channels = base_channel, kernel_size=(h_k, h_k), bn = True, same_padding=True)
        self.se      = SEBlock(base_channel)
        self.res2d   = ResidualBlock2D(base_channel, kernel_size=(h_k, h_k), padding=(h_p,h_p)) 
        self.res1d   = ResidualBlock1D(base_channel*4)
        self.avgpool = nn.AvgPool2d((self.n_features,1))
        self.gpool   = nn.AdaptiveAvgPool1d(1)
        self.fc      = nn.Linear(base_channel*4*8, 1)
        self._initialize_weights()

    def _initialize_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.Conv1d):
                nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.BatchNorm2d):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.BatchNorm1d):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.Linear):
                nn.init.normal_(m.weight, 0, 0.01)
                nn.init.constant_(m.bias, 0)
    
    def forward(self, input):
        """[forward]
        
        Args:
            input ([tensor],N,C,W,H): input features
        """
        # print("模型输入原始数据大小 ",input.shape)
        x = self.conv(input)
        # print("二维卷积层输出大小 ",x.shape)
        x = F.dropout(x, 0.1, training=self.training)
        z = self.se(x)
        # print("注意力输出大小 ",z.shape)
        # print("二维残差块输入大小 ",(x*z).shape)
        x = self.res2d(x*z)
        # print("二维残差块输出大小 ",x.shape)
        # x = F.dropout(x, 0.5, training=self.training)
        x = F.dropout(x, 0.5, training=self.training)
        x = self.avgpool(x)
        # print("平均池化层输出大小 ",x.shape)
        x = x.view(x.shape[0], x.shape[1], x.shape[3])
        # print("reshape结果 ",x.shape)
        x = self.res1d(x)
        # print("一维残差块输出大小 ",x.shape)
        x = F.dropout(x, 0.3, training=self.training)
        x = self.gpool(x)
        # print("自适应平均池化层输出 ",x.shape)
        x = x.view(x.shape[0], x.shape[1])
        # print("reshape结果 ",x.shape)
        x = self.fc(x)
        # print(x.shape)
        return x





class CNN_SHAP(nn.Module):
    def __init__(self):
        super(CNN_SHAP, self).__init__()
        h_k = 0
        self.n_features = 6
        h_p, h_k = 1, 3
        base_channel = 16
        self.conv    = Conv2d(in_channels = 1, out_channels = base_channel, kernel_size=(h_k, h_k), bn = True, same_padding=True)
        self.se      = SEBlock(base_channel)
        self.res2d   = ResidualBlock2D(base_channel, kernel_size=(h_k, h_k), padding=(h_p,h_p)) 
        self.res1d   = ResidualBlock1D(base_channel*4)
        self.avgpool = nn.AvgPool2d((self.n_features,1))
        self.gpool   = nn.AdaptiveAvgPool1d(1)
        self.fc      = nn.Linear(base_channel*4*8, 1)
        self._initialize_weights()

    def _initialize_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.Conv1d):
                nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.BatchNorm2d):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.BatchNorm1d):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.Linear):
                nn.init.normal_(m.weight, 0, 0.01)
                nn.init.constant_(m.bias, 0)
    
    def forward(self, input):
        """[forward]
        
        Args:
            input ([tensor],N,C,W,H): input features
        """
        # print("模型输入原始数据大小 ",input.shape)
        x = self.conv(input)
        # print("二维卷积层输出大小 ",x.shape)
        x = F.dropout(x, 0.1, training=self.training)
        z = self.se(x)
        # print("注意力输出大小 ",z.shape)
        # print("二维残差块输入大小 ",(x*z).shape)
        x = self.res2d(x*z)
        # print("二维残差块输出大小 ",x.shape)
        # x = F.dropout(x, 0.5, training=self.training)
        x = F.dropout(x, 0.5, training=self.training)
        x = self.avgpool(x)
        # print("平均池化层输出大小 ",x.shape)
        x = x.view(x.shape[0], x.shape[1], x.shape[3])
        # print("reshape结果 ",x.shape)
        x = self.res1d(x)
        # print("一维残差块输出大小 ",x.shape)
        x = F.dropout(x, 0.3, training=self.training)
        x = self.gpool(x)
        print("自适应平均池化层输出 ",x.shape)
        # fcl = self.gpool(x)
        x = x.view(x.shape[0], x.shape[1])

        print(x.shape)
        x = self.fc(x)
        # print(x.shape)
        return x
