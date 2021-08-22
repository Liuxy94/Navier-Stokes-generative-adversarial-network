import torch,copy,argparse,csv
from torch.autograd import Variable
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import numpy as np
import pandas as pd
import torch.nn as nn
import torch.nn.functional as F

# 固定随机种子
torch.manual_seed(15)
torch.cuda.manual_seed_all(15)
np.random.seed(15)
torch.backends.cudnn.deterministic = True





def main():
    ## ======================== 数据加载 ======================
    table = pd.read_csv("POD_coeffs_3900_new_grid_221_42.csv", header=None, dtype = float).values.astype(np.float16)
    data = copy.deepcopy(table[:,:-1].T)
    label = copy.deepcopy(table[:,1:].T)

    ## ======================== 数据集划分 ======================
    x_train, x_test, y_train, y_test = train_test_split(data,label, test_size=0.3 , random_state=2) # 70%用于训练，30%用于验证（测试）
    x_train = torch.from_numpy(x_train).type(torch.FloatTensor)
    y_train = torch.from_numpy(y_train).type(torch.FloatTensor)
    x_test = torch.from_numpy(x_test).type(torch.FloatTensor)
    y_test = torch.from_numpy(y_test).type(torch.FloatTensor)


    ## ======================== 搭建MLP网络 ======================
    # 网络G，用于预测下一时刻的特征向量
    G = torch.nn.Sequential(
        torch.nn.Linear(10+10,30),   # 特征向量+噪声向量
        torch.nn.ReLU(),   
        torch.nn.Linear(30,15), 
        torch.nn.ReLU(),   
        torch.nn.Linear(15,10),   
    )
    # 网络D，用于一组特征向量是网络G制造的还是真实标签，促使G的输出接近于真实标签
    D = torch.nn.Sequential(
        torch.nn.Linear(10,6),   
        torch.nn.ReLU(),   
        torch.nn.Linear(6,3), 
        torch.nn.ReLU(),   
        torch.nn.Linear(3,1),
        nn.Sigmoid()   
    )


    loss_train = []
    loss_val = []

    ## ========================模型验证======================
    def valid(model,data,criterion):
        model.eval()
        z = torch.from_numpy(np.random.randn(data.shape[0], 10)).float() # 随机噪声
        with torch.no_grad():   # 不记录梯度信息
            y_pred = model(torch.cat([z, data], dim=1))
            loss = criterion(y_pred, y_test)
            loss_val.append(loss.item())

    ## ======================== 训练模型 ======================
    loss_func_bce = nn.BCELoss()    # 用于鉴别器模型D（对抗损失），D的任务判断是否为标签，属于二分类任务，因此用交叉熵损失
    loss_func_reg = nn.MSELoss()    # 用于生成器G，计算模型预测和标签的差异，因为是回归任务所以用均方根损失
    opt_g = torch.optim.Adam(G.parameters(), lr=3e-4)   
    opt_d = torch.optim.Adam(D.parameters(), lr=1e-4)
    epochs = 1000
    for epoch in range(epochs):
        # 训练判别器D，目的是能够区分G的输出和真实标签
        for d in range(1):
            D.train()
            G.eval()
            # 前向传播
            z=torch.from_numpy(np.random.randn(x_train.shape[0], 10)).float()   # 随机噪声
            d_real = D(y_train)     #  （标签）输入判别器的结果
            y_gen = G(torch.cat([z, x_train], dim=1))      #  （噪声+x）输入生成器，得到预测结果
            d_gen = D(y_gen)        #  （预测结果）输入判别器的结果
            # 计算损失
            Dloss_real = loss_func_bce(d_real, torch.ones((x_train.shape[0],1))) # 对于（标签+x）组，判别器输出应趋向于全1
            Dloss_gen = loss_func_bce(d_gen, torch.zeros((x_train.shape[0],1)))  # 对于（预测结果+x）组，判别器输出应趋向于全0
            Dloss = Dloss_real + Dloss_gen
            # 反向传播（只对判别器参数进行更新）
            Dloss.backward()
            opt_d.step()
            opt_d.zero_grad()
            opt_g.zero_grad()
        # 训练生成器G，目的是G的输出能够欺骗D，让D以为G的输出就是真实标签
        for g in range(3):
            D.eval()
            G.train()
            # 前向传播
            z = torch.from_numpy(np.random.randn(x_train.shape[0], 10)).float() # 随机噪声
            y_gen = G(torch.cat([z, x_train], dim=1))
            d_gen = D(y_gen)  
            # 计算损失函数
            Gloss_adventure = 0.3 * loss_func_bce(d_gen, torch.ones((x_train.shape[0],1)))   # G的目的是，让D以为它的输出就是真实标签，因此G趋向于让d_gen等于1
            Gloss_regression = loss_func_reg(y_gen,y_train)
            Gloss = Gloss_regression + Gloss_adventure
            # 反向传播
            Gloss.backward()
            opt_g.step()
            opt_g.zero_grad()
            opt_d.zero_grad()
            loss_train.append(Gloss_regression.item())
            valid(G,x_test,loss_func_reg)
        D.eval()
        G.eval()
    


    ## ========================可视化训练过程======================
    x = [i for i in range(epochs*3)]    # 每一轮G都训练了3次，所以乘3
    plt.plot(x,loss_train, label='train')
    plt.plot(x,loss_val, label='val')
    plt.title('GAN_with_noise')
    plt.legend()
    plt.savefig(fname="result/GAN_with_noise.png")
    np.save('result/GAN_with_noise.npy',loss_train)   # 保存为npy文件,供不同方法对比
    print('GAN_with_noise——测试集均方误差：',loss_val[-1])
    plt.show()
    
if __name__=='__main__':
    main()