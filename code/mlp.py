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




def cal_score(model,data,label):
    """ 计算准确率 """
    y_pred = torch.max(F.softmax(model(data),dim=1), 1)[1].numpy()
    accuracy = (y_pred==label).sum()/len(label)
    return accuracy



def main():
    ## ========================数据加载======================
    table = pd.read_csv("POD_coeffs_3900_new_grid_221_42.csv", header=None, dtype = float).values.astype(np.float16)
    data = copy.deepcopy(table[:,:-1].T)
    label = copy.deepcopy(table[:,1:].T)

    ## ========================数据集划分======================
    x_train, x_test, y_train, y_test = train_test_split(data,label, test_size=0.3 , random_state=2) # 70%用于训练，30%用于验证（测试）
    x_train = torch.from_numpy(x_train).type(torch.FloatTensor)     # 数据格式转换成pytorch要求的Tensor类型
    y_train = torch.from_numpy(y_train).type(torch.FloatTensor)
    x_test = torch.from_numpy(x_test).type(torch.FloatTensor)
    y_test = torch.from_numpy(y_test).type(torch.FloatTensor)


    ## ========================搭建MLP网络======================
    model=torch.nn.Sequential(
        torch.nn.Linear(10,30),   
        torch.nn.ReLU(),   
        torch.nn.Linear(30,15), 
        torch.nn.ReLU(),   
        torch.nn.Linear(15,10),   
    )

    
    loss_train = []
    loss_val = []

    ## ========================模型验证======================
    def valid(model,data,criterion):
        model.eval()
        with torch.no_grad():   # 不记录梯度信息
            y_pred = model(data)
            loss = criterion(y_pred, y_test)
            loss_val.append(loss.item())

    ## ========================训练模型======================
    loss_fn = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=3e-4)
    epochs = 3000
    for e in range(epochs):
        y_pred = model(x_train)
        loss = loss_fn(y_pred, y_train)
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        loss_train.append(loss.item())
        valid(model,x_test,loss_fn)

    ## ========================可视化训练过程======================
    x = [i for i in range(epochs)]
    plt.plot(x,loss_train, label='train')
    plt.plot(x,loss_val, label='val')
    plt.title('MLP')
    plt.legend()
    plt.savefig(fname="result/mlp.png")
    np.save('result/mlp.npy',loss_train)   # 保存为npy文件,供不同方法对比
    print('MLP——测试集均方误差：',loss_val[-1])
    plt.show()
    
if __name__=='__main__':
    main()