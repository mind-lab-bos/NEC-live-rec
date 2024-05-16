import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

def plotAllCondition():
    mat_contents = sio.loadmat('allQs.mat')

    allQ_data = np.zeros([21,4,9])
    allQ_data[:,:,0:8] = mat_contents['aggregate']
    allQ_data[:,:,-1] = np.mean(allQ_data[:,:,0:8],2)
    q_names = mat_contents['questions']
    # print(allQ_data.shape)
    for i in range(9):
        qname = q_names[0,i][0]
        single_Q_data = np.squeeze(allQ_data[:,:,i])
        plt.subplot(3,3,i+1)
        plt.violinplot(single_Q_data,showextrema=False,showmeans=True)
        plt.title(qname)
        plt.yticks([1, 2, 3, 4, 5])
        plt.xticks(np.arange(1,5),['LF', 'LS', 'RF', 'RS'])
    # single_Q_data = np.squeeze(allQ_data[:,:,1])
    # plt.subplot(3,3,2)
    # plt.violinplot(single_Q_data,showextrema=False,showmedians=True)
    plt.show()
def plotLiveRec():
    mat_contents = sio.loadmat('allQs.mat')

    allQ_data = np.zeros([21,4,9])
    allQ_data[:,:,0:8] = mat_contents['aggregate']
    allQ_data[:,:,-1] = np.mean(allQ_data[:,:,0:8],2)
    q_names = mat_contents['questions']
    # print(allQ_data.shape)
    for i in range(9):
        qname = q_names[0,i][0]
        Live_Q_data = np.reshape(np.squeeze(allQ_data[:,0:2,i]),42)
        Rec_Q_data = np.reshape(np.squeeze(allQ_data[:,2:4,i]),42)
        print(Live_Q_data)
        plotData = np.vstack((Live_Q_data.T,Rec_Q_data.T)).T
        print(plotData.shape)
        plt.subplot(3,3,i+1)
        plt.violinplot(plotData,showextrema=False,showmeans=True)
        plt.title(qname)
        plt.yticks([1, 2, 3, 4, 5])
        plt.xticks(np.arange(1,3),['Live', 'Rec'])
    # # single_Q_data = np.squeeze(allQ_data[:,:,1])
    # # plt.subplot(3,3,2)
    # # plt.violinplot(single_Q_data,showextrema=False,showmedians=True)
    plt.show()

def main():
    plotAllCondition()
    # plotLiveRec()



if __name__ == "__main__":
    main()