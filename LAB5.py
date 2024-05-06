# 判断某一个矩阵是实Schur矩阵的方法
import numpy as np

def Power(A):
    eigenvalues = []
    eigen_vector = []
    x_0 = np.ones(A.shape[0])
    while True:
        x_1 = np.dot(A, x_0)

        # eigenvalue = np.linalg.norm(x_1, np.inf)
        # 求 x_1 最大分量所在的位置
        max_index = np.argmax(abs(x_1))
        eigenvalue = x_1[max_index]
        eigenvalues.append(eigenvalue)
        x_0 = x_1 / eigenvalue
        eigen_vector.append(x_0)
        if len(eigenvalues) >= 2:
            if abs(eigenvalues[-1] - eigenvalues[-2]) < 1e-5:
                break

    return eigenvalues[-1], eigen_vector[-1]

def Subject1():
    A = [[5,-4, 1], [-4, 6, -4], [1, -4, 7]]
    A = np.array(A)
    eigenvalue, eigen_vector = Power(A)
    print('eigenvalue for A:', eigenvalue)
    print('eigen_vector for A:', eigen_vector)
    
    B = [[25, -41, 10, -6], [-41, 68, -17, 10], [10, -17, 5, -3], [-6, 10, -3, 2]]
    B = np.array(B)
    eigenvalue, eigen_vector = Power(B)
    print('eigenvalue for B:', eigenvalue)
    print('eigen_vector for B:', eigen_vector)

eps = 1e-10

def isZero(x):
    return abs(x) < eps

def isSchur(A):
    n = A.shape[0]
    nowIndex = 0
    while nowIndex < n - 1:
        if isZero(A[nowIndex][nowIndex]):
            nowIndex += 1
            continue
        # 先找这个元素对应的分块矩阵的阶数
        if ~isZero(A[nowIndex][nowIndex + 1]) or ~isZero(A[nowIndex + 1][nowIndex]):
            targetIndex = nowIndex + 1
        else:
            targetIndex = nowIndex
        
        # 判断除去这个分块之后的行和列是否均为 0
        for i in range(nowIndex, targetIndex + 1):
            for j in range(targetIndex + 1, n):
                if ~isZero(A[i][j]):
                    return False
        for i in range(targetIndex + 1, n):
            for j in range(nowIndex, targetIndex + 1):
                if ~isZero(A[i][j]):
                    return False
        
        nowIndex = targetIndex + 1
    return True


def QR_Algorithm(A):
    n = A.shape[0]
    num = 0
    while not isSchur(A):
        Q, R = np.linalg.qr(A)
        A = np.dot(R, Q) 
        num += 1
        if num % 100000 == 0:
            print('Schur:', A, "num:", num)
   
    # 根据 A 的对角块求特征值
    eigenvalues = []
    nowIndex = 0
    while nowIndex < n:
        if nowIndex == n - 1 or (A[nowIndex][nowIndex + 1] == 0 and A[nowIndex + 1][nowIndex] == 0):
            eigenvalues.append(A[nowIndex][nowIndex])
            nowIndex += 1
        else:
            a = A[nowIndex][nowIndex]
            b = A[nowIndex][nowIndex + 1]
            c = A[nowIndex + 1][nowIndex]
            d = A[nowIndex + 1][nowIndex + 1]
            eigenvalues.append((a + d + np.sqrt((a + d) ** 2 - 4 * (a * d - b * c))) / 2)
            eigenvalues.append((a + d - np.sqrt((a + d) ** 2 - 4 * (a * d - b * c))) / 2)
            nowIndex += 2
    print('eigenvalues:', eigenvalues)

def updateQR(A):
    n = A.shape[0]
    while not isSchur(A):
        T = A - np.eye(n) * A[n-1][n-1]
        Q, R = np.linalg.qr(T)
        A = np.dot(R, Q) + np.eye(n) * A[n-1][n-1]
    # 根据 A 的对角块求特征值
    eigenvalues = []
    nowIndex = 0
    while nowIndex < n:
        if nowIndex == n - 1 or (A[nowIndex][nowIndex + 1] == 0 and A[nowIndex + 1][nowIndex] == 0):
            eigenvalues.append(A[nowIndex][nowIndex])
            nowIndex += 1
        else:
            a = A[nowIndex][nowIndex]
            b = A[nowIndex][nowIndex + 1]
            c = A[nowIndex + 1][nowIndex]
            d = A[nowIndex + 1][nowIndex + 1]
            eigenvalues.append((a + d + np.sqrt((a + d) ** 2 - 4 * (a * d - b * c))) / 2)
            eigenvalues.append((a + d - np.sqrt((a + d) ** 2 - 4 * (a * d - b * c))) / 2)
            nowIndex += 2
    print('eigenvalues:', eigenvalues)
def Subject2():
    A = [[0.5, 0.5, 0.5, 0.5], [0.5, 0.5, -0.5, -0.5], [0.5, -0.5, 0.5, -0.5], [0.5, -0.5, -0.5, 0.5]]
    A = np.array(A)
    QR_Algorithm(A)

def Subject3():
    A = [[0.5, 0.5, 0.5, 0.5], [0.5, 0.5, -0.5, -0.5], [0.5, -0.5, 0.5, -0.5], [0.5, -0.5, -0.5, 0.5]]
    A = np.array(A)
    updateQR(A)

if __name__ == '__main__':
    print('Subject1:')
    Subject1()
    print('Subject2:')
    # Subject2()
    print('Subject2 need to run for a long time, so I comment it')
    print('Subject3:')
    Subject3()
    
