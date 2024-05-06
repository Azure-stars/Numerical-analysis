# LAB 5 矩阵特征值计算 REPORT

> 姓名：郑友捷		学号：2021010771		班级：计 14

## 上机题 1

### 实验思路

- 幂法的实现参照课本上的实现即可
- 不采用瑞利商加速，只进行规格化操作。每次迭代过程中矩阵和向量点乘得到的结果作为本次迭代得到的特征向量，其绝对值最大的分量作为本次迭代得到的特征值。
- 注意规格化除去的分量可以是负数，只要保证其绝对值最大即可
- 当相邻两次的特征值绝对值之差小于 1e-5 时，结束循环



### 实验代码

```python
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
```

### 实验输出

```sh
Subject1:
eigenvalue for A: 12.254320584751564
eigen_vector for A: [-0.67401981  1.         -0.88955964]
eigenvalue for B: 98.52169772379699
eigen_vector for B: [-0.60397234  1.         -0.25113513  0.14895345]
```

## 上机题 2

### 实验思路

- QR 分解的算法思路参考课本上的实现，调用 python 中对应的函数即可。
- 需要自己手动实现判断矩阵是否为 Schur 矩阵，这里采用了较为朴素的形式。直接判断当前的每一个分块矩阵是否为 1 或者 2 阶。若不是则返回 False。
- 完成分解之后，由于结果为实 Schur 阵。因此特征值较容易求，分为一阶块和二阶块分类讨论即可。

### 实验代码

```python
# 判断为 0 的阈值
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
```



### 实验输出

运行代码，发现卡死。调试输出发现收敛速度极慢。


迭代了  8000000 次之后，A 的结果如下：

```sh
[ 0.49937142  0.49999838  0.50011893  0.5005106 ]
[ 0.49999838  0.50063102 -0.49948938 -0.49988055]
[ 0.50011893 -0.49948938  0.5003902  -0.50000107]
[ 0.5005106  -0.49988055 -0.50000107  0.49960736]
```

分析可知，A 本身就是一个正定阵。在进行 QR 分解时，若 R 为对角阵，可以设置 Q = A，R 为单位阵。此时 $A_{k+1} = RQ = A_{k}$。因此理论上会死循环。实际计算时由于误差的存在，会导致 $A_{k+1}$和 $A_{k}$不同，但是因为误差而导致的区别较小，误差积累速度很慢，很难收敛。



## 上机题 3

### 实验思路

- 原点位移采用课件上提到的 $s_k = A_k(n,n)$的计算方法
- 计算公式变为：$QR=A_k-s_kI_n, A_{k+1}=RQ + s_kI_n$。

### 实验代码

```python
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
```

### 实验输出

```sh
Subject3:
eigenvalues: [1.0, -0.9999999999999998, 1.0, 1.0]
```

### 实验分析

同一个矩阵，原始的 QR 算法几乎无法收敛，而带原点位移的QR算法收敛速度较快。可以看出原点位移的QR算法由于引入了偏移量，防止了参与 QR 分解的矩阵为正定阵导致的死循环问题，可以在一定程度上解决 QR 算法存在的无法收敛的问题。
