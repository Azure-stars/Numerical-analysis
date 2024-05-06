# LAB 3 线性方程组的直接解法 REPORT

> 姓名：郑友捷		学号：2021010771		班级：计 14

## 上机题 6

### 实验思路

- Cholesky 分解算法直接使用课件上提到的直接分解算法，通过传入对称正定矩阵可以直接得出其对应的下三角矩阵
- 对于 Hilbert 矩阵，其问题的不稳定性很高，且由于阶数较高、数值较小，很容易出现数值误差问题。据本人测算，使用 `numpy.float128` 位数值，不加任何数值处理技巧地使用直接分解算法时，对 14 阶 Hilbert 矩阵计算到最后会出现对角线元素为**复数**的情况，即对负数开方导致出现了问题。而这个负数的出现是**计算误差逐步累积**导致的结果。
- 为了解决这个问题，我和刘明道助教沟通之后，决定采用**python 的  sympy 库**，通过分数的方式完成整个问题的精确计算，将问题从数值计算转化为精确计算。因此不会存在计算上的误差，即截断和舍入误差均为 0。而这个方式虽然改变了问题的稳定性，但是**对 Hilbert 矩阵和整个问题的病态性（敏感性）却没有本质上的影响。仍然可以通过添加扰动计算误差的形式来直观估计问题的敏感性**。



### 实验代码

```python
# Cholesky 分解算法
import numpy as np
import sympy as sp

def symbolic_cholesky(A):
    n = A.shape[0]
    L = sp.zeros(n, n)

    for i in range(n):
        for j in range(i + 1):
            if i == j:
                L[i, i] = sp.sqrt(A[i, i] - sum(L[i, k]**2 for k in range(i)))
            else:
                L[i, j] = (A[i, j] - sum(L[i, k]*L[j, k] for k in range(j))) / L[j, j]

    return L

def Hilbert(n):
    return sp.Matrix([[sp.Rational(1, i + j + 1) for j in range(n)] for i in range(n)])

# 求解 L L^T x = b
def solveTriangular(L, b):
    # 先令 y = L^T x
    # 求解 L y = b
    n = len(b)
    y = sp.Matrix([0 for i in range(n)])
    for i in range(n):
        y[i] = (b[i] - sum(L[i, j]*y[j] for j in range(i))) / L[i, i]

    # 求解 L^T x = y
    L_T = sp.transpose(L)
    x = sp.Matrix([0 for i in range(n)])
    for i in range(n-1, -1, -1):
        x[i] = (y[i] - sum(L_T[i, j]*x[j] for j in range(i+1, n))) / L_T[i, i]
    return x

def queryHilbert(H, n):
    L = symbolic_cholesky(H)
    # H = np.array(H).astype(np.float64)
    # L = np.array(L).astype(np.float64)
    x = sp.Matrix([1 for i in range(n)])
    b = H * x

    x_answer = solveTriangular(L, b)

    # 计算残差
    r = b - H * x_answer
    r_inf = max(abs(r))
    # 计算误差
    e = x - x_answer
    # 计算误差的无穷范数
    e_inf = max(abs(e))
    print("n =", n)
    print("Before adding perturbation:")
    print("Infinite norm of residual error is", r_inf)
    print("Infinite norm of error", e_inf)

    # 对 b 沿着绝对值最大的方向加一个扰动，值为 1e-7

    # 找到绝对值最大的元素
    max_index = np.argmax(abs(b))
    if b[max_index] > 0:
        b[max_index] += 1e-7
    else:
        b[max_index] -= 1e-7

    x_answer = solveTriangular(L, b)
    # 计算残差
    r = b - H * x_answer
    r_inf = max(abs(r))
    # 计算误差
    e = x - x_answer
    # 计算误差的无穷范数
    e_inf = max(abs(e))
    print("After adding perturbation:")
    print("Infinite norm of residual error is", r_inf)
    print("Infinite norm of error", e_inf)


def queryRoot(n):
    H = Hilbert(n)
    queryHilbert(H, n)

def Subject6():
    n = 8
    queryRoot(n)
    n = 10
    queryRoot(n)
    n = 12
    queryRoot(n)
    n = 14
    queryRoot(n)

```

### 实验输出

```sh
n = 8
Before adding perturbation:
Infinite norm of residual error is 0
Infinite norm of error 0
After adding perturbation:
Infinite norm of residual error is 1.11022302462516e-16
Infinite norm of error 0.0288287221496351
n = 10
Before adding perturbation:
Infinite norm of residual error is 0
Infinite norm of error 0
After adding perturbation:
Infinite norm of residual error is 4.44089209850063e-16
Infinite norm of error 0.961800623796269
n = 12
Before adding perturbation:
Infinite norm of residual error is 0
Infinite norm of error 0
After adding perturbation:
Infinite norm of residual error is 8.88178419700125e-16
Infinite norm of error 32.9843004869018
n = 14
Before adding perturbation:
Infinite norm of residual error is 0
Infinite norm of error 0
After adding perturbation:
Infinite norm of residual error is 8.43769498715119e-15
Infinite norm of error 572.367493793506
```



### 实验分析

- 可以看出，由于使用了 sympy 进行精确计算，因此在不添加扰动前，计算结果的误差均为 0。（关于为什么使用 sympy 库进行计算，详见`实验思路`部分）
- 当添加了扰动之后，可以看出 Hilbert 矩阵对应的计算残差和误差都出现了变化，随着矩阵阶数的增大而逐渐增大。其中残差的变化不是很明显，而误差变化十分明显。因此可以较为明显地看出 Hilbert 矩阵存在的病态性问题。
- 从条件数来看，$cond(H_n)=𝑂(\frac{(1+\sqrt{2})^{4n}}{\sqrt{n}}$。随着 𝑛 增大，求解 $𝐻^n𝑥=𝑏$ 将会变得十分敏感。另外计算误差的累计也会导致上述提到的出现复数的问题。

