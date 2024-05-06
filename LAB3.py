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

if __name__ == "__main__":
    Subject6()