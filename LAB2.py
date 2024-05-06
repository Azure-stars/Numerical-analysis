# 编程实现牛顿法与牛顿下山法求解方程
# 认为迭代判停绝对误差小于 1e-4，且迭代次数不超过 1000000 次
def Newton(f, df, x0, eps):
    xk = x0
    ans = []
    # 如果超过 1000000 次迭代，认为无法收敛
    Num = 1000000
    while True:
        x_k1 = xk - f(xk) / df(xk)
        if abs(x_k1 - xk) < eps:
            break
        ans.append(x_k1)
        if len(ans) > Num:
            print("Error: Newton Method failed to converge")
            return None
        xk = x_k1
    print("Newton Method Solution:")
    for index, i in enumerate(ans):
        print("Iteration i: {} val: {}".format(index, i))

    return x_k1

# 牛顿下山因子，认为有 100 项
NewtonDownhillFactor = [1/x for x in range(1, 101)]

def NewtonDownhill(f, df, x0, eps):
    xk = x0
    ans = []
    while True:
        factor_index = 0
        s = f(xk) / df(xk)
        x_k1 = xk - s
        while (abs(f(x_k1)) >= abs(f(xk))):
            factor_index += 1
            if (factor_index >= 100):
                print("Error: Newton Downhill Method failed to converge")
                return None
            x_k1 = xk - s * NewtonDownhillFactor[factor_index]
        ans.append((x_k1, factor_index))
        
        if abs(x_k1 - xk) < eps:
            break
        xk = x_k1
    print("Newton Downhill Method Solution:")
    for index, i in enumerate(ans):
        print("Iteration i: {} val: {} factor: {}".format(index, i[0], i[1]) )
    return x_k1

def f(x):
    return x**3 - 2*x + 2

def f_prime(x):
    return 3*x**2 - 2

def g(x):
    return -x**3 + 5*x

def g_prime(x):
    return -3*x**2 + 5

def std_answer():
    import scipy
    from scipy.optimize import fsolve

    # 之所以起始为 - 1，是为了能够收敛
    print(fsolve(f, -1, fprime=f_prime))

    # 求解 g(x) = 0 的解
    print(fsolve(g, 1.35))

def Subject2():
    print("Newton Method for fx:")
    ans = Newton(f, f_prime, 0, 1e-4)
    print("Answer X: ", ans)

    print("Newton Downhill Method for fx:")
    ans = NewtonDownhill(f, f_prime, 0, 1e-4)
    print("Answer X: ", ans)

    print("Newton Method for gx:")
    ans = Newton(g, g_prime, 1.35, 1e-4)
    print("Answer X: ", ans)

    print("Newton Downhill Method for gx:")
    ans = NewtonDownhill(g, g_prime, 1.35, 1e-4)
    print("Answer X: ", ans)


fzero_eps = 1e-5


def fzerotx(f, a, b):
    fa = f(a)
    fb = f(b)
    if (fa * fb > 0):
        print("Error: f(a) * f(b) > 0")
        return None
    c = a
    fc = fa
    d = b - c
    e = d
    while (fb != 0):
        if (fa * fb > 0):
            a = c
            fa = fc
            d = b - c
            e = d
        if (abs(fa) < abs(fb)):
            c = b
            b = a
            a = c
            fc = fb
            fb = fa
            fa = fc
        m = 0.5 * (a - b)
        tol = 2.0 * fzero_eps * max(abs(b), 1.0)
        if (abs(m) <= tol or fb == 0):
            break
        if (abs(e) < tol or abs(fc) <= abs(fb)):
            d = m
            e = m
        else:
            s = fb / fc
            if (a == c):
                p = 2.0 * m * s
                q = 1.0 - s
            else:
                q = fc / fa
                r = fb / fa
                p = s * (2.0 * m * q * (q - r) - (b - c) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            if (p > 0):
                q = -q
            p = abs(p)
            if (2.0 * p < 3.0 * m * q - abs(tol * q) and p < abs(0.5 * e * q)):
                e = d
                d = p / q
            else:
                d = m
                e = m
        c = b
        fc = fb
        if (abs(d) > tol):
            b += d
        else:
            b += tol if m > 0 else -tol
        fb = f(b)
    return b

def getRootForBessel():
    import scipy
    from scipy.special import jv
    # 第一类的零阶贝塞尔函数
    f = lambda x: jv(0, x)

    # 求其前 10 个正的零点
    root_range = [(2.0,2.5),(5.5,6.0),(8.5,9.0),(11.5,12.0),(14.5,15.0),(18.0,18.5),(21.0,21.5),(24.0,24.5),(27.0,27.5),(30.5,31.0)]

    answer = []
    for i in range(10):
        print("Root", i+1, end=' ')
        val = fzerotx(f, root_range[i][0], root_range[i][1])
        print(val)
        answer.append(val)

    # 绘制 bessel 曲线，同时标记这十个零点
    import numpy
    import matplotlib.pyplot as plt
    # x 需要加入这十个零点
    # 0 到 32 取 0, 0.5, 1, 1.5 等类型，共 65 个点
    answer_index = 0
    x = []
    slot_index = []
    for i in range(65):
        x.append(i * 0.5)
        if (answer_index < 10 and i * 0.5 < answer[answer_index] and (i+1) * 0.5 >= answer[answer_index]):
            slot_index.append(len(x))
            x.append(answer[answer_index])
            answer_index += 1
    x = numpy.array(x)
    y = jv(0, x)
    plt.plot(x, y)
    plt.scatter(answer, [0]*10, color='red')
    for i in range(10):
        plt.text(answer[i], 0, f"Root{i+1}={round(x[slot_index[i]],2)}", fontsize=5)
    plt.xlabel('x')
    plt.ylabel('J0(x)')
    plt.title('Bessel Function J0')
    plt.savefig('./figures/LAB2.png')
    # plt.show()

if __name__ == "__main__":
    print("Subject2:")
    Subject2()
    print("Subject3:")
    getRootForBessel()
    print("The result has been saved to ./figures/LAB2.png")