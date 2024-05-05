# 绘制差商误差、总误差限、截断误差和舍入误差四个变量和步长 h 的关系图
# 取 x = 1
# f(x) = sin(x)
# f'(x) = cos(x)
# 差商 f'(x) = (f(x+h) - f(x)) / h
# 截断误差 = M * h / 2, M = max|f''(x)|, x = 1
# 舍入误差 = 2 * u / h, u 为一次计算误差，设置为 1e-16
# 总误差限 = 截断误差 + 舍入误差

import numpy
import matplotlib.pyplot as plt

def Subject1():
    # 枚举不同的 h，计算差商误差、总误差限、截断误差和舍入误差
    h_list = [1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1]

    x = 1

    f = numpy.sin
    f_prime = numpy.cos

    f_prime_diff = numpy.zeros(len(h_list))
    f_truncation_diff = numpy.zeros(len(h_list))
    f_roundoff_diff = numpy.zeros(len(h_list))
    f_total_diff = numpy.zeros(len(h_list))

    for i in range(len(h_list)):
        h = h_list[i]
        f_prime_diff[i] = abs((f(x+h) - f(x)) / h - f_prime(x))
        f_truncation_diff[i] = h / 2
        f_roundoff_diff[i] = 2 * 1e-16 / h
        f_total_diff[i] = f_truncation_diff[i] + f_roundoff_diff[i]

    # 绘制差商误差、总误差限、截断误差和舍入误差四个变量和步长 h 的关系图
    # 四条曲线在同一张图上绘制

    plt.plot(h_list, f_prime_diff, label='Actual Error')
    plt.plot(h_list, f_truncation_diff, label='Truncation Error')
    plt.plot(h_list, f_roundoff_diff, label='Roundoff Error')
    plt.plot(h_list, f_total_diff, label='Total Error')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('h')
    plt.ylabel('Error')
    plt.legend()
    # plt.show()
    plt.savefig('./figures/LAB1.png')

def Subject3():
    # 使用单精度浮点数计算无穷级数
    ans_32 = numpy.float32(0)
    # 计算无穷级数，并且判断当 n 等于几时，无穷级数的和不再变化
    n = 1
    while True:
        temp_ans = ans_32 + numpy.float32(1 / n)
        if ans_32 == temp_ans:
            break
        ans_32 = temp_ans
        n += 1
    print("单精度所用次数: {} ".format(n))
    print("单精度答案:")
    print(ans_32)
    # 使用双精度浮点数计算无穷级数
    # 计算前 n 项
    ans_64 = numpy.float64(0)
    for i in range(1, n + 1):
        ans_64 += numpy.float64(1 / i)
    print("双精度计算相同次数答案 {}".format(ans_64))


if __name__ == '__main__':
    Subject1()
    Subject3()