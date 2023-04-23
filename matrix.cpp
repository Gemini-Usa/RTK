#include <format>
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "function.h"

/// <summary>
/// 向量相加 c=a+b
/// </summary>
/// <param name="a">传入向量1</param>
/// <param name="b">传入向量2</param>
/// <param name="c">传出向量</param>
/// <param name="m">向量维数</param>
void VectorAdd(const double* a, const double* b, double* c, int m)
{
	for (int i = 0; i < m; i++)
	{
		c[i] = a[i] + b[i];
	}
}

/// <summary>
/// 矩阵相减 c=a-b
/// </summary>
/// <param name="a">传入向量1</param>
/// <param name="b">传入向量2</param>
/// <param name="c">传出向量</param>
/// <param name="m">向量维数</param>
void VectorSub(const double* a, const double* b, double* c, int m)
{
	for (int i = 0; i < m; i++)
	{
		c[i] = a[i] - b[i];
	}
}

/// <summary>
/// 向量数乘 b=scalar*a
/// </summary>
/// <param name="a">传入向量</param>
/// <param name="b">传出向量</param>
/// <param name="scalar">传入标量</param>
/// <param name="m">向量维数</param>
void VectorScalarMult(const double* a, double* b, double scalar, int m)
{
	for (int i = 0; i < m; i++)
	{
		b[i] = scalar * a[i];
	}
}

/// <summary>
/// 向量点乘 a·b or aTb
/// </summary>
/// <param name="a">传入向量1</param>
/// <param name="b">传入向量2</param>
/// <param name="m">向量维数</param>
/// <returns>a·b or aTb</returns>
double VectorInnerProduct(const double* a, const double* b, int m)
{
	double res = 0.0;
	for (int i = 0; i < m; i++)
	{
		res += a[i] * b[i];
	}
	return res;
}

/// <summary>
/// 向量(欧式)范数 ||a||
/// </summary>
/// <param name="a">传入向量</param>
/// <param name="m">向量维数</param>
/// <returns>||a||</returns>
double Norm(const double* a, int m)
{
	double res = 0.0;
	for (int i = 0; i < m; i++)
	{
		res += a[i] * a[i];
	}
	return sqrt(res);
}

/// <summary>
/// 向量归一化 b=a/||a||
/// </summary>
/// <param name="a">传入向量</param>
/// <param name="b">传出向量</param>
/// <param name="m">向量维数</param>
void VectorNormalize(const double* a, double* b, int m)
{
	double norm = Norm(a, m);
	for (int i = 0; i < m; i++)
	{
		b[i] = a[i] / norm;
	}
}

/// <summary>
/// 向量外积 c=a×b
/// </summary>
/// <param name="a">传入向量1</param>
/// <param name="b">传入向量2</param>
/// <param name="c">传出向量</param>
void VectorOuterProduct(const double* a, const double* b, double* c)
{
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

/// <summary>
/// 矩阵相加 c=a+b
/// </summary>
/// <param name="a">传入矩阵1∈R{m×n}</param>
/// <param name="b">传入矩阵2∈R{m×n}</param>
/// <param name="c">传出向量∈R{m×n}</param>
/// <param name="m">矩阵行数</param>
/// <param name="n">矩阵列数</param>
void MatrixAdd(const double* a, const double* b, double* c, int m, int n)
{
	for (int i = 0; i < m * n; i++)
	{
		c[i] = a[i] + b[i];
	}
}

/// <summary>
/// 矩阵相减 c=a-b
/// </summary>
/// <param name="a">传入矩阵1∈R{m×n}</param>
/// <param name="b">传入矩阵2∈R{m×n}</param>
/// <param name="c">传出矩阵∈R{m×n}</param>
/// <param name="m">矩阵行数</param>
/// <param name="n">矩阵列数</param>
void MatrixSub(const double* a, const double* b, double* c, int m, int n)
{
	for (int i = 0; i < m * n; i++)
	{
		c[i] = a[i] - b[i];
	}
}

/// <summary>
/// 矩阵乘法 c=ab
/// </summary>
/// <param name="a">传入矩阵1∈R{m×n}</param>
/// <param name="b">传入矩阵2∈R{n×s}</param>
/// <param name="c">传出矩阵∈R{m×s}</param>
/// <param name="m">矩阵1行数</param>
/// <param name="n">矩阵1列数,矩阵2行数</param>
/// <param name="s">矩阵2列数</param>
void MatrixMultiply(const double* a, const double* b, double* c, int m, int n, int s)
{
	int num = 0;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < s; j++)
		{
			double sum = 0;
			for (int k = 0; k < n; k++)
			{
				sum = sum + a[i * n + k] * b[k * s + j];
			}
			c[num] = sum;
			num++;
		}
	}
}

/// <summary>
/// 矩阵转置 b=aT
/// </summary>
/// <param name="a">传入矩阵∈R{m×n}</param>
/// <param name="b">传出矩阵∈R{n×m}</param>
/// <param name="m">矩阵行数</param>
/// <param name="n">矩阵列数</param>
void MatrixTranspose(const double* a, double* b, int m, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			b[i * m + j] = a[j * n + i];
		}
	}
}

/// <summary>
/// 矩阵增广 b=[a I]
/// </summary>
/// <param name="a">传入矩阵∈R{m×m}</param>
/// <param name="b">传出矩阵∈R{m×2m}</param>
/// <param name="m">矩阵行数,矩阵列数</param>
void MatrixExtend(const double* a, double* b, int m)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < 2 * m; j++)
		{
			if (j < m)
			{
				b[2 * m * i + j] = a[m * i + j];
			}
			else
			{
				if (j == m + i)
				{
					b[2 * m * i + j] = 1.0;
				}
				else
				{
					b[2 * m * i + j] = 0.0;
				}
			}
		}
	}
}

/// <summary>
/// 矩阵行交换
/// </summary>
/// <param name="a">传入矩阵∈R{m×n}</param>
/// <param name="m">矩阵行数</param>
/// <param name="n">矩阵列数</param>
/// <param name="row1">交换行1</param>
/// <param name="row2">交换行2</param>
void MatrixRowExchange(double* a, int m, int n, int row1, int row2)
{
	double temp = 0.0;
	for (int i = 0; i < n; i++)
	{
		temp = a[(row2 - 1) * n + i];
		a[(row2 - 1) * n + i] = a[(row1 - 1) * n + i];
		a[(row1 - 1) * n + i] = temp;
	}
}

/// <summary>
/// 矩阵求逆 b=a^-1
/// </summary>
/// <param name="a">传入矩阵∈R{m×m}</param>
/// <param name="b">传出矩阵∈R{m×m}</param>
/// <param name="m">矩阵行数,矩阵列数</param>
void MatrixInverse(const double* a, double* b, int m)
{
	double* ext = new double[m * 2 * m];
	MatrixExtend(a, ext, m);
	for (int i = 0; i < m; i++)
	{
		if (fabs(ext[i * 2 * m + i]) < 1E-9)
		{
			int j;
			for (j = 0; j < m; j++)
			{
				if (fabs(ext[j * 2 * m + i]) > 1E-9)
				{
					MatrixRowExchange(ext, m, 2 * m, i, j);
					break;
				}
			}
			if (j >= m)
			{
				b = nullptr;
				return;
			}
		}
	}
	for (int i = 0; i < m; i++)
	{
		double pivot = ext[i * 2 * m + i];
		for (int j = 0; j < 2 * m; j++)
		{
			ext[i * 2 * m + j] /= pivot;
		}
		for (int k = 0; k < m; k++)
		{
			if (k == i) continue;
			double factor = ext[k * 2 * m + i];
			for (int l = 0; l < m * 2; l++)
			{
				ext[k * 2 * m + l] -= ext[i * 2 * m + l] * factor;
			}
		}
	}
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < 2 * m; j++)
		{
			if (j >= m)
			{
				b[i * m + j - m] = ext[i * m * 2 + j];
			}
		}
	}
	delete[] ext;
	ext = nullptr;
}

/// <summary>
/// 打印矩阵:向命令行打印矩阵a
/// </summary>
/// <param name="a">打印矩阵∈R{m×n}</param>
/// <param name="m">矩阵行数</param>
/// <param name="n">矩阵列数</param>
void MatrixPrint(const double* a, int m, int n)
{
	for (int i = 0; i < m * n; i++) {
        std::cout << std::format("{:8.4f},\t", a[i]);
        if ((i + 1) % n == 0) std::cout << std::endl;
    }
}