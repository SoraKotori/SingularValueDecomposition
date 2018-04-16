#pragma once
#include <cmath>
#include <iostream>
#include <vector>

#define MAX(x,y) ((x)>(y)?(x):(y))
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

using namespace std;

class SVD
{
public:
	SVD(int InputMatRows, int InputMatCols) : MatRows(InputMatRows), _ColumnCount(InputMatCols),
		U_Size(InputMatRows), W_Size_Rows(InputMatRows),
		W_Size_Cols(InputMatCols), V_Size(InputMatCols)
	{
		rv1 = new double[_ColumnCount];
		W = new float[_ColumnCount];

		U_Mat = new float*[MatRows];
		for (int row = 0; row < MatRows; ++row)
			U_Mat[row] = new float[MatRows];

		W_Mat = new float*[MatRows];
		for (int row = 0; row < MatRows; ++row)
			W_Mat[row] = new float[_ColumnCount];

		V_Mat = new float*[_ColumnCount];
		for (int row = 0; row < _ColumnCount; ++row)
			V_Mat[row] = new float[_ColumnCount];
	}

	~SVD()
	{
		delete[] rv1;
		delete[] W;

		for (int row = 0; row < MatRows; ++row)
			delete[] U_Mat[row];
		delete[] U_Mat;

		for (int row = 0; row < MatRows; ++row)
			delete[] W_Mat[row];
		delete[] W_Mat;

		for (int row = 0; row < _ColumnCount; ++row)
			delete[] V_Mat[row];
		delete[] V_Mat;
	}

	int Decompose(float **InputMat)
	{
		int flag, its, j, jj, k, l, nm;
		double c, f, h, s, x, y, z;
		double anorm = 0.0, g = 0.0, scale = 0.0;

		if (MatRows < _ColumnCount)
		{
			fprintf(stderr, "#rows must be > #cols \n");
			return 0;
		}

		/* Householder reduction to bidiagonal form */
		for (decltype(_ColumnCount) _ColumnIndex = 0; _ColumnIndex < _ColumnCount; _ColumnIndex++)
		{
			/* left-hand reduction */
			l = _ColumnIndex + 1;
			rv1[_ColumnIndex] = scale * g;
			g = s = scale = 0.0;

			if (_ColumnIndex < MatRows)
			{
				for (k = _ColumnIndex; k < MatRows; k++)
				{
					scale += fabs((double)InputMat[k][_ColumnIndex]);
				}

				if (scale)
				{
					for (k = _ColumnIndex; k < MatRows; k++)
					{
						InputMat[k][_ColumnIndex] = (float)((double)InputMat[k][_ColumnIndex] / scale);
						s += ((double)InputMat[k][_ColumnIndex] * (double)InputMat[k][_ColumnIndex]);
					}

					f = (double)InputMat[_ColumnIndex][_ColumnIndex];
					g = -SIGN(sqrt(s), f);
					h = f * g - s;
					InputMat[_ColumnIndex][_ColumnIndex] = (float)(f - g);

					if (_ColumnIndex != _ColumnCount - 1)
					{
						for (j = l; j < _ColumnCount; j++)
						{
							for (s = 0.0, k = _ColumnIndex; k < MatRows; k++)
							{
								s += ((double)InputMat[k][_ColumnIndex] * (double)InputMat[k][j]);
							}

							f = s / h;
							for (k = _ColumnIndex; k < MatRows; k++)
							{
								InputMat[k][j] += (float)(f * (double)InputMat[k][_ColumnIndex]);
							}
						}
					}
					for (k = _ColumnIndex; k < MatRows; k++)
					{
						InputMat[k][_ColumnIndex] = (float)((double)InputMat[k][_ColumnIndex] * scale);
					}
				}
			}
			W[_ColumnIndex] = (float)(scale * g);

			/* right-hand reduction */
			g = s = scale = 0.0;
			if (_ColumnIndex < MatRows && _ColumnIndex != _ColumnCount - 1)
			{
				for (k = l; k < _ColumnCount; k++)
				{
					scale += fabs((double)InputMat[_ColumnIndex][k]);
				}

				if (scale)
				{
					for (k = l; k < _ColumnCount; k++)
					{
						InputMat[_ColumnIndex][k] = (float)((double)InputMat[_ColumnIndex][k] / scale);
						s += ((double)InputMat[_ColumnIndex][k] * (double)InputMat[_ColumnIndex][k]);
					}
				
					f = (double)InputMat[_ColumnIndex][l];
					g = -SIGN(sqrt(s), f);
					h = f * g - s;
					InputMat[_ColumnIndex][l] = (float)(f - g);
					
					for (k = l; k < _ColumnCount; k++)
					{
						rv1[k] = (double)InputMat[_ColumnIndex][k] / h;
					}

					if (_ColumnIndex != MatRows - 1)
					{
						for (j = l; j < MatRows; j++)
						{
							for (s = 0.0, k = l; k < _ColumnCount; k++)
							{
								s += ((double)InputMat[j][k] * (double)InputMat[_ColumnIndex][k]);
							}

							for (k = l; k < _ColumnCount; k++)
							{
								InputMat[j][k] += (float)(s * rv1[k]);
							}
						}
					}

					for (k = l; k < _ColumnCount; k++)
					{
						InputMat[_ColumnIndex][k] = (float)((double)InputMat[_ColumnIndex][k] * scale);
					}
				}
			}
			anorm = MAX(anorm, (fabs((double)W[_ColumnIndex]) + fabs(rv1[_ColumnIndex])));
		}

		/* accumulate the right-hand transformation */
		for (decltype(_ColumnCount) _ColumnIndex = _ColumnCount - 1; _ColumnIndex >= 0; _ColumnIndex--)
		{
			if (_ColumnIndex < _ColumnCount - 1)
			{
				if (g)
				{
					for (j = l; j < _ColumnCount; j++)
						V_Mat[j][_ColumnIndex] = (float)(((double)InputMat[_ColumnIndex][j] / (double)InputMat[_ColumnIndex][l]) / g);
					/* double division to avoid underflow */
					for (j = l; j < _ColumnCount; j++)
					{
						for (s = 0.0, k = l; k < _ColumnCount; k++)
							s += ((double)InputMat[_ColumnIndex][k] * (double)V_Mat[k][j]);
						for (k = l; k < _ColumnCount; k++)
							V_Mat[k][j] += (float)(s * (double)V_Mat[k][_ColumnIndex]);
					}
				}
				for (j = l; j < _ColumnCount; j++)
					V_Mat[_ColumnIndex][j] = V_Mat[j][_ColumnIndex] = 0.0;
			}
			V_Mat[_ColumnIndex][_ColumnIndex] = 1.0;
			g = rv1[_ColumnIndex];
			l = _ColumnIndex;
		}

		/* accumulate the left-hand transformation */
		for (decltype(_ColumnCount) _ColumnIndex = _ColumnCount - 1; _ColumnIndex >= 0; _ColumnIndex--)
		{
			l = _ColumnIndex + 1;
			g = (double)W[_ColumnIndex];
			if (_ColumnIndex < _ColumnCount - 1)
				for (j = l; j < _ColumnCount; j++)
					InputMat[_ColumnIndex][j] = 0.0;
			if (g)
			{
				g = 1.0 / g;
				if (_ColumnIndex != _ColumnCount - 1)
				{
					for (j = l; j < _ColumnCount; j++)
					{
						for (s = 0.0, k = l; k < MatRows; k++)
							s += ((double)InputMat[k][_ColumnIndex] * (double)InputMat[k][j]);
						f = (s / (double)InputMat[_ColumnIndex][_ColumnIndex]) * g;
						for (k = _ColumnIndex; k < MatRows; k++)
							InputMat[k][j] += (float)(f * (double)InputMat[k][_ColumnIndex]);
					}
				}
				for (j = _ColumnIndex; j < MatRows; j++)
					InputMat[j][_ColumnIndex] = (float)((double)InputMat[j][_ColumnIndex] * g);
			}
			else
			{
				for (j = _ColumnIndex; j < MatRows; j++)
					InputMat[j][_ColumnIndex] = 0.0;
			}
			++InputMat[_ColumnIndex][_ColumnIndex];
		}

		/* diagonalize the bidiagonal form */
		for (k = _ColumnCount - 1; k >= 0; k--)
		{                             /* loop over singular values */
			for (its = 0; its < 30; its++)
			{                         /* loop over allowed iterations */
				flag = 1;
				for (l = k; l >= 0; l--)
				{                     /* test for splitting */
					nm = l - 1;
					if (fabs(rv1[l]) + anorm == anorm)
					{
						flag = 0;
						break;
					}
					if (fabs((double)W[nm]) + anorm == anorm)
						break;
				}
				if (flag)
				{
					c = 0.0;
					s = 1.0;
					for (decltype(_ColumnCount) _ColumnIndex = l; _ColumnIndex <= k; _ColumnIndex++)
					{
						f = s * rv1[_ColumnIndex];
						if (fabs(f) + anorm != anorm)
						{
							g = (double)W[_ColumnIndex];
							h = PYTHAG(f, g);
							W[_ColumnIndex] = (float)h;
							h = 1.0 / h;
							c = g * h;
							s = (-f * h);
							for (j = 0; j < MatRows; j++)
							{
								y = (double)InputMat[j][nm];
								z = (double)InputMat[j][_ColumnIndex];
								InputMat[j][nm] = (float)(y * c + z * s);
								InputMat[j][_ColumnIndex] = (float)(z * c - y * s);
							}
						}
					}
				}
				z = (double)W[k];
				if (l == k)
				{                  /* convergence */
					if (z < 0.0)
					{              /* make singular value nonnegative */
						W[k] = (float)(-z);
						for (j = 0; j < _ColumnCount; j++)
							V_Mat[j][k] = (-V_Mat[j][k]);
					}
					break;
				}
				if (its >= 30) {
					fprintf(stderr, "No convergence after 30,000! iterations \n");
					return(0);
				}

				/* shift from bottom 2 x 2 minor */
				x = (double)W[l];
				nm = k - 1;
				y = (double)W[nm];
				g = rv1[nm];
				h = rv1[k];
				f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
				g = PYTHAG(f, 1.0);
				f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

				/* next QR transformation */
				c = s = 1.0;
				for (j = l; j <= nm; j++)
				{
                    decltype(_ColumnCount) _ColumnIndex = j + 1;
					g = rv1[_ColumnIndex];
					y = (double)W[_ColumnIndex];
					h = s * g;
					g = c * g;
					z = PYTHAG(f, h);
					rv1[j] = z;
					c = f / z;
					s = h / z;
					f = x * c + g * s;
					g = g * c - x * s;
					h = y * s;
					y = y * c;
					for (jj = 0; jj < _ColumnCount; jj++)
					{
						x = (double)V_Mat[jj][j];
						z = (double)V_Mat[jj][_ColumnIndex];
						V_Mat[jj][j] = (float)(x * c + z * s);
						V_Mat[jj][_ColumnIndex] = (float)(z * c - x * s);
					}
					z = PYTHAG(f, h);
					W[j] = (float)z;
					if (z)
					{
						z = 1.0 / z;
						c = f * z;
						s = h * z;
					}
					f = (c * g) + (s * y);
					x = (c * y) - (s * g);
					for (jj = 0; jj < MatRows; jj++)
					{
						y = (double)InputMat[jj][j];
						z = (double)InputMat[jj][_ColumnIndex];
						InputMat[jj][j] = (float)(y * c + z * s);
						InputMat[jj][_ColumnIndex] = (float)(z * c - y * s);
					}
				}
				rv1[l] = 0.0;
				rv1[k] = f;
				W[k] = (float)x;
			}
		}

		ReturnLeftOrth(InputMat);
		ReturnSingularMat();

		return 1;
	}

	int U_Size = 0, V_Size = 0;
	int W_Size_Rows = 0, W_Size_Cols = 0;

	float **U_Mat = nullptr;
	float **V_Mat = nullptr;
	float **W_Mat = nullptr;

private:
	int MatRows = 0, _ColumnCount = 0;
	double *rv1 = nullptr;

	float *W = nullptr;

	double PYTHAG(double a, double b)
	{
		double at = fabs(a), bt = fabs(b), ct, result;

		if (at > bt) { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
		else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
		else result = 0.0;
		return result;
	}

	void ReturnLeftOrth(float **InputMat)
	{
		for (int row = 0; row < MatRows; ++row)
		{
			for (int col = 0; col < MatRows; ++col)
			{
				if (col < _ColumnCount)
					U_Mat[row][col] = InputMat[row][col];
				else
					U_Mat[row][col] = 0;
			}
		}
	}

	void ReturnSingularMat()
	{
		for (int row = 0; row < MatRows; ++row)
		{
			for (int col = 0; col < _ColumnCount; ++col)
			{
				if (row == col)
					W_Mat[row][col] = W[col];
				else
					W_Mat[row][col] = 0;
			}
		}
	}
};

namespace svd
{
	template<typename _Type>
	using VectorType = std::vector<_Type>;

	template<typename _Type>
    class SingularValueDecomposition
    {
    public:
        SingularValueDecomposition() = default;
        ~SingularValueDecomposition() = default;

    private:
		void HouseholderTransformations()
		{

		}
    };
}