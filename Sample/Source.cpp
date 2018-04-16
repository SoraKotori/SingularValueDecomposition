#include "SingularValueDecomposition.h"
#include <cstdlib>
#include <iostream>

using namespace std;

//int main(void)
//{
//
//    return EXIT_SUCCESS;
//}

int main(void)
{
	int Rows = 4, Cols = 3;

	float **Mat = new float*[Rows];
	for (int row = 0; row < Rows; ++row)
		Mat[row] = new float[Cols];

	int aa = 0;
	cout << "H =>" << endl;
	for (int row = 0; row < Rows; ++row)
	{
		for (int col = 0; col < Cols; ++col)
		{
			aa += 1;
			Mat[row][col] = (float)aa;
			printf("%lf, ", Mat[row][col]);
		}
		cout << endl;
	}
	cout << endl;

	//Implementation
	SVD Svd(Rows, Cols);
	Svd.Decompose(Mat);

	cout << "U =>" << endl;
	for (int row = 0; row < Rows; ++row)
	{
		for (int col = 0; col < Rows; ++col)
		{
			printf("%lf, ", Svd.U_Mat[row][col]);
		}
		printf("\n");
	}
	cout << endl;

	cout << "W =>" << endl;
	for (int row = 0; row < Rows; ++row)
	{
		for (int col = 0; col < Cols; ++col)
		{
			printf("%lf, ", Svd.W_Mat[row][col]);
		}
		printf("\n");
	}
	cout << endl;

	cout << "V(not Vt) =>" << endl;
	for (int row = 0; row < Cols; ++row)
	{
		for (int col = 0; col < Cols; ++col)
			printf("%lf, ", Svd.V_Mat[row][col]);
		printf("\n");
	}
	cout << endl;

	return 0;
}
