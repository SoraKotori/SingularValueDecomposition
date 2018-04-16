#include "mkl.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <functional>
#include <vector>

void print_matrix(lapack_int m, lapack_int n, float* a)
{
    for (lapack_int i = 0; i < m; ++i)
    {
        for (lapack_int j = 0; j < n; ++j)
        {
            printf("% 6.2f", a[i * n + j]);
        }
        printf("\n");
    }
}

template<int matrix_layout, char jobu, char jobvt>
class IntelSingleSingularValueDecomposition
{
public:
    IntelSingleSingularValueDecomposition() = default;
    ~IntelSingleSingularValueDecomposition() = default;

    IntelSingleSingularValueDecomposition(lapack_int m, lapack_int n)
    {

    }

    void Compute()
    {

    }
private:
    lapack_int m;
    lapack_int n;
};


lapack_int intelsvd(lapack_int m, lapack_int n, float *a, float *s, float *u, float *vt, float *superb)
{
    const int matrix_layout = LAPACK_ROW_MAJOR;
    const char jobu = 'A';
    const char jobvt = 'A';

    const lapack_int lda = n;
    const lapack_int ldu = m;
    const lapack_int ldvt = n;

    lapack_int info = LAPACKE_sgesvd(
        matrix_layout, jobu, jobvt, m, n,
        a, lda,
        s,
        u, ldu,
        vt, ldvt,
        superb);

    return info;
}

lapack_int intelbrd(lapack_int m, lapack_int n, float* a, float* d, float* e, float* tauq, float* taup)
{
    const int matrix_layout = LAPACK_ROW_MAJOR;
    const lapack_int lda = n;

    lapack_int info = LAPACKE_sgebrd(
        matrix_layout, m, n,
        a, lda,
        d,
        e,
        tauq,
        taup);

    return info;
}

lapack_int intelgbr(char vect, lapack_int m, lapack_int n, float* a, const float* tau)
{
    const int matrix_layout = LAPACK_ROW_MAJOR;
    const lapack_int lda = n;

    lapack_int info = LAPACKE_sorgbr(
        matrix_layout, vect,
        vect == 'Q' ? m : n, vect == 'Q' ? m : n, vect == 'Q' ? n : m,
        a, lda,
        tau);

    return info;
}

lapack_int intelgbr2(lapack_int m, lapack_int n, float* a, float* tauq, float* taup, float *vt, float *u)
{
    if (m > n)
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                u[i * m + j] = a[i * n + j];
            }
        }

        lapack_int info = intelgbr('Q', m, n, u, tauq);
        if (0 != info)
        {
            return info;
        }

        info = intelgbr('P', m, n, a, taup);
        if (0 != info)
        {
            return info;
        }

        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                vt[i * n + j] = a[i * n + j];
            }
        }

        return info;
    }
    else
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                vt[i * n + j] = a[i * n + j];
            }
        }

        lapack_int info = intelgbr('P', m, n, vt, taup);
        if (0 != info)
        {
            return info;
        }

        info = intelgbr('Q', m, n, a, tauq);
        if (0 != info)
        {
            return info;
        }

        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < m; ++j)
            {
                u[i * m + j] = a[i * n + j];
            }
        }

        return info;
    }
}

lapack_int intelsqr(lapack_int m, lapack_int n, float* d, float* e, float *vt, float *u)
{
    const int matrix_layout = LAPACK_ROW_MAJOR;
    const lapack_int ncvt = n;
    const lapack_int nru = m;
    const lapack_int ldvt = n;
    const lapack_int ldu = m;

    lapack_int info = LAPACKE_sbdsqr(
        matrix_layout, m >= n ? 'U' : 'L', std::min(m, n),
        ncvt, nru, 0,
        d,
        e,
        vt, ldvt,
        u, ldu,
        nullptr, 0);

    return info;
}

lapack_int myintelsvd(lapack_int m, lapack_int n, float *a, float *s, float *u, float *vt, float *superb)
{
    lapack_int info;
    return info;
}

namespace svd
{
    template<typename _Type>
    using VectorType = std::vector<_Type>;

    template<typename _Type, typename _EngineType = std::default_random_engine>
    class SingularValueDecomposition
    {
    public:
        using SizeType = typename VectorType<_Type>::size_type;

        SingularValueDecomposition() = default;
        ~SingularValueDecomposition() = default;

        SingularValueDecomposition(SizeType m, SizeType n) :
            u(m),
            v(m),
            P(m, VectorType<_Type>(m)),
            P(n, VectorType<_Type>(n)),
            _ColumnCount(n)
        {
        }

    private:
        VectorType<_Type> u;
        VectorType<_Type> v;

        VectorType<VectorType<_Type>> P;
        VectorType<VectorType<_Type>> Q;

        SizeType _ColumnCount;
        _EngineType _Engine;

        template<typename _ForwardIterator>
        void GenerateUnitVector(_ForwardIterator&& _First, _ForwardIterator&& _Last)
        {
            _Type _Sum = 0;
            std::for_each(_First, _Last, [&](auto&& _Value)
            {
                _Value = _Engine();
                _Sum += std::pow(_Value, 2);
            });

            std::for_each(_First, _Last, [_Sqrt = std::sqrt(_Sum)](auto&& _Value)
            {
                _Value /= _Sqrt;
            });
        }

        template<typename _MatrixType>
        void Lanczos(_MatrixType&& _A)
        {
            _Type _Beta = 0;
            GenerateUnitVector(std::begin(Q[0]), std::end(Q[0]));


        }

        template<typename _MatrixType>
        void HouseholderTransformations(_MatrixType&& _A)
        {
            auto RowCount = u.size();
            for (decltype(_ColumnCount) _ColumnIndex = 0; _ColumnIndex < _ColumnCount; ++_ColumnIndex)
            {
                std::fill_n(std::begin(u), _ColumnIndex, _Type(0));
                std::fill_n(std::begin(v), _ColumnIndex, _Type(0));

                _Type _SumSquare = 0;
                for (decltype(RowCount) RowIndex = _ColumnIndex; RowIndex < RowCount; ++RowIndex)
                {
                    u[RowIndex] = _A[RowIndex][_ColumnIndex];
                    _SumSquare += std::pow(u[RowIndex], 2);
                }
                auto _Alpha = u[_ColumnIndex] < 0 ? std::sqrt(_SumSquare) : -std::sqrt(_SumSquare);

                _SumSquare = 0;
                for (decltype(RowCount) RowIndex = _ColumnIndex; RowIndex < RowCount; ++RowIndex)
                {
                    v[RowIndex] = _ColumnIndex == RowIndex ? u[RowIndex] + _Alpha : u[RowIndex];
                    _SumSquare += std::pow(v[RowIndex], 2);
                }

                _SumSquare = std::sqrt(_SumSquare);
                if (_SumSquare < 0.0000000001)
                {
                    continue;
                }

                for (decltype(RowCount) RowIndex = _ColumnIndex; RowIndex < RowCount; ++RowIndex)
                {
                    v[RowIndex] /= _SumSquare;
                }

            }
        }
    };
}

void IntelCombination()
{
    float d[std::min(m, n)];
    float e[std::min(m, n) - 1];
    float tauq[std::min(m, n)];
    float taup[std::min(m, n)];
    lapack_int info = intelbrd(m, n, a, d, e, tauq, taup);

    print_matrix(m, n, a);
    printf("\n");
    print_matrix(std::min(m, n), 1, d);
    printf("\n");
    print_matrix(std::min(m, n) - 1, 1, e);
    printf("\n");
    print_matrix(std::min(m, n), 1, tauq);
    printf("\n");
    print_matrix(std::min(m, n), 1, taup);
    printf("\n");

    printf("-------------------------\n");

    float vt[n * n];
    float u[m * m];
    info = intelgbr2(m, n, a, tauq, taup, vt, u);

    print_matrix(m, n, a);
    printf("\n");
    print_matrix(std::min(m, n), 1, tauq);
    printf("\n");
    print_matrix(std::min(m, n), 1, taup);
    printf("\n");
    print_matrix(n, n, vt);
    printf("\n");
    print_matrix(m, m, u);
    printf("\n");

    printf("-------------------------\n");

    info = intelsqr(m, n, d, e, vt, u);
    print_matrix(m, n, a);
    printf("\n");
    print_matrix(std::min(m, n), 1, d);
    printf("\n");
    print_matrix(std::min(m, n) - 1, 1, e);
    printf("\n");
    print_matrix(n, n, vt);
    printf("\n");
    print_matrix(m, m, u);
    printf("\n");
}

int main(void)
{
    const lapack_int m = 4;
    const lapack_int n = 5;

    float a[m * n] =
    {
        1, 0, 0, 0, 2,
        0, 0, 3, 0, 0,
        0, 0, 0, 0, 0,
        0, 4, 0, 0, 0
    };
    print_matrix(m, n, a);
    printf("\n");

    float s[std::min(m, n)];
    float u[m * m];
    float vt[n * n];
    float superb[std::min(m, n) - 1];
    lapack_int info = intelsvd(m, n, a, s, u, vt, superb);

    print_matrix(m, m, u);
    printf("\n");
    print_matrix(std::min(m, n), 1, s);
    printf("\n");
    print_matrix(n, n, vt);
    printf("\n");
    print_matrix(std::min(m, n) - 1, 1, superb);
    printf("\n");

    printf("-------------------------\n");


    return 0;
}