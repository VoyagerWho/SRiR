#include "LUDecomposition.h"

LUDecomposition::LUDecomposition(const Matrix& A)
    :u(A), l(Matrix(A.size()))
{
    unsigned n = u.size();
    for(unsigned i=0; i<n; ++i)
    {
        l[i][i]=1.0;
    }
    for(unsigned k=0; k<n-1; ++k)
    {
        for(unsigned i=k+1; i<n; ++i)
            l[i][k]=u[i][k]/u[k][k];
        for(unsigned j=k+1; j<n; ++j)
        {
            for(unsigned i=k+1; i<n; ++i)
                u[i][j] = u[i][j]-l[i][k]*u[k][j];
        }
    }
    for(unsigned i=0; i<n-1; ++i)
    {
        for(unsigned j=i+1; j<n; ++j)
        {
            u[j][i] = 0;
        }
    }
}

LUDecomposition::~LUDecomposition()
{
    //dtor
}

Matrix LUDecomposition::getL() const
{
    return l;
}

Matrix LUDecomposition::getU() const
{
    return u;
}

Matrix LUDecomposition::getLU() const
{
    Matrix lu(u);
    for(unsigned i=0; i<lu.size()-1; ++i)
    {
        for(unsigned j=i+1; j<lu.size(); ++j)
        {
            lu[j][i] = l[j][i];
        }
    }
    return lu;
}
