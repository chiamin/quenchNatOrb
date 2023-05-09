# ifndef __FOURIERTRANSF_H_CMC__
# define __FOURIERTRANSF_H_CMC__
#include <cmath>
#include "ContainerUtility.h"

CMatrix Fourier_matrix (int N)
{
    CMatrix M (N,N);
    for(int kn = 0; kn < N; kn++)
        for(int xn = 0; xn < N; xn++)
            M(kn,xn) = std::exp(-1i * 2 * M_PI * kn * xn / Real(N));
    M *= 1./std::sqrt(N);
    return M;
}

Real Fourier_k (int kn, int N)
{
    // In unit of pi
    Real k = 2 * kn / Real(N);
    // Shift k to be in range (-pi,pi) instead of (0,2pi)
    if (k > 1.)
         k -= 2.;
    return k;
}

// Find the dominate frequency in unit of pi
int find_kn (const CVector& fx, const vector<int>& exclude_kn)
{
    int N = fx.size();
    CMatrix U = Fourier_matrix (N);
    CVector fk = U * fx;
    vector<Real> ws, kns;
    for(int kn = 0; kn < N; kn++)
    {
        if (!iut::in_vector(exclude_kn, kn))
        {
            ws.push_back (std::norm(fk(kn)));
            kns.push_back (kn);
        }
    }
    auto ind = iut::max_index (ws);
    return kns.at(ind);
}
#endif
