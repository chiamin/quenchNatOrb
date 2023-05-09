#include "itensor/all.h"
using namespace itensor;
using namespace std;

//---------------- Lazy ways -------------------------
template <typename SiteType, typename NumType>
Mat<NumType> exact_corr_by_MPO (const MPS& psi)
{
    int N = length (psi);
    auto sites = SiteType (siteInds(psi));
    auto corr = Mat<NumType> (N,N);
    for(int i = 1; i <= N; i++)
    for(int j = 1; j <= N; j++)
    {
        AutoMPO ampo (sites);
        ampo += 1.,"Cdag",i,"C",j;
        auto C = toMPO (ampo);
        if constexpr (is_same_v <NumType, Real>)
            corr(i-1,j-1) = inner (psi, C, psi);
        else
            corr(i-1,j-1) = innerC (psi, C, psi);
    }
    return corr;
}

template <typename SiteType, typename NumType>
Mat<NumType> exact_corr_by_inner (const MPS& psi)
{
    int N = length (psi);
    auto sites = SiteType (siteInds(psi));
    auto corr = Mat<NumType> (N,N);
    auto phi = psi;
    for(int i = 1; i <= N; i++)
    for(int j = 1; j <= N; j++)
    {
        phi.ref(j) *= sites.op("C",j);
        phi.ref(j).noPrime("Site");
        phi.ref(i) *= sites.op("Cdag",i);
        phi.ref(i).noPrime("Site");
        if constexpr (is_same_v <NumType, Real>)
            corr(i-1,j-1) = inner (psi, phi);
        else
            corr(i-1,j-1) = innerC (psi, phi);
    }
    return corr;
}

//---------------- More efficient way -------------------------
ITensor multSite (ITensor const& op, ITensor A)
{
    A *= op;
    A.noPrime();
    return A;
}

ITensor multSite (ITensor const& op1, ITensor const& op2, ITensor A)
{
    A = multSite (op2, A);
    A = multSite (op1, A);
    return A;
}

template <typename SiteType>
void MeasureCorrOne (const SiteType& sp, MPS& psi, int a, int b, CMatrix& corr, const string& op1, const string& op2)
// Update the correlation matrix of "shell" elements (a,b)...(b,b)...(b,a)
// spin_str can be "up" or "dn"
{
    if (a > b) {
        cout << "Error: MeasureCorrOne: a > b" << endl
             << "       a, b = " << a << ", " << b << endl;
        throw;
    }
    psi.position(b);

    ITensor A = multSite (sp.op(op1,b), sp.op(op2,b), psi(b));
    corr(b-a,b-a) = eltC (A * dag(psi(b)));

    if (a == b) return;

    // 
    ITensor R2 = multSite (sp.op(op2,b), psi(b));
    ITensor R1 = multSite (sp.op(op1,b), psi(b));
    R2 *= prime(dag(psi(b)), rightLinkIndex(psi,b-1));
    R1 *= prime(dag(psi(b)), rightLinkIndex(psi,b-1));
    for(int j = b-1; j >= a; j--)
    {
        auto L1 = multSite (sp.op(op1,j), sp.op("F",j), psi(j));
        auto L2 = multSite (sp.op("F",j), sp.op(op2,j), psi(j));
        L1 *= prime(dag(psi(j)), rightLinkIndex(psi,j));
        L2 *= prime(dag(psi(j)), rightLinkIndex(psi,j));
        corr(j-a,b-a) = eltC(L1*R2);
        corr(b-a,j-a) = eltC(L2*R1);
        R1 =  R1 * prime(dag(psi(j)),"Link") * multSite(psi(j), sp.op("F",j));
        R2 =  R2 * prime(dag(psi(j)),"Link") * multSite(psi(j), sp.op("F",j));
    }
}

template <typename SiteType>
CMatrix MeasureCorrBlock (const SiteType& sites, MPS& psi, int a, int b, const string& op1, const string& op2)
{
    auto corr = CMatrix (b-a+1, b-a+1);
    for(int i = a; i <= b; i++)
        MeasureCorrOne (sites, psi, a, i, corr, op1, op2);
    return corr;
}

template <typename SiteType>
CMatrix MeasureCorr (const SiteType& sites, MPS psi, const string& op1, const string& op2)
{
    return MeasureCorrBlock (sites, psi, 1, length(psi), op1, op2);
}
