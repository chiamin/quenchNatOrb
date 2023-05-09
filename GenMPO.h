#ifndef __GENMPO_H_CMC__
#define __GENMPO_H_CMC__
#include "itensor/mps/mpo.h"
using namespace itensor;
using namespace std;

tuple <AutoMPO,int,int>
t_mu_V_ampo
(const SiteSet& sites, int L_lead, int L_device,
 Real t_leadL, Real t_leadR, Real t_device, Real t_contactL, Real t_contactR,
 Real mu_leadL, Real mu_leadR, const vector<Real>& mu_device,
 Real V_leadL, Real V_leadR, Real V_device, Real V_contactL, Real V_contactR)
{
    auto N = length (sites);
    int first_site_device = L_lead + 1;
    int last_site_device = L_lead + L_device;
    if (N != 2*L_lead+L_device)
    {
        cout << "Error: size not match" << endl;
        cout << "N, L_lead, L_device = " << N << " " << L_lead << " " << L_device << endl;
        throw;
    }
    AutoMPO ampo (sites);
    for(int i = 1; i <= N; ++i)
    {
        Real mui, ti, Vi;
        // Set mu
        if      (i < first_site_device)  mui = mu_leadL;
        else if (i > last_site_device)   mui = mu_leadR;
        else
        {
            mui = mu_device.at(i-first_site_device);
        }
        // Set t and V
        if (i == first_site_device-1)
        {
            ti = t_contactL;
            Vi = V_contactL;
        }
        else if (i == last_site_device)
        {
            ti = t_contactR;
            Vi = V_contactR;
        }
        else if (i < first_site_device-1)
        {
            ti = t_leadL;
            Vi = V_leadL;
        }
        else if (i > last_site_device)
        {
            ti = t_leadR;
            Vi = V_leadR;
        }
        else
        {
            ti = t_device;
            Vi = V_device;
        }

        if (ti != 0. and i != N)
        {
            ampo += -ti,"Cdag",i,"C",i+1;
            ampo += -ti,"Cdag",i+1,"C",i;
//cout << "t: " << i << " " << i+1 << " " << ti << endl;
        }
        if (mui != 0.)
        {
            ampo += -mui,"N",i;
//cout << "mu: " << i << " " << mui << endl;
        }
        if (Vi != 0. and i != N)
        {
            ampo += Vi,"N",i,"N",i+1;
//cout << "V: " << i << " " << i+1 << " " << Vi << endl;
        }
    }
    return {ampo, first_site_device, last_site_device};
}

tuple <AutoMPO,int,int>
t_mu_V_ampo
(const SiteSet& sites, int L_lead, int L_device,
 Real t_leadL, Real t_leadR, Real t_device, Real t_contactL, Real t_contactR,
 Real mu_leadL, Real mu_leadR, Real mu_device,
 Real V_leadL, Real V_leadR, Real V_device, Real V_contactL, Real V_contactR)
{
    vector<Real> mus (L_device, mu_device);
    return t_mu_V_ampo (sites, L_lead, L_device, t_leadL, t_leadR, t_device, t_contactL, t_contactR,
                        mu_leadL, mu_leadR, mus, V_leadL, V_leadR, V_device, V_contactL, V_contactR);
}

Matrix get_tij
(int L_lead, int L_device,
 Real t_leadL, Real t_leadR, Real t_device, Real t_contactL, Real t_contactR,
 Real mu_leadL, Real mu_leadR, Real mu_device)
{
    auto N = 2*L_lead+L_device;
    int first_site_device = L_lead + 1;
    int last_site_device = L_lead + L_device;

    Matrix tij (N,N);
    for(int i = 1; i <= N; ++i)
    {
        Real mui, ti, Vi;
        // Set mu
        if (i < first_site_device)
            mui = mu_leadL;
        else if (i > last_site_device)
            mui = mu_leadR;
        else
            mui = mu_device;
        tij(i-1,i-1) = mui;

        // Set t
        if (i != N)
        {
            if (i == first_site_device-1)
                ti = t_contactL;
            else if (i == last_site_device)
                ti = t_contactR;
            else if (i < first_site_device-1)
                ti = t_leadL;
            else if (i > last_site_device)
                ti = t_leadR;
            else
                ti = t_device;
            tij(i-1,i) = ti;
            tij(i,i-1) = ti;
        }
    }
    return tij;
}

template <typename MatrixType>
void add_t (AutoMPO& ampo, const MatrixType& tij)
{
    int N = ncols(tij);
    for(int i = 1; i <= N; i++)
    {
        ampo += -tij(i-1,i-1),"N",i;
        for(int j = 1; j < i; j++)
        {
            ampo += -tij(i-1,j-1),"Cdag",i,"C",j;
            ampo += -iut::conj(tij(i-1,j-1)),"Cdag",j,"C",i;
        }
    }
}

template <typename MatrixType>
void add_t_rot (AutoMPO& ampo, const Matrix& tij, const MatrixType& U)
{
    auto tij_rot = conj(transpose(U)) * tij * U;
    add_t (ampo, tij_rot);
}

void add_tij (AutoMPO& ampo, int L_lead, int L_device,
              Real t_leadL, Real t_leadR, Real t_device, Real t_contactL, Real t_contactR,
              Real mu_leadL, Real mu_leadR, Real mu_device, const CMatrix& U)
{
    auto tij = get_tij (L_lead, L_device, t_leadL, t_leadR, t_device, t_contactL, t_contactR, mu_leadL, mu_leadR, mu_device);
    add_t_rot (ampo, tij, U);
}

Matrix get_Vij (int L_lead, int L_device, Real V_leadL, Real V_leadR, Real V_device, Real V_contactL, Real V_contactR)
{
    auto N = 2*L_lead+L_device;
    int first_site_device = L_lead + 1;
    int last_site_device = L_lead + L_device;

    Matrix V (N,N);
    for(int i = 1; i < N; i++)
    {
        // Set V
        if (i == first_site_device-1)
            V(i-1,i) = V_contactL;
        else if (i == last_site_device)
            V(i-1,i) = V_contactR;
        else if (i < first_site_device-1)
            V(i-1,i) = V_leadL;
        else if (i > last_site_device)
            V(i-1,i) = V_leadR;
        else
            V(i-1,i) = V_device;
    }
    return V;
}

tuple<ITensor,ITensor,ITensor> get_decomposed_Vij (int L_lead, int L_device, Real V_leadL, Real V_leadR, Real V_device, Real V_contactL, Real V_contactR)
{
    auto Vij = get_Vij (L_lead, L_device, V_leadL, V_leadR, V_device, V_contactL, V_contactR);

    int N = ncols(Vij);
    Index ii(N), kk(N);
    auto V = matrixITensor (Vij, ii, prime(ii));
    auto[uu,S,vv] = svd(V,{ii},{"Cutoff",1e-12});
    return {uu, S, vv};
}

MPO V_MPO_rot (const Fermion& sites, const ITensor& u_V, const ITensor& s_V, const ITensor& v_V, const CMatrix U)
{
    auto ii = noPrime (uniqueIndex(u_V, s_V));
    auto kk = Index (dim(ii));
    auto UU = matrixITensor (U, ii, prime(kk));
    auto Udag = noPrime(dag(UU));

    auto iu = commonIndex (u_V, s_V);
    auto iv = commonIndex (v_V, s_V);
    int N = length(sites);
    int D = dim(s_V.inds()(1));
    vector<MPO> mpos;
    mpos.reserve (D);

    for(int i = 1; i <= D; i++)
    {
        auto s = elt (s_V, i, i);
        auto u = u_V * setElt(iu=i);
        auto v = v_V * setElt(iv=i);

        cout << i << ": V singular = " << s << endl;

        u /= UU;
        u *= Udag;
        v /= prime(UU,ii);
        v *= prime(Udag,ii);

// Check hermitian
for(int n = 1; n <= dim(u.inds()(1)); n++)
for(int m = 1; m <= dim(u.inds()(1)); m++)
{
    if (std::abs(eltC(u,n,m)-iut::conj(eltC(u,m,n))) > 1e-12 or std::abs(eltC(v,n,m)-iut::conj(eltC(v,m,n))) > 1e-12)
    {
        cout << "gg" << endl;
        cout << eltC(u,n,m) << " " << eltC(u,m,n) << endl;
        cout << eltC(v,n,m) << " " << eltC(v,m,n) << endl;
        exit(0);
    }
}


        AutoMPO ampo1 (sites),
                ampo2 (sites);
        cout << "construct tmp mpo" << endl;
        for(int j = 1; j <= N; j++)
        {
            // Diagonal terms
            auto coef1 = eltC (u, kk=j, prime(kk)=j);
            ampo1 += coef1,"N",j;

            auto coef2 = eltC (v, kk=j, prime(kk)=j);
            ampo2 += coef2,"N",j;

            // Off-diagonal temrs
            for(int k = 1; k < j; k++)
            {
                coef1 = eltC (u, kk=j, prime(kk)=k);
                if (std::abs(coef1) > 1e-12)
                {
                    ampo1 += coef1,"Cdag",j,"C",k;
                    ampo1 += iut::conj(coef1),"Cdag",k,"C",j;
                }

                coef2 = eltC (v, kk=j, prime(kk)=k);
                if (std::abs(coef2) > 1e-12)
                {
                    ampo2 += coef2,"Cdag",j,"C",k;
                    ampo2 += iut::conj(coef2),"Cdag",k,"C",j;
                }
            }
        }
        auto mpo1 = toMPO (ampo1);
        auto mpo2 = toMPO (ampo2);

        cout << "multiply mpo" << endl;
        auto mpo = s * nmultMPO (prime(mpo1), mpo2);

        mpo.mapPrime(2,1);
        mpos.push_back (mpo);
    }
    cout << "sum mpo" << endl;
    auto MPO = sum (mpos);

    return MPO;
}

// Return the coefficient of V_k1_k2_k3_k4 for C_k1^dagger C_k2^dagger C_k3 C_k4
Cplx Vkkkk (const Matrix& Vij, const CMatrix U, int k1, int k2, int k3, int k4)
{
    int N = ncols(Vij);
    Index ii(N), kk(N);
    Index kkpr = prime(kk),
          kkpr2 = prime(kk,2),
          kkpr3 = prime(kk,3);
    auto V = matrixITensor (Vij, ii, prime(ii));

    auto UU = matrixITensor (U, ii, kk);
    auto Udag = dag(UU);

    auto U1 = Udag;
    auto U2 = prime(UU,1,kk);
    auto U3 = prime(Udag,2,kk);
    U3.prime(ii);
    auto U4 = prime(UU,3,kk);
    U4.prime(ii);

/*cout << U1 << endl << U2 << endl;
cout << U1/U2 << endl;
check_unitary (U1/U2);
check_unitary (U3/U4);*/


    U1 *= setElt(kk=k1);
    U2 *= setElt(kkpr=k2);
    U3 *= setElt(kkpr2=k3);
    U4 *= setElt(kkpr3=k4);

    auto Uip = U1 / U2;
    auto Ui = U3 / U4;
    auto C = Uip * V * Ui;

    return eltC(C);
}

void add_V (AutoMPO& ampo, const Matrix& Vij, const CMatrix& U)
{
    int N = ncols(Vij);
    for(int k1 = 1; k1 <= N; k1++)
    for(int k2 = 1; k2 <= N; k2++)
    for(int k3 = 1; k3 <= N; k3++)
    for(int k4 = 1; k4 <= N; k4++)
    {
        auto c = Vkkkk (Vij, U, k1, k2, k3, k4);
        if (true)//abs(c) > 1e-12)
        {
            //cout << "*** " << k1 << " " << k2 << " " << k3 << " " << k4 << ": " << c << endl;
            ampo += c,"Cdag",k1,"C",k2,"Cdag",k3,"C",k4;
        }
    }
}

void add_VNN (AutoMPO& ampo, int L_lead, int L_device, Real V_leadL, Real V_leadR, Real V_device, Real V_contactL, Real V_contactR, const CMatrix& U)
{
    auto Vij = get_Vij (L_lead, L_device, V_leadL, V_leadR, V_device, V_contactL, V_contactR);
cout << Vij << endl;
    add_V (ampo, Vij, U);
}
#endif
