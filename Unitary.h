#ifndef __UNITARY_H_CMC__
#define __UNITARY_H_CMC__
#include "itensor/all.h"
using namespace itensor;
using namespace std;
/*
MPO make_zeroMPO (const Fermion& sites)
{
    int N = length(sites);
    MPO p (N);
    Index ilinkpre;
    for(int k = 1; k <= N; k++)
    {
        auto ii = sites(k);
        ITensor A (dag(ii));
        ITensor Ap (prime(ii));

        A.set (1, 1.);
        Ap.set (1, 1.);

        auto ilink = Index(QN({"Pf",0,-2}),1,Out,"Link");
        if (k == 1)
            p.ref(k) = (A * Ap) * setElt(ilink=1);
        else if (k == N)
            p.ref(k) = (A * Ap) * setElt(ilinkpre=1);
        else
            p.ref(k) = (A * Ap) * setElt(ilinkpre=1) * setElt(ilink=1);
        ilinkpre = dag(ilink);
    }
    return p;
}

MPO projectMPO (const Fermion& sites, int i, int j, MPO zeroMPO)
{
    if (i == j)
    {
        auto ii = sites(i);
        auto A = ITensor (zeroMPO(i).inds());
        auto ilinks = linkInds (zeroMPO, i);
        if (order(ilinks) == 1)
            A.set (ii=2, prime(ii)=2, ilinks(1)=1, 1.);
        else
            A.set (ii=2, prime(ii)=2, ilinks(1)=1, ilinks(2)=1, 1.);
        zeroMPO.ref(i) = A;
    }
    else
    {
        auto ii = sites(i);
        auto A = ITensor (zeroMPO(i).inds());
        auto ilinks = linkInds (zeroMPO, i);
        if (order(ilinks) == 1)
            A.set (ii=2, prime(ii)=1, ilinks(1)=1, 1.);
        else
            A.set (ii=2, prime(ii)=1, ilinks(1)=1, ilinks(2)=1, 1.);
        zeroMPO.ref(i) = A;

        ii = sites(j);
        A = ITensor (zeroMPO(j).inds());
        ilinks = linkInds (zeroMPO, j);
        if (order(ilinks) == 1)
            A.set (ii=1, prime(ii)=2, ilinks(1)=1, 1.);
        else
            A.set (ii=1, prime(ii)=2, ilinks(1)=1, ilinks(2)=1, 1.);
        zeroMPO.ref(j) = A;
    }
    return zeroMPO;
}

template <typename VecType>
MPO projector1_old (const Fermion& sites, const VecType& phi, int i, const MPO& zeroMPO)
{
    int N = length(sites);
    vector<MPO> ps (N);
    for(int j = 1; j <= N; j++)
    {
        auto pj = projectMPO (sites, i, j, zeroMPO);
        ps.at(j-1) = phi(j-1) * pj;
    }
    auto p = sum (ps, {"Cutoff",1E-12});
    p.orthogonalize();
    return p;
}

MPO projectMPO (const Fermion& sites, int i, int j)
{
    int N = length(sites);
    MPO p (sites);
    Index ilinkpre;
    for(int k = 1; k <= N; k++)
    {
        auto ii = sites(k);
        ITensor A (dag(ii));
        ITensor Ap (prime(ii));

        int st = (k == i ? 2 : 1);
        int stp = (k == j ? 2 : 1);
        A.set (st, 1.);
        Ap.set (stp, 1.);

        auto ilink = Index(QN({"Pf",0,-2}),1,Out,"Link");
        if (k == 1)
            p.ref(k) = (A * Ap) * setElt(ilink=1);
        else if (k == N)
            p.ref(k) = (A * Ap) * setElt(ilinkpre=1);
        else
            p.ref(k) = (A * Ap) * setElt(ilinkpre=1) * setElt(ilink=1);
        ilinkpre = dag(ilink);
    }
    p.orthogonalize();
    return p;
}

template <typename VecType>
MPO projector1_old2 (const Fermion& sites, const VecType& phi, int i)
{
    MPO p;
    int N = length(sites);
    for(int j = 1; j <= N; j++)
    {
        auto pj = projectMPO (sites, i, j);
        pj *= phi(j-1);
        if (!p)
            p = pj;
        else
            p = sum (p, pj, {"Cutoff",1E-12});
    }
    p.orthogonalize();
    return p;
}

*/
MPO projector1_tmp (const Fermion& sites)
{
    int N = length(sites);
    MPO p (N);
    Index ilinkpre;
    for(int k = 1; k <= N; k++)
    {
        // Define empty tensor
        auto ii = sites(k);
        ITensor A;
        auto ilink = Index(QN({"Pf",0,-2}),1,QN({"Pf",1,-2}),1,Out,"Link");
        if (k == 1)
            A = ITensor (dag(ii), prime(ii), ilink);
        else if (k == N)
            A = ITensor (dag(ii), prime(ii), ilinkpre);
        else
            A = ITensor (dag(ii), prime(ii), ilink, ilinkpre);
        ilinkpre = dag(ilink);

        int iv = 1; // Set ket states all 0
        if (k == 1)
        {
            A.set (iv,1,1, 1.);
        }
        else if (k == N)
        {
            A.set (iv,1,2, 1.);
        }
        else
        {
            A.set (iv,1,1,1, 1.);
            A.set (iv,1,2,2, 1.);
        }
        p.ref(k) = A;
    }
    return p;
}

template <typename VecType>
MPO projector1_from_tmp (const VecType& phi, int i, MPO tmpMPO)
{
    int N = length(tmpMPO);

    tmpMPO.ref(i).store() = 0;
    if (i == 1)
    {
        tmpMPO.ref(i).set (2,1,1, 1.);
    }
    else if (i == N)
    {
        tmpMPO.ref(i).set (2,1,2, 1.);
    }
    else
    {
        tmpMPO.ref(i).set (2,1,1,1, 1.);
        tmpMPO.ref(i).set (2,1,2,2, 1.);
    }

    for(int k = 1; k <= N; k++)
    {
        auto coef = phi(k-1);
        int iv = (k == i ? 2 : 1);
        if (k == 1)
        {
            tmpMPO.ref(k).set (iv,2,2, coef);
        }
        else if (k == N)
        {
            tmpMPO.ref(k).set (iv,2,1, coef);
        }
        else
        {
            tmpMPO.ref(k).set (iv,2,2,1, coef);
        }
    }
    return tmpMPO;
}


MPO projector_out1 (const Fermion& sites, int i)
{
    int N = length(sites);
    MPO p (sites);  // identity
    p.ref(i).store() = 0;   // clear storage

    // |0> -> |0>
    // |1> -> |0>
    if (i == 1 or i == length(sites))
    {
        p.ref(i).set (1,1,1, 1.);
        p.ref(i).set (2,1,1, 1.);
    }
    else
    {
        p.ref(i).set (1,1,1,1, 1.);
        p.ref(i).set (2,1,1,1, 1.);
    }
    return p;
}

template <typename VecType>
MPO projector1 (const Fermion& sites, const VecType& phi, int i)
{
    int N = length(sites);
    MPO p (N);
    Index ilinkpre;
    for(int k = 1; k <= N; k++)
    {
        // Define empty tensor
        auto ii = sites(k);
        ITensor A;
        auto ilink = Index(QN({"Pf",0,-2}),1,QN({"Pf",1,-2}),1,Out,"Link");
        if (k == 1)
            A = ITensor (dag(ii), prime(ii), ilink);
        else if (k == N)
            A = ITensor (dag(ii), prime(ii), ilinkpre);
        else
            A = ITensor (dag(ii), prime(ii), ilink, ilinkpre);
        ilinkpre = dag(ilink);

        auto coef = phi(k-1);
        int iv = (k == i ? 2 : 1);
        if (k == 1)
        {
            A.set (iv,1,1, 1.);
            A.set (iv,2,2, coef);
        }
        else if (k == N)
        {
            A.set (iv,1,2, 1.);
            A.set (iv,2,1, coef);
        }
        else
        {
            A.set (iv,1,1,1, 1.);
            A.set (iv,2,2,1, coef);
            A.set (iv,1,2,2, 1.);
        }
        p.ref(k) = A;
    }
    return p;
}


template <typename MatType>
MPO unitaryMPO (const Fermion& sites, const MatType& U)
{
    auto tmpMPO = projector1_tmp (sites);

    int N = length(sites);
    vector<MPO> ps (N);
    for(int i = 1; i <= N; i++)
    {
        auto pi = projector1_from_tmp (column (U,i-1), i, tmpMPO);
    cout << "dim = " << maxLinkDim(pi) << endl;
        ps.at(i-1) = pi;
    }
cout << "sum" << endl;
    auto p = sum (ps, {"Cutoff",1E-12});
cout << "sum end" << endl;
    p.orthogonalize();
    return p;
}

template <typename MatType>
MPO unitaryMPO_backup (const Fermion& sites, const MatType& U)
{
    MPO p;
    int N = length(sites);
    for(int i = 1; i <= N; i++)
    {
        auto pi = projector1 (sites, column (U,i-1), i);
        if (!p)
            p = pi;
        else
            p = sum (p, pi, {"Cutoff",1E-12});
cout << "dim = " << maxLinkDim(pi) << endl;
    }
    p.orthogonalize();
    return p;
}

void printbig (const MPS& mps)
{
    ITensor A = mps(1);
    for(int i = 2; i <= length(mps); i++)
        A *= mps(i);
    PrintData(A);
}

void testk (const Fermion& sites)
{
    auto pout = projector_out1 (sites, 1);

    InitState st (sites);
    for(int i = 1; i <= length(sites); i++)
        st.set (i,"Occ");
    auto mps = MPS (st);

    auto phi = applyMPO (pout, mps);
    printbig (phi);
}


#endif
