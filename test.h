#ifndef __TEST_H_CMC__
#define __TEST_H_CMC__
#include "Unitary.h"
template <typename VecType>
void test_project1 (const Fermion& sites, const VecType& eigvecs)
{
    auto phi = CVector (column (eigvecs,0));
    auto tmpMPO = projector1_tmp (sites);
    auto p = projector1_from_tmp (phi, 1, tmpMPO);

    InitState st (sites);
    st.set (1, "Occ");
    auto prod_state = MPS (st);

    // ---------------------------------------
    auto rot_state = applyMPO (p, prod_state);
    ITensor A = rot_state(1);
    for(int i = 2; i <= length(prod_state); i++)
        A *= rot_state(i);
    PrintData(A);

    cout << phi << endl;
}

template <typename VecType>
void test_uMPO (const Fermion& sites, const VecType& eigvecs)
{
    InitState st (sites);
    st.set (2, "Occ");
    auto prod_state = MPS (st);

    auto phi = CVector (column (eigvecs,1));
    auto uMPO = unitaryMPO (sites, eigvecs);


    // ---------------------------------------
    auto rot_state = applyMPO (uMPO, prod_state);
    ITensor A = rot_state(1);
    for(int i = 2; i <= length(prod_state); i++)
        A *= rot_state(i);
    PrintData(A);

    cout << phi << endl;

    // --------------------------------------
    rot_state = applyMPO (dag(uMPO), rot_state);

    A = rot_state(1);
    for(int i = 2; i <= length(prod_state); i++)
        A *= rot_state(i);
    PrintData(A);
}
#endif
