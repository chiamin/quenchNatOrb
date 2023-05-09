#ifndef __ROTATEBASIS_H_CMC__
#define __ROTATEBASIS_H_CMC__
#include "itensor/all.h"
#include "FourierTransf.h"
using namespace itensor;
using namespace std;

void basis_analysis (const MPS& psi)
{
    auto sites = Fermion (iut::siteInds_byTag (psi,"Site"));
    auto corr = MeasureCorr (sites, psi, "Cdag", "C");

    Vector eigvals;
    CMatrix eigvecs;
    diagHermitian (corr, eigvecs, eigvals);

    Real tr_corr = sumels(eigvals);
    cout << "tot N = " << tr_corr << endl;

    vector<int> kn_obtained;
    int N = eigvals.size();
    ofstream ofs ("k.dat");
    ofs << scientific << setprecision(18) << endl;
    for(int i = 0; i < N; i++)
    {
        auto phi = CVector (column (eigvecs,i));
        auto kn = find_kn (phi, kn_obtained);
        if (kn > 1) kn -= 2;
        
        kn_obtained.push_back (kn);
        ofs << Fourier_k (kn,N) << " " << eigvals(i) << endl;
    }
}

void basis_analysis2 (const MPS& psi)
{
    auto sites = Fermion (iut::siteInds_byTag (psi,"Site"));
    auto corr = MeasureCorr (sites, psi, "Cdag", "C");
    using MatType = decltype(corr);

    int N = length(psi);
    for(int i = 2; i <= N; i++)
    {
        auto subcorr = MatType (subMatrix (corr,1,i,1,i));

        Vector eigvals;
        CMatrix eigvecs;
        diagHermitian (subcorr, eigvecs, eigvals);

        cout << "i = " << i << endl;
        cout << eigvals << endl;
    }
}

void project_to_product_orbitals (MPS& psi, const Vector& occs, Real zero_crit)
{
    MPS psi0 = psi;
    for(int i = 1; i <= occs.size(); i++)
    {
        auto ii = findIndex (psi(i), "Site");
        auto si = FermionSite (ii);
        auto occ = occs(i-1);
        if (abs(occ) < zero_crit)
        {
            auto proj = setElt (si.state("Emp"));
            psi.ref(i) *= dag(proj);
            psi.ref(i) *= proj;
            cout << "site " << i << ": " << occ << " -> 0" << endl;
        }
        else if (abs(occ-1.) < zero_crit)
        {
            auto proj = setElt (si.state("Occ"));
            psi.ref(i) *= dag(proj);
            psi.ref(i) *= proj;
            cout << "site " << i << ": " << occ << " -> 1" << endl;
        }
        else
        {
            cout << "site " << i << ": " << occ << endl;
        }
    }
    psi.orthogonalize();
    cout << "norm = " << norm(psi) << endl;    // norm should be close to 1
    cout << "overlap = " << innerC(psi0, psi) << endl;   // overlap should be close to 1
    psi.normalize();
}

vector<int> get_sort_index_by_k (const CMatrix& eigvecs)
{
    int N = ncols(eigvecs);
    vector<int> kn_obtained;
    vector<Real> abs_k;
    for(int i = 0; i < N; i++)
    {
        auto phi = CVector (column (eigvecs,i));
        auto kn = find_kn (phi, kn_obtained);
        kn_obtained.push_back (kn);
        auto k = Fourier_k (kn, N);
        abs_k.push_back (std::abs(k));
    }
    auto inds = iut::get_sort_indices (abs_k);
    return inds;
}

tuple<CMatrix,Vector> sort_orbitals (const CMatrix& orbitals, const Vector& occs)
{
    auto inds = get_sort_index_by_k (orbitals);

    int N = occs.size();
    CMatrix orbs_re (N,N);
    Vector occs_re (N);
    for(int i = 0; i < inds.size(); i++)
    {
        int j = inds.at(i);
        occs_re(i) = occs(j);
        column (orbs_re,i) &= column (orbitals,j);
    }
    return {orbs_re, occs_re};
}

template <typename MatType>
bool check_unitary (const MatType& U)
{
    auto I = conj(transpose(U)) * U;
    int N = ncols(U);
    MatType Iexact (N,N);
    for(int i = 0; i < N; i++)
        Iexact(i,i) = 1;
    auto zero = I - Iexact;
    if (norm(zero) > 1e-10)
    {
        cout << "check unitary matrix failed " << norm(zero) << endl;
        return false;
    }
    return true;
}

void check_unitary (const ITensor& A)
{
    if (order(A) != 2)
    {
        cout << "Error: check_unitary only supports rank-2 tensor";
        cout << A << endl;
        throw;
    }
    auto Adag = dag(A);
    Adag.prime(Adag.inds()(2));
    auto I = A * Adag;
    auto Iexact = iut::Identity(I.inds()(1), I.inds()(2));
    auto zero = I - Iexact;
    if (norm(zero) > 1e-9)
    {
        cout << "check hermitian tensor failed " << norm(zero) << endl;
        throw;
    }
}

inline CMatrix Rot_onebody (Cplx a, Cplx b)
// One-body rotational matrix
//     [ a ]   [ c ]
// R * [ b ] = [ 0 ]
{
    CMatrix rot(2,2);
    rot(0,0) = std::conj(a);
    rot(0,1) = std::conj(b);
    rot(1,0) = -b;
    rot(1,1) = a;
    Real c = std::sqrt (std::norm(a) + std::norm(b));
    rot *= 1./c;
    return rot;
}

// Make a many-body gate for a single-particle-basis rotation.
// <s1> and <s2> must be for spinless fermion.
// <rot> must be a 2x2 unitary matrix
template <typename MatType>
ITensor RotGate (const Index& s1, const Index& s2, const MatType& rot)
{
    assert (check_unitary(rot));
    assert (ncols(rot) == 2 and nrows(rot) == 2);

    Index sP1 = prime(s1);
    Index sP2 = prime(s2);
    ITensor res (sP1, sP2, dag(s1), dag(s2));

    // Basis:
    // |1,0>: particle at site 1
    // |0,1>: particle at site 2
    //
    // <rot> basis:
    // 0 -> |1,0>
    // 1 -> |0,1>
    // 
    // gate basis:
    // The first two indices are for rows;
    //     last                      columns.
    // 2,1 -> |1,0>
    // 1,2 -> |0,1>
    // 1,1 -> |0,0>
    // 2,2 -> |1,1>
    //
    // mapping from gate indices to <rot> indices:
    // 2,1 -> 0
    // 1,2 -> 1
    //
    res.set (1,1, 1,1, 1.0);
    res.set (2,1, 2,1, rot(0,0));
    res.set (2,1, 1,2, rot(0,1));
    res.set (1,2, 2,1, rot(1,0));
    res.set (1,2, 1,2, rot(1,1));
    res.set (2,2, 2,2, 1.0);

    return res;
}

// Singe-site phase operator for a particular single-particle orbital
template <typename NumType>
ITensor PhaseGate (const Index& si, NumType phase)
{
    // Check <phase> is a phase
    assert (std::abs(norm(phase)) - 1. < 1e-10);

    ITensor res (dag(si), prime(si));
    res.set (1,1,1.0);
    res.set (2,2,phase);
    return res;
}

vector<CMatrix> Get_rots (CVector orbit, int ibeg=1)
// Rotation gates from the last gate to the first gate
// Suppose the dimension of the orbital is N
// The orbital elements are assumed to be zero for i < ibeg, so we consider only a subvector from ibeg to N.
// The ortibal after rotation will become
//
//   [ 0 ]
//   | 0 |   --> 0 by assumption
//   | : |
//   | 1 |   --> ibeg; 1 by rotate 
//   | 0 |
//   | : |   --> 0 by rotations
//   [ 0 ]
//
// In this explanation, the index starts from 1,
// but in Vector object, the index starts from 0.
{
    vector<CMatrix> rots;

    // If the last site, only the phase needed to be fixed.
    // Store the last matrix element to a 1x1 matrix
    if (ibeg == orbit.size())
    {
        CMatrix rot (1,1);
        rot(0,0) = iut::conj(orbit(ibeg-1));
        rots.push_back (rot);
        return rots;
    }

    // Otherwise, store the rotation matrices
    int jend = orbit.size()-1;  // starting from the ending 2x2 block
    for(int j = jend; j >= ibeg; j--)
    {
        Cplx a = orbit(j-1),
             b = orbit(j);
        CMatrix rot = Rot_onebody (a, b);
        rots.push_back (rot);

        subVector (orbit,j-1,j+1) &= rot * subVector (orbit,j-1,j+1);
    }
    return rots;
}

void print_psi (const MPS& psi)
{
    ITensor A(1.);
    for(int i = 1; i <= length(psi); i++)
        A *= psi(i);
    //PrintData(A);
    cout << "110  " << eltC(A,2,2,1) << endl;
    cout << "101  " << eltC(A,2,1,2) << endl;
    cout << "011  " << eltC(A,1,2,2) << endl;
}

void print_psi4 (const MPS& psi)
{
    ITensor A(1.);
    for(int i = 1; i <= length(psi); i++)
        A *= psi(i);
    //PrintData(A);
    cout << "1100  " << eltC(A,2,2,1,1) << endl;
    cout << "1010  " << eltC(A,2,1,2,1) << endl;
    cout << "1001  " << eltC(A,2,1,1,2) << endl;
    cout << "0110  " << eltC(A,1,2,2,1) << endl;
    cout << "0101  " << eltC(A,1,2,1,2) << endl;
    cout << "0011  " << eltC(A,1,1,2,2) << endl;
}

// Apply the gates from the end to the front
// Suppose the length of MPS is N
// For two-site gates, the 1st gate will be applied on (N-1,N),
//                         2nd                         (N-2,N-1),
// and so on.
void apply_rotation (MPS& psi, const vector<CMatrix>& rots, Args const& args)
{
    int Nrot = rots.size();
    int N = length(psi);
    for(int ir = 0; ir < Nrot; ir++)
    {
        auto rot = rots.at(ir);
        // one-site gate
        if (ncols(rot) == 1)
        {
            // One-site gate can happen only on the last site
            assert (ir == 0 and rots.size() == 1);

            // Save the orthogonality center so that we can reset it later
            auto leftlim = leftLim( psi);
            auto rightlim = rightLim (psi);

            // Apply the gate
            int i = N - ir;
            auto is = findIndex (psi(N), "Site");
            auto rot_gate = PhaseGate (is, rot(0,0));

            psi.ref(N) *= rot_gate;
            psi.ref(N).noPrime("Site");

            // Restore the orthogonality center,
            // since a single-site unitary gate will not change the orthogonality center
            psi.leftLim (leftlim);
            psi.rightLim (rightlim);
        }
        // two-site gates
        else
        {
            int i = N-1 - ir;
            psi.position(i);
            auto is1 = findIndex (psi(i), "Site");
            auto is2 = findIndex (psi(i+1), "Site");
            auto rot_gate = RotGate (is1, is2, rot);
            applyGate (rot_gate, psi, {args,"Fromleft",false});
        }
    }
}

void apply_rotation (CMatrix& orbits, const vector<CMatrix>& rots)
{
    int N = ncols (orbits);
    for(int ir = 0; ir < rots.size(); ir++)
    {
        CMatrix rot = rots.at(ir);
        int rsize = ncols(rot);         // size of the gate; can be 1 or 2
        int j = N - rsize - ir;         // for ir=0, it should apply on the final 2 sites or 1 site
        subMatrix (orbits,j,j+rsize,0,N) &= rot * CMatrix(subMatrix (orbits,j,j+rsize,0,N));
    }
}

// Rotate the basis of <psi> based on its correlation matrix
// <psi> will be changed after the function call
template <typename SiteType>
void rotate_basis (const SiteType& sites, MPS& psi, CMatrix orbs)
{
    int N = ncols(orbs);

    // Rotate to the orbital basis one by one
    for(int i = 1; i <= N; i++)
    {
        // Get the ith orbital
        auto orb = CVector (column (orbs,i-1));
        // Get the corresponding rotation matrices
        auto rots = Get_rots (orb, i);
        // Apply the rotations to psi
        apply_rotation (psi, rots, {"Cutoff",1e-10});
        // Apply the rotations to the orbitals, including the ones that are not yet rotated
        apply_rotation (orbs, rots);
        // Orthogonalize psi
        psi.orthogonalize({"Cutoff",1e-10});
        cout << "max dim " << i << " = " << maxLinkDim(psi) << endl;
    }

    // Test: <orbs> after the rotation, orbs should be an identity matrix
    auto test = orbs;
    for(int i = 0; i < N; i++)
        test(i,i) -= 1.;
    if (norm(test) > 1e-10)
    {
        cout << "Error: rotation failed" << endl;
        cout << orbs << endl;
        throw;
    }

    // Test: Check that <psi> is normalized after rotation
    cout << "test norm: " << norm(psi) << endl;
}
#endif
