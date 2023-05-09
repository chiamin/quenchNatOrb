#include "itensor/all.h"
#include "ReadInput.h"
#include "IUtility.h"
#include "MyObserver.h"
#include "MPSUtility.h"
#include "tdvp.h"
#include "TDVPObserver.h"
#include "Corr.h"
#include "GenMPO.h"
#include "RotateBasis.h"
#include "Unitary.h"
#include "test.h"
using namespace itensor;
using namespace std;

MPS make_initstate (const Fermion& sites, int Np)
{
    int N = length(sites);
    int Npar = Np;
    InitState init (sites);

    for(int i = 1; i <= N; i++)
    {
        string state;
        if (i % 2 == 0 && Np-- > 0)
            state = "Occ";
        else
            state = "Emp";
        init.set (i, state);
        cout << i << ": " << state << endl;
    }
    if (Np > 0)
    {
        for(int i = 1; i <= N; i += 2)
            if (Np-- > 0)
                init.set (i,"Occ");
    }
    auto psi = MPS (init);
    return psi;
}

AutoMPO set_H (const SiteSet& sites, int L_lead, int L_device, Real t, Real mu, Real V, bool periodic=false)
{
    int N = length(sites);
    AutoMPO ampo (sites);
    for(int i = 1; i <= N; ++i)
    {
        if (i != N)
        {
            ampo += -t,"Cdag",i,"C",i+1;
            ampo += -t,"Cdag",i+1,"C",i;
            cout << "H t " << i << " " << t << endl;
        }
        if (mu != 0.)
            ampo += -mu,"N",i;
    }
    if (periodic)
    {
        ampo += -t,"Cdag",N,"C",1;
        ampo += -t,"Cdag",1,"C",N;
        cout << "H t " << N << " " << t << endl;
    }

    int Ri = L_lead + L_device;
    for(int i = L_lead+1; i <= Ri; i++)
    for(int j = L_lead+1; j < i; j++)
    {
        ampo += V,"N",i,"N",j;
    }
    return ampo;
}

template <typename MPSType>
ITensor print_wf (const MPSType& psi)
{
    ITensor pp (1.);
    vector<Index> iis;
    for(int i = 1; i <= length(psi); i++)
    {
        pp *= psi(i);
        auto is = findIndex (psi(i), "Site,0");
        iis.push_back (is);
        if constexpr (is_same_v <MPO, MPSType>)
            iis.push_back (prime(is));
    }
    pp.permute (iis);
    PrintData(pp);
    return pp;
}

MPS test_rotate (const Fermion& sites, MPS psi, const CMatrix& U, const Vector& occs, Real zero_crit)
{
    MPS psi0 = psi;

/*    CMatrix U;
    Vector occs;

    // Compute the correlation matrix
    CMatrix corr = MeasureCorr (sites, psi, "Cdag", "C");
    diagHermitian (corr, U, occs);*/

    cout << "occupasion" << endl;
    for(int i = 0; i < occs.size(); i++)
        cout << i << " " << occs(i) << endl;

    // Rotate the basis
    rotate_basis (sites, psi, U);
    // For each site, if the occupation number ~ 0 or 1, project <psi> to a |0> or |1> state
    project_to_product_orbitals (psi, occs, zero_crit);
    // Check bond dimension
    cout << "bond dimensions: " << "site original rotated" << endl;
    for(int i = 1; i < length(psi); i++)
        cout << i << ": " << dim(iut::rightIndex(psi0, i)) << " " << dim(iut::rightIndex(psi, i)) << endl;

    return psi;
}

void check_hermitian (const MPO& mpo)
{
    ITensor A (1.);
    for(int i = 1; i <= length(mpo); i++)
        A *= mpo(i);
    auto Adag = dag(A);
    Adag.prime();
    Adag.mapPrime(2,0);
    auto zero = A - Adag;
    if (norm(zero) > 1e-9)
    {
        cout << "check hermitian failed. " << norm(zero) << endl;
        throw;
    }
}

// Test the energies are the same before and after rotation
void test_rotate_energy (const Fermion& sites, MPS psi, MPO H, Real zero_crit,
                         int L_lead, int L_device,
                         Real t_lead, Real t_device, Real t_contactL, Real t_contactR,
                         Real mu_leadL, Real mu_leadR, Real mu_device,
                         Real V_lead, Real V_device, Real V_contact)
{
    // Compute the correlation matrix and natural orbitals
    CMatrix U;
    Vector occs;
    CMatrix corr = MeasureCorr (sites, psi, "Cdag", "C");

    diagHermitian (corr, U, occs);

    // Rotate psi
    auto psi_rot = test_rotate (sites, psi, U, occs, zero_crit);

    // Backup H
    auto H0 = H;

    // Interaction
    auto[u,s,v] = get_decomposed_Vij (L_lead, L_device, V_lead, V_lead, V_device, V_contact, V_contact);
    auto HV = V_MPO_rot (sites, u, s, v, U);
    cout << "HV bond dim = " << maxLinkDim(HV) << endl;

    // Kinetic energy
    AutoMPO ampor (sites);
    add_tij (ampor, L_lead, L_device,
             t_lead, t_lead, t_device, t_contactL, t_contactR,
             mu_leadL, mu_leadR, mu_device, U);
    auto Hk = toMPO(ampor);
    cout << "Hk bond dim = " << maxLinkDim(Hk) << endl;

    // Sum the two parts
    H = sum (Hk, HV);

    H.orthogonalize({"Cutoff",1e-12});
    cout << "MPO bond dim after rotation = " << maxLinkDim(H) << endl;

    // Print energies
    cout << "E " << innerC(psi,H0,psi) << " " << innerC(psi_rot,H,psi_rot) << endl;
}

MPS psi_rotate (const Fermion& sites, MPS psi, const CMatrix& U, const Vector& occs, Real zero_crit)
{
    cout << "occupasion" << endl;
    for(int i = 0; i < occs.size(); i++)
        cout << i << " " << occs(i) << endl;

    // Rotate the basis
    rotate_basis (sites, psi, U);
    // For each site, if the occupation number ~ 0 or 1, project <psi> to a |0> or |1> state
    project_to_product_orbitals (psi, occs, zero_crit);
    // Check bond dimension
    cout << "bond dimensions: " << endl;
    for(int i = 1; i < length(psi); i++)
        cout << i << ": " << dim(iut::rightIndex(psi, i)) << endl;

    return psi;
}

MPO H_rot (const Fermion& sites, const CMatrix& U,
           int L_lead, int L_device,
           Real t_lead, Real t_device, Real t_contactL, Real t_contactR,
           Real mu_leadL, Real mu_leadR, Real mu_device,
           Real V_lead, Real V_device, Real V_contact)
{
    // Interaction
    auto[u,s,v] = get_decomposed_Vij (L_lead, L_device, V_lead, V_lead, V_device, V_contact, V_contact);
    auto HV = V_MPO_rot (sites, u, s, v, U);
    cout << "HV bond dim = " << maxLinkDim(HV) << endl;

    // Kinetic energy
    AutoMPO ampor (sites);
    add_tij (ampor, L_lead, L_device,
             t_lead, t_lead, t_device, t_contactL, t_contactR,
             mu_leadL, mu_leadR, mu_device, U);
    auto Hk = toMPO(ampor);
    cout << "Hk bond dim = " << maxLinkDim(Hk) << endl;

    // Sum the two parts
    auto H = sum (Hk, HV);

    H.orthogonalize({"Cutoff",1e-12});
    cout << "MPO bond dim after rotation = " << maxLinkDim(H) << endl;
    return H;
}


CMatrix random_U (int N)
{
    auto cc = randomMatC (N,N);
    auto corr = cc + conj(transpose(cc));
    CMatrix U;
    Vector occs;
    diagHermitian (corr, U, occs);
    return U;
}

void print_occ (const Fermion& sites, const MPS& psi)
{
    // Compute the correlation matrix and natural orbitals
    CMatrix U;
    Vector occs;
    CMatrix corr = MeasureCorr (sites, psi, "Cdag", "C");
    diagHermitian (corr, U, occs);
    for(int i = 0; i < length(psi); i++)
        cout << occs(i) << endl;
}

void test (const Fermion& sites, MPS psi)
{
    int N = length(sites);
    auto U1 = random_U (N);
    auto U2 = random_U (N);
    auto U12 = U1 * U2;


    Vector occs (N);
    for(int i = 0; i < N; i++)
        occs(i) = 0.5;

    // U1|psi>
    auto psi1 = psi_rotate (sites, psi, U1, occs, 0);
    // U2.U1|psi>
    auto psi2 = psi_rotate (sites, psi1, U2, occs, 0);
    // U12|psi>
    auto psi12 = psi_rotate (sites, psi, U12, occs, 0);


    CMatrix Ut;
    Vector occs2, occs12;

    auto corr2 = MeasureCorr (sites, psi2, "Cdag", "C");
    diagHermitian (corr2, Ut, occs2);

    auto corr12 = MeasureCorr (sites, psi12, "Cdag", "C");
    diagHermitian (corr12, Ut, occs12);

    for(int i = 0; i < N; i++)
        cout << "occ " << occs2(i) << " " << occs12(i) << " " << (occs2(i)-occs12(i)) << endl;
}

void test2 (const Fermion& sites, MPO Ht, int time_steps, Sweeps sweepst, Real dt)
{
    int N = length(sites);
    auto U1 = random_U (N);

    Vector occs_tmp (N);
    for(int i = 0; i < N; i++)
        occs_tmp(i) = 0.5;

    auto psi = make_initstate (sites, N/2);
    auto psi1 = psi_rotate (sites, psi, U1, occs_tmp, 0);

    // Evolve psi
    Args args_tdvp  = {"Quiet",true,"NumCenter",2,"DoNormalize",true,"Truncate",true};
    auto obs = TDVPObserver (sites, psi);
    for(int step = 1; step <= time_steps; step++)
    {
        cout << "step = " << step << endl;
        tdvp (psi, Ht, 1_i*dt, sweepst, obs, args_tdvp);
    }
    // Evolve psi1
    for(int step = 1; step <= time_steps; step++)
    {
        cout << "step = " << step << endl;
        tdvp (psi1, Ht, 1_i*dt, sweepst, obs, args_tdvp);
    }

    CMatrix Ut;
    Vector occs, occs1;

    auto corr = MeasureCorr (sites, psi, "Cdag", "C");
    diagHermitian (corr, Ut, occs);

    auto corr1 = MeasureCorr (sites, psi1, "Cdag", "C");
    diagHermitian (corr1, Ut, occs1);

    for(int i = 0; i < N; i++)
        cout << "occ " << occs(i) << " " << occs1(i) << endl;

    int gg; cin>>gg;
}

CMatrix identity (int N)
{
    auto I = CMatrix (N,N);
    for(int i = 0; i < N; i++)
        I(i,i) = 1.;
    return I;
}

void test3 (const Fermion& sites, int time_steps, Sweeps sweepst, Real dt, MPS psi,
            int L_lead, int L_device,
            Real t_lead, Real t_device, Real t_contactL, Real t_contactR,
            Real mu_leadL, Real mu_leadR, Real mu_device,
            Real V_lead, Real V_device, Real V_contact)
{
    int N = length(sites);
    auto U = random_U (N);

    Vector occs_tmp (N);
    for(int i = 0; i < N; i++)
        occs_tmp(i) = 0.5;

    // H
    auto [ampot, idev_first, idev_last] = t_mu_V_ampo (sites, L_lead, L_device,
                                                      t_lead, t_lead, t_device, t_contactL, t_contactR,
                                                      mu_leadL, mu_leadR, mu_device,
                                                      V_lead, V_lead, V_device, V_contact, V_contact);
    auto H = toMPO (ampot);

    // rotate psi
    auto psi2 = psi_rotate (sites, psi, U, occs_tmp, 0);
    // rotate H
    auto H2 = H_rot (sites, U, L_lead, L_device, t_lead, t_device, t_contactL, t_contactR,
                     mu_leadL, mu_leadR, mu_device, V_lead, V_device, V_contact);


    // Evolve psi
    Args args_tdvp  = {"Quiet",true,"NumCenter",2,"DoNormalize",true,"Truncate",true};
    auto obs = TDVPObserver (sites, psi);
    for(int step = 1; step <= time_steps; step++)
    {
        cout << "step = " << step << endl;
        tdvp (psi, H, 1_i*dt, sweepst, obs, args_tdvp);
    }
    // Evolve psi2
    for(int step = 1; step <= time_steps; step++)
    {
        cout << "step = " << step << endl;
        tdvp (psi2, H2, 1_i*dt, sweepst, obs, args_tdvp);
    }


    auto en = real(innerC(psi,H,psi));
    auto en2 = real(innerC(psi2,H2,psi2));
    cout << "en " << en << " " << en2 << " " << (en-en2) << endl;

    CMatrix Ut;
    Vector occs, occs2;

    auto corr = MeasureCorr (sites, psi, "Cdag", "C");
    diagHermitian (corr, Ut, occs);

    auto corr2 = MeasureCorr (sites, psi2, "Cdag", "C");
    diagHermitian (corr2, Ut, occs2);

    for(int i = 0; i < N; i++)
        cout << "occ " << occs(i) << " " << occs2(i) << " " << (occs(i)-occs2(i)) << endl;

    int gg; cin>>gg;
}

tuple<MPS,MPO,CMatrix> rotate_psi_H (const Fermion& sites, const MPS& psi, const MPO& H, const CMatrix& U, Real zero_crit,
                                     int L_lead, int L_device,
                                     Real t_lead, Real t_device, Real t_contactL, Real t_contactR,
                                     Real mu_leadL, Real mu_leadR, Real mu_device,
                                     Real V_lead, Real V_device, Real V_contact)
{
    // Compute the correlation matrix and natural orbitals
    CMatrix Ut;
    Vector occs;
    CMatrix corr = MeasureCorr (sites, psi, "Cdag", "C");
    diagHermitian (corr, Ut, occs);

    // Update U
    auto Ur = U * Ut;

    // Update psi
    auto psir = psi_rotate (sites, psi, Ut, occs, zero_crit);

    // Update H
    auto Hr = H_rot (sites, Ur, L_lead, L_device, t_lead, t_device, t_contactL, t_contactR,
               mu_leadL, mu_leadR, mu_device, V_lead, V_device, V_contact);

    cout << "E " << innerC(psi,H,psi) << " " << innerC(psir,Hr,psir) << endl;
    return {psir, Hr, Ur};
}

int main(int argc, char* argv[])
{

    string infile = argv[1];
    InputGroup input (infile,"basic");

    auto L_lead   = input.getInt("L_lead");
    auto L_device   = input.getInt("L_device");
    auto t_lead     = input.getReal("t_lead");
    auto t_device   = input.getReal("t_device");
    auto t_contactL  = input.getReal("t_contactL");
    auto t_contactR  = input.getReal("t_contactR");
    auto mu_leadL   = input.getReal("mu_leadL");
    auto mu_leadR   = input.getReal("mu_leadR");
    auto mu_device  = input.getReal("mu_device");
    auto V_lead     = input.getReal("V_lead");
    auto V_device   = input.getReal("V_device");
    auto V_contact  = input.getReal("V_contact");

    auto do_write   = input.getYesNo("write_to_file");
    auto out_dir    = input.getString("outdir",".");
    auto out_minm   = input.getInt("out_minm",0);
    auto ConserveQNs = input.getYesNo("ConserveQNs",false);
    auto ConserveNf  = input.getYesNo("ConserveNf",false);
    auto WriteDim   = input.getInt("WriteDim",-1);
    auto sweeps     = iut::Read_sweeps (infile, "DMRG");

    auto mu_biasL      = input.getReal("mu_biasL");
    auto mu_biasS      = input.getReal("mu_biasS");
    auto mu_biasR      = input.getReal("mu_biasR");
    auto dt            = input.getReal("dt");
    auto time_steps    = input.getInt("time_steps");
    auto NumCenter     = input.getInt("NumCenter");
    auto Truncate      = input.getYesNo("Truncate");
    auto sweepst       = iut::Read_sweeps (infile, "TDVP");

    cout << "device site = " << L_lead+1 << " " << (L_lead + L_device) << endl;

    // Site set
    using SitesType = Fermion;
    int L = 2*L_lead + L_device;
    auto sites = SitesType (L, {"ConserveQNs",ConserveQNs,"ConserveNf",ConserveNf});

    // Initialze MPS
    MPS psi;
    int Np = L / 2;
    psi = make_initstate (sites, Np);
    psi.position(1);

    // Make MPO    
    auto [ampo, idev_first, idev_last] = t_mu_V_ampo (sites, L_lead, L_device,
                                                      t_lead, t_lead, t_device, t_contactL, t_contactR,
                                                      mu_leadL, mu_leadR, mu_device,
                                                      V_lead, V_lead, V_device, V_contact, V_contact);
    auto H = toMPO (ampo);
    cout << "MPO dim = " << maxLinkDim(H) << endl;

int pause;
auto zero_crit = input.getReal("zero_crit");   // critiria for close to 0 or 1
auto rotate_inteval = input.getInt("rotate_inteval");   // critiria for close to 0 or 1



    // DMRG
    MyObserver<SitesType> myobs (sites, psi, {"Write",do_write,"out_dir",out_dir,"out_minm",out_minm});
    dmrg (psi, H, sweeps, myobs, {"WriteDim",WriteDim});


    // Time evolution H
    AutoMPO ampot (sites);
    tie (ampot, idev_first, idev_last) = t_mu_V_ampo (sites, L_lead, L_device,
                                                      t_lead, t_lead, t_device, t_contactL, t_contactR,
                                                      mu_leadL+mu_biasL, mu_leadR+mu_biasR, mu_device+mu_biasS,
                                                      V_lead, V_lead, V_device, V_contact, V_contact);
    auto Ht = toMPO (ampot);


    // -- rotate --
    auto U = identity (L);
/*    if (rotate_inteval != 0)
    {
        tie(psi,Ht,U) = rotate_psi_H (sites, psi, Ht, U, zero_crit,
                                      L_lead, L_device, t_lead, t_device, t_contactL, t_contactR,
                                      mu_leadL+mu_biasL, mu_leadR+mu_biasR, mu_device+mu_biasS, V_lead, V_device, V_contact);
    }*/
//test2 (sites, Ht, rotate_inteval, sweepst, dt);

    /*auto tij = get_tij (L_lead, L_device,
                        t_lead, t_lead, t_device, t_contactL, t_contactR,
                        mu_leadL+mu_biasL, mu_leadR+mu_biasR, mu_device+mu_biasS);
    add_t_rot (ampot, tij, U);
    auto Ht = toMPO(ampot);
    cout << "Ht MPO bond dim after rotation = " << maxLinkDim(Ht) << endl;*/


    // Current
    /*int N = length (psi);
    vector<int> spec_links = {L_lead, L_lead+L_device};
    vector<MPO> JMPOs (N);
    for(int i = 1; i < L; i++)
    {
        AutoMPO ampoj (sites);
        ampoj += 1.,"Cdag",i,"C",i+1;
        JMPOs.at(i) = toMPO(ampoj);
    }*/

    // TDVP
    Args args_tdvp  = {"Quiet",true,"NumCenter",NumCenter,"DoNormalize",true,"Truncate",Truncate};
    auto obs = TDVPObserver (sites, psi);
    for(int step = 1; step <= time_steps; step++)
    {
        cout << "step = " << step << endl;
        tdvp (psi, Ht, 1_i*dt, sweepst, obs, args_tdvp);


        if (step % rotate_inteval == 0)
        {
            tie(psi,Ht,U) = 
                            rotate_psi_H (sites, psi, Ht, U, zero_crit,
                                          L_lead, L_device, t_lead, t_device, t_contactL, t_contactR,
                                          mu_leadL+mu_biasL, mu_leadR+mu_biasR, mu_device+mu_biasS, V_lead, V_device, V_contact);

            //test (sites, psi);
            //test3 (sites, rotate_inteval, sweepst, dt, psi,
            //       L_lead, L_device, t_lead, t_device, t_contactL, t_contactR,
            //       mu_leadL+mu_biasL, mu_leadR+mu_biasR, mu_device+mu_biasS, V_lead, V_device, V_contact);
/*{
    CMatrix Ut;
    Vector occs;
    CMatrix corr = MeasureCorr (sites, psi, "Cdag", "C");
    diagHermitian (corr, Ut, occs);
    for(int i = 0; i < L; i++)
        cout << "occ " << occs(i) << endl;
}*/
            //cout << "pause" << endl; int pp; cin>>pp;
        }

        // Measure currents
        /*for(int j = 1; j < L; j++)
        {
            auto Jtmp = innerC (psi, JMPOs.at(j), psi);
            auto J = -2. * imag(Jtmp);
            cout << "\tcurrent " << j << " " << j+1 << " = " << J << endl;
            //cout << "\t*current spec " << j << " " << j+1 << " = " << J << endl;
        }*/
        
        {
            // measure the correlations
            CMatrix corr = MeasureCorr (sites, psi, "Cdag", "C");
            // Compute currents from correlations
            auto corr_realspace = U * corr * conj(transpose(U));
            for(int j = 1; j < L; j++)
            {
                auto J = -2. * imag(corr_realspace(j-1,j));
                cout << "\tcurrent " << j << " " << j+1 << " = " << J << endl;
            }
        }
    }
    return 0;
}
