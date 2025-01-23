#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <iostream>
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "branches.hpp"
#include "constants.hpp"
#include "physics.hpp"
#include <TRandom3.h>
#include <cmath>
#include <vector>
class Reaction
{
protected:
        std::shared_ptr<Branches12> _data;

        double _beam_energy = 10.6041;

        std::unique_ptr<TLorentzVector> _beam;
        std::unique_ptr<TLorentzVector> _elec;
        std::unique_ptr<TLorentzVector> _gamma;
        std::unique_ptr<TLorentzVector> _target;

        std::unique_ptr<TLorentzVector> _Energy_loss_uncorr_pim;
        std::unique_ptr<TLorentzVector> _Energy_loss_uncorr_pip;
        std::unique_ptr<TLorentzVector> _Energy_loss_uncorr_prot;

        std::vector<std::shared_ptr<TLorentzVector>> _photons;

        std::unique_ptr<TLorentzVector> _other;
        std::unique_ptr<TLorentzVector> _neutron;

        std::unique_ptr<TLorentzVector> _elecUnSmear;
        std::unique_ptr<TLorentzVector> _protUnSmear;
        std::unique_ptr<TLorentzVector> _pipUnSmear;
        std::unique_ptr<TLorentzVector> _pimUnSmear;

        std::vector<std::unique_ptr<TLorentzVector>> _prot;
        std::vector<std::unique_ptr<TLorentzVector>> _pip;
        std::vector<std::unique_ptr<TLorentzVector>> _pim;

        std::vector<int> _prot_indices;
        std::vector<int> _pip_indices;
        std::vector<int> _pim_indices;

        std::unique_ptr<TLorentzVector> _boosted_pim_measured;

        float _weight = NAN;

        bool _is_eff_corrected = false;

        bool _is_boosted = false;

        bool _is_boosted_swapped = false;

        bool _hasE = false;
        bool _hasP = false;
        bool _hasPip = false;
        bool _hasPim = false;
        bool _hasOther = false;
        bool _hasNeutron = false;

        short _numPart = 0;
        short _numProt = 0;
        short _numPip = 0;
        short _numPim = 0;
        short _numPos = 0;
        short _numNeg = 0;
        short _numNeutral = 0;
        short _numOther = 0;
        short _numPhoton = 0;
        short _sector = -1;

        bool _is_FD_Prot = false;
        bool _is_CD_Prot = false;

        bool _is_FD_Pip = false;
        bool _is_CD_Pip = false;

        bool _is_FD_Pim = false;
        bool _is_CD_Pim = false;
        bool _is_FD = false;
        bool _is_CD = false;

        float _MM_mPim = NAN;
        float _MM2_mPim = NAN;
        float _MM_mpip = NAN;
        float _MM2_mpip = NAN;
        float _MM_mprot = NAN;
        float _MM2_mprot = NAN;
        float _MM2_exclusive = NAN;
        float _MM_exclusive = NAN;
        float _W = NAN;
        float _Q2 = NAN;
        float _MM2_mPim_swapped = NAN;

        float _inv_Ppip = NAN;
        float _inv_Ppim = NAN;
        float _inv_pip_pim = NAN;

        float _phi_gamma = NAN;
        float _phi_prot = NAN;
        float _phi_pip = NAN;
        float _phi_pim = NAN;

        float _alpha_ppip_pipim = NAN;
        float _alpha_pippim_pipf = NAN;
        float _alpha_ppim_pipip = NAN;

        float _alpha_ppip_pipim_swapped = NAN;
        float _alpha_pippim_pipf_swapped = NAN;
        float _alpha_ppim_pipip_swapped = NAN;

        float _beam_theta = NAN;
        float _elec_theta = NAN;
        float _E_elec = NAN;
        float _pim_theta_measured = NAN;
        short _pim_sec = -9999;
        float _theta_e = NAN;
        float _P_elec = NAN;
        float _elec_status;
        float _prot_status = NAN;
        float _pip_status = NAN;
        float _pim_status = NAN;

        int _sectorElec = -1;
        int _sectorPim = -1;
        int _sectorPip = -1;
        int _sectorProt = -1;
        // elec mom corr
        double _cx = NAN;
        double _cy = NAN;
        double _cz = NAN;

        double _px_prime_elec = NAN;
        double _py_prime_elec = NAN;
        double _pz_prime_elec = NAN;

        // energy loss corr

        double _px_prime_prot_E = NAN;
        double _py_prime_prot_E = NAN;
        double _pz_prime_prot_E = NAN;

        double _prot_mom = NAN;
        double _prot_mom_uncorr = NAN;
        float _E_corr_val_prot = NAN;
        double _prot_theta_uncorr = NAN;
        float _prot_phi_uncorr = NAN;
        double _prot_mom_tmt = NAN;
        double _prot_mom_prime = NAN;
        double _prot_theta = NAN;
        double _prot_phi = NAN;

        double _px_prime_pip_E = NAN;
        double _py_prime_pip_E = NAN;
        double _pz_prime_pip_E = NAN;

        double _pip_mom = NAN;
        double _pip_mom_tmt = NAN;
        double _pip_mom_uncorr = NAN;
        float _E_corr_val_pip = NAN;
        double _pip_theta_uncorr = NAN;
        double _pip_phi_uncorr = NAN;
        double _pip_mom_prime = NAN;

        double _pip_theta = NAN;
        double _pip_phi = NAN;

        double _px_prime_pim_E = NAN;
        double _py_prime_pim_E = NAN;
        double _pz_prime_pim_E = NAN;

        double _pim_mom = NAN;
        double _pim_mom_prime = NAN;
        double _pim_mom_tmt = NAN;
        double _pim_mom_uncorr = NAN;
        float _E_corr_val_pim = NAN;
        double _pim_theta_uncorr = NAN;
        double _pim_phi_uncorr = NAN;
        double _pim_theta = NAN;
        double _pim_phi = NAN;

        double _px_prime_prot_mom = NAN;
        double _py_prime_prot_mom = NAN;
        double _pz_prime_prot_mom = NAN;

        double _px_prime_pip_mom = NAN;
        double _py_prime_pip_mom = NAN;
        double _pz_prime_pip_mom = NAN;

        double _px_prime_pim_mom = NAN;
        double _py_prime_pim_mom = NAN;
        double _pz_prime_pim_mom = NAN;

        float _excl_Energy = NAN;

        double fe = NAN;
        double fpro = NAN;
        double fpip = NAN;
        double fpim = NAN;

        float _pr_p = NAN;
        float _pr_th = NAN;
        float _pr_ph_eff = NAN;

        float _pip_p = NAN;
        float _pip_th = NAN;
        float _pip_ph_eff = NAN;

        float _pim_p = NAN;
        float _pim_th = NAN;
        float _pim_ph_eff = NAN;

        float _eff_corr_fact_mPim = NAN;
        float _eff_corr_fact_Excl = NAN;
        //
        static const int CD_SEC = 3;
        static const int FD_SEC = 6;

        void SetElec();

        double _elec_mom_corrected = NAN;
        double _elec_mom = NAN;

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

public:
        Reaction() {};
        Reaction(const std::shared_ptr<Branches12> &data, float beam_energy);
        ~Reaction();

        inline bool mc()
        {
                return _mc;
        }

        void SetProton(int i);
        void SetPip(int i);
        void SetPim(int i);
        void SetOther(int i);
        void SetNeutron(int i);

        const std::vector<std::unique_ptr<TLorentzVector>> &GetProtons() const { return _prot; }
        const std::vector<std::unique_ptr<TLorentzVector>> &GetPips() const { return _pip; }
        const std::vector<std::unique_ptr<TLorentzVector>> &GetPims() const { return _pim; }

        const std::vector<int> &GetProtonIndices() const { return _prot_indices; }
        const std::vector<int> &GetPipIndices() const { return _pip_indices; }
        const std::vector<int> &GetPimIndices() const { return _pim_indices; }

        void CalcMissMassPim(const TLorentzVector &prot, const TLorentzVector &pip);
        void CalcMissMassPip(const TLorentzVector &prot, const TLorentzVector &pim);
        void CalcMissMassProt(const TLorentzVector &pip, const TLorentzVector &pim);
        void CalcMissMassExcl(const TLorentzVector &prot, const TLorentzVector &pip, const TLorentzVector &pim);
        float MM_mPim();
        float MM2_mPim();
        float MM2_mpip();
        float MM2_mprot();
        float MM2_exclusive();
        float MM_exclusive();
        inline float weight()
        {
                if (_mc)
                        return _data->mc_weight();

                else
                        return 1.0;
        }

        /// smearing fx's function
        void SmearingFunc(int part_id, int status_part, double p, double theta, double phi, double w_val, double &pNew, double &thetaNew,
                          double &phiNew)
        {
                //
                const double pS1 = 0.0184291 - 0.0110083 * theta + 0.00227667 * pow(theta, 2) - 0.000140152 * pow(theta, 3) +
                                   3.07424e-6 * pow(theta, 4);
                const double pR = 0.02 * sqrt(pow(pS1 * p, 2) + pow(0.02 * theta, 2));
                const double thetaR = 2.5 * sqrt(pow((0.004 * theta + 0.1) * (pow(p, 2) + 0.13957 * 0.13957) / pow(p, 2), 2));
                const double phiS1 = 0.85 - 0.015 * theta;
                const double phiS2 = 0.17 - 0.003 * theta;
                const double phiR = 3.5 * sqrt(pow(phiS1 * sqrt(pow(p, 2) + 0.13957 * 0.13957) / pow(p, 2), 2) + pow(phiS2, 2));

                // Generate new values
                if (part_id == ELECTRON)
                {
                        phiNew = phi + phiR * gRandom->Gaus(0, 1) * 0.4;
                        thetaNew = theta + thetaR * gRandom->Gaus(0, 1) * 0.4;
                        pNew = p + pR * gRandom->Gaus(0, 1) * p * 0.4;
                }
                else if (part_id == PROTON)
                {
                        double fact_cd = 0;
                        double fact_fd = 0;
                        double fact_cd1 = 0;
                        double fact_fd1 = 0;
                        if (status_part > 4000)
                        {
                                fact_cd = (0.000821) * pow(p, 3) + (-0.016500) * pow(p, 2) + (0.103611) * p + (1.393237);
                                fact_cd1 = 1.0;

                                phiNew = phi + 1 / (fact_cd)*phiR * gRandom->Gaus(0, 1);
                                thetaNew = theta + 1 / (fact_cd)*thetaR * gRandom->Gaus(0, 1);
                                pNew = p + 1 / (fact_cd)*pR * gRandom->Gaus(0, 1) * p;
                        }
                        else if (status_part <= 4000)
                        {
                                fact_fd = (0.000264) * pow(p, 3) + (-0.006454) * pow(p, 2) + (0.032683) * p + (1.658142);
                                fact_fd1 = 1.0;

                                phiNew = phi + 1 / (fact_fd * fact_fd1) * phiR * gRandom->Gaus(0, 1);
                                thetaNew = theta + 1 / (fact_fd * fact_fd1) * thetaR * gRandom->Gaus(0, 1);
                                pNew = p + 1 / (fact_fd * fact_fd1) * pR * gRandom->Gaus(0, 1) * p;
                        }
                }

                else if (part_id == PIP)
                {
                        double fact_cd = 0;
                        double fact_fd = 0;
                        double fact_cd1 = 0;
                        double fact_fd1 = 0;
                        if (status_part > 4000)
                        {
                                fact_cd = (0.000981) * pow(p, 3) + (-0.016882) * pow(p, 2) + (0.046752) * p + (1.720426);
                                fact_cd1 = 1.0;

                                phiNew = phi + 1 / (fact_cd * fact_cd1) * phiR * gRandom->Gaus(0, 1);
                                thetaNew = theta + 1 / (fact_cd * fact_cd1) * thetaR * gRandom->Gaus(0, 1);
                                pNew = p + 1 / (fact_cd * fact_cd1) * pR * gRandom->Gaus(0, 1) * p;
                        }
                        else if (status_part <= 4000)
                        {
                                fact_fd = (0.000085) * pow(p, 3) + (-0.003096) * pow(p, 2) + (0.023553) * p + (1.509910);
                                fact_fd1 = 1.0;
                                phiNew = phi + 1 / (fact_fd * fact_fd1) * phiR * gRandom->Gaus(0, 1);
                                thetaNew = theta + 1 / (fact_fd * fact_fd1) * thetaR * gRandom->Gaus(0, 1);
                                pNew = p + 1 / (fact_fd * fact_fd1) * pR * gRandom->Gaus(0, 1) * p;
                        }
                }

                else if (part_id == PIM)
                {
                        double fact_cd = 0;
                        double fact_fd = 0;
                        double fact_cd1 = 0;
                        double fact_fd1 = 0;
                        if (status_part > 4000)
                        {
                                fact_cd = (-0.001788) * pow(p, 3) + (0.025796) * pow(p, 2) + (-0.136577) * p + (2.007917);
                                fact_cd1 = 1.0;

                                phiNew = phi + 1 / (fact_cd * fact_cd1) * phiR * gRandom->Gaus(0, 1);
                                thetaNew = theta + 1 / (fact_cd * fact_cd1) * thetaR * gRandom->Gaus(0, 1);
                                pNew = p + 1 / (fact_cd * fact_cd1) * pR * gRandom->Gaus(0, 1) * p;
                        }
                        else if (status_part <= 4000)
                        {
                                fact_fd = (0.000760) * pow(p, 3) + (-0.021295) * pow(p, 2) + (0.171180) * p + (1.238299);
                                fact_fd1 = 1.0;
                                phiNew = phi + 1 / (fact_fd * fact_fd1) * phiR * gRandom->Gaus(0, 1);
                                thetaNew = theta + 1 / (fact_fd * fact_fd1) * thetaR * gRandom->Gaus(0, 1);
                                pNew = p + 1 / (fact_fd * fact_fd1) * pR * gRandom->Gaus(0, 1) * p;
                        }
                }
        }

        ///////////////////////////////////////////////////////
        void invMassPpip(const TLorentzVector &prot, const TLorentzVector &pip);
        void invMassPpim(const TLorentzVector &prot, const TLorentzVector &pip);
        void invMasspippim(const TLorentzVector &prot, const TLorentzVector &pip);

        float inv_Ppip();
        float inv_Ppim();
        float inv_pip_pim();

        ///////////// related to lab system ////////////
        float prot_momentum(const TLorentzVector &prot);
        float prot_theta_lab(const TLorentzVector &prot);
        float prot_Phi_lab(const TLorentzVector &prot);
        float prot_momT(const TLorentzVector &prot);

        float pip_momentum(const TLorentzVector &pip);
        float pip_theta_lab(const TLorentzVector &pip);
        float pip_Phi_lab(const TLorentzVector &pip);
        float pip_momT(const TLorentzVector &pip);

        float pim_momentum(const TLorentzVector &prot, const TLorentzVector &pip);
        float pim_momentum_measured(const TLorentzVector &pim);
        float pim_theta_lab(const TLorentzVector &prot, const TLorentzVector &pip);
        float pim_theta_lab_measured(const TLorentzVector &pim);
        float pim_E(const TLorentzVector &prot, const TLorentzVector &pip);
        float pim_E_measured(const TLorentzVector &prot);
        float pim_Phi_lab(const TLorentzVector &prot, const TLorentzVector &pip);
        float pim_Phi_lab_measured(const TLorentzVector &pim);
        float pim_momT(const TLorentzVector &pim);

        float elec_momentum();
        float theta_elec();
        float Phi_elec();

        float Energy_excl();

        inline float W()
        {
                return _W;
        }
        inline float Q2()
        {
                return _Q2;
        }

        inline short sec()
        {
                return _data->dc_sec(0);
        }
        inline int cc_tot_nphe()
        {
                return _data->cc_nphe_tot(0);
        }
        inline int vz()
        {
                return _data->vz(0);
        }
        inline int det()
        {
                return abs(_data->status(0) / 1000);
        }

        inline bool Fixed_MM_cut()
        {
                // exp
                return (Reaction::MM2_mPim() < 0.08 && Reaction::MM2_mPim() > -0.06);
        }

        inline bool TwoPion_missingPim()
        {
                bool _channelTwoPi = true;
                _channelTwoPi &= ((_numProt >= 1 && _numPip >= 1) && (_hasE && _hasP && _hasPip));
                return _channelTwoPi;
        }
        inline bool TwoPion_exclusive()
        {
                bool _channelTwoPi_excl = true;

                _channelTwoPi_excl &= ((_numProt >= 1 && _numPip >= 1 && _numPim >= 1) && (_hasE && _hasP && _hasPip && _hasPim));
                return _channelTwoPi_excl;
        }
        inline bool TwoPion_missingPip()
        {
                bool _channelTwoPi_mpip = true;

                _channelTwoPi_mpip &= ((_numProt >= 1 && _numPim >= 1) &&
                                       (_hasE && _hasP && _hasPim));
                return _channelTwoPi_mpip;
        }
        inline bool TwoPion_missingProt()
        {
                bool _channelTwoPi_mprot = true;
                _channelTwoPi_mprot &= ((_numPip >= 1 && _numPim >= 1) && (_hasE && _hasPip && _hasPim));
                return _channelTwoPi_mprot;
        }
        inline bool inclusive()
        {
                return (_hasE);
        }

        const TLorentzVector &e_mu()
        {
                return *_beam;
        }
        const TLorentzVector &e_mu_prime()
        {
                return *_elec;
        }
        const TLorentzVector &gamma()
        {
                return *_gamma;
        }
};

class MCReaction : public Reaction
{
private:
        float _weight_mc = NAN;
        float _W_mc = NAN;
        float _Q2_mc = NAN;

        std::unique_ptr<TLorentzVector> _beam_mc;
        std::unique_ptr<TLorentzVector> _elec_mc;
        std::unique_ptr<TLorentzVector> _gamma_mc;
        std::unique_ptr<TLorentzVector> _other_mc;

        std::vector<std::unique_ptr<TLorentzVector>> _prot_mc;
        std::vector<std::unique_ptr<TLorentzVector>> _pip_mc;
        std::vector<std::unique_ptr<TLorentzVector>> _pim_mc;

        std::vector<int> _prot_mc_indices;
        std::vector<int> _pip_mc_indices;
        std::vector<int> _pim_mc_indices;

        bool _is_boosted_mc = false;

        float _MM_mc = NAN;
        float _MM2_mc = NAN;
        float _MM2_exclusive_mc = NAN;

        float _alpha_ppip_pipim_mc = NAN;
        float _alpha_pippim_pipf_mc = NAN;
        float _alpha_ppim_pipip_mc = NAN;

public:
        MCReaction(const std::shared_ptr<Branches12> &data, float beam_energy);
        void SetMCElec();
        inline float weight()
        {
                return _data->mc_weight();
                // return 1.0;
        }
        inline float W_mc()
        {
                return _W_mc;
        }
        inline float Q2_mc()
        {
                return _Q2_mc;
        }

        void SetMCProton(int i);
        void SetMCPip(int i);
        void SetMCPim(int i);
        void SetMCOther(int i);

        const std::vector<std::unique_ptr<TLorentzVector>> &GetMcProtons() const { return _prot_mc; }
        const std::vector<std::unique_ptr<TLorentzVector>> &GetMcPips() const { return _pip_mc; }
        const std::vector<std::unique_ptr<TLorentzVector>> &GetMcPims() const { return _pim_mc; }

        const std::vector<int> &GetProtonMcIndices() const { return _prot_mc_indices; }
        const std::vector<int> &GetPipMcIndices() const { return _pip_mc_indices; }
        const std::vector<int> &GetPimMcIndices() const { return _pim_mc_indices; }

        float pim_mom_mc_gen();
        float pip_mom_mc_gen();
        float prot_mom_mc_gen();

        float MCinv_Ppip();
        float MCinv_Ppim();
        float MCinv_pip_pim();
        //

        float MCprot_theta_lab();
        float MCpip_theta_lab();
        float MCpim_theta_lab();

        float pim_phi_mc_gen();
        float pip_phi_mc_gen();
        float prot_phi_mc_gen();
};

#endif