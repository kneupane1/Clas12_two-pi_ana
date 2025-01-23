
#ifndef HIST_H_GUARD
#define HIST_H_GUARD
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TThread.h"
#include "colors.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "deltat.hpp"
#include "reaction.hpp"
#include <mutex>

using namespace std;

using TH2D_ptr = std::shared_ptr<TH2D>;
using TH1D_ptr = std::shared_ptr<TH1D>;
using THnSparse_ptr = std::shared_ptr<THnSparse>;
using TGraph_ptr = std::shared_ptr<TGraph>;

class Histogram
{
protected:
    std::shared_ptr<TFile> RootOutputFile;
    std::shared_ptr<TCanvas> def;

    int bins = 500;
    double p_min = 0.0;
    double p_max = 6.0;
    double Dt_max = 10.0;
    double Dt_min = -Dt_max;
    double q2_min = 0.0;
    double q2_max = 15.0;

    double w_max = 3.5;
    double w_min = 1.0;

    double zero = 0.0;

    static const short particle_num = 2; // 1-Pi 2-P
    std::string particle_name[particle_num] = {"pi", "Prot"};
    static const short charge_num = 2; // 0-pos 1-neg
    std::string charge_name[charge_num] = {"positive", "negative"};
    static const short with_id_num = 2; // 0-without 1-with
    std::string id_name[with_id_num] = {"withoutID", "withID"};

    static const short num_sectors = 6;
    std::string sec_name[num_sectors] = {"1", "2", "3", "4", "5", "6"};
    static const short CUTS = 4;
    enum cuts
    {
        before_any_cuts,
        with_one_cut,
        outside_one_cut,
        after_all_cuts
    };
    std::mutex mutex;

    TH1D_ptr inv_mass_pPip;
    TH1D_ptr inv_mass_pPim;
    TH1D_ptr inv_mass_pipPim;

    TH1D_ptr W_hist;
    TH1D_ptr Q2_hist;
    TH2D_ptr W_vs_q2;
    TH1D_ptr W_P2pi_hist;

    TH1D_ptr W_thrown;
    TH2D_ptr W_vs_Q2_thrown;
    TH1D_ptr Q2_thrown;

    TH1D_ptr vz_position[CUTS];
    TH2D_ptr pcal_sec[CUTS];
    TH2D_ptr pcal_hx_hy_sec[CUTS];
    TH2D_ptr dcr1_sec[CUTS];
    TH2D_ptr dcr2_sec[CUTS];
    TH2D_ptr dcr3_sec[CUTS];

    // EC Sampling Fraction
    TH2D_ptr EC_sampling_fraction[CUTS];
    TH2D_ptr ECin_sf_vs_PCAL_sf[CUTS];
    TH1D_ptr momentum[CUTS];

    //// Hadron pid
    TH1D_ptr prot_Delta_vz_cut_fd[CUTS];
    TH1D_ptr prot_Chi2pid_cut_fd[CUTS];
    TH1D_ptr pip_Delta_vz_cut_fd[CUTS];
    TH1D_ptr pip_Chi2pid_cut_fd[CUTS];
    TH1D_ptr pim_Delta_vz_cut[CUTS];
    TH1D_ptr pim_Chi2pid_cut[CUTS];

    TH1D_ptr prot_Delta_vz_cut_cd[CUTS];
    TH1D_ptr prot_Chi2pid_cut_cd[CUTS];
    TH1D_ptr pip_Delta_vz_cut_cd[CUTS];
    TH1D_ptr pip_Chi2pid_cut_cd[CUTS];

    TH2D_ptr dcr1_sec_prot[CUTS];
    TH2D_ptr dcr2_sec_prot[CUTS];
    TH2D_ptr dcr3_sec_prot[CUTS];
    TH2D_ptr dcr1_sec_pip[CUTS];
    TH2D_ptr dcr2_sec_pip[CUTS];
    TH2D_ptr dcr3_sec_pip[CUTS];
    TH2D_ptr dcr1_sec_pim[CUTS];
    TH2D_ptr dcr2_sec_pim[CUTS];
    TH2D_ptr dcr3_sec_pim[CUTS];

    TH2D_ptr W_vs_q2_sec[num_sectors];
    TH1D_ptr W_sec[num_sectors];

    //////////////////////////   electron pid cuts ///////////////
    TH1D_ptr htcc_nphe_sec[CUTS][num_sectors];
    TH1D_ptr elec_Chi2pid_sec[CUTS][num_sectors];
    TH1D_ptr vz_sec[CUTS][num_sectors];
    TH2D_ptr SF_VS_MOM[CUTS][num_sectors];

    TH1D_ptr MM2_twoPi_excl;
    TH1D_ptr MM_twoPi_excl;
    TH1D_ptr missing_Energy_hist;

    TH1D_ptr MM2_twoPi_mPim;
    TH1D_ptr MM_twoPi_mPim;
    TH1D_ptr MM2_twoPi_missingPip;
    TH1D_ptr MM2_twoPi_missingProt;

    // Mom vs Beta
    TH2D_ptr momvsbeta_hist[particle_num][charge_num][with_id_num];
    TH2D_ptr momvsbeta_hist_prot[2][3];
    TH2D_ptr momvsbeta_hist_pip[2][3];
    TH2D_ptr momvsbeta_hist_pim[2][3];

    // Delta T
    TH2D_ptr delta_t_hist[3][2][3]; //

public:
    Histogram(const std::string &output_file);

    ~Histogram();
    // Constructor

    // W and Q^2
    void makeHists_sector();
    void Fill_WvsQ2(const std::shared_ptr<Reaction> &_e);
    void Fill_WvsQ2_twoPi_thrown(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<MCReaction> &_e);
    void Write_WvsQ2();

    void Fill_W_vs_Q2_thrown();
    void Fill_inv_mass_hist();

    // P and E
    // ecectron cuts
    void makeHists_electron_cuts();
    void FillHists_electron_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e);
    void FillHists_electron_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e);

    void FillHists_prot_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i);
    void FillHists_prot_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i, const TLorentzVector &prot);

    void FillHists_pip_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i);
    void FillHists_pip_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i, const TLorentzVector &pip);
    void FillHists_pim_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i);
    void FillHists_pim_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i);

    void Write_Electron_cuts();
    void Write_Proton_cuts();
    void Write_Pip_cuts();
    void Write_Pim_cuts();

    void makeHists_MomVsBeta();
    void Fill_MomVsBeta(const std::shared_ptr<Branches12> &data, int part, const std::shared_ptr<Reaction> &_e);
    void Write_MomVsBeta();

    // Delta T
    void makeHists_deltat();

    void Fill_deltat_before_cut(const std::shared_ptr<Branches12> &data,
                                const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_prot_after_cut(const std::shared_ptr<Branches12> &data,
                                    const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_pip_after_cut(const std::shared_ptr<Branches12> &data,
                                   const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Fill_deltat_pim_after_cut(const std::shared_ptr<Branches12> &data,
                                   const std::shared_ptr<Delta_T> &dt, int part, const std::shared_ptr<Reaction> &_e);
    void Write_deltat();

    void Write();
};

#endif
