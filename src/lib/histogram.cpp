
#include "histogram.hpp"

Histogram::Histogram(const std::string &output_file)
{
        RootOutputFile = std::make_shared<TFile>(output_file.c_str(), "RECREATE");
        def = std::make_shared<TCanvas>("def");

        inv_mass_pPip = std::make_shared<TH1D>("pPip_mass", "Prot-Pip mass", bins, 1.0, 2.25);
        inv_mass_pPim = std::make_shared<TH1D>("pPim_mass", "Prot-Pim mass", bins, 0.75, 2.25);
        inv_mass_pipPim = std::make_shared<TH1D>("pipPim_mass", "Pip-Pim mass", bins, 0.0, 1.5);

        W_hist = std::make_shared<TH1D>("W", "W Reconstructed", bins, w_min, w_max);
        Q2_hist = std::make_shared<TH1D>("Q2", "Q^{2} reconstructed", bins, zero, q2_max);
        W_vs_q2 = std::make_shared<TH2D>("W_vs_q2", "W_vs_q2 reconstructed", bins, w_min, w_max,
                                         bins, q2_min, q2_max);

        W_thrown = std::make_shared<TH1D>("W_thrown", "W thrown", bins, w_min, w_max);
        Q2_thrown = std::make_shared<TH1D>("Q2_thrown", "Q^{2} thrown", bins, q2_min, q2_max);

        W_vs_Q2_thrown =
            std::make_shared<TH2D>("W_vs_q2_thrown", "W_vs_q2_thrown", bins, w_min,
                                   w_max, bins, q2_min, q2_max);

        MM2_twoPi_excl = std::make_shared<TH1D>("MMSQ_excl", "MMSQ excl", bins,
                                                -0.1, 0.1);
        MM_twoPi_excl = std::make_shared<TH1D>("MM_excl", "MM excl", bins,
                                               -0.25, 0.25);

        missing_Energy_hist = std::make_shared<TH1D>("e#pi^{+}#pi^{-}pX", "missing Energy", 500, -0.6, 0.6);

        MM_twoPi_mPim = std::make_shared<TH1D>("MM_e#pi^{+}pX", "Missing Mass: expecting #pi^{-}", bins,
                                               -0.4, 0.4);
        MM2_twoPi_mPim = std::make_shared<TH1D>("MMSQ_e#pi^{+}pX", "MMSQ: expecting #pi^{-}", bins,
                                                -0.4, 0.4);
        MM2_twoPi_missingPip = std::make_shared<TH1D>(
            "e#pi^{-}pX", "MMSQ: expecting #pi^{+}", bins, -0.4, 0.4);

        MM2_twoPi_missingProt = std::make_shared<TH1D>(
            "e#pi^{-}#pi^{+}X", "MMSQ: expecting proton", bins, 0.6, 1.3);

        makeHists_sector();
        makeHists_deltat();
        makeHists_MomVsBeta();
        makeHists_electron_cuts();
}

Histogram::~Histogram()
{
        this->Write();
}

void Histogram::Write()
{
        std::cout << GREEN << "Writting" << DEF << std::endl;

        std::cerr << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
        TDirectory *WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
        WvsQ2_folder->cd();
        Write_WvsQ2();

        std::cerr << BOLDBLUE << "Write_MomVsBeta()" << DEF << std::endl;
        TDirectory *Write_MomVsBeta_folder = RootOutputFile->mkdir("Mom Vs Beta");
        Write_MomVsBeta_folder->cd();
        Write_MomVsBeta();

        std::cerr << BOLDBLUE << "Write_Electron_cuts()" << DEF << std::endl;
        TDirectory *Electron_Cuts = RootOutputFile->mkdir("Electron_Cuts");
        Electron_Cuts->cd();
        Write_Electron_cuts();

        std::cerr << BOLDBLUE << "Write_Proton_cuts()" << DEF << std::endl;
        TDirectory *Proton_Cuts = RootOutputFile->mkdir("Proton_Cuts");
        Proton_Cuts->cd();
        Write_Proton_cuts();

        std::cerr << BOLDBLUE << "Write_Pip_cuts()" << DEF << std::endl;
        TDirectory *Pip_Cuts = RootOutputFile->mkdir("Pip_Cuts");
        Pip_Cuts->cd();
        Write_Pip_cuts();

        std::cerr << BOLDBLUE << "Write_Pim_cuts()" << DEF << std::endl;
        TDirectory *Pim_Cuts = RootOutputFile->mkdir("Pim_Cuts");
        Pim_Cuts->cd();
        Write_Pim_cuts();

        std::cerr << BOLDBLUE << "Write_deltat()" << DEF << std::endl;
        TDirectory *Write_deltat_folder = RootOutputFile->mkdir("Delta_t");
        Write_deltat_folder->cd();
        Write_deltat();

        std::cerr << BOLDBLUE << "Done Writing!!!" << DEF << std::endl;
}

void Histogram::Fill_WvsQ2(const std::shared_ptr<Reaction> &_e)
{
        short sec = _e->sec();
        TThread::Lock(); // Lock the thread to ensure exclusive access to the histograms
        {
                W_vs_q2->Fill(_e->W(), _e->Q2(), _e->weight());
                W_hist->Fill(_e->W(), _e->weight());
                Q2_hist->Fill(_e->Q2(), _e->weight());

                inv_mass_pPip->Fill(_e->inv_Ppip(), _e->weight());
                inv_mass_pPim->Fill(_e->inv_Ppim(), _e->weight());
                inv_mass_pipPim->Fill(_e->inv_pip_pim(), _e->weight());
                MM_twoPi_mPim->Fill(_e->MM_mPim(), _e->weight());
                MM2_twoPi_mPim->Fill(_e->MM2_mPim(), _e->weight());

                TThread::UnLock(); // Unlock after the operation
        }

        if (_e->TwoPion_exclusive())
        {
                missing_Energy_hist->Fill(_e->Energy_excl(), _e->weight());
                MM2_twoPi_excl->Fill(_e->MM2_exclusive(), _e->weight());
                MM_twoPi_excl->Fill(_e->MM_exclusive(), _e->weight());
        }
        if (_e->TwoPion_missingPip())
        {
                MM2_twoPi_missingPip->Fill(_e->MM2_mpip(), _e->weight());
        }
        if (_e->TwoPion_missingProt())
        {
                MM2_twoPi_missingProt->Fill(_e->MM2_mprot(), _e->weight());
        }

        if (sec > 0 && sec <= 6)
        {
                W_vs_q2_sec[sec - 1]->Fill(_e->W(), _e->Q2(), _e->weight());
                W_sec[sec - 1]->Fill(_e->W(), _e->weight());
        }
}
void Histogram::Fill_WvsQ2_twoPi_thrown(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<MCReaction> &_e)
{
        // short sec = _e->sec();

        W_vs_Q2_thrown->Fill(_e->W_mc(), _e->Q2_mc(), _e->weight());
        W_thrown->Fill(_e->W_mc(), _e->weight());
        Q2_thrown->Fill(_e->Q2_mc(), _e->weight());
}
void Histogram::Write_WvsQ2()
{
        W_thrown->SetXTitle("W_thrown (GeV)");
        if (W_thrown->GetEntries())
                W_thrown->Write();

        Q2_thrown->SetTitle("Q^{2} thrown ");
        Q2_thrown->SetXTitle("Q^{2} thrown (GeV^{2})");
        if (W_vs_Q2_thrown->GetEntries())
                Q2_thrown->Write();

        W_vs_Q2_thrown->SetTitle("Q^{2} versus W Thrown");
        W_vs_Q2_thrown->SetXTitle("W (GeV)");
        W_vs_Q2_thrown->SetYTitle("Q^{2} (GeV^{2})");
        W_vs_Q2_thrown->SetOption("COLZ1");
        if (W_vs_Q2_thrown->GetEntries())
                W_vs_Q2_thrown->Write();

        inv_mass_pPip->SetTitle("Inv Mass #Delta^{++} e(p,p'pi+X)e' events");
        inv_mass_pPip->SetXTitle("Mass (GeV)");
        inv_mass_pPip->Write();

        inv_mass_pPim->SetTitle("Inv Mass #Delta^{0} e(p,p'pi+X)e' events");
        inv_mass_pPim->SetXTitle("Mass (GeV)");
        inv_mass_pPim->Write();

        inv_mass_pipPim->SetTitle("Inv Mass #rho^{0} e(p,p'pi+X)e' events");
        inv_mass_pipPim->SetXTitle("Mass (GeV)");
        inv_mass_pipPim->Write();

        W_vs_q2->SetTitle("Q^{2} versus W Rec");
        W_vs_q2->SetXTitle("W (GeV)");
        W_vs_q2->SetYTitle("Q^{2} (GeV^{2})");
        W_vs_q2->SetOption("COLZ1");
        if (W_vs_q2->GetEntries())
                W_vs_q2->Write();

        W_hist->SetXTitle("W Distribution reconstructed");
        W_hist->SetXTitle("W (GeV)");
        if (W_hist->GetEntries())
                W_hist->Write();

        Q2_hist->SetXTitle("Q^{2} Distribution reconstructed");
        Q2_hist->SetXTitle("Q^{2} (GeV^{2})");
        if (Q2_hist->GetEntries())
                Q2_hist->Write();

        auto mmsq_4_topology = RootOutputFile->mkdir("mmsq_4_topology");
        mmsq_4_topology->cd();

        // excl

        MM2_twoPi_excl->SetTitle("MMSQ Excl e(p, p',pi+pi-X) Events");
        MM2_twoPi_excl->SetXTitle("MMSQ (GeV^{2}) ");
        if (MM2_twoPi_excl->GetEntries())
                MM2_twoPi_excl->Write();

        MM_twoPi_excl->SetTitle("Missing Mass Excl e(p, p',pi+pi-X) Events");
        MM_twoPi_excl->SetXTitle("MM (GeV) ");
        if (MM_twoPi_excl->GetEntries())
                MM_twoPi_excl->Write();

        // mPim
        MM2_twoPi_mPim->SetXTitle("MMSQ mPim (GeV^{2})");
        if (MM2_twoPi_mPim->GetEntries())
                MM2_twoPi_mPim->Write();

        MM_twoPi_mPim->SetXTitle("MM mPim (GeV^{2})");
        if (MM_twoPi_mPim->GetEntries())
                MM_twoPi_mPim->Write();

        // mPip (excl events)

        MM2_twoPi_missingPip->SetTitle("MMSQ mPip e(p, p',pi+pi-X) Events ");
        MM2_twoPi_missingPip->SetXTitle("MMSQ(GeV^{2}) ");

        if (MM2_twoPi_missingPip->GetEntries())
                MM2_twoPi_missingPip->Write();

        MM2_twoPi_missingProt->SetTitle("MMSQ mProt e(p, p',pi+pi-X) Events ");
        MM2_twoPi_missingProt->SetXTitle("MMSQ (GeV^{2}) ");
        if (MM2_twoPi_missingProt->GetEntries())
                MM2_twoPi_missingProt->Write();

        missing_Energy_hist->SetTitle("Missing energy e(p, p',pi+pi-X) Events");
        missing_Energy_hist->SetXTitle("E (GeV)");
        if (missing_Energy_hist->GetEntries())
                missing_Energy_hist->Write();

        auto wvsq2_sec = RootOutputFile->mkdir("wvsq2_sec");
        wvsq2_sec->cd();
        for (short i = 0; i < num_sectors; i++)
        {
                W_vs_q2_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
                W_vs_q2_sec[i]->SetXTitle("W (GeV)");
                W_vs_q2_sec[i]->SetOption("COLZ1");
                W_vs_q2_sec[i]->Write();
        }
        auto w_sec = RootOutputFile->mkdir("w_sec");
        w_sec->cd();
        for (short i = 0; i < num_sectors; i++)
        {
                W_sec[i]->SetXTitle("W (GeV)");
                W_sec[i]->Write();
        }

        ///////////////////////// numner of photoelectrons in HTCC //////////////////////////
        auto nphe_htcc_sec = RootOutputFile->mkdir("nphe_htcc_sec");
        nphe_htcc_sec->cd();

        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;

                for (short i = 0; i < num_sectors; i++)
                {
                        htcc_nphe_sec[c][i]->SetXTitle("Nphe");
                        htcc_nphe_sec[c][i]->SetYTitle("Count");

                        // if (htcc_nphe_sec[c][i]->GetEntries())
                        htcc_nphe_sec[c][i]->Write();
                }
        }
        ///////////////////////// //////////////////////////
        auto Vz_sec = RootOutputFile->mkdir("Vz_sec");
        Vz_sec->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;

                for (short i = 0; i < num_sectors; i++)
                {
                        vz_sec[c][i]->SetXTitle("Vz (cm)");
                        if (vz_sec[c][i]->GetEntries())
                                vz_sec[c][i]->Write();
                }
        }

        /////////////////////////  //////////////////////////
        auto SF_VS_MOM_sec = RootOutputFile->mkdir("SF_VS_MOM_sec");
        SF_VS_MOM_sec->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;

                for (short i = 0; i < num_sectors; i++)
                {
                        SF_VS_MOM[c][i]->SetOption("COLZ1");
                        SF_VS_MOM[c][i]->SetYTitle("SF ");
                        SF_VS_MOM[c][i]->SetXTitle("MOM (GeV)");
                        // if (SF_VS_MOM[i]->GetEntries())
                        SF_VS_MOM[c][i]->Write();
                }
        }
        /////////////////////////  //////////////////////////
        auto chi2pid_elec_sec = RootOutputFile->mkdir("chi2pid_elec_sec");
        chi2pid_elec_sec->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;

                for (short i = 0; i < num_sectors; i++)
                {
                        elec_Chi2pid_sec[c][i]->SetXTitle("Chi2pid");
                        elec_Chi2pid_sec[c][i]->SetYTitle("Count");

                        if (elec_Chi2pid_sec[c][i]->GetEntries())
                                elec_Chi2pid_sec[c][i]->Write();
                }
        }
}

void Histogram::makeHists_electron_cuts()
{
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;
                auto type = cut.second.c_str();
                EC_sampling_fraction[c] =
                    std::make_shared<TH2D>(Form("EC_sampling_fraction%s", type),
                                           Form("EC_sampling_fraction%s", type), bins, 1.5,
                                           10.0, bins, zero, 0.5);

                ECin_sf_vs_PCAL_sf[c] = std::make_shared<TH2D>(
                    Form("ECin_sf_vs_PCAL_sf%s", type), Form("ECin_sf_vs_PCAL_sf%s", type), bins,
                    -0., 0.25, bins, 0, 0.25);
                vz_position[c] = std::make_shared<TH1D>(Form("vz_position%s", type),
                                                        Form("vz_position%s", type), bins, -15, 15);
                momentum[c] = std::make_shared<TH1D>(Form("mom%s", type), Form("mom%s", type), bins, 0, 10);

                pcal_sec[c] =
                    std::make_shared<TH2D>(Form("pcal_sec%s", type), Form("pcal_sec%s", type), bins, -420, 420, bins, -420, 420);

                pcal_hx_hy_sec[c] =
                    std::make_shared<TH2D>(Form("pcal_hx_hy_sec%s", type), Form("pcal_hx_hy_sec%s", type), bins, -420, 420, bins, -420, 420);
                dcr1_sec[c] =
                    std::make_shared<TH2D>(Form("dcr1_sec%s", type), Form("dcr1_sec%s", type), bins, -180, 180, bins, -180, 180);
                dcr2_sec[c] =
                    std::make_shared<TH2D>(Form("dcr2_sec%s", type), Form("dcr2_sec%s", type), bins, -270, 270, bins, -270, 270);
                dcr3_sec[c] =
                    std::make_shared<TH2D>(Form("dcr3_sec%s", type), Form("dcr3_sec%s", type), bins, -320, 320, bins, -320, 320);

                prot_Delta_vz_cut_fd[c] = std::make_shared<TH1D>(Form("fd_prot_dvz_position%s", type),
                                                                 Form("fd_prot_dvz_position%s", type), bins, -40, 40);

                prot_Chi2pid_cut_fd[c] = std::make_shared<TH1D>(Form("fd_prot_chi2pid%s", type),
                                                                Form("fd_prot_chi2pid%s", type), bins, -20, 20);

                prot_Delta_vz_cut_cd[c] = std::make_shared<TH1D>(Form("cd_prot_dvz_position%s", type),
                                                                 Form("cd_prot_dvz_position%s", type), bins, -40, 40);

                prot_Chi2pid_cut_cd[c] = std::make_shared<TH1D>(Form("cd_prot_chi2pid%s", type),
                                                                Form("cd_prot_chi2pid%s", type), bins, -20, 20);

                pip_Delta_vz_cut_fd[c] = std::make_shared<TH1D>(Form("fd_pip_dvz_position%s", type),
                                                                Form("fd_pip_dvz_position%s", type), bins, -40, 40);

                pip_Chi2pid_cut_fd[c] = std::make_shared<TH1D>(Form("fd_pip_chi2pid%s", type),
                                                               Form("fd_pip_chi2pid%s", type), bins, -20, 20);

                pip_Delta_vz_cut_cd[c] = std::make_shared<TH1D>(Form("cd_pip_dvz_position%s", type),
                                                                Form("cd_pip_dvz_position%s", type), bins, -40, 40);

                pip_Chi2pid_cut_cd[c] = std::make_shared<TH1D>(Form("cd_pip_chi2pid%s", type),
                                                               Form("cd_pip_chi2pid%s", type), bins, -20, 20);

                pim_Delta_vz_cut[c] = std::make_shared<TH1D>(Form("pim_dvz_position%s", type),
                                                             Form("pim_dvz_position%s", type), bins, -40, 40);

                pim_Chi2pid_cut[c] = std::make_shared<TH1D>(Form("pim_chi2pid%s", type),
                                                            Form("pim_chi2pid%s", type), bins, -20, 20);

                dcr1_sec_prot[c] =
                    std::make_shared<TH2D>(Form("dcr1_sec_prot%s", type), Form("dcr1_sec_prot%s", type), bins, -180, 180, bins, -180, 180);
                dcr2_sec_prot[c] =
                    std::make_shared<TH2D>(Form("dcr2_sec_prot%s", type), Form("dcr2_sec_prot%s", type), bins, -270, 270, bins, -270, 270);
                dcr3_sec_prot[c] =
                    std::make_shared<TH2D>(Form("dcr3_sec_prot%s", type), Form("dcr3_sec_prot%s", type), bins, -420, 420, bins, -420, 420);

                dcr1_sec_pip[c] =
                    std::make_shared<TH2D>(Form("dcr1_sec_pip%s", type), Form("dcr1_sec_pip%s", type), bins, -180, 180, bins, -180, 180);
                dcr2_sec_pip[c] =
                    std::make_shared<TH2D>(Form("dcr2_sec_pip%s", type), Form("dcr2_sec_pip%s", type), bins, -270, 270, bins, -270, 270);
                dcr3_sec_pip[c] =
                    std::make_shared<TH2D>(Form("dcr3_sec_pip%s", type), Form("dcr3_sec_pip%s", type), bins, -420, 420, bins, -420, 420);

                dcr1_sec_pim[c] =
                    std::make_shared<TH2D>(Form("dcr1_sec_pim%s", type), Form("dcr1_sec_pim%s", type), bins, -180, 180, bins, -180, 180);
                dcr2_sec_pim[c] =
                    std::make_shared<TH2D>(Form("dcr2_sec_pim%s", type), Form("dcr2_sec_pim%s", type), bins, -270, 270, bins, -270, 270);
                dcr3_sec_pim[c] =
                    std::make_shared<TH2D>(Form("dcr3_sec_pim%s", type), Form("dcr3_sec_pim%s", type), bins, -360, 360, bins, -360, 360);
        }
}

void Histogram::FillHists_electron_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e)
{

        auto elec_cuts = std::make_shared<Pass2_Cuts>(_d);

        int sec = _e->sec();
        int ith_part = 0;
        if (_e->W() > 1.0 && _e->W() <= 3.0 && _e->Q2() > 1.5 && _e->Q2() <= 12.0)
        {

                vz_position[before_any_cuts]->Fill(_d->vz(0), _e->weight());
                if (_d->vz(0) > -10 && _d->vz(0) < 5)
                        vz_position[with_one_cut]->Fill(_d->vz(0), _e->weight());
                else
                        vz_position[outside_one_cut]->Fill(_d->vz(0), _e->weight());

                //// pcal
                pcal_sec[before_any_cuts]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0), _e->weight());
                if (elec_cuts->EC_hit_position_fiducial_cut_homogeneous())
                        pcal_sec[with_one_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0), _e->weight());
                else
                        pcal_sec[outside_one_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0), _e->weight());
                //// pcal hx hy
                pcal_hx_hy_sec[before_any_cuts]->Fill(_d->ec_pcal_hx(0), _d->ec_pcal_hy(0), _e->weight());
                if (elec_cuts->PCAL_fiducial_cut_HX_HY())
                        pcal_hx_hy_sec[with_one_cut]->Fill(_d->ec_pcal_hx(0), _d->ec_pcal_hy(0), _e->weight());
                else
                        pcal_hx_hy_sec[outside_one_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0), _e->weight());

                // dc
                dcr1_sec[before_any_cuts]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0), _e->weight());
                dcr2_sec[before_any_cuts]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0), _e->weight());
                dcr3_sec[before_any_cuts]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0), _e->weight());

                if (elec_cuts->DC_fiducial_cut_XY(ith_part, 0))
                {
                        dcr1_sec[with_one_cut]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0), _e->weight());
                        dcr2_sec[with_one_cut]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0), _e->weight());
                        dcr3_sec[with_one_cut]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0), _e->weight());
                }
                else
                {
                        dcr1_sec[outside_one_cut]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0), _e->weight());
                        dcr2_sec[outside_one_cut]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0), _e->weight());
                        dcr3_sec[outside_one_cut]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0), _e->weight());
                }

                momentum[before_any_cuts]->Fill(_d->p(0), _e->weight());
                if (_d->p(0) > 1.50)
                        momentum[with_one_cut]->Fill(_d->p(0), _e->weight());
                else
                        momentum[outside_one_cut]->Fill(_d->p(0), _e->weight());

                EC_sampling_fraction[before_any_cuts]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0), _e->weight());
                if (elec_cuts->EC_sampling_fraction_cut())

                        EC_sampling_fraction[with_one_cut]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0), _e->weight());
                else
                        EC_sampling_fraction[outside_one_cut]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0), _e->weight());

                //
                ECin_sf_vs_PCAL_sf[before_any_cuts]->Fill(((_d->ec_pcal_energy(0) / _d->p(0))), (_d->ec_ecin_energy(0) / _d->p(0)), _e->weight());
                if (elec_cuts->EC_sampling_fraction_cut())

                        ECin_sf_vs_PCAL_sf[with_one_cut]->Fill(((_d->ec_pcal_energy(0) / _d->p(0))), (_d->ec_ecin_energy(0) / _d->p(0)), _e->weight());
                else
                        ECin_sf_vs_PCAL_sf[outside_one_cut]->Fill(((_d->ec_pcal_energy(0) / _d->p(0))), (_d->ec_ecin_energy(0) / _d->p(0)), _e->weight());
                /////////
                if (sec > 0 && sec <= 6)
                {
                        /// nphe
                        htcc_nphe_sec[before_any_cuts][sec - 1]->Fill(_d->cc_htcc_nphe(0), _e->weight());
                        if (_d->cc_htcc_nphe(0) > 2)
                                htcc_nphe_sec[with_one_cut][sec - 1]->Fill(_d->cc_htcc_nphe(0), _e->weight());
                        else
                                htcc_nphe_sec[outside_one_cut][sec - 1]->Fill(_d->cc_htcc_nphe(0), _e->weight());
                        // chi2pid
                        elec_Chi2pid_sec[before_any_cuts][sec - 1]->Fill(_d->chi2pid(0), _e->weight());
                        if (_d->chi2pid(0) < 3)
                                elec_Chi2pid_sec[with_one_cut][sec - 1]->Fill(_d->chi2pid(0), _e->weight());
                        else
                                elec_Chi2pid_sec[outside_one_cut][sec - 1]->Fill(_d->chi2pid(0), _e->weight());
                        // vz cut
                        vz_sec[before_any_cuts][sec - 1]->Fill(_d->vz(0), _e->weight());
                        if (_d->vz(0) > -10 && _d->vz(0) < 5)
                                vz_sec[with_one_cut][sec - 1]->Fill(_d->vz(0), _e->weight());
                        else
                                vz_sec[outside_one_cut][sec - 1]->Fill(_d->vz(0), _e->weight());
                        ////////// EC/PCAL
                        SF_VS_MOM[before_any_cuts][sec - 1]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0), _e->weight());
                        if (elec_cuts->EC_sampling_fraction_cut())
                                SF_VS_MOM[with_one_cut][sec - 1]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0), _e->weight());
                        else
                                SF_VS_MOM[outside_one_cut][sec - 1]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0), _e->weight());
                }
        }
}

void Histogram::FillHists_electron_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e)
{
        int sec = _e->sec();

        if (_e->W() > 1.0 && _e->W() <= 3.0 && _e->Q2() > 1.5 && _e->Q2() <= 12.0)
        {

                vz_position[after_all_cuts]->Fill(_d->vz(0), _e->weight());
                pcal_sec[after_all_cuts]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0), _e->weight());
                pcal_hx_hy_sec[after_all_cuts]->Fill(_d->ec_pcal_hx(0), _d->ec_pcal_hy(0), _e->weight());
                dcr1_sec[after_all_cuts]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0), _e->weight());
                dcr2_sec[after_all_cuts]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0), _e->weight());
                dcr3_sec[after_all_cuts]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0), _e->weight());
                EC_sampling_fraction[after_all_cuts]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0), _e->weight());
                ECin_sf_vs_PCAL_sf[after_all_cuts]->Fill(((_d->ec_pcal_energy(0) / _d->p(0))), (_d->ec_ecin_energy(0) / _d->p(0)), _e->weight());

                momentum[after_all_cuts]->Fill(_d->p(0), _e->weight());

                if (sec > 0 && sec <= 6)
                {
                        // htcc_nph
                        htcc_nphe_sec[after_all_cuts][sec - 1]->Fill(_d->cc_htcc_nphe(0), _e->weight());
                        elec_Chi2pid_sec[after_all_cuts][sec - 1]->Fill(_d->chi2pid(0), _e->weight());
                        vz_sec[after_all_cuts][sec - 1]->Fill(_d->vz(0), _e->weight());
                        SF_VS_MOM[after_all_cuts][sec - 1]->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0), _e->weight());
                }
        }
}

void Histogram::FillHists_prot_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i)
{
        auto _cuts = std::make_shared<Pass2_Cuts>(_d);
        if (_e->W() > 1.0 && _e->W() <= 3.0 && _e->Q2() > 1.5 && _e->Q2() <= 12.0)
        {
                if (abs(_d->status(i)) < 4000)
                {
                        prot_Delta_vz_cut_fd[before_any_cuts]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());
                        prot_Chi2pid_cut_fd[before_any_cuts]->Fill(_d->chi2pid(i), _e->weight());

                        dcr1_sec_prot[before_any_cuts]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i), _e->weight());
                        dcr2_sec_prot[before_any_cuts]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i), _e->weight());
                        dcr3_sec_prot[before_any_cuts]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i), _e->weight());
                        if (_cuts->Hadron_Delta_vz_cut(i))
                                prot_Delta_vz_cut_fd[with_one_cut]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());
                        else
                                prot_Delta_vz_cut_fd[outside_one_cut]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());

                        if (_cuts->Hadron_Chi2pid_cut(i))
                                prot_Chi2pid_cut_fd[with_one_cut]->Fill(_d->chi2pid(i), _e->weight());
                        else
                                prot_Chi2pid_cut_fd[outside_one_cut]->Fill(_d->chi2pid(i), _e->weight());

                        if (_cuts->DC_fiducial_cut_XY(i, 1))
                        {
                                dcr1_sec_prot[with_one_cut]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i), _e->weight());
                                dcr2_sec_prot[with_one_cut]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i), _e->weight());
                                dcr3_sec_prot[with_one_cut]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i), _e->weight());
                        }
                        else
                        {

                                dcr1_sec_prot[outside_one_cut]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i), _e->weight());
                                dcr2_sec_prot[outside_one_cut]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i), _e->weight());
                                dcr3_sec_prot[outside_one_cut]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i), _e->weight());
                        }
                }
                else if (abs(_d->status(i)) >= 4000)
                {
                        prot_Delta_vz_cut_cd[before_any_cuts]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());
                        prot_Chi2pid_cut_cd[before_any_cuts]->Fill(_d->chi2pid(i), _e->weight());
                        if (_cuts->Hadron_Delta_vz_cut(i))
                                prot_Delta_vz_cut_cd[with_one_cut]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());
                        else
                                prot_Delta_vz_cut_cd[outside_one_cut]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());

                        if (_cuts->Hadron_Chi2pid_cut(i))
                                prot_Chi2pid_cut_cd[with_one_cut]->Fill(_d->chi2pid(i), _e->weight());
                        else
                                prot_Chi2pid_cut_cd[outside_one_cut]->Fill(_d->chi2pid(i), _e->weight());
                }
        }
}

void Histogram::FillHists_prot_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i, const TLorentzVector &prot)
{
        auto _cuts = std::make_shared<Pass2_Cuts>(_d);

        if (_e->W() > 1.0 && _e->W() <= 3.0 && _e->Q2() > 1.5 && _e->Q2() <= 12.0)
        {

                auto dt = std::make_shared<Delta_T>(_d);
                if (abs(_d->status(i)) < 4000)
                // if (dt->isCtof() == false)
                {
                        dcr1_sec_prot[after_all_cuts]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i), _e->weight());
                        dcr2_sec_prot[after_all_cuts]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i), _e->weight());
                        dcr3_sec_prot[after_all_cuts]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i), _e->weight());
                        prot_Delta_vz_cut_fd[after_all_cuts]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());
                        prot_Chi2pid_cut_fd[after_all_cuts]->Fill(_d->chi2pid(i), _e->weight());
                }
                else if (abs(_d->status(i)) > 4000)
                {

                        prot_Delta_vz_cut_cd[after_all_cuts]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());
                        prot_Chi2pid_cut_cd[after_all_cuts]->Fill(_d->chi2pid(i), _e->weight());
                }
        }
}
void Histogram::FillHists_pip_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i)
{
        auto _cuts = std::make_shared<Pass2_Cuts>(_d);

        if (_e->W() > 1.0 && _e->W() <= 3.0 && _e->Q2() > 1.5 && _e->Q2() <= 12.0)
        {

                if (abs(_d->status(i)) < 4000)
                {
                        pip_Delta_vz_cut_fd[before_any_cuts]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());
                        pip_Chi2pid_cut_fd[before_any_cuts]->Fill(_d->chi2pid(i), _e->weight());

                        dcr1_sec_pip[before_any_cuts]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i), _e->weight());
                        dcr2_sec_pip[before_any_cuts]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i), _e->weight());
                        dcr3_sec_pip[before_any_cuts]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i), _e->weight());

                        if (_cuts->Hadron_Delta_vz_cut(i))
                                pip_Delta_vz_cut_fd[with_one_cut]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());
                        else
                                pip_Delta_vz_cut_fd[outside_one_cut]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());

                        if (_cuts->Hadron_Chi2pid_cut(i))
                                pip_Chi2pid_cut_fd[with_one_cut]->Fill(_d->chi2pid(i), _e->weight());
                        else
                                pip_Chi2pid_cut_fd[outside_one_cut]->Fill(_d->chi2pid(i), _e->weight());
                }
                else if (abs(_d->status(i)) > 4000)
                {
                        pip_Delta_vz_cut_cd[before_any_cuts]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());
                        pip_Chi2pid_cut_cd[before_any_cuts]->Fill(_d->chi2pid(i), _e->weight());

                        if (_cuts->Hadron_Delta_vz_cut(i))
                                pip_Delta_vz_cut_cd[with_one_cut]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());
                        else
                                pip_Delta_vz_cut_cd[outside_one_cut]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());

                        if (_cuts->Hadron_Chi2pid_cut(i))
                                pip_Chi2pid_cut_cd[with_one_cut]->Fill(_d->chi2pid(i), _e->weight());
                        else
                                pip_Chi2pid_cut_cd[outside_one_cut]->Fill(_d->chi2pid(i), _e->weight());
                }

                if (_cuts->DC_fiducial_cut_XY(i, 2))
                {
                        dcr1_sec_pip[with_one_cut]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i), _e->weight());
                        dcr2_sec_pip[with_one_cut]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i), _e->weight());
                        dcr3_sec_pip[with_one_cut]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i), _e->weight());
                }
                else
                {
                        dcr1_sec_pip[outside_one_cut]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i), _e->weight());
                        dcr2_sec_pip[outside_one_cut]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i), _e->weight());
                        dcr3_sec_pip[outside_one_cut]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i), _e->weight());
                }
        }
}
void Histogram::FillHists_pip_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i, const TLorentzVector &pip)
{
        auto _cuts = std::make_shared<Pass2_Cuts>(_d);

        if (_e->W() > 1.0 && _e->W() <= 3.0 && _e->Q2() > 1.5 && _e->Q2() <= 12.0) // && _cuts->HadronsCuts(i))
        {
                auto dt = std::make_shared<Delta_T>(_d);

                if (abs(_d->status(i)) < 4000)
                // if (dt->isCtof() == false)
                {
                        pip_Delta_vz_cut_fd[after_all_cuts]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());
                        pip_Chi2pid_cut_fd[after_all_cuts]->Fill(_d->chi2pid(i), _e->weight());
                        dcr1_sec_pip[after_all_cuts]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i), _e->weight());
                        dcr2_sec_pip[after_all_cuts]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i), _e->weight());
                        dcr3_sec_pip[after_all_cuts]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i), _e->weight());
                }
                else if (abs(_d->status(i)) > 4000)
                // else if (dt->isCtof() == true)
                {
                        pip_Delta_vz_cut_cd[after_all_cuts]->Fill((_d->vz(i) - _d->vz(0)), _e->weight());
                        pip_Chi2pid_cut_cd[after_all_cuts]->Fill(_d->chi2pid(i), _e->weight());
                }
        }
}
void Histogram::FillHists_pim_pid_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i)
{
        pim_Delta_vz_cut[before_any_cuts]->Fill((_d->vz(0) - _d->vz(i)), _e->weight());
        pim_Chi2pid_cut[before_any_cuts]->Fill(_d->chi2pid(i), _e->weight());
        dcr1_sec_pim[before_any_cuts]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i), _e->weight());
        dcr2_sec_pim[before_any_cuts]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i), _e->weight());
        dcr3_sec_pim[before_any_cuts]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i), _e->weight());
}
void Histogram::FillHists_pim_pid_with_cuts(const std::shared_ptr<Branches12> &_d, const std::shared_ptr<Reaction> &_e, int i)
{
        pim_Delta_vz_cut[after_all_cuts]->Fill((_d->vz(0) - _d->vz(i)), _e->weight());
        pim_Chi2pid_cut[after_all_cuts]->Fill(_d->chi2pid(i), _e->weight());
        dcr1_sec_pim[after_all_cuts]->Fill(_d->dc_r1_x(i), _d->dc_r1_y(i), _e->weight());
        dcr2_sec_pim[after_all_cuts]->Fill(_d->dc_r2_x(i), _d->dc_r2_y(i), _e->weight());
        dcr3_sec_pim[after_all_cuts]->Fill(_d->dc_r3_x(i), _d->dc_r3_y(i), _e->weight());
}
void Histogram::Write_Electron_cuts()
{
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;

                vz_position[c]->SetXTitle("vz (cm)");
                if (vz_position[c]->GetEntries())
                        vz_position[c]->Write();
                pcal_sec[c]->SetXTitle("x (cm)");
                pcal_sec[c]->SetYTitle("y (cm)");
                pcal_sec[c]->SetOption("COLZ1");
                if (pcal_sec[c]->GetEntries())
                        pcal_sec[c]->Write();

                pcal_hx_hy_sec[c]->SetXTitle("x (cm)");
                pcal_hx_hy_sec[c]->SetYTitle("y (cm)");
                pcal_hx_hy_sec[c]->SetOption("COLZ1");
                if (pcal_hx_hy_sec[c]->GetEntries())
                        pcal_hx_hy_sec[c]->Write();

                dcr1_sec[c]->SetXTitle("x (cm)");
                dcr1_sec[c]->SetYTitle("y (cm)");
                dcr1_sec[c]->SetOption("COLZ1");
                if (dcr1_sec[c]->GetEntries())
                        dcr1_sec[c]->Write();

                dcr2_sec[c]->SetXTitle("x (cm)");
                dcr2_sec[c]->SetYTitle("y (cm)");
                dcr2_sec[c]->SetOption("COLZ1");
                if (dcr2_sec[c]->GetEntries())
                        dcr2_sec[c]->Write();

                dcr3_sec[c]->SetXTitle("x (cm)");
                dcr3_sec[c]->SetYTitle("y (cm)");
                dcr3_sec[c]->SetOption("COLZ1");
                if (dcr3_sec[c]->GetEntries())
                        dcr3_sec[c]->Write();

                EC_sampling_fraction[c]->SetXTitle("Momentum (GeV)");
                EC_sampling_fraction[c]->SetYTitle("Sampling Fraction");
                EC_sampling_fraction[c]->SetOption("COLZ1");
                if (EC_sampling_fraction[c]->GetEntries())
                        EC_sampling_fraction[c]->Write();

                ECin_sf_vs_PCAL_sf[c]->SetYTitle("ECin SF (GeV)");
                ECin_sf_vs_PCAL_sf[c]->SetXTitle("(0.2 - PCAL) SF");
                ECin_sf_vs_PCAL_sf[c]->SetOption("COLZ1");
                if (ECin_sf_vs_PCAL_sf[c]->GetEntries())
                        ECin_sf_vs_PCAL_sf[c]->Write();

                momentum[c]->SetXTitle("Momentum (GeV)");
                if (momentum[c]->GetEntries())
                        momentum[c]->Write();
        }
}

void Histogram::Write_Proton_cuts()
{
        auto proton_cuts_fd = RootOutputFile->mkdir("proton_cuts_fd");
        proton_cuts_fd->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;
                prot_Delta_vz_cut_fd[c]->SetXTitle("#Delta vz (cm)");
                if (prot_Delta_vz_cut_fd[c]->GetEntries())
                        prot_Delta_vz_cut_fd[c]->Write();

                prot_Chi2pid_cut_fd[c]->SetXTitle("Chi2pid");
                if (prot_Chi2pid_cut_fd[c]->GetEntries())
                        prot_Chi2pid_cut_fd[c]->Write();

                dcr1_sec_prot[c]->SetXTitle("x (cm)");
                dcr1_sec_prot[c]->SetYTitle("y (cm)");
                dcr1_sec_prot[c]->SetOption("COLZ1");
                if (dcr1_sec_prot[c]->GetEntries())
                        dcr1_sec_prot[c]->Write();

                dcr2_sec_prot[c]->SetXTitle("x (cm)");
                dcr2_sec_prot[c]->SetYTitle("y (cm)");
                dcr2_sec_prot[c]->SetOption("COLZ1");
                if (dcr2_sec_prot[c]->GetEntries())
                        dcr2_sec_prot[c]->Write();

                dcr3_sec_prot[c]->SetXTitle("x (cm)");
                dcr3_sec_prot[c]->SetYTitle("y (cm)");
                dcr3_sec_prot[c]->SetOption("COLZ1");
                if (dcr3_sec_prot[c]->GetEntries())
                        dcr3_sec_prot[c]->Write();
        }
        auto proton_cuts_cd = RootOutputFile->mkdir("proton_cuts_cd");
        proton_cuts_cd->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;
                prot_Delta_vz_cut_cd[c]->SetXTitle("#Delta vz (cm)");
                if (prot_Delta_vz_cut_cd[c]->GetEntries())
                        prot_Delta_vz_cut_cd[c]->Write();

                prot_Chi2pid_cut_cd[c]->SetXTitle("Chi2pid");
                if (prot_Chi2pid_cut_cd[c]->GetEntries())
                        prot_Chi2pid_cut_cd[c]->Write();
        }
}
void Histogram::Write_Pip_cuts()
{
        auto pip_cuts_fd = RootOutputFile->mkdir("pip_cuts_fd");
        pip_cuts_fd->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;
                pip_Delta_vz_cut_fd[c]->SetXTitle("#Delta vz (cm)");
                if (pip_Delta_vz_cut_fd[c]->GetEntries())
                        pip_Delta_vz_cut_fd[c]->Write();

                pip_Chi2pid_cut_fd[c]->SetXTitle("Chi2pid");
                if (pip_Chi2pid_cut_fd[c]->GetEntries())
                        pip_Chi2pid_cut_fd[c]->Write();

                dcr1_sec_pip[c]->SetXTitle("x (cm)");
                dcr1_sec_pip[c]->SetYTitle("y (cm)");
                dcr1_sec_pip[c]->SetOption("COLZ1");
                if (dcr1_sec_pip[c]->GetEntries())
                        dcr1_sec_pip[c]->Write();

                dcr2_sec_pip[c]->SetXTitle("x (cm)");
                dcr2_sec_pip[c]->SetYTitle("y (cm)");
                dcr2_sec_pip[c]->SetOption("COLZ1");
                if (dcr2_sec_pip[c]->GetEntries())
                        dcr2_sec_pip[c]->Write();

                dcr3_sec_pip[c]->SetXTitle("x (cm)");
                dcr3_sec_pip[c]->SetYTitle("y (cm)");
                dcr3_sec_pip[c]->SetOption("COLZ1");
                if (dcr3_sec_pip[c]->GetEntries())
                        dcr3_sec_pip[c]->Write();
        }
        auto pip_cuts_cd = RootOutputFile->mkdir("pip_cuts_cd");
        pip_cuts_cd->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;
                pip_Delta_vz_cut_cd[c]->SetXTitle("#Delta vz (cm)");
                if (pip_Delta_vz_cut_cd[c]->GetEntries())
                        pip_Delta_vz_cut_cd[c]->Write();

                pip_Chi2pid_cut_cd[c]->SetXTitle("Chi2pid");
                if (pip_Chi2pid_cut_cd[c]->GetEntries())
                        pip_Chi2pid_cut_cd[c]->Write();
        }
}

void Histogram::Write_Pim_cuts()
{

        auto pim_cuts = RootOutputFile->mkdir("pim_cuts");
        pim_cuts->cd();
        for (auto &&cut : before_after_cut)
        {
                int c = cut.first;
                pim_Delta_vz_cut[c]->SetXTitle("#Delta vz (cm)");
                if (pim_Delta_vz_cut[c]->GetEntries())
                        pim_Delta_vz_cut[c]->Write();

                pim_Chi2pid_cut[c]->SetXTitle("Chi2pid");
                if (pim_Chi2pid_cut[c]->GetEntries())
                        pim_Chi2pid_cut[c]->Write();

                dcr1_sec_pim[c]->SetXTitle("x (cm)");
                dcr1_sec_pim[c]->SetYTitle("y (cm)");
                dcr1_sec_pim[c]->SetOption("COLZ1");
                if (dcr1_sec_pim[c]->GetEntries())
                        dcr1_sec_pim[c]->Write();

                dcr2_sec_pim[c]->SetXTitle("x (cm)");
                dcr2_sec_pim[c]->SetYTitle("y (cm)");
                dcr2_sec_pim[c]->SetOption("COLZ1");
                if (dcr2_sec_pim[c]->GetEntries())
                        dcr2_sec_pim[c]->Write();

                dcr3_sec_pim[c]->SetXTitle("x (cm)");
                dcr3_sec_pim[c]->SetYTitle("y (cm)");
                dcr3_sec_pim[c]->SetOption("COLZ1");
                if (dcr3_sec_pim[c]->GetEntries())
                        dcr3_sec_pim[c]->Write();
        }
}

void Histogram::makeHists_sector()
{

        for (short i = 0; i < num_sectors; i++)
        {

                W_vs_q2_sec[i] = std::make_shared<TH2D>(
                    Form("wvsq2_sec_%d", i + 1), Form("W vs Q^{2} Sector: %d", i + 1), bins,
                    w_min, w_max, bins, zero, q2_max);

                W_sec[i] =
                    std::make_shared<TH1D>(Form("w_sec_%d", i + 1),
                                           Form("W Sector: %d", i + 1), bins, w_min, w_max);

                for (auto &&cut : before_after_cut)
                {
                        int c = cut.first;
                        auto type = cut.second.c_str();
                        // std::cout << " cut number = " << c << "  cut type is : " << type << std::endl;

                        htcc_nphe_sec[c][i] = std::make_shared<TH1D>(Form("htcc_nphe_sec%d%s", i + 1, type),
                                                                     Form("htcc nphe sec: %d%s", i + 1, type),
                                                                     200, 0, 70);
                        elec_Chi2pid_sec[c][i] = std::make_shared<TH1D>(Form("elec_Chi2pid_sec%d%s", i + 1, type),
                                                                        Form("elec_Chi2pid sec: %d%s", i + 1, type),
                                                                        200, -10, 5);
                        vz_sec[c][i] = std::make_shared<TH1D>(Form("vz_sec%d%s", i + 1, type),
                                                              Form("vz sec: %d%s", i + 1, type),
                                                              200, -20, 20);

                        SF_VS_MOM[c][i] = std::make_shared<TH2D>(
                            Form("SF_VS_MOM_%d%s", i + 1, type), Form("SF_VS_MOM_%d%s", i + 1, type), bins,
                            0, 10.0, bins, 0, 0.5);
                }
        }
}

void Histogram::makeHists_deltat()
{

        {
                std::string tof = "";
                static const short p_num = 3; // 0-P 1-Pip 2-Pim
                std::string p_name[p_num] = {"P", "pip", "Pim"};
                static const short cut_num = 2; // 0-no cuts, 1- with cuts
                std::string cut_name[cut_num] = {"before cut", "after cut"};

                for (short p = 0; p < p_num; p++)
                {
                        for (short c = 0; c < cut_num; c++)
                        {

                                tof = "both";
                                delta_t_hist[p][c][0] =
                                    std::make_shared<TH2D>(Form("delta_t_%s_%s_%s", p_name[p].c_str(),
                                                                cut_name[c].c_str(), tof.c_str()),
                                                           Form("#Deltat %s %s %s", p_name[p].c_str(),
                                                                cut_name[c].c_str(), tof.c_str()),
                                                           bins, p_min, p_max, bins, Dt_min, Dt_max);

                                tof = "ftof";
                                delta_t_hist[p][c][1] =
                                    std::make_shared<TH2D>(Form("delta_t_%s_%s_%s", p_name[p].c_str(),
                                                                cut_name[c].c_str(), tof.c_str()),
                                                           Form("#Deltat %s %s %s", p_name[p].c_str(),
                                                                cut_name[c].c_str(), tof.c_str()),
                                                           bins, p_min, p_max, bins, Dt_min, Dt_max);
                                tof = "ctof";
                                delta_t_hist[p][c][2] =
                                    std::make_shared<TH2D>(Form("delta_t_%s_%s_%s", p_name[p].c_str(),
                                                                cut_name[c].c_str(), tof.c_str()),
                                                           Form("#Deltat %s %s %s", p_name[p].c_str(),
                                                                cut_name[c].c_str(), tof.c_str()),
                                                           bins, p_min, p_max, bins, -6, 6);
                        }
                }
        }
}

void Histogram::Fill_deltat_before_cut(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt, int part,
                                       const std::shared_ptr<Reaction> &_e)
{
        if (_e->W() > 1.0 && _e->W() <= 3.0 && _e->Q2() > 1.5 && _e->Q2() <= 12.0)
        {
                int charge = data->charge(part);
                bool cd_part = (data->status(part) > 4000);
                bool fd_part = (data->status(part) > 2000 && data->status(part) < 4000);

                int pid = data->pid(part);
                float mom = data->p(part);
                float time[] = {NAN, NAN, NAN};
                float time_cd[] = {NAN, NAN, NAN};
                float time_fd[] = {NAN, NAN, NAN};

                if (fd_part)
                {

                        if (charge == 1)
                        {
                                // if (pid == PIP)
                                {
                                        time[1] = dt->dt_Pi();
                                        time_fd[1] = dt->dt_Pi();
                                        delta_t_hist[1][0][1]->Fill(mom, time_fd[1], _e->weight());
                                }

                                // else if (pid == PROTON)
                                {

                                        time[0] = dt->dt_P();
                                        time_fd[0] = dt->dt_P();
                                        delta_t_hist[0][0][1]->Fill(mom, time_fd[0], _e->weight());
                                }
                        }
                        else if (charge == -1)
                        {
                                // if (pid == PIM)
                                if (pid != ELECTRON)
                                {

                                        time[2] = dt->dt_Pi();
                                        time_fd[2] = dt->dt_Pi();
                                        delta_t_hist[2][0][1]->Fill(mom, time_fd[2], _e->weight());
                                }
                        }
                }
                else if (cd_part)
                {
                        if (charge == 1)
                        {
                                // if (pid == PROTON)
                                {

                                        time[0] = dt->dt_P();
                                        time_cd[0] = dt->dt_P();
                                        delta_t_hist[0][0][2]->Fill(mom, time_cd[0], _e->weight());
                                }

                                // else if (pid == PIP)
                                {
                                        time[1] = dt->dt_Pi();
                                        time_cd[1] = dt->dt_Pi();
                                        delta_t_hist[1][0][2]->Fill(mom, time_cd[1], _e->weight());
                                }
                        }
                        else if (charge == -1)
                        {
                                // if (pid == PIM)
                                if (pid != ELECTRON)
                                {
                                        time[2] = dt->dt_Pi();
                                        time_cd[2] = dt->dt_Pi();
                                        delta_t_hist[2][0][2]->Fill(mom, time_cd[2], _e->weight());
                                }
                        }
                }
                delta_t_hist[0][0][0]->Fill(mom, time[0], _e->weight());
                delta_t_hist[1][0][0]->Fill(mom, time[1], _e->weight());
                delta_t_hist[2][0][0]->Fill(mom, time[2], _e->weight());
        }
}
void Histogram::Fill_deltat_prot_after_cut(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt, int part,
                                           const std::shared_ptr<Reaction> &_e)
{
        if (_e->W() > 1.0 && _e->W() <= 3.0 && _e->Q2() > 1.5 && _e->Q2() <= 12.0)
        {
                int charge = data->charge(part);
                bool cd_part = (data->status(part) > 4000);
                bool fd_part = (data->status(part) > 2000 && data->status(part) < 4000);

                int pid = data->pid(part);
                float mom = data->p(part);
                float time = NAN;
                float time_cd = NAN;
                float time_fd = NAN;

                // if (pid == PROTON)
                {

                        if (charge == 1)

                        {
                                if (fd_part)
                                {
                                        time = dt->dt_P();
                                        time_fd = dt->dt_P();
                                        delta_t_hist[0][1][1]->Fill(mom, time_fd, _e->weight());
                                }
                                else if (cd_part)
                                {

                                        time = dt->dt_P();
                                        time_cd = dt->dt_P();
                                        delta_t_hist[0][1][2]->Fill(mom, time_cd, _e->weight());
                                }

                                delta_t_hist[0][1][0]->Fill(mom, time, _e->weight());
                        }
                }
        }
}
void Histogram::Fill_deltat_pip_after_cut(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt, int part,
                                          const std::shared_ptr<Reaction> &_e)
{
        if (_e->W() > 1.0 && _e->W() <= 3.0 && _e->Q2() > 1.5 && _e->Q2() <= 12.0)
        {
                int charge = data->charge(part);
                bool cd_part = (data->status(part) > 4000);
                bool fd_part = (data->status(part) > 2000 && data->status(part) < 4000);

                int pid = data->pid(part);
                float mom = data->p(part);
                float time = NAN;
                float time_cd = NAN;
                float time_fd = NAN;

                // if (pid == PIP)
                {

                        if (charge == 1)
                        {
                                if (fd_part)
                                {
                                        time = dt->dt_Pi();
                                        time_fd = dt->dt_Pi();
                                        delta_t_hist[1][1][1]->Fill(mom, time_fd, _e->weight());
                                }
                                else if (cd_part)
                                {

                                        time = dt->dt_Pi();
                                        time_cd = dt->dt_Pi();
                                        delta_t_hist[1][1][2]->Fill(mom, time_cd, _e->weight());
                                }

                                delta_t_hist[1][1][0]->Fill(mom, time, _e->weight());
                        }
                }
        }
}
void Histogram::Fill_deltat_pim_after_cut(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt, int part,
                                          const std::shared_ptr<Reaction> &_e)
{
        if (_e->W() > 1.0 && _e->W() <= 3.0 && _e->Q2() > 1.5 && _e->Q2() <= 12.0)
        {
                int charge = data->charge(part);
                bool cd_part = (data->status(part) > 4000);
                bool fd_part = (data->status(part) > 2000 && data->status(part) < 4000);

                int pid = data->pid(part);
                float mom = data->p(part);
                float time = NAN;
                float time_cd = NAN;
                float time_fd = NAN;

                // if (pid == PIM)
                {

                        if (charge == NEGATIVE)
                        {
                                if (fd_part)
                                {
                                        time = dt->dt_Pi();
                                        time_fd = dt->dt_Pi();
                                        delta_t_hist[2][1][1]->Fill(mom, time_fd, _e->weight());
                                }
                                else if (cd_part)
                                {

                                        time = dt->dt_Pi();
                                        time_cd = dt->dt_Pi();
                                        delta_t_hist[2][1][2]->Fill(mom, time_cd, _e->weight());
                                }

                                delta_t_hist[2][1][0]->Fill(mom, time, _e->weight());
                        }
                }
        }
}
void Histogram::Write_deltat()
{

        for (short i = 0; i < 3; i++)
        {
                for (short j = 0; j < 2; j++)
                {
                        for (short k = 0; k < 3; k++)
                        {
                                delta_t_hist[i][j][k]->SetXTitle("Momentum (GeV)");
                                delta_t_hist[i][j][k]->SetYTitle("#Deltat");
                                delta_t_hist[i][j][k]->SetOption("COLZ1");
                                // if (delta_t_hist[i][j][k]->GetEntries() > 1)
                                delta_t_hist[i][j][k]->Write();
                        }
                }
        }
}
void Histogram::makeHists_MomVsBeta()
{
        for (short p = 0; p < particle_num; p++)
        {
                for (short c = 0; c < charge_num; c++)
                {
                        for (short i = 0; i < with_id_num; i++)
                        {
                                momvsbeta_hist[p][c][i] = std::make_shared<TH2D>(
                                    Form("mom_vs_beta_%s_%s_%s", particle_name[p].c_str(), charge_name[c].c_str(), id_name[i].c_str()),
                                    Form("Momentum vs #beta %s %s %s", particle_name[p].c_str(), charge_name[c].c_str(), id_name[i].c_str()),
                                    bins, p_min, p_max, bins, zero, 1.2);
                        }
                }
        }
}
void Histogram::Fill_MomVsBeta(const std::shared_ptr<Branches12> &data, int part, const std::shared_ptr<Reaction> &_e)
{
        float beta = data->beta(part);
        float mom = data->p(part);
        int charge = data->charge(part);
        int pid = data->pid(part);
        if (beta != 0)
        {
                if (charge == NEGATIVE && pid != ELECTRON)
                {
                        momvsbeta_hist[0][1][0]->Fill(mom, beta, _e->weight());
                        if (pid == PIM)
                                momvsbeta_hist[0][1][1]->Fill(mom, beta, _e->weight());
                }
                else if (charge == POSITIVE)
                {
                        momvsbeta_hist[1][0][0]->Fill(mom, beta, _e->weight());
                        momvsbeta_hist[0][0][0]->Fill(mom, beta, _e->weight());

                        if (pid == PIP)
                                momvsbeta_hist[0][0][1]->Fill(mom, beta, _e->weight());

                        if (pid == PROTON)
                                momvsbeta_hist[1][0][1]->Fill(mom, beta, _e->weight());
                }
        }
}

void Histogram::Write_MomVsBeta()
{

        for (short p = 0; p < particle_num; p++)
        {
                for (short c = 0; c < charge_num; c++)
                {
                        for (short i = 0; i < with_id_num; i++)
                        {
                                momvsbeta_hist[p][c][i]->SetXTitle("Momentum (GeV)");
                                momvsbeta_hist[p][c][i]->SetYTitle("#beta");
                                momvsbeta_hist[p][c][i]->SetOption("COLZ1");
                                if (momvsbeta_hist[p][c][i]->GetEntries() > 1)
                                        momvsbeta_hist[p][c][i]->Write();
                        }
                }
        }
}
