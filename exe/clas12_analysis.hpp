#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "histogram.hpp"
#include "reaction.hpp"

template <class CutType>
size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<Histogram> &_hists, int thread_id)
{
        // Get the number of events in this thread
        size_t num_of_events = (int)_chain->GetEntries();

        float beam_energy = 10.6041;

        if (std::is_same<CutType, Pass2_Cuts>::value)
        {
                if (thread_id == 0)
                        std::cout << BLUE << "Using Pass2 RGA Cuts" << DEF << std::endl;

                beam_energy = 10.6041;
        }

        if (getenv("CLAS12_E") != NULL)
                beam_energy = atof(getenv("CLAS12_E"));

        // Print some information for each thread
        std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
                  << num_of_events << " Events " << DEF << "===============\n";

        // // // Make a data object which all the branches can be accessed from
        auto data = std::make_shared<Branches12>(_chain, _mc);

        // Total number of events "Processed"
        float no_of_events = 0;
        int prot = 0, pip = 0, pim = 0, elec = 0;

        // For each event
        for (size_t current_event = 0; current_event < num_of_events; current_event++)
        {

                // Get current event
                _chain->GetEntry(current_event);
                // If we are the 0th thread print the progress of the thread every 1000 events
                if (thread_id == 0 && current_event % 1000 == 0)
                        std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;

                ///////// This part is only for mc events /////////////
                if (_mc)
                {
                        if (data->mc_weight() <= 0)
                                if (data->mc_npart() < 1 || data->mc_weight() <= 0)
                                        continue;

                        // Make a reaction class from the data given
                        auto mc_event = std::make_shared<MCReaction>(data, beam_energy);

                        for (int part = 0; part < data->mc_npart(); part++)
                        {
                                // Check particle ID's and fill the reaction class
                                if (data->mc_pid(part) == PIP)
                                {
                                        mc_event->SetMCPip(part);
                                }
                                if (data->mc_pid(part) == PROTON)
                                {
                                        mc_event->SetMCProton(part);
                                }
                                if (data->mc_pid(part) == PIM)
                                {
                                        mc_event->SetMCPim(part);
                                        // } else {
                                        //   mc_event->SetMCOther(part);
                                }
                        }

                        if (mc_event->W_mc() < 3.0 && mc_event->W_mc() > 1.0 && mc_event->Q2_mc() < 12.0 && mc_event->Q2_mc() > 1.0)
                        {
                                _hists->Fill_WvsQ2_twoPi_thrown(data, mc_event);
                        }
                }

                ///////// This part is only for Rec events both mc and exp/////////////
                // Class reaction
                auto event = std::make_shared<Reaction>(data, beam_energy);
                no_of_events++;

                //// one basic cut to plot before electron-cuts
                if (data->charge(0) == NEGATIVE)
                        _hists->FillHists_electron_cuts(data, event);

                /// Class cut
                auto cuts = std::make_unique<Pass2_Cuts>(data);

                if (!cuts->ElectronCuts())
                        continue;

                // If we pass electron cuts the event is processed
                elec++;

                /// Class hist
                _hists->FillHists_electron_with_cuts(data, event);

                int prot_idx = -1;
                int pip_idx = -1;

                // CLass Deltat
                auto dt = std::make_shared<Delta_T>(data);
                auto dt_proton = std::make_shared<Delta_T>(data);
                auto dt_pip = std::make_shared<Delta_T>(data);

                // ///////////////////////////////////  Particle loop  ///////////////////////////

                for (int part = 1; part < data->gpart(); part++)
                {
                        dt->dt_calc(part);

                        if ((data->charge(part) != ELECTRON) && (data->charge(part) != 0))
                        {
                                _hists->Fill_MomVsBeta(data, part, event);
                                _hists->Fill_deltat_before_cut(data, dt, part, event);

                                if (data->charge(part) == POSITIVE)
                                {
                                        _hists->FillHists_prot_pid_cuts(data, event, part);
                                        _hists->FillHists_pip_pid_cuts(data, event, part);

                                        if (cuts->IsProton(part))
                                        {

                                                prot++;
                                                event->SetProton(part);
                                                prot_idx++;
                                                _hists->Fill_deltat_prot_after_cut(data, dt, part, event);
                                                _hists->FillHists_prot_pid_with_cuts(data, event, part, *event->GetProtons()[prot_idx]);
                                        }

                                        else if (cuts->IsPip(part))
                                        {

                                                pip++;
                                                event->SetPip(part);
                                                pip_idx++;
                                                _hists->Fill_deltat_pip_after_cut(data, dt, part, event);
                                                _hists->FillHists_pip_pid_with_cuts(data, event, part, *event->GetPips()[pip_idx]);
                                        }
                                }
                                else
                                {
                                        _hists->FillHists_pim_pid_cuts(data, event, part);

                                        if (cuts->IsPim(part))
                                        {
                                                event->SetPim(part);
                                                pim++;
                                                _hists->Fill_deltat_pim_after_cut(data, dt, part, event);
                                                _hists->FillHists_pim_pid_with_cuts(data, event, part);
                                        }
                                }
                        }
                }

                ///////////////////////////////////  Now start cutting on our kinamatics and selecting twoPion events for analysis //////////////
                if (event->W() > 1. && event->W() <= 3.0 && event->Q2() <= 12.0 && event->Q2() >= 1.)
                {
                        {
                                // _hists->Fill_WvsQ2(event);

                                if (event->TwoPion_exclusive())
                                {
                                        for (size_t i = 0; i < event->GetProtons().size(); ++i)
                                        {
                                                for (size_t j = 0; j < event->GetPips().size(); ++j)
                                                {
                                                        for (size_t k = 0; k < event->GetPims().size(); ++k)
                                                        {
                                                                event->CalcMissMassExcl(*event->GetProtons()[i], *event->GetPips()[j], *event->GetPims()[k]);
                                                        }
                                                }
                                        }
                                }

                                if (event->TwoPion_missingPip())
                                {
                                        for (size_t i = 0; i < event->GetProtons().size(); ++i)
                                        {
                                                for (size_t k = 0; k < event->GetPims().size(); ++k)
                                                {
                                                        event->CalcMissMassPip(*event->GetProtons()[i], *event->GetPims()[k]);
                                                }
                                        }
                                }

                                if (event->TwoPion_missingProt())
                                {
                                        for (size_t j = 0; j < event->GetPips().size(); ++j)
                                        {
                                                for (size_t k = 0; k < event->GetPims().size(); ++k)
                                                {
                                                        event->CalcMissMassProt(*event->GetPips()[j], *event->GetPims()[k]);
                                                }
                                        }
                                }
                                // }

                                if (event->TwoPion_missingPim())
                                { // // Retrieve the number of protons and pions in the event
                                        for (size_t i = 0; i < event->GetProtons().size(); ++i)
                                        {
                                                for (size_t j = 0; j < event->GetPips().size(); ++j)
                                                {
                                                        {
                                                                event->CalcMissMassPim(*event->GetProtons()[i], *event->GetPips()[j]);
                                                                event->invMassPpip(*event->GetProtons()[i], *event->GetPips()[j]);
                                                                event->invMassPpim(*event->GetProtons()[i], *event->GetPips()[j]);
                                                                event->invMasspippim(*event->GetProtons()[i], *event->GetPips()[j]);
                                                        }
                                                }
                                        }
                                }
                        }

                        _hists->Fill_WvsQ2(event);
                }
        }
        // std::cout.precision(3);

        std::cout << "Percent = " << 100.0 * num_of_events / num_of_events << std::endl;
        std::cout << "   no of total events  " << num_of_events << std::endl;
        std::cout << " elec " << elec << " prot " << prot << " pip " << pip << " pim " << pim << '\n';
        // Return the total number of events
        return num_of_events;
}
#endif
