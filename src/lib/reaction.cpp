#include "reaction.hpp"

Reaction::Reaction(const std::shared_ptr<Branches12> &data, float beam_energy)
{
        _data = data;
        _beam = std::make_unique<TLorentzVector>();
        _beam_energy = 10.6041; //
        // _beam_energy = atof(getenv("CLAS12_E"));

        _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

        _gamma = std::make_unique<TLorentzVector>();
        _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
        _elecUnSmear = std::make_unique<TLorentzVector>();

        _elec = std::make_unique<TLorentzVector>();
        this->SetElec();

        _prot = std::vector<std::unique_ptr<TLorentzVector>>(); // Initialize _protons as an empty vector
        _pip = std::vector<std::unique_ptr<TLorentzVector>>();
        _pim = std::vector<std::unique_ptr<TLorentzVector>>();

        _prot_indices = std::vector<int>();
        _pip_indices = std::vector<int>();
        _pim_indices = std::vector<int>();

        _Energy_loss_uncorr_prot = std::make_unique<TLorentzVector>();
        _Energy_loss_uncorr_pip = std::make_unique<TLorentzVector>();
        _Energy_loss_uncorr_pim = std::make_unique<TLorentzVector>();

        _other = std::make_unique<TLorentzVector>();
        _neutron = std::make_unique<TLorentzVector>();

        _protUnSmear = std::make_unique<TLorentzVector>();
        _pipUnSmear = std::make_unique<TLorentzVector>();
        _pimUnSmear = std::make_unique<TLorentzVector>();
}

Reaction::~Reaction() {}

void Reaction::SetElec()
{
        _numPart++;
        _hasE = true;
        _sectorElec = _data->dc_sec(0);
        _elec_status = abs(_data->status(0));

        if (_mc)
        {

                _elecUnSmear->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);

                double W_unsmear = physics::W_calc(*_beam, *_elecUnSmear);

                double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear;

                pUnSmear = _elecUnSmear->P();

                thetaUnSmear = _elecUnSmear->Theta() * 180 / PI;

                if (_elecUnSmear->Phi() > 0)
                        phiUnSmear = _elecUnSmear->Phi() * 180 / PI;
                else if (_elecUnSmear->Phi() < 0)
                        phiUnSmear = (_elecUnSmear->Phi() + 2 * PI) * 180 / PI;

                ////////////////////////////////////////////////////////////////

                // Generate new values
                Reaction::SmearingFunc(ELECTRON, _elec_status, pUnSmear, thetaUnSmear, phiUnSmear, W_unsmear, pSmear, thetaSmear, phiSmear);

                _pxPrimeSmear = _elecUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
                _pyPrimeSmear = _elecUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
                _pzPrimeSmear =
                    _elecUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

                _elec->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_E); // smeared
                // _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);  // unsmeared

                *_gamma += *_beam - *_elec;

                // // // W and Q2 for sim data
                _W = physics::W_calc(*_beam, *_elec);
                _Q2 = physics::Q2_calc(*_beam, *_elec);

                _elec_mom = _elec->P();
                _E_elec = _elec->E();
                _theta_e = _elec->Theta() * 180 / PI;
        }
        else
        {
                _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);
                *_gamma += *_beam - *_elec;
                // // Exp W and Q2 here
                _W = physics::W_calc(*_beam, *_elec);
                _Q2 = physics::Q2_calc(*_beam, *_elec);
                _elec_mom = _elec->P();
                _E_elec = _elec->E();
                _theta_e = _elec->Theta() * 180 / PI;
        }
}

void Reaction::SetProton(int i)
{
        _numProt++;
        _numPos++;
        _hasP = true;
        auto proton = std::make_unique<TLorentzVector>();

        _prot_status = abs(_data->status(i));
        _Energy_loss_uncorr_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
        // _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);

        _sectorProt = _data->dc_sec(i);
        _prot_mom_uncorr = _Energy_loss_uncorr_prot->P();
        _prot_theta_uncorr = _Energy_loss_uncorr_prot->Theta() * 180 / PI;
        if (_Energy_loss_uncorr_prot->Phi() > 0)
                _prot_phi_uncorr = _Energy_loss_uncorr_prot->Phi() * 180 / PI;
        else if (_Energy_loss_uncorr_prot->Phi() < 0)
                _prot_phi_uncorr = (_Energy_loss_uncorr_prot->Phi() + 2 * PI) * 180 / PI;

        if (abs(_prot_status) >= 4000)
        {
                _prot_mom_tmt = _prot_mom_uncorr;
        }
        if ((2000 < abs(_prot_status)) && (abs(_prot_status) < 4000))
        {

                // // these are Andrey's corrections
                if (_prot_theta_uncorr < 27)
                {
                        _prot_mom_tmt = _prot_mom_uncorr + exp(-2.739 - 3.932 * _prot_theta_uncorr) + 0.002907;
                }
                else
                {
                        _prot_mom_tmt = _prot_mom_uncorr + exp(-1.2 - 4.228 * _prot_mom_uncorr) + 0.007502;
                }
        }

        _px_prime_prot_E = _data->px(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
        _py_prime_prot_E = _data->py(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
        _pz_prime_prot_E = _data->pz(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));

        if (_mc)
        {
                // /////////////////// SMEARING PART ////////////////////////////////////////////////////////////////////////////

                _protUnSmear->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P); // energy loss corrected

                //////////////////////////////////////////////////////////////
                double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, W_unsmear, pSmear, thetaSmear,
                    phiSmear;

                pUnSmear = _protUnSmear->P();

                thetaUnSmear = _protUnSmear->Theta() * 180 / PI;

                if (_protUnSmear->Phi() > 0)
                        phiUnSmear = _protUnSmear->Phi() * 180 / PI;
                else if (_protUnSmear->Phi() < 0)
                        phiUnSmear = (_protUnSmear->Phi() + 2 * PI) * 180 / PI;

                // Generate new values

                Reaction::SmearingFunc(PROTON, _prot_status, pUnSmear, thetaUnSmear, phiUnSmear, W_unsmear, pSmear, thetaSmear, phiSmear);

                _pxPrimeSmear = _protUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
                _pyPrimeSmear = _protUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
                _pzPrimeSmear =
                    _protUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);
                proton->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_P);

                _prot.push_back(std::move(proton));
                _prot_indices.push_back(i);
        }
        else
        {
                proton->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E,
                                MASS_P);
                _prot.push_back(std::move(proton));
                _prot_indices.push_back(i);
        }
}

void Reaction::SetPip(int i)
{
        _numPip++;
        _numPos++;
        _hasPip = true;
        auto pip = std::make_unique<TLorentzVector>();

        _pip_status = abs(_data->status(i));
        _sectorPip = _data->dc_sec(i);

        // /////////////////////////////////     SMEARING PART  /////////////////////////////
        if (_mc)
        {
                _pipUnSmear->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);

                double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, W_unsmear, pSmear, thetaSmear,
                    phiSmear;

                pUnSmear = _pipUnSmear->P();

                thetaUnSmear = _pipUnSmear->Theta() * 180 / PI;

                if (_pipUnSmear->Phi() > 0)
                        phiUnSmear = _pipUnSmear->Phi() * 180 / PI;
                else if (_pipUnSmear->Phi() < 0)
                        phiUnSmear = (_pipUnSmear->Phi() + 2 * PI) * 180 / PI;

                // Generate new values
                Reaction::SmearingFunc(PIP, _pip_status, pUnSmear, thetaUnSmear, phiUnSmear, W_unsmear, pSmear, thetaSmear, phiSmear);

                _pxPrimeSmear = _pipUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
                _pyPrimeSmear = _pipUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
                _pzPrimeSmear = _pipUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

                // _pipSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIP);  // smeared
                // _pip->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIP); // smeared

                pip->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIP); // smeared
                _pip.push_back(std::move(pip));
                _pip_indices.push_back(i); // Store the index
        }
        else
        {

                /////// vector method for many pip
                pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
                _pip.push_back(std::move(pip)); // Add pip to the vector
                _pip_indices.push_back(i);      // Store the index}
        }
}

void Reaction::SetPim(int i)
{
        _numPim++;
        _numNeg++;
        _hasPim = true;
        auto pim = std::make_unique<TLorentzVector>();

        _pim_status = abs(_data->status(i));
        _sectorPim = _data->dc_sec(i);

        _pim_mom_uncorr = _Energy_loss_uncorr_pim->P();
        _pim_theta_uncorr = _Energy_loss_uncorr_pim->Theta() * 180 / PI;
        if (_Energy_loss_uncorr_pim->Phi() > 0)
                _pim_phi_uncorr = _Energy_loss_uncorr_pim->Phi() * 180 / PI;
        else if (_Energy_loss_uncorr_pim->Phi() < 0)
                _pim_phi_uncorr = (_Energy_loss_uncorr_pim->Phi() + 2 * PI) * 180 / PI;

        // /////////////////////////////////     SMEARING PART  /////////////////////////////
        if (_mc)
        {
                _pimUnSmear->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);

                double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, W_unsmear, pSmear, thetaSmear,
                    phiSmear;

                pUnSmear = _pimUnSmear->P();

                thetaUnSmear = _pimUnSmear->Theta() * 180 / PI;

                if (_pimUnSmear->Phi() > 0)
                        phiUnSmear = _pimUnSmear->Phi() * 180 / PI;
                else if (_pimUnSmear->Phi() < 0)
                        phiUnSmear = (_pimUnSmear->Phi() + 2 * PI) * 180 / PI;

                // Generate new values
                Reaction::SmearingFunc(PIM, _pim_status, pUnSmear, thetaUnSmear, phiUnSmear, W_unsmear, pSmear, thetaSmear, phiSmear);

                _pxPrimeSmear = _pimUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
                _pyPrimeSmear = _pimUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                                sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
                _pzPrimeSmear = _pimUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

                pim->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIM); // smeared
                _pim.push_back(std::move(pim));                                      // Add pim to the vector
                _pim_indices.push_back(i);                                           // Store the index
        }
        else
        {
                pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
                _pim.push_back(std::move(pim)); // Add pim to the vector
                _pim_indices.push_back(i);      // Store the index
        }
}
void Reaction::SetNeutron(int i)
{
        _numNeutral++;
        _hasNeutron = true;
        _neutron->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_N);
}

void Reaction::SetOther(int i)
{
        if (_data->pid(i) == NEUTRON)
        {
                SetNeutron(i);
        }
        else
        {
                _numOther++;
                _hasOther = true;
                _other->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), mass[_data->pid(i)]);
        }
}

/////////////////// new added ////////////////
void Reaction::CalcMissMassPim(const TLorentzVector &prot, const TLorentzVector &pip)
{
        auto mm_mpim = std::make_unique<TLorentzVector>();

        *mm_mpim += (*_gamma + *_target);
        *mm_mpim -= prot;
        *mm_mpim -= pip;

        _MM_mPim = mm_mpim->M();
        _MM2_mPim = mm_mpim->M2();
}

void Reaction::CalcMissMassExcl(const TLorentzVector &prot, const TLorentzVector &pip, const TLorentzVector &pim)
{
        auto mm_excl = std::make_unique<TLorentzVector>();

        if (TwoPion_exclusive())
        {

                *mm_excl += (*_gamma + *_target);
                *mm_excl -= prot;
                *mm_excl -= pip;
                *mm_excl -= pim;
                _MM_exclusive = mm_excl->M();
                _MM2_exclusive = mm_excl->M2();
                _excl_Energy = mm_excl->E();
        }
}

void Reaction::CalcMissMassPip(const TLorentzVector &prot, const TLorentzVector &pim)
{
        auto mm_mpip = std::make_unique<TLorentzVector>();

        if (TwoPion_missingPip())
        {
                *mm_mpip += (*_gamma + *_target);
                *mm_mpip -= prot;
                *mm_mpip -= pim;
                _MM2_mpip = mm_mpip->M2();
        }
}

void Reaction::CalcMissMassProt(const TLorentzVector &pip, const TLorentzVector &pim)
{
        auto mm_mprot = std::make_unique<TLorentzVector>();

        if (TwoPion_missingProt())
        {
                *mm_mprot += (*_gamma + *_target);
                *mm_mprot -= pip;
                *mm_mprot -= pim;
                _MM2_mprot = mm_mprot->M2();
        }
}

float Reaction::MM_mPim()
{
        return _MM_mPim;
}

float Reaction::MM2_mPim()
{
        return _MM2_mPim;
}

float Reaction::MM2_exclusive()
{
        return _MM2_exclusive;
}
float Reaction::MM_exclusive()
{
        return _MM_exclusive;
}
float Reaction::MM2_mpip()
{
        return _MM2_mpip;
}
float Reaction::MM2_mprot()
{
        return _MM2_mprot;
}
float Reaction::Energy_excl()
{
        return _excl_Energy;
}

float Reaction::elec_momentum()
{
        if (TwoPion_missingPim())
                return _elec->P();
        else
                return NAN;
}

float Reaction::prot_momentum(const TLorentzVector &prot)
{
        // if (TwoPion_missingPim())
        if (_hasP)
                return prot.P();
        // return _prot->P();
        else
                return NAN;
}

float Reaction::pip_momentum(const TLorentzVector &pip)
{
        // if (TwoPion_missingPim())
        if (_hasPip)
                return pip.P();
        // return _pip->P();
        else
                return NAN;
}
float Reaction::pim_momentum(const TLorentzVector &prot, const TLorentzVector &pip)
{
        if (TwoPion_missingPim())
        {
                auto missingpim_ = std::make_unique<TLorentzVector>();
                // *missingpim_ += *_gamma + *_target - *_prot - *_pip;
                *missingpim_ += *_gamma + *_target - prot - pip;
                return missingpim_->P();
        }
        else
                return NAN;
}

float Reaction::pim_E(const TLorentzVector &prot, const TLorentzVector &pip)
{
        if (TwoPion_missingPim())
        {
                auto missingpim_ = std::make_unique<TLorentzVector>();
                // *missingpim_ += *_gamma + *_target - *_prot - *_pip;
                *missingpim_ += *_gamma + *_target - prot - pip;
                return missingpim_->E();
        }
        else
                return NAN;
}
float Reaction::pim_momentum_measured(const TLorentzVector &pim)
{
        if (TwoPion_exclusive())
                return pim.P();
        else
                return NAN;
}
float Reaction::pim_E_measured(const TLorentzVector &pim)
{
        if (TwoPion_exclusive())
                return pim.E();
        else
                return NAN;
}
float Reaction::theta_elec()
{ /// lab theta mattrai hunchha electron ko case ma
        // if (TwoPion_missingPim())
        if (_elec->Theta() > -500)
                return _elec->Theta() * 180.0 / PI;
        else
                return NAN;
}
float Reaction::Phi_elec()
{
        if (_elec->Phi() > -500)
        {
                if (_elec->Phi() > 0)
                        return _elec->Phi() * 180 / PI;
                else if (_elec->Phi() < 0)
                        return (_elec->Phi() + 2 * PI) * 180 / PI;
                else
                        return NAN;
        }
        else
                return NAN;
}
float Reaction::prot_theta_lab(const TLorentzVector &prot)
{
        // if (TwoPion_missingPim())
        if (_hasP)
                return prot.Theta() * 180.0 / PI;
        else
                return NAN;
}
float Reaction::pip_theta_lab(const TLorentzVector &pip)
{
        // if (TwoPion_missingPim())
        if (_hasPip)
                return pip.Theta() * 180.0 / PI;
        else
                return NAN;
}
float Reaction::pim_theta_lab(const TLorentzVector &prot, const TLorentzVector &pip)
{
        if (TwoPion_missingPim())
        {
                auto missingpim_ = std::make_unique<TLorentzVector>();
                *missingpim_ += *_gamma + *_target - prot - pip;
                return missingpim_->Theta() * 180.0 / PI;
        }
        else
                return NAN;
}
float Reaction::pim_theta_lab_measured(const TLorentzVector &pim)
{ /////////////////////////////////////work here

        if (TwoPion_exclusive())
                return pim.Theta() * 180.0 / PI;
        else
                return NAN;
}
float Reaction::prot_Phi_lab(const TLorentzVector &prot)
{
        // if (TwoPion_missingPim())
        // {
        if (_hasP)
        {
                // if (prot.Phi() > 0)
                return prot.Phi() * 180 / PI;
                // else if (prot.Phi() < 0)
                // return (prot.Phi() + 2 * PI) * 180 / PI;
                // else return NAN;
        }
        else
                return NAN;
}
float Reaction::prot_momT(const TLorentzVector &prot)
{
        // if (TwoPion_missingPim())
        if (_hasP)
                return prot.Perp();
        else
                return NAN;
}
float Reaction::pip_momT(const TLorentzVector &pip)
{
        // if (TwoPion_missingPim())
        if (_hasPip)
                return pip.Perp();
        else
                return NAN;
}
float Reaction::pim_momT(const TLorentzVector &pim)
{
        // if (TwoPion_missingPim())
        return pim.Perp();
        // else
        //         return NAN;
}
float Reaction::pip_Phi_lab(const TLorentzVector &pip)
{

        if (_hasPip)
        {
                // if (pip.Phi() > 0)
                return pip.Phi() * 180 / PI;
                // else if (pip.Phi() <= 0)
                //         return (pip.Phi() + 2 * PI) * 180 / PI;
        }
        else
                return NAN;
}
float Reaction::pim_Phi_lab(const TLorentzVector &prot, const TLorentzVector &pip)
{
        auto missingpim_ = std::make_unique<TLorentzVector>();
        *missingpim_ += *_gamma + *_target - prot - pip;
        if (TwoPion_missingPim())
        {
                // if (missingpim_->Phi() > 0)
                return missingpim_->Phi() * 180 / PI;
                // else if (missingpim_->Phi() <= 0)
                //         return (missingpim_->Phi() + 2 * PI) * 180 / PI;
        }
        else
                return NAN;
}

float Reaction::pim_Phi_lab_measured(const TLorentzVector &pim)
{
        if (TwoPion_exclusive())
        {
                // if (pim.Phi() > 0)
                return pim.Phi() * 180 / PI;
                // else if (pim.Phi() <= 0)
                //         return (pim.Phi() + 2 * PI) * 180 / PI;
        }
        else
                return NAN;
}

// // // // Calculate invariant masse
void Reaction::invMassPpim(const TLorentzVector &prot, const TLorentzVector &pip)
{
        auto delta0 = std::make_unique<TLorentzVector>();
        auto missingpim_ = std::make_unique<TLorentzVector>();
        *missingpim_ += *_gamma + *_target - prot - pip;
        *delta0 += prot + *missingpim_;
        _inv_Ppim = delta0->M();
}
void Reaction::invMasspippim(const TLorentzVector &prot, const TLorentzVector &pip)
{
        auto rho0 = std::make_unique<TLorentzVector>();
        auto missingpim_ = std::make_unique<TLorentzVector>();
        *missingpim_ += *_gamma + *_target - prot - pip;
        *rho0 += pip + *missingpim_;
        _inv_pip_pim = rho0->M();
}
void Reaction::invMassPpip(const TLorentzVector &prot, const TLorentzVector &pip)
{
        auto deltaPP = std::make_unique<TLorentzVector>();
        *deltaPP += prot + pip;
        _inv_Ppip = deltaPP->M();
}

float Reaction::inv_Ppip()
{
        return _inv_Ppip;
}

float Reaction::inv_Ppim()
{
        return _inv_Ppim;
}
float Reaction::inv_pip_pim()
{
        return _inv_pip_pim;
}

///////////////////////////////////////////////////////  MC CLAS //////////////////////////////////////////////////

MCReaction::MCReaction(const std::shared_ptr<Branches12> &data,
                       float beam_enrgy)
{
        _data = data;
        if (!_data->mc())
                _data->mc_branches();
        _beam_mc = std::make_unique<TLorentzVector>();
        _beam_energy = 10.6041;
        // _beam_energy = atof(getenv("CLAS12_E"));
        _weight_mc = _data->mc_weight();
        _beam_mc->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);
        _gamma_mc = std::make_unique<TLorentzVector>();
        _target = std::make_unique<TLorentzVector>();
        _target->SetXYZM(0.0, 0.0, 0.0, MASS_P);
        _elec_mc = std::make_unique<TLorentzVector>();
        this->SetMCElec();

        _other_mc = std::make_unique<TLorentzVector>();
        //_neutron = std::make_unique<TLorentzVector>();

        _prot_mc = std::vector<std::unique_ptr<TLorentzVector>>(); // Initialize _protons as an empty vector
        _pip_mc = std::vector<std::unique_ptr<TLorentzVector>>();
        _pim_mc = std::vector<std::unique_ptr<TLorentzVector>>();

        _prot_mc_indices = std::vector<int>();
        _pip_mc_indices = std::vector<int>();
        _pim_mc_indices = std::vector<int>();
}
void MCReaction::SetMCElec()
{
        _elec_mc->SetXYZM(_data->mc_px(0), _data->mc_py(0), _data->mc_pz(0), MASS_E);
        *_gamma_mc += *_beam_mc - *_elec_mc;

        // Can calculate W and Q2 for MC thrown data here
        _W_mc = physics::W_calc(*_beam_mc, *_elec_mc);
        _Q2_mc = physics::Q2_calc(*_beam_mc, *_elec_mc);
}

void MCReaction::SetMCProton(int i)
{
        auto proton_mc = std::make_unique<TLorentzVector>();
        proton_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_P);
        _prot_mc.push_back(std::move(proton_mc));
        _prot_mc_indices.push_back(i);
}

void MCReaction::SetMCPip(int i)
{
        auto pionp_mc = std::make_unique<TLorentzVector>();
        pionp_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIP);
        _pip_mc.push_back(std::move(pionp_mc));
        _pip_mc_indices.push_back(i);
}

void MCReaction::SetMCPim(int i)
{
        auto pionm_mc = std::make_unique<TLorentzVector>();
        pionm_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIM);
        _pim_mc.push_back(std::move(pionm_mc));
        _pim_mc_indices.push_back(i);
}
// void MCReaction::SetMCOther(int i) {
//   _other_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i),
//   mass[_data->pid(i)]);
// }

float MCReaction::pim_mom_mc_gen()
{
        TLorentzVector *pim_mc = _pim_mc[0].get();
        return pim_mc->P();
}
float MCReaction::pip_mom_mc_gen()
{
        TLorentzVector *pip_mc = _pip_mc[0].get();
        return pip_mc->P();
}
float MCReaction::prot_mom_mc_gen()
{
        TLorentzVector *prot_mc = _prot_mc[0].get();
        return prot_mc->P();
}

float MCReaction::MCpim_theta_lab()
{
        TLorentzVector *pim_mc = _pim_mc[0].get();
        return pim_mc->Theta() * 180 / PI;
}
float MCReaction::MCpip_theta_lab()
{
        TLorentzVector *pip_mc = _pip_mc[0].get();
        return pip_mc->Theta() * 180 / PI;
}
float MCReaction::MCprot_theta_lab()
{
        TLorentzVector *prot_mc = _prot_mc[0].get();
        return prot_mc->Theta() * 180 / PI;
}

float MCReaction::pim_phi_mc_gen()
{
        TLorentzVector *pim_mc = _pim_mc[0].get();
        if (pim_mc->Phi() >= 0)
                return (pim_mc->Phi() * 180 / PI);
        else if (pim_mc->Phi() < 0)
                return ((pim_mc->Phi() + 2 * PI) * 180 / PI);
        else
                return NAN;
}
float MCReaction::pip_phi_mc_gen()
{
        TLorentzVector *pip_mc = _pip_mc[0].get();
        if (pip_mc->Phi() >= 0)
                return (pip_mc->Phi() * 180 / PI);
        else if (pip_mc->Phi() < 0)
                return ((pip_mc->Phi() + 2 * PI) * 180 / PI);
        else
                return NAN;
}
float MCReaction::prot_phi_mc_gen()
{
        TLorentzVector *prot_mc = _prot_mc[0].get();
        if (prot_mc->Phi() >= 0)
                return (prot_mc->Phi() * 180 / PI);
        else if (prot_mc->Phi() < 0)
                return ((prot_mc->Phi() + 2 * PI) * 180 / PI);
        else
                return NAN;
}

// // // // // Calculate invariant masse
float MCReaction::MCinv_Ppip()
{
        auto deltaPPMC = std::make_unique<TLorentzVector>();
        TLorentzVector *prot_mc = _prot_mc[0].get();
        TLorentzVector *pip_mc = _pip_mc[0].get();
        *deltaPPMC += *prot_mc + *pip_mc;
        return deltaPPMC->M();
}
float MCReaction::MCinv_Ppim()
{
        auto delta0MC = std::make_unique<TLorentzVector>();
        TLorentzVector *prot_mc = _prot_mc[0].get();
        TLorentzVector *pim_mc = _pim_mc[0].get();
        *delta0MC += *prot_mc;
        *delta0MC += *pim_mc;
        return delta0MC->M();
}
float MCReaction::MCinv_pip_pim()
{
        auto rho0MC = std::make_unique<TLorentzVector>();
        TLorentzVector *pip_mc = _pip_mc[0].get();
        TLorentzVector *pim_mc = _pim_mc[0].get();
        *rho0MC += *pip_mc + *pim_mc;
        return rho0MC->M();
}
