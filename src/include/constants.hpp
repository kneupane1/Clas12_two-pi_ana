
#ifndef CONSTANTS_H_GUARD
#define CONSTANTS_H_GUARD
#include "TMath.h"
#include <TROOT.h>
#include <map>
#include <unordered_map>

////////////////////////////////////////////
// static const bool _mc = true;
// static const bool _mc = false;
extern bool _mc; // Declare as external variable

///////////////////////////////////////////

static const int MAX_PARTS = 100;
static const int N_SIGMA = 3;
static const double PI = TMath::Pi();
static const float INV_SQRT_2PI = TMath::Sqrt(2 * TMath::Pi());
static const double DEG2RAD = PI / 180.0;
static const double RAD2DEG = 180.0 / PI;
static const int POSITIVE = 1;
static const int NEGATIVE = -1;

static const double c_special_units = 29.9792458; // speed of light in centimeters per nanosecond (cm/ns).

// particle codes, usually PDG codes, but always those used in BOS
static const int PROTON = 2212;
static const int NEUTRON = 2112;
static const int PIP = 211;
static const int PIM = -211;
static const int PI0 = 111;
static const int KP = 321;
static const int KM = -321;
static const int PHOTON = 22;
static const int ELECTRON = 11;

// PDG particle masses in GeV/c2
static const double MASS_P = 0.93827203;
static const double MASS_N = 0.93956556;
static const double MASS_E = 0.000511;
static const double MASS_PIP = 0.13957018;
static const double MASS_PIM = 0.13957018;
static const double MASS_PI0 = 0.1349766;
static const double MASS_KP = 0.493677;
static const double MASS_KM = 0.493677;
static const double MASS_G = 0.0;
static const double MASS_OMEGA = 0.78265;

static std::unordered_map<int, double> mass = {
    {PROTON, MASS_P}, {-PROTON, MASS_P}, {NEUTRON, MASS_N}, {PIP, MASS_PIP}, {PIM, MASS_PIM}, {PI0, MASS_PI0}, {KP, MASS_KP}, {KM, MASS_KM}, {PHOTON, MASS_G}, {ELECTRON, MASS_E}, {-ELECTRON, MASS_E}};

static std::map<int, std::string> before_after_cut = {{0, "_before_any_cuts"}, {1, "_with_one_cut"}, {2, "_outside_one_cut"}, {3, "_after_all_cuts"}};

#endif
