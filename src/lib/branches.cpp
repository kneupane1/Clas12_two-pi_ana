
#include "branches.hpp"

Branches12::Branches12(const std::shared_ptr<TChain> &tree)
{
        _tree = tree;
        Branches12::init();
}

Branches12::Branches12(const std::shared_ptr<TChain> &tree, bool mc)
{
        _tree = tree;
        _is_mc = mc;
        Branches12::init();
        if (_is_mc)
                Branches12::initMC();
}

bool Branches12::mc()
{
        return _is_mc;
}

void Branches12::mc_branches()
{
        _is_mc = true;
        Branches12::initMC();
}

void Branches12::init()
{
        _pid = 0;
        _p = 0;
        _p2 = 0;
        _px = 0;
        _py = 0;
        _pz = 0;
        _vx = 0;
        _vy = 0;
        _vz = 0;
        //_vt = 0;
        _charge = 0;
        _beta = 0;
        _chi2pid = 0;
        _status = 0;
        _dc_sec = 0;
        _dc_r1_x = 0;
        _dc_r1_y = 0;
        _dc_r1_z = 0;
        _dc_r2_x = 0;
        _dc_r2_y = 0;
        _dc_r2_z = 0;
        _dc_r3_x = 0;
        _dc_r3_y = 0;
        _dc_r3_z = 0;
        _cvt_x = 0;
        _cvt_y = 0;
        _cvt_z = 0;
        _ec_tot_energy = 0;
        _ec_pcal_energy = 0;
        _ec_pcal_sec = 0;
        _ec_pcal_time = 0;
        _ec_pcal_path = 0;
        _ec_pcal_x = 0;
        _ec_pcal_y = 0;
        _ec_pcal_z = 0;
        _ec_pcal_hx = 0;
        _ec_pcal_hy = 0;
        _ec_pcal_hz = 0;
        _ec_pcal_lu = 0;
        _ec_pcal_lv = 0;
        _ec_pcal_lw = 0;

        _ec_ecin_energy = 0;
        _ec_ecin_sec = 0;
        _ec_ecin_time = 0;
        _ec_ecin_path = 0;
        _ec_ecin_x = 0;
        _ec_ecin_y = 0;
        _ec_ecin_z = 0;
        _ec_ecin_hx = 0;
        _ec_ecin_hy = 0;
        _ec_ecin_hz = 0;
        _ec_ecin_lu = 0;
        _ec_ecin_lv = 0;
        _ec_ecin_lw = 0;

        _ec_ecout_energy = 0;
        _ec_ecout_sec = 0;
        _ec_ecout_time = 0;
        _ec_ecout_path = 0;
        _ec_ecout_x = 0;
        _ec_ecout_y = 0;
        _ec_ecout_z = 0;
        _ec_ecout_hx = 0;
        _ec_ecout_hy = 0;
        _ec_ecout_hz = 0;
        _ec_ecout_lu = 0;
        _ec_ecout_lv = 0;
        _ec_ecout_lw = 0;

        _cc_nphe_tot = 0;
        _cc_ltcc_sec = 0;
        _cc_ltcc_nphe = 0;

        _cc_htcc_sec = 0;
        _cc_htcc_nphe = 0;

        _cc_rich_sec = 0;
        _cc_rich_nphe = 0;

        _sc_ftof_1a_sec = 0;
        _sc_ftof_1a_time = 0;
        _sc_ftof_1a_path = 0;

        _sc_ftof_1b_sec = 0;
        _sc_ftof_1b_time = 0;
        _sc_ftof_1b_path = 0;

        _sc_ftof_2_sec = 0;
        _sc_ftof_2_time = 0;
        _sc_ftof_2_path = 0;

        _sc_ctof_time = 0;
        _sc_ctof_path = 0;

        // Set branch addresses and branch pointers
        if (!_tree)
                return;

        _tree->SetMakeClass(1);

        _tree->SetBranchAddress("run", &_run, &b_run);
        _tree->SetBranchAddress("event", &_event, &b_event);
        _tree->SetBranchAddress("beamCharge", &_beamCharge, &b_beamCharge);
        _tree->SetBranchAddress("liveTime", &_liveTime, &b_liveTime);
        _tree->SetBranchAddress("startTime", &_startTime, &b_startTime);
        _tree->SetBranchAddress("helicity", &_helicity, &b_helicity);
        _tree->SetBranchAddress("pid", &_pid, &b_pid);
        _tree->SetBranchAddress("p", &_p, &b_p);
        _tree->SetBranchAddress("p2", &_p2, &b_p2);
        _tree->SetBranchAddress("px", &_px, &b_px);
        _tree->SetBranchAddress("py", &_py, &b_py);
        _tree->SetBranchAddress("pz", &_pz, &b_pz);
        _tree->SetBranchAddress("vx", &_vx, &b_vx);
        _tree->SetBranchAddress("vy", &_vy, &b_vy);
        _tree->SetBranchAddress("vz", &_vz, &b_vz);
        //_tree->SetBranchAddress("vt", &_vt, &b_vt);
        _tree->SetBranchAddress("charge", &_charge, &b_charge);
        _tree->SetBranchAddress("beta", &_beta, &b_beta);
        _tree->SetBranchAddress("chi2pid", &_chi2pid, &b_chi2pid);
        _tree->SetBranchAddress("status", &_status, &b_status);
        _tree->SetBranchAddress("dc_sec", &_dc_sec, &b_dc_sec);
        _tree->SetBranchAddress("dc_r1_x", &_dc_r1_x, &b_dc_r1_x);
        _tree->SetBranchAddress("dc_r1_y", &_dc_r1_y, &b_dc_r1_y);
        _tree->SetBranchAddress("dc_r1_z", &_dc_r1_z, &b_dc_r1_z);
        _tree->SetBranchAddress("dc_r2_x", &_dc_r2_x, &b_dc_r2_x);
        _tree->SetBranchAddress("dc_r2_y", &_dc_r2_y, &b_dc_r2_y);
        _tree->SetBranchAddress("dc_r2_z", &_dc_r2_z, &b_dc_r2_z);
        _tree->SetBranchAddress("dc_r3_x", &_dc_r3_x, &b_dc_r3_x);
        _tree->SetBranchAddress("dc_r3_y", &_dc_r3_y, &b_dc_r3_y);
        _tree->SetBranchAddress("dc_r3_z", &_dc_r3_z, &b_dc_r3_z);
        _tree->SetBranchAddress("cvt_x", &_cvt_x, &b_cvt_x);
        _tree->SetBranchAddress("cvt_y", &_cvt_y, &b_cvt_y);
        _tree->SetBranchAddress("cvt_z", &_cvt_z, &b_cvt_z);
        _tree->SetBranchAddress("ec_tot_energy", &_ec_tot_energy, &b_ec_tot_energy);
        _tree->SetBranchAddress("ec_pcal_energy", &_ec_pcal_energy, &b_ec_pcal_energy);
        _tree->SetBranchAddress("ec_pcal_sec", &_ec_pcal_sec, &b_ec_pcal_sec);
        _tree->SetBranchAddress("ec_pcal_time", &_ec_pcal_time, &b_ec_pcal_time);
        _tree->SetBranchAddress("ec_pcal_path", &_ec_pcal_path, &b_ec_pcal_path);
        _tree->SetBranchAddress("ec_pcal_x", &_ec_pcal_x, &b_ec_pcal_x);
        _tree->SetBranchAddress("ec_pcal_y", &_ec_pcal_y, &b_ec_pcal_y);
        _tree->SetBranchAddress("ec_pcal_z", &_ec_pcal_z, &b_ec_pcal_z);
        _tree->SetBranchAddress("ec_pcal_hx", &_ec_pcal_hx, &b_ec_pcal_hx);
        _tree->SetBranchAddress("ec_pcal_hy", &_ec_pcal_hy, &b_ec_pcal_hy);
        _tree->SetBranchAddress("ec_pcal_hz", &_ec_pcal_hz, &b_ec_pcal_hz);
        _tree->SetBranchAddress("ec_pcal_lu", &_ec_pcal_lu, &b_ec_pcal_lu);
        _tree->SetBranchAddress("ec_pcal_lv", &_ec_pcal_lv, &b_ec_pcal_lv);
        _tree->SetBranchAddress("ec_pcal_lw", &_ec_pcal_lw, &b_ec_pcal_lw);

        _tree->SetBranchAddress("ec_ecin_energy", &_ec_ecin_energy, &b_ec_ecin_energy);
        _tree->SetBranchAddress("ec_ecin_sec", &_ec_ecin_sec, &b_ec_ecin_sec);
        _tree->SetBranchAddress("ec_ecin_time", &_ec_ecin_time, &b_ec_ecin_time);
        _tree->SetBranchAddress("ec_ecin_path", &_ec_ecin_path, &b_ec_ecin_path);
        _tree->SetBranchAddress("ec_ecin_x", &_ec_ecin_x, &b_ec_ecin_x);
        _tree->SetBranchAddress("ec_ecin_y", &_ec_ecin_y, &b_ec_ecin_y);
        _tree->SetBranchAddress("ec_ecin_z", &_ec_ecin_z, &b_ec_ecin_z);
        _tree->SetBranchAddress("ec_ecin_hx", &_ec_ecin_hx, &b_ec_ecin_hx);
        _tree->SetBranchAddress("ec_ecin_hy", &_ec_ecin_hy, &b_ec_ecin_hy);
        _tree->SetBranchAddress("ec_ecin_hz", &_ec_ecin_hz, &b_ec_ecin_hz);
        _tree->SetBranchAddress("ec_ecin_lu", &_ec_ecin_lu, &b_ec_ecin_lu);
        _tree->SetBranchAddress("ec_ecin_lv", &_ec_ecin_lv, &b_ec_ecin_lv);
        _tree->SetBranchAddress("ec_ecin_lw", &_ec_ecin_lw, &b_ec_ecin_lw);

        _tree->SetBranchAddress("ec_ecout_energy", &_ec_ecout_energy, &b_ec_ecout_energy);
        _tree->SetBranchAddress("ec_ecout_sec", &_ec_ecout_sec, &b_ec_ecout_sec);
        _tree->SetBranchAddress("ec_ecout_time", &_ec_ecout_time, &b_ec_ecout_time);
        _tree->SetBranchAddress("ec_ecout_path", &_ec_ecout_path, &b_ec_ecout_path);
        _tree->SetBranchAddress("ec_ecout_x", &_ec_ecout_x, &b_ec_ecout_x);
        _tree->SetBranchAddress("ec_ecout_y", &_ec_ecout_y, &b_ec_ecout_y);
        _tree->SetBranchAddress("ec_ecout_z", &_ec_ecout_z, &b_ec_ecout_z);
        _tree->SetBranchAddress("ec_ecout_hx", &_ec_ecout_hx, &b_ec_ecout_hx);
        _tree->SetBranchAddress("ec_ecout_hy", &_ec_ecout_hy, &b_ec_ecout_hy);
        _tree->SetBranchAddress("ec_ecout_hz", &_ec_ecout_hz, &b_ec_ecout_hz);
        _tree->SetBranchAddress("ec_ecout_lu", &_ec_ecout_lu, &b_ec_ecout_lu);
        _tree->SetBranchAddress("ec_ecout_lv", &_ec_ecout_lv, &b_ec_ecout_lv);
        _tree->SetBranchAddress("ec_ecout_lw", &_ec_ecout_lw, &b_ec_ecout_lw);

        _tree->SetBranchAddress("cc_nphe_tot", &_cc_nphe_tot, &b_cc_nphe_tot);

        _tree->SetBranchAddress("cc_htcc_nphe", &_cc_htcc_nphe, &b_cc_htcc_nphe);

        _tree->SetBranchAddress("sc_ftof_1a_sec", &_sc_ftof_1a_sec, &b_sc_ftof_1a_sec);
        _tree->SetBranchAddress("sc_ftof_1a_time", &_sc_ftof_1a_time, &b_sc_ftof_1a_time);
        _tree->SetBranchAddress("sc_ftof_1a_path", &_sc_ftof_1a_path, &b_sc_ftof_1a_path);

        _tree->SetBranchAddress("sc_ftof_1b_sec", &_sc_ftof_1b_sec, &b_sc_ftof_1b_sec);
        _tree->SetBranchAddress("sc_ftof_1b_time", &_sc_ftof_1b_time, &b_sc_ftof_1b_time);
        _tree->SetBranchAddress("sc_ftof_1b_path", &_sc_ftof_1b_path, &b_sc_ftof_1b_path);

        _tree->SetBranchAddress("sc_ftof_2_sec", &_sc_ftof_2_sec, &b_sc_ftof_2_sec);
        _tree->SetBranchAddress("sc_ftof_2_time", &_sc_ftof_2_time, &b_sc_ftof_2_time);
        _tree->SetBranchAddress("sc_ftof_2_path", &_sc_ftof_2_path, &b_sc_ftof_2_path);

        _tree->SetBranchAddress("sc_ctof_time", &_sc_ctof_time, &b_sc_ctof_time);
        _tree->SetBranchAddress("sc_ctof_path", &_sc_ctof_path, &b_sc_ctof_path);
}

void Branches12::initMC()
{
        _mc_run = 0;
        _mc_event = 0;
        _mc_type = 0;
        _mc_helicity = 0;
        _mc_weight = 0;
        _mc_npart = 0;
        _mc_pid = 0;
        _mc_helicity = 0;
        _mc_px = 0;
        _mc_py = 0;
        _mc_pz = 0;
        _mc_vx = 0;
        _mc_vy = 0;
        _mc_vz = 0;
        _mc_vt = 0;

        _tree->SetBranchAddress("mc_helicity", &_mc_helicity);
        _tree->SetBranchAddress("mc_weight", &_mc_weight);
        _tree->SetBranchAddress("mc_npart", &_mc_npart);
        _tree->SetBranchAddress("mc_pid", &_mc_pid);
        _tree->SetBranchAddress("mc_px", &_mc_px);
        _tree->SetBranchAddress("mc_py", &_mc_py);
        _tree->SetBranchAddress("mc_pz", &_mc_pz);
        _tree->SetBranchAddress("mc_vx", &_mc_vx);
        _tree->SetBranchAddress("mc_vy", &_mc_vy);
        _tree->SetBranchAddress("mc_vz", &_mc_vz);
        _tree->SetBranchAddress("mc_vt", &_mc_vt);
}
// Add these getter methods
int Branches12::getRun() { return _run; }
int Branches12::getEvent() { return _event; }

int Branches12::gpart()
{
        return _pid->size();
}
int Branches12::pid(int i)
{
        if (i >= _pid->size())
                return -9999;
        else
                return _pid->at(i);
}
float Branches12::p(int i)
{
        if (i >= _pid->size())
                return NAN;
        else
                return _p->at(i);
}
float Branches12::p2(int i)
{
        if (i >= _pid->size())
                return NAN;
        else
                return _p2->at(i);
}
float Branches12::px(int i)
{
        if (i >= _pid->size())
                return NAN;
        else
                return _px->at(i);
}
float Branches12::py(int i)
{
        if (i >= _pid->size())
                return NAN;
        else
                return _py->at(i);
}
float Branches12::pz(int i)
{
        if (i >= _pid->size())
                return NAN;
        else
                return _pz->at(i);
}
float Branches12::vx(int i)
{
        if (i >= _pid->size())
                return NAN;
        else
                return _vx->at(i);
}
float Branches12::vy(int i)
{
        if (i >= _pid->size())
                return NAN;
        else
                return _vy->at(i);
}
float Branches12::vz(int i)
{
        if (i >= _pid->size())
                return NAN;
        else
                return _vz->at(i);
}
int Branches12::charge(int i)
{
        if (i >= _pid->size())
                return -9999;
        else
                return _charge->at(i);
}
float Branches12::beta(int i)
{
        if (i >= _pid->size())
                return NAN;
        else
                return _beta->at(i);
}

float Branches12::chi2pid(int i)
{
        if (i >= _pid->size())
                return NAN;
        else
                return _chi2pid->at(i);
}

int Branches12::status(int i)
{
        if (i >= _pid->size())
                return -9999;
        else
                return _status->at(i);
}

int Branches12::dc_sec(int i)
{
        if (i >= _dc_sec->size())
                return -9999;
        else
                return _dc_sec->at(i);
}
float Branches12::dc_r1_x(int i)
{
        if (i >= _dc_sec->size())
                return NAN;
        else
                return _dc_r1_x->at(i);
}
float Branches12::dc_r1_y(int i)
{
        if (i >= _dc_sec->size())
                return NAN;
        else
                return _dc_r1_y->at(i);
}
float Branches12::dc_r1_z(int i)
{
        if (i >= _dc_sec->size())
                return NAN;
        else
                return _dc_r1_z->at(i);
}
float Branches12::dc_r2_x(int i)
{
        if (i >= _dc_sec->size())
                return NAN;
        else
                return _dc_r2_x->at(i);
}
float Branches12::dc_r2_y(int i)
{
        if (i >= _dc_sec->size())
                return NAN;
        else
                return _dc_r2_y->at(i);
}
float Branches12::dc_r2_z(int i)
{
        if (i >= _dc_sec->size())
                return NAN;
        else
                return _dc_r2_z->at(i);
}
float Branches12::dc_r3_x(int i)
{
        if (i >= _dc_sec->size())
                return NAN;
        else
                return _dc_r3_x->at(i);
}
float Branches12::dc_r3_y(int i)
{
        if (i >= _dc_sec->size())
                return NAN;
        else
                return _dc_r3_y->at(i);
}
float Branches12::dc_r3_z(int i)
{
        if (i >= _dc_sec->size())
                return NAN;
        else
                return _dc_r3_z->at(i);
}

int Branches12::sc_ftof_1a_sec(int i)
{
        if (i >= _sc_ftof_1b_sec->size())
                return -9999;
        else
                return _sc_ftof_1a_sec->at(i);
}
float Branches12::sc_ftof_1a_time(int i)
{
        if (i >= _sc_ftof_1b_sec->size())
                return NAN;
        else
                return _sc_ftof_1a_time->at(i);
}
float Branches12::sc_ftof_1a_path(int i)
{
        if (i >= _sc_ftof_1b_sec->size())
                return NAN;
        else
                return _sc_ftof_1a_path->at(i);
}

int Branches12::sc_ftof_1b_sec(int i)
{
        if (i >= _sc_ftof_1b_sec->size())
                return -9999;
        else
                return _sc_ftof_1b_sec->at(i);
}
float Branches12::sc_ftof_1b_time(int i)
{
        if (i >= _sc_ftof_1b_sec->size())
                return NAN;
        else
                return _sc_ftof_1b_time->at(i);
}
float Branches12::sc_ftof_1b_path(int i)
{
        if (i >= _sc_ftof_1b_sec->size())
                return NAN;
        else
                return _sc_ftof_1b_path->at(i);
}

int Branches12::sc_ftof_2_sec(int i)
{
        if (i >= _sc_ftof_1b_sec->size())
                return -9999;
        else
                return _sc_ftof_2_sec->at(i);
}
float Branches12::sc_ftof_2_time(int i)
{
        if (i >= _sc_ftof_1b_sec->size())
                return NAN;
        else
                return _sc_ftof_2_time->at(i);
}
float Branches12::sc_ftof_2_path(int i)
{
        if (i >= _sc_ftof_1b_sec->size())
                return NAN;
        else
                return _sc_ftof_2_path->at(i);
}

float Branches12::sc_ctof_time(int i)
{
        if (i >= _sc_ftof_1b_sec->size())
                return NAN;
        else
                return _sc_ctof_time->at(i);
}
float Branches12::sc_ctof_path(int i)
{
        if (i >= _sc_ftof_1b_sec->size())
                return NAN;
        else
                return _sc_ctof_path->at(i);
}

int Branches12::mc_run()
{
        return _mc_run;
}
int Branches12::mc_event()
{
        return _mc_event;
}
int Branches12::mc_type()
{
        return _mc_type;
}
int Branches12::mc_helicity()
{
        return _mc_helicity;
}
float Branches12::mc_weight()
{
        return _mc_weight;
}
int Branches12::mc_npart()
{
        return _mc_npart;
}

int Branches12::mc_pid(int i)
{
        if (i >= _mc_npart)
                return -9999;
        else
                return _mc_pid->at(i);
}
float Branches12::mc_px(int i)
{
        if (i >= _mc_npart)
                return NAN;
        else
                return _mc_px->at(i);
}
float Branches12::mc_py(int i)
{
        if (i >= _mc_npart)
                return NAN;
        else
                return _mc_py->at(i);
}
float Branches12::mc_pz(int i)
{
        if (i >= _mc_npart)
                return NAN;
        else
                return _mc_pz->at(i);
}
float Branches12::mc_vx(int i)
{
        if (i >= _mc_npart)
                return NAN;
        else
                return _mc_vx->at(i);
}
float Branches12::mc_vy(int i)
{
        if (i >= _mc_npart)
                return NAN;
        else
                return _mc_vy->at(i);
}
float Branches12::mc_vz(int i)
{
        if (i >= _mc_npart)
                return NAN;
        else
                return _mc_vz->at(i);
}
float Branches12::mc_vt(int i)
{
        if (i >= _mc_npart)
                return NAN;
        else
                return _mc_vt->at(i);
}

float Branches12::ec_tot_energy(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_tot_energy->at(i);
}
float Branches12::ec_pcal_energy(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_pcal_energy->at(i);
}
int Branches12::ec_pcal_sec(int i)
{
        if (i >= _ec_ecin_sec->size())
                return -9999;
        else
                return _ec_pcal_sec->at(i);
};
float Branches12::ec_pcal_time(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_pcal_time->at(i);
}
float Branches12::ec_pcal_path(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_pcal_path->at(i);
}
float Branches12::ec_pcal_x(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_pcal_x->at(i);
}
float Branches12::ec_pcal_y(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_pcal_y->at(i);
}
float Branches12::ec_pcal_z(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_pcal_z->at(i);
}
float Branches12::ec_pcal_hx(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_pcal_hx->at(i);
}
float Branches12::ec_pcal_hy(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_pcal_hy->at(i);
}
float Branches12::ec_pcal_hz(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_pcal_hz->at(i);
}
float Branches12::ec_pcal_lu(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_pcal_lu->at(i);
}
float Branches12::ec_pcal_lv(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_pcal_lv->at(i);
}
float Branches12::ec_pcal_lw(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_pcal_lw->at(i);
}

float Branches12::ec_ecin_energy(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecin_energy->at(i);
}
int Branches12::ec_ecin_sec(int i)
{
        if (i >= _ec_ecin_sec->size())
                return -9999;
        else
                return _ec_ecin_sec->at(i);
}
float Branches12::ec_ecin_time(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecin_time->at(i);
}
float Branches12::ec_ecin_path(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecin_path->at(i);
}
float Branches12::ec_ecin_x(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecin_x->at(i);
}
float Branches12::ec_ecin_y(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecin_y->at(i);
}
float Branches12::ec_ecin_z(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecin_z->at(i);
}
float Branches12::ec_ecin_hx(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecin_hx->at(i);
}
float Branches12::ec_ecin_hy(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecin_hy->at(i);
}
float Branches12::ec_ecin_hz(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecin_hz->at(i);
}
float Branches12::ec_ecin_lu(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecin_lu->at(i);
}
float Branches12::ec_ecin_lv(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecin_lv->at(i);
}
float Branches12::ec_ecin_lw(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecin_lw->at(i);
}

float Branches12::ec_ecout_energy(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecout_energy->at(i);
}
int Branches12::ec_ecout_sec(int i)
{
        if (i >= _ec_ecin_sec->size())
                return -9999;
        else
                return _ec_ecout_sec->at(i);
}
float Branches12::ec_ecout_time(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecout_time->at(i);
}
float Branches12::ec_ecout_path(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecout_path->at(i);
}
float Branches12::ec_ecout_x(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecout_x->at(i);
}
float Branches12::ec_ecout_y(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecout_y->at(i);
}
float Branches12::ec_ecout_z(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecout_z->at(i);
}
float Branches12::ec_ecout_hx(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecout_hx->at(i);
}
float Branches12::ec_ecout_hy(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecout_hy->at(i);
}
float Branches12::ec_ecout_hz(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecout_hz->at(i);
}
float Branches12::ec_ecout_lu(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecout_lu->at(i);
}
float Branches12::ec_ecout_lv(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecout_lv->at(i);
}
float Branches12::ec_ecout_lw(int i)
{
        if (i >= _ec_ecin_sec->size())
                return NAN;
        else
                return _ec_ecout_lw->at(i);
}

float Branches12::cvt_x(int i)
{
        return _cvt_x->at(i);
}
float Branches12::cvt_y(int i)
{
        return _cvt_y->at(i);
}
float Branches12::cvt_z(int i)
{
        return _cvt_z->at(i);
}

float Branches12::cc_nphe_tot(int i)
{
        return _cc_nphe_tot->at(i);
}
int Branches12::cc_ltcc_sec(int i)
{
        return _cc_ltcc_sec->at(i);
}
float Branches12::cc_ltcc_nphe(int i)
{
        return _cc_ltcc_nphe->at(i);
}

int Branches12::cc_htcc_sec(int i)
{
        return _cc_htcc_sec->at(i);
}
float Branches12::cc_htcc_nphe(int i)
{
        return _cc_htcc_nphe->at(i);
}

int Branches12::cc_rich_sec(int i)
{
        return _cc_rich_sec->at(i);
}
float Branches12::cc_rich_nphe(int i)
{
        return _cc_rich_nphe->at(i);
}
