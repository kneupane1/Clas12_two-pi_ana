
#include "deltat.hpp"

Delta_T::Delta_T(std::shared_ptr<Branches12> data)
{
  _data = data;
  if (!std::isnan(_data->sc_ftof_1b_time(0)))
  {
    _sc_t_v = _data->sc_ftof_1b_time(0);
    _sc_r_v = _data->sc_ftof_1b_path(0);
  }
  else if (!std::isnan(_data->sc_ftof_1a_time(0)))
  {
    _sc_t_v = _data->sc_ftof_1a_time(0);
    _sc_r_v = _data->sc_ftof_1a_path(0);
  }
  else if (!std::isnan(_data->sc_ftof_2_time(0)))
  {
    _sc_t_v = _data->sc_ftof_2_time(0);
    _sc_r_v = _data->sc_ftof_2_path(0);
  }
  _vertex = _vertex_time(_sc_t_v, _sc_r_v, 1.0);

  if (!std::isnan(_data->sc_ctof_time(0)))
  {
    _ctof = true;
    _ctof_t_v = _data->sc_ctof_time(0);
    _ctof_r_v = _data->sc_ctof_path(0);
    _ctof_vertex = _vertex_time(_ctof_t_v, _ctof_r_v, 1.0);
  }
  else
  {
    _ctof_vertex = _vertex;
  }
}

Delta_T::~Delta_T() {}

float Delta_T::_vertex_time(float sc_time, float sc_pathlength,
                            float relatavistic_beta)
{
  if (std::isnan(sc_time) || std::isnan(sc_pathlength))
    return NAN;
  return sc_time - sc_pathlength / (relatavistic_beta * c_special_units);
}

float Delta_T::_deltat(int pid)
{

  _beta = 1.0 / sqrt(1.0 + (mass[pid] / _momentum) * (mass[pid] / _momentum));

  // std::cout << "  status of ctof particle at  _deltat :  " << _data->status(pid) << "  _sc_t  is " << _sc_t << std::endl;

  if (_sc_t == _sc_t && _sc_r == _sc_r)
  {
    return _vertex - _vertex_time(_sc_t, _sc_r, _beta);
  }
  else
  {
    return NAN;
  }
}

float Delta_T::_ctof_deltat(int pid)
{
  _beta = 1.0 / sqrt(1.0 + (mass[pid] / _momentum) * (mass[pid] / _momentum));
  float _v = _vertex;
  if (_ctof)
    _v = _ctof_vertex;
  if (_ctof_t == _ctof_t && _ctof_r == _ctof_r)
  {
    return _v - _vertex_time(_ctof_t, _ctof_r, _beta);
  }
  else
  {
    return NAN;
  }
}

void Delta_T::dt_calc(int i)
{
  if (i > _data->gpart())
    return;
  _momentum = _data->p(i);
  if (!std::isnan(_data->sc_ftof_1b_time(i)))
  {
    _sc_t = _data->sc_ftof_1b_time(i);
    _sc_r = _data->sc_ftof_1b_path(i);
  }
  else if (!std::isnan(_data->sc_ftof_1a_time(i)))
  {

    _sc_t = _data->sc_ftof_1a_time(i);
    _sc_r = _data->sc_ftof_1a_path(i);
  }
  else if (!std::isnan(_data->sc_ftof_2_time(i)))
  {
    // std::cout << "  status of ftof2 particle :  " << _data->status(i) << std::endl;

    _sc_t = _data->sc_ftof_2_time(i);
    _sc_r = _data->sc_ftof_2_path(i);
  }

  else if (!std::isnan(_data->sc_ctof_time(i)))
  {
    _ctof = true;
    _sc_t = _data->sc_ctof_time(i);
    _sc_r = _data->sc_ctof_path(i);
    // std::cout << "  status of ctof particle at dt_calc  1  :  " << _data->status(i) << "   _sc_t   is " << _sc_t << std::endl;
  }
  // // // // findings: no difference when this run alone or  together with _sc_t = _data->sc_ctof_time(i); part
  if (!std::isnan(_data->sc_ctof_time(i)))
  {
    _ctof_t = _data->sc_ctof_time(i);
    _ctof_r = _data->sc_ctof_path(i);
  }
}
// void Delta_T::dt_calc_ctof(int i) {
//   if (!std::isnan(_data->sc_ctof_time(i))) {
//     _ctof_t = _data->sc_ctof_time(i);
//     _ctof_r = _data->sc_ctof_path(i);
//   }
// }

float Delta_T::dt_E(int i)
{
  this->dt_calc(i);
  return _deltat(ELECTRON);
}
float Delta_T::dt_P(int i)
{
  // std::cout << "  status of ctof proton at dt_calc 2 :  " << _data->status(i) << std::endl;

  this->dt_calc(i);
  return _deltat(PROTON);
}
float Delta_T::dt_Pi(int i)
{
  this->dt_calc(i);
  return _deltat(PIP);
}
float Delta_T::dt_K(int i)
{
  this->dt_calc(i);
  return _deltat(KP);
}

float Delta_T::dt_E() { return _deltat(ELECTRON); }
float Delta_T::dt_P() { return _deltat(PROTON); }
float Delta_T::dt_Pi() { return _deltat(PIP); }
float Delta_T::dt_K() { return _deltat(KP); }
float Delta_T::dt(int pid) { return _deltat(pid); }

float Delta_T::dt_ctof_E(int i)
{
  this->dt_calc(i);
  return _ctof_deltat(ELECTRON);
}
float Delta_T::dt_ctof_P(int i)
{
  this->dt_calc(i);
  return _ctof_deltat(PROTON);
}
float Delta_T::dt_ctof_Pi(int i)
{
  this->dt_calc(i);
  return _ctof_deltat(PIP);
}
float Delta_T::dt_ctof_K(int i)
{
  this->dt_calc(i);
  return _ctof_deltat(KP);
}

float Delta_T::dt_ctof_E() { return _ctof_deltat(ELECTRON); }
float Delta_T::dt_ctof_P() { return _ctof_deltat(PROTON); }
float Delta_T::dt_ctof_Pi() { return _ctof_deltat(PIP); }
float Delta_T::dt_ctof_K() { return _ctof_deltat(KP); }
float Delta_T::dt_ctof(int pid) { return _ctof_deltat(pid); }

float Delta_T::momentum() { return _momentum; }

bool Delta_T::ctof_particle(int i)
{
  this->dt_calc(i);
  return _ctof_particle;
}
