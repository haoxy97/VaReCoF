#include "dynamics.hh"
#include "slatec.hh"
#include "math.hh"

#include <cctype>
#include <cstdlib>
#include <cmath>
#include <iomanip>

namespace Dynamics {

  ostream* logout = &cout;
  int debug = 1;

  double Propagator::rel_tol = 1.e-5; // relative integration error tolerance
  double Propagator::abs_tol = 1.e-5; // absolute integration error tolerance
  double Propagator::time_step = 1.;  // time-integration step
  const double Propagator::time_tol = 0.1;  // time tolerance
  const double Propagator::ener_tol = 1.e-5; // energy tolerance

  bool State::is_init = false;
  vector<Atom> State::system;

  // atomic numbers
  enum  {
    DB_END =-1,
    HYDROGEN = 1,
    CARBON   = 6,
    NITROGEN = 7,
    OXYGEN   = 8
  };

}

/***********************************************************
 ********************** Atom **************************
 ***********************************************************/

const pair<int, const char*> Dynamics::Atom::eldb [] = {
  pair<int, const char*>(HYDROGEN, "H"), 
  pair<int, const char*>(CARBON,   "C"), 
  pair<int, const char*>(NITROGEN, "N"), 
  pair<int, const char*>(OXYGEN,   "O"),
  pair<int, const char*>(DB_END, "") // should be the last one
};

int Dynamics::Atom::name2num(const string& n)
{
  static const char funame [] = "Dynamics::Atom::name2num: ";

  for(const pair<int, const char*>* el = eldb; el->first != DB_END; ++el)
    if(el->second == n)
      return el->first;
    
  *logout << funame << "element name " << n << " is not in the database\n";
  throw Error();
}

const char* Dynamics::Atom::num2name(int n)
{
  static const char funame [] = "Dynamics::Atom::num2name: ";

  for(const pair<int, const char*>* el = eldb; el->first != DB_END; ++el)
    if(el->first == n)
      return el->second;
    
  *logout << funame << "element number " << n << " is not in the database\n";
  throw Error();
}

Dynamics::Atom::Atom(istream& from)
{
  static const char funame [] = "Dynamics::Atom::Atom: ";

  from >> _name >> isot;
  if(!from) {
    *logout << funame << "format error\n";
    throw Form_Err();
  }

  if(isdigit(_name[0])) {
    num   = atoi(_name.c_str());
    _name = num2name(num);
  }
  else
    num = name2num(_name);
}

Dynamics::Atom::Atom(int n, int i)
  : num(n), isot(i)
{
  static const char funame [] = "Dynamics::Atom::Atom: ";
  _name = num2name(num);
}

Dynamics::Atom::Atom(const string& n, int i)
  : _name(n), isot(i)
{
  static const char funame [] = "Dynamics::Atom::Atom: ";
  num = name2num(_name);
}

double Dynamics::Atom::mass () const
{
  static const char funame [] = "Dynamics::Atom::mass: ";

  switch (num) {
  case HYDROGEN:
    switch (isot) {
    case 1: return Phys_const::amu * 1.007825;
    case 2: return Phys_const::amu * 2.014;
    case 3: return Phys_const::amu * 3.01605;
    default: 
      *logout << funame << "unknown isotope " << isot 
	      << " of " << name() << "\n";
      throw Error();
    }
  case CARBON:
    switch (isot) {
    case 12: return Phys_const::amu * 12.0;
    case 13: return Phys_const::amu * 13.00335;
    default:
      *logout << funame << "unknown isotope " << isot 
	      << " of " << name() << "\n";
      throw Error();
    }
  case NITROGEN:
    switch (isot) {
    case 14: return Phys_const::amu * 14.00307;
    case 15: return Phys_const::amu * 15.00011;
    default:
      *logout << funame << "unknown isotope " << isot 
	   << " of " << name() << "\n";
      throw Error();
    }
  case OXYGEN:
    switch (isot) {
    case 16: return Phys_const::amu * 15.99491;
    case 17: return Phys_const::amu * 17.0;
    case 18: return Phys_const::amu * 18.0;
    default:
      *logout << funame << "unknown isotope " << isot 
	      << " of " << name() << "\n";
      throw Error();
    }
  default:
    *logout << funame << "unknown element number " 
	    << num << "\n";
    throw Error();
  }
}

/***********************************************************
 ************************ State *****************************
 ***********************************************************/

Dynamics::State::State (double* rd)
{
  if(!is_init) {
    *logout << "Dynamics::State::State: system is not initialized\n";
    throw Error();
  }

  if(rd) {
    _data = rd;
    is_ref = true;
  }
  else {
    _data = new double[6*system.size()];
    is_ref = false;
  }
}

Dynamics::State::State (const State& s) 
  : _data(new double[6*system.size()]), is_ref(false)
{
  int n = 6*system.size();
  for(int i = 0; i < n; ++i)
    _data[i] = s._data[i];
}

Dynamics::State::~State ()
{
  if(!is_ref)
    delete[] _data;
}

Dynamics::State::State& Dynamics::State::operator= (const State& s)
{
  int n = 6*system.size();
  for(int i = 0; i < n; ++i)
    _data[i] = s._data[i];
  return *this;
}

Dynamics::State::State& Dynamics::State::operator= (const double* s)
{
  int n = 6*system.size();
  for(int i = 0; i < n; ++i)
    _data[i] = s[i];
  return *this;
}

double Dynamics::State::kinetic_energy () const
{
  double res = 0.;
  for(int at = 0; at < system.size(); ++at) {
    const double* dp = vel(at);
    res += system[at].mass() * vdot(dp, dp, 3);
  }
  return res/2.;
}

void Dynamics::State::set_gauss_cluster ()
{
  Gauss::cluster.resize(system.size());
  for(int at = 0; at < system.size(); ++at)
    Gauss::cluster[at].num = system[at].number();
}

void Dynamics::State::set_gauss_pos ()
{
  for(int at = 0; at < size(); ++at)
    for(int i = 0; i < 3; ++i)
      Gauss::cluster[at].pos[i] = pos(at)[i];
}

Gauss::Method* Dynamics::State::find_method() const
{
  Gauss::Method* res = 0;
  set_gauss_pos();
  try {
    res = Gauss::min_pot();
  } catch(Gauss::Err) {}
  return res;
}

ostream& Dynamics::operator<< (ostream& to, const State& s)
{
  for(int at = 0; at < State::size(); ++at) {
    to << State::name(at) << " ";
    for(int i = 0; i < 3; ++i)
      to << s.pos(at)[i] << " ";
    for(int i = 0; i < 3; ++i)
      to << s.vel(at)[i] << " ";
    to << "\n";
  }
  
  return to;
}

istream& Dynamics::operator>> (istream& from, State& s)
{
  string stemp;
  for(int at = 0; at < State::size(); ++at) {
    from >> stemp;
    for(int i = 0; i < 3; ++i)
      from >> s.pos(at)[i];
    for(int i = 0; i < 3; ++i)
      from >> s.vel(at)[i];
  }
  
  return from;
}

/***********************************************************
 ********************** Probase *************************
 ***********************************************************/

ostream& Dynamics::operator<< (ostream& to, const Probase& p)
{
  to << p.time() << "\n";
  to << static_cast<const Dynamics::State&>(p);
  return to;
}

istream& Dynamics::operator>> (istream& from, Probase& p)
{
  from >> p.write_time();
  from >> static_cast<Dynamics::State&>(p);
  return from;
}

/**********************************************************
 ********************* Postate ****************************
 **********************************************************/

ostream& Dynamics::operator<< (ostream& to, const Postate& p)
{
  to << p.energy << " ";
  to << static_cast<const Dynamics::Probase&>(p);
  return to;
}

istream& Dynamics::operator>> (istream& from, Postate& p)
{
  from >> p.energy;
  from >> static_cast<Dynamics::Probase&>(p);
  return from;
}

/***********************************************************
 ********************** Propagator *************************
 ***********************************************************/

Dynamics::Propagator::Propagator(int id, double* s, double t) 
  : Probase(s, t), lrw(130 + 126 * size()),rwork(new double[lrw]), 
  liw(51), iwork(new int[liw]), idid(0), traj_id(id), _method(0)
{
  // ddeabm subroutine parameters setting
  info [0] = 0;     // start new trajectory
  info [1] = 0;     // rel_tol and abs_tol are scalars
  info [2] = 0;     // the solution only at the end point
  info [3] = 0;     // the integration can go beyond final point

  traj.reserve(max_backup_size);
  backup.reserve(max_backup_size);
}

Dynamics::Propagator::Propagator(int id, const State& s, double t) 
  : Probase(s, t), lrw(130 + 126 * size()),rwork(new double[lrw]), 
  liw(51), iwork(new int[liw]), idid(0), traj_id(id), _method(0)
{
  // ddeabm subroutine parameters setting
  info [0] = 0;     // start new trajectory
  info [1] = 0;     // rel_tol and abs_tol are scalars
  info [2] = 0;     // the solution only at the end point
  info [3] = 0;     // the integration can go beyond final point

  traj.reserve(max_backup_size);
  backup.reserve(max_backup_size);
}

Dynamics::Propagator::Propagator(int id, const Probase& pb) 
  : Probase(pb), lrw(130 + 126 * size()),rwork(new double[lrw]), 
  liw(51), iwork(new int[liw]), idid(0), traj_id(id), _method(0)
{
  // ddeabm subroutine parameters setting
  info [0] = 0;     // start new trajectory
  info [1] = 0;     // rel_tol and abs_tol are scalars
  info [2] = 0;     // the solution only at the end point
  info [3] = 0;     // the integration can go beyond final point

  traj.reserve(max_backup_size);
  backup.reserve(max_backup_size);
}

/*
Dynamics::Propagator::Propagator(const Propagator& in) 
  : Probase(in), lrw(in.lrw), rwork(new double[lrw]), 
  liw(in.liw), iwork(new int[liw]), idid(in.idid), _method(in._method)
{
  for(int i = 0; i < 15; ++i)
    info[i] = in.info[i];
  for(int i = 0; i < lrw; ++i)
    rwork[i] = in.rwork[i];
  for(int i = 0; i < liw; ++i)
    iwork[i] = in.iwork[i];

  traj = in.traj;
  backup = in.backup;
  traj.reserve(max_backup_size);
  backup.reserve(max_backup_size);
}

Dynamics::Propagator& Dynamics::Propagator::operator=(const Propogator& in) 
{
  Probase::operator=(in); 
  idid = in.idid;
  for(int i = 0; i < 15; ++i)
    info[i] = in.info[i];
  for(int i = 0; i < lrw; ++i)
    rwork[i] = in.rwork[i];
  for(int i = 0; i < liw; ++i)
    iwork[i] = in.iwork[i];

  _method = in._method;

  traj = in.traj;
  backup = in.backup;
  traj.reserve(max_backup_size);
  backup.reserve(max_backup_size);

  return *this;
}
*/

Dynamics::Propagator::~Propagator()
{
  delete[] rwork;
  delete[] iwork;

  // save trajectory & backup data
  ofstream to(traj_name(), ios::app);
  for(int i = 0; i < traj.size(); ++i)
    to << traj[i] << "\n";
  to.close();
  to.open(back_name(), ios::app);
  for(int i = 0; i < backup.size(); ++i)
    to << backup[i] << "\n";
}

Dynamics::Propagator& Dynamics::Propagator::operator= (const Probase& in) 
{
  Probase::operator=(in); 

  info [0] = 0;     // start new trajectory
  info [1] = 0;     // rel_tol and abs_tol are scalars
  info [2] = 0;     // the solution only at the end point
  info [3] = 0;     // the integration can go beyond final point

  return *this;
}

const char* Dynamics::Propagator::traj_name () const
{
  static char name [99];
  sprintf(name, "traj_%i.dat", traj_id);
  return name;
}

const char* Dynamics::Propagator::back_name () const
{
  static char name [99];
  sprintf(name, "traj_%i.back", traj_id);
  return name;
}

namespace Dynamics {

  extern "C" void set_dvd (const double&, const double*, 
			   double*, void*, void*);

  // parameters transfered between set_dvd and calling function
  struct dvd_par {
    Postate post;
    jmp_buf jmp;
    Gauss::Method* method;

    dvd_par() : method(0) {}
  };

}

extern "C" void Dynamics::set_dvd (const double& time, const double* x, 
				   double* dx, void* rpar, void* ipar)
{
  static const char funame [] = " Dynamics::set_dvd: ";

  dvd_par& par = *static_cast<dvd_par*>(rpar);

  // find the energy and forces
  Postate dv(const_cast<double*>(x), time);
  dv.set_gauss_pos();
  par.method->init();
  if(!par.method->apply(Gauss::DIRECT | Gauss::FORCE | Gauss::READ)) {
    *logout << funame << "potential calculation failed, long jump initiated\n";
    longjmp(par.jmp, Propagator::EPOT);
  }
  dv.energy = par.method->rel_energy();

  // set return value for dynamic variables, time, and energy
  par.post = dv;

  // find dynamic variables derivatives
  State dvd(dx);
  try {
    par.method->read_forces();
  } catch(Gauss::Err) {
    *logout << funame << "reading forces failed, long jump initiated\n";
    longjmp(par.jmp, Propagator::EPOT);
  }

  for(int at = 0; at < State::size(); ++at)
    for(int i = 0; i < 3; ++i) {
      dvd.pos(at)[i] = dv.vel(at)[i];
      dvd.vel(at)[i] = Gauss::cluster[at].force[i]/State::mass(at);
    }
}

void Dynamics::Propagator::push_back(const Postate& post)
{
  // clean up
  while(backup.size() && backup.rbegin()->time() > post.time() - time_tol)
    backup.pop_back();

  // resize backup
  if(backup.size() >= max_backup_size) {
    int itemp = backup.size() - max_backup_size + max_backup_size / 10;
    ofstream to(back_name(), ios::app);
    for(int i = 0; i < itemp; ++i)
      to << backup[i] << "\n";
    backup.erase(backup.begin(), backup.begin() + itemp);
  }
  // push back
  backup.push_back(post);
}

double Dynamics::Propagator::run (double fin_time)
{
  static const char* funame = "Dynamics::Propagator::run: ";

  double dtemp;
  int itemp;

  const int step_num = (int)ceil((fin_time - time())/time_step);

  // set method
  if(!method()) {
    *logout << funame << "method is not initialized\n";
    throw Error();
  }
  
  // resize trajectory container
  itemp = traj.size() + step_num - max_backup_size;
  if(itemp > 0) {
    // output
    // ...
    ofstream to(traj_name(), ios::app);
    for(int i = 0; i < itemp; ++i)
      to << traj[i] << "\n";
    traj.erase(traj.begin(), traj.begin() + itemp);
  }

  dvd_par par; // parameters transfered between set_dvd and calling function
  par.method = method();

  Gauss::Method* new_method = 0;
  // force calculation failure
  switch(setjmp(par.jmp)) {
  case EPOT:
    *this = *traj.rbegin();
    traj.pop_back();
    *logout << funame << "long jump acknowledged, time = "
	    << time() << "\n";
    new_method = backup.rbegin()->find_method();
    if(!new_method || new_method->rel_energy() > 
       backup.rbegin()->energy - ener_tol) {
      *logout << funame << "cannot find a new method\n";
      throw Pot_Err();
    }
    _method = new_method;
    step_back();
    return time();

  default:
    *logout << funame << "unknown long jump code\n";
    throw Error();
  }

  const int dv_size = size() * 6;
  for(int step_ind = 0; step_ind < step_num; ++step_ind) {//main cycle
    traj.push_back(*this);
    ddeabm_(set_dvd, dv_size, write_time(), data(), time() + time_step, info, 
	    rel_tol, abs_tol, idid, rwork, lrw, iwork, liw, 
	    &par, 0);
    info[0] = 1;

    if(idid < -1) {// integrator failure
      *this = *traj.rbegin(); 
      traj.pop_back();
      *logout << funame << "Integration error, idid = " << idid 
	      << " time = " << time() << "\n";
      new_method = backup.rbegin()->find_method();
      if(!new_method || new_method->rel_energy() > 
	 backup.rbegin()->energy - ener_tol) {
	*logout << funame << "cannot find a new method\n";
	throw Int_Err();
      }
      _method = new_method;
      step_back();
      return time();
    }

    push_back(par.post);

    if(backup.size() == 1)
      continue;
    
    double ener_incr = backup.rbegin()->energy 
	+ backup.rbegin()->kinetic_energy()
	- (backup.rbegin() + 1)->energy 
	- (backup.rbegin() + 1)->kinetic_energy();

    if(ener_incr < - ener_tol) {// energy DOWN
      *logout << funame << "energy jump DOWN " 
	      << ener_incr / Phys_const::kcal 
	      << " kcal/mol, time = " << time() << "\n";
      
      backup.rbegin()->energy -= ener_incr;
      step_back();
      return time();
    }

    if(ener_incr > ener_tol) {// energy DOWN
      *logout << funame << "energy jump UP " 
	      << ener_incr / Phys_const::kcal 
	      << " kcal/mol, time = " << time() << "\n";
    }
  }// main cycle

  method()->backup();
  new_method = backup.rbegin()->find_method();
  if(!new_method || new_method->rel_energy() > 
       backup.rbegin()->energy - ener_tol) {
    method()->restore();
    return time();
  }
  _method = new_method;
  step_back();

  return time();
}

void Dynamics::Propagator::step_back ()
{
  static const char funame [] = " Dynamics::Propagator::step_back: ";

  double dtemp;
  int itemp;

  double curr_ener = method()->rel_energy() - backup.rbegin()->energy;
  if(curr_ener > -ener_tol)
    return;
  
  bool is_fail = false;
  double prev_ener;
  vector<Postate>::reverse_iterator curr;
  for(curr = backup.rbegin() + 1; curr != backup.rend(); ++curr) {
    method()->backup();
    prev_ener = curr_ener;
    curr->set_gauss_pos();
    method()->init();
    if(!method()->apply(Gauss::DIRECT | Gauss::FORCE | Gauss::READ)) {
      is_fail = true;
      break;
    }
    curr_ener = method()->rel_energy() - curr->energy;
    if(curr_ener > -ener_tol)
      break;
  }

  if(curr == backup.rend()) {// backup exhausted
    jump_back(backup.begin()->time());
    *logout << funame << "backup scanned up to the beginning, "
      "no intersection found, time = " << time() << "\n";
    return;
  }
	 
  if(curr == backup.rbegin() + 1) {// no intersection found
    *logout << funame << "energy jump disregarded\n";
    method()->restore();
    return;
  }
	  
  if(is_fail) {
    method()->restore();
    jump_back((curr-1)->time());
    *logout << funame << "method " << method()->name() 
	    << " failed, restart time = " << time() << "\n";
    return;
  }

  method()->restore();
  dtemp = (curr_ener * (curr - 1)->time() - prev_ener * curr->time())
    / (curr_ener - prev_ener);
  jump_back(dtemp);

}

void Dynamics::Propagator::jump_back (double new_time)
{
  static const char funame [] = " Dynamics::Propagator::jump_back: ";

  double dtemp;
  int itemp;

  int step_ind = traj.size() - 1 + 
    (int)floor((new_time - traj.rbegin()->time())/time_step);
  if(step_ind < 0) {
    *logout << funame << "new time is too far in the past\n";
    throw Error();
  }
  else if(step_ind >= traj.size()) {
    *logout << funame << "new time is in the future\n";
    throw Error();
  }

  *this = traj[step_ind];
  traj.resize(step_ind);

  set_gauss_pos();
  if(!method()->apply(Gauss::DIRECT | Gauss::FORCE | Gauss::READ)) {
    *logout << funame << "new method failed\n";
    throw Error();
  }

  push_back(*this);
  backup.rbegin()->energy = method()->rel_energy();

  *logout << funame << "new time = " << time() << "\n";
}

void Dynamics::Propagator::propagate (double time_range)
{
  static const char* funame = "Dynamics::Propagator::propagate: ";

  //static const int max_fail_count = 5;
  //int fail_count = 0;

  double dtemp;
  int itemp;

  // set method
  if(!method() && !(_method = find_method())) {
    *logout << funame << "method initialization failed\n";
    throw Pot_Err();
  }

  double start_time = time();
  const double finish_time = time() + time_range;
  while(run(finish_time) < finish_time - time_tol){}
}
