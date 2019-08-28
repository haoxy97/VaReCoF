#ifndef MY_DYNAMICS
#define MY_DYNAMICS

#include <iostream>
#include <vector>
#include <csetjmp>

#include "global.hh"
#include "gauss.hh"

namespace Dynamics {

  using namespace std;

  extern int debug;
  extern ostream* logout;

  class Error {};
  class Form_Err: public Error {}; // format error
  class Pot_Err: public Error {};  // potential calculation error
  class Int_Err: public Error {};  // integration error
  class Range_Err: public Error {}; // indexing error
  class Run_Err: public Error {};

  class Atom {

    int num;
    int isot;
    string _name;

    static const pair<int, const char*> eldb []; // elements database
    static int         name2num (const string&);
    static const char* num2name (int);

  public:
  
    Atom (istream&);
    Atom (int, int);
    Atom (const string&, int);

    double      mass () const;
    const char* name () const { return _name.c_str(); }
    int       number () const { return num; }

  };

  /**********************************************************
   ******************** State *******************************
   **********************************************************/

  class State {

    double* _data;
    bool is_ref;

    static bool is_init;
    static vector<Atom> system;

  public:

    State(double* = 0);
    State (const State&);
    ~State ();
    State& operator= (const State&);
    State& operator= (const double*);

    static void init (const vector<Atom>& s) { system = s; is_init = true; }
    static int size () { return system.size(); }
    static double      mass (int at) { return system[at].mass(); }
    static const char* name (int at) { return system[at].name(); }

    double* pos(int at) { return _data + 3 * at; }
    double* vel(int at) { return _data + 3 * (at + system.size()); }
    const double* pos(int at) const {return _data + 3 * at;}
    const double* vel(int at) const {return _data + 3 * (at + system.size());}
    double* data () { return _data; }
    const double* data () const { return _data; }

    double kinetic_energy () const;

    // g98 stuff
    static void set_gauss_cluster ();
    void set_gauss_pos () const;
    Gauss::Method* find_method () const;

  };
    
  ostream& operator<< (ostream&, const State&);
  istream& operator>> (istream&, State&);


  /**********************************************************
   *********************** Probase **************************
   **********************************************************/

  class Probase : public State {
 
    double _time;
    
  public:

    Probase (double* =0, double =0.);
    Probase (const State&, double =0.);

    double time () const { return _time; }
    double& write_time() { return _time; }
  };

  inline Probase::Probase (double* s, double t) 
    : State(s), _time(t) {}
  inline Probase::Probase (const State& s, double t) 
    : State(s), _time(t) {}

  ostream& operator<< (ostream&, const Probase&);
  istream& operator>> (istream&, Probase&);

  /**********************************************************
   ********************* Postate ****************************
   **********************************************************/

  struct Postate : public Probase {

    double energy;

    Postate (double* =0, double =0.);
    Postate (const State&, double =0.);
  };

  inline Postate::Postate (double* s, double t) : Probase(s, t) {}
  inline Postate::Postate (const State& s, double t) : Probase(s,t) {}

  ostream& operator<< (ostream&, const Postate&);
  istream& operator>> (istream&, Postate&);

  /**********************************************************
   ******************** Propagator **************************
   **********************************************************/

  class Propagator : public Probase {

    const int lrw;
    double* rwork;
    const int liw;
    int* iwork;
    int idid;       // reports the result attributes
    int info [15];

    Gauss::Method* _method;

    vector<Postate> backup;
    vector<Probase> traj;

    const int traj_id;

    static const int max_backup_size = 10000;
    static const double time_tol; //time tolerance
    static const double ener_tol; // energy tolerance
 
    Gauss::Method* method () const { return _method; }

    const char* traj_name() const;
    const char* back_name() const;

    void push_back (const Postate&);
    void jump_back (double);
    void step_back ();
    double run     (double);

    Propagator (const Propagator&);
    Propagator& operator= (const Propagator&);

  public:

    // error codes
    enum {
      EPOT = 1, // potential calculation failed
      EJDN,     // energy jump down
      EJUP,     // energy jump up
      EINT,     // integration error
      EMTH      // method is not set
    };

    static double rel_tol; // relative error
    static double abs_tol; // absolute error
    static double time_step;    // time-integration step

    explicit Propagator(int, double* =0, double =0.);
    explicit Propagator(int, const State&, double =0.);
    explicit Propagator(int, const Probase&);

    ~Propagator ();
    Propagator& operator= (const Probase&);

    void propagate (double);

    friend ostream& operator<< (ostream&, const Propagator&);
    friend istream& operator>> (istream&, Propagator&);
  };

  ostream& operator<< (ostream&, const Propagator&);
  istream& operator>> (istream&, Propagator&);
}

#endif
