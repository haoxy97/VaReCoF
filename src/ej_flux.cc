#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <sstream>

#include "input.hh"
#include "units.hh"
#include "rotd.hh"
#include "math.hh"
#include "tmatrix.hh"
#include "slatec.hh"
#include "error.hh"
#include "divsur.hh"
#include "interflux.hh"

int main (int argc, char* argv[])
{
  const char funame [] = "ej_flux: ";

  if (argc != 2)
    error("usage: ej_flux inp_file");

  // temporary variables
  int itemp;
  double dtemp;
  std::string token, line, comment, key;

  // input 
  const InputKey multi_key  = "MultiInputFile";// input file for multi
  std::string multi_input_file;

  std::map<std::pair<int, int>, double> ener_shift;
  const InputKey corr_key = "EnergyCorrection(kcal/mol)";// energy correction
  //std::map<int, double> ecorr;

  const InputKey ener_key = "EnergyGrid(1/cm)";// new energy grid
  Grid new_ener_grid;

  const InputKey amom_key = "AngularMomentumGrid(a.u.)";// new angular momentum grid
  Grid  new_amom_grid;

  const InputKey  out_key = "OutputFile";// ej_flux output file
  std::string out_file;

  const InputKey face_key = "Face"; 

  int pes_val = 0;
  const InputKey pes_key = "ElectronicSurface"; // electronic state weights

  std::ifstream from(argv[1]);
  if(!from) {
    std::cout << funame << "cannot open " << argv[1] << " file\n";
    exit(1);
  }

  double min_val, step;
  int ds, ds_size, nstep, size;
  int ds_face = -1;
  while(from >> token) {
    if(multi_key == token) {
      from >> multi_input_file;
    }
    else if(key == corr_key) {
      from >> size;
      if(!from) {
	std::cout << funame << token << "failed to read number of entries\n";
	exit(1);
      }
      std::getline(from, comment);
      for(int i = 0; i < size; ++i) {
	std::getline(from, line);
	std::istringstream ist(line);
	int sur, face;
	ist >> sur >> face >> dtemp;
	if(!ist) {
	  std::cout << funame << token
		    << ": entry format: sur# face# shift\n";
	  exit(1);
	}
	ener_shift[std::make_pair(sur, face)] = dtemp * Phys_const::kcal;
      }
      if(!from) {
	std::cout << funame << token << ": wrong number of entries\n";
	exit(1);
      }
    }
    else if(ener_key == token) {// energy grid
      from >> min_val >> nstep >> step;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	exit(1);
      }
      new_ener_grid.set(Phys_const::incm*min_val, nstep, Phys_const::incm*step);
    }
    else if(amom_key == token) {//angular momentum grid
      from >> min_val >> nstep >> step;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	exit(1);
      }
      new_amom_grid.set(min_val, nstep, step);
    }
    else if(pes_key == token) {// surface
      from >> pes_val;
      if(!from) {
	std::cout << funame << token << " reading failed\n";
	exit(1);
      }
    }
    else if(out_key == token) {// output file name
      from >> out_file;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	exit(1);
      }
    }
    else if (face_key == token) {// read face
      from >> ds_face;
      if(!from) {
	std::cout << funame << "cannot read " << token << "\n";
	exit(1);
      }
    }
    else {
      std::cout << funame << "unknown keyword " << token << ". Known keywords are:\n";
      InputKey::show_all();
      exit(1);
    }
  }
  from.close();
  from.clear();

  // reading multi input
  std::map <string, Read> input_data;
  typedef map<string, Read>::iterator input_data_it;
 
  std::vector<double> temperature;
  input_data ["tmpr_grid"]   = Read(temperature);
  std::vector<double> ener_grid;
  input_data ["ener_grid"]   = Read(ener_grid);
  std::vector<double> amom_grid;
  input_data ["amom_grid"]  = Read(amom_grid);

  std::string flux_file; // flux output file
  input_data ["flux_file"]  = Read(flux_file, "flux.dat");      
  std::string mol_spec_file; // molecular specifications file name
  input_data ["mol_spec_file"] = Read(mol_spec_file, "structure.inp");

  int pes_size;
  input_data ["pes_size"] = Read(pes_size, 1);

  from.open(multi_input_file.c_str());
  if(!from) {
    std::cout << funame << "cannot open " << argv[1] << " file\n";
    exit(1);
  }

  while(from >> key) {
    input_data_it i = input_data.find(key);
    if (i == input_data.end())
      getline(from, comment);
    else
      from >> i->second;
  }
  from.close();
  from.clear();

  // check if all parameters were initialized
  bool is_init = true;
  for(input_data_it it = input_data.begin(); it != input_data.end(); ++it)
    if(!it->second.is_init()) {
      std::cout << funame << it->first << " is not initialized\n";
      is_init = false;
    }
  if(!is_init)
    exit(1);

  // checking for default values
  std::cout << "Default parameters:\n";
  for(input_data_it it = input_data.begin(); it != input_data.end(); ++it)
    if(it->second.is_default())
      std::cout << it->first << "  " << it->second << "\n";
  std::cout << "\n";

  // read molecular structure from the file  
  mol_init(mol_spec_file);

  Dynvar temp_dv;
  Array_2<double> ej_flux(new_ener_grid.size(), amom_grid.size());
  Array_2<double> total_ej_flux(new_ener_grid.size(), amom_grid.size());
  Matrix<double> min_flux(new_ener_grid.size(), amom_grid.size());
  Matrix<double> total_min_flux(new_ener_grid.size(), amom_grid.size());
  double vmin;

  for(int en = 0; en < new_ener_grid.size(); ++en)
    for(int am = 0; am < amom_grid.size(); ++am)
      total_min_flux[en][am] = min_flux[en][am] = -1.;

  std::map<std::pair<int, int>, double>::const_iterator esp;
  // shifted grid and flux
  Array<double> shift_grid(ener_grid.size());
  Array<double> shift_flux(ener_grid.size());

  // read from flux file
  from.open(flux_file.c_str());
  if(!from) {
    std::cout << funame << "cannot open " << flux_file.c_str() << " file\n";
    exit(1);
  }

  while (from >> ds >> ds_size) {//input cycle
    std::cout << "ds = " << ds << std::endl;
    total_ej_flux.init();
    ej_flux.init();

    // read ej-resolved flux
    for(int face = 0; face < ds_size; ++face) {//face cycle
      
      from >> itemp;// sampling numbers
      std::getline(from, token);

      for(int pes = 0; pes < pes_size; ++pes)
	if(pes != pes_val)
	  from >> dtemp;
	else
	  from >> vmin;
      
      for(int pes = 0; pes < pes_size; ++pes)
	from >> temp_dv;// minimal configurations

      for(int pes = 0; pes < pes_size; ++pes)
	for(int i = 0; i < temperature.size(); ++i)// thermal flux
	  from >> dtemp >> dtemp;

      for(int pes = 0; pes < pes_size; ++pes)
	for(int i = 0; i < ener_grid.size(); ++i)// microcannonical flux
	  from >> dtemp >> dtemp;

      for(int en = 0; en < ener_grid.size(); ++en)
	shift_grid[en] = ener_grid[en];

      esp = ener_shift.find(std::make_pair(ds, face));
      if(esp != ener_shift.end()) {
	//std::cout << "energy shift = " << esp->second / Phys_const::kcal << "\n";
	vmin += esp->second;
	for(int en = 0; en < ener_grid.size(); ++en)
	  shift_grid[en] +=  esp->second;
      }

      for(int pes = 0; pes < pes_size; ++pes)// ej_resolved flux
	  for(int am = 0; am < amom_grid.size(); ++am)
	    if(pes != pes_val)
	      for(int en = 0; en < ener_grid.size(); ++en)
		from >> dtemp >> dtemp;
	    else {
	      for(int en = 0; en < ener_grid.size(); ++en)
		from >> shift_flux[en] >> dtemp;
	      //std::cout << "am = " << am << std::endl;
	      InterFlux inter_flux(shift_grid, shift_flux, vmin);
	      for(int en = 0; en < new_ener_grid.size(); ++en) {
		//std::cout << "en = " << en << std::endl;
		dtemp = inter_flux(new_ener_grid[en]);
		total_ej_flux(en, am) += dtemp;
		if(face == ds_face)
		  ej_flux(en, am) = dtemp;
	      }
	    }
    }// face cycle
    std::cout << "e-interpolation done" << std::endl;

    if(!from) {
      std::cout << funame << "flux file " << flux_file << " is corrupted\n";
      exit(1);
    }

    // minimal flux
    for(int am = 0; am < amom_grid.size(); ++am)
      for(int en = 0; en < new_ener_grid.size(); ++en) {
	dtemp = total_ej_flux(en, am);
	if(total_min_flux[en][am] < 0. || dtemp < total_min_flux[en][am]) {
	  total_min_flux[en][am] = dtemp;
	  if(ds_face < 0)
	    min_flux[en][am] = dtemp;
	  else
	    min_flux[en][am] = ej_flux(en, am);
	}
      }
  }//input cycle
  from.close();
  from.clear();

  // check that all flux values vere initialized
  for(int am = 0; am < amom_grid.size(); ++am)
    for(int en = 0; en < new_ener_grid.size(); ++en)
      if(min_flux[en][am] < 0.) {
	std::cout << funame << "flux with energy = "
		  << new_ener_grid[en]/Phys_const::incm	 << "cm^-1 not initialized\n";
	exit(1);
      }

  // normalization
  double nfac = sqrt(2.*M_PI);
  dtemp =sqrt(mol_array[0]->mass() * mol_array[1]->mass()
	      /(mol_array[0]->mass() + mol_array[1]->mass()));
  nfac *= dtemp * dtemp * dtemp;
  nfac *= mol_array[0]->stat_sum() * mol_array[1]->stat_sum();

  std::cout << "normalization factor = " << nfac << std::endl;

  for(int am = 0; am < amom_grid.size(); ++am) {
    dtemp = amom_grid[am] * nfac;
    for(int en = 0; en < new_ener_grid.size(); ++en)
      min_flux[en][am] *= dtemp;
  }

  
  // angular momentum interpolation
  Array<double> mod_amom_grid(amom_grid.size());
  for(int am = 0; am < amom_grid.size(); ++am)
    mod_amom_grid[am] = amom_grid[am];

  Matrix<double> flux_res(new_ener_grid.size(), new_amom_grid.size());
  for(int en = 0; en < new_ener_grid.size(); ++en) {
    Spline sp(mod_amom_grid.data(), min_flux[en], mod_amom_grid.size());
    for(int am = 0; am < new_amom_grid.size(); ++am) {
      flux_res[en][am] = sp.fit(new_amom_grid[am], 0);
      if(flux_res[en][am] < 0.1)
	flux_res[en][am] = 0.;
    }
  }

  // make sure that the flux as a function of e is monothonic
  for(int am = 0; am < new_amom_grid.size(); ++am) {
    dtemp = flux_res[new_ener_grid.size() - 1][am];
    for(int en = new_ener_grid.size() - 1; en > 0; --en) {
      double& fl = flux_res[en - 1][am];
      if(fl > dtemp)
	fl = dtemp;
      dtemp = fl;
    }
  }

  int dof_num = 3;
  for (int i = 0; i < 2; ++i)
    switch(mol_array[i]->type()) {
    case ATOM:
      break;
    case LINEAR:
      dof_num += 2;
      break;
    case NONLINEAR:
      dof_num += 3;
    }

  std::cout << "Rate, 10^11 cm^3/sec:\n";
  nfac  /= 612.6 * 4. * M_PI * new_amom_grid.step() * new_ener_grid.step();
  double tfac;
  for(int t = 0; t < temperature.size(); ++t) {
    tfac = nfac * pow(temperature[t],  double(dof_num) / 2.);
    dtemp = 0.;
    for(int am = 0; am < new_amom_grid.size(); ++am)
      for(int en = 0; en < new_ener_grid.size(); ++en)
	dtemp += new_amom_grid[am] * exp(-new_ener_grid[en]/temperature[t]) 
	  * flux_res[en][am];
    dtemp /=  tfac;
    std::cout << temperature[t] / Phys_const::kelv << " "
	      << dtemp << "\n";
  }

  // output
  ofstream ej_out(out_file.c_str());
  for(int am = 0; am < new_amom_grid.size(); ++am)
    for(int en = 0; en < new_ener_grid.size(); ++en)
      ej_out << new_ener_grid[en] / Phys_const::incm << "  " << new_amom_grid[am] << "  "
	     << flux_res[en][am] << "\n";
}
