#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "math.hh"
#include "tmatrix.hh"
#include "error.hh"

int main(int argc, char* argv []) {

  // temporary variables
  int itemp;
  double dtemp;

  if (argc != 2)
    error("usage: interpolate data_file");

  ifstream from;
  from.open(argv[1]);
  if(!from) {
    cout << "cannot open " << argv[1] << " file\n";
    exit(1);
  }

  int n;
  from >> n;
  Array<double> x(n), y(n);
  for(int i = 0; i < n; ++i) {
    from >> x[i] >> y[i];
    if(i && x[i] <= x[i-1]) {
      std::cout << "x should be in ascending order\n";
      exit(1);
    }
  }
  if(!from) {
    std::cout << "cannot read interpolated data\n";
    exit(1);
  }
  

  std::vector<double> x0;
  while(from >> dtemp)
    x0.push_back(dtemp);

  Spline sp(x.data(), y.data(), x.size());
  for(int i = 0; i < x0.size(); ++i) {
    if(x[0] < x0[i] && x0[i] < x[x.size()-1])
      std::cout << sp.fit(x0[i], 0) << "\n";
    else
      std::cout << "x0 = " << x0[i] << " is out of range\n";
  }

  return 0;
}
