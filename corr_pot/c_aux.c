#include <stdio.h>
#include <math.h>

/************** Vector product function needed by sjk_pot ************/

void cross_ (const double* v1, const double* v2, double* res)
{
  res[0] = v1[1] * v2[2] - v1[2] * v2[1];
  res[1] = v1[2] * v2[0] - v1[0] * v2[2];
  res[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

/******************* distance between two vectors ********************/

double vdistance_ (const double* vec1, const double* vec2, int dim)
{
  double dist2 = 0.0, temp;
  int i;

  for (i = 0; i < dim; ++i) {
    temp = vec1[i] - vec2[i];
    dist2 += temp * temp;
  }
  return sqrt(dist2);
}

/************ Tak-San Ho's OH + H potential ***************/

void kernel3d_(const char*);
void pes_h2x__(const int*, const double*, const double*, 
	       const double*, double*, double*, double*);

double oh_h (const double* r, const double* rpar, const int* ipar)
{
  const char h2o_dir [] = "/tcghome/sjk/rotd/src/sjk_pot/h2o/COH21Ap";

  const int nfrag = 2;
  const int ndim = 3;
  const int natom = 100;

  const int xtype = 2;
  const double infty = 30.;

  static int is_first = 1;
  static double ref_en;
  static double oh1;

  double oh2, hh;
  double v2, v3;
  double rr[3][3];
  double ener;
  int i, ifrag, iatom, idim;

  for(i = 0; i < 3; ++i) {
    ifrag = i < 2 ? 0 : 1;
    iatom = i < 2 ? i : 0;
    for(idim = 0; idim < 3; ++idim)
      rr[i][idim] = *(ifrag + nfrag * iatom + nfrag * natom * idim + r);
  }

  if(is_first) {
    is_first = 0;
    oh1 = vdistance_(rr[0], rr[1], 3);
    kernel3d_(h2o_dir);
    pes_h2x__(&xtype, &oh1, &infty, &infty, &ref_en, &v2, &v3);
  }

  oh2 = vdistance_(rr[0], rr[2], 3);
  hh  = vdistance_(rr[1], rr[2], 3);

  pes_h2x__(&xtype, &oh1, &oh2, &hh, &ener, &v2, &v3);
  return ener - ref_en;
}

/*********** Kendrick O + O + H potential *************************/

void dimpot_(const double*, double*);

double o2_h_kendrick (const double* r, const double* rpar, const int* ipar)
{
  const int nfrag = 2;
  const int ndim = 3;
  const int natom = 100;

  const double infinity = 30.;

  static double ref_en;
  static double dist [3];
  static int is_first = 1;

  int i, ifrag, iatom, idim;
  double rr[3][3];
  double ener[2];

  for(i=0; i < 3; ++i) {
    ifrag = i < 2 ? 0 : 1;
    iatom = i < 2 ? i : 0;
    for(idim = 0; idim < 3; ++idim)
      rr[i][idim] = *(ifrag + nfrag * iatom + nfrag * natom * idim + r);
  }

  if(is_first) {
    is_first = 0;
    dist[0] = vdistance_(rr[0], rr[1], 3);
    dist[1] = infinity;
    dist[2] = infinity;
    dimpot_(dist, ener);
    ref_en = ener[0];
  }

  dist[1] = vdistance_(rr[0], rr[2], 3);
  dist[2] = vdistance_(rr[1], rr[2], 3);

  dimpot_(dist, ener);
  return ener[0] - ref_en;
}

double oh_o_kendrick (const double* r, const double* rpar, const int* ipar)
{
  const int nfrag = 2;
  const int ndim = 3;
  const int natom = 100;

  const double infinity = 30.;

  static double ref_en;
  static double dist [3];
  static int is_first = 1;

  int i, ifrag, iatom, idim;
  double rr[3][3];
  double ener[2];

  for(i=0; i < 3; ++i) {
    ifrag = i < 2 ? 0 : 1;
    iatom = i < 2 ? i : 0;
    for(idim = 0; idim < 3; ++idim)
      rr[i][idim] = *(ifrag + nfrag * iatom + nfrag * natom * idim + r);
  }

  if(is_first) {
    is_first = 0;
    dist[1] = vdistance_(rr[0], rr[1], 3);
    dist[0] = infinity;
    dist[2] = infinity;
    dimpot_(dist, ener);
    ref_en = ener[0];
  }

  dist[0] = vdistance_(rr[0], rr[2], 3);
  dist[2] = vdistance_(rr[1], rr[2], 3);

  dimpot_(dist, ener);
  return ener[0] - ref_en;
}

/***************** Varandas O + O + H potential *********************/

void ho2sur_(const double*, double*);

double o2_h_varandas (const double* r, const double* rpar, const int* ipar)
{
  const int nfrag = 2;
  const int ndim = 3;
  const int natom = 100;

  const double infinity = 30.;

  static double ref_en;
  static double dist [3];
  static int is_first = 1;

  int i, ifrag, iatom, idim;
  double rr[3][3];
  double ener;

  for(i=0; i < 3; ++i) {
    ifrag = i < 2 ? 0 : 1;
    iatom = i < 2 ? i : 0;
    for(idim = 0; idim < 3; ++idim)
      rr[i][idim] = *(ifrag + nfrag * iatom + nfrag * natom * idim + r);
  }

  if(is_first) {
    is_first = 0;
    dist[0] = vdistance_(rr[0], rr[1], 3);
    dist[1] = infinity;
    dist[2] = infinity;
    ho2sur_(dist, &ref_en);
  }

  dist[1] = vdistance_(rr[0], rr[2], 3);
  dist[2] = vdistance_(rr[1], rr[2], 3);

  ho2sur_(dist, &ener);
  return ener - ref_en;
}

double oh_o_varandas (const double* r, const double* rpar, const int* ipar)
{
  const int nfrag = 2;
  const int ndim = 3;
  const int natom = 100;

  const double infinity = 30.;

  static double ref_en;
  static double dist [3];
  static int is_first = 1;

  int i, ifrag, iatom, idim;
  double rr[3][3];
  double ener;

  for(i=0; i < 3; ++i) {
    ifrag = i < 2 ? 0 : 1;
    iatom = i < 2 ? i : 0;
    for(idim = 0; idim < 3; ++idim)
      rr[i][idim] = *(ifrag + nfrag * iatom + nfrag * natom * idim + r);
  }

  if(is_first) {
    is_first = 0;
    dist[1] = vdistance_(rr[0], rr[1], 3);
    dist[0] = infinity;
    dist[2] = infinity;
    ho2sur_(dist, &ref_en);
  }

  dist[0] = vdistance_(rr[0], rr[2], 3);
  dist[2] = vdistance_(rr[1], rr[2], 3);

  ho2sur_(dist, &ener);
  return ener - ref_en;
}

/************************* CH + H potential *************************/

void surfin_(const int*);
void surf_(const double*, double*, double*, 
	   double*, double*, int*, const int*);

double ch_h (const double* r, const double* rpar, const int* ipar)
{
  const int zero = 0;

  const int nfrag = 2;
  const int natom = 100;
  const int ndim = 3;

  double res, dea, deb, hab;
  int flag;
  int i, idim, iatom, ifrag;
  double rr[3][3];

  static int is_first = 1;
  static double ref_en; 

  for(i=0; i < 3; ++i) {
    ifrag = i < 2 ? 0 : 1;
    iatom = i < 2 ? i : 0;
    for(idim = 0; idim < 3; ++idim)
      rr[i][idim] = *(ifrag + nfrag * iatom + nfrag * natom * idim + r);
  }

  if(is_first) {
    is_first = 0;
    /*    printf("ch_h: initializing potential\n");*/
    surfin_(&zero);

    rr[2][0] += 100.;
    surf_(rr[0], &ref_en, &dea, &deb, &hab, &flag, &zero);
    rr[2][0] -= 100.;
  }
  
  surf_(rr[0], &res, &dea, &deb, &hab, &flag, &zero);
  res -= ref_en;
  return res;
}
