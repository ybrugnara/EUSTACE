#ifndef SEGCLUST_SEG_MEAN_H
#define SEGCLUST_SEG_MEAN_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <iostream>
//#include <fstream>
#define Min(i, j)  ((i) < (j) ? (i) : (j))
#define Max(i, j)  ((i) > (j) ? (i) : (j))
namespace segclust
{

void compute_segmentation_mean(gsl_vector *x, int Kmax, int lmin, int lmax, bool vh, gsl_vector **J_est, gsl_matrix **t_est);

class Segmentation_mean
{
  public:
  int      _lengthx;     // length of the data
  int      _lmin;        // minimum size for a segment
  int      _lmax;        // maximum size for a segment
  bool     _vh;          // variance homogeneity 
  double  *_x;           // data
  int      _Kmax;        // maximum number of segments
  double **_Cost;        // Cost matrix
  double **_D;           // matrix of the lengths of the minimum distances
  int    **_Breaks;      // possible breakpoints
  int    **_BestBreaks;  // best breakpoints
  Segmentation_mean(int n, int k, int minlength, int maxlength, bool varh);
  void Init(double *Data, bool varh);
  void DynProg(int k);
  void Terminate();
  friend std::ostream & operator << (std::ostream &s, const Segmentation_mean &Seg_mean);
  ~Segmentation_mean();
};
}
#endif
