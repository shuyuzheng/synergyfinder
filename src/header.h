#ifdef _OPENMP
#include <omp.h>
#endif
#define USE_FC_LEN_T
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#ifndef FCONE
#define FCONE
#endif
//#include <R_ext/Applic.h>
#include <R_ext/Utils.h>
#include <time.h>

#define MINF -1.0e15
#define EPS DBL_EPSILON

// ///////////////////////////////////
// //  From schlather.c
// //
// void schlatherfull(int *covmod, double *data, double *dist, int *nSite, int *nObs,
//                    int *dim, int *weighted, double *weights, double *locs,
//                    double *scales, double *shapes, double *nugget, double *range,
//                    double *smooth, double *smooth2, int *fitmarge, double *dns);
// void schlatherdsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs, int *dim,
//                       int *weighted, double *weights, double *locdsgnmat, double *locpenmat,
//                       int *nloccoeff, int *npparloc, double *locpenalty, double *scaledsgnmat,
//                       double *scalepenmat, int *nscalecoeff, int *npparscale,
//                       double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
//                       int *nshapecoeff, int *npparshape, double *shapepenalty, int *usetempcov,
//                       double *tempdsgnmatloc, double *temppenmatloc, int *ntempcoeffloc,
//                       int *nppartempcoeffloc, double *temppenaltyloc, double *tempdsgnmatscale,
//                       double *temppenmatscale, int *ntempcoeffscale, int *nppartempcoeffscale,
//                       double *temppenaltyscale, double *tempdsgnmatshape, double *temppenmatshape,
//                       int *ntempcoeffshape, int *nppartempcoeffshape, double *temppenaltyshape,
//                       double *loccoeff, double *scalecoeff, double *shapecoeff,
//                       double *tempcoeffloc, double *tempcoeffscale, double *tempcoeffshape,
//                       double *nugget, double *range, double *smooth, double *smooth2, double *dns);
// 
// ///////////////////////////////////
// //  From schlatherind.c
// //
// void schlatherindfull(int *covmod, double *data, double *dist, int *nSite,
//                       int *nObs, int *dim, int *weighted, double *weights,
//                       double *locs, double *scales, double *shapes,
//                       double *alpha, double *nugget, double *range, double *smooth,
//                       double *smooth2, int *fitmarge,double *dns);
// void schlatherinddsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs, int *dim,
//                          int *weighted, double *weights, double *locdsgnmat, double *locpenmat,
//                          int *nloccoeff, int *npparloc, double *locpenalty, double *scaledsgnmat,
//                          double *scalepenmat, int *nscalecoeff, int *npparscale,
//                          double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
//                          int *nshapecoeff, int *npparshape, double *shapepenalty, int *usetempcov,
//                          double *tempdsgnmatloc, double *temppenmatloc, int *ntempcoeffloc,
//                          int *nppartempcoeffloc, double *temppenaltyloc, double *tempdsgnmatscale,
//                          double *temppenmatscale, int *ntempcoeffscale, int *nppartempcoeffscale,
//                          double *temppenaltyscale, double *tempdsgnmatshape, double *temppenmatshape,
//                          int *ntempcoeffshape, int *nppartempcoeffshape, double *temppenaltyshape,
//                          double *loccoeff, double *scalecoeff, double *shapecoeff,
//                          double *tempcoeffloc, double *tempcoeffscale, double *tempcoeffshape,
//                          double *alpha, double *nugget, double *range, double *smooth,
//                          double *smooth2, double *dns);
// 
// ///////////////////////////////////
// //  From geomgauss.c
// //
// void geomgaussfull(int *covmod, double *data, double *dist, int *nSite, int *nObs, int *dim,
//                    int *weighted, double *weights, double *locs, double *scales, double *shapes,
//                    double *sigma2, double *sigma2Bound, double *nugget, double *range,
//                    double *smooth, double *smooth2, int *fitmarge,double *dns);
// void geomgaussdsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs, int *dim,
//                       int *weighted, double *weights, double *locdsgnmat, double *locpenmat,
//                       int *nloccoeff, int *npparloc, double *locpenalty, double *scaledsgnmat,
//                       double *scalepenmat, int *nscalecoeff, int *npparscale,
//                       double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
//                       int *nshapecoeff, int *npparshape, double *shapepenalty, int *usetempcov,
//                       double *tempdsgnmatloc, double *temppenmatloc, int *ntempcoeffloc,
//                       int *nppartempcoeffloc, double *temppenaltyloc, double *tempdsgnmatscale,
//                       double *temppenmatscale, int *ntempcoeffscale, int *nppartempcoeffscale,
//                       double *temppenaltyscale, double *tempdsgnmatshape, double *temppenmatshape,
//                       int *ntempcoeffshape, int *nppartempcoeffshape, double *temppenaltyshape,
//                       double *loccoeff, double *scalecoeff, double *shapecoeff,
//                       double *tempcoeffloc, double *tempcoeffscale, double *tempcoeffshape,
//                       double *sigma2, double *sigma2Bound, double *nugget, double *range,
//                       double *smooth, double *smooth2, double *dns);
// 
// ///////////////////////////////////
// //  From brownResnick.c
// //
// void brownresnickfull(double *data, double *dist, int *nSite, int *nObs, int *weighted,
//                       double *weights, double *locs, double *scales, double *shapes,
//                       double *range, double *smooth, int *fitmarge, double *dns);
// void brownresnickdsgnmat(double *data, double *dist, int *nSite, int *nObs, int *weighted,
//                          double *weights, double *locdsgnmat, double *locpenmat, int *nloccoeff,
//                          int *npparloc, double *locpenalty, double *scaledsgnmat, double *scalepenmat,
//                          int *nscalecoeff, int *npparscale, double *scalepenalty, double *shapedsgnmat,
//                          double *shapepenmat, int *nshapecoeff, int *npparshape, double *shapepenalty,
//                          int *usetempcov, double *tempdsgnmatloc, double *temppenmatloc,
//                          int *ntempcoeffloc, int *nppartempcoeffloc, double *temppenaltyloc,
//                          double *tempdsgnmatscale, double *temppenmatscale, int *ntempcoeffscale,
//                          int *nppartempcoeffscale, double *temppenaltyscale, double *tempdsgnmatshape,
//                          double *temppenmatshape, int *ntempcoeffshape, int *nppartempcoeffshape,
//                          double *temppenaltyshape, double *loccoeff, double *scalecoeff,
//                          double *shapecoeff, double *tempcoeffloc, double *tempcoeffscale,
//                          double *tempcoeffshape,double *range, double *smooth, double *dns);
// 
// ///////////////////////////////////
// //  From smith.c
// //
// void smithfull(double *data, double *distVec, int *nSite, int *nObs, int *weighted, double *weights,
//                double *locs, double *scales, double *shapes, double *cov11, double *cov12,
//                double *cov22, int *fitmarge, double *dns);
// void smithdsgnmat(double *data, double *distVec, int *nSite, int *nObs, int *weighted,
//                   double *weights, double *locdsgnmat, double *locpenmat, int *nloccoeff,
//                   int *npparloc, double *locpenalty, double *scaledsgnmat,
//                   double *scalepenmat, int *nscalecoeff, int *npparscale,
//                   double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
//                   int *nshapecoeff, int *npparshape, double *shapepenalty,
//                   int *usetempcov, double *tempdsgnmatloc, double *temppenmatloc,
//                   int *ntempcoeffloc, int *nppartempcoeffloc, double *temppenaltyloc,
//                   double *tempdsgnmatscale, double *temppenmatscale, int *ntempcoeffscale,
//                   int *nppartempcoeffscale, double *temppenaltyscale, double *tempdsgnmatshape,
//                   double *temppenmatshape, int *ntempcoeffshape, int *nppartempcoeffshape,
//                   double *temppenaltyshape, double *loccoeff, double *scalecoeff,
//                   double *shapecoeff, double *tempcoeffloc, double *tempcoeffscale,
//                   double *tempcoeffshape, double *cov11, double *cov12, double *cov22,
//                   double *dns);
// 
// ///////////////////////////////////
// //  From smith3d.c
// //
// void smithfull3d(double *data, double *distVec, int *nSite, int *nObs, int *weighted,
//                  double *weights, double *locs, double *scales, double *shapes,
//                  double *cov11, double *cov12, double *cov13, double *cov22,
//                  double *cov23, double *cov33, int *fitmarge, double *dns);
// void smithdsgnmat(double *data, double *distVec, int *nSite, int *nObs, int *weighted,
//                   double *weights, double *locdsgnmat, double *locpenmat, int *nloccoeff,
//                   int *npparloc, double *locpenalty, double *scaledsgnmat,
//                   double *scalepenmat, int *nscalecoeff, int *npparscale,
//                   double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
//                   int *nshapecoeff, int *npparshape, double *shapepenalty,
//                   int *usetempcov, double *tempdsgnmatloc, double *temppenmatloc,
//                   int *ntempcoeffloc, int *nppartempcoeffloc, double *temppenaltyloc,
//                   double *tempdsgnmatscale, double *temppenmatscale, int *ntempcoeffscale,
//                   int *nppartempcoeffscale, double *temppenaltyscale, double *tempdsgnmatshape,
//                   double *temppenmatshape, int *ntempcoeffshape, int *nppartempcoeffshape,
//                   double *temppenaltyshape, double *loccoeff, double *scalecoeff,
//                   double *shapecoeff, double *tempcoeffloc, double *tempcoeffscale,
//                   double *tempcoeffshape, double *cov11, double *cov12, double *cov22,
//                   double *dns);
// 
// ///////////////////////////////////
// //  From extremalt.c
// //
// void extremaltfull(int *covmod, double *data, double *dist, int *nSite, int *nObs,
//                    int *dim, int *weighted, double *weights, double *locs, double *scales,
//                    double *shapes, double *nugget, double *range, double *smooth, double *smooth2,
//                    double *df, int *fitmarge, double *dns);
// void extremaltdsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs, int *dim,
//                       int *weighted, double *weights, double *locdsgnmat, double *locpenmat,
//                       int *nloccoeff, int *npparloc, double *locpenalty, double *scaledsgnmat,
//                       double *scalepenmat, int *nscalecoeff, int *npparscale,
//                       double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
//                       int *nshapecoeff, int *npparshape, double *shapepenalty, int *usetempcov,
//                       double *tempdsgnmatloc, double *temppenmatloc, int *ntempcoeffloc,
//                       int *nppartempcoeffloc, double *temppenaltyloc, double *tempdsgnmatscale,
//                       double *temppenmatscale, int *ntempcoeffscale, int *nppartempcoeffscale,
//                       double *temppenaltyscale, double *tempdsgnmatshape, double *temppenmatshape,
//                       int *ntempcoeffshape, int *nppartempcoeffshape, double *temppenaltyshape,
//                       double *loccoeff, double *scalecoeff, double *shapecoeff,
//                       double *tempcoeffloc, double *tempcoeffscale, double *tempcoeffshape,
//                       double *nugget, double *range, double *smooth, double *smooth2, double *df,
//                       double *dns);
// 
///////////////////////////////////
//  From utils.c
//
void distance(double *coord, int *nDim, int *nSite, int *vec, double *dist);
void distance2orig(double *coord, int n, int dim, double *dist, int grid);
double gev2frech(double *data, int nObs, int nSite, double *locs,
                 double *scales, double *shapes, double *jac, double *frech);
double gev2frechTrend(double *data, int nObs, int nSite, double *locs, double *scales,
                      double *shapes, double *trendlocs, double *trendscales,
                      double *trendshapes,double *jac, double *frech);
double dsgnmat2Param(double *locdsgnmat, double *scaledsgnmat,
                     double *shapedsgnmat, double *loccoeff,
                     double *scalecoeff, double *shapecoeff,
                     int nSite, int nloccoeff, int nscalecoeff,
                     int nshapecoeff, double *locs, double *scales,
                     double *shapes);
double dsgnmat2Param2(double *locdsgnmat, double *scaledsgnmat, double *shapedsgnmat,
                      double *tempdsgnmatLoc, double *tempdsgnmatScale, double *tempdsgnmatShape,
                      double *loccoeff, double *scalecoeff, double *shapecoeff,
                      double *tempcoeffLoc, double *tempcoeffScale, double *tempcoeffShape,
                      int nSite, int nObs, int *usetempcov, int nloccoeff, int nscalecoeff,
                      int nshapecoeff, int ntempcoeffLoc, int ntempcoeffScale,
                      int ntempcoeffShape, double *locs, double *scales, double *shapes);
void dsgnmat2temptrend(double *dsgnmatloc, double *dsgnmatscale, double *dsgnmatshape,
                       double *loccoeff, double *scalecoeff, double *shapecoeff, int nSite,
                       int nObs, int *usetempcov, int nloccoeff, int nscalecoeff,
                       int nshapecoeff, double *trendlocs, double *trendscales,
                       double *trendshapes);
void dsgnmat2Alpha(double *alphadsgnmat, double *alphacoeff,
                   int nSite, int nalphacoeff, double *alphas);
void dsgnmat2Sigma2(double *sigma2dsgnmat, double *sigma2coeff,
                    int nSite, int nsigma2coeff, double *sigma2);
void gev(double *prob, int *n, double *locs, double *scales, double *shapes,
         double *quant);
double gev2unif(double *data, int nObs, int nSite, double *locs,
                double *scales, double *shapes, double *unif,
                double *logdens);
double gev2unifTrend(double *data, int nObs, int nSite, double *locs,
                     double *scales, double *shapes, double *trendlocs,
                     double *trendscales, double *trendshapes, double *unif,
                     double *logdens);
int getCurrentPair(int site1, int site2, int nSite);
void getSiteIndex(int currentPair, int nSite, int *site1, int *site2);

// ///////////////////////////////////
// //  From univllik.c
// //
// void gevlik(double *data, int *n, double *loc, double *scale,
//             double *shape, double *dns);
// void gpdlik(double *exceed, int *n, double *thresh, double *scale,
//             double *shape, double *dns);
// 
// ///////////////////////////////////
// //  From covariance.c
// //
double whittleMatern(double *dist, int n, double nugget, double sill, double range,
                     double smooth, double *rho);
double cauchy(double *dist, int n, double nugget, double sill, double range,
              double smooth, double *rho);
double caugen(double *dist, int n, double nugget, double sill, double range,
              double smooth, double smooth2, double *rho);
double powerExp(double *dist, int n, double nugget, double sill, double range,
                double smooth, double *rho);
double bessel(double *dist, int n, int dim, double nugget, double sill,
              double range, double smooth, double *rho);
double mahalDistFct(double *distVec, int n, double *cov11,
                    double *cov12, double *cov22, double *mahal);
double mahalDistFct3d(double *distVec, int n, double *cov11,
                      double *cov12, double *cov13, double *cov22,
                      double *cov23, double *cov33, double *mahal);
double geomCovariance(double *dist, int n, int dim, int covmod,
                      double sigma2, double sigma2Bound, double nugget,
                      double range, double smooth, double smooth2,
                      double *rho);
double brownResnick(double *dist, int n, double range, double smooth,
                    double *rho);
double fbm(double *coord, double *dist, int dim, int nSite, double sill, double range,
           double smooth, double *rho);

// ///////////////////////////////////
// //  From mcmc.c
// //
// SEXP gibbs(SEXP n, SEXP np, SEXP thin, SEXP init,
//            SEXP psd, SEXP f, SEXP rho);
// 
// ///////////////////////////////////
// //  From pairwiselik.c
// //
// double lpliksmith(double *data, double *rho, double *jac,
//                   int nObs, int nSite);
// double lplikschlather(double *data, double *rho, double *jac,
//                       int nObs, int nSite);
// double lplikschlatherind(double *data, double alpha, double *rho,
//                          double *jac, int nObs, int nSite);
// double lplikextremalt(double *data, double *rho, double df, double *jac,
//                       int nObs, int nSite);
// 
// ///////////////////////////////////
// //  From weightedPairwiselik.c
// //
// double wlplikschlather(double *data, double *rho, double *jac,
//                        int nObs, int nSite, double *weights);
// double wlpliksmith(double *data, double *mahalDist, double *jac,
//                    int nObs, int nSite, double *weights);
// double wlplikschlatherind(double *data, double alpha, double *rho,
//                           double *jac, int nObs, int nSite, double *weights);
// double wlplikextremalt(double *data, double *rho, double df, double *jac,
//                        int nObs, int nSite, double *weights);
// 
// ///////////////////////////////////
// //  From penalizations.c
// //
// double penalization(double *penmat, double *beta, double pencoeff, int n,
//                     int nppar);
// double penalization2(double *penmat, double *beta, double pencoeff, int n,
//                      int nppar);
// 
// 
// ///////////////////////////////////
// //  From extcoeff.c
// //
// void extCoeffSmith(double *frech, int *nObs, int *nSite,
//                    double *extCoeff);
// void extCoeffST(double *frech, double *xBar, double *z, double *theta,
//                 int *nObs, double *dns);
// 
// ///////////////////////////////////
// //  From fitcovmat.c
// //
// void fitcovmat2d(double *cov11, double *cov12, double *cov22,
//                  int *nPairs, double *dist, double *extcoeff,
//                  double *weights, double *ans);
// void fitcovmat3d(double *cov11, double *cov12, double *cov13,
//                  double *cov22, double *cov23, double *cov33,
//                  int *nPairs, double *dist, double *extcoeff,
//                  double *weights, double *ans);
// void fitcovariance(int *covmod, double *nugget, double *range, double *smooth,
//                    double *smooth2, int *nPairs, int *dim, double *distVec,
//                    double *extcoeff, double *weights, double *ans);
// void fittcovariance(int *covmod, double *nugget, double *range, double *smooth,
//                     double *smooth2, double *DoF, int *nPairs, int *dim, double *dist,
//                     double *extcoeff, double *weights, double *ans);
// void fiticovariance(int *covmod, double *alpha, double *nugget, double *range,
//                     double *smooth, double *smooth2, int *nPairs, int *dim,
//                     double *dist, double *extcoeff, double *weights, double *ans);
// void fitgcovariance(int *covmod, double *sigma2, double *sigma2Bound, double *nugget,
//                     double *range, double *smooth, double *smooth2, int *nPairs,
//                     int *dim, double *dist, double *extcoeff, double *weights,
//                     double *ans);
// void fitbrcovariance(double *range, double *smooth, int *nPairs,
//                      double *dist, double *extcoeff, double *weights,
//                      double *ans);
// 
// ///////////////////////////////////
// //  From spatgevlik.c
// //
// void spatgevlik(double *data, double *covariables, int *nSite, int *nObs,
//                 double *locdsgnmat, double *locpenmat, int *nloccoeff,
//                 int *npparloc, double *locpenalty, double *scaledsgnmat,
//                 double *scalepenmat, int *nscalecoeff, int *npparscale,
//                 double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
//                 int *nshapecoeff, int *npparshape, double *shapepenalty,
//                 int *usetempcov, double *tempdsgnmatLoc, double *temppenmatLoc,
//                 int *ntempcoeffLoc, int *nppartempcoeffLoc, double *temppenaltyLoc,
//                 double *tempdsgnmatScale, double *temppenmatScale, int *ntempcoeffScale,
//                 int *nppartempcoeffScale, double *temppenaltyScale, double *tempdsgnmatShape,
//                 double *temppenmatShape, int *ntempcoeffShape, int *nppartempcoeffShape,
//                 double *temppenaltyShape, double *loccoeff, double *scalecoeff,
//                 double *shapecoeff, double *tempcoeffLoc, double *tempcoeffScale,
//                 double *tempcoeffShape, double *dns);
// 
// ///////////////////////////////////
// //  From madogram.c
// //
// void madogram(double *data, int *nObs, int *nSite, double *mado);
// void variogram(double *data, int *nObs, int *nSite, double *vario);
// void lmadogram(double *data, int *nObs, int *nSite, double *lambda,
//                int *nLambda, double *lmado);
// 
// ///////////////////////////////////
// //  From simsmith.c
// //
// void rsmith1d(double *coord, double *center, double *edge, int *nObs,
//               int *nSites, double *var, double *ans);
// void rsmith2d(double *coord, double *center, double *edge, int *nObs,
//               int *nSites, int *grid, double *cov11, double *cov12,
//               double *cov22, double *ans);
// 
// ///////////////////////////////////
// //  From direct.c
// //
// void buildcovmat(int *nSite, int *grid, int *covmod, double *coord, int *dim,
//                  double *nugget, double *sill, double *range, double *smooth,
//                  double *covMat);
// void direct(int *n, int *nSite, int *grid, int *covmod, double *coord, int *dim,
//             double *nugget, double *sill, double *range, double *smooth,
//             double *ans);
// 
// ///////////////////////////////////
// //  From randomlines.c
// //
// void vandercorput(int *n, double *coord);
// void rotation(double *coord, int *n, double *u, double *v, double *w,
//               double *angle);
// 
// ///////////////////////////////////
// //  From turningbands.c
// //
// void tbm(int *nobs, int *nsite, int *dim, int *covmod, int *grid,
//          double *coord, double *nugget, double *sill, double *range,
//          double *smooth, int *nlines, double *ans);
// 
// ///////////////////////////////////
// //  From simschlather.c
// //
// void rschlathertbm(double *coord, int *nObs, int *nSites, int *dim,
//                    int *covmod, int *grid, double *nugget, double *range,
//                    double *smooth, double *uBound, int *nlines,
//                    double *ans);
// void rschlatherdirect(double *coord, int *nObs, int *nSites, int *dim,
//                       int *covmod, int *grid, double *nugget, double *range,
//                       double *smooth, double *uBound, double *ans);
// void rschlathercirc(int *nObs, int *ngrid, double *steps, int *dim,
//                     int *covmod, double *nugget, double *range,
//                     double *smooth, double *uBound, double *ans);
// void tbmcore(int *nsite, int *neffSite, int *dim, int *covmod,
//              int *grid, double *coord, double *nugget, double *sill,
//              double *range, double *smooth, int *nlines, double *lines,
//              double *ans);
// void circcore(double *rho, double *a, double *ia, int m, int halfM, int mdag,
//               int mdagbar, int ngrid, int nbar, double isqrtMbar, double nugget,
//               double *ans);
// 
// ///////////////////////////////////
// //  From simgeometric.c
// //
// void rgeomtbm(double *coord, int *nObs, int *nSite, int *dim,
//               int *covmod, int *grid, double *sigma2, double *nugget,
//               double *range, double *smooth, double *uBound, int *nlines,
//               double *ans);
// void rgeomdirect(double *coord, int *nObs, int *nSite, int *dim,
//                  int *covmod, int *grid, double *sigma2, double *nugget,
//                  double *range, double *smooth, double *uBound,
//                  double *ans);
// void rgeomcirc(int *nObs, int *ngrid, double *steps, int *dim,
//                int *covmod, double *sigma2, double *nugget, double *range,
//                double *smooth, double *uBound, double *ans);
// 
// ///////////////////////////////////
// //  From simBrownResnick.c
// //
// void rbrowndirect(double *coord, double *bounds, int *nObs, int *nSite,
//                   int *dim, int *grid, double *range, double *smooth,
//                   double *uBound, int *method, int *maxSim, int *nPP,
//                   int *idxsubOrig, int *nsubOrig, double *ans);
// 
// ///////////////////////////////////
// //  From simextremalt.c
// //
// void rextremalttbm(double *coord, int *nObs, int *nSite, int *dim,
//                    int *covmod, int *grid, double *nugget, double *range,
//                    double *smooth, double *DoF, double *uBound, int *nlines,
//                    double *ans);
// void rextremaltdirect(double *coord, int *nObs, int *nSite, int *dim,
//                       int *covmod, int *grid, double *nugget, double *range,
//                       double *smooth, double *DoF, double *uBound, double *ans);
// void rextremaltcirc(int *nObs, int *ngrid, double *steps, int *dim,
//                     int *covmod, double *nugget, double *range,
//                     double *smooth, double *DoF, double *uBound, double *ans);
// 
// ///////////////////////////////////
// //  From gpdproc.c
// //
// void gpdprocfull(double *data, double *distVec, int *nSite,
//                  int *nObs, double *excRates, double *threshs, double *scales,
//                  double *shapes, double *cov11, double *cov12,
//                  double *cov22, int *fitmarge, double *dns);
// double gpd2ugpd(double *data, int nObs, int nSite, double *excRates,
//                 double *threshs, double *scales, double *shapes,
//                 double *jac, double *ugpd);
// double lpliksmithgpd(double *data, double *mahalDist, double *jac,
//                      double *excRates, int nObs, int nSite);
// 
// ///////////////////////////////////
// //  From standardErrors.c
// //
// void smithstderr(double *data, double *distVec, int *nSite, int *nObs, double *locdsgnmat,
//                  int *nloccoeff, double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
//                  int *nshapecoeff, double *tempdsgnmatloc, int *ntemploccoeff,
//                  double *tempdsgnmatscale, int *ntempscalecoeff, double *tempdsgnmatshape,
//                  int *ntempshapecoeff,double *loccoeff, double *scalecoeff, double *shapecoeff,
//                  double *temploccoeff, double *tempscalecoeff, double *tempshapecoeff,
//                  double *cov11, double *cov12, double *cov22, int *fitmarge, int *usetempcov,
//                  double *weights, double *hess, double *grad);
// void smithstderr3d(double *data, double *distVec, int *nSite, int *nObs, double *locdsgnmat,
//                    int *nloccoeff, double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
//                    int *nshapecoeff, double *tempdsgnmatloc, int *ntemploccoeff,
//                    double *tempdsgnmatscale, int *ntempscalecoeff, double *tempdsgnmatshape,
//                    int *ntempshapecoeff, double *loccoeff, double *scalecoeff, double *shapecoeff,
//                    double *temploccoeff, double *tempscalecoeff, double *tempshapecoeff, double *cov11,
//                    double *cov12, double *cov13, double *cov22, double *cov23, double *cov33, int *fitmarge,
//                    int *usetempcov, double *weights, double *hess, double *grad);
// void schlatherstderr(int *covmod, double *data, double *dist, int *nSite, int *nObs,
//                      double *locdsgnmat, int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
//                      double *shapedsgnmat, int *nshapecoeff, double *tempdsgnmatloc,
//                      int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
//                      double *tempdsgnmatshape, int *ntempshapecoeff, double *loccoeff,
//                      double *scalecoeff, double *shapecoeff, double *temploccoeff,
//                      double *tempscalecoeff, double *tempshapecoeff, double *nugget, double *range,
//                      double *smooth, double *smooth2, int *fitmarge, int *usetempcov, double *weights,
//                      double *hess, double *grad);
// void schlatherindstderr(int *covmod, double *data, double *dist, int *nSite, int *nObs,
//                         double *locdsgnmat, int *nloccoeff, double *scaledsgnmat,
//                         int *nscalecoeff, double *shapedsgnmat,	int *nshapecoeff,
//                         double *tempdsgnmatloc, int *ntemploccoeff, double *tempdsgnmatscale,
//                         int *ntempscalecoeff, double *tempdsgnmatshape, int *ntempshapecoeff,
//                         double *loccoeff, double *scalecoeff, double *shapecoeff,
//                         double *temploccoeff, double *tempscalecoeff, double *tempshapecoeff,
//                         double *alpha, double *nugget, double *range, double *smooth,
//                         double *smooth2, int *fitmarge, int *usetempcov, double *weights, double *hess,
//                         double *grad);
// void geomgaussstderr(int *covmod, double *data, double *dist, int *nSite, int *nObs,
//                      double *locdsgnmat, int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
//                      double *shapedsgnmat, int *nshapecoeff,  double *tempdsgnmatloc,
//                      int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
//                      double *tempdsgnmatshape, int *ntempshapecoeff, double *loccoeff,
//                      double *scalecoeff, double *shapecoeff, double *temploccoeff,
//                      double *tempscalecoeff, double *tempshapecoeff, double *sigma2,
//                      double *nugget, double *range, double *smooth, double *smooth2,
//                      int *fitmarge, int *usetempcov, double *weights, double *hess, double *grad);
// void brownresnickstderr(double *data, double *dist, int *nSite, int *nObs, double *locdsgnmat,
//                         int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
//                         double *shapedsgnmat, int *nshapecoeff, double *tempdsgnmatloc,
//                         int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
//                         double *tempdsgnmatshape, int *ntempshapecoeff, double *loccoeff,
//                         double *scalecoeff, double *shapecoeff, double *temploccoeff,
//                         double *tempscalecoeff, double *tempshapecoeff, double *range,
//                         double *smooth, int *fitmarge, int *usetempcov, double *weights,
//                         double *hess, double *grad);
// void spatgevstderr(double *data, int *nSite, int *nObs, double *locdsgnmat,
//                    int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
//                    double *shapedsgnmat, int *nshapecoeff, double *tempdsgnmatloc,
//                    int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
//                    double *tempdsgnmatshape, int *ntempshapecoeff,  double *loccoeff,
//                    double *scalecoeff, double *shapecoeff, double *temploccoeff,
//                    double *tempscalecoeff, double *tempshapecoeff, int *usetempcov,
//                    double *hess, double *grad);
// void extremaltstderr(int *covmod, double *data, double *dist, int *nSite, int *nObs,
//                      double *locdsgnmat, int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
//                      double *shapedsgnmat, int *nshapecoeff, double *tempdsgnmatloc,
//                      int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
//                      double *tempdsgnmatshape, int *ntempshapecoeff, double *loccoeff,
//                      double *scalecoeff, double *shapecoeff, double *temploccoeff,
//                      double *tempscalecoeff, double *tempshapecoeff, double *nugget, double *range,
//                      double *smooth, double *smooth2, double *df, int *fitmarge, int *usetempcov,
//                      double *weights, double *hess, double *grad);
// 
// ///////////////////////////////////
// //  From standardErrorsCommonPart.c
// //
// void marginalPartSmith(int *start, int *nObs, int *nSite, double *data, double *frech,
//                        double *mahalDist, double *locs, double *scales, double *shapes,
//                        double *trendlocs, double *trendscales, double *trendshapes,
//                        int *nloccoeff, int *nscalecoeff, int *nshapecoeff, int *ntemploccoeff,
//                        int *ntempscalecoeff, int *ntempshapecoeff, double *locdsgnmat,
//                        double *scaledsgnmat, double *shapedsgnmat, double *tempdsgnmatloc,
//                        double *tempdsgnmatscale, double *tempdsgnmatshape, double *weights,
//                        double *hess, double *grad);
// void marginalPartSchlat(int *start, int *nObs, int *nSite, double *data, double *frech,
//                         double *rho, double *locs, double *scales, double *shapes,
//                         double *trendlocs, double *trendscales, double *trendshapes,
//                         int *nloccoeff, int *nscalecoeff, int *nshapecoeff, int *ntemploccoeff,
//                         int *ntempscalecoeff, int *ntempshapecoeff, double *locdsgnmat,
//                         double *scaledsgnmat, double *shapedsgnmat, double *tempdsgnmatloc,
//                         double *tempdsgnmatscale, double *tempdsgnmatshape, double *weights,
//                         double *hess, double *grad);
// void marginalPartiSchlat(int *start, int *nObs, int *nSite, double *data, double *frech,
//                          double *alpha, double *rho, double *locs, double *scales, double *shapes,
//                          double *trendlocs, double *trendscales, double *trendshapes,
//                          int *nloccoeff, int *nscalecoeff, int *nshapecoeff, int *ntemploccoeff,
//                          int *ntempscalecoeff, int *ntempshapecoeff, double *locdsgnmat,
//                          double *scaledsgnmat, double *shapedsgnmat, double *tempdsgnmatloc,
//                          double *tempdsgnmatscale, double *tempdsgnmatshape, double *weights,
//                          double *hess, double *grad);
// void marginalPartExtremalt(int *start, int *nObs, int *nSite, double *data, double *frech,
//                            double *df, double *rho, double *locs, double *scales, double *shapes,
//                            double *trendlocs, double *trendscales, double *trendshapes,
//                            int *nloccoeff, int *nscalecoeff, int *nshapecoeff, int *ntemploccoeff,
//                            int *ntempscalecoeff, int *ntempshapecoeff, double *locdsgnmat,
//                            double *scaledsgnmat, double *shapedsgnmat, double *tempdsgnmatloc,
//                            double *tempdsgnmatscale, double *tempdsgnmatshape, double *weights,
//                            double *hess, double *grad);
// 
// ///////////////////////////////////
// //  From circulant.c
// //
// void circemb(int *nsim, int *ngrid, double *steps, int *dim, int *covmod,
//              double *nugget, double *sill, double *range, double *smooth,
//              double *ans);
// 
// ///////////////////////////////////
// //  From latentVariable.c
// //
// void latentgev(int *n, double *data, int *nSite, int *nObs, int *covmod,
//                int *dim, double *distMat, double *dsgnMat, int *nBeta, double *beta,
//                double *sills, double *ranges, double *smooths, double *gevParams,
//                double *hyperSill, double *hyperRange, double *hyperSmooth,
//                double *hyperBetaMean, double *hyperBetaIcov, double *propGev,
//                double *propRanges, double *propSmooths, double *mcLoc,
//                double *mcScale, double *mcShape, double *accRates,
//                double *extRates, int *thin, int *burnin);
// void DIC(int *n, int *nSite, int *nObs, double *data, double *chainLoc,
//          double *chainScale, double *chainShape, double *postLoc,
//          double *postScale, double *postShape, double *dic, double *effNpar,
//          double *dbar);
// 
///////////////////////////////////
//  From kriging.c
//
void skriging(int *nSite, int *nSiteKrig, int *covmod, int *dim,
              double *icovMat, double *coord, double *coordKrig, double *obs,
              double *sill, double *range, double *smooth, double *smooth2,
              double *weights);

// ///////////////////////////////////
// //  From copula.c
// //
// void copula(int *copula, int *covmod, double *dist, double *data, int *nSite, int *nObs,
//             int *dim, int * fitmarge, double *locdsgnmat,  double *locpenmat, int *nloccoeff,
//             int *npparloc, double *locpenalty, double *scaledsgnmat, double *scalepenmat,
//             int *nscalecoeff, int *npparscale, double *scalepenalty,
//             double *shapedsgnmat, double *shapepenmat, int *nshapecoeff, int *npparshape,
//             double *shapepenalty, int *usetempcov, double *tempdsgnmatloc,
//             double *temppenmatloc, int *ntempcoeffloc, int *nppartempcoeffloc,
//             double *temppenaltyloc, double *tempdsgnmatscale, double *temppenmatscale,
//             int *ntempcoeffscale, int *nppartempcoeffscale, double *temppenaltyscale,
//             double *tempdsgnmatshape, double *temppenmatshape, int *ntempcoeffshape,
//             int *nppartempcoeffshape, double *temppenaltyshape, double *loccoeff,
//             double *scalecoeff, double *shapecoeff, double *tempcoeffloc,
//             double *tempcoeffscale, double *tempcoeffshape, double *DoF, double *nugget, double *range,
//             double *smooth, double *smooth2, double *dns);
// double gaussianCopula(double *unif, double sd, double *covMat, int nObs, int nSite);
// double studentCopula(double *data, double DoF, double *covMat, int nObs,
//                      int nSite);
// 
// ///////////////////////////////////
// //  From maxLinear.c
// //
// void rcondMaxLin(double *data, double *dsgnMat, int *p, int *nSite, int *nSim,
//                  double *Z);
// void maxLinear(int *nSim, double *dsgnMat, double *Z, int *nSite, int *p,
//                int *grid, double *sim);
// void maxLinDsgnMat(double *coord, double *grid, int *nSite, int *p,
//                    double *areaPixel, int *dim, double *param, double *dsgnMat);
// 
// ///////////////////////////////////
// //  From condsimMaxStab.c
// //
// // Utilities functions
// void getSubMatrix(double *mat, int *dim, int *nr, int *rows, int *nc,
//                   int *cols, double *submat);
// void listAllPartOfASet(int *n, int *nPart, int *allPart, int *allSize);
// void stirling2ndKind(int *n, int *k, double *ans);
// void bell(int *n, int *ans);
// void gettau(int *nCond, int *part, int *set, int *tau);
// void gettaubar(int *nCond, int *part, int *set, int *taubar);
// void standardize(double *quant, double *cov, double *mean, int *n);
// void pmvnorm(double *bounds, int *n, double *cor, double *prob,
//              double *err, int *nMc);
// void pmvt(double *bounds, int *n, double *DoF, double *mu, double *scaleMat,
//           double *prob, double *err, int *nMc);
// void convert2rightformat(int *partition, int *n, int *size);
// void sampleDiscreteDist(int *n, double *prob, int *ans);
// void validPart(int *partition, int *n, int *valid);
// void getBounds(int *partition, int *n, int *idx, int *lbound, int *ubound);
// int getPartSize(int *partition, int *n);
// 
// // Functions for conditional simulations of Schlather processes
// void getStartingPartitionSC(int *nsim, int *n, double *covChol, int *startPart);
// void gibbsForPartSC(int *nchain, int *nthin, int *burnin, int *nCond,
//                     int *currentPart, double *cov, double *y, int *chain,
//                     double *timings);
// void getParametersSC(int *tau, int *taubar, int *ntau, int *ntaubar, double *cov,
//                      double *y, double *DoF, double *mu, double *scaleMat);
// void getfvaluesSC(double *y, int *n, int *ntau, int *tau, double *cov,
//                   double *f);
// void computeprobaSC(double *DoF, double *mu, double *scaleMat, double *y,
//                     int *ntaubar, int *taubar, double *prob);
// void computeWeightsSC(int *nCond, double *y, int *nPart, int *allPart,
//                       int *allSize, double *cov, double *weights);
// double computeWeightForOneSetSC(int *idxSet, int *nCond, int *partition, double *cov,
//                                 double *y);
// void condsimschlather(int *nsim, int *n, int *nCond, int *allPart, double *cov,
//                       double *y, double *sim, double *subextfct,
//                       double *extfct, double *timings);
// 
// // Functions for conditional simulations of Brown-Resnick processes
// void computeWeightsBR(int *nCond, double *y, int *nPart, int *allPart, int *allSize,
//                       double *cov, double *sigma2, double *covChol, double *ham,
//                       double *mean1, double *weights);
// void computeprobaBR(double *icovChol, double *mu, double *y, int *n,
//                     int *r, int *nCond, int *taubar, double *prob);
// void getfvaluesBR(double *y, double *sigma2, double *covjchol, int *r,
//                   double *f);
// void buildJ(int *tau, int *n, int *r, double *J);
// void buildJtilde(int *tau, int *n, int *r, double *Jtilde);
// void getParametersBR(double *J, double *Jtilde, int *n, int *nr,
//                      double *covChol, double *ham, double *mean1,
//                      double *ytilde, double *iBchol, double *mu);
// void condsimbrown(int *nsim, int *n, int *nCond, int *allPart, double *covChol,
//                   double *sigma2, double *ham, double *mean1, double *ytilde, double *coord,
//                   double *range, double *smooth, int *dim,  double *sim, double *subextfct,
//                   double *extfct, double *timings);
// void condsimbrown2(int *nsim, int *n, int *nCond, int *allPart, double *covChol,
//                    double *sigma2, double *ham, double *mean1, double *ytilde, double *sim,
//                    double *coord, double *range, double *smooth, double *xlim);
// double computeWeightForOneSetBR(int *idxSet, int *nCond, int *partition, double *cov,
//                                 double *sigma2, double *covChol, double *ham, double *mean1,
//                                 double *y);
// void getStartingPartitionBR(int *nSim, int *n, double *coord, double *range, double *smooth,
//                             int *startPart);
// void gibbsForPartBR(int *nchain, int *nthin, int *burnin, int *nCond,
//                     int *currentPart, double *cov, double *sigma2,
//                     double *covChol, double *ham, double *mean1, double *y,
//                     int *chain, double *timing);
// 
// // Functions for conditional simulations of extremal-t processes
// void computeWeightsExtt(int *nCond, double *y, int *nPart, int *allPart,
//                         int *allSize, double *cov, double *nu, double *weights);
// void computeprobaExtt(double *nu, double *DoF, double *mu, double *scaleMat, double *y,
//                       int *ntaubar, int *taubar, double *prob);
// void getfvaluesExtt(double *y, int *n, int *ntau, int *tau, double *cov,
//                     double *nu, double *f);
// void getParametersExtt(int *tau, int *taubar, int *ntau, int *ntaubar, double *cov,
//                        double *y, double *nu, double *DoF, double *mu, double *scaleMat);
// void condsimextt(int *nsim, int *n, int *nCond, int *allPart, double *nu,
//                  double *cov, double *y, double *sim, double *subextfct,
//                  double *extfct, double *timings);
// void gibbsForPartExtt(int *nchain, int *nthin, int *burnin, int *nCond,
//                       int *currentPart, double *nu, double *cov, double *y, int *chain,
//                       double *timings);
// double computeWeightForOneSetExtt(int *idxSet, int *nCond, int *partition,
//                                   double *nu, double *cov, double *y);
// void totoExtt(int *nsim, int *n, double *nu, double *covChol, double *ans);
// void getStartingPartitionExtt(int *nsim, int *n, double *nu, double *covChol,
//                               int *startPart);
// 
// ///////////////////////////////////
// //  From fft.c
// //
// //static void fftmx(double *a, double *b, int ntot, int n, int nspan, int isn,
// //		  int m, int kt, double *at, double *ck, double *bt, double *sk,
// //		  int *np, int *nfac);
// void fft_factor(int n, int *pmaxf, int *pmaxp);
// Rboolean fft_work(double *a, double *b, int nseg, int n, int nspn, int isn,
//                   double *work, int *iwork);
// 
// 
// ///////////////////////////////////
// //  From completedLogLik.c
// //
// void completellik(int *nObs, int *nSite, double *data, int *parts, int *partSizes,
//                   double *cov, double *weights);
// void pschlather(double *q, int *dim, double *cov, double *prob, int *nMC);
// 
// 
// 
// 
// ///////////////////////////////////
// //  From concurrence.c
// //
// void empiricalConcProb(double *data, int *nSite, int *nObs, int *blockSize,
//                        int *nBlock, double *concProb);
// void empiricalBootConcProb(double *data, int *nSite, int *nObs, int *blockSize,
//                            double *concProb);
// void concProbKendall(double *data, int *nSite, int *nObs, double *concProb,
//                      double *jackKnife, int *computeStdErr);
// 
// ///////////////////////////////////
// //  From maxStableExactSim.c
// //
// void rbrownexact(double *coord, int *nObs, int *nSite, int *dim,
//                  int *grid, double *range, double *smooth,
//                  double *ans);
// void rextremaltexact(double *coord, int *nObs, int *nSite, int *dim,
//                      int *covmod, int *grid, double *nugget, double *range,
//                      double *smooth, double *DoF, double *ans);
// void rschlatherexact(double *coord, int *nObs, int *nSite, int *dim,
//                      int *covmod, int *grid, double *nugget, double *range,
//                      double *smooth, double *ans);