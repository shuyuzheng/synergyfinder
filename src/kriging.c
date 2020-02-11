# include "header.h"

void skriging(int *nSite, int *nSiteKrig, int *covmod, int *dim,
              double *icovMat, double *coord, double *coordKrig, double *obs,
              double *sill, double *range, double *smooth, double *smooth2,
              double *weights){
  
  /* This function computes the kriging weights using simple kriging,
   i.e., the mean is supposed to be 0. */
  
  double zero = 0, one = 1,
    *dist = malloc(*nSite * *nSiteKrig * sizeof(double)),
    *covariances = malloc(*nSite * *nSiteKrig * sizeof(double));
    
    for (int i=0;i<(*nSite * *nSiteKrig);i++)
      dist[i] = 0;
    
    /* 1. Compute the distances between the kriging locations and the
     locations where we got data */
    for (int i=0;i<*nSiteKrig;i++){
      for (int j=0;j<*nSite;j++){
        for (int k=0;k<*dim;k++)
          dist[j + i * *nSite] += (coord[j + k * *nSite] - coordKrig[i + k * *nSiteKrig]) *
            (coord[j + k * *nSite] - coordKrig[i + k * *nSiteKrig]);
        
        dist[j + i * *nSite] = sqrt(dist[j + i * *nSite]);
        
      }
    }
    
    // 2. Compute the covariance from these distances
    switch(*covmod){
    case 1:
      whittleMatern(dist, *nSite * *nSiteKrig, zero, *sill, *range, *smooth,
                    covariances);
      break;
    case 2:
      cauchy(dist, *nSite * *nSiteKrig, zero, *sill, *range, *smooth,
             covariances);
      break;
    case 3:
      powerExp(dist, *nSite * *nSiteKrig, zero, *sill, *range, *smooth,
               covariances);
      break;
    case 4:
      bessel(dist, *nSite * *nSiteKrig, *dim, zero, *sill, *range, *smooth,
             covariances);
      break;
    case 5:
      caugen(dist, *nSite * *nSiteKrig, zero, *sill, *range, *smooth, *smooth2,
             covariances);
      break;
    }
    
    /* 3. Compute the kriging weights i.e. weights = icovMat %*%
     covariances */
    F77_CALL(dsymm)("L", "U", nSite, nSiteKrig, &one, icovMat, nSite,
             covariances, nSite, &zero, weights, nSite FCONE FCONE);
    
    free(dist); free(covariances);
    return;
}