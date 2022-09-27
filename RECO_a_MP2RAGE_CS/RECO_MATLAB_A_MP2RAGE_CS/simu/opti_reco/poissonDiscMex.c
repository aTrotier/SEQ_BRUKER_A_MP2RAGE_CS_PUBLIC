#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <malloc.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>		/* For time() */
#include "poissonDiscMex.h"

#define FAILURE 1
#define SUCCESS 0
#define MAXVIEWS 2049
#define MAXSLICES 512
#define MAXACQS MAXVIEWS*MAXSLICES

/*----------------------------------Start Mex interface-------------------
 For the rest delete DB_MSG --------*/
#include "mex.h"


// function to shutdown with ctrl + c
#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
#else
    extern bool utIsInterruptPending();
#endif

// global parameters
int Genpoint;
double AccTotale;
        
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* function [mask] = vdPoisMex(sx,sy,fovx,fovy,accelx,accely,calib,ellipse)
    function [mask] = poissonDiscMex(fovy,fovz,sky,skz,ry,rz,ncal,cutcorners, pp)*/
  int err;
  int i;
  int csmask[MAXACQS];
  double *mask;
    
  // parameters
  float fovy;
  float fovz;
  int sky;
  int skz;
  float ry;
  float rz;
  int ncal;
  int cutcorners;
  float pp;

  /*	Check Arguments 	*/

  if (nrhs != 9)
	mexErrMsgTxt("Incorrect Input Argument List");

   if (nlhs != 1)
	mexErrMsgTxt("Function has only one output");

  fovy = *mxGetPr(prhs[0]);
  fovz = *mxGetPr(prhs[1]);
  sky =  *mxGetPr(prhs[2]);
  skz =  *mxGetPr(prhs[3]);
  ry =  *mxGetPr(prhs[4]);
  rz =  *mxGetPr(prhs[5]);
  ncal =  *mxGetPr(prhs[6]);
  cutcorners =  *mxGetPr(prhs[7]);
  pp =   *mxGetPr(prhs[8]);
  
 // Warning i inverse the sky and skz because the output of matlab is also inverted
int genPoint =  genVDPoissonSampling(csmask,  fovz,  fovy,  skz, sky,  ry, rz,  ncal, cutcorners,  pp); 
   
  /*
  err = genVDPoissonSampling(csmask,  opvthick, yfov , csnsz, yviews,  redfacKz, 
			     redfacKy, centKy, doellipyz, pp );
  */

  plhs[0] = mxCreateDoubleMatrix(sky,skz,mxREAL);
  mask = mxGetPr(plhs[0]);
  
  for (i=0;i < sky*skz; i++)
  {
	mask[i] = (double)csmask[i];
  }
}

/*----------------------------- End of mex interface ------------
 The rest of code is suppose to be the same as the one in the sequence poissonDisc.c
 Except :
        - DB_MSG should be remove*/



/*
int genVDPoissonSampling(int *acqmask, float fovy, float fovz, int sky, 
			 int skz, float ry, float rz, int ncal, 
			 int cutcorners, float pp); 
int initListGrid2D(VarCFloatComplexGrid *, int, int);
int appendListGrid2D(VarCFloatComplexGrid *, int, int, float complex);
void freeListGrid2D(VarCFloatComplexGrid *);
int initRandomQueue(VarCFloatComplexArray *, int);
int pushRandomQueue(VarCFloatComplexArray *, float complex);
int popRandomQueue(VarCFloatComplexArray *, float complex *, long *);
void freeRandomQueue(VarCFloatComplexArray *);
long ran_seed;
*/


/* Set up an 2D array to mimic the python grid structure.  Each cell
   has space for some moderate amount of points, and the append
   function will increase that later if needed.  A bad return only
   happens if memory can't be allocated */
int initListGrid2D(VarCFloatComplexGrid *p, int w, int h) 
{
  int i;
  int initsz = 10;
  
  p->h = h;
  p->arrSz = w*h;
  /* These need to be initialized in case freeListGrid2D is called early */
  p->nele = NULL;
  p->sz = NULL;
  p->arrR = NULL;
  p->arrI = NULL;
  /* The number of actual points in each cell */
  p->nele =(int *) malloc(p->arrSz*sizeof(*p->nele)); 
  if (!p->nele) {
    freeListGrid2D(p);
    return FAILURE;
  }
  /* The max number of points in each cell */
  p->sz =(int *) malloc(p->arrSz*sizeof(*p->sz)); 
  if (!p->sz) {
    freeListGrid2D(p);
    return FAILURE;
  }
  /* Space for the grid */
  p->arrR =(float **) malloc(p->arrSz*sizeof(*p->arrR));
  p->arrI =(float **) malloc(p->arrSz*sizeof(*p->arrI));
  if (!p->arrR || !p->arrI) {
    freeListGrid2D(p);
    return FAILURE;
  }
  /* Also in case freeListGrid2D gets called */
  for (i = 0; i < p->arrSz; i++) 
  {
      p->arrR[i] = NULL;
      p->arrI[i] = NULL;
  }
  
  for (i = 0; i < p->arrSz; i++) {
    /* Start each out with initsz spaces */
    p->arrI[i] = (float *) malloc(initsz*sizeof(**p->arrI));
    p->arrR[i] = (float *) malloc(initsz*sizeof(**p->arrR));
    
    if (!p->arrI[i] || !p->arrR[i]) {
      freeListGrid2D(p);
      return FAILURE;
    }
    p->sz[i] = initsz;
    p->nele[i] = 0;
  }
  return SUCCESS;
}

/* Put the point "pt" into the grid at spot (wc, hc).  Failure only
   occurs if memory allocation fails as bad point locations are just
   ignored */
int appendListGrid2D(VarCFloatComplexGrid *p, int wc, int hc, float complex pt) 
{
  /* Cell index */
  int idx = wc*p->h + hc;
  /* If out of range just return */
  if (idx < 0 || idx >= p->arrSz) return SUCCESS;
  /* If not enough space, resize first */
  if (p->nele[idx] + 1 > p->sz[idx]) {
    float  *tmpR, *tmpI;
    tmpR = (float  *) realloc(p->arrR[idx], 2*p->sz[idx]*sizeof(**p->arrR));
    tmpI = (float  *) realloc(p->arrI[idx], 2*p->sz[idx]*sizeof(**p->arrI));
    /* Didn't work! */
    if (!tmpR || !tmpI) return FAILURE;
    p->arrR[idx] = tmpR;
    p->arrI[idx] = tmpI;
    p->sz[idx] *= 2;
  }

  p->arrR[idx][p->nele[idx]] = creal(pt);
  p->arrI[idx][p->nele[idx]] = cimag(pt);
  p->nele[idx] += 1;

  return SUCCESS;
}

void freeListGrid2D(VarCFloatComplexGrid *p) 
{
  if (p->nele) {
    free(p->nele);
    p->nele = NULL;
  }
  if (p->sz) {
    free(p->sz);
    p->sz = NULL;
  }
  if (p->arrR) {
    int i;
    for (i = 0; i < p->arrSz; i++) {
      if (p->arrR[i]) {
	free(p->arrR[i]);
	p->arrR[i] = NULL;
      }
    }
    free(p->arrR);  
    p->arrR = NULL;
  }
  
    if (p->arrI) {
    int i;
    for (i = 0; i < p->arrSz; i++) {
      if (p->arrI[i]) {
	free(p->arrI[i]);
	p->arrI[i] = NULL;
      }
    }
    free(p->arrI);  
    p->arrI = NULL;
  }
}

/* Set up an array to mimic the python Random Queue construct */
int initRandomQueue(VarCFloatComplexArray *p, int sz) 
{
  /* Allocate the array */
  p->sz = sz;
  p->nele = 0;
  p->arrR = (float *) malloc(sz*sizeof(*p->arrR));
  p->arrI = (float *) malloc(sz*sizeof(*p->arrI));
  if (!p->arrR || !p->arrI) return FAILURE;

  return SUCCESS;
}

int pushRandomQueue(VarCFloatComplexArray *p, float complex pt) 
{
  /* If not enough space, resize first */
  if (p->nele + 1 > p->sz) {
    float  *tmpR, *tmpI;
    tmpR = (float *) realloc(p->arrR, 2*p->sz*sizeof(*p->arrR));
    tmpI = (float *) realloc(p->arrI, 2*p->sz*sizeof(*p->arrI));
    
    /* Didn't work! */
    if (!tmpR || !tmpI) {
      free(p->arrR);
      p->arrR = NULL;
      free(p->arrI);
      p->arrI = NULL;
      
      return FAILURE;
    }
    p->arrR = tmpR;
    p->arrI = tmpI;
    p->sz *= 2;
  }
  p->arrR[p->nele] = creal(pt);
  p->arrI[p->nele] = cimag(pt);
  p->nele += 1;
  
  return SUCCESS;
}
/* This only fails if there are no points in the queue */
int popRandomQueue(VarCFloatComplexArray *p, float complex *pt, long *seed) 
{
  if (p->nele == 0) return FAILURE;
  else {
    int idx = (int) floor(p->nele*ran0(seed));
    /* Just in case */
    if (idx == p->nele) idx--;
    *pt = p->arrR[idx]+p->arrI[idx]*I ;
    p->arrR[idx] = p->arrR[p->nele - 1];
    p->arrI[idx] = p->arrI[p->nele - 1];
    p->nele--;
  }
  
  return SUCCESS;
}
void freeRandomQueue(VarCFloatComplexArray *p) 
{
  if (p->arrR || p->arrI) {
    free(p->arrR);
    free(p->arrI);
    p->arrR = NULL;
    p->arrI = NULL;
  }
}

/* Already in mfast */
void epic_error(int blah, char *str, int foo, int bar) 
{
    blah=foo;
    foo=bar;
  fputs(str, stderr);
  fflush(stderr);
  return;
}

/* Numrec C random number generator, returns values between 0 and 1 */
float ran0(long *idum)
{
  long k;
  float ans;
  long mask = 123459876;
  long iq = 127773;
  long ia = 16807;
  long im = 2147483647;
  long ir = 2836;
  float am;
  
  am = 1.0 / im;
  *idum ^= mask;
  k = (*idum) / iq;
  *idum = ia*(*idum - k*iq) - ir*k;
  if (*idum < 0) *idum += im;
  ans = am*(*idum);
  *idum ^= mask;
  return ans;
}






/*
int main(int argc, char *argv[])
{
  int err;
  int i;
  int csmask[MAXACQS];
  float yfov = 256;
  float opvthick = 256;
  int yviews = 256*3;
  int csnsz = 256*3;
  float redfacKy = 2*3;
  float redfacKz = 2*3;
  int centKy = 10;
  int doellipyz = 0;
  FILE *filep;
  
  ran_seed = (long) time(0);
  err = genVDPoissonSampling(csmask, yfov, opvthick, yviews, csnsz, redfacKy, 
			     redfacKz, centKy, doellipyz);

  // The csmask table as a disp image:
  //   Number of rows (-h) in mask is yviews
  //   Number of columns (-w) in mask is nsz 
  if ((filep = fopen("vdpoisson_lg.img", "w")) == NULL) {
    printf("Can't open file!\n");
    return FAILURE;
  }
  for (i = 0; i < yviews*csnsz; i++) {
    short blah = (short) csmask[i];
    fwrite(&blah, sizeof(short), 1, filep);
  }
  fclose(filep);

  exit(0);
}  */

/* Generates a variable density resolution Poisson disc sampling
   pattern. From M. Lustig

   Output: acqmask - sampling pattern array. The array size has to be
                     at least sky*skz, and is indexed by acqmask[i*skz + j], 
		     with 0 <= i < sky, 0 <= j < skz. 
		     The row/i index reference ky values, while the
		     column/j index reference kz values

   Input: fovy - FOV in Y, in mm
          fovz - FOV in Z, in mm
	  sky - number of ky values
          skz - number of kz values, has to be > 1
	  ry - acceleration rate in Y with no calibration area
	  rz - acceleration rate in Z with no calibration area
	  ncal - size of the calibration area in ky and kz
	  cutcorners - remove the ky, kz corners
          pp - polynomial order of variable density (0 = uniform)
*/
/* Stop here, check, can I do 2d? */

int genVDPoissonSampling(int *acqmask, float fovy, float fovz, int sky, 
			 int skz, float ry, float rz, int ncal, 
			 int cutcorners, float pp) 
{

  const char* myname = "genVDPoissonSampling:";
  char tmpstr[80];
  int i, j;
  int masksum;
  int cellDim = 2;
  Genpoint = 0;
  int numneighpts = 30;
  float rsy = 1;
  float rsz = 1;
  float mr_top = 400;
  float mr_bot = 0;
  float res_y, res_z;
  float r0;
  float tol = 0.05;
  float grid_max;
  float gridCsz;
  int grid_w, grid_h;
  float pW, pH;
  VarCFloatComplexGrid grid;
  VarCFloatComplexArray proclist;

  short *mask = NULL;
  short *automask = NULL;
  float *R = NULL;
  float complex *r_grid = NULL;
    //long ran_seed = (long) time(0);
  long ran_seed = 10000;

  /* 2D check */
  if (skz == 1) rz = 1;

  /* Setup */
  r0 = (sky > skz) ? (float)ncal/(float)sky : (float)ncal/(float)skz;
  /* Very unlikely, but protect against divide by 0 errors below */
  if (r0 == 1) {
    if (sky == ncal)
      sprintf(tmpstr, "%s The number of in-plane PE's must be > %d\n",
	      myname, ncal);
    if (skz == ncal)
      sprintf(tmpstr, "%s The number of slice PE's must be > %d\n",
	      myname, ncal);
    epic_error(0, tmpstr, 0, 0);
    return FAILURE;
  }

  /* Allocate arrays */
  mask = (short*) malloc(sky*skz*sizeof(*mask));
  automask = (short*) malloc(sky*skz*sizeof(*automask));
  R = (float*) malloc(sky*skz*sizeof(*R));
  r_grid = (float complex* ) malloc(sky*skz*sizeof(*r_grid));

  /* Clear out the acquisition array */
  for (i = 0; i < sky*skz; i++) acqmask[i] = 0;

  /* Adjust the scaling for the R matrix */
  res_y = fovy/sky;
  res_z = fovz/skz;
  

  if (res_y > res_z) rsy = res_z/res_y;
  else               rsz = res_y/res_z;

  
 /* 
  for(i = int(sky/2 - ncal/2); i < int(sky/2 + ncal/2); i++)
        {
            for(j = int(skz/2 - ncal/2); j < int(skz/2 + ncal/2); j++)
            {                
                // Add the autocalibration lines and counts line.

                    automask[i*skz + j]=1;
            }
        }
  */
  
  /* Calculate the masks and the R array */
  masksum = 0;
  for (i = 0; i < sky; i++) {
    float y = -1 + 2.0*i/(sky - 1);
    for (j = 0; j < skz; j++) {
      float z = (skz == 1) ? 0 : -1 + 2.0*j/(skz - 1);
      /* Autocalibration lines */

      /*automask[i*skz + j] = 
	(sqrt(rsy*rsy*y*y) < r0 && sqrt(rsz*rsz*z*z) < r0) ? 1 : 0;
       */
      
      if(i>=floor(sky/2 - ncal/2) && i<floor(sky/2 + ncal/2) && j>=floor(skz/2 - ncal/2) && j<floor(skz/2 + ncal/2) )
      {
          automask[i*skz + j] = 1;
      }
      else
      {
          automask[i*skz + j] = 0;
      }
      
      /* Mask */
      if (cutcorners) mask[i*skz + j] = (sqrt(y*y + z*z) <= 1) ? 1 : 0;
      else            mask[i*skz + j] = 1;
      masksum += mask[i*skz + j];
      /* R */
      R[i*skz + j] = powf((rsy*rsy*y*y + rsz*rsz*z*z),pp/2.0);
    }
  }

  /* Calculate the parameters that give the desired acceleration
     numerically using bisection */
  int iter=0;
  while (1) {
    float est_accl = 0;
    float dsum = 0;
    float mr = mr_bot/2.0 + mr_top/2.0;
    float rgridR, rgridI;
    grid_max = 0;
    for (i = 0; i < sky; i++) {
      for (j = 0; j < skz; j++) {
	int idx = i*skz + j;
	float rrval = ((R[i*skz + j] - r0)*(mr - 1.0))/(1.0 - r0);
	if (rrval < 0) rrval = 0;
	rgridR = rrval + 1;
	rgridI = rrval*rz/ry + 1;
	/* Save the value for the sampling density array */
	r_grid[idx] = rgridR + I*rgridI;
	dsum += mask[i*skz + j]/(rgridR*rgridI);
	/* Keep track of the max/min grid values */
	if (rgridR > grid_max) grid_max = rgridR;
	if (rgridI > grid_max) grid_max = rgridI;
      }
    }
    float est_accl_old = est_accl;
    est_accl = 1.24*1.24*masksum/dsum;
    /* All done */
    //if (fabs(est_accl - ry*rz) < tol) break;
    if (fabs(est_accl - ry*rz) < tol || abs(est_accl_old-est_accl) < tol/10.0) break;
    
    if(iter == 200000)
    {
    //mexPrintf("Reach 200000 iteration withour converging =  %d \n",iter);
    sprintf(tmpstr, "Reach 200000 iteration withour converging\n");
    break;
    }
    
    /* Not done, so try different parameters */
    if (est_accl < ry*rz) mr_bot = mr;
    else                     mr_top = mr;
    
    iter++;
  }

  /* The "sample_poisson_ellipse" section of Miki's script */
  gridCsz = grid_max/sqrt(2.0);
  grid_w = (int) ceil(sky/gridCsz);
  grid_h = (int) ceil(skz/gridCsz);

  /* Initialize the grid */
  if (initListGrid2D(&grid, grid_w, grid_h) == FAILURE) {
    sprintf(tmpstr, "%s Can't allocate memory for ListGrid2D!\n", myname);
    epic_error(0, tmpstr, 0, 0);
    if (mask) free(mask);
    if (automask) free(automask);
    if (R) free(R);
    if (r_grid) free(r_grid);
    return FAILURE;
  }
  /* Initialize the random queue array */
  if (initRandomQueue(&proclist, sky*skz) == FAILURE) {
    sprintf(tmpstr, "%s Can't allocate memory for Random Queue!\n", myname);
    epic_error(0, tmpstr, 0, 0);
    freeListGrid2D(&grid);
    if (mask) free(mask);
    if (automask) free(automask);
    if (R) free(R);
    if (r_grid) free(r_grid);
    return FAILURE;
  }
  /* Generate the first point */
  pW = sky*ran0(&ran_seed);
  pH = (skz == 1) ? 0 : skz*ran0(&ran_seed);
  /* Put the point in the queue */
  if (pushRandomQueue(&proclist, pW + I*pH) == FAILURE) {
    sprintf(tmpstr, "%s Can't allocate memory in RQ push!\n", myname);
    epic_error(0, tmpstr, 0, 0);
    freeListGrid2D(&grid);
    freeRandomQueue(&proclist);
    if (mask) free(mask);
    if (automask) free(automask);
    if (R) free(R);
    if (r_grid) free(r_grid);
    return FAILURE;
  }
  /* Put the point in the grid */
  if (appendListGrid2D(&grid, (int) (pW/gridCsz), 
		       (int) (pH/gridCsz), pW + I*pH) == FAILURE) {
    sprintf(tmpstr, "%s ListGrid2D append failed!\n", myname);
    epic_error(0, tmpstr, 0, 0);
    freeListGrid2D(&grid);
    freeRandomQueue(&proclist);
    if (mask) free(mask);
    if (automask) free(automask);
    if (R) free(R);
    if (r_grid) free(r_grid);
    return FAILURE;
  }
  /* And use it */
  acqmask[(int) pW *skz + (int) pH] = mask[(int) pW *skz + (int) pH];
  
  /* Generate other points from points in the queue */
  do {
    float rW, rH;
    float complex p;
    /* Get a point. This fails only if there are no more points */
    if (popRandomQueue(&proclist, &p, &ran_seed) == FAILURE) break;
    pW = creal(p);
    pH = cimag(p);
    rW = creal(r_grid[(int) pW*skz + (int) pH]);
    rH = cimag(r_grid[(int) pW*skz + (int) pH]);

    /* Generate a number of points around points already in the
       sample, and then check if they are too close to other
       points. Typically numneighpts = 30 is sufficient. The larger
       numneighpts the slower the algorithm, but the more sample
       points are produced */
    for (i = 0; i < numneighpts; i++) {
      int cW, cH;
      /* Generate a point randomly selected around p, between r and
	 2*r units away.  Note, r >= 1 (or should be!) */
      float ratio = rH/rW;
      float rr = rW*(1 + ran0(&ran_seed));
      float rt = 2*3.14159265358979323846*ran0(&ran_seed);
      /* The generated point q */
      float qW = rr*sin(rt)*ratio + pW;
      float qH = (skz == 1) ? 0 : rr*cos(rt) + pH;
      int in_neighborhood = 0;
      /* If inside the rectangle continue */
      if (0 <= qW && qW < sky && 0 <= qH && qH < skz) {
	int cell, l;
	/* This is the "in_neighbourhood(q, r)" condition */
	for (cW = -cellDim; cW <= cellDim; cW++) {
	  int cellidxW = (int) (qW/gridCsz) + cW;
	  if (cellidxW < 0 || cellidxW >= grid_w) continue;
	  for (cH = -cellDim; cH <= cellDim; cH++) {
	    int cellidxH = (int) (qH/gridCsz) + cH;
	    if (cellidxH < 0 || cellidxH >= grid_h) continue;
	    /* Cell index in the grid array */
	    cell = cellidxW*grid_h + cellidxH;
	    /* Number of points in the cell */
	    for (l = 0; l < grid.nele[cell]; l++) {
	      float cptW = (grid.arrR[cell][l]);
	      float cptH = (grid.arrI[cell][l]);
	      /* To keep the point q, the following has to be false
		 for all points l in all cells */
	      if ((qW - cptW)*(qW - cptW) + 
		  (qH - cptH)*(qH - cptH)/(ratio*ratio) <= rW*rW) {
		in_neighborhood = 1;
		break;
	      }
	    }
	    if (in_neighborhood) break;
	  }
	  if (in_neighborhood) break;
	}
	/* Add the point */
	if (!in_neighborhood) {
	  /* Put the point in the queue */
	  if (pushRandomQueue(&proclist, qW + I*qH) == FAILURE) {
	    sprintf(tmpstr, "%s Can't allocate memory in RQ push!\n", myname);
	    epic_error(0, tmpstr, 0, 0);
	    freeListGrid2D(&grid);
	    freeRandomQueue(&proclist);
	    if (mask) free(mask);
	    if (automask) free(automask);
	    if (R) free(R);
	    if (r_grid) free(r_grid);
	    return FAILURE;
	  }
	  /* Put the point in the grid */
	  if (appendListGrid2D(&grid, (int) (qW/gridCsz), 
			       (int) (qH/gridCsz), qW + I*qH) == FAILURE) {
	    sprintf(tmpstr, "%s ListGrid2D append failed!\n", myname);
	    epic_error(0, tmpstr, 0, 0);
	    freeListGrid2D(&grid);
	    freeRandomQueue(&proclist);
	    if (mask) free(mask);
	    if (automask) free(automask);
	    if (R) free(R);
	    if (r_grid) free(r_grid);
	    return FAILURE;
	  }
	  acqmask[(int) qW*skz + (int) qH] = mask[(int) qW*skz + (int) qH];
	}
      }
    }
  } while (proclist.nele);

  /* Add the autocalibration lines and count all the points.  Genpoint
     isn't actually used, but it's easy to throw in here */
  Genpoint = 0;
  for (i = 0; i < sky; i++) {
    for (j = 0; j < skz; j++) {
      int idx = i*skz + j;
      if (automask[idx] && !acqmask[idx]) acqmask[idx] = 1;
      Genpoint += acqmask[idx];
    }
  }
  
  AccTotale = masksum / (float) Genpoint;
  
  freeListGrid2D(&grid);
  freeRandomQueue(&proclist);
  if (mask) free(mask);
  if (automask) free(automask);
  if (R) free(R);
  if (r_grid) free(r_grid);
  
  return SUCCESS;
  
}
/*
    printf("Processing point (%f, %f), r = (%f, %f)\n", pW, pH, rW, rH);


  printf("ListGrid:\n");
  for (i = 0; i < grid_w; i++) {
    for (j = 0; j < grid_h; j++) {
      int idx = i*grid_h + j;
      if (grid.nele[idx]) {
	int k;
	printf("Index %d has %d elements\n", idx, grid.nele[idx]);
	printf("Which are:\n");
	for (k = 0; k < grid.nele[idx]; k++)
	  printf("(%f, %f)\n", creal(grid.arr[idx][k]), 
		 cimag(grid.arr[idx][k]));
      }
    }
  }
  

  printf("masksum = %d/%d\n", masksum, sky*skz);
  printf("min, max = %f, %f\n", grid_min, grid_max);
  printf("grid cell = %f, max2 = %f, grid_w/h = %d/%d\n", 
	 gridCsz, grid_max2, grid_w, grid_h);


      printf("Got point p (%f, %f)\n", creal(p), cimag(p));


	if (i == 50 && j == 42) {
	  printf("R[%d][%d] = %f, %f, %f\n", i, j, R[i*skz + j], r0, mr);
	  printf("RR[%d][%d] = %f\n", i, j, rrval);
	  printf("D[%d][%d] = %f\n", i, j, 
		 mask[i*skz + j]/((rrval + 1)*(rrval*ratio_accl + 1)));
	}
    printf("sum mask = %d\n", masksum);
    printf("sum D = %f\n", dsum);
    printf("est_accel/tgt_accl = %f/%f\n\n", est_accl, tgt_accl);

  for (i = 0; i < sky; i++) {
    for (j = 0; j < skz; j++) {
      int idx = i*skz + j;
      if (i == 50)
	printf("rgrid[%d][%d] = (%f, %f)\n", i, j, creal(r_grid[idx]), 
	       cimag(r_grid[idx]));
    }
  }

*/
