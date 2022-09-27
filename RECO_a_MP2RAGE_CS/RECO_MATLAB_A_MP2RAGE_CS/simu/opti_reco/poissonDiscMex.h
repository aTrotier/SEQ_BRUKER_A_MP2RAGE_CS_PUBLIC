#ifndef StructPoisson
#define StructPoisson

typedef struct VarCFloatComplexArray{
    int sz;
    int nele;
    float  *arrR;
    float  *arrI;
}VarCFloatComplexArray;



typedef struct {
  int h;
  int arrSz;
  int *nele;
  int *sz;
  float  **arrR;
  float  **arrI;
}VarCFloatComplexGrid;

float ran0(long *);
void epic_error(int, char *, int, int);

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

#endif