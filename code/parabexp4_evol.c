#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..N+1 */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..N+1 */


#define M (100) //a tipusok szama
#define EXT_LIM (1E-5) // ez a kihalasi limit
#define ADDCONC (2E-5) // ennyivel adjuk hozza az ujat
#define COUNT (3E-5) // ettol a relativ koncentraciotol kezdve szamoljuk be az egyuttelok koze

#define N (2*M) // number of equations
#define RTOL  RCONST(1.0e-8)   /*scalar relative tolerance*/
#define T0    RCONST(0.0)      /*initial time*/
#define NOUT  (2E7)          /*number of output times */

#define TN (20) // a tipusok dimenzionkenti szama
#define STN (/*(TN + 1) **/ (TN + 1) * (TN + 1))


/******* RANDOM NUMBER GENERATOR **********/
#define AA 471
#define BB 1586
#define CC 6988
#define DD 9689
#define MMM 16383
#define RIMAX 2147483648.0        /* = 2^31 */
#define RandomInteger (++nd, ra[nd & MMM] = ra[(nd-AA) & MMM] ^ ra[(nd-BB) & MMM] ^ ra[(nd-CC) & MMM] ^ ra[(nd-DD) & MMM])
void seed(long seed);
static long ra[MMM+1], nd;

void seed(long seed)
{
 int  i;
 if(seed<0) { puts("SEED error."); exit(1); }
 ra[0]= (long) fmod(16807.0*(double)seed, 2147483647.0);
 for(i=1; i<=MMM; i++)
 {
  ra[i] = (long)fmod( 16807.0 * (double) ra[i-1], 2147483647.0);
 }
}

long randl(long num)      /* random number between 0 and num-1 */
{
 return(RandomInteger % num);
}

double randd(void)
{
 return((double) RandomInteger / RIMAX);
}

/********* END OF RANDOM NUMBER GENERATOR ********/

void seed(long seed);

realtype mass, a[M+1], b[M+1], c[M+1];


static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);
static int check_flag(void *flagvalue, char *funcname, int opt);
double r;
realtype tout;
N_Vector y;
void *cvode_mem;
int freetyp[STN + 1], numfreetyp, coexnum;
int lookup[M + 1], present[M + 1];

typedef struct _IND
{
  double a, b, c;
}PARAMETERS;

PARAMETERS param[STN + 1];

void shuffle(int *array, size_t n)
{
  if(n>1)
  {
    size_t i;
    for(i=0;i<n-1;i++)
    {
      size_t j=i+randl(n-i);
      int t=array[j];
      array[j]=array[i];
      array[i]=t;
    }
  }
}

double lambda(double b, double c, double r)
{
  double x;
  x = 0.5*(-1.0 * (b + c * r) + sqrt(b * b + c * c * r * r + 6.0 * b * c * r));
  return(x);
}

void initialize_set()
{
  int i, na, nb, nc, num = 1, rr;
  
  //for(na = 0; na <= TN; na++)
  //{
    for(nb = 0; nb <= TN; nb++)
    {
      for(nc = 0; nc <= TN; nc++)
      {
        param[num].a = 50.0/* + 50.0 / TN * (double)na*/;
        param[num].b = 5.0 + 5.0 / TN * (double)nb;
        param[num].c = 1.0 + 4.0 / TN * (double)nc;
        num++;
      }
    }
  //}
  for(i = 1; i <= M; i++)
  {
    lookup[i] = -999;
    present[i] = 0; 
  }
  
  rr = randl(STN) + 1;  
  na = 1;
  for(i = 1; i <= STN; i++)
  {  
    if(i == rr) continue;
    freetyp[na] = i;
    na++;
  }
  numfreetyp = STN - 1;
  
  a[1] = param[rr].a;
  b[1] = param[rr].b;
  c[1] = param[rr].c; 

  lookup[1] = rr;
  present[1] = 1;
}

void print_out(int j)
{
  int i, cou = 0;
  double parsum[M + 1], sum, Phi;
  
  Phi = 0.0;
  for(i = 1; i <= M; i++)
    Phi += r * c[i] * Ith(y, 2 * i - 1);
  
  for(i = 1; i <= M; i++)
    parsum [i] = 0.0;
  sum = 0.0;
  for(i = 1; i <= M; i++)
  {
    parsum[i] = Ith(y, 2* i - 1) + 2.0 * Ith(y, 2* i);
    sum += parsum[i];
    if(parsum[i] > COUNT)
      cou++;
  }

  printf("%lf\t%lf\t%d\t%lf\t", tout, r, cou, a[1]);
  if(j>0)
  {
    for(i = 1; i <= M; i++)
      printf("%le\t", parsum[i] / sum);
  }
  printf("\n");
  fflush(stdout);

}

void test()
{
  int i;
  /*double summ;
  
  for(i=1; i <= 5; i++)
    printf("a[%d]=%lf\tb[%d]=%lf\tc[%d]=%lf\n",i, a[i], i, b[i], i, c[i]);
  summ=0.0;
  for(i=1;i<=M;i++)
    summ+=Ith(y, 2 * i - 1) + 2.0 * Ith(y, 2 * i);
  for(i=1; i <= 5; i++)
    printf("%d freq=%lf\n",i, (Ith(y, 2 * i - 1) + 2.0 * Ith(y, 2 * i))/summ);
  for(i=1; i <= 5; i++)  
    printf("%d. lu:=%d\n",i,lookup[i]);
  printf("numfreetyp=%d\n", numfreetyp);
  for(i = 1; i <= numfreetyp; i++)
    printf("freetyp %d.:=%d\n", i, freetyp[i]);*/
  printf("coexnum=%d\n",coexnum);
  printf("present: ");  
  for(i = 1; i <= M; i++)
    printf("%d ", present[i]);
  printf("\n\n");
  fflush(stdout);
  
}

void Add_new()
{
  int i, rr, newtyp, na;
  
  coexnum = 0;
  for(i = 1; i <= M; i++)
    if((Ith(y, 2* i -1) + 2.0 * Ith(y, 2* i)) >  EXT_LIM)
      coexnum++;  
/*  
  printf("#> ADD ELOTT:\n");
  test();
*/  
  rr = randl(numfreetyp) + 1;
  newtyp = freetyp[rr];
  
  coexnum++;
    
  a[coexnum] = param[newtyp].a;
  b[coexnum] = param[newtyp].b;
  c[coexnum] = param[newtyp].c;
  
  lookup[coexnum]  = newtyp;
  present[coexnum] = 1;
  
  Ith(y, 2 * coexnum - 1) = ADDCONC;
  CVodeReInit(cvode_mem,tout,y);
  
//  printf("#Adding type %d as species %d\n", newtyp, coexnum);
  
  na = 1;
  for(i = 1; i <= numfreetyp; i++)
  {
     if(i == rr) continue;
     if(i > rr)
      freetyp[na] = freetyp[na + 1];
    na++;
  }
  freetyp[numfreetyp] = -999;
  numfreetyp--;
/*  
  printf("#> ADD UTAN:\n");
  test();
*/  
}


void Remove(int w)
{
  int i, typ;

/*  
  printf("#> REMOVE ELOTT:\n");
  test();
*/
  typ = lookup[w];
  numfreetyp++;
  freetyp[numfreetyp] = typ;
//  printf("#Removing type %d as species %d\n", typ, w);

  for(i = w; i < M; i++)
  {
    Ith(y, 2*i-1) = Ith(y, 2*(i+1)-1);
    Ith(y, 2*i) = Ith(y, 2*(i+1));
    a[i] = a[i + 1];
    b[i] = b[i + 1];
    c[i] = c[i + 1];
    lookup[i] = lookup[i + 1];
  }
  Ith(y, 2*M - 1) = 0.0;
  Ith(y, 2*M) = 0.0;
  lookup[M] = -999;
  for(i=1; i<=M; i++)
    if(present[i] == 0) break;
  present[i - 1] = 0;
  CVodeReInit(cvode_mem,tout,y);
  
/*  
  printf("#> REMOVE UTAN:\n");
  test();
*/
}

int main()
{
  int i, ind, flag, iout, coexnum;
  realtype reltol, t;
  N_Vector abstol;
  double acurrent;

  seed(40128);// 549190
  mass = 2.0;
  
  initialize_set();

  y = abstol = NULL;
  cvode_mem = NULL;

  y = N_VNew_Serial(N);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  abstol = N_VNew_Serial(N); 
  if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);

  for(i = 1; i <= M; i++)
  {
    Ith(y, 2 * i - 1) = 0.0;
    Ith(y, 2 * i) = 0.0;
  }
  Ith(y, 1) = (double)mass/10.0;
  
  for (i = 1; i <= N; i++)
    Ith(abstol, i) = RCONST(1.0e-6);

  reltol = RTOL;

  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
  flag = CVodeInit(cvode_mem, f, T0, y);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);
  flag = CVDense(cvode_mem, N);
  if (check_flag(&flag, "CVDense", 1)) return(1);

  iout = 1;
  tout = 0.1;

  double min=0.1, max=50.0;
  
  while(1)
  {
    
    r = min + (max-min) / 2.0 + (max-min) / 2.0 * sin(6.28*tout/200000.0);
    if(tout<=1.0E6) acurrent = 10.0 + 90.0 * tout/1.0E6;
    else acurrent  = 100.0 - 90.0 * (tout-1.0E6)/1.0E6;
    
    for(i=1; i <= M; i++)
      a[i]=acurrent;
    
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    ind = 0;
    for(i = 1; i <= N; i++) // negativok nullazasa
    {
      if(Ith(y,i) < 0.0)
      {
        Ith(y,i) = 0.0;
        ind = 1;
      }
    }
    if(ind != 0)
      CVodeReInit(cvode_mem,tout,y);
    
    for(i = 1; i <= M; i++) // kihalas
      if((present[i] == 1) && ((Ith(y, 2* i - 1) + 2.0 * Ith(y, 2* i)) <=  EXT_LIM))
        Remove(i);
    
    if(iout%500 == 0)
      print_out(0);
    
    if((iout%1000)==0)
    {
      coexnum = 0;
      for(i = 1; i <= M; i++)
        if((Ith(y, 2* i -1) + 2.0 * Ith(y, 2* i)) >  EXT_LIM)
          coexnum++;
      if(coexnum == M)
      {
        printf("System full!\n");
        exit(1);
      }      
      Add_new();
    }
    
    if (check_flag(&flag, "CVode", 1)) break;
    if (flag == CV_SUCCESS)
    {
      iout++;
      tout += 0.1;
    }
  
    if (iout > NOUT) break;
  }// itt az integralas vege


  N_VDestroy_Serial(y);
  N_VDestroy_Serial(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  return(0);
}

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  int i;
    
  realtype xx[M+1];
  realtype yy[M+1];
  realtype Phi = 0.0;

  for(i = 1; i <= M; i++)
  {
    xx[i] = Ith(y, 2 * i - 1);
    yy[i] = Ith(y, 2 * i);
  }
  
  for(i = 1; i <= M; i++)
    Phi += r * c[i] * xx[i];
  
  for(i = 1; i <= M; i++)
  {  
    Ith(ydot, 2 * i - 1) = -2.0 * a[i] * pow(xx[i], 2.0) + 2.0 * b[i] * yy[i] - r * c[i] * xx[i] - xx[i] * Phi / mass;
    Ith(ydot, 2 * i) = a[i] * pow(xx[i], 2.0) - b[i] * yy[i] + r * c[i] * xx[i] - yy[i] * Phi / mass; 
  }
  
  return(0);
}

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
