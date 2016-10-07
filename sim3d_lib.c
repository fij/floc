// ==================================== model description =========================================
//
// ----- geometry (including: initialization of coordinates) -----
//
// . At simulation time t = 0: 
//   - N point-like particles are placed randomly inside a cube with size L x L x L
//     such that the minimum distance between any two particles is at least  0.6 * L * N ^ { - 1 / 3 }
// . The cube has periodic boundary conditions in all three directions
// . Book-keeping of the particles: The simulated volume is divided into rectangular book-keeping grid cells (little cubes) of equal size:
//   - the minimum number of grid cells is used such that the side length of each grid cell (cube) remains above R
// . For each particle the coordinates pfx,pfy,pfz store the three integer coordinates of its grid cell in the grid of book-keeping cells
// . For each book-keeping grid cell 
//   - the data structures BP and PN list the particles in that grid cell
//   - the three coordinates of the grid cell (pfx,pfy,pfz) are mapped to a single grid cell index by BXYZ2I
//   - the three coordinates of the grid cell (pfx,pfy,pfz) are mapped to the number of particles in the grid cell by BXYZ2N
//
// ----- rules of motion ------
//
// . Each particle has unit mass, thus, the sum of forces acting on the particle equals the acceleration of the particle.
// . There are three forces acting on a particle:
//   (a) self-propelling
//       - the particle adjusts, linearly with time constant Tau, the magnitude of its velocity to a preferred magnitude, V0
//       - note that the particle does not change the direction of its own velocity
//   (b) radial repulsion between any two particles
//       - if the distance of the two (point-like) particles is r, then the magnitude of the repulsion is   ( 1 / r^2 ) - ( 1 / R^2 )
//       - this interaction has a cutoff at the interaction radius R
//   (c) noise 
//       - the magnitude of the noise vector depends on the adaptive integration time step, dt
//       - it is added at each update to each particle's velocity as a vector with random direction and length S * sqrt(dt)
//
// -------- numerical integration: midpoint method with fixed time step ---------
//
// (1) With Euler's method compute positions and velocities at the midpoint
// (2) With the midpoint method -- using the initial and midpoint values -- compute the final positions and velocities
//
// -------- initialization of velocity vectors --------
//
// At simulation time = 0 :
//
// . the speed of each particle is the preferred speed, V_0
// . IF _START_STATE = 0, then the direction of motion of each particle is a _different_ random unit vector selected from the full 4 PI solid angle
//   ELSE the direction of motion of each particle is the _same_ random unit vector
//
// -------- output --------
//
// IF startState = 0
// THEN: saving data at simulation time values of  t_i = 10 ** ( i / n ) for i = 0, 1, 2, ...
// ELSE:
// - saving at these times,
// - saving at simulation time        t_warmUp = N * L / V0
// - saving at times            t_i + t_warmUp
//
// ===============================================================================================================
//
// -------- implementation of time lag --------
//
// The sum of forces acting on a particle at time t  are those that were computed at time  t - tau
// Note that the sum of forces includes all three effects affecting the movement of a particle:
// - self-propelling
// - interactions
// - noise
//
// The program sets the time lag (tau) as a multiple of the integration time step (dt):  tau = MTL * dt  (MTL: Multiplier for Time Lag)
//
// Because of the time lag, the program stores the values of all variables from the current time point and MTL previous time points
// The variable IC (Index Circular) is the index of the current time step in the circular arrays storing these data
//
// Initialization of the past values of all variables: all (including midpoint values) are set to the current (t=0) values
//
// From the simulation time step t=0 to the simulation time step t=MTL the values of the variables are saved at position IC=t
// At t=MTL+1 the counter jumps to the beginning of the circular array and it becomes IC=0. 
//
// In summary:
// - for the simulation time step t, values of the variables are saved at position IC = t modulo (MTL+1) in the circular data arrays
// - data saved tau simulation time ago are at the index IC_tau = (t-tau) modulo (MTL+1) = ( IC - MTL ) modulo (MTL+1)
// - the next position where data will be saved is IC_next = ( IC + 1 ) modulo (MTL+1) = IC_tau
//
// ===============================================================================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifdef XWIN
  #include "Xlibext-20120911.c"
  #define XFillCircle(disp, win, gc, x, y, rad) XFillArc(disp, win, gc, x-rad, y-rad, 2*rad, 2*rad, 0, 64*360)
  #define XDrawCircle(disp, win, gc, x, y, rad) XDrawArc(disp, win, gc, x-rad, y-rad, 2*rad, 2*rad, 0, 64*360)
#endif

// === definitions ===

#define SQR(a)   ( (a) * (a) )
#define MAX(a,b) ( (a) > (b) ? (a) : (b) )
#define MIN(a,b) ( (a) < (b) ? (a) : (b) )
#define ABS(a)   ( (a) > 0.0 ? (a) : ((-1.0)*a) )
#define PI 3.14159265
#define MY_STRLEN 1000 // default string length
#define CHAIN_END -1 // this number tells that the chain of particles listed for the current field of the board ends at the current item
typedef unsigned long int uli;
typedef unsigned short ush;

// === parameters ===
int _RND_SEED; // random seed number
int _N; // number of particles
int _START_STATE; // 0: disordered, all others: ordered (see "initialization of velocity vectors" above) 
double _V0; // preferred magnitude of the velocity of every particle
double _TAU; // time constant for adjusting the magnitude of the velocity to the preferred magnitude
double _S; // magnitude of noise vector
double _DT; // time step of a full simulation update
int _MTL; // Multiplier for Time Lag = length of time lag (tau) / length of simulation update step (dt)
int _IC; // Index of Current data in the arrays saving current and past data (see explanation above at "implementation of time lag"
int _ICN; // Index Circular Next: the index of the next time step, which is the same as the index of the time step tau time ago
// ---- all variables are saved for the current time step and MTL past time steps ----
double *  _X, *  _Y, *  _Z, *  _X_M, *  _Y_M, *  _Z_M; // coordinates of the particle at the current time point and at the midpoint
double * _VX, * _VY, * _VZ, * _VX_M, * _VY_M, * _VZ_M, * _V, * _V_M; // velocity components and speed at the current time point and the midpoint
// EX,EY,EZ are the (x,y,z) coordinates (projections) of the unit vector pointing in the direction of motion (_M suffix: at the midpoint)
double * _EX, * _EY, * _EZ, * _EX_M, * _EY_M, * _EZ_M; 
double ** _FX_SUM, ** _FY_SUM, ** _FZ_SUM, ** _FX_SUM_M, ** _FY_SUM_M, ** _FZ_SUM_M; // components of the sum of forces acting on a particle now and at the midpoint

// ---- X11 display ----
// dimensions (in pixels) when visualizing with the X server (X11): info field width and height, X11 margin
int _X11InfoFieldWidth = 280; int _X11InfoFieldHeight = 180; int _X11Margin = 20;
int _X11DisplayAreaSize = 500; // horizontal and vertical size of the area displaying the particles
double _X11_PIC_MAGN, _X11_OBJ_MAGN; // magnification for drawing entire X11 image, additional magnification of a single object on an X11 image
int _X11GraphFieldHeight = 80; int _X11LineHeight = 25; char * _X11FontName = "10x20";
int _X11GraphLenNow; // length of X11 graph now, i.e., number of points in the graph now
// parameters set in the intializing function
int _X11WinWidth; int _X11GraphFieldWidth; int _X11GraphFieldUpEnd; int _X11WinHeight;
double * _X11Graph; // the graph displaying data in the X11 window
int _NPC = 20; // NPC: number of particle colors
int * _ColorCode; // color codes, from 0 to NPC-1: particle colors, NPC: background color, NPC+1: info color, NPC+2: particle boundary color
int _X11_DrawMethod; // which method should be used for drawing the particles

// ----- wall clock time, simulation time, visualization and saving time frequencies -----
uli _WT_NOW; // current wall clock time (WT: wall clock time), number of seconds elapsed since the epoch (Jan 1, 1970)
uli _WT_MAX; // stop simulation when this number of wall clock seconds has been used
uli _WT_START; // wall clock time (physical time) when the program is started: seconds elapased since the epoch (Jan 1, 1970)
double _DRAW_FREQ_ST; // draw after this amount of simulation time (ST: Simulation Time)
int _SAVE_FREQ_N; // saveFreqN = n: determines the simulation time values at which output should be written (see above at "output")
int _DRAW_SLEEP_WT; // number of seconds to sleep after drawing (WT: Wall clock Time)
double _ST_NOW, _ST_PREV; // current and previous simulation time
double _ST_MAX; // maximum simulation time for which the the program should be run
double _WARM_UP_TIME; // IF the simulation is started from the ordered state, THEN saving starts from this time value (see above at "output")
int _IS_FIRST_UPDATE; // 1: we are at the first update step, 0: we are not

// ----- geometry and book-keeping -----
double _L; // the simulated field is a cube with side length L
double _R; // interaction radius: two particles interact only if their distance is below R
           // the simulated field is divided into BN x BN x BN book-keeping fields of equal size, the size (BL) of each cell is the smallest number >= R
double _BL; // size of one grid cell, the letter "B" indicates book-keeping, this is the smallest real-valued number that is 
            // (i) greater or equal than R and (ii) multiplied by an integer gives the size of the simulated field, L
int _BN; // the simulated field is divided into  BN x BN  book-keeping fields, "B" as book-keeping
int *** _BXYZ2I; // BXYZ2I[ix][iy][iz] maps the three coordinates of a grid cell to a single grid cell index
// ---- all variables are saved for the current time step and MTL past time steps ----
int * _PFX, * _PFY, * _PFZ, * _PFX_M, * _PFY_M, * _PFZ_M; // the indexes of the book-keeping field of the i. particle now and at the midpoint
int *** _BP, *** _BP_M; // BP[ix][iy][iz] is the index of the 1st particle in the chain of particles in the current (ix,iy,iz) field of the board
                        // the CHAIN_END value tells: there are no particles in (ix,iy,iz), _M suffix: same thing at the midpoint of the simulation update
int * _PN, * _PN_M; // PN[i] (Particle Next): is the index of the particle following the <i>th particle in the chain of particles in the same field
                    // the CHAIN_END value tells: no more particles in this field, _M: at the midpoint of the simulation update
int *** _BXYZ2N, *** _BXYZ2N_M; // BXYZ2N[ix][iy][iz]: number of particles in the grid cell (ix,iy,iz) now and at the midpoint of the simulation update

// === function definitions ===

void CoordsToGridCellIndexes( double x, double y, double z, int gridCellNum, double gridCellNum_over_fieldSize, int * ix, int * iy, int * iz )
{
  // particle -> field index book-keeping; ix and iy are the indexes of the book-keeping field where the coordinate pair (x,y) belongs
  // fieldSize: size of the simulated file (in both directions)
  // gridCellNum: number of grid cells (in both directions)
  //
  // IF the particle falls beyond the first/last valid grid cell (due to rounding maybe), THEN put it in the first/last valid grid cell
  * ix = (int) floor( x * gridCellNum_over_fieldSize );
  * iy = (int) floor( y * gridCellNum_over_fieldSize );
  * iz = (int) floor( z * gridCellNum_over_fieldSize );
  //
  // rounding errors may put the coordinates outside the valid range: correct them, if necessary
  if( gridCellNum <= *ix ){ -- *ix; } // if( gridCellNum <= *ix ){ -- *ix; }else if( 0 > *ix ){ ++ *ix; }  
  if( gridCellNum <= *iy ){ -- *iy; } // if( gridCellNum <= *iy ){ -- *iy; }else if( 0 > *iy ){ ++ *iy; }
  if( gridCellNum <= *iz ){ -- *iz; }
}

// ---------------------------------------

void RemoveParticleFromBookKeepingTable( int i, int pfx_i, int pfy_i, int pfz_i, int **** bp, int ** pn, int **** bxyz2n )
{
  // removing the i. particle from the "grid cell -> particle index list" book-keeping

  // IF the i. particle is the first one in the chain of particles saved for its grid cell,
  // THEN make the table's grid cell point to the chain value of the i. particle:
  //   - If the i. particle is the last one in the chain of particles saved for its grid cell,
  //     then this chain value will be CHAIN_END and the table will receive this CHAIN_END value.
  //     This will mean that after the removal of the i. particle as the only from its cell, its cell will contain no particles.
  //   - If the i. particle is followed by at least one other particle in the list of
  //     particles in its grid cell, then after removing the i. particle from the cell,
  //     the first particle in the list will be the particle that was following it.
  if( i == (*bp)[ pfx_i ][ pfy_i ][ pfz_i ] ){ (*bp)[ pfx_i ][ pfy_i ][ pfz_i ] = (*pn)[ i ]; }
  else{
      // IF the i. particle is not the first in the list of particles listed for its grid cell,
      // THEN first find the particle that is before it in the list: This will be the j. particle.
      int j = (*bp)[ pfx_i ][ pfy_i ][ pfz_i ]; while( i != (*pn)[ j ] ){ j = (*pn)[ j ];
 }
      // Now the particle before the j. particle should not point to the i. particle,
      // but to the particle (or CHAIN_END flag) that the i. particle points to.
      (*pn)[ j ] = (*pn)[ i ];
  }

  // decrease by 1 the counter counting the number of particles in the current book-keeping grid cell of the particle
  --(*bxyz2n)[ pfx_i ][ pfy_i ][ pfz_i ];
}

// ---------------------------------------

void InsertParticleIntoBookKeepingTable( int i, int pfx_i, int pfy_i, int pfz_i, int **** bp, int ** pn, int **** bxyz2n )
{
  // field index -> list of the indexes of particles in this field
  //
  // IF the i. particle is the first in its book-keeping field,
  // THEN add the i. particle to the list of particles in this field and close the chain of particles in its book-keeping field
  if( CHAIN_END == (*bp)[pfx_i][pfy_i][pfz_i] ){ (*bp)[pfx_i][pfy_i][pfz_i] = i; (*pn)[i] = CHAIN_END; }
  else{
    // IF the current particle is not the first in its book-keeping field, 
    // THEN 
    // (1) find the last particle in its book-keeping field,
    int j = (*bp)[pfx_i][pfy_i][pfz_i]; while( CHAIN_END != (*pn)[j] ){ j = (*pn)[j]; }
    // (2) append the current particle to the chain after the last particle and
    (*pn)[ j ] = i;
    // (3) the i. particle has now become the last in the list of particles saved for its grid cells, thus, the chain ends after it
    (*pn)[ i ] = CHAIN_END;
  }

  // increment the counter that counts the number of particles in the current grid cell
  ++ (*bxyz2n)[ pfx_i ][ pfy_i ][ pfz_i ];
}

// ---------------------------------------

// generate a 3d unit (unit length) vector pointing in a direction distributed evenly in the full 4 PI solid angle
void Generate3dVector_randomDir_unitLength( double * coord_x, double * coord_y, double * coord_z )
{
  // (1) generate random coordinates in the (+/-1,+/-1,+/-1) cube
  // (2) keep vector only if its length is between 0.1 and 1
  // (3) normalize the vector to unit length
  double sumSqr; do{ 
      * coord_x = 1.0 - 2.0 * rand() / ( 1.0 + RAND_MAX );
      * coord_y = 1.0 - 2.0 * rand() / ( 1.0 + RAND_MAX );
      * coord_z = 1.0 - 2.0 * rand() / ( 1.0 + RAND_MAX );
      sumSqr = SQR( * coord_x ) + SQR( * coord_y ) + SQR( * coord_z );
  }while( sumSqr < .01 || sumSqr >= 1 );

  // 1.0 / the length of the vector
  double one_over_length_of_vector = pow( sumSqr, -0.5 );

  // normalize the vector to unit length
  * coord_x *= one_over_length_of_vector;
  * coord_y *= one_over_length_of_vector;
  * coord_z *= one_over_length_of_vector;
}

// ---------------------------------------

// allocate memory to a 2-dimensional tensor of integers (a matrix) with sizes nA and nB
void Alloc_int2tensor( int *** t, int nA, int nB )
{
  * t = (int**)calloc( nA, sizeof(int*) );
  int iA; for( iA=0; iA<nA; ++iA ){
      ( * t )[iA] = (int*)calloc( nB, sizeof(int) );
  }
}

// ---------------------------------------

// allocate memory to a 2-dimensional tensor of doubles (a matrix) with sizes nA and nB
void Alloc_d2tensor( double *** t, int nA, int nB )
{
  * t = (double**)calloc( nA, sizeof(double*) );
  int iA; for( iA=0; iA<nA; ++iA ){
      ( * t )[iA] = (double*)calloc( nB, sizeof(double) );
  }
}

// ---------------------------------------

// allocate memory to a 2-dimensional tensor of integers (a matrix) with sizes nA and nB
// fill up the tensor with the value 'val'
void AllocInit_int2tensor( int *** t, int nA, int nB, int val )
{
  * t = (int**)calloc( nA, sizeof(int*) );
  int iA; for( iA=0; iA<nA; ++iA ){
      ( * t )[iA] = (int*)calloc( nB, sizeof(int) );
      int iB; for( iB=0; iB<nB; ++iB ){
	  ( * t )[iA][iB] = val;
      }
  }
}

// ---------------------------------------

// allocate memory to a 3-dimensional tensor of integers with sizes nA, nB and nC
// fill up the tensor with the value 'val'
void AllocInit_int3tensor( int **** t, int nA, int nB, int nC, int val )
{
  * t = (int***)calloc( nA, sizeof(int**) );
  int iA; for( iA=0; iA<nA; ++iA ){
      ( * t )[iA] = (int**)calloc( nB, sizeof(int*) );
      int iB; for( iB=0; iB<nB; ++iB ){
	  ( * t )[iA][iB] = (int*)calloc( nC, sizeof(int) );
	  int iC; for( iC=0; iC<nC; ++iC ){
	      ( * t )[iA][iB][iC] = val;
	  }
      }
  }
}

// ---------------------------------------

int nextNumberWithModulo( int number, int modulo_base ){ return ( number + 1 ) % modulo_base; }

// ---------------------------------------

void Init( int argc, char * argv[], int * rndSeed, int * n, int * startState, double * v0,
	   double * warmUpTime, double * tau, double * s, double * dt, int * mtl, int * ic, int * icn, int * is_first_update,
	   double * l, double * r, double * bl, int * bn, int **** bxyz2i, 
	   int ** pfx, int ** pfy, int ** pfz, int ** pfx_m, int ** pfy_m, int ** pfz_m,
	   int **** bp, int **** bp_m, int ** pn, int ** pn_m, int **** bxyz2n, int **** bxyz2n_m,
	   double * stMax, double * stNow, double * stPrev, double * drawFreqST, int * saveFreqN, int * drawSleepWT,
	   uli * wtMax, uli * wtStart, uli *wtNow,
	   double ** x,    double ** y,    double ** z,    double  ** vx,      double  ** vy,      double  ** vz,      double ** v,
	   double ** x_m,  double ** y_m,  double ** z_m,  double  ** vx_m,    double  ** vy_m,    double  ** vz_m,    double ** v_m,
	   double ** ex,   double ** ey,   double ** ez,   double *** fxSum,   double *** fySum,   double *** fzSum, 
	   double ** ex_m, double ** ey_m, double ** ez_m, double *** fxSum_m, double *** fySum_m, double *** fzSum_m,
	   int * x11drawMethod, int x11displayAreaSize, int * x11winWidth, int x11infoFieldWidth, double * x11picMagn, double * x11objMagn, int x11margin,
	   int * x11graphFieldWidth, int * x11graphFieldUpEnd, int x11infoFieldHeight, int * x11winHeight, int x11graphFieldHeight,
	   int npc, int ** colorCode, double ** x11graph, char * x11fontName, int * x11graphLenNow )
{
  // ------------- initialize parameters -----------------
  // IF the number of arguments is incorrect, THEN write to stderr how to use the program
  if( 18 != argc ){

    fprintf(stderr,"\n\tUsage: %s \\\n",argv[0]);
    fprintf(stderr,  "\t       \\\n");

    fprintf(stderr,  "\t       <seed number of the random number generator> \\\n");
    fprintf(stderr,  "\t       <particle number> \\\n");
    fprintf(stderr,  "\t       <start state (0: disordered, 1: ordered) -- see header of sim3d_lib.c for details> \\\n");
    fprintf(stderr,  "\t       <side length of the simulated cubic volume> \\\n");
    fprintf(stderr,  "\t       <interaction radius, two particles interact only if their distance is below this value> \\\n");
    fprintf(stderr,  "\t       \\\n");

    fprintf(stderr,  "\t       <v0: preferred velocity (length of velocity vector) for each particle> \\\n");
    fprintf(stderr,  "\t       <tau: time constant for adjusting the length of the velocity vector to the preferred value, v0> \\\n");
    fprintf(stderr,  "\t       <noise amplitude (set to zero before \'warmUpTime\' (see sim3d.c)> \\\n");
    fprintf(stderr,  "\t       <time step of one full simulation update (with the midpoint method)> \\\n");
    fprintf(stderr,  "\t       <an _integer_ Multiplier for the Time Lag (mtl) = length of time lag (tau) / length of simulation update step (dt)> \\\n");
    fprintf(stderr,  "\t       \\\n");

    fprintf(stderr,  "\t       <drawing method of particles in the X11 window (0: arrow, 1: disc)> \\\n");
    fprintf(stderr,  "\t       <magnification ratio for the size of an object within the X11 image, a real number> \\\n");
    fprintf(stderr,  "\t       <maximum simulation time (simulation stops when reaching this simulation time), negative number means: no limit> \\\n");
    fprintf(stderr,  "\t       <draw at simulation time intervals of this length> \\\n");
    fprintf(stderr,  "\t       <n (integer): save at simulation time values of t = 10 ** ( i / n ) for i = 0, 1, 2, ... > \\\n");
    fprintf(stderr,  "\t       <if showing visual output: sleep (i.e., do nothing) this number of wall clock seconds (real time) after drawing> \\\n");
    fprintf(stderr,  "\t       <stop simulation when this number of wall clock seconds has been used, negative number means: no limit>\n");
    fprintf(stderr,  "\n");

    exit(0);
  }
  // IF the number of parameters is correct, THEN save the parameters
  int iArg=0;
  //
  * rndSeed    = atoi(argv[++iArg]);
  * n          = atoi(argv[++iArg]);
  * startState = atoi(argv[++iArg]);
  * l          = atof(argv[++iArg]);
  * r          = atof(argv[++iArg]);
  //
  * v0     = atof(argv[++iArg]);
  * tau    = atof(argv[++iArg]);
  * s      = atof(argv[++iArg]);
  * dt     = atof(argv[++iArg]);
  * mtl    = atoi(argv[++iArg]);
  //
  * x11drawMethod = atof(argv[++iArg]);
  * x11objMagn    = atof(argv[++iArg]);
  * stMax         = atof(argv[++iArg]);
  * drawFreqST    = atof(argv[++iArg]);
  * saveFreqN     = atoi(argv[++iArg]);
  * drawSleepWT   = atoi(argv[++iArg]);
  * wtMax    = (uli)atoi(argv[++iArg]);
  //
  // print start message, set start and current wall clock time (wt), set random seed, set simulation time and previous simulation time
  fprintf(stderr,"Initializing ... "); fflush(stderr); * wtStart = (uli)time(NULL); * wtNow = * wtStart; srand( * rndSeed ); * stNow = 0.0; * stPrev = 0.0;


  // ------------ initialize the coordinate-based book-keeping data structures of the particles -----------------
  // set size of the cells of the book-keeping grid and the size (both height and width) of one grid cell
  * bn = (int)floor( * l / * r ); * bl = * l / * bn;
  // allocate memory to the changing book-keeping arrays and fill them up with the values provided in the last argument
  AllocInit_int3tensor( bp,   *bn, *bn, *bn, CHAIN_END ); AllocInit_int3tensor( bxyz2n,   *bn, *bn, *bn, 0 );
  AllocInit_int3tensor( bp_m, *bn, *bn, *bn, CHAIN_END ); AllocInit_int3tensor( bxyz2n_m, *bn, *bn, *bn, 0 );
  // allocate memory to the constant table converting between the two types of grid cell indexing
  int ix; for(ix=0;ix<*bn;++ix){
      * bxyz2i = (int***)calloc( * bn, sizeof(int**) ); 
      int iy; for(iy=0;iy<*bn;++iy){
	  ( * bxyz2i )[ix] = (int**)calloc( * bn, sizeof(int*) );
	  int iz; for(iz=0;iz<*bn;++iz){
	      ( * bxyz2i )[ix][iy] = (int*)calloc( * bn, sizeof(int) );
	      // mapping the three coordinates of a grid cell to a single grid cell index
	      ( * bxyz2i )[ix][iy][iz] = ( ( ix * (*bn) ) + iy ) * (*bn) + iz;
	  }
      }
  }


  // allocate memory to the book-keeping showing for each particle the next particle in its book-keeping grid cell, _m: at the midpoint of the sim.update
  *pn = (int*)calloc(*n,sizeof(int)); *pn_m = (int*)calloc(*n,sizeof(int));
  // the initial values: all particles are chain end particles of their book-keeping chains
  int i; for(i=0; i<*n; ++i){ (*pn)[i] = CHAIN_END; (*pn_m)[i] = CHAIN_END; }
  // allocate memory to the arrays storing the book-keeping grid cell index by particle index
  // these book-keepings (pfx,pfy,pfz) will be initialized when the particles are placed
  *pfx   = (int*)calloc(*n,sizeof(int)); *pfy   = (int*)calloc(*n,sizeof(int)); *pfz   = (int*)calloc(*n,sizeof(int));
  *pfx_m = (int*)calloc(*n,sizeof(int)); *pfy_m = (int*)calloc(*n,sizeof(int)); *pfz_m = (int*)calloc(*n,sizeof(int)); 


  // ----------- coordinates, velocities ---------------
  // allocate memory to particle coordinates, components of the velocity, direction of motion
  *x    = (double*)calloc(*n,sizeof(double)); *y    = (double*)calloc(*n,sizeof(double)); *z    = (double*)calloc(*n,sizeof(double)); 
  *x_m  = (double*)calloc(*n,sizeof(double)); *y_m  = (double*)calloc(*n,sizeof(double)); *z_m  = (double*)calloc(*n,sizeof(double)); 
  *vx   = (double*)calloc(*n,sizeof(double)); *vy   = (double*)calloc(*n,sizeof(double)); *vz   = (double*)calloc(*n,sizeof(double));
  *vx_m = (double*)calloc(*n,sizeof(double)); *vy_m = (double*)calloc(*n,sizeof(double)); *vz_m = (double*)calloc(*n,sizeof(double)); 
  *ex   = (double*)calloc(*n,sizeof(double)); *ey   = (double*)calloc(*n,sizeof(double)); *ez   = (double*)calloc(*n,sizeof(double));
  *ex_m = (double*)calloc(*n,sizeof(double)); *ey_m = (double*)calloc(*n,sizeof(double)); *ez_m = (double*)calloc(*n,sizeof(double)); 
  *v    = (double*)calloc(*n,sizeof(double));
  *v_m  = (double*)calloc(*n,sizeof(double));


  // -------- allocate memory to the current and past values of force component sums --------
  Alloc_d2tensor( fxSum,   1 + *mtl, *n ); Alloc_d2tensor( fySum,   1 + *mtl, *n ); Alloc_d2tensor( fzSum,   1 + *mtl, *n );
  Alloc_d2tensor( fxSum_m, 1 + *mtl, *n ); Alloc_d2tensor( fySum_m, 1 + *mtl, *n ); Alloc_d2tensor( fzSum_m, 1 + *mtl, *n );


  // -------- indexes for handling the time lag --------
  // Index of Current data in the arrays saving current and past data (see explanation above at "implementation of time lag"
  * ic = 0;
  // Index Circular Next, index of the next time step in the circular arrays,
  //                      which is the same as the index of the time step that was tau time ago
  * icn = nextNumberWithModulo( * ic, * mtl + 1 );
  // set that we are before the first completed simulation update
  * is_first_update = 1;


  // sequentially place the particles at random positions: the distance between any two particles has to be above the minimum value tested below at the (TT) sign
  /*int i;*/ for( i=0; i < *n; ++i ){
      // (1) find a random point not too close to other already placed points
      int coords_ok; do{ coords_ok = 1;
	  //
	  // random coordinates inside the L x L square
	  (*x)[i] = *l * ( rand() / ( 1.0 + RAND_MAX ) );
	  (*y)[i] = *l * ( rand() / ( 1.0 + RAND_MAX ) );
	  (*z)[i] = *l * ( rand() / ( 1.0 + RAND_MAX ) );
	  // set the (pfx[i],pfy[i]) indexes of the grid cell where the (x,y) coordinates of the i. particle belong
	  CoordsToGridCellIndexes( (*x)[i], (*y)[i], (*z)[i], *bn, ( 1.0 * *bn ) / ( *l ), *pfx + i, *pfy + i, *pfz + i );
	  //
	  // - only particles in the current and the neighboring 26 grid cells can be closer than R to the current particle 
	  // - so loop through the list of particles in the current book-keeping grid cell and the neighboring 26 book-keeping grid cells
	  // NOTE: for listing neighboring grid cells consider the periodic boundary conditions
	  // - (dx,dy,dz): relative position of the investigated grid cell from the current grid cell (dx, dy, and dz can be independently -1, 0, +1)
	  // - also: compute the min. squared distance of already placed particles (in the 27 grid cells) and the attempted coordinates of the i. particle
	  //         default value of minimum squared distance: above the square of the side length of the simulated volume (a cube)
	  double minDistSqr = 2 * SQR(*l); 
	  int dx, dy, dz; for( dx = -1; coords_ok && dx <= 1; ++dx ){ for( dy = -1; coords_ok && dy <= 1; ++dy ){ for( dz = -1; coords_ok && dz <= 1; ++dz ){
	      //
	      // (jx,jy,jz): the position of the investigated other grid cell that is at a relative position (dx,dy,dz) from the current grid cell
	      // xo: due to the periodic boundaries in the book-keeping the relative position of the neighboring cell 
	      //     is ( dx + xo * BN , dy + yo * BN ) from the current cell  ( xo can be = -1, 0, +1 )
	      // yo,zo: similar to xo
	      int jx = (*pfx)[i] + dx; int xo = 0; if( jx < 0 ){ jx += * bn; --xo; }else if( jx > * bn - 1 ){ jx -= * bn; ++xo; }
	      int jy = (*pfy)[i] + dy; int yo = 0; if( jy < 0 ){ jy += * bn; --yo; }else if( jy > * bn - 1 ){ jy -= * bn; ++yo; }
	      int jz = (*pfz)[i] + dz; int zo = 0; if( jz < 0 ){ jz += * bn; --zo; }else if( jz > * bn - 1 ){ jz -= * bn; ++zo; }
	      //
	      // . loop through the list of particles already placed into the (jx,jy,jz) book-keeping grid cell
	      //   in the (jx,jy,jz) grid cell the index of the first already placed particle is "j"
	      // . NOTE: if we see the CHAIN_END marker right at the start, THEN there are no particles in the (jx,jy,jz) grid cell
	      int j = (*bp)[jx][jy][jz]; while( CHAIN_END != j ){
		  // IF the square of the distance between the i. and j. particles is below the minimum saved so far, THEN set the minimum to this new value
		  minDistSqr = MIN( minDistSqr,   SQR( (*x)[j] + xo * *l - (*x)[i] ) 
				                + SQR( (*y)[j] + yo * *l - (*y)[i] )
				                + SQR( (*z)[j] + zo * *l - (*z)[i] )
				  );
		  // step to next particle or the "CHAIN_END" marker of the current book-keeping field 
		  j = (*pn)[j];
	      }
	  }}}
	  // (TT) IF the minimum distance is below  0.6 * L * N ^ { -1/3 },
	  //      THEN it is too low and the currently found coordinates are not OK
	  if( sqrt( minDistSqr ) < 0.6 * (*l) * pow((*n),(-1.0)/3.0) ){ coords_ok = 0; }
	  // IF the currently selected coordinates are not OK, THEN reject these coordinates and try to find new coordinates
      }while( ! coords_ok );
      //
      // IF we have found good coordinates for the i. particle,
      // THEN insert into the "field index -> particle list" book-keeping (called the book-keeping "table")
      InsertParticleIntoBookKeepingTable( i, (*pfx)[i], (*pfy)[i], (*pfz)[i], bp, pn, bxyz2n );
  }

  // initializing each particle's velocity vector:
  // - the speed (i.e., the magnitude of the velocity) of each particle is v0
  // - the direction of the velocity:
  //   IF startState = 0, THEN the direction of motion of each particle is a _different_ random unit vector selected from the full 4 PI solid angle
  switch ( * startState ) {
  case 0:
      // for each particle: set _different_ random direction for the velocity
      // generate a 3d unit (i.e., length=1) vector pointing in a direction distributed evenly in the 4 PI solid angle
      for( i=0; i < *n; ++i ){ Generate3dVector_randomDir_unitLength( *ex + i, *ey + i, *ez + i ); }
      break;
  // all other cases of startState:
  default:
      // generate a 3d unit (i.e., length=1) vector pointing in a direction distributed evenly in the 4 PI solid angle
      // save this unit vector as the direction of motion of the first particle (i.e., the particle with the index zero)
      Generate3dVector_randomDir_unitLength( *ex + 0, *ey + 0, *ez + 0 );
      // for each particle from 1 to n-1: set _same_ random direction for the velocity: copy from the first particle, i.e., the particle with index 0
      /*int i;*/ for(i=1; i<*n; ++i){ (*ex)[i] = (*ex)[0]; (*ey)[i] = (*ey)[0]; (*ez)[i] = (*ez)[0]; }
      break;
  }
  // for each particle: set speed and velocity vector
  for( i=0; i < *n; ++i ){
      // set the particle's speed to the preferred speed
      ( *v)[i] = *v0;
      // set the velocity vector using the vector's magnitude (i.e., the speed) and the direction of the velocity vector
      (*vx)[i] = (*v)[i] * (*ex)[i];
      (*vy)[i] = (*v)[i] * (*ey)[i];
      (*vz)[i] = (*v)[i] * (*ez)[i];
  }

  // ------- initializing other variables --------
  *warmUpTime = 1.0 * *n * *l / *v0; // simulation "warm up time" when starting from the ordered state, see in the header of this file at "output"


  // ------- set variables for displaying the simulation on the X11 output, start X11 output ------
#ifdef XWIN
  // set magnification of X11 image
  * x11picMagn = x11displayAreaSize / * l;
  // set width and height of X11 window, set width and upper coordinate of graph field
  * x11winWidth = x11infoFieldWidth + x11displayAreaSize + x11margin;
  * x11graphFieldWidth = *x11winWidth - 2 * x11margin;
  * x11graphFieldUpEnd = MAX( x11infoFieldHeight, x11displayAreaSize ) + 2 * x11margin;
  * x11winHeight = *x11graphFieldUpEnd + x11graphFieldHeight + 2 * x11margin;
 
  // initialize X11 display, npc: number of particle colors
  XColor * xcol = (XColor*)calloc( 3 + npc, sizeof(XColor) ); // list of data structures for particle colors and other colors
  * colorCode = (int*)calloc( 3 + npc, sizeof(int) );
  g_win( "sim", "sim", 0, 0, * x11winWidth, * x11winHeight, 4 );
  g_font( x11fontName );
  * x11graph = (double*)calloc( * x11graphFieldWidth, sizeof(double) ); //for(i=0; i < * x11graphFieldWidth; ++i) { ( * x11graph )[i] = 0.0; }

  // set data structures for particle colors (shades of red)
  for(i=0;i<npc;++i){ xcol[2+i].flags = DoRed|DoGreen|DoBlue; xcol[i].red = 0; xcol[i].green = (ush)floor(65535*i/(npc-1)); xcol[i].blue = 0; }
  // set data structures for background (black), info color (wheat) and particle boundary (dark red)
  xcol[  npc].flags = DoRed | DoGreen | DoBlue; xcol[  npc].red = 0;     xcol[  npc].green = 0;     xcol[  npc].blue = 0;
  xcol[1+npc].flags = DoRed | DoGreen | DoBlue; xcol[1+npc].red = 63567; xcol[1+npc].green = 57015; xcol[1+npc].blue = 45875;
  xcol[2+npc].flags = DoRed | DoGreen | DoBlue; xcol[2+npc].red = 0;     xcol[2+npc].green = 30000; xcol[2+npc].blue = 0;

  // allocate the colors
  for(i=0; i<3+npc; ++i) {
      // try to allocate color and check, if successful
      if( !XAllocColor(display,cmap,xcol+i) ) { fprintf( stderr, "Error: couldn't allocate the %d. color.\n",(1+i)); exit(0); }
      // IF the color allocation was successful, THEN save the color code
      (*colorCode)[i] = xcol[i].pixel;
  }

  // set the displayed length of the efficiency graph to zero
  * x11graphLenNow = 0;
#endif // ifdef XWIN


  // --------- data output -----------
  // write data header to stdout
  fprintf( stdout, "# Simulation time\n#\tEfficiency\n#\t\tDirection of motion (three coordinates of unit vector)\n\n" ); fflush( stdout );

  // --------- send message: finished init ----------
  fprintf(stderr,"done.\n"); fflush(stderr);
}

// ------------------------------------------------------------------

void EfficiencyAndDirectionOfMotion( int n, double * vx, double * vy, double * vz, double v0, double * eff, double * ex, double * ey, double * ez )
{
  // vectorial sum of the velocities of all particles
  double vxSum = 0.0, vySum = 0.0, vzSum = 0.0; int i; for( i=0; i<n; ++i ){ vxSum += vx[i]; vySum += vy[i]; vzSum += vz[i]; }
  // divide the length of the sum vector by the number of particles and the preferred velocity
  double len = sqrt( SQR(vxSum) + SQR(vySum) + SQR(vzSum) ); * eff = len / ( n * v0 );
  // the unit vector pointing in the direction of the vectorial sum of all movements
  * ex = vxSum / len; * ey = vySum / len; * ez = vzSum / len;
}

// ------------------------------------------------------------------

#ifdef XWIN
void XPutPic( int * colorCode, int nParticleColors, int x11displayAreaSize, int x11winWidth, int x11winHeight, int n, double simTimeNow,
	      double ** x, double ** y, double ** z, double ** vx, double ** vy, double ** vz, double v0, double ** ex, double ** ey, double ** ez,
	      double tau, double s, double dt, double l, double r,
	      int x11drawMethod, int x11infoFieldWidth, int x11margin, double x11picMagn, double x11objMagn, int drawSleepWT,
	      double ** x11graph, int x11graphFieldWidth, int x11graphFieldHeight, int * x11graphLenNow, int x11graphFieldUpEnd, int x11lineHeight )
{
  // 1, cleaning the window: filling the entire window with the background color
  // 2, drawing particles
  // 3, drawing the boundary of the simulated area
  // 4, updating and drawing the efficiency graph
  // 5, write info in the left info area
  // 6, showing it all (flushing)


  // 1, set color to black, fill entire window with this color
  XSetForeground( display, gc, colorCode[nParticleColors] ); XFillRectangle( display, pix1, gc, 0, 0, x11winWidth, x11winHeight );


  // 2, set color to brightest particle color, draw particles
  XSetForeground( display, gc, colorCode[ nParticleColors - 1 ] );
  // variables used for drawing each particle (lxm: left x margin, uym: up y margin, pm: picture magnification, om: object magnification)
  int lxm = x11infoFieldWidth; int uym = x11margin; double pm = x11picMagn; double om = x11objMagn; int i; 
  // method of drawing particles depends on this switch
  switch( x11drawMethod ){
  case 0: { 
      // set cosine and sine of 0.4 (radian)
      double cos04 = cos( 0.4 ); double sin04 = sin( 0.4 );
      // draw all particles
      for(i=0; i<n; ++i){
	  // set particle color based on the "z" coordinate of the particle, use only 80% of the particle colors from brightest to 20% of darkest
	  int iz = (int)floor ( nParticleColors * .8 * (*z)[i]/l ); if( iz == nParticleColors ){ --iz; }; XSetForeground( display, gc, colorCode[ iz ] );
	  // the body of the arrow
	  XDrawLine( display, pix1, gc,
		     (int) round ( lxm + pm *   (*x)[i] ),
		     (int) round ( uym + pm *   (*y)[i] ),
		     (int) round ( lxm + pm * ( (*x)[i] + om * (*vx)[i] ) ),
		     (int) round ( uym + pm * ( (*y)[i] + om * (*vy)[i] ) ) );
	  // the head of the arrow (only on one side)
	  XDrawLine( display, pix1, gc,
		     (int) round ( lxm + pm * ( (*x)[i] + om * (*vx)[i] ) ),
		     (int) round ( uym + pm * ( (*y)[i] + om * (*vy)[i] ) ),
		     (int) round ( lxm + pm * ( (*x)[i] + om * (*vx)[i] - om * 0.4 * ( cos04 * (*vx)[i] - sin04 * (*vy)[i] ) ) ),
	             (int) round ( uym + pm * ( (*y)[i] + om * (*vy)[i] - om * 0.4 * ( sin04 * (*vx)[i] + cos04 * (*vy)[i] ) ) ) );
      }
      break;
  }
    /*
  case 1: {
      // fill circles with the dimmest particle color
      XSetForeground( display, gc, colorCode[ 1 ] );
      for(i=0; i<n; ++i){ XFillCircle(display, pix1, gc, (int) round ( lxm + pm * x[i] ), (int) round ( uym + pm * y[i] ), (int)round ( pm * .5 * r ) ); }
      // draw circles with the 2nd brightest particle color
      XSetForeground( display, gc, colorCode[ nParticleColors - 2 ] );
      for(i=0; i<n; ++i){ XDrawCircle(display, pix1, gc, (int) round ( lxm + pm * x[i] ), (int) round ( uym + pm * y[i] ), (int)round ( pm * .5 * r ) ); }
      break;
  }
    */
  default: {
      fprintf(stderr,"\n\tError, no draw method with this code: %d\n\n",x11drawMethod);fflush(stderr);exit(-1);
  }
  }


  // 3, set the color to the info color and draw the boundary of the simulated circular area with dots
  XSetForeground( display, gc, colorCode[ 1 + nParticleColors ] );
  // rxm: right x margin, lym: low y margin
  int rxm = (int) round ( lxm + x11displayAreaSize ); int lym = (int) round ( uym + x11displayAreaSize );
  // draw dots at 3 pixels distance from each other, use the "lxm" and "uym" variables from above (see before the section where particles were drawn)
  for( i = 0; 4*i <= x11displayAreaSize; ++i ){
      XDrawPoint( display, pix1, gc, (int) round ( lxm + 4.0 * i ), uym                           ); // upper ("up" on screen) boundary
      XDrawPoint( display, pix1, gc,               lxm            , (int) round ( uym + 4.0 * i ) ); // left boundary
      XDrawPoint( display, pix1, gc, (int) round ( lxm + 4.0 * i ), lym                           ); // lower ("low" on screen) boundary
      XDrawPoint( display, pix1, gc,               rxm            , (int) round ( uym + 4.0 * i ) ); // left boundary
  }


  // 4, update and draw the efficiency graph
  //
  // 4a, update the efficiency graph
  //
  // IF the efficiency graph has reached its maximum length,
  // THEN shift all values to the left, i.e., toward the front of the list
  if( *x11graphLenNow == x11graphFieldWidth ){ int i; for( i = 0 ; i <= x11graphFieldWidth - 2 ; ++i ){ (*x11graph)[i] = (*x11graph)[i+1]; }}
  // ELSE: increment the length of the graph
  else{ ++ *x11graphLenNow; }
  //
  // save the current simulation time and the current efficiency value at the current position of the efficiency graph
  double eff, dir_x,dir_y,dir_z; EfficiencyAndDirectionOfMotion( n, *vx, *vy, *vz, v0, & eff, & dir_x, & dir_y, & dir_z ); (*x11graph)[ *x11graphLenNow - 1 ] = eff;
  //
  // 4b, set the color to the info color and draw the efficiency graph
  XSetForeground( display, gc, colorCode[ 1 + nParticleColors ] );
  // draw the frames (axes) for the displayed efficiency graph
  int ig; for(ig=0; 4.0 * ig < x11graphFieldWidth; ++ig) { // two horizontal dotted lines
      XDrawPoint( display, pix1, gc, x11margin + 4 * ig, x11graphFieldUpEnd                       );
      XDrawPoint( display, pix1, gc, x11margin + 4 * ig, x11graphFieldUpEnd + x11graphFieldHeight );
  }
  for(ig=0; 4.0 * ig < x11graphFieldHeight; ++ig) { // two vertical dotted lines
      XDrawPoint( display, pix1, gc, x11margin,                      x11graphFieldUpEnd + 4 * ig );
      XDrawPoint( display, pix1, gc, x11margin + x11graphFieldWidth, x11graphFieldUpEnd + 4 * ig );
  }
  // draw the points of the efficiency graph
  for(ig=1; ig < *x11graphLenNow; ++ig){
      XDrawLine( display, pix1, gc,
		 x11margin + ( ig - 1 ),
		 x11graphFieldUpEnd + (int) x11graphFieldHeight * ( 1.0 - (*x11graph)[ ig - 1 ] ),
		 x11margin +   ig,
		 x11graphFieldUpEnd + (int) x11graphFieldHeight * ( 1.0 - (*x11graph)[ ig     ] ) );
  }


  // 5, set color to info color, write info in the left info area, dispHi: vertical coordinate of the current displayed string, dispStr: currently displayed string
  XSetForeground( display, gc, colorCode[ 1 + nParticleColors ] );
  // Write label of the graph
  char dispStr[MY_STRLEN]; sprintf( dispStr, "E(t)" );
  XDrawString( display, pix1, gc, x11margin, x11graphFieldUpEnd - .2 * x11lineHeight, dispStr, (signed int)strlen(dispStr) );
  // Side length of the simulated rectangular field
  int dispHi = x11margin + x11lineHeight; sprintf( dispStr, "L = %.2f", l );
  XDrawString( display, pix1, gc, x11margin, dispHi, dispStr, (signed int)strlen(dispStr) );
  // Number of particles
  dispHi += x11lineHeight; sprintf( dispStr, "N = %d", n );
  XDrawString( display, pix1, gc, x11margin, dispHi, dispStr, (signed int)strlen(dispStr) );
  // Preferred magnitude of the velocity
  dispHi += x11lineHeight; sprintf( dispStr, "V0 = %g", v0 );
  XDrawString( display, pix1, gc, x11margin, dispHi, dispStr, (signed int)strlen(dispStr) );
  // Time constant used by a particle to adjust its magnitude of velocity to the preferred magnitude of velocity
  dispHi += x11lineHeight; sprintf( dispStr, "Tau = %g", tau );
  XDrawString( display, pix1, gc, x11margin, dispHi, dispStr, (signed int)strlen(dispStr) );
  // Interaction cutoff radius
  dispHi += x11lineHeight; sprintf( dispStr, "R = %g", r );
  XDrawString( display, pix1, gc, x11margin, dispHi, dispStr, (signed int)strlen(dispStr) );
  // Noise amplitude
  dispHi += x11lineHeight; sprintf( dispStr, "S = %g", s );
  XDrawString( display, pix1, gc, x11margin, dispHi, dispStr, (signed int)strlen(dispStr) );
  // Default simulation time step
  dispHi += x11lineHeight; sprintf( dispStr, "dt = %g", dt );
  XDrawString( display, pix1, gc, x11margin, dispHi, dispStr, (signed int)strlen(dispStr) );
  // Current simulation time
  dispHi += x11lineHeight; sprintf( dispStr, "t = %.2f", simTimeNow );
  XDrawString( display, pix1, gc, x11margin, dispHi, dispStr, (signed int)strlen(dispStr) );
  // Current efficiency
  dispHi += x11lineHeight; sprintf( dispStr, "E = %.2f", (*x11graph)[*x11graphLenNow-1] );
  XDrawString( display, pix1, gc, x11margin, dispHi, dispStr, (signed int)strlen(dispStr) );


  // 6, flush (show) the entire image that we have drawn
  h_show( x11winWidth, x11winHeight );
  // if requested, wait after drawing
  if( 0 < drawSleepWT ){ sleep( drawSleepWT ); }
}
#endif

// ------------------------------------------------------------------

void WriteOutputData( int n, double * x, double * y, double * z, double * vx, double * vy, double * vz, double v0, double stNow )
{
  // compute the efficiency of the motion and write it to stdout
  double eff, ex, ey, ez; EfficiencyAndDirectionOfMotion( n, vx, vy, vz, v0, & eff, & ex, & ey, & ez );
  fprintf(stdout,"%g\t%g\t%g\t%g\t%g\n",stNow,eff,ex,ey,ez); fflush(stdout);
}

// ------------------------------------------------------------------

// The distance between the centers of two particles is  "dist"
// The interaction force acts between the two particles if  "dist" <= "r"  where "r" is the cutoff radius and  cutoff_rSqr = r * r
void AddInteractionsToForces( int bn, int **** bxyz2n, int **** bxyz2i, int **** bp, int ** pn, double ** x, double ** y, double ** z,
			      double cutoff_rSqr, double ** fxSum, double ** fySum, double ** fzSum, double l )
{
  // 2b, repulsion among interacting particles
  // loop through the list of the non-empty book-keeping grid cells
  // bxyz2n[ix][iy][iz]: number of particles in grid cell (ix,iy,iz) has to be >0
  int ix; for(ix=0; ix<bn; ++ix){ int iy; for(iy=0; iy<bn; ++iy){ int iz; for(iz=0; iz<bn; ++iz){ if( (*bxyz2n)[ix][iy][iz] ){
      // (dx,dy,dz): relative position of the other grid cell from the current grid cell
      int dx; for(dx=-1; dx<=1; ++dx){ int dy; for(dy=-1; dy<=1; ++dy){ int dz; for(dz=-1; dz<=1; ++dz){
	  // IF the investigated other grid cell is the same as the current grid cell, THEN proceed only if this grid cell contains at least two particles 
	  if( ! dx && ! dy && ! dz ){ // if( 0 == dx && 0 == dy && 0 == dz ){
	      if( 2 <= (*bxyz2n)[ix][iy][iz] ){
		  // loop through the list of particles in the current grid cell
		  // j1: the index of the first particle in the current grid cell
		  int j1 = (*bp)[ix][iy][iz]; while( CHAIN_END != j1 ){
		      // j2: the index of the other particle in the current grid cell
		      // NOTE: each pair of particles is considered only once
		      int j2 = (*pn)[j1]; while( CHAIN_END != j2 ){
			  // signed (x,y,z) components of the distance from particle j2 to j1
			  double signed_x_dist_from_j2_to_j1 = (*x)[j1] - (*x)[j2];
			  double signed_y_dist_from_j2_to_j1 = (*y)[j1] - (*y)[j2];
			  double signed_z_dist_from_j2_to_j1 = (*z)[j1] - (*z)[j2];
			  // square of the two particles' distance
			  double distSqr = SQR( signed_x_dist_from_j2_to_j1 ) + SQR( signed_y_dist_from_j2_to_j1 ) + SQR( signed_z_dist_from_j2_to_j1 );
			  // the two particles interact only if the square of their distance is below the square of the interaction cutoff radius
			  if( distSqr < cutoff_rSqr ){
			      // r is the distance between the two particles (both are thought to have zero size)
			      // F = (1/r^2) is the magnitude of the repulsion force between the two particles
			      // Therefore, F / r  =  1 / r^3 =  ( r^2 ) ^ (-3/2)
			      double f_over_r = pow( distSqr, -1.5 );
			      double fx = signed_x_dist_from_j2_to_j1 * f_over_r;
			      double fy = signed_y_dist_from_j2_to_j1 * f_over_r;
			      double fz = signed_z_dist_from_j2_to_j1 * f_over_r;
			      // increment the components of the sum of forces acting on particle j1
			      (*fxSum)[j1] += fx;
			      (*fySum)[j1] += fy;
			      (*fzSum)[j1] += fz;
			      // values with opposite signs act on particle j2
			      (*fxSum)[j2] -= fx;
			      (*fySum)[j2] -= fy;
			      (*fzSum)[j2] -= fz;
			  }
			  // step to next particle or the "CHAIN_END" marker of the current book-keeping field 
			  j2 = (*pn)[j2];
		      }
		      // step to next particle or the "CHAIN_END" marker of the current book-keeping field 
		      j1 = (*pn)[j1];
		  }
	      }
	  }
	  // IF the investigated other grid cell is different from the current grid cell
	  else{
 	      // (jx,jy,jz): the coordinates of the investigated other grid cell that is at a relative position (dx,dy,dz) from the current grid cell
	      // xo: due to the periodic boundaries in the book-keeping the relative position of the neighboring cell 
	      //     is ( dx + xo * BN , dy + yo * BN , dz + zo * BN ) from the current cell  ( xo can be = -1, 0, +1 )
	      // yo: similar to xo
	      int jx = ix + dx; int xo = 0; if( jx < 0 ){ jx += bn; --xo; }else if( jx > bn-1 ){ jx -= bn; ++xo; }
	      int jy = iy + dy; int yo = 0; if( jy < 0 ){ jy += bn; --yo; }else if( jy > bn-1 ){ jy -= bn; ++yo; }
	      int jz = iz + dz; int zo = 0; if( jz < 0 ){ jz += bn; --zo; }else if( jz > bn-1 ){ jz -= bn; ++zo; }
	      //
	      // proceed only:
	      // IF  the grid cell index of the other grid cell is larger than the current grid cell's index
	      //     - note that this condition ensures that each pair of particles is considered only once
              if(    ( (*bxyz2i)[jx][jy][jz] > (*bxyz2i)[ix][iy][iz] )
	      // AND the other grid cell contains at least one particle
		  && (*bxyz2n)[jx][jy][jz]
		)
	      {
		  // loop through the list of particles in the _current_ grid cell, k1: the index of the particle in the _current_ grid cell
 		  int k1 = (*bp)[ix][iy][iz]; while( CHAIN_END != k1 ){
		      // loop through the list of particles in the _other_ grid cell, k2: the index of the particle in the _other_ grid cell
		      int k2 = (*bp)[jx][jy][jz]; while( CHAIN_END != k2 ){
			  // signed (x,y,z) components of the distance from particle k2 to k1
			  double signed_x_dist_from_k2_to_k1 = (*x)[k1] - ( (*x)[k2] + xo * l );
			  double signed_y_dist_from_k2_to_k1 = (*y)[k1] - ( (*y)[k2] + yo * l );
			  double signed_z_dist_from_k2_to_k1 = (*z)[k1] - ( (*z)[k2] + zo * l );
			  // square of the two particles' distance
			  double distSqr = SQR( signed_x_dist_from_k2_to_k1 ) + SQR( signed_y_dist_from_k2_to_k1 ) + SQR( signed_z_dist_from_k2_to_k1 );
			  // the two particles interact only if the square of their distance is below the square of the interaction cutoff radius
			  if( distSqr < cutoff_rSqr ){
			      // r is the distance between the two particles (both are thought to have zero size)
			      // F = (1/r^2) is the magnitude of the repulsion force between the two particles
			      // Therefore, F / r  =  1 / r^3 =  ( r^2 ) ^ (-3/2)
			      double f_over_r = pow( distSqr, -1.5 );
			      double fx = signed_x_dist_from_k2_to_k1 * f_over_r;
			      double fy = signed_y_dist_from_k2_to_k1 * f_over_r;
			      double fz = signed_z_dist_from_k2_to_k1 * f_over_r;
			      // increment the sum of force components acting on particle j1
			      (*fxSum)[k1] += fx;
			      (*fySum)[k1] += fy;
			      (*fzSum)[k1] += fz;
			      // values with opposite signs act on particle j2
			      (*fxSum)[k2] -= fx;
			      (*fySum)[k2] -= fy;
			      (*fzSum)[k2] -= fz;
			  }
			  // step to next particle or the "CHAIN_END" marker of the current book-keeping field 
			  k2 = (*pn)[k2];
		      }
		      // step to next particle or the "CHAIN_END" marker of the current book-keeping field 
		      k1 = (*pn)[k1];
		  }
	      }
	  }
      }}}}
   }}}
}

// ------------------------------------------------------------------

void UpdateSim_MidpMeth(
			 double dt, int mtl, int * ic, int * icn, int * is_first_update, int n, double v0, double tau, double s, 
			 double l, double cutoff_rSqr, double * stNow, double * stPrev, uli * wtNow,

			 double ** x,    double ** y,    double ** z,    double ** vx,      double ** vy,      double ** vz,      double ** v,
			 double ** x_m,  double ** y_m,  double ** z_m,  double ** vx_m,    double ** vy_m,    double ** vz_m,    double ** v_m,
			 double ** ex,   double ** ey,   double ** ez,   double *** fxSum,   double *** fySum,   double *** fzSum, 
			 double ** ex_m, double ** ey_m, double ** ez_m, double *** fxSum_m, double *** fySum_m, double *** fzSum_m,

			 double bl, int bn, int **** bxyz2i,

			 int ** pfx,   int ** pfy,   int ** pfz,   int **** bp,   int ** pn,   int **** bxyz2n, 
			 int ** pfx_m, int ** pfy_m, int ** pfz_m, int **** bp_m, int ** pn_m, int **** bxyz2n_m 
		      )
{
  // ---------------------------------------
  // NOTE:
  // - Forces for the current time are computed in (1) and (3) below,
  //   but the forces from tau time ago are applied for updating velocities in (2) and (4)
  // - If we are at the first simulation update, then the forces computed in (1) are copied into all past values,
  //   similarly the midpoint force values computed in (3) are copied into all past midpoint values
  //
  // 1, compute the self-propelling and interaction forces acting on each particle
  //    1.a, set all force sums to zero
  //    1.b, self-propelling forces: loop through the list of particles
  //    1.c, particle-particle interactions: loop through the list of non-empty book-keeping grid cells
  //         - a particle interacts with all other particles in the 27 nearby book-keeping grid cells
  //           (the grid cell itself and the 26 neighboring grid cells)
  //         - for each interacting particle-particle pair: increment the force components of both interacting particles
  //         - for each book-keeping grid cell:
  //           . consider interactions inside the grid cell itself only if it contains at least 2 particles
  //           . consider interactions with a neighboring grid cell only if both grid cells contain at least 1 particle
  //           . to compute an interaction between two particles only once:
  //             consider the neighboring grid cell only if its grid cell index is above the current grid cell's index
  //    1.d, add noise for a time interval of dt/2
  //         - noise amplitude is S * sqrt(dt/2)
  //    1.e, IF we are at the first simulation update, THEN copy the computed forces to all past values
  //
  // 2, based on the current coordinates and velocities
  //    2a, compute the coordinates and velocities at the midpoint: dt/2 time after the current simulation time, t
  //    2b, IF a particle has left the simulated field in either direction, THEN put it back using the periodic boundary conditions
  //    2c, particle -> grid cell index book-keeping: compute grid cell indexes at the midpoint
  //    2d, grid cell index -> particle book-keeping: set book-keeping arrays at the midpoint
  //
  // 3, compute forces at the midpoint
  //    3a, set force sums to zero
  //    3b, self-propelling
  //    3c, interactions
  //    3d, noise, magnitude is S * sqrt(dt)
  //    3e, IF we are at the first simulation update, THEN copy the computed midpoint forces values to all past midpoint values
  //
  // 4, change the initial coordinate and velocity vectors using the velocities and forces at the midpoint
  //
  // 5, update simulation time, update current index and next index in the circular arrays, set is_first_update to 0
  // ---------------------------------------

  // -------- Compute the forces acting on each particle at the current time step --------

  // 1a, for each particle: initialize the components of the sum of forces acting on the particle
  int i; for(i=0; i<n; ++i){ (*fxSum)[*ic][i] = 0.0; (*fySum)[*ic][i] = 0.0; (*fzSum)[*ic][i] = 0.0; };

  // 1b, self-propelling force with a magnitude of  q = (v_0-v)/tau
  // IF the magnitude of the particle's current velocity is not v0, THEN it has a non-zero self-propelling force
  double one_over_tau = 1.0 / tau;
  for(i=0; i<n; ++i){
      // compute the magnitude of the self-propelling force
      double q = ( v0 - (*v)[i] ) * one_over_tau;
      // add components of the self-propelling force to the components of the sum of forces acting on the <i>th particle
      (*fxSum)[*ic][i] += q * (*ex)[i];
      (*fySum)[*ic][i] += q * (*ey)[i];
      (*fzSum)[*ic][i] += q * (*ez)[i];
  }
    
  // 1c, add interactions to the forces
  AddInteractionsToForces( bn, bxyz2n, bxyz2i, bp, pn, x, y, z, cutoff_rSqr, (*fxSum)+*ic, (*fySum)+*ic, (*fzSum)+*ic, l );

  // 1d, IF the noise amplitude is non-zero, THEN add noise to the forces
  //     NOTE:  amplitude of noise is   sqrt ( dt / 2 )
  if( s > 0.0 ){ 
      // noise amplitude/magnitude (length of the noise vector)
      double noiseMagn = s * sqrt( .5 * dt );
      for(i=0; i<n; ++i){
	  // generate a unit vector pointing in a direction distributed evenly in the 4 PI solid angle
	  double noiseDir_x, noiseDir_y, noiseDir_z;
	  Generate3dVector_randomDir_unitLength( & noiseDir_x, & noiseDir_y, & noiseDir_z );
	  // add this unit vector multiplied by 'noiseMagn' to the sum of the forces acting on the <i>th particle
	  (*fxSum)[*ic][i] += noiseMagn * noiseDir_x;
	  (*fySum)[*ic][i] += noiseMagn * noiseDir_y;
	  (*fzSum)[*ic][i] += noiseMagn * noiseDir_z;
      }
  }

  // 1e, IF we are at the first simulation update, THEN copy current values into all past values
  if( * is_first_update ){
      int it; for(it=1; it<=mtl; ++it){ // at the first update ic = 0, so the range of all other values is it=1..mtl
	  (*fxSum)[it][i] = (*fxSum)[*ic][i]; 
	  (*fySum)[it][i] = (*fySum)[*ic][i]; 
	  (*fzSum)[it][i] = (*fzSum)[*ic][i]; 
      }
  }


  // -------- Use the forces from tau time ago to update the velocity components and coordinates of all particles --------

  // 2a, compute positions, velocities and directions at the midpoint
  for(i=0; i<n; ++i){
      // position (coordinates) of the <i>th particle at the midpoint
      (* x_m)[i] = (* x)[i] + .5 * dt * (*vx   )[i];
      (* y_m)[i] = (* y)[i] + .5 * dt * (*vy   )[i];
      (* z_m)[i] = (* z)[i] + .5 * dt * (*vz   )[i];
      // compute velocity vector at the midpoint based on the sum of forces from the past
      (*vx_m)[i] = (*vx)[i] + .5 * dt * (*fxSum)[*icn][i];
      (*vy_m)[i] = (*vy)[i] + .5 * dt * (*fySum)[*icn][i];
      (*vz_m)[i] = (*vz)[i] + .5 * dt * (*fzSum)[*icn][i];
      // speed at the midpoint
      (*v_m)[i] = sqrt( SQR( (*vx_m)[i] ) + SQR( (*vy_m)[i] ) + SQR( (*vz_m)[i] ) );
      double one_over_v_m = 1.0 / (*v_m)[i];
      // direction of the velocity vector at the midpoint
      (*ex_m)[i] = (*vx_m)[i] * one_over_v_m;
      (*ey_m)[i] = (*vy_m)[i] * one_over_v_m;
      (*ez_m)[i] = (*vz_m)[i] * one_over_v_m;
  }
  // 2b, particles leaving the simulated volume should be re-inserted using the periodic boundary conditions
  //     NOTE: the interval [0,L) is closed from below and open from above
  for(i=0; i<n; ++i){
      if( (*x_m)[i] < 0 ){ (*x_m)[i] += l; }else if( (*x_m)[i] >= l ){ (*x_m)[i] -= l; }
      if( (*y_m)[i] < 0 ){ (*y_m)[i] += l; }else if( (*y_m)[i] >= l ){ (*y_m)[i] -= l; }
      if( (*z_m)[i] < 0 ){ (*z_m)[i] += l; }else if( (*z_m)[i] >= l ){ (*z_m)[i] -= l; }
  }  
  // 2c, compute grid cell indexes at the midpoint, bn_over_l = bn/l
  double bn_over_l = 1.0*bn/l; for(i=0; i<n; ++i){ CoordsToGridCellIndexes( (*x_m)[i], (*y_m)[i], (*z_m)[i], bn, bn_over_l, *pfx_m + i, *pfy_m + i, *pfz_m + i ); }
  // 2d, set book-keeping data structures at the midpoint
  // default: identical to original values
  for(i=0; i<n; ++i){ (*pn_m)[i] = (*pn)[i]; }; 
  int ix; for(ix=0; ix<bn; ++ix){ int iy; for(iy=0; iy<bn; ++iy){ int iz; for(iz=0; iz<bn; ++iz){
      (*bp_m)[ix][iy][iz] = (*bp)[ix][iy][iz]; (*bxyz2n_m)[ix][iy][iz] = (*bxyz2n)[ix][iy][iz]; }}}
  // check all particles: IF there is a change, THEN change the book-keeping table accoringly
  for(i=0; i<n; ++i){  
      // IF the current particle has moved to a different book-keeping cell, THEN change its book-keeping tables at the midpoint
      if( (*pfx_m)[i] != (*pfx)[i] || (*pfy_m)[i] != (*pfy)[i] || (*pfz_m)[i] != (*pfz)[i] ){
	  // in the midpoint book-keeping (bp_m, pn_m, bxyz2n_m) remove the particle from its original grid cell (pfx[i], pfy[i], pfz[i])
	  RemoveParticleFromBookKeepingTable( i, (*pfx  )[i], (*pfy  )[i], (*pfz  )[i], bp_m, pn_m, bxyz2n_m );
	  // in the midpoint book-keeping (bp_m, pn_m, bxyz2n_m) insert the particle into its grid cell at the midpoint (pfx_m[i], pfy_m[i], pfz_m[i])
          InsertParticleIntoBookKeepingTable( i, (*pfx_m)[i], (*pfy_m)[i], (*pfz_m)[i], bp_m, pn_m, bxyz2n_m );
      }
  }

  // ----------------------------------------

  // 3a, for each particle: initialize the components of the sum of forces acting on the particle
  for(i=0; i<n; ++i){ (*fxSum_m)[*ic][i] = 0.0; (*fySum_m)[*ic][i] = 0.0; (*fzSum_m)[*ic][i] = 0.0; };

  // 3b, self-propelling force
  // IF the magnitude of the particle's current velocity is not v0, THEN it has a non-zero self-propelling force
  /*double one_over_tau = 1.0 / tau;*/
  for(i=0; i<n; ++i){
      // compute the magnitude of the self-propelling force
      double q_m = ( v0 - (*v_m)[i] ) * one_over_tau;
      // add components of the self-propelling force to the components of the sum of forces acting on the <i>th particle
      (*fxSum_m)[*ic][i] += q_m * (*ex_m)[i];
      (*fySum_m)[*ic][i] += q_m * (*ey_m)[i];
      (*fzSum_m)[*ic][i] += q_m * (*ez_m)[i]; 
  }

  // 3c, add interactions to the forces
  AddInteractionsToForces( bn, bxyz2n_m, bxyz2i, bp_m, pn_m, x_m, y_m, z_m, cutoff_rSqr, (*fxSum_m)+*ic, (*fySum_m)+*ic, (*fzSum_m)+*ic, l );

  // 3d, IF the noise amplitude is non-zero, THEN add noise to the forces
  //     NOTE: amplitude of noise is   sqrt ( dt )
  if( s > 0.0 ){
      // noise amplitude/magnitude (length of the noise vector)
      double noiseMagn = s * sqrt( dt );
      for(i=0; i<n; ++i){
	  // generate a unit vector pointing in a direction distributed evenly in the 4 PI solid angle
	  double noiseDir_x, noiseDir_y, noiseDir_z;
	  Generate3dVector_randomDir_unitLength( & noiseDir_x, & noiseDir_y, & noiseDir_z );
	  // add this unit vector multiplied by 'noiseMagn' to the sum of the forces acting on the <i>th particle
	  (*fxSum_m)[*ic][i] += noiseMagn * noiseDir_x;
	  (*fySum_m)[*ic][i] += noiseMagn * noiseDir_y;
	  (*fzSum_m)[*ic][i] += noiseMagn * noiseDir_z;
      }
  }

  // 3e, IF we are at the first simulation update, THEN copy current values into all past values
  if( * is_first_update ){
      int it; for(it=1; it<=mtl; ++it){ // at the first update ic = 0, so the range of all other values is it=1..mtl
	  (*fxSum_m)[it][i] = (*fxSum_m)[*ic][i]; 
	  (*fySum_m)[it][i] = (*fySum_m)[*ic][i]; 
	  (*fzSum_m)[it][i] = (*fzSum_m)[*ic][i]; 
      }
  }

  // ----------------------------------------

  // 4a, use values at the midpoint to compute positions, velocities and directions at the final point
  for(i=0; i<n; ++i){
      // position at the final point of the simulation update
      (* x)[i] = (* x_m)[i] + dt * (*   vx_m)[i];
      (* y)[i] = (* y_m)[i] + dt * (*   vy_m)[i];
      (* z)[i] = (* z_m)[i] + dt * (*   vy_m)[i];
      // velocity vector at the final point of the simulation update based on the sum of forces in the past
      (*vx)[i] = (*vx_m)[i] + dt * (*fxSum_m)[*icn][i];
      (*vy)[i] = (*vy_m)[i] + dt * (*fySum_m)[*icn][i];
      (*vz)[i] = (*vz_m)[i] + dt * (*fzSum_m)[*icn][i];
      // speed at the final point of the simulation update
      ( *v)[i] = sqrt( SQR( (*vx)[i] ) + SQR( (*vy)[i] ) + SQR( (*vz)[i] ) );
      double one_over_v = 1.0 / (*v)[i];
      // direction of the velocity at the final point
      (*ex)[i] = (*vx)[i] * one_over_v;
      (*ey)[i] = (*vy)[i] * one_over_v;
      (*ez)[i] = (*vz)[i] * one_over_v;
  }
  
  // 4b, particles leaving the simulated volume should be re-inserted using the periodic boundary conditions
  //     NOTE: the interval [0,L) is closed from below and open from above
  for(i=0; i<n; ++i){
      if( (*x)[i] < 0 ){ (*x)[i] += l; }else if( (*x)[i] >= l ){ (*x)[i] -= l; }
      if( (*y)[i] < 0 ){ (*y)[i] += l; }else if( (*y)[i] >= l ){ (*y)[i] -= l; }
      if( (*z)[i] < 0 ){ (*z)[i] += l; }else if( (*z)[i] >= l ){ (*z)[i] -= l; }
  }  

  // 4c, compute grid cell indexes at the final point, bn_over_l = bn/l
  /*double bn_over_l = 1.0*bn/l;*/ for(i=0; i<n; ++i){ CoordsToGridCellIndexes( (*x)[i], (*y)[i], (*z)[i], bn, bn_over_l, *pfx + i, *pfy + i, *pfz + i ); }

  // 4d, set book-keeping data structures at the final point
  //     default: identical to values at the midpoint
  for(i=0; i<n; ++i){ (*pn)[i] = (*pn_m)[i]; }; 
  /*int ix;*/ for(ix=0; ix<bn; ++ix){ int iy; for(iy=0; iy<bn; ++iy){ int iz; for(iz=0; iz<bn; ++iz){
      (*bp)[ix][iy][iz] = (*bp_m)[ix][iy][iz]; (*bxyz2n)[ix][iy][iz] = (*bxyz2n_m)[ix][iy][iz]; }}}
  // check all particles: IF there is a change, THEN change the book-keeping table accoringly
  for(i=0; i<n; ++i){
      // IF the current particle has moved to a different book-keeping cell, THEN change its book-keeping tables at the final time point
      if( (*pfx)[i] != (*pfx_m)[i] || (*pfy)[i] != (*pfy_m)[i] || (*pfz)[i] != (*pfz_m)[i] ){
	  // in the book-keeping of the final time point (bp, pn, bxyz2n) remove the particle from its grid cell at the midpoint (pfx_m[i], pfy_m[i], pfz_m[i])
	  RemoveParticleFromBookKeepingTable( i, (*pfx_m)[i], (*pfy_m)[i], (*pfz_m)[i], bp, pn, bxyz2n );
	  // in the book-keeping of the final time point (bp, pn, bxyz2n) insert the particle into its grid cell at the final time point (pfx[i], pfy[i], pfz[i])
	  InsertParticleIntoBookKeepingTable( i, (*pfx  )[i], (*pfy  )[i], (*pfz  )[i], bp, pn, bxyz2n );
      }
  }

  // ----------------------------------------

  // 5, update counters
  // previous and current simulation time, wall clock time (wt)
  * stPrev = * stNow; * stNow += dt; * wtNow = (uli)time(NULL);
  // current index in the circular arrays saving past data, next index (which is also the index of the time step tau time ago)
  * ic = * icn; * icn = nextNumberWithModulo( * ic, mtl + 1 );
  // even if this has been the first simulation update, the next will not be the first
  * is_first_update = 0;
}
