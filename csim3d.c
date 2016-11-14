// sim ( = simulation )
//
// ROUGH DOC:
// - continuous model of interacting particles in 3d
//
// DETAILED DOC:
// - read sim3d_lib.c (starting at the top, proceeding towards the bottom):
//   = global constants, parameters and variables
//   = at the beginning of each function definition
//
// COMPILING:
//  - see the file compile3d.pl in the current directory
//  - the two compiled files are: 'sim3d' (demo simulation)
//    and 'csim3d' (number-crunching, "c" stands for "crunching")
//
// RUNNING:
// - start either 'sim3d' or 'csim3d' with zero arguments to see the requested list of arguments
// - see the README file for examples
//
// ==================================================================

#include "sim3d_lib.c" // include library file

int main( int ArgC, char * ArgV[] )
{
  // initializing
  Init( ArgC, ArgV, & _RND_SEED, & _N, & _START_STATE, & _V0, & _WARM_UP_TIME, & _TAU, & _S, 
	& _DT, & _MTL, & _IC, & _ICN, & _IS_FIRST_UPDATE, & _L, & _R, & _BL, & _BN, & _BXYZ2I,
	& _PFX, & _PFY, & _PFZ, & _PFX_M, & _PFY_M, & _PFZ_M, & _BP, & _BP_M, & _PN, & _PN_M, & _BXYZ2N, & _BXYZ2N_M,
	& _ST_MAX, & _ST_NOW, & _ST_PREV, & _DRAW_FREQ_ST, & _SAVE_FREQ_N, & _DRAW_SLEEP_WT, & _WT_MAX, & _WT_START, & _WT_NOW,
        & _X,    & _Y,    & _Z,    & _VX,       & _VY,       & _VZ,       & _V,
	& _X_M,  & _Y_M,  & _Z_M,  & _VX_M,     & _VY_M,     & _VZ_M,     & _V_M,
	& _EX,   & _EY,   & _EZ,   & _FX_SUM,   & _FY_SUM,   & _FZ_SUM, 
	& _EX_M, & _EY_M, & _EZ_M, & _FX_SUM_M, & _FY_SUM_M, & _FZ_SUM_M,
	& _X11_DrawMethod, _X11DisplayAreaSize, & _X11WinWidth, _X11InfoFieldWidth, & _X11_PIC_MAGN, & _X11_OBJ_MAGN, _X11Margin,
        & _X11GraphFieldWidth, & _X11GraphFieldUpEnd, _X11InfoFieldHeight, & _X11WinHeight, _X11GraphFieldHeight,
        _NPC, & _ColorCode, & _X11Graph, _X11FontName, & _X11GraphLenNow );

  // loop: drawing, saving and updating the simulation
  // exit condition: having reached at least one of the maximum allowed values
  // negative threshold values mean: no threshold
  while(    ( 0 > _ST_MAX      || _ST_NOW             < _ST_MAX      )
	 && ( 0 > _WT_MAX      || _WT_NOW - _WT_START < _WT_MAX      ) )
  {
#ifdef XWIN
      // drawing data: (i) at the beginning, (ii) at the requested simulation time interval
      if( 0.0 == _ST_NOW || floor( _ST_NOW / _DRAW_FREQ_ST ) > floor( _ST_PREV / _DRAW_FREQ_ST ) ){
	  XPutPic( _ColorCode, _NPC, _X11DisplayAreaSize, _X11WinWidth, _X11WinHeight, _N, _ST_NOW, 
		   & _X, & _Y, & _Z, & _VX, & _VY, & _VZ, _V0, & _EX, & _EY, & _EZ,
		   _TAU, _S, _DT, _L, _R, _X11_DrawMethod, _X11InfoFieldWidth, _X11Margin, _X11_PIC_MAGN, _X11_OBJ_MAGN, _DRAW_SLEEP_WT,
		   & _X11Graph, _X11GraphFieldWidth, _X11GraphFieldHeight, & _X11GraphLenNow, _X11GraphFieldUpEnd, _X11LineHeight );
      }
#endif
      // . with n=SaveFreqN, saving data to STDOUT: 
      // . save at simulation time values of t_i = 1 * 10 ** ( i / n ) for i = 0, 1, 2, ... , saving first at simulation time = N * update time step
      if(    0.0 == _ST_NOW
	  || (    ( _ST_NOW >= _N * _DT )
	                 // IF the simulation is started from the disordered state, 
	       && (    (    ( 0 == _START_STATE )
			    // THEN save at the selected t_i = 1 * 10 ** ( i / n ) simulation times
			 && (   floor( _SAVE_FREQ_N * log( _ST_PREV ) / log( 10.0 ) )
			      < floor( _SAVE_FREQ_N * log( _ST_NOW  ) / log( 10.0 ) ) 
			    )
		       )
		    ||
		       // IF the simulation is started from the ordered state, 
		       (    ( 0 != _START_STATE )
			    // THEN after time = 0
			    // save at t_warmUp (warm-up time)
			 && (    ( _WARM_UP_TIME <= _ST_NOW && _ST_NOW < _WARM_UP_TIME + 2.0 * _DT )
			    // AND after that at the  t_i - t_warmUp  simulation times
			      || (    ( _ST_NOW >= _WARM_UP_TIME + _N * _DT )
				   && (   floor( _SAVE_FREQ_N * log( _ST_PREV - _WARM_UP_TIME ) / log( 10.0 ) )
					< floor( _SAVE_FREQ_N * log( _ST_NOW  - _WARM_UP_TIME ) / log( 10.0 ) ) )
				 )
			    )
		       )
		  )
	     )
	)
      {
	  WriteOutputData( _N, _X, _Y, _Z, _VX, _VY, _VZ, _V0, _ST_NOW );
      }
      // Update the simulation with the midpoint method
      // Note: IF the simulation is started from the ordered state, THEN before reaching warmUpTime the applied noise is zero
      UpdateSim_MidpMeth(
			  _DT, _MTL, & _IC, & _ICN, & _IS_FIRST_UPDATE, _N, _V0, _TAU, ( ( ( 0 != _START_STATE ) && ( _ST_NOW < _WARM_UP_TIME ) ) ? 0.0 : _S ),
			  _L, SQR( _R ), & _ST_NOW, & _ST_PREV, & _WT_NOW,

			  & _X,    & _Y,    & _Z,    & _VX,       & _VY,       & _VZ,       & _V, 
			  & _X_M,  & _Y_M,  & _Z_M,  & _VX_M,     & _VY_M,     & _VZ_M,     & _V_M, 
			  & _EX,   & _EY,   & _EZ,   & _FX_SUM,   & _FY_SUM,   & _FZ_SUM, 
			  & _EX_M, & _EY_M, & _EZ_M, & _FX_SUM_M, & _FY_SUM_M, & _FZ_SUM_M, 

			  _BL, _BN, & _BXYZ2I,

			  & _PFX,   & _PFY,   & _PFZ,   & _BP,   & _PN,   & _BXYZ2N, 
			  & _PFX_M, & _PFY_M, & _PFZ_M, & _BP_M, & _PN_M, & _BXYZ2N_M 
			);
  }

  // write data to stdout at the end of the simulation
  WriteOutputData( _N, _X, _Y, _Z, _VX, _VY, _VZ, _V0, _ST_NOW );

  // done
  return 1;
}
