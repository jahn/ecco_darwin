C $Header: /u/gcmpack/MITgcm_contrib/ecco_darwin/v4_llc270_JAMES_budget/code_darwin/PTRACERS_SIZE.h,v 1.1 2020/02/29 01:06:09 dcarroll Exp $
C $Name:  $

#ifdef ALLOW_PTRACERS

CBOP
C    !ROUTINE: PTRACERS_SIZE.h
C    !INTERFACE:
C #include PTRACERS_SIZE.h
 
C    !DESCRIPTION:
C Contains passive tracer array size (number of tracers).

C PTRACERS_num defines how many passive tracers are allocated/exist.
C  and is set here (default 1)
C
C     Number of tracers
      INTEGER PTRACERS_num
      PARAMETER(PTRACERS_num = 39 )

#ifdef ALLOW_AUTODIFF_TAMC
      INTEGER    iptrkey
      INTEGER    maxpass
      PARAMETER( maxpass     = PTRACERS_num + 2 )
#endif

CEOP
#endif /* ALLOW_PTRACERS */

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
