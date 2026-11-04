#include "OFFLINE_OPTIONS.h"

CBOP
C     !ROUTINE: OFFLINE_CALC_DIFF

C     !INTERFACE:
      SUBROUTINE OFFLINE_CALC_DIFF(
     I        bi, bj, iMin, iMax, jMin, jMax, kArg, kSize,
     U        KappaRx,
     I        myThid )

C     !DESCRIPTION:
C     *===========================================================*
C     | SUBROUTINE OFFLINE_CALC_DIFF                              |
C     | o Add contribution to net diffusivity read in offline pkg |
C     *===========================================================*

C     !USES:
      IMPLICIT NONE
C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
c#include "DYNVARS.h"
c#include "GRID.h"
#include "OFFLINE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     bi, bj    :: current tile indices
C     iMin,iMax :: Range of points for which calculation is done
C     jMin,jMax :: Range of points for which calculation is done
C     kArg      :: = 0 -> do the k-loop here and treat all levels
C                  > 0 -> k-loop is done outside and treat only level k=kArg
C     kSize     :: 3rd Dimension of the vertical diffusivity array KappaRx
C     KappaRx   :: vertical diffusivity array
C     myThid    :: my Thread Id number
      INTEGER bi, bj, iMin, iMax, jMin, jMax, kArg, kSize
      _RL KappaRx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,kSize)
      INTEGER myThid

#ifdef ALLOW_OFFLINE
C     !LOCAL VARIABLES:
C     i, j, k   :: Loop counters
      INTEGER i, j, k
CEOP

C--   Add Vertical diffusivity contribution from OFFLINE
      IF ( offlineDiffKrFile .NE. ' ' ) THEN
       IF ( kArg .EQ. 0 ) THEN
C-     do all levels :
        DO k=1,MIN(Nr,kSize)
         DO j=jMin,jMax
          DO i=iMin,iMax
            KappaRx(i,j,k) = KappaRx(i,j,k)
     &                  +( offlineDiffKr(i,j,k,bi,bj)
     &                     - diffKrNrS(k) )
          ENDDO
         ENDDO
        ENDDO
       ELSE
C-     do level k=kArg only :
        k = MIN(kArg,kSize)
        DO j=jMin,jMax
         DO i=iMin,iMax
           KappaRx(i,j,k) = KappaRx(i,j,k)
     &                 +( offlineDiffKr(i,j,kArg,bi,bj)
     &                    - diffKrNrS(kArg) )
         ENDDO
        ENDDO
       ENDIF
      ENDIF

#endif /* ALLOW_OFFLINE */

      RETURN
      END
