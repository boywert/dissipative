#include "../../build/arepoconfig.h"
#ifdef SGCHEM
      SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
      IMPLICIT NONE
#include "cool.h"
      integer NEQ, IPAR
      REAL T, Y(NEQ), RPAR
      integer ML, MU, NROWPD
      REAL PD(NROWPD,NEQ)
      return
      end
#endif
