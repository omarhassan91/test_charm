#if ENBONDA_EXPAND == 1

!-------------------------------------------------------------------
! Do VDWOPT expansion of code
! ##EXPAND  V dswit dfswit dshft dvips gvswit .when. EXPAND (expand_vdwopt)
! ##PASS1   dswit
  IF (DSWIT) THEN

#define ENBONDA_DSWIT 1
#include "enbonda_expand.inc"
#undef ENBONDA_DSWIT
     
! ##PASS2   dfswit
  ELSE IF (DFSWIT) THEN

#define ENBONDA_DFSWIT 1
#include "enbonda_expand.inc"
#undef ENBONDA_DFSWIT

! ##PASS3   dshft
  ELSE IF (DSHFT) THEN

#define ENBONDA_DSHFT 1
#include "enbonda_expand.inc"
#undef ENBONDA_DSHFT

#if KEY_NBIPS == 1
! ##PASS4   dvips   .when. NBIPS
  ELSE IF (DVIPS) THEN

#define ENBONDA_DVIPS 1
#include "enbonda_expand.inc"
#undef ENBONDA_DVIPS

#endif /* KEY_NBIPS */

! ##PASSA   gvswit
  ELSE IF (GVSWIT) THEN

#define ENBONDA_GVSWIT 1
#include "enbonda_expand.inc"
#undef ENBONDA_GVSWIT

! ##PASS5
  ELSE

#include "enbonda_expand.inc"

! ##EXFIN
  ENDIF
! ##EXEND
! ##ENDEX    (expand_vdwopt)

#else /* ENBONDA_EXPAND */

#define ENBONDA_V_FLAG 1
#define ENBONDA_DSWIT 1
#define ENBONDA_DFSWIT 1
#define ENBONDA_DSHFT 1
#define ENBONDA_DVIPS 1
#define ENBONDA_GVSWIT 1
  
#include "enbonda_expand.inc"

#undef ENBONDA_V_FLAG
#undef ENBONDA_DSWIT
#undef ENBONDA_DFSWIT
#undef ENBONDA_DSHFT
#undef ENBONDA_DVIPS
#undef ENBONDA_GVSWIT
  
#endif /* ENBONDA_EXPAND */
