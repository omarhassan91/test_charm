#define PASTE(a) a
#define ADD_BLOCK_SUFFIX(name) PASTE(name)ENBAEXP_BLOCK_SUFFIX
#define ADD_ELECOPT_SUFFIX(name) PASTE(name)ENBAEXP_ELECOPT_SUFFIX
#define ADD_VDWOPT_SUFFIX(name) PASTE(name)ENBAEXP_VDWOPT_SUFFIX

CALL ADD_BLOCK_SUFFIX(ADD_ELECOPT_SUFFIX(ADD_VDWOPT_SUFFIX(ENBA)))(ENB,EEL,ELECFG,FIRST,NATOMX, &
     CGX,JNBL,INBL, &
     CCNBA,CCNBB,CCNBC,CCNBD, &
#if ENBAEXP_E0 == 1
     LEWLD,CSHIFT,CFSWIT,RSHFT,GESWIT,   & ! E0
#endif  /* ENBAEXP_E0 == 1 */
#if ENBAEXP_V0 == 1
     LVSW,LVFSW,LVSH,GVSWIT,             & ! V0
#endif  /* ENBAEXP_V0 == 1 */
#if KEY_NBIPS == 1
     LVIPSX,LEIPSX,                      & ! NBIPS
#endif  /* KEY_NBIPS */
     IACNBX)

! fool the cpp into expanding the function name
! by putting the quotes on separate lines
IF(PRNLEV > 6) WRITE(OUTU,125) &
  " &
  ADD_BLOCK_SUFFIX(ADD_ELECOPT_SUFFIX(ADD_VDWOPT_SUFFIX(ENBA))) &
  "
