module coordc
  use chm_kinds
  use dimens_fcm
  implicit none
  character(len=*),private,parameter :: file_name   ="coordc_ltm.src"
  !
  !     The Comparison Coordinates
  !
  !     Purpose:
  !     Maintaining an alternate set of coordinates for analysis and
  !     for coordinate manipulations (with normal modes).
  !     Variable  Purpose
  !
  !     XCOMP         X component
  !     YCOMP         Y component
  !     ZCOMP         Z component
  !     WCOMP         Weighting or temp factor array
  !     second comparison set:
  !     XCOMP2         X component
  !     YCOMP2         Y component
  !     ZCOMP2         Z component
  !
  !
  !
  real(chm_real),dimension(:),allocatable,save :: XCOMP,YCOMP,ZCOMP,WCOMP
#if KEY_COMP2==1
  real(chm_real),dimension(:),allocatable,save :: XCOMP2,YCOMP2,ZCOMP2,WCOMP2 
#endif

contains
  subroutine allocate_coordc()
    use memory
    character(len=*),parameter :: routine_name="allocate_coordc"
    call chmalloc(file_name,routine_name,'xcomp ',maxa,crl=xcomp)
    call chmalloc(file_name,routine_name,'ycomp ',maxa,crl=ycomp)
    call chmalloc(file_name,routine_name,'zcomp ',maxa,crl=zcomp)
    call chmalloc(file_name,routine_name,'wcomp ',maxa,crl=wcomp)
#if KEY_COMP2==1
    call chmalloc(file_name,routine_name,'xcomp2',maxa,crl=xcomp2) 
#endif
#if KEY_COMP2==1
    call chmalloc(file_name,routine_name,'ycomp2',maxa,crl=ycomp2) 
#endif
#if KEY_COMP2==1
    call chmalloc(file_name,routine_name,'zcomp2',maxa,crl=zcomp2) 
#endif
#if KEY_COMP2==1
    call chmalloc(file_name,routine_name,'wcomp2',maxa,crl=wcomp2) 
#endif
    return
  end subroutine allocate_coordc

end module coordc

