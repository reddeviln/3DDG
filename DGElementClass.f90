MODULE DGElementClass
  USE NodalDGStorageClass

  IMPLICIT NONE

  TYPE DGElement
     REAL(KIND=RP)                                :: delta_x,xL,xR,yL,yR,zL,zR
     INTEGER                                      :: nEqn
     REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:,:)       :: FstarR,FstarL,GstarR,GstarL,HstarR,HstarL
     REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:,:,:) :: Q,Q_dot
  END TYPE DGElement

  !
CONTAINS
  !
  SUBROUTINE ConstructDGElement(this,nEqn,xL,xR,yL,yR,zL,zR,N)
    IMPLICIT NONE
    INTEGER        , INTENT(IN)      :: N
    TYPE(DGElement), INTENT(INOUT)   :: this
    REAL(KIND=RP)  , INTENT(IN)      :: xL,xR,yL,yR,zL,zR
    INTEGER        , INTENT(IN)      :: nEqn

    ALLOCATE(this%FstarR(0:N,0:N,nEqn),this%FstarL(0:n,0:n,nEqn),this%GstarR(0:n,0:n,nEqn)&
         ,this%GstarL(0:n,0:n,nEqn),this%HstarR(0:n,0:n,nEqn),this%HstarL(0:n,0:n,nEqn))
    ALLOCATE(this%Q(0:n,0:n,0:n,nEqn),this%Q_dot(0:n,0:n,0:n,nEqn))
    this%nEqn       = nEqn
    this%xL         = xL
    this%xR         = xR
    this%yL         = yL
    this%yR         = yR
    this%zL         = zL
    this%zR         = zR
    this%delta_x    = xR-xL
    this%FstarL     = 0.0_RP
    this%FstarR     = 0.0_RP
    this%GstarL     = 0.0_RP
    this%GstarR     = 0.0_RP
    this%HstarL     = 0.0_RP
    this%HstarR     = 0.0_RP
    this%Q          = 0.0_RP
    this%Q_dot      = 0.0_RP

  end SUBROUTINE ConstructDGElement

end MODULE DGELEMENTCLASS
