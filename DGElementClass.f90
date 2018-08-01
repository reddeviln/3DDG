MODULE DGElementClass
  USE NodalDGStorageClass

  IMPLICIT NONE

  TYPE DGElement
     REAL(KIND=RP)                                :: delta_x,xL,xR,yL,yR,zL,zR,lambdamax
     INTEGER                                      :: nEqn
     REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:,:)   :: FstarR,FstarL&
          &,GstarR,GstarL,HstarR,HstarL,QLx,QRx,QLy,QRy,QRz,QLz
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
    this%QLx        = 0.0_RP
    this%QRx        = 0.0_RP
    this%QLy        = 0.0_RP
    this%QRy        = 0.0_RP
    this%QLz        = 0.0_RP
    this%QRz        = 0.0_RP
    this%Q_dot      = 0.0_RP
    this%lambdamax  = 0.0_RP
  END SUBROUTINE ConstructDGElement
  
!////////////////////////////////////////////////////////////////////

  SUBROUTINE LocalTimeDerivative(this,DG)
    IMPLICIT NONE
    TYPE(DGElement)      ,INTENT(INOUT)  :: this
    TYPE(NodalDGStorage) ,INTENT(IN)     :: DG
    
    ! Local Variables

    INTEGER                                   :: j,N,nEqn
    REAL(KIND=RP),DIMENSION(0:DG%N,0:DG%N,0:DG%N,this%nEqn) ::F_prime,G_prime&
         &,H_prime
    
    nEqn = this%nEqn
    N    = DG%N
    
    
    CALL SystemDGDerivative(this%FstarR,this%FstarL,this%GstarR,this&
         &%GstarL,this%HstarR,this%HstarL,F_prime,G_prime&
         &,H_prime,DG%D,DG%weights,DG%l_at_one,DG%l_at_minus_one,this%Q&
         &,nEqn,N)
    this%Q_dot=(-8.0_RP/this%delta_x**3)*(F_prime+G_prime+H_prime)

  END SUBROUTINE LocalTimeDerivative
  
!////////////////////////

  SUBROUTINE SystemDGDerivative(FR,FL,GR,GL,HR,HL,Fprime,Gprime&
       &,Hprime,D,weights,l_one,l_minus_one,Q,nEqn,N)
    IMPLICIT NONE
    INTEGER                                  ,INTENT(IN) :: nEqn,N
    REAL(KIND=RP),DIMENSION(0:N,0:N,nEqn)    ,INTENT(IN) :: FR,FL,GL,GR,HL,HR
    REAL(KIND=RP),DIMENSION(0:N,0:N,0:N,nEqn),INTENT(IN) :: Q
    REAL(KIND=RP),DIMENSION(0:N,0:N,0:N,nEqn),INTENT(OUT):: Fprime,Gprime,Hprime
    REAL(KIND=RP),DIMENSION(0:N,0:N)         ,INTENT(IN) :: D
    REAL(KIND=RP),DIMENSION(0:N)             ,INTENT(IN) :: weights,l_one,l_minus_one

    !Local Variables

    INTEGER :: i,j,k,l
    !initializing output variables
    Fprime = 0.0_RP
    Gprime = 0.0_RP
    Hprime = 0.0_RP

    !compute volume Terms

    CALL computeVolumePI(D,Q,1,Fprime,N,nEqn)
    CALL computeVolumePI(D,Q,2,Gprime,N,nEqn)
    CALL computeVolumePI(D,Q,3,Hprime,N,nEqn)

    !compute surface Terms

    DO l=0,N
       DO k=0,N
          DO j=0,N
             DO i=1,nEqn
                Fprime(l,k,j,i)=Fprime(l,k,j,i)+(FR(k,j,i)*l_one(l)+&
                     FL(k,j,i)*l_minus_one(l))/weights(l)
                Gprime(l,k,j,i)=Gprime(l,k,j,i)+(GR(l,j,i)*l_one(k)+&
                     GL(l,j,i)*l_minus_one(k))/weights(k)
                Hprime(l,k,j,i)=Hprime(l,k,j,i)+(HR(l,k,i)*l_one(j)+&
                     HL(l,k,i)*l_minus_one(j))/weights(j)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE SystemDGDerivative

  !/////////////////////

  SUBROUTINE DestructElement(this)
    IMPLICIT NONE
    TYPE(DGElement),INTENT(INOUT) :: this

    this%nEqn       = -1
    this%xL         = 0.0_RP
    this%xR         = 0.0_RP
    this%yL         = 0.0_RP
    this%yR         = 0.0_RP
    this%zL         = 0.0_RP
    this%zR         = 0.0_RP
    this%delta_x    = 0.0_RP
    DEALLOCATE(this%FstarR,this%FstarL,this%GstarR,this%GstarL,this&
         &%HstarR,this%HstarL,this%Q,this%Q_dot,this%QLx,this%QRx&
         &,this%QLy,this%QRy,this%QRz,this%QLz)

  END SUBROUTINE DestructElement
  
END MODULE DGELEMENTCLASS
