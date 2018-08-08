MODULE DGElementClass
  USE NodalDGStorageClass
  USE FluxRoutines

  IMPLICIT NONE

  TYPE DGElement
     REAL(KIND=RP)                                :: delta_x,xL,xR,yL,yR,zL,zR
     REAL(KIND=RP),DIMENSION(6)                   :: lambdamax
     INTEGER                                      :: nEqn
     REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:,:)   :: FstarR,FstarL&
          &,GstarR,GstarL,HstarR,HstarL,QLx,QRx,QLy,QRy,QRz,QLz
     REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:,:,:) :: Q,Q_dot,res
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
    ALLOCATE(this%Q(0:n,0:n,0:n,nEqn),this%Q_dot(0:n,0:n,0:n,nEqn)&
         &,this%QLx(0:n,0:n,nEqn),this%QLy(0:n,0:n,nEqn),this&
         &%QLz(0:n,0:n,nEqn),this%QRx(0:n,0:n,nEqn),this%QRy(0:n,0:n&
         &,nEqn),this%QRz(0:n,0:n,nEqn))
    ALLOCATE(this%res(0:N,0:N,0:N,nEqn))
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

  SUBROUTINE LocalTimeDerivative(this,DG,t)
    IMPLICIT NONE
    TYPE(DGElement)      ,INTENT(INOUT)  :: this
    TYPE(NodalDGStorage) ,INTENT(IN)     :: DG
    REAL(KIND=RP)        ,INTENT(IN)     :: t
    
    ! Local Variables

    INTEGER                                   :: j,N,nEqn
    REAL(KIND=RP),DIMENSION(0:DG%N,0:DG%N,0:DG%N,this%nEqn) ::F_prime,G_prime&
         &,H_prime
    
    nEqn = this%nEqn
    N    = DG%N
    
    CALL SystemDGDerivative(this,this%FstarR,this%FstarL,this%GstarR,this&
         &%GstarL,this%HstarR,this%HstarL,F_prime,G_prime&
         &,H_prime,DG%D,DG%weights,DG%l_at_one,DG%l_at_minus_one,this%Q&
         &,nEqn,N)
    this%Q_dot=(-8.0_RP/this%delta_x**3)*(F_prime+G_prime+H_prime)+this%res
       
  END SUBROUTINE LocalTimeDerivative
  
  !////////////////////////


  SUBROUTINE SystemDGDerivative(this,FR,FL,GR,GL,HR,HL,Fprime,Gprime&
       &,Hprime,D,weights,l_one,l_minus_one,Q,nEqn,N)
    IMPLICIT NONE
    TYPE(DGElement)                          ,INTENT(IN) :: this
    INTEGER                                  ,INTENT(IN) :: nEqn,N
    REAL(KIND=RP),DIMENSION(0:N,0:N,nEqn)    ,INTENT(IN) :: FR,FL,GL,GR,HL,HR
    REAL(KIND=RP),DIMENSION(0:N,0:N,0:N,nEqn),INTENT(IN) :: Q
    REAL(KIND=RP),DIMENSION(0:N,0:N,0:N,nEqn),INTENT(OUT):: Fprime,Gprime,Hprime
    REAL(KIND=RP),DIMENSION(0:N,0:N)         ,INTENT(IN) :: D
    REAL(KIND=RP),DIMENSION(0:N)             ,INTENT(IN) :: weights,l_one,l_minus_one

    !Local Variables

    INTEGER :: i,j,k,l
    REAL(KIND=RP),DIMENSION(0:N,0:N,nEqn) :: FanaL,FanaR,GanaL,GanaR&
         &,HanaL,HanaR
    !initializing output variables
    Fprime = 0.0_RP
    Gprime = 0.0_RP
    Hprime = 0.0_RP
    !compute volume Terms
    CALL computeVolumePI(D,Q,1,Fprime,N,nEqn)
    CALL computeVolumePI(D,Q,2,Gprime,N,nEqn)
    CALL computeVolumePI(D,Q,3,Hprime,N,nEqn)
    !compute surface Terms
    CALL EulerAnalyticFlux(this%QLx,FanaL,1,N,nEqn)
    CALL EulerAnalyticFlux(this%QRx,FanaR,1,N,nEqn)
    CALL EulerAnalyticFlux(this%QLy,GanaL,2,N,nEqn)
    CALL EulerAnalyticFlux(this%QRy,GanaR,2,N,nEqn)
    CALL EulerAnalyticFlux(this%QLz,HanaL,3,N,nEqn)
    CALL EulerAnalyticFlux(this%QRz,HanaR,3,N,nEqn)
    DO l=0,N
       DO k=0,N
          DO i=1,nEqn
             Fprime(0,k,l,i)=Fprime(0,k,l,i)+(FL(k,l,i)+FanaL(k,l,i))
             Fprime(N,k,l,i)=Fprime(N,k,l,i)+(FR(k,l,i)-FanaR(k,l,i))
             Gprime(k,0,l,i)=Gprime(k,0,l,i)+(GL(k,l,i)+GanaL(k,l,i))
             Gprime(k,N,l,i)=Gprime(k,N,l,i)+(GR(k,l,i)-GanaR(k,l,i))
             Hprime(k,l,0,i)=Hprime(k,l,0,i)+(HL(k,l,i)+HanaL(k,l,i))
             Hprime(k,l,N,i)=Hprime(k,l,N,i)+(HR(k,l,i)-HanaR(k,l,i))
          ENDDO
       ENDDO
    ENDDO
    Fprime=Fprime*this%delta_x**2/4.0_RP
    Gprime=Gprime*this%delta_x**2/4.0_RP
    Hprime=Hprime*this%delta_x**2/4.0_RP
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
         &,this%QLy,this%QRy,this%QRz,this%QLz,this%res)

  END SUBROUTINE DestructElement
  
END MODULE DGELEMENTCLASS
