MODULE TimeIntegration
  USE DGMeshClass
  IMPLICIT NONE
  REAL(KIND=RP) :: CFL,lambdaMaxtime
  REAL(KIND=RP),DIMENSION(5) :: a=(/0.0_rp, -567301805773.0_rp/1357537059087.0_rp,&
      -2404267990393.0_rp/2016746695238.0_rp, -3550918686646.0_rp/2091501179385.0_rp,&
      -1275806237668.0_rp/842570457699.0_rp/),b=(/0.0_rp, 1432997174477.0_rp/9575080441755.0_rp,&
      2526269341429.0_rp/6820363962896.0_rp, 2006345519317.0_rp/3224310063776.0_rp,&
2802321613138.0_rp/2924317926251.0_rp /),c=(/1432997174477.0_rp/9575080441755.0_rp, 5161836677717.0_rp/13612068292357.0_rp,&
      1720146321549.0_rp/2090206949498.0_rp, 3134564353537.0_rp/4481467310338.0_rp,&
2277821191437.0_rp/14882151754819.0_rp /)
CONTAINS

  SUBROUTINE ConstructSimulation(this,NQ,N,nEqn,xmin,xmax,ymin,ymax&
       &,zmin,zmax,CFLVal)
    IMPLICIT NONE
    TYPE(DGMesh) ,INTENT(INOUT) :: this
    INTEGER      ,INTENT(IN)    :: NQ,N,nEqn
    REAL(KIND=RP),INTENT(IN)    :: xmin,xmax,ymin,ymax,zmin,zmax,CFLVal

    !local variables

    REAL(KIND=RP)                 :: dx
    INTEGER                       :: i,K
    REAL(KIND=RP),DIMENSION(0:NQ) :: nodes

    K=NQ**3
    dx=(xmax-xmin)/real(NQ,kind=RP)
    DO i=0,NQ
       nodes(i)=i*dx+xmin
    END DO
    CFL=CFLVal
    CALL ConstructMesh3D(this,nodes,K,N,nEqn,0.0_RP)
  END SUBROUTINE ConstructSimulation

  !///////////////////////////////////////////////////////////

  SUBROUTINE getLambdaMaxglobally(this)
    IMPLICIT NONE

    TYPE(DGMesh),INTENT(INOUT) :: this

    !local
    INTEGER                                                        :: i
    REAL(KIND=RP),DIMENSION(this%K,0:this%DG%N,0:this%DG%N,0:this%DG%N) ::&
         & p,c
    REAL(KIND=RP),DIMENSION(this%K)    :: maxatElement                                

    DO i=1,this%K
       p(i,:,:,:)=(gamma-1.0_RP)*(this%e(i)%Q_dot(:,:,:,5)-0.5_RP&
            &*(this%e(i)%Q_dot(:,:,:,2)**2+this%e(i)%Q_dot(:,:,:,3)**2&
            &+this%e(i)%Q_dot(:,:,:,4)**2)/this%e(i)%Q_dot(:,:,:,1)
    END DO
    DO i=1,this%K
       IF(any(gamma*p(i,:,:,:)/this%e(i)%Q_dot(:,:,:,1)<=0)) THEN
          print*,"pressure/density negativ in element ",i
          exit(1)
       ELSE
          c(i,:,:,:)=sqrt(gamma*p(i,:,:,:)/this%e(i)%Q_dot(:,:,:,1))
       ENDIF
    END DO
    DO i=1,this%K
       maxatElement(i)=max(maxval(abs(this%e(i)%Q_dot(:,:,:,2)/this&
            &%e(i)%Q_dot(:,:,:,1)+c)),maxval(abs(this%e(i)%Q_dot(:,:,:&
            &,3)/this%e(i)%Q_dot(:,:,:,1)+c)),maxval(abs(this%e(i)&
            &%Q_dot(:,:,:,4)/this%e(i)%Q_dot(:,:,:,1)+c)),maxval(abs(this%e(i)%Q_dot(:,:,:,2)/this&
            &%e(i)%Q_dot(:,:,:,1)-c)),maxval(abs(this%e(i)%Q_dot(:,:,:&
            &,3)/this%e(i)%Q_dot(:,:,:,1)-c)),maxval(abs(this%e(i)&
            &%Q_dot(:,:,:,4)/this%e(i)%Q_dot(:,:,:,1)-c))
       this%e(i)%lambdamax=maxatElement(i)
    END DO
    this%lambdamax=maxval(maxatElement)
  END SUBROUTINE getLambdaMaxglobally

  !///////////////////////////////////////////////////////

  SUBROUTINE getRungeKuttaStep(this,t)
    IMPLICIT NONE

    TYPE(DGMesh) ,INTENT(INOUT) :: this
    REAL(KIND=RP),INTENT(IN) :: t

    !local
    INTEGER :: step,i
    REAL(KIND=RP) :: dt
    TYPE(DGMesh) :: previous
    REAL(KIND=RP),DIMENSION(1:this%K,0:this%DG%N,0:this%DG%N,0:this&
         &%DG%N,this%e(1)%nEqn) :: g
    CALL getLambdaMaxglobally(this)
    
    dt=CFL/(3.0_RP*this%lambdamax)*this%e(1)%delta_x/(2*real(this%DG&
         &%N,kind=rp)+1)
    previous=this
    CALL GlobalTimeDerivative(this,t)
    
    DO i=1,K
       g(i,:,:,:,:)=this%e(i)%Q_dot
    END DO
    DO step=1,5
       g=a(step)*g
       CALL GlobalTimeDerivative(this,t+b(step)*dt)
       DO i=1,K
          g(i,:,:,:,:)=g(i,:,:,:,:)+this%e(i)%Q_dot
       END DO
       DO i=1,K
          previous%e(i)%Q_dot=previous&e(i)%Q_dot+c(step)*dt*g(i,:,:,:&
               &,:)
       END DO
       this=previous
    END DO
  END SUBROUTINE getRungeKuttaStep
END MODULE TimeIntegration

    
    
       

    
  
       

    

    
    
    
  
