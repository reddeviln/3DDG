MODULE TimeIntegration
  USE DGMeshClass
  IMPLICIT NONE
  REAL(KIND=RP),DIMENSION(:,:,:,:,:),ALLOCATABLE :: points
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
    print*,"Constructing Mesh..."
    CALL ConstructMesh3D(this,nodes,NQ,N,nEqn,0.0_RP)
    print*,"Imposing initial condition ..."
    CALL InitialCondition(this)
  END SUBROUTINE ConstructSimulation

  !///////////////////////////////////////////////////////////

  SUBROUTINE InitialCondition(this)
    IMPLICIT NONE
    TYPE(DGMesh),INTENT(INOUT) :: this
    !localVariables
    REAL(KIND=RP),DIMENSION(0:this%DG%N) :: x,w
    REAL(KIND=RP),DIMENSION(1:this%K,0:this%DG%N,0:this%DG%N,0:this&
         &%DG%N) :: Qplot
    INTEGER :: i,j,l
       !u(:,:,:,:,:,:,1)=2.0_RP+SIN(pi*(xyz(:,:,:,:,:,:,1)+xyz(:,:,:,:,:,:,2)+xyz(:,:,:,:,:,:,3)))/10.0_RP
   !u(:,:,:,:,:,:,2)=u(:,:,:,:,:,:,1)
   !u(:,:,:,:,:,:,3)=u(:,:,:,:,:,:,1)
   !u(:,:,:,:,:,:,4)=u(:,:,:,:,:,:,1)
    !u(:,:,:,:,:,:,5)=u(:,:,:,:,:,:,1)*u(:,:,:,:,:,:,1)
    ALLOCATE(points(1:this%K,0:this%DG%N,0:this%DG%N,0:this%DG%N,3))
    CALL LegendreGaussLobattoNodesandWeights(this%DG%N,x,w)
    DO i=1,this%K
       DO j=0,this%DG%N
          DO l=0,this%DG%N
             points(i,:,j,l,1)=(this%e(i)%xR-this%e(i)%xL)*0.5_RP*x+(this&
                  &%e(i)%xR+this%e(i)%xL)*0.5_RP
          END DO
       END DO
       DO j=0,this%DG%N
          DO l=0,this%DG%N
             points(i,j,:,l,2)=(this%e(i)%yR-this%e(i)%yL)*0.5_RP*x+(this&
                  &%e(i)%yR+this%e(i)%yL)*0.5_RP
          END DO
       END DO
       DO j=0,this%DG%N
          DO l=0,this%DG%N
             points(i,j,l,:,3)=(this%e(i)%zR-this%e(i)%zL)*0.5_RP*x+(this&
                  &%e(i)%zR+this%e(i)%zL)*0.5_RP
          END DO
       END DO
       this%e(i)%Q(:,:,:,1)=2.0_RP+SIN(pi*(points(i,:,:,:,1)+points(i,:,:&
            &,:,2)+points(i,:,:,:,3)))/10.0_RP
       this%e(i)%Q(:,:,:,2)=this%e(i)%Q(:,:,:,1)
       this%e(i)%Q(:,:,:,3)=this%e(i)%Q(:,:,:,1)
       this%e(i)%Q(:,:,:,4)=this%e(i)%Q(:,:,:,1)
       this%e(i)%Q(:,:,:,5)=this%e(i)%Q(:,:,:,1)*this%e(i)%Q(:,:,:,1)
    END DO
    DO i=1,this%K
       Qplot(i,:,:,:)=this%e(i)%Q(:,:,:,1)
    END DO
    
    OPEN(file='Plots/initial.tec',unit=15)
    CALL ExportToTecplot_3D(points(:,:,:,:,1),points(:,:,:,:,2)&
         &,points(:,:,:,:,3),Qplot,this%DG%N,this%K&
         &,15,"rho")
    CLOSE(15)
  END SUBROUTINE InitialCondition
  
!//////////////////////////////////////////////////////////////
  
  SUBROUTINE getLambdaMaxGlobally(this)
    IMPLICIT NONE

    TYPE(DGMesh),INTENT(INOUT) :: this

    !local
    INTEGER                                                        :: i
    REAL(KIND=RP),DIMENSION(this%K,0:this%DG%N,0:this%DG%N,0:this%DG%N) ::&
         & p,c
    REAL(KIND=RP),DIMENSION(this%K)    :: maxatElement                                

    DO i=1,this%K
       p(i,:,:,:)=(gamma-1.0_RP)*(this%e(i)%Q(:,:,:,5)-0.5_RP&
            &*(this%e(i)%Q(:,:,:,2)**2+this%e(i)%Q(:,:,:,3)**2&
            &+this%e(i)%Q(:,:,:,4)**2)/this%e(i)%Q(:,:,:,1))
    END DO
    DO i=1,this%K
       IF(any(gamma*p(i,:,:,:)/this%e(i)%Q(:,:,:,1)<=0)) THEN
          print*,"pressure/density negativ in element ",i,minval(p(i&
               &,:,:,:)),minval(this%e(i)%Q(:,:,:,1))
          CALL EXIT(1)
       ELSE
          c(i,:,:,:)=sqrt(gamma*p(i,:,:,:)/this%e(i)%Q(:,:,:,1))
       ENDIF
    END DO
    DO i=1,this%K
       maxatElement(i)=max(maxval(abs(this%e(i)%Q(:,:,:,2)/this&
            &%e(i)%Q(:,:,:,1)+c(i,:,:,:))),maxval(abs(this%e(i)%Q(:,:,:&
            &,3)/this%e(i)%Q(:,:,:,1)+c(i,:,:,:))),maxval(abs(this%e(i)&
            &%Q(:,:,:,4)/this%e(i)%Q(:,:,:,1)+c(i,:,:,:))),maxval(abs(this%e(i)%Q(:,:,:,2)/this&
            &%e(i)%Q(:,:,:,1)-c(i,:,:,:))),maxval(abs(this%e(i)%Q(:,:,:&
            &,3)/this%e(i)%Q(:,:,:,1)-c(i,:,:,:))),maxval(abs(this%e(i)&
            &%Q(:,:,:,4)/this%e(i)%Q(:,:,:,1)-c(i,:,:,:))))
       this%e(i)%lambdamax=maxatElement(i)
    END DO
    this%lambdamax=maxval(maxatElement)
  END SUBROUTINE getLambdaMaxglobally

  !///////////////////////////////////////////////////////

  SUBROUTINE getEulerResidual(this,t)
    IMPLICIT NONE
    TYPE(DGMesh) :: this
    REAL(KIND=RP)   :: t
    !local variables
    REAL(KIND=RP) :: c1,c2,c3,c4,c5
    INTEGER :: i,j,l,m
    c1=pi/10.0_RP
    c2=-0.2_RP*pi+0.05_RP*pi*(1.0_RP+5.0_RP*gamma)
    c3=pi/100.0_RP*(gamma-1.0_RP)
    c4=0.05_RP*(-16.0_RP*pi+pi*(9.0_RP+15.0_RP*gamma))
    c5=0.01_RP*(3.0_RP*pi*gamma-2.0_RP*pi)
    DO i=1,this%K
       DO m=0,this%dg%n
          DO j=0,this%dg%n
             DO l=0,this%dg%n
                this%e(i)%res(l,j,m,1)=c1*cos(pi*(points(i,l,j,m,1)&
                     &+points(i,l,j,m,2)+points(i,l,j,m,3)-2.0_RP*t))
                this%e(i)%res(l,j,m,2)=c2*cos(pi*(points(i,l,j,m,1)&
                     &+points(i,l,j,m,2)+points(i,l,j,m,3)-2.0_RP*t))&
                     &+c3*cos(2.0_RP*pi*(points(i,l,j,m,1)&
                     &+points(i,l,j,m,2)+points(i,l,j,m,3)-2.0_RP*t))
                this%e(i)%res(l,j,m,3)=this%e(i)%res(l,j,m,2)
                this%e(i)%res(l,j,m,4)=this%e(i)%res(l,j,m,2)
                this%e(i)%res(l,j,m,5)=c4*cos(pi*(points(i,l,j,m,1)&
                     &+points(i,l,j,m,2)+points(i,l,j,m,3)-2.0_RP*t))&
                     &+c5*cos(2.0_RP*pi*(points(i,l,j,m,1)&
                     &+points(i,l,j,m,2)+points(i,l,j,m,3)-2.0_RP*t))
             END DO
          END DO
       END DO
    END DO
  END SUBROUTINE getEulerResidual
  
 !//////////////////////////////////////////

  SUBROUTINE getRungeKuttaStep(this,t,tend)
    IMPLICIT NONE

    TYPE(DGMesh) ,INTENT(INOUT) :: this
    REAL(KIND=RP),INTENT(INOUT) :: t
    REAL(KIND=RP),INTENT(IN)    :: tend

    !local
    INTEGER :: step,i
    REAL(KIND=RP) :: dt
    TYPE(DGMesh) :: previous
    REAL(KIND=RP),DIMENSION(1:this%K,0:this%DG%N,0:this%DG%N,0:this&
         &%DG%N,this%e(1)%nEqn) :: g    
    dt=CFL/(3.0_RP*this%lambdamax)*this%e(1)%delta_x/(2*real(this%DG&
         &%N,kind=rp)+1)
    dt=min(dt,tend-t)
    print*,"dt"
    print*,dt
    print*,"Constructing global time Derivative..."
    CALL getEulerResidual(this,t)
    CALL GlobalTimeDerivative(this,t)
    DO i=1,this%K
       g(i,:,:,:,:)=this%e(i)%Q_dot
    END DO
    DO step=1,5
       print*,"RK stage" ,step
       g=a(step)*g
       CALL getEulerResidual(this,t+b(step)*dt)
       CALL GlobalTimeDerivative(this,t+b(step)*dt)
       DO i=1,this%K
          g(i,:,:,:,:)=g(i,:,:,:,:)+this%e(i)%Q_dot
       END DO
       DO i=1,this%K
          this%e(i)%Q=this%e(i)%Q+c(step)*dt*g(i,:,:,:&
               &,:)
       END DO
    END DO
    t=t+dt
    print*,"Finished with RK step"
  END SUBROUTINE getRungeKuttaStep

  !///////////////////////////////////////////////////////

  SUBROUTINE preparePlot(this,fname,var)
    IMPLICIT NONE
    TYPE(DGMesh)     ,INTENT(IN) :: this
    CHARACTER(len=*) ,INTENT(IN) :: fname
    INTEGER          ,INTENT(IN) :: var
    !local
    REAL(KIND=RP),DIMENSION(1:this%K,0:this%DG%N,0:this%DG%N,0:this&
         &%DG%N) ::Qplot
    INTEGER :: i
    DO i=1,this%K
       Qplot(i,:,:,:)=this%e(i)%Q(:,:,:,var)
    END DO
    
    OPEN(file=fname,unit=15)
    CALL ExportToTecplot_3D(points(:,:,:,:,1),points(:,:,:,:,2),points(:,:,:,:,3),Qplot,this%DG%N,this%K,15,"var")
    CLOSE(15)
  END SUBROUTINE preparePlot
    
END MODULE TimeIntegration
    

    
    
    
  
