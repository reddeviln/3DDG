PROGRAM Driver
  USE TimeIntegration
  IMPLICIT NONE
  REAL(KIND=RP) :: muVal=0.2_RP,gammaVal=1.4_RP,PrVal=0.7_RP,CFLVal=0.5_RP,t=0.0_RP,tend=1.0_RP
  REAL(KIND=RP) :: xmin=-1.0_RP,xmax=1.0_RP,ymin=-1.0_RP,ymax=1.0_RP&
       &,zmin=-1.0_RP,zmax=1.0_RP
  TYPE(DGMesh)  :: Simulation
  INTEGER,PARAMETER       :: NQ=3,N=3,nEqn=5
  INTEGER       :: i,m=0,j=0
  CHARACTER(len=3) :: numChar
  CHARACTER(len=21) :: fname ='Plots/Movies/UXXX.tec'
  REAL(KIND=RP),DIMENSION(1:nq**3,0:N,0:N,0:N) :: Qplot
  print*,"THIS IS A DGCODE BY NILS DORNBUSCH inspired by 'FORTRAN NOTES'&
       &"," by ANDREW WINTERS"
  print*,"Starting everything up..."
  CALL setPhysPar(muVal,gammaVal,PrVal)
  print*,"physical constants set successfully"
  CALL ConstructSimulation(Simulation,NQ,N,nEqn,xmin,xmax,ymin,ymax&
       &,zmin,zmax,CFLVal)
  print*,"finished setting up. Computing LambdaMax now"
  DO WHILE (tend-t>epsilon(t))
     CALL getLambdaMaxGlobally(Simulation)
     print*,"t"
     print*,t
     print*,"sum Energy 1st element"
     print*,sum(Simulation%e(1)%Q(:,:,:,5))
     CALL getRungeKuttaStep(Simulation,t,tend)
     IF(MODULO(j,5).EQ.0) THEN
        m=m+1
        WRITE(numChar,'(i3)')m
        IF (m.GE.100) THEN
           fName(15:17) = numChar
        ELSEif(m.GE.10) then
           fName(15:15)    = "0"
           fName(16:17)  = numChar(2:3)
        ELSE
           fName(15:16) = "00"
           fName(17:17)=numChar(3:3)
        END IF
        CALL preparePlot(Simulation,fName,1)
     END IF
     j=j+1
  END DO
  CALL preparePlot(Simulation,"Plots/endTime.tec",1)
  print*,"Done without a problem! See you next time."
END PROGRAM Driver

  
