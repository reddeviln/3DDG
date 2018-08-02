PROGRAM Driver
  USE TimeIntegration
  IMPLICIT NONE
  REAL(KIND=RP) :: muVal=0.2_RP,gammaVal=1.4_RP,PrVal=0.7_RP,CFLVal=0.5_RP,t=0.0_RP,tend=0.5_RP
  REAL(KIND=RP) :: xmin=-1.0_RP,xmax=1.0_RP,ymin=-1.0_RP,ymax=1.0_RP&
       &,zmin=-1.0_RP,zmax=1.0_RP
  TYPE(DGMesh)  :: Simulation
  INTEGER       :: NQ=3,N=3,nEqn=5
  print*,"THIS IS DGCODE BY NILS DORNBUSCH inspired by 'FORTRAN NOTES'&
       &"," by ANDREW WINTERS"
  print*,"Starting everything up..."
  CALL setPhysPar(muVal,gammaVal,PrVal)
  print*,"physical constants set succesfully"
  CALL ConstructSimulation(Simulation,NQ,N,nEqn,xmin,xmax,ymin,ymax&
       &,zmin,zmax,CFLVal)
  print*,"finished setting up. Computing LambdaMax now"
  DO WHILE (tend-t>epsilon(t))
     CALL getLambdaMaxGlobally(Simulation)
     print*,"t"
     print*,t
     print*,"sum Energy 2nd element"
     print*,sum(Simulation%e(2)%Q_dot(:,:,:,5))
     CALL getRungeKuttaStep(Simulation,t,tend)
  END DO
  
  
END PROGRAM Driver

  
