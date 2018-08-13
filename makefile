DGEulerDebug : Diverses.f90 DGtoolbox.f90 FluxRoutines.f90 NodalDGStorage.f90 DGElementClass.f90 DGMeshClass.f90 TimeIntegration.f90 Driver.f90
		gfortran -g -fbacktrace -o DGEulerDebug -fbounds-check -fcheck=all -ffpe-trap=invalid,denormal,zero,underflow,overflow Diverses.f90 DGtoolbox.f90 FluxRoutines.f90 NodalDGStorage.f90 DGElementClass.f90 DGMeshClass.f90 TimeIntegration.f90 Driver.f90
DGEuler : Diverses.f90 DGtoolbox.f90 FluxRoutines.f90 NodalDGStorage.f90 DGElementClass.f90 DGMeshClass.f90 TimeIntegration.f90 Driver.f90
		gfortran -O3 -o DGEuler Diverses.f90 DGtoolbox.f90 FluxRoutines.f90 NodalDGStorage.f90 DGElementClass.f90 DGMeshClass.f90 TimeIntegration.f90 Driver.f90
