project(DGEuler)
cmake_minimum_required(VERSION 3.10)
enable_language(Fortran)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG} -fbounds-check -ffpe-trap=invalid,denormal,zero -fbacktrace")
add_executable(DGEuler Driver.f90 TimeIntegration.f90 DGMeshClass.f90 DGElementClass.f90 NodalDGStorage.f90 FluxRoutines.f90 DGtoolbox.f90 Diverses.f90)

install(TARGETS DGEuler DESTINATION bin)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/bin/DGEuler.exe DESTINATION /Build)
