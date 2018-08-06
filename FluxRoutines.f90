MODULE FluxRoutines
  USE Diverses
  IMPLICIT NONE
  !Contains many different flux functions and Riemann solvers
  REAL(KIND=RP)         :: mu,g=9.812_RP,R=287.058_RP,gamma,Pr

CONTAINS

  SUBROUTINE setPhysPar(muVal,gammaVal,PrVal)
    IMPLICIT NONE
    REAL(KIND=RP),INTENT(IN) ::muVal,gammaVal,PrVal

    mu=muVal
    gamma=gammaVal
    Pr=PrVal
  END SUBROUTINE setPhysPar

  !//////////////////////////////////////////////////////
  
  SUBROUTINE computeVolumePI(D,Q,dir,Fout,N,nEqn)
    IMPLICIT NONE
    REAL(KIND=RP),DIMENSION(0:N,0:N)         ,INTENT(IN) :: D
    REAL(KIND=RP),DIMENSION(0:N,0:N,0:N,nEqn),INTENT(IN) :: Q
    INTEGER                                  ,INTENT(IN) :: dir,N,nEqn
    REAL(KIND=RP),DIMENSION(0:N,0:N,0:N,nEqn),INTENT(OUT):: Fout

    !Local Variables
    INTEGER                                   :: i,j,k,l
    REAL(KIND=RP),DIMENSION(nEqn,0:N,0:N,0:N) :: Foutc
    REAL(KIND=RP),DIMENSION(0:N,nEqn)         :: Fsharp
    SELECT CASE (dir)
    CASE(1)
       DO k=0,N
          DO j=0,N
             DO i=0,N
                CALL computeEulerPI(Q(i,j,k,:),Q(:,j,k,:),N,nEqn&
                     &,dir,Fsharp)                
                DO l=1,nEqn
                   Foutc(l,i,j,k)=2.0_RP*dot_product(D(i,:),Fsharp(:,l))
                   
                END DO
                
             END DO
          END DO   
       END DO
       
    CASE(2)   
       DO k=0,N
          DO j=0,N
             DO i=0,N
                CALL computeEulerPI(Q(i,j,k,:),Q(i,:,k,:),N,nEqn&
                     &,dir,Fsharp)
                DO l=1,nEqn
                   Foutc(l,i,j,k)=2.0_RP*dot_product(D(j,:),Fsharp(:,l))
                END DO
             END DO
          END DO
       END DO
    CASE(3)
       
       DO k=0,N
          DO j=0,N
             DO i=0,N
                CALL computeEulerPI(Q(i,j,k,:),Q(i,j,:,:),N,nEqn&
                     &,dir,Fsharp)
                DO l=1,nEqn
                   Foutc(l,i,j,k)=2.0_RP*dot_product(D(k,:),Fsharp(:,l))
                END DO
             END DO
          END DO
       END DO
    CASE DEFAULT
       print*,"specify 1,2 or 3 for flux direction"
    END SELECT
    
    DO i=1,nEqn
       Fout(:,:,:,i)=Foutc(i,:,:,:)
    END DO
    
    
  END SUBROUTINE computeVolumePI

  !////////////////////////////////////////////////////

  SUBROUTINE computeEulerPI(QL,QR,N,nEqn,dir,Fsharp)
    IMPLICIT NONE

    INTEGER                           ,INTENT(IN) :: N, nEqn,dir
    REAL(KIND=RP),DIMENSION(nEqn)     ,INTENT(IN) :: QL
    REAL(KIND=RP),DIMENSION(0:N,nEqn) ,INTENT(IN) :: QR
    REAL(KIND=RP),DIMENSION(0:N,nEqn) ,INTENT(OUT):: Fsharp

    !Local Variables

    INTEGER                            :: i,j
    REAL(KIND=RP)                      :: p1,h1
    REAL(KIND=RP),DIMENSION(0:N)       :: p2,h2
    
    p1=(gamma-1.0_RP)*(QL(5)-0.5_RP*(QL(2)*QL(2)+QL(3)*QL(3)+QL(4)*QL(4))/QL(1))
    p2=(gamma-1.0_RP)*(QR(:,5)-0.5_RP*(QR(:,2)*QR(:,2)+QR(:,3)*QR(:&
         &,3)+QR(:,4)*QR(:,4))/QR(:,1))
    h1=(QL(5)+p1)/QL(1)
    h2=(QR(:,5)+p2)/QR(:,1)
    SELECT CASE(dir)
    CASE(1)
       
       DO i=0,N
          Fsharp(i,1)=0.25_RP*(QL(1)+QR(i,1))*(QL(2)/QL(1)+QR(i,2)/QR(i,1))
          Fsharp(i,2)=Fsharp(i,1)*0.5_RP*(QL(2)/QL(1)+QR(i,2)/QR(i,1))+0.5_RP&
               &*(p1+p2(i))
          Fsharp(i,3)=Fsharp(i,1)*0.5_RP*(QL(3)/QL(1)+QR(i,3)/QR(i,1))
          Fsharp(i,4)=Fsharp(i,1)*0.5_RP*(QL(4)/QL(1)+QR(i,4)/QR(i,1))
          Fsharp(i,5)=Fsharp(i,1)*0.5_RP*(h1+h2(i))
       END DO
!
    CASE(2)
       DO i=0,N
          Fsharp(i,1)=0.25_RP*(QL(1)+QR(i,1))*(QL(3)/QL(1)+QR(i,3)/QR(i,1))
          Fsharp(i,3)=Fsharp(i,1)*0.5_RP*(QL(3)/QL(1)+QR(i,3)/QR(i,1))+0.5_RP&
               &*(p1+p2(i))
          Fsharp(i,2)=Fsharp(i,1)*0.5_RP*(QL(2)/QL(1)+QR(i,2)/QR(i,1))
          Fsharp(i,4)=Fsharp(i,1)*0.5_RP*(QL(4)/QL(1)+QR(i,4)/QR(i,1))
          Fsharp(i,5)=Fsharp(i,1)*0.5_RP*(h1+h2(i))
       END DO
!
    CASE(3)

       DO i=0,N
          Fsharp(i,1)=0.25_RP*(QL(1)+QR(i,1))*(QL(4)/QL(1)+QR(i,4)/QR(i,1))
          Fsharp(i,4)=Fsharp(i,1)*0.5_RP*(QL(4)/QL(1)+QR(i,4)/QR(i,1))+0.5_RP&
               &*(p1+p2(i))
          Fsharp(i,3)=Fsharp(i,1)*0.5_RP*(QL(3)/QL(1)+QR(i,3)/QR(i,1))
          Fsharp(i,2)=Fsharp(i,1)*0.5_RP*(QL(2)/QL(1)+QR(i,2)/QR(i,1))
          Fsharp(i,5)=Fsharp(i,1)*0.5_RP*(h1+h2(i))
       END DO
    END SELECT
  END SUBROUTINE computeEulerPI

  !////////////////////////////////////////////////////////////////////
  
  SUBROUTINE computeEulerPIinterf(QL,QR,N,nEqn,dir,Fsharp)
    IMPLICIT NONE
!does the same as computeEulerPI but only takes one node as an argument
    INTEGER                           ,INTENT(IN) :: N, nEqn,dir
    REAL(KIND=RP),DIMENSION(nEqn)     ,INTENT(IN) :: QL,QR
    REAL(KIND=RP),DIMENSION(nEqn) ,INTENT(OUT):: Fsharp

    !Local Variables

    INTEGER                            :: i,j
    REAL(KIND=RP)                      :: p1,h1,p2,h2
    
    p1=(gamma-1.0_RP)*(QL(5)-0.5_RP*(QL(2)*QL(2)+QL(3)*QL(3)+QL(4)*QL(4))/QL(1))
    p2=(gamma-1.0_RP)*(QR(5)-0.5_RP*(QR(2)*QR(2)+QR(3)*QR(3)+QR(4)*QR(4))/QR(1))
    h1=(QL(5)+p1)/QL(1)
    h2=(QL(5)+p2)/QL(1)
    SELECT CASE(dir)
    CASE(1)
       
       Fsharp(1)=0.25_RP*(QL(1)+QR(1))*(QL(2)/QL(1)+QR(2)/QR(1))
       Fsharp(2)=Fsharp(1)*0.5_RP*(QL(2)/QL(1)+QR(2)/QR(1))+0.5_RP&
            &*(p1+p2)
       Fsharp(3)=Fsharp(1)*0.5_RP*(QL(3)/QL(1)+QR(3)/QR(1))
       Fsharp(4)=Fsharp(1)*0.5_RP*(QL(4)/QL(1)+QR(4)/QR(1))
       Fsharp(5)=Fsharp(1)*0.5_RP*(h1+h2)
!
    CASE(2)
       
       Fsharp(1)=0.25_RP*(QL(1)+QR(1))*(QL(3)/QL(1)+QR(3)/QR(1))
       Fsharp(3)=Fsharp(1)*0.5_RP*(QL(3)/QL(1)+QR(3)/QR(1))+0.5_RP&
            &*(p1+p2)
       Fsharp(2)=Fsharp(1)*0.5_RP*(QL(2)/QL(1)+QR(2)/QR(1))
       Fsharp(4)=Fsharp(1)*0.5_RP*(QL(4)/QL(1)+QR(4)/QR(1))
       Fsharp(5)=Fsharp(1)*0.5_RP*(h1+h2)
!
    CASE(3)
       
       Fsharp(1)=0.25_RP*(QL(1)+QR(1))*(QL(4)/QL(1)+QR(4)/QR(1))
       Fsharp(4)=Fsharp(1)*0.5_RP*(QL(4)/QL(1)+QR(4)/QR(1))+0.5_RP&
            &*(p1+p2)
       Fsharp(3)=Fsharp(1)*0.5_RP*(QL(3)/QL(1)+QR(3)/QR(1))
       Fsharp(2)=Fsharp(1)*0.5_RP*(QL(2)/QL(1)+QR(2)/QR(1))
       Fsharp(5)=Fsharp(1)*0.5_RP*(h1+h2)

    END SELECT
  END SUBROUTINE computeEulerPIinterf
  
 !////////////////////////////////////////////////////////////////////

  SUBROUTINE RiemannSolver(QR,QL,Fout,nEqn,N,dir,lambdamax)
    IMPLICIT NONE
    REAL(KIND=RP),DIMENSION(0:N,0:N,nEqn),INTENT(IN)  :: QR,QL
    REAL(KIND=RP),DIMENSION(0:N,0:N,nEqn),INTENT(OUT) :: Fout
    INTEGER                              ,INTENT(IN)  :: nEqn,N,dir
    REAL(KIND=RP)                        ,INTENT(IN)  :: lambdamax

    !local variables

    REAL(KIND=RP),DIMENSION(0:N,0:N,nEqn)   :: Fsharp
    INTEGER                                 :: i,j

    DO i=0,N
       DO j=0,N
          CALL computeEulerPIinterf(QL(i,j,:),QR(i,j,:),N,nEqn,dir&
               &,Fsharp(i,j,:))
       END DO
    END DO
    
    Fout=Fsharp-0.5_RP*lambdamax*(-QR+QL)
  END SUBROUTINE RiemannSolver

  !////////////////////////////////////////////////////////////////////

  SUBROUTINE EulerAnalyticFlux(Q,Fout,dir,N,nEqn)
    IMPLICIT NONE
    INTEGER                              ,INTENT(IN) :: N, nEqn, dir
    REAL(KIND=RP),DIMENSION(0:N,0:N,nEqn),INTENT(IN) :: Q
    REAL(KIND=RP),DIMENSION(0:N,0:N,nEqn),INTENT(OUT):: Fout

    !local variables
    REAL(KIND=RP),DIMENSION(0:n,0:n) :: p,h
    p=(gamma-1.0_RP)*(Q(:,:,5)-0.5_RP*(Q(:,:,2)*Q(:,:,2)+Q(:,:,3)*Q(:&
         &,:,3)+Q(:,:,4)*Q(:,:,4))/Q(:,:,1))
    h=(Q(:,:,5)+p)/Q(:,:,1)
    SELECT CASE(dir)
    CASE(1)
       Fout(:,:,1)=Q(:,:,2)
       Fout(:,:,2)=Q(:,:,2)**2/Q(:,:,1)+p
       Fout(:,:,3)=Q(:,:,2)*Q(:,:,3)/Q(:,:,1)
       Fout(:,:,4)=Q(:,:,2)*Q(:,:,4)/Q(:,:,1)
       Fout(:,:,5)=Q(:,:,2)*h
    CASE(2)
       Fout(:,:,1)=Q(:,:,3)
       Fout(:,:,2)=Q(:,:,3)*Q(:,:,2)/Q(:,:,1)
       Fout(:,:,3)=Q(:,:,3)**2/Q(:,:,1)+p
       Fout(:,:,4)=Q(:,:,3)*Q(:,:,4)/Q(:,:,1)
       Fout(:,:,5)=Q(:,:,3)*h
    CASE(3)
       Fout(:,:,1)=Q(:,:,4)
       Fout(:,:,2)=Q(:,:,4)*Q(:,:,2)/Q(:,:,1)
       Fout(:,:,3)=Q(:,:,4)*Q(:,:,3)/Q(:,:,1)
       Fout(:,:,4)=Q(:,:,4)**2/Q(:,:,1)+p
       Fout(:,:,5)=Q(:,:,4)*h
    END SELECT
  END SUBROUTINE EulerAnalyticFlux
  
END MODULE FluxRoutines

       
       
       
       
    
