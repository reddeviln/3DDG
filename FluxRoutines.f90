MODULE FluxRoutines
  IMPLICIT NONE
  !Contains many different flux functions and Riemann solvers
  REAL(KIND=RP)         :: mu,g=9.81,R,gamma

CONTAINS

  SUBROUTINE computeVolumePI(D,Q,dir,Fout,N,nEqn)
    IMPLICIT NONE
    REAL(KIND=RP),DIMENSION(0:N,0:N)         ,INTENT(IN) :: D
    REAL(KIND=RP),DIMENSION(0:N,0:N,0:N,nEqn),INTENT(IN) :: Q
    INTEGER                                  ,INTENT(IN) :: dir,N,nEqn
    REAL(KIND=RP),DIMENSION(0:N,0:N,0:N,nEqn),INTENT(OUT):: Fout

    !Local Variables
    INTEGER                                   :: i,j,k,l
    REAL(KIND=RP),DIMENSION(nEqn,0:N,0:N,0:N) :: Foutc
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
                   Foutc(l,i,j,k)=2.0_RP*dot_product(D(i,:),Fsharp(:,l))
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
                   Foutc(l,i,j,k)=2.0_RP*dot_product(D(i,:),Fsharp(:,l))
                END DO
             END DO
          END DO
       END DO
    CASE DEFAULT
       print*,"specify 1,2 or 3 for flux direction"
    END SELECT
    Fout=cshift(Foutc,-1,1)
    
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
    REAL(KIND=RP),DIMENSION(nEqn)      :: p1,h1
    REAL(KIND=RP),DIMENSION(0:N,nEqn)  :: p2,h2
    
    p1=(gamma-1.0_RP)*(QL(5)-0.5_RP*(QL(2)*QL(2)+QL(3)*QL(3)+QL(4)*QL(4))/QL(1))
    p2=(gamma-1.0_RP)*(QR(:,5)-0.5_RP*(QR(:,2)*QR(:,2)+QR(:,3)*QR(:&
         &,3)+QR(:,4)*QR(:,4))/QR(:,1))
    h1=(QL(5)+p1)/QL(1)
    h2=(QL(:,5)+p2)/QL(:,1)
    SELECT CASE(dir)
    CASE(1)
       
       DO i=0,N
          Fsharp(i,1)=0.25_RP*(QL(1)+QR(i,1))*(QL(2)/QL(1)+QR(i,2)/QR(i,1))
          Fsharp(i,2)=Fsharp(i,1)*0.5_RP*(QL(2)/QL(1)+QR(i,2)/QR(i,1))+0.5_RP&
               &*(p1+p2)
          Fsharp(i,3)=Fsharp(i,1)*0.5_RP*(QL(3)/QL(1)+QR(i,3)/QR(i,1))
          Fsharp(i,4)=Fsharp(i,1)*0.5_RP*(QL(4)/QL(1)+QR(i,4)/QR(i,1))
          Fsharp(i,5)=Fsharp(i,1)*0.5_RP*(h1+h2)
       END DO
!
    CASE(2)
       DO i=0,N
          Fsharp(i,1)=0.25_RP*(QL(1)+QR(i,1))*(QL(3)/QL(1)+QR(i,3)/QR(i,1))
          Fsharp(i,3)=Fsharp(i,1)*0.5_RP*(QL(3)/QL(1)+QR(i,3)/QR(i,1))+0.5_RP&
               &*(p1+p2)
          Fsharp(i,2)=Fsharp(i,1)*0.5_RP*(QL(2)/QL(1)+QR(i,2)/QR(i,1))
          Fsharp(i,4)=Fsharp(i,1)*0.5_RP*(QL(4)/QL(1)+QR(i,4)/QR(i,1))
          Fsharp(i,5)=Fsharp(i,1)*0.5_RP*(h1+h2)
       END DO
!
    CASE(3)

       DO i=0,N
          Fsharp(i,1)=0.25_RP*(QL(1)+QR(i,1))*(QL(4)/QL(1)+QR(i,4)/QR(i,1))
          Fsharp(i,4)=Fsharp(i,1)*0.5_RP*(QL(4)/QL(1)+QR(i,4)/QR(i,1))+0.5_RP&
               &*(p1+p2)
          Fsharp(i,3)=Fsharp(i,1)*0.5_RP*(QL(3)/QL(1)+QR(i,3)/QR(i,1))
          Fsharp(i,2)=Fsharp(i,1)*0.5_RP*(QL(2)/QL(1)+QR(i,2)/QR(i,1))
          Fsharp(i,5)=Fsharp(i,1)*0.5_RP*(h1+h2)
       END DO
    END SELECT
  END SUBROUTINE computeEulerPI

 !////////////////////////////////////////////////////////////////////

  SUBROUTINE RiemannSolver(QR,QL,Fout,nEqn,N,dir,lambdamax)
    IMPLICIT NONE
    REAL(KIND=RP),DIMENSION(0:N,0:N,nEqn),INTENT(IN)  :: QR,QL
    REAL(KIND=RP),DIMENSION(0:N,0:N,nEqn),INTENT(OUT) :: Fout
    INTEGER                              ,INTENT(IN)  :: nEqn,N,dir
    REAL(KIND=RP)                        ,INTENT(IN)  :: lambdamax

    !local variables

    REAL(KIND=RP),DIMENSION(0:N,0:N,nEqn)             :: Fsharp
!///////////////////////////////////////////////////
!@TODO dimension issue!
!/////////////////////////////////////////////////
    CALL computeEulerPI(QL,QR,N,nEqn,dir,Fsharp)
    Fout=Fsharp-0.5_RP*lambdamax*(QR-QL)
  END SUBROUTINE RiemannSolver
END MODULE FluxRoutines

       
       
       
       
    
