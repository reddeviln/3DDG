MODULE DGMeshClass
  USE DGElementClass
  USE FluxRoutines
  IMPLICIT NONE

  TYPE NodePointers
     INTEGER :: eLeft, eRight
  END TYPE NodePointers
  
  TYPE DGMesh
     INTEGER                                     :: K
     REAL(KIND=RP)                               :: lambdamax
     TYPE(NodalDGStorage)                        :: DG
     TYPE(DGElement)   ,ALLOCATABLE,DIMENSION(:) :: e
     TYPE(NodePointers),ALLOCATABLE,DIMENSION(:) :: px,py,pz
  END TYPE DGMesh

CONTAINS

  SUBROUTINE ConstructMesh3D(this,x_nodes,NQ,N,nEqn,lambdamax)
    IMPLICIT NONE
    TYPE(DGMesh)                    ,INTENT(INOUT)  :: this
    INTEGER                         ,INTENT(IN)     :: NQ,N,nEqn
    REAL(KIND=RP),DIMENSION(0:NQ)   ,INTENT(IN)     :: x_nodes
    REAL(KIND=RP)                   ,INTENT(IN)     :: lambdamax
    
    !Local Variables

    INTEGER :: i,j,l,K

    k=nq**3
    
    this%K = K
    this%lambdamax=lambdamax
    ALLOCATE(this%e(K),this%px(K),this%py(K),this%pz(K))
    CALL ConstructNodalDGStorage(this%DG,N)
    print*,"Constructing the Elements..."
    DO l=1,NQ
       DO j=1,NQ
          DO i=1,NQ
             CALL ConstructDGElement(this%e(i+(j-1)*NQ+(l-1)*NQ**2),nEqn,x_nodes(i-1),x_nodes(i)&
                  &,x_nodes(j-1),x_nodes(j),x_nodes(l-1),x_nodes(l),N)
          END DO
       END DO
    END DO
    print*,"Setting the element neighbors.."
    !setting the  element pointers
    this%px%eLeft=0
    this%px%eRight=0
    this%py%eLeft=0
    this%py%eRight=0
    this%pz%eLeft=0
    this%pz%eRight=0
    !x-neighbors
    DO l=1,nq
       DO j=1,nq
          DO i=1+(j-1)*nq+(l-1)*nq**2,nq+nq*(j&
               &-1)+(l-1)*nq**2-1
             this%px(i)%eLeft = i
             this%px(i)%eRight= i+1
          ENDDO
       ENDDO
    ENDDO
    DO j=1,nq**2
        l=j*nq
        this%px(l)%eLeft = l
        this%px(l)%eRight= l-nq+1
    ENDDO
    !y-neighbors
    DO j= 1,nq
       DO i=(1)+(j-1)*nq**2,j*nq**2-nq
          this%py(i)%eLeft = i
          this%py(i)%eRight= i+nq
       ENDDO
    ENDDO
    DO j=1,nq
       DO i=1+j*nq**2-nq,j*nq**2
          this%py(i)%eLeft = i
          this%py(i)%eRight= i-nq**2+nq
       ENDDO
    ENDDO
    
    !z-neighbors
    DO i=1,K-nq**2
       this%pz(i)%eLeft =i
       this%pz(i)%eRight=i+nq**2
    ENDDO
    DO i=K-nq**2+1,K
       this%pz(i)%eLeft = i
       this%pz(i)%eRight= i-k+nq**2
    ENDDO
    IF(any(this%px%eLeft==0).OR. any(this%py%eLeft==0).OR. any(this%pz%eLeft&
         ==0)) THEN
       print*,"Initialization failed. Pointers were not set correctly!"
       CALL exit(1)
    END IF
    

  END SUBROUTINE ConstructMesh3D

  !//////////////////////////////////////

  SUBROUTINE GlobalTimeDerivative(this,t)
    IMPLICIT NONE
    TYPE(DGMESH) ,INTENT(INOUT)                     :: this
    REAL(KIND=RP),INTENT(IN)                        :: t

    !local variables
    
    REAL(KIND=RP),DIMENSION(0:this%DG%N,0:this%DG%N,&
         & this%e(1)%nEqn) :: F,G,H
    INTEGER                                         :: i,j,idL,idR&
         &,nEqn,N

    N    = this%DG%N
    nEqn = this%e(1)%nEqn

    !Solve Riemann problem to get numerical flux at the interface
    !Set boundary data
    DO i=1,this%K
       this%e(i)%QLx=this%e(i)%Q(0,:,:,:)
       this%e(i)%QRx=this%e(i)%Q(N,:,:,:)
       this%e(i)%QLy=this%e(i)%Q(:,0,:,:)
       this%e(i)%QRy=this%e(i)%Q(:,N,:,:)
       this%e(i)%QLz=this%e(i)%Q(:,:,0,:)
       this%e(i)%QRz=this%e(i)%Q(:,:,N,:)
    END DO
    CALL getLambdaMaxLocally(this)
       !x-direction
    DO j=1,this%K
       idL=this%px(j)%eLeft
       idR=this%px(j)%eRight
       CALL RiemannSolver(this%e(idL)%QRx,this%e(idR)%QLx,F,this%e(idR)&
            &%nEqn,N,1,max(this%e(idL)%lambdamax(2),this%e(idR)%lambdamax(1)))
       this%e(idR)%FstarL = -F
       this%e(idL)%FstarR = F
    ENDDO

    !y-direction
    DO j=1,this%K
       idL=this%py(j)%eLeft
       idR=this%py(j)%eRight
       CALL RiemannSolver(this%e(idL)%QRy,this%e(idR)%QLy,G,this%e(idR)&
            &%nEqn,N,2,max(this%e(idL)%lambdamax(4),this%e(idR)%lambdamax(3)))
       this%e(idR)%GstarL = -G
       this%e(idL)%GstarR = G
    ENDDO

    !z-direction
    DO j=1,this%K
       idL=this%pz(j)%eLeft
       idR=this%pz(j)%eRight
       CALL RiemannSolver(this%e(idL)%QRz,this%e(idR)%QLz,H,this%e(idR)&
            &%nEqn,N,3,max(this%e(idL)%lambdamax(6),this%e(idR)%lambdamax(5)))
       this%e(idR)%HstarL = -H
       this%e(idL)%HstarR = H
    ENDDO
    !Compute local time derivative on each element
    DO j=1,this%K
       CALL LocalTimeDerivative(this%e(j),this%DG,t)
    ENDDO
    
  END SUBROUTINE GlobalTimeDerivative
     
  !//////////////////////////////////////

  SUBROUTINE getLambdaMaxLocally(this)
    IMPLICIT NONE
    TYPE(DGMesh),INTENT(INOUT) :: this

    !local variables

    REAL(KIND=RP),DIMENSION(this%K,0:this%dg%n,0:this%dg%n) :: pLx&
         &,pRx,pLy,pRy,pLz,pRz,cLx,cRx,cLy,cRy,cLz,cRz
    INTEGER :: i
    DO i=1,this%K
       pLx(i,:,:)=(gamma-1.0_RP)*(this%e(i)%QLx(:,:,5)-0.5_RP&
            &*(this%e(i)%QLx(:,:,2)**2+this%e(i)%QLx(:,:,3)**2&
            &+this%e(i)%QLx(:,:,4)**2)/this%e(i)%QLx(:,:,1))
       pRx(i,:,:)=(gamma-1.0_RP)*(this%e(i)%QRx(:,:,5)-0.5_RP&
            &*(this%e(i)%QRx(:,:,2)**2+this%e(i)%QRx(:,:,3)**2&
            &+this%e(i)%QRx(:,:,4)**2)/this%e(i)%QRx(:,:,1))
       pLy(i,:,:)=(gamma-1.0_RP)*(this%e(i)%QLy(:,:,5)-0.5_RP&
            &*(this%e(i)%QLy(:,:,2)**2+this%e(i)%QLy(:,:,3)**2&
            &+this%e(i)%QLy(:,:,4)**2)/this%e(i)%QLy(:,:,1))
       pRy(i,:,:)=(gamma-1.0_RP)*(this%e(i)%QRy(:,:,5)-0.5_RP&
            &*(this%e(i)%QRy(:,:,2)**2+this%e(i)%QRy(:,:,3)**2&
            &+this%e(i)%QRy(:,:,4)**2)/this%e(i)%QRy(:,:,1))
       pLz(i,:,:)=(gamma-1.0_RP)*(this%e(i)%QLz(:,:,5)-0.5_RP&
            &*(this%e(i)%QLz(:,:,2)**2+this%e(i)%QLz(:,:,3)**2&
            &+this%e(i)%QLz(:,:,4)**2)/this%e(i)%QLz(:,:,1))
       pRz(i,:,:)=(gamma-1.0_RP)*(this%e(i)%QRz(:,:,5)-0.5_RP&
            &*(this%e(i)%QRz(:,:,2)**2+this%e(i)%QRz(:,:,3)**2&
            &+this%e(i)%QRz(:,:,4)**2)/this%e(i)%QRz(:,:,1))
    END DO
    
    DO i=1,this%K
       IF(ANY(pLx(i,:,:)<0.OR.this%e(i)%QLx(:,:,1)<0).OR.ANY(pRx(i,:&
            &,:) <0.OR.this%e(i)%QRx(:,:,1)<0).OR.ANY(pLy(i,:,:)&
            &<0.OR.this%e(i)%QLy(:,:,1)<0).OR.ANY(pRy(i,:,:)<0.OR.this%e(i)%QRy(:,:&
            &,1)<0) .OR. ANY(pLz(i,:,:)<0.OR.this%e(i)%QLz(:,:,1)<0) .OR.&
            & ANY(pRz(i,:,:)<0.OR. this%e(i)%QRz(:,:,1)<0)) THEN
          print*,"pressure/density negativ in Element",i
          CALL exit(2);
       END IF
       cLx(i,:,:)=sqrt(gamma*pLx(i,:,:)/this%e(i)%QLx(:,:,1))
       cRx(i,:,:)=sqrt(gamma*pRx(i,:,:)/this%e(i)%QRx(:,:,1))
       cLy(i,:,:)=sqrt(gamma*pLy(i,:,:)/this%e(i)%QLy(:,:,1))
       cRy(i,:,:)=sqrt(gamma*pRy(i,:,:)/this%e(i)%QRy(:,:,1))
       cLz(i,:,:)=sqrt(gamma*pLz(i,:,:)/this%e(i)%QLz(:,:,1))
       cRz(i,:,:)=sqrt(gamma*pRz(i,:,:)/this%e(i)%QRz(:,:,1))
    END DO
    DO i=1,this%K
       this%e(i)%lambdamax(1)=max(maxval(abs(this%e(i)%QLx(:,:,2)/this&
            &%e(i)%QLx(:,:,1)+cLx(i,:,:))),maxval(abs(this%e(i)%QLx(:,:&
            &,3)/this%e(i)%QLx(:,:,1)+cLx(i,:,:))),maxval(abs(this%e(i)&
            &%QLx(:,:,4)/this%e(i)%QLx(:,:,1)+cLx(i,:,:))),maxval(abs(this%e(i)%QLx(:,:,2)/this&
            &%e(i)%QLx(:,:,1)-cLx(i,:,:))),maxval(abs(this%e(i)%QLx(:,:&
            &,3)/this%e(i)%QLx(:,:,1)-cLx(i,:,:))),maxval(abs(this%e(i)&
            &%QLx(:,:,4)/this%e(i)%QLx(:,:,1)-cLx(i,:,:))))
       this%e(i)%lambdamax(2)=max(maxval(abs(this%e(i)%QRx(:,:,2)/this&
            &%e(i)%QRx(:,:,1)+cRx(i,:,:))),maxval(abs(this%e(i)%QRx(:,:&
            &,3)/this%e(i)%QRx(:,:,1)+cRx(i,:,:))),maxval(abs(this%e(i)&
            &%QRx(:,:,4)/this%e(i)%QRx(:,:,1)+cRx(i,:,:))),maxval(abs(this%e(i)%QRx(:,:,2)/this&
            &%e(i)%QRx(:,:,1)-cRx(i,:,:))),maxval(abs(this%e(i)%QRx(:,:&
            &,3)/this%e(i)%QRx(:,:,1)-cRx(i,:,:))),maxval(abs(this%e(i)&
            &%QRx(:,:,4)/this%e(i)%QRx(:,:,1)-cRx(i,:,:))))
       this%e(i)%lambdamax(3)=max(maxval(abs(this%e(i)%QLy(:,:,2)/this&
            &%e(i)%QLy(:,:,1)+cLy(i,:,:))),maxval(abs(this%e(i)%QLy(:,:&
            &,3)/this%e(i)%QLy(:,:,1)+cLy(i,:,:))),maxval(abs(this%e(i)&
            &%QLy(:,:,4)/this%e(i)%QLy(:,:,1)+cLy(i,:,:))),maxval(abs(this%e(i)%QLy(:,:,2)/this&
            &%e(i)%QLy(:,:,1)-cLy(i,:,:))),maxval(abs(this%e(i)%QLy(:,:&
            &,3)/this%e(i)%QLy(:,:,1)-cLy(i,:,:))),maxval(abs(this%e(i)&
            &%QLy(:,:,4)/this%e(i)%QLy(:,:,1)-cLy(i,:,:))))
       this%e(i)%lambdamax(4)=max(maxval(abs(this%e(i)%QRy(:,:,2)/this&
            &%e(i)%QRy(:,:,1)+cRy(i,:,:))),maxval(abs(this%e(i)%QRy(:,:&
            &,3)/this%e(i)%QRy(:,:,1)+cRy(i,:,:))),maxval(abs(this%e(i)&
            &%QRy(:,:,4)/this%e(i)%QRy(:,:,1)+cRy(i,:,:))),maxval(abs(this%e(i)%QRy(:,:,2)/this&
            &%e(i)%QRy(:,:,1)-cRy(i,:,:))),maxval(abs(this%e(i)%QRy(:,:&
            &,3)/this%e(i)%QRy(:,:,1)-cRy(i,:,:))),maxval(abs(this%e(i)&
            &%QRy(:,:,4)/this%e(i)%QRy(:,:,1)-cRy(i,:,:))))
       this%e(i)%lambdamax(5)=max(maxval(abs(this%e(i)%QLz(:,:,2)/this&
            &%e(i)%QLz(:,:,1)+cLz(i,:,:))),maxval(abs(this%e(i)%QLz(:,:&
            &,3)/this%e(i)%QLz(:,:,1)+cLz(i,:,:))),maxval(abs(this%e(i)&
            &%QLz(:,:,4)/this%e(i)%QLz(:,:,1)+cLz(i,:,:))),maxval(abs(this%e(i)%QLz(:,:,2)/this&
            &%e(i)%QLz(:,:,1)-cLz(i,:,:))),maxval(abs(this%e(i)%QLz(:,:&
            &,3)/this%e(i)%QLz(:,:,1)-cLz(i,:,:))),maxval(abs(this%e(i)&
            &%QLz(:,:,4)/this%e(i)%QLz(:,:,1)-cLz(i,:,:))))
       this%e(i)%lambdamax(6)=max(maxval(abs(this%e(i)%QRz(:,:,2)/this&
            &%e(i)%QRz(:,:,1)+cRz(i,:,:))),maxval(abs(this%e(i)%QRz(:,:&
            &,3)/this%e(i)%QRz(:,:,1)+cRz(i,:,:))),maxval(abs(this%e(i)&
            &%QRz(:,:,4)/this%e(i)%QRz(:,:,1)+cRz(i,:,:))),maxval(abs(this%e(i)%QRz(:,:,2)/this&
            &%e(i)%QRz(:,:,1)-cRz(i,:,:))),maxval(abs(this%e(i)%QRz(:,:&
            &,3)/this%e(i)%QRz(:,:,1)-cRz(i,:,:))),maxval(abs(this%e(i)&
            &%QRz(:,:,4)/this%e(i)%QRz(:,:,1)-cRz(i,:,:))))
    END DO
    
  END SUBROUTINE getLambdaMaxLocally
  
  
  !///////////////////////////////////////
  
  SUBROUTINE DestructDGMesh(this)
    TYPE(DGMesh),INTENT(INOUT) :: this
    INTEGER                    :: i
    
    DO i=1,this%K
       CALL DestructElement(this%e(i))
    ENDDO
    CALL DestructNodalDGStorage(this%DG)
    DEALLOCATE(this%e,this%px,this%py,this%pz)
    this%K=-1
  END SUBROUTINE DestructDGMesh
  
 
END MODULE DGMeshClass
