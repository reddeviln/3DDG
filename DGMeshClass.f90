MODULE DGMeshClass
  USE DGElementClass
  IMPLICIT NONE

  TYPE NodePointers
     INTEGER :: eLeft, eRight
  END TYPE NodePointers
  
  TYPE DGMesh
     INTEGER                                     :: K
     TYPE(NodalDGStorage)                        :: DG
     TYPE(DGElement)   ,ALLOCATABLE,DIMENSION(:) :: e
     TYPE(NodePointers),ALLOCATABLE,DIMENSION(:) :: px,py,pz
  END TYPE DGMesh

CONTAINS

  SUBROUTINE ConstructMesh3D(this,x_nodes,K,N,nEqn)
    IMPLICIT NONE
    TYPE(DGMesh)                  ,INTENT(INOUT)  :: this
    INTEGER                       ,INTENT(IN)     :: K,N,nEqn
    REAL(KIND=RP),DIMENSION(0:K)  ,INTENT(IN)     :: x_nodes
    
    !Local Variables

    INTEGER :: i,j,l,nq

    nq=k**(1/3)

    this%K = K
    ALLOCATE(this%e(K),this%px(0:K),this%py(0:K),this%pz(0:K))
    CALL ConstructNodalDGStorage(this%DG,N)
    DO i=1,K
       CALL ConstructDGElement(this%e(i),nEqn,x_nodes(i-1),x_nodes(i)&
            &,x_nodes(i-1),x_nodes(i),x_nodes(i-1),x_nodes(i),N)
    ENDDO
    !setting the neighboring element pointers
    this%px%eLeft=0
    this%px%eRight=0
    !x-neighbors
    DO l=1,nq
       DO j=1,nq
          DO i=2+(j-1)*nq+(l-1)*nq**2,nq+nq*(j&
               &-1)+(l-1)*nq**2-1
             this%px(i)%eLeft = i-1
             this%px(i)%eRight= i+1
          ENDDO
       ENDDO
    ENDDO
    DO j=1,nq**2
        i=1+(j-1)*nq
        this%px(i)%eLeft = i+nq-1
        this%px(i)%eRight= i+1
        l=j*nq
        this%px(l)%eLeft = i-1
        this%px(l)%eRight= i-nq+1
    ENDDO
    !y-neighbors
    DO j= 1,nq
       DO i=(nq+1)+(j-1)*nq**2,(3*nq)+(j-1)*nq**2
          this%py(i)%eLeft = i+nq
          this%py(i)%eRight= i-nq
       ENDDO
    ENDDO
    DO j=1,nq
       DO i=1+(j-1)*nq**2,nq+(j-1)*nq**2
          this%py(i)%eLeft = i+nq
          this%py(i)%eRight= i+nq**2-nq
       ENDDO
       DO l=(j)*nq**2-nq+1,j*nq**2
          this%py(i)%eLeft = i-nq**2+nq
          this%py(i)%eRight= i-nq
       ENDDO
    ENDDO
    
    !z-neighbors
    DO i=nq**2+1,K-nq**2
       this%pz(i)%eLeft =i-nq**2
       this%pz(i)%eRight=i+nq**2
    ENDDO
    DO i=1,nq**2
       this%pz(i)%eLeft = i+k-nq**2
       this%pz(i)%eRight= i+nq**2
    ENDDO
    DO i=K-nq**2+1,K
       this%pz(i)%eLeft = i-nq**2
       this%pz(i)%eRight= i-k+nq**2
    ENDDO
  END SUBROUTINE ConstructMesh3D

  !//////////////////////////////////////

  SUBROUTINE DestructDGMesh(this)
    TYPE(DGMesh),INTENT(INOUT) :: this
    INTEGER                    :: i
    
    DO i=1,this%K
       CALL DestructDGElement(this%e(i))
    ENDDO
    CALL DestructNodalDGStorage(this%DG)
    DEALLOCATE(this%e,this%px,this%py,this%pz)
    this%K=-1
  END SUBROUTINE DestructDGMesh
  
 
END MODULE DGMeshClass
