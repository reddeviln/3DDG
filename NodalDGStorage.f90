MODULE NodalDGStorageClass
  USE DGToolbox
  IMPLICIT NONE

  !
  TYPE NodalDGStorage
     INTEGER                                  :: N
     REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)   :: l_at_minus_one
     REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)   :: l_at_one
     REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)   :: weights
     REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:) :: D
     REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:) :: D_hat
  END TYPE NodalDGStorage
  !
CONTAINS

  SUBROUTINE ConstructNodalDGStorage(this,N)
    IMPLICIT NONE
    
    INTEGER             ,INTENT(IN)   :: N
    TYPE(NodalDGStorage),INTENT(OUT)  :: this
    !Local Variables
    REAL(KIND=RP),ALLOCATABLE,DIMENSION(0:N) :: x, BaryWeights !nodes and weights of interp.
    INTEGER                                  :: i,j

    ALLOCATE(this%l_at_minus_one(0:N),this%l_at_one(0:N),this%weights(0:N))
    ALLOCATE(this%D(0:N,0:N),this%D_hat(0:N,0:N)
    this%l_at_minus_one   = 0.0_RP
    this%l_at_one         = 0.0_RP
    this%weights          = 0.0_RP
    this%D                = 0.0_RP
    this%N                = N

    call LegendreGaussLobattoNodesAndWeights(N,x,this%weights)
    BaryWeights=barzweight(x,N)
    do i=0,n
       this%l_at_minus_one(i)=lagrangevalue(x,i,N,-1.0_RP)
       this%l_at_one(i)=lagrangevalue(x,i,N,1.0_RP)
    end do
    this%D=baryzdiffmatrix(x,this%weights,N)
    do j= 0,N
       do i= 0,N
          this%D_hat(i,j)= -this%D(j,i)*(this%weights(j)/this%weights(i))
       enddo
    enddo
  end SUBROUTINE ConstructNodalDGStorage
  !
  !/////////////////////////////////////////////////////////////////////////////////////
  !

  SUBROUTINE DestructNodalDGStorage(this)
    IMPLCIT NONE
    TYPE(NODALDGStorage) :: this

    this%N=-1
    DEALLOCATE(this%l_at_minus_one,this%l_at_one,this%weights,this%D_hat,this%D)
  end SUBROUTINE DestructNodalDGStorage
  !
end MODULE NodalDGStorageClass

