module userInput
  implicit none
  public
  double precision, parameter :: &
    E = 200.0d9, &
    ! E = 1.4d6, &
    v = 0.3d0, &
    rho = 7850, &
    Lx = 10.0d0, &
    Ly = 1.0d0, &
    Lz = 0.1d0, &
    ! totalLoad = -rho*2.0d0, &
    totalLoad = -6.0d7
  ! pointLoad = -200000.0d0
  integer, parameter :: &
    degEl = 2, &
    flag = 0, &
    nElx = 40, &
    nEly = 4, &
    nDofPerNode = 2

  double precision, parameter :: &
    alpha = 0.25d0, & !Newmark Parameters
    beta = 0.5d0, &
    tStart = 0.0d0, &
    tEnd = 3.0d0, &
    dt = 0.01d0

  integer, parameter ::          noOfLoadSteps = 1

contains

  subroutine setupBC(dofBC)
    implicit none
    integer :: i!, nNodes, nDof, cnt
    integer, allocatable, dimension(:), intent(out) :: dofBC

    allocate (dofBC((nEly*degEl + 1)*nDofPerNode))

    ! nNodes = (nElx*degEl + 1)*(nEly*degEl + 1) - flag*(nElx*nEly)*(degEl - 1)**2 !total no of global nodes
    ! nDof = nDofPerNode*nNodes !total DOF including BC DOFs

    ! ! cnt = 0
    ! do i = 1, nDof, (nElx*degEl + 1)*nDofPerNode
    !    cnt = cnt + 1
    !    dofBC(cnt) = i !Global DOF where BC is specified on primary vars
    !    cnt = cnt + 1
    !    dofBC(cnt) = i + 1
    ! end do

    do i = 1, (nEly*degEl + 1)*nDofPerNode
      dofBC(i) = i !Global DOF where BC is specified on primary vars
    end do

  end subroutine setupBC

end module userInput
