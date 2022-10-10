module userInput
  implicit none
  public
  double precision, parameter :: &
    E = 200.0d9, &
    v = 0.4d0, &
    rho = 1.0d3, &
    Lx = 0.35d0, &
    Ly = 0.02d0, &
    Lz = 1.0d0, &
    totalLoad = -rho*2.0d0, &
    pointLoad = -200000.0d0

  integer, parameter :: &
    degEl = 2, &
    flag = 0, &
    nElx = 35, &
    nEly = 2, &
    nDofPerNode = 2

  double precision, parameter :: &
    alpha = 0.25d0, & !Newmark Parameters
    beta = 0.5d0, &
    tStart = 0.0d0, &
    tEnd = 2.0d0, &
    dt = 0.02d0

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
