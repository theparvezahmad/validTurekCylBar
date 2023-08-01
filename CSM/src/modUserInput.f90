module userInput
  implicit none
  public
  double precision, parameter :: &
    E = 1.4d6, &
    v = 0.4d0, &
    rho = 1.0d3, &
    Lx = 0.02d0, &
    Ly = 0.35d0, &
    Lz = 1.0d0, &
    totalLoad = rho*2.0d0

  integer, parameter :: &
    degEl = 2, &
    flag = 0, &
    nElx = 16, &
    nEly = 280, &
    nDofPerNode = 2

  double precision, parameter :: &
    alpha = 0.25d0, & !Newmark Parameters
    beta = 0.5d0, &
    tStart = 0.0d0, &
    tEnd = 10.0d0, &
    dt = 0.02d0

  integer, parameter ::          noOfLoadSteps = 1

contains

  subroutine setupBC(dofBC)
    implicit none
    integer :: i
    integer, allocatable, dimension(:), intent(out) :: dofBC

    allocate (dofBC((nElx*degEl + 1)*nDofPerNode))

    do i = 1, (nElx*degEl + 1)*nDofPerNode
      dofBC(i) = i !Global DOF where BC is specified on primary vars
    end do

  end subroutine setupBC

end module userInput
