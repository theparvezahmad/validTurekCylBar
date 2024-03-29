module input
  implicit none
  public

  character(len=10), parameter:: case = "C4"
  !C0: No solid shape update
  integer, parameter:: &
    dim = 2, &
    chanH_ = 205, &
    q = 9, &
    totTime_ = 200000, &
    noOfSnaps = 401, &
    dispFreq = 1, &
    avgSpan = 2

  double precision, parameter:: &
    rhoF_ = 1.0d0, &
    rhoF = 1000.0d0, &
    chanL = 2.5d0, & !Length of channel
    chanH = 0.41d0, & !Height of channel
    chanW = 1.0d0, & !Width of channel(into the plane)
    barL = 0.35d0, &
    barH = 0.02d0, &
    uMean = 0.2d0, &
    nu = 0.001d0, &
    dia = 0.1d0, &
    xc = 0.2d0, &
    yc = 0.2d0

  double precision, parameter :: &
    E = rhoF*uMean**2.0d0*3.5d4, &
    ! E = uMean**2.0d0*1.4d6, &
    ! E = 1.4d6, &
    v = 0.4d0, &
    rhoS = 1.0d3, &
    ! Lx = 0.02d0, &
    ! Ly = 0.35d0, &
    ! Lz = 1.0d0
    totalLoad = -rhoS*2.0d0

  integer, parameter :: &
    degEl = 2, &
    flag = 0, &
    nElx = 35, &
    nEly = 2, &
    nDofPerNode = 2

  double precision, parameter :: &
    alpha = 0.25d0, & !Newmark Parameters
    beta = 0.5d0
  ! tStart = 0.0d0, &
  ! tEnd = 10.0d0, &
  ! dt = 0.02d0
  integer, parameter ::          noOfLoadSteps = 1

  integer, allocatable, dimension(:), protected :: dofBC

contains

  subroutine setupBC()
    implicit none
    integer :: i
    ! integer, allocatable, dimension(:), intent(out) :: dofBC

    allocate (dofBC((nEly*degEl + 1)*nDofPerNode))

    do i = 1, (nEly*degEl + 1)*nDofPerNode
      dofBC(i) = i !Global DOF where BC is specified on primary vars
    end do

  end subroutine setupBC

end module input
