module input
   implicit none
   ! public

   integer, parameter:: &
      dim = 2, &
      chanH_ = 164, &
      q = 9, &
      time_ = 500000, &
      noOfSnaps = 3, &
      dispFreq = 100, &
      avgSpan = 2

   double precision, parameter:: &
      rhoF_ = 1.0d0, &
      rhoF = 1000.0d0, &
      chanL = 2.5d0, & !Length of channel
      chanH = 0.41d0, & !Width of channel
      barL = 0.35d0, &
      barH = 0.02d0, &
      uMean = 2.0d0, &
      nu = 0.001d0, &
      dia = 0.1d0, &
      xc = 0.2d0, &
      yc = 0.2d0

   double precision, parameter :: &
      ! E = uMean**2.0d0*1.4d6, &
      E = 1.4d6, &
      v = 0.4d0, &
      rhoS = 1.0d3, &
      ! Lx = 0.02d0, &
      ! Ly = 0.35d0, &
      ! Lz = 1.0d0
      totalLoad = rhoS*2.0d0

   integer, parameter :: &
      degEl = 2, &
      flag = 0, &
      nElx = 140, &
      nEly = 8, &
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

      allocate (dofBC((nEly*degEl + 1)*nDofPerNode))

      do i = 1, (nEly*degEl + 1)*nDofPerNode
         dofBC(i) = i !Global DOF where BC is specified on primary vars
      end do

   end subroutine setupBC

end module input
