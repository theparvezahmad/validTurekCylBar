program turekFSI
   use input
   use math
   use fem
   implicit none
   !Real units are in SI units
   !Underscore Variables are in LBM units

   integer, parameter:: &
      chanH_ = 205, &
      q = 9, &
      time_ = 500000, &
      noOfSnaps = 3, &
      dispFreq = 100

   double precision, parameter:: &
      rhoF_ = 1.0d0, &
      rhoF = 1000.0d0, &
      chanL = 2.5d0, & !Length of channel
      chanH = 0.41d0, & !Width of channel
      barL = 0.35d0, &
      barH = 0.02d0, &
      uMean = 0.2d0, &
      nu = 0.001d0, &
      dia = 0.1d0, &
      xc = 0.2d0, &
      yc = 0.2d0

   double precision, parameter:: &
      d0 = 0.0d0, &
      pi = 4.0d0*datan(1.0d0), &
      one36th = 1.0d0/36.0d0

   integer::nx, ny
   double precision:: nu_, uMean_, uPara_, uParaRamp_, dia_, xc_, yc_, chanL_, barL_, barH_, ii, jj
   double precision:: Clen, Crho, Ct, Cnu, CVel, CFor, tau, t
   integer:: i, j, a, a1, t_, ia, ja, ex(0:q - 1), ey(0:q - 1), kb(0:q - 1), solnumber
   integer, allocatable, dimension(:, :)::isn
   double precision:: tmp1, tmp2, tmp3, rhoAvg, feq, fx_t, fy_t, Cd, Cl
   double precision:: fx(2), fy(2)
   double precision:: wt(0:q - 1)
   double precision, allocatable, dimension(:, :, :):: f, ft
   double precision, allocatable, dimension(:, :):: ux, uy, rho
   logical::isCyl, isBar
   character(len=30):: filename

   integer:: cnt, k
   integer :: nEl, nNodesPerEl, nNodes, nDofPerEl, nDof, nDofBC
   integer, allocatable, dimension(:)::bottomEl, topEl, leftEl, rightEl
   integer, allocatable, dimension(:) :: dofBC
   integer, allocatable, dimension(:, :) :: dofMapBC
   double precision, allocatable, dimension(:, :, :) :: DeltaEl
   integer:: coupleRange
   double precision, allocatable, dimension(:) :: DeltaG
   double precision, allocatable, dimension(:, :)::bounDofTopEl, bounDofBottomEl, bounDofRightEl
   double precision, allocatable, dimension(:, :)::refBounTopEl, refBounBottomEl, refBounRightEl
   double precision, allocatable, dimension(:, :)::bounDispTopEl, bounDispBottomEl, bounDispRightEl, bounCoord
   double precision:: lenElx_, lenEly_
   double precision, allocatable, dimension(:)::xIntersect
   double precision, allocatable, dimension(:)::sortedA, uniqSortedA

   solnumber = 0

   write (*, *) '======================================================'
   write (*, *) 'Program started at :', dateTime()
!===Conversion Factors===
   Clen = chanH/chanH_
   Crho = rhoF/rhoF_
   Ct = dia/uMean*0.002d0 !Time Period of fundamental node
   Cnu = Clen**2.0d0/Ct
   CVel = Clen/Ct
   ! CFor = Crho*Clen**4.0d0*Ct**(-2.0d0)
   CFor = Crho*Clen**3.0d0*Ct**(-2.0d0)

!===Other LBM parameters===
   chanL_ = chanL/Clen
   dia_ = dia/Clen
   xc_ = xc/Clen
   yc_ = yc/Clen
   barL_ = barL/Clen
   barH_ = barH/Clen
   nx = int(chanL_)
   ny = chanH_
   nu_ = nu/Cnu
   uMean_ = uMean/CVel
   tau = 3*nu_ + 0.5d0

   write (*, *) 'Clen = ', Clen
   write (*, *) 'Crho = ', Crho
   write (*, *) 'Ct   = ', Ct
   write (*, *) 'CVel = ', CVel
   write (*, *) 'tau  = ', tau
   write (*, *) 'uMax_ = ', 1.5d0*uMean_
   write (*, *) 'Re_SI = ', uMean*dia/nu
   write (*, *) 'Re_LBM = ', uMean_*dia_/nu_
   ! write (*, *) 'Aborted by User for Confirmation'; stop
!----------------------------------------------------------------------
   allocate (isn(nx + 2, ny + 2))
   allocate (ux(nx + 2, ny + 2))
   allocate (uy(nx + 2, ny + 2))
   allocate (rho(nx + 2, ny + 2))
   allocate (f(0:q - 1, nx + 2, ny + 2))
   allocate (ft(0:q - 1, nx + 2, ny + 2))
!----------------------------------------------------------------------
   ex(0) = 0; ey(0) = 0; 
   ex(1) = 1; ey(1) = 0; 
   ex(2) = 0; ey(2) = 1; 
   ex(3) = -1; ey(3) = 0; 
   ex(4) = 0; ey(4) = -1; 
   ex(5) = 1; ey(5) = 1; 
   ex(6) = -1; ey(6) = 1; 
   ex(7) = -1; ey(7) = -1; 
   ex(8) = 1; ey(8) = -1; 
!----------------------------------------------------------------------
   do a = 0, q - 1
      if (a .eq. 0) wt(a) = 16*one36th
      if (a .ge. 1 .and. a .le. 4) wt(a) = 4*one36th
      if (a .ge. 5 .and. a .le. 8) wt(a) = one36th
   end do
!----------------------------------------------------------------------
   do a = 0, q - 1
      do a1 = a, q - 1
         if (ex(a) + ex(a1) .eq. 0 .and. ey(a) + ey(a1) .eq. 0) then
            kb(a) = a1
            kb(a1) = a
         end if
      end do
   end do
!----------------------------------------------------------------------
   do i = 1, nx + 2
      do j = 1, ny + 2
         uPara_ = d0! 6*uMean_*(ny - (j - 1.5))*(j - 1.5)/ny**2;
         do a = 0, q - 1
            tmp1 = uPara_*ex(a)
            tmp2 = uPara_*uPara_
            f(a, i, j) = wt(a)*rhoF_*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
         end do
      end do
   end do
!----------------------------------------------------------------------
   fx(1) = d0
   fy(1) = d0
   open (unit=10, file="../output/tRhoCdCl.dat")
   write (10, *) "Variables=timeLBM,timeReal,rho,Cd,Cl"
   ! open (unit=11, file="VelSig.dat")

   !-------------------------Working Start---------------------------
   !Assume displacement field X is known numbered from left-bottom spanning height first
   nEl = nElx*nEly
   nNodesPerEl = (degEl + 1)**2 - flag*(degEl - 1)**2 !nodes per element
   nNodes = (nElx*degEl + 1)*(nEly*degEl + 1) - flag*(nElx*nEly)*(degEl - 1)**2 !total no of global nodes
   nDofPerEl = nDofPerNode*nNodesPerEl
   nDof = nDofPerNode*nNodes !total DOF including BC DOFs
   lenElx_ = barL_/nElx
   lenEly_ = barH_/nEly

   call setupBC(dofBC)
   nDofBC = nDof - size(dofBC)

   call setupElemMap(dofMapBC, coupleRange)!check

   allocate (bottomEl(nElx))
   allocate (topEl(nElx))
   allocate (leftEl(nEly))
   allocate (rightEl(nEly))

   allocate (DeltaG(nDofBC))
   DeltaG = 0.0d0 !Put arbitrary values here
   DeltaG(55:56) = [0.0d0, 5.0d0]
   DeltaG(43:44) = [0.0d0, 5.0d0]
   DeltaG(83:84) = [0.0d0, -5.0d0]
   DeltaG(71:72) = [0.0d0, -5.0d0]
   ! DeltaG(111:112) = [2.0d0, 2.0d0]

   ! cnt = 1
   ! k = 1
   ! do while (k .le. (nElx*degEl + 1)*nDofPerNode)
   !    leftNodesU(cnt) = k
   !    leftNodesV(cnt) = k + 1
   !    k = k + 2
   !    cnt = cnt + 1
   ! end do

   ! topNodesU(1) = 0.0d0
   ! topNodesV(1) = 0.0d0
   cnt = 1
   do k = 1, nEl, nEly
      bottomEl(cnt) = k
      cnt = cnt + 1
   end do

   cnt = 1
   do k = nEly, nEl, nEly
      topEl(cnt) = k
      cnt = cnt + 1
   end do

   cnt = 1
   do k = 1, nEly
      leftEl(cnt) = k
      cnt = cnt + 1
   end do

   cnt = 1
   do k = nEly*(nElx - 1) + 1, nEl
      rightEl(cnt) = k
      cnt = cnt + 1
   end do

   allocate (bounDofTopEl(0:size(topEl), 2))
   allocate (bounDofBottomEl(0:size(bottomEl), 2))
   allocate (bounDofRightEl(0:size(rightEl), 2))

   allocate (refBounTopEl(0:size(topEl), 2))
   allocate (refBounBottomEl(0:size(bottomEl), 2))
   allocate (refBounRightEl(0:size(rightEl), 2))

   allocate (bounDispTopEl(0:size(topEl), 2))
   allocate (bounDispBottomEl(0:size(bottomEl), 2))
   allocate (bounDispRightEl(0:size(rightEl), 2))

   allocate (bounCoord(size(topEl) + size(bottomEl) + size(rightEl) + 2, 2))

   refBounTopEl(0, 1:2) = [xc_ + 0.5d0*dia_, yc_ + 0.5d0*barH_]
   do k = 1, size(topEl)
      refBounTopEl(k, 1) = refBounTopEl(0, 1) + k*lenElx_
      refBounTopEl(k, 2) = refBounTopEl(0, 2)
   end do

   refBounBottomEl(0, 1:2) = [xc_ + 0.5d0*dia_, yc_ - 0.5d0*barH_]
   do k = 1, size(bottomEl)
      refBounBottomEl(k, 1) = refBounBottomEl(0, 1) + k*lenElx_
      refBounBottomEl(k, 2) = refBounBottomEl(0, 2)
   end do

   refBounRightEl(0, 1:2) = [xc_ + 0.5d0*dia_ + barL_, yc_ - 0.5d0*barH_]
   do k = 1, size(rightEl)
      refBounRightEl(k, 1) = refBounRightEl(0, 1)
      refBounRightEl(k, 2) = refBounRightEl(0, 2) + k*lenEly_
   end do
   ! cnt = 2
   ! k = (nElx*degEl + 1)*nDofPerNode - 1
   ! do while (k .le. nDofBC)
   !    topNodesU(cnt) = k
   !    topNodesV(cnt) = k + 1
   !    k = k + (nElx*degEl + 1)*nDofPerNode
   !    cnt = cnt + 1
   ! end do

   ! cnt = 2
   ! k = 1
   ! do while (k .le. nDofBC)
   !    bottomNodesU(cnt) = k
   !    bottomNodesV(cnt) = k + 1
   !    k = k + (nElx*degEl + 1)*nDofPerNode
   !    cnt = cnt + 1
   ! end do

   ! cnt = 1
   ! k = nDofBC - (nElx*degEl + 1)*nDofPerNode + 1
   ! do while (k .le. nDofBC)
   !    rightNodesU(cnt) = k
   !    rightNodesV(cnt) = k + 1
   !    k = k + 2
   !    cnt = cnt + 1
   ! end do

   allocate (DeltaEl(nNodesPerEl, nDofPerNode, nEl))
   DeltaEl = 0.0d0

   !Back-Mapping Global unknown vector to Element-level vector
   do i = 1, nEl !spans element
      do j = 1, nDofPerEl !spans DOF
         if (dofMapBC(i, j) == 0) then !BC DOFs
            cycle ! If BC dof is zero, no action needed. If non-zero, fill DeltaEl
         else if (modulo(j, 2) == 0) then !because v is stored in even places and u in odd
            DeltaEl(ceiling(j/2.0d0), 2, i) = DeltaG(dofMapBC(i, j)) !Stores only v for ith element
         else
            DeltaEl(ceiling(j/2.0d0), 1, i) = DeltaG(dofMapBC(i, j)) !Stores only u for ith element
         end if
      end do
   end do
   !----------------------Working Ends----------------------------
   !--------------------------------------------------------------
   do i = 2, nx + 1
      do j = 2, ny + 2

         isBar = .false.
         isCyl = .false.

         ii = i - 1.5
         jj = j - 1.5

         isCyl = ((ii - xc_)**2.0 + (jj - yc_)**2.0)**0.5 .le. 0.5*dia_

         bounDofTopEl(0, 1:2) = 0.0d0
         do k = 1, size(topEl)
            bounDofTopEl(k, 1) = dotProd(psiN(1.0d0, 1.0d0), DeltaEl(:, 1, topEl(k)))
            bounDofTopEl(k, 2) = dotProd(psiN(1.0d0, 1.0d0), DeltaEl(:, 2, topEl(k)))
         end do

         bounDofBottomEl(0, 1:2) = 0.0d0
         do k = 1, size(bottomEl)
            bounDofBottomEl(k, 1) = dotProd(psiN(1.0d0, -1.0d0), DeltaEl(:, 1, bottomEl(k)))
            bounDofBottomEl(k, 2) = dotProd(psiN(1.0d0, -1.0d0), DeltaEl(:, 2, bottomEl(k)))
         end do

         bounDofRightEl(0, 1) = dotProd(psiN(1.0d0, -1.0d0), DeltaEl(:, 1, rightEl(1)))
         bounDofRightEl(0, 2) = dotProd(psiN(1.0d0, -1.0d0), DeltaEl(:, 2, rightEl(1)))
         do k = 1, size(rightEl)
            bounDofRightEl(k, 1) = dotProd(psiN(1.0d0, 1.0d0), DeltaEl(:, 1, rightEl(k)))
            bounDofRightEl(k, 2) = dotProd(psiN(1.0d0, 1.0d0), DeltaEl(:, 2, rightEl(k)))
         end do
         !------------------------------------------------------------
         ! write (*, *) linPieceWiseYonX(-3.0d0, reshape([0.0d0, 1.0d0, 2.0d0, 3.0d0, 10.0d0, 20.0d0, 15.0d0, 0.0d0], [4, 2]))
         ! write (*, *) linPieceWiseXonY(2.8d0, reshape([10.0d0, 20.0d0, 15.0d0, 0.0d0, 0.0d0, 1.0d0, 2.0d0, 3.0d0], [4, 2]))
         ! xIntersect = linPieceWiseIntersect([1.0d0,20.0d0],&
         ! reshape([0.0d0, 1.0d0, 2.0d0, 3.0d0,0.0d0,0.0d0, 10.0d0, 20.0d0, 15.0d0, 0.0d0,0.0d0,10.0d0], [6, 2]))
         ! xIntersect = linPieceWiseIntersect([10.0d0,5.0d0],&
         ! reshape([0.0d0,0.0d0, 1.0d0, 2.0d0, 3.0d0,3.0d0,0.0d0,0.0d0,10.0d0, 10.0d0, 10.0d0, 10.0d0, 0.0d0,0.0d0], [7, 2]))

         ! uniqSortedA = sortAscendUnique([-5.0d0,0.0d0,1.0d0,3.0d0,10.0d0,5.8d0,2.0d0,5.0d0])
         ! uniqSortedA = sortAscendUnique([0.0d0,0.0d0,0.0d0,0.0d0,3.2d0,0.0d0,0.0d0,0.0d0])
         ! uniqSortedA = sortAscendUnique(xIntersect)

         bounDispTopEl = refBounTopEl + bounDofTopEl
         bounDispBottomEl = refBounBottomEl + bounDofBottomEl
         bounDispRightEl = refBounRightEl + bounDofRightEl

         bounCoord(1:size(topEl) + 1, :) = bounDispTopEl
         bounCoord(size(topEl) + size(rightEl) + 1:size(topEl) + 1:-1, :) = bounDispRightEl
         bounCoord(size(topEl) + size(rightEl) + size(bottomEl) + 1:size(topEl) + size(rightEl) + 1:-1, :) = bounDispBottomEl
         bounCoord(size(topEl) + size(rightEl) + size(bottomEl) + 2, :) = bounDispTopEl(0, :) !to close the polygon

         xIntersect = linPieceWiseIntersect([ii, jj], bounCoord)
         uniqSortedA = sortAscendUnique(xIntersect)

         do k = 1, size(uniqSortedA)
            if (ii .le. uniqSortedA(k)) then
               if (mod(k, 2) == 0) then
                  isBar = .true.
               end if
               exit
            end if
         end do

         if (isCyl .or. isBar) then
            isn(i, j) = 1
            ! elseif (j == 1 .or. j == ny + 2) then
            !    isn(i, j) = 2
         else
            isn(i, j) = 0
         end if

      end do
   end do

   isn(:, 1) = 2
   isn(:, ny + 2) = 2

   open (unit=12, file="region.dat")

   write (12, *) "Variables=x,y,region"
   write (12, '(2(a,I5))') "Zone I=", nx, ",J=", ny

   do j = 2, ny + 1
      do i = 2, nx + 1
         write (12, '(2(2X,I5),2X,I3)') i, j, isn(i, j)
      end do
      write (12, *)
   end do
   close (12)
!----------------------------------------------------------------------
   do t_ = 0, time_

      t = t_*Ct

      rhoAvg = d0
      do i = 2, nx + 1
         do j = 2, ny + 1
            tmp1 = d0
            tmp2 = d0
            tmp3 = d0
            do a = 0, q - 1
               tmp1 = tmp1 + f(a, i, j)
               tmp2 = tmp2 + f(a, i, j)*ex(a)
               tmp3 = tmp3 + f(a, i, j)*ey(a)
            end do

            rho(i, j) = tmp1
            ux(i, j) = tmp2/tmp1
            uy(i, j) = tmp3/tmp1
            rhoAvg = rhoAvg + tmp1
         end do
      end do

      rhoAvg = rhoAvg/(nx*ny)
      if (rhoAvg .gt. 10.0d0) then
         write (*, *) 'Code Diverged'
         stop
      end if
!----------------------------------------------------------------------
      do i = 2, nx + 1
         do j = 2, ny + 1
            do a = 0, q - 1
               tmp1 = ux(i, j)*ex(a) + uy(i, j)*ey(a)
               tmp2 = ux(i, j)*ux(i, j) + uy(i, j)*uy(i, j)
               feq = wt(a)*rho(i, j)*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
               ft(a, i, j) = f(a, i, j) - (f(a, i, j) - feq)/tau !collision
            end do
         end do
      end do
!----------------------------------------------------------------------
      do i = 2, nx + 1 !Streaming post-collision
         do j = 2, ny + 1
            do a = 0, q - 1
               ia = i + ex(a)
               ja = j + ey(a)
               !if (ia<1 )        { ia = nx  }
               !if (ia>nx)        { ia = 1   }

               f(a, ia, ja) = ft(a, i, j)
            end do
         end do
      end do
!----------------------------------------------------------------------
      ! do i = 2, nx + 1
      !    j = 2
      !    f(5, i, j) = f(8, i, j - 1)
      !    f(2, i, j) = f(4, i, j - 1)
      !    f(6, i, j) = f(7, i, j - 1)
      ! end do

      ! do i = 2, nx + 1
      !    j = ny + 1
      !    f(8, i, j) = f(5, i, j + 1)
      !    f(4, i, j) = f(2, i, j + 1)
      !    f(7, i, j) = f(6, i, j + 1)
      ! end do

      do j = 2, ny + 1
         i = 2
         uPara_ = 6.0d0*uMean_*(ny - (j - 1.5))*(j - 1.5)/ny**2.0d0; 
         if (t .lt. 2.0d0) then
            uParaRamp_ = uPara_*(1 - cos(pi*t/2.0d0))/2.0d0
         else
            uParaRamp_ = uPara_
         end if
         rho(i, j) = (f(0, i, j) + f(2, i, j) + f(4, i, j) + 2*(f(6, i, j) + f(3, i, j) + f(7, i, j)))/(1 - uParaRamp_)
         f(1, i, j) = f(3, i, j) + ((2.0/3.0)*rho(i, j)*uParaRamp_)
         f(5, i, j) = f(7, i, j) - (0.5*(f(2, i, j) - f(4, i, j))) + ((1.0/6.0)*rho(i, j)*uParaRamp_)
         f(8, i, j) = f(6, i, j) + (0.5*(f(2, i, j) - f(4, i, j))) + ((1.0/6.0)*rho(i, j)*uParaRamp_)
      end do

      do j = 2, ny + 1
         i = nx + 1
         ux(i, j) = (f(0, i, j) + f(2, i, j) + f(4, i, j) + 2*(f(1, i, j) + f(5, i, j) + f(8, i, j)))/rhoF_ - 1
         f(3, i, j) = f(1, i, j) - ((2.0/3.0)*rhoF_*ux(i, j))
         f(6, i, j) = f(8, i, j) - 0.5*(f(2, i, j) - f(4, i, j)) - ((1.0/6.0)*rhoF_*ux(i, j))
         f(7, i, j) = f(5, i, j) + 0.5*(f(2, i, j) - f(4, i, j)) - ((1.0/6.0)*rhoF_*ux(i, j))
      end do
!----------------------------------------------------------------------
      fx(2) = d0
      fy(2) = d0
      do i = 2, nx + 1 !BC
         do j = 2, ny + 1
            if (isn(i, j) .eq. 0) then

               do a = 0, q - 1

                  ia = i + ex(a)
                  ja = j + ey(a)

                  if (isn(ia, ja) .eq. 1) then !structure
                     f(kb(a), i, j) = ft(a, i, j)
                     f(a, ia, ja) = ft(kb(a), ia, ja)
                     fx(2) = fx(2) + ex(a)*2.0*(-ft(kb(a), ia, ja) + ft(a, i, j))
                     fy(2) = fy(2) + ey(a)*2.0*(-ft(kb(a), ia, ja) + ft(a, i, j))
                  end if

                  if (isn(ia, ja) .eq. 2) then !wall
                     f(kb(a), i, j) = ft(a, i, j)
                  end if

               end do
            end if
         end do
      end do
!----------------------------------------------------------------------
      fx_t = 0.5*(fx(1) + fx(2))
      fy_t = 0.5*(fy(1) + fy(2))
      Cd = fx_t*CFor!/(0.5*rhoF_*uMean_*uMean_*dia_)
      Cl = fy_t*CFor!/(0.5*rhoF_*uMean_*uMean_*dia_)
!----------------------------------------------------------------------
      if (mod(t_, dispFreq) .eq. 0) then
         write (10, '(I8,4(3X,F10.6))') t_, t, rhoAvg, Cd, Cl
         !write (*, '(I8,4(3X,F10.6))') ts, ts*Ct, rhoAvg, Cd, Cl

         !write (11, '(I8,6(3X,F10.6))') ts, ux(150, 201), uy(150, 201), ux(200, 250), uy(200, 250), ux(250, 201), uy(250, 201)
      end if
!----------------------------------------------------------------------
      fx(1) = fx(2)
      fy(1) = fy(2)
!----------------------------------------------------------------------
      if (t_ .le. time_ .and. mod(t_, (time_/(noOfSnaps - 1))) .eq. 0) then

         solnumber = solnumber + 1
         write (filename, '(a,i3.3,a)') "../output/snap", solnumber, ".dat"
         open (unit=12, file=filename)

         write (12, *) "Variables=x,y,u,v,rho,region"
         write (12, '(2(a,I5))') "Zone I=", nx, ",J=", ny

         do j = 2, ny + 1
            do i = 2, nx + 1
               write (12, '(2(2X,I5),3(2X,E12.4),2X,I3)') i, j, ux(i, j), uy(i, j), rho(i, j), isn(i, j)
            end do
            write (12, *)
         end do
         close (12)
         write (*, '(a,I3,a,I8)') "snap", solnumber, " recorded at LBM time %d", t_
      end if
!----------------------------------------------------------------------
   end do!Time loop Ends

   close (10)
   ! close (11)

   write (*, *) 'Program ended at :', dateTime()
   write (*, *) '======================================================'

contains

   subroutine setupElemMap(dofMapBC_, coupleRange_)
      implicit none
      ! Arguments declarations
      integer, allocatable, dimension(:, :), intent(out) :: dofMapBC_
      integer, intent(out) :: coupleRange_
      integer, allocatable, dimension(:, :) :: dofMap
      integer, allocatable, dimension(:) :: cntS
      integer :: iEl, nNd, cnt, i, j, k, p, m, n
      !
      !global degEl flag nEl nElx nEly dofBC nDofPerNode nDofPerEl
      !Node Map relating global DOF to local DOF for all elements
      !m2f: dofMap=zeros(nEl,nDofPerEl)
      allocate (dofMap(nEl, nDofPerEl))
      allocate (dofMapBC_(nEl, nDofPerEl))

      cnt = 1
      allocate (cntS(nEl))
      cntS = 0

      do n = 0, nElx - 1
         if (n > 0 .and. n < nElx) then
            cnt = cnt - (nEly*degEl + 1)*nDofPerNode !to account for common nodes along y
         end if
         do i = 0, degEl
            do m = 0, nEly - 1
               if (m > 0 .and. m < nEly) then
                  cnt = cnt - nDofPerNode !to account for common nodes along x
               end if
               do j = 0, degEl
                  iEl = (m + 1) + nEly*n !global element number
                  if (flag == 1 .and. i > 0 .and. i < degEl .and. j > 0 .and. j < degEl) then
                     cntS(iEl) = cntS(iEl) + 1 !count interior nodes for a given element
                     cycle !skip interior node for serendipity element
                  end if
                  nNd = (j + 1) + (degEl + 1)*i - cntS(iEl) !local node number for a given element
                  !dofMap(iEl,[2*nNd-1 2*nNd])=[ cnt,cnt+1 ] !allot global nodal dof
                  dofMap(iEl, 2*nNd - 1) = cnt
                  dofMap(iEl, 2*nNd) = cnt + 1
                  cnt = cnt + nDofPerNode
               end do
            end do
         end do
      end do

      coupleRange_ = dofMap(1, nDofPerEl) - dofMap(1, 1)

      !Replacing BC nodes with zero and renumbering the DOF map(works even for unordered dofBC)
      dofMapBC_ = dofMap
      do k = 1, size(dofBC)!for each BC node
         do i = 1, nEl!for each element
            do j = 1, nDofPerEl!for each DOF
               if (dofMapBC_(i, j) == dofBC(k)) then
                  dofMapBC_(i, j) = 0 !set BC DOF to zero
               end if
               if (dofMapBC_(i, j) > dofBC(k)) then
                  dofMapBC_(i, j) = dofMapBC_(i, j) - 1 !decrement global DOF numbering
               end if
            end do
         end do

         do p = 1, size(dofBC)!helps work with unordered dofBC
            if (dofBC(p) > dofBC(k)) then
               dofBC(p) = dofBC(p) - 1
            end if
         end do

      end do

   end subroutine setupElemMap

end program turekFSI