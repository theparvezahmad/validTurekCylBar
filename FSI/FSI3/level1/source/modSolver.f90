module solver

   use input
   use math
   use fem
   use lbm
   implicit none

   double precision, allocatable, dimension(:, :):: ci
   double precision, allocatable, dimension(:):: wi
   integer, allocatable, dimension(:)::kb

   integer::nx, ny
   double precision:: Clen, Crho, Ct, Cnu, CVel, CFor, tau, invTau
   double precision:: nu_, uMean_, xc_, yc_, dia_, chanL_, barL_, barH_

   double precision, allocatable, dimension(:, :, :):: f

   integer :: nEl, nNodesPerEl, nNodes, nDofPerEl, nDof, nDofBC
   double precision:: lenElx, lenEly, lenElx_, lenEly_

   ! integer:: cnt, k, m, n, p
   integer, allocatable, dimension(:, :) :: dofMapBC
   integer, allocatable, dimension(:) :: dofBC
   integer:: coupleRange

   integer, allocatable, dimension(:)::bottomEl, topEl, leftEl, rightEl
   double precision, allocatable, dimension(:, :)::refBounTopEl, refBounBottomEl, refBounRightEl

contains

   subroutine compFSIresponse()
      implicit none

      integer, allocatable, dimension(:, :)::isn
      double precision, allocatable, dimension(:, :, :):: ft
      double precision, allocatable, dimension(:, :):: ux, uy, rho
      !====================================
      ! integer::a, ia, ja, i, j, k, m, n, p, cnt, t_
      integer:: cFSinteract
      ! double precision:: t, tmp1, tmp2, tmp3, uPara_, uParaRamp_, ii, jj
      double precision:: Cd, Cl, rhoSum, Fx(avgSpan), Fy(avgSpan)
      ! double precision:: Fx(avgSpan), Fy(avgSpan), FxLocal, FyLocal
      ! double precision:: wi(0:q - 1)

      double precision, allocatable, dimension(:, :, :) :: DeltaEl
      integer:: nUdiag, nLdiag, ldab
      ! double precision, allocatable, dimension(:) :: DeltaG
      ! integer, allocatable, dimension(:)::hash

      double precision, allocatable, dimension(:, :)::bounDofTopEl, bounDofBottomEl, bounDofRightEl
      double precision, allocatable, dimension(:, :)::bounDispTopEl, bounDispBottomEl, bounDispRightEl, bounCoord

      double precision, allocatable, dimension(:, :)::surfForce, uniqSurfForce

      double precision, allocatable, dimension(:, :):: PSItPointForce

      ! double precision, allocatable, dimension(:) :: tempSum
      !========================
      double precision :: t, dt, totTime!, tStart_, tEnd_, dt_
      integer:: t_
      integer, dimension(2) :: probeDof
      external :: DGBSV
      !double precision, allocatable, dimension(:,:), intent(out) :: soln

      double precision, allocatable, dimension(:, :):: KgB, MgB, Mtemp
      double precision, allocatable, dimension(:):: Fg
      double precision :: a0, a1, a2, a3, a4, a5, a6, a7
      double precision, allocatable, dimension(:) :: X, Xd, Xdd, XOld, XdOld, XddOld
      double precision, allocatable, dimension(:) :: Xi, Xidd, Xii, Xti

      double precision :: tol, err
      double precision, allocatable, dimension(:, :):: stepLoad
      ! double precision::stepLoad
      integer, allocatable, dimension(:)::pivot
      integer :: iLoad, cntIter, ok, cnt

      character(len=30):: filename
      integer:: solnumber
      solnumber = 0

      allocate (ft(0:q - 1, nx + 2, ny + 2))
      allocate (ux(nx + 2, ny + 2))
      allocate (uy(nx + 2, ny + 2))
      allocate (rho(nx + 2, ny + 2))
      allocate (isn(nx + 2, ny + 2))

      allocate (bounDofTopEl(0:degEl*size(topEl), 2))
      allocate (bounDofBottomEl(0:degEl*size(bottomEl), 2))
      allocate (bounDofRightEl(0:degEl*size(rightEl), 2))

      allocate (bounDispTopEl(0:degEl*size(topEl), 2))
      allocate (bounDispBottomEl(0:degEl*size(bottomEl), 2))
      allocate (bounDispRightEl(0:degEl*size(rightEl), 2))

      allocate (bounCoord(degEl*(size(topEl) + size(bottomEl) + size(rightEl)) + 2, 2))

      dt = Ct
      a0 = 1.0d0/(alpha*dt**2.0d0)
      a1 = beta/(alpha*dt)
      a2 = 1.0d0/(alpha*dt)
      a3 = (1.0d0/(2.0d0*alpha)) - 1.0d0
      a4 = (beta/alpha) - 1.0d0
      a5 = (dt/2.0d0)*((beta/alpha) - 2.0d0)
      a6 = dt*(1.0d0 - beta)
      a7 = beta*dt

      allocate (Mtemp(nDofBC, nDofBC))

      allocate (pivot(nDofBC))
      allocate (X(nDofBC), Xd(nDofBC), Xdd(nDofBC))
      allocate (XOld(nDofBC), XdOld(nDofBC), XddOld(nDofBC))
      allocate (Xi(nDofBC), Xii(nDofBC), Xti(nDofBC), Xidd(nDof))

      allocate (PSItPointForce(nDofPerEl, nEl))
      allocate (stepLoad(nDofPerEl, nEl))

      XOld = 0.0d0
      XdOld = 0.0d0
      XddOld = 0.0d0

      !cnt=1
      !nCycle=11
      !totSnaps=(tEnd-tStart)/dt+1

      ! allocate(DeltaConvG(noOfLoadSteps,3))
      ! DeltaConvG = 0.0d0

      !totT0=totalLoad !Traction(per unit area)

      X = 0.0d0
      ! t = tStart_

      Fx = d0
      Fy = d0
      open (unit=10, file="../output/Fluid_tRhoCdCl.dat")
      write (10, *) "Variables=timeLBM,timeReal,rho,Cd,Cl"
      open (UNIT=11, file='../output/Solid_tUxUy.dat')
      write (11, *) "Variables=timeReal,ux,uy"

      nUdiag = coupleRange
      nLdiag = coupleRange
      ldab = 2*nLdiag + nUdiag + 1
      ! topLeftPt = nDofBC - (nElx*degEl + 1)*nDofPerNode + 1
      probeDof = nDofBC - (nEly/2*degEl + 1)*nDofPerNode + [1, 2]
      !----------------------------------------------------------------------

      t = 0.0d0!tStart
      t_ = 0
      totTime = totTime_*Ct

      do while (t .lt. totTime)

         t_ = t_ + 1
         t = t_*dt

         call calcMacroVarLBM()

         ! rhoAvg = sum(rho)/(nx*ny)
         if (rhoSum/(nx*ny) .gt. 10.0d0) then
            write (*, *) 'Code Diverged'
            stop
         end if

         call collide()
         call stream()
         call applyInletOutletBC2()
         !----------------------------------------------------------------------
         call applyObjWallBC_calcForceObj()!(cFSinteract, surfForce)
         ! allocate(surfForce(6,4))
         ! surfForce = reshape([1, 3, 5, 6,&
         !                      1, 3, 3, 9,&
         !                      1, 3, 8, 7,&
         !                      1, 3, 4, 5,&
         !                      1, 3, 10, 11,&
         !                      1, 3, 20, 30], [6, 4], order=[2, 1])
         call addDuplicateFields()!(surfForce, uniqSurfForce)
         deallocate (surfForce)
         ! do i = 1, size(uniqSurfForce,1)
         !    uniqSurfForce(i,3)=1
         !    uniqSurfForce(i,4)=3*i
         !    ! write(*,*) uniqSurfForce(i,1),uniqSurfForce(i,2)
         ! end do

         call distSurfForce2Elem()!(uniqSurfForce, PSItPointForce)
         deallocate (uniqSurfForce)
         ! call calcMgKgFg(PSItPointForce)
         ! write(*,*) tempSum
         ! write(*,*) sum(tempSum),sum(uniqSurfForce(:,3))
         ! stop
         !==========================compDynRes=========================
         do iLoad = 1, noOfLoadSteps
            !loadVec=linspace(totalLoad/noOfLoadSteps,totalLoad,noOfLoadSteps)
            !stepLoad=loadVec(iLoad)
            ! stepLoad = totalLoad*iLoad/noOfLoadSteps
            stepLoad = PSItPointForce*iLoad/noOfLoadSteps

            !Xi=zeros(nDofBC,1);
            call calcMgKgFg(X, stepLoad, MgB, KgB, Fg)
            !call calcMgKgFg(X, stepLoad)
            !Xi=(a0*Mg + Kg)\(Fg + Mg*(a2*XdOld + a3*XddOld))
            Xi = Fg + mulMatBVec(MgB, a2*XdOld + a3*XddOld, nUdiag, nLdiag)
            Mtemp = a0*MgB + KgB
            !call DGESV(nDofBC, 1, Mtemp, nDofBC, pivot, Xi, nDofBC, ok)
            call DGBSV(nDofBC, nLdiag, nUdiag, 1, Mtemp, ldab, pivot, Xi, nDofBC, ok)

            err = 1
            tol = 1e-8
            cntIter = 0

            do while (err > tol)
               cntIter = cntIter + 1
               Xidd = a0*Xi - a2*XdOld - a3*XddOld
               Xti = X + Xi
               call calcMgKgFg(Xti, stepLoad, MgB, KgB, Fg)
               !Xii=(a0*Mg + Kg)\(Fg - Mg*Xidd)
               Xii = Fg - mulMatBVec(MgB, Xidd, nUdiag, nLdiag)
               Mtemp = a0*MgB + KgB
               !call DGESV(nDofBC, 1, Mtemp, nDofBC, pivot, Xii, nDofBC, ok)
               call DGBSV(nDofBC, nLdiag, nUdiag, 1, Mtemp, ldab, pivot, Xii, nDofBC, ok)
               Xi = Xi + Xii
               err = norm2(Xii)/norm2(Xi + X)
               if (cntIter == 10) exit
            end do !iter loop ends

            if (modulo(int(t/dt), 20) == 0) then
               write (*, *) t, err, cntIter!, dateTime()
            end if

            ! Res=(Mg*Xidd + Kg*Xii - Fg);!Residual
            ! DeltaConvG(iLoad,1:3)=[Xii(topLeftPt) Res(topLeftPt) cntIter]; !Converged solution for the last load step

         end do !iLoad loop ends

         Xdd = a0*Xi - a2*XdOld - a3*XddOld
         Xd = XdOld + a6*XddOld + a7*Xdd
         X = X + Xi

         XOld = X
         XdOld = Xd
         XddOld = Xdd

         ! !Analytical solution
         ! Xa=inv(V)*inv(Mg_)*Fg_./Freq.^2.*(1-cos(Freq*t));!constant load
         ! Xa=V*Xa;

         !Recording time history
         ! soln(i,1)=t !/funT;
         ! do iProbeDof=1,size(probeDof)
         !   soln(i,iProbeDof+1)=X(probeDof(iProbeDof))
         ! enddo
         !==========================compDynRes=========================
         !----------------------------------------------------------------------
         !----------------------------------------------------------------------
         if (mod(t_, dispFreq) .eq. 0) then
            write (10, '(I8,4(3X,F10.6))') t_, t, rhoSum/(nx*ny), Cd, Cl
            !write (*, '(I8,4(3X,F10.6))') ts, ts*Ct, rhoAvg, Cd, Cl

            !write (11, '(I8,6(3X,F10.6))') ts, ux(150, 201), uy(150, 201), ux(200, 250), uy(200, 250), ux(250, 201), uy(250, 201)
         end if
         !----------------------------------------------------------------------
         if (t_ .le. totTime_ .and. mod(t_, (totTime_/(noOfSnaps - 1))) .eq. 0) then

            solnumber = solnumber + 1
            write (filename, '(a,i3.3,a)') "../output/snap", solnumber, ".dat"
            call writeSoln(filename)
            write (*, '(a,I3,a,I8)') "snap", solnumber, " recorded at LBM time %d", t_
         end if
         !----------------------------------------------------------------------
      end do!Time loop Ends

      close (10)
      close (11)

   contains
      subroutine calcMacroVarLBM()
         implicit none

         integer:: a, i, j
         double precision::tmp1, tmp2, tmp3, tmp4
         rhoSum = 0.0d0

         do j = 2, ny + 1
            do i = 2, nx + 1
               tmp1 = d0
               tmp2 = d0
               tmp3 = d0
               do a = 0, q - 1
                  tmp1 = tmp1 + f(a, i, j)
                  tmp2 = tmp2 + f(a, i, j)*ci(a, 1)
                  tmp3 = tmp3 + f(a, i, j)*ci(a, 2)
               end do

               rho(i, j) = tmp1
               tmp4 = 1.0d0/tmp1
               ux(i, j) = tmp2*tmp4!tmp2/tmp1
               uy(i, j) = tmp3*tmp4!tmp3/tmp1
               rhoSum = rhoSum + tmp1
            end do
         end do
      end subroutine calcMacroVarLBM

      subroutine collide()
         implicit none

         integer::a, i, j
         double precision::tmp1, tmp2, feq

         do j = 2, ny + 1
            do i = 2, nx + 1
               do a = 0, q - 1
                  tmp1 = ux(i, j)*ci(a, 1) + uy(i, j)*ci(a, 2)
                  tmp2 = ux(i, j)*ux(i, j) + uy(i, j)*uy(i, j)
                  feq = wi(a)*rho(i, j)*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
                  ft(a, i, j) = f(a, i, j) - (f(a, i, j) - feq)*invTau !collision
               end do
            end do
         end do
      end subroutine collide

      subroutine stream()
         implicit none

         integer::a, i, j, ia, ja

         do j = 2, ny + 1
            do i = 2, nx + 1 !Streaming post-collision
               do a = 0, q - 1
                  ia = i + int(ci(a, 1))
                  ja = j + int(ci(a, 2))
                  !if (ia<1 )        { ia = nx  }
                  !if (ia>nx)        { ia = 1   }

                  f(a, ia, ja) = ft(a, i, j)
               end do
            end do
         end do
      end subroutine stream

      subroutine applyInletOutletBC2()
         implicit none

         integer::i, j
         double precision::uPara_, uParaRamp_, tmp1, tmp2
         tmp1 = 1.0d0/(ny**2.0d0)

         do j = 2, ny + 1
            i = 2
            uPara_ = 6.0d0*uMean_*(ny - (j - 1.5))*(j - 1.5)*tmp1; 
            if (t .lt. 2.0d0) then
               uParaRamp_ = uPara_*(1 - cos(pi*t*0.5d0))*0.5d0
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

      end subroutine applyInletOutletBC2

      subroutine detectDeformedShape()!(isn)
         implicit none

         ! integer, dimension(:, :), intent(out) :: isn
         integer::i, j, k, m
         logical::isCyl, isBar
         double precision::ii, jj
         double precision, allocatable, dimension(:)::xIntersect, uniqSortedA

         call mapGlobal2Local(X, DeltaEl)
         isn = 0
         do i = 2, nx + 1
            do j = 2, ny + 2

               isBar = .false.
               isCyl = .false.

               ii = i! - 1.5
               jj = j! - 1.5

               isCyl = ((ii - xc_)**2.0 + (jj - yc_)**2.0)**0.5 .le. 0.5*dia_

               bounDofTopEl(0, 1:2) = 0.0d0
               do k = 1, size(topEl)
                  do m = 0, degEl - 1
                     bounDofTopEl(degEl*k - m, 1) = dotProd(psiN(1.0d0 - dble(2*m)/degEl, 1.0d0), DeltaEl(:, 1, topEl(k)))
                     bounDofTopEl(degEl*k - m, 2) = dotProd(psiN(1.0d0 - dble(2*m)/degEl, 1.0d0), DeltaEl(:, 2, topEl(k)))
                  end do
                  ! bounDofTopEl(k, 1) = dotProd(psiN(1.0d0, 1.0d0), DeltaEl(:, 1, topEl(k)))
                  ! bounDofTopEl(k, 2) = dotProd(psiN(1.0d0, 1.0d0), DeltaEl(:, 2, topEl(k)))
               end do

               bounDofBottomEl(0, 1:2) = 0.0d0
               do k = 1, size(bottomEl)
                  do m = 0, degEl - 1
                     bounDofBottomEl(degEl*k - m, 1) = dotProd(psiN(1.0d0 - dble(2*m)/degEl, -1.0d0), DeltaEl(:, 1, bottomEl(k)))
                     bounDofBottomEl(degEl*k - m, 2) = dotProd(psiN(1.0d0 - dble(2*m)/degEl, -1.0d0), DeltaEl(:, 2, bottomEl(k)))
                  end do
                  ! bounDofBottomEl(k, 1) = dotProd(psiN(1.0d0, -1.0d0), DeltaEl(:, 1, bottomEl(k)))
                  ! bounDofBottomEl(k, 2) = dotProd(psiN(1.0d0, -1.0d0), DeltaEl(:, 2, bottomEl(k)))
               end do

               bounDofRightEl(0, 1) = dotProd(psiN(1.0d0, -1.0d0), DeltaEl(:, 1, rightEl(1)))
               bounDofRightEl(0, 2) = dotProd(psiN(1.0d0, -1.0d0), DeltaEl(:, 2, rightEl(1)))
               do k = 1, size(rightEl)
                  do m = 0, degEl - 1
                     bounDofRightEl(degEl*k - m, 1) = dotProd(psiN(1.0d0, 1.0d0 - dble(2*m)/degEl), DeltaEl(:, 1, rightEl(k)))
                     bounDofRightEl(degEl*k - m, 2) = dotProd(psiN(1.0d0, 1.0d0 - dble(2*m)/degEl), DeltaEl(:, 2, rightEl(k)))
                  end do
                  ! bounDofRightEl(k, 1) = dotProd(psiN(1.0d0, 1.0d0), DeltaEl(:, 1, rightEl(k)))
                  ! bounDofRightEl(k, 2) = dotProd(psiN(1.0d0, 1.0d0), DeltaEl(:, 2, rightEl(k)))
               end do
               !------------------------------------------------------------
               ! write (*, *) linPieceWiseYonX(-3.0d0, reshape([0.0d0, 1.0d0, 2.0d0, 3.0d0, 10.0d0, 20.0d0, 15.0d0, 0.0d0], [4, 2]))
               ! write (*, *) linPieceWiseXonY(2.8d0, reshape([10.0d0, 20.0d0, 15.0d0, 0.0d0, 0.0d0, 1.0d0, 2.0d0, 3.0d0], [4, 2]))
               ! xIntersect = linPieceWiseIntersect([1.0d0,20.0d0],&
               ! reshape([0.0d0, 1.0d0, 2.0d0, 3.0d0,0.0d0,0.0d0, 10.0d0, 20.0d0, 15.0d0, 0.0d0,0.0d0,10.0d0], [6, 2]))
               ! xIntersect = quadPieceWiseIntersect([1.5d0,5.0d0],&
               ! reshape([0.0d0,0.0d0, 1.0d0, 2.0d0, 3.0d0,3.0d0,0.0d0,0.0d0,10.0d0, 10.0d0, 10.0d0, 10.0d0, 0.0d0,0.0d0], [7, 2]))
               ! write(*,*) quadPieceWiseXonY(1.5,&
               ! reshape([0.0d0,0.0d0, 1.0d0, 2.0d0, 3.0d0,3.0d0,0.0d0,0.0d0,10.0d0, 10.0d0, 10.0d0, 10.0d0, 0.0d0,0.0d0], [7, 2]))

               ! uniqSortedA = sortAscendUnique([-5.0d0,0.0d0,1.0d0,3.0d0,10.0d0,5.8d0,2.0d0,5.0d0])
               ! uniqSortedA = sortAscendUnique([0.0d0,0.0d0,0.0d0,0.0d0,3.2d0,0.0d0,0.0d0,0.0d0])
               ! uniqSortedA = sortAscendUnique(xIntersect)

               bounDispTopEl = refBounTopEl + bounDofTopEl
               bounDispBottomEl = refBounBottomEl + bounDofBottomEl
               bounDispRightEl = refBounRightEl + bounDofRightEl

               bounCoord(1:degEl*size(topEl) + 1, :) = bounDispTopEl
               bounCoord(degEl*(size(topEl) + size(rightEl)) + 1:degEl*size(topEl) + 1:-1, :) = bounDispRightEl
               bounCoord(degEl*(size(topEl) + size(rightEl) + size(bottomEl)) + 1:degEl*(size(topEl) + size(rightEl)) + 1:-1, :) = &
                  bounDispBottomEl
               bounCoord(degEl*(size(topEl) + size(rightEl) + size(bottomEl)) + 2, :) = bounDispTopEl(0, :) !to close the polygon
               ! bounCoord(degEl*(size(topEl)+size(rightEl)+size(bottomEl))+3,:)=0.5d0*(bounDispTopEl(0,:) + bounDispBottomEl(0,:)) !to close the polygon

               xIntersect = linPieceWiseIntersect([ii, jj], bounCoord)
               ! xIntersect=quadPieceWiseIntersect([ii,jj],bounCoord)
               uniqSortedA = sortAscendUnique(xIntersect)

               do k = 1, size(uniqSortedA)
                  if (ii .le. uniqSortedA(k)) then
                     if (mod(k, 2) == 0) then
                        isBar = .true.
                     end if
                     exit
                  end if
               end do

               if (isCyl) then
                  isn(i, j) = 1
               end if

               if (isBar) then
                  isn(i, j) = 2
               end if

            end do
         end do

         isn(:, 1) = 3
         isn(:, ny + 2) = 3

         open (unit=12, file="../output/region.dat")

         write (12, *) "Variables=x,y,region"
         write (12, '(2(a,I5))') "Zone I=", nx, ",J=", ny

         do j = 2, ny + 1
            do i = 2, nx + 1
               write (12, '(2(2X,I5),2X,I3)') i, j, isn(i, j)
            end do
            write (12, *)
         end do
         close (12)
         ! stop
      end subroutine detectDeformedShape

      subroutine applyObjWallBC_calcForceObj()
         implicit none

         ! double precision, allocatable, dimension(:, :) :: surfForce
         integer:: i, j, a, ia, ja!, cFSinteract
         ! double precision, allocatable, dimension(:, :) :: surfForce_
         double precision:: tmp1, tmp2, Fx_t, Fy_t

         ! allocate (surfForce_(int(2*barL_ + barH_)*4, 4))

         Fx(avgSpan) = d0
         Fy(avgSpan) = d0
         cFSinteract = 0

         do j = 2, ny + 1
            do i = 2, nx + 1 !BC
               if (isn(i, j) .eq. 0) then

                  do a = 0, q - 1

                     ia = i + int(ci(a, 1))
                     ja = j + int(ci(a, 2))

                     if (isn(ia, ja) .eq. 1) then !cylinder
                        f(kb(a), i, j) = ft(a, i, j)
                        f(a, ia, ja) = ft(kb(a), ia, ja)

                        tmp1 = ci(a, 1)*2.0*(-ft(kb(a), ia, ja) + ft(a, i, j))
                        Fx(avgSpan) = Fx(avgSpan) + tmp1

                        tmp2 = ci(a, 2)*2.0*(-ft(kb(a), ia, ja) + ft(a, i, j))
                        Fy(avgSpan) = Fy(avgSpan) + tmp2
                     end if

                     if (isn(ia, ja) .eq. 2) then !elastic bar
                        f(kb(a), i, j) = ft(a, i, j)
                        f(a, ia, ja) = ft(kb(a), ia, ja)

                        tmp1 = ci(a, 1)*2.0*(-ft(kb(a), ia, ja) + ft(a, i, j))
                        Fx(avgSpan) = Fx(avgSpan) + tmp1

                        tmp2 = ci(a, 2)*2.0*(-ft(kb(a), ia, ja) + ft(a, i, j))
                        Fy(avgSpan) = Fy(avgSpan) + tmp2

                        ! cFSinteract = cFSinteract + 1
                        ! ! surfForce_(cFSinteract, :) = [0.5*(i + ia) - 1.5d0, 0.5*(j + ja) - 1.5d0, tmp1, tmp2]
                        ! surfForce_(cFSinteract, :) = [0.5d0*(i + ia), 0.5d0*(j + ja), tmp1, tmp2]
                     end if

                     if (isn(ia, ja) .eq. 3) then !wall
                        f(kb(a), i, j) = ft(a, i, j)
                     end if

                  end do
               end if
            end do
         end do

         ! allocate (surfForce(cFSinteract, 4))
         ! surfForce = surfForce_(1:cFSinteract, :)

         Fx_t = sum(Fx)/min(t_ + 1, avgSpan)
         Fy_t = sum(Fy)/min(t_ + 1, avgSpan)
         Cd = Fx_t*CFor!/(0.5*rhoF_*uMean_*uMean_*dia_)
         Cl = Fy_t*CFor!/(0.5*rhoF_*uMean_*uMean_*dia_)
         !----------------------------------------------------------------------
         Fx(1:avgSpan - 1) = Fx(2:avgSpan)
         Fy(1:avgSpan - 1) = Fy(2:avgSpan)

      end subroutine applyObjWallBC_calcForceObj

      subroutine addDuplicateFields()!(surfForce, uniqSurfForce)
         implicit none

         ! double precision, dimension(:, :), intent(in) :: surfForce
         ! double precision, allocatable, dimension(:, :), intent(out) :: uniqSurfForce

         double precision, dimension(size(surfForce, 1), size(surfForce, 2)) :: uniqSurfForce_
         integer, allocatable, dimension(:)::hash
         integer:: i, k, j, cntUniq, cntHash

         cntUniq = 0
         cntHash = 0
         allocate (hash(size(surfForce, 1) - 1))

         outer: do i = 1, size(surfForce, 1)

            do k = 1, cntHash
               if (i .eq. hash(k)) then
                  cycle outer
               end if
            end do

            cntUniq = cntUniq + 1
            uniqSurfForce_(cntUniq, :) = surfForce(i, :)
            inner: do j = i, size(surfForce, 1)

               if (j == i) then
                  cycle
               end if

               do k = 1, cntHash
                  if (j .eq. hash(k)) then
                     cycle inner
                  end if
               end do

               if (surfForce(i, 1) .eq. surfForce(j, 1) .and. surfForce(i, 2) .eq. surfForce(j, 2)) then
                  cntHash = cntHash + 1
                  hash(cntHash) = j
                  uniqSurfForce_(cntUniq, 3:4) = uniqSurfForce_(cntUniq, 3:4) + surfForce(j, 3:4)
               end if

            end do inner
         end do outer

         allocate (uniqSurfForce(cntUniq, size(surfForce, 2)))
         uniqSurfForce = uniqSurfForce_(1:cntUniq, :)

      end subroutine addDuplicateFields

      subroutine distSurfForce2Elem()!(uniqSurfForce, PSItPointForce)
         implicit none
         ! input: uniqSurfForce(N*4) N(=no of unique force interaction points) X [x y Fx Fy]
         ! output: PSItPointForce(nEl) point force integral for each element
         ! double precision, dimension(:, :), intent(in) :: uniqSurfForce
         ! double precision, allocatable, dimension(:, :), intent(out) :: PSItPointForce
         ! double precision,allocatable,dimension(:),intent(out) :: tempSum
         integer::iUniqSurfForce, k, i
         double precision:: x1, x2, y1, y2, x, y, xNat, yNat
         double precision, allocatable, dimension(:, :) :: PSI, PSIt
         double precision, allocatable, dimension(:) :: psiNEval

         allocate (PSI(2, nDofPerEl))
         allocate (PSIt(nDofPerEl, 2))
         allocate (psiNEval(nNodesPerEl))

         PSItPointForce = d0
         ! allocate(tempSum(nEl))
         ! tempSum=d0
         ! write(*,*) size(tempSum)

         do i = 1, size(topEl)
            x1 = bounDispTopEl(degEl*(i - 1), 1)
            x2 = bounDispTopEl(degEl*i, 1)
            y1 = bounDispTopEl(degEl*(i - 1), 2)
            y2 = bounDispTopEl(degEl*i, 2)

            do iUniqSurfForce = 1, size(uniqSurfForce, 1)
               x = uniqSurfForce(iUniqSurfForce, 1)
               y = uniqSurfForce(iUniqSurfForce, 2)

               ! if ( x .ge. x1 .and. x .lt. x2 .and. y .gt. (min(y1,y2)-0.5d0) .and. y .lt. (max(y1,y2)+0.5d0)) then
               if (x .ge. x1 .and. x .lt. x2 .and. y .gt. (y1 - 0.5d0) .and. y .lt. (y2 + 0.5d0)) then

                  xNat = (2.0d0*x - (x1 + x2))/(x2 - x1)
                  psiNEval = psiN(xNat, 1.0d0)

                  do k = 1, nNodesPerEl
                     PSI(1, 2*k - 1) = psiNEval(k)
                     PSI(2, 2*k) = psiNEval(k)
                  end do
                  PSIt = transpose(PSI)

                  PSItPointForce(:, topEl(i)) = PSItPointForce(:, topEl(i)) + mulMatVec(PSIt, uniqSurfForce(iUniqSurfForce, 3:4))
                  ! tempSum(topEl(i)) = tempSum(topEl(i)) + uniqSurfForce(iUniqSurfForce,3)
               end if
            end do
         end do

         do i = 1, size(bottomEl)
            x1 = bounDispBottomEl(degEl*(i - 1), 1)
            x2 = bounDispBottomEl(degEl*i, 1)
            y1 = bounDispBottomEl(degEl*(i - 1), 2)
            y2 = bounDispBottomEl(degEl*i, 2)

            do iUniqSurfForce = 1, size(uniqSurfForce, 1)
               x = uniqSurfForce(iUniqSurfForce, 1)
               y = uniqSurfForce(iUniqSurfForce, 2)

               ! if ( x .ge. x1 .and. x .lt. x2 .and. y .gt. (min(y1,y2)-0.5d0) .and. y .lt. (max(y1,y2)+0.5d0)) then
               if (x .ge. x1 .and. x .lt. x2 .and. y .gt. (y1 - 0.5d0) .and. y .lt. (y2 + 0.5d0)) then

                  xNat = (2.0d0*x - (x1 + x2))/(x2 - x1)
                  psiNEval = psiN(xNat, -1.0d0)

                  do k = 1, nNodesPerEl
                     PSI(1, 2*k - 1) = psiNEval(k)
                     PSI(2, 2*k) = psiNEval(k)
                  end do
                  PSIt = transpose(PSI)

               PSItPointForce(:, bottomEl(i)) = PSItPointForce(:, bottomEl(i)) + mulMatVec(PSIt, uniqSurfForce(iUniqSurfForce, 3:4))
                  ! tempSum(bottomEl(i)) = tempSum(bottomEl(i)) + uniqSurfForce(iUniqSurfForce,3)
               end if
            end do
         end do

         do i = 1, size(rightEl)
            x1 = bounDispRightEl(degEl*(i - 1), 1)
            x2 = bounDispRightEl(degEl*i, 1)
            y1 = bounDispRightEl(degEl*(i - 1), 2)
            y2 = bounDispRightEl(degEl*i, 2)

            do iUniqSurfForce = 1, size(uniqSurfForce, 1)
               x = uniqSurfForce(iUniqSurfForce, 1)
               y = uniqSurfForce(iUniqSurfForce, 2)

               ! if ( x .gt. (min(x1,x2)-0.5d0) .and. x .lt. (max(x1,x2)+0.5d0) .and. y .ge. y1 .and. y .lt. y2) then
               if (x .gt. (x1 - 0.5d0) .and. x .lt. (x2 + 0.5d0) .and. y .ge. y1 .and. y .lt. y2) then

                  yNat = (2.0d0*y - (y1 + y2))/(y2 - y1)
                  psiNEval = psiN(1.0d0, yNat)

                  do k = 1, nNodesPerEl
                     PSI(1, 2*k - 1) = psiNEval(k)
                     PSI(2, 2*k) = psiNEval(k)
                  end do
                  PSIt = transpose(PSI)

                 PSItPointForce(:, rightEl(i)) = PSItPointForce(:, rightEl(i)) + mulMatVec(PSIt, uniqSurfForce(iUniqSurfForce, 3:4))
                  ! tempSum(rightEl(i)) = tempSum(rightEl(i)) + uniqSurfForce(iUniqSurfForce,3)
               end if
            end do
         end do

      end subroutine distSurfForce2Elem

      subroutine calcMgKgFg(DeltaG, stepLoad, MgB, KgB, Fg)
         !subroutine calcMgKgFg(DeltaG, stepLoad)
         implicit none

         double precision, dimension(:), intent(in) :: DeltaG
         double precision, dimension(:, :), intent(in) :: stepLoad
         ! double precision:: stepLoad
         double precision, allocatable, dimension(:, :), intent(out) :: KgB, MgB
         double precision, allocatable, dimension(:), intent(out) :: Fg

         double precision, allocatable, dimension(:, :) :: BL, BL0, BLt, BLtCBL, BLuv
         double precision, allocatable, dimension(:, :) :: BNL, BNLt, BNLtCBNL
         double precision, allocatable, dimension(:, :, :) :: DeltaEl ! m2f:check dim(nNod
         double precision, allocatable, dimension(:) :: psiNGP, dpsiNdxGP, dpsiNdyGP
         double precision, allocatable, dimension(:) :: F01, F02, Fe, fPSI, BLtS01
         double precision, allocatable, dimension(:, :) :: Ke, KL, KNL, Me ! m2f:check dim(nDof
         double precision, allocatable, dimension(:, :) :: PSI, PSIt, PSItPSI
         double precision :: dudx, dudy, dvdx, dvdy
         double precision, dimension(3, 3) :: C
         double precision, dimension(3) :: Pt, Wt
         double precision, dimension(3) :: E01, S01 ! m2f:check dim(3,1)
         double precision, dimension(4, 4) :: S_Mat

         double precision :: tmp
         integer :: i, j, k, iEl, ii, jj

         allocate (MgB(ldab, nDofBC))
         allocate (KgB(ldab, nDofBC))
         allocate (Fg(nDofBC))

         allocate (BL(3, nDofPerEl))
         allocate (BL0(3, nDofPerEl))
         allocate (BLuv(3, nDofPerEl))
         allocate (BLt(nDofPerEl, 3))
         allocate (BLtCBL(nDofPerEl, nDofPerEl))
         allocate (KL(nDofPerEl, nDofPerEl))
         allocate (BNL(4, nDofPerEl))
         allocate (BNLt(nDofPerEl, 4))
         allocate (BNLtCBNL(nDofPerEl, nDofPerEl))
         allocate (KNL(nDofPerEl, nDofPerEl))
         allocate (Ke(nDofPerEl, nDofPerEl))

         allocate (BLtS01(nDofPerEl))
         allocate (F01(nDofPerEl))
         allocate (F02(nDofPerEl))
         allocate (fPSI(nDofPerEl))
         allocate (Fe(nDofPerEl))

         allocate (PSI(2, nDofPerEl))
         allocate (PSIt(nDofPerEl, 2))
         allocate (PSItPSI(nDofPerEl, nDofPerEl))
         allocate (Me(nDofPerEl, nDofPerEl))

         allocate (psiNGP(nNodesPerEl))
         allocate (dpsiNdxGP(nNodesPerEl))
         allocate (dpsiNdyGP(nNodesPerEl))

         allocate (DeltaEl(nNodesPerEl, nDofPerNode, nEl))

         Wt = [5.0d0/9.0d0, 8.0d0/9.0d0, 5.0d0/9.0d0] !Gauss Weights
         Pt = [-sqrt(3.0d0/5.0d0), 0.0d0, sqrt(3.0d0/5.0d0)] !Gauss Points
         tmp = 0.25d0*lenElx*lenEly

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
         ! !Constant shear BC on top surface for Quadratic Lagrange
         ! for i=1:nElx
         ! F02(nDofPerEl-nDofPerNode*(degEl+1)+1:nDofPerEl,nEl-i+1)=Lz*a/6*T0*[1,0,4,0,1,0]';
         ! end
         ! ! F02(13:18,nEl-3)=Lz*a/6*T0*[1,0,4,0,1,0]';
         ! ! F02(13:18,nEl-2)=Lz*a/6*T0*[1,0,4,0,1,0]';
         ! ! F02(13:18,nEl-1)=Lz*a/6*T0*[1,0,4,0,1,0]';
         ! ! F02(13:18,nEl )=Lz*a/6*T0*[1,0,4,0,1,0]';

         KgB = 0.0d0
         MgB = 0.0d0
         Fg = 0.0d0

         do iEl = 1, nEl ! spans element
            KL = 0.0d0
            KNL = 0.0d0
            F01 = 0.0d0
            F02 = 0.0d0
            Me = 0.0d0
            !iEl=1;
            do j = 1, 3 ! spans x- Gauss points
               do i = 1, 3 ! spans y- Gauss points
                  ! psiGP =psi (XXi(Pt(i)),YEta(Pt(j)));
                  ! dpsidxGP=dpsidx(XXi(Pt(i)),YEta(Pt(j)));
                  ! dpsidyGP=dpsidy(XXi(Pt(i)),YEta(Pt(j)));

                  psiNGP = psiN(Pt(i), Pt(j)) !*GP refers to vars post evaluation at Gauss Points
                  dpsiNdxGP = dpsiNdXi(Pt(i), Pt(j))*2/lenElx !=dpsiNdXi*dXidx
                  dpsiNdyGP = dpsiNdEta(Pt(i), Pt(j))*2/lenEly !=dpsiNdEta*dEtady

                  dudx = dotProd(dpsiNdxGP, DeltaEl(:, 1, iEl)) !Uses converged u from last load step
                  dudy = dotProd(dpsiNdyGP, DeltaEl(:, 1, iEl))
                  dvdx = dotProd(dpsiNdxGP, DeltaEl(:, 2, iEl))
                  dvdy = dotProd(dpsiNdyGP, DeltaEl(:, 2, iEl))

                  PSI = 0.0d0
                  BL0 = 0.0d0
                  BLuv = 0.0d0
                  BNL = 0.0d0

                  do k = 1, nNodesPerEl
                     PSI(1, 2*k - 1) = psiNGP(k)
                     PSI(2, 2*k) = psiNGP(k)

                     BL0(1, 2*k - 1) = dpsiNdxGP(k)
                     BL0(2, 2*k) = dpsiNdyGP(k)
                     BL0(3, 2*k - 1) = dpsiNdyGP(k)
                     BL0(3, 2*k) = dpsiNdxGP(k)

                     BLuv(1, 2*k - 1) = dudx*dpsiNdxGP(k)
                     BLuv(2, 2*k - 1) = dudy*dpsiNdyGP(k)
                     BLuv(3, 2*k - 1) = dudx*dpsiNdyGP(k) + dudy*dpsiNdxGP(k)

                     BLuv(1, 2*k) = dvdx*dpsiNdxGP(k)
                     BLuv(2, 2*k) = dvdy*dpsiNdyGP(k)
                     BLuv(3, 2*k) = dvdx*dpsiNdyGP(k) + dvdy*dpsiNdxGP(k)

                     BNL(1, 2*k - 1) = dpsiNdxGP(k)
                     BNL(2, 2*k - 1) = dpsiNdyGP(k)
                     BNL(3, 2*k) = dpsiNdxGP(k)
                     BNL(4, 2*k) = dpsiNdyGP(k)
                  end do
                  BL = BL0 + BLuv

                  E01(1) = dudx + 0.0d0 + 0.5d0*(dudx*dudx + dvdx*dvdx)
                  E01(2) = dvdy + 0.0d0 + 0.5d0*(dudy*dudy + dvdy*dvdy)
                  E01(3) = dudy + dvdx + (dudx*dudy + dvdx*dvdy)

                  !C=E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2];
                  C = E/((1 + v)*(1 - 2*v))*reshape([1 - v, v, 0.0d0, v, 1 - v, 0.0d0, 0.0d0, 0.0d0, (1 - 2*v)/2], &
                                                    [3, 3], order=[2, 1])

                  S01 = mulMatVec(C, E01)

                  S_Mat = reshape([S01(1), S01(3), 0.0d0, 0.0d0, S01(3), S01(2), 0.0d0, 0.0d0, &
                                   0.0d0, 0.0d0, S01(1), S01(3), 0.0d0, 0.0d0, S01(3), S01(2)], [4, 4], order=[2, 1])

                  BLt = transpose(BL)
                  BLtCBL = mulMat(BLt, mulMat(C, BL)) !Integrand for KL
                  KL = KL + tmp*Wt(j)*Wt(i)*BLtCBL !Summing up all Gauss Points

                  BNLt = transpose(BNL)
                  BNLtCBNL = mulMat(BNLt, mulMat(S_Mat, BNL)) !Integrand for KNL
                  KNL = KNL + tmp*Wt(j)*Wt(i)*BNLtCBNL !Summing up all Gauss Points

                  BLtS01 = mulMatVec(BLt, S01)
                  F01 = F01 + tmp*Wt(j)*Wt(i)*BLtS01 !Summing up all Gauss Points

                  PSIt = transpose(PSI)
                  PSItPSI = rhoS*mulMat(PSIt, PSI) !Integrand for Me
                  Me = Me + tmp*Wt(j)*Wt(i)*PSItPSI !Summing up all Gauss Points

                  ! fPSI = mulMatVec(PSIt, [0.0d0, stepLoad])
                  fPSI = stepLoad(:, iEl)
                  F02 = F02 + tmp*Wt(j)*Wt(i)*fPSI
               end do
            end do

            Ke = KL + KNL
            Fe = F02 - F01

            !Assembling
            !for iEl=1:nEl !for each element
            do i = 1, nDofPerEl !for i in Ke_ij
               do j = 1, nDofPerEl !for j in Ke_ij

                  if (dofMapBC(iEl, i) == 0 .or. dofMapBC(iEl, j) == 0) then
                     cycle !skip Ke_ij for BC dof
                  end if

                  !Kg(dofMapBC(iEl, i), dofMapBC(iEl, j)) = Kg(dofMapBC(iEl, i), dofMapBC(iEl, j)) + Ke(i, j)
                  !Mg(dofMapBC(iEl, i), dofMapBC(iEl, j)) = Mg(dofMapBC(iEl, i), dofMapBC(iEl, j)) + Me(i, j)

                  ii = dofMapBC(iEl, i)
                  jj = dofMapBC(iEl, j)
                  if (i .ge. max(1, j - nUdiag) .and. i .le. min(nDofBC, j + nLdiag)) then
                     KgB(nLdiag + nUdiag + 1 + ii - jj, jj) = KgB(nLdiag + nUdiag + 1 + ii - jj, jj) + Ke(i, j)
                     MgB(nLdiag + nUdiag + 1 + ii - jj, jj) = MgB(nLdiag + nUdiag + 1 + ii - jj, jj) + Me(i, j)
                  end if

               end do

               if (dofMapBC(iEl, i) == 0) then
                  cycle !skip Qe_i for BC dof
               end if

               Fg(dofMapBC(iEl, i)) = Fg(dofMapBC(iEl, i)) + Fe(i)
            end do

         end do

      end subroutine calcMgKgFg

      subroutine writeSoln(filename)
         implicit none

         integer::i, j
         character(len=*), intent(in) :: filename

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

      end subroutine writeSoln

      subroutine mapGlobal2Local(DeltaG, DeltaEl)
         implicit none

         double precision, dimension(:), intent(in) :: DeltaG
         double precision, allocatable, dimension(:, :, :), intent(out) :: DeltaEl
         integer::i, j

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

      end subroutine mapGlobal2Local

   end subroutine compFSIresponse

   subroutine setupLBMvars()
      implicit none

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
      xc_ = xc/Clen + 1.5d0
      yc_ = yc/Clen + 1.5d0
      barL_ = barL/Clen
      barH_ = barH/Clen
      nx = int(chanL_)
      ny = chanH_
      nu_ = nu/Cnu
      uMean_ = uMean/CVel
      tau = 3*nu_ + 0.5d0
      invTau = 1.0d0/tau

      write (*, *) 'Clen = ', Clen
      write (*, *) 'Crho = ', Crho
      write (*, *) 'Ct   = ', Ct
      write (*, *) 'CVel = ', CVel
      write (*, *) 'tau  = ', tau
      write (*, *) 'uMax_ = ', 1.5d0*uMean_
      write (*, *) 'Re_SI = ', uMean*dia/nu
      write (*, *) 'Re_LBM = ', uMean_*dia_/nu_

   end subroutine setupLBMvars

   subroutine setupFEMvars()
      implicit none

      nEl = nElx*nEly
      nNodesPerEl = (degEl + 1)**2 - flag*(degEl - 1)**2 !nodes per element
      nNodes = (nElx*degEl + 1)*(nEly*degEl + 1) - flag*(nElx*nEly)*(degEl - 1)**2 !total no of global nodes
      nDofPerEl = nDofPerNode*nNodesPerEl
      nDof = nDofPerNode*nNodes !total DOF including BC DOFs
      nDofBC = nDof - size(dofBC) !!total DOF excluding BC DOFs
      lenElx = barL/nElx
      lenEly = barH/nEly
      lenElx_ = barL_/nElx
      lenEly_ = barH_/nEly
   end subroutine setupFEMvars

   subroutine initProbDist()
      implicit none

      integer:: i, j, a
      double precision:: tmp1, tmp2, uPara_

      allocate (f(0:q - 1, nx + 2, ny + 2))

      do i = 1, nx + 2
         do j = 1, ny + 2
            uPara_ = d0! 6*uMean_*(ny - (j - 1.5))*(j - 1.5)/ny**2;
            do a = 0, q - 1
               tmp1 = uPara_*ci(a, 1)
               tmp2 = uPara_*uPara_
               f(a, i, j) = wi(a)*rhoF_*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2)
            end do
         end do
      end do

   end subroutine initProbDist

   subroutine setupElemMap(dofMapBC_, coupleRange_)
      implicit none
      ! Arguments declarations
      integer, allocatable, dimension(:, :), intent(out) :: dofMapBC_
      integer, intent(out) :: coupleRange_
      integer, allocatable, dimension(:, :) :: dofMap
      integer, allocatable, dimension(:) :: cntS
      integer :: iEl, nNd
      integer :: cnt, i, j, k, p, m, n
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

   subroutine demarcateElem()!(bottomEl, topEl, rightEl, leftEl)
      implicit none

      ! integer, allocatable, dimension(:), intent(out) :: bottomEl, topEl, rightEl, leftEl
      integer:: k, cnt

      allocate (bottomEl(nElx))
      allocate (topEl(nElx))
      allocate (leftEl(nEly))
      allocate (rightEl(nEly))

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

   end subroutine demarcateElem

   subroutine undeformedInterface()!(refBounTopEl, refBounBottomEl, refBounRightEl)
      implicit none

      ! double precision, allocatable, dimension(:, :), intent(out) :: refBounTopEl, refBounBottomEl, refBounRightEl
      integer::k, m

      allocate (refBounTopEl(0:degEl*size(topEl), 2))
      allocate (refBounBottomEl(0:degEl*size(bottomEl), 2))
      allocate (refBounRightEl(0:degEl*size(rightEl), 2))

      refBounTopEl(0, 1:2) = [xc_ + 0.5d0*dia_, yc_ + 0.5d0*barH_]
      do k = 1, size(topEl)
         do m = 0, degEl - 1
            refBounTopEl(degEl*k - m, 1) = refBounTopEl(0, 1) + (k - dble(m)/degEl)*lenElx_
            refBounTopEl(degEl*k - m, 2) = refBounTopEl(0, 2)
         end do

         ! refBounTopEl(2*k-1, 1) = refBounTopEl(0, 1) + (k-0.5d0)*lenElx_
         ! refBounTopEl(2*k  , 1) = refBounTopEl(0, 1) + k*lenElx_
         ! refBounTopEl(2*k-1, 2) = refBounTopEl(0, 2)
         ! refBounTopEl(2*k  , 2) = refBounTopEl(0, 2)
      end do

      refBounBottomEl(0, 1:2) = [xc_ + 0.5d0*dia_, yc_ - 0.5d0*barH_]
      do k = 1, size(bottomEl)
         do m = 0, degEl - 1
            refBounBottomEl(degEl*k - m, 1) = refBounBottomEl(0, 1) + (k - dble(m)/degEl)*lenElx_
            refBounBottomEl(degEl*k - m, 2) = refBounBottomEl(0, 2)
         end do
         ! refBounBottomEl(2*k-1, 1) = refBounBottomEl(0, 1) + (k-0.5d0)*lenElx_
         ! refBounBottomEl(2*k  , 1) = refBounBottomEl(0, 1) + k*lenElx_
         ! refBounBottomEl(2*k-1, 2) = refBounBottomEl(0, 2)
         ! refBounBottomEl(2*k  , 2) = refBounBottomEl(0, 2)
      end do

      refBounRightEl(0, 1:2) = [xc_ + 0.5d0*dia_ + barL_, yc_ - 0.5d0*barH_]
      do k = 1, size(rightEl)
         do m = 0, degEl - 1
            refBounRightEl(degEl*k - m, 1) = refBounRightEl(0, 1)
            refBounRightEl(degEl*k - m, 2) = refBounRightEl(0, 2) + (k - dble(m)/degEl)*lenEly_
         end do
         ! refBounRightEl(2*k-1, 1) = refBounRightEl(0, 1)
         ! refBounRightEl(2*k  , 1) = refBounRightEl(0, 1)
         ! refBounRightEl(2*k-1, 2) = refBounRightEl(0, 2) + (k-0.5d0)*lenEly_
         ! refBounRightEl(2*k  , 2) = refBounRightEl(0, 2) + k*lenEly_
      end do

   end subroutine undeformedInterface

   function dateTime()

      implicit none
      character(len=30)::dateTime
      character(len=10):: ampm
      integer:: d, h, m, n, s, y, mm, values(8)
      character(len=3), parameter, dimension(12) :: &
         month = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

      call date_and_time(values=values)

      y = values(1)
      m = values(2)
      d = values(3)
      h = values(5)
      n = values(6)
      s = values(7)
      mm = values(8)

      if (h < 12) then
         ampm = 'AM'
      elseif (h == 12) then
         if (n == 0 .and. s == 0) then
            ampm = 'Noon'
         else
            ampm = 'PM'
         end if
      else
         h = h - 12
         if (h < 12) then
            ampm = 'PM'
         elseif (h == 12) then
            if (n == 0 .and. s == 0) then
               ampm = 'Midnight'
            else
               ampm = 'AM'
            end if
         end if
      end if

      write (dateTime, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)') &
         d, trim(month(m)), y, h, ':', n, ':', s, '.', mm, trim(ampm)

   end function dateTime

end module solver
