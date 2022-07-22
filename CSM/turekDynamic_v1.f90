program main
   use userInput
   use shapeFun
   use mathOp

   implicit none
   !Implements banded solver from LAPACK
   !This code uses standard numbering system
   !==========================================================================
   !Solves nonlinear 2D elasticity(static) problem using FEM
   !Uses linear rectangular element
   !DOF: u_x, u_y (displacement along x and y in a plane)
   !BC: Bottom face is fixed, Top face has constant shear stress
   !==========================================================================
   ! Variable declarations
   double precision :: a, b
   integer :: nEl, nNodesPerEl, nNodes, nDofPerEl, nDof, nDofBC, coupleRange, nUdiag, nLdiag, ldab
   integer, allocatable, dimension(:) :: dofBC
   integer, allocatable, dimension(:, :) :: dofMapBC

   double precision, allocatable, dimension(:, :) :: KgB, MgB, Mtemp
   double precision, allocatable, dimension(:) :: Fg

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

   !double precision :: Fg
   !double precision :: Freq,funT
   !integer :: invMK
   !double precision :: Kg
   !double precision :: Mg
   integer :: probeDof(2), topLeftPt

   nEl = nElx*nEly
   a = Lx/nElx
   b = Ly/nEly
   nNodesPerEl = (degEl + 1)**2 - flag*(degEl - 1)**2 !nodes per element
   nNodes = (nElx*degEl + 1)*(nEly*degEl + 1) - flag*(nElx*nEly)*(degEl - 1)**2 !total no of global nodes
   nDofPerEl = nDofPerNode*nNodesPerEl
   nDof = nDofPerNode*nNodes !total DOF including BC DOFs

   call setupBC(dofBC)
   nDofBC = nDof - size(dofBC) !total DOF excluding BC DOFs

   call setupElemMap(dofMapBC, coupleRange)!check

   !allocate (Kg(nDofBC, nDofBC))
   !allocate (Mg(nDofBC, nDofBC))
   nUdiag = coupleRange
   nLdiag = coupleRange
   ldab = 2*nLdiag + nUdiag + 1
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

   write (*, *) nElx, nEly, coupleRange
   topLeftPt = nDofBC - (nElx*degEl + 1)*nDofPerNode + 1
   probeDof = nDofBC - (0.5*nElx*degEl + 1)*nDofPerNode + [1, 2]

   call compDynRes(tStart, tEnd, dt, probeDof)

   ! Xtx(cnt)=X(probePt(2)) !/L;
   ! Xty(cnt)=-X(probePt(1)) !/L;
   !Call calcMgKgFg() for DeltaG=0;
   ! call calcMgKgFg( zeros(nDofBC,1),0 , Mg,Kg )!check
   ! invMK=Mg\Kg
   ! ![V,~]=eig(invMK);
   ! Freq=sqrt(eig(invMK))
   ! funT=2*pi/min(Freq) !Fundamental time period

contains

   subroutine compDynRes(tStart_, tEnd_, dt_, probeDof_)
      implicit none
      double precision, intent(in) :: tStart_, tEnd_, dt_
      integer, dimension(:), intent(in) :: probeDof_
      !double precision, allocatable, dimension(:,:), intent(out) :: soln

      !double precision, allocatable, dimension(:, :):: Kg, Mg, Mtemp
      !double precision, allocatable, dimension(:):: Fg
      double precision :: a0, a1, a2, a3, a4, a5, a6, a7
      double precision, allocatable, dimension(:) :: X, Xd, Xdd, XOld, XdOld, XddOld
      double precision, allocatable, dimension(:) :: Xi, Xidd, Xii, Xti

      double precision :: stepLoad, tol, err, t
      integer, allocatable, dimension(:)::pivot
      integer :: i, iLoad, cntIter, ok
      !double precision, allocatable, dimension(:) :: tVec,Xtx,Xty

      !nPerT=200 !No of samples per fundamental time period
      a0 = 1.0d0/(alpha*dt_**2.0d0)
      a1 = beta/(alpha*dt_)
      a2 = 1.0d0/(alpha*dt_)
      a3 = (1.0d0/(2.0d0*alpha)) - 1.0d0
      a4 = (beta/alpha) - 1.0d0
      a5 = (dt_/2.0d0)*((beta/alpha) - 2.0d0)
      a6 = dt_*(1.0d0 - beta)
      a7 = beta*dt_

      allocate (Mtemp(nDofBC, nDofBC))

      allocate (pivot(nDofBC))
      allocate (X(nDofBC), Xd(nDofBC), Xdd(nDofBC))
      allocate (XOld(nDofBC), XdOld(nDofBC), XddOld(nDofBC))
      allocate (Xi(nDofBC), Xii(nDofBC), Xti(nDofBC), Xidd(nDof))

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
      t = 0.0d0

      open (UNIT=10, file='dynamic.dat')

      do while (t .le. tEnd_)

         do iLoad = 1, noOfLoadSteps
            !loadVec=linspace(totalLoad/noOfLoadSteps,totalLoad,noOfLoadSteps)
            !stepLoad=loadVec(iLoad)
            stepLoad = totalLoad*iLoad/noOfLoadSteps

            !Xi=zeros(nDofBC,1);
            !call calcMgKgFg(X, stepLoad, MgB, KgB, Fg)
            call calcMgKgFg(X, stepLoad)
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
               call calcMgKgFg(Xti, stepLoad)
               !Xii=(a0*Mg + Kg)\(Fg - Mg*Xidd)
               Xii = Fg - mulMatBVec(MgB, Xidd, nUdiag, nLdiag)
               Mtemp = a0*MgB + KgB
               !call DGESV(nDofBC, 1, Mtemp, nDofBC, pivot, Xii, nDofBC, ok)
               call DGBSV(nDofBC, nLdiag, nUdiag, 1, Mtemp, ldab, pivot, Xii, nDofBC, ok)
               Xi = Xi + Xii
               err = norm2(Xii)/norm2(Xi + X)
            end do !iter loop ends

            !if (modulo(t/dt,20.0d0)==0.0d0) then
            write (*, *) t, err, cntIter
            !disp([num2str(t),' ',num2str(err),' ',num2str(cntIter)])
            !end if

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
         t = t + dt_
         write (10, '(3(E12.4,2X))') t, X(probeDof_(1)), X(probeDof_(2))

      end do
      close (10)

   end subroutine compDynRes

   !subroutine calcMgKgFg(DeltaG, stepLoad, MgB, KgB, Fg)
   subroutine calcMgKgFg(DeltaG, stepLoad)
      implicit none
      double precision, dimension(:), intent(in) :: DeltaG
      double precision, intent(in) :: stepLoad
      !double precision, dimension(:, :), intent(out) :: KgB, MgB
      !double precision, dimension(:), intent(out) :: Fg

      double precision :: tmp
      integer :: i, j, k, iEl, ii, jj

      Wt = [5.0d0/9.0d0, 8.0d0/9.0d0, 5.0d0/9.0d0] !Gauss Weights
      Pt = [-sqrt(3.0d0/5.0d0), 0.0d0, sqrt(3.0d0/5.0d0)] !Gauss Points
      tmp = 0.25d0*a*b

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
               dpsiNdxGP = dpsiNdXi(Pt(i), Pt(j))*2/a !=dpsiNdXi*dXidx
               dpsiNdyGP = dpsiNdEta(Pt(i), Pt(j))*2/b !=dpsiNdEta*dEtady

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
               BLtCBL = Lz*mulMat(BLt, mulMat(C, BL)) !Integrand for KL
               KL = KL + tmp*Wt(j)*Wt(i)*BLtCBL !Summing up all Gauss Points

               BNLt = transpose(BNL)
               BNLtCBNL = Lz*mulMat(BNLt, mulMat(S_Mat, BNL)) !Integrand for KNL
               KNL = KNL + tmp*Wt(j)*Wt(i)*BNLtCBNL !Summing up all Gauss Points

               BLtS01 = Lz*mulMatVec(BLt, S01)
               F01 = F01 + tmp*Wt(j)*Wt(i)*BLtS01 !Summing up all Gauss Points

               PSIt = transpose(PSI)
               PSItPSI = rho*Lz*mulMat(PSIt, PSI) !Integrand for Me
               Me = Me + tmp*Wt(j)*Wt(i)*PSItPSI !Summing up all Gauss Points

               fPSI = Lz*mulMatVec(PSIt, [stepLoad, 0.0d0])
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
         !end
      end do
   end subroutine calcMgKgFg

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

      do n = 0, nEly - 1
         if (n > 0 .and. n < nEly) then
            cnt = cnt - (nElx*degEl + 1)*nDofPerNode !to account for common nodes along y
         end if
         do j = 0, degEl
            do m = 0, nElx - 1
               if (m > 0 .and. m < nElx) then
                  cnt = cnt - nDofPerNode !to account for common nodes along x
               end if
               do i = 0, degEl
                  iEl = (m + 1) + nElx*n !global element number
                  if (flag == 1 .and. i > 0 .and. i < degEl .and. j > 0 .and. j < degEl) then
                     cntS(iEl) = cntS(iEl) + 1 !count interior nodes for a given element
                     cycle !skip interior node for serendipity element
                  end if
                  nNd = (i + 1) + (degEl + 1)*j - cntS(iEl) !local node number for a given element
                  !dofMap(iEl,[2*nNd-1 2*nNd])=[ cnt,cnt+1 ] !allot global nodal dof
                  dofMap(iEl, 2*nNd - 1) = cnt
                  dofMap(iEl, 2*nNd) = cnt + 1
                  cnt = cnt + nDofPerNode
               end do
            end do
         end do
      end do

      coupleRange_ = dofMap(1, nDofPerEl) - dofMap(1, 1)

      !Replacing BC nodes wi !Replacing BC nodes with zero and renumbering the DOF map(works even for unordered dofBC)
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

end program main
