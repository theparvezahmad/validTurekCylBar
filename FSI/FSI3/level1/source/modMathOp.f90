module mathOp
   implicit none

   public

contains

   function dotProd(A, B)
      implicit none

      double precision, dimension(:), intent(in)  :: A
      double precision, dimension(:), intent(in)  :: B
      double precision :: dotProd
      integer :: i

      if (size(A) /= size(B)) then
         write (*, *) 'Incompatible Vectors'
      else
         dotProd = 0.0
         do i = 1, size(A)
            dotProd = dotProd + A(i)*B(i)
         end do
      end if
   end function dotProd

   function mulMat(A, B)
      implicit none

      double precision, dimension(:, :), intent(in)  :: A
      double precision, dimension(:, :), intent(in)  :: B
      double precision, dimension(size(A, 1), size(B, 2)) :: mulMat
      integer :: i, j, k

      if (size(A, 2) /= size(B, 1)) then
         write (*, *) 'Incompatible Martices'
      else
         do i = 1, size(A, 1)
            do j = 1, size(B, 2)
               mulMat(i, j) = 0.0
               do k = 1, size(A, 2)
                  mulMat(i, j) = mulMat(i, j) + A(i, k)*B(k, j)
               end do
            end do
         end do
      end if

   end function mulMat

   function mulMatVec(A, B)
      implicit none

      double precision, dimension(:, :), intent(in)  :: A
      double precision, dimension(:), intent(in)  :: B
      double precision, dimension(size(A, 1)) :: mulMatVec
      integer :: i, k

      if (size(A, 2) /= size(B)) then
         write (*, *) 'Incompatible Martix-Vector Multiplication'
      else
         do i = 1, size(A, 1)
            mulMatVec(i) = 0.0d0
            do k = 1, size(A, 2)
               mulMatVec(i) = mulMatVec(i) + A(i, k)*B(k)
            end do
         end do
      end if

   end function mulMatVec

   function mulMatBVec(AB, b, nUdiag, nLdiag)
      implicit none

      double precision, dimension(:, :), intent(in)  :: AB
      double precision, dimension(:), intent(in)  :: b
      integer, intent(in):: nUdiag, nLdiag
      double precision, dimension(size(b)) :: mulMatBVec
      integer :: i, j, n

      if (size(AB, 2) /= size(b)) then
         write (*, *) 'Incompatible Band Martix-Vector Multiplication'
      else
         n = size(AB, 2)
         do i = 1, n
            mulMatBVec(i) = 0.0d0
            do j = 1, n
               if (i .ge. max(1, j - nUdiag) .and. i .le. min(n, j + nLdiag)) then
                  mulMatBVec(i) = mulMatBVec(i) + AB(nLdiag + nUdiag + 1 + i - j, j)*b(j)
               end if
            end do
         end do
      end if

   end function mulMatBVec

   function transMat(M)
      implicit none
      double precision, dimension(:, :), intent(in)  :: M
      double precision, dimension(size(M, 2), size(M, 1)) :: transMat
      integer :: i, j

      do i = 1, size(M, 1)
         do j = 1, size(M, 2)
            transMat(j, i) = M(i, j)
         end do
      end do

   end function transMat

   function linPieceWiseYonX(x, listOfCoord) result(y)
      implicit none

      double precision, intent(in)::x
      double precision, dimension(:, :), intent(in):: listOfCoord
      double precision::y, x1, x2, y1, y2
      integer:: i

      if (x .lt. minval(listOfCoord(:, 1)) .or. x .gt. maxval(listOfCoord(:, 1))) then
         write (*, *) "Function linPieceWiseYonX Aborted"
         write (*, *) "Error: Requested evaluation point out of data range"
         stop
      end if

      do i = 1, size(listOfCoord, 1) - 1

         x1 = listOfCoord(i, 1)
         x2 = listOfCoord(i + 1, 1)
         y1 = listOfCoord(i, 2)
         y2 = listOfCoord(i + 1, 2)

         if (x .ge. x1 .and. x .le. x2) then
            y = dotProd([1 - (x - x1)/(x2 - x1), (x - x1)/(x2 - x1)], [y1, y2])
            exit
         end if

      end do
   end function linPieceWiseYonX

   function linPieceWiseXonY(y, listOfCoord) result(x)
      implicit none

      double precision, intent(in)::y
      double precision, dimension(:, :), intent(in):: listOfCoord
      double precision::x, x1, x2, y1, y2
      integer:: i

      if (y .lt. minval(listOfCoord(:, 2)) .or. y .gt. maxval(listOfCoord(:, 2))) then
         write (*, *) "Function linPieceWiseXonY Aborted"
         write (*, *) "Error: Requested evaluation point out of data range"
         stop
      end if

      do i = 1, size(listOfCoord, 1) - 1

         x1 = listOfCoord(i, 1)
         x2 = listOfCoord(i + 1, 1)
         y1 = listOfCoord(i, 2)
         y2 = listOfCoord(i + 1, 2)

         if (y .ge. y1 .and. y .le. y2) then
            x = dotProd([1 - (y - y1)/(y2 - y1), (y - y1)/(y2 - y1)], [x1, x2])
            exit
         end if

      end do

   end function linPieceWiseXonY

   function linPieceWiseIntersect(pt, listOfCoord) result(xIntersect)
      !Accepts list of (x,y), joins them using piece-wise lines
      !Outputs intersection of a horizontal line passing through pt
      implicit none

      double precision, dimension(:), intent(in)::pt
      double precision, dimension(:, :), intent(in):: listOfCoord
      double precision::x, y, x1, x2, y1, y2
      double precision, dimension(20):: xIntersect_!, yIntersect_
      double precision, allocatable, dimension(:):: xIntersect!, yIntersect
      integer:: i, cntX!, cntY

      cntX = 0
      ! cntY=0

      x = pt(1)
      y = pt(2)
      do i = 1, size(listOfCoord, 1) - 1

         ! if(i==size(listOfCoord,1)) then
         !    ip1=1
         ! else
         !    ip1=i+1
         ! end if

         x1 = listOfCoord(i, 1)
         x2 = listOfCoord(i + 1, 1)
         y1 = listOfCoord(i, 2)
         y2 = listOfCoord(i + 1, 2)

         ! if (x .ge. x1 .and. x .le. x2 .and. x1 /= x2) then
         !    cntY = cntY + 1
         !    yIntersect_(cntY) = dotProd([1 - (x - x1)/(x2 - x1), (x - x1)/(x2 - x1)], [y1, y2])
         ! end if

         if (y .ge. min(y1, y2) .and. y .le. max(y1, y2) .and. y1 /= y2) then
            cntX = cntX + 1
            xIntersect_(cntX) = dotProd([1 - (y - y1)/(y2 - y1), (y - y1)/(y2 - y1)], [x1, x2])
         end if

      end do

      allocate (xIntersect(cntX))
      ! allocate(yIntersect(cntY))

      xIntersect = xIntersect_(1:cntX)
      ! yIntersect=yIntersect_(1:cntY)

   end function linPieceWiseIntersect

   function sortAscendUnique(A) result(uniqSortedA)
      !Accepts an array A
      !Outputs sorted array in ascending order without duplication
      implicit none

      double precision, dimension(:), intent(in)::A
      double precision, allocatable, dimension(:)::uniqSortedA
      double precision, dimension(size(A))::sortedA, uniqSortedA_
      double precision::tmp
      integer::i, cnt
      logical::isUnsorted

      cnt = 0
      isUnsorted = .true.
      sortedA = A
      do while (isUnsorted)
         isUnsorted = .false.
         do i = 1, size(sortedA) - 1
            if (sortedA(i) .gt. sortedA(i + 1)) then
               isUnsorted = .true.
               cnt = cnt + 1
               tmp = sortedA(i)
               sortedA(i) = sortedA(i + 1)
               sortedA(i + 1) = tmp
            end if
         end do
      end do

      cnt = 1
      uniqSortedA_(1) = sortedA(1)
      do i = 1, size(sortedA) - 1
         if (sortedA(i) /= sortedA(i + 1)) then
            cnt = cnt + 1
            uniqSortedA_(cnt) = sortedA(i + 1)
         end if
      end do

      allocate (uniqSortedA(cnt))
      uniqSortedA = uniqSortedA_(1:cnt)

   end function sortAscendUnique

   function quadPieceWiseIntersect(pt, listOfCoord) result(xIntersect)
      !Accepts list of (x,y), joins them using piece-wise lines
      !Outputs intersection of a horizontal line passing through pt
      implicit none

      double precision, dimension(:), intent(in)::pt
      double precision, dimension(:, :), intent(in):: listOfCoord
      double precision::x, y, x1, x2, x3, y1, y2, y3, l2
      double precision, dimension(20):: xIntersect_!, yIntersect_
      double precision, allocatable, dimension(:):: xIntersect!, yIntersect
      integer:: i, cntX!, cntY

      cntX = 0
      ! cntY=0

      x = pt(1)
      y = pt(2)
      do i = 1, size(listOfCoord, 1) - 2, 2

         x1 = listOfCoord(i, 1)
         x2 = listOfCoord(i + 1, 1)
         x3 = listOfCoord(i + 2, 1)
         y1 = listOfCoord(i, 2)
         y2 = listOfCoord(i + 1, 2)
         y3 = listOfCoord(i + 2, 2)
         l2 = (y3 - y1)**2.0d0

         ! if (x .ge. x1 .and. x .le. x2 .and. x1 /= x2) then
         !    cntY = cntY + 1
         !    yIntersect_(cntY) = dotProd([1 - (x - x1)/(x2 - x1), (x - x1)/(x2 - x1)], [y1, y2])
         ! end if

         if (y .ge. min(y1, y3) .and. y .le. max(y1, y3) .and. y1 /= y3) then
            cntX = cntX + 1
    xIntersect_(cntX) = dotProd([2.0d0*(y2 - y)*(y3 - y)/l2, -4.0d0*(y1 - y)*(y3 - y)/l2, 2.0d0*(y1 - y)*(y2 - y)/l2], [x1, x2, x3])
         end if

      end do

      allocate (xIntersect(cntX))
      ! allocate(yIntersect(cntY))

      xIntersect = xIntersect_(1:cntX)
      ! yIntersect=yIntersect_(1:cntY)

   end function quadPieceWiseIntersect

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

   subroutine addDuplicateFields(surfForce, uniqSurfForce)
      implicit none
      double precision, dimension(:, :), intent(in) :: surfForce
      double precision, allocatable, dimension(:, :), intent(out) :: uniqSurfForce
      
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

end module mathOp
