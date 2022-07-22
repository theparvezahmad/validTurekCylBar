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

end module mathOp
