module lbm
   implicit none
   public

   double precision, allocatable, dimension(:, :), protected:: ci
   double precision, allocatable, dimension(:), protected:: wi
   integer, allocatable, dimension(:), protected::kb

contains

   subroutine oppoVector()!(ci, kb)
      implicit none

      ! double precision, dimension(0:, :), intent(in) :: ci
      ! integer, allocatable, dimension(:), intent(out) :: kb
      integer:: a, a1

      allocate (kb(0:size(ci, 1) - 1))

      do a = 0, size(ci, 1) - 1
         do a1 = a, size(ci, 1) - 1
            if (ci(a, 1) + ci(a1, 1) .eq. 0.0d0 .and. ci(a, 2) + ci(a1, 2) .eq. 0.0d0) then
               kb(a) = a1
               kb(a1) = a
            end if
         end do
      end do

   end subroutine oppoVector

   subroutine setupD2Q9()!(ci, wi, minusVector)
      implicit none
      integer, parameter::d = 2, q = 9
      integer:: a
      ! double precision, allocatable, dimension(:, :), intent(out) :: ci
      ! double precision, allocatable, dimension(:), intent(out) :: wi
      ! integer, allocatable, dimension(:), intent(out):: minusVector
      allocate (ci(0:q - 1, d))
      allocate (wi(0:q - 1))

      ci(0, :) = [0.0d0, 0.0d0]
      ci(1, :) = [1.0d0, 0.0d0]
      ci(2, :) = [0.0d0, 1.0d0]
      ci(3, :) = [-1.0d0, 0.0d0]
      ci(4, :) = [0.0d0, -1.0d0]
      ci(5, :) = [1.0d0, 1.0d0]
      ci(6, :) = [-1.0d0, 1.0d0]
      ci(7, :) = [-1.0d0, -1.0d0]
      ci(8, :) = [1.0d0, -1.0d0]

      do a = 0, q - 1
         if (a .eq. 0) wi(a) = 4.0d0/9.0d0
         if (a .ge. 1 .and. a .le. 4) wi(a) = 1.0d0/9.0d0
         if (a .ge. 5 .and. a .le. 8) wi(a) = 1.0d0/36.0d0
      end do

      call oppoVector()!(ci, minusVector)

   end subroutine setupD2Q9

   subroutine setupD3Q15()!(ci, wi)
      implicit none
      integer, parameter::d = 3, q = 15
      integer:: a
      ! integer, allocatable, dimension(:, :), intent(out) :: ci
      ! double precision, allocatable, dimension(:), intent(out) :: wi
      allocate (ci(0:q - 1, d))
      allocate (wi(0:q - 1))

      ci(0, :) = [0, 0, 0]
      ci(1, :) = [1, 0, 0]
      ci(2, :) = [-1, 0, 0]
      ci(3, :) = [0, 1, 0]
      ci(4, :) = [0, -1, 0]
      ci(5, :) = [0, 0, 1]
      ci(6, :) = [0, 0, -1]
      ci(7, :) = [1, 1, 1]
      ci(8, :) = [-1, -1, -1]
      ci(9, :) = [1, 1, -1]
      ci(10, :) = [-1, -1, 1]
      ci(11, :) = [1, -1, 1]
      ci(12, :) = [-1, 1, -1]
      ci(13, :) = [-1, 1, 1]
      ci(14, :) = [1, -1, -1]

      do a = 0, q - 1
         if (a .eq. 0) wi(a) = 2.0d0/9.0d0
         if (a .ge. 1 .and. a .le. 6) wi(a) = 1.0d0/9.0d0
         if (a .ge. 7 .and. a .le. 14) wi(a) = 1.0d0/72.0d0
      end do

   end subroutine setupD3Q15

   subroutine setupD3Q19()!(ci, wi)
      implicit none
      integer, parameter::d = 3, q = 19
      integer:: a
      ! integer, allocatable, dimension(:, :), intent(out) :: ci
      ! double precision, allocatable, dimension(:), intent(out) :: wi
      allocate (ci(0:q - 1, d))
      allocate (wi(0:q - 1))

      ci(0, :) = [0, 0, 0]
      ci(1, :) = [1, 0, 0]
      ci(2, :) = [-1, 0, 0]
      ci(3, :) = [0, 1, 0]
      ci(4, :) = [0, -1, 0]
      ci(5, :) = [0, 0, 1]
      ci(6, :) = [0, 0, -1]
      ci(7, :) = [1, 1, 0]
      ci(8, :) = [-1, -1, 0]
      ci(9, :) = [1, 0, 1]
      ci(10, :) = [-1, 0, -1]
      ci(11, :) = [0, 1, 1]
      ci(12, :) = [0, -1, -1]
      ci(13, :) = [1, -1, 0]
      ci(14, :) = [-1, 1, 0]
      ci(15, :) = [1, 0, -1]
      ci(16, :) = [-1, 0, 1]
      ci(17, :) = [0, 1, -1]
      ci(18, :) = [0, -1, 1]

      do a = 0, q - 1
         if (a .eq. 0) wi(a) = 1.0d0/3.0d0
         if (a .ge. 1 .and. a .le. 6) wi(a) = 1.0d0/18.0d0
         if (a .ge. 7 .and. a .le. 14) wi(a) = 1.0d0/36.0d0
      end do

   end subroutine setupD3Q19

end module lbm
