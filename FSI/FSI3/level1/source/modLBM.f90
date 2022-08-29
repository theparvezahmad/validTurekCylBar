module lbm
   implicit none

contains

   subroutine oppoVector(vector, kb)
      implicit none

      integer, dimension(:, 0:), intent(in) :: vector
      integer, allocatable, dimension(:), intent(out) :: kb
      integer:: a, a1

      allocate (kb(0:size(vector, 2) - 1))

      do a = 0, size(vector, 2) - 1
         do a1 = a, size(vector, 2) - 1
            if (vector(1, a) + vector(1, a1) .eq. 0 .and. vector(2, a) + vector(2, a1) .eq. 0) then
               kb(a) = a1
               kb(a1) = a
            end if
         end do
      end do

   end subroutine oppoVector

   subroutine setupD2Q9(vector, weight, minusVector)
      implicit none
      integer, parameter::d = 2, q = 9
      integer:: a
      integer, allocatable, dimension(:, :), intent(out) :: vector
      double precision, allocatable, dimension(:), intent(out) :: weight
      integer, allocatable, dimension(:), intent(out):: minusVector
      allocate (vector(d, 0:q - 1))
      allocate (weight(0:q - 1))

      vector(:, 0) = [0, 0]
      vector(:, 1) = [1, 0]
      vector(:, 2) = [0, 1]
      vector(:, 3) = [-1, 0]
      vector(:, 4) = [0, -1]
      vector(:, 5) = [1, 1]
      vector(:, 6) = [-1, 1]
      vector(:, 7) = [-1, -1]
      vector(:, 8) = [1, -1]

      do a = 0, q - 1
         if (a .eq. 0) weight(a) = 4.0d0/9.0d0
         if (a .ge. 1 .and. a .le. 4) weight(a) = 1.0d0/9.0d0
         if (a .ge. 5 .and. a .le. 8) weight(a) = 1.0d0/36.0d0
      end do

      call oppoVector(vector, minusVector)

   end subroutine setupD2Q9

   subroutine setupD3Q15(vector, weight)
      implicit none
      integer, parameter::d = 3, q = 15
      integer:: a
      integer, allocatable, dimension(:, :), intent(out) :: vector
      double precision, allocatable, dimension(:), intent(out) :: weight
      allocate (vector(d, 0:q - 1))
      allocate (weight(0:q - 1))

      vector(:, 0) = [0, 0, 0]
      vector(:, 1) = [1, 0, 0]
      vector(:, 2) = [-1, 0, 0]
      vector(:, 3) = [0, 1, 0]
      vector(:, 4) = [0, -1, 0]
      vector(:, 5) = [0, 0, 1]
      vector(:, 6) = [0, 0, -1]
      vector(:, 7) = [1, 1, 1]
      vector(:, 8) = [-1, -1, -1]
      vector(:, 9) = [1, 1, -1]
      vector(:, 10) = [-1, -1, 1]
      vector(:, 11) = [1, -1, 1]
      vector(:, 12) = [-1, 1, -1]
      vector(:, 13) = [-1, 1, 1]
      vector(:, 14) = [1, -1, -1]

      do a = 0, q - 1
         if (a .eq. 0) weight(a) = 2.0d0/9.0d0
         if (a .ge. 1 .and. a .le. 6) weight(a) = 1.0d0/9.0d0
         if (a .ge. 7 .and. a .le. 14) weight(a) = 1.0d0/72.0d0
      end do

   end subroutine setupD3Q15

   subroutine setupD3Q19(vector, weight)
      implicit none
      integer, parameter::d = 3, q = 19
      integer:: a
      integer, allocatable, dimension(:, :), intent(out) :: vector
      double precision, allocatable, dimension(:), intent(out) :: weight
      allocate (vector(d, 0:q - 1))
      allocate (weight(0:q - 1))

      vector(:, 0) = [0, 0, 0]
      vector(:, 1) = [1, 0, 0]
      vector(:, 2) = [-1, 0, 0]
      vector(:, 3) = [0, 1, 0]
      vector(:, 4) = [0, -1, 0]
      vector(:, 5) = [0, 0, 1]
      vector(:, 6) = [0, 0, -1]
      vector(:, 7) = [1, 1, 0]
      vector(:, 8) = [-1, -1, 0]
      vector(:, 9) = [1, 0, 1]
      vector(:, 10) = [-1, 0, -1]
      vector(:, 11) = [0, 1, 1]
      vector(:, 12) = [0, -1, -1]
      vector(:, 13) = [1, -1, 0]
      vector(:, 14) = [-1, 1, 0]
      vector(:, 15) = [1, 0, -1]
      vector(:, 16) = [-1, 0, 1]
      vector(:, 17) = [0, 1, -1]
      vector(:, 18) = [0, -1, 1]

      do a = 0, q - 1
         if (a .eq. 0) weight(a) = 1.0d0/3.0d0
         if (a .ge. 1 .and. a .le. 6) weight(a) = 1.0d0/18.0d0
         if (a .ge. 7 .and. a .le. 14) weight(a) = 1.0d0/36.0d0
      end do

   end subroutine setupD3Q19

end module lbm
