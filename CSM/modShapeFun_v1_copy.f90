module shapeFun
   use userInput, only: degEl, flag
   implicit none

   integer, parameter, private::n = (degEl + 1)**2 - flag*(degEl - 1)**2 !nodes per element

contains

   function psiN(eta, xi)
      implicit none

      double precision, allocatable, dimension(:)::psiN
      double precision, intent(in) :: eta, xi
      allocate (psiN(n))

      select case (n)
      case (4)
         psiN = [((eta - 1)*(xi - 1))/4, &
                 -((eta - 1)*(xi + 1))/4, &
                 -((eta + 1)*(xi - 1))/4, &
                 ((eta + 1)*(xi + 1))/4]
      case (8)
         psiN = [-((eta - 1)*(xi - 1)*(eta + xi + 1))/4, &
                 ((eta - 1)*(xi - 1)*(xi + 1))/2, &
                 ((eta - 1)*(xi + 1)*(eta - xi + 1))/4, &
                 ((eta - 1)*(eta + 1)*(xi - 1))/2, &
                 -((eta - 1)*(eta + 1)*(xi + 1))/2, &
                 ((eta + 1)*(xi - 1)*(xi - eta + 1))/4, &
                 -((eta + 1)*(xi - 1)*(xi + 1))/2, &
                 ((eta + 1)*(xi + 1)*(eta + xi - 1))/4]
      case (9)
         psiN = [(eta*xi*(eta - 1)*(xi - 1))/4, &
                 -(eta*(eta - 1)*(xi - 1)*(xi + 1))/2, &
                 (eta*xi*(eta - 1)*(xi + 1))/4, &
                 -(xi*(eta - 1)*(eta + 1)*(xi - 1))/2, &
                 (eta - 1)*(eta + 1)*(xi - 1)*(xi + 1), &
                 -(xi*(eta - 1)*(eta + 1)*(xi + 1))/2, &
                 (eta*xi*(eta + 1)*(xi - 1))/4, &
                 -(eta*(eta + 1)*(xi - 1)*(xi + 1))/2, &
                 (eta*xi*(eta + 1)*(xi + 1))/4]
      case (12)
         psiN = [((eta - 1)*(xi - 1)*(9*eta**2 + 9*xi**2 - 10))/32, &
                 -(9*(3*xi - 1)*(eta - 1)*(xi - 1)*(xi + 1))/32, &
                 (9*(3*xi + 1)*(eta - 1)*(xi - 1)*(xi + 1))/32, &
                 -((eta - 1)*(xi + 1)*(9*eta**2 + 9*xi**2 - 10))/32, &
                 -(9*(3*eta - 1)*(eta - 1)*(eta + 1)*(xi - 1))/32, &
                 (9*(3*eta - 1)*(eta - 1)*(eta + 1)*(xi + 1))/32, &
                 (9*(3*eta + 1)*(eta - 1)*(eta + 1)*(xi - 1))/32, &
                 -(9*(3*eta + 1)*(eta - 1)*(eta + 1)*(xi + 1))/32, &
                 -((eta + 1)*(xi - 1)*(9*eta**2 + 9*xi**2 - 10))/32, &
                 (9*(3*xi - 1)*(eta + 1)*(xi - 1)*(xi + 1))/32, &
                 -(9*(3*xi + 1)*(eta + 1)*(xi - 1)*(xi + 1))/32, &
                 ((eta + 1)*(xi + 1)*(9*eta**2 + 9*xi**2 - 10))/32]
      case (16)
         psiN = [((3*eta - 1)*(3*eta + 1)*(3*xi - 1)*(3*xi + 1)*(eta - 1)*(xi - 1))/256, &
                 -(9*(3*eta - 1)*(3*eta + 1)*(3*xi - 1)*(eta - 1)*(xi - 1)*(xi + 1))/256, &
                 (9*(3*eta - 1)*(3*eta + 1)*(3*xi + 1)*(eta - 1)*(xi - 1)*(xi + 1))/256, &
                 -((3*eta - 1)*(3*eta + 1)*(3*xi - 1)*(3*xi + 1)*(eta - 1)*(xi + 1))/256, &
                 -(9*(3*eta - 1)*(3*xi - 1)*(3*xi + 1)*(eta - 1)*(eta + 1)*(xi - 1))/256, &
                 (81*(3*eta - 1)*(3*xi - 1)*(eta - 1)*(eta + 1)*(xi - 1)*(xi + 1))/256, &
                 -(81*(3*eta - 1)*(3*xi + 1)*(eta - 1)*(eta + 1)*(xi - 1)*(xi + 1))/256, &
                 (9*(3*eta - 1)*(3*xi - 1)*(3*xi + 1)*(eta - 1)*(eta + 1)*(xi + 1))/256, &
                 (9*(3*eta + 1)*(3*xi - 1)*(3*xi + 1)*(eta - 1)*(eta + 1)*(xi - 1))/256, &
                 -(81*(3*eta + 1)*(3*xi - 1)*(eta - 1)*(eta + 1)*(xi - 1)*(xi + 1))/256, &
                 (81*(3*eta + 1)*(3*xi + 1)*(eta - 1)*(eta + 1)*(xi - 1)*(xi + 1))/256, &
                 -(9*(3*eta + 1)*(3*xi - 1)*(3*xi + 1)*(eta - 1)*(eta + 1)*(xi + 1))/256, &
                 -((3*eta - 1)*(3*eta + 1)*(3*xi - 1)*(3*xi + 1)*(eta + 1)*(xi - 1))/256, &
                 (9*(3*eta - 1)*(3*eta + 1)*(3*xi - 1)*(eta + 1)*(xi - 1)*(xi + 1))/256, &
                 -(9*(3*eta - 1)*(3*eta + 1)*(3*xi + 1)*(eta + 1)*(xi - 1)*(xi + 1))/256, &
                 ((3*eta - 1)*(3*eta + 1)*(3*xi - 1)*(3*xi + 1)*(eta + 1)*(xi + 1))/256]
      case default
         write (*, *) 'Selected Element Type is Unavailable'
         stop

      end select

   end function psiN

   function dpsiNdXi(eta, xi)
      implicit none

      double precision, allocatable, dimension(:)::dpsiNdXi
      double precision, intent(in) :: eta, xi
      allocate (dpsiNdXi(n))

      select case (n)
      case (4)
         dpsiNdXi = [eta/4 - 1/4, &
                     1/4 - eta/4, &
                     -eta/4 - 1/4, &
                     eta/4 + 1/4]
      case (8)
         dpsiNdXi = [-((eta + 2*xi)*(eta - 1))/4, &
                     xi*(eta - 1), &
                     ((eta - 2*xi)*(eta - 1))/4, &
                     ((eta - 1)*(eta + 1))/2, &
                     -((eta - 1)*(eta + 1))/2, &
                     -((eta - 2*xi)*(eta + 1))/4, &
                     -xi*(eta + 1), &
                     ((eta + 2*xi)*(eta + 1))/4]
      case (9)
         dpsiNdXi = [(eta*(2*xi - 1)*(eta - 1))/4, &
                     -eta*xi*(eta - 1), &
                     (eta*(2*xi + 1)*(eta - 1))/4, &
                     -((2*xi - 1)*(eta - 1)*(eta + 1))/2, &
                     2*xi*(eta - 1)*(eta + 1), &
                     -((2*xi + 1)*(eta - 1)*(eta + 1))/2, &
                     (eta*(2*xi - 1)*(eta + 1))/4, &
                     -eta*xi*(eta + 1), &
                     (eta*(2*xi + 1)*(eta + 1))/4]
      case (12)
         dpsiNdXi = [-((eta - 1)*(-9*eta**2 - 27*xi**2 + 18*xi + 10))/32, &
                     (9*(eta - 1)*(-9*xi**2 + 2*xi + 3))/32, &
                     (9*(eta - 1)*(9*xi**2 + 2*xi - 3))/32, &
                     -((eta - 1)*(9*eta**2 + 27*xi**2 + 18*xi - 10))/32, &
                     -(9*(3*eta - 1)*(eta - 1)*(eta + 1))/32, &
                     (9*(3*eta - 1)*(eta - 1)*(eta + 1))/32, &
                     (9*(3*eta + 1)*(eta - 1)*(eta + 1))/32, &
                     -(9*(3*eta + 1)*(eta - 1)*(eta + 1))/32, &
                     ((eta + 1)*(-9*eta**2 - 27*xi**2 + 18*xi + 10))/32, &
                     -(9*(eta + 1)*(-9*xi**2 + 2*xi + 3))/32, &
                     -(9*(eta + 1)*(9*xi**2 + 2*xi - 3))/32, &
                     ((eta + 1)*(9*eta**2 + 27*xi**2 + 18*xi - 10))/32]
      case (16)
         dpsiNdXi = [-((3*eta - 1)*(3*eta + 1)*(eta - 1)*(-27*xi**2 + 18*xi + 1))/256, &
                     (9*(3*eta - 1)*(3*eta + 1)*(eta - 1)*(-9*xi**2 + 2*xi + 3))/256, &
                     (9*(3*eta - 1)*(3*eta + 1)*(eta - 1)*(9*xi**2 + 2*xi - 3))/256, &
                     -((3*eta - 1)*(3*eta + 1)*(eta - 1)*(27*xi**2 + 18*xi - 1))/256, &
                     (9*(3*eta - 1)*(eta - 1)*(eta + 1)*(-27*xi**2 + 18*xi + 1))/256, &
                     -(81*(3*eta - 1)*(eta - 1)*(eta + 1)*(-9*xi**2 + 2*xi + 3))/256, &
                     -(81*(3*eta - 1)*(eta - 1)*(eta + 1)*(9*xi**2 + 2*xi - 3))/256, &
                     (9*(3*eta - 1)*(eta - 1)*(eta + 1)*(27*xi**2 + 18*xi - 1))/256, &
                     -(9*(3*eta + 1)*(eta - 1)*(eta + 1)*(-27*xi**2 + 18*xi + 1))/256, &
                     (81*(3*eta + 1)*(eta - 1)*(eta + 1)*(-9*xi**2 + 2*xi + 3))/256, &
                     (81*(3*eta + 1)*(eta - 1)*(eta + 1)*(9*xi**2 + 2*xi - 3))/256, &
                     -(9*(3*eta + 1)*(eta - 1)*(eta + 1)*(27*xi**2 + 18*xi - 1))/256, &
                     ((3*eta - 1)*(3*eta + 1)*(eta + 1)*(-27*xi**2 + 18*xi + 1))/256, &
                     -(9*(3*eta - 1)*(3*eta + 1)*(eta + 1)*(-9*xi**2 + 2*xi + 3))/256, &
                     -(9*(3*eta - 1)*(3*eta + 1)*(eta + 1)*(9*xi**2 + 2*xi - 3))/256, &
                     ((3*eta - 1)*(3*eta + 1)*(eta + 1)*(27*xi**2 + 18*xi - 1))/256]
      case default
         write (*, *) 'Selected Element Type is Unavailable'
         stop

      end select

   end function dpsiNdXi

   function dpsiNdEta(eta, xi)
      implicit none

      double precision, allocatable, dimension(:)::dpsiNdEta
      double precision, intent(in) :: eta, xi
      allocate (dpsiNdEta(n))

      select case (n)
      case (4)
         dpsiNdEta = [xi/4 - 1/4, &
                      -xi/4 - 1/4, &
                      1/4 - xi/4, &
                      xi/4 + 1/4]
      case (8)
         dpsiNdEta = [-((2*eta + xi)*(xi - 1))/4, &
                      ((xi - 1)*(xi + 1))/2, &
                      ((xi + 1)*(2*eta - xi))/4, &
                      eta*(xi - 1), &
                      -eta*(xi + 1), &
                      -((xi - 1)*(2*eta - xi))/4, &
                      -((xi - 1)*(xi + 1))/2, &
                      ((2*eta + xi)*(xi + 1))/4]
      case (9)
         dpsiNdEta = [(xi*(2*eta - 1)*(xi - 1))/4, &
                      -((2*eta - 1)*(xi - 1)*(xi + 1))/2, &
                      (xi*(2*eta - 1)*(xi + 1))/4, &
                      -eta*xi*(xi - 1), &
                      2*eta*(xi - 1)*(xi + 1), &
                      -eta*xi*(xi + 1), &
                      (xi*(2*eta + 1)*(xi - 1))/4, &
                      -((2*eta + 1)*(xi - 1)*(xi + 1))/2, &
                      (xi*(2*eta + 1)*(xi + 1))/4]
      case (12)
         dpsiNdEta = [-((xi - 1)*(-27*eta**2 + 18*eta - 9*xi**2 + 10))/32, &
                      -(9*(3*xi - 1)*(xi - 1)*(xi + 1))/32, &
                      (9*(3*xi + 1)*(xi - 1)*(xi + 1))/32, &
                      ((xi + 1)*(-27*eta**2 + 18*eta - 9*xi**2 + 10))/32, &
                      (9*(xi - 1)*(-9*eta**2 + 2*eta + 3))/32, &
                      -(9*(xi + 1)*(-9*eta**2 + 2*eta + 3))/32, &
                      (9*(xi - 1)*(9*eta**2 + 2*eta - 3))/32, &
                      -(9*(xi + 1)*(9*eta**2 + 2*eta - 3))/32, &
                      -((xi - 1)*(27*eta**2 + 18*eta + 9*xi**2 - 10))/32, &
                      (9*(3*xi - 1)*(xi - 1)*(xi + 1))/32, &
                      -(9*(3*xi + 1)*(xi - 1)*(xi + 1))/32, &
                      ((xi + 1)*(27*eta**2 + 18*eta + 9*xi**2 - 10))/32]

      case (16)
         dpsiNdEta = [-((3*xi - 1)*(3*xi + 1)*(xi - 1)*(-27*eta**2 + 18*eta + 1))/256, &
                      (9*(3*xi - 1)*(xi - 1)*(xi + 1)*(-27*eta**2 + 18*eta + 1))/256, &
                      -(9*(3*xi + 1)*(xi - 1)*(xi + 1)*(-27*eta**2 + 18*eta + 1))/256, &
                      ((3*xi - 1)*(3*xi + 1)*(xi + 1)*(-27*eta**2 + 18*eta + 1))/256, &
                      (9*(3*xi - 1)*(3*xi + 1)*(xi - 1)*(-9*eta**2 + 2*eta + 3))/256, &
                      -(81*(3*xi - 1)*(xi - 1)*(xi + 1)*(-9*eta**2 + 2*eta + 3))/256, &
                      (81*(3*xi + 1)*(xi - 1)*(xi + 1)*(-9*eta**2 + 2*eta + 3))/256, &
                      -(9*(3*xi - 1)*(3*xi + 1)*(xi + 1)*(-9*eta**2 + 2*eta + 3))/256, &
                      (9*(3*xi - 1)*(3*xi + 1)*(xi - 1)*(9*eta**2 + 2*eta - 3))/256, &
                      -(81*(3*xi - 1)*(xi - 1)*(xi + 1)*(9*eta**2 + 2*eta - 3))/256, &
                      (81*(3*xi + 1)*(xi - 1)*(xi + 1)*(9*eta**2 + 2*eta - 3))/256, &
                      -(9*(3*xi - 1)*(3*xi + 1)*(xi + 1)*(9*eta**2 + 2*eta - 3))/256, &
                      -((3*xi - 1)*(3*xi + 1)*(xi - 1)*(27*eta**2 + 18*eta - 1))/256, &
                      (9*(3*xi - 1)*(xi - 1)*(xi + 1)*(27*eta**2 + 18*eta - 1))/256, &
                      -(9*(3*xi + 1)*(xi - 1)*(xi + 1)*(27*eta**2 + 18*eta - 1))/256, &
                      ((3*xi - 1)*(3*xi + 1)*(xi + 1)*(27*eta**2 + 18*eta - 1))/256]
      case default
         write (*, *) 'Selected Element Type is Unavailable'
         stop

      end select

   end function dpsiNdEta

end module shapeFun
