module shapeFun
   use userInput, only: degEl, flag
   implicit none

   integer, parameter, private::n = (degEl + 1)**2 - flag*(degEl - 1)**2 !nodes per element
   double precision, parameter:: &
      d1 = 1.0d0, &
      d2 = 2.0d0, &
      d3 = 3.0d0, &
      d4 = 4.0d0, &
      d5 = 5.0d0, &
      d6 = 6.0d0, &
      d8 = 8.0d0, &
      d9 = 9.0d0, &
      d10 = 10.0d0, &
      d12 = 12.0d0, &
      d16 = 16.0d0, &
      d18 = 18.0d0, &
      d27 = 27.0d0, &
      d32 = 32.0d0, &
      d36 = 36.0d0, &
      d81 = 81.0d0, &
      d256 = 256.0d0

contains

   function psiN(eta, xi)
      implicit none

      double precision, allocatable, dimension(:)::psiN
      double precision, intent(in) :: eta, xi
      allocate (psiN(n))

      select case (n)
      case (4)
         psiN = [((eta - d1)*(xi - d1))/d4, &
                 -((eta - d1)*(xi + d1))/d4, &
                 -((eta + d1)*(xi - d1))/d4, &
                 ((eta + d1)*(xi + d1))/d4]
      case (8)
         psiN = [-((eta - d1)*(xi - d1)*(eta + xi + d1))/d4, &
                 ((eta - d1)*(xi - d1)*(xi + d1))/d2, &
                 ((eta - d1)*(xi + d1)*(eta - xi + d1))/d4, &
                 ((eta - d1)*(eta + d1)*(xi - d1))/d2, &
                 -((eta - d1)*(eta + d1)*(xi + d1))/d2, &
                 ((eta + d1)*(xi - d1)*(xi - eta + d1))/d4, &
                 -((eta + d1)*(xi - d1)*(xi + d1))/d2, &
                 ((eta + d1)*(xi + d1)*(eta + xi - d1))/d4]
      case (9)
         psiN = [(eta*xi*(eta - d1)*(xi - d1))/d4, &
                 -(eta*(eta - d1)*(xi - d1)*(xi + d1))/d2, &
                 (eta*xi*(eta - d1)*(xi + d1))/d4, &
                 -(xi*(eta - d1)*(eta + d1)*(xi - d1))/d2, &
                 (eta - d1)*(eta + d1)*(xi - d1)*(xi + d1), &
                 -(xi*(eta - d1)*(eta + d1)*(xi + d1))/d2, &
                 (eta*xi*(eta + d1)*(xi - d1))/d4, &
                 -(eta*(eta + d1)*(xi - d1)*(xi + d1))/d2, &
                 (eta*xi*(eta + d1)*(xi + d1))/d4]
      case (12)
         psiN = [((eta - d1)*(xi - d1)*(d9*eta**d2 + d9*xi**d2 - d10))/d32, &
                 -(d9*(d3*xi - d1)*(eta - d1)*(xi - d1)*(xi + d1))/d32, &
                 (d9*(d3*xi + d1)*(eta - d1)*(xi - d1)*(xi + d1))/d32, &
                 -((eta - d1)*(xi + d1)*(d9*eta**d2 + d9*xi**d2 - d10))/d32, &
                 -(d9*(d3*eta - d1)*(eta - d1)*(eta + d1)*(xi - d1))/d32, &
                 (d9*(d3*eta - d1)*(eta - d1)*(eta + d1)*(xi + d1))/d32, &
                 (d9*(d3*eta + d1)*(eta - d1)*(eta + d1)*(xi - d1))/d32, &
                 -(d9*(d3*eta + d1)*(eta - d1)*(eta + d1)*(xi + d1))/d32, &
                 -((eta + d1)*(xi - d1)*(d9*eta**d2 + d9*xi**d2 - d10))/d32, &
                 (d9*(d3*xi - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d32, &
                 -(d9*(d3*xi + d1)*(eta + d1)*(xi - d1)*(xi + d1))/d32, &
                 ((eta + d1)*(xi + d1)*(d9*eta**d2 + d9*xi**d2 - d10))/d32]
      case (16)
         psiN = [((d3*eta - d1)*(d3*eta + d1)*(d3*xi - d1)*(d3*xi + d1)*(eta - d1)*(xi - d1))/d256, &
                 -(d9*(d3*eta - d1)*(d3*eta + d1)*(d3*xi - d1)*(eta - d1)*(xi - d1)*(xi + d1))/d256, &
                 (d9*(d3*eta - d1)*(d3*eta + d1)*(d3*xi + d1)*(eta - d1)*(xi - d1)*(xi + d1))/d256, &
                 -((d3*eta - d1)*(d3*eta + d1)*(d3*xi - d1)*(d3*xi + d1)*(eta - d1)*(xi + d1))/d256, &
                 -(d9*(d3*eta - d1)*(d3*xi - d1)*(d3*xi + d1)*(eta - d1)*(eta + d1)*(xi - d1))/d256, &
                 (d81*(d3*eta - d1)*(d3*xi - d1)*(eta - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d256, &
                 -(d81*(d3*eta - d1)*(d3*xi + d1)*(eta - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d256, &
                 (d9*(d3*eta - d1)*(d3*xi - d1)*(d3*xi + d1)*(eta - d1)*(eta + d1)*(xi + d1))/d256, &
                 (d9*(d3*eta + d1)*(d3*xi - d1)*(d3*xi + d1)*(eta - d1)*(eta + d1)*(xi - d1))/d256, &
                 -(d81*(d3*eta + d1)*(d3*xi - d1)*(eta - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d256, &
                 (d81*(d3*eta + d1)*(d3*xi + d1)*(eta - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d256, &
                 -(d9*(d3*eta + d1)*(d3*xi - d1)*(d3*xi + d1)*(eta - d1)*(eta + d1)*(xi + d1))/d256, &
                 -((d3*eta - d1)*(d3*eta + d1)*(d3*xi - d1)*(d3*xi + d1)*(eta + d1)*(xi - d1))/d256, &
                 (d9*(d3*eta - d1)*(d3*eta + d1)*(d3*xi - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d256, &
                 -(d9*(d3*eta - d1)*(d3*eta + d1)*(d3*xi + d1)*(eta + d1)*(xi - d1)*(xi + d1))/d256, &
                 ((d3*eta - d1)*(d3*eta + d1)*(d3*xi - d1)*(d3*xi + d1)*(eta + d1)*(xi + d1))/d256]
      case (25)
         psiN = [(eta*xi*(d2*eta - d1)*(d2*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(eta - d1)*(xi - d1))/d36, &
                 -(d2*eta*xi*(d2*eta - d1)*(d2*eta + d1)*(d2*xi - d1)*(eta - d1)*(xi - d1)*(xi + d1))/d9, &
                 (eta*(d2*eta - d1)*(d2*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(eta - d1)*(xi - d1)*(xi + d1))/d6, &
                 -(d2*eta*xi*(d2*eta - d1)*(d2*eta + d1)*(d2*xi + d1)*(eta - d1)*(xi - d1)*(xi + d1))/d9, &
                 (eta*xi*(d2*eta - d1)*(d2*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(eta - d1)*(xi + d1))/d36, &
                 -(d2*eta*xi*(d2*eta - d1)*(d2*xi - d1)*(d2*xi + d1)*(eta - d1)*(eta + d1)*(xi - d1))/d9, &
                 (d16*eta*xi*(d2*eta - d1)*(d2*xi - d1)*(eta - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d9, &
                 -(d4*eta*(d2*eta - d1)*(d2*xi - d1)*(d2*xi + d1)*(eta - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d3, &
                 (d16*eta*xi*(d2*eta - d1)*(d2*xi + d1)*(eta - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d9, &
                 -(d2*eta*xi*(d2*eta - d1)*(d2*xi - d1)*(d2*xi + d1)*(eta - d1)*(eta + d1)*(xi + d1))/d9, &
                 (xi*(d2*eta - d1)*(d2*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(eta - d1)*(eta + d1)*(xi - d1))/d6, &
                 -(d4*xi*(d2*eta - d1)*(d2*eta + d1)*(d2*xi - d1)*(eta - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d3, &
                 (d2*eta - d1)*(d2*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(eta - d1)*(eta + d1)*(xi - d1)*(xi + d1), &
                 -(d4*xi*(d2*eta - d1)*(d2*eta + d1)*(d2*xi + d1)*(eta - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d3, &
                 (xi*(d2*eta - d1)*(d2*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(eta - d1)*(eta + d1)*(xi + d1))/d6, &
                 -(d2*eta*xi*(d2*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(eta - d1)*(eta + d1)*(xi - d1))/d9, &
                 (d16*eta*xi*(d2*eta + d1)*(d2*xi - d1)*(eta - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d9, &
                 -(d4*eta*(d2*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(eta - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d3, &
                 (d16*eta*xi*(d2*eta + d1)*(d2*xi + d1)*(eta - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d9, &
                 -(d2*eta*xi*(d2*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(eta - d1)*(eta + d1)*(xi + d1))/d9, &
                 (eta*xi*(d2*eta - d1)*(d2*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(eta + d1)*(xi - d1))/d36, &
                 -(d2*eta*xi*(d2*eta - d1)*(d2*eta + d1)*(d2*xi - d1)*(eta + d1)*(xi - d1)*(xi + d1))/d9, &
                 (eta*(d2*eta - d1)*(d2*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(eta + d1)*(xi - d1)*(xi + d1))/d6, &
                 -(d2*eta*xi*(d2*eta - d1)*(d2*eta + d1)*(d2*xi + d1)*(eta + d1)*(xi - d1)*(xi + d1))/d9, &
                 (eta*xi*(d2*eta - d1)*(d2*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(eta + d1)*(xi + d1))/d36]
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
         dpsiNdXi = [eta/d4 - d1/d4, &
                     d1/d4 - eta/d4, &
                     -eta/d4 - d1/d4, &
                     eta/d4 + d1/d4]
      case (8)
         dpsiNdXi = [-((eta + d2*xi)*(eta - d1))/d4, &
                     xi*(eta - d1), &
                     ((eta - d2*xi)*(eta - d1))/d4, &
                     ((eta - d1)*(eta + d1))/d2, &
                     -((eta - d1)*(eta + d1))/d2, &
                     -((eta - d2*xi)*(eta + d1))/d4, &
                     -xi*(eta + d1), &
                     ((eta + d2*xi)*(eta + d1))/d4]
      case (9)
         dpsiNdXi = [(eta*(d2*xi - d1)*(eta - d1))/d4, &
                     -eta*xi*(eta - d1), &
                     (eta*(d2*xi + d1)*(eta - d1))/d4, &
                     -((d2*xi - d1)*(eta - d1)*(eta + d1))/d2, &
                     d2*xi*(eta - d1)*(eta + d1), &
                     -((d2*xi + d1)*(eta - d1)*(eta + d1))/d2, &
                     (eta*(d2*xi - d1)*(eta + d1))/d4, &
                     -eta*xi*(eta + d1), &
                     (eta*(d2*xi + d1)*(eta + d1))/d4]
      case (12)
         dpsiNdXi = [-((eta - d1)*(-d9*eta**d2 - d27*xi**d2 + d18*xi + d10))/d32, &
                     (d9*(eta - d1)*(-d9*xi**d2 + d2*xi + d3))/d32, &
                     (d9*(eta - d1)*(d9*xi**d2 + d2*xi - d3))/d32, &
                     -((eta - d1)*(d9*eta**d2 + d27*xi**d2 + d18*xi - d10))/d32, &
                     -(d9*(d3*eta - d1)*(eta - d1)*(eta + d1))/d32, &
                     (d9*(d3*eta - d1)*(eta - d1)*(eta + d1))/d32, &
                     (d9*(d3*eta + d1)*(eta - d1)*(eta + d1))/d32, &
                     -(d9*(d3*eta + d1)*(eta - d1)*(eta + d1))/d32, &
                     ((eta + d1)*(-d9*eta**d2 - d27*xi**d2 + d18*xi + d10))/d32, &
                     -(d9*(eta + d1)*(-d9*xi**d2 + d2*xi + d3))/d32, &
                     -(d9*(eta + d1)*(d9*xi**d2 + d2*xi - d3))/d32, &
                     ((eta + d1)*(d9*eta**d2 + d27*xi**d2 + d18*xi - d10))/d32]
      case (16)
         dpsiNdXi = [-((d3*eta - d1)*(d3*eta + d1)*(eta - d1)*(-d27*xi**d2 + d18*xi + d1))/d256, &
                     (d9*(d3*eta - d1)*(d3*eta + d1)*(eta - d1)*(-d9*xi**d2 + d2*xi + d3))/d256, &
                     (d9*(d3*eta - d1)*(d3*eta + d1)*(eta - d1)*(d9*xi**d2 + d2*xi - d3))/d256, &
                     -((d3*eta - d1)*(d3*eta + d1)*(eta - d1)*(d27*xi**d2 + d18*xi - d1))/d256, &
                     (d9*(d3*eta - d1)*(eta - d1)*(eta + d1)*(-d27*xi**d2 + d18*xi + d1))/d256, &
                     -(d81*(d3*eta - d1)*(eta - d1)*(eta + d1)*(-d9*xi**d2 + d2*xi + d3))/d256, &
                     -(d81*(d3*eta - d1)*(eta - d1)*(eta + d1)*(d9*xi**d2 + d2*xi - d3))/d256, &
                     (d9*(d3*eta - d1)*(eta - d1)*(eta + d1)*(d27*xi**d2 + d18*xi - d1))/d256, &
                     -(d9*(d3*eta + d1)*(eta - d1)*(eta + d1)*(-d27*xi**d2 + d18*xi + d1))/d256, &
                     (d81*(d3*eta + d1)*(eta - d1)*(eta + d1)*(-d9*xi**d2 + d2*xi + d3))/d256, &
                     (d81*(d3*eta + d1)*(eta - d1)*(eta + d1)*(d9*xi**d2 + d2*xi - d3))/d256, &
                     -(d9*(d3*eta + d1)*(eta - d1)*(eta + d1)*(d27*xi**d2 + d18*xi - d1))/d256, &
                     ((d3*eta - d1)*(d3*eta + d1)*(eta + d1)*(-d27*xi**d2 + d18*xi + d1))/d256, &
                     -(d9*(d3*eta - d1)*(d3*eta + d1)*(eta + d1)*(-d9*xi**d2 + d2*xi + d3))/d256, &
                     -(d9*(d3*eta - d1)*(d3*eta + d1)*(eta + d1)*(d9*xi**d2 + d2*xi - d3))/d256, &
                     ((d3*eta - d1)*(d3*eta + d1)*(eta + d1)*(d27*xi**d2 + d18*xi - d1))/d256]
      case (25)
         dpsiNdXi = [-(eta*(d2*eta - d1)*(d2*eta + d1)*(d4*xi - d1)*(eta - d1)*(-d4*xi**d2 + d2*xi + d1))/d36, &
                     (d2*eta*(d2*eta - d1)*(d2*eta + d1)*(eta - d1)*(-d8*xi**d3 + d3*xi**d2 + d4*xi - d1))/d9, &
                     (eta*xi*(d2*eta - d1)*(d2*eta + d1)*(d8*xi**d2 - d5)*(eta - d1))/d3, &
                     (d2*eta*(d2*eta - d1)*(d2*eta + d1)*(eta - d1)*(-d8*xi**d3 - d3*xi**d2 + d4*xi + d1))/d9, &
                     (eta*(d2*eta - d1)*(d2*eta + d1)*(d4*xi + d1)*(eta - d1)*(d4*xi**d2 + d2*xi - d1))/d36, &
                     (d2*eta*(d2*eta - d1)*(d4*xi - d1)*(eta - d1)*(eta + d1)*(-d4*xi**d2 + d2*xi + d1))/d9, &
                     -(d16*eta*(d2*eta - d1)*(eta - d1)*(eta + d1)*(-d8*xi**d3 + d3*xi**d2 + d4*xi - d1))/d9, &
                     -(d8*eta*xi*(d2*eta - d1)*(d8*xi**d2 - d5)*(eta - d1)*(eta + d1))/d3, &
                     -(d16*eta*(d2*eta - d1)*(eta - d1)*(eta + d1)*(-d8*xi**d3 - d3*xi**d2 + d4*xi + d1))/d9, &
                     -(d2*eta*(d2*eta - d1)*(d4*xi + d1)*(eta - d1)*(eta + d1)*(d4*xi**d2 + d2*xi - d1))/d9, &
                     -((d2*eta - d1)*(d2*eta + d1)*(d4*xi - d1)*(eta - d1)*(eta + d1)*(-d4*xi**d2 + d2*xi + d1))/d6, &
                     (d4*(d2*eta - d1)*(d2*eta + d1)*(eta - d1)*(eta + d1)*(-d8*xi**d3 + d3*xi**d2 + d4*xi - d1))/d3, &
                     d2*xi*(d2*eta - d1)*(d2*eta + d1)*(d8*xi**d2 - d5)*(eta - d1)*(eta + d1), &
                     (d4*(d2*eta - d1)*(d2*eta + d1)*(eta - d1)*(eta + d1)*(-d8*xi**d3 - d3*xi**d2 + d4*xi + d1))/d3, &
                     ((d2*eta - d1)*(d2*eta + d1)*(d4*xi + d1)*(eta - d1)*(eta + d1)*(d4*xi**d2 + d2*xi - d1))/d6, &
                     (d2*eta*(d2*eta + d1)*(d4*xi - d1)*(eta - d1)*(eta + d1)*(-d4*xi**d2 + d2*xi + d1))/d9, &
                     -(d16*eta*(d2*eta + d1)*(eta - d1)*(eta + d1)*(-d8*xi**d3 + d3*xi**d2 + d4*xi - d1))/d9, &
                     -(d8*eta*xi*(d2*eta + d1)*(d8*xi**d2 - d5)*(eta - d1)*(eta + d1))/d3, &
                     -(d16*eta*(d2*eta + d1)*(eta - d1)*(eta + d1)*(-d8*xi**d3 - d3*xi**d2 + d4*xi + d1))/d9, &
                     -(d2*eta*(d2*eta + d1)*(d4*xi + d1)*(eta - d1)*(eta + d1)*(d4*xi**d2 + d2*xi - d1))/d9, &
                     -(eta*(d2*eta - d1)*(d2*eta + d1)*(d4*xi - d1)*(eta + d1)*(-d4*xi**d2 + d2*xi + d1))/d36, &
                     (d2*eta*(d2*eta - d1)*(d2*eta + d1)*(eta + d1)*(-d8*xi**d3 + d3*xi**d2 + d4*xi - d1))/d9, &
                     (eta*xi*(d2*eta - d1)*(d2*eta + d1)*(d8*xi**d2 - d5)*(eta + d1))/d3, &
                     (d2*eta*(d2*eta - d1)*(d2*eta + d1)*(eta + d1)*(-d8*xi**d3 - d3*xi**d2 + d4*xi + d1))/d9, &
                     (eta*(d2*eta - d1)*(d2*eta + d1)*(d4*xi + d1)*(eta + d1)*(d4*xi**d2 + d2*xi - d1))/d36]
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
         dpsiNdEta = [xi/d4 - d1/d4, &
                      -xi/d4 - d1/d4, &
                      d1/d4 - xi/d4, &
                      xi/d4 + d1/d4]
      case (8)
         dpsiNdEta = [-((d2*eta + xi)*(xi - d1))/d4, &
                      ((xi - d1)*(xi + d1))/d2, &
                      ((xi + d1)*(d2*eta - xi))/d4, &
                      eta*(xi - d1), &
                      -eta*(xi + d1), &
                      -((xi - d1)*(d2*eta - xi))/d4, &
                      -((xi - d1)*(xi + d1))/d2, &
                      ((d2*eta + xi)*(xi + d1))/d4]
      case (9)
         dpsiNdEta = [(xi*(d2*eta - d1)*(xi - d1))/d4, &
                      -((d2*eta - d1)*(xi - d1)*(xi + d1))/d2, &
                      (xi*(d2*eta - d1)*(xi + d1))/d4, &
                      -eta*xi*(xi - d1), &
                      d2*eta*(xi - d1)*(xi + d1), &
                      -eta*xi*(xi + d1), &
                      (xi*(d2*eta + d1)*(xi - d1))/d4, &
                      -((d2*eta + d1)*(xi - d1)*(xi + d1))/d2, &
                      (xi*(d2*eta + d1)*(xi + d1))/d4]
      case (12)
         dpsiNdEta = [-((xi - d1)*(-d27*eta**d2 + d18*eta - d9*xi**d2 + d10))/d32, &
                      -(d9*(d3*xi - d1)*(xi - d1)*(xi + d1))/d32, &
                      (d9*(d3*xi + d1)*(xi - d1)*(xi + d1))/d32, &
                      ((xi + d1)*(-d27*eta**d2 + d18*eta - d9*xi**d2 + d10))/d32, &
                      (d9*(xi - d1)*(-d9*eta**d2 + d2*eta + d3))/d32, &
                      -(d9*(xi + d1)*(-d9*eta**d2 + d2*eta + d3))/d32, &
                      (d9*(xi - d1)*(d9*eta**d2 + d2*eta - d3))/d32, &
                      -(d9*(xi + d1)*(d9*eta**d2 + d2*eta - d3))/d32, &
                      -((xi - d1)*(d27*eta**d2 + d18*eta + d9*xi**d2 - d10))/d32, &
                      (d9*(d3*xi - d1)*(xi - d1)*(xi + d1))/d32, &
                      -(d9*(d3*xi + d1)*(xi - d1)*(xi + d1))/d32, &
                      ((xi + d1)*(d27*eta**d2 + d18*eta + d9*xi**d2 - d10))/d32]

      case (16)
         dpsiNdEta = [-((d3*xi - d1)*(d3*xi + d1)*(xi - d1)*(-d27*eta**d2 + d18*eta + d1))/d256, &
                      (d9*(d3*xi - d1)*(xi - d1)*(xi + d1)*(-d27*eta**d2 + d18*eta + d1))/d256, &
                      -(d9*(d3*xi + d1)*(xi - d1)*(xi + d1)*(-d27*eta**d2 + d18*eta + d1))/d256, &
                      ((d3*xi - d1)*(d3*xi + d1)*(xi + d1)*(-d27*eta**d2 + d18*eta + d1))/d256, &
                      (d9*(d3*xi - d1)*(d3*xi + d1)*(xi - d1)*(-d9*eta**d2 + d2*eta + d3))/d256, &
                      -(d81*(d3*xi - d1)*(xi - d1)*(xi + d1)*(-d9*eta**d2 + d2*eta + d3))/d256, &
                      (d81*(d3*xi + d1)*(xi - d1)*(xi + d1)*(-d9*eta**d2 + d2*eta + d3))/d256, &
                      -(d9*(d3*xi - d1)*(d3*xi + d1)*(xi + d1)*(-d9*eta**d2 + d2*eta + d3))/d256, &
                      (d9*(d3*xi - d1)*(d3*xi + d1)*(xi - d1)*(d9*eta**d2 + d2*eta - d3))/d256, &
                      -(d81*(d3*xi - d1)*(xi - d1)*(xi + d1)*(d9*eta**d2 + d2*eta - d3))/d256, &
                      (d81*(d3*xi + d1)*(xi - d1)*(xi + d1)*(d9*eta**d2 + d2*eta - d3))/d256, &
                      -(d9*(d3*xi - d1)*(d3*xi + d1)*(xi + d1)*(d9*eta**d2 + d2*eta - d3))/d256, &
                      -((d3*xi - d1)*(d3*xi + d1)*(xi - d1)*(d27*eta**d2 + d18*eta - d1))/d256, &
                      (d9*(d3*xi - d1)*(xi - d1)*(xi + d1)*(d27*eta**d2 + d18*eta - d1))/d256, &
                      -(d9*(d3*xi + d1)*(xi - d1)*(xi + d1)*(d27*eta**d2 + d18*eta - d1))/d256, &
                      ((d3*xi - d1)*(d3*xi + d1)*(xi + d1)*(d27*eta**d2 + d18*eta - d1))/d256]
      case (25)
         dpsiNdEta = [-(xi*(d4*eta - d1)*(d2*xi - d1)*(d2*xi + d1)*(xi - d1)*(-d4*eta**d2 + d2*eta + d1))/d36, &
                      (d2*xi*(d4*eta - d1)*(d2*xi - d1)*(xi - d1)*(xi + d1)*(-d4*eta**d2 + d2*eta + d1))/d9, &
                      -((d4*eta - d1)*(d2*xi - d1)*(d2*xi + d1)*(xi - d1)*(xi + d1)*(-d4*eta**d2 + d2*eta + d1))/d6, &
                      (d2*xi*(d4*eta - d1)*(d2*xi + d1)*(xi - d1)*(xi + d1)*(-d4*eta**d2 + d2*eta + d1))/d9, &
                      -(xi*(d4*eta - d1)*(d2*xi - d1)*(d2*xi + d1)*(xi + d1)*(-d4*eta**d2 + d2*eta + d1))/d36, &
                      (d2*xi*(d2*xi - d1)*(d2*xi + d1)*(xi - d1)*(-d8*eta**d3 + d3*eta**d2 + d4*eta - d1))/d9, &
                      -(d16*xi*(d2*xi - d1)*(xi - d1)*(xi + d1)*(-d8*eta**d3 + d3*eta**d2 + d4*eta - d1))/d9, &
                      (d4*(d2*xi - d1)*(d2*xi + d1)*(xi - d1)*(xi + d1)*(-d8*eta**d3 + d3*eta**d2 + d4*eta - d1))/d3, &
                      -(d16*xi*(d2*xi + d1)*(xi - d1)*(xi + d1)*(-d8*eta**d3 + d3*eta**d2 + d4*eta - d1))/d9, &
                      (d2*xi*(d2*xi - d1)*(d2*xi + d1)*(xi + d1)*(-d8*eta**d3 + d3*eta**d2 + d4*eta - d1))/d9, &
                      (eta*xi*(d2*xi - d1)*(d2*xi + d1)*(d8*eta**d2 - d5)*(xi - d1))/d3, &
                      -(d8*eta*xi*(d2*xi - d1)*(d8*eta**d2 - d5)*(xi - d1)*(xi + d1))/d3, &
                      d2*eta*(d2*xi - d1)*(d2*xi + d1)*(d8*eta**d2 - d5)*(xi - d1)*(xi + d1), &
                      -(d8*eta*xi*(d2*xi + d1)*(d8*eta**d2 - d5)*(xi - d1)*(xi + d1))/d3, &
                      (eta*xi*(d2*xi - d1)*(d2*xi + d1)*(d8*eta**d2 - d5)*(xi + d1))/d3, &
                      (d2*xi*(d2*xi - d1)*(d2*xi + d1)*(xi - d1)*(-d8*eta**d3 - d3*eta**d2 + d4*eta + d1))/d9, &
                      -(d16*xi*(d2*xi - d1)*(xi - d1)*(xi + d1)*(-d8*eta**d3 - d3*eta**d2 + d4*eta + d1))/d9, &
                      (d4*(d2*xi - d1)*(d2*xi + d1)*(xi - d1)*(xi + d1)*(-d8*eta**d3 - d3*eta**d2 + d4*eta + d1))/d3, &
                      -(d16*xi*(d2*xi + d1)*(xi - d1)*(xi + d1)*(-d8*eta**d3 - d3*eta**d2 + d4*eta + d1))/d9, &
                      (d2*xi*(d2*xi - d1)*(d2*xi + d1)*(xi + d1)*(-d8*eta**d3 - d3*eta**d2 + d4*eta + d1))/d9, &
                      (xi*(d4*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(xi - d1)*(d4*eta**d2 + d2*eta - d1))/d36, &
                      -(d2*xi*(d4*eta + d1)*(d2*xi - d1)*(xi - d1)*(xi + d1)*(d4*eta**d2 + d2*eta - d1))/d9, &
                      ((d4*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(xi - d1)*(xi + d1)*(d4*eta**d2 + d2*eta - d1))/d6, &
                      -(d2*xi*(d4*eta + d1)*(d2*xi + d1)*(xi - d1)*(xi + d1)*(d4*eta**d2 + d2*eta - d1))/d9, &
                      (xi*(d4*eta + d1)*(d2*xi - d1)*(d2*xi + d1)*(xi + d1)*(d4*eta**d2 + d2*eta - d1))/d36]
      case default
         write (*, *) 'Selected Element Type is Unavailable'
         stop

      end select

   end function dpsiNdEta

end module shapeFun
