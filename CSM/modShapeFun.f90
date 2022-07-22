module shapeFun
  use userInput, only: degEl,flag
  implicit none

  integer,parameter,private::n=(degEl+1)**2 - flag*(degEl-1)**2 !nodes per element

  contains

  function psiN(eta,xi)
    implicit none

    double precision,dimension(n)::psiN
    double precision,intent(in) :: eta,xi

    psiN =[ (eta*xi*(eta - 1)*(xi - 1))/4,&
    -(eta*(eta - 1)*(xi - 1)*(xi + 1))/2,&
    (eta*xi*(eta - 1)*(xi + 1))/4,&
    -(xi*(eta - 1)*(eta + 1)*(xi - 1))/2,&
    (eta - 1)*(eta + 1)*(xi - 1)*(xi + 1),&
    -(xi*(eta - 1)*(eta + 1)*(xi + 1))/2,&
    (eta*xi*(eta + 1)*(xi - 1))/4,&
    -(eta*(eta + 1)*(xi - 1)*(xi + 1))/2,&
    (eta*xi*(eta + 1)*(xi + 1))/4]
  end function psiN

  function dpsiNdXi(eta,xi)
    implicit none

    double precision,dimension(n)::dpsiNdXi
    double precision,intent(in) :: eta,xi

    dpsiNdXi =[ (eta*xi*(eta - 1))/4 + (eta*(eta - 1)*(xi - 1))/4,&
    - (eta*(eta - 1)*(xi - 1))/2 - (eta*(eta - 1)*(xi + 1))/2,&
    (eta*xi*(eta - 1))/4 + (eta*(eta - 1)*(xi + 1))/4,&
    - (xi*(eta - 1)*(eta + 1))/2 - ((eta - 1)*(eta + 1)*(xi - 1))/2,&
    (eta - 1)*(eta + 1)*(xi - 1) + (eta - 1)*(eta + 1)*(xi + 1),&
    - (xi*(eta - 1)*(eta + 1))/2 - ((eta - 1)*(eta + 1)*(xi + 1))/2,&
    (eta*xi*(eta + 1))/4 + (eta*(eta + 1)*(xi - 1))/4,&
    - (eta*(eta + 1)*(xi - 1))/2 - (eta*(eta + 1)*(xi + 1))/2,&
    (eta*xi*(eta + 1))/4 + (eta*(eta + 1)*(xi + 1))/4]
  end function dpsiNdXi

  function dpsiNdEta(eta,xi)
    implicit none

    double precision,dimension(n)::dpsiNdEta
    double precision,intent(in) :: eta,xi

    dpsiNdEta=[ (eta*xi*(xi - 1))/4 + (xi*(eta - 1)*(xi - 1))/4,&
    - (eta*(xi - 1)*(xi + 1))/2 - ((eta - 1)*(xi - 1)*(xi + 1))/2,&
    (eta*xi*(xi + 1))/4 + (xi*(eta - 1)*(xi + 1))/4,&
    - (xi*(eta - 1)*(xi - 1))/2 - (xi*(eta + 1)*(xi - 1))/2,&
    (eta - 1)*(xi - 1)*(xi + 1) + (eta + 1)*(xi - 1)*(xi + 1),&
    - (xi*(eta - 1)*(xi + 1))/2 - (xi*(eta + 1)*(xi + 1))/2,&
    (eta*xi*(xi - 1))/4 + (xi*(eta + 1)*(xi - 1))/4,&
    - (eta*(xi - 1)*(xi + 1))/2 - ((eta + 1)*(xi - 1)*(xi + 1))/2,&
    (eta*xi*(xi + 1))/4 + (xi*(eta + 1)*(xi + 1))/4]
  end function dpsiNdEta

end module shapeFun
