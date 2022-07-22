program try
   implicit none
   integer:: cnt, i, odd(5), even(5)
   cnt = 1
   do i = 1, 10
      odd(cnt) = i
      i = i + 1
      even(cnt) = i + 1
      cnt = cnt + 1
   end do
end program try
