module qsort
!  use real_precision  
 implicit none
 contains
! quicksort.f -*-f90-*-
! Author: t-nissie, some tweaks by 1AdAstra1
! License: GPLv3
! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea
!!
  recursive subroutine quicksort(a)
    real  :: a(:)
    real  :: x, t
    integer :: first = 1, last
    integer ::  i, j    
    last = size(a, 1)
    x = a( (first+last) / 2 )
    i = first
    j = last    
    do
       do while (a(i) < x)
          i=i+1
       end do
       do while (x < a(j))
          j=j-1
       end do
       if (i >= j) exit
       t = a(i);  a(i) = a(j);  a(j) = t
       i=i+1
       j=j-1
    end do    

    if (first < i - 1) call quicksort(a(first : i - 1))
    if (j + 1 < last)  call quicksort(a(j + 1 : last))
  end subroutine quicksort
end module qsort
