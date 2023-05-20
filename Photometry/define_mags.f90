MODULE define_mags
 implicit none
 !define the magnitude variables
 !these are the arrays for quantities which are in the fit
 ! region
 real, dimension(:), allocatable :: xin,yin,magin
 real, dimension(:), allocatable :: xout,yout,magout
 real, dimension(:), allocatable :: dist
 logical, dimension(:),allocatable :: recovered
 !these arrays are all of the stars which are read in
 real, dimension(:), allocatable :: xti,yti,xto,yto
END MODULE define_mags