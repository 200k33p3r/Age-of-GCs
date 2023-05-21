
MODULE define_errors
  implicit none
  ! define the variables needed for the distribution function of  photometric errors

  real, dimension(2), parameter :: SGBartStar = [-13.77,-13.446]  !magnitude of SGB in the
  !artifical stars
  integer, parameter :: nfilter = 2   !number of filters for the HESS diagram
  
  INTEGER :: Nbinrad,Nbinmag,nmaxstars
  integer,dimension(:,:,:), allocatable :: nhist ![binrad,binmag,filter]
  ! number of stars in each error histogram
  
  real, dimension(:,:,:), allocatable :: completness ![binrad,binmag,filter]
  !copmletenes fraction in each error histogram
  
  real, dimension(:,:), allocatable :: binrad,binmag ! [n, filter]
  !bin boundaries for the various error histograms
  
  real, dimension(:,:,:,:),allocatable ::  histMagerr ! [n,binrad,binmag,filter]
  real, dimension(:,:,:,:), allocatable :: cdfy_photerrors ![n,binrad,binmag,filter]
  
END MODULE define_errors
