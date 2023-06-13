program RadialDensity  
  use qsort  
  implicit none  !  program to read in ACS treasury photometry file and
  ! find the center of the cluster, by looking at the
  ! average x value of stars within 200 pixels (10'') of a
  ! guess center (by eye), and then output the
  ! distance of each star from the center of the cluster
  ! I use VV for F606W and VI for F606W-F814W filters
  ! just look at the stars used in the fit   !  Variable Declarations
  real, parameter :: xguess = 3000
  real, parameter :: yguess = 3000
  real, parameter :: radius = 200  !200 pixels from guess center
  real, parameter :: Vmax = 20.779
  real, parameter :: Vmin = 13.797
  real, parameter :: VImin = 0.463
  real, parameter :: VImax = 0.916
  integer :: Nstars,nfit,i,j,id,nage
  REAL, dimension(:),allocatable :: xx,yy,fit_id
  REAL, dimension(:),allocatable :: xt,yt,vt,vit,dist2,idlist
  real :: t1,t2,t3,t4,rad2,xcen,ycen,area  
  character(len=40) :: filename
  character(len=10) :: cNstars
  integer :: n_args  !read in Sample Size  from command line
  n_args = command_argument_count()
  if (n_args /= 2) then
     write(*,*)'Usage: ./a.out num_stars  filename'
     stop
  end if
  call get_command_argument(1,cNstars)
  cNstars = trim(cNstars)
  read(cNstars,*) Nstars
  call get_command_argument(2,filename)
  filename = trim(filename)  
  allocate (xt(Nstars),yt(Nstars),vt(Nstars), &
       vit(Nstars), dist2(Nstars), idlist(Nstars)  )  !  Read problem dat and count number of stars used for center fit
  nfit = 0
  nage = 0
  rad2 = radius**2
  open (unit=4, file = filename, status='old')
  read(4,*)    !skip header line
  do i = 1,Nstars
     read(4,*) id, vt(i), t1, t2, t3, vit(i), t4, xt(i), yt(i)
     dist2(i) = (xt(i) - xguess)**2 + (yt(i) -yguess)**2
     idlist(i) = id
     if( dist2(i) <= rad2 .and. vt(i) <= Vmax .and. &
          vt(i) >= Vmin .and. vit(i) >= VImin .and. &
          vit(i) <= VImax ) nfit = nfit + 1
     if( vt(i) <= Vmax .and.  vt(i) >= Vmin  &
         .and. vit(i) >= VImin .and. vit(i) <= VImax ) nage = nage + 1
    enddo
   close(4)
   write(*,*)'Number of stars close to center:', nfit
   write(*,*)'Number of stars for age fit:', nage   
   allocate (xx(nfit),yy(nfit))
   j = 0
   do i = 1,Nstars
      if( dist2(i) <= rad2 .and. vt(i) <= Vmax .and. &
          vt(i) >= Vmin .and. vit(i) >= VImin .and. &
          vit(i) <= VImax ) then
         j = j + 1
         xx(j) = xt(i)
         yy(j) = yt(i)
      endif
   enddo   
   write(*,*)'j, nfit:', j,nfit
   xcen = sum(xx)/nfit
   ycen = sum(yy)/nfit
   write(*,*)'Average x,y:', xcen,ycen
   deallocate(xx,yy, dist2)
   allocate (dist2(nage),fit_id(nage))
   j = 0
   do i = 1, Nstars
      if( vt(i) <= Vmax .and.  vt(i) >= Vmin  &
           .and. vit(i) >= VImin .and. vit(i) <= VImax )  then
         j = j + 1
         dist2(j) = sqrt( (xt(i)-xcen)**2 + (yt(i) - ycen)**2)
         fit_id(j) = idlist(i)
      endif
   enddo   
   call quicksort(dist2)
   open(unit=20, file ='Distance.dat')
   write(20,'(A1,I7)') "#", nage
   write(*,*)'nage=',nage   
   do i = 1,nage
!      write(*,*) i
      write(20,*) dist2(i)
   enddo   !for fun, calculate the stellar density as a function of radius
!   do i = 1, nage
!      !pixscale is 0.05 arcseconds per pixel
!      area = 3.14156* ( (dist2(i)*0.05)**2 )
!      write(*,*) dist2(i)*0.05, real(i)/area
!   enddo   deallocate (xt,yt,vt,vit,dist2) 
   close(unit=20)
   open(unit=21, file='fit_id.dat')
   do i = 1,nage
!      write(*,*) i
      write(21,*) fit_id(i)
   enddo
   close(unit=21)
END program RadialDensity
