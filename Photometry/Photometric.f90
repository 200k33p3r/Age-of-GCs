program PhotometricProperties
  implicit none  
  !  program to read in ACS treasury artificial star file and
  ! characterize the completenss, and photometric uncerainty
  ! as a function of distance from center of the cluster
  ! and F606W magnitude for stars around the main sequence
  ! turn-off (defined by Vmin, Vmax, VImin and VImax).
  ! I use VV for F606W and VI for F606W-F814W filters
  ! assumses a guess center of cluster is at  3000,3000   !  Variable Declarations
  real, parameter :: Vmax = -11.95
  real, parameter :: Vmin = -15.9493
  real, parameter :: VImin = -0.893
  real, parameter :: VImax = 1.255
  real, parameter :: Vshort_long = -14.18 !magnitude between short/long
                     !exposure
  integer, parameter :: Nmagbins = 12 !totalnumber of mag bins
  integer, parameter :: Nradbins = 10 !number of radial bins  
  integer :: Nstars,nfit,i,j,id
  REAL, dimension(:),allocatable :: xin,yin,vin,viin,  &
       xout,yout,vout,viout,dist
  real, dimension(Nradbins,Nmagbins) ::  completeness
  real, dimension(Nradbins) :: binrad
  real, dimension(Nmagbins) :: binmag
  logical, dimension(:),allocatable :: recovered
  REAL, dimension(:),allocatable :: xti,yti,vti,viti, xto,yto, &
       vto, vito  !temp arrays
  real :: ii_in,ii_out,xcenter,ycenter
  integer :: i1,i2,i3,i4  
  character(len=40) :: filename
  character(len=10) :: cNstars
  integer :: n_args  !read in Sample Size & filename from command line
  n_args = command_argument_count()
  if (n_args /= 2) then
     write(*,*)'Usage: ./a.out num_stars filename'
     stop
  end if
  call get_command_argument(1,cNstars)
  cNstars = trim(cNstars)
  read(cNstars,*) Nstars
  call get_command_argument(2,filename)
  filename = trim(filename)  
  allocate (xti(Nstars),yti(Nstars),vti(Nstars), &
       viti(Nstars), xto(Nstars),yto(Nstars), vto(Nstars), &
       vito(Nstars) )  !Read in data and count number of stars in the MSTO fit region
  nfit = 0
  open (unit=4, file = filename, status='old')
  read(4,*)    !skip header line
  do i = 1,Nstars
     read(4,*) xti(i),yti(i),vti(i),ii_in, &
          xto(i),yto(i),vto(i),ii_out
     viti(i) = vti(i) - ii_in
     vito(i) = vto(i) - ii_out
     if(  vti(i) <= Vmax .and. vti(i) >= Vmin .and. &
          viti(i) >= VImin .and. viti(i) <= VImax )   nfit = nfit + 1
   enddo
   close(4)
!   write(*,*)'Number of stars in MSTO region:', nfit   
   allocate (xin(nfit),yin(nfit),vin(nfit),viin(nfit), &
        xout(nfit),yout(nfit),vout(nfit),viout(nfit), &
        dist(nfit), recovered(nfit) )   
   j = 0
   do i = 1,Nstars
      if(  vti(i) <= Vmax .and. vti(i) >= Vmin .and. &
           viti(i) >= VImin .and. viti(i) <= VImax ) then
         j = j + 1
         xin(j) = xti(i)
         yin(j) = yti(i)
         vin(j) = vti(i)
         viin(j)= viti(i)
         xout(j) = xto(i)
         yout(j) = yto(i)
         vout(j) = vto(i)
         viout(j)= vito(i)
      endif
   enddo
!   write(*,*)'nfit,j:', nfit, j   
   deallocate (xti,yti,vti,viti,xto,yto,vto,vito)   
!find out which stars are bad/missing
   call FindMissing(xin,yin,xout,yout,vin,viin,vout,viout,nfit,recovered)   
!find the center of the input stars, and distance from center
   call FindCenter(xin,yin,dist,xcenter,ycenter,nfit)   
! find completenss and magnitude error as a function of distance from the
   ! cluster center and magnitude.  The routine  sets bin boundaries, and
   ! outputs the results.
   call GetPhotProp(dist,vin,viin,vout,viout,recovered,completeness, &
        binrad,binmag, Vshort_long, Nradbins,Nmagbins,nfit) 
END program PhotometricProperties 
! -----------------------------------------------------------------
 subroutine GetPhotProp(dist,vv,viin,vout,viout, recovered,completeness, &
      binrad,binmag, Vcut, Nradbins,Nmagbins,npt)   
   implicit none
   integer :: npt,Nradbins,Nmagbins
   real, dimension(npt) :: dist, vv,viin,vout,viout
   logical, dimension(npt) :: recovered
   real, dimension(Nradbins,Nmagbins) ::  completeness
   real, dimension(Nradbins) :: binrad
   real, dimension(Nmagbins) :: binmag
   real,intent(IN) :: Vcut   
   integer :: i,j,indxrad,indxV , outfile,k
   real :: maxdist,width1,width2,minv,maxv,widthrad
   integer, dimension(Nradbins,Nmagbins) :: ngood,nbad
   real, dimension(Nradbins,Nmagbins) :: sigV,sigVI
   real, dimension(Nradbins,Nmagbins,800) :: histVerr,histVIerr
   character(len=30) :: filename1, filename2   
   maxdist = maxval(dist)
!   write(*,*)maxdist
   widthrad = maxdist/Nradbins   
   binrad(1) = widthrad*1.5
   binrad(Nradbins) = maxdist
   maxdist = maxdist- 3.0*widthrad
   widthrad = maxdist/Nradbins
   do i = 2,Nradbins-1
      binrad(i) = binrad(i-1) + widthrad
   enddo
!   do i = 1, Nradbins
!      write(*,*)'i,binrad:',i,binrad(i)
!   enddo   
   if (mod(Nmagbins,2) .ne. 0 ) then
      write(*,*)'Nmagbins must be an even number'
      stop
   endif
   minv = minval(vv)
   maxv = maxval(vv)
   width1 = (Vcut - minv)/(Nmagbins)*2.0
 !  write(*,*)'min,maxV,Vcut:', minv,maxv,Vcut
   do i = 1, Nmagbins/2
      binmag(i) = minv + i*width1
   enddo   
   width2 = (maxv - Vcut)/(Nmagbins)*2.0
   do i = Nmagbins/2+1,Nmagbins
      binmag(i) = Vcut + (i- Nmagbins/2)*width2
   enddo
!   do i = 1, Nmagbins
!      write(*,*)'i,magbin:',i,binmag(i)
!   enddo   
   ngood = 0
   nbad = 0
   sigV = 0.0
   sigVI = 0.0
   do i = 1, npt
      if (dist(i) <= binrad(1) ) then
         indxrad = 1
      elseif (dist(i) > binrad(Nradbins-1) ) then
         indxrad = Nradbins
      else
         indxrad = int( (dist(i)- binrad(1) )/widthrad) + 2
!         write(*,*)indxrad,dist(i),binrad(indxrad)
      endif
      if (vv(i) <= vcut ) then
         indxV = int( (vv(i) - minv)/width1 ) + 1
      else
         indxV = int( (vv(i) - vcut)/width2 ) +  Nmagbins/2 + 1
         if(indxV > Nmagbins ) indxV=Nmagbins
      endif
 !          write(*,*) vv(i),binmag(indxV),indxV
      if(recovered(i) ) then
         ngood(indxrad,indxV) = ngood(indxrad,indxV) + 1
!         write(*,*)i,indxrad,indxV,ngood(indxrad,indxV)
         sigV(indxrad,indxV) = sigV(indxrad,indxV) + abs (vv(i) - vout(i) )
         sigVI(indxrad,indxV)=sigVI(indxrad,indxV) + abs (viin(i) - viout(i))
!         indxh(indxrad,indxV) = indxh(indxrad,indxV) + 1
         histVerr(indxrad,indxV,ngood(indxrad,indxV)) = vv(i) - vout(i)
         histVIerr(indxrad,indxV,ngood(indxrad,indxV))= viin(i) - viout(i)
      else
         nbad(indxrad,indxV) = nbad(indxrad,indxV) + 1
      endif
   enddo   !write everything out, including the histograms of the photometric errors
   !for now, write to file numbers 21+
   outfile = 1
   completeness = ngood/(real(ngood + nbad) )
   write(*,*)Nradbins,Nmagbins
   do j = 1,Nmagbins
      do i = 1,Nradbins
         sigV(i,j) = sigV(i,j)/real(ngood(i,j))
         sigVI(i,j)=sigVI(i,j)/real(ngood(i,j))
         write(*,*) binrad(i),binmag(j), ngood(i,j), &
              completeness(i,j),nbad(i,j),sigV(i,j),sigVI(i,j)         
         write(filename1, '(A4,I0.2,A4)') "Verr", outfile,".dat"
         filename1 = trim(filename1)
         write(filename2, '(A5,I0.2,A4)') "VIerr", outfile,".dat"
         filename2 = trim(filename2)
 !        write(*,*)filename1, filename2
         open(unit=20, file=filename1)
         open(unit=21, file=filename2)
         write(20,'("#Npts     Radius    F606W    Completeness")' )
         write(20,'(A1,I5,F10.4,F12.6,F9.6)')"#",  &
              ngood(i,j),binrad(i),binmag(j), completeness(i,j)
         write(21,'("#Npts     Radius    F606W   Completeness")' )
         write(21,'(A1,I5,F10.4,F12.6,F9.6)')"#", &
              ngood(i,j),binrad(i),binmag(j), completeness(i,j)
         write(20,'("# Verr")')
         write(21,'("# VIerr")')
         do k = 1,ngood(i,j)
            write(20,*) histVerr(i,j,k)
            write(21,*) histVIerr(i,j,k)
         enddo
         close(unit=20)
         close(unit=21)
         outfile = outfile + 1
      enddo
      write(*,*)
   enddo   

   return
 end subroutine GetPhotProp
! -----------------------------------------------------------------
 subroutine FindMissing(xin,yin,xout,yout,vin,viin,vout,viout,npt,recovered)   !based upon Anderson, Jay et al. 2008, consider a star recovered if position is within 0.5 pix and
   ! magnidue with 0.75 mag of their input values   
   implicit none
   integer :: npt   
   real, dimension(npt) :: xin,yin,xout, yout, vin,viin,vout,viout
   logical, dimension(npt) :: recovered   
   real, parameter :: pix_miss = 0.5 !max shift between in/out coord
   real, parameter :: mag_miss = 0.75 !max shift between in/out magnitude
   real :: distance,iin,iout   
   integer :: i, j   
   j = 0
   do i = 1,npt
      distance = sqrt( (xin(i) - xout(i) )**2  + ( yin(i) - yout(i) )**2 )
      iin = vin(i) - viin(i)
      iout = vout(i)  - viout(i)
      if ( distance < pix_miss .and. &
            (abs(vin(i) - vout(i) ) < mag_miss ) .and. &
            (abs(iin - iout) < mag_miss ) ) then
         recovered(i) = .true.
         j = j + 1
      else
         recovered(i) = .false.
      endif
   enddo
!   write(*,*)'npt, ngood:', npt, j   
   return
 end subroutine FindMissing 
! ----------------------------------------------------------------- 
 subroutine FindCenter(x,y,dist, xcen,ycen,n)   
   implicit none
   integer :: n
   real, dimension(n) :: x,y,dist
   real :: xcen, ycen   
   real, parameter :: xguess = 3000
   real, parameter :: yguess = 3000
   real, parameter :: radius = 200  !200 pixels from guess center   
   real :: rad2,dist2, sumx,sumy
   integer i,ncen   
   ncen = 0
   xcen = 0.0
   ycen = 0.0
   do i = 1,n
      rad2 = radius**2
      dist2 = (x(i) - xguess)**2 + (y(i) -yguess)**2
      if (dist2 <= rad2 ) then
         ncen = ncen + 1
         xcen = xcen + x(i)
         ycen = ycen + y(i)
      endif
   enddo
   xcen = xcen/ncen
   ycen = ycen/ncen
!   write(*,*)'ncen,x,y:', ncen, xcen, ycen   
   do i = 1,n
      dist(i) = sqrt( (x(i) - xcen)**2 + (y(i) - ycen)**2 )
   enddo 
!  write(*,*)'min/max dist:', minval(dist),maxval(dist)
   return
 end subroutine FindCenter
