program PhotoPropI
  use define_mags
  
  implicit none
  ! program to read in ACS treasury artificial star file and
  ! characterize the completenss, and photometric uncerainty
  ! as a function of distance from center of the cluster, 
  ! F606W and F814 magnitudes for stars around the main sequence
  ! turn-off (defined by Vmin, Vmax, VImin and VImax).
  ! I use VV for F606W and II for F814W filters
  ! assumses a guess center of cluster is at  3000,3000
  ! go +/- 2.0 mag around the SGB point 
   !  Variable Declarations
  real, parameter :: Vmax = -11.7
  real, parameter :: Vmin = -15.683
  real, parameter :: IImax = -11.568
  real, parameter :: IImin = -15.737
  real, parameter,dimension(2) :: short_long = [-14.08, -14.14]
  !magnitude between short and long exposures, first F606W and then F814W!
!  real, parameter :: Ishort_long = -13.9
  integer, parameter :: Nmagbins = 8 !totalnumber of mag bins (must be even)
  integer, parameter :: Nradbins= 10 !number of radial bins
  integer :: Nstars,nfitI,nfitV,nfit,i,j,id,k
  real, dimension(:), allocatable :: vti,iiti, vto,iito
  real :: ii_in,ii_out,xcenter,ycenter
  real :: i1,i2,i3,i4
  logical :: L_AS
 
  character(len=40) :: filename
  character(len=10) :: cNstars
  character(len=10) :: cL_AS
  integer :: n_args
  
  !read in Sample Size & filename from command line
  n_args = command_argument_count()
  if (n_args /= 3) then
     write(*,*)'Usage: ./a.out num_stars filename L_AS'
     stop
  end if
  call get_command_argument(1,cNstars)
  cNstars = trim(cNstars)
  read(cNstars,*) Nstars
  call get_command_argument(2,filename)
  filename = trim(filename)
  call get_command_argument(3,cL_AS)
  cL_AS = trim(cL_AS)
  if (cL_AS == 'True') then
      L_AS = .True.
  else
      L_AS = .False.
  endif
  allocate (xti(Nstars),yti(Nstars),vti(Nstars), &
       iiti(Nstars), xto(Nstars),yto(Nstars), vto(Nstars), &
       iito(Nstars) )
  !Read in data and count number of stars in the MSTO fit region
  nfitV = 0
  nfitI = 0
  open (unit=4, file = filename, status='old')
  read(4,*)    !skip header line
  do i = 1,Nstars
     read(4,*) xti(i),yti(i),vti(i),iiti(i),&
          xto(i),yto(i),vto(i),iito(i)
     if(  vti(i) <= Vmax .and. vti(i) >= Vmin ) nfitV=nfitV+1.
     if(  iiti(i) <= IImax .and. iiti(i) >= IImin ) nfitI=nfitI+1.
   enddo
   close(4)
   write(*,*)'Number of stars in MSTO region:', nfitV,nfitI
   ! treat V and I separately, so do this twice
   do k = 1, 2
      if(k == 1 ) then
         call instars(vti,vto,Vmax,Vmin,nstars,nfitV)
!         do i = 1, nfitV
!            write(*,*)i,magin(i),magout(i)
!         enddo
      else
         call instars(iiti,iito,IImax,IImin,nstars,nfitI)
!         do i = 1, nfitI
!            write(*,*)i,magin(i),magout(i)
!         enddo
      endif
      !find out which stars are bad/missing from output file
      call FindMissing
      !find the center of the input stars, and distance from center  
      call FindCenter(xcenter,ycenter)
      write(*,*)'Main xcen,ycen', xcenter, ycenter
! find completenss and magnitude error as a function of distance from the 
! cluster center and magnitude.  The routine  sets bin boundaries, and
! outputs the results.
      call GetPhotProp(short_long(k), Nradbins,Nmagbins,k, L_AS)
         
      deallocate (xin,yin,magin,xout,yout, magout, &
           dist, recovered )
   enddo
 END program PhotoPropI
 !------------------------------------------------------------------
 subroutine instars(magi,mago,magmax,magmin,nstars,Nfit)
 
   use define_mags
   implicit none
   integer :: Nfit,nstars
   real, dimension(nstars) :: magi,mago
   real :: magmax, magmin
   integer :: i,j
   allocate(xin(Nfit),yin(Nfit),magin(Nfit), &
        xout(Nfit),yout(Nfit), magout(Nfit), &
        dist(Nfit), recovered(Nfit) )
   write(*,*) 'instars N:', Nfit,nstars
   j = 0 
   do i = 1,nstars
      if(  magi(i)  <= magmax .and. magi(i) >= magmin ) then 
         j = j + 1
         xin(j) = xti(i)
         yin(j) = yti(i)
         magin(j) = magi(i)
         xout(j) = xto(i)
         yout(j) = yto(i)
         magout(j) = mago(i)
      endif
   enddo
   return
 end subroutine instars
 ! -----------------------------------------------------------------
 subroutine GetPhotProp(vcut, Nradbins,Nmagbins,kin,L_AS)
   use define_mags
   implicit none
   real,intent(IN) :: Vcut
   integer, intent(IN) :: kin
   integer :: Nradbins,Nmagbins
   real, dimension(Nradbins,Nmagbins) ::  completeness
   real, dimension(Nradbins) :: binrad
   real, dimension(Nmagbins) :: binmag
   integer :: i,j,indxrad,indxV , outfile,k, npt
   real :: maxdist,width1,width2,minv,maxv,widthrad
   integer, dimension(Nradbins,Nmagbins) :: ngood,nbad
   real, dimension(Nradbins,Nmagbins) :: sigV
   real, dimension(Nradbins,Nmagbins,20000) :: histVerr
   character(len=30) :: filename1
   character(len=16) :: filename2
   character(len=40),parameter,dimension(2) :: startfile =["Verr", "Ierr"]
   logical :: L_AS
   npt = size(dist)
   write(*,*)'in GetPhotProp npt:', npt
   maxdist = maxval(dist)
!   write(*,*)maxdist
   if (Nradbins .ne. 2 ) then
      widthrad = maxdist/Nradbins
      binrad(1) = widthrad*1.5
      binrad(Nradbins) = maxdist
      maxdist = maxdist- 3.0*widthrad
      widthrad = maxdist/Nradbins
      do i = 2,Nradbins-1
         binrad(i) = binrad(i-1) + widthrad
      enddo
   else
      widthrad = maxdist/Nradbins
      binrad(1) = widthrad
      binrad(2) = maxdist
   endif
!   do i = 1, Nradbins
!      write(*,*)'i,binrad:',i,binrad(i)
!   enddo
   minv = minval(magin)
   maxv = maxval(magin)
   !when L_AS is False, we want to use only 1 bin for short exposure
   !as there is very few bright stars
   if (L_AS) then
      if (mod(Nmagbins,2) .ne. 0 ) then
         write(*,*)'Nmagbins must be an even number'
         stop
      endif
      width1 = (Vcut - minv)/(Nmagbins)*2.0
   !  write(*,*)'min,maxV,Vcut:', minv,maxv,Vcut
      do i = 1, Nmagbins/2
         binmag(i) = minv + i*width1
      enddo
      width2 = (maxv - Vcut)/(Nmagbins)*2.0
      do i = Nmagbins/2+1,Nmagbins
         binmag(i) = Vcut + (i- Nmagbins/2)*width2
      enddo
   else
      width1 = (Vcut - minv)
      binmag(1) = minv
      width2 = (maxv - Vcut)/(Nmagbins - 1)
      do i = 2, Nmagbins
         binmag(i) = binmag(1) + (i-1)*width2
      enddo
   endif
!   do i = 1, Nmagbins
!      write(*,*)'i,magbin:',i,binmag(i)
!   enddo
   ngood = 0
   nbad = 0
   sigV = 0.0
   do i = 1, npt
      if (dist(i) <= binrad(1) ) then
         indxrad = 1
      elseif (dist(i) > binrad(Nradbins-1) ) then
         indxrad = Nradbins
      else
         indxrad = int( (dist(i)- binrad(1) )/widthrad) + 2
!         write(*,*)indxrad,dist(i),binrad(indxrad)
      endif
      if (magin(i) <= vcut ) then
         if (L_AS) then
            indxV = int( (magin(i) - minv)/width1 ) + 1
         else
            indxV = 1
         endif
      else
         if (L_AS) then
            indxV = int( (magin(i) - vcut)/width2 ) +  Nmagbins/2 + 1
            if(indxV > Nmagbins ) indxV=Nmagbins 
         else
            indxV = int( (magin(i) - vcut)/width2 ) +  2
            if(indxV > Nmagbins ) indxV=Nmagbins
         endif
      endif
!           write(*,*) magin(i),binmag(indxV),indxV
      if(recovered(i) ) then
         ngood(indxrad,indxV) = ngood(indxrad,indxV) + 1
!         write(*,*)i,indxrad,indxV,ngood(indxrad,indxV)
         sigV(indxrad,indxV) = sigV(indxrad,indxV) + &
              abs (magin(i) - magout(i) )
!         indxh(indxrad,indxV) = indxh(indxrad,indxV) + 1
         histVerr(indxrad,indxV,ngood(indxrad,indxV)) =  &
              magin(i) - magout(i)
      else
         nbad(indxrad,indxV) = nbad(indxrad,indxV) + 1
      endif
   enddo
!write everything out, including the histograms of the photometric errors
   !for now, write to file numbers 21+
   outfile = 1
   write(filename2,'(A4,A12)')  startfile(kin), "Boundary.dat"
   write(*,*) "Filename:", filename2
   open(unit=30, file=filename2)
   completeness = ngood/(real(ngood + nbad) )
   write(30,'("#      Nradbin      Nmagbin")' ) 
   write(30,*)Nradbins,Nmagbins
   write(30,'("# outboundRad     minMagBound            Nstar  Completenss")')
   do j = 1,Nmagbins
      do i = 1,Nradbins
         sigV(i,j) = sigV(i,j)/real(ngood(i,j))
         write(30,*) binrad(i),binmag(j), ngood(i,j), &
              completeness(i,j)
         
         write(filename1, '(A4,I0.3,A4)') startfile(kin), outfile,".dat"
         filename1 = trim(filename1)
 !        write(*,*)filename1
         open(unit=20, file=filename1)
         write(20,'("#Npts     Radius    Mag      Completeness")' )
         write(20,'(A1,I5,F10.4,F12.6,F9.6)')"#",  &
              ngood(i,j),binrad(i),binmag(j), completeness(i,j)
         write(20,'(A1, 1x,A4)') "#", startfile(kin)
         do k = 1,ngood(i,j) 
            write(20,*) histVerr(i,j,k)
         enddo
         close(unit=20)
         outfile = outfile + 1
      enddo
!      write(*,*)
   enddo
   close(unit=30) 
   return
 end subroutine GetPhotProp
! -----------------------------------------------------------------
 subroutine FindMissing
   use define_mags
   !based upon Anderson, Jay et al. 2008, consider a star recovered if position is within 0.5 pix and
   ! magnidue with 0.75 mag of their input values
   implicit none
   
   real, parameter :: pix_miss = 0.5 !max shift between in/out coord
   real, parameter :: mag_miss = 0.75 !max shift between in/out magnitude
   
   integer :: i, j, npt
   real :: distance 
   npt=size(xin)
   write(*,*)'in FindMissing, npt=', npt
   j = 0 
   do i = 1,npt
      distance = sqrt( (xin(i) - xout(i) )**2  + ( yin(i) - yout(i) )**2 )
      if ( distance < pix_miss .and. &
            (abs(magin(i) - magout(i) ) < mag_miss )  ) then
         recovered(i) = .true.
         j = j + 1
      else
         recovered(i) = .false.
      endif
   enddo
   write(*,*)'npt, ngood:', npt, j
   return
 end subroutine FindMissing
 ! -----------------------------------------------------------------
 subroutine FindCenter(xcen,ycen)
   use define_mags
   
   implicit none
   
   real :: xcen, ycen
   real, parameter :: xguess = 3000
   real, parameter :: yguess = 3000
   real, parameter :: radius = 200  !200 pixels from guess center
   real :: rad2,dist2, sumx,sumy
   integer i,ncen, n
   n = size(xin)
   ncen = 0
   xcen = 0.0
   ycen = 0.0
   do i = 1,n
      rad2 = radius**2
      dist2 = (xin(i) - xguess)**2 + (yin(i) -yguess)**2
      if (dist2 <= rad2 ) then
         ncen = ncen + 1
         xcen = xcen + xin(i)
         ycen = ycen + yin(i)
      endif
   enddo
   xcen = xcen/ncen
   ycen = ycen/ncen
   write(*,*)'ncen,x,y:', ncen, xcen, ycen
   do i = 1,n
      dist(i) = sqrt( (xin(i) - xcen)**2 + (yin(i) - ycen)**2 )
   enddo
   write(*,*)'min/max dist:', minval(dist),maxval(dist)
   return
 end subroutine FindCenter