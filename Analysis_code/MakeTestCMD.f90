program MakeFakeCMD
  ! program to generate a test CMD +/- 2.5 mag from the subgiant branch
  ! from an isochrone, using the observed location of points and artifical
  ! star testss 
  ! note that right now, the magnitude of the SGB in the artifical stars is hardwired
  ! into the code, set in module define_error
  ! the slope of the power law present day mass function (called IMFslope in the code)
  ! is read in from the command line.  Paust et al. 2010 found -1.23 for M92, where
  ! Saltpeter slope is -2.35, but this is a global mass function.  Figure 2 shows the
  ! HST ACS GC treasury program mass fuction.  Taking points with 0.274 < M < 0.71 one
  ! finds a slope of -1.02, so let's use that value (see MassFuncCut.dat for PDMF data)
  !binary fraction is also read in on the command line.  Milone et al. 2012 characterize
  ! the photometric binary fraction with q> 0.5 (mass ratio) and found for M92 a binary
  ! fraction of 2% and it looks pretty flat as a function of mass (see figures 24, 25 & 31)
  ! code assumes a flat binary ratio betwen q= 0.5 to 1.0

  !for now, maximum number of stars in the error historgams is hardwired into the code
  ! see get_phot_errors
  
!  use real_precision
!  use qsort
  use random
  use define_errors

  implicit none

  integer, parameter :: nmaxisoage = 100 !maximum number of isochrone ages

  real, parameter :: ten = 10.0, negpfour= -0.4, ntwopfive =-2.5
  
  real :: xx,Viso,Iiso
  
  REAL, dimension(:),allocatable :: dist !input distribution of distances
  real, dimension(:),allocatable :: rdist  !randomly generated distance from center
  real, dimension(:),allocatable :: mass   !randomly drawn from powerlaw IMF
  real, dimension(:),allocatable :: VV    ! F606W mag of simulated star
  real, dimension(:),allocatable :: VI,II    ! F814W colour of simulated star
  real, dimension(:),allocatable :: cdfy !culumative distribution function
  real, dimension(:),allocatable :: isomass,isoVI,isoV,isoI
  real :: temp,minmass,maxmass,IMFslope,m1,m2,powmass,age,Vsgb,Isgb, &
       Verr,Ierr, binaryfraction,secondarymass,Vsecondary,CMDage,  &
       VIsecondary,Isecondary,Iprimary,Imag,Vtemp,Itemp,VItemp

  integer :: nstar  !number of fake stars to generate

  
  integer :: nobs !number of oberved stars in the CMD
  integer :: i,j,k,n_args,niso,endf,indxmin,indxmax,intage, &
       indxmcnum,numrad,nummag, &
       fitindxmax,fitindxmin,npt

  character(len=12) ::cIMFslope,cbinfrac,cnstar,cCMDage
  character(len=150) :: filename
  character(len=1) :: c1
  character(len=4) :: chkstr="cmd."    ! extract MC number from isochorne filename
  character(len=5) :: mcnumber,cage
  character(len=14) ::outfile

  n_args = command_argument_count()
  if (n_args /= 5) then
     write(*,*)'Usage: ./a.out  filename MassFunction_power_law_exponent Binary_Fraction Nstars CMDage'
     call get_command_argument(1,filename)
     filename = trim(filename)
     call get_command_argument(2,cIMFslope)
     cIMFslope = trim(cIMFslope)
     read(cIMFslope,*) IMFslope
     write(*,*) "filename:", filename
     write(*,*) "imfslope:", IMFslope
     stop
  end if

  call get_command_argument(1,filename)
  filename = trim(filename)
  call get_command_argument(2,cIMFslope)
  cIMFslope = trim(cIMFslope)
  read(cIMFslope,*) IMFslope
  call get_command_argument(3,cbinfrac)
  cbinfrac = trim(cbinfrac)
  read(cbinfrac,*) binaryfraction
  call get_command_argument(4,cnstar)
  cnstar = trim(cnstar)
  read(cnstar,*) nstar
  call get_command_argument(5,cCMDage)
  cCMDage = trim(cCMDage)
  read(cCMDage,*) CMDage
 
  write(*,*)"MassSlope, BinaryFrac, nstari, CMDage:", IMFslope, binaryfraction,nstar, CMDage

  allocate (rdist(nstar),mass(nstar),VV(nstar),II(nstar),VI(nstar) )

  indxmcnum = index(filename,chkstr)
!  write(*,*)chkstr, indxmcnum
  mcnumber = filename(indxmcnum+4:indxmcnum+9)
!  write(*,*)"number: ", mcnumber
!  write(*,*)'done'
  
  
  call init_ranseed()

  !read in the distrubtion of distances (from cluster center)
  open(unit=10,file='./inputfiles/Distance.dat', status='old')
  read(10,*)c1,nobs
  allocate ( dist(nobs),cdfy(nobs) )
  do i = 1,nobs
     read(10,*) dist(i)
     cdfy(i) = real(i)/nobs
  enddo
  close(10)

  !read in the photometric errors as a function of distance and magnitude
  !see module define_errors for the global variable declarations

  call get_phot_errors

!  write(*,'("# indx    distance    mass")' )
  ! calculate the radial distibution of points
  do i = 1,nstar
     call random_number(xx)
!     write(*,*)'xx:', xx
     do while (xx < cdfy(1) )
        call random_number(xx)
     enddo
     call lininterp(cdfy,dist,xx,rdist(i),nobs)
  enddo

  open(unit=10, file=filename, status='old')
  !MY 04/22 modifed the code to generate sCMD for a given age of isochrone
  do j = 1, nmaxisoage
     read(10,*,iostat=endf)
     if (endf < 0) exit     !end of isochrone file, we are done
     read(10, '(1x,I3,F9.6,F8.4,F7.0)') niso,temp,m1,age
     niso = niso
     if ( abs( age - CMDage ) > 0.01 ) then
        read(10,*)
        do i = 1,niso
           read(10,*)
        enddo
        read(10,*)
        read(10,*)
     else        
 !    write(*,*)'age=',age
        read(10,*)
        allocate (isomass(niso), isoV(niso),isoVI(niso),isoI(niso) )
        do i = 1,niso
           read(10,*)k, isomass(i),isoV(i),isoVI(i)
           isoI(i) = isoV(i) - isoVI(i)
        enddo
    !    if (isomass(niso) < 0.000000001) then
    !       isomass(niso) = isomass(niso - 1)
    !       isoV(niso) = isoV(niso - 1)
    !       isoVI(niso) = isoVI(niso - 1)
    !       isoI(niso) = isoV(i) - isoVI(i)
    !    endif
        read(10,*)
        read(10,*)
   !     read(10,*)
        call findFitPoints(isoV,isoVI,niso,Vsgb,Isgb, indxmin,indxmax,fitindxmin,fitindxmax)
   !     write(*,*)'in main:',indxmin,indxmax,isomass(indxmin),isomass(indxmax), Vsgb

        minmass = isomass(fitindxmin)
        maxmass = isomass(fitindxmax)
   !     write(*,*)'age, min/max mass:',age, minmass,maxmass
        m1 = maxmass**(IMFslope+1.) - minmass**(IMFslope+1.)
        m2 =  minmass**(IMFslope+1.)
        powmass = (1.0/(IMFslope+1.) )
        
        npt = 0
        
        do i = 1, nstar
           call random_number(xx)
           mass(i) = ( m1*xx + m2 )**powmass
           !interpolate in isochrone mass to get V and VI of simulated star
           call lininterp(isomass(fitindxmin:fitindxmax),isoV(fitindxmin:fitindxmax), &
                mass(i),Vtemp,fitindxmax-fitindxmin)
           call lininterp(isomass(fitindxmin:fitindxmax),isoVI(fitindxmin:fitindxmax), &
                mass(i),VItemp,fitindxmax-fitindxmin)
           Itemp = Vtemp - VItemp
   !        write(*,'(I6,2x,F11.4,4F10.5)')i,rdist(i),mass(i),VV(i),VI(i),II(i)
           
    !make a binary if needed
           call random_number(xx)
           if(xx <= binaryfraction) then
              !assume flat secondary mass distrubtion with q=0.5 to 1.0)
              call random_number(xx)
              secondarymass =(0.5 + 0.5*xx )*mass(i)
              !get magnitude and color of binary star
              call lininterp(isomass(1:indxmax),isoV(1:indxmax), &
                    secondarymass,Vsecondary,indxmax)
              call lininterp(isomass(1:indxmax),isoVI(1:indxmax), &
                   secondarymass,VIsecondary,indxmax)
   !            write(*,*)'binary:',secondarymass, Vsecondary,VIsecondary
   !            write(*,*)'primary:', mass(i),VV(i),VI(I)
   !           Iprimary = VV(i) - VI(i)
              !combine the two magnitudes to get one simulated star
              Vtemp = ntwopfive*log10( &
                   ten**(negpfour*Vtemp) + ten**(negpfour*Vsecondary) )
              Isecondary = Vsecondary - VIsecondary
              Itemp = ntwopfive*log10(  &
                   ten**(negpfour*Itemp) + ten**(negpfour*Isecondary) )
              VItemp =Vtemp- Itemp
   !            write(*,*)'Imag:', Iprimary, Isecondary, Imag
   !            write(*,*)'combined:', VV(i),VI(i)
            endif
    !end binary
          
           call findIndxForErrors(Vtemp,rdist(i), Vsgb, 1, numrad,nummag)
           call random_number(xx)
           if (xx < completness(numrad,nummag,1) ) then 
              call getError(numrad,nummag,1, Verr)
              if (Vtemp + Verr >= isoV(fitindxmax) .and.  &
                   Vtemp + Verr <= isoV(fitindxmin) )   then
                 call findIndxForErrors(Itemp,rdist(i), Isgb, 2, numrad,nummag)
                 call random_number(xx)
                 if (xx < completness(numrad,nummag,2) ) then 
                    call getError(numrad,nummag,2, Ierr)
   !                 write(j+20,*) rdist(i),mass(i),Vtemp, &
   !                      Itemp,Vtemp + Verr, Itemp+Ierr
   !                 call lininterp(isoI(indxmin:indxmax), &
   !                      isoV(indxmin:indxmax), &
   !                      Itemp + Ierr,Viso,indxmax-indxmin)
   !                 if (Vtemp +  Verr <= Viso + 0.1 .and.  &
   !                      Vtemp +  Verr >= Viso - 0.1) then
                    call lininterp(isoI(indxmin:indxmax), &
                         isoV(indxmin:indxmax), &
                         Itemp + Ierr,Viso,indxmax-indxmin)
                    call lininterp(isoV(indxmin:indxmax), &
                         isoI(indxmin:indxmax), &
                         Vtemp + Verr,Iiso,indxmax-indxmin)
                    if (Vtemp +  Verr <= Viso + 0.08 .and. &
                         Vtemp +  Verr >= Viso - 0.08 .and. &
                         Itemp + Ierr <= Iiso + 0.08 .and. &
                         Itemp + Ierr >= Iiso - 0.08) then
                       npt = npt + 1
                       VV(npt) = Vtemp +  Verr
                       II(npt) = Itemp + Ierr
                       VI(npt) = VV(npt) - II(npt)
                    end if
   !                 call lininterp(isoI(indxmin:indxmax), &
   !                      isoV(indxmin:indxmax), &
   !                      II(npt),Viso,indxmax-indxmin)
   !                 VIsolist(npt) = Viso
   !                 endif
                 endif
              endif
           endif
        enddo

        call wrtout(vv(1:npt),ii(1:npt),age,IMFslope, &
             binaryfraction,npt,mcnumber)

     
 
        deallocate (isomass,isoV, isoVI,isoI) 
     endif
  enddo
  close(10)

  stop
end program MakeFakeCMD

!---------------------------------------------------------------------
subroutine wrtout(vv,ii,age,IMFslope, binaryfraction, npt,mcnumber)
  implicit none

  integer :: npt
  real, dimension(npt) :: vv,ii
  real :: age, IMFslope, binaryfraction
  character(len=5) :: mcnumber


  integer :: nhx,nhy !number of bins in x and y 
  !working variables 
  real :: minx,maxx,miny,maxy
  integer j,indxx,indxy

  real :: vi(npt),x(npt),y(npt)


  integer :: intage,i
  character(len=5) :: cage
  character(len=14) ::outfile

  intage = nint(age)
  write(cage,"(I0.5)") intage
  cage = trim(cage)
  outfile = "mc"//mcnumber//".a"//cage
  write(*,*)"age: ",cage," outfile: ", outfile
!  open(unit=25,file=outfile,form='unformatted')
  open(unit=25,file=outfile)
!  write(25, '(A29)' ) '# Obs_F606W         Obs_F814W'


!  write(*,*)'widthX,widthV:', dx,dy
  vi = vv - ii
!  x = vi
!  y = vv
  write(25,'("#Total_stars  MassSlope BinaryFraction")' )
  write(25,'(A1,1x,i11,3x,2F9.4)' )  &
       "#", npt, IMFslope, binaryfraction
  write(25,'("#  F606W            F606-F814 " )' )
!  do j=1,nhy
!     do i = 1,nhx
!        write(25,*) xhist(i),yhist(j),hist(i,j)
!     enddo
     !comment out the line below to remove the blank lines
!     write(25,*)
!  enddo
  
  do i=1,npt
     write(25,*) vi(i),vv(i)
  enddo
  close(25)
  return
end subroutine wrtout


!---------------------------------------------------------------------
subroutine  getError(numrad,nummag,ifilter, Verr)
  use define_errors
  implicit none
  
  ! using the artifical star error distribution function for a given radius and magnitude
  ! bin, calculate a random error in F606W (Verr) and F606W-F814W (VIerr)

  integer :: numrad, nummag,ifilter
  real :: Verr

  real :: xx

!  write(*,*)'In GetError:', numrad,nummag
  call random_number(xx)
!  write(*,*)'xx:', xx
  do while (xx < cdfy_photerrors(1,numrad,nummag,ifilter) )
     call random_number(xx)
  enddo
  call lininterp(cdfy_photerrors(:,numrad,nummag,ifilter), &
       histMagerr(:,numrad,nummag,ifilter),xx,Verr,nhist(numrad,nummag,ifilter) )
  
end subroutine getError

!----------------------------------------------------------------------------
subroutine  findIndxForErrors(VV,distance, Vsgb, ifilter, numrad,nummag)
  use define_errors
  implicit none

  real :: VV,distance, Vsgb
  integer :: ifilter
  integer :: numrad, nummag

  integer :: i
  real :: deltamag

  !search through the artifical star radial and magnitude bins to find the
  ! bin of our target star, which is located at distance with VV mag

  if (distance > binrad(Nbinrad-1, ifilter) ) then
     numrad = Nbinrad
!     write(*,*)'HitHigh:', numrad, binrad(Nbinrad,ifilter), distance
  else
     do i = 1, Nbinrad-1
        if (distance <= binrad(i,ifilter) ) then
           numrad = i
            exit
        endif
     enddo
  endif
  
!  write(*,*)'hit: ',  numrad,binrad(numrad,ifilter), distance

  !magnitudes have different zero-points, find the offset from the known
  !subgiant branch magnitudes
!  write(*,*)'ifilter, SGBmag:', ifilter,Vsgb,SGBartStar(ifilter) 
  deltamag = Vsgb - SGBartStar(ifilter)

  !ie to put artfical stars into the isochrone zero point we have
  ! V_artstars_in_iso =V_artstars +  deltaMag

  if (vv > binmag(Nbinmag -1,ifilter )  + deltaMag  ) then
     nummag = Nbinmag
!     write(*,*)'HitHigh:', nummag, binmag(nummag,ifilter)+deltaMag,  &
!           binmag(nummag-1,ifilter)+deltaMag, vv
     return
  endif

  do i = 1,Nbinmag-1
     if (vv <= binmag(i,ifilter) + deltaMag  ) then
        nummag = i
!        write(*,*)'Hit:', nummag, binmag(nummag,ifilter)+deltaMag,  &
!              binmag(nummag-1,ifilter)+deltaMag, vv
        return
     endif
  enddo

  write(*,*)'did not find the magnitude'
  write(*,*) nummag, binmag(i,ifilter)+deltaMag, vv
  stop

end subroutine findIndxForErrors


!----------------------------------------------------------------------------
subroutine findFitPoints(VV,VI,niso,Vsgb,Isgb,indxmin,indxmax, fitindxmin,fitindxmax)
  implicit none

  !find the location of the subgiant branch, and then the indexes to the
  ! fitpoints, which are +/- 2.5 mag from the subgiant branch, plus the index
  ! of the points which are +/- 3.25 mag from the subgiant branch, which are
  ! used to generate the initial simulated fake CMD

  integer :: niso, indxmin, indxmax,fitindxmin,fitindxmax
  real, dimension(niso) :: VV, VI

  real ::  Vsgb, Isgb

  real, parameter  :: sgbDeltaColour = 0.05  !difference in color between MSTO & SGB
  real, parameter :: deltaMag = 2.0   !fit region is defined to be within +/-deltaMag of the SGB
  integer :: i,ito,nmax,locate
  real :: minVI,VIsgb,hitmag

  minVI = 50.0

  !define the main sequence turn-off as the bluest point
  do i = 1,niso
     if(vi(i) < minVI ) then
        minVI = vi(i)
        ito = i
     endif
  enddo

  !find subgiant branch magnitude, taking advantage of the fact we
  ! know it is close the the turnoff, so don't need to search the entire array
  VIsgb = vi(ito) + sgbDeltaColour
  nmax = min(100,niso-ito)
!  write(*,*)'ito,nmax,npt:', ito, nmax
  call lininterp(VI(ito:nmax+ito), VV(ito:nmax+ito), VIsgb, Vsgb, nmax )
  Isgb = Vsgb - VIsgb
!  write(*,*)"TO,SGB: ", minVI,VIsgb, Vsgb,Isgb

  !find the location of the faint end of the fit region
  hitmag = Vsgb + deltaMag
  fitindxmin = locate (VV(1:ito), hitmag, ito )
  !find the location of the bright end of the fit region
  hitmag = Vsgb - deltaMag
  fitindxmax = locate(VV(ito:niso),hitmag,niso-ito)
  fitindxmax = ito + fitindxmax
!  write(*,*)indxmin, Vsgb + deltaMag, VV(indxmin), indxmax, Vsgb - deltaMag,VV(indxmax)
  !  write(*,*)indxmin,indxmax,VV(indxmin), VV(indxmax)

  !find the location of the point 1 mag fainter than the fit region
  hitmag = Vsgb + deltaMag + 1
  !hitmag = Vsgb + deltaMag
  indxmin = locate (VV(1:fitindxmin), hitmag, fitindxmin ) - 1
  !find the location of the point 1 mag brighter than the fit region
  hitmag = Vsgb - deltaMag - 1
  !hitmag = Vsgb - deltaMag
  indxmax = locate(VV(fitindxmax:niso),hitmag,niso-fitindxmax) + 1
  indxmax = fitindxmax + indxmax

!   write(*,*)indxmin, Vsgb + deltaMag+1.0, VV(indxmin), indxmax, Vsgb - deltaMag -1.0,VV(indxmax)

  
  return
end subroutine findFitPoints


!----------------------------------------------------------------------------
  subroutine get_phot_errors

     use define_errors

     implicit none

  !read in the binned artificial star photometric errors.  Note that these are
  !kept in many different files.... ./inputfiles/Verr??s.dat .and ./inputfiles/VIerr??s.dat
  !files are already sorted.  the input files ./inputfiles/VerrBopundary.dat and
  ! IerrBoundary.dat give us
  ! the bin boundaries and number of points in each bin 

 
  integer :: i,j,k,iopen,ifilter
 
  real :: deltamag
  real :: tr,tv,tcomp
  integer :: numt,numx

  character(len=29) :: filename
  character(len=25) :: filename1
  character(len=17),parameter,dimension(nfilter) ::  &
       startfile =["./inputfiles/Verr", "./inputfiles/Ierr"]


  write(filename,'(A17,A12)') startfile(1), "Boundary.dat"
  open(unit=30, file=filename)
  read(30,*)
  read(30,*) Nbinrad,Nbinmag
  read(30,*)
!  rewind(30)

  nmaxstars = 750
  allocate ( nhist(Nbinrad,Nbinmag,nfilter), completness (Nbinrad,Nbinmag,nfilter) , &
       binrad(Nbinrad,nfilter), binmag(Nbinmag,nfilter) )
  allocate( cdfy_photerrors(nmaxstars,Nbinrad,Nbinmag,nfilter) )
  allocate ( histMagerr(nmaxstars,Nbinrad,Nbinmag, nfilter) )

  
  do ifilter = 1,2
     if (ifilter == 2) then 
        write(filename,'(A17,A12)') startfile(ifilter), "Boundary.dat"
!        write(*,*) "Filename: ", filename
        open(unit=30, file=filename)
        read(30,*)
        read(30,*) Nbinrad,Nbinmag
     !  write(*,*) 'Nbins: ', Nbinrad,Nbinmag
        read(30,*)
     endif
 
     do j=1,Nbinmag
        do i = 1,Nbinrad
!        write(*,*)'before read:', i,j
           read(30,*) binrad(i,ifilter),binmag(j,ifilter), nhist(i,j,ifilter), &
                completness(i,j,ifilter)
           !        write(*,*)i,j,binrad(i),nhist(i,j)
        enddo
     enddo
     close(30)

!     write(*,*)'radbin 1, magbin 4:', binrad(1,ifilter),binmag(4,ifilter),nhist(1,4,ifilter)
  
!     write(*,*)'radbin 2, magbin 3:', binrad(2,ifilter),binmag(3,ifilter), nhist(2,3,ifilter)

     nmaxstars = maxval(nhist)
!     write(*,*) 'Maximum number of stars in a bin: ', nmaxstars
     do k = 1,Nbinmag
        do j = 1,Nbinrad
           iopen = (k-1)*10 + j
!        write(*,*) 'iopen:', iopen
           write(filename1, '(A17,I0.3,A5)') startfile(ifilter), iopen,"s.dat"
!           filename1 = trim(filename1)

!           write(*,*)"reading file: ", filename1
           open(unit=11, file=filename1)
           read(11,*)
           read(11,'(3X,I3, 1x, F9.4, 2x, F10.6, 1x, F8.6)') numt,tr,tv,tcomp
           read(11,*)
           if ( numt /= nhist(j,k,ifilter)  ) then
              write(*,*)'something odd: ', numt, nhist(j,k,ifilter),j,k,ifilter
              stop
           endif
           if(   abs( binrad(j,ifilter) - tr)> 0.0001  &
                .or. abs(binmag(k,ifilter) - tv) > 0.001 ) then
              write(*,*)'something odd: ' 
              write(*,*)binrad(j,ifilter),tr,binrad(j,ifilter) - tr, j       
              write(*,*)binmag(k,ifilter), tv,binmag(k,ifilter) - tv,k
              stop
           endif
!        write(*,*)'j,k,nhist:', j,k,nhist(j,k,ifilter)
           do i = 1, nhist(j,k,ifilter)
              cdfy_photerrors(i,j,k,ifilter) = real(i)/nhist(j,k,ifilter)
              read(11,*)histMagerr(i,j,k,ifilter)
           enddo
           close(11)
        enddo
     enddo
  enddo

end subroutine get_phot_errors

  

!----------------------------------------------------------------------------
subroutine lininterp(x,y,xval,yval,npts)
  implicit none

  !subroutine for linear interpolation in x,y, assuming they are sorted
  integer :: npts
  real :: x(npts),y(npts)
  real  :: xval
  real :: yval

  integer :: i,locate
  logical :: ascnd


!  write(*,*)'Linterp:',x(1), xval
  ascnd = (x(npts) >= x(1))
  if (ascnd) then
     if (xval <= x(1) ) then
        yval = y(1)
!     yval = (y(1) - 0.0_R8)/x(1)*xval +0.0_R8
!     write(*,'(4F12.5)')xval,x(1),y(1),yval
        return
     elseif (xval >= x(npts) ) then
        yval = y(npts)
        return
     else
        i=locate(x,xval,npts)
        yval = (y(i+1)-y(i))/(x(i+1)-x(i)) * (xval-x(i)) + y(i) 
!     write(*,'(i6,6F12.5)')i,xval,x(i+1),x(i),y(i+1),y(i),yval
     end if
  else
     if (xval <= x(npts) ) then
        yval = y(npts)
!     yval = (y(1) - 0.0_R8)/x(1)*xval +0.0_R8
!     write(*,'(4F12.5)')xval,x(1),y(1),yval
        return
     elseif (xval >= x(1) ) then
        yval = y(1)
        return
     else
        i=locate(x,xval,npts)
        yval = (y(i)-y(i-1))/(x(i)-x(i-1)) * (xval-x(i-1)) + y(i-1) 
!     write(*,'(i6,6F12.5)')i,xval,x(i+1),x(i),y(i+1),y(i),yval
     end if
  end if

  return

end subroutine lininterp

!-----------------------------------------------------------------------
integer function locate(xx,x,n)
    ! Locate a value in a sorted array

  implicit none
  integer :: n
  real, dimension(n), intent(in) :: xx
  real, intent(in) :: x
  integer :: jl,jm,ju
  logical :: ascnd

    ascnd = (xx(n) >= xx(1))
    jl=0
    ju=n+1
    do
       if (ju-jl <= 1) exit
       jm=(ju+jl)/2
       if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
       else
          ju=jm
       end if
    end do

    if (x == xx(1)) then
       locate = 1
    else if (x == xx(n)) then
       locate = n-1
    else if(ascnd.and. (x > xx(n) .or. x < xx(1))) then
       locate = -1
    else if(.not.ascnd.and. (x < xx(n) .or. x > xx(1))) then
       locate = -1
    else
       locate = jl
    end if
  end function locate

!---------------------------------------------------------------------


  
