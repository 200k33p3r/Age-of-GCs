
program test_iso
  implicit none

  ! put isochrones into sigma space and 
  ! find best fit compared to calibrating stars
  ! uses the fact that input isochrone is ordered from
  ! faintest to brightest


   !  Variable Declarations
  integer, parameter :: nstar = 2   !number of calibrating stars
  integer, parameter :: nmaxiso = 100 !maxnumber of ages in isochrone 

  real :: calF606(nstar), calF68(nstar),sigF606(nstar),sigF68(nstar)

  real, dimension (:), allocatable :: f606w,f68,sigmaDist, xs, ys
  real, dimension (:), allocatable :: xx,yy

  integer :: i,j,npts,jminarr(1),nall,ieof,ii, &
       jmin,nage
  real :: x1,x2, slope,intcpt,mindist,sumdist, &
       isoage,pslope, minchi2,chi2(nmaxiso)
  character(len=256) :: filename

  character(len=11) :: chkstr="feh230isol."    ! use to extract MC number from isochorne filename
  character(len=5) :: mcnumber
  character(len=21) :: varfile
  real, dimension(20) :: varnums
  integer :: varcolour
  
  integer :: n_args, indxmcnum
 
  data calF606 / 5.78667, 6.04062 /!absolute F606W mag of calibraters
  data calF68  / 0.566, 0.601 / !F606W-F814W color of calibraters
  data sigF606 / 0.00258611, 0.00373387 / !uncertainty in abs. mag
  data sigF68 / 0.0023, 0.0054 / !uncertainty in color of calibrators
  
!  read in file name  from command line
  n_args = command_argument_count()
  if (n_args /= 1) then
     write(*,*)'Usage: ./iso_sigma filename'
     stop
  end if

  call get_command_argument(1,filename)
  filename = trim(filename)

!figure out the MC index number and find the corrosponding var file and read in the variables
!  write(*,*) filename
  indxmcnum = index(filename,chkstr)
!  write(*,*)chkstr, indxmcnum
  mcnumber = filename(indxmcnum+11:indxmcnum+16)
!  write(*,*)"number: ", mcnumber

  varfile = "./var/varfeh230."//mcnumber
!  write(*,*) "varfile: ", varfile

  open(unit=25, file=varfile)
  do i = 1,20
     read(25,*)varnums(i)
!     write(*,*)varnums(i)
  enddo
  read(25,*)varcolour
!  write(*,*)varcolour
  close(25)

 !read in isochrone of specified age

  nage = 0 
  open (unit=4, file = filename, status='old')
  do 
     read(4,*,iostat=ieof)    !skip header lines
     if (ieof < 0 ) exit
     read(4,'(1x,I3,17X,F6.0)' ) nall,isoage
!     write(*,*) nall,isoage
     nage = nage + 1
     read(4,*)
     allocate (xx(nall), yy(nall) )
     do i = 1, nall
        read(4,'(107X,F6.3,1x,F6.3)' ) yy(i),xx(i)
           !           write(*,*)f606w(i),f68(i)
     enddo
     read(4,*)
     read(4,*)
  !trim isofile to the main sequence
     npts = 0
     do i = 1,nall
 !    write(*,*)'yy:',yy(i)
        if (yy(i) < 8.0 .and. yy(i) > 5.0 ) npts = npts + 1
        if (yy(i) < 5.0) exit
     enddo

!        write(*,*)'age,npts = ',isoage,npts

     allocate (f606w(npts), f68(npts), sigmaDist(npts), &
             xs(npts), ys(npts) )

     i = 0
     do j=1,nall
           !     write(*,*)'here:',yy(j)
        if (yy(j) < 8.0 .and. yy(j) > 5.0)  then
           i = i + 1
           f606w(i) = yy(j)
           f68(i) = xx(j)
!           if(abs(isoage-16500) < 0.01) write(*,*)j,f606w(i),f68(i)
        endif
        if (yy(j) < 5.0) exit
     enddo

     sumdist = 0.0
     do i = 1,nstar
     !put into sigma space
        call sigma(f606w,f68,calF606(i),calF68(i),sigF606(i), &
                sigF68(i),xs,ys,sigmaDist,npts)
        
 !near the sigma space origin, approximate by a straight line and
 ! find slope and intercept.
        jminarr = minloc(sigmaDist)
        jmin = jminarr(1)
!           write(*,*)'age, jmin:',isoage,jmin,jminarr(1),npts
        if (jmin < 2 .or. jmin > npts - 1 ) then
           write(*,*) trim(filename) ,isoage, "bad Jmin", jmin,npts
           stop
        endif
        slope = (ys(jmin + 1) - ys(jmin - 1 ) ) / &
             (xs(jmin + 1) - xs(jmin - 1) )
        intcpt = ys(jmin + 1) - slope*xs(jmin + 1 )
 !now that we have an equation for our the isoline in sigma space
 !geometric formulae for the minimum distance from the origin
 ! to a straight line
        mindist = abs( intcpt) / sqrt(slope*slope + 1.0)
!           write(*,*)isoage,i,mindist
        sumdist = sumdist + mindist
     enddo
!     write(*,*)isoage,sumdist,MCrunNumber,filename
     chi2(nage) = sumdist
  
     deallocate (f606w,f68,sigmaDist,xs,ys)
     deallocate (xx,yy)

  enddo
  !  write(*,*)(chi2(i),i=1,nage)
  minchi2 = minval(chi2(1:nage))
  write(*,'(F10.3,1x,20F9.6,I2,1x,A5)')minchi2,(varnums(i), i=1,20), varcolour, mcnumber

  close(4)
  
END program test_iso

subroutine sigma(yy,xx,y0,x0,sigy0,sigx0,xsig,ysig,sigDist,npts)

  implicit none
  real,dimension(npts) :: yy,xx,xsig,ysig,sigDist
  real :: y0,x0,sigy0,sigx0
  integer :: npts

  !convert to sigma space
  ysig = (y0 - yy)/sigy0
  xsig = (x0 - xx)/sigx0
  sigDist = sqrt(xsig**2 + ysig**2)
 
!  do i = 1,npts
!      ysig(i) = (y0 - yy(i) )/sigy0
!      xsig(i) = (x0 - xx(i) )/sigx0
!      sigDist(i) = sqrt(xsig(i)**2 + ysig(i)**2) 
!     write(*,*)xsig(i),ysig(i),sigDist(i)
!  enddo

end subroutine sigma
