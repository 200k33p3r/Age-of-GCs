      program lowmass
c this subroutine reads in an exsiting isochrone file, and replaces
c the lines below a specified mass by single points from input 
c track files.  Used for lower main sequence MC runs
C     CODE DOES NOT DOE THE MASS FUNCTION, SO JUST SET THESE TO 0.0
C    Isochrone ages in the input file must be ordered low to high  

      implicit none
      integer maxiso,maxpt,maxtrk,nbad,endf,maxptrk
      parameter (maxtrk = 20)     ! maximum number of input low mass tracks
      parameter (maxiso = 90)     !maximum number of imput isochrone ages
      parameter (maxpt = 800)     !maximum number of points on a given isochrone
      parameter (maxptrk = 3000)  !maximum number of points on a given track 
      character*80 hd1
      character*128 hd2
      character*1 cext
      integer npts,npts2,itrk,KTTAU,ntrk,niso,numages
      real*8 minmass,maxage
      real*8 isoage(maxiso),xx,feht,Yprim,alphafe,cmixlen,FGRY,FGRZ,
     $      oversh,sstandard(7),alexcoef,opalcoef2,talphacoef,
     $     plascoef,cocoef,ageneg,teffl(maxtrk,maxiso),
     $     gl(maxtrk,maxiso),bl(maxtrk,maxiso),rl(maxtrk,maxiso)
      real*8 mixlen,ovsh,age,yy,zz,zeff,feh,afe,smass(maxtrk),xenv,zenv,
     $     afe2,cmixl
      real*8 age1(maxptrk),tlog(maxptrk),glog(maxptrk),blog(maxptrk),
     $     rlog(maxptrk),ycore(maxptrk)
      real*8 mass(maxpt),logg(maxpt),logteff(maxpt),logl(maxpt),
     $     num0(maxpt),num235(maxpt),num400(maxpt),vv(maxpt),
     $     vi(maxpt),f606(maxpt),f606f814(maxpt)
      real*8 solbol,filter(13),colr(maxtrk,maxiso,8),vvv(maxtrk,maxiso),
     $        boll(maxtrk,maxiso)
      integer i,j,k,nhit,eep,itable,ii,ntrkpt
      logical LPHXCOL,LVCSCOL,lbad(maxtrk)

      character*10 :: cniso,cntrk,cminmass
      integer n_args
      data itrk/20/
c set some defaults


!     get number of input tracks, number of isochrone ages and
! the minimum mass which should NOT be replaced by the low mass tracks
      n_args = command_argument_count()
      if (n_args /= 3) then
         write(*,*)'./lowmass2 num_inpt_trks, num_iso_ages, minmass'
         stop
      end if
      call get_command_argument(1,cntrk)
      cntrk = trim(cntrk)
      read(cntrk,*) ntrk
      call get_command_argument(2,cniso)
      cniso = trim(cniso)
      read(cniso,*) niso
      call get_command_argument(3,cminmass)
      cminmass = trim(cminmass)
      read(cminmass,*) minmass


!      write(*,*)'ntrk,numages=',ntrk,niso
!      write(*,*)'minmass:', minmass

c******* MC work, read in the color table from the var file.

c     write(*,*)'PHX,VC:', LPHXCOL,LVCSCOL
c      if ( (lphxcol .and. lvcscol) .or.
c     $     (.not.lphxcol .and. .not. lvcscol) ) then
c         write(*,*)'CODE STOPPED:COLOR TABLE SELECTION:',lphxcol,lvcscol 
c         stop
c      endif

      open(unit=75,status='old')
      do i = 1,20
         read(75,*)
      enddo
      read(75,337)itable
 337  format(i2)
      close(75)
      lphxcol=.false.
      lvcscol=.false.
!      write(*,*)'itable:',itable
      if(itable.eq.0) then
         LPHXCOL=.TRUE.
         write(*,*)'using phoenix'
      elseif(itable.eq.1) then
         LVCSCOL=.TRUE.
         write(*,*)'using V&C'
      else
         write(*,*)'incorrect colour table in var file'
         stop
      endif
c***********end of MC changes
      

c     get the ages of the isochrones in the input file

      open(unit=11,status='old')
      do i=1,niso
         read(11,200)hd1
 200     format(a80)
         read(11,210)npts,mixlen,ovsh,isoage(i),yy,zz,zeff,feht,alphafe
 210     format(1x,I3,F9.6,F8.4,F7.0,F8.4,2E12.4,2F6.2)
c         write(*,*)isoage(i)
c         write(*,*)npts
         read(11,230)hd2
 230     format(a128)
         do j=1,npts
            read(11,250)
         enddo
         read(11,*)
         read(11,*)
      enddo
      rewind(11)
      close(11)
!      write(*,*)'read in isochrone ages'
      do i = 2,niso
         if(isoage(i-1) > isoage(i) ) then
            write(*,*)'CODE STOPPED: isochrone ages not increasing'
            write(*,*)(isoage(j),j=1,niso)
            stop
         endif
      enddo
      maxage = isoage(niso)
      write(*,*)'feht,alphafe,maxage:',feht,alphafe,maxage

c initialize the colour tables      
      SOLBOL=4.75D0
      IF(LPHXCOL) THEN 
         itable = 4
         call color_init(itable,feht,alphafe)
         itable = 5
         call color_init(itable,feht,alphafe)
      endif
      IF(LVCSCOL) THEN
         itable = 8
         call color_init(itable,feht,alphafe)      
      endif

      nbad=0
      do i=1,ntrk
         itrk=itrk+1
!ensure that track file includes all requeste ages
         call trkcheck(alphafe,maxage, itrk, nbad,lbad(i) )
         if (lbad(i) ) then
            goto 190        !skip track and go to end of  loop
         endif

c ****** MC changes start
c         read(itrk,*)
c         read(itrk,*)
c         read(itrk,*)
c ****** end MC changes

C read in the standard track header line
         READ(itrk,55)SMASS(i), XENV,ZENV,AFE2,CMIXL
 55      FORMAT(4X,F5.3,3x,E10.4,3x,E10.4,18x,F5.2,6x,F6.4)
         READ(itrk,*)
        
         read(itrk,*)age1(1),tlog(1),glog(1),blog(1),rlog(1),ycore(1)
c        write(*,*)smass(i),xenv,zenv,afe2,cmixl
c        write(*,*)age1,tlog,glog,blog,rlog,ycore1
c  skip part of pre-main sequence
         do while (age1(1).le.1.5D8)
            read(itrk,*)age1(1),tlog(1),glog(1),blog(1),rlog(1),ycore(1)
         enddo

         ageneg=age1(1)/1.0D6

         ii = 2
         do j = 2,maxptrk
            read(itrk,*,iostat=endf)age1(j),tlog(j),glog(j),blog(j),
     $           rlog(j),ycore(j)
            if(endf < 0 ) exit
            age1(j)=age1(j)/1.0D6 - ageneg
            ii = ii + 1
         enddo
         ntrkpt = ii
c         write(*,*) 'ntrkpt=',ntrkpt

c     go through a do a linear interpolation as we hit each age
         ii = 2
         k = 1
c         write(*,*)'ii,age(ii-1):', ii, age1(ii-1)/1.0D6,isoage(niso)
         do while (ii <= ntrkpt .and. k <= niso )
            if(age1(ii) > isoage(k) ) then
c               write(*,*)'hitage,', k,isoage(k),age1(ii)
               call linint(isoage(k),age1(ii-1),age1(ii),tlog(ii-1),
     $              tlog(ii),teffl(i,k) )
               call linint(isoage(k),age1(ii-1),age1(ii),glog(ii-1),
     $              glog(ii),gl(i,k) )
               call linint(isoage(k),age1(ii-1),age1(ii),blog(ii-1),
     $              blog(ii),bl(i,k) )
               call linint(isoage(k),age1(ii-1),age1(ii),rlog(ii-1),
     $              rlog(ii),rl(i,k) )
               boll(i,k) = SOLBOL-(2.5*bl(i,k))
               if(lphxcol) then
                  itable=4
                  call get_mags(itable,feht,alphafe,boll(i,k),gl(i,k),
     $                 teffl(i,k),filter)
!                  write(*,*)bl(i,k),gl(i,k),teffl(i,k)
                  vvv(i,k) = filter(3)
                  colr(i,k,2) = filter(3) - filter(5)
                  itable = 5
                  call get_mags(itable,feht,alphafe,boll(i,k),gl(i,k),
     $                 teffl(i,k),filter)
                  colr(i,k,3) = filter(4)
                  colr(i,k,4) = filter(4) - filter(7)
               else
                  itable = 8
                  call get_mags(itable,feht,alphafe,boll(i,k),gl(i,k),
     $                 teffl(i,k),filter)
                  vvv(i,k) = filter(2)
                  colr(i,k,2) = filter(2) - filter(3)
                  colr(i,k,3) = filter(4)
                  colr(i,k,4) = filter(4) - filter(5)
               endif
c     write(*,*)teffl(i,k),gl(i,k),bl(i,k),rl(i,k)
               k = k + 1
            else
               ii = ii +1 
            endif
         enddo
c         write(*,*)'k,niso,age',k,niso, isoage(k),age1(ii) 
 190     continue
      enddo

c finsihed reading in the tracks, now read in the isochrones 
c and replace the points

      open(unit=11,status='old')
      open(unit=12)
      do i=1,niso
         read(11,200)hd1
         write(12,200)hd1
         read(11,210)npts,mixlen,ovsh,isoage(i),yy,zz,zeff,feh,afe
         read(11,230)hd2
         do j=1,npts-1
            read(11,250)eep,cext,mass(j),logg(j),logteff(j),logl(j),
     $           num0(j),num235(j),num400(j),vv(j),vi(j),
     $           f606(j),f606f814(j)
 250        format(I3,A1,F9.6,2F10.6,F11.7,3E16.9,8F7.3)
            if(mass(j).lt.minmass) nhit=j
         enddo
         read(11,251)eep,cext,mass(npts),logg(npts),logteff(npts),
     $           logl(npts),num0(npts),num235(npts),vv(npts),
     $           vi(npts),f606(npts),f606f814(npts)
 251     format(I3,A1,F9.6,2F10.6,F11.7,2E16.9,16x,8F7.3)
         num400(npts)=0.000000000E0
         read(11,*)
         read(11,*)
c         write(*,*)npts,nhit,ntrk,npts-nhit+ntrk
         npts2 = npts-nhit + ntrk - nbad
         write(12,220)npts2,mixlen,ovsh,isoage(i),yy,zz,zeff,feh,afe
 220     format('#',I3,F9.6,F8.4,F7.0,F8.4,2E12.4,2F6.2)
         write(12,230)hd2
         do j=1,ntrk
            if(.not. lbad(j) )
     $          write(12,250)j,' ',smass(j),gl(j,i),teffl(j,i),bl(j,i),
     $       0.0,0.0,0.0,vvv(j,i),colr(j,i,2),colr(j,i,3),colr(j,i,4)
         enddo
         k=ntrk+1
         do j=nhit+1,npts
            write(12,250)k,' ',mass(j),logg(j),logteff(j),
     $           logl(j),num0(j),num235(j),num400(j),vv(j),vi(j),
     $           f606(j),f606f814(j)
            k=k+1
         enddo
         write(12,*)
         write(12,*)
      enddo


      close(11)
      close(12)
      stop
      end

C-----------------------------------------------------------
      subroutine linint(xhit,x1,x2,y1,y2,y)
      implicit none
      real*8 xhit,x1,x2,y1,y2,y
c do simple linear interpolation 
      real*8 wx
      wx=(xhit-x1)/(x2 - x1)
      y=(1.0d0-wx)*y1 + wx*y2
      return
      end


c----------------------------------------------------------------------------
c color-teff transforms: 
c NCT=1 Vandenberg & Clem - BVRI
c     4 PHOENIX           - UBVRI & 2MASS JHK
c     5    "              - HST ACS
c     6    "              - HST WFPC2
c     7    "              - SDSS ugriz
c     8 GCTreasury        - V & C/PHOENIX: V, I and ACS F606W, F814W
c Indices for these filters are:
c (Index) | 1      2      3      4      5      6      7      8
c ----------------------------------------------------------------
c         | U      B      V      R      I      J      H      K
c (ACS)   | F435W  F475W  F555W  F606W  F625W  F775W  F814W  F850LP
c (WFPC2) | F336W  F439W  F450W  F555W  F606W  F791W  F814W  F850LW  
c         | u      g      r      i      z
c GCTrsry | B      V      I      F606W  F814W  F606W  F814W (ACS/WFPC2)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine color_init(nct,feh,afe)
      
      real*8 feh,feha,afe,t0,g0,bmv,vmr,vmi,bc,sby,sm1,sc1 
      integer nct

      feha=feh+5.0d-1*afe   !V&C, GW do not account for [a/Fe] variation
      if(nct.eq.1.or.nct.eq.8) CALL BVRI(0,FEHA,G0,T0,BMV,VMR,VMI,BC) !BVRI
      if(nct.eq.4) call phx_color_init(feh,afe,1) !UBVRIJHKs-PHX
      if(nct.eq.5) call phx_color_init(feh,afe,2) !HST-ACSWF 
      if(nct.eq.6) call phx_color_init(feh,afe,3) !HST-WFPC2
      if(nct.eq.7) call phx_color_init(feh,afe,4) !SDSS ugriz
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_mags(nct,feh,afe,mbol,grav,teff,filter)

      real*8 filter(13),color(13),feh,afe,grav,teff,feha,mv
      real*8 bmv,vmr,vmi,bc,mbol
      real*8 f606w_acs,f814w_acs,f606w_wfpc2,f814w_wfpc2
      real clrs(9)
      integer ncol,icol(9),nct,j

      data icol/4,4,8,8,8,8,5,7,8/

!adjust [Fe/H] for semi-empirical color-Teff transformations
      feha=feh+5.0d-1*afe

!***************************
!for Vandenberg & Clem BVRI
!***************************
      if(nct.eq.1.or.nct.eq.8)then
         call bvri(1,feha,grav,teff,bmv,vmr,vmi,bc)
         mv=mbol-bc
         filter(1)=bmv+mv       !B from B-V
         filter(2)=mv           !V
         filter(3)=mv-vmr       !R from V-R
         filter(4)=mv-vmi       !I from V-I
         if(nct.eq.8)then
!***************************
!     for the GC Treasury project
!***************************
!     V&C versions of HST mags: filters 1, 2, 3 = B, V, I
            filter(3)=filter(4) !don't keep R, change to I
! filters 4, 5, 6, 7 = ACS F606W + F814W, WFPC2 F606W + F814W
            filter(4)=f606w_acs(filter(2),filter(3)) !F606W (ACS)
            filter(5)=f814w_acs(filter(2),filter(3)) !F814W (ACS)
            filter(6)=f606w_wfpc2(filter(4),filter(5)) !F606W (WFPC2)
            filter(7)=f814w_wfpc2(filter(4),filter(5)) !F814W (WFPC2)
         endif
!***************************
!for PHOENIX synthetic
!***************************
      elseif(nct.gt.3.and.nct.lt.8)then !PHOENIX
         call phx_color_interp(teff,grav,color,nct)
c         write(*,*)mbol,color(1)
         mv=mbol+color(1)
         if(nct.eq.4)then       !UBVRIJHK
            filter(1)=color(2)+color(3)+mv !U from U-B and B-V
            filter(2)=color(3)+mv !B from B-V
            filter(3)=mv        !V
            do j=4,icol(nct)
               filter(j)=mv-color(j) ![filter] from V-[filter]
            enddo
c            write(*,*)'Filter:',filter(1),filter(2),filter(3),filter(4)
         elseif(nct.eq.5.or.nct.eq.6)then !HST ACS/WFPC2
            do j=1,icol(nct)
               filter(j)=mv-color(j+1) ![filter] from V-[filter]
            enddo
         elseif(nct.eq.7)then   !SDSS
            filter(1)=color(2)+mv !g: here "mv" is really "mg"
            do j=2,icol(nct)
               filter(j)=filter(j-1)-color(j)
            enddo
         endif
!anything else...
      else
         stop "INVALID COLOR-TEFF TRANSFORMATION SPECIFIED"
      endif

      return
      end
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function f606w_acs(v,i)
c     Sirianni et al 2005 HST-ACS transformation equation from V and I to
c     F606W (Vega mag system)
      real*8 v,i
      if(v-i.ge.4.0d-1) then
         f606w_acs=v-2.6331d1+2.6398d1-3.4d-1*(v-i)+3.8d-2*(v-i)*(v-i)
      else
         f606w_acs=v-2.6394d1+2.6398d1-1.53d-1*(v-i)-9.6d-2*(v-i)*(v-i)
      endif
       
      return
      end
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function f606w_wfpc2(v,i)
c     Sirianni et al 2005 HST-ACS transformation equation from ACS to WFPC2
c     Here V=F606W and I=F814W
      real*8 v,i
 
      f606w_wfpc2= v + 1.6d-2*(v-i) - 3.0d-2*(v-i)*(v-i) !- 4.168d0
      return
      end
  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function f814w_acs(v,i)
c     Sirianni et al 2005 HST-ACS transformation equation from V and I to
c     F814W (Vega mag system)
      real*8 v,i
      if(v-i.ge.1.0d-1)then
         f814w_acs=i-2.5496d1+2.5501d1+1.4d-2*(v-i)-1.5d-2*(v-i)*(v-i)
      else
         f814w_acs=i-2.5489d1+2.5501d1-4.1d-2*(v-i)+9.3d-2*(v-i)*(v-i)
      endif
       
      return
      end
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function f814w_wfpc2(v,i)
c     Sirianni et al 2005 HST-ACS transformation equation from ACS to WFPC2
c     Here V=F606W and I=F814W
      real*8 v,i
      f814w_wfpc2= i + 2.3d-2*(v-i) - 2.0d-3*(v-i)*(v-i) !- 4.601d0
      return
      end
 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                       END OF ALL_COLOR FILE                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       PHOENIX COLOR INTERPOLATION IN TEFF, LOG G, [Fe/H],& [a/Fe]      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                        C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

ccccccccccccccccccccccc subroutine phx_color_init cccccccccccccccccccccccc
c this subroutine determines which filter set will be used, reads in     c
c the appropriate tables, and stores them in arrays in a common block    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine phx_color_init(feh,afe,filter)

      implicit none
      integer filter, iz, iafe, nz, nzafe
      parameter(nz=10,nzafe=7)
c      parameter(nz=7,nzafe=7)
      real*8 feh, afe, z_0(nz), z_afe(nzafe)
      data z_0/-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5/
c      data   z_0/-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5/
c      data z_afe/-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5/
c filter is an integer representing the choice of filter set
c options include:
c     1. J-C UBVRI and 2MASS JHK
c     2. HST-ACS WFC
c     3. HST-WFPC2
c     4. SDSS ugriz

c set [a/Fe] index parameter
      iafe=0
      if(afe.eq.-0.2d0) iafe=1
      if(afe.eq. 0.0d0) iafe=2
      if(afe.eq. 0.2d0) iafe=3
      if(afe.eq. 0.4d0) iafe=4
      if(afe.eq. 0.6d0) iafe=5
      if(afe.eq. 0.8d0) iafe=6
      write(*,*)'phoenix',feh,afe,filter
      if(iafe.eq.0) then
         write(*,*)'afe=',afe
         write(*,*)"PROBLEM WITH [a/Fe] IN ISOCHRONE FILE"
         stop
      endif

c locate [Fe/H] in the array of available values
      call hunt(z_0,nz,feh,iz)

c check to make sure [Fe/H] contained within available values
c note that the alpha-enhanced colors begin at -2.5 (not -4)
c$$$      if(iz.eq.0) stop '[Fe/H] out of bounds!'
c$$$      if(iz.eq.nz) stop '[Fe/H] out of bounds!'
      if(iafe.eq.2.and.iz.lt.2) iz=2
      if(iafe.ne.2.and.iz.lt.5) iz=5
      if(iz.gt.nz-2) iz=nz-2
      if(iafe.eq.6.and.iz.gt.nz-3) iz=nz-3

c read in the required tables for a range of [Fe/H] at fixed [a/Fe]    
      call phx_color_read(iz,iafe,filter)

c interpolate in [Fe/H], only need once
      call z_interp(feh)

c all set for interpolation in T_eff and log(g)      
      return
      end


ccccccccccccccccccccccccc subroutine phx_color_read ccccccccccccccccccccccc
c reads in color tables based on specified Z, [a/Fe], filter set          c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine phx_color_read(iz,iafe,filter)

c put in a hack here, so that ACS tables (filter=2) are 
c saved to a separate common  block.... this means we only will have 
c to read in things once for V,I F606W, F814W colours

      implicit none
      integer iz,iafe,filter,nmax,ifile,nt,it,nz,nfltr
      parameter(nt=66,nz=10,nfltr=4,nmax=1000,ifile=21)
c      parameter(nt=66,nz=7,nfltr=4,nmax=1000,ifile=21)
      integer icol(nfltr),i,j,k
      integer ncol,itnum,isize,itemp
      integer ncolA,itnumA,isizeA,itempA
      character filename(4)*23,zfile(nz)*5,afile(6)*3,suffix(nfltr)*9
      real*8 zl,fz,coltbl,teff,ggl
      real*8 coltblA,teffA,gglA
      data zfile/'Zm4d0','Zm3d5','Zm3d0','Zm2d5','Zm2d0','Zm1d5',
     *     'Zm1d0','Zm0d5','Zp0d0','Zp0d5'/
      data afile/'am2','ap0','ap2','ap4','ap6','ap8'/
      data suffix/'STD_2MASS','HST_ACSWF','HST_WFPC2','SDSSugriz'/
      data icol/8,8,8,5/
      common/phx/coltbl(4,nmax,13),teff(4,nmax),ggl(4,nmax),
     *     itemp(4,nt),itnum(4),isize(4),ncol
      common/phxA/coltblA(4,nmax,13),teffA(4,nmax),gglA(4,nmax),
     *     itempA(4,nt),itnumA(4),isizeA(4),ncolA
      common/z/zl(4),fz(4)
c      data zl/-4.0,-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5/

      ncol=icol(filter)
      ncolA=icol(filter)   !NOT USIN THIS RIGHT NOW.....
      do i=1,4
c set up filenames based on the input variables
         filename(i)(1:4)="phx/"
         filename(i)(5:9)=zfile(iz+i-2)
         filename(i)(10:12)=afile(iafe)
         filename(i)(13:13)="."
         do j=1,nfltr
            if(filter.eq.j) filename(i)(14:23)=suffix(j)
         enddo
c open table for reading
         write(*,*)filename(i)
         open(ifile,file=filename(i),status='old')
         read(ifile,*)
c read data in table
         do j=1,nmax
c***BCC 1/2007 my  hack -- ASSUMES Z is the same in each of the files
            if(filter.ne.2) then
               read(ifile,*,end=2) 
     *           teff(i,j),ggl(i,j),zl(i),(coltbl(i,j,k),k=1,ncol)
               if(j.ge.2.and.teff(i,j).eq.teff(i,j-1).and.
     *              ggl(i,j).eq.ggl(i,j-1)) then
c print warning message if duplicate lines exist in table
                  write(*,*) "DUPLICATE LINE IN COLOR TABLE"
                  write(*,*)zl(i),teff(i,j),teff(i,j-1),
     $                 ggl(i,j),ggl(i,j-1)
            endif
            elseif(filter.eq.2) then
               read(ifile,*,end=2) 
     *           teffA(i,j),gglA(i,j),zl(i),(coltblA(i,j,k),k=1,ncol)
               if(j.ge.2.and.teffA(i,j).eq.teffA(i,j-1).and.
     *              gglA(i,j).eq.gglA(i,j-1)) then
c print warning message if duplicate lines exist in table
                  write(*,*) "DUPLICATE LINE IN COLOR TABLE"
                  write(*,*)zl(i),teffa(i,j),teffa(i,j-1),
     $                 ggla(i,j),ggla(i,j-1)
               endif
            else
               write(*,*)'something wrong here'
               stop
            endif


         enddo
c close data file
 2       close(ifile)
c isize counts number of data points, itnum counts number of distinct T's,
c itemp gives location of each T value in the j-array
         if(filter.ne.2) then 
            isize(i)=j-1
            itemp(i,1)=1
            itnum(i)=1
            do j=2,isize(i)
               if(teff(i,j).gt.teff(i,j-1)) then
                  itnum(i)=itnum(i)+1
                  itemp(i,itnum(i))=j
               endif
            enddo
         else

C***BCC more of the hack

            isizeA(i)=j-1
            itempA(i,1)=1
            itnumA(i)=1
            do j=2,isizeA(i)
               if(teffA(i,j).gt.teffA(i,j-1)) then
                  itnumA(i)=itnumA(i)+1
                  itempA(i,itnumA(i))=j
               endif
            enddo
         endif



      enddo
      
      return
      end

ccccccccccccccccccccccc subroutine phx_color_interp ccccccccccccccccccccccc
c     interpolates the mags and colors based on T, g, [Fe/H], [a/Fe]      c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine phx_color_interp(tl,gl,color,nct)
      implicit none
      integer nmax,nt,i,j,k,iii,jj,inewt,inewg,tint,gint,nct
      integer ncol,itnum,isize,itemp
      integer ncolA,itnumA,isizeA,itempA
      parameter(nmax=1000,nt=66)
      real*8 zl,fz,fg(4),ft(4),coltbl,tl,gl,teff,ggl
      real*8 coltblA,teffA,gglA
      real*8 cln,t,qt(4),qg(4),colr(4,13),col(4,13,4),color(13)
      common/phx/coltbl(4,nmax,13),teff(4,nmax),ggl(4,nmax),itemp(4,nt),
     *     itnum(4),isize(4),ncol
      common/phxA/coltblA(4,nmax,13),teffA(4,nmax),gglA(4,nmax),
     *     itempA(4,nt),itnumA(4),isizeA(4),ncolA
      common/z/zl(4),fz(4)
      
      cln=dlog(1.0d1)



      if(nct.ne.5) then
      do iii=1,4
c     locate T in the Teff array
         t=dexp(cln*tl)
         do i=1,itnum(iii)-1
            if(t.ge.teff(iii,itemp(iii,i)).and.
     *           t.lt.teff(iii,itemp(iii,i+1))) inewt=i
         enddo
         if(inewt.lt.2) inewt=2
         if(inewt.gt.itnum(iii)-2) inewt=itnum(iii)-2

c     find interpolation coeff.'s in T
         do i=1,4
            qt(i)=teff(iii,itemp(iii,inewt+i-2))
         enddo
         call interp(qt,ft,t,4)

c     locate Log G in the Log G array for each T
         do j=1,4
            jj=inewt+j-2
            do k=itemp(iii,jj),itemp(iii,jj+1)-1
               if(gl.ge.ggl(iii,k).and.gl.lt.ggl(iii,k+1)) inewg=k
            enddo
            if(inewg.lt.2) inewg=2
            if(inewg.gt.itemp(iii,jj+1)-2) inewg=itemp(iii,jj+1)-2
            do k=1,4
               qg(k)=ggl(iii,inewg+k-2)
            enddo

c     find interpolation coefficients in Log G
            call interp(qg,fg,gl,4)

c     g-interpolation
            do k=1,ncol
               col(iii,k,j) = fg(1)*coltbl(iii,inewg-1,k)
     *                      + fg(2)*coltbl(iii,inewg  ,k)
     *                      + fg(3)*coltbl(iii,inewg+1,k)
     *                      + fg(4)*coltbl(iii,inewg+2,k)
            enddo               !k-loop
         enddo                  !j-loop

c     T-interpolation
         do i=1,ncol
            colr(iii,i)=ft(1)*col(iii,i,1)+ft(2)*col(iii,i,2)
     *           +ft(3)*col(iii,i,3)+ft(4)*col(iii,i,4)
         enddo                  !i-loop
        
      enddo                     !iii-loop

      else
c do the interpolation for the ACS filters, which are in different common block
      do iii=1,4
c     locate T in the Teff array
         t=dexp(cln*tl)
         do i=1,itnumA(iii)-1
            if(t.ge.teffA(iii,itempA(iii,i)).and.
     *           t.lt.teffA(iii,itempA(iii,i+1))) inewt=i
         enddo
         if(inewt.lt.2) inewt=2
         if(inewt.gt.itnumA(iii)-2) inewt=itnumA(iii)-2

c     find interpolation coeff.'s in T
         do i=1,4
            qt(i)=teffA(iii,itempA(iii,inewt+i-2))
         enddo
         call interp(qt,ft,t,4)

c     locate Log G in the Log G array for each T
         do j=1,4
            jj=inewt+j-2
            do k=itempA(iii,jj),itempA(iii,jj+1)-1
               if(gl.ge.gglA(iii,k).and.gl.lt.gglA(iii,k+1)) inewg=k
            enddo
            if(inewg.lt.2) inewg=2
            if(inewg.gt.itempA(iii,jj+1)-2) inewg=itempA(iii,jj+1)-2
            do k=1,4
               qg(k)=gglA(iii,inewg+k-2)
            enddo

c     find interpolation coefficients in Log G
            call interp(qg,fg,gl,4)

c     g-interpolation
            do k=1,ncol
               col(iii,k,j) = fg(1)*coltblA(iii,inewg-1,k)
     *                      + fg(2)*coltblA(iii,inewg  ,k)
     *                      + fg(3)*coltblA(iii,inewg+1,k)
     *                      + fg(4)*coltblA(iii,inewg+2,k)
            enddo               !k-loop
         enddo                  !j-loop

c     T-interpolation
         do i=1,ncol
            colr(iii,i)=ft(1)*col(iii,i,1)+ft(2)*col(iii,i,2)
     *           +ft(3)*col(iii,i,3)+ft(4)*col(iii,i,4)
         enddo                  !i-loop
        
      enddo                     !iii-loop
      endif



c     Z-interpolation
      do i=1,ncol
         color(i)=fz(1)*colr(1,i) + fz(2)*colr(2,i)
     *        + fz(3)*colr(3,i) + fz(4)*colr(4,i)
      enddo

      return
      end

ccccccccccccccccccccccccc subroutine z_interp cccccccccccccccccccccccccccc
c      Get interpolation coefficients in [Fe/H], only needed once.       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine z_interp(feh)
      implicit none
      real*8 zl,fz,feh
      common/z/zl(4),fz(4)

      call interp(zl,fz,feh,4)

      return
      end

**********************************************************************
*                                HUNT                                *
**********************************************************************
      subroutine hunt(xx,n,x,jlo)
C Subroutine HUNT searches an array of length N to return a value
C of JLO such that X is located between array elements JLO and JLO+1.
C This search is based on the last value of JLO returned. JLO = 0 or
C JLO = N is returned to indicate that X is out of range of the array.
C This routine is taken from Numerical Recipes, pp.91-92
      implicit none
      integer n, jlo, jhi, inc, jm
      logical lascnd
      real*8 xx(n), x
      
      lascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
C Input guess not useful; go immediately to bisection.
         jlo=0
         jhi=n+1
         goto 3
      endif
C Set up the hunting increment. 
      inc=1
C Hunt up.
      if(x.ge.xx(jlo).eqv.lascnd)then
 1       jhi=jlo+inc
         if(jhi.gt.n)then
C Done hunting, since off end of the table.
            jhi=n+1
         else if(x.ge.xx(jhi).eqv.lascnd)then
C Not done hunting.
            jlo=jhi
            inc=inc+inc
            goto 1
C Done hunting.
         endif
      else
C Hunt down.
         jhi=jlo
 2       jlo=jhi-inc
         if(jlo.lt.1) then
C Done hunting, since off end of table.
            jlo=0
         else if(x.lt.xx(jlo).eqv.lascnd)then
C Not done hunting.
            jhi=jlo
            inc=inc+inc
            goto 2
C Done hunting.
         endif
      endif
C Hunt is done; begin the final bisection phase.
 3    if(jhi-jlo.eq.1) return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.lascnd)then
         jlo=jm
      else
         jhi=jm
      endif
      goto 3
      
      return
      end
      
C******************************************************************************
C               N-PT. Lagrangian Interpolation Coefficients
C******************************************************************************
      subroutine interp(a,b,x,n)
c {a} are the tabulated values for use in interpolation
c {b} are coefficients of the interpolating polynomial
c  x  is the abscissa to be interpolated
c  n  is the number of points to be used, interpolating polynomial
c     has order n-1 
      implicit none
      integer i,j,n
      real*8 a(n),b(n),x
      do i=1,n
         b(i)=1.0d0
         do j=1,n
            if(j.ne.i) b(i)=b(i)*(x-a(j))/(a(i)-a(j))
         enddo
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                           END of PHX_COLOR.F                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c BEGIN VRSUB SUBROUTINES: comprises all the useful SR's in the V-R code c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c BEGIN VANDENBERG & CLEM COLOR TRANSFORM SUBROUTINES                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE BVRI(MODE,FE,GV,TEFF,BMV,VMR,VMI,BCV)
C ----------------------------------------------------------------------
C *** THIS READS, AND INTERPOLATES IN, bvrilo.data AND bvrihi.data (SEE
C     VANDENBERG & CLEM 2003, AJ, 126, IN PRESS).            
C *** SET MODE=0 TO READ THE TABLES AND TO INTERPOLATE THEM TO THE
C                DESIRED VALUE OF [Fe/H]  (= FE IN THE ARGUMENT LIST)   
C *** SET MODE=1 TO INTERPOLATE FOR THE COLORS AT THE DESIRED VALUES OF
C                [Fe/H], log g, and log Teff (GV AND TEFF REPRESENT THE
C                LAST TWO OF THESE QUANTITIES IN THE ARGUMENT LIST).
C                (NOTE THAT THE [Fe/H] INTERPOLATION IS CARRIED OUT ONLY
C                WHEN THE SUBROUTINE IS CALLED WITH MODE=0 OR MODE=-1.)
C *** SET MODE=-1 IF THE INPUT TABLES, WHICH HAVE ALREADY BEEN READ
C                USING MODE=0, ARE TO BE RE-INTERPOLATED TO BE
C                CONSISTENT WITH A NEW VALUE OF [Fe/H].        
C *** THE OUTPUT CONSISTS OF THE B-V, V-R, AND V-I COLORS, ALONG WITH
C     THE BOLOMETRIC CORRECTION TO V (ON THE SCALE WHERE THE SUN HAS
C     M_bol = 4.75 and M_V = 4.84).  THESE QUANTITIES ARE REPRESENTED,
C     IN TURN, BY BMV, VMR, VMI, AND BCV IN THE ARGUMENT LIST. (B AND V
C     ARE ON THE JOHNSON SYSTEM, R AND I ON THE COUSINS SYSTEM.)
C *** MODE MUST BE SET TO ZERO THE FIRST TIME THAT THIS SUBROUTINE IS
C     CALLED, AS IT IS ONLY WHEN MODE=0 THAT BVRILO.DATA AND BVRIHI.DATA
C     ARE READ.  ONCE THE TABLES HAVE BEEN INTERPOLATED TO A PARTICULAR
C     VALUE OF [Fe/H] (USING MODE=0 OR MODE=-1), INTERPOLATIONS FOR THE
C     COLORS APPROPRIATE TO INPUT VALUES OF [Fe/H], log g, and log Teff
C     MAY BE CARRIED OUT *ANY NUMBER OF TIMES* USING MODE=1 (THE ONLY
C     VALUE OF MODE FOR WHICH COLOR INFORMATION IS OBTAINED).
C *** NOTE THAT THE BVRILO.DATA AND BVRIHI.DATA FILES ARE "ATTACHED" TO
C     UNITS 10 AND 11, RESPECTIVELY.  WHENEVER AN [Fe/H] INTERPOLATION 
C     IS CARRIED OUT, A MESSAGE IS SENT TO UNIT 6 ("ATTACHED" TO THE  
C     TERMINAL) TO INDICATE THE VALUE OF [FE/H] THAT APPLIES TO
C     SUBSEQUENT INTERPOLATIONS FOR log g and log Teff. 
C *** CHKBVRI.FOR (ALSO ON DISK) MAY BE USED TO TEST THE INTERPOLATION
C     SUBROUTINE
C *** NOTE, IN THE OPEN STATEMENTS BELOW, THAT THE NAMES OF THE INPUT
C *** DATA FILES ARE WRITTEN IN LOWER-CASE TEXT
C ----------------------------------------------------------------------
*23456789012345678901234567890123456789012345678901234567890123456789012
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(34,7,6,4),B(34,7,4),BV(11,12,6),VR(11,12,6), 
     1   VI(11,12,6),BC(11,12,6),C(11,12),D(11,12),E(11,12),F(11,12), 
     2   H(4,4),TB(34),TT(34),GG(7),TA(11),T(11),G(12),FEH(6),P(4),Q(4),    
     3   R(4),S(4),AG(4),AT(4),X(11),Y(11)
      CHARACTER*19 NAMCOL
      CHARACTER*10 NAMFE                                       
      INTEGER ILO, IHI
      PARAMETER(ILO=17,IHI=18)
      SAVE
 1000 FORMAT(I3,13X,I3,13X,I3,13X,I3)                                   
 1001 FORMAT(13F6.0)                                                    
 1002 FORMAT(20F4.1)                                                    
 1003 FORMAT(13F6.3)                                                    
 1004 FORMAT('     Color grid interpolated to [FE/H] =',F6.2)           
 1005 FORMAT(' log g =',F6.3,'  log Teff =',F7.4,' ARE OUTSIDE THE',          
     1   1X,'RANGE OF THE COLOR TABLES')                          
 1006 FORMAT(A10,F5.2,A19)
 1007 FORMAT(1X,' **** INPUT DATA FILE DOES NOT EXIST **** ')      
      IF(MODE) 11,3,16                                                  
C *** WHEN MODE=0 THE COLOR TRANSFORMATION TABLES ARE READ
    3 TBND=LOG10(5000.)                                                
      GBND=2.
      OPEN(UNIT=ILO,FILE='vdb/bvrilo.data',ERR=49,STATUS='OLD')
      OPEN(UNIT=IHI,FILE='vdb/bvrihi.data',ERR=49,STATUS='OLD')      
      READ(ILO,1000) NT,NG,NFE,NDX                                       
      READ(ILO,1001) (TA(I),I=1,NT)                                      
      READ(ILO,1002) (G(I),I=1,NG)                                       
      DO 4 I=1,NT                                                       
    4 T(I)=LOG10(TA(I))                                                
      DO 5 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                    
      DO 5 J=1,NG                                                       
      READ(ILO,1003) (BV(I,J,K),I=1,NT)                                  
    5 CONTINUE                                                          
      DO 6 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                 
      DO 6 J=1,NG                                                       
      READ(ILO,1003) (VR(I,J,K),I=1,NT)                                  
    6 CONTINUE                                                          
      DO 7 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                   
      DO 7 J=1,NG                                                       
      READ(ILO,1003) (VI(I,J,K),I=1,NT)                                  
    7 CONTINUE                                                          
      DO 8 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                   
      DO 8 J=1,NG
      READ(ILO,1003) (BC(I,J,K),I=1,NT)
    8 CONTINUE                                             
      READ(IHI,1000) NTT,NGG,NFE,NDXX                                    
      READ(IHI,1001) (TB(I),I=1,NTT)                                     
      READ(IHI,1002) (GG(I),I=1,NGG)                                     
      DO 9 I=1,NTT                                                      
    9 TT(I)=LOG10(TB(I))                                               
      DO 10 L=1,NDXX                                                    
      DO 10 K=1,NFE                                                     
      READ(IHI,1006) NAMFE,FEH(K),NAMCOL
      DO 10 J=1,NGG                                                     
      READ(IHI,1003) (A(I,J,K,L),I=1,NTT)
   10 CONTINUE
      CLOSE(UNIT=ILO,STATUS='KEEP')
      CLOSE(UNIT=IHI,STATUS='KEEP')         
C *** WHEN MODE=0 OR MODE=-1, COLOR TRANSFORMATION TABLES ARE
C *** CREATED FOR THE INPUT [Fe/H] VALUE USING LINEAR INTERPOLATION
   11 DO 12 M=2,NFE                                                     
      K=M-1                                                             
      IF(FE.LE.FEH(M)) GO TO 13                                         
   12 CONTINUE                                                          
   13 M=K+1                                                             
      SLOPE=(FE-FEH(K))/(FEH(M)-FEH(K))                                 
      DO 14 J=1,NG                                                      
      DO 14 I=1,NT                                                      
      C(I,J)=BV(I,J,K)+SLOPE*(BV(I,J,M)-BV(I,J,K))                      
      D(I,J)=VR(I,J,K)+SLOPE*(VR(I,J,M)-VR(I,J,K))                      
      E(I,J)=VI(I,J,K)+SLOPE*(VI(I,J,M)-VI(I,J,K))                      
      F(I,J)=BC(I,J,K)+SLOPE*(BC(I,J,M)-BC(I,J,K))                      
   14 CONTINUE                                                          
      DO 15 L=1,NDXX                                                    
      DO 15 J=1,NGG                                                     
      DO 15 I=1,NTT                                                     
      B(I,J,L)=A(I,J,K,L)+SLOPE*(A(I,J,M,L)-A(I,J,K,L))                 
   15 CONTINUE                                                          
      !WRITE(0,1004) FE                                                  
C *** WHEN MODE=0 OR MODE=-1, CONTROL RETURNS TO THE CALLING
C *** PROGRAM ONCE THE [Fe/H] INTERPOLATION IS CARRIED OUT (I.E., NO
C *** INTERPOLATIONS ARE PERFORMED FOR INPUT VALUES OF GV AND TEFF)
      GO TO 50
C *** WHEN MODE=1, INTERPOLATIONS ARE CARRIED OUT FOR THE INPUT VALUES
C *** OF [Fe/H], log g, AND log Teff (= FE, GV, AND TEFF IN THE ARGUMENT
C *** LIST)                                    
   16 IF(TEFF.LE.TBND) GO TO 18                                         
      IF(GV.GE.GBND) GO TO 40                                           
      IF(TEFF.LE.T(NT)) GO TO 18                                        
C      WRITE(0,1005) GV,TEFF                                             
C *** EXECUTION HALTS WITH A "STOP30" CODE IF THE INPUT TEMPERATURE IS
C *** OUTSIDE THE RANGE OF THE TABLES ON THE HIGH SIDE
C      STOP30                                                            
      BCV=-99.999d0
      BMV=-99.999d0
      VMR=-99.999d0
      VMI=-99.999d0
      RETURN
C *** THE NEXT SECTION ASSUMES THAT THE LOW-TEMPERATURE TABLES ARE THE
C *** RELEVANT ONES TO USE
   18 NGM=NG-1                                                          
      DO 19 I=3,NGM                                                     
      MG=I-2                                                            
      IF(GV.LE.G(I)) GO TO 20                                           
   19 CONTINUE                                                          
      IF(GV.LE.G(NG)) GO TO 20
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log g CONSIDERED IN THE
C *** TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT TWO
C *** LINES SHOULD BE ACTIVATED                                          
C     WRITE(6,1005) GV,TEFF                                             
C     STOP31                                                            
   20 R(1)=G(MG)                                                        
      R(2)=G(MG+1)                                                      
      R(3)=G(MG+2)                                                      
      R(4)=G(MG+3)                                                      
      CALL LGRAN4(R,AG,GV)                                              
      NTM=NT-1                                                          
      DO 25 I=3,NTM                                                     
      MT=I-2                                                            
      IF(TEFF.LE.T(I)) GO TO 30                                         
   25 CONTINUE                                                          
      IF(TEFF.LE.T(NT)) GO TO 30                                        
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log Teff CONSIDERED IN
C *** THE TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT
C *** TWO LINES SHOULD BE ACTIVATED.
C     WRITE(6,1005) GV,TEFF                                             
C     STOP32                                                            
   30 R(1)=T(MT)                                                        
      R(2)=T(MT+1)                                                      
      R(3)=T(MT+2)                                                      
      R(4)=T(MT+3)                                                      
      CALL LGRAN4(R,AT,TEFF)                                            
      L=MG-1                                                            
      DO 35 I=1,4                                                       
      L=L+1                                                             
      P(I)=AT(1)*C(MT,L)+AT(2)*C(MT+1,L)+AT(3)*C(MT+2,L)+AT(4)*C(MT+3,L)
      Q(I)=AT(1)*D(MT,L)+AT(2)*D(MT+1,L)+AT(3)*D(MT+2,L)+AT(4)*D(MT+3,L)
      R(I)=AT(1)*E(MT,L)+AT(2)*E(MT+1,L)+AT(3)*E(MT+2,L)+AT(4)*E(MT+3,L)
   35 S(I)=AT(1)*F(MT,L)+AT(2)*F(MT+1,L)+AT(3)*F(MT+2,L)+AT(4)*F(MT+3,L)
C *** 4-POINT LAGRANGIAN INTERPOLATION IS USED TO DERIVE THE VALUES OF
C *** B-V, V-R, V-I, AND BC_V CORRESPONDING TO THE INPUT VALUES OF
C *** [Fe/H], log g, and log Teff.
      BMV=AG(1)*P(1)+AG(2)*P(2)+AG(3)*P(3)+AG(4)*P(4)                   
      VMR=AG(1)*Q(1)+AG(2)*Q(2)+AG(3)*Q(3)+AG(4)*Q(4)                   
      VMI=AG(1)*R(1)+AG(2)*R(2)+AG(3)*R(3)+AG(4)*R(4)                   
      BCV=AG(1)*S(1)+AG(2)*S(2)+AG(3)*S(3)+AG(4)*S(4)                  
C *** SPLINE INTERPOLATION IS USED TO FIND THE VALUE OF BC_V AT LOW
C *** TEMPERATURES.
      DO 38 I=1,NT                                                      
      X(I)=T(I)                                                         
   38 Y(I)=AG(1)*F(I,MG)+AG(2)*F(I,MG+1)+AG(3)*F(I,MG+2)+AG(4)*F(I,MG+3)
      CALL INTEP(TEFF,BCV,X,Y,NT,IER)                                  
      GO TO 50                                                          
C *** THE NEXT SECTION ASSUMES THAT THE HIGH-TEMPERATURE TABLES ARE THE
C *** RELEVANT ONES TO USE.
   40 IF(TEFF.LE.TT(NTT)) GO TO 42                                      
      WRITE(0,1005) GV,TEFF                                             
      STOP33                                                            
   42 NGM=NGG-1                                                         
      DO 43 I=3,NGM                                                     
      MG=I-2                                                            
      IF(GV.LE.GG(I)) GO TO 44                                          
   43 CONTINUE                                                          
      IF(GV.LE.GG(NGG)) GO TO 44                                         
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log g CONSIDERED IN THE
C *** TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT TWO
C *** LINES SHOULD BE ACTIVATED.
C     WRITE(6,1005) GV,TEFF                                             
C     STOP34                                                            
   44 R(1)=GG(MG)                                                       
      R(2)=GG(MG+1)                                                     
      R(3)=GG(MG+2)                                                     
      R(4)=GG(MG+3)                                                     
      CALL LGRAN4(R,AG,GV)                                              
      NTM=NTT-1                                                         
      DO 45 I=3,NTM                                                     
      MT=I-2                                                            
      IF(TEFF.LE.TT(I)) GO TO 46                                        
   45 CONTINUE                                                          
      IF(TEFF.LE.TT(NTT)) GO TO 46                                      
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log Teff CONSIDERED IN
C *** THE TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT
C *** TWO LINES SHOULD SHOULD BE ACTIVATED.
C     WRITE(6,1005) GV,TEFF                                             
C     STOP35                                                            
   46 R(1)=TT(MT)                                                       
      R(2)=TT(MT+1)                                                     
      R(3)=TT(MT+2)                                                     
      R(4)=TT(MT+3)                                                     
      CALL LGRAN4(R,AT,TEFF)                                            
      L=MG-1                                                            
      DO 48 I=1,4                                                       
      L=L+1                                                             
      DO 48 K=1,NDXX                                                    
      H(I,K)=AT(1)*B(MT,L,K)+AT(2)*B(MT+1,L,K)+AT(3)*B(MT+2,L,K)+       
     1   AT(4)*B(MT+3,L,K)                                              
   48 CONTINUE                                                          
C *** 4-POINT LAGRANGIAN INTERPOLATION IS USED TO DERIVE THE VALUES OF
C *** B-V, V-R, V-I, AND BC_V CORRESPONDING TO THE INPUT VALUES OF
C *** [Fe/H], log g, and log Teff.
      BMV=AG(1)*H(1,1)+AG(2)*H(2,1)+AG(3)*H(3,1)+AG(4)*H(4,1)           
      VMR=AG(1)*H(1,2)+AG(2)*H(2,2)+AG(3)*H(3,2)+AG(4)*H(4,2)           
      VMI=AG(1)*H(1,3)+AG(2)*H(2,3)+AG(3)*H(3,3)+AG(4)*H(4,3)           
      BCV=AG(1)*H(1,4)+AG(2)*H(2,4)+AG(3)*H(3,4)+AG(4)*H(4,4)
      GO TO 50
   49 WRITE(0,1007)
      STOP          
   50 RETURN                                                            
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
             
      SUBROUTINE UVBY(MODE,FE,GV,TEFF,SBY,SM1,SC1,BCV)           
C ----------------------------------------------------------------------
C *** THIS READS, AND INTERPOLATES IN, UVBYLO.DATA AND UVBYHI.DATA (SEE
C     CLEM, VANDENBERG, GRUNDAHL, & BELL, AJ, IN PRESS).            
C *** SET MODE=0 TO READ THE TABLES AND TO INTERPOLATE THEM TO THE
C                DESIRED VALUE OF [Fe/H]  (= FE IN THE ARGUMENT LIST)   
C *** SET MODE=1 TO INTERPOLATE FOR THE COLORS AT THE DESIRED VALUES OF
C                [Fe/H], log g, and log Teff (GV AND TEFF REPRESENT THE
C                LAST TWO OF THESE QUANTITIES IN THE ARGUMENT LIST).
C                (NOTE THAT THE [Fe/H] INTERPOLATION IS CARRIED OUT ONLY
C                WHEN THE SUBROUTINE IS CALLED WITH MODE=0 OR MODE=-1.)
C *** SET MODE=-1 IF THE INPUT TABLES, WHICH HAVE ALREADY BEEN READ
C                USING MODE=0, ARE TO BE RE-INTERPOLATED TO BE
C                CONSISTENT WITH A NEW VALUE OF [Fe/H].        
C *** THE OUTPUT CONSISTS OF THE b-y, m1, AND c1 INDICES, ALONG 
C     WITH THE BOLOMETRIC CORRECTION TO V (ON THE SCALE WHERE THE SUN 
C     HAS M_bol = 4.75 and M_V = 4.84).  THESE QUANTITIES ARE 
C     REPRESENTED, IN TURN, BY SBY, SM1, SC1, AND BCV IN THE 
C     ARGUMENT LIST.   NOTE THAT IT IS ASSUMED THE BOLOMETRIC 
C     CORRECTIONS TO STROMGREN y ARE IDENTICAL TO THOSE IN JOHNSON V.
C *** MODE MUST BE SET TO ZERO THE FIRST TIME THAT THIS SUBROUTINE IS
C     CALLED, AS IT IS ONLY WHEN MODE=0 THAT UVBYLO.DATA AND UVBYHI.DATA
C     ARE READ.  ONCE THE TABLES HAVE BEEN INTERPOLATED TO A PARTICULAR
C     VALUE OF [Fe/H] (USING MODE=0 OR MODE=-1), INTERPOLATIONS FOR THE
C     COLORS APPROPRIATE TO INPUT VALUES OF [Fe/H], log g, and log Teff
C     MAY BE CARRIED OUT *ANY NUMBER OF TIMES* USING MODE=1 (THE ONLY
C     VALUE OF MODE FOR WHICH COLOR INFORMATION IS OBTAINED).
C *** NOTE THAT THE UVBYLO.DATA AND UVBYHI.DATA FILES ARE "ATTACHED" TO
C     UNITS 12 AND 13, RESPECTIVELY.  WHENEVER AN [Fe/H] INTERPOLATION 
C     IS CARRIED OUT, A MESSAGE IS SENT TO UNIT 6 ("ATTACHED" TO THE  
C     TERMINAL) TO INDICATE THE VALUE OF [FE/H] THAT APPLIES TO
C     SUBSEQUENT INTERPOLATIONS FOR log g and log Teff. 
C *** CHKUVBY.FOR (ALSO ON DISK) MAY BE USED TO TEST THE INTERPOLATION
C     SUBROUTINE
C ----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 M1(13,12,8)
      DIMENSION A(51,7,8,4),B(51,7,4),BY(13,12,8),C1(13,12,8),
     1  BC(13,12,8),C(13,12),D(13,12),E(13,12),F(13,12),H(4,4),TB(51),
     2  TT(51),GG(7),TA(13),T(13),G(12),FEH(8),O(4),P(4),Q(4),R(4),
     3  AG(4),AT(4),X(13),Y(13)
      CHARACTER*19 NAMCOL
      CHARACTER*10 NAMFE
      INTEGER ILO,IHI
      PARAMETER(ILO=19,IHI=20)
      SAVE
 1000 FORMAT(I3,13X,I3,13X,I3,13X,I3)                                   
 1001 FORMAT(13F6.0)                                                    
 1002 FORMAT(20F4.1)                                                    
 1003 FORMAT(13F6.3)                                                    
 1004 FORMAT('     Color grid interpolated to [FE/H] =',F6.2)           
 1005 FORMAT(1X,7HLOG G =,F6.3,11H LOG TEFF =,F7.4,8H OUTSIDE,          
     1   1X,11HCOLOR TABLE)                                             
 1006 FORMAT(A10,F5.2,A19)
 1007 FORMAT(1X,' **** INPUT DATA FILE DOES NOT EXIST **** ')      
      IF(MODE) 11,3,16                                                  
C *** WHEN MODE=0 THE COLOR TRANSFORMATION TABLES ARE READ
    3 TBND=LOG10(5500.)                                                
      GBND=2.
      OPEN(UNIT=ILO,ERR=49,FILE='uvbylo.data',STATUS='OLD')
      OPEN(UNIT=IHI,ERR=49,FILE='uvbyhi.data',STATUS='OLD')      
      READ(ILO,1000) NT,NG,NFE,NDX                                       
      READ(ILO,1001) (TA(I),I=1,NT)                                      
      READ(ILO,1002) (G(I),I=1,NG)                                       
      DO 4 I=1,NT                                                       
    4 T(I)=LOG10(TA(I))                                                
      DO 5 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                    
      DO 5 J=1,NG                                                       
      READ(ILO,1003) (BY(I,J,K),I=1,NT)                                  
    5 CONTINUE                                                          
      DO 6 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                 
      DO 6 J=1,NG                                                       
      READ(ILO,1003) (M1(I,J,K),I=1,NT)                                  
    6 CONTINUE                                                          
      DO 7 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                   
      DO 7 J=1,NG                                                       
      READ(ILO,1003) (C1(I,J,K),I=1,NT)                                  
    7 CONTINUE                                                          
      DO 8 K=1,NFE                                                      
      READ(ILO,1006) NAMFE,FEH(K),NAMCOL                                   
      DO 8 J=1,NG
      READ(ILO,1003) (BC(I,J,K),I=1,NT)
    8 CONTINUE                                             
      READ(IHI,1000) NTT,NGG,NFE,NDXX                                    
      READ(IHI,1001) (TB(I),I=1,NTT)                                     
      READ(IHI,1002) (GG(I),I=1,NGG)                                     
      DO 9 I=1,NTT                                                      
    9 TT(I)=LOG10(TB(I))                                               
      DO 10 L=1,NDXX                                                    
      DO 10 K=1,NFE                                                     
      READ(IHI,1006) NAMFE,FEH(K),NAMCOL
      DO 10 J=1,NGG                                                     
      READ(IHI,1003) (A(I,J,K,L),I=1,NTT)
   10 CONTINUE
      CLOSE(UNIT=ILO,STATUS='KEEP')
      CLOSE(UNIT=IHI,STATUS='KEEP')         
C *** WHEN MODE=0 OR MODE=-1, COLOR TRANSFORMATION TABLES ARE
C *** CREATED FOR THE INPUT [Fe/H] VALUE USING LINEAR INTERPOLATION
   11 DO 12 M=2,NFE                                                     
      K=M-1                                                             
      IF(FE.LE.FEH(M)) GO TO 13                                         
   12 CONTINUE                                                          
   13 M=K+1                                                             
      SLOPE=(FE-FEH(K))/(FEH(M)-FEH(K))                                 
      DO 14 J=1,NG                                                      
      DO 14 I=1,NT                                                      
      C(I,J)=BY(I,J,K)+SLOPE*(BY(I,J,M)-BY(I,J,K))
      D(I,J)=M1(I,J,K)+SLOPE*(M1(I,J,M)-M1(I,J,K))
      E(I,J)=C1(I,J,K)+SLOPE*(C1(I,J,M)-C1(I,J,K))
      F(I,J)=BC(I,J,K)+SLOPE*(BC(I,J,M)-BC(I,J,K))
   14 CONTINUE                                                          
      DO 15 L=1,NDXX                                                    
      DO 15 J=1,NGG                                                     
      DO 15 I=1,NTT                                                     
      B(I,J,L)=A(I,J,K,L)+SLOPE*(A(I,J,M,L)-A(I,J,K,L))                 
   15 CONTINUE                                                          
      !WRITE(0,1004) FE  
C *** WHEN MODE=0 OR MODE=-1, CONTROL RETURNS TO THE CALLING
C *** PROGRAM ONCE THE [Fe/H] INTERPOLATION IS CARRIED OUT (I.E., NO
C *** INTERPOLATIONS ARE PERFORMED FOR INPUT VALUES OF GV AND TEFF)
      GO TO 50
C *** WHEN MODE=1, INTERPOLATIONS ARE CARRIED OUT FOR THE INPUT VALUES
C *** OF [Fe/H], log g, AND log Teff (= FE, GV, AND TEFF IN THE ARGUMENT
C *** LIST)
   16 IF(TEFF.LE.TBND) GO TO 18                                         
      IF(GV.GE.GBND) GO TO 40                                           
      IF(TEFF.LE.T(NT)) GO TO 18                                        
C *** EXECUTION HALTS WITH A "STOP30" CODE IF THE INPUT TEMPERATURE IS
C *** OUTSIDE THE RANGE OF THE TABLES ON THE HIGH SIDE
C      WRITE(0,1005) GV,TEFF                                             
C      STOP30                                                            
      BCV=-99.999d0
      BMV=-99.999d0
      VMR=-99.999d0
      VMI=-99.999d0
      RETURN
C *** THE NEXT SECTION ASSUMES THAT THE LOW-TEMPERATURE TABLES ARE THE
C *** RELEVANT ONES TO USE
   18 NGM=NG-1                                                          
      DO 19 I=3,NGM                                                     
      MG=I-2                                                            
      IF(GV.LE.G(I)) GO TO 20                                           
   19 CONTINUE                                                          
      IF(GV.LE.G(NG)) GO TO 20                                          
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log g CONSIDERED IN THE 
C *** TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT TWO
C *** LINES SHOULD BE ACTIVATED   
C     WRITE(6,1005) GV,TEFF
C     STOP31
   20 R(1)=G(MG)                                                        
      R(2)=G(MG+1)                                                      
      R(3)=G(MG+2)                                                      
      R(4)=G(MG+3)                                                      
      CALL LGRAN4(R,AG,GV)                                              
      NTM=NT-1                                                          
      DO 25 I=3,NTM                                                     
      MT=I-2                                                            
      IF(TEFF.LE.T(I)) GO TO 30                                         
   25 CONTINUE                                                          
      IF(TEFF.LE.T(NT)) GO TO 30                                        
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log Teff CONSIDERED IN
C *** THE TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT
C *** TWO LINES SHOULD BE ACTIVATED.
C     WRITE(6,1005) GV,TEFF
C     STOP32
   30 R(1)=T(MT)                                                        
      R(2)=T(MT+1)                                                      
      R(3)=T(MT+2)                                                      
      R(4)=T(MT+3)                                                      
      CALL LGRAN4(R,AT,TEFF)                                            
      L=MG-1                                                            
      DO 35 I=1,4                                                       
      L=L+1                                                             
      O(I)=AT(1)*C(MT,L)+AT(2)*C(MT+1,L)+AT(3)*C(MT+2,L)+AT(4)*C(MT+3,L)
      P(I)=AT(1)*D(MT,L)+AT(2)*D(MT+1,L)+AT(3)*D(MT+2,L)+AT(4)*D(MT+3,L)
      Q(I)=AT(1)*E(MT,L)+AT(2)*E(MT+1,L)+AT(3)*E(MT+2,L)+AT(4)*E(MT+3,L)
   35 R(I)=AT(1)*F(MT,L)+AT(2)*F(MT+1,L)+AT(3)*F(MT+2,L)+AT(4)*F(MT+3,L)
C *** 4-POINT LAGRANGIAN INTERPOLATION IS USED TO DERIVE THE VALUES OF
C *** b-y, m1, c1, AND BC_V CORRESPONDING TO THE INPUT VALUES OF
C *** [Fe/H], log g, and log Teff.
      SBY=AG(1)*O(1)+AG(2)*O(2)+AG(3)*O(3)+AG(4)*O(4)
      SM1=AG(1)*P(1)+AG(2)*P(2)+AG(3)*P(3)+AG(4)*P(4)                   
      SC1=AG(1)*Q(1)+AG(2)*Q(2)+AG(3)*Q(3)+AG(4)*Q(4)                   
      BCV=AG(1)*R(1)+AG(2)*R(2)+AG(3)*R(3)+AG(4)*R(4)                  
C *** SPLINE INTERPOLATION IS USED TO FIND THE VALUE OF BC_V AT LOW
C *** TEMPERATURES.
      DO 38 I=1,NT                                                      
      X(I)=T(I)                                                         
   38 Y(I)=AG(1)*F(I,MG)+AG(2)*F(I,MG+1)+AG(3)*F(I,MG+2)+AG(4)*F(I,MG+3)
      CALL INTEP(TEFF,BCV,X,Y,NT,IER)                                  
      GO TO 50                                                          
C *** THE NEXT SECTION ASSUMES THAT THE HIGH-TEMPERATURE TABLES ARE THE
C *** RELEVANT ONES TO USE.
   40 IF(TEFF.LE.TT(NTT)) GO TO 42                                      
      WRITE(0,1005) GV,TEFF                                             
      STOP33                                                            
   42 NGM=NGG-1                                                         
      DO 43 I=3,NGM                                                     
      MG=I-2                                                            
      IF(GV.LE.GG(I)) GO TO 44                                          
   43 CONTINUE                                                          
      IF(GV.LE.GG(NGG)) GO TO 44                                        
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log g CONSIDERED IN THE
C *** TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT TWO
C *** LINES SHOULD BE ACTIVATED.
C     WRITE(6,1005) GV,TEFF   
C     STOP34
   44 R(1)=GG(MG)                                                       
      R(2)=GG(MG+1)                                                     
      R(3)=GG(MG+2)                                                     
      R(4)=GG(MG+3)                                                     
      CALL LGRAN4(R,AG,GV)                                              
      NTM=NTT-1                                                         
      DO 45 I=3,NTM                                                     
      MT=I-2                                                            
      IF(TEFF.LE.TT(I)) GO TO 46                                        
   45 CONTINUE                                                          
      IF(TEFF.LE.TT(NTT)) GO TO 46                                      
C *** SOME EXTRAPOLATION OUTSIDE THE RANGE OF log Teff CONSIDERED IN
C *** THE TABLES IS PERMITTED.  TO AVOID ANY EXTRAPOLATIONS, THE NEXT
C *** TWO LINES SHOULD SHOULD BE ACTIVATED.
C     WRITE(6,1005) GV,TEFF
C     STOP35
   46 R(1)=TT(MT)                                                       
      R(2)=TT(MT+1)                                                     
      R(3)=TT(MT+2)                                                     
      R(4)=TT(MT+3)                                                     
      CALL LGRAN4(R,AT,TEFF)                                            
      L=MG-1                                                            
      DO 48 I=1,4                                                       
      L=L+1                                                             
      DO 48 K=1,NDXX                                                    
      H(I,K)=AT(1)*B(MT,L,K)+AT(2)*B(MT+1,L,K)+AT(3)*B(MT+2,L,K)+       
     1   AT(4)*B(MT+3,L,K)                                              
   48 CONTINUE                                                          
C *** 4-POINT LAGRANGIAN INTERPOLATION IS USED TO DERIVE THE VALUES OF
C *** b-y, m1, c1, AND BC_V CORRESPONDING TO THE INPUT VALUES OF
C *** [Fe/H], log g, and log Teff.
      SBY=AG(1)*H(1,1)+AG(2)*H(2,1)+AG(3)*H(3,1)+AG(4)*H(4,1)           
      SM1=AG(1)*H(1,2)+AG(2)*H(2,2)+AG(3)*H(3,2)+AG(4)*H(4,2)           
      SC1=AG(1)*H(1,3)+AG(2)*H(2,3)+AG(3)*H(3,3)+AG(4)*H(4,3)           
      BCV=AG(1)*H(1,4)+AG(2)*H(2,4)+AG(3)*H(3,4)+AG(4)*H(4,4)
      GO TO 50
   49 WRITE(0,1007)
      STOP          
   50 RETURN                                                            
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE LGRAN4(X,A,XX)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      SAVE
      DIMENSION X(4),A(4)
      R1=(X(1)-X(2))*(X(1)-X(3))*(X(1)-X(4))
      R2=(X(2)-X(1))*(X(2)-X(3))*(X(2)-X(4))
      R3=(X(3)-X(1))*(X(3)-X(2))*(X(3)-X(4))
      R4=(X(4)-X(1))*(X(4)-X(2))*(X(4)-X(3))
      A(1)=((XX-X(2))*(XX-X(3))*(XX-X(4)))/R1
      A(2)=((XX-X(1))*(XX-X(3))*(XX-X(4)))/R2
      A(3)=((XX-X(1))*(XX-X(2))*(XX-X(4)))/R3
      A(4)=((XX-X(1))*(XX-X(2))*(XX-X(3)))/R4
      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE INTEP(XP,P,X,F,N,IER)
C *** PURPOSE:  To interpolate a function value P for a given argument
C     XP from a table of N values (X,F).  This is a spline interpolation
C     scheme based on Hermite polynomials.  The source is U.S. Airforce
C     Surveys in Geophysics No 272.
C *** USAGE:  For random values of XP
C               CALL INTEP(XP,P,X,F,N,IER)
C     or after the first call to INTEP with monotonically increasing or
C     decreasing values of XP consistent with the X vector
C               CALL EINTEP(XP,P,X,F,N,IER)
C     DESCRIPTION OF PARAMETERS:
C     XP  - the chosen argument value
C     P   - the resultant interpolated value
C     X   - the vector of independent values
C     F   - the vector of dependent values
C     N   - the number of points in the (X,F) vectors
C     IER - the resultant error parameter (set = 2 if XP is beyond
C           either extreme of X; in which case P is set to the value of
C           F at that extremum)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*8 LP1,LP2,L1,L2
      DIMENSION F(*),X(*)
      IER=1
      IO=1
      IUP=0
      IF(X(2).LT.X(1)) IUP=1
      N1=N-1
      IF((XP.GE.X(N).AND.IUP.EQ.0).OR.(XP.LE.X(N).AND.IUP.EQ.1)) THEN
*    5  P=F(N)
       P=F(N)
*       GO TO 6
       IER=2
       RETURN
      ELSE IF((XP.LE.X(1).AND.IUP.EQ.0).OR.
     1   (XP.GE.X(1).AND.IUP.EQ.1)) THEN
       P=F(1)
*    6  IER=2
       IER=2
       RETURN
      ENDIF
      ENTRY EINTEP(XP,P,X,F,N,IER)
    8 DO 1 I=IO,N
      IF(XP.LT.X(I).AND.IUP.EQ.0) GO TO 2
      IF(XP.GT.X(I).AND.IUP.EQ.1) GO TO 2
    1 CONTINUE
      P=F(N)
      IER=2
      RETURN
*      GO TO 5
    2 I=I-1
      IF(I.EQ.IO-1) GO TO 4
      IO=I+1
      LP1=1./(X(I)-X(I+1))
      LP2=1./(X(I+1)-X(I))
      IF(I.EQ.1) FP1=(F(2)-F(1))/(X(2)-X(1))
      IF(I.EQ.1) GO TO 3
      FP1=(F(I+1)-F(I-1))/(X(I+1)-X(I-1))
    3 IF(I.GE.N1) FP2=(F(N)-F(N-1))/(X(N)-X(N-1))
      IF(I.EQ.N1) GO TO 4
      FP2=(F(I+2)-F(I))/(X(I+2)-X(I))
    4 XPI1=XP-X(I+1)
      XPI=XP-X(I)
      L1=XPI1*LP1
      L2=XPI*LP2
      P=F(I)*(1.-2.*LP1*XPI)*L1*L1+F(I+1)*(1.-2.*LP2*XPI1)*L2*L2+
     1   FP2*XPI1*L2*L2+FP1*XPI*L1*L1
      RETURN
      END      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c END OF VANDENBERG & CLEM COLOR TRANSFORM SUBROUTINES
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c BEGIN AKIMA SPLINE SUBROUTINES
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE AKM_CSPL (n,x,f,c)

*  DESCRIPTION
*  Subroutine akm_cspl is based on the method of interpolation described
*  by Hiroshi Akima (1970 Journal of the Association for Computing Machinery,
*  Vol. 17, pp. 580-602) and also by Carl de Boor (1978 Applied Mathematical
*  Sciences, Vol. 27, "A Practical Guide to Splines").

*  This routine computes the Akima spline interpolant with the "not-a-knot
*  end conditions which are useful when there is no information about the
*  derivatives at the endpoints. In this case P(1)=P(2) and P(n-2)=P(n-1);
*  in other words, x(2) and x(n-1) are not knots even though they define
*  the endpoints of a polynomial piece --- the spline must pass through
*  these two points. Thus, it is a requirement that f''' be continuous
*  across x(2) and x(n-1). 

*  DIMENSIONS
*  The internal work arrays dd, w, and s are restricted to 2000 elements 
*  by the parameter nmx; consequently the data arrays x and f are also
*  limited to 2000 elements. The spline coefficients are contained in
*  the array c, a 4x2000 2 dimensional array.

*  CALL PARAMETERS
*     n      - number of data pairs
*     x(n)   - array, independent variable
*     f(n)   - array, dependent variable
*     c(4,n) - array, spline coefficients

      implicit double precision (a-h, o-z)
      parameter (nmx = 2000)
      parameter (zero = 0.0d0, two = 2.0d0, three = 3.0d0)
      dimension x(n),f(n),c(4,n) 
      dimension dd(nmx),w(nmx),s(nmx)

      nm1 = n - 1
      nm2 = n - 2
      nm3 = n - 3

*  Create two additional points outside the region of the spline by fitting
*  a quadratic to the three adjacent data points. This is a special feature
*  of the Akima spline.

      call akm_ends (x(1),x(2),x(3),f(1),f(2),f(3),f0,fm1)
      call akm_ends (x(n),x(nm1),x(nm2),f(n),f(nm1),f(nm2),fnp1,fnp2)

*  Compute the divided differences

      ddm1 = (f0 - fm1)/(x(2) - x(1))
      dd0 = (f(1) - f0)/(x(3) - x(2))

      do i = 1, nm1
        ip1 = i + 1
        dd(i) = (f(ip1) - f(i))/(x(ip1) - x(i))
      end do

      dd(n) = (f(n) - fnp1)/(x(nm2) - x(nm1))
      ddnp1 = (fnp1 - fnp2)/(x(nm1) - x(n))

*  Compute the Akima weights

      w0 = abs (dd0 - ddm1)
      w(1) = abs (dd(1) - dd0)

      do i = 2, nm1
        im1 = i - 1
        w(i) = abs(dd(i) - dd(im1))
      end do

      w(n) = abs (dd(n) - dd(nm1))
      wnp1 = abs (ddnp1 - dd(n))

*  Compute Akima slopes at interior knots

      if (w(2).eq.zero.and.w0.eq.zero) then
        s(1) = 5.0d-1*(dd(1) + dd0)
      else
        s(1) = (w(2)*dd0 + w0*dd(1))/(w0 + w(2))
      end if

      do i = 2, nm1
        im1 = i - 1
        ip1 = i + 1
        if (w(ip1).eq.zero.and.w(im1).eq.zero) then
          s(i) = 5.0d-1*(dd(i) + dd(im1))
        else
          s(i) = (w(ip1)*dd(im1) + w(im1)*dd(i))/(w(im1) + w(ip1))
        end if
      end do

      if (wnp1.eq.zero.and.w(nm1).eq.zero) then
        s(n) = 5.0d-1*(dd(n) + dd(nm1))
      else
        s(n) = (wnp1*dd(nm1) + w(nm1)*dd(n))/(w(nm1) + wnp1)
      end if

*  Spline coefficients for all the polynomial pieces

      do i = 1, nm1
        ip1 = i + 1
        dx = x(ip1) - x(i)
        c(1,i) = f(i)
        c(2,i) = s(i)
        c(3,i) = (three*dd(i) - two*s(i) - s(ip1))/dx
        c(4,i) = (s(ip1) + s(i) - two*dd(i))/dx/dx
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE AKM_ENDS(x1,x2,x3,f1,f2,f3,f00,f01)
      
*  DESCRIPTION
*  Subroutine akm_ends is based on the method of interpolation described
*  by Hiroshi Akima (1970 Journal of the Association for Computing Machinery,
*  Vol. 17, pp. 580-602).

*  This routine is required by the routine AKM_CSPL. It sets up the two 
*  additional external points required by the Akima method of constructing
*  a spline.

*  CALL PARAMETERS
*  (x1,f1), (x2,f2) and (x3,f3) are the known data pairs at the ends of the
*  region of interpolation. f00 and f01 are the computed values of the
*  interpolating polynomial external to the region of interpolation.

      implicit double precision (a-h,o-z)

      g0 = f1
      df21 = f2 - f1
      df31 = f3 - f1
      dx21 = x2 - x1
      dx31 = x3 - x1
      dx32 = x3 - x2
      den = dx21*dx32*dx31

      g1 = (df21*dx31*dx31 - df31*dx21*dx21)/den
      g2 = (df31*dx21 - df21*dx31)/den

      f00 = g0 - g1*dx32 + g2*dx32*dx32
      f01 = g0 - g1*dx31 + g2*dx31*dx31

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE AKM_EVAL (n,x,c,xp,fp,dfp)
      
*  DESCRIPTION
*  This routine is used to evaluate the spline defined by the
*  routine AKM_CSPL and its derivative at an interpolating point.

*  CALL PARAMETERS
*     n      - number of data pairs
*     x(n)   - array, independent variable
*     c(4,n) - array, spline coefficients
*     xp     - interpolating point
*     fp     - value of the spline at xp
*     dfp    - first derivative of the spline at xp   
      
      implicit double precision (a-h, o-z)
      dimension x(n), c(4,n)
      logical more

*  Statement function defines cubic spline interpolation

      csval(dx,a,b,e,d) = a +dx*(b + dx*(e + dx*d))

*  Statement function defines cubic spline derivative

      csder(dx,b,e,d) = b + dx*(2.0d0*e + dx*3.0d0*d)

*  Check direction of spline

      if (x(n).gt.x(1)) then

*  Check to see that x(1) < xp < x(n)

      if (xp.lt.x(1)) then
        write(*,'('' ***Test point is less than x(1)***'')')
        stop
      else if (xp.gt.x(n)) then
        write(*,'('' ***Test point is greater than x(n)***'')')
        stop
      end if

*  Find the polynomial piece containing the test point

      i = 2
      more = .true.
      do while (more)
        if (xp.lt.x(i).or.i.eq.n) then
          more = .false.
          nl = i - 1
        else 
          i = i + 1
        end if
      end do

      else

*  The spline is backwards
*  Check to see that x(n) < xp < x(1)

      if (xp.gt.x(1)) then
        write(*,'('' ***Test point is greater than x(1)***'')')
        stop
      else if (xp.lt.x(n)) then
        write(*,'('' ***Test point is less than x(n)***'')')
        stop
      end if

*  Find the polynomial piece containing the test point

      i = 2
      more = .true.
      do while (more)
        if (xp.gt.x(i).or.i.eq.n) then
          more = .false.
          nl = i - 1
        else 
          i = i + 1
        end if
      end do

      end if

*  Evaluate spline at the test point xp

      dx = xp - x(nl)
      fp = csval(dx,c(1,nl),c(2,nl),c(3,nl),c(4,nl))

*  Evaluate the derivative at the test point

      dfp = csder(dx,c(2,nl),c(3,nl),c(4,nl))
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   
      SUBROUTINE AKM_INT (a,b,n,x,c,sum)

*  DESCRIPTION
*  This routine is used to integrate the spline defined by the
*  routine AKM_CSPL.

*  CALL PARAMETERS
*     a      - lower bound
*     b      - upper bound
*     n      - number of data pairs
*     x(n)   - array, independent variable
*     c(4,n) - array, spline coefficients
*     sum     - value of the integral over the interval [a,b]

      implicit double precision (a-h, o-z)
      dimension x(n),c(4,n)

*  Statement function defined by cubic spline integration

      csint(dx,a,b,g,d)=abs(dx*(a+dx*(b/2.0d0+dx*(g/3.0d0+dx*d/4.0d0))))

*  Check direction of integration

      if (x(1).lt.x(n)) then

*  Forward integration
*  Check limits of integration

      if (a.lt.x(1)) then
        write(*,'('' ***Lower bound is less than x(1)***'')')
        stop
      else if (b.gt.x(n)) then
        write(*,'('' ***Upper bound is greater than x(n)***'')')
        stop
      end if

*  Find the polynomial piece containing the lower bound

      i = 2
      do while (a.gt.x(i))
        i = i + 1
      end do
      nl = i - 1

*  Find the polynomial piece containing the upper bound

      i = n - 1
      do while (b.lt.x(i))
        i = i - 1
      end do
      nu = i

      else

*  Backward integration
*  Check limits of integration

      if (a.gt.x(1)) then
        write(*,'('' ***Upper bound is greater than x(1)***'')')
        stop
      else if (b.lt.x(n)) then
        write(*,'('' ***Lower bound is less than x(n)***'')')
        stop
      end if

*  Find the polynomial piece containing the upper bound

      i = 2
      do while (a.lt.x(i))
        i = i + 1
      end do
      nl = i - 1

*  Find the polynomial piece containing the lower bound

      i = n - 1
      do while (b.gt.x(i))
        i = i - 1
      end do
      nu = i

      end if

*  Initialize the sum

      sum = 0.0d0

*  Subtract the portion of the 1st polynomial piece outside the lower bound

      dx = a - x(nl)
      if (dx.ne.0.0d0) then
        sum = sum - csint(dx,c(1,nl),c(2,nl),c(3,nl),c(4,nl))
      end if

*  Integrate from x(nl) to x(nu)

      if (nu.gt.nl) then
        nt = nu - 1
        do i = nl, nt
          dx = x(i+1) - x(i)
          sum = sum + csint(dx,c(1,i),c(2,i),c(3,i),c(4,i))
        end do
      end if

*  Add the portion of the last polynomial piece within the upper bound

      dx = b - x(nu)
      if (dx.ne.0.0d0) then
        sum = sum + csint(dx,c(1,nu),c(2,nu),c(3,nu),c(4,nu))
      end if

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE AKM_BISECTION(n,c,x,x1,x2,froot,xroot)

*  DESCRIPTION
*  This routine is used to invert the spline defined by the
*  routine AKM_CSPL in order to find a root by the bisection method.
*  The routine AKM_EVAL is called within this routine.

*  CALL PARAMETERS
*     n      - number of data pairs
*     c(4,n) - array, spline coefficients
*     x(n)   - array, independent variable
*     x1     - estimated lower limit of region to be searched
*     x2     - estimated upper limit of region to be searched
*     froot  - value of the spline at the location of the root
*     xroot  - calculated value of the root

      implicit double precision (a-h, o-z)
      dimension x(*), c(4,n)
      logical more

*  This routine must be preceded by a call to akm_cspl(n,x,f,c) to set up
*  the spline which is used to find the root.

      xa = x1
      xb = x2

      call akm_eval(n,x,c,xa,fa,dfdx)
      fa = froot - fa
      call akm_eval(n,x,c,xb,fb,dfdx)
      fb = froot - fb
      kl = 0

      do while (fa*fb.ge.0.0d0.and.kl.lt.20) 

        kl = kl + 1

        if (abs(fa).lt.abs(fb)) then
          xa = xa + 1.6d0*(xa - xb)
          xa = max(xa,x(1))
          call akm_eval(n,x,c,xa,fa,dfdx)
          fa = froot - fa
        else
          xb = xb + 1.6d0*(xb - xa)
          xb = min(xb,x(n))
          call akm_eval(n,x,c,xb,fb,dfdx)
          fb = froot - fb
        end if

       write(*,'(i4,1pd14.7,0pf7.3,1pd14.7,0pf7.3)') kl,xa,fa,xb,fb
      end do

      if (kl.eq.20) then
        write(*,'('' *** akm_bisection failed to find the region'')')
        write(*,'(''     containing the root after 20 iterations.'')')
        stop
      end if

      diff = abs(xb - xa)
      more = .true.
      kl = 0

      do while (more)

        kl = kl + 1
        xab = 5.0d-1*(xa + xb)
        call akm_eval(n,x,c,xab,fb,dfdx)
        fb = froot - fb

        if (fa*fb.le.0.0d0) then
          xb = xab
        else
          xa = xab
          fa = fb
        end if

      diff = abs(xb - xa)/abs(5.0d-1*(xa + xb))
      more = (diff.gt.1.0d-8.and.kl.lt.40)

      end do

      if (kl.eq.40) then
        write(*,'('' ***akm_bisection root did not converge within'')')
        write(*,'(''    one in 1.0d-8 after 40 iterations.'')')
        write(*,'(''    Converged within '',1pd12.5)') diff
      end if

      xroot = 5.0d-1*(xa + xb)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c END AKIMA SPLINE SUBROUTINES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c BEGIN UTILITY SUBROUTINES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE PARSE (FILESPEC,FIELD,OUTSPEC,LENGTH)
      CHARACTER FILESPEC*(*),FIELD*(*),OUTSPEC*(*)

      IF (FIELD.EQ.'DEVICE') THEN

        NC1 = 1
        NC2 = INDEX(filespec,':')
        OUTSPEC = FILESPEC(NC1:NC2)
        LENGTH = NC2 - NC1 + 1

      ELSE IF (FIELD.EQ.'DIRECTORY') THEN

        NC1 = index(filespec,'[')
        NC2 = index(filespec,']')
        OUTSPEC = FILESPEC(NC1:NC2)
        LENGTH = NC2 - NC1 + 1

      ELSE IF (FIELD.EQ.'NAME') THEN

        NC1 = index(filespec,']')
        NC2 = index(filespec,'.')
        OUTSPEC = FILESPEC(NC1+1:NC2-1)
        LENGTH = NC2 - NC1 - 1

      ELSE IF (FIELD.EQ.'TYPE') THEN

        NC1 = index(filespec,'.')
        NC2 = index(filespec,';')
        OUTSPEC = FILESPEC(NC1:NC2-1)
        LENGTH = NC2 - NC1

      ELSE IF (FIELD.EQ.'VERSION') THEN

        NC1 = index(filespec,';')
        NC2 = index(filespec,' ')
        OUTSPEC = FILESPEC(NC1+1:NC2-1)
        LENGTH = NC2 - NC1 - 1

      ELSE

        OUTSPEC = 'INVALID FILE SPECIFICATION'
        LENGTH = 26

      END IF

      RETURN
      END
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine str_trim(nl,str,nch)
      character*(*) str
      
      nch = nl
      do while (str(nch:nch).eq.' '.and.nch.gt.0)
        nch = nch - 1
      end do
      
      return
      end    
        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        SUBROUTINE OPEN_FILE(NUNIT,FSTAT,FNME,ASK,PROMPT)
        INTEGER NUNIT
        CHARACTER*(*) FSTAT,FNME,PROMPT
        CHARACTER*1 ANS
        LOGICAL ASK, MORE

        MORE = .TRUE.

        DO WHILE (MORE)
          IF (ASK) THEN
           WRITE(*,'(T25,A,$)') PROMPT
            READ(*,'(A)') FNME
          END IF
          
          NCH=1
          DO WHILE (FNME(NCH:NCH).NE.' ')
            NCH = NCH + 1
          END DO
          NCH = NCH - 1

          IF (FSTAT.EQ.'IN'.OR.FSTAT.EQ.'in') THEN
            OPEN (NUNIT,IOSTAT=IERR,FILE=FNME,STATUS='OLD')
            IF (IERR.EQ.0) THEN
              MORE=.FALSE.
            ELSE
              WRITE(*,*)'*** ',FNME(1:NCH),' DOES NOT EXIST***'
            END IF

          ELSE IF(FSTAT.EQ.'OUT'.OR.FSTAT.EQ.'out') THEN
            OPEN (NUNIT,IOSTAT=IERR,FILE=FNME,STATUS='NEW')
            IF (IERR.EQ.0) THEN
              MORE=.FALSE.
            ELSE
              WRITE(*,*)'***',FNME(1:NCH),' ALREADY EXISTS***'
              WRITE(*,'(T25,''DO YOU WANT TO OVERWRITE? '',$)')
              READ(*,'(A)') ANS
              IF (ANS.EQ.'Y'.OR.ANS.EQ.'y') THEN
                OPEN(NUNIT,FILE=FNME,STATUS='OLD')
                CLOSE(NUNIT,STATUS='DELETE')
                OPEN(NUNIT,FILE=FNME,STATUS='NEW')
                MORE = .FALSE.
              ELSE
                WRITE(*,'(T25,''New file name: '',$)')
                READ(*,'(A)') FNME
              END IF
            END IF
          END IF

       END DO

       RETURN
       END

*23456789012345678901234567890123456789012345678901234567890123456789012
        SUBROUTINE IOFILE(NUNIT,FSTAT,DNME,DASK,DPROM,FNME,ASK,PROM)
        INTEGER NUNIT
        CHARACTER*(*) FSTAT,FNME,PROM,DPROM,DNME
        CHARACTER*40 TNME
        CHARACTER*1 ANS
        LOGICAL ASK, MORE, DASK, TASK

        IF (DASK) THEN
          LP = LEN(DPROM)
          LD = LEN(DNME)
          WRITE(*,'(T25,A,A,A,A,$)') DPROM(1:LP),'[',DNME(1:LD),']: '
          READ(*,'(A)') TNME
          IF (TNME.EQ.' ') THEN
            FNME = DNME
          ELSE
            FNME = TNME
          END IF
        END IF
        
        MORE = .TRUE.
        TASK = ASK

        DO WHILE (MORE)
          IF (TASK) THEN
           WRITE(*,'(T25,A,$)') PROM
            READ(*,'(A)') FNME
          END IF
          
          NCH=1
          DO WHILE (FNME(NCH:NCH).NE.' ')
            NCH = NCH + 1
          END DO
          NCH = NCH - 1

          IF (FSTAT.EQ.'IN'.or.fstat.eq.'in') THEN
            OPEN (NUNIT,iostat=ierr,FILE=FNME,STATUS='OLD')
            if (ierr.eq.0) then
              more=.false.
            else
              WRITE(*,*)'***',FNME(1:NCH),' DOES NOT EXIST***'
              if (.not.ask) task = .true.
            end if

          ELSE IF(FSTAT.EQ.'OUT'.or.fstat.eq.'out') THEN
            OPEN (NUNIT,iostat=ierr,FILE=FNME,STATUS='NEW')
            if (ierr.eq.0) then
              more=.false.
            else
              WRITE(*,*)'***',FNME(1:NCH),' ALREADY EXISTS***'
              write(*,'(''       Do you want to OVERWRITE? '',$)')
              read(*,'(a)') ans
              if (ans.eq.'y'.or.ans.eq.'Y') then
                open(nunit,file=fnme,status='old')
                close(nunit,status='delete')
                open(nunit,file=fnme,status='new')
                more = .false.
              else         
                WRITE(*,'(T25,''New file name: '',$)')
                READ(*,'(A)') FNME
              end if
            end if

          ELSE
            OPEN(NUNIT,iostat=ierr,STATUS='SCRATCH')
            if (ierr.eq.0) then
              more=.false.
            else
              WRITE(*,*)'***SCRATCH FILE ALREADY EXISTS***'
              write(*,'(''       Do you want to OVERWRITE? '',$)')
              read(*,'(a)') ans
              if (ans.eq.'y'.or.ans.eq.'Y') then
                open(nunit,file=fnme,status='old')
                close(nunit,status='delete')
                open(nunit,file=fnme,status='new')
                more = .false.
              end if
            end if

          END IF

       end do

       return
       END

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      LOGICAL FUNCTION YES(QUESTION)

*  Prompts question on 'stdout', checks for valid answer on 'stdin'
*  Returns true if answer valid and affirmative. Else false.

      CHARACTER QUESTION*(*),LETTER*1
      LOGICAL MORE

      MORE = .TRUE.
      DO WHILE (MORE)
         WRITE(*,'(T25,A,'' [Y/N] '',$)') QUESTION
         READ (*,'(A)') LETTER
         IF (LETTER.EQ.'Y'.OR.LETTER.EQ.'y') THEN
            YES = .TRUE.
            MORE = .FALSE.
         ELSE IF (LETTER.EQ.'N'.OR.LETTER.EQ.'n') THEN
            YES = .FALSE.
            MORE = .FALSE.
         ELSE
            WRITE(*,'('' **** INVALID ANSWER **** '')')
         END IF
      END DO

      RETURN
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c END UTILITY SUBROUTINES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c $$$$$$$$$$$$$$$$$$$$$$ END VRSUB SUBROUTINES $$$$$$$$$$$$$$$$$$$$$$$$$$
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine trkcheck(alphafe,maxage, itrk, nbad,lbad)

      implicit none
      real*8 alphafe,maxage
      integer itrk,nbad
      logical lbad

      real*8 smass,xenv,zenv,afe2,cmixl,age
      integer endf

      lbad = .false.
      
      open(unit=itrk,status='old')

c check to make sure the *.iso file contains a good evolutionary track
c it must extend to at least MAXAGE
c ****** MC changes start
C         READ(itrk,*)
C         READ(itrk,*)
c         read(itrk,*)
c     ****** end MC changes
c         read(itrk,'(a80)')hd1
c         write(*,'(a80)')hd1
         READ(itrk,55)smass, XENV,ZENV,AFE2,CMIXL
c **** MC changed format
 55      FORMAT(4X,F5.3,3x,E10.4,3x,E10.4,18x,F5.2,6x,F6.4)
c     55      FORMAT(4X,F5.3,3x,E10.4,3x,E10.4,18x,F5.2,6x,F6.4)
c***  end of MC changes
!         write(*,*)'mass=',smass
         if(abs(afe2 - alphafe) > 0.001) then
            write(*,*)'CODESTOPPED Incorrect [alpha/Fe], itrk=',
     $           itrk, afe2,alphafe
            stop
         endif

         READ(itrk,*)
         read(itrk,*)age
         do while (age <= maxage*1.0D6)
            read(itrk,*,iostat=endf)age
            if(endf < 0) then
               write(*,*) 'input file', itrk, 'is bad.'
               close(itrk)
               nbad = nbad + 1
               lbad = .true.
               write(*,*)'itrk,nbad:',itrk,nbad
               return
            endif
         enddo
         rewind(itrk)
         return
         end
