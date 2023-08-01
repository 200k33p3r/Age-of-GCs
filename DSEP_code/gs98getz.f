cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                    gs98getz.f
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This program uses the GS98 Space Sci. Rev. 85, 161 mixture 
c input [Fe/H] and [a/Fe] on the command line and output total Z 
c for stellar evolution; 
c due to diffusion and that we don't know X, this is not exact
c


      subroutine gs98getz(feh,alpha,sumz)
      implicit none
      real*8 awt(23), a(23), b(23), b1(23),sum1,cln,sum,sumam,sumz
      real*8 feh,alpha,erra(23),errb(23),errb1(23),errsum,errsumz
      integer i
      character*2 cname(23)  
      character*20 cfeh,cafe
      data awt/ 1.0079, 4.0026, 12.0111, 14.0067, 15.9994, 20.179,
     1 22.98977, 24.305, 26.98154, 28.0855, 30.97376, 32.06, 35.453,
     * 39.948,40.08,
     2 47.9, 51.996, 54.938, 55.847, 58.70,6.939,9.0122,10.811 /
      data cname/'H','He','C','N','O','Ne','Na','Mg','Al','Si','P',
     $           'S','Cl','Ar','Ca','Ti','Cr','Mn','Fe','Ni',
     $           'Li','Be','B'/
C Element identification:
C 1=H; 2=He; 3=C; 4=N; 5=O; 6=Ne; 7=NA; 8=Mg; 9=Al; 10=Si; 
C 11=P; 12=S; 13=Cl; 14=Ar; 15=Ca; 16=Ti; 17=Cr; 18=Mn; 19=Fe; 20=Ni
C 21=Li; 22=Be; 23=B
C GS98 abundances
C allow for possiblilty of alpha enhancement

      CALL SETUP(feh,alpha,a,erra)
      cln = log(10.0D0)
      sum = 0.0D0
      errsum = 0.0D0
      do i=1, 23
	 b(i) = exp(cln*a(i))
         errb(i) = b(i)*cln*erra(i)
	 b1(i) = b(i)
	 sum = sum + b(i)
         errsum = errsum +errb(i)**2 
      end do
      errsum = sqrt(errsum)
      b1(1) = 0.0
      b1(2) = 0.0
      do i=3, 23
	 sum1 = sum1 + b1(i)
      end do
      sumam = 0.0
      do i=1, 23
	 b(i) = b(i) / sum
	 sumam = sumam + b(i)*awt(i)
	 b1(i) = b1(i) / sum1
         errb1(i) = errb(i)/sum1
         errb(i) = errb(i)/sum
      end do
      do i=1, 23
	 a(i) = b(i) * awt(i)/ sumam
         erra(i) = errb(i)*awt(i)/sumam
      end do
      sumz=0.0d0
      errsumz=0.0d0
      do i=3, 23
	 sumz = sumz + a(i) 
         errsumz = errsumz + erra(i)**2
      end do
      errsumz = sqrt(errsumz)
      write(*,140)sumz
  140 format(F12.9)
      end subroutine gs98getz
C
C
      subroutine setup(feh,alpha,a,erra)
      implicit none
      real*8 a(23), erra(23), feh,alpha
C put in logarithm form numbers from
C GS98
C allow for possiblilty of alpha enhancement
C Element identification:
C 1=H; 2=He; 3=C; 4=N; 5=O; 6=Ne; 7=NA; 8=Mg; 9=Al; 10=Si; 
C 11=P; 12=S; 13=Cl; 14=Ar; 15=Ca; 16=Ti; 17=Cr; 18=Mn; 19=Fe; 20=Ni
C 21=Li; 22=Be; 23=B
      a(1) = 12.0
      a(2) = 10.93
      a(3) = 8.52 + feh            
      a(4) = 7.92 + feh            
      a(5) = 8.83 + alpha + feh    
      a(6) = 8.08 + alpha + feh
      a(7) = 6.32 + feh
      a(8) = 7.58 + alpha + feh
      a(9) = 6.49 + feh
      a(10) = 7.56 + alpha + feh
      a(11) = 5.50 + feh
      a(12) = 7.22 + alpha + feh
      a(13) = 5.28 + feh
      a(14) = 6.40 + feh
      a(15) = 6.35 + alpha + feh
      a(16) = 4.98 + alpha + feh
      a(17) = 5.69 + feh
      a(18) = 5.51 - alpha + feh
      a(19) = 7.50 + feh 
      a(20) = 6.25 + feh
C add in li, be and B
      a(21) = 3.31
      a(22) = 1.42
      a(23) = 2.79  
C put in errors
      erra(1) = 0.0
      erra(2) = 0.0
      erra(3) = 0.05
      erra(4) = 0.07
      erra(5) = 0.07
      erra(6) = 0.06
      erra(7) = 0.03
      erra(8) = 0.05
      erra(9) = 0.07
      erra(10) = 0.05
      erra(11) = 0.04
      erra(12) = 0.06
      erra(13) = 0.3
      erra(14) = 0.14
      erra(15) = 0.03
      erra(16) = 0.02
      erra(17) = 0.03
      erra(18) = 0.03
      erra(19) = 0.01
      erra(20) = 0.04
      erra(21) = 0.1
      erra(22) = 0.1
      erra(23) = 0.04
      return
      end subroutine setup
