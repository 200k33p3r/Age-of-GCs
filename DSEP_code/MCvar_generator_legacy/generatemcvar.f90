
program GenerateMCvar
  !generate monte carlo input files for DSEP

  use real_precision
  use random
  
  implicit none
  integer, parameter :: numuniform =  9   !number of variables with a uniform distribution

  real(dp),parameter :: minbound = 0.01_dp  !minimum value for some gaussian variables

  
  integer i, kttau,mcstart,mcend
  real(dp) :: Yprim, alphafe, cmixla, fgry, fgrz, sstandard(7), alphae
  real(dp) :: alexcoef, opalcoef2, talphacoef, plascoef, cocoef,feh,ttau
  real(dp) :: fehmean,fehsig,afemean,afesig,alphac,dydz, sumz
  real(dp),dimension(numuniform) :: uniform

  real(dp) :: onethird, twothird

  character(len=7) fileprefix
  character(len=11) :: filestart
  character(len=17) :: outfile

  onethird = 1.0_dp/3.0_dp
  twothird = 2.0_dp/3.0_dp

  !read in composition, and the starting/ending MC numbers from the command line

  call get_variables(fehmean,fehsig,afemean,afesig,mcstart,mcend)  

!  write(*,'(A,4F9.3)') '[Fe/H] sigma_[Fe/H] [alpha/Fe], sigma_[alpha_Fe]', &
!       fehmean, fehsig,afemean,afesig
!  write(*,'(A,2I6)') 'MCstartnum, MCendnum:', mcstart,mcend

  if (fehmean < 0.0d0) then
     fileprefix = 'varfeh-'
  else
     fileprefix = 'varfeh+'
  endif
  write(filestart, '(A7,i3.3,A1)' ) fileprefix,nint(abs(fehmean*100._dp)),"."
  
!  write(*,*) 'filestart: ', filestart
  
! initialize random number generator
  call init_ranseed
      
  do i = mcstart,mcend
     call random_number(uniform)
    
!  generate values for variables with uniform distributions 
     FGRY = 0.5_dp + 0.8_dp*uniform(2)        !helium diffusion
     FGRZ = 0.5_dp + 0.8_dp*uniform(3)        !heavy element diffusion
     alexcoef = 0.7_dp + 0.6_dp*uniform(4)    !low temperature opacities
     ALPHAE =  0.2_dp*uniform(5)              !convetive envelope overshoot
     ALPHAC =  0.2_dp*uniform(8)              !convetive core overshoot
     CMIXLA = 1.0_dp + 1.5_dp*uniform(6)      !mixing length
     TTAU = uniform(7)                        !surface boundry condition
     if(TTAU <= onethird) then
        KTTAU = 0
     elseif(TTAU > onethird .and. TTAU <= twothird)  then
        KTTAU = 1
     else
        KTTAU = 5
     endif
!
! for [a/Fe], only have tables in steps of 0.2dex, so need to round to the 
! nearest tabulated value

     alphafe = afesig*random_normal() + afemean
!     write(*,*) alphafe
     if (alphafe < -0.10_dp ) then
        alphafe = -0.2_dp
     elseif (alphafe < 0.10_dp) then
        alphafe = 0.0_dp
     elseif (alphafe < 0.30_dp) then
        alphafe  = 0.20_dp
     elseif (alphafe < 0.50_dp ) then
        alphafe = 0.40_dp
     elseif (alphafe < 0.70_dp ) then
        alphafe = 0.60_dp
     else
        alphafe = 0.80_dp
     endif
!     write(*,*)'j,alphafe:',j,alphafe

     feh = fehsig*random_normal() + fehmean 

!  generate values for variables with gaussian distributions,
     !  enforcing a hard lower bound at minbound

     !high temperature opacities 
     opalcoef2 = 0.03_dp*random_normal() + 1.0_dp       
     do while (opalcoef2 < minbound )
        opalcoef2 = 0.03_dp*random_normal() + 1.0_dp
     enddo

     !triple-alpha nuclear reaction rate (He fusion)
     talphacoef = 0.15_dp*random_normal() + 1.0_dp
     do while (talphacoef < minbound )
        talphacoef = 0.15_dp*random_normal() + 1.0_dp
     enddo

     !plasma neutrino loses 
     plascoef = 0.05_dp*random_normal() + 1.0_dp
     do while (plascoef < minbound )
        plascoef = 0.05_dp*random_normal() + 1.0_dp
     enddo

     !conducitve opacities 
     cocoef = 0.20_dp*random_normal() + 1.0_dp
     do while (cocoef < minbound )
        cocoef = 0.2_dp*random_normal() + 1.0_dp
     enddo

     !pp->h2 + e +nu
     sstandard(1) = 0.0098_dp*random_normal() + 0.997543_dp
     do while (sstandard(1) < minbound)
        sstandard(1) = 0.0098_dp*random_normal() + 0.997543_dp
     enddo

     !he3+he3->he4+p+p
     sstandard(2) =0.0971_dp*random_normal() + 1.02913_dp
     do while (sstandard(2) < minbound)
        sstandard(2) =0.0971_dp*random_normal() + 1.02913_dp
     enddo

     !he3+he4->be7+g
     sstandard(3) =0.0556_dp*random_normal() + 1.03704_dp
     do while (sstandard(3) < minbound)
        sstandard(3) =0.0556_dp*random_normal() + 1.03704_dp   
     enddo

     !c12+p->n13+g
     sstandard(4) = 0.3448_dp*random_normal() + 0.965517_dp
     do while (sstandard(4) < minbound)
        sstandard(4) = 0.3448_dp*random_normal() + 0.965517_dp
     enddo

     !c13+p->n14+g
     sstandard(5) = 0.2182_dp*random_normal() + 1.47273_dp
     do while (sstandard(5) < minbound)
        sstandard(5) = 0.2182_dp*random_normal() + 1.47273_dp
     enddo

     !n14+p->o15+g
     sstandard(6) = 0.0331_dp*random_normal() + 0.478916_dp
     do while (sstandard(6) < minbound)
        sstandard(6) = 0.0331_dp*random_normal() + 0.478916_dp
     enddo

     !n16+p->f17+g
     sstandard(7) = 0.0851_dp*random_normal() + 1.12766_dp
     do while (sstandard(7) < minbound)
        sstandard(7) = 0.0851_dp*random_normal() + 1.12766_dp
     enddo

     call gs98getz(feh, alphafe, sumz)
     dydz = 1 + 1.500_dp*uniform(9)
     Yprim = 0.244_dp + 0.005_dp*uniform(1) + sumz*dydz

     write(outfile,'(A11,i0)' ) filestart, i
!     write(*,*)'outfile: ', outfile
     
     open(unit=24, file=outfile)

     write(24, 80) feh,   'FeH' 
     write(24, 80) Yprim, 'Yprim'
     write(24, 80) alphafe, 'alphafe'
     write(24, 80) CMIXLA, 'CMIXLA'
     write(24, 80) FGRY, 'FGRY'
     write(24, 80) FGRZ, 'FGRZ'
     write(24, 90) KTTAU, 'KTTAU'
     write(24, 80) ALPHAE, 'ALPHAE'
     write(24, 80) ALPHAC, 'ALPHAC'
     write(24, 80) SSTANDARD(1), 'SSTANDARD(1) -- PP'
     write(24, 80) SSTANDARD(2), 'SSTANDARD(2) -- He3+He3'
     write(24, 80) SSTANDARD(3), 'SSTANDARD(3) -- He3+He4'
     write(24, 80) SSTANDARD(4), 'SSTANDARD(4) -- P+C12'
     write(24, 80) SSTANDARD(5), 'SSTANDARD(5) -- P+C13'
     write(24, 80) SSTANDARD(6), 'SSTANDARD(6) -- P+N14'
     write(24, 80) SSTANDARD(7), 'SSTANDARD(7) -- P+O16'
     write(24, 80) alexcoef, 'alexcoef'
     write(24, 80) opalcoef2, 'opalcoef2'
     write(24, 80) talphacoef, 'talphacoef'
     write(24, 80) plascoef, 'plascoef'
     write(24, 80) cocoef, 'cocoef'
80   format(F9.6, 5X, 30A)
90   format(I2, 12X, 30A)

     close(24)


  enddo
  stop
end program GenerateMCvar

subroutine get_variables(feh,fehsig,afe,afesig,mcstart,mcend)
 
  use real_precision

  implicit none

  real(dp) :: feh,fehsig,afe,afesig
  integer :: mcstart,mcend

  
!  character(len=30) :: filename
  character(len=10) :: cfeh,cfehsig,cafe,cafesig, cmcstart,cmcend
  integer :: n_args

  
  !read in Sample Size  from command line
  n_args = command_argument_count()
  if (n_args /= 6) then
     write(*,*) &
    'Usage: ./a.out [Fe/H] sigma_[Fe/H] [alpha/Fe]  sigma_[alpha/Fe] startMCnum endMCnum '
     stop
  end if

  call get_command_argument(1,cfeh)
  cfeh = trim(cfeh)
  read(cfeh,*) feh
  call get_command_argument(2,cfehsig)
  cfehsig = trim(cfehsig)
  read(cfehsig,*) fehsig
  call get_command_argument(3,cafe)
  cafe = trim(cafe)
  read(cafe,*) afe
  call get_command_argument(4,cafesig)
  cafesig = trim(cafesig)
  read(cafesig,*) afesig
  call get_command_argument(5,cmcstart)
  cmcstart = trim(cmcstart)
  read(cmcstart,*) mcstart
  call get_command_argument(6,cmcend)
  cmcend = trim(cmcend)
  read(cmcend,*) mcend

  if(mcend < mcstart) then
     write(*,*)'ERROR: MC ending number must be larger than MC starting number'
     write(*,*)'MCstart, MCend:', mcstart,mcend
     stop
  endif



  return
end subroutine get_variables


	

