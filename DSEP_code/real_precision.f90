!*******************************************************************************
!>
!  Numeric kinds.

    module real_precision

    use iso_fortran_env, only: real64,int64

    private

    integer,parameter,public :: dp = real64  !! default real kind
    integer,parameter,public :: ip = int64   !! default integer kind

    end module real_precision
!*******************************************************************************

