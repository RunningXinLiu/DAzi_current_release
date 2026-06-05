subroutine caldespersion(nx, ny, nz, vel, pv, iwave, igr, kmax, periods, depz, minthk)
    implicit none
    integer, intent(in) :: nx, ny, nz, iwave, igr, kmax
    real, intent(in) :: vel(nx,ny,nz), depz(nz), minthk
    real*8, intent(in) :: periods(*)
    real*8, intent(out) :: pv(nx*ny,*)

    stop 'caldespersion is not enabled in active parallel FMM ISO backend'
end subroutine caldespersion
