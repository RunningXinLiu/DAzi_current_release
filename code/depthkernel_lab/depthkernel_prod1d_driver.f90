program depthkernel_prod1d_driver
  implicit none
  interface
    subroutine depthkernel_column_phase_fd(nz, vs, kmax, periods, iwave, igr, depz, minthk, pv)
      integer, intent(in) :: nz, kmax, iwave, igr
      real, intent(in) :: vs(nz), depz(nz), minthk
      real(kind=8), intent(in) :: periods(*)
      real(kind=8), intent(out) :: pv(*)
    end subroutine depthkernel_column_phase_fd

    subroutine depthkernel_column_depth_fd(nz, vs, idep, kmax, periods, iwave, igr, depz, minthk, &
        dlnvs, dlnvp, dlnrho, sen_vs, sen_vp, sen_rho)
      integer, intent(in) :: nz, idep, kmax, iwave, igr
      real, intent(in) :: vs(nz), depz(nz), minthk
      real, intent(in) :: dlnvs, dlnvp, dlnrho
      real(kind=8), intent(in) :: periods(*)
      real(kind=8), intent(out) :: sen_vs(*), sen_vp(*), sen_rho(*)
    end subroutine depthkernel_column_depth_fd

    subroutine depthkernel_column_fd(nz, vs, kmax, periods, iwave, igr, depz, minthk, &
        dlnvs, dlnvp, dlnrho, pv, sen_vs, sen_vp, sen_rho)
      integer, intent(in) :: nz, kmax, iwave, igr
      real, intent(in) :: vs(nz), depz(nz), minthk
      real, intent(in) :: dlnvs, dlnvp, dlnrho
      real(kind=8), intent(in) :: periods(*)
      real(kind=8), intent(out) :: pv(*), sen_vs(kmax,*), sen_vp(kmax,*), sen_rho(kmax,*)
    end subroutine depthkernel_column_fd
  end interface

  integer, parameter :: max_nz = 200, max_k = 160
  integer :: nz, kmax, iwave, igr, ios
  real :: depz(max_nz), vs(max_nz), minthk
  real(kind=8) :: periods(max_k), pv(max_k)
  real(kind=8) :: sen_vs(max_k,max_nz), sen_vp(max_k,max_nz), sen_rho(max_k,max_nz)
  character(len=256) :: infile, outfile
  character(len=256) :: arg
  integer :: argc, unit_in, unit_out
  real :: t0, t1
  real :: dlnvs, dlnvp, dlnrho

  argc = command_argument_count()
  if (argc >= 1) then
    call get_command_argument(1, infile)
  else
    infile = 'sample_column.in'
  endif
  if (argc >= 2) then
    call get_command_argument(2, outfile)
  else
    outfile = 'prod1d_kernel.out'
  endif

  unit_in = 10
  open(unit_in, file=trim(infile), status='old', action='read', iostat=ios)
  if (ios /= 0) stop 'Cannot open input file.'

  call skip_comments(unit_in)
  read(unit_in,*) nz, kmax
  if (nz > max_nz .or. kmax > max_k) stop 'Input exceeds compiled limits.'

  call skip_comments(unit_in)
  read(unit_in,*) depz(1:nz)
  call skip_comments(unit_in)
  read(unit_in,*) vs(1:nz)
  call skip_comments(unit_in)
  read(unit_in,*) periods(1:kmax)
  call skip_comments(unit_in)
  read(unit_in,*) minthk, iwave, igr
  close(unit_in)

  dlnvs = 0.01
  dlnvp = 0.01
  dlnrho = 0.01
  if (argc >= 3) then
    call get_command_argument(3, arg)
    read(arg,*) dlnvs
    dlnvp = dlnvs
    dlnrho = dlnvs
  endif
  if (argc >= 4) then
    call get_command_argument(4, arg)
    read(arg,*) minthk
  endif
  if (argc >= 5) then
    call get_command_argument(5, arg)
    read(arg,*) dlnvp
  endif
  if (argc >= 6) then
    call get_command_argument(6, arg)
    read(arg,*) dlnrho
  endif

  call cpu_time(t0)
  call depthkernel_column_fd(nz, vs, kmax, periods, iwave, igr, depz, minthk, &
      dlnvs, dlnvp, dlnrho, pv, sen_vs, sen_vp, sen_rho)
  call cpu_time(t1)

  unit_out = 20
  open(unit_out, file=trim(outfile), status='replace', action='write')
  write(unit_out,'(a)') '# depthkernel_prod1d_driver output'
  write(unit_out,'(a,i0,1x,a,i0,1x,a,i0,1x,a,i0)') '# nz ', nz, 'kmax ', kmax, 'iwave ', iwave, 'igr ', igr
  write(unit_out,'(a,4f12.6)') '# minthk_dlnvs_dlnvp_dlnrho ', minthk, dlnvs, dlnvp, dlnrho
  write(unit_out,'(a,f12.5)') '# cpu_seconds ', t1 - t0
  write(unit_out,'(a)') '# periods_and_phase_or_group_velocity'
  write(unit_out,'(a)') '# period_s velocity_km_s'
  call write_period_table(unit_out, kmax, periods, pv)
  write(unit_out,'(a)') '# kernels: depth_index depth_km period_s dC_dVs dC_dVp dC_drho'
  call write_kernel_table(unit_out, nz, kmax, depz, periods, sen_vs, sen_vp, sen_rho)
  close(unit_out)

  write(*,'(a)') 'Wrote '//trim(outfile)
  write(*,'(a,f10.4)') 'CPU seconds: ', t1 - t0

contains

  subroutine skip_comments(unit_id)
    integer, intent(in) :: unit_id
    character(len=512) :: line
    integer :: ios_local
    do
      read(unit_id,'(a)',iostat=ios_local) line
      if (ios_local /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#') cycle
      backspace(unit_id)
      exit
    enddo
  end subroutine skip_comments

  subroutine write_period_table(unit_id, nk, t, c)
    integer, intent(in) :: unit_id, nk
    real(kind=8), intent(in) :: t(nk), c(nk)
    integer :: i
    do i = 1, nk
      write(unit_id,'(2f16.8)') t(i), c(i)
    enddo
  end subroutine write_period_table

  subroutine write_kernel_table(unit_id, nz_local, nk, dep, t, kvs, kvp, krho)
    integer, intent(in) :: unit_id, nz_local, nk
    real, intent(in) :: dep(nz_local)
    real(kind=8), intent(in) :: t(nk), kvs(nk,nz_local), kvp(nk,nz_local), krho(nk,nz_local)
    integer :: iz_local, ik
    do iz_local = 1, nz_local
      do ik = 1, nk
        write(unit_id,'(i5,2f16.8,3es18.8)') iz_local, dep(iz_local), t(ik), &
            kvs(ik,iz_local), kvp(ik,iz_local), krho(ik,iz_local)
      enddo
    enddo
  end subroutine write_kernel_table

end program depthkernel_prod1d_driver

subroutine refineGrid2LayerMdl(minthk0, mmax, dep, vp, vs, rho, &
    rmax, rdep, rvp, rvs, rrho, rthk)
  implicit none
  integer, parameter :: NL = 200
  integer, intent(in) :: mmax
  integer, intent(out) :: rmax
  real, intent(in) :: minthk0
  real, intent(in) :: dep(*), vp(*), vs(*), rho(*)
  real, intent(out) :: rdep(NL), rvp(NL), rvs(NL), rrho(NL), rthk(NL)
  integer :: nsublay(NL)
  real :: thk, minthk, newthk, initdep
  integer :: i, j, k

  k = 0
  initdep = 0.0
  do i = 1, mmax-1
    thk = dep(i+1) - dep(i)
    minthk = thk/minthk0
    nsublay(i) = int((thk + 1.0e-4)/minthk) + 1
    newthk = thk/nsublay(i)
    do j = 1, nsublay(i)
      k = k + 1
      if (k > NL) stop 'Too many refined layers.'
      rthk(k) = newthk
      rdep(k) = initdep + rthk(k)
      initdep = rdep(k)
      rvp(k) = vp(i) + (2*j - 1)*(vp(i+1) - vp(i))/(2*nsublay(i))
      rvs(k) = vs(i) + (2*j - 1)*(vs(i+1) - vs(i))/(2*nsublay(i))
      rrho(k) = rho(i) + (2*j - 1)*(rho(i+1) - rho(i))/(2*nsublay(i))
    enddo
  enddo

  k = k + 1
  if (k > NL) stop 'Too many refined layers.'
  rthk(k) = 0.0
  rvp(k) = vp(mmax)
  rvs(k) = vs(mmax)
  rrho(k) = rho(mmax)
  rdep(k) = dep(mmax)
  rmax = k
end subroutine refineGrid2LayerMdl
