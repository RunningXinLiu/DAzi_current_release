program depthkernel_ti_worker
  implicit none

  integer, parameter :: max_nx = 20000
  integer, parameter :: max_ny = 20000
  integer, parameter :: max_nz = 200
  integer, parameter :: max_k = 160
  integer :: nx, ny, nz, kmax, iwave, igr
  integer :: col_start, col_count, col_end
  integer :: argc, ios, unit_id
  integer :: post, local_col, ii, jj
  real :: minthk
  real, allocatable :: vel(:,:,:), depz(:)
  real(kind=8), allocatable :: periods(:), pv_part(:,:)
  real(kind=4), allocatable :: lsen_part(:,:,:)
  character(len=1024) :: vel_bin, depz_bin, periods_bin, outdir
  character(len=64) :: arg

  argc = command_argument_count()
  if (argc /= 13) then
    write(*,'(a)') 'usage: depthkernel_ti_worker VEL_BIN DEPZ_BIN PERIODS_BIN'
    write(*,'(a)') '       NX NY NZ KMAX MINTHK IWAVE IGR COL_START COL_COUNT OUTDIR'
    stop 2
  endif

  call get_command_argument(1, vel_bin)
  call get_command_argument(2, depz_bin)
  call get_command_argument(3, periods_bin)
  call get_command_argument(4, arg); read(arg,*) nx
  call get_command_argument(5, arg); read(arg,*) ny
  call get_command_argument(6, arg); read(arg,*) nz
  call get_command_argument(7, arg); read(arg,*) kmax
  call get_command_argument(8, arg); read(arg,*) minthk
  call get_command_argument(9, arg); read(arg,*) iwave
  call get_command_argument(10, arg); read(arg,*) igr
  call get_command_argument(11, arg); read(arg,*) col_start
  call get_command_argument(12, arg); read(arg,*) col_count
  call get_command_argument(13, outdir)

  if (nx <= 0 .or. ny <= 0 .or. nz <= 1 .or. kmax <= 0) stop 'Invalid dimensions.'
  if (nz > max_nz .or. kmax > max_k) stop 'Worker compiled limits exceeded.'
  if (col_start < 1 .or. col_count < 1 .or. col_start + col_count - 1 > nx*ny) stop 'Invalid column range.'
  if (nx > max_nx .or. ny > max_ny) stop 'Invalid horizontal dimensions.'

  allocate(vel(nx,ny,nz), depz(nz), periods(kmax), pv_part(col_count,kmax), &
      lsen_part(col_count,kmax,nz-1))

  open(newunit=unit_id, file=trim(vel_bin), form='unformatted', access='stream', &
      status='old', action='read', iostat=ios)
  if (ios /= 0) stop 'Cannot open vel bin.'
  read(unit_id) vel
  close(unit_id)

  open(newunit=unit_id, file=trim(depz_bin), form='unformatted', access='stream', &
      status='old', action='read', iostat=ios)
  if (ios /= 0) stop 'Cannot open depz bin.'
  read(unit_id) depz
  close(unit_id)

  open(newunit=unit_id, file=trim(periods_bin), form='unformatted', access='stream', &
      status='old', action='read', iostat=ios)
  if (ios /= 0) stop 'Cannot open periods bin.'
  read(unit_id) periods
  close(unit_id)

  pv_part = 0.0d0
  lsen_part = 0.0
  col_end = col_start + col_count - 1
  do post = col_start, col_end
    local_col = post - col_start + 1
    ii = mod(post - 1, nx) + 1
    jj = (post - 1) / nx + 1
    call depthkernel_ti_column(nz, vel(ii,jj,1:nz), kmax, periods, iwave, igr, depz, minthk, &
        pv_part(local_col,1:kmax), lsen_part(local_col,1:kmax,1:nz-1))
  enddo

  call system('mkdir -p ' // trim(outdir))
  open(newunit=unit_id, file=trim(outdir)//'/pv_part.bin', &
      form='unformatted', access='stream', status='replace')
  write(unit_id) pv_part
  close(unit_id)
  open(newunit=unit_id, file=trim(outdir)//'/lsen_gsc_part.bin', &
      form='unformatted', access='stream', status='replace')
  write(unit_id) lsen_part
  close(unit_id)

  deallocate(vel, depz, periods, pv_part, lsen_part)
end program depthkernel_ti_worker

subroutine depthkernel_ti_column(nz, vsz, kmaxRc, tRc, iwave, igr, depz, minthk, pv_col, lsen_col)
  implicit none
  integer, parameter :: NL = 200
  integer, parameter :: NP=160
  integer, intent(in) :: nz, kmaxRc, iwave, igr
  real, intent(in) :: vsz(nz), depz(nz), minthk
  real(kind=8), intent(in) :: tRc(kmaxRc)
  real(kind=8), intent(out) :: pv_col(kmaxRc)
  real(kind=4), intent(out) :: lsen_col(kmaxRc,nz-1)

  real :: vpz(nz), rhoz(nz)
  integer :: mmax, iflsph, mode, rmax
  integer :: i, j, k, jjj
  real(kind=8) :: cgRc(NP)
  real :: rdep(NL), rvp(NL), rvs(NL), rrho(NL), rthk(NL)
  real :: t_in(kmaxRc), cp_in(kmaxRc)
  real :: TA_in(NL), TC_in(NL), TF_in(NL)
  real :: TL_in(NL), TN_in(NL), TRho_in(NL)
  real :: qp(NL), qs(NL), etap(NL)
  real :: etas(NL), frefp(NL), frefs(NL)
  real(kind=4) :: dcdah(NP,NL), dcdn(NP,NL), dcdbv(NP,NL)
  real(kind=4) :: dcR_dL, dcR_dA
  integer :: nsublay(NL)

  if (kmaxRc > NP .or. nz > NL) stop 'depthkernel_ti_column limits exceeded.'
  mmax = nz
  iflsph = 1
  mode = 1
  lsen_col = 0.0

  do k = 1, nz
    vpz(k) = 0.9409 + 2.0947*vsz(k) - 0.8206*vsz(k)**2 + &
        0.2683*vsz(k)**3 - 0.0251*vsz(k)**4
    rhoz(k) = 1.6612*vpz(k) - 0.4721*vpz(k)**2 + &
        0.0671*vpz(k)**3 - 0.0043*vpz(k)**4 + &
        0.000106*vpz(k)**5
  enddo

  call refineLayerMdl_worker(minthk, mmax, depz, vpz, vsz, rhoz, rmax, rdep, rvp, rvs, rrho, rthk, nsublay)

  call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, kmaxRc, tRc, cgRc)
  pv_col(1:kmaxRc) = cgRc(1:kmaxRc)

  do i = 1, rmax
    TA_in(i) = rrho(i)*rvp(i)**2
    TC_in(i) = TA_in(i)
    TL_in(i) = rrho(i)*rvs(i)**2
    TN_in(i) = TL_in(i)
    TF_in(i) = 1.0*(TA_in(i) - 2 * TL_in(i))
    TRho_in(i) = rrho(i)
  enddo
  qp(1:rmax) = 150.0
  qs(1:rmax) = 50.0
  etap(1:rmax) = 0.00
  etas(1:rmax) = 0.00
  frefp(1:rmax) = 1.00
  frefs(1:rmax) = 1.00

  cp_in(1:kmaxRc) = sngl(cgRc(1:kmaxRc))
  t_in(1:kmaxRc) = sngl(tRc(1:kmaxRc))

  call tregn96(rmax, rthk, TA_in, TC_in, TF_in, TL_in, TN_in, TRho_in, &
      qp, qs, etap, etas, frefp, frefs, &
      kmaxRc, t_in, cp_in(1:kmaxRc), dcdah, dcdbv, dcdn)

  do i = 1, kmaxRc
    k = 0
    do j = 1, nz-1
      do jjj = 1, nsublay(j)
        k = k + 1
        dcR_dA = 0.5/(rrho(k)*rvp(k))*dcdah(i,k) - TF_in(k)/((TA_in(k)-2.0*TL_in(k))**2)*dcdn(i,k)
        dcR_dL = 0.5/(rrho(k)*rvs(k))*dcdbv(i,k) + 2.0*TF_in(k)/((TA_in(k)-2.0*TL_in(k))**2)*dcdn(i,k)
        lsen_col(i,j) = lsen_col(i,j) + dcR_dA*TA_in(k) + dcR_dL*TL_in(k)
      enddo
    enddo
  enddo
end subroutine depthkernel_ti_column

subroutine refineLayerMdl_worker(minthk0, mmax, dep, vp, vs, rho, rmax, rdep, rvp, rvs, rrho, rthk, nsublay)
  implicit none
  integer, parameter :: NL = 200
  integer, intent(in) :: mmax
  integer, intent(out) :: rmax
  real, intent(in) :: minthk0
  real, intent(in) :: dep(*), vp(*), vs(*), rho(*)
  real, intent(out) :: rdep(NL), rvp(NL), rvs(NL), rrho(NL), rthk(NL)
  integer, intent(out) :: nsublay(NL)
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
end subroutine refineLayerMdl_worker
