program surfdisp96_forward_batch
  implicit none

  integer, parameter :: max_k = 128
  integer :: nmodels, nlay, kmax, iflsph, iwave, igr, mode
  integer :: imodel, iper, unit_in, unit_out, ios
  real, allocatable :: rthk(:), rvp(:), rvs(:), rrho(:)
  real(kind=8) :: periods(max_k), curve(max_k)
  character(len=512) :: infile, outfile

  if (command_argument_count() >= 1) then
    call get_command_argument(1, infile)
  else
    infile = 'surfdisp96_forward_batch.in'
  endif
  if (command_argument_count() >= 2) then
    call get_command_argument(2, outfile)
  else
    outfile = 'surfdisp96_forward_batch.out'
  endif

  unit_in = 10
  open(unit_in, file=trim(infile), status='old', action='read', iostat=ios)
  if (ios /= 0) stop 'Cannot open input file.'

  call skip_comments(unit_in)
  read(unit_in,*) nmodels, nlay, kmax, iflsph, iwave, igr
  if (kmax > max_k) stop 'kmax exceeds compiled max_k.'

  allocate(rthk(nlay), rvp(nlay), rvs(nlay), rrho(nlay))

  call skip_comments(unit_in)
  read(unit_in,*) periods(1:kmax)

  unit_out = 20
  open(unit_out, file=trim(outfile), status='replace', action='write')
  write(unit_out,'(a)') '# model_index period_index period_s velocity_km_s'

  mode = 1
  do imodel = 1, nmodels
    call skip_comments(unit_in)
    read(unit_in,*) rthk(1:nlay)
    call skip_comments(unit_in)
    read(unit_in,*) rvp(1:nlay)
    call skip_comments(unit_in)
    read(unit_in,*) rvs(1:nlay)
    call skip_comments(unit_in)
    read(unit_in,*) rrho(1:nlay)

    call surfdisp96(rthk, rvp, rvs, rrho, nlay, iflsph, iwave, mode, igr, kmax, periods, curve)
    do iper = 1, kmax
      write(unit_out,'(i8,1x,i8,1x,f16.8,1x,f18.10)') imodel, iper, periods(iper), curve(iper)
    enddo
  enddo

  close(unit_out)
  close(unit_in)
  deallocate(rthk, rvp, rvs, rrho)
  write(*,'(a)') 'Wrote '//trim(outfile)

contains

  subroutine skip_comments(unit_id)
    integer, intent(in) :: unit_id
    character(len=4096) :: line
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

end program surfdisp96_forward_batch
