subroutine scipy_lsmr_solve(tag, m, n, leniw, lenrw, iw, rw, b, damp, &
                            atol, btol, conlim, itnlim, x, istop, itn, &
                            anorm, acond, rnorm, arnorm, xnorm)
    implicit none
    character(len=*), intent(in) :: tag
    integer, intent(in) :: m, n, itnlim
    integer(kind=8), intent(in) :: leniw, lenrw
    integer, intent(in) :: iw(*)
    real, intent(in) :: rw(*), b(*)
    real, intent(in) :: damp, atol, btol, conlim
    real, intent(out) :: x(*)
    integer, intent(out) :: istop, itn
    real, intent(out) :: anorm, acond, rnorm, arnorm, xnorm

    character(len=2048) :: python_cmd, script_path, base_dir, outdir
    character(len=2048) :: row_file, col_file, rw_file, b_file
    character(len=4096) :: cmd
    character(len=1024) :: env_value
    character(len=64) :: m_s, n_s, nar_s, damp_s, atol_s, btol_s
    character(len=64) :: conlim_s, maxiter_s
    integer :: status, unit_id, cmd_status, nnz_i
    integer :: meta_i(2)
    real(kind=8) :: meta_r(5)
    logical :: file_exists, log_open, keep_binaries

    if (leniw < 2_8 * lenrw) then
        write(*,'(A)') 'ERROR: scipy_lsmr_solve received leniw < 2*lenrw'
        stop 9301
    endif
    if (lenrw > int(huge(nnz_i), kind=8)) then
        write(*,'(A)') 'ERROR: scipy_lsmr_solve nnz exceeds INTEGER(4) capacity'
        stop 9302
    endif
    nnz_i = int(lenrw, kind=4)

    python_cmd = '/opt/miniconda3/envs/seisloc/bin/python'
    script_path = '/Users/liuxin/Desktop/DSurf_test/DAzi_vorotomo_handoff_20260515/current_release_20260515/code/depthkernel_lab/scipy_lsmr_triplet_solver.py'
    base_dir = '.dazi_lsmr'

    env_value = ''
    call get_environment_variable('DAZI_PYTHON', env_value, status=status)
    if (status == 0) python_cmd = trim(adjustl(env_value))
    env_value = ''
    call get_environment_variable('DAZI_LSMR_SCRIPT', env_value, status=status)
    if (status == 0) script_path = trim(adjustl(env_value))
    env_value = ''
    call get_environment_variable('DAZI_LSMR_DIR', env_value, status=status)
    if (status == 0) base_dir = trim(adjustl(env_value))
    keep_binaries = .false.
    env_value = ''
    call get_environment_variable('DAZI_LSMR_KEEP_BINARIES', env_value, status=status)
    if (status == 0) then
        env_value = trim(adjustl(env_value))
        keep_binaries = env_value == '1' .or. env_value == 'true' .or. &
                        env_value == 'TRUE' .or. env_value == 'yes' .or. env_value == 'YES'
    endif

    outdir = trim(adjustl(base_dir)) // '/' // trim(adjustl(tag))
    call system('mkdir -p ' // trim(adjustl(outdir)))

    row_file = trim(adjustl(outdir)) // '/row.bin'
    col_file = trim(adjustl(outdir)) // '/col.bin'
    rw_file = trim(adjustl(outdir)) // '/rw.bin'
    b_file = trim(adjustl(outdir)) // '/b.bin'

    open(newunit=unit_id, file=trim(row_file), form='unformatted', &
         access='stream', status='replace')
    write(unit_id) iw(1:nnz_i)
    close(unit_id)

    open(newunit=unit_id, file=trim(col_file), form='unformatted', &
         access='stream', status='replace')
    write(unit_id) iw(nnz_i + 1:2 * nnz_i)
    close(unit_id)

    open(newunit=unit_id, file=trim(rw_file), form='unformatted', &
         access='stream', status='replace')
    write(unit_id) rw(1:nnz_i)
    close(unit_id)

    open(newunit=unit_id, file=trim(b_file), form='unformatted', &
         access='stream', status='replace')
    write(unit_id) b(1:m)
    close(unit_id)

    write(m_s,'(I0)') m
    write(n_s,'(I0)') n
    write(nar_s,'(I0)') nnz_i
    write(damp_s,'(ES24.16E3)') damp
    write(atol_s,'(ES24.16E3)') atol
    write(btol_s,'(ES24.16E3)') btol
    write(conlim_s,'(ES24.16E3)') conlim
    write(maxiter_s,'(I0)') itnlim

    cmd = trim(adjustl(python_cmd)) // ' ' // trim(adjustl(script_path)) // &
          ' --row-bin ' // trim(row_file) // &
          ' --col-bin ' // trim(col_file) // &
          ' --rw-bin ' // trim(rw_file) // &
          ' --b-bin ' // trim(b_file) // &
          ' --outdir ' // trim(adjustl(outdir)) // &
          ' --m ' // trim(adjustl(m_s)) // &
          ' --n ' // trim(adjustl(n_s)) // &
          ' --nar ' // trim(adjustl(nar_s)) // &
          ' --damp ' // trim(adjustl(damp_s)) // &
          ' --atol ' // trim(adjustl(atol_s)) // &
          ' --btol ' // trim(adjustl(btol_s)) // &
          ' --conlim ' // trim(adjustl(conlim_s)) // &
          ' --maxiter ' // trim(adjustl(maxiter_s))

    inquire(unit=66, opened=log_open)
    write(*,'(A,A)') 'LSMR backend: scipy, tag=', trim(adjustl(tag))
    write(*,'(A,A)') 'LSMR backend command: ', trim(cmd)
    if (log_open) then
        write(66,'(A,A)') 'LSMR backend: scipy, tag=', trim(adjustl(tag))
        write(66,'(A,A)') 'LSMR backend command: ', trim(cmd)
    endif

    call execute_command_line(trim(cmd), exitstat=cmd_status)
    if (cmd_status /= 0) then
        write(*,'(A,I0)') 'ERROR: SciPy LSMR backend failed with exit status ', cmd_status
        if (log_open) write(66,'(A,I0)') 'ERROR: SciPy LSMR backend failed with exit status ', cmd_status
        stop 9303
    endif

    inquire(file=trim(adjustl(outdir)) // '/x_lsmr.bin', exist=file_exists)
    if (.not. file_exists) then
        write(*,'(A)') 'ERROR: missing SciPy LSMR output x_lsmr.bin'
        if (log_open) write(66,'(A)') 'ERROR: missing SciPy LSMR output x_lsmr.bin'
        stop 9304
    endif
    open(newunit=unit_id, file=trim(adjustl(outdir)) // '/x_lsmr.bin', &
         form='unformatted', access='stream', status='old')
    read(unit_id) x(1:n)
    close(unit_id)

    inquire(file=trim(adjustl(outdir)) // '/meta_i.bin', exist=file_exists)
    if (.not. file_exists) then
        write(*,'(A)') 'ERROR: missing SciPy LSMR output meta_i.bin'
        if (log_open) write(66,'(A)') 'ERROR: missing SciPy LSMR output meta_i.bin'
        stop 9305
    endif
    open(newunit=unit_id, file=trim(adjustl(outdir)) // '/meta_i.bin', &
         form='unformatted', access='stream', status='old')
    read(unit_id) meta_i
    close(unit_id)

    inquire(file=trim(adjustl(outdir)) // '/meta_r.bin', exist=file_exists)
    if (.not. file_exists) then
        write(*,'(A)') 'ERROR: missing SciPy LSMR output meta_r.bin'
        if (log_open) write(66,'(A)') 'ERROR: missing SciPy LSMR output meta_r.bin'
        stop 9306
    endif
    open(newunit=unit_id, file=trim(adjustl(outdir)) // '/meta_r.bin', &
         form='unformatted', access='stream', status='old')
    read(unit_id) meta_r
    close(unit_id)

    istop = meta_i(1)
    itn = meta_i(2)
    anorm = real(meta_r(1))
    acond = real(meta_r(2))
    rnorm = real(meta_r(3))
    arnorm = real(meta_r(4))
    xnorm = real(meta_r(5))
    if (.not. keep_binaries) then
        call execute_command_line('rm -f ' // trim(row_file) // ' ' // trim(col_file) // ' ' // &
                                  trim(rw_file) // ' ' // trim(b_file))
    endif
end subroutine scipy_lsmr_solve
