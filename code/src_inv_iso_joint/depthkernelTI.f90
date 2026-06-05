! use tregn96 from cps to calculate the dcdL and dcdA
subroutine depthkernelTI(nx, ny, nz, vel, pvRc, iwave, igr, kmaxRc, tRc, depz, minthk, Lsen_Gsc)
    use omp_lib
    implicit none

    integer nx, ny, nz
    real vel(nx, ny, nz)
    integer iwave,igr
    real minthk
    real depz(nz)
    integer kmaxRc
    real*8 tRc(kmaxRc)
    ! output
    real*8 pvRc(nx*ny, kmaxRc)
    ! parameter list
    real vpz(nz),vsz(nz),rhoz(nz)
    integer mmax,iflsph,mode,rmax
    integer ii,jj,k,i,j,jjj
    integer, parameter:: NL=200
    integer, parameter:: NP=160

    real*8 cgRc(NP)
    real rdep(NL),rvp(NL),rvs(NL),rrho(NL),rthk(NL)
    ! for tregn96
    real t_in(kmaxRc), cp_in(kmaxRc)
    real TA_in(NL), TC_in(NL), TF_in(NL)
    real TL_in(NL), TN_in(NL), TRho_in(NL)
    real qp(NL), qs(NL), etap(NL)
    real etas(NL), frefp(NL), frefs(NL)

    real*4 dcdah(NP,NL),dcdn(NP,NL)
    real*4 dcdbv(NP,NL)
    real*4 dcR_dL, dcR_dA
    real*4 Lsen_Gsc(nx*ny, kmaxRc, nz-1)
    integer nsublay(NL), post
    integer rmax_min, rmax_max, rmax_sum
    integer ti_threads_observed
    real*8 timer_refine_phase, timer_tregn, timer_splice
    real*8 timer_tmp_start, timer_tmp_end
    logical log_open
    integer env_status
    character(len=128) backend_name
    character(len=1024) dump_dir

    mmax=nz
    iflsph=1
    mode=1
    pvRc=0.0
    Lsen_Gsc=0.0
    rmax_min=NL
    rmax_max=0
    rmax_sum=0
    ti_threads_observed=1
    timer_refine_phase=0.0d0
    timer_tregn=0.0d0
    timer_splice=0.0d0
    backend_name=''
    call get_environment_variable('DAZI_DEPTHKERNEL_TI_BACKEND', backend_name)
    if(len_trim(backend_name) == 0) backend_name='default_fortran_tregn96'
    inquire(unit=66, opened=log_open)
    write(6,'(a,a)') 'Depth kernel TI backend: ', trim(backend_name)
    write(6,'(a,i0,a,i0,a,i0,a,i0)') 'Depth kernel TI grid/periods: nx=', nx, ' ny=', ny, &
        ' nz=', nz, ' kmax=', kmaxRc
    write(6,'(a,i0)') 'Depth kernel TI OpenMP max threads: ', omp_get_max_threads()
	    if(log_open) then
	        write(66,'(a,a)') 'Depth kernel TI backend: ', trim(backend_name)
	        write(66,'(a,i0,a,i0,a,i0,a,i0)') 'Depth kernel TI grid/periods: nx=', nx, ' ny=', ny, &
	            ' nz=', nz, ' kmax=', kmaxRc
	        write(66,'(a,i0)') 'Depth kernel TI OpenMP max threads: ', omp_get_max_threads()
	    endif
	    if(trim(adjustl(backend_name)) == 'process_batch') then
	        call depthkernelTI_process_batch(nx,ny,nz,vel,pvRc,iwave,igr,kmaxRc,tRc,depz,minthk,Lsen_Gsc)
	        dump_dir = ''
	        call get_environment_variable('DAZI_DEPTHKERNEL_TI_DUMP_DIR', dump_dir, status=env_status)
	        if(env_status == 0) then
	            if(log_open) write(66,'(a,a)') 'Depth kernel TI dump dir: ', trim(adjustl(dump_dir))
	            call depthkernelTI_dump_arrays(nx,ny,nz,kmaxRc,pvRc,Lsen_Gsc,trim(adjustl(dump_dir)))
	        endif
	        return
	    else if(trim(adjustl(backend_name)) /= 'default_fortran_tregn96') then
	        write(6,'(a,a)') 'ERROR: unknown Depth kernel TI backend: ', trim(adjustl(backend_name))
	        if(log_open) write(66,'(a,a)') 'ERROR: unknown Depth kernel TI backend: ', trim(adjustl(backend_name))
	        stop 9201
	    endif
	    ! write(6, *) ' depth kernel parallel:'
    !$omp parallel &
    !$omp default(private) &
    !$omp shared(nx,ny,nz,vel,minthk,mmax,depz,kmaxRc) &
    !$omp shared(tRc,pvRc,iflsph,iwave,mode,igr,Lsen_Gsc) &
    !$omp reduction(+:timer_refine_phase,timer_tregn,timer_splice,rmax_sum) &
    !$omp reduction(min:rmax_min) reduction(max:rmax_max,ti_threads_observed)
    !$omp do
    do jj=1,ny
        do ii=1,nx
            ti_threads_observed = max(ti_threads_observed, omp_get_num_threads())
            timer_tmp_start = omp_get_wtime()
            post=ii+(jj-1)*nx
            vsz(1:nz)=vel(ii,jj,1:nz)
            ! some other emperical relationship maybe better,
            do k=1,nz
                vpz(k)=0.9409 + 2.0947*vsz(k) - 0.8206*vsz(k)**2+ &
                0.2683*vsz(k)**3 - 0.0251*vsz(k)**4
                rhoz(k)=1.6612*vpz(k) - 0.4721*vpz(k)**2 + &
                0.0671*vpz(k)**3 - 0.0043*vpz(k)**4 + &
                0.000106*vpz(k)**5
            enddo
            ! change from refineGrid2LayerMdl into refineLayerMdl
            ! call refineGrid2LayerMdl(minthk, mmax, depz, vpz, vsz, rhoz, rmax, rdep, &
            ! rvp, rvs, rrho, rthk)
            call refineLayerMdl(minthk, mmax, depz, vpz, vsz, rhoz, rmax, rdep, &
            rvp, rvs, rrho, rthk, nsublay)
            rmax_min = min(rmax_min, rmax)
            rmax_max = max(rmax_max, rmax)
            rmax_sum = rmax_sum + rmax

            call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, kmaxRc, &
            tRc, cgRc)
            pvRc(ii+(jj-1)*nx,1:kmaxRc)=cgRc(1:kmaxRc)
            timer_tmp_end = omp_get_wtime()
            timer_refine_phase = timer_refine_phase + timer_tmp_end - timer_tmp_start
            timer_tmp_start = timer_tmp_end
            !print*,cgRc(1:kmaxRc)
            !------------------------------------------------------------------!
            do i = 1, rmax
                TA_in(i)=rrho(i)*rvp(i)**2
                TC_in(i)=TA_in(i)
                TL_in(i)=rrho(i)*rvs(i)**2
                TN_in(i)=TL_in(i)
                TF_in(i)=1.0*(TA_in(i) - 2 * TL_in(i))
                TRho_in(i)=rrho(i)
            enddo
            qp(1:rmax)=150.0
            qs(1:rmax)=50.0
            etap(1:rmax)=0.00
            etas(1:rmax)=0.00
            frefp(1:rmax)=1.00
            frefs(1:rmax)=1.00

            cp_in(1:kmaxRc)=sngl(cgRc(1:kmaxRc))
            t_in(1:kmaxRc)=sngl(tRc(1:kmaxRc))

            ! ! write(6, *)'tregn96'
            call tregn96(rmax, rthk, TA_in, TC_in, TF_in, TL_in, TN_in, TRho_in, &
            qp, qs, etap, etas, frefp, frefs,  &
            kmaxRc, t_in, cp_in(1:kmaxRc),&
            dcdah, dcdbv, dcdn)
            timer_tmp_end = omp_get_wtime()
            timer_tregn = timer_tregn + timer_tmp_end - timer_tmp_start
            timer_tmp_start = timer_tmp_end
            !
            ! ! write(*,*)"nsublay:", nsublay(1:nz)
            do i=1,kmaxRc  ! period
                k=0
                do j=1,nz-1                ! inversion layer
                    do jjj=1,nsublay(j)    ! refined layer k-th in jth inversion layer
                        k=k+1
                        dcR_dA = 0.5/(rrho(k)*rvp(k))*dcdah(i, k) - TF_in(k)/((TA_in(k)-2.0*TL_in(k))**2)*dcdn(i,k)
                        dcR_dL = 0.5/(rrho(k)*rvs(k))*dcdbv(i, k) + 2.0*TF_in(k)/((TA_in(k)-2.0*TL_in(k))**2)*dcdn(i,k)
                        Lsen_Gsc(post,i,j)=Lsen_Gsc(post,i,j)+dcR_dA*TA_in(k)+dcR_dL*TL_in(k)
                    enddo
                enddo
            enddo
            timer_tmp_end = omp_get_wtime()
            timer_splice = timer_splice + timer_tmp_end - timer_tmp_start

        enddo
    enddo
    !$omp end do
    !$omp end parallel
    write(6,'(a)') '  depthkernelTI internal timing:'
    write(6,'(a,f12.3,a)') '    refine+surfdisp96 phase velocity = ', timer_refine_phase, ' s'
    write(6,'(a,f12.3,a)') '    tregn96 sensitivity kernel = ', timer_tregn, ' s'
    write(6,'(a,f12.3,a)') '    splice refined kernel = ', timer_splice, ' s'
    write(6,'(a,i0,a,i0,a,f10.2)') '    refined sublayers rmax min/max/avg = ', rmax_min, '/', &
        rmax_max, '/', dble(rmax_sum)/dble(nx*ny)
    write(6,'(a,i0)') '    observed worker threads = ', ti_threads_observed
    inquire(unit=66, opened=log_open)
    if(log_open) then
        write(66,'(a)') '  depthkernelTI internal timing:'
        write(66,'(a,f12.3,a)') '    refine+surfdisp96 phase velocity = ', timer_refine_phase, ' s'
        write(66,'(a,f12.3,a)') '    tregn96 sensitivity kernel = ', timer_tregn, ' s'
        write(66,'(a,f12.3,a)') '    splice refined kernel = ', timer_splice, ' s'
        write(66,'(a,i0,a,i0,a,f10.2)') '    refined sublayers rmax min/max/avg = ', rmax_min, '/', &
            rmax_max, '/', dble(rmax_sum)/dble(nx*ny)
        write(66,'(a,i0)') '    observed worker threads = ', ti_threads_observed
    endif
    dump_dir = ''
    call get_environment_variable('DAZI_DEPTHKERNEL_TI_DUMP_DIR', dump_dir, status=env_status)
    if(env_status == 0) then
        if(log_open) write(66,'(a,a)') 'Depth kernel TI dump dir: ', trim(adjustl(dump_dir))
        call depthkernelTI_dump_arrays(nx,ny,nz,kmaxRc,pvRc,Lsen_Gsc,trim(adjustl(dump_dir)))
    endif
end subroutine

subroutine depthkernelTI_dump_arrays(nx,ny,nz,kmaxRc,pvRc,Lsen_Gsc,dump_dir)
    implicit none
    integer nx, ny, nz, kmaxRc
    real*8 pvRc(nx*ny, kmaxRc)
    real*4 Lsen_Gsc(nx*ny, kmaxRc, nz-1)
    character(len=*) dump_dir
    integer unit_id

    if(len_trim(dump_dir) == 0) return
    call system('mkdir -p ' // trim(adjustl(dump_dir)))
    open(newunit=unit_id,file=trim(adjustl(dump_dir))//'/pv_ti.bin', &
        form='unformatted',access='stream',status='replace')
    write(unit_id) pvRc
    close(unit_id)
    open(newunit=unit_id,file=trim(adjustl(dump_dir))//'/lsen_gsc_ti.bin', &
        form='unformatted',access='stream',status='replace')
    write(unit_id) Lsen_Gsc
    close(unit_id)
end subroutine depthkernelTI_dump_arrays

subroutine depthkernelTI_process_batch(nx,ny,nz,vel,pvRc,iwave,igr,kmaxRc,tRc,depz,minthk,Lsen_Gsc)
    implicit none
    integer nx, ny, nz, iwave, igr, kmaxRc
    real vel(nx,ny,nz), depz(nz), minthk
    real*8 tRc(kmaxRc), pvRc(nx*ny,kmaxRc)
    real*4 Lsen_Gsc(nx*ny,kmaxRc,nz-1)
    character(len=1024) python_cmd, script_path, worker_path, outdir, cmd, env_value
    character(len=32) nx_s, ny_s, nz_s, kmax_s, minthk_s, iwave_s, igr_s
    character(len=32) workers_s, tile_s
    integer status, cmd_status, unit_id, workers, tile_columns
    logical log_open, file_exists

    python_cmd = '/opt/miniconda3/envs/dispa/bin/python'
    script_path = '/Users/liuxin/Desktop/DSurf_test/DAzi_vorotomo_migration/depthkernel_lab/depthkernel_ti_process_batch.py'
    worker_path = '/Users/liuxin/Desktop/DSurf_test/DAzi_vorotomo_migration/depthkernel_lab/depthkernel_ti_worker'
    outdir = '.dazi_ti_batch'
    workers = 4
    tile_columns = 64

    env_value = ''
    call get_environment_variable('DAZI_PYTHON', env_value, status=status)
    if(status == 0 .and. len_trim(env_value) > 0) python_cmd = env_value
    env_value = ''
    call get_environment_variable('DAZI_DEPTHKERNEL_TI_BATCH_SCRIPT', env_value, status=status)
    if(status == 0 .and. len_trim(env_value) > 0) script_path = env_value
    env_value = ''
    call get_environment_variable('DAZI_DEPTHKERNEL_TI_WORKER', env_value, status=status)
    if(status == 0 .and. len_trim(env_value) > 0) worker_path = env_value
    env_value = ''
    call get_environment_variable('DAZI_DEPTHKERNEL_TI_BATCH_DIR', env_value, status=status)
    if(status == 0 .and. len_trim(env_value) > 0) outdir = env_value
    env_value = ''
    call get_environment_variable('DAZI_DEPTHKERNEL_TI_WORKERS', env_value, status=status)
    if(status == 0 .and. len_trim(env_value) > 0) read(env_value,*) workers
    env_value = ''
    call get_environment_variable('DAZI_DEPTHKERNEL_TI_TILE_COLUMNS', env_value, status=status)
    if(status == 0 .and. len_trim(env_value) > 0) read(env_value,*) tile_columns
    if(workers < 1) workers = 1
    if(tile_columns < 1) tile_columns = 1

    call system('mkdir -p ' // trim(adjustl(outdir)))
    open(newunit=unit_id,file=trim(adjustl(outdir))//'/vel_ti.bin', &
        form='unformatted',access='stream',status='replace')
    write(unit_id) vel
    close(unit_id)
    open(newunit=unit_id,file=trim(adjustl(outdir))//'/depz.bin', &
        form='unformatted',access='stream',status='replace')
    write(unit_id) depz
    close(unit_id)
    open(newunit=unit_id,file=trim(adjustl(outdir))//'/tRc.bin', &
        form='unformatted',access='stream',status='replace')
    write(unit_id) tRc
    close(unit_id)

    write(nx_s,'(I0)') nx
    write(ny_s,'(I0)') ny
    write(nz_s,'(I0)') nz
    write(kmax_s,'(I0)') kmaxRc
    write(minthk_s,'(F12.6)') minthk
    write(iwave_s,'(I0)') iwave
    write(igr_s,'(I0)') igr
    write(workers_s,'(I0)') workers
    write(tile_s,'(I0)') tile_columns

    cmd = trim(adjustl(python_cmd))//' '//trim(adjustl(script_path))// &
        ' --vel-bin '//trim(adjustl(outdir))//'/vel_ti.bin'// &
        ' --depz-bin '//trim(adjustl(outdir))//'/depz.bin'// &
        ' --periods-bin '//trim(adjustl(outdir))//'/tRc.bin'// &
        ' --nx '//trim(adjustl(nx_s))// &
        ' --ny '//trim(adjustl(ny_s))// &
        ' --nz '//trim(adjustl(nz_s))// &
        ' --kmax '//trim(adjustl(kmax_s))// &
        ' --minthk '//trim(adjustl(minthk_s))// &
        ' --iwave '//trim(adjustl(iwave_s))// &
        ' --igr '//trim(adjustl(igr_s))// &
        ' --workers '//trim(adjustl(workers_s))// &
        ' --tile-columns '//trim(adjustl(tile_s))// &
        ' --worker '//trim(adjustl(worker_path))// &
        ' --outdir '//trim(adjustl(outdir))

    inquire(unit=66, opened=log_open)
    write(6,'(a,a)') 'Depth kernel TI batch Python: ', trim(adjustl(python_cmd))
    write(6,'(a,a)') 'Depth kernel TI batch script: ', trim(adjustl(script_path))
    write(6,'(a,a)') 'Depth kernel TI worker: ', trim(adjustl(worker_path))
    write(6,'(a,a)') 'Depth kernel TI batch outdir: ', trim(adjustl(outdir))
    write(6,'(a,i0)') 'Depth kernel TI batch workers: ', workers
    write(6,'(a,i0)') 'Depth kernel TI batch tile columns: ', tile_columns
    write(6,'(a,a)') 'Depth kernel TI backend command: ', trim(cmd)
    if(log_open) then
        write(66,'(a,a)') 'Depth kernel TI batch Python: ', trim(adjustl(python_cmd))
        write(66,'(a,a)') 'Depth kernel TI batch script: ', trim(adjustl(script_path))
        write(66,'(a,a)') 'Depth kernel TI worker: ', trim(adjustl(worker_path))
        write(66,'(a,a)') 'Depth kernel TI batch outdir: ', trim(adjustl(outdir))
        write(66,'(a,i0)') 'Depth kernel TI batch workers: ', workers
        write(66,'(a,i0)') 'Depth kernel TI batch tile columns: ', tile_columns
        write(66,'(a,a)') 'Depth kernel TI backend command: ', trim(cmd)
    endif

    call execute_command_line(trim(cmd), exitstat=cmd_status)
    if(cmd_status /= 0) then
        write(6,'(a,i0)') 'ERROR: depth kernel TI backend command failed with exit status ', cmd_status
        if(log_open) write(66,'(a,i0)') 'ERROR: depth kernel TI backend command failed with exit status ', cmd_status
        stop 9202
    endif

    inquire(file=trim(adjustl(outdir))//'/pv_ti_batch.bin', exist=file_exists)
    if(.not. file_exists) then
        write(6,'(a)') 'ERROR: missing TI batch output pv_ti_batch.bin'
        if(log_open) write(66,'(a)') 'ERROR: missing TI batch output pv_ti_batch.bin'
        stop 9203
    endif
    inquire(file=trim(adjustl(outdir))//'/lsen_gsc_ti_batch.bin', exist=file_exists)
    if(.not. file_exists) then
        write(6,'(a)') 'ERROR: missing TI batch output lsen_gsc_ti_batch.bin'
        if(log_open) write(66,'(a)') 'ERROR: missing TI batch output lsen_gsc_ti_batch.bin'
        stop 9204
    endif

    open(newunit=unit_id,file=trim(adjustl(outdir))//'/pv_ti_batch.bin', &
        form='unformatted',access='stream',status='old',action='read')
    read(unit_id) pvRc
    close(unit_id)
    open(newunit=unit_id,file=trim(adjustl(outdir))//'/lsen_gsc_ti_batch.bin', &
        form='unformatted',access='stream',status='old',action='read')
    read(unit_id) Lsen_Gsc
    close(unit_id)
end subroutine depthkernelTI_process_batch

! do i=1, kmaxRc
! 	do j=1, rmax
! 		dcR_dA(i, j)= 0.5/(rrho(j)*rvp(j))*dcdah(i, j) - TF_in(j)/((TA_in(j)-2.0*TL_in(j))**2)*dcdn(i,j)
! 		dcR_dL(i, j)= 0.5/(rrho(j)*rvs(j))*dcdbv(i, j) + 2.0*TF_in(j)/((TA_in(j)-2.0*TL_in(j))**2)*dcdn(i,j)
! 	enddo
! enddo
