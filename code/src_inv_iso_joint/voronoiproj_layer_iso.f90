subroutine voronoiproj_layer_iso(leniw,lenrw,colg,nrow,rw,dres,goxd,dvxd,gozd,dvzd,depz,&
                                nx,ny,nz,nd,ncell,acell,nzrand,spfra,hvratio,damp,iproj,iter,dv,&
                                rnorm_out, acond_out, itn_out, istop_out)
        use lsmrModule, only:lsmr

        implicit none
        integer(kind=8), intent(in) :: leniw, lenrw
        integer nx,ny,nz
    !      integer iw(leniw)
        integer colg(lenrw),nrow(nd)
        real depz(nz)
        real rw(lenrw)
        integer ncells,acells,ncell,acell
        real dv(*),dres(*)
        real goxd,gozd,dvxd,dvzd
        real damp
        real(kind=8) :: spfra
        real :: hvratio, cmb
        integer ndim,nd
        integer iproj,iter,nzrand
        ! Output arguments for LSMR quality metrics
        real, intent(out) :: rnorm_out, acond_out
        integer, intent(out) :: itn_out, istop_out

        real,parameter:: radius = 6371.0,ftol = 0.0001,pi = 3.141592654
        integer ii,jj,ix,iy,iz
        real,dimension(:),allocatable:: grow,gcol,subrow,dis,dws,xunknown
        real,dimension(:),allocatable:: lat,lon,rad,theta,phi,rrad,xpts,ypts,zpts
        real,dimension(:),allocatable :: rw_p,rwgp,norm
        real,dimension(:),allocatable:: radnz
        integer,dimension(:),allocatable:: iw_p,row,col,iwgp,colgp
        integer idx
        integer(kind=8) :: maxnar, nzid, i
        integer :: iseed(4)
        integer :: iseed_theta(4), iseed_phi(4), iseed_rad(4), iseed_idx(4)
        integer :: iseed_theta_l(4), iseed_phi_l(4)
        real xs,ys,zs
        character(len=32) :: citer
        character(len=64) :: lsmr_backend, lsmr_tag
        character(len=32) :: gp_backend, gp_spfra_env
        character(len=16) :: timer_env
        logical :: timer_enabled, timer_log_open
        integer :: timer_status, clock_rate
        integer :: t_start, t_cells, t_map, t_gp, t_lsmr0, t_lsmr1, t_end
        real rx

        real atol,btol
        real conlim
        integer istop
        integer itnlim
        real acond
        real anorm
        real arnorm
        real rnorm
        real xnorm
        integer localSize,nout,itn
        integer env_status
        integer gp_env_status, gp_spfra_status, grid_id
        integer(kind=8) :: leniw_p, lenrw_p, leniwgp, lenrwgp
        integer(kind=8) :: start
        real(kind=8) :: gp_spfra_eff, gp_spfra_value
        logical :: use_adaptive_cells

        ! eps scaled by typical magnitude of dis (units: km^2)
        real, parameter :: eps_tie = 1.0e-12
        ! integer time_begin,time_end ! Ying

        timer_env = ''
        call get_environment_variable('DAZI_VORO_TIMER', timer_env, status=timer_status)
        timer_enabled = (timer_status == 0 .and. &
                        (trim(adjustl(timer_env)) == '1' .or. trim(adjustl(timer_env)) == 'T' .or. &
                         trim(adjustl(timer_env)) == 'true'))
        if (timer_enabled) call system_clock(t_start, clock_rate)

        allocate(lat(nx-2),lon(ny-2),rad(nz-1))
        ndim = (nx-2)*(ny-2)*(nz-1)

        do ii = 1,nx-2
        lat(ii) = (goxd-(ii-1)*dvxd)*pi/180
        enddo

        do ii = 1,ny-2
        lon(ii) = (gozd+(ii-1)*dvzd)*pi/180
        enddo

        !cmb = radius - depz(nz-1)*hvratio
        do ii = 1,nz-1
        rad(ii) = radius-depz(ii)*hvratio
        !rad(ii) = cmb+depz(ii)*hvratio
        enddo

        ! nzrand = int((nz-1)/2)
        ! acells = int(ncell/2)
        acells = acell*nzrand

        use_adaptive_cells = (acell > 0)
        if(use_adaptive_cells) then
            ! write(*,*) 'Adding ISO cells!'
            ncells = ncell*nzrand+acell*nzrand
        else
            if (iproj == 1 .and. iter == 1) write(*,*) 'No adaptive cells requested.'
            ncells = ncell*nzrand
        endif
        allocate(theta(ncells),phi(ncells),rrad(ncells),norm(ncells))
        allocate(xpts(ncells),ypts(ncells),zpts(ncells),dis(ncells),xunknown(ncells))
        allocate(rw_p(ndim))
        allocate(iw_p(2*ndim+1),row(ndim),col(ndim),dws(ndim))

        ! Deterministic, OMP-safe RNG streams per (iproj,iter)
        iseed(1:3) = (/38,62,346/)
        iseed(4)   = 2*(iproj+iter)+1
        iseed(4) = mod(iseed(4), 4095)
        if (iseed(4) == 0) iseed(4) = 1
        iseed_theta = iseed
        iseed_phi   = iseed
        iseed_rad   = iseed
        iseed_idx   = iseed
        iseed_theta(4) = iseed_theta(4) + 11
        iseed_phi(4)   = iseed_phi(4)   + 23
        iseed_rad(4)   = iseed_rad(4)   + 37
        iseed_idx(4)   = iseed_idx(4)   + 59
        iseed_theta(4) = mod(iseed_theta(4), 4095); if (iseed_theta(4)==0) iseed_theta(4)=1
        iseed_phi(4)   = mod(iseed_phi(4),   4095); if (iseed_phi(4)==0)   iseed_phi(4)=1
        iseed_rad(4)   = mod(iseed_rad(4),   4095); if (iseed_rad(4)==0)   iseed_rad(4)=1
        iseed_idx(4)   = mod(iseed_idx(4),   4095); if (iseed_idx(4)==0)   iseed_idx(4)=1
        ! Normal cells: independent random (lat,lon) per layer (no replication across layers)
        do ii = 1, nzrand
            ! Make a deterministic, layer-specific seed derived from the base stream
            iseed_theta_l = iseed_theta
            iseed_phi_l   = iseed_phi

            ! Offset the 4th seed component by layer index to decorrelate layers
            iseed_theta_l(4) = iseed_theta_l(4) + 97 * ii
            iseed_phi_l(4)   = iseed_phi_l(4)   + 193 * ii

            iseed_theta_l(4) = mod(iseed_theta_l(4), 4095); if (iseed_theta_l(4) == 0) iseed_theta_l(4) = 1
            iseed_phi_l(4)   = mod(iseed_phi_l(4),   4095); if (iseed_phi_l(4)   == 0) iseed_phi_l(4)   = 1

            call slarnv(1, iseed_theta_l, ncell, theta((ii-1)*ncell+1:ii*ncell))
            theta((ii-1)*ncell+1:ii*ncell) = (gozd + theta((ii-1)*ncell+1:ii*ncell) * (ny-3) * dvzd) * pi/180

            call slarnv(1, iseed_phi_l, ncell, phi((ii-1)*ncell+1:ii*ncell))
            phi((ii-1)*ncell+1:ii*ncell) = pi/2 - (goxd - phi((ii-1)*ncell+1:ii*ncell) * (nx-3) * dvxd) * pi/180
        enddo

        allocate(radnz(nzrand))
        call slarnv(1, iseed_rad, nzrand, radnz)
        radnz = radnz * depz(nz-1)

        do ii = 1, nzrand
            do jj = 1, ncell
                rrad((ii-1)*ncell + jj) = radius - radnz(ii) * hvratio
            enddo
        enddo

        if(iproj==1 .and. iter == 1) then
            open(77,file='pvpoints.dat')
            do ii=1,nzrand
            do jj=1,ncell
            write(77,'(3f10.4)') rrad((ii-1)*ncell+jj),phi((ii-1)*ncell+jj),theta((ii-1)*ncell+jj)
            enddo
            enddo
            close(77)
            !print*,goxd,gozd,dvxd,dvzd,nx,ny
            !stop
        endif
        !call slarnv(1,iseed,ncell,rrad)
        !rrad = radius-rrad*depz(nz-1)*hvratio

        ! adaptive cells based on dws, assume 1/2 of all ncells are used
        ! as adaptive cells
        if (use_adaptive_cells) then
        dws = 0
        ! write(*,*) 'Adding ISO cells!'
        do i = 1,lenrw
            dws(colg(i)) = dws(colg(i))+abs(rw(i))
        enddo
 
        do ii = 1, acell
            call random_index_lapack(idx, dws, iseed_idx)
            dws(idx) = 0.0  ! reduce the probability of being selected again
            !col(idx) = (iz-1)*(nx-2)*(ny-2)+(iy-1)*(nx-2)+ix
            ! the ix, iy is wrong, bug fixed by Ying
            !ix = mod(idx,nx-2)
            !iy = idx/(nx-2)
            iz = (idx-1)/((nx-2)*(ny-2))+1  ! Ying
            iy = (idx-1+(1-iz)*(nx-2)*(ny-2))/(nx-2)+1  ! Ying
            ix = idx+(nx-2)*(1-iy)+(nx-2)*(ny-2)*(1-iz)  ! Ying
            do jj = 1,nzrand
                theta(nzrand*ncell+(ii-1)*nzrand+jj) = (gozd+(iy-1)*dvzd)*pi/180
                phi(nzrand*ncell+(ii-1)*nzrand+jj) = pi/2-(goxd-(ix-1)*dvxd)*pi/180
                rrad(nzrand*ncell+(ii-1)*nzrand+jj) = radius-radnz(jj)*hvratio ! because the depz(0) is not the boundary
                ! rrad(nzrand*ncell+(ii-1)*nzrand+jj) = radius-depz(iz)*hvratio ! because the depz(0) is not the boundary
            enddo
        enddo
    endif

    xpts = rrad*sin(phi)*cos(theta)
    ypts = rrad*sin(phi)*sin(theta)
    zpts = rrad*cos(phi)

    if (iproj == 1 .and. use_adaptive_cells .and. iter == 1) then
        open(77,file='pvpoints_acell.dat')
        do ii=1,ncells
        write(77,'(3f10.4)') rrad(ii),phi(ii),theta(ii)
        enddo
        close(77)
        !print*,goxd,gozd,dvxd,dvzd,nx,ny
        !stop
        ! ! Ying output points for plot
        ! open(36,file='vorosites.txt')
        ! do ii = 1,ncells
        ! write(36, '(3f10.2)') xpts(ii), ypts(ii), zpts(ii)
        ! enddo
        ! close(36)

    endif

    if (iproj == 1) then
        write(citer,'(I0)') iter
        open(37, file='vorosites_latlon_iter'//trim(citer)//'.txt')
        do ii = 1, ncells
            write(37, '(3f10.2)') theta(ii)*180/pi, (pi/2-phi(ii))*180/pi, (radius-rrad(ii))/hvratio
        enddo
        close(37)
    endif

        if (iproj == 1 .and. use_adaptive_cells .and. iter == 1) then
        open(38,file='col.txt')
    endif
    if (timer_enabled) call system_clock(t_cells)
    idx = 0
    do iz = 1,nz-1
        do iy = 1,ny-2
            do ix = 1,nx-2
            xs = rad(iz)*sin(pi/2-lat(ix))*cos(lon(iy))
            ys = rad(iz)*sin(pi/2-lat(ix))*sin(lon(iy))
            zs = rad(iz)*cos(pi/2-lat(ix))
            dis =  (xpts-xs)**2+(ypts-ys)**2+(zpts-zs)**2
            ! Deterministic tie-break to stabilize minloc when distances are nearly equal
            ! eps scaled by typical magnitude of dis (units: km^2)
            ! real, parameter :: eps_tie = 1.0e-12
            do jj = 1, ncells
                dis(jj) = dis(jj) + eps_tie * real(jj)
            enddo
            idx = idx+1
            col(idx) = (iz-1)*(nx-2)*(ny-2)+(iy-1)*(nx-2)+ix
            row(idx) = minloc(dis,1)

            if (iproj == 1 .and.  use_adaptive_cells .and. iter == 1) then
                write(38,*) ix, iy, iz, col(idx), lat(ix)*180/pi, lon(iy)*180/pi, dws(col(idx))
            endif

            enddo
        enddo
    enddo
    if (iproj == 1 .and. use_adaptive_cells .and. iter == 1) then
        close(38)
    endif



    rw_p = 1.0
    leniw_p = 2_8*int(ndim,kind=8)
    lenrw_p = int(ndim,kind=8)
    iw_p(1:ndim) = row
    iw_p(ndim+1:2*ndim) = col
    if (timer_enabled) call system_clock(t_map)

    gp_spfra_eff = spfra
    gp_spfra_env = ''
    call get_environment_variable('DAZI_VORO_GP_SPFRA', gp_spfra_env, status=gp_spfra_status)
    if (gp_spfra_status == 0) then
        read(gp_spfra_env,*,iostat=gp_spfra_status) gp_spfra_value
        if (gp_spfra_status == 0 .and. gp_spfra_value > 0.0_8) gp_spfra_eff = gp_spfra_value
    endif

    allocate(grow(ndim),gcol(nd),subrow(ncells))
    ! maxnar = int(0.6*nd*ncells)
    maxnar = int(gp_spfra_eff * real(nd,kind=8) * real(ncells,kind=8), kind=8)  ! projection-only sparse allocation
    if((iproj==1)) then
        write(*,*) 'Max Number of sparse matrix (GP) non-zero element ',maxnar
        write(*,'(A,F8.3)') 'Vorotomo GP sparse fraction: ', gp_spfra_eff
    endif
    allocate(iwgp(2_8*maxnar), colgp(maxnar), rwgp(maxnar))

    gp_backend = 'direct'
    call get_environment_variable('DAZI_VORO_GP_BACKEND', gp_backend, status=gp_env_status)
    if (iproj == 1 .and. trim(adjustl(gp_backend)) == 'direct') then
        write(*,'(A)') 'Vorotomo GP backend: direct sparse row projection'
    endif

    ! call system_clock(time_begin) 
    nzid = 0_8
    if (trim(adjustl(gp_backend)) == 'direct') then
        start = 0_8
        do ii = 1,nd
            subrow = 0
            do ix = 1,nrow(ii)
                grid_id = colg(start+ix)
                if (grid_id < 1 .or. grid_id > ndim) then
                    stop 'DAZI_VORO_GP_BACKEND=direct: grid column out of range'
                endif
                subrow(row(grid_id)) = subrow(row(grid_id)) + rw(start+ix)
            enddo
            do ix = 1,ncells
                if(abs(subrow(ix))>ftol) then
                    if (nzid >= maxnar) then
                        write(6,*) 'Voronoi projection part sparsefrac is too small, please increase your sparsefrac'
                        stop 'nzid > maxnar, please increase your sparsefrac'
                    endif
                    nzid = nzid +1_8
                    rwgp(nzid) = subrow(ix)
                    iwgp(nzid) = ii
                    colgp(nzid) = ix
                endif
            enddo
            start = start + int(nrow(ii), kind=8)
        enddo
    else
        start = 0_8
        do ii = 1,nd
            grow = 0
            do ix = 1,nrow(ii)
                grow(colg(start+ix)) = rw(start+ix)
            enddo
            !gcol = 0
            !gcol(ii) = 1.0 
            !call aprod(2,nd,ndim,grow,gcol,leniw,lenrw,iw,rw)
            subrow = 0

            call aprod(1,ncells,ndim,grow,subrow,leniw_p,lenrw_p,iw_p,rw_p)
            do ix = 1,ncells
                if(abs(subrow(ix))>ftol) then
                    if (nzid >= maxnar) then
                        write(6,*) 'Voronoi projection part sparsefrac is too small, please increase your sparsefrac'
                        stop 'nzid > maxnar, please increase your sparsefrac'
                    endif
                    nzid = nzid +1_8
                    rwgp(nzid) = subrow(ix)
                    iwgp(nzid) = ii
                    colgp(nzid) = ix
                endif
            enddo
            start = start + int(nrow(ii), kind=8)
        enddo
    endif

    if (nzid > maxnar) then
        write(6,*) 'Voronoi projection part sparsefrac is too small, please increase your sparsefrac'
        stop 'nzid > maxnar, please increase your sparsefrac'
    endif
    leniwgp = 2_8*nzid
    lenrwgp = nzid

    ! ---- Pack into LSMR triplet header layout ----
    ! iwgp(1:nnz)         = row indices (already stored during assembly as iwgp(1+k))
    ! iwgp(nnz+1:2*nnz)   = col indices
    if (nzid > int(huge(iwgp(1)), kind=8)) then
        stop 'nzid exceeds INTEGER(4) capacity for iwgp header'
    end if
    iwgp(int(nzid,kind=4)+1 : 2*int(nzid,kind=4)) = colgp(1:int(nzid,kind=4))
    if (timer_enabled) call system_clock(t_gp)

    ! call system_clock(time_end) 
    ! if(iproj==1) print *, "time for projecting grid to cells:",(time_end-time_begin)/1000

    !norm = 0
    !do ii=1,nzid
    !norm(iwgp(ii+nzid)) = norm(iwgp(ii+nzid))+rwgp(ii)**2
    !enddo

    !do ii =1,ncells
    !norm(ii) = sqrt(norm(ii)/nd+0.01)
    !enddo

    !do ii =1,nzid
    !rwgp(ii) = rwgp(ii)/norm(iwgp(ii+nzid))
    !enddo

    ! conlim = 100
    ! itnlim = 400
    ! atol = 1e-6!1e-3/((dvxd+dvzd)*111.19/2.0*0.1) !1e-2
    ! btol = 1e-6!1e-3/(dvxd*nx*111.19/3.0)!1e-3
    conlim = 50
    itnlim = 100
    atol = 1e-3/((dvxd+dvzd)*111.19/2.0*0.1) !1e-2
    btol = 1e-3/(dvxd*nx*111.19/3.0)!1e-3
    istop = 0
    anorm = 0.0
    acond = 0.0
    arnorm = 0.0
    xnorm = 0.0
    ! localSize = 10!int(ncells/4)
    localSize = int(ncells/4)
    !damp = dampvel
    ! using lsmr to solve for the projection coefficients
    !print*, 'LSMR beginning ...'

    nout = -1
    !nout = 36
    !open(nout,file='lsmrout_sub.txt')

    lsmr_backend = ''
    call get_environment_variable('DAZI_VORO_LSMR_BACKEND', lsmr_backend, status=env_status)
    if (env_status /= 0) then
        call get_environment_variable('DAZI_LSMR_BACKEND', lsmr_backend, status=env_status)
    endif

    !print*, 'sth is wrong...'
    ! call system_clock(time_begin) 
    if (timer_enabled) call system_clock(t_lsmr0)
    if (trim(adjustl(lsmr_backend)) == 'scipy') then
        write(lsmr_tag,'("voro_iso_i",I0,"_r",I0)') iter, iproj
        call scipy_lsmr_solve(lsmr_tag, nd, ncells, leniwgp, lenrwgp, iwgp, rwgp, dres, damp, &
        atol, btol, conlim, itnlim, xunknown, istop, itn, anorm, acond, rnorm, arnorm, xnorm)
    else
        call LSMR(nd, ncells, leniwgp, lenrwgp,iwgp,rwgp,dres,damp,&
        atol, btol, conlim, itnlim, localSize,nout,&
        xunknown, istop, itn, anorm, acond,rnorm, arnorm, xnorm)
    endif
    if (timer_enabled) call system_clock(t_lsmr1)
    ! call system_clock(time_end) 
    ! if(iproj==1) print *, "time for lsqr:",(time_end-time_begin)/1000
    !close(nout)
    !do ii = 1,ncells
    !  xunknown(ii) = xunknown(ii)/norm(ii)
    !enddo

    !norm = (norm**2-0.01)*nd
    !do ii = 1,ncells
    !if (norm(ii)<0.01) then
    !    call random_number(rx)
    !    xunknown(ii) = xunknown(ii)+rx-0.5
    !endif
    !enddo

    ! call system_clock(time_begin) 
    dv(1:ndim) = 0
    call aprod(2,ncells,ndim,dv,xunknown,leniw_p,lenrw_p,iw_p,rw_p)
    if (timer_enabled) then
        call system_clock(t_end)
        inquire(unit=66, opened=timer_log_open)
        write(*,'(A,I0,A,I0,A,5(F10.4,A))') 'Voro timer iso iter=', iter, ' real=', iproj, &
            ' cells=', real(t_cells-t_start)/real(clock_rate), 's map=', &
            real(t_map-t_cells)/real(clock_rate), 's gp=', &
            real(t_gp-t_map)/real(clock_rate), 's lsmr=', &
            real(t_lsmr1-t_lsmr0)/real(clock_rate), 's back=', &
            real(t_end-t_lsmr1)/real(clock_rate), 's'
        if (timer_log_open) then
            write(66,'(A,I0,A,I0,A,5(F10.4,A))') 'Voro timer iso iter=', iter, ' real=', iproj, &
                ' cells=', real(t_cells-t_start)/real(clock_rate), 's map=', &
                real(t_map-t_cells)/real(clock_rate), 's gp=', &
                real(t_gp-t_map)/real(clock_rate), 's lsmr=', &
                real(t_lsmr1-t_lsmr0)/real(clock_rate), 's back=', &
                real(t_end-t_lsmr1)/real(clock_rate), 's'
        endif
    endif
    ! call system_clock(time_end) 
    ! if(iproj==1) print *, "time for solution:",(time_end-time_begin)/1000
    rnorm_out = rnorm
    acond_out = acond
    itn_out = itn
    istop_out = istop
    deallocate(grow,gcol,subrow)
    deallocate(theta,phi,rrad,dws,norm)
    deallocate(xpts,ypts,zpts,dis,xunknown)
    deallocate(iw_p,rw_p,row,col)
    deallocate(lat,lon,rad)
    deallocate(iwgp,colgp,rwgp)
    deallocate(radnz)

contains

subroutine random_index_lapack(idx, weights, iseed)
    implicit none
    integer, intent(out) :: idx
    real,    intent(in)  :: weights(:)
    integer, intent(inout) :: iseed(4)

    real :: x, wsum, prob
    real :: tmp(1)
    integer :: n

    n = size(weights)
    wsum = sum(weights)

    ! Draw x ~ U(0,1) from LAPACK RNG (thread-safe when seed is private)
    call slarnv(1, iseed, 1, tmp)
    x = tmp(1)
    if (x < 0.0) x = 0.0
    if (x > 1.0) x = 1.0

    if (wsum <= 0.0) then
        ! Fallback: uniform index when weights are all zero
        idx = 1 + int(x * n)
        if (idx < 1) idx = 1
        if (idx > n) idx = n
        return
    endif

    prob = 0.0
    do idx = 1, n
        prob = prob + weights(idx) / wsum
        if (x <= prob) exit
    enddo

    if (idx < 1) idx = 1
    if (idx > n) idx = n
end subroutine random_index_lapack


end subroutine voronoiproj_layer_iso
