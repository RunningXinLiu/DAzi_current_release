subroutine voronoiproj_layer_iso(leniw,lenrw,colg,nrow,rw,dres,goxd,dvxd,gozd,dvzd,depz,&
                                nx,ny,nz,nd,ncell,acell,nzrand,spfra,hvratio,damp,iproj,iter,dv,&
                                rnorm_out, acond_out, itn_out, istop_out)
        use lsmrModule, only:lsmr

        implicit none
        integer leniw,lenrw
        integer nx,ny,nz
    !      integer iw(leniw)
        integer colg(lenrw),nrow(nd)
        real depz(nz)
        real rw(lenrw)
        integer ncells,acells,ncell,acell,athresh
        real dv(*),dres(*)
        real goxd,gozd,dvxd,dvzd
        real damp
        real hvratio,spfra,cmb
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
        integer maxnar,nzid
        integer iseed(4)
        integer :: iseed_theta(4), iseed_phi(4), iseed_rad(4), iseed_idx(4)
        real xs,ys,zs
        character(len=32) :: citer
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
        integer leniw_p,lenrw_p,leniwgp,lenrwgp
        integer start

        integer :: ncell_r, ncell_a, ncell_per
        integer, allocatable :: ause_layer(:)
        real,    allocatable :: S_layer(:)
        real :: Smax
        integer :: off, nxy, idx2d
        real, parameter :: gamma_cell = 0.5

        ! --- AUTO sizing (adapts to nd and grid extent) ---
        integer :: N_total_iso, N_rand_total, N_adapt_total
        real :: kappa_cell, rand_frac
        integer, parameter :: nxy_ref = 525      ! (23-2)*(27-2)
        integer, parameter :: nzrand_ref = 8     ! reference depth slices
        integer, parameter :: Ntot_ref = 360     ! 45 per layer * 8 layers
        integer, parameter :: Nd_ref  = 64594

        ! Per-layer effective area for random-site activation (shallow denser)
        integer, allocatable :: ruse_layer(:)
        real,    allocatable :: A_layer(:)
        real :: Amax

        ! eps scaled by typical magnitude of dis (units: km^2)
        real, parameter :: eps_tie = 1.0e-12
        ! integer time_begin,time_end ! Ying

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

        ! --- Voronoi site layout (layer-major):
        athresh = 10

        ! AUTO mode: if caller passes ncell<=0 or acell<=0, choose counts from nd and grid size.
        rand_frac = 0.67
        kappa_cell = real(Ntot_ref) / (sqrt(real(Nd_ref)) * sqrt(real(nxy_ref)/real(nxy_ref)) * sqrt(real(nzrand_ref)/real(nzrand_ref)))
        ! kappa_cell evaluates to ~1.42 for the reference case.

        nxy = (nx-2) * (ny-2)
        if (ncell <= 0 .or. acell <= 0) then
            N_total_iso = int( kappa_cell * sqrt(real(nd)) * sqrt(real(nxy)/real(nxy_ref)) * sqrt(real(max(1,nzrand))/real(nzrand_ref)) + 0.5 )
            ! Clamp to a safe range (per-layer 12..120)
            if (N_total_iso < 12*max(1,nzrand))  N_total_iso = 12*max(1,nzrand)
            if (N_total_iso > 120*max(1,nzrand)) N_total_iso = 120*max(1,nzrand)

            N_rand_total  = int(rand_frac * real(N_total_iso) + 0.5)
            N_adapt_total = N_total_iso - N_rand_total

            ncell_r = max(10, int(real(N_rand_total)/real(max(1,nzrand)) + 0.5))
            ncell_a = max(0,  int(real(N_adapt_total)/real(max(1,nzrand)) + 0.5))
        else
            ncell_r = ncell
            if (acell > athresh) then
                ncell_a = acell
            else
                ncell_a = 0
                if (iproj==1) write(*,*) 'Too few adaptive cells, not adding them!'
            endif
        endif

        ncell_per = ncell_r + ncell_a
        ncells = ncell_per * nzrand

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
        ! Normal cells: same (lat,lon) pattern replicated across nzrand layers
        call slarnv(1, iseed_theta, ncell_r, theta(1:ncell_r))
        theta(1:ncell_r) = (gozd + theta(1:ncell_r) * (ny-3) * dvzd) * pi/180
        do ii = 1, nzrand-1
            theta(ii*ncell_r+1:(ii+1)*ncell_r) = theta(1:ncell_r)
        enddo

        call slarnv(1, iseed_phi, ncell_r, phi(1:ncell_r))
        phi(1:ncell_r) = pi/2 - (goxd - phi(1:ncell_r) * (nx-3) * dvxd) * pi/180
        do ii = 1, nzrand-1
            phi(ii*ncell_r+1:(ii+1)*ncell_r) = phi(1:ncell_r)
        enddo

        allocate(radnz(nzrand))
        call slarnv(1, iseed_rad, nzrand, radnz)
        radnz = radnz * depz(nz-1)

        do ii = 1, nzrand
            do jj = 1, ncell_r
                rrad((ii-1)*ncell_per + jj) = radius - radnz(ii) * hvratio
            enddo
        enddo

        allocate(S_layer(nz-1), ause_layer(nz-1), A_layer(nz-1), ruse_layer(nz-1))

        S_layer = 0.0
        A_layer = 0.0
        ruse_layer = ncell_r

        if(iproj==1) then
            write(*,*) 'Voronoi cell counts: ncell_r=', ncell_r, ', ncell_a=', ncell_a, ', total per layer=', ncell_per, ', nzrand=', nzrand
        endif

        ! adaptive cells based on dws, assume 1/2 of all ncells are used
        ! as adaptive cells

        if (acell > athresh) then
            dws = 0
            ! write(*,*) 'Adding ISO cells!'
            do ii = 1,lenrw
                dws(colg(ii)) = dws(colg(ii))+abs(rw(ii))
            enddo

            do iz = 1, nz-1
                off = (iz-1) * (nx-2) * (ny-2)
                S_layer(iz) = sum( dws(off+1:off+(nx-2)*(ny-2)) )
                ! Effective illuminated area: count points above 20% of layer max
                if (maxval(dws(off+1:off+nxy)) > 0.0) then
                    A_layer(iz) = real( count( dws(off+1:off+nxy) > 0.2*maxval(dws(off+1:off+nxy)) ) )
                else
                    A_layer(iz) = 0.0
                endif
            enddo
        else
            S_layer = 0.0
            A_layer = 0.0
        endif

        Smax = maxval(S_layer)
        Amax = maxval(A_layer)

        do iz = 1, nz-1
            ! Adaptive active count
            if (ncell_a == 0 .or. Smax <= 0.0) then
                ause_layer(iz) = 0
            else
                ause_layer(iz) = int( real(ncell_a) * (max(0.0, S_layer(iz)/Smax))**gamma_cell + 0.5 )
                if (ause_layer(iz) < 0) ause_layer(iz) = 0
                if (ause_layer(iz) > ncell_a) ause_layer(iz) = ncell_a
            endif

            ! Random active count (nested subset) for shallow dense / deep sparse
            if (Amax > 0.0) then
                ruse_layer(iz) = int( real(ncell_r) * (max(0.0, A_layer(iz)/Amax))**gamma_cell + 0.5 )
            else
                ruse_layer(iz) = ncell_r
            endif
            if (ruse_layer(iz) < 8) ruse_layer(iz) = 8
            if (ruse_layer(iz) > ncell_r) ruse_layer(iz) = ncell_r
        enddo

        if (iproj == 1) then
            if (ncell <= 0 .or. acell <= 0) then
                write(*,*) 'AUTO cell sizing enabled.'
                write(*,*) 'AUTO totals: N_total_iso=', N_total_iso, ' N_rand_total=', N_rand_total, ' N_adapt_total=', N_adapt_total
            endif
            write(*,*) 'Per-layer (max) Voronoi sites: random=', ncell_r, ' adaptive_max=', ncell_a, '  ncell_per=', ncell_per
            write(*,*) 'Active sites per depth layer: iz  ruse  ause  S_layer  A_layer'
            do iz = 1, nz-1
                write(*,'(I4,2X,I6,2X,I6,2X,1PE12.4,2X,1PE12.4)') iz, ruse_layer(iz), ause_layer(iz), S_layer(iz), A_layer(iz)
            enddo
        endif

        do ii = 1, ncell_a
            call random_index_lapack(idx, dws, iseed_idx)
            dws(idx) = 0.0  ! reduce the probability of being selected again
            !col(idx) = (iz-1)*(nx-2)*(ny-2)+(iy-1)*(nx-2)+ix
            ! the ix, iy is wrong, bug fixed by Ying
            !ix = mod(idx,nx-2)
            !iy = idx/(nx-2)
            iz = (idx-1)/((nx-2)*(ny-2))+1  ! Ying
            iy = (idx-1+(1-iz)*(nx-2)*(ny-2))/(nx-2)+1  ! Ying
            ix = idx+(nx-2)*(1-iy)+(nx-2)*(ny-2)*(1-iz)  ! Ying
            do jj = 1, nzrand
                theta(nzrand*ncell_r+(ii-1)*nzrand+jj) = (gozd+(iy-1)*dvzd)*pi/180
                phi(nzrand*ncell_r+(ii-1)*nzrand+jj) = pi/2-(goxd-(ix-1)*dvxd)*pi/180
                rrad(nzrand*ncell_r+(ii-1)*nzrand+jj) = radius-radnz(jj)*hvratio ! because the depz(0) is not the boundary
                ! rrad(nzrand*ncell+(ii-1)*nzrand+jj) = radius-depz(iz)*hvratio ! because the depz(0) is not the boundary
            enddo
        enddo

        xpts = rrad*sin(phi)*cos(theta)
        ypts = rrad*sin(phi)*sin(theta)
        zpts = rrad*cos(phi)

        if (iproj == 1 .and. acell > athresh .and. iter == 1) then
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

        if (iproj == 1 .and. acell > athresh .and. iter == 1) then
            open(38,file='col.txt')
        endif
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

                if (iz <= nz-1) then
                    integer :: j0
                    real, parameter :: huge_penalty = 1.0e30
                    j0 = (iz-1) * ncell_per

                    ! Disable inactive random sites (keep nested subset 1..ruse_layer)
                    if (ruse_layer(iz) < ncell_r) then
                        do jj = j0 + ruse_layer(iz) + 1, j0 + ncell_r
                            dis(jj) = dis(jj) + huge_penalty
                        enddo
                    endif

                    ! Disable inactive adaptive sites (keep nested subset 1..ause_layer)
                    if (ncell_a > 0) then
                        if (ause_layer(iz) < ncell_a) then
                            do jj = j0 + ncell_r + ause_layer(iz) + 1, j0 + ncell_per
                                dis(jj) = dis(jj) + huge_penalty
                            enddo
                        endif
                    endif
                endif

                row(idx) = minloc(dis,1)

                if (iproj == 1 .and.  acell > athresh .and. iter == 1) then
                    write(38,*) ix, iy, iz, col(idx), lat(ix)*180/pi, lon(iy)*180/pi, dws(col(idx))
                endif

                enddo
            enddo
        enddo
        if (iproj == 1 .and. acell > athresh .and. iter == 1) then
            close(38)
        endif



        rw_p = 1.0
        leniw_p = 2*ndim+1
        lenrw_p = ndim
        iw_p(1) = ndim
        iw_p(2:ndim+1) = row
        iw_p(ndim+2:2*ndim+1) = col

        allocate(grow(ndim),gcol(nd),subrow(ncells))
        ! maxnar = int(0.6*nd*ncells)
        maxnar = int(spfra*nd*ncells)  ! *3 for because of anisotropic case need more
        if((iproj==1)) then
            write(*,*) 'Max Number of sparse matrix (GP) non-zero element ',maxnar
        endif
        allocate(iwgp(maxnar*2+1),colgp(maxnar),rwgp(maxnar))

        ! call system_clock(time_begin) 
        nzid = 0
        do ii = 1,nd
            grow = 0
            start = sum(nrow(1:ii-1))
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
                    if(nzid .ge. maxnar) then
                        write(6,*) 'Voronoi projection part sparsefrac is too small, please increase your sparsefrac'
                        stop 'nzid > maxnar, please increase your sparsefrac'
                    endif
                    nzid = nzid+1
                    rwgp(nzid) = subrow(ix)
                    iwgp(1+nzid) = ii
                    colgp(nzid) = ix
                endif
            enddo
        enddo

        if (nzid.gt.maxnar) then
            write(6,*) 'Voronoi projection part sparsefrac is too small, please increase your sparsefrac'
            stop 'nzid > maxnar, please increase your sparsefrac'
        endif
        leniwgp = nzid*2+1
        lenrwgp = nzid
        iwgp(1) = lenrwgp
        iwgp(nzid+2:nzid*2+1) = colgp(1:nzid)
        ! call system_clock(time_end) 
        ! if(iproj==1) print *, "time for projecting grid to cells:",(time_end-time_begin)/1000

        !norm = 0
        !do ii=1,nzid
        !norm(iwgp(1+ii+nzid)) = norm(iwgp(1+ii+nzid))+rwgp(ii)**2
        !enddo

        !do ii =1,ncells
        !norm(ii) = sqrt(norm(ii)/nd+0.01)
        !enddo

        !do ii =1,nzid
        !rwgp(ii) = rwgp(ii)/norm(iwgp(1+ii+nzid))
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

        !print*, 'sth is wrong...'
        ! call system_clock(time_begin) 
        call LSMR(nd, ncells, leniwgp, lenrwgp,iwgp,rwgp,dres,damp,&
        atol, btol, conlim, itnlim, localSize,nout,&
        xunknown, istop, itn, anorm, acond,rnorm, arnorm, xnorm)
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
        if (allocated(S_layer)) deallocate(S_layer)
        if (allocated(ause_layer)) deallocate(ause_layer)
        if (allocated(A_layer)) deallocate(A_layer)
        if (allocated(ruse_layer)) deallocate(ruse_layer)

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
