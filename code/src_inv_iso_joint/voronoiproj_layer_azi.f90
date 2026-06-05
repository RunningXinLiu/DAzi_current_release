subroutine voronoiproj_layer_azi(leniw,lenrw,colg,nrow,rw,dres,goxd,dvxd,gozd,dvzd,depz,&
                               nx,ny,nz,nd,ncell,acell,nzrand,spfra,hvratio,iaratio,damp,iproj,iter,dv,&
                               rnorm_out, acond_out, itn_out, istop_out)
      ! This subroutine builds the Voronoi projection matrix for both Vs and azimuthal anisotropy
      ! and solves for the projection coefficients using LSMR.
      ! The cell for Vs and Gc Gs are built separately to allow different complexity.
      use lsmrModule, only:lsmr

      implicit none
      integer(kind=8), intent(in) :: leniw, lenrw
      integer(kind=8) :: lenrw_iso
      integer nx,ny,nz
!      integer iw(leniw)
      integer colg(lenrw),nrow(nd)
      real depz(nz)
      real rw(lenrw)
      integer ncells,acells,ncell,acell
      integer :: ncells_v, ncells_a
      integer :: ncell_a, acell_a
      real dv(*),dres(*)
      real goxd,gozd,dvxd,dvzd
      real damp
      real(kind=8) :: spfra
      real :: hvratio, cmb, iaratio
      integer ndim,nd,ndim0
      integer iproj,iter,nzrand
      ! Output arguments for LSMR quality metrics
      real, intent(out) :: rnorm_out, acond_out
      integer, intent(out) :: itn_out, istop_out
      real,parameter:: radius = 6371.0,ftol = 0.0001,pi = 3.141592654
      integer ii,jj,ix,iy,iz
      real,dimension(:),allocatable:: lat,lon,rad
      real,dimension(:),allocatable:: theta_v,phi_v,rrad_v,xpts_v,ypts_v,zpts_v
      real,dimension(:),allocatable:: theta_a,phi_a,rrad_a,xpts_a,ypts_a,zpts_a
      real,dimension(:),allocatable :: rw_p,rwgp,norm
      real,dimension(:),allocatable:: radnz
      integer,dimension(:),allocatable:: iw_p,row_v,row_a,col,iwgp,colgp
      integer idx

      integer :: nunk
      integer :: base_a
      integer,dimension(:),allocatable :: ncell_v_l, acell_v_l
      integer,dimension(:),allocatable :: ncell_a_l, acell_a_l
      integer,dimension(:),allocatable :: off_nv_l, off_av_l
      integer,dimension(:),allocatable :: off_na_l, off_aa_l
      integer :: cnt_layer, layer_start
      real,dimension(:),allocatable :: dis
      real,dimension(:),allocatable :: xunknown
      real,dimension(:),allocatable :: grow, gcol, subrow
      real,dimension(:),allocatable :: dws_v, dws_a

      real,dimension(:),allocatable :: rw_p_3
      real,dimension(:),allocatable:: rw_iso, dws_iso
      integer(kind=8),dimension(:),allocatable :: idx_iso
      integer,dimension(:),allocatable :: colg_iso
      integer,dimension(:),allocatable:: iw_p_3

      integer(kind=8) :: maxnar,nzid,i
      integer iseed(4)
      integer :: iseed_theta(4), iseed_phi(4), iseed_rad(4), iseed_idx(4)
      integer :: iseed_theta_l(4), iseed_phi_l(4)
      integer :: iseed_theta_a(4), iseed_phi_a(4)
      integer :: iseed_theta_a_l(4), iseed_phi_a_l(4)
      real xs,ys,zs
      character(len=32) :: citer
      character(len=64) :: lsmr_backend, lsmr_tag
      character(len=32) :: gp_backend, gp_spfra_env
      character(len=64) :: azi_vs_mode
      character(len=64) :: taper_mode
      character(len=16) :: timer_env
      logical :: timer_enabled, timer_log_open
      logical :: freeze_vs
      logical :: taper_enabled, taper_vs, taper_ani, taper_deep_suppress, taper_yunnan
      integer :: timer_status, clock_rate
      integer :: t_start, t_cells, t_map, t_gp, t_lsmr0, t_lsmr1, t_end
      real rx
      real :: zfrac, f_vs, f_ani

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
      integer azi_vs_status
      integer taper_status
      integer gp_env_status, gp_spfra_status, model_id, grid_id, cell_id
      integer(kind=8) :: leniw_p,lenrw_p,leniwgp,lenrwgp
      integer(kind=8) :: leniw_p_3,lenrw_p_3
      integer(kind=8) :: start
      real(kind=8) :: gp_spfra_eff, gp_spfra_value
      logical :: use_adaptive_cells
      ! eps scaled by typical magnitude of dis (units: km^2)
      real, parameter :: eps_tie = 1.0e-12

      timer_env = ''
      call get_environment_variable('DAZI_VORO_TIMER', timer_env, status=timer_status)
      timer_enabled = (timer_status == 0 .and. &
                      (trim(adjustl(timer_env)) == '1' .or. trim(adjustl(timer_env)) == 'T' .or. &
                       trim(adjustl(timer_env)) == 'true'))
      if (timer_enabled) call system_clock(t_start, clock_rate)
      azi_vs_mode = 'joint'
      call get_environment_variable('DAZI_AZI_VS_MODE', azi_vs_mode, status=azi_vs_status)
      if (azi_vs_status /= 0 .or. len_trim(azi_vs_mode) == 0) azi_vs_mode = 'joint'
      freeze_vs = (trim(adjustl(azi_vs_mode)) == 'freeze' .or. trim(adjustl(azi_vs_mode)) == 'frozen' .or. &
                   trim(adjustl(azi_vs_mode)) == 'ani_only' .or. trim(adjustl(azi_vs_mode)) == 'residual')
      if (iproj == 1 .and. freeze_vs) then
        write(*,'(A)') 'Vorotomo AZI projection: frozen Vs columns; solving Gc/Gs residual only.'
      endif
      taper_mode = 'none'
      call get_environment_variable('DAZI_VORO_DEPTH_TAPER', taper_mode, status=taper_status)
      if (taper_status /= 0 .or. len_trim(taper_mode) == 0) taper_mode = 'none'
      taper_enabled = (trim(adjustl(taper_mode)) /= 'none' .and. trim(adjustl(taper_mode)) /= '0')
      taper_vs = (index(trim(adjustl(taper_mode)), 'vs_ani') > 0 .or. &
                  index(trim(adjustl(taper_mode)), 'vs_only') > 0)
      taper_ani = (index(trim(adjustl(taper_mode)), 'ani') > 0)
      taper_deep_suppress = (index(trim(adjustl(taper_mode)), 'deep_suppressed') > 0 .or. &
                             index(trim(adjustl(taper_mode)), 'deep_suppress') > 0)
      taper_yunnan = (index(trim(adjustl(taper_mode)), 'yunnan') > 0)
      if (iproj == 1 .and. taper_enabled .and. iter == 1) then
        write(*,'(A,A)') 'Vorotomo depth taper mode: ', trim(adjustl(taper_mode))
      endif

      allocate(lat(nx-2),lon(ny-2),rad(nz-1))
      ndim0 = (nx-2)*(ny-2)*(nz-1) !this is the dimension of one variable
      ndim = (nx-2)*(ny-2)*(nz-1)*3

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
      use_adaptive_cells = (acell > 0)

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

      allocate(radnz(nzrand))
      call slarnv(1, iseed_rad, nzrand, radnz)
      radnz = radnz * depz(nz-1)

      ! -------------------------------
      ! Separate Voronoi parameterizations:
      !   Vs uses (ncell, acell)
      !   Azimuthal anisotropy (Gc,Gs) uses reduced complexity
      ! -------------------------------
      ncell_a = max(1, int(iaratio*real(ncell)))
      acell_a = max(0, int(iaratio*real(acell)))

      allocate(ncell_v_l(nzrand), acell_v_l(nzrand))
      allocate(ncell_a_l(nzrand), acell_a_l(nzrand))
      allocate(off_nv_l(nzrand), off_av_l(nzrand))
      allocate(off_na_l(nzrand), off_aa_l(nzrand))

      do ii = 1, nzrand
        zfrac = radnz(ii) / max(1.0, depz(nz-1))

        ! Depth taper: shallow bands get more cells, weakly resolved deep
        ! bands get fewer.  The default thresholds are shallow and were used
        ! for the SRL synthetic tests.  The yunnan modes use relative depth
        ! breaks for the 0--120 km Yunnan model so the taper does not suppress
        ! all anisotropy below 12 km.
        if (taper_yunnan) then
          if (zfrac <= 0.20) then
            f_vs = 1.20
            f_ani = 1.20
          elseif (zfrac <= 0.45) then
            f_vs = 1.00
            f_ani = 0.85
          elseif (zfrac <= 0.70) then
            f_vs = 0.75
            f_ani = 0.45
          else
            f_vs = 0.55
            f_ani = 0.25
          endif
        elseif (zfrac <= 4.0/max(1.0, depz(nz-1))) then
          f_vs = 1.40
          f_ani = 1.80
        elseif (zfrac <= 8.0/max(1.0, depz(nz-1))) then
          f_vs = 1.15
          f_ani = 1.20
        elseif (zfrac <= 12.0/max(1.0, depz(nz-1))) then
          f_vs = 0.80
          f_ani = 0.55
        else
          f_vs = 0.45
          f_ani = 0.15
        endif
        if (taper_deep_suppress) then
          if (taper_yunnan .and. zfrac > 0.70) then
            f_ani = min(f_ani, 0.15)
          elseif ((.not. taper_yunnan) .and. zfrac > 8.0/max(1.0, depz(nz-1))) then
            f_ani = min(f_ani, 0.20)
          endif
        endif

        if (taper_enabled .and. taper_vs) then
          ncell_v_l(ii) = max(3, int(real(ncell)*f_vs + 0.5))
          acell_v_l(ii) = max(0, int(real(acell)*f_vs + 0.5))
        else
          ncell_v_l(ii) = ncell
          acell_v_l(ii) = acell
        endif

        if (taper_enabled .and. taper_ani) then
          ncell_a_l(ii) = max(1, int(real(ncell_a)*f_ani + 0.5))
          acell_a_l(ii) = max(0, int(real(acell_a)*f_ani + 0.5))
        else
          ncell_a_l(ii) = ncell_a
          acell_a_l(ii) = acell_a
        endif
        if (.not. use_adaptive_cells) then
          acell_v_l(ii) = 0
          acell_a_l(ii) = 0
        endif
      enddo

      ! Vs Voronoi
      ncells_v = 0
      do ii = 1, nzrand
        off_nv_l(ii) = ncells_v + 1
        ncells_v = ncells_v + ncell_v_l(ii)
      enddo
      do ii = 1, nzrand
        off_av_l(ii) = ncells_v + 1
        ncells_v = ncells_v + acell_v_l(ii)
      enddo

      ! Anisotropy Voronoi (reduced complexity; add adaptive cells whenever available)
      ncells_a = 0
      do ii = 1, nzrand
        off_na_l(ii) = ncells_a + 1
        ncells_a = ncells_a + ncell_a_l(ii)
      enddo
      do ii = 1, nzrand
        off_aa_l(ii) = ncells_a + 1
        ncells_a = ncells_a + acell_a_l(ii)
      enddo

      if (iproj == 1 .and. iter == 1 .and. taper_enabled) then
        open(41, file='vorotomo_depth_taper_counts_iter1.txt')
        write(41,'(A)') '# mode depth_km vs_normal vs_adaptive ani_normal ani_adaptive'
        do ii = 1, nzrand
          write(41,'(A,1X,F8.3,4(1X,I8))') trim(adjustl(taper_mode)), radnz(ii), &
               ncell_v_l(ii), acell_v_l(ii), ncell_a_l(ii), acell_a_l(ii)
        enddo
        close(41)
        write(*,'(A,I0,A,I0,A,I0,A,I0)') 'Depth-taper Voronoi totals: Vs normal=', sum(ncell_v_l), &
             ' Vs adaptive=', sum(acell_v_l), ' ANI normal=', sum(ncell_a_l), ' ANI adaptive=', sum(acell_a_l)
      endif

      ! Total unknowns: Vs uses its own cells, Gc and Gs share a (coarser) cell set
      ncells = ncells_v

      allocate(theta_v(ncells_v),phi_v(ncells_v),rrad_v(ncells_v))
      allocate(xpts_v(ncells_v),ypts_v(ncells_v),zpts_v(ncells_v))

      allocate(theta_a(ncells_a),phi_a(ncells_a),rrad_a(ncells_a))
      allocate(xpts_a(ncells_a),ypts_a(ncells_a),zpts_a(ncells_a))

      allocate(dis(max(ncells_v,ncells_a)))
      allocate(xunknown(ncells_v + 2*ncells_a))
      allocate(norm(ncells_v + 2*ncells_a))

      ! Mapping from regular grid points to Voronoi cells (Vs and aniso)
      allocate(row_v(ndim0), row_a(ndim0), col(ndim0))

      ! -------------------------------
      ! Vs normal cells: independent random (lat,lon) per layer (no replication)
      do ii = 1, nzrand
          cnt_layer = ncell_v_l(ii)
          layer_start = off_nv_l(ii)
          iseed_theta_l = iseed_theta
          iseed_phi_l   = iseed_phi

          ! decorrelate per layer deterministically
          iseed_theta_l(4) = iseed_theta_l(4) + 97 * ii
          iseed_phi_l(4)   = iseed_phi_l(4)   + 193 * ii
          iseed_theta_l(4) = mod(iseed_theta_l(4), 4095); if (iseed_theta_l(4)==0) iseed_theta_l(4)=1
          iseed_phi_l(4)   = mod(iseed_phi_l(4),   4095); if (iseed_phi_l(4)==0)   iseed_phi_l(4)=1

          call slarnv(1, iseed_theta_l, cnt_layer, theta_v(layer_start:layer_start+cnt_layer-1))
          theta_v(layer_start:layer_start+cnt_layer-1) = (gozd + theta_v(layer_start:layer_start+cnt_layer-1) * (ny-3) * dvzd) * pi/180

          call slarnv(1, iseed_phi_l, cnt_layer, phi_v(layer_start:layer_start+cnt_layer-1))
          phi_v(layer_start:layer_start+cnt_layer-1) = pi/2 - (goxd - phi_v(layer_start:layer_start+cnt_layer-1) * (nx-3) * dvxd) * pi/180
          rrad_v(layer_start:layer_start+cnt_layer-1) = radius - radnz(ii) * hvratio
      enddo

      ! -------------------------------
      ! Anisotropy normal cells: reduced complexity, independent random (lat,lon) per layer
      ! Use separate seeds to decorrelate from Vs (do NOT mutate iseed_theta/iseed_phi)
      iseed_theta_a = iseed_theta
      iseed_phi_a   = iseed_phi
      iseed_theta_a(4) = mod(iseed_theta_a(4) + 101, 4095); if (iseed_theta_a(4)==0) iseed_theta_a(4)=1
      iseed_phi_a(4)   = mod(iseed_phi_a(4)   + 103, 4095); if (iseed_phi_a(4)==0)   iseed_phi_a(4)=1

      do ii = 1, nzrand
          cnt_layer = ncell_a_l(ii)
          layer_start = off_na_l(ii)
          iseed_theta_a_l = iseed_theta_a
          iseed_phi_a_l   = iseed_phi_a

          ! decorrelate per layer deterministically
          iseed_theta_a_l(4) = iseed_theta_a_l(4) + 149 * ii
          iseed_phi_a_l(4)   = iseed_phi_a_l(4)   + 251 * ii
          iseed_theta_a_l(4) = mod(iseed_theta_a_l(4), 4095); if (iseed_theta_a_l(4)==0) iseed_theta_a_l(4)=1
          iseed_phi_a_l(4)   = mod(iseed_phi_a_l(4),   4095); if (iseed_phi_a_l(4)==0)   iseed_phi_a_l(4)=1

          call slarnv(1, iseed_theta_a_l, cnt_layer, theta_a(layer_start:layer_start+cnt_layer-1))
          theta_a(layer_start:layer_start+cnt_layer-1) = (gozd + theta_a(layer_start:layer_start+cnt_layer-1) * (ny-3) * dvzd) * pi/180

          call slarnv(1, iseed_phi_a_l, cnt_layer, phi_a(layer_start:layer_start+cnt_layer-1))
          phi_a(layer_start:layer_start+cnt_layer-1) = pi/2 - (goxd - phi_a(layer_start:layer_start+cnt_layer-1) * (nx-3) * dvxd) * pi/180
          rrad_a(layer_start:layer_start+cnt_layer-1) = radius - radnz(ii) * hvratio
      enddo

      !separate the big G matrix to Giso and Gc Gs-------------------!
      call range_idx_1d(idx_iso, colg, 0, ndim0+1)
      lenrw_iso = int(size(idx_iso), kind=8)
      if (lenrw_iso <= 0_8) then
          write(*,*) 'ERROR: empty iso-slice in voronoiproj_layer_azi (lenrw_iso=0). Check colg indexing.'
          stop
      endif
      allocate(colg_iso(int(lenrw_iso,kind=4)), rw_iso(int(lenrw_iso,kind=4)))
      colg_iso = colg(int(idx_iso,kind=4))
      rw_iso   = rw(int(idx_iso,kind=4))

      allocate(dws_iso(ndim0))

      ! adaptive cells based on data support (computed from iso-slice weights)
      dws_iso = 0.0
      do i = 1, lenrw_iso
        dws_iso(colg_iso(i)) = dws_iso(colg_iso(i)) + abs(rw_iso(i))
      enddo

      ! -------- Vs adaptive cells
      if (use_adaptive_cells) then
        if (allocated(dws_v)) deallocate(dws_v)
        allocate(dws_v(size(dws_iso)))
        dws_v = dws_iso

        do jj = 1, nzrand
          do ii = 1, acell_v_l(jj)
            call random_index_lapack(idx, dws_v, iseed_idx)
            dws_v(idx) = 0.0

            iz = (idx-1)/((nx-2)*(ny-2)) + 1
            iy = (idx-1 + (1-iz)*(nx-2)*(ny-2)) / (nx-2) + 1
            ix = idx + (nx-2)*(1-iy) + (nx-2)*(ny-2)*(1-iz)

            layer_start = off_av_l(jj) + ii - 1
            theta_v(layer_start) = (gozd + (iy-1)*dvzd) * pi/180
            phi_v(layer_start)   = pi/2 - (goxd - (ix-1)*dvxd) * pi/180
            rrad_v(layer_start)  = radius - radnz(jj) * hvratio
          enddo
        enddo
        deallocate(dws_v)
      endif

      ! -------- Anisotropy adaptive cells (always added when available)
      if (sum(acell_a_l) > 0) then
        if (allocated(dws_a)) deallocate(dws_a)
        allocate(dws_a(size(dws_iso)))
        dws_a = dws_iso

        do jj = 1, nzrand
          do ii = 1, acell_a_l(jj)
            call random_index_lapack(idx, dws_a, iseed_idx)
            dws_a(idx) = 0.0

            iz = (idx-1)/((nx-2)*(ny-2)) + 1
            iy = (idx-1 + (1-iz)*(nx-2)*(ny-2)) / (nx-2) + 1
            ix = idx + (nx-2)*(1-iy) + (nx-2)*(ny-2)*(1-iz)

            layer_start = off_aa_l(jj) + ii - 1
            theta_a(layer_start) = (gozd + (iy-1)*dvzd) * pi/180
            phi_a(layer_start)   = pi/2 - (goxd - (ix-1)*dvxd) * pi/180
            rrad_a(layer_start)  = radius - radnz(jj) * hvratio
          enddo
        enddo
        deallocate(dws_a)
      endif

      xpts_v = rrad_v*sin(phi_v)*cos(theta_v)
      ypts_v = rrad_v*sin(phi_v)*sin(theta_v)
      zpts_v = rrad_v*cos(phi_v)

      xpts_a = rrad_a*sin(phi_a)*cos(theta_a)
      ypts_a = rrad_a*sin(phi_a)*sin(theta_a)
      zpts_a = rrad_a*cos(phi_a)

      if (iproj == 1 .and. use_adaptive_cells .and. iter == 1) then
        open(77,file='pvpoints_acell.dat')
        do ii=1,ncells_v
          write(77,'(3f10.4)') rrad_v(ii),phi_v(ii),theta_v(ii)
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

        open(37,file='vorosites_latlon.txt')
        do ii = 1,ncells_v
        write(37, '(3f10.2)') theta_v(ii)*180/pi, (pi/2-phi_v(ii))*180/pi, (radius-rrad_v(ii))/hvratio
        !write(36, '(3f10.2)') atan(ypts(ii)/xpts(ii))*180/pi, acos(zpts(ii)/sqrt(xpts(ii)**2+ypts(ii)**2+zpts(ii)**2))*180/pi, &
        !    (radius-sqrt(xpts(ii)**2+ypts(ii)**2+zpts(ii)**2))/hvratio
        enddo
        close(37)
      endif

      if (iproj == 1) then
        write(citer,'(I0)') iter
        open(37, file='vorosites_latlon_vs_iter'//trim(citer)//'.txt')
        do ii = 1, ncells_v
            write(37, '(3f10.2)') theta_v(ii)*180/pi, (pi/2-phi_v(ii))*180/pi, (radius-rrad_v(ii))/hvratio
        enddo
        close(37)

        open(38, file='vorosites_latlon_aniso_iter'//trim(citer)//'.txt')
        do ii = 1, ncells_a
            write(38, '(3f10.2)') theta_a(ii)*180/pi, (pi/2-phi_a(ii))*180/pi, (radius-rrad_a(ii))/hvratio
        enddo
        close(38)
      endif

      if (iproj == 1 .and. use_adaptive_cells .and. iter == 1) then
        open(39,file='col.txt')
      endif
      if (timer_enabled) call system_clock(t_cells)
      idx = 0
      do iz = 1,nz-1
        do iy = 1,ny-2
          do ix = 1,nx-2
            idx = idx + 1
            xs = rad(iz)*sin(pi/2-lat(ix))*cos(lon(iy))
            ys = rad(iz)*sin(pi/2-lat(ix))*sin(lon(iy))
            zs = rad(iz)*cos(pi/2-lat(ix))

            ! Vs cells
            dis(1:ncells_v) = (xpts_v-xs)**2 + (ypts_v-ys)**2 + (zpts_v-zs)**2
            do jj = 1, ncells_v
              dis(jj) = dis(jj) + eps_tie * real(jj)
            enddo
            row_v(idx) = minloc(dis(1:ncells_v), 1)
            if (row_v(idx) < 1) row_v(idx) = 1
            if (row_v(idx) > ncells_v) row_v(idx) = ncells_v

            ! Anisotropy cells
            dis(1:ncells_a) = (xpts_a-xs)**2 + (ypts_a-ys)**2 + (zpts_a-zs)**2
            do jj = 1, ncells_a
              dis(jj) = dis(jj) + eps_tie * real(jj)
            enddo
            row_a(idx) = minloc(dis(1:ncells_a), 1)
            if (row_a(idx) < 1) row_a(idx) = 1
            if (row_a(idx) > ncells_a) row_a(idx) = ncells_a

            col(idx) = (iz-1)*(nx-2)*(ny-2)+(iy-1)*(nx-2)+ix
            if (iproj == 1 .and. use_adaptive_cells .and. iter == 1) then
              write(39,*) ix, iy, iz, col(idx), lat(ix)*180/pi, lon(iy)*180/pi, dws_iso(col(idx))
            endif
          enddo
        enddo
      enddo
      if (iproj == 1 .and. use_adaptive_cells .and. iter == 1) then
        close(39)
      endif
      ! Allocate mapping arrays and projection structures before use

      allocate(rw_p(ndim0))
      allocate(iw_p(2*ndim0+1))

      allocate(rw_p_3(ndim))
      allocate(iw_p_3(2*ndim+1))

      rw_p = 1.0
      leniw_p = 2_8*int(ndim0,kind=8)
      lenrw_p = int(ndim0,kind=8)
      iw_p(1:ndim0) = row_v
      iw_p(ndim0+1:2*ndim0) = col
      !---------------------------------------------------------------!
      !compose a bigger projection matrix 
      ! P is ndim0 * ncells matrix
      ! Big P is ndim (ndim0*3) * (3*ncells)
      !|P 0 0|
      !|0 P 0|
      !|0 0 P|
      rw_p_3 = 1.0
      leniw_p_3 = 2_8*int(ndim,kind=8)
      lenrw_p_3 = int(ndim,kind=8)

      ! Big projection matrix: [P_v 0 0; 0 P_a 0; 0 0 P_a]
      iw_p_3(1:ndim0) = row_v
      iw_p_3(ndim0+1:2*ndim0) = row_a + ncells_v
      iw_p_3(2*ndim0+1:3*ndim0) = row_a + ncells_v + ncells_a

      iw_p_3(3*ndim0+1:4*ndim0) = col
      iw_p_3(4*ndim0+1:5*ndim0) = col + ndim0
      iw_p_3(5*ndim0+1:6*ndim0) = col + 2*ndim0
      !---------------------------------------------------------------!
      if (timer_enabled) call system_clock(t_map)

      ! Set nunk only once, at top
      nunk = ncells_v + 2*ncells_a
      gp_spfra_eff = spfra
      gp_spfra_env = ''
      call get_environment_variable('DAZI_VORO_GP_SPFRA', gp_spfra_env, status=gp_spfra_status)
      if (gp_spfra_status == 0) then
        read(gp_spfra_env,*,iostat=gp_spfra_status) gp_spfra_value
        if (gp_spfra_status == 0 .and. gp_spfra_value > 0.0_8) gp_spfra_eff = gp_spfra_value
      endif

      allocate(grow(ndim),gcol(nd),subrow(nunk))
      ! maxnar = int(0.6*nd*ncells)
      maxnar = int(gp_spfra_eff * real(nd,kind=8) * real(nunk,kind=8), kind=8)
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

      nzid = 0
      if (trim(adjustl(gp_backend)) == 'direct') then
        start = 0_8
        do ii = 1,nd
        subrow = 0
        do ix = 1,nrow(ii)
          model_id = colg(start+ix)
          if (model_id >= 1 .and. model_id <= ndim0) then
            if (freeze_vs) cycle
            cell_id = row_v(model_id)
          elseif (model_id > ndim0 .and. model_id <= 2*ndim0) then
            grid_id = model_id - ndim0
            cell_id = ncells_v + row_a(grid_id)
          elseif (model_id > 2*ndim0 .and. model_id <= ndim) then
            grid_id = model_id - 2*ndim0
            cell_id = ncells_v + ncells_a + row_a(grid_id)
          else
            stop 'DAZI_VORO_GP_BACKEND=direct: model column out of range'
          endif
          subrow(cell_id) = subrow(cell_id) + rw(start+ix)
        enddo
        do ix = 1, nunk
        if(abs(subrow(ix))>ftol) then
          if(nzid .ge. maxnar) then
            write(6,*) 'Voronoi projection part sparsefrac is too small, please increase your sparsefrac'
            stop 'nzid > maxnar, please increase your sparsefrac'
          endif
          nzid = nzid+1
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
        if (freeze_vs) grow(1:ndim0) = 0.0
        !gcol = 0
        !gcol(ii) = 1.0 
        !call aprod(2,nd,ndim,grow,gcol,leniw,lenrw,iw,rw)
        subrow = 0
        call aprod(1,nunk,ndim,grow,subrow,leniw_p_3,lenrw_p_3,iw_p_3,rw_p_3)
        do ix = 1, nunk
        if(abs(subrow(ix))>ftol) then
          if(nzid .ge. maxnar) then
            write(6,*) 'Voronoi projection part sparsefrac is too small, please increase your sparsefrac'
            stop 'nzid > maxnar, please increase your sparsefrac'
          endif
          nzid = nzid+1
          rwgp(nzid) = subrow(ix)
          iwgp(nzid) = ii
          colgp(nzid) = ix
        endif
        enddo
        start = start + int(nrow(ii), kind=8)
        enddo
      endif
      if (nzid.gt.maxnar) then
        write(6,*) 'Voronoi projection part sparsefrac is too small, please increase your sparsefrac'
        stop 'nzid > maxnar, please increase your sparsefrac'
      endif
      leniwgp = 2_8*nzid
      lenrwgp = nzid

      ! ---- Pack into LSMR triplet header layout ----
      ! iwgp(1:nnz)         = row indices (already stored during assembly as iwgp(1+k))
      ! iwgp(nnz+1:2*nnz)   = col indices
      iwgp(int(nzid,kind=4)+1 : 2*int(nzid,kind=4)) = colgp(1:int(nzid,kind=4))
      if (timer_enabled) call system_clock(t_gp)

      ! norm = 0
      ! do ii=1,nzid
      ! norm(iwgp(ii+nzid)) = norm(iwgp(ii+nzid))+rwgp(ii)**2
      ! enddo
    
      ! do ii =1,3*ncells
      ! norm(ii) = sqrt(norm(ii)/nd+0.01)
      ! enddo
    
      ! do ii =1,nzid
      ! rwgp(ii) = rwgp(ii)/norm(iwgp(1+ii+nzid))
      ! enddo

      conlim = 50
      itnlim = 100
      atol = 1e-3/((dvxd+dvzd)*111.19/2.0*0.1) !1e-2
      btol = 1e-3/(dvxd*nx*111.19/3.0)!1e-3
      istop = 0
      anorm = 0.0
      acond = 0.0
      arnorm = 0.0
      xnorm = 0.0
      localSize = max(1, int(nunk/4))
      !damp = dampvel
      ! using lsmr to solve for the projection coefficients
      !print*, 'LSMR beginning ...'

      nout = -1
      !nout = 36
      !open(nout,file='lsmrout_sub.txt')
      xunknown = 0.0

      lsmr_backend = ''
      call get_environment_variable('DAZI_VORO_LSMR_BACKEND', lsmr_backend, status=env_status)
      if (env_status /= 0) then
        call get_environment_variable('DAZI_LSMR_BACKEND', lsmr_backend, status=env_status)
      endif

      if (trim(adjustl(lsmr_backend)) == 'scipy') then
        if (timer_enabled) call system_clock(t_lsmr0)
        write(lsmr_tag,'("voro_azi_i",I0,"_r",I0)') iter, iproj
        call scipy_lsmr_solve(lsmr_tag, nd, nunk, leniwgp, lenrwgp, iwgp, rwgp, dres, damp, &
        atol, btol, conlim, itnlim, xunknown, istop, itn, anorm, acond, rnorm, arnorm, xnorm)
      else
        if (timer_enabled) call system_clock(t_lsmr0)
        call LSMR(nd, nunk, leniwgp, lenrwgp,iwgp,rwgp,dres,damp,&
        atol, btol, conlim, itnlim, localSize,nout,&
        xunknown, istop, itn, anorm, acond,rnorm, arnorm, xnorm)
      endif
      if (timer_enabled) call system_clock(t_lsmr1)
      !close(nout)
      ! do ii = 1,ncells
      !   xunknown(ii) = xunknown(ii)/norm(ii)
      ! enddo
      ! norm = (norm**2-0.01)*nd
      ! do ii = 1,3*ncells
      ! if (norm(ii)<0.01) then
      !     call random_number(rx)
      !     xunknown(ii) = xunknown(ii)+rx-0.5
      ! endif
      ! enddo

      dv(1:ndim) = 0
      call aprod(2,nunk,ndim,dv,xunknown,leniw_p_3,lenrw_p_3,iw_p_3,rw_p_3)
      if (freeze_vs) dv(1:ndim0) = 0.0
      if (timer_enabled) then
        call system_clock(t_end)
        inquire(unit=66, opened=timer_log_open)
        write(*,'(A,I0,A,I0,A,5(F10.4,A))') 'Voro timer azi iter=', iter, ' real=', iproj, &
          ' cells=', real(t_cells-t_start)/real(clock_rate), 's map=', &
          real(t_map-t_cells)/real(clock_rate), 's gp=', &
          real(t_gp-t_map)/real(clock_rate), 's lsmr=', &
          real(t_lsmr1-t_lsmr0)/real(clock_rate), 's back=', &
          real(t_end-t_lsmr1)/real(clock_rate), 's'
        if (timer_log_open) then
          write(66,'(A,I0,A,I0,A,5(F10.4,A))') 'Voro timer azi iter=', iter, ' real=', iproj, &
            ' cells=', real(t_cells-t_start)/real(clock_rate), 's map=', &
            real(t_map-t_cells)/real(clock_rate), 's gp=', &
            real(t_gp-t_map)/real(clock_rate), 's lsmr=', &
            real(t_lsmr1-t_lsmr0)/real(clock_rate), 's back=', &
            real(t_end-t_lsmr1)/real(clock_rate), 's'
        endif
      endif
      rnorm_out = rnorm
      acond_out = acond
      itn_out = itn
      istop_out = istop
      deallocate(grow,gcol,subrow)
      deallocate(theta_v,phi_v,rrad_v)
      deallocate(theta_a,phi_a,rrad_a)
      deallocate(xpts_v,ypts_v,zpts_v)
      deallocate(xpts_a,ypts_a,zpts_a)
      deallocate(dis,xunknown)
      if (allocated(iw_p)) deallocate(iw_p)
      if (allocated(row_v)) deallocate(row_v)
      if (allocated(row_a)) deallocate(row_a)
      if (allocated(col)) deallocate(col)
      deallocate(lat,lon,rad)
      deallocate(iwgp,colgp,rwgp)
      deallocate(rw_p_3)
      if (allocated(rw_p)) deallocate(rw_p)
      deallocate(iw_p_3)
      deallocate(colg_iso,idx_iso)
      deallocate(dws_iso)
      if (allocated(dws_v)) deallocate(dws_v)
      if (allocated(dws_a)) deallocate(dws_a)
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

subroutine range_idx_1d( idx, arr, value1, value2 )
    !no multi-modal, only for 1d array index
    implicit none
    integer(kind=8),dimension(:),allocatable :: idx(:)
    integer, intent(in) :: arr(:)
    integer(kind=8) :: i, count1
    integer :: value1, value2

    allocate(idx(size(arr)))
    count1 = 0_8
    do i = 1_8, int(size(arr), kind=8)
        if ( arr(int(i,kind=4)) < value2 .and. arr(int(i,kind=4)) > value1 ) then
            count1 = count1 + 1_8
            idx(count1) = i
        endif
    enddo

    if (count1 > 0_8) then
        idx = idx(1:int(count1,kind=4))
    else
        deallocate(idx)
        allocate(idx(0))
    endif
end subroutine range_idx_1d

end subroutine
