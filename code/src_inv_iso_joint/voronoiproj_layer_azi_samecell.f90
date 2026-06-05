subroutine voronoiproj_layer_azi(leniw,lenrw,colg,nrow,rw,dres,goxd,dvxd,gozd,dvzd,depz,&
                           nx,ny,nz,nd,ncell,acell,nzrand,spfra,hvratio,damp,iproj,iter,dv,&
                           rnorm_out, acond_out, itn_out, istop_out)
      use lsmrModule, only:lsmr

      implicit none
      integer leniw,lenrw
      integer lenrw_iso
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
      integer ndim,nd,ndim0
      integer iproj,iter,nzrand
      ! Output arguments for LSMR quality metrics
      real, intent(out) :: rnorm_out, acond_out
      integer, intent(out) :: itn_out, istop_out
      real,parameter:: radius = 6371.0,ftol = 0.0001,pi = 3.141592654
      integer ii,jj,ix,iy,iz
      real,dimension(:),allocatable:: grow,gcol,subrow,dis,xunknown
      real,dimension(:),allocatable:: lat,lon,rad,theta,phi,rrad,xpts,ypts,zpts
      real,dimension(:),allocatable :: rw_p,rwgp,norm
      real,dimension(:),allocatable:: radnz
      integer,dimension(:),allocatable:: iw_p,row,col,iwgp,colgp
      integer idx

      real,dimension(:),allocatable :: rw_p_3
      real,dimension(:),allocatable:: rw_iso, dws_iso
      integer,dimension(:),allocatable:: idx_iso, colg_iso
      integer,dimension(:),allocatable:: iw_p_3

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
      integer leniw_p_3,lenrw_p_3
      integer start
      ! eps scaled by typical magnitude of dis (units: km^2)
      real, parameter :: eps_tie = 1.0e-12

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
      acells = acell*nzrand
      athresh = 10
      if(acell>athresh) then
        ! write(*,*) 'Adding AZI cells!'
        ncells = ncell*nzrand+acell*nzrand
      else
        write(*,*) 'Too few adaptive cells, not adding them!'
        ncells = ncell*nzrand
      endif

      allocate(theta(ncells),phi(ncells),rrad(ncells),norm(3*ncells))
      allocate(xpts(ncells),ypts(ncells),zpts(ncells),dis(ncells),xunknown(3*ncells))
      allocate(rw_p(ndim0)) ! the projection matrix
      allocate(iw_p(2*ndim0+1),row(ndim0),col(ndim0),dws_iso(ndim0))
      allocate(rw_p_3(ndim)) ! the projection matrix
      allocate(iw_p_3(2*ndim+1))
      ! allocate(colg_iso(ndim0))

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
      iseed_idx(4)   = mod(iseed_idx(4),   4095); if (iseed_idx(4)==0)   iseed_idx(4)=1      ! Normal cells: same (lat,lon) pattern replicated across nzrand layers
      call slarnv(1, iseed_theta, ncell, theta(1:ncell))
      theta(1:ncell) = (gozd + theta(1:ncell) * (ny-3) * dvzd) * pi/180
      do ii = 1, nzrand-1
          theta(ii*ncell+1:(ii+1)*ncell) = theta(1:ncell)
      enddo

      call slarnv(1, iseed_phi, ncell, phi(1:ncell))
      phi(1:ncell) = pi/2 - (goxd - phi(1:ncell) * (nx-3) * dvxd) * pi/180
      do ii = 1, nzrand-1
          phi(ii*ncell+1:(ii+1)*ncell) = phi(1:ncell)
      enddo

      allocate(radnz(nzrand))
      call slarnv(1, iseed_rad, nzrand, radnz)
      radnz = radnz * depz(nz-1)

      do ii = 1, nzrand
          do jj = 1, ncell
              rrad((ii-1)*ncell + jj) = radius - radnz(ii) * hvratio
          enddo
      enddo

      !separate the big G matrix to Giso and Gc Gs-------------------!
      call range_idx_1d(idx_iso, colg, 0, ndim0+1)
      lenrw_iso = size(idx_iso)
      if (lenrw_iso <= 0) then
          write(*,*) 'ERROR: empty iso-slice in voronoiproj_layer_azi (lenrw_iso=0). Check colg indexing.'
          stop
      endif
      allocate(colg_iso(lenrw_iso), rw_iso(lenrw_iso))
      colg_iso = colg(idx_iso)
      rw_iso   = rw(idx_iso)

      ! write(*,*) 'iproj=',iproj,'iter=',iter

      ! write(*,*) 'Here 4!'
      ! adaptive cells based on dws, assume 1/2 of all ncells are used
      ! as adaptive cells
      if(acell >athresh) then
        ! write(*,*) 'Adding AZI cells!'
        dws_iso = 0
        do ii = 1,lenrw_iso
          dws_iso(colg_iso(ii)) = dws_iso(colg_iso(ii))+abs(rw_iso(ii))
        enddo

        do ii = 1,acell
          call random_index_lapack(idx, dws_iso, iseed_idx)
          dws_iso(idx) = 0.0  ! reduce the probability of being selected again
          ! ix = mod(idx,nx-2)
          ! iy = idx/(nx-2)
          ! iz = idx/((nx-2)*(ny-2))

          iz = (idx-1)/((nx-2)*(ny-2))+1  ! Ying
          iy = (idx-1+(1-iz)*(nx-2)*(ny-2))/(nx-2)+1  ! Ying
          ix = idx+(nx-2)*(1-iy)+(nx-2)*(ny-2)*(1-iz)  ! Ying

          do jj = 1,nzrand
          theta(nzrand*ncell+(ii-1)*nzrand+jj) = (gozd+(iy-1)*dvzd)*pi/180
          phi(nzrand*ncell+(ii-1)*nzrand+jj) = pi/2-(goxd-(ix-1)*dvxd)*pi/180
          ! rrad(nzrand*ncell+(ii-1)*nzrand+jj) = radius-depz(iz+1)*hvratio
          rrad(nzrand*ncell+(ii-1)*nzrand+jj) = radius-radnz(jj)*hvratio ! because the depz(0) is not the boundary
          enddo
        enddo
      endif

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

        open(37,file='vorosites_latlon.txt')
        do ii = 1,ncells
        write(37, '(3f10.2)') theta(ii)*180/pi, (pi/2-phi(ii))*180/pi, (radius-rrad(ii))/hvratio
        !write(36, '(3f10.2)') atan(ypts(ii)/xpts(ii))*180/pi, acos(zpts(ii)/sqrt(xpts(ii)**2+ypts(ii)**2+zpts(ii)**2))*180/pi, &
        !    (radius-sqrt(xpts(ii)**2+ypts(ii)**2+zpts(ii)**2))/hvratio
        enddo
        close(37)
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
                        ! eps_tie scaled by typical magnitude of dis (units: km^2)
            do jj = 1, ncells
                dis(jj) = dis(jj) + eps_tie * real(jj)
            enddo
            idx = idx+1
            col(idx) = (iz-1)*(nx-2)*(ny-2)+(iy-1)*(nx-2)+ix
            row(idx) = minloc(dis,1)
            if (row(idx) < 1) row(idx) = 1
            if (row(idx) > ncells) row(idx) = ncells
            if (iproj == 1 .and. acell > athresh .and. iter == 1) then
              write(38,*) ix, iy, iz, col(idx), lat(ix)*180/pi, lon(iy)*180/pi, dws_iso(col(idx))
            endif
          enddo
        enddo
      enddo
      if (iproj == 1 .and. acell > athresh .and. iter == 1) then
        close(38)
      endif

! write(*,*) 'Here 5!'
      rw_p = 1.0
      leniw_p = 2*ndim0+1
      lenrw_p = ndim0
      iw_p(1) = ndim0
      iw_p(2:ndim0+1) = row
      iw_p(ndim0+2:2*ndim0+1) = col
      !---------------------------------------------------------------!
      !compose a bigger projection matrix 
      ! P is ndim0 * ncells matrix
      ! Big P is ndim (ndim0*3) * (3*ncells)
      !|P 0 0|
      !|0 P 0|
      !|0 0 P|
      rw_p_3 = 1.0
      leniw_p_3 = 2*ndim+1
      lenrw_p_3 = ndim
      iw_p_3(1) = ndim

      ! iw_p_3(2:ndim0+1) = row
      ! iw_p_3(ndim0+2:2*ndim0+1) = row + ndim0
      ! iw_p_3(2*ndim0+2:3*ndim0+1) = row + 2*ndim0
      ! iw_p_3(3*ndim0+2:4*ndim0+1) = col
      ! iw_p_3(4*ndim0+2:5*ndim0+1) = col + ncells
      ! iw_p_3(5*ndim0+2:6*ndim0+1) = col + 2*ncells

      iw_p_3(2:ndim0+1) = row
      iw_p_3(ndim0+2:2*ndim0+1) = row + ncells
      iw_p_3(2*ndim0+2:3*ndim0+1) = row + 2*ncells
      iw_p_3(3*ndim0+2:4*ndim0+1) = col
      iw_p_3(4*ndim0+2:5*ndim0+1) = col + ndim0
      iw_p_3(5*ndim0+2:6*ndim0+1) = col + 2*ndim0
      !---------------------------------------------------------------!

      allocate(grow(ndim),gcol(nd),subrow(3*ncells))
      ! maxnar = int(0.6*nd*ncells)
      maxnar = int(spfra*nd*ncells*3)
      if((iproj==1)) then
        write(*,*) 'Max Number of sparse matrix (GP) non-zero element ',maxnar
      endif
      allocate(iwgp(maxnar*2+1),colgp(maxnar),rwgp(maxnar))

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
      call aprod(1,3*ncells,ndim,grow,subrow,leniw_p_3,lenrw_p_3,iw_p_3,rw_p_3)
      do ix = 1,3*ncells
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

      ! norm = 0
      ! do ii=1,nzid
      ! norm(iwgp(1+ii+nzid)) = norm(iwgp(1+ii+nzid))+rwgp(ii)**2
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
      localSize = int(ncells/4)
      !damp = dampvel
      ! using lsmr to solve for the projection coefficients
      !print*, 'LSMR beginning ...'

      nout = -1
      !nout = 36
      !open(nout,file='lsmrout_sub.txt')

      call LSMR(nd, 3*ncells, leniwgp, lenrwgp,iwgp,rwgp,dres,damp,&
      atol, btol, conlim, itnlim, localSize,nout,&
      xunknown, istop, itn, anorm, acond,rnorm, arnorm, xnorm)
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
      call aprod(2,3*ncells,ndim,dv,xunknown,leniw_p_3,lenrw_p_3,iw_p_3,rw_p_3)
      rnorm_out = rnorm
      acond_out = acond
      itn_out = itn
      istop_out = istop
      deallocate(grow,gcol,subrow)
      deallocate(theta,phi,rrad,norm)
      deallocate(xpts,ypts,zpts,dis,xunknown)
      deallocate(iw_p,rw_p,row,col)
      deallocate(lat,lon,rad)
      deallocate(iwgp,colgp,rwgp)
      deallocate(rw_p_3)
      deallocate(iw_p_3)
      deallocate(colg_iso,idx_iso)
      deallocate(dws_iso)
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
    integer,dimension(:),allocatable :: idx(:)
    integer, intent(in) :: arr(:)
    integer :: i,count1
    integer :: value1, value2
    allocate(idx(size( arr )))
    count1 = 0
    do i = 1, size( arr )
        if ( arr(i) < value2 .and. arr(i) > value1) then
            count1 = count1+1
            idx(count1) = i
            ! write(*,'(i5)') idx
        endif
    enddo
    if (count1 > 0) then
        idx = idx(1:count1)
    else
        deallocate(idx)
        allocate(idx(0))
    endif
end subroutine range_idx_1d

end subroutine
