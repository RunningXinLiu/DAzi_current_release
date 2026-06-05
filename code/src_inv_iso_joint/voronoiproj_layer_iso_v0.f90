subroutine voronoiproj_layer_iso(leniw,lenrw,colg,nrow,rw,dres,goxd,dvxd,gozd,dvzd,depz,&
                       nx,ny,nz,nd,ncell,spfra,hvratio,damp,iproj,iter,dv)
      use lsmrModule, only:lsmr

      implicit none
      integer leniw,lenrw
      integer nx,ny,nz
!      integer iw(leniw)
      integer colg(lenrw),nrow(nd)
      real depz(nz)
      real rw(lenrw)
      integer ncell,acells,ncells,nzrand
      real dv(*),dres(*)
      real goxd,gozd,dvxd,dvzd
      real damp
      real hvratio,spfra,cmb
      integer ndim,nd
      integer iproj,iter

      real,parameter:: radius = 6371.0,ftol = 0.0001,pi = 3.141592654
      integer ii,ix,iy,iz,jj
      real,dimension(:),allocatable:: grow,gcol,subrow,dis,dws,xunknown
      real,dimension(:),allocatable:: lat,lon,rad,theta,phi,rrad,xpts,ypts,zpts,rad1,rad2
      real,dimension(:),allocatable :: rw_p,rwgp,norm
      integer,dimension(:),allocatable:: iw_p,row,col,iwgp,colgp,radnz
      integer idx
      integer maxnar,nzid
      integer iseed(4)
      real xs,ys,zs
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
      integer start,athresh

      integer time_begin,time_end ! Ying

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

      nzrand = int((nz-1)/2)
      ! Ying
      ! acells = int(ncell/3)
      acells = int(ncell/2)

      !==================Manipulating the voronoi cell depth============================!
      iseed(1:3) = (/38,62,346/)
      iseed(4) = 2*(iproj+iter)+1
      write(*,*) 'iproj=',iproj,'iter=',iter
      allocate(rad1(nzrand),rad2(nzrand),radnz(nzrand))
      call slarnv(1,iseed,nzrand,rad2)
      rad2 = rad2*(nz-1)
      radnz = ceiling(rad2)
      rad2 = depz(radnz)
      call unique_1d(rad1,rad2)
      nzrand = size(rad1)
      !==============================================!
      athresh = 10
      if(acells>athresh) then
          ncells = ncell*nzrand+acells*nzrand
      else
          write(*,*) 'Too few adaptive cells, not adding them!'
          ncells = ncell*nzrand
      endif
      allocate(theta(ncells),phi(ncells),rrad(ncells),norm(ncells))
      allocate(xpts(ncells),ypts(ncells),zpts(ncells),dis(ncells),xunknown(ncells))
      allocate(rw_p(ndim))
      allocate(iw_p(2*ndim+1),row(ndim),col(ndim),dws(ndim))


      call slarnv(1,iseed,ncell,theta)
      theta = (gozd+theta*(ny-3)*dvzd)*pi/180
      do ii = 1,nzrand-1
      theta(ii*ncell+1:(ii+1)*ncell) = theta(1:ncell)
      enddo
      call slarnv(1,iseed,ncell,phi)
      phi = pi/2-(goxd-phi*(nx-3)*dvxd)*pi/180
      do ii = 1,nzrand-1
      phi(ii*ncell+1:(ii+1)*ncell) = phi(1:ncell)
      enddo

      do ii = 1,nzrand
      do jj = 1,ncell
      !rrad((ii-1)*ncell+jj) = radius - (depz((ii-1)*2+1)+depz(ii*2))/2.0*hvratio
      rrad((ii-1)*ncell+jj) = radius - rad1(ii)*hvratio!(depz((ii-1)*2+1)+depz(ii*2))/2.0*hvratio
      enddo
      enddo

    if(iproj==1) then
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
      if (acells > athresh) then
      write(*,*),'Adjusting cells!'
      dws = 0
      do ii = 1,lenrw
      dws(colg(ii)) = dws(colg(ii))+abs(rw(ii))
      enddo
      !acells = int(ncell/2.0)
      do ii = 1, acells
      call random_index(idx,dws) 
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
      rrad(nzrand*ncell+(ii-1)*nzrand+jj) = radius-rad1(jj)*hvratio ! because the depz(0) is not the boundary
      ! rrad(nzrand*ncell+(ii-1)*nzrand+jj) = radius-depz(iz)*hvratio ! because the depz(0) is not the boundary
      enddo
      enddo
      endif

    if(iproj==1) then
        open(77,file='pvpoints_acell.dat')
        do ii=1,ncells
        write(77,'(3f10.4)') rrad(ii),phi(ii),theta(ii)
        enddo
        close(77)
        !print*,goxd,gozd,dvxd,dvzd,nx,ny
        !stop
    endif
      xpts = rrad*sin(phi)*cos(theta)
      ypts = rrad*sin(phi)*sin(theta)
      zpts = rrad*cos(phi)

      ! Ying output points for plot
      open(36,file='vorosites.txt')
      do ii = 1,ncells
      write(36, '(3f10.2)') xpts(ii), ypts(ii), zpts(ii)
      enddo
      close(36)

      open(37,file='vorosites_latlon.txt')
      do ii = 1,ncells
      write(37, '(3f10.2)') theta(ii)*180/pi, (pi/2-phi(ii))*180/pi, (radius-rrad(ii))/hvratio
      !write(36, '(3f10.2)') atan(ypts(ii)/xpts(ii))*180/pi, acos(zpts(ii)/sqrt(xpts(ii)**2+ypts(ii)**2+zpts(ii)**2))*180/pi, &
      !    (radius-sqrt(xpts(ii)**2+ypts(ii)**2+zpts(ii)**2))/hvratio
      enddo
      close(37)

      open(38,file='col.txt')
      idx = 0
      do iz = 1,nz-1
      do iy = 1,ny-2
      do ix = 1,nx-2
      xs = rad(iz)*sin(pi/2-lat(ix))*cos(lon(iy))
      ys = rad(iz)*sin(pi/2-lat(ix))*sin(lon(iy))
      zs = rad(iz)*cos(pi/2-lat(ix))
      dis =  (xpts-xs)**2+(ypts-ys)**2+(zpts-zs)**2
      idx = idx+1
      col(idx) = (iz-1)*(nx-2)*(ny-2)+(iy-1)*(nx-2)+ix
      row(idx) = minloc(dis,1)
      write(38,*) ix, iy, iz, col(idx), lat(ix)*180/pi, lon(iy)*180/pi, dws(idx)
      enddo
      enddo
      enddo
      close(38)

      rw_p = 1.0
      leniw_p = 2*ndim+1
      lenrw_p = ndim
      iw_p(1) = ndim
      iw_p(2:ndim+1) = row
      iw_p(ndim+2:2*ndim+1) = col

      allocate(grow(ndim),gcol(nd),subrow(ncells))
      maxnar = int(spfra*nd*ncells)
      allocate(iwgp(maxnar*2+1),colgp(maxnar),rwgp(maxnar))

      call system_clock(time_begin) 
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
          nzid = nzid+1
          rwgp(nzid) = subrow(ix)
          iwgp(1+nzid) = ii
          colgp(nzid) = ix
      endif
      enddo
      enddo
      leniwgp = nzid*2+1
      lenrwgp = nzid
      iwgp(1) = lenrwgp
      iwgp(nzid+2:nzid*2+1) = colgp(1:nzid)
      call system_clock(time_end) 
      if(iproj==1) print *, "time for projecting grid to cells:",(time_end-time_begin)/1000

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

      conlim = 100
      itnlim = 400
      atol = 1e-6!1e-3/((dvxd+dvzd)*111.19/2.0*0.1) !1e-2
      btol = 1e-6!1e-3/(dvxd*nx*111.19/3.0)!1e-3
      istop = 0
      anorm = 0.0
      acond = 0.0
      arnorm = 0.0
      xnorm = 0.0
      localSize = 10!int(ncells/4)
      !damp = dampvel
      ! using lsmr to solve for the projection coefficients
      !print*, 'LSMR beginning ...'

      nout = -1
      !nout = 36
      !open(nout,file='lsmrout_sub.txt')

      !print*, 'sth is wrong...'
      call system_clock(time_begin) 
      call LSMR(nd, ncells, leniwgp, lenrwgp,iwgp,rwgp,dres,damp,&
      atol, btol, conlim, itnlim, localSize,nout,&
      xunknown, istop, itn, anorm, acond,rnorm, arnorm, xnorm)
      call system_clock(time_end) 
      if(iproj==1) print *, "time for lsqr:",(time_end-time_begin)/1000
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

      call system_clock(time_begin) 
      dv(1:ndim) = 0
      call aprod(2,ncells,ndim,dv,xunknown,leniw_p,lenrw_p,iw_p,rw_p)
      call system_clock(time_end) 
      if(iproj==1) print *, "time for solution:",(time_end-time_begin)/1000
      deallocate(grow,gcol,subrow)
      deallocate(theta,phi,rrad,dws,norm)
      deallocate(xpts,ypts,zpts,dis,xunknown)
      deallocate(iw_p,rw_p,row,col)
      deallocate(lat,lon,rad)
      deallocate(iwgp,colgp,rwgp)
      deallocate(rad1,rad2,radnz)

contains
subroutine random_index( idx, weights )
    integer :: idx
    real, intent(in) :: weights(:)

    real x, wsum, prob

    wsum = sum( weights )

    call random_number( x )

    prob = 0
    do idx = 1, size( weights )
        prob = prob + weights( idx ) / wsum   !! 0 < prob < 1
        if ( x <= prob ) exit
    enddo
end subroutine random_index

subroutine unique_1d(arr1d_unique, arr1d )
    integer :: ii,jj
    real, intent(in) :: arr1d(:)
    real,dimension(:),allocatable:: arr1d_unique
    integer,dimension(:),allocatable:: mask1d
    allocate(mask1d(size( arr1d )))
    mask1d = 1
    do ii = 1, size( arr1d )
      do jj = ii+1, size(arr1d)
        if(arr1d(ii)==arr1d(jj)) then
          mask1d(ii)=0
        endif
      enddo
    enddo
    arr1d_unique = pack(arr1d, mask1d>0)
    deallocate(mask1d)

end subroutine unique_1d

end subroutine
