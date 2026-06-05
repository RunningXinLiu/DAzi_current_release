subroutine voronoiproj(leniw,lenrw,colg,nrow,rw,dres,goxd,dvxd,gozd,dvzd,depz,&
                       nx,ny,nz,nd,ncell,hvratio,damp,iproj,iter,dv)
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
      real hvratio,cmb
      integer ndim,nd
      integer iproj,iter
      integer,parameter:: idist=1

      real,parameter:: radius = 6371.0,ftol = 0.0001,pi = 3.141592654
      integer ii,ix,iy,iz,jj
      real,dimension(:),allocatable:: grow,subrow,dis,dws,xunknown!,gcol
      real,dimension(:),allocatable:: lat,lon,rad,theta,thetalayer,phi,philayer,rrad,xpts,ypts,zpts
      real,dimension(:),allocatable :: rw_p,rwgp,norm
      integer,dimension(:),allocatable:: iw_p,row,col,iwgp,colgp
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
      integer start

      allocate(lat(nx-2),lon(ny-2),rad(nz-1))
      ndim = (nx-2)*(ny-2)*(nz-1)

      nzrand = 10!int((nz-1)/2)
      acells = int(ncell/3)
      if(acells>300) then
            ncells = ncell*nzrand+acells*nzrand
      else
          ncells = ncell*nzrand
      endif

      allocate(theta(ncells),phi(ncells),thetalayer(ncell),philayer(ncell),rrad(ncells),norm(ncells))
      allocate(xpts(ncells),ypts(ncells),zpts(ncells),dis(ncells),xunknown(ncells))
      allocate(rw_p(ndim))
      allocate(iw_p(2*ndim+1),row(ndim),col(ndim),dws(ndim))

      iseed(1:3) = (/38,62,346/)
      iseed(4) = 2*(iproj+iter)+1
      call slarnv(idist,iseed,ncell,thetalayer)
      !print*,iproj,idist,minval(thetalayer),maxval(thetalayer)
      thetalayer = (gozd+thetalayer*(ny-3)*dvzd)*pi/180
      !print*,iproj,minval(thetalayer),maxval(thetalayer)
      do ii = 1,nzrand
      theta((ii-1)*ncell+1:ii*ncell) = thetalayer(1:ncell)
      enddo
      call slarnv(idist,iseed,ncell,philayer)
      !print*,iproj,minval(philayer),maxval(philayer)
      philayer = pi/2-(goxd-philayer*(nx-3)*dvxd)*pi/180
      !print*,iproj,minval(philayer),maxval(philayer)
      do ii = 1,nzrand
      phi((ii-1)*ncell+1:ii*ncell) = philayer(1:ncell)
      enddo
      call slarnv(idist,iseed,nzrand,rad)
      rad = rad*depz(nz-1)
      do ii = 1,nzrand
      do jj = 1,ncell
      !rrad((ii-1)*ncell+jj) = radius - (depz((ii-1)*2+1)+depz(ii*2))/2.0*hvratio
      rrad((ii-1)*ncell+jj) = radius - rad(ii)*hvratio!(depz((ii-1)*2+1)+depz(ii*2))/2.0*hvratio
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
      if (acells > 300) then
      dws = 0
      do ii = 1,lenrw
      dws(colg(ii)) = dws(colg(ii))+abs(rw(ii))
      enddo
      !acells = int(ncell/2.0)
      do ii = 1, acells
      call random_index(idx,dws) 
      ix = mod(idx,nx-2)
      iy = idx/(nx-2)
      !iz = idx/((nx-2)*(ny-2))
      do jj = 1,nzrand
      theta(nzrand*ncell+(ii-1)*nzrand+jj) = (gozd+(ix+1)*dvzd)*pi/180
      phi(nzrand*ncell+(ii-1)*nzrand+jj) = pi/2-(goxd-(iy+1)*dvxd)*pi/180
      !rrad(nzrand*ncell+(ii-1)*nzrand+jj) = radius-depz(iz+1)*hvratio
      rrad(nzrand*ncell+(ii-1)*nzrand+jj) = radius - (depz((jj-1)*2+1)+depz(jj*2))/2.0*hvratio
      enddo
      enddo
      endif

      xpts = rrad*sin(phi)*cos(theta)
      ypts = rrad*sin(phi)*sin(theta)
      zpts = rrad*cos(phi)


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
      enddo
      enddo
      enddo
      rw_p = 1.0
      leniw_p = 2*ndim+1
      lenrw_p = ndim
      iw_p(1) = ndim
      iw_p(2:ndim+1) = row
      iw_p(ndim+2:2*ndim+1) = col

      allocate(grow(ndim),subrow(ncells))
      maxnar = int(0.6*nd*ncells)
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
      
      call aprod(1,ncells,ndim,grow,subrow,leniw_p,lenrw_p,iw_p,rw_p)
    !   write(*,*) ' Doing vorotomo !'
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
      itnlim = 100
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
      call LSMR(nd, ncells, leniwgp, lenrwgp,iwgp,rwgp,dres,damp,&
      atol, btol, conlim, itnlim, localSize,nout,&
      xunknown, istop, itn, anorm, acond,rnorm, arnorm, xnorm)
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

      dv(1:ndim) = 0
      call aprod(2,ncells,ndim,dv,xunknown,leniw_p,lenrw_p,iw_p,rw_p)
      deallocate(grow,subrow)
      deallocate(theta,phi,thetalayer,philayer,rrad,dws,norm)
      deallocate(xpts,ypts,zpts,dis,xunknown)
      deallocate(iw_p,rw_p,row,col)
      deallocate(lat,lon,rad)
      deallocate(iwgp,colgp,rwgp)

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


end subroutine
