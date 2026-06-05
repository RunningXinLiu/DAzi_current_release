subroutine voronoiproj_3d_azi(leniw,lenrw,colg,nrow,rw,dres,goxd,dvxd,gozd,dvzd,depz,&
  nx,ny,nz,nd,ncells,acells,spfra,hvratio,damp,iproj,iter,dv)

use lsmrModule, only:lsmr

implicit none
integer(kind=8) leniw,lenrw
integer(kind=8) lenrw_iso
integer nx,ny,nz
!      integer iw(leniw)
integer colg(lenrw)
integer nrow(nd)

real depz(nz)
real rw(lenrw)
integer ncells,acells,athresh
integer allcells ! total voronoi cells
real dv(*),dres(*)
real goxd,gozd,dvxd,dvzd
real damp
real hvratio
real(kind=8) spfra
integer ndim,nd,ndim0
integer iproj,iter


real,parameter:: radius = 6371.0,ftol = 0.1,pi = 3.141592654
integer ii,jj,ix,iy,iz
real,dimension(:),allocatable:: grow,gcol,subrow,dis,xunknown
real,dimension(:),allocatable:: lat,lon,rad,theta,phi,rrad,xpts,ypts,zpts
real,dimension(:),allocatable :: rw_p,rwgp,norm
real,dimension(:),allocatable:: radnz
real,dimension(:),allocatable:: rad1,rad2
integer,dimension(:),allocatable:: row,col,colgp
integer,dimension(:),allocatable:: iw_p,iwgp,iw_p_3
integer idx

real,dimension(:),allocatable :: rw_p_3
real,dimension(:),allocatable:: rw_iso, dws_iso
integer(kind=8),dimension(:),allocatable:: idx_iso
integer,dimension(:),allocatable:: colg_iso


integer(kind=8) maxnar
integer iseed(4)
real xs,ys,zs

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
integer(kind=8) leniw_p,lenrw_p
integer(kind=8) leniwgp,lenrwgp
integer(kind=8) leniw_p_3,lenrw_p_3
integer(kind=8) start,nzid,i
! integer(kind=16) start,nzid
! integer*8 start,nzid
! integer, parameter :: LargeInt_K = selected_int_kind (18)
! integer (kind=LargeInt_K) :: start,nzid

allocate(lat(nx-2),lon(ny-2),rad(nz-1))
ndim0 = (nx-2)*(ny-2)*(nz-1) !this is the dimession of one variable
ndim = (nx-2)*(ny-2)*(nz-1)*3

! the regular grid
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

athresh = 100
! if the acell is too few, dont use acell
if(acells>athresh) then
  allcells = ncells+acells
else
  allcells = ncells
  write(*,*) 'The number of adaptive cell is too few! Not adding adaptive cell.'
endif

allocate(theta(allcells),phi(allcells),rrad(allcells),norm(3*allcells))
allocate(xpts(allcells),ypts(allcells),zpts(allcells),dis(allcells),xunknown(3*allcells))
allocate(rw_p(ndim0)) ! the projection matrix
allocate(iw_p(2*ndim0),row(ndim0),col(ndim0),dws_iso(ndim0))
allocate(rw_p_3(ndim)) ! the projection matrix
allocate(iw_p_3(2*ndim))
! allocate(colg_iso(ndim0))

!separate the big G matrix to Giso and Gc Gs-------------------!
call range_idx_1d(idx_iso,colg,0,ndim0+1)
lenrw_iso = size(idx_iso)
allocate(colg_iso(lenrw_iso))
colg_iso = colg(idx_iso)
rw_iso = rw(idx_iso)

iseed(1:3) = (/38,62,346/)
iseed(4) = 2*(iproj+iter)+1

call slarnv(1,iseed,ncells,theta)
theta(1:ncells) = (gozd+theta(1:ncells)*(ny-3)*dvzd)*pi/180
call slarnv(1,iseed,ncells,phi)
phi(1:ncells) = pi/2-(goxd-phi(1:ncells)*(nx-3)*dvxd)*pi/180
call slarnv(1,iseed,ncells,rrad)
rrad(1:ncells) = radius-rrad(1:ncells)*depz(nz-1)*hvratio

if(iproj==1) then
  open(77,file='pvpoints.dat')
  do ii=1,ncells
    write(77,'(3f10.4)') rrad(ii),phi(ii),theta(ii)
  enddo
  close(77)
  !print*,goxd,gozd,dvxd,dvzd,nx,ny
  !stop
endif


! write(*,*) 'Here 4!'
! adaptive cells based on dws, assume 1/2 of all ncells are used
! add adaptive cells
if(acells>athresh) then
  write(*,*)'Adjusting acells!'
  dws_iso = 0
  do i = 1,lenrw_iso
    dws_iso(colg_iso(i)) = dws_iso(colg_iso(i))+abs(rw_iso(i))
  enddo
  ! acells = int(ncells/2.0)
  do ii = 1,acells
    call random_index(idx,dws_iso) 

    ! ix = mod(idx,nx-2)
    ! iy = idx/(nx-2)
    ! iz = idx/((nx-2)*(ny-2))

    iz = (idx-1)/((nx-2)*(ny-2))+1  ! Ying
    iy = (idx-1+(1-iz)*(nx-2)*(ny-2))/(nx-2)+1  ! Ying
    ix = idx+(nx-2)*(1-iy)+(nx-2)*(ny-2)*(1-iz)  ! Ying

    theta(ncells+ii) = (gozd+(iy-1)*dvzd)*pi/180
    phi(ncells+ii) = pi/2-(goxd-(ix-1)*dvxd)*pi/180
    rrad(ncells+ii) = radius-depz(iz)*hvratio ! because the depz(0) is not the boundary
  enddo
endif

if(iproj==1) then
open(77,file='pvpoints_acell.dat')
do ii=1,allcells
write(77,'(3f10.4)') rrad(ii),phi(ii),theta(ii)
enddo
close(77)
!print*,goxd,gozd,dvxd,dvzd,nx,ny
!stop
endif

xpts = rrad*sin(phi)*cos(theta)
ypts = rrad*sin(phi)*sin(theta)
zpts = rrad*cos(phi)

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

! write(*,*) 'Here 5!'
rw_p = 1.0
leniw_p = 2*ndim0
lenrw_p = ndim0
iw_p(1:ndim0) = row
iw_p(ndim0+1:2*ndim0) = col
!---------------------------------------------------------------!
!compose a bigger projection matrix 
! P is ndim0 * ncells matrix
! Big P is ndim (ndim0*3) * (3*ncells)
!|P 0 0|
!|0 P 0|
!|0 0 P|
rw_p_3 = 1.0
leniw_p_3 = 2*ndim
lenrw_p_3 = ndim
! iw_p_3(2:ndim0+1) = row
! iw_p_3(ndim0+2:2*ndim0+1) = row + ndim0
! iw_p_3(2*ndim0+2:3*ndim0+1) = row + 2*ndim0
! iw_p_3(3*ndim0+2:4*ndim0+1) = col
! iw_p_3(4*ndim0+2:5*ndim0+1) = col + ncells
! iw_p_3(5*ndim0+2:6*ndim0+1) = col + 2*ncells

iw_p_3(1:ndim0) = row
iw_p_3(ndim0+1:2*ndim0) = row + allcells
iw_p_3(2*ndim0+1:3*ndim0) = row + 2*allcells
iw_p_3(3*ndim0+1:4*ndim0) = col
iw_p_3(4*ndim0+1:5*ndim0) = col + ndim0
iw_p_3(5*ndim0+1:6*ndim0) = col + 2*ndim0
!---------------------------------------------------------------!

allocate(grow(ndim),gcol(nd),subrow(3*allcells))
! Just use the same spafrac as in the main function
! if(INT8(nd)*3*ncells.lt.100000000) then
! maxnar = 0.6*INT8(nd)*3*ncells
! else
! maxnar = 0.6*INT8(nd)*3*ncells
! endif
maxnar = INT8(spfra*INT8(nd)*3*allcells)

if((iproj==1)) then
write(*,*) 'Max Number of sparse matrix (GP) non-zero element ',maxnar
endif
allocate(iwgp(maxnar*2),colgp(maxnar),rwgp(maxnar))

nzid = 0
do ii = 1,nd
grow = 0
start = 0
! use a do loop to realize sum, to avoid initialize nrow as int8
do jj = 1,ii-1
  start = start + nrow(jj)
enddo
! start = sum(nrow(1:ii-1)) ! make nrow int8 because sum return the same data type as the array

! if((iproj==1).and.(ii==nd)) then
! if((iproj==1)) then
!   write(*,*) 'start is ',start,'Number of non-zero of G is ',lenrw,', left ',lenrw-start,'at iter',iter,', proj ',iproj
! endif
do ix = 1,nrow(ii)
grow(colg(start+ix)) = rw(start+ix)
if(start+ix.gt.lenrw) stop 'clog index exceed range!'
enddo
if((iproj==1)) then
! if((iproj==1).and.(ii == nd)) then
! write(*,'(a,I14,a,I14,a,I14,a,I4,a,I4)') 'start is ',start+ix,'. Number of non-zero of G is ',lenrw,', left '&
! ,lenrw-start-ix,', at iter',iter,', proj ',iproj
endif
!gcol = 0
!gcol(ii) = 1.0 
!call aprod(2,nd,ndim,grow,gcol,leniw,lenrw,iw,rw)
subrow = 0
call aprod(1,3*allcells,ndim,grow,subrow,leniw_p_3,lenrw_p_3,iw_p_3,rw_p_3)
do ix = 1,3*allcells
if(abs(subrow(ix))>ftol) then
nzid = nzid+1
rwgp(nzid) = subrow(ix)
iwgp(nzid) = ii
colgp(nzid) = ix
endif
enddo

enddo

leniwgp = nzid*2
lenrwgp = nzid
iwgp(nzid+1:nzid*2) = colgp(1:nzid)
if(iproj==1) then
write(*,*) 'Number of element for index (lenrwgp) of GP is ',lenrwgp
write(*,*) 'Max number of sparse matrix (GP) non-zero element ',maxnar
if (lenrwgp > maxnar) stop 'increase sparsity fraction(spfra) of voronoi projection!'
endif
! norm = 0
! do ii=1,nzid
! norm(iwgp(ii+nzid)) = norm(iwgp(ii+nzid))+rwgp(ii)**2
! enddo

! do ii =1,3*ncells
! norm(ii) = sqrt(norm(ii)/nd+0.01)
! enddo

! do ii =1,nzid
! rwgp(ii) = rwgp(ii)/norm(iwgp(ii+nzid))
! enddo
! write(*,*) 'Here 7!'
conlim = 50
itnlim = 100
atol = 1e-3/((dvxd+dvzd)*111.19/2.0*0.1) !1e-2
btol = 1e-3/(dvxd*nx*111.19/3.0)!1e-3
istop = 0
anorm = 0.0
acond = 0.0
arnorm = 0.0
xnorm = 0.0
localSize = int(allcells/4)
!damp = dampvel
! using lsmr to solve for the projection coefficients
!print*, 'LSMR beginning ...'

nout = -1
!nout = 36
!open(nout,file='lsmrout_sub.txt')

call LSMR(nd, 3*allcells, leniwgp, lenrwgp,iwgp,rwgp,dres,damp,&
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
! write(*,*) 'Here 8!'
dv(1:ndim) = 0
call aprod(2,3*allcells,ndim,dv,xunknown,leniw_p_3,lenrw_p_3,iw_p_3,rw_p_3)

deallocate(grow,gcol,subrow)
deallocate(theta,phi,rrad,norm)
deallocate(xpts,ypts,zpts,dis,xunknown)
deallocate(iw_p,rw_p,row,col)
deallocate(lat,lon,rad)
deallocate(iwgp,colgp,rwgp)
deallocate(rw_p_3)
deallocate(iw_p_3)
deallocate(colg_iso,idx_iso,rw_iso)
deallocate(dws_iso)


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

subroutine range_idx_1d( idx, arr, value1, value2 )
!no multi-modal, only for 1d array index
implicit none
integer(kind=8),dimension(:),allocatable :: idx(:)
integer, intent(in) :: arr(:)
integer(kind=8) :: i,count1
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
idx=idx(1:count1)
end subroutine range_idx_1d

end subroutine
