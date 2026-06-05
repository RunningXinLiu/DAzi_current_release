module fmm_backend
    implicit none

contains

    subroutine CalSurfGParallelIso(nx, ny, nz, nparpi, vels, iw, rw, col, dsurf, &
        goxdf, gozdf, dvxdf, dvzdf, kmaxRc, tRc, periods, depz, minthk, &
        scxf, sczf, rcxf, rczf, nrc1, nsrcsurf1, wavetype, igrt, &
        kmax, nsrcsurf, nrcf, nar, maxnar, writepath, dall, tRcV, nthreads)
        implicit none
        integer, intent(in) :: nx, ny, nz, nparpi, kmaxRc, kmax, nsrcsurf, nrcf
        integer, intent(in) :: writepath, dall, nthreads
        integer(kind=8), intent(inout) :: nar
        integer(kind=8), intent(in) :: maxnar
        real, intent(in) :: vels(nx,ny,nz)
        real, intent(inout) :: rw(*), dsurf(*)
        integer, intent(inout) :: iw(*), col(*)
        real, intent(in) :: goxdf, gozdf, dvxdf, dvzdf, minthk
        real*8, intent(in) :: tRc(*)
        integer, intent(in) :: periods(nsrcsurf,kmax)
        integer, intent(in) :: nrc1(nsrcsurf,kmax), nsrcsurf1(kmax)
        integer, intent(in) :: wavetype(nsrcsurf,kmax), igrt(nsrcsurf,kmax)
        real, intent(in) :: depz(nz), scxf(nsrcsurf,kmax), sczf(nsrcsurf,kmax)
        real, intent(in) :: rcxf(nrcf,nsrcsurf,kmax), rczf(nrcf,nsrcsurf,kmax)
        real*8, intent(inout) :: tRcV((nx-2)*(ny-2),kmaxRc)

        integer, allocatable :: countsurf(:,:,:)
        real*8 :: dummy_periods(1)
        integer :: checkstat, row_id, iper, isrc, irec

        interface
            subroutine CalSurfG_parallel(nx,ny,nz,nparpi,vels,iw,rw,col,dsurf, &
                goxdf,gozdf,dvxdf,dvzdf,kmaxRc,kmaxRg,kmaxLc,kmaxLg, &
                tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk, &
                scxf,sczf,rcxf,rczf,nrc1,nsrcsurf1,kmax,nsrcsurf,nrcf, &
                nar,maxnar,writepath,dall,countsurf,nthreads,use_aniso,Lsen_Gsc,rmax,tRcV)
                implicit none
                integer :: nx,ny,nz,nparpi,kmaxRc,kmaxRg,kmaxLc,kmaxLg
                integer :: kmax,nsrcsurf,nrcf,writepath,dall,nthreads
                logical, optional :: use_aniso
                integer, optional :: rmax
                integer(kind=8) :: nar,maxnar
                real :: vels(nx,ny,nz),rw(*),dsurf(*)
                integer :: iw(*)
                integer :: col(*),countsurf(nrcf,nsrcsurf,kmax)
                real :: goxdf,gozdf,dvxdf,dvzdf,minthk
                real*8 :: tRc(*),tRg(*),tLc(*),tLg(*)
                real, optional :: Lsen_Gsc(nx*ny,kmaxRc,nz-1)
                real*8, optional :: tRcV((nx-2)*(ny-2),kmaxRc)
                integer :: wavetype(nsrcsurf,kmax),igrt(nsrcsurf,kmax)
                integer :: periods(nsrcsurf,kmax),nrc1(nsrcsurf,kmax),nsrcsurf1(kmax)
                real :: depz(nz),scxf(nsrcsurf,kmax),sczf(nsrcsurf,kmax)
                real :: rcxf(nrcf,nsrcsurf,kmax),rczf(nrcf,nsrcsurf,kmax)
            end subroutine CalSurfG_parallel
        end interface

        if (maxnar < 1_8) stop 'CalSurfGParallelIso: maxnar must be positive'
        allocate(countsurf(nrcf, nsrcsurf, kmax), stat=checkstat)
        if (checkstat /= 0) stop 'CalSurfGParallelIso: allocate countsurf failed'

        countsurf = 0
        dummy_periods = 0.0_8

        row_id = 0
        do iper = 1, kmax
            do isrc = 1, nsrcsurf1(iper)
                do irec = 1, nrc1(isrc,iper)
                    row_id = row_id + 1
                    countsurf(irec,isrc,iper) = row_id
                enddo
            enddo
        enddo
        if (row_id /= dall) stop 'CalSurfGParallelIso: countsurf row mapping does not match dall'

        call CalSurfG_parallel(nx, ny, nz, nparpi, vels, iw, rw, col, dsurf, &
            goxdf, gozdf, dvxdf, dvzdf, kmaxRc, 0, 0, 0, &
            tRc, dummy_periods, dummy_periods, dummy_periods, &
            wavetype, igrt, periods, depz, minthk, &
            scxf, sczf, rcxf, rczf, nrc1, nsrcsurf1, kmax, nsrcsurf, nrcf, &
            nar, maxnar, writepath, dall, countsurf, nthreads, use_aniso=.false., tRcV=tRcV)

        deallocate(countsurf)
    end subroutine CalSurfGParallelIso

    subroutine CalSurfGParallelAniso(nx, ny, nz, nparpi, vels, iw, rw, col, dsurf, &
        Lsen_Gsc, rmax, tRcV, goxdf, gozdf, dvxdf, dvzdf, kmaxRc, tRc, periods, depz, minthk, &
        scxf, sczf, rcxf, rczf, nrc1, nsrcsurf1, kmax, nsrcsurf, nrcf, &
        nar, maxnar, writepath, dall, nthreads)
        implicit none
        integer, intent(in) :: nx, ny, nz, nparpi, rmax, kmaxRc, kmax, nsrcsurf, nrcf
        integer, intent(in) :: writepath, dall, nthreads
        integer(kind=8), intent(inout) :: nar
        integer(kind=8), intent(in) :: maxnar
        real, intent(in) :: vels(nx,ny,nz)
        real, intent(inout) :: rw(*), dsurf(*), Lsen_Gsc(nx*ny,kmaxRc,nz-1)
        real*8, intent(inout) :: tRcV((nx-2)*(ny-2),kmaxRc)
        integer, intent(inout) :: iw(*), col(*)
        real, intent(in) :: goxdf, gozdf, dvxdf, dvzdf, minthk
        real*8, intent(in) :: tRc(*)
        integer, intent(in) :: periods(nsrcsurf,kmax)
        integer, intent(in) :: nrc1(nsrcsurf,kmax), nsrcsurf1(kmax)
        real, intent(in) :: depz(nz), scxf(nsrcsurf,kmax), sczf(nsrcsurf,kmax)
        real, intent(in) :: rcxf(nrcf,nsrcsurf,kmax), rczf(nrcf,nsrcsurf,kmax)

        integer, allocatable :: countsurf(:,:,:), wavetype(:,:), igrt(:,:)
        real*8 :: dummy_periods(1)
        integer :: checkstat, row_id, iper, isrc, irec

        interface
            subroutine CalSurfG_parallel(nx,ny,nz,nparpi,vels,iw,rw,col,dsurf, &
                goxdf,gozdf,dvxdf,dvzdf,kmaxRc,kmaxRg,kmaxLc,kmaxLg, &
                tRc,tRg,tLc,tLg,wavetype,igrt,periods,depz,minthk, &
                scxf,sczf,rcxf,rczf,nrc1,nsrcsurf1,kmax,nsrcsurf,nrcf, &
                nar,maxnar,writepath,dall,countsurf,nthreads,use_aniso,Lsen_Gsc,rmax,tRcV)
                implicit none
                integer :: nx,ny,nz,nparpi,kmaxRc,kmaxRg,kmaxLc,kmaxLg
                integer :: kmax,nsrcsurf,nrcf,writepath,dall,nthreads
                logical, optional :: use_aniso
                integer, optional :: rmax
                integer(kind=8) :: nar,maxnar
                real :: vels(nx,ny,nz),rw(*),dsurf(*)
                integer :: iw(*)
                integer :: col(*),countsurf(nrcf,nsrcsurf,kmax)
                real :: goxdf,gozdf,dvxdf,dvzdf,minthk
                real*8 :: tRc(*),tRg(*),tLc(*),tLg(*)
                real, optional :: Lsen_Gsc(nx*ny,kmaxRc,nz-1)
                real*8, optional :: tRcV((nx-2)*(ny-2),kmaxRc)
                integer :: wavetype(nsrcsurf,kmax),igrt(nsrcsurf,kmax)
                integer :: periods(nsrcsurf,kmax),nrc1(nsrcsurf,kmax),nsrcsurf1(kmax)
                real :: depz(nz),scxf(nsrcsurf,kmax),sczf(nsrcsurf,kmax)
                real :: rcxf(nrcf,nsrcsurf,kmax),rczf(nrcf,nsrcsurf,kmax)
            end subroutine CalSurfG_parallel
        end interface

        if (maxnar < 1_8) stop 'CalSurfGParallelAniso: maxnar must be positive'
        allocate(countsurf(nrcf, nsrcsurf, kmax), stat=checkstat)
        if (checkstat /= 0) stop 'CalSurfGParallelAniso: allocate countsurf failed'
        allocate(wavetype(nsrcsurf,kmax), igrt(nsrcsurf,kmax), stat=checkstat)
        if (checkstat /= 0) stop 'CalSurfGParallelAniso: allocate ray type arrays failed'

        countsurf = 0
        wavetype = 2
        igrt = 0
        dummy_periods = 0.0_8

        row_id = 0
        do iper = 1, kmax
            do isrc = 1, nsrcsurf1(iper)
                do irec = 1, nrc1(isrc,iper)
                    row_id = row_id + 1
                    countsurf(irec,isrc,iper) = row_id
                enddo
            enddo
        enddo
        if (row_id /= dall) stop 'CalSurfGParallelAniso: countsurf row mapping does not match dall'

        call CalSurfG_parallel(nx, ny, nz, nparpi, vels, iw, rw, col, dsurf, &
            goxdf, gozdf, dvxdf, dvzdf, kmaxRc, 0, 0, 0, &
            tRc, dummy_periods, dummy_periods, dummy_periods, &
            wavetype, igrt, periods, depz, minthk, &
            scxf, sczf, rcxf, rczf, nrc1, nsrcsurf1, kmax, nsrcsurf, nrcf, &
            nar, maxnar, writepath, dall, countsurf, nthreads, .true., Lsen_Gsc, rmax, tRcV)

        deallocate(countsurf, wavetype, igrt)
    end subroutine CalSurfGParallelAniso

end module fmm_backend
