! Direct Inversion for 3-D Azimuthal Anisotropy
! V1.1 (2017) based on Mineos to calculate kernel
! V2.0 (2019) based on CPS-tregn96 to calcualte kernel
! Copyright:
!    Author: Chuanming Liu (USTC)
!     Email: Chuanming.liu@colorado.edu
! Reference:
! Liu, C., Yao, H., Yang, H., Shen, W., Fang, H., Hu, S., Qiao, L., 2019. Direct
! inversion for three-dimensional shear wavespeed azimuthal anisotropy based on
! surface-wave ray tracing: methodology and application to Yunnan, southwest China.
! Journal of Geophysics Research: Solid Earth. doi: 10.1029/2018JB016920

program SurfAniso
        use lsmrModule, only:lsmr
        use lsmrblasInterface, only : dnrm2
        use omp_lib
        use, intrinsic :: ieee_arithmetic
        use fmm_backend, only: CalSurfGParallelIso, CalSurfGParallelAniso
        implicit none
        character inputfile*100
        character logfile*100
        ! character outmodel*100
        character outmodel_iso*100
        character outmodel_azi*100
        character outsyn*100
        logical ex
        character dummy*40
        character datafile*80
        integer nx,ny,nz
        real goxd,gozd
        real dvxd,dvzd
        integer nsrc,nrc
        real weightVs,weight0,weight1,weightGcs,weight
        real damp,dampAA,dampVs
        real minthk
        integer kmax,kmaxRc
        real*8,dimension(:),allocatable:: tRc
        real,dimension(:),allocatable:: depz
        integer itn
        integer nout
        integer localSize
        real mean,std_devs,balances,balanceb
        integer msurf
        real,dimension(:),allocatable:: obst,dsyn,cbst,dist
        real,dimension(:),allocatable:: pvall
        real sta1_lat,sta1_lon,sta2_lat,sta2_lon
        real dist1
        integer dall
        integer istep
        real,parameter :: pi=3.1415926535898
        integer checkstat
        integer ii,jj,kk
        real, dimension (:,:), allocatable :: scxf,sczf
        real, dimension (:,:,:), allocatable :: rcxf,rczf
        integer,dimension(:,:),allocatable::wavetype,igrt,nrc1
        integer,dimension(:),allocatable::nsrc1
        integer,dimension(:,:),allocatable::periods
        real,dimension(:),allocatable::rw
        integer, dimension(:),allocatable::iw,col,nrow
        real,dimension(:),allocatable::dv,norm,dvsub,dvstd,dvstdall,dvall
        real,dimension(:),allocatable :: dv_mean_accum
        real,dimension(:),allocatable :: dv_mean2_accum
        real,dimension(:),allocatable :: dv_var_within_accum
        integer :: nuq_iter
        real,dimension(:,:,:),allocatable::vsf
        character strf
        integer veltp,wavetp
        real velvalue
        integer knum,knumo,err
        integer istep1,istep2
        integer period
        integer knumi,srcnum,count1
        integer HorizonType,VerticalType
        character line*200
        integer iter,maxiter
        ! integer*4/8
        integer(kind=8) maxnar
        real acond,anorm,arnorm,rnorm,xnorm
        character str1
        real atol,btol
        real conlim
        integer istop
        integer itnlim
        integer(kind=8) lenrw,leniw
        integer(kind=8) nar
        integer(kind=8) nar_tmp,nars
        integer count3,nvz,nvx
        integer m,maxvp,n
        integer(kind=8) i,j,k
        real(kind=8) spfra
        real noiselevel
        integer ifsyn
        integer writepath
        real averdws
        real maxnorm
        real maxborder
        logical use_sparse_g

        character SenGTruefile*200, SenGfile*200
        integer lay
        integer maxm
        real, dimension(:), allocatable:: sigmaT,resSigma,datweight
        real, dimension(:), allocatable:: Tdata,resbst,fwdT,fwdTvs,fwdTaa,Tref,RefTaa,resbst_iso

        real meandeltaT,stddeltaT
        real,dimension(:,:),allocatable:: GVs,GGs,GGc
        real,dimension(:,:,:), allocatable :: Lsen_Gsc
        real gaussian
        real synT
        !real,dimension(:,:,:),allocatable::vsRela,vsAbs,vsreal
        real,dimension(:,:,:),allocatable:: gcf,gsf
        integer(kind=8) nar1,nar2
        real vsref
        real res2Norm,Mnorm2,Enorm2,Vsnorm2,Gcnorm2,Gsnorm2
        real resNSigma2Norm
        real Enorm2Lame,Vsnorm2Lame,Gcnorm2Lame,Gsnorm2Lame
        real MwNorm2
        ! ...existing code...

        !For Poisson Voronoi inverison
        integer iproj,vorotomo,ncell,acell,nrealizations,idx
        integer nzrand
        real hvratio,iaratio
        real,dimension(:),allocatable:: dvsub_local
        integer(kind=8) :: nnz_check

        ! for voronoi inversion misfit metric
        integer :: n_nan, n_clip
        real    :: mf_min, mf_max
        integer :: maxvp1,maxvp2
        integer,dimension(:),allocatable:: itn_voro, istop_voro
        real    :: rnorm_out, acond_out
        integer :: itn_out, istop_out
        real,dimension(:),allocatable:: misfit_voro, acond_voro
        ! for sorted misfit
        real, allocatable :: misfit_sorted(:)
        integer, allocatable :: idx_sorted(:)
        ! for weighted misfit
        real :: eps, sumw
        integer :: n_keep, n_good
        integer :: i_q1, i_q3, k0, kmin
        real    :: q1, q3, iqr, alpha
        real    :: csum
        integer :: nkeep2
        integer :: fid_misfit
        integer :: fid_misfit_voro
        real, allocatable :: w_rank(:)       ! weights in rank order (sorted)
        real, allocatable :: w_by_real(:)    ! weights mapped back to realization id (optional)
        real, parameter :: MISFIT_BAD = 1.0e20
        real*8,dimension(:,:),allocatable:: tRcV
        ! ---- hard-gate thresholds for voronoi realizations (minimal, gfortran-safe) ----
        real,    parameter :: COND_MAX  = 500.0
        integer, parameter :: ITN_MIN   = 2
        real sumObs
        integer Refwritepath,realwritepath
        real VariGc,VariGs,VariVs
        integer count4,para
        real LambdaVs,LambdaGc,LambdaGs
        real PreRes,meanAbs
        integer iter_mod,rmax
        character(len=2) :: id='00'
        character(len=60) :: filename
        real*8 :: startT,endT
        logical:: isTikh
        real :: summ, threshold0,thresholdVs,threshold
        real :: pertV
        logical :: iso_inv, iso_mod
        real Minvel,MaxVel
        real,dimension(:),allocatable:: GcInv, GsInv, VsInv
        real mindVs, maxdVs, meadVs
        real minGs, maxGs, meaGs
        real minGc, maxGc, meaGc
        real lameGcs, lameVs
        real mean2
        integer Nre, NreVs, narVs

        integer output_intermediate
        integer :: ncost
        real    :: data_misfit_norm, data_misfit_norm2
        real    :: model_l2_norm, model_l2_norm2
        real    :: mae_res, std_res, rms_res
        real*8  :: timer_iter_start, timer_iter_end
        real*8  :: timer_section_start, timer_section_end
        real*8  :: timer_g_matrix, timer_residual_weight
        real*8  :: timer_solver, timer_post_update, timer_raypath_output
        character(len=64) :: lsmr_backend, lsmr_tag
        integer :: lsmr_env_status
        character(len=64) :: ray_backend, fmm_threads_text
        integer :: ray_env_status, fmm_env_status, fmm_threads, read_status
        character(len=64) :: azi_vs_mode
        integer :: azi_vs_mode_status
        logical :: azi_freeze_vs
        integer :: residual_m, residual_n
        ! ---- defaults to avoid uninitialized voronoi controls when vorotomo==0 ----
        vorotomo      = 0
        ncell         = 0
        acell         = 0
        nrealizations = 1
        nzrand        = 0
        iaratio       = 1.0
        azi_vs_mode   = 'joint'
        azi_freeze_vs = .false.
        ! -------------------------------------------------------------------------

        ! open output for LSMR
        startT=OMP_get_wtime()
        nout=36
        open(nout,file='lsmr.txt')
        ! terminal output
        write(*,*)
        write(*,*) '                       DAzimSurfTomo'
        write(*,*)

        ! read contral file
        if (iargc() < 1) then
            write(*,*) 'input file [para.in (Default)]:'
            read(*,'(a)') inputfile
            if (len_trim(inputfile) <=1 ) then
                inputfile='para.in'
            else
                inputfile=inputfile(1:len_trim(inputfile))
            endif
        else
            call getarg(1, inputfile)
        endif
        inquire(file=inputfile, exist=ex)
        if (.not. ex)   stop 'unable to open the inputfile'

        open(10,file=inputfile,status='old',action='read')
        read(10,'(a30)') dummy
        read(10,'(a30)') dummy
        read(10,'(a30)') dummy
        ! read(10,*) SenGfile
        read(10,*) datafile
        read(10,*) nx,ny,nz
        read(10,*) goxd,gozd
        read(10,*) dvxd,dvzd
        read(10,*) minthk
        read(10,*) Minvel,Maxvel
        read(10,*) nsrc
        read(10,*) spfra
        read(10,*) maxiter
        read(10,*) iso_mod
        read(10,'(a30)') dummy
        read(10,*) weightVs
        read(10,*) weightGcs
        read(10,*) damp
        ! read(10,*) isTikh
        ! read(10,*) LcorrXY,LcorrZZ
        ! write(*,*) 'integrated  sensitivity kernel for Vs and Gcs file based on MOD  '
        ! write(*,'(a)')SenGfile
        write(*,*)'input Rayleigh wave phase velocity data file:'
        write(*,'(a)') datafile
        write(*,*)  'model origin:latitude,longitue'
        write(*,'(2f10.4)') goxd,gozd
        write(*,*) 'grid spacing:latitude,longitue'
        write(*,'(2f10.4)') dvxd,dvzd
        write(*,*) 'model dimension:nx,ny,nz'
        write(*,'(3i5)') nx,ny,nz
        write(*,*)'depth refined interval layer '
        write(*,'(f8.1)')minthk
        write(*,*)'smoothing for dVsv '
        write(*,'(f8.1)')weightVs
        write(*,*)'smoothing for Gc, Gs '
        write(*,'(f8.1)')weightGcs
        write(*,*)'damping'
        write(*,'(f8.1)')damp

        ! write(*,*)' Regularization Type: (T) 1st order Tikhonov ;(F) Gaussian'
        ! write(*,*) 'correlation length: XY, ZZ (km)'
        ! write(*,'(f8.1)')LcorrXY, LcorrZZ

        if (nz.LE.1)  stop 'error nz value.'
        read(10,'(a30)')dummy
        read(10,*) kmaxRc
        write(*,*) 'number of period'
        write(*,'(i6)') kmaxRc

        if(kmaxRc.gt.0)then
            allocate(tRc(kmaxRc),stat=checkstat)
            if (checkstat > 0) stop 'error allocating RP'
            read(10,*)(tRc(i),i=1,kmaxRc)
        else
            stop 'Can only deal with Rayleigh wave phase velocity data!'
        endif

        ! read in output control parameters
        read(10,'(a30)') dummy
        read(10,*) output_intermediate
        read(10,*) Refwritepath
        ! read in vorotomo parameters
        read(10,'(a30)') dummy
        read(10,*) vorotomo
        write(logfile,'(a,a)') trim(inputfile),'_inv.log'
        !open(66,file=logfile,action='write')
        open(66, file=logfile)
        if (vorotomo /= 0) then
            read(10,*) ncell ! normal voronoi cell number per layer, around (nx*ny)/4
            read(10,*) acell ! adaptive voronoi cell number per layer, around (nx*ny)/10
            read(10,*) nrealizations ! number of voronoi realizations, usually 100-200
            read(10,*) nzrand ! number of random layers for voronoi cell selection, usually nz/2

            if (.not. iso_mod) then
                read(10,*) iaratio ! iaratio = cell number ratio for isotropy vs anisotropy, usually 0.5-0.7
            endif
            write(*,*) 'Using vorotomo method!'
            write(*,'(a,i5,a,i5,a,i5)') 'Number of normal voronoi cells per layer: ',ncell, &
            ' | Number of adaptive voronoi cells per layer: ',acell, &
            ' | Number of voronoi realizations: ',nrealizations
            write(66,*) 'Using vorotomo method!'
            write(66,'(a,i5,a,i5,a,i5)') 'Number of normal voronoi cells per layer: ',ncell, &
            ' | Number of adaptive voronoi cells per layer: ',acell, &
            ' | Number of voronoi realizations: ',nrealizations
            if (.not. iso_mod) then
                write(*,'(a,f5.2)') 'Isotropy/ anisotropy cell number ratio: ',iaratio
                write(66,'(a,f5.2)') 'Isotropy/ anisotropy cell number ratio: ',iaratio
            endif
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! A) Dimensionality reduction ratio: maxvp / allcells
        !     • It is recommended to be around 4–10 for iso, and 8–15 for joint (to reduce dimensions without being too coarse)

        ! B) Average number of observations per cell: dall / allcells
        !     • It is recommended to be > 15–20 for iso, and > 20–30 for joint
        ! C) Number of realizations:
        !     • It is recommended to be around 50–100 for iso, and 100–200 for joint (Increase if you want to do Uncertainty Quantification)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !======= output control parameters =======
        ! write(logfile,'(a,a)') trim(inputfile),'_inv.log'
        !open(66,file=logfile,action='write')
        ! open(66, file=logfile)
        write(66,*)
        write(66,*) '                  DAzimSurfTomo'
        write(66,*)
        write(66,*) 'model origin:latitude,longitue'
        write(66,'(2f10.4)') goxd,gozd
        write(66,*) 'grid spacing:latitude,longitue'
        write(66,'(2f10.4)') dvxd,dvzd
        write(66,*) 'model dimension:nx,ny,nz'
        write(66,'(3i5)') nx,ny,nz

        write(*,*)'Rayleigh wave phase velocity used,periods:(s)'
        write(*,'(50f6.1)')(tRc(i),i=1,kmaxRc)
        write(66,*)'Rayleigh wave phase velocity used,periods:(s)'
        write(66,'(50f6.1)')(tRc(i),i=1,kmaxRc)
        close(10)   ! reading input file finished

        nrc=nsrc
        kmax=kmaxRc
        lay=nz-1
        !-----------------------------------------------------------------------!
        ! read measurements
        inquire(file=datafile, exist=ex)
        if (.not. ex) then
            write(66,'(a)') 'unable to open the datafile'
            close(66)
            stop 'unable to open the datafile'
        endif
        write(*,*) 'begin load data file.....'

        open(unit=87,file=datafile,status='old')

        allocate(scxf(nsrc,kmax),sczf(nsrc,kmax),&
        rcxf(nrc,nsrc,kmax),rczf(nrc,nsrc,kmax),stat=checkstat)
        if(checkstat > 0) stop 'scxf/rcxf allocation failed'
        scxf = 0.0
        sczf = 0.0
        rcxf = 0.0
        rczf = 0.0

        allocate(periods(nsrc,kmax),wavetype(nsrc,kmax),nrc1(nsrc,kmax),nsrc1(kmax),&
        igrt(nsrc,kmax),stat=checkstat)
        if(checkstat > 0) stop 'period/source-count allocation failed'
        periods = 0
        wavetype = 0
        igrt = 0
        nrc1 = 0
        nsrc1 = 0

        allocate(obst(nrc*nsrc*kmax),dist(nrc*nsrc*kmax),stat=checkstat)
        if(checkstat > 0) stop 'observation array allocation failed'
        obst=0
        dist=0

        allocate(pvall(nrc*nsrc*kmax),stat=checkstat)
        if(checkstat > 0) stop 'phase velocity array allocation failed'

        istep=0
        istep2=0
        dall=0
        knumo=12345
        knum=0
        istep1=0
        pvall=0
        do
            read(87,'(a)',iostat=err) line
            if(err.eq.0) then
                if(line(1:1).eq.'#') then
                    read(line,*) str1, sta1_lat, sta1_lon, period, wavetp, veltp
                    if(wavetp.eq.2.and.veltp.eq.0) knum=period
                    if(wavetp.eq.2.and.veltp.eq.1) stop 'can not deal with Rayleigh wave group data'
                    if(wavetp.eq.1.and.veltp.eq.0) stop 'can not deal with Love wave phase data'
                    if(wavetp.eq.1.and.veltp.eq.1) stop 'can not deal with Love wave group data'
                    if(knum.lt.1 .or. knum.gt.kmax) stop 'period index out of range in data header'
                    if(knum.ne.knumo) then
                        istep=0
                        istep2=istep2+1
                    endif
                    istep=istep+1
                    if(istep.gt.nsrc) stop 'source count exceeds nsrc for this period'
                    istep1=0
                    sta1_lat=(90.0-sta1_lat)*pi/180.0
                    sta1_lon=sta1_lon*pi/180.0
                    scxf(istep,knum)=sta1_lat
                    sczf(istep,knum)=sta1_lon
                    periods(istep,knum)=period
                    wavetype(istep,knum)=wavetp
                    igrt(istep,knum)=veltp
                    nsrc1(knum)=istep
                    knumo=knum
                else
                    read(line,*) sta2_lat,sta2_lon,velvalue
                    if(knum.lt.1 .or. knum.gt.kmax .or. istep.lt.1) stop 'receiver row encountered before a valid source header'
                    istep1=istep1+1
                    if(istep1.gt.nrc) stop 'receiver count exceeds nrc for this source'
                    dall=dall+1
                    if(dall.gt.nrc*nsrc*kmax) stop 'data count exceeds allocated observation arrays'
                    sta2_lat=(90.0-sta2_lat)*pi/180.0
                    sta2_lon=sta2_lon*pi/180.0
                    rcxf(istep1,istep,knum)=sta2_lat
                    rczf(istep1,istep,knum)=sta2_lon
                    call delsph(sta1_lat,sta1_lon,sta2_lat,sta2_lon,dist1)
                    dist(dall)=dist1
                    obst(dall)=dist1/velvalue
                    pvall(dall)=velvalue
                    nrc1(istep,knum)=istep1
                endif
            else
                exit
            endif
        enddo
        close(87)
        write(*,'(a,i7)') ' Number of all measurements', dall
        !----------------------------------------------------------------------!
        !  Initialization
        allocate(depz(nz), stat=checkstat)
        allocate(vsf(nx,ny,nz), stat=checkstat)

        ray_backend = 'serial'
        call get_environment_variable('DAZI_RAY_BACKEND', ray_backend, status=ray_env_status)
        if (ray_env_status /= 0 .or. len_trim(ray_backend) == 0) ray_backend = 'serial'
        use_sparse_g = (trim(adjustl(ray_backend)) == 'fmm_parallel')

        maxvp=(nx-2)*(ny-2)*(nz-1)
        if (iso_mod) then
            maxm = maxvp
            ! maxnar: upper bound for number of nonzeros in sparse G (kind=8 / int64).
            ! Use REAL(kind=8) in the product to avoid intermediate overflow in default INTEGER arithmetic.
            maxnar = int(spfra * real(dall,kind=8) * real(nx,kind=8) * real(ny,kind=8) * real(nz,kind=8), kind=8)
        else
            maxm = (nx-2)*(ny-2)*(nz-1)*3
            ! maxnar: upper bound for number of nonzeros in sparse G (kind=8 / int64).
            ! Use REAL(kind=8) in the product to avoid intermediate overflow in default INTEGER arithmetic.
            maxnar = int(spfra * real(dall,kind=8) * real(nx,kind=8) * real(ny,kind=8) * real(nz,kind=8) * 3.0_8, kind=8)
        endif
        if (maxnar < 0_8) stop 'number overflow, decrease your sparsefrac'
        ! allocate(dv(maxm), stat=checkstat)
                allocate(dv(maxm),dvsub(maxm),dvstdall(maxm),dvstd(maxm),dvall(maxm*max(1,nrealizations)))
        allocate(dv_mean_accum(maxm), dv_mean2_accum(maxm), dv_var_within_accum(maxm))
        dvstdall              = 0.0
        dv_mean_accum         = 0.0
        dv_mean2_accum        = 0.0
        dv_var_within_accum   = 0.0
        nuq_iter              = 0
        ! allocate(dv(maxm),dvstdall(maxm),dvstd(maxm),dvall(maxm*nrealizations), stat=checkstat)
!        allocate(dvall(maxvp*nrealizations),stats=checkstat)
        ! allocate(dv(maxm), stat=checkstat)
        allocate(GcInv(maxvp), GsInv(maxvp), VsInv(maxvp), stat=checkstat)
        allocate(rw(maxnar), stat=checkstat)
        allocate(iw(2*maxnar), stat=checkstat)
        ! allocate(col(maxnar), stat=checkstat)
        allocate(col(maxnar),nrow(dall), stat=checkstat)
        if(checkstat > 0)then
                write(6,*)'error with allocate:  integer iw'
        endif
        allocate(norm(maxvp), stat=checkstat)
        !allocate(cbst(dall+maxm*maxm),dsyn(dall),stat=checkstat)
        allocate(cbst(dall+maxm*3), dsyn(dall), stat=checkstat)
        allocate(sigmaT(dall), Tref(dall), datweight(dall), stat=checkstat)
        allocate(Tdata(dall), resbst(dall), fwdT(dall), fwdTvs(dall), fwdTaa(dall), RefTaa(dall), stat=checkstat)
        allocate(resbst_iso(dall), stat=checkstat)
        allocate(resSigma(dall), stat=checkstat)
        if (iso_mod) then
            allocate(Lsen_Gsc(1,1,1), stat=checkstat)
        else
            allocate(Lsen_Gsc(nx*ny,kmaxRc,nz-1), stat=checkstat)
        endif
        if(checkstat > 0) stop 'Lsen_Gsc allocation failed'
        if (use_sparse_g) then
            ! FMM backend writes sparse triplets directly; dense G diagnostic
            ! matrices are legacy-only and are huge for real data.
            allocate(GVs(1,1), GGc(1,1), GGs(1,1), stat=checkstat)
        elseif (iso_mod) then
            allocate(GVs(dall,maxvp), GGc(1,1), GGs(1,1), stat=checkstat)
        else
            allocate(GVs(dall,maxvp), GGc(dall,maxvp), GGs(dall,maxvp), stat=checkstat)
        endif
        if(checkstat > 0) stop 'G sensitivity array allocation failed'
        if (iso_mod) then
            allocate(gcf(1,1,1), gsf(1,1,1), stat=checkstat)
        else
            allocate(gcf(nx-2,ny-2,nz-1),gsf(nx-2,ny-2,nz-1), stat=checkstat)
        endif
        if(checkstat > 0) stop 'anisotropic model array allocation failed'
        allocate( tRcV((nx-2)*(ny-2),kmaxRc), stat=checkstat)
        !----------------------------------------------------------------------!
        !  read reference isotropic model: vsf
        open(11, file='MOD', status='old')
        vsf=0
        read(11,*) (depz(i),i=1,nz)
        do k = 1,nz
            do j = 1,ny
                read(11,*)(vsf(i,j,k),i=1,nx)
            enddo
        enddo
        close(11)
        write(*,*) ' grid points in depth direction:(km)'
        write(*,'(50f7.2)') depz

        call  CalRmax(nz,depz,minthk,rmax)
        fmm_threads = omp_get_max_threads()
        fmm_threads_text = ''
        call get_environment_variable('DAZI_FMM_THREADS', fmm_threads_text, status=fmm_env_status)
        if (fmm_env_status == 0 .and. len_trim(fmm_threads_text) > 0) then
            read(fmm_threads_text, *, iostat=read_status) fmm_threads
            if (read_status /= 0 .or. fmm_threads < 1) fmm_threads = omp_get_max_threads()
        endif
        write(6,'(A,A)') 'Ray backend: ', trim(adjustl(ray_backend))
        write(66,'(A,A)') 'Ray backend: ', trim(adjustl(ray_backend))
        if (trim(adjustl(ray_backend)) == 'fmm_parallel') then
            write(6,'(A,I0)') 'FMM ray-tracing OpenMP threads: ', fmm_threads
            write(66,'(A,I0)') 'FMM ray-tracing OpenMP threads: ', fmm_threads
        endif
        azi_vs_mode = 'joint'
        call get_environment_variable('DAZI_AZI_VS_MODE', azi_vs_mode, status=azi_vs_mode_status)
        if (azi_vs_mode_status /= 0 .or. len_trim(azi_vs_mode) == 0) azi_vs_mode = 'joint'
        azi_vs_mode = adjustl(azi_vs_mode)
        azi_freeze_vs = (trim(azi_vs_mode) == 'freeze' .or. trim(azi_vs_mode) == 'frozen' .or. &
                         trim(azi_vs_mode) == 'ani_only' .or. trim(azi_vs_mode) == 'residual')
        if (.not. iso_mod) then
            if (azi_freeze_vs) then
                write(6,'(A)')  'Anisotropic Vs mode: freeze (solve residual with Gc/Gs only; Vs is not updated)'
                write(66,'(A)') 'Anisotropic Vs mode: freeze (solve residual with Gc/Gs only; Vs is not updated)'
            else
                write(6,'(A)')  'Anisotropic Vs mode: joint (solve dVs, Gc, and Gs together)'
                write(66,'(A)') 'Anisotropic Vs mode: joint (solve dVs, Gc, and Gs together)'
            endif
        endif
        !----------------------------------------------------------------------!
        !                     iteration part
        !----------------------------------------------------------------------!
        RefTaa=0
        ! iso_mod=.false.
        open(34, file='IterVel.out')
        open(67,file='costVSiter.dat')
        ncost = 68
        open(ncost, file='cost_terms.dat', status='replace', action='write')
        write(ncost,'(a)') '# iter  data_misfit_L2  data_misfit_L2_sq  model_L2  model_L2_sq  mae_res  std_res  rms_res'
        if (vorotomo /= 0) then
            fid_misfit = 99
            open(fid_misfit, file='misfit_sorted.dat', status='unknown')
            fid_misfit_voro = 98
            open(fid_misfit_voro, file='misfit_voro_raw_sanitized.dat', status='unknown')
        endif
        dvstdall = 0

        do iter=1, maxiter
            timer_iter_start = omp_get_wtime()
            timer_g_matrix = 0.0d0
            timer_residual_weight = 0.0d0
            timer_solver = 0.0d0
            timer_post_update = 0.0d0
            timer_raypath_output = 0.0d0
            iter_mod=mod(iter, 2)
            write(6,*)  ' -----------------------------------------------------------'
            write(66,*) ' -----------------------------------------------------------'

            if (iso_mod) then
                iso_inv=.true.
                write(66,*)iter,'th iteration, invert for isotropic Vs para.'
                write(6,*) iter,'th iteration, invert for isotropic Vs para.'
                ! maxm =  maxvp
            else
                iso_inv=.false.
                write(66,*)iter,'th iteration, invert for dVs, Gc, Gs '
                write(6,*) iter,'th iteration, invert for dVs, Gc, Gs '
                ! maxm =  maxvp*2 ! joint
                ! maxm =  maxvp*3
            endif
            write(6,*)  ' -----------------------------------------------------------'
            write(66,*) ' -----------------------------------------------------------'
            !----------------------------------------------------------------------!
            ! compute G matrix based on the sensitivity kernel and ray-tracing
            ! forward calculation of traveltime misfit based on the reference iso- or aniso- model.
            ! Refwritepath = 0 ! switch of write ray path
            dsyn = 0
            GGc = 0
            GGs = 0
            GVs = 0
            tRcV = 0
            iw = 0
            rw = 0.0
            col = 0
            nar = 0
            Lsen_Gsc=0.0
            timer_section_start = omp_get_wtime()
            if (iso_inv) then
                write(66,*) 'Calculating isotropic sensitivity kernel matrix G_iso.'
                write(6,*) 'Calculating isotropic sensitivity kernel matrix G_iso.'
                if (trim(adjustl(ray_backend)) == 'fmm_parallel') then
                    call CalSurfGParallelIso(nx, ny, nz, maxvp, vsf, iw, rw, col, dsyn, &
                    goxd, gozd, dvxd, dvzd, kmaxRc, tRc, periods, depz, minthk, &
                    scxf, sczf, rcxf, rczf, nrc1, nsrc1, wavetype, igrt, &
                    kmax, nsrc, nrc, nar, maxnar, 0, dall, tRcV, fmm_threads)
                else
                    call CalSurfG(nx,ny,nz,maxvp,vsf,iw,rw,col,dsyn,&
                    GVs, dall,&
                    goxd, gozd, dvxd, dvzd, kmaxRc, tRc, periods, depz, minthk,&
                    scxf, sczf, rcxf, rczf, nrc1, nsrc1, kmax, nsrc, nrc, nar, 0)  ! No output raypath for the every iteration
                endif
            else
                write(66,*) 'Calculating azimuthal anisotropic sensitivity kernel matrix G_ani.'
                write(6,*) 'Calculating azimuthal anisotropic sensitivity kernel matrix G_ani.'
                if (trim(adjustl(ray_backend)) == 'fmm_parallel') then
                    call CalSurfGParallelAniso(nx, ny, nz, maxvp, vsf, iw, rw, col, dsyn, &
                    Lsen_Gsc, rmax, tRcV, goxd, gozd, dvxd, dvzd, kmaxRc, tRc, periods, depz, minthk, &
                    scxf, sczf, rcxf, rczf, nrc1, nsrc1, kmax, nsrc, nrc, nar, maxnar, 0, dall, fmm_threads)
                else
                    call CalSurfGAnisoJoint(nx, ny, nz, maxvp, vsf, iw, rw, col, dsyn, &
                    GVs, GGc, GGs, Lsen_Gsc, dall, rmax, tRcV, &
                    goxd, gozd, dvxd, dvzd, kmaxRc, tRc, periods,depz, minthk,&
                    scxf, sczf, rcxf, rczf, nrc1, nsrc1, kmax, nsrc, nrc, nar, 0) ! No output raypath for the every iteration
                endif
            endif
            timer_section_end = omp_get_wtime()
            timer_g_matrix = timer_section_end - timer_section_start
            write(*,*) 'Number of non-zero element (nar) of G is ',nar
            write(*,*) 'Max Number of sparse matrix non-zero element ',maxnar
            if (nar.gt.maxnar) then
                write(6,*) 'nar > maxnar, please increase your sparsefrac'
                stop 'nar > maxnar, please increase your sparsefrac'
            endif
            write(*,*) ' Finish G matrix calculation.'
            !----------------------------------------------------------------------!
            ! output the corresponding isotropic phase velocity map of the MOD.
            if (iter.eq.1)then
                open(77,file='period_phaseVMOD.dat')
                call WTPeriodPhaseV(nx, ny, gozd, goxd, dvzd, dvxd, kmaxRc, tRc, tRcV, 77)
            end if
            timer_section_start = omp_get_wtime()
            !----------------------------------------------------------------------!
            ! calculate data residual based on the reference iso- or aniso- model.
            ! cbst: travel time residual.
            ! dsyn: predicated isotropic travel time
            ! Tdata: traveltime residual, used in the inversion
            cbst=0
            Tdata=0
            sigmaT=0
            Tref=0
            ! For both joint and iso inversion. Data fitting: obs(i) - Tref-iso(i)
            do i=1,dall
                Tref(i)=dsyn(i)
                cbst(i)=obst(i) - Tref(i)
                Tdata(i)=cbst(i)
            enddo

            ! Statistic of the data residual based on the reference model.
            mean=sum(cbst(1:dall))/dall
            std_devs=sqrt(sum((cbst(1:dall)-mean)**2)/dall)
            meanAbs=sum(abs(cbst(1:dall)))/dall

            write(6 ,'(a, f12.4,a,f10.2,a,f10.2,a)') '  Before Inversion: abs mean, std, RMS of Res:', meanAbs,' s ',&
            std_devs,' s ', dnrm2(dall,cbst,1)/sqrt(real(dall)),' s'
            write(66,'(a, f12.4,a,f10.2,a,f10.2,a)') '  Before Inversion: abs mean, std, RMS of Res:', meanAbs,' s ',&
            std_devs,' s ', dnrm2(dall,cbst,1)/sqrt(real(dall)),' s'

            !----------------------------------------------------------------------!
            ! data weight
            ! datweight=0.0
            ! if (iso_inv) then
            !     thresholdVs=0.2
            !     do i = 1,dall
            !     datweight(i) = 1.0
            !     if(abs(cbst(i)) > thresholdVs)  datweight(i) = exp(-(abs(cbst(i))-thresholdVs))
            !     cbst(i) = cbst(i)*datweight(i)
            !     enddo
            ! else
            !     ! thresholdVs=0.5
            !     threshold=std_devs*2.
            !     do i = 1,dall
            !     datweight(i) = 1.0
            !     if(abs(cbst(i)) > threshold)  datweight(i) = exp(-(abs(cbst(i))-threshold))
            !     cbst(i) = cbst(i)*datweight(i)
            !     enddo
            ! For Yunnan.
            call CalDdatSigma(dall, obst, cbst, sigmaT, meandeltaT)
            do i=1,dall
                datweight(i)=1/sigmaT(i)
                cbst(i)=cbst(i)*datweight(i)
            enddo

            do i = 1, nar
                rw(i) = rw(i) * datweight(iw(i))
            end do

            meanAbs=sum(abs(cbst(1:dall)))/dall
            write(6 ,'(a, f8.3, a, f8.3,a, f7.3, a)') '  mean data weight:',sum(datweight(1:dall))/dall, &
            ' |  abs data mean with weight:',meanAbs,'s  |  dt/t0:',  meandeltaT*100, ' %'
            write(66,'(a, f8.3, a, f8.3,a, f7.3, a)') '  mean data weight:',sum(datweight(1:dall))/dall, &
            ' |  abs data mean with weight:',meanAbs,'s  |  dt/t0:',  meandeltaT*100, ' %'
            if (iso_mod) then
                norm=0
                do i=1,nar
                norm(col(i))=norm(col(i))+abs(rw(i))
                enddo
            endif

            ! write out res for 1-th and final inversion
            if ((iter .eq. 1) .or. (iter .eq. maxiter)) then
                write(id,'(I2.2)') iter
                filename = 'Traveltime_use_'//TRIM(id)//'th.dat'
                open(88, file=filename)
                write(88,'(a)') 'Dist(km)        T_obs(s)       T_forward(s)         Res(s)        weight        W_Res'
                do i=1,dall
                    write(88,*) dist(i), obst(i), dsyn(i), Tdata(i), datweight(i), cbst(i)
                enddo
                close(88)
            endif
            timer_section_end = omp_get_wtime()
            timer_residual_weight = timer_section_end - timer_section_start
            !----------------------------------------------------------------------!
            ! Setting of iteration.
            !----------------------------------------------------------------------!
            ! A*x = b
            ! m       input      m, the number of rows in A. (data)
            ! n       input      n, the number of columns in A. (model)
            ! count3: number of increase of row number.
            ! nar: number of G which value is not zero, including augmented matrix.
            !      if Gn*m does not contain zeros,  nar =n*m
            ! nar1: nar before regularization
            ! rw(i) gives value of the i-th non-zero value of G
            ! iw(i+1) means this value is in the iw(i+1)-th row.
            ! col(i) means this value is in the col(i)-th column. = iw(i+1+nar)
            !----------------------------------------------------------------------!

    !---------------------------------------------------------------------------------------!
    !---------------------------------------------------------------------------------------!
            timer_section_start = omp_get_wtime()
            if (vorotomo /= 0) then
                if (.not. allocated(misfit_voro)) then
                    allocate(misfit_voro(nrealizations))
                    allocate(acond_voro(nrealizations))
                    allocate(itn_voro(nrealizations))
                    allocate(istop_voro(nrealizations))
                    allocate(misfit_sorted(nrealizations))
                    allocate(idx_sorted(nrealizations))
                end if
                nrow = 0
                ! hvratio = 1.0 !dvxd*(nx-3)*111.19/depz(nz-1)
                maxborder = max(dvxd*(nx-3),dvzd*(ny-3))
                hvratio = maxborder*111.19/depz(nz-1)
                dv = 0
                dvstd = 0
                dvsub = 0
                dvall = 0
                leniw = 2_8*nar
                lenrw = nar
                iw(nar+1:2*nar) = col(1:nar)

                ! Build per-row nnz counts for the DATA block (rows 1..dall).
                ! IMPORTANT: iw(1:nar) must contain row indices in [1..dall] here.
                do i = 1, nar
                    if (iw(i) < 1 .or. iw(i) > dall) then
                        write(*,'(a,1x,i0,1x,a,1x,i0)') 'ERROR: iw(i) row index out of [1..dall] at i=', i, ' row=', iw(i)
                        stop 'iw row index out of range while building nrow (vorotomo)'
                    end if
                    nrow(iw(i)) = nrow(iw(i)) + 1
                end do
                ! Detect integer overflow/wrap: if it happens, nrow could become negative.
                do i = 1, dall
                    if (nrow(i) < 0) stop 'nrow overflow/wrap detected!'
                end do
                ! ------------------------------------------------------------
                ! Sanity checks for CSR-like row structure of G
                ! nrow(r) counts how many nonzeros are in data row r (r=1..dall)
                ! The total nnz implied by nrow should equal nar.
                ! We accumulate in INT8 to avoid overflow when nar is large.
                ! ------------------------------------------------------------
                
                nnz_check = 0_8
                do i = 1, dall
                    nnz_check = nnz_check + int(nrow(i), kind=8)
                end do
                if (nnz_check /= nar) then
                    write(*,'(a,1x,i0,1x,a,1x,i0)') 'WARNING: sum(nrow)=', nnz_check, ' != nar=', nar
                    write(66,'(a,1x,i0,1x,a,1x,i0)') 'WARNING: sum(nrow)=', nnz_check, ' != nar=', nar
                endif
                ! ------------------------------------------------------------
                ! Diagnostic: verify (iw, col, rw) triplets are consistent.
                ! - iw(1:nar) must be row ids in [1..dall]
                ! - col(1:nar) must be model ids in [1..maxm]
                ! - Prefer (iw) non-decreasing for CSR-like access; if not,
                !   voronoi assembly that assumes row-contiguity may be wrong.
                ! ------------------------------------------------------------
                if (iter == 1) then
                    call check_triplets(iw, col, rw, nar, dall, maxm, 66)
                endif

                do i =1, nrealizations
                    misfit_voro(i) = 0.0
                    acond_voro(i) = 0.0
                    itn_voro(i) = 0
                    istop_voro(i) = 0
                enddo

                write(*,'(a,i0)') 'OMP max threads: ', omp_get_max_threads()
                write(66,'(a,i0)') 'OMP max threads: ', omp_get_max_threads()
                write(*,'(a)') 'Seed rule: iseed(1:3)=(/38,62,346/), iseed(4)=2*(iproj+iter)+1 (mod 4095),' // &
                ' with stream offsets {11,23,37,59}'
                write(66,'(a)') 'Seed rule: iseed(1:3)=(/38,62,346/), iseed(4)=2*(iproj+iter)+1 (mod 4095),' // &
                ' with stream offsets {11,23,37,59}'
                write(*,*) "Doing parallel Poisson Voronoi inversion with ", nrealizations, " realizations."
                write(66,*) "Doing parallel Poisson Voronoi inversion with ", nrealizations, " realizations."
                !$omp parallel default(none) &
                !$omp shared(iso_inv,leniw,lenrw,col,nrow,rw,cbst,goxd,dvxd,gozd,dvzd,depz,maxvp,maxm) &
                !$omp shared(nx,ny,nz,dall,ncell,acell,nzrand,spfra,hvratio,iaratio,damp,nrealizations,iter,dvall) &
                !$omp shared(misfit_voro,acond_voro,itn_voro,istop_voro) &
                !$omp private(iproj,rnorm_out,acond_out,itn_out,istop_out,dvsub_local)
                allocate(dvsub_local(maxm))
                !$omp do schedule(static)
                do iproj = 1,nrealizations
                    dvsub_local = 0.0
                    if (iso_inv) then
                        ! write(*,*) "  Voronoi ISO projection realization ", iproj
                        call voronoiproj_layer_iso(leniw,lenrw,col,nrow,rw,cbst,goxd,dvxd,gozd,dvzd,depz,&
                                nx,ny,nz,dall,ncell,acell,nzrand,spfra,hvratio,damp,iproj,iter,dvsub_local,&
                                rnorm_out, acond_out, itn_out, istop_out)
                        dvall((iproj-1)*maxvp+1:iproj*maxvp) = dvsub_local(1:maxvp)
                        misfit_voro(iproj) = rnorm_out
                        acond_voro(iproj)  = acond_out
                        itn_voro(iproj)    = itn_out
                        istop_voro(iproj)  = istop_out
                    else
                        ! write(*,*) "  Voronoi AZI projection realization ", iproj
                        call voronoiproj_layer_azi(leniw,lenrw,col,nrow,rw,cbst,goxd,dvxd,gozd,dvzd,depz,&
                                nx,ny,nz,dall,ncell,acell,nzrand,spfra,hvratio,iaratio,damp,iproj,iter,dvsub_local,&
                                rnorm_out, acond_out, itn_out, istop_out)
                        dvall((iproj-1)*maxm+1:iproj*maxm) = dvsub_local(1:maxm)
                        misfit_voro(iproj) = rnorm_out
                        acond_voro(iproj)  = acond_out
                        itn_voro(iproj)    = itn_out
                        istop_voro(iproj)  = istop_out
                    endif
                enddo
                !$omp end do
                if (allocated(dvsub_local)) then
                    deallocate(dvsub_local)
                endif
                !$omp end parallel
                
                ! ============================================================
                ! Weighted ensemble fusion (misfit-sorted) + logging
                ! Drop-in replacement block for Main_Jt.f90
                !
                ! Assumptions:
                !   - misfit_voro(1:nrealizations) contains per-realization misfit (smaller is better)
                !   - dvall stores per-realization solution vectors in contiguous blocks:
                !       iso_inv=.true.  -> block size = maxvp
                !       iso_inv=.false. -> block size = maxm
                !   - idx_sorted(1:nrealizations) is available (integer)
                !   - median(x, n) function exists and returns median of x(1:n)
                !
                ! Notes:
                !   - Weights are defined on the SORTED RANK i (i=1 is best misfit).
                !     Therefore, use idx_sorted(i) when accessing dvall blocks.
                !   - For file output, we print: rank, realization_id, misfit_sorted, weight_rank
                ! ============================================================
                ! ------------------------------
                ! Misfit health check (before sorting)
                ! ------------------------------

                ! ---- log (rank, real_id, misfit, weight) ----
                write(fid_misfit_voro,'(a,i0)') 'Iteration(raw) ', iter
                write(fid_misfit_voro,'(a)') 'real_id        misfit           acond          itn  istop'
                do i=1,nrealizations
                    write(fid_misfit_voro,'(i7,1x,es14.6,1x,es14.6,1x,i6,1x,i6)') &
                        i, misfit_voro(i), acond_voro(i), itn_voro(i), istop_voro(i)
                end do
                write(fid_misfit_voro,'(a)') ''
                ! ---- basic stats ----
                ! n_nan  = count(misfit_voro /= misfit_voro)   ! NaN check: NaN /= NaN is true
                ! n_nan = count(ieee_is_nan(misfit_voro))
                n_nan = count(misfit_voro >= MISFIT_BAD*0.5) ! when ieee_is_nan is unavailable
                mf_min = minval(misfit_voro)
                mf_max = maxval(misfit_voro)
                n_clip = count(misfit_voro >= MISFIT_BAD*0.5)

                write(*,'(a,i0,a,i0,a,es12.4,a,es12.4,a,i0)') &
                    'Iter=', iter, '  NaN=', n_nan, '  min=', mf_min, '  max=', mf_max, '  clip=', n_clip
                ! ---- prepare sorted misfit and indices ----
                ! ---- sanitize misfit before sorting (important under OMP / LSMR failures) ----
                do i = 1, nrealizations
                    if (.not.(misfit_voro(i) >= 0.0 .and. misfit_voro(i) < MISFIT_BAD)) misfit_voro(i) = MISFIT_BAD

                    ! mark bad terminations (istop=4,5) and (istop=3 and damp==0)
                    ! istop=4,5: solver did not converge
                    ! istop=3 and damp==0: solver did not converge
                    ! istop=3 and damp\=0: solver converged under the damping threshold
                    if (istop_voro(i) == 4 .or. istop_voro(i) == 5) misfit_voro(i) = MISFIT_BAD
                    if (abs(damp) < 1e-8 .and. istop_voro(i) == 3) misfit_voro(i) = MISFIT_BAD
                    ! ---- additional hard gate (minimal): drop numerically ill/degenerate solutions ----
                    if (acond_voro(i) > COND_MAX) misfit_voro(i) = MISFIT_BAD
                    if (itn_voro(i)   < ITN_MIN)  misfit_voro(i) = MISFIT_BAD

                    ! ! ---- normalize misfit ----
                    ! misfit_voro(i) = misfit_voro(i)/sqrt(real(dall))
                    ! !
                enddo

                write(fid_misfit_voro,'(a,i0)') 'Iteration(sanitized) ', iter
                write(fid_misfit_voro,'(a)') 'real_id        misfit           acond          itn  istop'
                do i=1,nrealizations
                    write(fid_misfit_voro,'(i7,1x,es14.6,1x,es14.6,1x,i6,1x,i6)') &
                        i, misfit_voro(i), acond_voro(i), itn_voro(i), istop_voro(i)
                end do
                write(fid_misfit_voro,'(a)') ''

                ! ---- sort misfit and get realization indices ----
                misfit_sorted = misfit_voro
                call index_sort(misfit_sorted, idx_sorted)

                ! ! ---- choose top-k to keep ----
                ! n_keep = max(1, int(real(nrealizations) * 0.9))
                ! if (n_keep > nrealizations) n_keep = nrealizations
                ! ---- adaptive choose top-k to keep ----
                ! Step-1: keep only numerically-valid solutions
                n_good = count(misfit_voro < MISFIT_BAD*0.5)
                if (n_good < 1) n_good = 1
                ! Step-2: robust tail cut on sorted misfit among good solutions
                ! Use IQR gate: misfit <= median + alpha * IQR
                ! (misfit_sorted is already sorted ascending at this point)
                alpha = 2.0
                kmin  = max(5, int(0.05*real(nrealizations)))
                kmin  = min(kmin, n_good)

                i_q1 = max(1, int(0.25*real(n_good)))
                i_q3 = max(1, int(0.75*real(n_good)))
                q1   = misfit_sorted(i_q1)
                q3   = misfit_sorted(i_q3)
                iqr  = max(1.0e-12, q3 - q1)

                k0 = 0
                do i = 1, n_good
                    if (misfit_sorted(i) <= median(misfit_sorted(1:n_good), n_good) + alpha*iqr) then
                        k0 = k0 + 1
                    else
                        exit
                    end if
                end do
                k0 = max(kmin, min(k0, n_good))

                ! Step-3: final n_keep determined by cumulative weight coverage
                n_keep = k0
                ! n_keep = max(1, min(n_good, int(real(nrealizations)*0.9)))

                ! ---- build rank weights (allocate per-iter is OK) ----
                if (allocated(w_rank)) deallocate(w_rank)
                allocate(w_rank(nrealizations))
                w_rank = 0.0

                eps = median(misfit_sorted(1:n_keep), n_keep)**2
                do i = 1, n_keep
                    w_rank(i) = 1.0 / (eps + misfit_sorted(i)**2)
                end do
                sumw = sum(w_rank(1:n_keep))
                if (sumw <= 0.0) then
                    w_rank(1:n_keep) = 1.0 / real(n_keep)
                else
                    w_rank(1:n_keep) = w_rank(1:n_keep) / sumw
                end if
                ! ---- refine n_keep by cumulative weight coverage (95%) ----
                csum = 0.0
                nkeep2 = n_keep
                do i = 1, n_keep
                    csum = csum + w_rank(i)
                    if (csum >= 0.95) then
                        nkeep2 = i
                        exit
                    end if
                end do
                if (nkeep2 < kmin) nkeep2 = kmin
                if (nkeep2 > n_keep) nkeep2 = n_keep
                n_keep = nkeep2

                ! renormalize after truncation
                sumw = sum(w_rank(1:n_keep))
                if (sumw <= 0.0) then
                    w_rank(1:n_keep) = 1.0 / real(n_keep)
                else
                    w_rank(1:n_keep) = w_rank(1:n_keep) / sumw
                end if
                ! ---- map to realization id weights (recommended) ----
                if (allocated(w_by_real)) deallocate(w_by_real)
                allocate(w_by_real(nrealizations))
                w_by_real = 0.0
                do i = 1, n_keep
                    idx = idx_sorted(i)
                    w_by_real(idx) = w_rank(i)
                end do

                ! ---- log (rank, real_id, misfit, weight) ----
                write(fid_misfit, '(a,i0)') 'Iteration ', iter
                write(fid_misfit, '(a)') 'rank  real_id        misfit           weight'
                do i = 1, nrealizations
                    idx = idx_sorted(i)
                    write(fid_misfit, '(i6,1x,i7,1x,es14.6,1x,es14.6)') i, idx, misfit_sorted(i), w_rank(i)
                end do
                write(fid_misfit, '(a)') ''

                ! ---- weighted mean/std ----
                if (iso_inv) then
                    dv(1:maxvp) = 0.0
                    dvstd(1:maxvp) = 0.0
                    do i = 1, n_keep
                        idx = idx_sorted(i)
                        dv(1:maxvp) = dv(1:maxvp) + w_rank(i) * dvall((idx-1)*maxvp+1 : idx*maxvp)
                    end do
                    do i = 1, n_keep
                        idx = idx_sorted(i)
                        dvstd(1:maxvp) = dvstd(1:maxvp) + w_rank(i) * &
                            (dvall((idx-1)*maxvp+1 : idx*maxvp) - dv(1:maxvp))**2
                    end do
                    dvstd(1:maxvp) = sqrt(dvstd(1:maxvp))
                else
                    dv(1:maxm) = 0.0
                    dvstd(1:maxm) = 0.0
                    do i = 1, n_keep
                        idx = idx_sorted(i)
                        dv(1:maxm) = dv(1:maxm) + w_rank(i) * dvall((idx-1)*maxm+1 : idx*maxm)
                    end do
                    do i = 1, n_keep
                        idx = idx_sorted(i)
                        dvstd(1:maxm) = dvstd(1:maxm) + w_rank(i) * &
                            (dvall((idx-1)*maxm+1 : idx*maxm) - dv(1:maxm))**2
                    end do
                    dvstd(1:maxm) = sqrt(dvstd(1:maxm))
                end if
                ! ============================================================
                !!!!!!!!!!!!!!!!!!!!!!!
                ! ! calculate mean
                ! do iproj = 1,nrealizations
                !     if (iso_inv) then
                !         dvsub(1:maxvp) = dvall((iproj-1)*maxvp+1:iproj*maxvp)!:,iproj)
                !     else
                !         dvsub(1:maxm) = dvall((iproj-1)*maxm+1:iproj*maxm)!:,iproj)
                !     endif
                !     dv = dv+dvsub
                !     ! dvstd = dvstd+dvsub**2
                ! enddo
                ! dv = dv/nrealizations
                ! ! dvstd = sqrt(dvstd/nrealizations-dv**2)
                ! ! dvstdall = dvstdall+dvstd
                ! ! calculate standard derivation
                ! do iproj = 1,nrealizations
                !     if (iso_inv) then
                !         dvsub(1:maxvp) = dvall((iproj-1)*maxvp+1:iproj*maxvp)!:,iproj)
                !     else
                !         dvsub(1:maxm) = dvall((iproj-1)*maxm+1:iproj*maxm)!:,iproj)
                !     endif
                !     dvstd = dvstd+(dvsub - dv)**2
                ! enddo
                ! ! std of each iteration
                ! dvstd = sqrt(dvstd/(nrealizations - 1)) 
                !!!!!!!!!!!!!!!!!!!!!!!
                ! std sum of all iterations

                ! --- UQ accumulation (law of total variance) ---
                ! dv is the (possibly misfit-weighted) ensemble mean for this iteration
                ! dvstd is the (possibly misfit-weighted) ensemble std for this iteration
                dv_mean_accum        = dv_mean_accum        + dv
                dv_mean2_accum       = dv_mean2_accum       + dv**2
                dv_var_within_accum  = dv_var_within_accum  + dvstd**2
                nuq_iter             = nuq_iter + 1
                ! update running dvstdall so intermediate outputs are meaningful
                if (nuq_iter > 0) then
                    dvstdall = sqrt( max(0.0, (dv_var_within_accum/nuq_iter) + &
                                        (dv_mean2_accum/nuq_iter - (dv_mean_accum/nuq_iter)**2)) )
                endif
            else ! not vorotomo, conventional tikhonov regularization
                !----------------------------------------------------------------------!
                nar1=nar
                lameGcs=weightGcs
                ! lameVs=dnrm2(dall,cbst,1)**2/dall*weightVs 
                lameVs=weightVs
                count3=0
                if (azi_freeze_vs .and. (.not. iso_inv)) then
                    do i=1,nar
                        if (col(i) >= 1 .and. col(i) <= maxvp) rw(i) = 0.0
                    enddo
                    write(* ,'(a)') '  Frozen-Vs ANI mode: zeroed dVs sensitivity columns before regularization.'
                    write(66,'(a)') '  Frozen-Vs ANI mode: zeroed dVs sensitivity columns before regularization.'
                endif
                if (iso_inv) then
                    call TikhonovRegularization(nx, ny, nz, maxvp, dall, nar, rw, iw, col, count3, iso_inv, lameGcs, lameVs, maxnar)
                else
                    call        TikhRegul_joint(nx, ny, nz, maxvp, dall, nar, rw, iw, col, narVs, count3, lameGcs, lameVs, maxnar)
                endif

                forall(i=1:count3)
                cbst(dall+i)=0
                end forall
                !----------------------------------------------------------------------!
                !----------------------------------------------------------------------!
                if (nar > maxnar) stop 'increase sparsity fraction(spfra)'
                write(*  ,'(a,3f8.2)') '  damp,  lamebda Gsc, lamebda Vs: ', damp, lameGcs, lameVs
                write(66 ,'(a,3f8.2)') '  damp,  lamebda Gsc, lamebda Vs: ', damp, lameGcs, lameVs

                m = dall + count3 ! data number, rows number of G
                n = maxm          ! model number, columns of G
                do i=1,nar
                iw(nar+i)=col(i)
                enddo
                leniw = 2*nar
                lenrw = nar
                dv = 0
                ! Control of LSQR
                ! atol:   estimate of error to G, accurate to about 4 digits, set btol = 1.0e-4. (default AA: 1e-4)
                ! btol:   estimate of error to data, accurate to about 4 digits, set btol = 1.0e-4. (default AA: 1e-4)
                ! conlim: the apparent condition number of the matrix Abar, controls the amp of Gc,s
                !         conlim and damp may be used separately or together to regularize ill-conditioned systems.
                ! itnlim: an upper limit on the number of iterations.
                if (iso_inv) then
                    atol = 1e-3
                    btol = 1e-3
                    conlim = 1200
                    itnlim = 1000
                    localSize = n/4
                else
                    atol = 1e-5
                    btol = 1e-4
                    conlim = 200
                    itnlim = 500
                    localSize = 10
                endif

                istop = 0
                anorm = 0.0
                acond = 0.0
                arnorm = 0.0
                xnorm = 0.0

                lsmr_backend = ''
                call get_environment_variable('DAZI_MAIN_LSMR_BACKEND', lsmr_backend, status=lsmr_env_status)
                if (lsmr_env_status /= 0) then
                    call get_environment_variable('DAZI_LSMR_BACKEND', lsmr_backend, status=lsmr_env_status)
                endif

                if (trim(adjustl(lsmr_backend)) == 'scipy') then
                    write(lsmr_tag,'("main_i",I0)') iter
                    call scipy_lsmr_solve(lsmr_tag, m, n, leniw, lenrw, iw, rw, cbst, damp, &
                    atol, btol, conlim, itnlim, dv, istop, itn, anorm, acond, rnorm, arnorm, xnorm)
                else
                    call LSMR(m, n, leniw, lenrw, iw, rw, cbst, damp,&
                    atol, btol, conlim, itnlim, localSize, nout,&
                    dv, istop, itn, anorm, acond, rnorm, arnorm, xnorm)
                endif
                if(istop==3) THEN
                    write(* , '(a)') '  istop = 3, large condition number, LSMR failed'
                    write(66, '(a)') '  istop = 3, large condition number, LSMR failed'
                endif
                write(*,'(a)') '  Finish LSMR.......'
                write(*, '(a, i7)')  '  itn=               ',itn
                write(66,'(a, i7)')  '  itn=               ',itn
                write(*,'(a, f7.1)') '  L2 norm of A=      ',anorm
                write(*, '(a, f7.1)')'  Condition NO. of A=',acond
                write(66,'(a, f7.1)')'  Condition NO. of A=',acond
                write(*,'(a, f7.1)') '  rnorm=             ',rnorm
                write(*,'(a, f7.1)') '  arnorm=            ',arnorm
                write(*,'(a, f7.3)') '  norm of dv =       ',xnorm
            endif
            timer_section_end = omp_get_wtime()
            timer_solver = timer_section_end - timer_section_start
            if (azi_freeze_vs .and. (.not. iso_inv)) then
                dv(1:maxvp) = 0.0
                dvstd(1:maxvp) = 0.0
                write(* ,'(a)') '  Frozen-Vs ANI mode: forced dVs update to zero.'
                write(66,'(a)') '  Frozen-Vs ANI mode: forced dVs update to zero.'
            endif
            !----------------------------------------------------------------------!
            ! Cost terms for BOTH regular and Voronoi inversions
            !   data misfit: || W (G dv - Tdata) ||_2  (using existing helper routines)
            !   model norm : || dv ||_2
            resbst   = 0.0
            timer_section_start = omp_get_wtime()
            resbst   = 0.0
            resSigma = 0.0
            fwdTvs   = 0.0
            fwdTaa   = 0.0

            if (vorotomo /= 0) then
                residual_m = dall
                residual_n = maxm
            else
                residual_m = m
                residual_n = n
            endif

            if (trim(adjustl(ray_backend)) == 'fmm_parallel') then
                call CalSparseReslNorm(residual_m, residual_n, maxvp, dall, leniw, lenrw, iw, rw, dv(1:maxm), &
                                    datweight, Tdata, fwdTvs, fwdTaa, resbst, iso_inv)
            elseif (iso_inv) then
                call CalVsReslNorm(maxvp, dall, GVs, dv(1:maxvp), datweight, Tdata, fwdTvs, resbst)
            else
                VsInv = dv(1:maxvp)
                GcInv = dv(1+maxvp:2*maxvp)
                GsInv = dv(1+2*maxvp:3*maxvp)
                call CalReslNormJoint(maxvp, dall, GVs, GGc, GGs, VsInv, GcInv, GsInv, &
                                    datweight, Tdata, fwdTvs, fwdTaa, resbst)
            endif

            resSigma(1:dall) = resbst(1:dall) * datweight(1:dall)
            data_misfit_norm  = dnrm2(dall, resSigma, 1)
            data_misfit_norm2 = data_misfit_norm * data_misfit_norm

            model_l2_norm  = dnrm2(maxm, dv, 1)
            model_l2_norm2 = model_l2_norm * model_l2_norm
            ! ---- unweighted residual summary (matches the log line "After Inversion") ----
            mean = sum(resbst(1:dall))/real(dall)
            mae_res = sum(abs(resbst(1:dall))) / real(dall)
            std_res = sqrt( sum((resbst(1:dall) - mean)**2) / real(dall) )
            rms_res = dnrm2(dall, resbst, 1) / sqrt(real(dall))

            write(ncost,'(i6, 7(1x, es16.6))') iter, data_misfit_norm, data_misfit_norm2, &
                                  model_l2_norm, model_l2_norm2, mae_res, std_res, rms_res
            !----------------------------------------------------------------------!
            ! Based on perturbation, construct isotropic or anisotropic velocity model
            ! For vsRela will be used in subroutine: CalSurfGAniso, it will be important to get a right value.
            ! vsf initially comes from MOD, then can be updated after each iteration.
            if (iso_inv) then
                do k=1,nz-1
                do j=1,ny-2
                do i=1,nx-2
                pertV = dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Clip the perturbation to [-0.5, 0.5] km/s
                if (pertV.ge.0.500) pertV=0.500
                if (pertV.le.-0.500) pertV=-0.500
                if (abs(pertV)<1e-5) pertV=0.0
                dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i) = pertV
                vsf(i+1,j+1,k)=vsf(i+1,j+1,k)+dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)

                if(vsf(i+1,j+1,k).lt.Minvel) vsf(i+1,j+1,k)=Minvel
                if(vsf(i+1,j+1,k).gt.Maxvel) vsf(i+1,j+1,k)=Maxvel
                enddo
                enddo
                enddo
            else
                gcf=0
                gsf=0
                do k=1,nz-1
                do j=1,ny-2
                do i=1,nx-2

                pertV = dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! Clip the perturbation to [-0.5, 0.5] km/s
                if (pertV.ge. 0.5)   pertV=0.5
                if (pertV.le.-0.5)   pertV=-0.5
                if (abs(pertV)<1e-5) pertV=0.0
                dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i) = pertV
                vsf(i+1,j+1,k)=vsf(i+1,j+1,k)+dv((k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)

                if(vsf(i+1,j+1,k).lt.Minvel) vsf(i+1,j+1,k)=Minvel
                if(vsf(i+1,j+1,k).gt.Maxvel) vsf(i+1,j+1,k)=Maxvel

                gcf(i,j,k)=dv(maxvp  +(k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)
                gsf(i,j,k)=dv(maxvp*2+(k-1)*(nx-2)*(ny-2)+(j-1)*(nx-2)+i)
                enddo
                enddo
                enddo
            endif
            !----------------------------------------------------------------------!
            ! result analysis
            if (iso_inv)then
                mindVs=minval(dv(1:maxvp))
                maxdVs=maxval(dv(1:maxvp))
                meadVs=sum(abs(dv(1:maxvp)))/maxvp
                write(6 ,'(a,3f10.4)')  '  min  max and abs mean  dVs (km/s)', mindVs, maxdVs, meadVs
                write(66,'(a,3f10.4)')  '  min  max and abs mean  dVs (km/s)', mindVs, maxdVs, meadVs
                do k=1,nz-1
                VariVs=sum(abs(dv((k-1)*(nx-2)*(ny-2)+1:k*(nx-2)*(ny-2) )))/((nx-2)*(ny-2))
                write(66,'(a,f5.1,a,f5.1,a,f10.4)') '  Z ',depz(k),' - ',depz(k+1),&
                ' km  abs mean dVs (km/s)',VariVs
                write(6,'(a,f5.1,a,f5.1,a,f10.4)')  '  Z ',depz(k),' - ',depz(k+1),&
                ' km  abs mean dVs (km/s)',VariVs
                enddo
            else
                mindVs=minval(dv(1:maxvp))
                maxdVs=maxval(dv(1:maxvp))
                meadVs=sum(abs(dv(1:maxvp)))/maxvp
                write(6 ,'(a,3f10.4)')  '  min  max and abs mean  dVs (km/s)', mindVs, maxdVs, meadVs
                write(66,'(a,3f10.4)')  '  min  max and abs mean  dVs (km/s)', mindVs, maxdVs, meadVs
                minGc=minval(dv(maxvp+1:maxvp+maxvp))*100
                maxGc=maxval(dv(maxvp+1:maxvp+maxvp))*100
                meaGc=sum(abs(dv(maxvp+1:maxvp+maxvp)))/maxvp*100
                minGs=minval(dv(maxvp*2+1:maxvp+maxvp*2))*100
                maxGs=maxval(dv(maxvp*2+1:maxvp+maxvp*2))*100
                meaGs=sum(abs(dv(maxvp*2+1:maxvp+maxvp*2)))/maxvp*100
                write(6 ,'(a,3f10.4)')  '  min  max and abs mean   Gc/L (%) ', minGc, maxGc, meaGc
                write(66,'(a,3f10.4)')  '  min  max and abs mean   Gc/L (%) ', minGc, maxGc, meaGc
                write(6 ,'(a,3f10.4)')  '  min  max and abs mean   Gs/L (%) ', minGs, maxGs, meaGs
                write(66,'(a,3f10.4)')  '  min  max and abs mean   Gs/L (%) ', minGs, maxGs, meaGs
                do k=1, nz-1
                VariVs=sum(abs(dv((k-1)*(nx-2)*(ny-2)+1:k*(nx-2)*(ny-2) )))/((nx-2)*(ny-2))
                VariGc=sum(abs(gcf(1:nx-2,1:ny-2,k)))/((nx-2)*(ny-2))
                VariGs=sum(abs(gsf(1:nx-2,1:ny-2,k)))/((nx-2)*(ny-2))
                write(66,'(a, f5.1, a, f5.1, a, 2f10.3, f9.4)')'  Z ',depz(k),' - ',depz(k+1),' km  Abs Mean &
                Gc (%)  Gs (%)   dVs (km/s)', VariGc*100, VariGs*100, VariVs
                write(6 ,'(a, f5.1, a, f5.1, a, 2f10.3, f9.4)')'  Z ',depz(k),' - ',depz(k+1),' km  Abs Mean &
                Gc (%)  Gs (%)   dVs (km/s)', VariGc*100, VariGs*100, VariVs
                enddo
            endif
            
            !----------------------------------------------------------------------!
            if (vorotomo == 0) then
                ! Calculate ||Lm||2
                Mnorm2 = 0
                MwNorm2 = 0
                if (iso_inv) then
                    Nre=nar-nar1
                    call Calmodel2Norm(nar1, nar, maxvp, Nre, rw, col, dv, lameVs)
                else
                    Nre=nar-nar1
                    NreVs=narVs-nar1
                    call Calmodel2NormJoint(nar1, nar, maxvp, NreVs, Nre, rw, col, dv, lameGcs, lameVs)
                endif
                !----------------------------------------------------------------------!
                ! Calculate ||W(Gm-d)||2
                ! first forward calculate the traveltime.
                ! Tdata: inversion used data: tobs- tref(iso)
                ! Calculate deltaT=deltaTvs+deltaTaa
                ! Here restT=deltaT_in-deltaT_out
                resbst=0
                fwdTvs=0
                fwdTaa=0
                if (trim(adjustl(ray_backend)) == 'fmm_parallel') then
                    call CalSparseReslNorm(m, n, maxvp, dall, leniw, lenrw, iw, rw, dv(1:maxm), &
                                        datweight, Tdata, fwdTvs, fwdTaa, resbst, iso_inv)
                elseif (iso_inv) then
                    call CalVsReslNorm(maxvp,dall, GVs, dv(1:maxvp), datweight, Tdata, fwdTvs, resbst)
                else
                    VsInv=dv(1:maxvp)
                    GcInv=dv(1+maxvp:2*maxvp)
                    GsInv=dv(1+maxvp*2:3*maxvp)
                    call  CalReslNormJoint(maxvp, dall, GVs, GGc, GGs, VsInv, GcInv, GsInv, datweight, Tdata, &
                    fwdTvs, fwdTaa, resbst)
                    ! update T_AA. useless.
                    RefTaa=0
                    forall(i=1:dall)
                        RefTaa(i)=fwdTaa(i)
                    end forall
                endif
            endif
            !----------------------------------------------------------------------!
            !  Residual Traveltime Analsysi
            !  dsyn: ref-iso traveltime from G-matrix based on vsf model.
            !  Tdata: inversion used data: tobs- tref(iso)
            !  resbst: residual, without weight(res2Norm with weight)
            if ((iter .eq. 1) .or. (iter .eq. maxiter)) then
                write(id,'(I2.2)') iter
                filename = 'Traveltime_statis_'//TRIM(id)//'th.dat'
                open(88,file=TRIM(filename))
                if (iso_mod) then
                    write(88,'(7a)')'   Dist(km)   T_obs(s)  T_ref_iso   Res(in)   dT(dvs)   Res(out)'
                    do i=1,dall
                        write(88,'(3f10.3, 3e12.3)') dist(i), obst(i), dsyn(i), Tdata(i), fwdTvs(i), resbst(i)
                    enddo
                else
                    write(88,'(7a)')'          Dist(km)       T_obs(s)        T_ref-iso        Res(in)   &
                    dT(aa)        dT(dvs)        Res(out)'
                    do i=1,dall
                        write(88,'(3f10.4, 4e12.3)') dist(i), obst(i), dsyn(i), Tdata(i), fwdTvs(i), fwdTaa(i), resbst(i)
                    enddo
                endif
                close(88)
            endif
            !----------------------------------------------------------------------!
            ! Statistic for the  data residual.
            mean = sum(resbst(1:dall))/dall
            meanAbs=sum(abs(resbst(1:dall)))/dall
            std_devs = sqrt(sum((resbst(1:dall)-mean)**2)/dall)

            write(6 ,'(a,f12.4,a,f10.2,a,f10.2,a)')  '  After Inversion: abs mean, std, RMS of Res :', meanAbs,' s ',&
            std_devs,' s ', dnrm2(dall,resbst,1)/sqrt(real(dall)),' s'
            write(66,'(a,f12.4,a,f10.2,a,f10.2,a)')  '  After Inversion: abs mean, std, RMS of Res :', meanAbs,' s ',&
            std_devs,' s ', dnrm2(dall,resbst,1)/sqrt(real(dall)),' s'

            write(67,'(3f12.4)') mean*1000,1000*std_devs,dnrm2(dall,cbst,1)/sqrt(real(dall))
            ! mean=sum(abs(RefTaa(1:dall)))/dall
            ! mean2=sum(abs(fwdTvs(1:dall)))/dall
            ! write(66,'(a,f12.4,a,f12.4)')'  ABS Mean Taa:',mean,'s    ABS Mean T(dv)', mean2
            !----------------------------------------------------------------------!
            if (iso_mod) then
                write(34,*)',OUTPUT S VELOCITY AT ITERATION',iter
                do k=1,nz
                    do j=1,ny
                        write(34,'(100f7.3)') (vsf(i,j,k),i=1,nx)
                    enddo
                enddo
                write(34,*)',OUTPUT DWS AT ITERATION',iter
                do k=1,nz-1
                    do j=2,ny-1
                        write(34,'(100f10.3)') (norm((k-1)*(ny-2)*(nx-2)+(j-2)*(nx-2)+i-1),i=2,nx-1)
                    enddo
                enddo
            endif
            !----------------------------------------------------------------------!
            write(66,'(a)') ' '
            write(6 ,'(a)') '  '

            ! output result of each iteration
            if (output_intermediate /= 0) then ! output all intermediate results?
                if (iter < maxiter) then ! dont output the final iteration results here
                    write(outmodel_iso,'(a,a,i3.3)') 'DSurfTomo.inv','.iter',iter
                    write(outmodel_azi,'(a,a,i3.3)') 'Gc_Gs_model.inv','.iter',iter
                    open(64,file=outmodel_iso)
                    open(74,file=outmodel_azi)
                    if (vorotomo /= 0) then
                        call writeVsmodel_voro(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,vsf,dvstd,dvstdall,64)
                        if (.not. iso_mod) then
                            call writeAzimuthal_voro(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,gcf,gsf,vsf,dvstd,dvstdall,74)
                        else
                            ! isotropic mode: no Gc/Gs uncertainty vectors exist; keep file consistent but without voro std columns
                            call writeAzimuthal_voro_iso(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,vsf,dvstd,dvstdall,74)
                        endif
                    else
                        call writeVsmodel(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,vsf,64)
                        if (.not. iso_mod) then
                            call writeAzimuthal(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,gcf,gsf,vsf,74)
                        else
                            call writeAzimuthal_iso(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,vsf,74)
                        endif
                    endif
                    close(64)
                    close(74)
                endif
            endif

            ! output the raypath only for the last iteration
            timer_section_end = omp_get_wtime()
            timer_post_update = timer_section_end - timer_section_start

            timer_section_start = omp_get_wtime()
            if (iter == maxiter .and. Refwritepath == 1) then
                if (iso_inv) then
                    write(66,*) 'Calculating isotropic raypaths.'
                    write(6,*) 'Calculating isotropic raypaths.'
                    if (trim(adjustl(ray_backend)) == 'fmm_parallel') then
                        call CalSurfGParallelIso(nx, ny, nz, maxvp, vsf, iw, rw, col, dsyn, &
                        goxd, gozd, dvxd, dvzd, kmaxRc, tRc, periods, depz, minthk, &
                        scxf, sczf, rcxf, rczf, nrc1, nsrc1, wavetype, igrt, &
                        kmax, nsrc, nrc, nar, maxnar, Refwritepath, dall, tRcV, fmm_threads)
                    else
                        call CalSurfG(nx,ny,nz,maxvp,vsf,iw,rw,col,dsyn,&
                        GVs, dall,&
                        goxd, gozd, dvxd, dvzd, kmaxRc, tRc, periods, depz, minthk,&
                        scxf, sczf, rcxf, rczf, nrc1, nsrc1, kmax, nsrc, nrc, nar, Refwritepath)
                    endif
                else
                    write(66,*) 'Calculating azimuthal anisotropic raypaths.'
                    write(6,*) 'Calculating azimuthal anisotropic raypaths.'
                    if (trim(adjustl(ray_backend)) == 'fmm_parallel') then
                        call CalSurfGParallelAniso(nx, ny, nz, maxvp, vsf, iw, rw, col, dsyn, &
                        Lsen_Gsc, rmax, tRcV, goxd, gozd, dvxd, dvzd, kmaxRc, tRc, periods, depz, minthk, &
                        scxf, sczf, rcxf, rczf, nrc1, nsrc1, kmax, nsrc, nrc, nar, maxnar, &
                        Refwritepath, dall, fmm_threads)
                    else
                        call CalSurfGAnisoJoint(nx, ny, nz, maxvp, vsf, iw, rw, col, dsyn, &
                        GVs, GGc, GGs, Lsen_Gsc, dall, rmax, tRcV, &
                        goxd, gozd, dvxd, dvzd, kmaxRc, tRc, periods,depz, minthk,&
                        scxf, sczf, rcxf, rczf, nrc1, nsrc1, kmax, nsrc, nrc, nar, Refwritepath)
                    endif
                endif
            endif
            timer_section_end = omp_get_wtime()
            timer_raypath_output = timer_section_end - timer_section_start

            timer_iter_end = omp_get_wtime()
            write(*,'(a,i0,a)') '  Timing summary for iteration ', iter, ':'
            write(*,'(a,f12.3,a)') '    G build incl. kernels/FMM/ray/G rows = ', timer_g_matrix, ' s'
            write(*,'(a,f12.3,a)') '    Residual + weighting + traveltime I/O = ', timer_residual_weight, ' s'
            write(*,'(a,f12.3,a)') '    Solver/projection/ensemble = ', timer_solver, ' s'
            write(*,'(a,f12.3,a)') '    Post-update/residual/output = ', timer_post_update, ' s'
            write(*,'(a,f12.3,a)') '    Final raypath output = ', timer_raypath_output, ' s'
            write(*,'(a,f12.3,a)') '    Iteration total = ', timer_iter_end - timer_iter_start, ' s'
            write(66,'(a,i0,a)') '  Timing summary for iteration ', iter, ':'
            write(66,'(a,f12.3,a)') '    G build incl. kernels/FMM/ray/G rows = ', timer_g_matrix, ' s'
            write(66,'(a,f12.3,a)') '    Residual + weighting + traveltime I/O = ', timer_residual_weight, ' s'
            write(66,'(a,f12.3,a)') '    Solver/projection/ensemble = ', timer_solver, ' s'
            write(66,'(a,f12.3,a)') '    Post-update/residual/output = ', timer_post_update, ' s'
            write(66,'(a,f12.3,a)') '    Final raypath output = ', timer_raypath_output, ' s'
            write(66,'(a,f12.3,a)') '    Iteration total = ', timer_iter_end - timer_iter_start, ' s'

        enddo ! iteration
        
        if (nuq_iter > 0) then
            dv_mean_accum        = dv_mean_accum  / nuq_iter
            dv_mean2_accum       = dv_mean2_accum / nuq_iter
            dv_var_within_accum  = dv_var_within_accum / nuq_iter

            dvstdall = sqrt( max(0.0, dv_var_within_accum + &
                                (dv_mean2_accum - dv_mean_accum**2)) )
        else
            dvstdall = 0.0
        endif

        close(ncost)
        close(67)
        ! if (iso_inv) then
        ! output the updated model, could be used as initial MOD for the next iteration.
        filename = 'MOD_Ref'
        open(11, file=filename)
        do k=1,nz
            write(11,'(f7.1)',advance='no') depz(k)
        enddo
        do k = 1,nz
            do j = 1,ny
                do i=1,nx
                    if (i.eq.1) then
                        write(11,'(/f8.4)',advance='no')vsf(i,j,k)
                    else
                        write(11,'(f8.4)',advance='no')vsf(i,j,k)
                    endif
                enddo
            enddo
        enddo
        close(11)
        ! endif
        !----------------------------------------------------------------------!
        ! output inversion result.
        open(63,file='DSurfTomo.inv')
        open(73,file='Gc_Gs_model.inv')
        
        if (vorotomo /= 0) then
            call writeVsmodel_voro(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,vsf,dvstd,dvstdall,63)
            if (.not. iso_mod) then
                call writeAzimuthal_voro(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,gcf,gsf,vsf,dvstd,dvstdall,73)
            else
                call writeAzimuthal_voro_iso(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,vsf,dvstd,dvstdall,73)
            endif
        else
            call writeVsmodel(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,vsf,63)
            if (.not. iso_mod) then
                call writeAzimuthal(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,gcf,gsf,vsf,73)
            else
                call writeAzimuthal_iso(nx,ny,nz,gozd,goxd,dvzd,dvxd,depz,vsf,73)
            endif
        endif
        close(63)
        close(73)
        !----------------------------------------------------------------------!
        open(77,file='phaseV_FWD.dat')
        call WTPeriodPhaseV(nx,ny,gozd,goxd,dvzd,dvxd,kmaxRc,tRc,tRcV,77)
        !----------------------------------------------------------------------!
        if (.not. iso_mod) then
            write(*,'(a)')'  Begin forward calculate period azimuthal A1, A2.'
            open(42,file='period_Azm_tomo.inv',status='replace',action='write')
            call FwdAzimuthalAniMap(nx,ny,nz,maxvp,&
                goxd,gozd,dvxd,dvzd,kmaxRc,tRc,&
                gcf,gsf,Lsen_Gsc,tRcV)
        else
            write(*,'(a)')'  Skip period azimuthal A1/A2 forward in isotropic mode.'
            open(42,file='period_Azm_tomo.inv',status='replace',action='write')
            close(42)
        endif
        !----------------------------------------------------------------------!
        ! output result model
        ! note: for the input Vs model is grid model (nz point)
        ! BUT, the inversion kernel and result is layer (nz-1 layer) model
        ! for Gc,Gs, the input and result are both layer models, so I set the lower boundary depth as the output
        ! for Vs, input is grid model, while the inversion result is layered model. 
        !  I set the mid depth between the points as the depth index, and the output model uses the average Vs.
        write(66,*)'  -----------------------------------------------------------'
        write(*,*)  '  Program finishes successfully'
        write(66,*) '  Program finishes successfully'

        write(*,*)  '  Output inverted shear velocity model: Vs_model_Syn.rela  Vs_model_Syn.abs'
        write(66,*) '  Output inverted shear velocity model: Vs_model_Syn.rela  Vs_model_Syn.abs'

        endT=OMP_get_wtime()
        write(* ,'(a, f13.1, a)') '   All time cost= ',endT-startT,"s"
        write(66,'(a, f13.1, a)') '   All time cost= ',endT-startT,"s"

        close(nout) ! close lsmr.txt

        close(66)   ! close surf_tomo.log
        close(34)   ! IterVel.out
        if (vorotomo /= 0) then
            close(fid_misfit)
            close(fid_misfit_voro)
        endif
        deallocate(obst)
        deallocate(dsyn)
        deallocate(dist)
        deallocate(pvall)
        deallocate(depz)
        deallocate(norm)
        deallocate(scxf,sczf)
        deallocate(rcxf,rczf)
        deallocate(wavetype,igrt,nrc1)
        deallocate(nsrc1,periods)
        deallocate(rw)
        deallocate(iw,col)
        deallocate(cbst)
        ! deallocate(dv)
        deallocate(dv,dvstd,dvstdall,dvsub,dvall)
        deallocate(dv_mean_accum, dv_mean2_accum, dv_var_within_accum)
        deallocate(vsf)
        deallocate(sigmaT,Tdata,resbst,fwdT,resSigma,resbst_iso)
        deallocate(fwdTvs,fwdTaa)
        deallocate(RefTaa,Tref)
        deallocate(GVs,GGs,GGc,Lsen_Gsc)
        deallocate(nrow)

        deallocate(gcf,gsf)
        if(kmaxRc.gt.0) then
            deallocate(tRc)
        endif
        deallocate(tRcV)
        ! voro misfit
        if (allocated(w_rank)) deallocate(w_rank)
        if (allocated(w_by_real)) deallocate(w_by_real)
        if (allocated(idx_sorted)) deallocate(idx_sorted)
        if (allocated(misfit_voro)) deallocate(misfit_voro)
        if (allocated(acond_voro)) deallocate(acond_voro)
        if (allocated(itn_voro)) deallocate(itn_voro)
        if (allocated(istop_voro)) deallocate(istop_voro)
        if (allocated(misfit_sorted)) deallocate(misfit_sorted)
        
contains
    ! Internal subroutine for triplet diagnostics for vorotomo
    subroutine check_triplets(iw, col, rw, nar, dall, nmodel, fid)
        use, intrinsic :: ieee_arithmetic
        implicit none
        integer, intent(in) :: dall, nmodel, fid
        integer(kind=8), intent(in) :: nar
        integer, intent(in) :: iw(:), col(:)
        real, intent(in) :: rw(:)

        integer(kind=8) :: k
        integer :: r, c
        integer :: bad_row, bad_col, bad_nan
        integer(kind=8) :: breaks
        integer :: last_r, last_c
        logical :: nondecreasing

        integer, allocatable :: cnt(:)
        integer, allocatable :: first_idx(:), last_idx(:)
        integer(kind=8) :: nnz_sum
        integer :: srow(5), si

        bad_row = 0
        bad_col = 0
        bad_nan = 0
        breaks  = 0_8
        nondecreasing = .true.

        last_r = iw(1)
        last_c = col(1)

        allocate(cnt(dall))
        allocate(first_idx(dall))
        allocate(last_idx(dall))
        cnt = 0
        first_idx = 0
        last_idx = 0

        do k = 1_8, nar
            r = iw(k)
            c = col(k)

            if (r < 1 .or. r > dall) bad_row = bad_row + 1
            if (c < 1 .or. c > nmodel) bad_col = bad_col + 1
            if (.not. ieee_is_finite(real(rw(k),kind=8))) bad_nan = bad_nan + 1

            if (k > 1_8) then
                if (r < last_r) then
                    nondecreasing = .false.
                    breaks = breaks + 1_8
                else if (r == last_r) then
                    if (c < last_c) breaks = breaks + 1_8
                end if
                last_r = r
                last_c = c
            end if

            if (r >= 1 .and. r <= dall) then
                cnt(r) = cnt(r) + 1
                if (first_idx(r) == 0) first_idx(r) = int(k)
                last_idx(r) = int(k)
            end if
        end do

        nnz_sum = 0_8
        do r = 1, dall
            nnz_sum = nnz_sum + int(cnt(r), kind=8)
        end do

        write(*,'(a)') '--- [VORO CHECK] Triplet diagnostics ---'
        write(*,'(a,1x,i0)') '  nar =', nar
        write(*,'(a,1x,i0)') '  dall=', dall
        write(*,'(a,1x,i0)') '  nmodel(maxm)=', nmodel
        write(*,'(a,1x,i0)') '  bad_row_count =', bad_row
        write(*,'(a,1x,i0)') '  bad_col_count =', bad_col
        write(*,'(a,1x,i0)') '  nonfinite_rw_count =', bad_nan
        write(*,'(a,1x,l1,1x,a,1x,i0)') '  iw nondecreasing =', nondecreasing, '  breaks=', breaks
        write(*,'(a,1x,i0)') '  sum(cnt) =', nnz_sum

        write(fid,'(a)') '--- [VORO CHECK] Triplet diagnostics ---'
        write(fid,'(a,1x,i0)') '  nar =', nar
        write(fid,'(a,1x,i0)') '  dall=', dall
        write(fid,'(a,1x,i0)') '  nmodel(maxm)=', nmodel
        write(fid,'(a,1x,i0)') '  bad_row_count =', bad_row
        write(fid,'(a,1x,i0)') '  bad_col_count =', bad_col
        write(fid,'(a,1x,i0)') '  nonfinite_rw_count =', bad_nan
        write(fid,'(a,1x,l1,1x,a,1x,i0)') '  iw nondecreasing =', nondecreasing, '  breaks=', breaks
        write(fid,'(a,1x,i0)') '  sum(cnt) =', nnz_sum

        if (bad_row > 0 .or. bad_col > 0 .or. bad_nan > 0) then
            stop 'invalid sparse triplets detected before Voronoi projection'
        end if
        if (.not. nondecreasing) then
            stop 'sparse triplets are not row-contiguous before Voronoi projection'
        end if

        ! sample a few rows to see if entries are contiguous (CSR-like)
        srow(1) = 1
        srow(2) = max(1, dall/4)
        srow(3) = max(1, dall/2)
        srow(4) = max(1, 3*dall/4)
        srow(5) = dall

        write(*,'(a)') '  sample rows: row  first_idx  last_idx  cnt'
        write(fid,'(a)') '  sample rows: row  first_idx  last_idx  cnt'
        do si = 1, 5
            r = srow(si)
            write(*,'(i8,3x,i10,3x,i10,3x,i10)') r, first_idx(r), last_idx(r), cnt(r)
            write(fid,'(i8,3x,i10,3x,i10,3x,i10)') r, first_idx(r), last_idx(r), cnt(r)
        end do

        write(*,'(a)') '--- [VORO CHECK] End ---'
        write(fid,'(a)') '--- [VORO CHECK] End ---'

        deallocate(cnt)
        deallocate(first_idx)
        deallocate(last_idx)
    end subroutine check_triplets

    ! Internal subroutine for index sort (ascending)
    subroutine index_sort(arr, idx)
        implicit none
        real, intent(inout) :: arr(:)
        integer, intent(out) :: idx(:)
        integer :: i, j, n
        real :: temp
        integer :: itemp

        n = size(arr)
        do i = 1, n
            idx(i) = i
        end do
        ! NOTE: O(n^2) selection-style sort is fine here because n = nrealizations
        ! is typically small (order 10^1–10^2). If you ever push nrealizations
        ! to thousands, replace this with an O(n log n) sort.
        do i = 1, n-1
            do j = i+1, n
                if (arr(j) < arr(i)) then
                    temp = arr(i)
                    arr(i) = arr(j)
                    arr(j) = temp
                    itemp = idx(i)
                    idx(i) = idx(j)
                    idx(j) = itemp
                end if
            end do
        end do
    end subroutine index_sort

    ! Median function for sorted array
    real function median(arr, n)
        real, intent(in) :: arr(:)
        integer, intent(in) :: n
        integer :: mid
        if (n <= 0) then
            median = 0.0
            return
        end if
        if (mod(n,2) == 1) then
            mid = (n+1)/2
            median = arr(mid)
        else
            mid = n/2
            median = 0.5*(arr(mid) + arr(mid+1))
        end if
    end function median

end program SurfAniso



subroutine writeVsmodel(nx, ny, nz, gozd, goxd, dvzd, dvxd, depz, vsf, Idout)
    integer nx,ny,nz
    real goxd,gozd
    real dvxd,dvzd
    real depz(nz)
    ! real vsAbs(nx-2,ny-2,nz-1), vsRela(nx-2,ny-2,nz-1)
    integer Idout
    real:: vsf(nx,ny,nz)
    integer:: k,j,i
    real :: vsref
    ! Because the (gozd,goxd) point is NE point of the inversion region.---> NE move 1 point.
    do k=1,nz
          do j=1,ny
                do i=1,nx
                    write(Idout,'(4f10.4)') gozd+(j-2)*dvzd,goxd-(i-2)*dvxd,depz(k),vsf(i,j,k)
                enddo
            enddo
        enddo
end subroutine writeVsmodel

subroutine writeVsmodel_voro(nx, ny, nz, gozd, goxd, dvzd, dvxd, depz, vsf, dvstd, dvstdall, Idout)
    implicit none
    integer, intent(in) :: nx, ny, nz, Idout
    real,    intent(in) :: goxd, gozd, dvxd, dvzd
    real,    intent(in) :: depz(nz)
    real,    intent(in) :: vsf(nx,ny,nz)
    real,    intent(in) :: dvstd((nx-2)*(ny-2)*(nz-1))
    real,    intent(in) :: dvstdall((nx-2)*(ny-2)*(nz-1))

    integer :: i, j, k
    integer :: maxvp, idx_cell
    real    :: std_out, stdall_out

    maxvp = (nx-2)*(ny-2)*(nz-1)

    do k = 1, nz
        do j = 1, ny
            do i = 1, nx
                std_out    = 0.0
                stdall_out = 0.0

                ! only valid for inversion cells: i=2..nx-1, j=2..ny-1, k=1..nz-1
                if (k <= nz-1 .and. i>=2 .and. i<=nx-1 .and. j>=2 .and. j<=ny-1) then
                    idx_cell = (k-1)*(nx-2)*(ny-2) + (j-2)*(nx-2) + (i-1)   ! 1..maxvp
                    std_out    = dvstd(idx_cell)
                    stdall_out = dvstdall(idx_cell)
                end if

                write(Idout,'(6f10.4)') gozd+(j-2)*dvzd, goxd-(i-2)*dvxd, depz(k), vsf(i,j,k), std_out, stdall_out
            end do
        end do
    end do
end subroutine writeVsmodel_voro

subroutine writeAzimuthal(nx, ny, nz, gozd, goxd, dvzd, dvxd, depz, Gc, Gs, vsf, Idout)
    implicit none
    integer, intent(in) :: nx, ny, nz, Idout
    real,    intent(in) :: goxd, gozd, dvxd, dvzd
    real,    intent(in) :: depz(nz)
    real,    intent(in) :: Gc(nx-2,ny-2,nz-1), Gs(nx-2,ny-2,nz-1)
    real,    intent(in) :: vsf(nx,ny,nz)

    integer :: i, j, k
    real    :: cosTmp, sinTmp
    real    :: vsref
    real    :: AzimAmp, AzimAng
    real, parameter :: pi = 3.14159265358979323846

    do k = 1, nz-1
        do j = 1, ny-2
            do i = 1, nx-2
                cosTmp  = Gc(i,j,k)
                sinTmp  = Gs(i,j,k)
                AzimAmp = 0.5*sqrt(cosTmp*cosTmp + sinTmp*sinTmp)
                AzimAng = atan2(sinTmp, cosTmp) / pi * 180.0
                if (AzimAng < 0.0) AzimAng = AzimAng + 360.0
                AzimAng = 0.5 * AzimAng

                vsref = 0.5*(vsf(i+1,j+1,k) + vsf(i+1,j+1,k+1))

                write(Idout,'(8f10.4)') gozd+(j-1)*dvzd, goxd-(i-1)*dvxd, depz(k+1), vsref, &
                                    AzimAng, AzimAmp, Gc(i,j,k)*100.0, Gs(i,j,k)*100.0
            end do
        end do
    end do
end subroutine writeAzimuthal

subroutine writeAzimuthal_iso(nx, ny, nz, gozd, goxd, dvzd, dvxd, depz, vsf, Idout)
    implicit none
    integer, intent(in) :: nx, ny, nz, Idout
    real,    intent(in) :: goxd, gozd, dvxd, dvzd
    real,    intent(in) :: depz(nz)
    real,    intent(in) :: vsf(nx,ny,nz)

    integer :: i, j, k
    real    :: vsref

    do k = 1, nz-1
        do j = 1, ny-2
            do i = 1, nx-2
                vsref = 0.5*(vsf(i+1,j+1,k) + vsf(i+1,j+1,k+1))
                write(Idout,'(8f10.4)') gozd+(j-1)*dvzd, goxd-(i-1)*dvxd, depz(k+1), vsref, &
                    0.0, 0.0, 0.0, 0.0
            end do
        end do
    end do
end subroutine writeAzimuthal_iso

subroutine writeAzimuthal_voro(nx, ny, nz, gozd, goxd, dvzd, dvxd, depz, Gc, Gs, vsf, dvstd, dvstdall, Idout)
    implicit none
    integer, intent(in) :: nx, ny, nz, Idout
    real,    intent(in) :: goxd, gozd, dvxd, dvzd
    real,    intent(in) :: depz(nz)
    real,    intent(in) :: Gc(nx-2,ny-2,nz-1), Gs(nx-2,ny-2,nz-1)
    real,    intent(in) :: vsf(nx,ny,nz)
    real,    intent(in) :: dvstd((nx-2)*(ny-2)*(nz-1)*3)
    real,    intent(in) :: dvstdall((nx-2)*(ny-2)*(nz-1)*3)

    integer :: i, j, k, idx
    integer :: maxvp
    real    :: cosTmp, sinTmp, vsref
    real    :: AzimAmp, AzimAng
    real, parameter :: pi = 3.14159265358979323846

    maxvp = (nx-2)*(ny-2)*(nz-1)

    do k = 1, nz-1
        do j = 1, ny-2
            do i = 1, nx-2
                cosTmp  = Gc(i,j,k)
                sinTmp  = Gs(i,j,k)
                AzimAmp = 0.5*sqrt(cosTmp*cosTmp + sinTmp*sinTmp)
                AzimAng = atan2(sinTmp, cosTmp) / pi * 180.0
                if (AzimAng < 0.0) AzimAng = AzimAng + 360.0
                AzimAng = 0.5 * AzimAng

                vsref = 0.5*(vsf(i+1,j+1,k) + vsf(i+1,j+1,k+1))
                idx   = (k-1)*(nx-2)*(ny-2) + (j-1)*(nx-2) + i   ! 1..maxvp

                write(Idout,'(14f10.4)') gozd+(j-1)*dvzd, goxd-(i-1)*dvxd, depz(k+1), vsref, &
                    AzimAng, AzimAmp, Gc(i,j,k)*100.0, Gs(i,j,k)*100.0, &
                    dvstd(idx), dvstd(idx+maxvp)*100.0, dvstd(idx+2*maxvp)*100.0, &
                    dvstdall(idx), dvstdall(idx+maxvp)*100.0, dvstdall(idx+2*maxvp)*100.0
            end do
        end do
    end do
end subroutine writeAzimuthal_voro

subroutine writeAzimuthal_voro_iso(nx, ny, nz, gozd, goxd, dvzd, dvxd, depz, vsf, dvstd, dvstdall, Idout)
    implicit none
    integer, intent(in) :: nx, ny, nz, Idout
    real,    intent(in) :: goxd, gozd, dvxd, dvzd
    real,    intent(in) :: depz(nz)
    real,    intent(in) :: vsf(nx,ny,nz)
    real,    intent(in) :: dvstd((nx-2)*(ny-2)*(nz-1))
    real,    intent(in) :: dvstdall((nx-2)*(ny-2)*(nz-1))

    integer :: i, j, k, idx
    integer :: maxvp
    real    :: vsref
    real    :: AzimAmp, AzimAng
    real    :: Gc_pct, Gs_pct
    real    :: std_vs, std_gc, std_gs
    real    :: stdall_vs, stdall_gc, stdall_gs

    maxvp = (nx-2)*(ny-2)*(nz-1)

    ! Isotropic mode: no meaningful azimuth (angle/amplitude) and no Gc/Gs.
    ! Keep 14-column format consistent with writeAzimuthal_voro:
    ! lon lat depth vsref AzimAng AzimAmp Gc% Gs% std(Vs) std(Gc%) std(Gs%) stdall(Vs) stdall(Gc%) stdall(Gs%)

    do k = 1, nz-1
        do j = 1, ny-2
            do i = 1, nx-2
                vsref = 0.5*(vsf(i+1,j+1,k) + vsf(i+1,j+1,k+1))
                idx   = (k-1)*(nx-2)*(ny-2) + (j-1)*(nx-2) + i   ! 1..maxvp

                AzimAng = 0.0
                AzimAmp = 0.0
                Gc_pct  = 0.0
                Gs_pct  = 0.0

                std_vs     = dvstd(idx)
                std_gc     = 0.0
                std_gs     = 0.0
                stdall_vs  = dvstdall(idx)
                stdall_gc  = 0.0
                stdall_gs  = 0.0

                write(Idout,'(14f10.4)') gozd+(j-1)*dvzd, goxd-(i-1)*dvxd, depz(k+1), vsref, &
                    AzimAng, AzimAmp, Gc_pct, Gs_pct, &
                    std_vs, std_gc, std_gs, &
                    stdall_vs, stdall_gc, stdall_gs
            end do
        end do
    end do
end subroutine writeAzimuthal_voro_iso

subroutine WTPeriodPhaseV(nx,ny,gozd,goxd,dvzd,dvxd,kmaxRc,tRc,pvRc,Idout)
        implicit none
        integer nx,ny
        real goxd,gozd
        real dvxd,dvzd
        integer kmaxRc
        real*8 tRc(kmaxRc)
        real*8 pvRc((nx-2)*(ny-2),kmaxRc)
        integer Idout
        integer:: tt,k,jj,ii
        ! note:  pvRc((nx-2)*(ny-2),kmaxRc) ----different pvRc in
        do tt=1,kmaxRc
            do jj=1,ny-2
                do ii=1,nx-2
                    write(Idout,'(5f10.4)') gozd+(jj-1)*dvzd,goxd-(ii-1)*dvxd,tRc(tt),pvRc((jj-1)*(nx-2)+ii,tt)
                enddo
            enddo
        enddo
        close(Idout)
end subroutine WTPeriodPhaseV

subroutine CalRmax(nz, depz, minthk0, rmax)
        implicit none
        integer nz
        real depz(nz)
        real minthk0, minthk
        integer rmax
        integer k,i
        integer NL
        parameter (NL=200)
        integer nsublay(NL)
        real thk
        rmax=0
        do i=1, nz-1
            thk = depz(i+1)-depz(i)
            minthk = thk/minthk0
            nsublay(i) = int((thk+1.0e-4)/minthk) + 1
        enddo
        rmax=sum(nsublay(1:nz-1))
end subroutine CalRmax
