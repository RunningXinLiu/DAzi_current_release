subroutine depthkernel_column_phase_fd(nz, vs, kmax, periods, iwave, igr, depz, minthk, pv)
    implicit none
    integer, parameter :: NL = 200, NP = 160
    integer, intent(in) :: nz, kmax, iwave, igr
    real, intent(in) :: vs(nz), depz(nz), minthk
    real*8, intent(in) :: periods(*)
    real*8, intent(out) :: pv(*)

    integer :: rmax, iflsph, mode, nk_surf
    real :: vp(nz), rho(nz)
    real :: rdep(NL), rvp(NL), rvs(NL), rrho(NL), rthk(NL)
    real*8 :: cg(NP)

    if (kmax > NP .or. nz > NL) stop 'depthkernel_column_phase_fd limits exceeded.'

    iflsph = 1
    mode = 1
    call depthkernel_empirical_vp_rho(nz, vs, vp, rho)
    call refineGrid2LayerMdl(minthk, nz, depz, vp, vs, rho, rmax, rdep, rvp, rvs, rrho, rthk)
    nk_surf = kmax
    call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, nk_surf, periods, cg)
    pv(1:kmax) = cg(1:kmax)
end subroutine depthkernel_column_phase_fd

subroutine depthkernel_column_depth_fd(nz, vs, idep, kmax, periods, iwave, igr, depz, minthk, &
    dlnvs, dlnvp, dlnrho, sen_vs, sen_vp, sen_rho)
    implicit none
    integer, parameter :: NL = 200, NP = 160
    integer, intent(in) :: nz, idep, kmax, iwave, igr
    real, intent(in) :: vs(nz), depz(nz), minthk
    real, intent(in) :: dlnvs, dlnvp, dlnrho
    real*8, intent(in) :: periods(*)
    real*8, intent(out) :: sen_vs(*), sen_vp(*), sen_rho(*)

    integer :: i, nn, rmax, iflsph, mode, nk_surf
    real :: vp(nz), rho(nz)
    real :: depm(NL), vpm(NL), vsm(NL), rhom(NL)
    real :: rdep(NL), rvp(NL), rvs(NL), rrho(NL), rthk(NL)
    real*8 :: cg1(NP), cg2(NP)

    if (kmax > NP .or. nz > NL) stop 'depthkernel_column_depth_fd limits exceeded.'
    if (idep < 1 .or. idep > nz) stop 'depthkernel_column_depth_fd invalid depth index.'

    iflsph = 1
    mode = 1
    sen_vs(1:kmax) = 0.0d0
    sen_vp(1:kmax) = 0.0d0
    sen_rho(1:kmax) = 0.0d0
    call depthkernel_empirical_vp_rho(nz, vs, vp, rho)

    do i = 1, nz-1
        depm(i) = depz(i)
        vsm(i) = vs(i)
        vpm(i) = vp(i)
        rhom(i) = rho(i)
    enddo
    depm(nz) = depz(nz)
    vsm(nz) = vs(nz)
    vpm(nz) = vp(nz)
    rhom(nz) = rho(nz)

    vsm(idep) = vs(idep) - 0.5*dlnvs*vs(idep)
    call refineGrid2LayerMdl(minthk, nz, depm, vpm, vsm, rhom, rmax, rdep, rvp, rvs, rrho, rthk)
    nk_surf = kmax
    call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, nk_surf, periods, cg1)
    vsm(idep) = vs(idep) + 0.5*dlnvs*vs(idep)
    call refineGrid2LayerMdl(minthk, nz, depm, vpm, vsm, rhom, rmax, rdep, rvp, rvs, rrho, rthk)
    nk_surf = kmax
    call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, nk_surf, periods, cg2)
    vsm(idep) = vs(idep)
    do nn = 1, kmax
        sen_vs(nn) = (cg2(nn) - cg1(nn))/(dlnvs*vs(idep))
    enddo

    vpm(idep) = vp(idep) - 0.5*dlnvp*vp(idep)
    call refineGrid2LayerMdl(minthk, nz, depm, vpm, vsm, rhom, rmax, rdep, rvp, rvs, rrho, rthk)
    nk_surf = kmax
    call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, nk_surf, periods, cg1)
    vpm(idep) = vp(idep) + 0.5*dlnvp*vp(idep)
    call refineGrid2LayerMdl(minthk, nz, depm, vpm, vsm, rhom, rmax, rdep, rvp, rvs, rrho, rthk)
    nk_surf = kmax
    call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, nk_surf, periods, cg2)
    vpm(idep) = vp(idep)
    do nn = 1, kmax
        sen_vp(nn) = (cg2(nn) - cg1(nn))/(dlnvp*vp(idep))
    enddo

    rhom(idep) = rho(idep) - 0.5*dlnrho*rho(idep)
    call refineGrid2LayerMdl(minthk, nz, depm, vpm, vsm, rhom, rmax, rdep, rvp, rvs, rrho, rthk)
    nk_surf = kmax
    call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, nk_surf, periods, cg1)
    rhom(idep) = rho(idep) + 0.5*dlnrho*rho(idep)
    call refineGrid2LayerMdl(minthk, nz, depm, vpm, vsm, rhom, rmax, rdep, rvp, rvs, rrho, rthk)
    nk_surf = kmax
    call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, nk_surf, periods, cg2)
    rhom(idep) = rho(idep)
    do nn = 1, kmax
        sen_rho(nn) = (cg2(nn) - cg1(nn))/(dlnrho*rho(idep))
    enddo
end subroutine depthkernel_column_depth_fd

subroutine depthkernel_column_fd(nz, vs, kmax, periods, iwave, igr, depz, minthk, &
    dlnvs, dlnvp, dlnrho, pv, sen_vs, sen_vp, sen_rho)
    implicit none
    integer, parameter :: NL = 200, NP = 160
    integer, intent(in) :: nz, kmax, iwave, igr
    real, intent(in) :: vs(nz), depz(nz), minthk
    real, intent(in) :: dlnvs, dlnvp, dlnrho
    real*8, intent(in) :: periods(*)
    real*8, intent(out) :: pv(*), sen_vs(kmax,*), sen_vp(kmax,*), sen_rho(kmax,*)

    integer :: i, nn, rmax, iflsph, mode, nk_surf
    real :: vp(nz), rho(nz)
    real :: depm(NL), vpm(NL), vsm(NL), rhom(NL)
    real :: rdep(NL), rvp(NL), rvs(NL), rrho(NL), rthk(NL)
    real*8 :: cg0(NP), cg1(NP), cg2(NP)

    if (kmax > NP .or. nz > NL) stop 'depthkernel_column_fd limits exceeded.'

    iflsph = 1
    mode = 1

    call depthkernel_empirical_vp_rho(nz, vs, vp, rho)

    call refineGrid2LayerMdl(minthk, nz, depz, vp, vs, rho, rmax, rdep, rvp, rvs, rrho, rthk)
    nk_surf = kmax
    call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, nk_surf, periods, cg0)
    pv(1:kmax) = cg0(1:kmax)

    do i = 1, nz-1
        depm(i) = depz(i)
        vsm(i) = vs(i)
        vpm(i) = vp(i)
        rhom(i) = rho(i)
    enddo
    depm(nz) = depz(nz)
    vsm(nz) = vs(nz)
    vpm(nz) = vp(nz)
    rhom(nz) = rho(nz)

    do i = 1, nz
        sen_vs(1:kmax,i) = 0.0d0
        sen_vp(1:kmax,i) = 0.0d0
        sen_rho(1:kmax,i) = 0.0d0

        vsm(i) = vs(i) - 0.5*dlnvs*vs(i)
        call refineGrid2LayerMdl(minthk, nz, depm, vpm, vsm, rhom, rmax, rdep, rvp, rvs, rrho, rthk)
        nk_surf = kmax
        call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, nk_surf, periods, cg1)
        vsm(i) = vs(i) + 0.5*dlnvs*vs(i)
        call refineGrid2LayerMdl(minthk, nz, depm, vpm, vsm, rhom, rmax, rdep, rvp, rvs, rrho, rthk)
        nk_surf = kmax
        call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, nk_surf, periods, cg2)
        vsm(i) = vs(i)
        do nn = 1, kmax
            sen_vs(nn,i) = (cg2(nn) - cg1(nn))/(dlnvs*vs(i))
        enddo

        vpm(i) = vp(i) - 0.5*dlnvp*vp(i)
        call refineGrid2LayerMdl(minthk, nz, depm, vpm, vsm, rhom, rmax, rdep, rvp, rvs, rrho, rthk)
        nk_surf = kmax
        call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, nk_surf, periods, cg1)
        vpm(i) = vp(i) + 0.5*dlnvp*vp(i)
        call refineGrid2LayerMdl(minthk, nz, depm, vpm, vsm, rhom, rmax, rdep, rvp, rvs, rrho, rthk)
        nk_surf = kmax
        call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, nk_surf, periods, cg2)
        vpm(i) = vp(i)
        do nn = 1, kmax
            sen_vp(nn,i) = (cg2(nn) - cg1(nn))/(dlnvp*vp(i))
        enddo

        rhom(i) = rho(i) - 0.5*dlnrho*rho(i)
        call refineGrid2LayerMdl(minthk, nz, depm, vpm, vsm, rhom, rmax, rdep, rvp, rvs, rrho, rthk)
        nk_surf = kmax
        call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, nk_surf, periods, cg1)
        rhom(i) = rho(i) + 0.5*dlnrho*rho(i)
        call refineGrid2LayerMdl(minthk, nz, depm, vpm, vsm, rhom, rmax, rdep, rvp, rvs, rrho, rthk)
        nk_surf = kmax
        call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, nk_surf, periods, cg2)
        rhom(i) = rho(i)
        do nn = 1, kmax
            sen_rho(nn,i) = (cg2(nn) - cg1(nn))/(dlnrho*rho(i))
        enddo
    enddo
end subroutine depthkernel_column_fd

subroutine depthkernel_empirical_vp_rho(nz, vs, vp, rho)
    implicit none
    integer, intent(in) :: nz
    real, intent(in) :: vs(nz)
    real, intent(out) :: vp(nz), rho(nz)
    integer :: k

    do k = 1, nz
        vp(k) = 0.9409 + 2.0947*vs(k) - 0.8206*vs(k)**2 + &
            0.2683*vs(k)**3 - 0.0251*vs(k)**4
        rho(k) = 1.6612*vp(k) - 0.4721*vp(k)**2 + &
            0.0671*vp(k)**3 - 0.0043*vp(k)**4 + 0.000106*vp(k)**5
    enddo
end subroutine depthkernel_empirical_vp_rho
