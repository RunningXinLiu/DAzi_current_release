! ADDING REGULARIZATION TERM
subroutine TikhonovRegularization(nx,ny,nz,maxvp,dall,nar,rw,iw,col,count3,iso_inv,weightGcs,weightVs,maxnar)
    use sparse_index_utils, only: check_sparse_capacity, require_default_int_range
    implicit none
    integer nx,ny,nz,maxvp
    integer dall
    integer(kind=8) nar,maxnar
    real rw(*)

    integer col(*),iw(*)
    integer count3
    real weightGcs, weightVs
    real weight
    logical iso_inv
    integer i,j,k
    integer nvx,nvz,sc
    !-------------------------------------------------------!
    nvz=ny-2
    nvx=nx-2
    count3=0
    if (iso_inv) then
        call require_default_int_range('TikhonovRegularization maximum row index', &
            int(dall,kind=8)+int(maxvp,kind=8))
        call require_default_int_range('TikhonovRegularization maximum column index', int(maxvp,kind=8))
    else
        call require_default_int_range('TikhonovRegularization maximum row index', &
            int(dall,kind=8)+2_8*int(maxvp,kind=8))
        call require_default_int_range('TikhonovRegularization maximum column index', 2_8*int(maxvp,kind=8))
    endif
    !----------------------------------------------------------!
    if (iso_inv) then
        weight=weightVs
        do k=1,nz-1
        do j=1,nvz
        do i=1,nvx
        if(i==1.or.i==nvx.or.j==1.or.j==nvz.or.k==1.or.k==nz-1)then
        call check_sparse_capacity(nar, 1_8, maxnar, 'TikhonovRegularization iso boundary')
        count3=count3+1
        col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i
        rw(nar+1)=2.0*weight
        iw(nar+1)=dall+count3

        nar=nar+1
        else
        call check_sparse_capacity(nar, 7_8, maxnar, 'TikhonovRegularization iso interior')
        count3=count3+1
        col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i
        rw(nar+1)=6.0*weight
        iw(nar+1)=dall+count3
        rw(nar+2)=-1.0*weight
        iw(nar+2)=dall+count3
        col(nar+2)=(k-1)*nvz*nvx+(j-1)*nvx+i-1
        rw(nar+3)=-1.0*weight
        iw(nar+3)=dall+count3
        col(nar+3)=(k-1)*nvz*nvx+(j-1)*nvx+i+1
        rw(nar+4)=-1.0*weight
        iw(nar+4)=dall+count3
        col(nar+4)=(k-1)*nvz*nvx+(j-2)*nvx+i
        rw(nar+5)=-1.0*weight
        iw(nar+5)=dall+count3
        col(nar+5)=(k-1)*nvz*nvx+j*nvx+i
        rw(nar+6)=-1.0*weight
        iw(nar+6)=dall+count3
        col(nar+6)=(k-2)*nvz*nvx+(j-1)*nvx+i
        rw(nar+7)=-1.0*weight
        iw(nar+7)=dall+count3
        col(nar+7)=k*nvz*nvx+(j-1)*nvx+i

        nar=nar+7
        endif
        enddo
        enddo
        enddo
    else
        weight=weightGcs
        do sc=1,2
        do k=1,nz-1
        do j=1,nvz
        do i=1,nvx
        if(i==1.or.i==nvx.or.j==1.or.j==nvz.or.k==1.or.k==nz-1)then
        call check_sparse_capacity(nar, 1_8, maxnar, 'TikhonovRegularization aniso boundary')
        count3=count3+1
        col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i+(sc-1)*maxvp
        rw(nar+1)=2.0*weight
        iw(nar+1)=dall+count3

        nar=nar+1
        else
        call check_sparse_capacity(nar, 7_8, maxnar, 'TikhonovRegularization aniso interior')
        count3=count3+1
        col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i+(sc-1)*maxvp
        rw(nar+1)=6.0*weight
        iw(nar+1)=dall+count3
        rw(nar+2)=-1.0*weight
        iw(nar+2)=dall+count3
        col(nar+2)=(k-1)*nvz*nvx+(j-1)*nvx+i-1+(sc-1)*maxvp
        rw(nar+3)=-1.0*weight
        iw(nar+3)=dall+count3
        col(nar+3)=(k-1)*nvz*nvx+(j-1)*nvx+i+1+(sc-1)*maxvp
        rw(nar+4)=-1.0*weight
        iw(nar+4)=dall+count3
        col(nar+4)=(k-1)*nvz*nvx+(j-2)*nvx+i+(sc-1)*maxvp
        rw(nar+5)=-1.0*weight
        iw(nar+5)=dall+count3
        col(nar+5)=(k-1)*nvz*nvx+j*nvx+i+(sc-1)*maxvp
        rw(nar+6)=-1.0*weight
        iw(nar+6)=dall+count3
        col(nar+6)=(k-2)*nvz*nvx+(j-1)*nvx+i+(sc-1)*maxvp
        rw(nar+7)=-1.0*weight
        iw(nar+7)=dall+count3
        col(nar+7)=k*nvz*nvx+(j-1)*nvx+i+(sc-1)*maxvp

        nar=nar+7
        endif
        enddo
        enddo
        enddo
        enddo
    endif
end subroutine


subroutine TikhRegul_joint(nx,ny,nz,maxvp,dall,nar,rw,iw,col,narVs,count3,weightGcs,weightVs,maxnar)
    use sparse_index_utils, only: check_sparse_capacity, require_default_int_range
    ! list
    ! count3: increase of rows, all rows= dall+count3
    ! nar   : all nonzero number of G
    ! narVs: nonzero number of G after Vs regularization
    implicit none
    integer nx,ny,nz,maxvp
    integer dall
    integer(kind=8) nar,maxnar
    real rw(*)
    integer col(*),iw(*)
    integer count3,narVs
    real weightGcs, weightVs
    real weight
    integer i,j,k
    integer nvx,nvz,sc
    !-------------------------------------------------------!
    nvz=ny-2
    nvx=nx-2
    count3=0
    call require_default_int_range('TikhRegul_joint maximum row index', &
        int(dall,kind=8)+3_8*int(maxvp,kind=8))
    call require_default_int_range('TikhRegul_joint maximum column index', 3_8*int(maxvp,kind=8))
    !----------------------------------------------------------!

    weight=weightVs
    do k=1,nz-1
    do j=1,nvz
    do i=1,nvx
    if(i==1.or.i==nvx.or.j==1.or.j==nvz.or.k==1.or.k==nz-1)then
        call check_sparse_capacity(nar, 1_8, maxnar, 'TikhRegul_joint Vs boundary')
        count3=count3+1
        col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i
        rw(nar+1)=2.0*weight
        iw(nar+1)=dall+count3
        nar=nar+1
    else
        call check_sparse_capacity(nar, 7_8, maxnar, 'TikhRegul_joint Vs interior')
        count3=count3+1
        col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i
        rw(nar+1)=6.0*weight
        iw(nar+1)=dall+count3
        rw(nar+2)=-1.0*weight
        iw(nar+2)=dall+count3
        col(nar+2)=(k-1)*nvz*nvx+(j-1)*nvx+i-1
        rw(nar+3)=-1.0*weight
        iw(nar+3)=dall+count3
        col(nar+3)=(k-1)*nvz*nvx+(j-1)*nvx+i+1
        rw(nar+4)=-1.0*weight
        iw(nar+4)=dall+count3
        col(nar+4)=(k-1)*nvz*nvx+(j-2)*nvx+i
        rw(nar+5)=-1.0*weight
        iw(nar+5)=dall+count3
        col(nar+5)=(k-1)*nvz*nvx+j*nvx+i
        rw(nar+6)=-1.0*weight
        iw(nar+6)=dall+count3
        col(nar+6)=(k-2)*nvz*nvx+(j-1)*nvx+i
        rw(nar+7)=-1.0*weight
        iw(nar+7)=dall+count3
        col(nar+7)=k*nvz*nvx+(j-1)*nvx+i
        nar=nar+7
    endif
    enddo
    enddo
    enddo
    narVs=nar

    weight=weightGcs
    do sc=1,2
    do k=1,nz-1
    do j=1,nvz
    do i=1,nvx
    if(i==1.or.i==nvx.or.j==1.or.j==nvz.or.k==1.or.k==nz-1)then
        call check_sparse_capacity(nar, 1_8, maxnar, 'TikhRegul_joint Gcs boundary')
        count3=count3+1
        col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i+(sc)*maxvp
        rw(nar+1)=2.0*weight
        iw(nar+1)=dall+count3
        nar=nar+1
    else
        call check_sparse_capacity(nar, 7_8, maxnar, 'TikhRegul_joint Gcs interior')
        count3=count3+1
        col(nar+1)=(k-1)*nvz*nvx+(j-1)*nvx+i+(sc)*maxvp
        rw(nar+1)=6.0*weight
        iw(nar+1)=dall+count3
        rw(nar+2)=-1.0*weight
        iw(nar+2)=dall+count3
        col(nar+2)=(k-1)*nvz*nvx+(j-1)*nvx+i-1+(sc)*maxvp
        rw(nar+3)=-1.0*weight
        iw(nar+3)=dall+count3
        col(nar+3)=(k-1)*nvz*nvx+(j-1)*nvx+i+1+(sc)*maxvp
        rw(nar+4)=-1.0*weight
        iw(nar+4)=dall+count3
        col(nar+4)=(k-1)*nvz*nvx+(j-2)*nvx+i+(sc)*maxvp
        rw(nar+5)=-1.0*weight
        iw(nar+5)=dall+count3
        col(nar+5)=(k-1)*nvz*nvx+j*nvx+i+(sc)*maxvp
        rw(nar+6)=-1.0*weight
        iw(nar+6)=dall+count3
        col(nar+6)=(k-2)*nvz*nvx+(j-1)*nvx+i+(sc)*maxvp
        rw(nar+7)=-1.0*weight
        iw(nar+7)=dall+count3
        col(nar+7)=k*nvz*nvx+(j-1)*nvx+i+(sc)*maxvp
        nar=nar+7
    endif
    enddo
    enddo
    enddo
    enddo

end subroutine
