program SurfAniso
        implicit none
        integer nd,ncells
        integer(kind=8) maxnar,nd1,ncells1
        integer,dimension(:),allocatable:: rwgp,colgp
        integer(kind=8),dimension(:),allocatable:: iwgp
        
        nd = 220649344
        ncells = 100

        nd1 = nd
        ncells1 = ncells

        maxnar = 0.1*long(nd)*ncells
        ! maxnar = 0.1*nd1*ncells1
        write(*,*) 'maxnar is ',maxnar
        allocate(iwgp(maxnar*2+1),colgp(maxnar),rwgp(maxnar))
        write(*,*) 'I am here 0'

        deallocate(iwgp,colgp,rwgp)
end program