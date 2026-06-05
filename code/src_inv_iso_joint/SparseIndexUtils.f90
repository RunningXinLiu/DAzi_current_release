module sparse_index_utils
    implicit none

contains

    subroutine require_default_int_range(name, value)
        character(len=*), intent(in) :: name
        integer(kind=8), intent(in) :: value

        if (value >= 0_8 .and. value <= int(huge(0), kind=8)) return

        write(6,*) 'ERROR: ', trim(name), ' exceeds default INTEGER range.'
        write(6,*) '       value=', value, ' max=', int(huge(0), kind=8)
        stop 'default INTEGER range exceeded'
    end subroutine require_default_int_range

    subroutine require_positive_default_int_range(name, value)
        character(len=*), intent(in) :: name
        integer(kind=8), intent(in) :: value

        if (value >= 1_8 .and. value <= int(huge(0), kind=8)) return

        write(6,*) 'ERROR: ', trim(name), ' is outside positive default INTEGER range.'
        write(6,*) '       value=', value, ' max=', int(huge(0), kind=8)
        stop 'positive default INTEGER range exceeded'
    end subroutine require_positive_default_int_range

    subroutine check_sparse_capacity(nar, needed, maxnar, context)
        integer(kind=8), intent(in) :: nar, needed, maxnar
        character(len=*), intent(in) :: context

        if (needed >= 0_8 .and. nar <= maxnar - needed) return

        write(6,*) 'ERROR: sparse G capacity exceeded in ', trim(context)
        write(6,*) '       current nar=', nar, ' needed additional nnz=', needed, ' maxnar=', maxnar
        write(6,*) '       increase sparse fraction (spfra)'
        stop 'sparse capacity exceeded'
    end subroutine check_sparse_capacity

end module sparse_index_utils
