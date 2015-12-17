module mbd_helper

implicit none

contains


function is_in(c, str) result(is)
    character(len=1), intent(in) :: c
    character(len=*), intent(in) :: str
    logical :: is

    integer :: i

    is = .false.
    do i = 1, len(str)
        if (c == str(i:i)) then
            is = .true.
            exit
        end if
    end do
end function is_in


function blanked(cs, str) result(str_out)
    character(len=*), intent(in) :: cs
    character(len=*), intent(in) :: str
    character(len=len(str)) :: str_out

    integer :: i

    str_out = str
    do i = 1, len(str)
        if (is_in(str(i:i), cs)) then
            str_out(i:i) = ' '
        end if
    end do
end function blanked


end module mbd_helper

