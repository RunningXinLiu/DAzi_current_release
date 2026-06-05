subroutine delsph(flat1,flon1,flat2,flon2,del)
use coord_utils, only: great_circle_distance_colat_rad
implicit none
real flat1,flat2
real flon1,flon2
real del

call great_circle_distance_colat_rad(flat1, flon1, flat2, flon2, del)
end subroutine
