module coord_utils
    implicit none
    integer, parameter :: coord_mode_spherical_angular = 1
    integer, parameter :: coord_mode_spherical_elevation = 2
    integer, parameter :: coord_mode_cartesian_xy = 3
    integer, parameter :: coord_mode_cartesian_xyz = 4
    integer, parameter :: coord_mode_default = coord_mode_spherical_angular
    real, parameter :: coord_pi = 3.1415926535898
    real, parameter :: earth_radius_km = 6371.0

contains

    pure function coord_mode_name(coord_mode) result(name)
        integer, intent(in) :: coord_mode
        character(len=32) :: name

        select case (coord_mode)
        case (coord_mode_spherical_angular)
            name = 'spherical_angular'
        case (coord_mode_spherical_elevation)
            name = 'spherical_elevation'
        case (coord_mode_cartesian_xy)
            name = 'cartesian_xy'
        case (coord_mode_cartesian_xyz)
            name = 'cartesian_xyz'
        case default
            name = 'unknown'
        end select
    end function coord_mode_name

    pure logical function coord_mode_is_implemented(coord_mode)
        integer, intent(in) :: coord_mode

        coord_mode_is_implemented = (coord_mode == coord_mode_spherical_angular)
    end function coord_mode_is_implemented

    pure logical function coord_mode_uses_elevation(coord_mode)
        integer, intent(in) :: coord_mode

        coord_mode_uses_elevation = (coord_mode == coord_mode_spherical_elevation .or. &
            coord_mode == coord_mode_cartesian_xyz)
    end function coord_mode_uses_elevation

    pure real function elevation_path_distance(horizontal_dist_km, dz_km)
        real, intent(in) :: horizontal_dist_km, dz_km

        elevation_path_distance = sqrt(horizontal_dist_km * horizontal_dist_km + dz_km * dz_km)
    end function elevation_path_distance

    pure real function elevation_path_scale(horizontal_dist_km, dz_km)
        real, intent(in) :: horizontal_dist_km, dz_km

        if (horizontal_dist_km <= 0.0) then
            elevation_path_scale = 1.0
        else
            elevation_path_scale = elevation_path_distance(horizontal_dist_km, dz_km) / horizontal_dist_km
        endif
    end function elevation_path_scale

    pure real function latdeg_to_colatrad(lat_deg)
        real, intent(in) :: lat_deg

        latdeg_to_colatrad = (90.0 - lat_deg) * coord_pi / 180.0
    end function latdeg_to_colatrad

    pure real function londeg_to_lonrad(lon_deg)
        real, intent(in) :: lon_deg

        londeg_to_lonrad = lon_deg * coord_pi / 180.0
    end function londeg_to_lonrad

    pure real function colatrad_to_latdeg(theta_rad)
        real, intent(in) :: theta_rad

        colatrad_to_latdeg = 90.0 - theta_rad * 180.0 / coord_pi
    end function colatrad_to_latdeg

    pure real function lonrad_to_londeg(phi_rad)
        real, intent(in) :: phi_rad

        lonrad_to_londeg = phi_rad * 180.0 / coord_pi
    end function lonrad_to_londeg

    pure subroutine great_circle_distance_colat_rad(theta1, phi1, theta2, phi2, dist_km)
        real, intent(in) :: theta1, phi1, theta2, phi2
        real, intent(out) :: dist_km
        real :: dtheta, dphi, lat1, lat2, a, c

        dtheta = theta2 - theta1
        dphi = phi2 - phi1
        lat1 = coord_pi / 2.0 - theta1
        lat2 = coord_pi / 2.0 - theta2
        a = sin(dtheta / 2.0) * sin(dtheta / 2.0) + &
            sin(dphi / 2.0) * sin(dphi / 2.0) * cos(lat1) * cos(lat2)
        a = min(1.0, max(0.0, a))
        c = 2.0 * atan2(sqrt(a), sqrt(1.0 - a))
        dist_km = earth_radius_km * c
    end subroutine great_circle_distance_colat_rad

end module coord_utils
