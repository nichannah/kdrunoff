
module read_netcdf_field_2d

   use netcdf

   implicit none

contains

   subroutine read_field_2d(field_out, file_name, field_name)

     real, intent(inout), dimension(:, :) :: field_out
     character(len=*), intent(in) :: file_name, field_name

     ! This will be the netCDF ID for the file and data variable.
     integer :: ncid, varid

     ! Open the file. NF90_NOWRITE tells netCDF we want read-only access
     ! to the file.
     call check( nf90_open(file_name, NF90_NOWRITE, ncid) )

     ! Get the varid of the data variable, based on its name.
     call check( nf90_inq_varid(ncid, field_name, varid) )

     ! Read the data.
     call check( nf90_get_var(ncid, varid, field_out) )

     ! Close the file, freeing all resources.
     call check( nf90_close(ncid) )

   end subroutine read_field_2d

   subroutine check(status)
     integer, intent ( in) :: status

     if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop "Stopped"
     end if
  end subroutine check

end module read_netcdf_field_2d

program test
  use kdrunoff
  use read_netcdf_field_2d, only : read_field_2d

  implicit none

  integer, parameter :: NX = 360
  integer, parameter :: NY = 300
  integer :: num_land_points
  real, dimension(:, :), allocatable :: x_t, y_t, mask, runoff
  logical, dimension(:, :), allocatable :: land_sea_mask
  real, dimension(:, :), allocatable :: old_runoff

  allocate(x_t(NX, NY), y_t(NX, NY), mask(NX, NY), runoff(NX, NY))
  allocate(land_sea_mask(NX, NY), old_runoff(NX, NY))

  call read_field_2d(x_t, 'test.nc', 'x_T')
  call read_field_2d(y_t, 'test.nc', 'y_T')
  call read_field_2d(mask, 'test.nc', 'wet')

  land_sea_mask(:, :) = .false.
  where (mask > 0.5) land_sea_mask = .true.

  call kdrunoff_init(land_sea_mask, x_t, y_t, num_land_points)

  ! Make a random runoff field.
  call random_number(runoff)
  old_runoff(:, :) = runoff(:, :)

  ! Make sure there is something on the land.
  if (sum((1 - mask(:, :))*runoff(:, :)) <= 0.0) then
    print*, 'No runoff of land.'
    stop 1
  endif

  call kdrunoff_remap(runoff, 1)

  ! Make sure there is nothing left of the land.
  if (sum((1 - mask(:, :))*runoff(:, :)) > 0.0) then
    print*, 'Runoff left on land.'
    stop 1
  endif

  ! Test that totals are the same.
  if (abs(sum(old_runoff) - sum(runoff)) > 1e15) then
    print*, 'Totals are different'
    print*, sum(old_runoff), sum(runoff)
    stop 1
  endif

  deallocate(x_t, y_t, mask, runoff, old_runoff)

  call kdrunoff_end()

end program test


