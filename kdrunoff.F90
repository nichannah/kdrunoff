
module read_netcdf_field_2d

   use netcdf

   implicit none

contains

   subroutine read_field_2d(field_out, field_name)

     real, intent(inout), dimension(:, :) :: field_out 
     character(len=*), intent(in) :: field_name

     ! This will be the netCDF ID for the file and data variable.
     integer :: ncid, varid

     ! Open the file. NF90_NOWRITE tells netCDF we want read-only access
     ! to the file.
     call check( nf90_open('grid_spec.nc', NF90_NOWRITE, ncid) )

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

program kdtree_test
  use kdtree2_module
  use read_netcdf_field_2d, only : read_field_2d

  implicit none

  integer, parameter :: NX = 360
  integer, parameter :: NY = 300
  integer, parameter :: D = 2
  real, dimension(:, :), allocatable :: x_t, y_t, mask, output

  integer :: i, j, n, m, nn
  integer :: lc, oc
  real(kdkind), dimension(:,:), allocatable :: ocean_points
  real(kdkind), dimension(:,:), allocatable :: ocean_indices
  real(kdkind), dimension(:,:), allocatable :: land_points
  real(kdkind), dimension(:,:), allocatable :: land_indices
  real(kdkind), allocatable :: query_vec(:)
  type(kdtree2), pointer :: tree

  type(kdtree2_result), allocatable :: results(:)

  allocate(x_t(NX, NY), y_t(NX, NY), mask(NX, NY), output(NX, NY))
  
  call read_field_2d(x_t, 'x_T')
  call read_field_2d(y_t, 'y_T')
  call read_field_2d(mask, 'wet')

  n = sum(mask)
  nn = 1

  allocate(ocean_points(D, n))
  allocate(ocean_indices(D, n))
  allocate(land_points(D, NX*NY - n))
  allocate(land_indices(D, NX*NY - n))
  allocate(query_vec(D))
  allocate(results(nn))

  print*, 'Making ocean and land points'

  oc = 1
  lc = 1
  do j=1,NY
    do i=1,NX
      if (mask(i, j) > 0.5) then
        ocean_points(1, oc) = x_t(i, j)
        ocean_points(2, oc) = y_t(i, j)
        ocean_indices(1, oc) = i
        ocean_indices(2, oc) = j
        oc = oc + 1
      else
        land_points(1, lc) = x_t(i, j)
        land_points(2, lc) = y_t(i, j)
        land_indices(1, lc) = i
        land_indices(2, lc) = j
        lc = lc + 1
      endif
    enddo
  enddo

  print*, 'ocean_points, land_points: ', size(ocean_points, 2), size(land_points, 2)
  tree => kdtree2_create(ocean_points, sort=.false., rearrange=.false.)

  print*, 'land_points(:, 10): ', land_points(:, 10)
  print*, 'land_indices(:, 10): ', land_indices(:, 10)
  query_vec(:) = land_points(:, 10)

  call kdtree2_n_nearest(tp=tree, qv=query_vec, nn=nn, results=results)

  print*, 'results(1)%idx: ', results(1)%idx
  print*, 'ocean_indices(:, idx): ', ocean_indices(:, results(1)%idx)

  !do i=1, size(land_points, 2)
  !  call kdtree2_n_nearest(distance, index)
  !enddo

  call kdtree2_destroy(tree)  

  deallocate(ocean_points)
  deallocate(land_points)

end program kdtree_test


