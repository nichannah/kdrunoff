
module kdrunoff

  use kdtree2_module

  implicit none
  private

  real(kdkind), dimension(:,:), allocatable :: ocean_points
  real(kdkind), dimension(:,:), allocatable :: ocean_indices
  real(kdkind), dimension(:,:), allocatable :: land_points
  real(kdkind), dimension(:,:), allocatable :: land_indices
  type(kdtree2), pointer :: tree

  integer, parameter :: D = 2

  public kdrunoff_init, kdrunoff_remap, kdrunoff_end

contains

  subroutine kdrunoff_init(land_sea_mask, x_t, y_t, num_land_points)
    real, dimension(:, :), intent(in) :: land_sea_mask, x_t, y_t
    integer, intent(out) :: num_land_points

    integer :: nx, ny, n, i, j, n_ocn, n_land

    ! Total number of ocean points.
    n = sum(land_sea_mask)
    nx = size(land_sea_mask, 1)
    ny = size(land_sea_mask, 2)

    num_land_points = nx*ny - n

    allocate(ocean_points(D, n))
    allocate(ocean_indices(D, n))
    allocate(land_points(D, num_land_points))
    allocate(land_indices(D, num_land_points))

    ! Make lists of ocean and land points, also indices to those points.
    n_ocn = 1
    n_land = 1
    do j=1,ny
      do i=1,nx
        if (land_sea_mask(i, j) > 0.5) then
          ocean_points(1, n_ocn) = x_t(i, j)
          ocean_points(2, n_ocn) = y_t(i, j)
          ocean_indices(1, n_ocn) = i
          ocean_indices(2, n_ocn) = j
          n_ocn = n_ocn + 1
        else
          land_points(1, n_land) = x_t(i, j)
          land_points(2, n_land) = y_t(i, j)
          land_indices(1, n_land) = i
          land_indices(2, n_land) = j
          n_land = n_land + 1
        endif
      enddo
    enddo

    ! Create kdtree data structure
    tree => kdtree2_create(ocean_points, sort=.false., rearrange=.false.)

  end subroutine kdrunoff_init

  subroutine kdrunoff_remap(runoff, rspread)
    real, dimension(:, :), intent(inout) :: runoff
    integer, intent(in) :: rspread

    integer :: n, i, j, nni, nnj
    real :: val
    type(kdtree2_result), allocatable :: results(:)

    allocate(results(rspread))

    do n=1, size(land_points, 2)
      i = land_indices(1, n)
      j = land_indices(2, n)
      val = runoff(i, j)
      if (val > 0.0) then
        call kdtree2_n_nearest(tp=tree, qv=land_points(:, i), nn=rspread, &
                               results=results)
        runoff(i, j) = runoff(i, j) - val
        nni = ocean_indices(1, results(1)%idx)
        nnj = ocean_indices(2, results(1)%idx)
        runoff(nni, nnj) = runoff(nni, nnj) + val
      endif
    enddo

    deallocate(results)
  end subroutine kdrunoff_remap

  subroutine kdrunoff_end()

    call kdtree2_destroy(tree)

    deallocate(ocean_points)
    deallocate(ocean_indices)
    deallocate(land_points)
    deallocate(land_indices)

  end subroutine kdrunoff_end

end module kdrunoff

