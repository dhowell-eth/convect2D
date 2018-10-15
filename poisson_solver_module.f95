module poissonSolver
  ! Model with various subroutines used in the "convect" model
  
  CONTAINS

  function iteration_2DPoisson(u,f,h,alpha,c) result(rms_residue)
    implicit none
    real,intent(inout) :: u(:,:)
    real,intent(in) :: f(:,:),h,alpha,c
    real :: res_terms,indv_residue,rms_residue
    integer :: i,j
    real :: del2_u,nx,ny
    real :: min_res,max_res
    real,allocatable :: test_d2u(:,:)
    ! Reset terms b/c of multiple function calls
    res_terms = 0.0
    indv_residue = 0.0
    ! Calculate Current Residue
      do i=2,size(u,1)-1
        do j=2,size(u,2)-1
          ! For each point:
          ! Get finite diff approximation

          del2_u = (u(i-1,j) + u(i+1,j)+ u(i,j-1)+u(i,j+1)-((4.0+c*(h**2.0))*u(i,j)))/(h**2.0)
          ! Calculate residue 
          indv_residue = del2_u-f(i,j)
          
          ! Update u with new estimate
          u(i,j) = u(i,j)+(alpha*indv_residue*((h**2.0)/(4.0+c*h**2.0)))

          ! Store residue term for RMS calculation
          res_terms = res_terms + (indv_residue**2.0)
        end do
      end do
    ! Calculate RMS residue
    
    nx = size(u,1)-2.0
    ny = size(u,2)-2.0
    
    rms_residue = sqrt(res_terms/(nx*ny))

    ! Reset terms b/c of multiple function calls
    res_terms = 0.0
    indv_residue = 0.0



  end function iteration_2DPoisson

! Computes rms of a 2d array
function rms_2D(a) result(rms)
  implicit none
  real,intent(in) :: a(:,:)
  real :: running_total,nx,ny,rms
  integer :: i,j
  running_total=0.0
  nx=0.0
  ny=0.0

  do i=2,size(a,1)-1
    do j=2,size(a,2)-1
      running_total = running_total + a(i,j)**2.0
    end do
  end do
  nx = float(size(a,1))
  ny = float(size(a,2))
  rms = sqrt(running_total/(nx*ny))

  running_total=0.0
end function rms_2D

  subroutine residue_2DPoisson(u,f,h,res,c)
      implicit none
      real, intent(in) :: u(:,:),f(:,:),h,c
      real, intent(inout) :: res(:,:)
      integer :: i,j
      real :: del2_u
      do i=2,size(u,1)-1
        do j=2,size(u,2)-1
          del2_u = (1/h**2.0) * ( u(i,j+1) + u(i,j-1) + u(i+1,j) + u(i-1,j) - ((4.0+c*(h**2.0))*u(i,j)) )
          res(i,j) = del2_u-f(i,j)
        end do
      end do
  end subroutine residue_2DPoisson

  subroutine restrict(fine,coarse)
    implicit none
    real,intent(in) :: fine(:,:)
    real,intent(inout) :: coarse(:,:)
    integer :: i,j

    do i=1,size(coarse,1)
      do j=1,size(coarse,2)
        coarse(i,j) = fine((i-1)+i,(j-1)+j)
      end do
    end do
  end subroutine restrict
!!

! Only works for grids, use another interpolation if you want a value off the grid
subroutine linear_interpolate(grid_cell)
  implicit none
  real,intent(inout) :: grid_cell(:,:)
  integer :: i

  ! Quick and dirty "interpolation" using averages
  grid_cell(1,2) = (grid_cell(1,3)+grid_cell(1,1))/2.0
  grid_cell(2,1) = (grid_cell(3,1)+grid_cell(1,1))/2.0
  grid_cell(2,3) = (grid_cell(1,3)+grid_cell(3,3))/2.0
  grid_cell(3,2) = (grid_cell(3,3)+grid_cell(3,1))/2.0
  grid_cell(2,2) = (grid_cell(1,1)+grid_cell(1,3)+grid_cell(3,1)+grid_cell(3,3))/4.0 

 end subroutine linear_interpolate   

subroutine prolongate(coarse,fine)
    implicit none
    real,intent(inout) :: fine(:,:)
    real,intent(in) :: coarse(:,:)
    integer :: i,j
  ! Copy every other point
    do i=1,size(coarse,1)
      do j=1,size(coarse,2)
        fine((i-1)+i,(j-1)+j) = coarse(i,j)
      end do
    end do
  ! Interpolate between coarse points (Using a simple linear interpolation)
    do i=1,size(fine,1)-2,2
      do j=1,size(fine,2)-2,2
        call linear_interpolate(fine(i:i+2,j:j+2))
      end do
    end do
end subroutine prolongate

  recursive function Vcycle_2DPoisson(u_f,rhs,h,c,use_temp_bc) result (resV)
    implicit none
    real resV
    real,intent(inout):: u_f(:,:)  ! arguments
    real,intent(in)   :: rhs(:,:),h,c
    logical ,intent(in) :: use_temp_bc
    integer         :: nx,ny,nxc,nyc, i,j  ! local variables
    real,allocatable:: res_c(:,:),corr_c(:,:),res_f(:,:),corr_f(:,:)
    real            :: alpha=0.7, res_rms

    nx=size(u_f,1); ny=size(u_f,2)  ! must be power of 2 plus 1
    nxc=1+(nx-1)/2; nyc=1+(ny-1)/2  ! coarse grid size

    !print*,"Inside function (c):",c

    if (min(nx,ny)>5) then  ! not the coarsest level

       allocate(res_f(nx,ny),corr_f(nx,ny), &
            corr_c(nxc,nyc),res_c(nxc,nyc))

       !---------- take 2 iterations on the fine grid--------------
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha,c) 
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha,c)
        !print*,"First two iters: ",res_rms
       !---------- restrict the residue to the coarse grid --------
       call residue_2DPoisson(u_f,rhs,h,res_f,c) 
       call restrict(res_f,res_c)

       !---------- solve for the coarse grid correction -----------
       corr_c = 0.  
       res_rms = Vcycle_2DPoisson(corr_c,res_c,h*2,c,use_temp_bc) ! *RECURSIVE CALL*

       !---- prolongate (interpolate) the correction to the fine grid 
       call prolongate(corr_c,corr_f)
       ! If switch is set, apply non-zero boundary conditions
       if (use_temp_bc) corr_f(:,1) = 1.0
       !---------- correct the fine-grid solution -----------------

       !print*,"corr_c:",MAXVAL(corr_c)
       !print*,res_rms
       !print*,"corr_f:",maxval(corr_f)
       u_f = u_f - corr_f  


       !---------- two more smoothing iterations on the fine grid---
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha,c)
       res_rms = iteration_2DPoisson(u_f,rhs,h,alpha,c)

       deallocate(res_f,corr_f,res_c,corr_c)

    else  

       !----- coarsest level (ny=5): iterate to get 'exact' solution

       do i = 1,100
          res_rms = iteration_2DPoisson(u_f,rhs,h,alpha,c)
       end do

    end if

    resV = res_rms   ! returns the rms. residue

  end function Vcycle_2DPoisson

end module poissonSolver