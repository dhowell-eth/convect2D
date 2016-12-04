! Written By: Dorran Howell
! 11/21/2016
! Numerical Modeling in Fortran
! HW 7
! Written for fortran95 syntax, compiled by gfortran

!! TODO add a use statement for the poisson solver module
module diffusion_2d
	use poissonSolver
	implicit none
	real :: pi = 3.14159265359

	type AdvectionGrid
		integer :: nx, ny
		real,allocatable :: T(:,:),phi(:,:),v_x(:,:),v_y(:,:)
		real :: h,ra
	end type AdvectionGrid

	contains

		!----------------------------------------------------------------------!
		subroutine set_boundary_conditions_2d(data,constant)
			! Description
	
			implicit none
			real, intent(inout) :: data(:,:)
			real, intent(in) :: constant
			integer :: ri,rf,ci,cf,i,j
			
			ri = 1
			rf = size(data,1)
			ci = 1
			cf = size(data,2)

			do i=ri,rf
				if (i == ri .or. i == rf) then
					do j=ci,cf
						data(i,j) = constant
					end do 
				else
					data(i,ci) = constant
					data(i,cf) = constant
				end if
			end do
		end subroutine

		!----------------------------------------------------------------------!
		function del_squared_2d(data,grid_spacing,boundary_value) result(output)
			! Calculate just the 2d gradient (this will have to be input into the diffusion equ
			! Uses centered finite differences
			! Return a 2D Array

			implicit none
			integer :: i,j
			real, intent(in) :: data(:,:),grid_spacing
			real, dimension(size(data,1),size(data,2)):: output
			real,optional,intent(in) :: boundary_value
			real :: boundary_constant

			! If no boundary condition provided, set equal to 0.0
			IF (.not. present(boundary_value)) then
				boundary_constant= 0.0
			else
				boundary_constant = boundary_value
			end if

			output = 0.0
			! Calculate 2nd derivative for input data array using centered method
			do i = 2,size(data,1)-1
				do j=2,size(data,2)-1
					output(i,j) = ((data(i-1,j) + data(i+1,j)+ & 
					data(i,j-1)+data(i,j+1)-(4.0*data(i,j))))/(grid_spacing**2.0)
				end do
			end do

		end function del_squared_2d

		!----------------------------------------------------------------------!
		! Routine for calculating the x (first dimension) derivative of a 2D Array
		! using finite differences.
		subroutine x_deriv_2D(data,output,grid_spacing)
			implicit none
			real,intent(in) :: data(:,:), grid_spacing
			real,intent(inout) :: output(:,:)
			integer :: i,j
			output=0.0
			do i=2,size(data,1)-1
				do j=2,size(data,2)-1
					output(i,j) = (data(i+1,j)-data(i-1,j))/(grid_spacing*2.0)
				end do
			end do
		end subroutine x_deriv_2D

		!----------------------------------------------------------------------!
		
		subroutine calculate_v_field(advection_grid)
			implicit none
			type(AdvectionGrid) :: advection_grid
			integer :: i,j

			! NOTE: Flipped X and Y here

			do i = 2,size(advection_grid%phi,1)-1
				do j=2,size(advection_grid%phi,2)-1
					advection_grid%v_x(i,j) = (advection_grid%phi(i,j+1) - advection_grid%phi(i,j-1))/(advection_grid%h*2.0)
					! TODO I think the - sign goes with v_y instead... [DONE]
					advection_grid%v_y(i,j) = -1*(advection_grid%phi(i+1,j) - advection_grid%phi(i-1,j))/(advection_grid%h*2.0)
				end do
			end do

			call set_boundary_conditions_2d(advection_grid%v_x,0.0)
			call set_boundary_conditions_2d(advection_grid%v_y,0.0)
		end subroutine calculate_v_field
		!----------------------------------------------------------------------!
		subroutine v_grad(advection_grid,output)
			implicit none
			integer :: i,j
			type(AdvectionGrid):: advection_grid
			real :: output(:,:),vx,vy,dTdx,dTdy

			do i = 2,size(advection_grid%T,1)-1
				do j =2,size(advection_grid%T,2)-1
					vx = advection_grid%v_x(i,j)
					vy = advection_grid%v_y(i,j)
					if (vx>0) then
						dTdx = (advection_grid%T(i,j)-advection_grid%T(i-1,j))/advection_grid%h
					else
						dTdx = (advection_grid%T(i+1,j)-advection_grid%T(i,j))/advection_grid%h
					end if
					if (vy>0) then
						dTdy = (advection_grid%T(i,j)-advection_grid%T(i,j-1))/advection_grid%h
					else
						dTdy = (advection_grid%T(i,j+1)-advection_grid%T(i,j))/advection_grid%h
					end if
					output(i,j) = vx*dTdx + vy*dTdy
				end do
			end do

			call set_boundary_conditions_2d(output,0.0)

		end subroutine v_grad
		!----------------------------------------------------------------------!
end module diffusion_2d