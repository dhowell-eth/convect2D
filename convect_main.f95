! Written By: Dorran Howell
! 11/21/2016
! Numerical Modeling in Fortran
! HW 7
! Written for fortran95 syntax, compiled by gfortran

program convect_main
	! This program simulates 2D Temperature convection/advection.
	! Adjust input parameters in the './convect_inputs.txt' file to modify the simulation
	! Only initial and final T profiles, as well as final stream function values, are written out (to CWD)

	use diffusion_2d
	implicit none
	real,allocatable :: init_temp(:,:),del2(:,:),vgrad(:,:),old_dtdx(:,:),w(:,:),s(:,:),vgrad_w(:,:),del2_w(:,:)
	integer :: nx=3,ny=3,i,j,centerx,centery
	character(len=50) :: outputfile_suffix,out_path,step_number,run_date,step_dir
	character(len=:),allocatable :: formatted_out_path,formatted_run_date
	real :: total_time,field_magnitude,diffusivity,dt,dt_dif,dt_adv,t,h,a_adv,a_dif,B,ra,placeholder,pr
	type(AdvectionGrid) :: this_advection_grid
	real :: rms_f=0.0,rms_residue=0.0,error_threshhold=0
	integer :: iter_counter,num_frames
	real :: dt_for_write
	real, allocatable :: frame_times(:)

	! Read in input parameters from a file
	namelist /inputs/ nx,ny,total_time,diffusivity,a_dif,a_adv,ra,pr,error_threshhold,num_frames
	open(1,file="convect_inputs.txt",status="old")
	read(1,inputs)
	close(1)

	! Print a welcome message
	print*,"______________________________________________________________"
	print*,"This program solves the 2D Convection-Advection Equation for a random and normal profile using a &
	centered finite difference approximation."
	print*,"Adjust parameters in ./convect_inputs.txt to modify the simulation."
	print*,"Executing..."

	! Allocate various arrays T, phi, V
	allocate(this_advection_grid%T(nx,ny),this_advection_grid%phi(nx,ny),this_advection_grid%v_x(nx,ny),this_advection_grid%v_y(nx,ny))
	allocate(init_temp(nx,ny),del2(nx,ny),vgrad(nx,ny))
	allocate(old_dtdx(nx,ny),w(nx,ny),s(nx,ny))
	allocate(this_advection_grid%w(nx,ny),del2_w(nx,ny),vgrad_w(nx,ny))

	! Initialize values for velocity arrays
	this_advection_grid%v_x = 0.0
	this_advection_grid%v_y = 0.0
	! Initialize omega 
	! Populate T field with random signal between 0:1
	this_advection_grid%w = 0.0
	do i=1,size(this_advection_grid%w,1)
		do j=1,size(this_advection_grid%w,2)
			this_advection_grid%w(i,j) = rand()
		end do
	end do
	! Set grid spacing (assuming dx = dy)
	h = 1.0/(float(ny)-1.0)
	this_advection_grid%h = h
	this_advection_grid%nx = nx
	this_advection_grid%ny = ny
	this_advection_grid%ra = ra

	! Populate T field with random signal between 0:1
	do i=1,size(this_advection_grid%T,1)
		do j=1,size(this_advection_grid%T,2)
			this_advection_grid%T(i,j) = rand()
		end do
	end do

	! Set boundary conditions on top and bottom
	this_advection_grid%T(:,1) = 1.0
	this_advection_grid%T(:,ny) = 0.0

	! Calculate diffusive dt (does not change with each iteration)
	dt_dif = a_dif*(this_advection_grid%h**2)/diffusivity
	


	! Figure out times to write out data
	dt_for_write = total_time/num_frames
	num_frames = num_frames-2.0
	allocate(frame_times(num_frames))
	do i=1,num_frames
		frame_times(i) = dt_for_write*float(i)
	end do

	!----------------------------------------------------------------------!
	! Run simulation

	t = 0.0
	iter_counter = 1
	!------------------------------------!
	! Write out initial T and phi fields
	open(42,file='T_field.dat')
	do i=1,size(this_advection_grid%T,1)
		write(42,*),t,this_advection_grid%T(i,:)
	end do
	open(43,file='phi_field.dat')
	do i=1,size(this_advection_grid%phi,1)
		write(43,*),t,this_advection_grid%phi(i,:)
	end do
	open(44,file='w_field.dat')
	do i=1,size(this_advection_grid%w,1)
		write(44,*),t,this_advection_grid%w(i,:)
	end do
	!------------------------------------!

	print*,num_frames,frame_times

	do
		!------------------------------------!
		! If this time is one you want in output file, write it.
		if ((t >= frame_times(iter_counter)) .and. (iter_counter <= num_frames)) THEN
			print*,"Writing..."
			! Write out T with tine
			do i=1,size(this_advection_grid%T,1)
				write(42,*),t,this_advection_grid%T(i,:)
			end do
			do i=1,size(this_advection_grid%phi,1)
				write(43,*),t,this_advection_grid%phi(i,:)
			end do
			do i=1,size(this_advection_grid%w,1)
				write(44,*),t,this_advection_grid%w(i,:)
			end do
				iter_counter = iter_counter + 1
		end if
		!------------------------------------!

		print*,"Percent Complete:",(t/total_time)*100,"%"

		!! --- Get inputs for solving for new T and W --- !!

		!------------------------------------!
		! Solve s/phi (streamfunction)

		s = 0.0
		! Get rms of the function
		rms_f = 0.0
		rms_f = rms_2D(this_advection_grid%w)
		rms_residue = rms_f*2.0

		do while (rms_residue/rms_f>error_threshhold)
			rms_residue = Vcycle_2DPoisson(s,this_advection_grid%w,this_advection_grid%h)
		end do	

		!-----------------!
		! Store streamfunction in grid

		this_advection_grid%phi = s

		! Set edge values for stream function
		call set_boundary_conditions_2d(this_advection_grid%phi,0.0)

		! Calculate velocities and set edge values
		call calculate_v_field(this_advection_grid)
		call set_boundary_conditions_2d(this_advection_grid%v_x,0.0)
		call set_boundary_conditions_2d(this_advection_grid%v_y,0.0)

		!-----------------!
		! Figure out time_step using minimum of advective and diffusive time steps
		
		dt_adv = a_adv*MIN(h/MAXVAL(ABS(this_advection_grid%v_x)),h/MAXVAL(ABS(this_advection_grid%v_y)))
		dt = MIN(dt_dif,dt_adv)

		!-----------------!
		! Get terms for calculating dT

		del2 = del_squared_2d(this_advection_grid%T,h)
		vgrad = 0.0
		call v_grad(this_advection_grid,vgrad)
		
		! Get terms for calculating dW
		del2_w = 0.0
		vgrad_w = 0.0
		del2_w = del_squared_2d(this_advection_grid%w,h)
		call v_grad_w(this_advection_grid,vgrad_w)

		!-----------------!

		! Calc dtDx for w time step
		old_dtdx=0.0
		call x_deriv_2D(this_advection_grid%T,old_dtdx,this_advection_grid%h)

		! Set BC for dtdx
		old_dtdx(1,:) = 0.0
		old_dtdx(this_advection_grid%nx,:) = 0.0

		!------------------------------------!	
		! Solve for new T and w

		! Get temperature at this time step
		this_advection_grid%T = this_advection_grid%T + dt*((diffusivity*del2)-vgrad)

		this_advection_grid%w = this_advection_grid%w + dt*(pr*del2_w - vgrad_w - ra*pr*old_dtdx)

		! Update boundaries
		this_advection_grid%T(1,:) = this_advection_grid%T(2,:)
		this_advection_grid%T(nx,:) = this_advection_grid%T(nx-1,:)
		!------------------------------------!

		! Increment time
		t = t + dt

		! Kill loop when it has lasted for the input duration
		if (t>total_time) exit

	end do
	!----------------------------------------------------------------------!

	! Write out final T and phi
	do i=1,size(this_advection_grid%T,1)
		write(42,*),t,this_advection_grid%T(i,:)
	end do
	close(42)
	do i=1,size(this_advection_grid%phi,1)
		write(43,*),t,this_advection_grid%phi(i,:)
	end do
	close(43)
	do i=1,size(this_advection_grid%w,1)
		write(44,*),t,this_advection_grid%w(i,:)
	end do
	close(44)

	! Print exit message
	print*,"Run complete!"
	print*,"______________________________________________________________"

	! Deallocate arrays
	deallocate(this_advection_grid%T,this_advection_grid%phi,this_advection_grid%v_x,this_advection_grid%v_y,init_temp,del2,vgrad)
	deallocate(old_dtdx,w,s,frame_times)
end program