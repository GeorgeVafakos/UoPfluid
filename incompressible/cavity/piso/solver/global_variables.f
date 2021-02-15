        !-----------------------------------------------------------------------------------------------------------------------------------------------------------
        !    ______    ______        _________                                                !     
        !      ||        ||            ||    \\     ___                             ____      !     Version:    2.0
        !      ||        ||            ||     ||   // \\ ==||                        ||       !     
        !      ||        ||            ||     ||  ||  //   ||                        ||       !     Creator:    George Vafakos
        !      ||        ||    ____    ||____//   ||       ||  ____  ____   ο    ____||       !     
        !      ||        ||   //  \\   ||       ==||==     ||   ||    ||  =||   //   ||       !     Site:       https://github.com/GeorgeVafakos/UoPfluid
        !      ||        ||  ||    ||  ||         ||       ||   ||    ||   ||  ||    ||       !
        !       \\______//    \\__// __||__     __||__   __||__  \\__//   _||_  \\___||_      !
        !                                                                                     !
        !-----------------------------------------------------------------------------------------------------------------------------------------------------------
        ! 
        ! 
        !   File name:          global_variables
        !   
        !   Type:               module
        ! 
        !   Description:        This is the file where all global variables are declared.
        ! 
        ! 


        module global_variables
            use Class_Axis

            ! Declare axis
            type (cartesian_axis) :: x
            type (cartesian_axis) :: y
            type (cartesian_axis) :: z
            type (global_axis) :: N

            ! Geometric variables
            integer i, j, k, ID
            real*8  x_begin, x_end, LengthX
            integer cells_x, NumberFaces_x, NumberNodes_x
            real*8  y_begin, y_end, LengthY
            integer cells_y, NumberFaces_y, NumberNodes_y
            real*8  z_begin, z_end, LengthZ
            integer cells_z, NumberFaces_z, NumberNodes_z
            integer TotalFaces, TotalNodes, TotalCells
            real*8, allocatable, dimension(:)     :: Dx, Dy, Dz, dxe, dxw, dyn, dys, dzt, dzb
            real*8, allocatable, dimension(:)     :: Dx_vec, Dy_vec, Dz_vec, dxe_vec, dxw_vec, dyn_vec, dys_vec, dzt_vec, dzb_vec
            real*8, allocatable, dimension(:,:)   :: Dx_2D, Dy_2D, Dz_2D, dxe_2D, dxw_2D, dyn_2D, dys_2D, dzt_2D, dzb_2D
            real*8, allocatable, dimension(:,:,:) :: Dx_3D, Dy_3D, Dz_3D, dxe_3D, dxw_3D, dyn_3D, dys_3D, dzt_3D, dzb_3D
            integer, allocatable, dimension(:)    :: glob_index
            integer, allocatable, dimension(:)    :: nodes_inner, nodes_T, nodes_N, nodes_E, nodes_P, nodes_W, nodes_S, nodes_B
            integer k_boole, NumberBoundaries
            logical boundary_node, flow_2D

            ! Read input variables
            character(len=20) :: advection_scheme
            character(len=20) :: line, start, dir_name
            character(len=50) :: X_stretch, Y_stretch, Z_stretch
            character(len=20), allocatable, dimension(:) :: boundary_names
            real*8, allocatable, dimension(:)   :: mat_inlet
            logical file_exists, end_loop, nondim_plots, adjustable_time_step
            integer start_time

            ! Flow parameters
            real*8 Dt, nu, rho, Re, r, V_char

            ! State variables
            real*8, allocatable, dimension(:)   :: Var_old

            ! Finite volume coefficients
            real*8, allocatable, dimension(:)   :: A, aP_time
            real*8, allocatable, dimension(:)   :: aT_lapl, aN_lapl, aE_lapl, aP_lapl, aW_lapl, aS_lapl, aB_lapl
            real*8, allocatable, dimension(:)   :: aT_conv, aN_conv, aE_conv, aP_conv, aW_conv, aS_conv, aB_conv
            real*8, allocatable, dimension(:)   :: aE_gradx, aP_gradx, aW_gradx, aN_grady, aP_grady, aS_grady, aT_gradz, aP_gradz, aB_gradz
            real*8, allocatable, dimension(:)   :: aE_divx, aP_divx, aW_divx, aN_divy, aP_divy, aS_divy, aT_divz, aP_divz, aB_divz

            ! Other parameters
            real*8  :: Time = 0, cputime_start, cputime_finish, time_cpu
            integer :: Iter_count = 0, PISO_counter
            real*8, allocatable, dimension(:)   :: Err_Vx, Err_Vy, Err_Vz, Err_p, Err_pr, Err_pr_v
            integer, allocatable, dimension(:)  :: Iterations
            integer print_results, print_to_screen
            integer :: unit = 10
            real*8 :: g_value = 9.81
            integer num_threads


        end module
