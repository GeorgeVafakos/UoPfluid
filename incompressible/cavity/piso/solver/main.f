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
        !   Solver:             cavity
        ! 
        !   Discription:        This is an incompressible solver for the solution of the Navier-Stokes equation, using the PISO algorithm.
        ! 
        !   Geometry:           Lid-driven cavity 2D incompressible flow.
        ! 
        !   Algorithms:         Jacobi, Gauss-Seidel, SOR
        !                       
        ! 


        program main

        use omp_lib
        use global_variables
        use Class_Advection
        use Class_Axis
        use Class_CFL
        use Class_Geometry
        use Class_Scalar_Equation
        use Class_Scalar_Variable
        use Class_Vector_Equation
        use Class_Vector_Variable
        use define_classes
        use subprograms_case
        use subprograms_default
        use subprograms_print

        implicit none

        ! Read input file
        include 'read_input.f'

        ! Generate mesh
        include 'mesh_generation.f'

        ! Allocate classes and variables
        call allocate_variables()

        ! Load initial conditions
        include 'initial_conditions.f'

        ! Calculate nondimensional numbers
        call calculate_nondim_numbers()

        ! Print simulation info in screen
        call print_sim_info_screen()

        ! Insert initial equation terms
        include 'equations_terms.f'

        ! Print initial conditions
        call update_sim_info()
        call print_initial_conditions()

        ! Initalize CPU timer
        cputime_start = omp_get_wtime()

        do
            ! Time step starts
            Time = Time + Dt
            Iter_count = Iter_count + 1

            ! Solve Navier-Stokes equation
            call NS_eqn%solve(domain, V, Vx_BC, Vy_BC, Vz_BC, 100, 1.0d-6, 1.0d-6, 1.0d-6, 1.0d-8)

            ! Create A matrix
            call create_A()

            ! PISO algorithm
            PISO_counter = 0
            do
                PISO_counter = PISO_counter + 1

                ! Create H matrix
                call create_H()

                ! Calculate pressure
                call Pres_eqn%solve(domain, p, p_BC, 1, 1.0d-6, 1.0d-8)

                ! Update velocity
                V%x(nodes_P) = H%x(nodes_P) - (1.0/rho)*gradx(p%field)/A
                V%y(nodes_P) = H%y(nodes_P) - (1.0/rho)*grady(p%field)/A
                V%z(nodes_P) = H%z(nodes_P) - (1.0/rho)*gradz(p%field)/A
                call V%BC_Adjust_x(domain, Vx_BC)
                call V%BC_Adjust_y(domain, Vy_BC)
                call V%BC_Adjust_z(domain, Vz_BC)

                ! Calculate divergence error
                call V%calc_mean_div_error

                if (abs(V%mean_div_error) <= 1.0d-3 .OR. PISO_counter==100) then
                    exit
                end if
            end do

            ! Calculate CFL condition or Adjust time step
            call Co%CFL_condition(V)

            ! Create Error arrays with increasing size
            call relative_residual(V, V_old, 1.0d-12)
            call relative_residual(p, p_old, 1.0d-12)
            call calculate_residual_arrays()

            ! Save results
            if (mod(Iter_count,print_results)==0) then
                call update_sim_info()
                call print_sim_results()
            end if

            ! Calculate elapsed CPU time
            time_cpu = omp_get_wtime() - cputime_start

            ! Print to screen
            if (mod(Iter_count,print_to_screen)==0) then
                print '(A,I17,A15,F10.4,A16,F9.4,A10,F12.8)', 'Iter.', Iter_count,  'Co_mean =', Co%mean, 'Co_max =', Co%max, 'Dt =', Dt
                print '(A,A21,ES10.2,A16,I9)', 'Velocity: x-axis', 'Error u =', V%e(1), 'Loops =', NS_eqn%counter(1)
                print '(A,A21,ES10.2,A16,I9)', 'Velocity: y-axis', 'Error v =', V%e(2), 'Loops =', NS_eqn%counter(2)
                if (flow_2D .eqv. .FALSE.)  then
                print '(A,A21,ES10.2,A16,I9)', 'Velocity: z-axis', 'Error w =', V%e(3), 'Loops =', NS_eqn%counter(3)
                end if
                print '(A,A29,ES10.2,A16,I9,A16,ES9.2)', 'Pressure', 'Error p =', p%e, 'PISO Loops =', PISO_counter, 'Div_Err = ', V%mean_div_error
                print '(A,F15.5,A4,/)', 'CPU time =',time_cpu,'sec'
            end if

            ! Start elapsed cpu time for the next time step
            cputime_start = omp_get_wtime()

            ! Check if fully developed flow
            if (V%e(1)<=1.0d-6 .AND. V%e(2)<=1.0d-6 .AND. V%e(3)<=1.0d-6 .AND. p%e<=1.0d-4 .AND. NS_eqn%counter(1)<=1 .AND. NS_eqn%counter(2)<=1 .AND. NS_eqn%counter(3)<=1 .AND. PISO_counter<=1) then
                print *, 'Fully developed flow'
                call print_sim_results_final()
                call update_sim_info()
                exit
            else
                call update_equation_terms()
            end if
        end do

        end program main
