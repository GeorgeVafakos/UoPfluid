        !-----------------------------------------------------------------------------------------------------------------------------------------------------------
        !                                                                                     !
        !    ______    ______        _________                                                !   Version: 2.0
        !      ||        ||            ||    \\     ___                             ____      !
        !      ||        ||            ||     ||   // \\ ==||                        ||       !
        !      ||        ||            ||     ||  ||  //   ||                        ||       !
        !      ||        ||    ____    ||____//   ||       ||  ____  ____   ο    ____||       !
        !      ||        ||   //  \\   ||       ==||==     ||   ||    ||  =||   //   ||       !
        !      ||        ||  ||    ||  ||         ||       ||   ||    ||   ||  ||    ||       !
        !       \\______//    \\__// __||__     __||__   __||__  \\__//   _||_  \\___||_      !
        !                                                                                     !
        !                                                                                     !
        !-----------------------------------------------------------------------------------------------------------------------------------------------------------
        ! 
        ! 
        ! 
        ! 
        !   Solver:             pipe_PISO
        ! 
        !   Discription:        This is a multidimensional solver for the diffusion equation. The present solver can solve the diffusion equation for a scalar 
        !                       variable, but it can easily be transformed for a vector variable as well.
        ! 
        !   Algorithms:         Jacobi, Gauss-Seidel, SOR
        !                       
        ! 
        ! 
        ! 


        program main

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
        call cpu_time(cputime_start)

        do
            ! Time step starts
            Time = Time + Dt
            Iter_count = Iter_count + 1

            ! Solve diffusion equation
            call NS_eqn%solve(domain, V, Vx_BC, Vy_BC, Vz_BC, 100, 1.0d-6, 1.0d-6, 1.0d-6, 1.0d-12)

            ! Create A matrix
            call create_A()

            ! PISO algorithm
            PISO_counter = 0
            do
                PISO_counter = PISO_counter + 1

                ! Create H matrix
                call create_H()

                ! Calculate pressure
                call Pres_eqn%solve(domain, p, p_BC, 1, 1.0d-6, 1.0d-12)

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

            ! Update induction equation terms
            call update_induction_terms()

            ! Solve induction equation
            call Ind_eqn%solve(domain, Bi, Bix_BC, Biy_BC, Biz_BC, 100, 1.0d-6, 1.0d-6, 1.0d-6, 1.0d-12)

            ! CVP algorithm
            dB%x = 0.0
            dB%y = 0.0
            dB%z = 0.0
            CVP_counter = 0
            do
                CVP_counter = CVP_counter + 1

                B%x = Bi%x + Bo%x
                B%y = Bi%y + Bo%y
                B%z = Bi%z + Bo%z

                ! Create G
                call create_G()

                ! Solve the magnetic field correction equation
                call dB_eqn%solve(domain, dB, dBx_BC, dBy_BC, dBz_BC, 1, 1.0d-6, 1.0d-6, 1.0d-6, 1.0d-12)

                ! Correct magnetic field
                Bi%x = Bi%x + dB%x
                Bi%y = Bi%y + dB%y
                Bi%z = Bi%z + dB%z
                call Bi%BC_Adjust_x(domain, Bix_BC)
                call Bi%BC_Adjust_y(domain, Biy_BC)
                call Bi%BC_Adjust_z(domain, Biz_BC)

                ! Calculate magnetic divergence error
                call Bi%calc_mean_div_error

                if (abs(Bi%mean_div_error) <= 1.0d-3 .OR. CVP_counter==100) then
                    exit
                end if
            end do

            ! Calculate total magnetic field
            B%x = Bi%x + Bo%x
            B%y = Bi%y + Bo%y
            B%z = Bi%z + Bo%z

            ! Calculate electric current density
            Je%x(nodes_P) = (1.0/mu)*curlx(Bi)
            Je%y(nodes_P) = (1.0/mu)*curly(Bi)
            Je%z(nodes_P) = (1.0/mu)*curlz(Bi)

            ! Update heat transport terms
            call update_heat_terms()

            ! Solve heat transport equation
            call Heat_eqn%solve(domain, T, T_BC, 100, 1.0d-6, 1.0d-12)

            ! Calculate CFL condition or Adjust time step
            call Co%CFL_condition(V)
            call Co_mag%CFL_condition(B)

            ! Create Error arrays with increasing size
            call relative_residual(V, V_old, 1.0d-12)
            call relative_residual(p, p_old, 1.0d-12)
            call relative_residual(Bi, Bi_old, 1.0d-12)
            call relative_residual(T, T_old, 1.0d-12)
            call calculate_residual_arrays()

            ! Save results
            if (mod(Iter_count,print_results)==0) then
                call update_sim_info()
                call print_sim_results()
            end if

            ! Calculate elapsed CPU time
            call cpu_time(cputime_finish)
            time_cpu = cputime_finish - cputime_start

            ! Print to screen
            print '(A,I17,A15,F10.4,A16,F9.4,A16,F9.4,A10,F12.8)', 'Iter.', Iter_count,  'Co_mean =', Co%mean, 'Co_max =', Co%max, 'Co_Mag_max =', Co_mag%max, 'Dt =', Dt
            print '(A,A21,ES10.2,A16,I9)', 'Velocity: x-axis', 'Error u =', V%e(1), 'Loops =', NS_eqn%counter(1)
            print '(A,A21,ES10.2,A16,I9)', 'Velocity: y-axis', 'Error v =', V%e(2), 'Loops =', NS_eqn%counter(2)
            print '(A,A21,ES10.2,A16,I9)', 'Velocity: z-axis', 'Error w =', V%e(3), 'Loops =', NS_eqn%counter(3)
            print '(A,A29,ES10.2,A16,I9,A16,ES9.2)', 'Pressure', 'Error p =', p%e, 'PISO Loops =', PISO_counter, 'Div_Err = ', V%mean_div_error
            print '(A,A15,ES10.2,A16,I9)', 'Magnetic Field: x-axis', 'Error Bx =', Bi%e(1), 'Loops =', Ind_eqn%counter(1)
            print '(A,A15,ES10.2,A16,I9)', 'Magnetic Field: y-axis', 'Error By =', Bi%e(2), 'Loops =', Ind_eqn%counter(2)
            print '(A,A15,ES10.2,A16,I9)', 'Magnetic Field: z-axis', 'Error Bz =', Bi%e(3), 'Loops =', Ind_eqn%counter(3)
            print '(A,A19,ES9.2,A16,I9,A15,ES9.2)', 'Magnetic Divergence', 'Mag_Div_Err = ', Bi%mean_div_error, 'CVP Loops =', CVP_counter
            print '(A,A26,ES10.2,A16,I9)', 'Temperature', 'Error T =', T%e, 'Loops =', Heat_eqn%counter
            print '(A,F10.5,A4,/)', 'CPU time =',time_cpu,'sec'

            ! Start elapsed cpu time for the next time step
            call cpu_time(cputime_start)

            ! Check if fully developed flow
            if (V%e(1)<=1.0d-6 .AND. V%e(2)<=1.0d-6 .AND. V%e(3)<=1.0d-6 .AND. p%e<=1.0d-6 .AND. Bi%e(1)<=1.0d-6 .AND. Bi%e(2)<=1.0d-6 .AND. Bi%e(3)<=1.0d-6 .AND. T%e<=1.0d-6 .AND. NS_eqn%counter(1)<=1 .AND. NS_eqn%counter(2)<=1 .AND. NS_eqn%counter(3)<=1 .AND. PISO_counter<=1 .AND. Ind_eqn%counter(1)<=1 .AND. Ind_eqn%counter(2)<=1 .AND. Ind_eqn%counter(3)<=1 .AND. CVP_counter<=1 .AND. Heat_eqn%counter<=1) then
                print *, 'Fully developed flow'
                call print_sim_results_final()
                call update_sim_info()
                exit
            else
                call update_equation_terms()
            end if
        end do

        end program main
