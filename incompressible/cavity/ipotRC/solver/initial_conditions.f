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
        !   File name:          initial_conditions
        ! 
        !   Type:               source
        ! 
        !   Description:        In this file the initial conditions of the variables and the starting values of parameters are declared.
        ! 
        ! 
        !-----------------------------------------------------------------------------------------------------------------------------------------------------------
        ! Initial Conditions
        !-----------------------------------------------------------------------------------------------------------------------------------------------------------
        ! Check if simulation_info file exists
        call chdir('../results')
        inquire(file = '.simulation_info', exist = file_exists)
        call chdir('../solver')

        ! Define number of CPU threads are used for the OpenMP library
        call omp_set_num_threads(num_threads)

        !------------------------------------------------------------------------------------------
        ! Define discretization scheme for each term
        !------------------------------------------------------------------------------------------
        if (advection_scheme=='CD') then
            NS_adv_term%create_terms => CD_scheme
        else if (advection_scheme=='upwind') then
            NS_adv_term%create_terms => upwind_scheme
        else if (advection_scheme == 'none') then
            NS_adv_term%create_terms => none
        end if

        !------------------------------------------------------------------------------------------
        ! Define solution scheme for each equation
        !------------------------------------------------------------------------------------------
        ! Navier-Stokes
        if (NS_scheme=='Jacobi') then
            NS_eqn%solve => jacobi_vector
        else if (NS_scheme=='Gauss_Seidel') then
            NS_eqn%solve => gauss_seidel_vector
        else if (NS_scheme == 'SOR') then
            NS_eqn%solve => SOR_vector
        end if

        ! Pressure Poisson
        if (Pres_scheme=='Jacobi') then
            Pres_eqn%solve => jacobi_scalar
        else if (Pres_scheme=='Gauss_Seidel') then
            Pres_eqn%solve => gauss_seidel_scalar
        else if (Pres_scheme == 'SOR') then
            Pres_eqn%solve => SOR_scalar
        end if

        ! Time step control
        if (adjustable_time_step .eqv. .TRUE.)  then
            Co%CFL_condition => calculate_Adj_Dt
        else
            Co%CFL_condition => calculate_CFL
        end if

        !------------------------------------------------------------------------------------------
        ! Prepare Boundary Conditions for Velocity
        !------------------------------------------------------------------------------------------
        ! Define which boundary condition will act on each boundary for the given variable
        call define_BC(domain, V, Vx_BC, Vy_BC, Vz_BC)

        ! if (V%typeBC(inlet)=='custom_inlet')  then
        !     call csvread(Vz_BC(inlet)%custom_inlet_file,mat_inlet)
        !     ID = 0
        !     do j = y%node_beg, y%node_end
        !         do i = x%node_beg, x%node_end
        !             if (i/=x%node_beg .AND. i/=x%node_end .AND. j/=y%node_beg .AND. j/=y%node_end) then
        !                 ID = ID + 1
        !                 Vz_BC(inlet)%value(ID) = mat_inlet(i+(j-1)*x%NumberNodes)
        !             end if
        !         end do
        !     end do
        !     Vz_BC(inlet)%value = Vz_BC(inlet)%value*Vz_BC(inlet)%custom_inlet_coeff
        ! end if

            ! ! Internal field initial conditions
            ! if (V%typeInit=='same_as_inlet') then
            !     V%x = 0.0
            !     V%y = 0.0
            !     V%z(nodes_P) = pack(spread(Vz_BC(inlet)%value,2,z%NumberCells),.TRUE.)
            ! else
            !     V%x = V%init(1)
            !     V%y = V%init(2)
            !     if (flow_2D .eqv. .FALSE.)  then
            !         V%z = V%init(3)
            !     end if
            ! end if
            
            ! ! Adjust boundary conditions
            ! call V%BC_Adjust_x(domain, Vx_BC)
            ! call V%BC_Adjust_y(domain, Vy_BC)
            ! if (flow_2D .eqv. .FALSE.)  then
            !     call V%BC_Adjust_z(domain, Vz_BC)
            ! end if

            ! ! Prepare for first time step
            ! V_old = V


        !------------------------------------------------------------------------------------------
        ! Prepare Boundary Conditions for Pressure
        !------------------------------------------------------------------------------------------
        ! Define which boundary condition will act on each boundary for the given variable
        call define_BC(domain, p, p_BC)

            ! ! Internal field initial conditions
            ! p%field = p%init

            ! ! Adjust boundary conditions
            ! call p%BC_Adjust(domain, p_BC)

            ! ! Prepare for first time step
            ! p_old = p
            
            


        !------------------------------------------------------------------------------------------
        ! Impose Initial Conditions According to Starting Point (latest_time or start_time)
        !------------------------------------------------------------------------------------------
        if (start == 'latest_time' .AND. file_exists .eqv. .TRUE.) then
            call chdir('../results')
            Time = obtain_variable('.simulation_info','last_saved_time')
            Iter_count = obtain_variable('.simulation_info','last_saved_iter')
            Dt = obtain_variable('.simulation_info','Dt')
            allocate(Err_Vx(Iter_count), Err_Vy(Iter_count), Err_p(Iter_count), Iterations(Iter_count) )
            if (flow_2D .eqv. .FALSE.)  then
                allocate(Err_Vz(Iter_count))
            end if

            ! Load initial conditions from last saved results
            call enter_last_result_dir()
            call csvread('Iterations.csv',Iterations)
            call csvread('u.csv',V%x)
            call csvread('v.csv',V%y)
            call csvread('p.csv',p%field)
            call csvread('Err_u.csv',Err_Vx)
            call csvread('Err_v.csv',Err_Vy)
            call csvread('Err_p.csv',Err_p)
            if (flow_2D .eqv. .FALSE.)  then
                call csvread('w.csv',V%z)
                call csvread('Err_w.csv',Err_Vz)
            end if
            call chdir('../../solver')

        else 
            
            !------------------------------------------------------------------------------------------
            ! Initial Boundary Conditions for Velocity
            !------------------------------------------------------------------------------------------
            ! Internal field initial conditions
            ! if (V%typeInit=='same_as_inlet') then
                V%x = 0.0
                V%y = 0.0
                V%z(nodes_P) = pack(spread(Vz_BC(inlet)%value,2,z%NumberCells),.TRUE.)
            ! else
                V%x = V%init(1)
                V%y = V%init(2)
                if (flow_2D .eqv. .FALSE.)  then
                    V%z = V%init(3)
                end if
            ! end if

             ! Adjust boundary conditions
            call V%BC_Adjust_x(domain, Vx_BC)
            call V%BC_Adjust_y(domain, Vy_BC)
            if (flow_2D .eqv. .FALSE.)  then
                call V%BC_Adjust_z(domain, Vz_BC)
            end if

            !------------------------------------------------------------------------------------------
            ! Initial Boundary Conditions for Pressure
            !------------------------------------------------------------------------------------------
            ! Internal field initial conditions
            p%field = p%init

             ! Adjust boundary conditions
            call p%BC_Adjust(domain, p_BC)

            !------------------------------------------------------------------------------------------
            ! Create Results directory
            !------------------------------------------------------------------------------------------
            call system("rm -rf ../results")
            call system("mkdir -p ../results")

        end if

            ! Prepare for first time step
            V_old = V
            p_old = p
             
