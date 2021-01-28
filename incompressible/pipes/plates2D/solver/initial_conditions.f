        !-----------------------------------------------------------------------------------------------------------------------------------------------------------
        !                                                                                     !
        !    ______    ______        _________                                                !   Version: 1.0
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
        !   File name:          initial_conditions
        ! 
        !   Type:               source
        ! 
        !   Description:        In this file the initial conditions of the variables and the starting values of parameters are declared.
        ! 
        ! 
        ! 
        !-----------------------------------------------------------------------------------------------------------------------------------------------------------
        ! Initial Conditions
        !-----------------------------------------------------------------------------------------------------------------------------------------------------------
        ! Check if input file exists
        call chdir('../results')
        inquire(file = '.simulation_info', exist = file_exists)
        call chdir('../solver')



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

            ! Navier-Stokes Maxwell tensions terms
            NS_max_tens%create_terms => CD_scheme

            ! Induction equation advection terms
            Ind_adv_term%create_terms => CD_scheme

            ! Induction equation Bi advection terms
            Ind_Bi_adv%create_terms => CD_scheme

            ! Induction equation Bo advection terms
            Ind_Bo_adv%create_terms => CD_scheme

            ! Heat transport equation convection terms
            Heat_conv_term%create_terms => CD_scheme

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

            ! Induction equation
            if (Ind_scheme=='Jacobi') then
                Ind_eqn%solve => jacobi_vector
            else if (Ind_scheme=='Gauss_Seidel') then
                Ind_eqn%solve => gauss_seidel_vector
            else if (Ind_scheme == 'SOR') then
                Ind_eqn%solve => SOR_vector
            end if

            ! Magnetic field correction equation
            if (dB_scheme=='Jacobi') then
                dB_eqn%solve => jacobi_vector
            else if (dB_scheme=='Gauss_Seidel') then
                dB_eqn%solve => gauss_seidel_vector
            else if (dB_scheme == 'SOR') then
                dB_eqn%solve => SOR_vector
            end if

            ! Heat transport Poisson
            if (Heat_scheme=='Jacobi') then
                Heat_eqn%solve => jacobi_scalar
            else if (Heat_scheme=='Gauss_Seidel') then
                Heat_eqn%solve => gauss_seidel_scalar
            else if (Heat_scheme == 'SOR') then
                Heat_eqn%solve => SOR_scalar
            end if

            ! Time step control
            if (adjustable_time_step .eqv. .TRUE.)  then
                Co%CFL_condition => calculate_Adj_Dt
            else
                Co%CFL_condition => calculate_CFL
            end if

            Co_mag%CFL_condition => calculate_CFL

            !------------------------------------------------------------------------------------------
            ! Prepare Boundary Conditions for Velocity
            !------------------------------------------------------------------------------------------
            ! Define which boundary condition will act on each boundary for the given variable
            call define_BC(domain, V, Vx_BC, Vy_BC, Vz_BC)

            if (V%typeBC(inlet)=='custom_inlet')  then
                call csvread(Vz_BC(inlet)%custom_inlet_file,mat_inlet)
                ID = 0
                do j = y%node_beg, y%node_end
                    do i = x%node_beg, x%node_end
                        if (i/=x%node_beg .AND. i/=x%node_end .AND. j/=y%node_beg .AND. j/=y%node_end) then
                            ID = ID + 1
                            Vz_BC(inlet)%value(ID) = mat_inlet(i+(j-1)*x%NumberNodes)
                        end if
                    end do
                end do
                Vz_BC(inlet)%value = Vz_BC(inlet)%value*Vz_BC(inlet)%custom_inlet_coeff
            end if

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
            ! Prepare Boundary Conditions for Induced Magnetic Field
            !------------------------------------------------------------------------------------------
            ! Define which boundary condition will act on each boundary for the given variable
            call define_BC(domain, Bi, Bix_BC, Biy_BC, Biz_BC)

            if (Bi%typeBC(inlet)=='custom_inlet')  then
                call csvread(Biz_BC(inlet)%custom_inlet_file,mat_inlet)
                ID = 0
                do j = y%node_beg, y%node_end
                    do i = x%node_beg, x%node_end
                        if (i/=x%node_beg .AND. i/=x%node_end .AND. j/=y%node_beg .AND. j/=y%node_end) then
                            ID = ID + 1
                            Biz_BC(inlet)%value(ID) = mat_inlet(i+(j-1)*x%NumberNodes)
                        end if
                    end do
                end do
                Biz_BC(inlet)%value = Biz_BC(inlet)%value*Biz_BC(inlet)%custom_inlet_coeff
            end if
            
            ! ! Adjust boundary conditions
            ! call Bi%BC_Adjust_x(domain, Bix_BC)
            ! call Bi%BC_Adjust_y(domain, Biy_BC)
            ! if (flow_2D .eqv. .FALSE.)  then
            !     call Bi%BC_Adjust_z(domain, Biz_BC)
            ! end if
            
            ! ! Prepare for first time step
            ! Bi_old = Bi

            !------------------------------------------------------------------------------------------
            ! Prepare Boundary Conditions for Outer Magnetic Field
            !------------------------------------------------------------------------------------------
            if (Bo%typeInit=='smooth_step') then
                if (flow_2D .eqv. .FALSE.)  then
                    Bo_function = (Bo_mag/2)*(1.0 + tanh(Bo_slope*(z%nodes-Bo_position)))
                    Bo%x = 0.0
                    Bo%y= pack(spread(Bo_function,1,x%NumberNodes*y%NumberNodes),.TRUE.)
                    Bo%z = 0.0
                else if (flow_2D .eqv. .TRUE.)  then
                    Bo_function = (Bo_mag/2)*(1.0 + tanh(Bo_slope*(y%nodes-Bo_position)))
                    Bo%x = 0.0
                    Bo%y= pack(spread(Bo_function,2,y%NumberNodes),.TRUE.)
                end if
            else
                Bo%x = Bo%init(1)
                Bo%y = Bo%init(2)
                if (flow_2D .eqv. .FALSE.)  Bo%z = Bo%init(3)
            end if

            ! ! Prepare for first time step
            ! B%x = Bi%x + Bo%x
            ! B%y = Bi%y + Bo%y
            ! B%z = Bi%z + Bo%z


            !------------------------------------------------------------------------------------------
            ! Prepare Boundary Conditions for Magnetic Field Correction
            !------------------------------------------------------------------------------------------
            ! Define which boundary condition will act on each boundary for the given variable
            call define_BC(domain, dB, dBx_BC, dBy_BC, dBz_BC)
            
            ! ! Adjust boundary conditions
            ! call dB%BC_Adjust_x(domain, dBx_BC)
            ! call dB%BC_Adjust_y(domain, dBy_BC)
            ! if (flow_2D .eqv. .FALSE.)  then
            !     call dB%BC_Adjust_z(domain, dBz_BC)
            ! end if


            !------------------------------------------------------------------------------------------
            ! Prepare Boundary Conditions for Temperature
            !------------------------------------------------------------------------------------------
            ! Define which boundary condition will act on each boundary for the given variable
            call define_BC(domain, T, T_BC)



        if (start == 'latest_time' .AND. file_exists .eqv. .TRUE.) then
            call chdir('../results')
            Time = obtain_variable('.simulation_info','last_saved_time')
            Iter_count = obtain_variable('.simulation_info','last_saved_iter')
            Dt = obtain_variable('.simulation_info','Dt')
            allocate(Err_Vx(Iter_count), Err_Vy(Iter_count), Err_p(Iter_count), Err_Bix(Iter_count), Err_Biy(Iter_count), Err_T(Iter_count), Iterations(Iter_count) )
            if (flow_2D .eqv. .FALSE.)  then
                allocate(Err_Vz(Iter_count), Err_Biz(Iter_count))
            end if

            ! Load initial conditions from last saved results
            call enter_last_result_dir()
            call csvread('Iterations.csv',Iterations)
            call csvread('u.csv',V%x)
            call csvread('v.csv',V%y)
            call csvread('p.csv',p%field)
            call csvread('Bix.csv',Bi%x)
            call csvread('Biy.csv',Bi%y)
            call csvread('T.csv',T%field)
            call csvread('Err_u.csv',Err_Vx)
            call csvread('Err_v.csv',Err_Vy)
            call csvread('Err_p.csv',Err_p)
            call csvread('Err_T.csv',Err_T)
            call csvread('Err_Bix.csv',Err_Bix)
            call csvread('Err_Biy.csv',Err_Biy)
            if (flow_2D .eqv. .FALSE.)  then
                call csvread('w.csv',V%z)
                call csvread('Biz.csv',Bi%z)
                call csvread('Err_w.csv',Err_Vz)
            end if
            call chdir('../../solver')

        else 
            
            !------------------------------------------------------------------------------------------
            ! Initial Boundary Conditions for Velocity
            !------------------------------------------------------------------------------------------
            ! Internal field initial conditions
            if (V%typeInit=='same_as_inlet') then
                V%x = 0.0
                V%y = 0.0
                V%z(nodes_P) = pack(spread(Vz_BC(inlet)%value,2,z%NumberCells),.TRUE.)
            else
                V%x = V%init(1)
                V%y = V%init(2)
                if (flow_2D .eqv. .FALSE.)  then
                    V%z = V%init(3)
                end if
            end if

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
            ! Initial Conditions for Induced Magnetic Field
            !------------------------------------------------------------------------------------------
            ! Internal field initial conditions
            Bi%x = Bi%init(1)
            Bi%y = Bi%init(2)
            if (flow_2D .eqv. .FALSE.)  then
                Bi%z = Bi%init(3)
            end if

            ! Adjust boundary conditions
            call Bi%BC_Adjust_x(domain, Bix_BC)
            call Bi%BC_Adjust_y(domain, Biy_BC)
            if (flow_2D .eqv. .FALSE.)  then
                call Bi%BC_Adjust_z(domain, Biz_BC)
            end if

            !------------------------------------------------------------------------------------------
            ! Initial Conditions for Magnetic Field Correction
            !------------------------------------------------------------------------------------------
            ! Internal field initial conditions
            dB%x = dB%init(1)
            dB%y = dB%init(2)
            if (flow_2D .eqv. .FALSE.)  then
                dB%z = dB%init(3)
            end if

            !------------------------------------------------------------------------------------------
            ! Initial Boundary Conditions for Temperature
            !------------------------------------------------------------------------------------------
            ! Internal field initial conditions
            T%field = p%init

             ! Adjust boundary conditions
            call T%BC_Adjust(domain, T_BC)

            !------------------------------------------------------------------------------------------
            ! Create Results directory
            !------------------------------------------------------------------------------------------
            call system("rm -rf ../results")
            call system("mkdir -p ../results")

        end if

            ! Prepare for first time step
            V_old = V
            p_old = p
            Bi_old = Bi
            T_old = T

            ! Prepare for first time step
            B%x = Bi%x + Bo%x
            B%y = Bi%y + Bo%y
            B%z = Bi%z + Bo%z
 