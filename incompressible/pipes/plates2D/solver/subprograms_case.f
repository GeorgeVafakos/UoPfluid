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
        !   File name:          subprograms_case
        !   
        !   Type:               module
        ! 
        !   Description:        In this file are stored all functions and subroutines that are used in the code.
        !
        !
        !   Containing Subprogramms             Type                Short Description
        !   =======================             ====                =================
        ! 
        !   - allocate_variables                subroutine          Allocates all arrays from variable and equation classes using the global index.
        !   - stretchFunct                      function            A algebraic stretching function where the user defines if the stretch will be 2-ways or 1-way.
        !   - meshgrid                          subroutine          Creates a 2D ar 3D mesh of the input arrays. The ouput arrays will be in 2 and 3 dimensions.
        !   - update_equation_terms             subroutine          Updates all equation terms and the variable values of the previous time step.
        !   - csvread                           subroutine          Reads an 1D array from a user defined file.
        !   - csvwrite                          subroutine          Writes a 1D or 2D array into a file.
        !   - convert2Dto1D                     function            Converts a 2D array into a 1D array, using the global index.
        !   - convert1Dto2D                     function            Converts a 1D array into a 2D array, using the global index.
        !   - enter_last_result_dir             subroutine          A subroutine that changes the working directory into the last saved result directory.
        !   - obtain_variable                   subroutine          Obtains the value of a variable from a file.
        !   - calculate_residual_arrays         subroutine          Calculates an increasing size array with the resudual of a variable in every time step.
        !
        !
        !
        !
        !

        module subprograms_case
            implicit none


        contains

            subroutine allocate_variables()
                use Class_Advection
                use Class_CFL
                use Class_Scalar_Equation
                use Class_Scalar_Variable
                use Class_Vector_Equation
                use Class_Vector_Variable
                use global_variables
                use define_classes

                call V%allocate
                call V_old%allocate
                call p%allocate
                call p_old%allocate
                call H%allocate

                call NS_eqn%allocate
                call Pres_eqn%allocate

                call NS_adv_term%allocate
                call Co%allocate

                call V%allocate_BCs(domain, Vx_BC)
                call V%allocate_BCs(domain, Vy_BC)
                call V%allocate_BCs(domain, Vz_BC)
                call p%allocate_BCs(domain, p_BC)

                allocate(mat_inlet(x%NumberNodes*y%NumberNodes))

            end subroutine


            subroutine update_equation_terms()
                use global_variables
                use Class_Geometry
                use Class_Scalar_Variable
                use Class_Scalar_Equation
                use Class_Vector_Variable
                use Class_Vector_Equation
                use define_classes
                use subprograms_default

                ! Save previous time step variables
                V_old = V
                p_old = p

                ! Update advection terms
                call NS_adv_term%create_terms(V)

                aP_time = (Dx*Dy/Dt)*(1.0-k_boole+k_boole*Dz)

                ! Navier-Stokes x-axis terms
                NS_eqn%x%aT = -nu*aT_lapl + NS_adv_term%aT
                NS_eqn%x%aN = -nu*aN_lapl + NS_adv_term%aN
                NS_eqn%x%aE = -nu*aE_lapl + NS_adv_term%aE
                NS_eqn%x%aP = -nu*aP_lapl + NS_adv_term%aP + aP_time
                NS_eqn%x%aW = -nu*aW_lapl + NS_adv_term%aW
                NS_eqn%x%aS = -nu*aS_lapl + NS_adv_term%aS
                NS_eqn%x%aB = -nu*aB_lapl + NS_adv_term%aB
                NS_eqn%x%B  =  aP_time*V%x(nodes_P) 

                ! Navier-Stokes y-axis terms
                NS_eqn%y%aT = -nu*aT_lapl + NS_adv_term%aT
                NS_eqn%y%aN = -nu*aN_lapl + NS_adv_term%aN
                NS_eqn%y%aE = -nu*aE_lapl + NS_adv_term%aE
                NS_eqn%y%aP = -nu*aP_lapl + NS_adv_term%aP + aP_time
                NS_eqn%y%aW = -nu*aW_lapl + NS_adv_term%aW
                NS_eqn%y%aS = -nu*aS_lapl + NS_adv_term%aS
                NS_eqn%y%aB = -nu*aB_lapl + NS_adv_term%aB
                NS_eqn%y%B  =  aP_time*V%y(nodes_P) 

                ! Navier-Stokes z-axis terms
                if (flow_2D .eqv. .FALSE.)  then
                    NS_eqn%z%aT = -nu*aT_lapl + NS_adv_term%aT
                    NS_eqn%z%aN = -nu*aN_lapl + NS_adv_term%aN
                    NS_eqn%z%aE = -nu*aE_lapl + NS_adv_term%aE
                    NS_eqn%z%aP = -nu*aP_lapl + NS_adv_term%aP + aP_time
                    NS_eqn%z%aW = -nu*aW_lapl + NS_adv_term%aW
                    NS_eqn%z%aS = -nu*aS_lapl + NS_adv_term%aS
                    NS_eqn%z%aB = -nu*aB_lapl + NS_adv_term%aB
                    NS_eqn%z%B  =  aP_time*V%z(nodes_P)
                end if
            end subroutine
            
            subroutine create_A()
                use global_variables
                use define_classes

                A = -nu*aP_lapl + NS_adv_term%aP + aP_time
            end subroutine


            subroutine create_H()
                use global_variables
                use define_classes
                use subprograms_default

                ! Create H
                H%x(nodes_P) = (-k_boole*NS_eqn%x%aT*V%x(nodes_T)-NS_eqn%x%aN*V%x(nodes_N)-NS_eqn%x%aE*V%x(nodes_E)-NS_eqn%x%aW*V%x(nodes_W)-NS_eqn%x%aS*V%x(nodes_S)-k_boole*NS_eqn%x%aB*V%x(nodes_B) + aP_time*V_old%x(nodes_P) )/A
                H%y(nodes_P) = (-k_boole*NS_eqn%y%aT*V%y(nodes_T)-NS_eqn%y%aN*V%y(nodes_N)-NS_eqn%y%aE*V%y(nodes_E)-NS_eqn%y%aW*V%y(nodes_W)-NS_eqn%y%aS*V%y(nodes_S)-k_boole*NS_eqn%y%aB*V%y(nodes_B) + aP_time*V_old%y(nodes_P) )/A
                H%z(nodes_P) = k_boole*(-NS_eqn%z%aT*V%z(nodes_T)-NS_eqn%z%aN*V%z(nodes_N)-NS_eqn%z%aE*V%z(nodes_E)-NS_eqn%z%aW*V%z(nodes_W)-NS_eqn%z%aS*V%z(nodes_S)-        NS_eqn%z%aB*V%z(nodes_B) + aP_time*V_old%z(nodes_P) )/A
                call H%BC_Adjust_x(domain, Vx_BC)
                call H%BC_Adjust_y(domain, Vy_BC)
                call H%BC_Adjust_z(domain, Vz_BC)

                ! Update Poisson equation B matrix (from the system Ax=B)
                Pres_eqn%B  = A*div(H)

            end subroutine


            subroutine calculate_residual_arrays()
                use global_variables
                use define_classes
                real*8, allocatable, dimension(:)   :: temp_err_Vx, temp_err_Vy, temp_err_Vz, temp_err_p
                integer, allocatable, dimension(:)  :: temp_iterations

                allocate(temp_iterations(Iter_count), temp_err_Vx(Iter_count), temp_err_Vy(Iter_count), temp_err_Vz(Iter_count), temp_err_p(Iter_count))
                temp_iterations(1:Iter_count-1) = Iterations(1:Iter_count-1)
                temp_err_Vx(1:Iter_count-1) = Err_Vx(1:Iter_count-1)
                temp_err_Vy(1:Iter_count-1) = Err_Vy(1:Iter_count-1)
                temp_err_Vz(1:Iter_count-1) = Err_Vz(1:Iter_count-1)
                temp_err_p(1:Iter_count-1) = Err_p(1:Iter_count-1)
                temp_iterations(Iter_count) = Iter_count
                temp_err_Vx(Iter_count) = V%e(1)
                temp_err_Vy(Iter_count) = V%e(2)
                temp_err_Vz(Iter_count) = V%e(3)
                temp_err_p(Iter_count) = p%e
                call move_alloc(temp_iterations,Iterations)
                call move_alloc(temp_err_Vx,Err_Vx)
                call move_alloc(temp_err_Vy,Err_Vy)
                call move_alloc(temp_err_Vz,Err_Vz)
                call move_alloc(temp_err_p,Err_p)
            end subroutine


            subroutine calculate_nondim_numbers()
                use global_variables
                use Class_Vector_Variable
                use define_classes
                use subprograms_default

                ! Chalacteristic Velocity
                if (flow_2D .eqv. .FALSE.)  then
                    if (V%typeBC(inlet)=='custom_inlet')  then
                        V_char = integrate(convert1Dto2D(V%z(1:x%NumberNodes*y%NumberNodes),x%NumberNodes,y%NumberNodes), x%nodes, y%nodes)/(x%Length*y%Length)
                    else
                        V_char = Vz_BC(inlet)%characteristic_valueBC
                    end if
                else
                    V_char = Vx_BC(inlet)%characteristic_valueBC
                end if

                Re = V_char*(y%Length/2.0)/nu

            end subroutine

        end module
