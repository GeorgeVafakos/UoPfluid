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
                call B%allocate
                call Bi%allocate
                call Bo%allocate
                call dB%allocate
                call G%allocate
                call Je%allocate
                call T%allocate
                call T_old%allocate

                call NS_eqn%allocate
                call Pres_eqn%allocate
                call Ind_eqn%allocate
                call dB_eqn%allocate
                call Heat_eqn%allocate

                call NS_adv_term%allocate
                call NS_max_tens%allocate
                call Ind_adv_term%allocate
                call Ind_Bi_adv%allocate
                call Ind_Bo_adv%allocate
                call Heat_conv_term%allocate
                call Co%allocate

                call V%allocate_BCs(domain, Vx_BC)
                call V%allocate_BCs(domain, Vy_BC)
                call V%allocate_BCs(domain, Vz_BC)
                call p%allocate_BCs(domain, p_BC)
                call Bi%allocate_BCs(domain, Bix_BC)
                call Bi%allocate_BCs(domain, Biy_BC)
                call Bi%allocate_BCs(domain, Biz_BC)
                call dB%allocate_BCs(domain, dBx_BC)
                call dB%allocate_BCs(domain, dBy_BC)
                call dB%allocate_BCs(domain, dBz_BC)
                call T%allocate_BCs(domain, T_BC)

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
                Bi_old = Bi
                T_old = T

                ! Update advection terms
                call NS_adv_term%create_terms(V)
                call Ind_adv_term%create_terms(V)
                call NS_max_tens%create_terms(B)

                aP_time = (Dx*Dy/Dt)*(1.0-k_boole+k_boole*Dz)

                ! Navier-Stokes x-axis terms
                NS_eqn%x%aT = -nu*aT_lapl + NS_adv_term%aT
                NS_eqn%x%aN = -nu*aN_lapl + NS_adv_term%aN
                NS_eqn%x%aE = -nu*aE_lapl + NS_adv_term%aE
                NS_eqn%x%aP = -nu*aP_lapl + NS_adv_term%aP + aP_time
                NS_eqn%x%aW = -nu*aW_lapl + NS_adv_term%aW
                NS_eqn%x%aS = -nu*aS_lapl + NS_adv_term%aS
                NS_eqn%x%aB = -nu*aB_lapl + NS_adv_term%aB
                NS_eqn%x%B  =  aP_time*V%x(nodes_P) - (1.0/rho)*gradx(p%field) + (1.0/rho)*crossx(Je,B)

                ! Navier-Stokes y-axis terms
                NS_eqn%y%aT = -nu*aT_lapl + NS_adv_term%aT
                NS_eqn%y%aN = -nu*aN_lapl + NS_adv_term%aN
                NS_eqn%y%aE = -nu*aE_lapl + NS_adv_term%aE
                NS_eqn%y%aP = -nu*aP_lapl + NS_adv_term%aP + aP_time
                NS_eqn%y%aW = -nu*aW_lapl + NS_adv_term%aW
                NS_eqn%y%aS = -nu*aS_lapl + NS_adv_term%aS
                NS_eqn%y%aB = -nu*aB_lapl + NS_adv_term%aB
                NS_eqn%y%B  =  aP_time*V%y(nodes_P) - (1.0/rho)*grady(p%field) + (1.0/rho)*crossy(Je,B)

                ! Navier-Stokes z-axis terms
                if (flow_2D .eqv. .FALSE.)  then
                    NS_eqn%z%aT = -nu*aT_lapl + NS_adv_term%aT
                    NS_eqn%z%aN = -nu*aN_lapl + NS_adv_term%aN
                    NS_eqn%z%aE = -nu*aE_lapl + NS_adv_term%aE
                    NS_eqn%z%aP = -nu*aP_lapl + NS_adv_term%aP + aP_time
                    NS_eqn%z%aW = -nu*aW_lapl + NS_adv_term%aW
                    NS_eqn%z%aS = -nu*aS_lapl + NS_adv_term%aS
                    NS_eqn%z%aB = -nu*aB_lapl + NS_adv_term%aB
                    NS_eqn%z%B  =  aP_time*V%z(nodes_P) - (1.0/rho)*gradz(p%field) + (1.0/rho)*crossz(Je,B)
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
                H%x(nodes_P) = (-k_boole*NS_eqn%x%aT*V%x(nodes_T)-NS_eqn%x%aN*V%x(nodes_N)-NS_eqn%x%aE*V%x(nodes_E)-NS_eqn%x%aW*V%x(nodes_W)-NS_eqn%x%aS*V%x(nodes_S)-k_boole*NS_eqn%x%aB*V%x(nodes_B) + aP_time*V_old%x(nodes_P)  +  (1.0/rho)*crossx(Je,B))/A
                H%y(nodes_P) = (-k_boole*NS_eqn%y%aT*V%y(nodes_T)-NS_eqn%y%aN*V%y(nodes_N)-NS_eqn%y%aE*V%y(nodes_E)-NS_eqn%y%aW*V%y(nodes_W)-NS_eqn%y%aS*V%y(nodes_S)-k_boole*NS_eqn%y%aB*V%y(nodes_B) + aP_time*V_old%y(nodes_P)  +  (1.0/rho)*crossy(Je,B))/A
                H%z(nodes_P) = k_boole*(-NS_eqn%z%aT*V%z(nodes_T)-NS_eqn%z%aN*V%z(nodes_N)-NS_eqn%z%aE*V%z(nodes_E)-NS_eqn%z%aW*V%z(nodes_W)-NS_eqn%z%aS*V%z(nodes_S)-        NS_eqn%z%aB*V%z(nodes_B) + aP_time*V_old%z(nodes_P)  +  (1.0/rho)*crossz(Je,B))/A
                call H%BC_Adjust_x(domain, Vx_BC)
                call H%BC_Adjust_y(domain, Vy_BC)
                call H%BC_Adjust_z(domain, Vz_BC)

                ! Update Poisson equation B matrix (from the system Ax=B)
                Pres_eqn%B  = A*div(H)

            end subroutine

            subroutine update_induction_terms()
                use global_variables
                use Class_Geometry
                use Class_Scalar_Variable
                use Class_Scalar_Equation
                use Class_Vector_Variable
                use Class_Vector_Equation
                use define_classes
                use subprograms_default

                ! Advection terms
                call Ind_adv_term%create_terms(V)
                call Ind_Bi_adv%create_terms(Bi)

                ! Induction equation x-axis terms
                Ind_eqn%x%aT = -(1.0/(mu*sigma))*aT_lapl + Ind_adv_term%aT
                Ind_eqn%x%aN = -(1.0/(mu*sigma))*aN_lapl + Ind_adv_term%aN
                Ind_eqn%x%aE = -(1.0/(mu*sigma))*aE_lapl + Ind_adv_term%aE
                Ind_eqn%x%aP = -(1.0/(mu*sigma))*aP_lapl + Ind_adv_term%aP + aP_time
                Ind_eqn%x%aW = -(1.0/(mu*sigma))*aW_lapl + Ind_adv_term%aW
                Ind_eqn%x%aS = -(1.0/(mu*sigma))*aS_lapl + Ind_adv_term%aS
                Ind_eqn%x%aB = -(1.0/(mu*sigma))*aB_lapl + Ind_adv_term%aB
                Ind_eqn%x%B  = aP_time*Bi%x(nodes_P) + flux(Ind_Bi_adv,V%x) + flux(Ind_Bo_adv,V%x) - flux(Ind_adv_term,Bo%x)

                ! Induction equation y-axis terms
                Ind_eqn%y%aT = -(1.0/(mu*sigma))*aT_lapl + Ind_adv_term%aT
                Ind_eqn%y%aN = -(1.0/(mu*sigma))*aN_lapl + Ind_adv_term%aN
                Ind_eqn%y%aE = -(1.0/(mu*sigma))*aE_lapl + Ind_adv_term%aE
                Ind_eqn%y%aP = -(1.0/(mu*sigma))*aP_lapl + Ind_adv_term%aP + aP_time
                Ind_eqn%y%aW = -(1.0/(mu*sigma))*aW_lapl + Ind_adv_term%aW
                Ind_eqn%y%aS = -(1.0/(mu*sigma))*aS_lapl + Ind_adv_term%aS
                Ind_eqn%y%aB = -(1.0/(mu*sigma))*aB_lapl + Ind_adv_term%aB
                Ind_eqn%y%B  = aP_time*Bi%y(nodes_P) + flux(Ind_Bi_adv,V%y) + flux(Ind_Bo_adv,V%y) - flux(Ind_adv_term,Bo%y)

                ! Induction equation z-axis terms
                if (flow_2D .eqv. .FALSE.)  then
                    Ind_eqn%z%aT = -(1.0/(mu*sigma))*aT_lapl + Ind_adv_term%aT
                    Ind_eqn%z%aN = -(1.0/(mu*sigma))*aN_lapl + Ind_adv_term%aN
                    Ind_eqn%z%aE = -(1.0/(mu*sigma))*aE_lapl + Ind_adv_term%aE
                    Ind_eqn%z%aP = -(1.0/(mu*sigma))*aP_lapl + Ind_adv_term%aP + aP_time
                    Ind_eqn%z%aW = -(1.0/(mu*sigma))*aW_lapl + Ind_adv_term%aW
                    Ind_eqn%z%aS = -(1.0/(mu*sigma))*aS_lapl + Ind_adv_term%aS
                    Ind_eqn%z%aB = -(1.0/(mu*sigma))*aB_lapl + Ind_adv_term%aB
                    Ind_eqn%z%B  = aP_time*Bi%z(nodes_P) + flux(Ind_Bi_adv,V%z) + flux(Ind_Bo_adv,V%z) - flux(Ind_adv_term,Bo%z)
                end if

            end subroutine


            subroutine create_G()
                use global_variables
                use define_classes
                use subprograms_default

                ! Calculate G variable and extrapolate to the boundaries
                G%field(nodes_P) = - div(Bi)
                call extrapolate(G%field, domain)

                ! Update Poisson equation B matrix (from the system Ax=B)
                dB_Eqn%x%B = gradx(G%field)
                dB_Eqn%y%B = grady(G%field)
                dB_Eqn%z%B = gradz(G%field)

            end subroutine


            subroutine update_heat_terms()
                use global_variables
                use Class_Geometry
                use Class_Scalar_Variable
                use Class_Scalar_Equation
                use Class_Vector_Variable
                use Class_Vector_Equation
                use define_classes
                use subprograms_default

                call Heat_conv_term%create_terms(V)

                Heat_eqn%aT = -alpha*aT_lapl + Heat_conv_term%aT
                Heat_eqn%aN = -alpha*aN_lapl + Heat_conv_term%aN
                Heat_eqn%aE = -alpha*aE_lapl + Heat_conv_term%aE
                Heat_eqn%aP = -alpha*aP_lapl + Heat_conv_term%aP + aP_time
                Heat_eqn%aW = -alpha*aW_lapl + Heat_conv_term%aW
                Heat_eqn%aS = -alpha*aS_lapl + Heat_conv_term%aS
                Heat_eqn%aB = -alpha*aB_lapl + Heat_conv_term%aB
                Heat_eqn%B  = aP_time*T%field(nodes_P)
            end subroutine


            subroutine calculate_residual_arrays()
                use global_variables
                use define_classes
                real*8, allocatable, dimension(:)   :: temp_err_Vx, temp_err_Vy, temp_err_Vz, temp_err_p, temp_err_Bix, temp_err_Biy, temp_err_Biz, temp_err_T
                integer, allocatable, dimension(:)  :: temp_iterations

                allocate(temp_iterations(Iter_count), temp_err_Vx(Iter_count), temp_err_Vy(Iter_count), temp_err_Vz(Iter_count), temp_err_p(Iter_count), temp_err_Bix(Iter_count), temp_err_Biy(Iter_count), temp_err_Biz(Iter_count), temp_err_T(Iter_count))
                temp_iterations(1:Iter_count-1) = Iterations(1:Iter_count-1)
                temp_err_Vx(1:Iter_count-1) = Err_Vx(1:Iter_count-1)
                temp_err_Vy(1:Iter_count-1) = Err_Vy(1:Iter_count-1)
                temp_err_Vz(1:Iter_count-1) = Err_Vz(1:Iter_count-1)
                temp_err_p(1:Iter_count-1) = Err_p(1:Iter_count-1)
                temp_err_Bix(1:Iter_count-1) = Err_Bix(1:Iter_count-1)
                temp_err_Biy(1:Iter_count-1) = Err_Biy(1:Iter_count-1)
                temp_err_Biz(1:Iter_count-1) = Err_Biz(1:Iter_count-1)
                temp_err_T(1:Iter_count-1) = Err_T(1:Iter_count-1)
                temp_iterations(Iter_count) = Iter_count
                temp_err_Vx(Iter_count) = V%e(1)
                temp_err_Vy(Iter_count) = V%e(2)
                temp_err_Vz(Iter_count) = V%e(3)
                temp_err_p(Iter_count) = p%e
                temp_err_Bix(Iter_count) = Bi%e(1)
                temp_err_Biy(Iter_count) = Bi%e(2)
                temp_err_Biz(Iter_count) = Bi%e(3)
                temp_err_T(Iter_count) = T%e
                call move_alloc(temp_iterations,Iterations)
                call move_alloc(temp_err_Vx,Err_Vx)
                call move_alloc(temp_err_Vy,Err_Vy)
                call move_alloc(temp_err_Vz,Err_Vz)
                call move_alloc(temp_err_p,Err_p)
                call move_alloc(temp_err_Bix,Err_Bix)
                call move_alloc(temp_err_Biy,Err_Biy)
                call move_alloc(temp_err_Biz,Err_Biz)
                call move_alloc(temp_err_T,Err_T)
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

                ! Characteristic Magnetic Field
                if (Bo%typeInit=='smooth_step')  then
                    B_char = Bo_mag
                else
                    B_char = Bo%init(2)
                end if

                Re = V_char*(y%Length/2.0)/nu
                Ha = B_char*(y%Length/2.0)*sqrt(sigma/(rho*nu))
                Rm = mu*sigma*V_char*(y%Length/2.0)

            end subroutine

        end module
