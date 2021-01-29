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
        !   File name:          Class_CFL
        !
        !   Type:               module / class
        !
        !   Description:        In this file the geometry and the grid is generated according to the input file parameters.
        !
        !


        module Class_CFL
            use global_variables
            use Class_Geometry
            use Class_Vector_Variable
            implicit none


            ! Create CFL class type
            type :: CFL_number
                real*8, allocatable, dimension(:) :: x
                real*8, allocatable, dimension(:) :: y
                real*8, allocatable, dimension(:) :: z
                real*8 :: x_max,  y_max,  z_max
                real*8 :: x_mean, y_mean, z_mean
                real*8 :: max, mean, fixed_value
                procedure(calculate_CFL), pointer :: CFL_condition => null()
            contains
                procedure :: allocate => allocate_CFL_number
            end type

            ! Declare CFL number for Navier-Stokes
            type (CFL_number) :: Co


        contains

            subroutine allocate_CFL_number(CFL)
                class (CFL_number) :: CFL

                allocate(CFL%x(TotalCells), CFL%y(TotalCells), CFL%z(TotalCells))
            end subroutine
            

            subroutine calculate_CFL(CFL, Var)
                class (CFL_number) :: CFL
                class (Vector_Variable) :: Var

                CFL%x= Var%x(nodes_P)*Dt/Dx
                CFL%x_max = maxval(CFL%x)
                CFL%x_mean = sum(CFL%x)/size(CFL%x)

                CFL%y= Var%y(nodes_P)*Dt/Dy
                CFL%y_max = maxval(CFL%y)
                CFL%y_mean = sum(CFL%y)/size(CFL%y)

                CFL%z= Var%z(nodes_P)*Dt/Dz
                CFL%z_max = maxval(CFL%z)
                CFL%z_mean = sum(CFL%z)/size(CFL%z)

                CFL%max  = max(CFL%x_max, CFL%y_max, k_boole*CFL%z_max)
                CFL%mean = (CFL%x_mean+CFL%y_mean+k_boole*CFL%z_mean)/(2+k_boole)
            end subroutine


            subroutine calculate_Adj_Dt(CFL, Var)
                class (CFL_number) :: CFL
                class (Vector_Variable) :: Var

                ! Dt = min(maxval(CFL%fixed_value*Dx/Var%x(nodes_P), mask=CFL%fixed_value*Dx/Var%x(nodes_P)==CFL%fixed_value*Dx/Var%x(nodes_P)), maxval(CFL%fixed_value*Dy/Var%y(nodes_P), mask=CFL%fixed_value*Dy/Var%y(nodes_P)==CFL%fixed_value*Dy/Var%y(nodes_P)))
                Dt = min(minval(CFL%fixed_value*Dx/abs(Var%x(nodes_P))), minval(CFL%fixed_value*Dy/abs(Var%y(nodes_P))))
                if (flow_2D .eqv. .FALSE.)    then
                    Dt = min(Dt, minval(CFL%fixed_value*Dz/abs(Var%z(nodes_P))))
                end if

                CFL%max = CFL%fixed_value
                CFL%mean = (sum(Var%x(nodes_P)*Dt/Dx)+sum(Var%y(nodes_P)*Dt/Dy)+k_boole*sum(Var%z(nodes_P)*Dt/Dz))/((2+k_boole)*TotalCells)
            end subroutine



        end module
