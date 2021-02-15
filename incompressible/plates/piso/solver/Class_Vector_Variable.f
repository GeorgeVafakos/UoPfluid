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
        !   File name:          Class_Vector_Variable
        ! 
        !   Type:               module / class
        ! 
        !   Description:        In this file the vector variable and their boundary condition class are created.
        ! 
        ! 



        module Class_Vector_Variable
            use global_variables
            use Class_Geometry
            implicit none

            ! Create variable class type
            type, public :: Vector_Variable
                real*8, allocatable, dimension(:) :: x
                real*8, allocatable, dimension(:) :: y
                real*8, allocatable, dimension(:) :: z
                real*8  e(3)
                real*8  init(3)
                real*8  mean_div_error
                character(len=20) :: typeInit
                character(len=20), allocatable, dimension(:) :: typeBC
            contains
                procedure :: allocate => allocate_vector_variable
                procedure :: allocate_BCs => allocate_BC_value_at_boundaries
                procedure :: define_typeBC => allocate_typeBC_charmatrix_for_the_input_file
                procedure :: BC_Adjust_x => implement_boundary_conditions_x
                procedure :: BC_Adjust_y => implement_boundary_conditions_y
                procedure :: BC_Adjust_z => implement_boundary_conditions_z
                procedure :: calc_mean_div_error => calculate_mean_divergence_error
            end type

            ! Create variable boundary conditions class type
            type :: Vector_Variable_BC
                real*8 :: characteristic_valueBC, coeff_robin, custom_inlet_coeff
                real*8, allocatable, dimension(:)  :: value
                character(len=20) :: custom_inlet_file
                procedure(dirichlet_vect), pointer :: boundary_condition => null()
            end type


        contains

            subroutine allocate_vector_variable(Var)
                class (Vector_Variable) :: Var

                allocate(Var%x(TotalNodes), Var%y(TotalNodes), Var%z(TotalNodes))
            end subroutine


            subroutine allocate_BC_value_at_boundaries(Var, domain, Var_BC)
                class (Vector_Variable) :: Var
                class (duct) :: domain
                class (Vector_Variable_BC), dimension(domain%NumberBoundaries) :: Var_BC

                do i = 1, domain%NumberBoundaries
                    allocate(Var_BC(i)%value(domain%NumberNodesBoundary(i)))
                end do
            end subroutine


            subroutine allocate_typeBC_charmatrix_for_the_input_file(Var, domain)
                class (Vector_Variable) :: Var
                class (duct) :: domain

                allocate(Var%typeBC(domain%NumberBoundaries))
            end subroutine


            subroutine implement_boundary_conditions_x(Var, domain, Var_BC)
                class (Vector_Variable) :: Var
                class (duct) :: domain
                class (Vector_Variable_BC), dimension(domain%NumberBoundaries) :: Var_BC

                do i = 1, domain%NumberBoundaries
                    call Var_BC(i)%boundary_condition(Var%x,domain%boundary(i))
                end do
            end subroutine


            subroutine implement_boundary_conditions_y(Var, domain, Var_BC)
                class (Vector_Variable) :: Var
                class (duct) :: domain
                class (Vector_Variable_BC), dimension(domain%NumberBoundaries) :: Var_BC

                do i = 1, domain%NumberBoundaries
                    call Var_BC(i)%boundary_condition(Var%y,domain%boundary(i))
                end do
            end subroutine


            subroutine implement_boundary_conditions_z(Var, domain, Var_BC)
                class (Vector_Variable) :: Var
                class (duct) :: domain
                class (Vector_Variable_BC), dimension(domain%NumberBoundaries) :: Var_BC

                do i = 1, domain%NumberBoundaries
                    call Var_BC(i)%boundary_condition(Var%z,domain%boundary(i))
                end do
            end subroutine


            subroutine neumann_vect(Var_BC, var, bound)
                class (Vector_Variable_BC) :: Var_BC
                class (boundary) :: bound
                real*8, dimension(:) :: var

                var(bound%nodes_bound) = (bound%nodes_dist1*bound%nodes_dist2*(bound%nodes_dist1+bound%nodes_dist2)*Var_BC%value + ((bound%nodes_dist1+bound%nodes_dist2)**2)*var(bound%nodes_next1) - (bound%nodes_dist1**2)*var(bound%nodes_next2))/(bound%nodes_dist2**2+2*bound%nodes_dist1*bound%nodes_dist2)
            end subroutine


            subroutine dirichlet_vect(Var_BC, var, bound)
                class (Vector_Variable_BC) :: Var_BC
                class (boundary) :: bound
                real*8, dimension(:) :: var

                var(bound%nodes_bound) = Var_BC%value
            end subroutine


            subroutine robin_vect(Var_BC, var, bound)
                class (Vector_Variable_BC) :: Var_BC
                class (boundary) :: bound
                real*8, dimension(:) :: var

                var(bound%nodes_bound) = (bound%nodes_dist1*bound%nodes_dist2*(bound%nodes_dist1+bound%nodes_dist2)*Var_BC%value + ((bound%nodes_dist1+bound%nodes_dist2)**2)*var(bound%nodes_next1) - (Var_BC%coeff_robin*bound%nodes_dist1*bound%nodes_dist2*(bound%nodes_dist1+bound%nodes_dist2)+bound%nodes_dist1**2)*var(bound%nodes_next2))/(bound%nodes_dist2**2+2*bound%nodes_dist1*bound%nodes_dist2)
            end subroutine


            subroutine calculate_mean_divergence_error(Var)
                class (Vector_Variable) :: Var
                real*8 :: Div_Err(TotalCells)

                Div_Err = aE_divx*Var%x(nodes_E)+aP_divx*Var%x(nodes_P)+aW_divx*Var%x(nodes_W) + aN_divy*Var%y(nodes_N)+aP_divy*Var%y(nodes_P)+aS_divy*Var%y(nodes_S) + k_boole*(aT_divz*Var%z(nodes_T)+aP_divz*Var%z(nodes_P)+aB_divz*Var%z(nodes_B))
                Var%mean_div_error = sum(Div_Err)/size(Div_Err)
            end subroutine


        end module
