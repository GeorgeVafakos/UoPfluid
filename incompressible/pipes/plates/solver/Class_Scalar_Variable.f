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
        !   File name:          Class_Scalar_Variable
        ! 
        !   Type:               source
        ! 
        !   Description:        In this file the scalar variable and their boundary condition class are created.
        ! 
        ! 



        module Class_Scalar_Variable
            use global_variables
            use Class_Geometry
            implicit none

            ! Create variable class type
            type, public :: Scalar_Variable
                real*8, allocatable, dimension(:) :: field
                real*8  e
                real*8  init
                character(len=20) :: typeInit
                character(len=20), allocatable, dimension(:) :: typeBC
            contains
                procedure :: allocate => allocate_scalar_variable
                procedure :: allocate_BCs => allocate_BC_value_at_boundaries
                procedure :: define_typeBC => allocate_typeBC_charmatrix_for_the_input_file
                procedure :: BC_Adjust => implement_boundary_conditions
            end type

            type :: Scalar_Variable_BC
                real*8 :: characteristic_valueBC, coeff_robin, custom_inlet_coeff
                real*8, allocatable, dimension(:) :: value
                character(len=20) :: custom_inlet_file
                procedure(dirichlet_scalar), pointer :: boundary_condition => null()
            end type


        contains

            subroutine allocate_scalar_variable(Var)
                class (Scalar_Variable) :: Var
                allocate(Var%field(TotalNodes))
            end subroutine


            subroutine allocate_BC_value_at_boundaries(Var, domain, Var_BC)
                class (Scalar_Variable) :: Var
                class (duct) :: domain
                class (Scalar_Variable_BC), dimension(domain%NumberBoundaries) :: Var_BC

                do i = 1, domain%NumberBoundaries
                    allocate(Var_BC(i)%value(domain%NumberNodesBoundary(i)))
                end do
            end subroutine


            subroutine allocate_typeBC_charmatrix_for_the_input_file(Var, domain)
                class (Scalar_Variable) :: Var
                class (duct) :: domain

                allocate(Var%typeBC(domain%NumberBoundaries))
            end subroutine


            subroutine implement_boundary_conditions(Var, domain, Var_BC)
                class (Scalar_Variable) :: Var
                class (duct) :: domain
                class (Scalar_Variable_BC), dimension(domain%NumberBoundaries) :: Var_BC

                do i = 1, domain%NumberBoundaries
                    call Var_BC(i)%boundary_condition(Var%field,domain%boundary(i))
                end do
                
            end subroutine


            subroutine neumann_scalar(Var_BC, var, bound)
                class (Scalar_Variable_BC) :: Var_BC
                class (boundary) :: bound
                real*8, dimension(:) :: var

                var(bound%nodes_bound) = (bound%nodes_dist1*bound%nodes_dist2*(bound%nodes_dist1+bound%nodes_dist2)*Var_BC%value + ((bound%nodes_dist1+bound%nodes_dist2)**2)*var(bound%nodes_next1) - (bound%nodes_dist1**2)*var(bound%nodes_next2))/(bound%nodes_dist2**2+2*bound%nodes_dist1*bound%nodes_dist2)

            end subroutine


            subroutine robin_scalar(Var_BC, var, bound)
                class (Scalar_Variable_BC) :: Var_BC
                class (boundary) :: bound
                real*8, dimension(:) :: var

                var(bound%nodes_bound) = (bound%nodes_dist1*bound%nodes_dist2*(bound%nodes_dist1+bound%nodes_dist2)*Var_BC%value + ((bound%nodes_dist1+bound%nodes_dist2)**2)*var(bound%nodes_next1) - (Var_BC%coeff_robin*bound%nodes_dist1*bound%nodes_dist2*(bound%nodes_dist1+bound%nodes_dist2)+bound%nodes_dist1**2)*var(bound%nodes_next2))/(bound%nodes_dist2**2+2*bound%nodes_dist1*bound%nodes_dist2)
            end subroutine


            subroutine dirichlet_scalar(Var_BC, var, bound)
                class (Scalar_Variable_BC) :: Var_BC
                class (boundary) :: bound
                real*8, dimension(:) :: var
                
                var(bound%nodes_bound) = Var_BC%value

            end subroutine


        end module
