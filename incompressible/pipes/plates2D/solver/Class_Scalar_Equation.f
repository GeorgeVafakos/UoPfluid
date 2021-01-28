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
        !   File name:          mesh_generation
        !
        !   Type:               module / class
        !
        !   Description:        In this file the geometry and the grid is generated according to the input file parameters.
        !
        !

        module Class_Scalar_Equation
            use global_variables
            use Class_Geometry
            use Class_Scalar_Variable

            implicit none

            ! Create equation class type
            type :: Scalar_Equation
                real*8, allocatable, dimension(:) :: aT
                real*8, allocatable, dimension(:) :: aN
                real*8, allocatable, dimension(:) :: aE
                real*8, allocatable, dimension(:) :: aP
                real*8, allocatable, dimension(:) :: aW
                real*8, allocatable, dimension(:) :: aS
                real*8, allocatable, dimension(:) :: aB
                real*8, allocatable, dimension(:) :: B
                integer counter
                procedure(jacobi_scalar), pointer :: solve => null()
            contains
                procedure :: allocate => allocate_discrete_vector_equation
            end type
            
            
        contains
            
            subroutine allocate_discrete_vector_equation(Eqn)
                class (Scalar_Equation) :: Eqn

                allocate(Eqn%aN(TotalCells), Eqn%aE(TotalCells), Eqn%aP(TotalCells), Eqn%aW(TotalCells), Eqn%aS(TotalCells), Eqn%B(TotalCells), Eqn%aT(TotalCells), Eqn%aB(TotalCells))
            end subroutine


            subroutine jacobi_scalar(Eqn, domain, Var, Var_BC, max_count, tol, small_numbr)
                class (Scalar_Equation) :: Eqn
                class (Scalar_Variable) :: Var
                class (duct) :: domain
                class (Scalar_Variable_BC), dimension(domain%NumberBoundaries) :: Var_BC
                integer max_count
                real*8  :: Error(TotalNodes)
                real*8  e, tol, small_numbr

                var_old = Var%field
                Eqn%counter = 0
                do 
                    Eqn%counter = Eqn%counter + 1

                    ! Calculate the internal field
                    do ID = N%cell_beg, N%cell_end
                        Var%field(glob_index(ID)) = (1.0/Eqn%aP(ID))*(Eqn%B(ID) - Eqn%aT(ID)*Var_old(glob_index(ID)+k_boole*x%NumberNodes*y%NumberNodes) - Eqn%aN(ID)*Var_old(glob_index(ID)+x%NumberNodes) - Eqn%aE(ID)*Var_old(glob_index(ID)+1) - Eqn%aW(ID)*Var_old(glob_index(ID)-1)  - Eqn%aS(ID)*Var_old(glob_index(ID)-x%NumberNodes) - Eqn%aB(ID)*Var_old(glob_index(ID)-k_boole*x%NumberNodes*y%NumberNodes) )
                    end do

                    ! Adjust boundary conditions
                    call Var%BC_Adjust(domain, Var_BC)

                    ! Error calculation (mean value)
                    Error = abs(Var%field-Var_old)/max(abs(Var%field),small_numbr)
                    e = sum(Error)/size(Error)

                    ! Stop condition
                    if (e <= tol .OR. e/=e .OR. Eqn%counter==max_count) then
                        exit
                    else
                        Var_old = Var%field
                    end if
                end do
            end subroutine


            subroutine gauss_seidel_scalar(Eqn, domain, Var, Var_BC, max_count, tol, small_numbr)
                class (Scalar_Equation) :: Eqn
                class (Scalar_Variable) :: Var
                class (duct) :: domain
                class (Scalar_Variable_BC), dimension(domain%NumberBoundaries) :: Var_BC
                integer max_count
                real*8  :: Error(TotalNodes)
                real*8  e, tol, small_numbr

                var_old = Var%field
                Eqn%counter = 0
                do 
                    Eqn%counter = Eqn%counter + 1

                    ! Calculate the internal field
                    do ID = N%cell_beg, N%cell_end
                        Var%field(glob_index(ID)) = (1.0/Eqn%aP(ID))*(Eqn%B(ID) - Eqn%aT(ID)*Var_old(glob_index(ID)+k_boole*x%NumberNodes*y%NumberNodes) - Eqn%aN(ID)*Var_old(glob_index(ID)+x%NumberNodes) - Eqn%aE(ID)*Var_old(glob_index(ID)+1) - Eqn%aW(ID)*Var%field(glob_index(ID)-1)  - Eqn%aS(ID)*Var%field(glob_index(ID)-x%NumberNodes) - Eqn%aB(ID)*Var%field(glob_index(ID)-k_boole*x%NumberNodes*y%NumberNodes) )
                    end do

                    ! Adjust boundary conditions
                    call Var%BC_Adjust(domain, Var_BC)

                    ! Error calculation (mean value)
                    Error = abs(Var%field-Var_old)/max(abs(Var%field),small_numbr)
                    e = sum(Error)/size(Error)

                    ! Stop condition
                    if (e <= tol .OR. e/=e .OR. Eqn%counter==max_count) then
                        exit
                    else
                        Var_old = Var%field
                    end if
                end do
            end subroutine


            subroutine SOR_scalar(Eqn, domain, Var, Var_BC, max_count, tol, small_numbr)
                class (Scalar_Equation) :: Eqn
                class (Scalar_Variable) :: Var
                class (duct) :: domain
                class (Scalar_Variable_BC), dimension(domain%NumberBoundaries) :: Var_BC
                integer max_count
                real*8  :: Error(TotalNodes)
                real*8  e, tol, small_numbr

                var_old = Var%field
                Eqn%counter = 0
                do 
                    Eqn%counter = Eqn%counter + 1

                    ! Calculate the internal field
                    do ID = N%cell_beg, N%cell_end
                        Var%field(glob_index(ID)) = (1.0-r)*Var%field(glob_index(ID)) + r*(1.0/Eqn%aP(ID))*(Eqn%B(ID) - Eqn%aT(ID)*Var_old(glob_index(ID)+k_boole*x%NumberNodes*y%NumberNodes) - Eqn%aN(ID)*Var_old(glob_index(ID)+x%NumberNodes) - Eqn%aE(ID)*Var_old(glob_index(ID)+1) - Eqn%aW(ID)*Var%field(glob_index(ID)-1)  - Eqn%aS(ID)*Var%field(glob_index(ID)-x%NumberNodes) - Eqn%aB(ID)*Var%field(glob_index(ID)-k_boole*x%NumberNodes*y%NumberNodes) )
                    end do

                    ! Adjust boundary conditions
                    call Var%BC_Adjust(domain, Var_BC)

                    ! Error calculation (mean value)
                    Error = abs(Var%field-Var_old)/max(abs(Var%field),small_numbr)
                    e = sum(Error)/size(Error)

                    ! Stop condition
                    if (e <= tol .OR. e/=e .OR. Eqn%counter==max_count) then
                        exit
                    else
                        Var_old = Var%field
                    end if
                end do
            end subroutine

        end module