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

        module Class_Vector_Equation
            use global_variables
            use Class_Geometry
            use Class_Vector_Variable
            use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

            implicit none

            ! Create equation class type
            type :: component
                real*8, allocatable, dimension(:) :: aT
                real*8, allocatable, dimension(:) :: aN
                real*8, allocatable, dimension(:) :: aE
                real*8, allocatable, dimension(:) :: aP
                real*8, allocatable, dimension(:) :: aW
                real*8, allocatable, dimension(:) :: aS
                real*8, allocatable, dimension(:) :: aB
                real*8, allocatable, dimension(:) :: B
            end type

            type, public :: Vector_Equation
                type (component) :: x
                type (component) :: y
                type (component) :: z
                integer :: counter(3)
                procedure(jacobi_vector), pointer :: solve => null()
            contains
                procedure :: allocate => allocate_discrete_vector_equation
            end type
            
            
        contains
            
            subroutine allocate_discrete_vector_equation(Eqn)
                class (Vector_Equation) :: Eqn

                allocate(Eqn%x%aT(TotalCells), Eqn%x%aN(TotalCells), Eqn%x%aE(TotalCells), Eqn%x%aP(TotalCells), Eqn%x%aW(TotalCells), Eqn%x%aS(TotalCells), Eqn%x%aB(TotalCells), Eqn%x%B(TotalCells))
                allocate(Eqn%y%aT(TotalCells), Eqn%y%aN(TotalCells), Eqn%y%aE(TotalCells), Eqn%y%aP(TotalCells), Eqn%y%aW(TotalCells), Eqn%y%aS(TotalCells), Eqn%y%aB(TotalCells), Eqn%y%B(TotalCells))
                allocate(Eqn%z%aT(TotalCells), Eqn%z%aN(TotalCells), Eqn%z%aE(TotalCells), Eqn%z%aP(TotalCells), Eqn%z%aW(TotalCells), Eqn%z%aS(TotalCells), Eqn%z%aB(TotalCells), Eqn%z%B(TotalCells))
            end subroutine


            subroutine jacobi_vector(Eqn, domain, Var, VarX_BC, VarY_BC, VarZ_BC, max_count, tol_x, tol_y, tol_z, small_numbr)
                class (Vector_Equation) :: Eqn
                class (Vector_Variable) :: Var
                class (duct) :: domain
                class (Vector_Variable_BC), dimension(domain%NumberBoundaries) :: VarX_BC, VarY_BC, VarZ_BC
                real*8 tol_x, tol_y, tol_z
                integer max_count
                real*8  :: Error(TotalNodes)
                real*8  e, small_numbr
                
                ! x-component of equation
                var_old = Var%x
                Eqn%counter(1) = 0
                do 
                    Eqn%counter(1) = Eqn%counter(1) + 1

                    ! Calculate the internal field
                    do ID = N%cell_beg, N%cell_end
                        Var%x(glob_index(ID)) = (1.0/Eqn%x%aP(ID))*(Eqn%x%B(ID)  - Eqn%x%aT(ID)*Var_old(glob_index(ID)+k_boole*x%NumberNodes*y%NumberNodes) - Eqn%x%aN(ID)*Var_old(glob_index(ID)+x%NumberNodes) - Eqn%x%aE(ID)*Var_old(glob_index(ID)+1) - Eqn%x%aW(ID)*Var_old(glob_index(ID)-1) - Eqn%x%aS(ID)*Var_old(glob_index(ID)-x%NumberNodes) - Eqn%x%aB(ID)*Var_old(glob_index(ID)-k_boole*x%NumberNodes*y%NumberNodes))
                    end do

                    ! Adjust boundary conditions
                    call Var%BC_Adjust_x(domain, VarX_BC)

                    ! Error calculation (mean value)
                    Error = abs(Var%x-Var_old)/max(abs(Var%x),small_numbr)
                    e = sum(Error)/size(Error)

                    ! Stop condition
                    if (e <= tol_x .OR. e/=e .OR. Eqn%counter(1)==max_count) then
                        exit
                    else
                        Var_old = Var%x
                    end if
                end do

                ! y-component of equation
                Var_old = Var%y
                Eqn%counter(2) = 0
                do 
                    Eqn%counter(2) = Eqn%counter(2) + 1

                    ! Calculate the internal field
                    do ID = N%cell_beg, N%cell_end
                        Var%y(glob_index(ID)) = (1.0/Eqn%y%aP(ID))*(Eqn%y%B(ID)  - Eqn%y%aT(ID)*Var_old(glob_index(ID)+k_boole*x%NumberNodes*y%NumberNodes) - Eqn%y%aN(ID)*Var_old(glob_index(ID)+x%NumberNodes) - Eqn%y%aE(ID)*Var_old(glob_index(ID)+1) - Eqn%y%aW(ID)*Var_old(glob_index(ID)-1) - Eqn%y%aS(ID)*Var_old(glob_index(ID)-x%NumberNodes) - Eqn%y%aB(ID)*Var_old(glob_index(ID)-k_boole*x%NumberNodes*y%NumberNodes))
                    end do
                    ! Adjust boundary conditions
                    call Var%BC_Adjust_y(domain, VarY_BC)

                    ! Error calculation (mean value)
                    Error = abs(Var%y-Var_old)/max(abs(Var%y),small_numbr)
                    e = sum(Error)/size(Error)

                    ! Stop condition
                    if (e <= tol_y .OR. e/=e .OR. Eqn%counter(2)==max_count) then
                        exit
                    else
                        Var_old = Var%y
                    end if
                end do

                ! z-component of equation
                Var_old = Var%z
                Eqn%counter(3) = 0
                do while (flow_2D .eqv. .FALSE.)
                    Eqn%counter(3) = Eqn%counter(3) + 1

                    ! Calculate the internal field
                    do ID = N%cell_beg, N%cell_end
                        Var%z(glob_index(ID)) = (1.0/Eqn%z%aP(ID))*(Eqn%z%B(ID)  - Eqn%z%aT(ID)*Var_old(glob_index(ID)+k_boole*x%NumberNodes*y%NumberNodes) - Eqn%z%aN(ID)*Var_old(glob_index(ID)+x%NumberNodes) - Eqn%z%aE(ID)*Var_old(glob_index(ID)+1) - Eqn%z%aW(ID)*Var_old(glob_index(ID)-1) - Eqn%z%aS(ID)*Var_old(glob_index(ID)-x%NumberNodes) - Eqn%z%aB(ID)*Var_old(glob_index(ID)-k_boole*x%NumberNodes*y%NumberNodes))
                    end do
                    ! Adjust boundary conditions
                    call Var%BC_Adjust_z(domain, VarZ_BC)

                    ! Error calculation (mean value)
                    Error = abs(Var%z-Var_old)/max(abs(Var%z),small_numbr)
                    e = sum(Error)/size(Error)

                    ! Stop condition
                    if (e <= tol_z .OR. e/=e .OR. Eqn%counter(3)==max_count) then
                        exit
                    else
                        Var_old = Var%z
                    end if
                end do

            end subroutine


            subroutine gauss_seidel_vector(Eqn, domain, Var, VarX_BC, VarY_BC, VarZ_BC, max_count, tol_x, tol_y, tol_z, small_numbr)
                class (Vector_Equation) :: Eqn
                class (Vector_Variable) :: Var
                class (duct) :: domain
                class (Vector_Variable_BC), dimension(domain%NumberBoundaries) :: VarX_BC, VarY_BC, VarZ_BC
                real*8 tol_x, tol_y, tol_z
                integer max_count
                real*8  :: Error(TotalNodes)
                real*8  e, small_numbr
                
                ! x-component of equation
                var_old = Var%x
                Eqn%counter(1) = 0
                do 
                    Eqn%counter(1) = Eqn%counter(1) + 1

                    ! Calculate the internal field
                    do ID = N%cell_beg, N%cell_end
                        Var%x(glob_index(ID)) = (1.0/Eqn%x%aP(ID))*(Eqn%x%B(ID)  - Eqn%x%aT(ID)*Var_old(glob_index(ID)+k_boole*x%NumberNodes*y%NumberNodes) - Eqn%x%aN(ID)*Var_old(glob_index(ID)+x%NumberNodes) - Eqn%x%aE(ID)*Var_old(glob_index(ID)+1) - Eqn%x%aW(ID)*Var%x(glob_index(ID)-1) - Eqn%x%aS(ID)*Var%x(glob_index(ID)-x%NumberNodes) - Eqn%x%aB(ID)*Var%x(glob_index(ID)-k_boole*x%NumberNodes*y%NumberNodes))
                    end do
                    ! Adjust boundary conditions
                    call Var%BC_Adjust_x(domain, VarX_BC)

                    ! Error calculation (mean value)
                    Error = abs(Var%x-Var_old)/max(abs(Var%x),small_numbr)
                    e = sum(Error)/size(Error)

                    ! Stop condition
                    if (e <= tol_x .OR. e/=e .OR. Eqn%counter(1)==max_count) then
                        exit
                    else
                        Var_old = Var%x
                    end if
                end do

                ! y-component of equation
                Var_old = Var%y
                Eqn%counter(2) = 0
                do 
                    Eqn%counter(2) = Eqn%counter(2) + 1

                    ! Calculate the internal field
                    do ID = N%cell_beg, N%cell_end
                        Var%y(glob_index(ID)) = (1.0/Eqn%y%aP(ID))*(Eqn%y%B(ID)  - Eqn%y%aT(ID)*Var_old(glob_index(ID)+k_boole*x%NumberNodes*y%NumberNodes) - Eqn%y%aN(ID)*Var_old(glob_index(ID)+x%NumberNodes) - Eqn%y%aE(ID)*Var_old(glob_index(ID)+1) - Eqn%y%aW(ID)*Var%y(glob_index(ID)-1) - Eqn%y%aS(ID)*Var%y(glob_index(ID)-x%NumberNodes) - Eqn%y%aB(ID)*Var%y(glob_index(ID)-k_boole*x%NumberNodes*y%NumberNodes))
                    end do
                    ! Adjust boundary conditions
                    call Var%BC_Adjust_y(domain, VarY_BC)

                    ! Error calculation (mean value)
                    Error = abs(Var%y-Var_old)/max(abs(Var%y),small_numbr)
                    e = sum(Error)/size(Error)

                    ! Stop condition
                    if (e <= tol_y .OR. e/=e .OR. Eqn%counter(2)==max_count) then
                        exit
                    else
                        Var_old = Var%y
                    end if
                end do

                ! z-component of equation
                Var_old = Var%z
                Eqn%counter(3) = 0
                do while (flow_2D .eqv. .FALSE.)
                    Eqn%counter(3) = Eqn%counter(3) + 1

                    ! Calculate the internal field
                    do ID = N%cell_beg, N%cell_end
                        Var%z(glob_index(ID)) = (1.0/Eqn%z%aP(ID))*(Eqn%z%B(ID)  - Eqn%z%aT(ID)*Var_old(glob_index(ID)+k_boole*x%NumberNodes*y%NumberNodes) - Eqn%z%aN(ID)*Var_old(glob_index(ID)+x%NumberNodes) - Eqn%z%aE(ID)*Var_old(glob_index(ID)+1) - Eqn%z%aW(ID)*Var%z(glob_index(ID)-1) - Eqn%z%aS(ID)*Var%z(glob_index(ID)-x%NumberNodes) - Eqn%z%aB(ID)*Var%z(glob_index(ID)-k_boole*x%NumberNodes*y%NumberNodes))
                    end do
                    ! Adjust boundary conditions
                    call Var%BC_Adjust_z(domain, VarZ_BC)

                    ! Error calculation (mean value)
                    Error = abs(Var%z-Var_old)/max(abs(Var%z),small_numbr)
                    e = sum(Error)/size(Error)

                    ! Stop condition
                    if (e <= tol_z .OR. e/=e .OR. Eqn%counter(3)==max_count) then
                        exit
                    else
                        Var_old = Var%z
                    end if
                end do

            end subroutine


            subroutine SOR_vector(Eqn, domain, Var, VarX_BC, VarY_BC, VarZ_BC, max_count, tol_x, tol_y, tol_z, small_numbr)
                class (Vector_Equation) :: Eqn
                class (Vector_Variable) :: Var
                class (duct) :: domain
                class (Vector_Variable_BC), dimension(domain%NumberBoundaries) :: VarX_BC, VarY_BC, VarZ_BC
                real*8 tol_x, tol_y, tol_z
                integer max_count
                real*8  :: Error(TotalNodes)
                real*8  e, small_numbr
                
                ! x-component of equation
                var_old = Var%x
                Eqn%counter(1) = 0
                do 
                    Eqn%counter(1) = Eqn%counter(1) + 1

                    ! Calculate the internal field
                    do ID = N%cell_beg, N%cell_end
                        Var%x(glob_index(ID)) = (1.0-r)*Var%x(glob_index(ID)) + r*(1.0/Eqn%x%aP(ID))*(Eqn%x%B(ID)  - Eqn%x%aT(ID)*Var_old(glob_index(ID)+k_boole*x%NumberNodes*y%NumberNodes) - Eqn%x%aN(ID)*Var_old(glob_index(ID)+x%NumberNodes) - Eqn%x%aE(ID)*Var_old(glob_index(ID)+1) - Eqn%x%aW(ID)*Var%x(glob_index(ID)-1) - Eqn%x%aS(ID)*Var%x(glob_index(ID)-x%NumberNodes) - Eqn%x%aB(ID)*Var%x(glob_index(ID)-k_boole*x%NumberNodes*y%NumberNodes))
                    end do
                    ! Adjust boundary conditions
                    call Var%BC_Adjust_x(domain, VarX_BC)

                    ! Error calculation (mean value)
                    Error = abs(Var%x-Var_old)/max(abs(Var%x),small_numbr)
                    e = sum(Error)/size(Error)

                    ! Stop condition
                    if (e <= tol_x .OR. e/=e .OR. Eqn%counter(1)==max_count) then
                        exit
                    else
                        Var_old = Var%x
                    end if
                end do

                ! y-component of equation
                Var_old = Var%y
                Eqn%counter(2) = 0
                do 
                    Eqn%counter(2) = Eqn%counter(2) + 1

                    ! Calculate the internal field
                    do ID = N%cell_beg, N%cell_end
                        Var%y(glob_index(ID)) = (1.0-r)*Var%y(glob_index(ID)) + r*(1.0/Eqn%y%aP(ID))*(Eqn%y%B(ID)  - Eqn%y%aT(ID)*Var_old(glob_index(ID)+k_boole*x%NumberNodes*y%NumberNodes) - Eqn%y%aN(ID)*Var_old(glob_index(ID)+x%NumberNodes) - Eqn%y%aE(ID)*Var_old(glob_index(ID)+1) - Eqn%y%aW(ID)*Var%y(glob_index(ID)-1) - Eqn%y%aS(ID)*Var%y(glob_index(ID)-x%NumberNodes) - Eqn%y%aB(ID)*Var%y(glob_index(ID)-k_boole*x%NumberNodes*y%NumberNodes))
                    end do
                    ! Adjust boundary conditions
                    call Var%BC_Adjust_y(domain, VarY_BC)

                    ! Error calculation (mean value)
                    Error = abs(Var%y-Var_old)/max(abs(Var%y),small_numbr)
                    e = sum(Error)/size(Error)

                    ! Stop condition
                    if (e <= tol_y .OR. e/=e .OR. Eqn%counter(2)==max_count) then
                        exit
                    else
                        Var_old = Var%y
                    end if
                end do

                ! z-component of equation
                Var_old = Var%z
                Eqn%counter(3) = 0
                do while (flow_2D .eqv. .FALSE.)
                    Eqn%counter(3) = Eqn%counter(3) + 1

                    ! Calculate the internal field
                    do ID = N%cell_beg, N%cell_end
                        Var%z(glob_index(ID)) = (1.0-r)*Var%z(glob_index(ID)) + r*(1.0/Eqn%z%aP(ID))*(Eqn%z%B(ID)  - Eqn%z%aT(ID)*Var_old(glob_index(ID)+k_boole*x%NumberNodes*y%NumberNodes) - Eqn%z%aN(ID)*Var_old(glob_index(ID)+x%NumberNodes) - Eqn%z%aE(ID)*Var_old(glob_index(ID)+1) - Eqn%z%aW(ID)*Var%z(glob_index(ID)-1) - Eqn%z%aS(ID)*Var%z(glob_index(ID)-x%NumberNodes) - Eqn%z%aB(ID)*Var%z(glob_index(ID)-k_boole*x%NumberNodes*y%NumberNodes))
                    end do
                    ! Adjust boundary conditions
                    call Var%BC_Adjust_z(domain, VarZ_BC)

                    ! Error calculation (mean value)
                    Error = abs(Var%z-Var_old)/max(abs(Var%z),small_numbr)
                    e = sum(Error)/size(Error)

                    ! Stop condition
                    if (e <= tol_z .OR. e/=e .OR. Eqn%counter(3)==max_count) then
                        exit
                    else
                        Var_old = Var%z
                    end if
                end do

            end subroutine

        end module