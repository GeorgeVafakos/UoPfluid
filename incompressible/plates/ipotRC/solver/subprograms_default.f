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
        !   File name:          subprograms_default
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


        module subprograms_default
            implicit none

            interface meshgrid
                module procedure meshgrid2D
                module procedure meshgrid3D
            end interface

            interface csvread
                module procedure csvread_real
                module procedure csvread_integer
            end interface

            interface csvwrite
                module procedure csvwrite1D_real
                module procedure csvwrite1D_integer
                module procedure csvwrite2D
            end interface

            interface define_BC
                module procedure define_vectorBC_type
                module procedure define_scalarBC_type
            end interface

            interface read_input_variable
                module procedure read_input_vector_variable
                module procedure read_input_scalar_variable
            end interface

            interface relative_residual
                module procedure relative_residual_vector
                module procedure relative_residual_scalar
            end interface

            interface integrate
                module procedure integrate1D
                module procedure integrate2D
            end interface


        contains

            subroutine meshgrid2D(X, Y, x_vec, y_vec)
                integer i,j
                real, intent(in) , dimension(:)   :: x_vec, y_vec
                real, intent(out), dimension(:,:) :: X, Y

                X = spread(x_vec, 1, size(y_vec) )
                Y = spread(y_vec, 2, size(x_vec) )

            end subroutine


            subroutine meshgrid3D(X, Y, Z, x_vec, y_vec, z_vec)
                integer i
                real, intent(in) , dimension(:)     :: x_vec, y_vec, z_vec
                real, intent(out), dimension(:,:,:) :: X, Y, Z
                  
                do i=1,size(z_vec)
                    X(:,:,i) = spread( x_vec, 2, size(y_vec) )
                    Y(:,:,i) = spread( y_vec, 1, size(x_vec) )
                end do
                do i=1,size(x_vec)
                    Z(i,:,:) = spread( z_vec, 1, size(y_vec) )
                end do
            end subroutine


            subroutine define_vectorBC_type(domain, Var, VarX_BC, VarY_BC, VarZ_BC)
                use Class_Geometry
                use Class_Vector_Variable

                type (duct) :: domain
                type (Vector_Variable) :: Var
                type (Vector_Variable_BC), dimension(domain%NumberBoundaries) :: VarX_BC, VarY_BC
                type (Vector_Variable_BC), optional, dimension(domain%NumberBoundaries) :: VarZ_BC

                do i = 1, domain%NumberBoundaries
                    if (Var%typeBC(i) == 'zero_gradient') then
                        VarX_BC(i)%boundary_condition => neumann_vect
                        VarY_BC(i)%boundary_condition => neumann_vect
                        if (present(VarZ_BC)) then
                            VarZ_BC(i)%boundary_condition => neumann_vect
                        end if
                    else if (Var%typeBC(i) == 'finite_cond_wall') then
                        VarX_BC(i)%boundary_condition => robin_vect
                        VarY_BC(i)%boundary_condition => robin_vect
                        if (present(VarZ_BC)) then
                            VarZ_BC(i)%boundary_condition => robin_vect
                        end if
                    else
                        VarX_BC(i)%boundary_condition => dirichlet_vect
                        VarY_BC(i)%boundary_condition => dirichlet_vect
                        if (present(VarZ_BC)) then
                            VarZ_BC(i)%boundary_condition => dirichlet_vect
                        end if
                    end if

                    if (Var%typeBC(i) == 'parabolic_inlet') then
                        VarX_BC(i)%value = (6.0*VarX_BC(i)%characteristic_valueBC/y%length**2)*(y%nodes(y%in_beg:y%in_end)-y%begin)*(y%length-(y%nodes(y%in_beg:y%in_end)-y%begin))
                        VarY_BC(i)%value = VarY_BC(i)%characteristic_valueBC
                    else
                        VarX_BC(i)%value = VarX_BC(i)%characteristic_valueBC
                        VarY_BC(i)%value = VarY_BC(i)%characteristic_valueBC
                    end if
                    if (present(VarZ_BC)) then
                        VarZ_BC(i)%value = VarZ_BC(i)%characteristic_valueBC
                    end if
                end do

            end subroutine


            subroutine define_scalarBC_type(domain, Var, Var_BC)
                use Class_Geometry
                use Class_Scalar_Variable

                type (duct) :: domain
                type (Scalar_Variable) :: Var
                type (Scalar_Variable_BC), dimension(domain%NumberBoundaries) :: Var_BC

                do i = 1, domain%NumberBoundaries
                    if (Var%typeBC(i) == 'zero_gradient') then
                        Var_BC(i)%boundary_condition => neumann_scalar
                    else if (Var%typeBC(i) == 'finite_cond_wall') then
                        Var_BC(i)%boundary_condition => robin_scalar
                    else
                        Var_BC(i)%boundary_condition => dirichlet_scalar
                    end if
                    Var_BC(i)%value = Var_BC(i)%characteristic_valueBC
                end do

            end subroutine


            subroutine read_input_vector_variable(domain, Var, Var_str, VarX_BC, VarY_BC, VarZ_BC)
                use global_variables
                use Class_Geometry
                use Class_Vector_Variable

                type (duct) :: domain
                type (Vector_Variable) :: Var
                type (Vector_Variable_BC), dimension(domain%NumberBoundaries) :: VarX_BC, VarY_BC
                type (Vector_Variable_BC), optional, dimension(domain%NumberBoundaries) :: VarZ_BC
                character(len=*) Var_str

                ! Check if input file exists
                inquire(file = 'input', exist = file_exists)

                ! Open input file
                if (file_exists .eqv. .TRUE.) then
                    open(unit = unit, file = 'input', status = 'old', action = 'read')
                else
                    stop 'Error: Input file does not exist!'
                end if

                ! Define varialbe typeBC matrix
                call Var%define_typeBC(domain)

                ! Read variables
                end_loop = .FALSE.
                do while (end_loop .eqv. .FALSE.)
                    read(unit,*) line
                    do i = 1, domain%NumberBoundaries
                        if (line==Var_str//'_'//boundary_names(i)) then
                            read(unit,*) Var%typeBC(i)
                            if (Var%typeBC(i)=='no_slip') then
                                VarX_BC(i)%characteristic_valueBC = 0.0
                                VarY_BC(i)%characteristic_valueBC = 0.0
                                if (flow_2D .eqv. .FALSE.)      VarZ_BC(i)%characteristic_valueBC = 0.0
                            else if (Var%typeBC(i)=='fixed_zero') then
                                VarX_BC(i)%characteristic_valueBC = 0.0
                                VarY_BC(i)%characteristic_valueBC = 0.0
                                if (flow_2D .eqv. .FALSE.)      VarZ_BC(i)%characteristic_valueBC = 0.0
                            else if (Var%typeBC(i)=='fixed') then
                                read(unit,*) VarX_BC(i)%characteristic_valueBC
                                read(unit,*) VarY_BC(i)%characteristic_valueBC
                                if (flow_2D .eqv. .FALSE.)  then
                                    read(unit,*) VarZ_BC(i)%characteristic_valueBC
                                    VarY_BC(i)%characteristic_valueBC = 0.0
                                end if
                            else if (Var%typeBC(i)=='parabolic_inlet')  then
                                read(unit,*) VarX_BC(i)%characteristic_valueBC
                            else if (Var%typeBC(i)=='custom_inlet') then
                                read(unit,*) VarZ_BC(i)%custom_inlet_file
                                read(unit,*) VarZ_BC(i)%custom_inlet_coeff
                                VarX_BC(i)%characteristic_valueBC = 0.0
                                VarY_BC(i)%characteristic_valueBC = 0.0
                            else if (Var%typeBC(i)=='finite_cond_wall') then
                                read(unit,*) VarX_BC(i)%coeff_robin
                                VarX_BC(i)%coeff_robin = 1.0/VarX_BC(i)%coeff_robin
                                VarY_BC(i)%coeff_robin = 1.0/VarX_BC(i)%coeff_robin
                                if (flow_2D .eqv. .FALSE.)  VarZ_BC(i)%coeff_robin = 1.0/VarX_BC(i)%coeff_robin
                                VarX_BC(i)%characteristic_valueBC = 0.0
                                VarY_BC(i)%characteristic_valueBC = 0.0
                                if (flow_2D .eqv. .FALSE.)  VarZ_BC(i)%characteristic_valueBC = 0.0
                            else if (Var%typeBC(i)=='zero_gradient') then
                                VarX_BC(i)%characteristic_valueBC = 0.0
                                VarY_BC(i)%characteristic_valueBC = 0.0
                                if (flow_2D .eqv. .FALSE.)      VarZ_BC(i)%characteristic_valueBC = 0.0
                            else
                                stop 'Error: Invalid Boundary Condition type!'
                            end if
                        end if
                    end do
                    if (line==Var_str//'_init') then
                        read(unit,*) Var%typeInit
                        if (Var%typeInit=='fixed_zero') then
                            Var%init(1) = 0.0
                            Var%init(2) = 0.0
                            if (flow_2D .eqv. .FALSE.)      Var%init(3) = 0.0
                        else if (Var%typeInit=='fixed') then
                            read(unit,*) Var%init(1)
                            read(unit,*) Var%init(2)
                            if (flow_2D .eqv. .FALSE.)      read(unit,*) Var%init(3)
                        else if (Var%typeInit=='same_as_inlet') then
                        else
                            stop 'Error: Invalid Initial Condition type!'
                        end if
                    end if
                    if (line=='end')                end_loop = .TRUE.
                end do

                ! Close file from memory
                close(unit)
            end subroutine


            subroutine read_input_scalar_variable(domain, Var, Var_str, Var_BC)
                use global_variables
                use Class_Geometry
                use Class_Scalar_Variable

                class (duct) :: domain
                type (Scalar_Variable) :: Var
                type (Scalar_Variable_BC), dimension(:) :: Var_BC
                character(len=*) Var_str

                ! Check if input file exists
                inquire(file = 'input', exist = file_exists)
                
                ! Open input file
                if (file_exists .eqv. .TRUE.) then
                    open(unit = unit, file = 'input', status = 'old', action = 'read')
                else
                    stop 'Error: Input file does not exist!'
                end if

                ! Define varialbe typeBC matrix
                call Var%define_typeBC(domain)

                ! Read variables
                end_loop = .FALSE.
                do while (end_loop .eqv. .FALSE.)
                    read(unit,*) line
                    do i = 1, domain%NumberBoundaries
                        if (line==Var_str//'_'//boundary_names(i)) then
                            read(unit,*) Var%typeBC(i)
                            if (Var%typeBC(i)=='no_slip') then
                                Var_BC(i)%characteristic_valueBC = 0.0
                            else if (Var%typeBC(i)=='fixed_zero') then
                                Var_BC(i)%characteristic_valueBC = 0.0
                            else if (Var%typeBC(i)=='fixed') then
                                read(unit,*) Var_BC(i)%characteristic_valueBC
                            else if (Var%typeBC(i)=='custom_inlet') then
                                read(unit,*) Var_BC(i)%custom_inlet_file
                                read(unit,*) Var_BC(i)%custom_inlet_coeff
                            else if (Var%typeBC(i)=='zero_gradient') then
                                Var_BC(i)%characteristic_valueBC = 0.0
                            else if (Var%typeBC(i)=='finite_cond_wall') then
                                read(unit,*) Var_BC(i)%coeff_robin
                                Var_BC(i)%coeff_robin = 1.0/Var_BC(i)%coeff_robin
                                Var_BC(i)%characteristic_valueBC = 0.0
                            else
                                stop 'Error: Invalid Boundary Condition type!'
                            end if
                        end if
                    end do
                    if (line==Var_str//'_init') then
                        read(unit,*) Var%typeInit
                        if (Var%typeInit=='fixed_zero') then
                            Var%init = 0.0
                        else if (Var%typeInit=='fixed') then
                            read(unit,*) Var%init
                        else
                            stop 'Error: Invalid Initial Condition type!'
                        end if
                    end if
                    if (line=='end')                end_loop = .TRUE.
                end do

                ! Close file from memory
                close(unit)

            end subroutine


            subroutine csvread_real(filename, A)
                integer i, n
                real, dimension(:)  :: A
                integer :: unit = 10
                character(len=*), intent(in) :: filename

                n = size(A)

                open(unit = unit, file = filename, status = 'old', action = 'read')
                do i = 1, n
                    read(unit,*) A(i)
                end do
                close(unit)
            end subroutine


            subroutine csvread_integer(filename, A)
                integer i, n
                integer, dimension(:)  :: A
                integer :: unit = 10
                character(len=*), intent(in) :: filename

                n = size(A)

                open(unit = unit, file = filename, status = 'old', action = 'read')
                do i = 1, n
                    read(unit,*) A(i)
                end do
                close(unit)
            end subroutine


            subroutine csvwrite1D_real(filename, A)
                integer i, n
                real, dimension(:)  :: A
                integer :: unit = 10
                character(len=*), intent(in) :: filename

                n = size(A)

                open(unit = unit, file = filename, status = 'replace', action = 'write')
                do i = 1, n
                    write(unit,'(F25.12)') A(i)
                end do
                close(unit)
            end subroutine


            subroutine csvwrite1D_integer(filename, A)
                integer i, n
                integer, dimension(:)  :: A
                integer :: unit = 10
                character(len=*), intent(in) :: filename

                n = size(A)

                open(unit = unit, file = filename, status = 'replace', action = 'write')
                do i = 1, n
                    write(unit,'(I15)') A(i)
                end do
                close(unit)
            end subroutine


           subroutine csvwrite2D(filename, A)
                integer col, row, i, j
                real, dimension(:,:)  :: A
                integer :: unit = 10
                character(len=*), intent(in) :: filename

                row = size(A,1)
                col = size(A,2)

                open(unit = unit, file = filename, status = 'replace', action = 'write')
                do j = 1, row
                    do i = 1, col
                        write(unit,'(F18.12)', advance='no') A(j,i)
                        if (i < col) write(unit,'(A)', advance='no') ','
                    end do
                    write (unit,*)
                end do
                close(unit)
            end subroutine


            function convert2Dto1D(Ar, I_max, J_max)    result(m)
                integer i_count, j_count, I_max, J_max
                real :: Ar(:,:), m(I_max*J_max)

                do j_count = 1, J_max
                    do i_count = 1, I_max
                        m(i_count+(j_count-1)*I_max) = Ar(i_count,j_count)
                    end do
                end do
            end function


            function convert1Dto2D(Ar, I_max, J_max)    result(m)
                integer i_count, j_count, I_max, J_max
                real :: Ar(:), m(I_max,J_max)

                do j_count = 1, J_max
                    do i_count = 1, I_max
                        m(i_count,j_count)= Ar(i_count+(j_count-1)*I_max)
                    end do
                end do
            end function


            function convert3Dto1D(Ar, I_max, J_max, K_max)    result(m)
                integer i_count, j_count, k_count, I_max, J_max, K_max
                real :: Ar(:,:,:), m(I_max*J_max*K_max)

                do k_count = 1, K_max
                    do j_count = 1, J_max
                        do i_count = 1, I_max
                            m(i_count+(j_count-1)*I_max+(k_count-1)*I_max*J_max) = Ar(i_count,j_count,k_count)
                        end do
                    end do
                end do
            end function


            function convert1Dto3D(Ar, I_max, J_max, K_max)    result(m)
                integer i_count, j_count, k_count, I_max, J_max, K_max
                real :: Ar(:), m(I_max,J_max,K_max)

                do k_count = 1, K_max
                    do j_count = 1, J_max
                        do i_count = 1, I_max
                            m(i_count,j_count,k_count)= Ar(i_count+(j_count-1)*I_max+(k_count-1)*I_max*J_max)
                        end do
                    end do
                end do
            end function


            subroutine enter_last_result_dir()
                use global_variables

                open(unit = unit, file = '.simulation_info', status = 'old', action = 'read')
                do
                    read(unit,*) line
                    if (line=='last_saved_iter')   then
                        read(unit,*) start_time
                        exit
                    end if
                end do
                close(unit)
                write(dir_name,'(I10)') start_time
                open(unit = unit, file = 'dirname.txt', status = 'replace', action = 'write')
                write(unit,'(A)', advance='no') '../results/'
                write(unit,'(A)', advance='no') adjustl(dir_name)
                close(unit)
                open(unit = unit, file = 'dirname.txt', status = 'old', action = 'read')
                read(unit,'(A)') dir_name
                close(unit)
                
                call system('rm dirname.txt')
                call chdir(dir_name)
            end subroutine


            function obtain_variable(filename, var_name)    result(var)     
                real var
                character(len=*) :: filename, var_name
                character(len=100) :: line
                logical end_loop
                integer :: unit = 10

                open(unit = unit, file = filename, status = 'old', action = 'read')
                end_loop = .FALSE.
                do while (end_loop .eqv. .FALSE.)
                    read(unit,*) line
                    if (line==var_name)   then
                        read(unit,*) var
                        exit
                    end if
                    if (line=='end')        end_loop = .TRUE.
                end do
                close(unit)
            end function


            function integrate1D(u, x)    result(R)
                real, dimension(:) :: u, x
                real R
                integer i, n

                n = size(x) - 1
                R = 0
                do i = 1, n
                    R = R + 0.5*(x(i+1)-x(i))*(u(i)+u(i+1))
                end do
            end function


            function integrate2D(u, x, y)    result(R)
                real, dimension(:,:) :: u
                real, dimension(:)   :: x, y
                real R
                integer i, j, i_end, j_end 
                real :: Ry(size(x)-1)

                i_end = size(x) - 1
                j_end = size(y) - 1
                R = 0
                Ry = 0
                do i = 1, i_end 
                    do j = 1, j_end 
                        ! Ry(i) = Ry(i) + 0.5*(y(j+1)-y(j))*(u(i,j)+u(i,j+1))
                        Ry(i) = Ry(i) + abs(y(j+1)-y(j))*(u(i,j))
                    end do
                end do
                do i = 1, i_end 
                    ! R = R + 0.5*(x(i+1)-x(i))*(Ry(i)+Ry(i+1))
                    R = R + abs(x(i+1)-x(i))*(Ry(i))
                end do
            end function


            function laplacian(var)       result(m)
                use global_variables
                use Class_Vector_Variable

                real, dimension(:) :: var
                real :: m(TotalCells)

                m = k_boole*aT_lapl*var(nodes_T)+aN_lapl*var(nodes_N)+aE_lapl*var(nodes_E)+aP_lapl*var(nodes_P)+aW_lapl*var(nodes_W)+aS_lapl*var(nodes_S)+k_boole*aB_lapl*var(nodes_B)
            end function


            function gradx(var)     result(m)
                use global_variables
                use Class_Scalar_Variable

                real, dimension(:) :: var
                real :: m(TotalCells)

                m = aE_gradx*var(nodes_E) + aP_gradx*var(nodes_P) + aW_gradx*var(nodes_W)
            end function


            function grady(var)     result(m)
                use global_variables
                use Class_Scalar_Variable

                real, dimension(:) :: var
                real :: m(TotalCells)

                m = aN_grady*var(nodes_N) + aP_grady*var(nodes_P) + aS_grady*var(nodes_S)
            end function


            function gradz(var)     result(m)
                use global_variables
                use Class_Scalar_Variable

                real, dimension(:) :: var
                real :: m(TotalCells)

                m = aT_gradz*var(nodes_T) + aP_gradz*var(nodes_P) + aB_gradz*var(nodes_B)
            end function


            function crossx(A_vect,B_vect)    result(m)
                use global_variables
                use Class_Vector_Variable

                type (Vector_variable) :: A_vect, B_vect
                real :: m(TotalCells)

                m = A_vect%y(Nodes_P)*B_vect%z(Nodes_P) - A_vect%z(Nodes_P)*B_vect%y(Nodes_P)
            end function


            function crossy(A_vect,B_vect)    result(m)
                use global_variables
                use Class_Vector_Variable

                type (Vector_variable) :: A_vect, B_vect
                real :: m(TotalCells)

                m = - A_vect%x(Nodes_P)*B_vect%z(Nodes_P) + A_vect%z(Nodes_P)*B_vect%x(Nodes_P)
            end function


            function crossz(A_vect,B_vect)    result(m)
                use global_variables
                use Class_Vector_Variable

                type (Vector_variable) :: A_vect, B_vect
                real :: m(TotalCells)

                m = A_vect%x(Nodes_P)*B_vect%y(Nodes_P) - A_vect%y(Nodes_P)*B_vect%x(Nodes_P)
            end function


            function curlx(A_vect)      result(m)
                use global_variables
                use Class_Vector_Variable

                type (Vector_variable) :: A_vect
                real :: m(TotalCells)

                m = grady(A_vect%z) - gradz(A_vect%y)
            end function


            function curly(A_vect)      result(m)
                use global_variables
                use Class_Vector_Variable

                type (Vector_variable) :: A_vect
                real :: m(TotalCells)

                m = gradz(A_vect%x) - gradx(A_vect%z)
            end function


            function curlz(A_vect)      result(m)
                use global_variables
                use Class_Vector_Variable

                type (Vector_variable) :: A_vect
                real :: m(TotalCells)

                m = gradx(A_vect%y) - grady(A_vect%x)
            end function


            function div(Var)       result(m)
                use global_variables
                use Class_Vector_Variable

                type (Vector_variable) :: Var
                real :: m(TotalCells)

                m = aE_divx*Var%x(nodes_E)+aP_divx*Var%x(nodes_P)+aW_divx*Var%x(nodes_W) + aN_divy*Var%y(nodes_N)+aP_divy*Var%y(nodes_P)+aS_divy*Var%y(nodes_S) + k_boole*(aT_divz*Var%z(nodes_T)+aP_divz*Var%z(nodes_P)+aB_divz*Var%z(nodes_B))
            end function


            function flux(Adv,var)      result(m)
                use global_variables
                use Class_Advection

                type (Advection) :: Adv
                real :: var(:) 
                real :: m(TotalCells)

                m = k_boole*Adv%aT*var(nodes_T) + Adv%aN*var(nodes_N) + Adv%aE*var(nodes_E) + Adv%aP*var(nodes_P) + Adv%aW*var(nodes_W) + Adv%aS*var(nodes_S) + k_boole*Adv%aB*var(nodes_B)

            end function
            

            subroutine relative_residual_vector(Var, Var0, small_num)
                use global_variables
                use Class_Vector_Variable
                use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

                type (Vector_Variable) :: Var, Var0
                real :: Error(TotalNodes)
                real res, small_num

                ! x-component error calculation
                Error = abs((Var%x-Var0%x))/abs(Var%x)
                where( Var%x<=small_num ) Error = 0.0
                where( ieee_is_finite(Error) ) 
                    Error = Error
                elsewhere
                    Error = 0.0
                end where
                Var%e(1) = sum(Error)/size(Error)

                ! y-component error calculation
                Error = abs((Var%y-Var0%y))/abs(Var%y)
                where( Var%y<=small_num ) Error = 0.0
                where( ieee_is_finite(Error) ) 
                    Error = Error
                elsewhere
                    Error = 0.0
                end where
                Var%e(2) = sum(Error)/size(Error)

                ! z-component error calculation
                if (flow_2D .eqv. .FALSE.)  then
                    Error = abs((Var%z-Var0%z))/abs(Var%z)
                    where( Var%z<=small_num ) Error = 0.0
                    where( ieee_is_finite(Error) ) 
                        Error = Error
                    elsewhere
                        Error = 0.0
                    end where
                    Var%e(3) = sum(Error)/size(Error)
                end if
            end subroutine


            subroutine relative_residual_scalar(Var, Var0, small_num)
                use global_variables
                use Class_Scalar_Variable
                use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

                type (Scalar_Variable) :: Var, Var0
                real :: Error(TotalNodes)
                real res, small_num

                ! Error calculation
                Error = abs((Var%field-Var0%field))/abs(Var%field)
                where( Var%field<=small_num ) Error = 0.0
                where( ieee_is_finite(Error) ) 
                    Error = Error
                elsewhere
                    Error = 0.0
                end where
                Var%e = sum(Error)/size(Error)

            end subroutine


            function mean_divergence_error(Var)     result(mean_val)
                use global_variables
                use Class_Vector_Variable

                type (Vector_Variable) :: Var
                real :: mean_val
                real :: B(TotalCells)

                ! B = aE_divx*Var%x(nodes_E)+aP_divx*Var%x(nodes_P)+aW_divx*Var%x(nodes_W) + aN_divy*Var%y(nodes_N)+aP_divy*Var%y(nodes_P)+aS_divy*Var%y(nodes_S) + k_boole*(aT_divz*Var%z(nodes_T)+aP_divz*Var%z(nodes_P)+aB_divz*Var%z(nodes_B))
                B = div(Var)
                mean_val = sum(B)/size(B)

            end function


            subroutine extrapolate(Var, domain)
                use global_variables
                use Class_Geometry

                type (duct) :: domain
                real, dimension(:) :: Var

                do i = 1, domain%NumberBoundaries
                    Var(domain%boundary(i)%nodes_bound) = Var(domain%boundary(i)%nodes_next2) + ((domain%boundary(i)%nodes_dist1+domain%boundary(i)%nodes_dist2)/(domain%boundary(i)%nodes_dist2))*(Var(domain%boundary(i)%nodes_next1)-Var(domain%boundary(i)%nodes_next2))
                end do

            end subroutine

    	end module
