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
        !   File name:          read_input
        !   
        !   Type:               source
        ! 
        !   Description:        This file reads all the input data and parameters from the 'input' text file.
        ! 
        ! 
        ! Check if the input file exists
        inquire(file = 'input', exist = file_exists)

        ! Open input file
        if (file_exists .eqv. .TRUE.) then
            open(unit = unit, file = 'input', status = 'old', action = 'read')
        else
            stop 'Error: Input file does not exist!'
        end if

        ! Read variables
        end_loop = .FALSE.
        do while (end_loop .eqv. .FALSE.)
            read(unit,*) line

            if (line=='flow_2D')                    read(unit,*) flow_2D
            if (line=='cells_x')                    read(unit,*) cells_x
            if (line=='cells_y')                    read(unit,*) cells_y
            if (line=='cells_z')                    read(unit,*) cells_z
            if (line=='x_begin')                    read(unit,*) x_begin
            if (line=='x_end')                      read(unit,*) x_end
            if (line=='y_begin')                    read(unit,*) y_begin
            if (line=='y_end')                      read(unit,*) y_end
            if (line=='z_begin')                    read(unit,*) z_begin
            if (line=='z_end')                      read(unit,*) z_end
            if (line=='nu')                         read(unit,*) nu
            if (line=='rho')                        read(unit,*) rho

            if (line=='X_stretch')  then
                read(unit,*) X_stretch
                if (X_stretch=='one_way_stretch')   then
                    read(unit,*) x%P
                    read(unit,*) x%Q
                    x%create_faces => one_way_stretch
                else if (X_stretch=='two_way_stretch')   then
                    read(unit,*) x%P
                    read(unit,*) x%Q
                    x%create_faces => two_way_stretch
                else if (X_stretch=='interior_point_stretch')   then
                    read(unit,*) x%P
                    read(unit,*) x%interior_point
                    x%create_faces => interior_point_stretch
                else
                    stop 'Error: Availiable X stretching functions: one_way_stretch, two_way_stretch or interior_point_stretch'
                end if
            end if

            if (line=='Y_stretch')  then
                read(unit,*) Y_stretch
                if (Y_stretch=='one_way_stretch')   then
                    read(unit,*) y%P
                    read(unit,*) y%Q
                    y%create_faces => one_way_stretch
                else if (Y_stretch=='two_way_stretch')   then
                    read(unit,*) y%P
                    read(unit,*) y%Q
                    y%create_faces => two_way_stretch
                else if (Y_stretch=='interior_point_stretch')   then
                    read(unit,*) y%P
                    read(unit,*) y%interior_point
                    y%create_faces => interior_point_stretch
                else
                    stop 'Error: Availiable Y stretching functions: one_way_stretch, two_way_stretch or interior_point_stretch'
                end if
            end if

            if (line=='Z_stretch')  then
                read(unit,*) Z_stretch
                if (Z_stretch=='one_way_stretch')   then
                    read(unit,*) z%P
                    read(unit,*) z%Q
                    z%create_faces => one_way_stretch
                else if (Z_stretch=='two_way_stretch')   then
                    read(unit,*) z%P
                    read(unit,*) z%Q
                    z%create_faces => two_way_stretch
                else if (Z_stretch=='interior_point_stretch')   then
                    read(unit,*) z%P
                    read(unit,*) z%interior_point
                    z%create_faces => interior_point_stretch
                else
                    stop 'Error: Availiable Z stretching functions: one_way_stretch, two_way_stretch or interior_point_stretch'
                end if
            end if

            if (line=='Dt')     then
                read(unit,*) Dt
                adjustable_time_step = .FALSE.
            end if

            if (line=='Co')     then
                read(unit,*) Co%fixed_value
                Dt = 0.001
                adjustable_time_step = .TRUE.
            end if

            if (line=='NS_scheme')  then
                read(unit,*) NS_scheme
                if (NS_scheme=='SOR')             read(unit,*) r
                if (NS_scheme/='Jacobi' .AND. NS_scheme/='Gauss_Seidel' .AND. NS_scheme/='SOR') stop 'Error: Availiable numerical schemes: Jacobi, Gauss_Seidel or SOR'
            end if

            if (line=='Pres_scheme')  then
                read(unit,*) Pres_scheme
                if (Pres_scheme=='SOR')             read(unit,*) r
                if (Pres_scheme/='Jacobi' .AND. Pres_scheme/='Gauss_Seidel' .AND. Pres_scheme/='SOR') stop 'Error: Availiable numerical schemes: Jacobi, Gauss_Seidel or SOR'
            end if

            if (line=='advection_scheme')  then
                read(unit,*) advection_scheme
                if (advection_scheme/='upwind' .AND. advection_scheme/='CD' .AND. advection_scheme/='none') stop 'Error: Availiable advection schemes: upwind, CD or none'
            end if

            if (line=='print_results')              read(unit,*) print_results
            if (line=='print_to_screen')            read(unit,*) print_to_screen
            if (line=='start')                      read(unit,*) start
            if (line=='end')                        end_loop = .TRUE.
        end do
        close(unit)


        !------------------------------------------------------------------------------------------
        ! Read Boundary Names
        !------------------------------------------------------------------------------------------
        call domain%allocate_bound_number
        allocate(boundary_names(domain%NumberBoundaries))

        ! Open input file
        if (file_exists .eqv. .TRUE.) then
            open(unit = unit, file = 'input', status = 'old', action = 'read')
        else
            stop 'Error: Input file does not exist!'
        end if

        ! Read variables
        end_loop = .FALSE.
        do while (end_loop .eqv. .FALSE.)
            read(unit,*) line

            if (line=='boundary_names') then
                do i = 1, domain%NumberBoundaries
                    read(unit,*) boundary_names(i)
                end do
            end if
            if (line=='end')                        end_loop = .TRUE.
        end do
        close(unit)

        
        !------------------------------------------------------------------------------------------
        ! Read Other Parameters
        !------------------------------------------------------------------------------------------
        ! Define a a boolean parameter for the solution schemes (2D=0 and 3D=1)
        if (flow_2D .eqv. .FALSE.)  then
            k_boole = 1
        else
            k_boole = 0
        end if        

        ! Read velocity boundary conditions
        allocate(Vx_BC(domain%NumberBoundaries))
        allocate(Vy_BC(domain%NumberBoundaries))
        allocate(Vz_BC(domain%NumberBoundaries))
        allocate(p_BC(domain%NumberBoundaries))
        call read_input_variable(domain, V, 'V', Vx_BC, Vy_BC, Vz_BC)
        call read_input_variable(domain, p, 'p', p_BC)
