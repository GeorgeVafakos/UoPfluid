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
        !   File name:          subprograms_print
        !   
        !   Type:               module
        ! 
        !   Description:        In this file are stored all the subroutines that are used to print to file, screen or to save data.
        ! 



        module subprograms_print
            implicit none

            interface print_screen_matrix
                module procedure print_screen_matrix1D_real
                module procedure print_screen_matrix1D_int
                module procedure print_screen_matrix2D_real
                module procedure print_screen_matrix2D_int
            end interface


        contains

            subroutine print_initial_conditions()
                use global_variables
                use subprograms_default
                use define_classes

                real*8 :: empty(1)=0.0

                write(dir_name,'(I10)') Iter_count
                call chdir('../results')
                call execute_command_line ('mkdir -p ../results' // dir_name )

                ! Print nodes
                call csvwrite('x_nodes.csv', x%nodes)
                call csvwrite('y_nodes.csv', y%nodes)
                call csvwrite('x_faces.csv', x%faces)
                call csvwrite('y_faces.csv', y%faces)
                
                ! Print variables
                call csvwrite('u.csv',V%x)
                call csvwrite('v.csv',V%y)
                call csvwrite('p.csv',p%field)

                ! Print residuals
                call csvwrite('Err_u.csv', empty)
                call csvwrite('Err_v.csv', empty)
                call csvwrite('Err_p.csv', empty)
                call csvwrite('Iterations.csv', empty)

                ! 3rd dimension
                if (flow_2D .eqv. .FALSE.)  then
                    call csvwrite('z_nodes.csv', z%nodes)
                    call csvwrite('w.csv',V%z)
                    call csvwrite('Err_w.csv', empty)

                end if
                
                call system('cp ../plots.py plots.py')
                call system('mv *.csv ' // dir_name)
                call system('mv *.py ' // dir_name)
                call system('cp .simulation_info ' // dir_name )
                
            end subroutine


            subroutine print_sim_results()
                use global_variables
                use subprograms_default
                use define_classes

                write(dir_name,'(I10)') Iter_count
                call chdir('../results')
                call execute_command_line ('mkdir -p ../results' // dir_name )

                ! Print nodes
                call csvwrite('x_nodes.csv', x%nodes)
                call csvwrite('y_nodes.csv', y%nodes)
                call csvwrite('x_faces.csv', x%faces)
                call csvwrite('y_faces.csv', y%faces)
                
                ! Print variables
                call csvwrite('u.csv',V%x)
                call csvwrite('v.csv',V%y)
                call csvwrite('p.csv',p%field)

                ! Print residuals
                call csvwrite('Err_u.csv', Err_Vx)
                call csvwrite('Err_v.csv', Err_Vy)
                call csvwrite('Err_p.csv', Err_p)
                call csvwrite('Iterations.csv', Iterations)

                ! 3rd dimension
                if (flow_2D .eqv. .FALSE.)  then
                    call csvwrite('z_nodes.csv', z%nodes)
                    call csvwrite('w.csv',V%z)
                    call csvwrite('Err_w.csv', Err_Vz)
                end if
                
                call system('cp ../plots.py plots.py')
                call system('mv *.csv ' // dir_name)
                call system('mv *.py ' // dir_name)
                call system('cp .simulation_info ' // dir_name )
            end subroutine


            subroutine print_sim_results_final()
                use global_variables
                use subprograms_default
                use define_classes

                call execute_command_line ('mkdir -p ../results/final')
                call chdir('../results/final')

                ! Print nodes
                call csvwrite('x_nodes.csv', x%nodes)
                call csvwrite('y_nodes.csv', y%nodes)
                call csvwrite('x_faces.csv', x%faces)
                call csvwrite('y_faces.csv', y%faces)
                
                ! Print variables
                call csvwrite('u.csv',V%x)
                call csvwrite('v.csv',V%y)
                call csvwrite('p.csv',p%field)

                ! Print residuals
                call csvwrite('Err_u.csv', Err_Vx)
                call csvwrite('Err_v.csv', Err_Vy)
                call csvwrite('Err_p.csv', Err_p)
                call csvwrite('Iterations.csv', Iterations)

                ! 3rd dimension
                if (flow_2D .eqv. .FALSE.)  then
                    call csvwrite('z_nodes.csv', z%nodes)
                    call csvwrite('w.csv',V%z)
                    call csvwrite('Err_w.csv', Err_Vz)

                end if

                call system('cp ../../plots.py plots.py')
            end subroutine


            subroutine print_screen_matrix1D_real(A)
                implicit none
                integer i
                real*8, dimension(:) :: A

                do i = 1, size(A)
                    write(*,*) A(i)
                end do
            end subroutine

            subroutine print_screen_matrix1D_int(A)
                implicit none
                integer i
                integer, dimension(:) :: A

                do i = 1, size(A)
                    write(*,*) A(i)
                end do
            end subroutine


            subroutine print_screen_matrix2D_real(A)
                implicit none
                integer i, j, row, col
                real(8), dimension(:,:) :: A

                row = size(A,1)
                col = size(A,2)

                do i = 1, col
                    do j = 1, row
                        write(*,'(F8.4)', advance="no") A(i,j)
                    end do
                    write(*,*)
                end do
            end subroutine


            subroutine print_screen_matrix2D_int(A)
                implicit none
                integer i, j, row, col
                integer, dimension(:,:) :: A

                row = size(A,1)
                col = size(A,2)

                do i = 1, col
                    do j = 1, row
                        write(*,'(I6)', advance="no") A(i,j)
                    end do
                    write(*,*)
                end do
            end subroutine


            subroutine print_UoPfluid_version_toscreen()
                write(*,*) ''
                write(*,*) '|----------------------------------------------------------------------------------------------------------------------------------------------'
                write(*,*) '|    ______    ______        _________                                                |                                                        '
                write(*,*) '|      ||        ||            ||    \\     ___                             ____      |   Version:  2.0                                       '
                write(*,*) '|      ||        ||            ||     ||   // \\ ==||                        ||       |                                                        '
                write(*,*) '|      ||        ||            ||     ||  ||  //   ||                        ||       |   Creator:  George Vafakos                            '
                write(*,*) '|      ||        ||    ____    ||____//   ||       ||  ____  ____   ο    ____||       |                                                        '
                write(*,*) '|      ||        ||   //  \\   ||       ==||==     ||   ||    ||  =||   //   ||       |   Site:     https://github.com/GeorgeVafakos/UoPfluid '
                write(*,*) '|      ||        ||  ||    ||  ||         ||       ||   ||    ||   ||  ||    ||       |                                                        '
                write(*,*) '|       \\______//    \\__// __||__     __||__   __||__  \\__//   _||_  \\___||_      |                                                        '
                write(*,*) '|                                                                                     |                                                        '
                write(*,*) '|----------------------------------------------------------------------------------------------------------------------------------------------'
                write(*,*) '|'
                write(*,*) '|'
            end subroutine


            subroutine print_UoPfluid_version_totextfile(unit)
                integer, intent(in) :: unit

                write(unit, '(A)') '#-----------------------------------------------------------------------------------------------------------------------------------------------------------'
                write(unit, '(A)') '#    ______    ______        _________                                                 !                                                                    '
                write(unit, '(A)') '#      ||        ||            ||    \\     ___                              ____      !   Version: 2.0                                                     '
                write(unit, '(A)') '#      ||        ||            ||     ||   // \\ ==||                         ||       !                                                                    '
                write(unit, '(A)') '#      ||        ||            ||     ||  ||  //   ||                         ||       !   Creator: George Vafakos                                          '
                write(unit, '(A)') '#      ||        ||    ____    ||____//   ||       ||  ____  ____    ο    ____||       !                                                                    '
                write(unit, '(A)') '#      ||        ||   //  \\   ||       ==||==     ||   ||    ||   =||   //   ||       !   Site:    https://github.com/GeorgeVafakos/UoPfluid               '
                write(unit, '(A)') '#      ||        ||  ||    ||  ||         ||       ||   ||    ||    ||  ||    ||       !                                                                    '
                write(unit, '(A)') '#       \\______//    \\__// __||__     __||__   __||__  \\___/|_  _||_  \\___||_      !                                                                    '
                write(unit, '(A)') '#                                                                                      !                                                                    '
                write(unit, '(A)') '#-----------------------------------------------------------------------------------------------------------------------------------------------------------'
                write(unit, '(A)') '# '
                write(unit, '(A)') '# '
            end subroutine


            subroutine print_UoPfluid_version_tosourcefile(unit)
                integer, intent(in) :: unit

                write(unit, '(A)') '        !-----------------------------------------------------------------------------------------------------------------------------------------------------------'
                write(unit, '(A)') '        !    ______    ______        _________                                                 !                                                                    '
                write(unit, '(A)') '        !      ||        ||            ||    \\     ___                              ____      !   Version: 2.0                                                     '
                write(unit, '(A)') '        !      ||        ||            ||     ||   // \\ ==||                         ||       !                                                                    '
                write(unit, '(A)') '        !      ||        ||            ||     ||  ||  //   ||                         ||       !   Creator: George Vafakos                                          '
                write(unit, '(A)') '        !      ||        ||    ____    ||____//   ||       ||  ____  ____    ο    ____||       !                                                                    '
                write(unit, '(A)') '        !      ||        ||   //  \\   ||       ==||==     ||   ||    ||   =||   //   ||       !   Site:    https://github.com/GeorgeVafakos/UoPfluid               '
                write(unit, '(A)') '        !      ||        ||  ||    ||  ||         ||       ||   ||    ||    ||  ||    ||       !                                                                    '
                write(unit, '(A)') '        !       \\______//    \\__// __||__     __||__   __||__  \\___/|_  _||_  \\___||_      !                                                                    '
                write(unit, '(A)') '        !                                                                                      !                                                                    '
                write(unit, '(A)') '        !-----------------------------------------------------------------------------------------------------------------------------------------------------------'
                write(unit, '(A)') '        ! '
                write(unit, '(A)') '        ! '
            end subroutine


            subroutine print_sim_info_screen()
                use global_variables

                call print_UoPfluid_version_toscreen()
                print '(A)',                 ' |    Solver:         mhd_pipe'
                print '(A)',                 ' | '
                print '(A)',                 ' | '
                print '(A)',                 ' |    Classes Created .................................... [OK]'
                print '(A)',                 ' | '
                print '(A)',                 ' |    Variables Loaded ................................... [OK]'
                print '(A)',                 ' | '
                print '(A)',                 ' |    Mesh Generated ..................................... [OK]'
                print '(A)',                 ' | '
                print '(A)',                 ' | '
                print '(A)',                 ' |    Mesh Info                                 Flow Parameters'
                print '(A)',                 ' |    =========                                 ==============='
                print '(A)',                 ' | '
                print '(A,I12,A22,F12.1)',    ' |    Total Cells:', TotalCells, 'Re =', Re
                print '(A,I12,A22,E12.1)',    ' |    Total Nodes:', TotalNodes, 'nu =', nu
                print '(A,I12,A22,E12.1)',    ' |    Total Faces:', TotalFaces, 'rho =', rho
                print '(A,I15,A22,E12.1)',    ' |    Cells x =', cells_x, 'Dt =', Dt
                print '(A,I15,A22,E12.1)',    ' |    Cells y =', cells_y
                print '(A,I15,A22,E12.1)',    ' |    Cells z =', k_boole*cells_z
                print '(A,F20.4,A22,E12.1)',   ' |    Px =', x%P
                print '(A,F20.4,A22,E12.1)',   ' |    Qx =', x%Q
                print '(A,F20.4,A22,E12.1)',   ' |    Py =', y%P
                print '(A,F20.4,A22,E12.1)',   ' |    Qy =', y%Q
                if (flow_2D .eqv. .FALSE.)  then
                print '(A,F20.4,A22,E12.1)',   ' |    Pz =', k_boole*z%P
                print '(A,F20.4,A22,F12.1)',   ' |    Qz =', k_boole*z%Q
                end if
                ! write(*,'(A,I10)')             ' |    Newton-Raphson Iter:', NewtRaph_Iter_y
                ! write(*,'(A,I9)')              ' |    Hartmann Layer Cells:', NumberCells_Layer_y
                ! write(*,'(A,I14)')             ' |    Number of nodes:', TotalNodes
                ! write(*,'(A,I14)')             ' |    Number of cells:', TotalCells
                ! write(*,'(A)')                 ' |'
                ! write(*,'(A)')                 ' |'
                write(*,*) 
                write(*,*) 
                write(*,*) 'Start'
                write(*,*) 
            end subroutine


            subroutine update_sim_info()
                use global_variables
                use Class_CFL

                call chdir('../results')
                open(unit = unit, file = '.simulation_info', status = 'replace', action = 'write')
                    call print_UoPfluid_version_totextfile(unit)
                    write(unit, '(A)') '#   File name:          simulation_info'
                    write(unit, '(A)') '# '
                    write(unit, '(A)') '#   Type:               text/hidden'
                    write(unit, '(A)') '# '
                    write(unit, '(A)') '#   Description:        In this file some crucial information about the current simulation are printed.'
                    write(unit, '(A)') '# '
                    write(unit, '(A)') '# '
                    write(unit, '(A)') 
                    write(unit, '(A)') '#  The latest saved time step (or final time step)'
                    write(unit, '(A)') 'last_saved_iter'
                    write(dir_name,'(I10)') Iter_count
                    write(unit, '(A)') trim(adjustl(dir_name))
                    write(unit, '(A)') 
                    write(unit, '(A)') 'last_saved_time'
                    write(dir_name,'(F18.4)') Time
                    write(unit, '(A)') trim(adjustl(dir_name))
                    write(unit, '(A)') 
                    write(unit, '(A)') '#  2D or 3D flow'
                    write(unit, '(A)') 'flow_2D'
                    if (flow_2D .eqv. .FALSE.)  then
                        write(unit,'(A)') 'false'
                    else
                        write(unit,'(A)') 'true'
                    end if
                    write(unit, '(A)')
                    write(unit, '(A)') '# Geometric info'
                    write(unit, '(A)') 'x_begin'
                    write(dir_name,'(F18.1)') x_begin
                    write(unit, '(A)') trim(adjustl(dir_name))
                    write(unit, '(A)') 'x_end'
                    write(dir_name,'(F18.1)') x_end
                    write(unit, '(A)') trim(adjustl(dir_name))
                    write(unit, '(A)') 'y_begin'
                    write(dir_name,'(F18.1)') y_begin
                    write(unit, '(A)') trim(adjustl(dir_name))
                    write(unit, '(A)') 'y_end'
                    write(dir_name,'(F18.1)') y_end
                    write(unit, '(A)') trim(adjustl(dir_name))
                    if (flow_2D .eqv. .FALSE.)  then
                        write(unit, '(A)') 'z_begin'
                        write(dir_name,'(F18.1)') z_begin
                        write(unit, '(A)') trim(adjustl(dir_name))
                        write(unit, '(A)') 'z_end'
                        write(dir_name,'(F18.1)') z_end
                        write(unit, '(A)') trim(adjustl(dir_name))
                    end if
                    write(unit, '(A)') 'TotalCells'
                    write(dir_name,'(I10)') TotalCells
                    write(unit, '(A)') trim(adjustl(dir_name))
                    write(unit, '(A)') 'TotalNodes'
                    write(dir_name,'(I10)') TotalNodes
                    write(unit, '(A)') trim(adjustl(dir_name))
                    write(unit, '(A)') 'TotalFaces'
                    write(dir_name,'(I10)') TotalFaces
                    write(unit, '(A)') trim(adjustl(dir_name))
                    write(unit, '(A)')
                    write(unit, '(A)') '# Flow parameters'
                    write(unit, '(A)') 'nu'
                    write(dir_name,'(F18.8)') nu
                    write(unit, '(A)') trim(adjustl(dir_name))
                    write(unit, '(A)') 'rho'
                    write(dir_name,'(F18.8)') rho
                    write(unit, '(A)') trim(adjustl(dir_name))
                    write(unit, '(A)') 'Dt'
                    write(dir_name,'(F18.8)') Dt
                    write(unit, '(A)') trim(adjustl(dir_name))
                    write(unit, '(A)') 'Co_max'
                    write(dir_name,'(F18.8)') Co%max
                    write(unit, '(A)') trim(adjustl(dir_name))
                    write(unit, '(A)')
                    write(unit, '(A)') '#------------------------------------------------------------------------------'
                    write(unit, '(A)') '# End of the file (the next line is necessary)'
                    write(unit, '(A)') 'end'
                    write(unit, '(A)') '#------------------------------------------------------------------------------'
                    write(unit, '(A)')
                close(unit)
            end subroutine

            end module
