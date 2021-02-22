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
        !   File name:          subprograms
        !   
        !   Type:               module
        ! 
        !   Description:        In this file are stored all functions and subroutines that are used in the code.
        ! 
        !                       
        ! 

        module subprograms
            implicit none

            interface csvwrite
                module procedure csvwrite1D
                module procedure csvwrite2D
            end interface

        contains
        

            subroutine meshgrid(X, Y, x_vec, y_vec)
                integer i,j
                real, intent(in) , dimension(:)   :: x_vec, y_vec
                real, intent(out), dimension(:,:) :: X, Y

                X = spread(x_vec, 2, size(y_vec) )
                Y = spread(y_vec, 1, size(x_vec) )

            end subroutine


            function csvread(filename, n) result(A)
                integer i, n
                real, dimension(n)  :: A
                integer :: unit = 10
                character(len=*), intent(in) :: filename

                open(unit = unit, file = filename, status = 'old', action = 'read')
                do i = 1, n
                    read(unit,*) A(i)
                end do
                close(unit)
            end function


            function count_size(filename1)   result(n1)
                character(len=*), intent(in) :: filename1
                integer :: unit1 = 10, n1, io
                real m1
                n1=0
                open(unit = unit1, file = filename1, status = 'old', action = 'read')
                do 
                    n1 = n1+1
                    read(unit1,*, iostat=io) m1
                    if (io > 0) then
                        print *, "Somwthing went wrong."
                        exit
                    else if (io<0) then
                        exit
                    end if
                end do
                n1 = n1-1
                close(unit1)
            end function


            subroutine csvwrite1D(filename, A)
                integer i, n
                real, dimension(:)  :: A
                character(len=*), intent(in) :: filename

                n = size(A)

                open(unit = 9, file = filename, status = 'replace', action = 'write')
                do i = 1, n
                    write(9,"(F25.12)") A(i)
                end do
                close(9)
            end subroutine


           subroutine csvwrite2D(filename, A)
                    integer col, row, i, j
                    real, dimension(:,:)  :: A
                    character(len=*), intent(in) :: filename

                    row = size(A,1)
                    col = size(A,2)

                    open(unit = 9, file = filename, status = 'replace', action = 'write')
                    do j = 1, row
                        do i = 1, col
                            write(9,"(F18.12)", advance="no") A(j,i)
                            if (i < col) write(9,"(A)", advance="no") ','
                         end do
                        write (9,*)
                    end do
                    close(9)
            end subroutine


            function convert1Dto2D(Ar, I_max, J_max)    result(m)
                integer i_count, j_count, I_max, J_max
                real :: Ar(:), m(I_max,J_max)

                do j_count = 1, J_max
                    do i_count = 1, I_max
                        m(i_count,j_count)= Ar(i_count+(j_count-1)*I_max)
                    end do
                end do
            end function

            function convert2Dto1D(Ar, I_max, J_max)    result(m)
                integer i_count, j_count, I_max, J_max
                real :: Ar(:,:), m(I_max*J_max)

                do j_count = 1, J_max
                    do i_count = 1, I_max
                        m(i_count+(j_count-1)*I_max) = Ar(i_count,j_count)
                    end do
                end do
            end function


            function laplacian(Ar)       result(B)
                use global_variables

                real, dimension(:,:) :: Ar
                real :: B(cells_y,cells_x)

                B = Dx/dyn*Ar(j_in_beg+1:j_in_end+1,i_in_beg:i_in_end)+Dy/dxe*Ar(j_in_beg:j_in_end,i_in_beg+1:i_in_end+1)-(Dx/dyn + Dy/dxe + Dy/dxw + Dx/dys)*Ar(j_in_beg:j_in_end,i_in_beg:i_in_end)+Dy/dxw*Ar(j_in_beg:j_in_end,i_in_beg-1:i_in_end-1)+Dx/dys*Ar(j_in_beg-1:j_in_end-1,i_in_beg:i_in_end)
            end function


            function gradx(Ar)       result(B)
                use global_variables

                real, dimension(:,:) :: Ar
                real :: B(cells_x,cells_y)

                B = 0.5*Dx*Dy/dxe*Ar(i_in_beg+1:i_in_end+1,j_in_beg:j_in_end)+0.5*Dx*Dy*(1.0/dxw - 1.0/dxe)*Ar(i_in_beg:i_in_end,j_in_beg:j_in_end)-0.5*Dx*Dy/dxw*Ar(i_in_beg-1:i_in_end-1,j_in_beg:j_in_end)

            end function


            function grady(Ar)       result(B)
                use global_variables

                real, dimension(:,:) :: Ar
                real :: B(cells_x,cells_y)

                B = 0.5*Dx*Dy/dyn*Ar(i_in_beg:i_in_end,j_in_beg+1:j_in_end+1)+0.5*Dx*Dy*(1.0/dys - 1.0/dyn)*Ar(i_in_beg:i_in_end,j_in_beg:j_in_end)-0.5*Dx*Dy/dys*Ar(i_in_beg:i_in_end,j_in_beg-1:j_in_end-1)

            end function



            function relative_residual(Ar, Ar0, small_num)      result(res)
                use global_variables
                use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

                real, dimension(:,:) :: Ar ,Ar0
                real :: Error(NumberNodes_y, NumberNodes_x)
                real res, small_num

                Error = abs((Ar-Ar0))/abs(Ar)
                where( Ar<=small_num ) Error = 0.0
                where( ieee_is_finite(Error) ) 
                    Error = Error
                elsewhere
                    Error = 0.0
                end where
                res = sum(Error)/size(Error)

            end function


            subroutine print_sim_info_screen()
                print '(A)', ''
                print '(A)', '|--------------------------------------------------------------------------------------------------------------------------------------------'
                print '(A)', '|    ______    ______        _________                                               | '
                print '(A)', '|      ||        ||            ||    \\     ___                             ____     |   Version:  2.0  '
                print '(A)', '|      ||        ||            ||     ||   // \\ ==||                        ||      | '
                print '(A)', '|      ||        ||            ||     ||  ||  //   ||                        ||      |   Creator:  George Vafakos '
                print '(A)', '|      ||        ||    ____    ||____//   ||       ||  ____  ____   ο    ____||      | '
                print '(A)', '|      ||        ||   //  \\   ||       ==||==     ||   ||    ||  =||   //   ||      |   Site:     https://github.com/GeorgeVafakos/UoPfluid '
                print '(A)', '|      ||        ||  ||    ||  ||         ||       ||   ||    ||   ||  ||    ||      | '
                print '(A)', '|       \\______//    \\__// __||__     __||__   __||__  \\__//   _||_  \\___||_     | '
                print '(A)', '|                                                                                    | '
                print '(A)', '|--------------------------------------------------------------------------------------------------------------------------------------------'
                print '(A)', '|'
                print '(A)', '|'
                print '(A)', '|    Tool Name:      Streamlines'
                print '(A)', '| '
                print '(A)', '|    Classes Created .................................... [OK]'
                print '(A)', '| '
                print '(A)', '|    Variables Loaded ................................... [OK]'
                print '(A)', '| '
                print '(A)', '|    Mesh Generated ..................................... [OK]'
                print '(A)', '| '
                print '(A)', '| '
                print '(A)', '|    Start '
                print '(A)', '| '
                print '(A)', '| '
                print '(A)', '| '
            end subroutine




        end module