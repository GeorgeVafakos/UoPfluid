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
        !   File name:          Class_Axis
        !
        !   Type:               module / class
        !
        !   Description:        In this file the axis class is created. 
        !
        !

        module Class_Axis
            implicit none

            type, public :: axis
                integer NumberNodes, NumberCells
                integer node_beg, node_end, cell_beg, cell_end
                procedure(one_way_stretch), pointer :: create_faces => null()
            contains
                procedure :: create_nodes => locate_cell_centers
                generic, public :: define => define_axis_index, define_global_index
                procedure, pass :: define_axis_index
                procedure, pass :: define_global_index
            end type

            type, extends (axis) :: cartesian_axis
                real*8 begin, end, length, stretching_point, perc_cells_Left
                integer NumberFaces, face_beg, face_end, in_end, in_beg
                real*8, allocatable, dimension(:)   :: faces, nodes
                real*8  P, Q, interior_point, mhd_Ha, mhd_LayerPoints
            end type cartesian_axis

            type, extends (axis) :: global_axis
            end type global_axis
            

        contains

            subroutine define_axis_index(ax, cells, ax_begin, ax_end)
                class (axis) :: ax
                integer ::  cells
                real*8  :: ax_begin, ax_end

                ax%NumberCells = cells

                select type (ax)
                    class is (cartesian_axis)

                        ax%begin = ax_begin
                        ax%end = ax_end

                        ax%NumberFaces = ax%NumberCells + 1                                         ! Number of faces in axis
                        ax%NumberNodes = ax%NumberCells + 2                                         ! Number of nodes in axis
                        ax%Length = ax%end - ax%begin                                               ! Length of axis
                        ax%face_beg = 1                                                             ! First face in axis
                        ax%face_end = ax%face_beg + (ax%NumberFaces-1)                              ! Last face in axis
                        ax%node_beg = 1                                                             ! First node in axis
                        ax%node_end = ax%node_beg + (ax%NumberNodes-1)                              ! Last node in axis
                        ax%in_beg = ax%node_beg + 1                                                 ! First internal node in axis
                        ax%in_end = ax%node_end - 1                                                 ! Last internal node in axis
                        ax%cell_beg = 1                                                             ! First cell in axis
                        ax%cell_end = ax%cell_beg + (ax%NumberCells-1)                              ! Last cell in axis

                    class default
                         stop 'define (axis index): unexpected type for axis object!'
                end select
            end subroutine


            subroutine define_global_index(global, ax1, ax2, ax3)
                class (axis) :: global
                class (axis) :: ax1, ax2
                class (axis), optional :: ax3

                select type (ax1)
                    class is (cartesian_axis)
                        select type (ax2)
                            class is (cartesian_axis)
                                global%NumberNodes = ax1%NumberNodes*ax2%NumberNodes
                                global%NumberCells = ax1%NumberCells*ax2%NumberCells
                                if (present(ax3)) then
                                    select type (ax3)
                                        class is (cartesian_axis)
                                            global%NumberNodes = ax1%NumberNodes*ax2%NumberNodes*ax3%NumberNodes
                                            global%NumberCells = ax1%NumberCells*ax2%NumberCells*ax3%NumberCells
                                        class default
                                            stop 'define_global_index: unexpected axis type'
                                    end select
                                else
                                    global%NumberNodes = ax1%NumberNodes*ax2%NumberNodes
                                    global%NumberCells = ax1%NumberCells*ax2%NumberCells
                                end if
                            class default
                                stop 'define_global_index: unexpected axis type'
                        end select

                    class default
                        stop 'define_global_index: unexpected axis type'
                end select

                global%node_beg = 1                                                                 ! First global node
                global%node_end = global%node_beg + (global%NumberNodes-1)                          ! Last global node
                global%cell_beg = 1                                                                 ! First global cell
                global%cell_end = global%cell_beg + (global%NumberCells-1)                          ! Last global cell

            end subroutine


            subroutine one_way_stretch(ax)
                class (axis) :: ax
                integer i
                real*8  n1, n2
                real*8, allocatable, dimension(:)  :: s, n

                select type (ax)
                    class is (cartesian_axis)

                        allocate(s(ax%NumberFaces), n(ax%NumberFaces))
                        n1 = 1.0
                        n2 = ax%NumberFaces
                        n = ((/(i, i=1,ax%NumberFaces)/) - n1)/(ax%NumberFaces - n1)
                        s = ax%P*n +(1.0-ax%P)*(1.0 - tanh(ax%Q*(1.0 - n))/tanh(ax%Q))
                        ax%faces = ax%begin + s*(ax%end - ax%begin)

                    class default
                        stop 'stretchFunct: unexpected axis type'
                end select
            end subroutine


            subroutine two_way_stretch(ax)
                class (axis) :: ax
                integer i
                real*8  n1, n2
                real*8, allocatable, dimension(:)  :: s, n

                select type (ax)
                    class is (cartesian_axis)

                        allocate(s(ax%NumberFaces), n(ax%NumberFaces))
                        n1 = 1.0
                        n2 = ax%NumberFaces
                        n = ((/(i, i=1,ax%NumberFaces)/) - n1)/(ax%NumberFaces - n1)
                        s = ax%P*n +(1.0-ax%P)*0.5*(1.0 + tanh(ax%Q*(n - 0.5))/tanh(0.5*ax%Q))
                        ax%faces = ax%begin + s*(ax%end - ax%begin)
                    class default
                        stop 'stretchFunct: unexpected axis type'
                end select
            end subroutine


            subroutine interior_point_stretch(ax)
                class (axis) :: ax
                integer i
                real*8  n1, b
                real*8, allocatable, dimension(:)  :: s, n

                select type (ax)
                    class is (cartesian_axis)

                        n1 = 1.0
                        n = ((/(i, i=1,ax%NumberFaces)/) - n1)/(ax%NumberFaces - n1)
                        b = (1/(2*ax%P))*log((1+(exp(ax%P)-1.0)*((ax%interior_point-ax%begin)/ax%length))/(1+(exp(-ax%P)-1.0)*((ax%interior_point-ax%begin)/ax%length)))
                        s = ((ax%interior_point-ax%begin)/ax%length)*(1+sinh(ax%P*(n-b))/sinh(ax%P*b))
                        ax%faces = ax%begin + s*(ax%end - ax%begin)

                    class default
                        stop 'stretchFunct: unexpected axis type'
                end select
            end subroutine


            subroutine locate_cell_centers(ax)
                class (axis) :: ax

                select type(ax)
                    class is (cartesian_axis)
                        allocate(ax%nodes(ax%NumberNodes))

                        ax%nodes(ax%node_beg) = ax%faces(ax%face_beg)
                        ax%nodes(ax%in_beg:ax%in_end) = (ax%faces(ax%face_beg:ax%face_end-1) + ax%faces(ax%face_beg+1:ax%face_end))/2.0
                        ax%nodes(ax%node_end) = ax%faces(ax%face_end)

                    class default
                        stop 'locate_midpoints: unexpected axis type'
                end select
            end subroutine

        end module
