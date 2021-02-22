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
        !   File name:          global_variables
        !   
        !   Type:               module
        ! 
        !   Description:        This is the file where all global variables are declared and used throughout the code.
        ! 
        ! 
        ! 


        module global_variables

            ! Geometric variables
            integer cells_x, NumberFaces_x, NumberNodes_x, cells_y, NumberFaces_y, NumberNodes_y, TotalNodes
            integer i, i_face_end, i_node_end, i_cell_end, i_in_end, i_in_beg, j, j_face_end, j_node_end, j_cell_end, j_in_end, j_in_beg
            integer :: i_face_beg=1, i_node_beg=1, i_cell_beg=1, j_face_beg=1, j_node_beg=1, j_cell_beg=1
            real, allocatable, dimension(:)   :: x_nodes,x_faces,Dx_vec,dxe_vec, dxw_vec, y_nodes,y_faces,Dy_vec,dyn_vec, dys_vec
            real, allocatable, dimension(:,:) :: Dx, Dy, dxe, dxw, dyn, dys
            real, allocatable, dimension(:,:) :: aN_lapl, aE_lapl, aP_lapl, aW_lapl, aS_lapl

            ! State variables
            real, allocatable, dimension(:,:) :: psi, psi0, u, v, omega
            real, allocatable, dimension(:)   :: u_vec, v_vec

            ! Other parameters
            real  cputime_start, cputime_finish, time_cpu
            integer  :: unit = 10, Iter_count=0, print_screen
            real r, e
            
        end module
