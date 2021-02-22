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
        !   Tools:            Streamline
        ! 
        !   Discription:      This is a software subroutine-like tool that can calculate the streamlines using a given velocity field.
        !                     No input file is required. Only the u, v velocities (1D global index) and the x,y nodes/faces in csv.
        ! 
        !   Algorithm:        SOR 
        ! 
        ! 
        ! 


        program main

        use global_variables
        use subprograms

        implicit none

        !------------------------------------------------------------------------------------------
        ! Mesh Generation
        !------------------------------------------------------------------------------------------
        ! x-axis
        NumberNodes_x = count_size('x_nodes.csv')
        NumberFaces_x = count_size('x_faces.csv')
        x_nodes = csvread('x_nodes.csv',NumberNodes_x)
        x_faces = csvread('x_faces.csv',NumberFaces_x)
        cells_x = NumberFaces_x - 1 
        i_face_end = i_face_beg + (NumberFaces_x-1) 
        i_node_end = i_node_beg + (NumberNodes_x-1) 
        i_in_beg = i_node_beg + 1 
        i_in_end = i_node_end - 1 
        i_cell_end = cells_x 
        Dx_vec = x_faces(i_face_beg+1:i_face_end) - x_faces(i_face_beg:i_face_end-1)
        dxe_vec = x_nodes(i_in_beg+1:i_in_end+1) - x_nodes(i_in_beg:i_in_end)
        dxw_vec = x_nodes(i_in_beg:i_in_end) - x_nodes(i_in_beg-1:i_in_end-1)

        ! y-axis
        NumberNodes_y = count_size('y_nodes.csv')
        NumberFaces_y = count_size('y_faces.csv')
        y_nodes = csvread('y_nodes.csv',NumberNodes_y)
        y_faces = csvread('y_faces.csv',NumberFaces_y)
        cells_y = NumberFaces_y - 1
        j_face_end = j_face_beg + (NumberFaces_y-1) 
        j_node_end = j_node_beg + (NumberNodes_y-1) 
        j_in_beg = j_node_beg + 1 
        j_in_end = j_node_end - 1 
        j_cell_end = cells_y 
        Dy_vec = y_faces(j_face_beg+1:j_face_end) - y_faces(j_face_beg:j_face_end-1)
        dyn_vec = y_nodes(j_in_beg+1:j_in_end+1) - y_nodes(j_in_beg:j_in_end)
        dys_vec = y_nodes(j_in_beg:j_in_end) - y_nodes(j_in_beg-1:j_in_end-1)

        ! Create grid
        TotalNodes = NumberNodes_x*NumberNodes_y
        allocate(Dx(cells_x,cells_y), Dy(cells_x,cells_y), dxe(cells_x,cells_y), dxw(cells_x,cells_y), dyn(cells_x,cells_y), dys(cells_x,cells_y))
        call meshgrid(Dx, Dy, Dx_vec, Dy_vec)
        call meshgrid(dxe, dyn, dxe_vec, dyn_vec)
        call meshgrid(dxw, dys, dxw_vec, dys_vec)

        !------------------------------------------------------------------------------------------
        ! Initial Conditions
        !------------------------------------------------------------------------------------------
        ! Print simulation info in screen
        call print_sim_info_screen()

        ! Load initial conditions
        allocate(psi(NumberNodes_x,NumberNodes_y), psi0(NumberNodes_x,NumberNodes_y), omega(NumberNodes_x,NumberNodes_y))
        psi0(i_node_beg:i_node_end,j_node_beg:j_node_end) = 0
        psi = psi0

        ! Load velocity field
        u_vec = csvread('u.csv',TotalNodes)
        v_vec = csvread('v.csv',TotalNodes)
        u = convert1Dto2D(u_vec, NumberNodes_x, NumberNodes_y)
        v = convert1Dto2D(v_vec, NumberNodes_x, NumberNodes_y)

        ! Calculate vorticity
        omega(i_in_beg:i_in_end,j_in_beg:j_in_end) = gradx(v) - grady(u)

        ! Calculate Laplace equation terms
        aN_lapl = Dx/dyn
        aE_lapl = Dy/dxe
        aP_lapl =-(Dy/dxe + Dy/dxw + Dx/dyn + Dx/dys)
        aW_lapl = Dy/dxw
        aS_lapl = Dx/dys


        ! Choose relacation for SOR algorithm
        r = 1.0

        !------------------------------------------------------------------------------------------
        ! Solve streamline equations
        !------------------------------------------------------------------------------------------
        do
            ! Time step starts
            Iter_count = Iter_count + 1

            call cpu_time(cputime_start)
            ! Solve Laplace equation with the SOR method
            do j = j_in_beg, j_in_end
                do i = i_in_beg, i_in_end
                   psi(i,j) = (1.0-r)*psi0(i,j) + r*(-omega(i,j) -aN_lapl(i-1,j-1)*psi0(i,j+1)-aE_lapl(i-1,j-1)*psi0(i+1,j)-aW_lapl(i-1,j-1)*psi(i-1,j)-aS_lapl(i-1,j-1)*psi(i,j-1))/aP_lapl(i-1,j-1)
               end do
            end do

            ! Left boundary conditions
            ! psi(i_node_beg, j_in_beg:j_in_end) = (dxw(i_cell_beg,:)*dxw(i_cell_beg+1,:)*(dxw(i_cell_beg,:)+dxw(i_cell_beg+1,:))*u(i_node_beg,j_in_beg:j_in_end) + ((dxw(i_cell_beg,:)+dxw(i_cell_beg+1,:))**2)*psi(i_node_beg+1,j_in_beg:j_in_end) - (dxw(i_cell_beg,:)**2)*psi(i_node_beg+2,j_in_beg:j_in_end))/(dxw(i_cell_beg+1,:)**2+2*dxw(i_cell_beg,:)*dxw(i_cell_beg+1,:))
            ! psi(i_node_beg, j_in_beg:j_in_end) = (((dxw(i_cell_beg,:)+dxw(i_cell_beg+1,:))**2)*psi(i_node_beg+1,j_in_beg:j_in_end) - (dxw(i_cell_beg,:)**2)*psi(i_node_beg+2,j_in_beg:j_in_end))/(dxw(i_cell_beg+1,:)**2+2*dxw(i_cell_beg,:)*dxw(i_cell_beg+1,:))
            psi(i_node_beg, j_in_beg:j_in_end) = 0

            ! Right boundary conditions
            ! psi(i_node_end, j_in_beg:j_in_end) = (dxe(i_cell_end,:)*dxe(i_cell_end-1,:)*(dxe(i_cell_end,:)+dxe(i_cell_end-1,:))*u(i_node_end,j_in_beg:j_in_end) + ((dxe(i_cell_end,:)+dxe(i_cell_end-1,:))**2)*psi(i_node_end-1,j_in_beg:j_in_end) - (dxe(i_cell_end,:)**2)*psi(i_node_end-2,j_in_beg:j_in_end))/(dxe(i_cell_end-1,:)**2+2*dxe(i_cell_end,:)*dxe(i_cell_end-1,:))
            ! psi(i_node_end, j_in_beg:j_in_end) = (((dxe(i_cell_end,:)+dxe(i_cell_end-1,:))**2)*psi(i_node_end-1,j_in_beg:j_in_end) - (dxe(i_cell_end,:)**2)*psi(i_node_end-2,j_in_beg:j_in_end))/(dxe(i_cell_end-1,:)**2+2*dxe(i_cell_end,:)*dxe(i_cell_end-1,:))
            psi(i_node_end, j_in_beg:j_in_end) = 0

            ! Bottom boundary conditions
            ! psi(i_in_beg:i_in_end, j_node_beg) = (dys(:,j_cell_beg)*dys(:,j_cell_beg+1)*(dys(:,j_cell_beg)+dys(:,j_cell_beg+1))*v(i_in_beg:i_in_end,j_node_beg) + ((dys(:,j_cell_beg)+dys(:,j_cell_beg+1))**2)*psi(i_in_beg:i_in_end,j_node_beg+1) - (dys(:,j_cell_beg)**2)*psi(i_in_beg:i_in_end,j_node_beg+2))/(dys(:,j_cell_beg+1)**2+2*dys(:,j_cell_beg)*dys(:,j_cell_beg+1))
            ! psi(i_in_beg:i_in_end, j_node_beg) = (((dys(:,j_cell_beg)+dys(:,j_cell_beg+1))**2)*psi(i_in_beg:i_in_end,j_node_beg+1) - (dys(:,j_cell_beg)**2)*psi(i_in_beg:i_in_end,j_node_beg+2))/(dys(:,j_cell_beg+1)**2+2*dys(:,j_cell_beg)*dys(:,j_cell_beg+1))
            psi(i_in_beg:i_in_end, j_node_beg) = 0

            ! Upper boundary conditions
            ! psi(i_in_beg:i_in_end, j_node_end) = (dyn(:,j_cell_end)*dyn(:,j_cell_end-1)*(dyn(:,j_cell_end)+dyn(:,j_cell_end-1))*v(i_in_beg:i_in_end,j_node_end) + ((dyn(:,j_cell_end)+dyn(:,j_cell_end-1))**2)*psi(i_in_beg:i_in_end,j_node_end-1) - (dyn(:,j_cell_end)**2)*psi(i_in_beg:i_in_end,j_node_end-2))/(dyn(:,j_cell_end-1)**2+2*dyn(:,j_cell_end)*dyn(:,j_cell_end-1))
            ! psi(i_in_beg:i_in_end, j_node_end) = (((dyn(:,j_cell_end)+dyn(:,j_cell_end-1))**2)*psi(i_in_beg:i_in_end,j_node_end-1) - (dyn(:,j_cell_end)**2)*psi(i_in_beg:i_in_end,j_node_end-2))/(dyn(:,j_cell_end-1)**2+2*dyn(:,j_cell_end)*dyn(:,j_cell_end-1))
            psi(i_in_beg:i_in_end, j_node_end) = 0

            ! Calculate errors
            e = relative_residual(psi, psi0, 1.0d-8)

            ! Calculate elapsed cpu time
            call cpu_time(cputime_finish)
            time_cpu = cputime_finish - cputime_start

            ! Print results in screen
            print_screen = 1
            if (mod(Iter_count,print_screen)==0) then
                write(*,'(A,I11,A15,ES9.2,A16,F9.5,A4)') 'Iter.', Iter_count, 'Error =', e, 'CPU time =',time_cpu,'sec'
            end if

            ! Check if fully developed flow
            if (e <= 1.0d-6) then
                write(*,'(A)') 'Fully developed flow'
                call csvwrite1D('psi.csv', convert2Dto1D(psi,NumberNodes_x,NumberNodes_y))
                exit
            else
                psi0 = psi
            end if
        end do


        end program main
        