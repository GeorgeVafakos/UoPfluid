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
        !   File name:          mesh_generation
        !
        !   Type:               source
        !
        !   Description:        In this file the geometry and the grid is generated according to the input file parameters.
        !
        !
        !------------------------------------------------------------------------------------------
        ! Define Axis
        !------------------------------------------------------------------------------------------
        
        ! Define axis indices
        call x%define(cells_x, x_begin, x_end)
        call y%define(cells_y, y_begin, y_end)
        if (flow_2D .eqv. .FALSE.)  then
            call z%define(cells_z, z_begin, z_end)
            call N%define(x, y, z)
        else
            call N%define(x, y)
        end if

        ! Define axis faces and nodes
        call x%create_faces
        call x%create_nodes
        call y%create_faces
        call y%create_nodes
        if (flow_2D .eqv. .FALSE.)  then
            call z%create_faces
            call z%create_nodes
        end if

        ! Define total faces, nodes and cells
        if (flow_2D .eqv. .FALSE.)  then
            TotalFaces = x%NumberFaces*y%NumberFaces*z%NumberFaces
            TotalNodes = x%NumberNodes*y%NumberNodes*z%NumberNodes
            TotalCells = x%NumberCells*y%NumberCells*z%NumberCells
        else
            TotalFaces = x%NumberFaces*y%NumberFaces
            TotalNodes = x%NumberNodes*y%NumberNodes
            TotalCells = x%NumberCells*y%NumberCells
        end if


        !------------------------------------------------------------------------------------------
        ! Create Mesh Distances
        !------------------------------------------------------------------------------------------
        Dx_vec = x%faces(x%face_beg+1:x%face_end) - x%faces(x%face_beg:x%face_end-1)
        Dy_vec = y%faces(y%face_beg+1:y%face_end) - y%faces(y%face_beg:y%face_end-1)

        dxe_vec = x%nodes(x%in_beg+1:x%in_end+1) - x%nodes(x%in_beg:x%in_end)
        dxw_vec = x%nodes(x%in_beg:x%in_end) - x%nodes(x%in_beg-1:x%in_end-1)
        
        dyn_vec = y%nodes(y%in_beg+1:y%in_end+1) - y%nodes(y%in_beg:y%in_end)
        dys_vec = y%nodes(y%in_beg:y%in_end) - y%nodes(y%in_beg-1:y%in_end-1)

        if (flow_2D .eqv. .FALSE.)  then
            Dz_vec = z%faces(z%face_beg+1:z%face_end) - z%faces(z%face_beg:z%face_end-1)
            dzt_vec = z%nodes(z%in_beg+1:z%in_end+1) - z%nodes(z%in_beg:z%in_end)
            dzb_vec = z%nodes(z%in_beg:z%in_end) - z%nodes(z%in_beg-1:z%in_end-1)
        end if 

        if (flow_2D .eqv. .TRUE.)   then
            allocate(Dx_2D(cells_x,cells_y), Dy_2D(cells_x,cells_y), Dz_2D(cells_x,cells_y), dxe_2D(cells_x,cells_y), dxw_2D(cells_x,cells_y), dyn_2D(cells_x,cells_y), dys_2D(cells_x,cells_y), dzt_2D(cells_x,cells_y), dzb_2D(cells_x,cells_y))
        else
            allocate(Dx_3D(cells_x,cells_y,cells_z), Dy_3D(cells_x,cells_y,cells_z), Dz_3D(cells_x,cells_y,cells_z), dxe_3D(cells_x,cells_y,cells_z), dxw_3D(cells_x,cells_y,cells_z), dyn_3D(cells_x,cells_y,cells_z), dys_3D(cells_x,cells_y,cells_z), dzt_3D(cells_x,cells_y,cells_z), dzb_3D(cells_x,cells_y,cells_z))
        end if


        ! Create a distance vector that consits of smaller identical vectors, placed one after another
        allocate(Dx(TotalCells), Dy(TotalCells), Dz(TotalCells), dxe(TotalCells), dxw(TotalCells), dyn(TotalCells), dys(TotalCells), dzt(TotalCells), dzb(TotalCells))
        if (flow_2D .eqv. .FALSE.)  then

            call meshgrid(Dx_3D, Dy_3D, Dz_3D, Dx_vec, Dy_vec, Dz_vec)
            call meshgrid(dxe_3D, dyn_3D, dzt_3D, dxe_vec, dyn_vec, dzt_vec)
            call meshgrid(dxw_3D, dys_3D, dzb_3D, dxw_vec, dys_vec, dzb_vec)

            ! Dx = pack(spread(Dx_vec,1,cells_y*cells_z),.TRUE.)
            ! Dy = pack(spread(Dy_vec,1,cells_x*cells_z),.TRUE.)
            ! Dz = pack(spread(Dz_vec,1,cells_x*cells_y),.TRUE.)
            
            ! dxe = pack(spread(dxe_vec,2,cells_y*cells_z),.TRUE.)
            ! dxw = pack(spread(dxw_vec,2,cells_y*cells_z),.TRUE.)

            ! dyn = pack(spread(dyn_vec,2,cells_x*cells_z),.TRUE.)
            ! dys = pack(spread(dys_vec,2,cells_x*cells_z),.TRUE.)

            ! dzt = pack(spread(dzt_vec,1,cells_x*cells_y),.TRUE.)
            ! dzb = pack(spread(dzb_vec,1,cells_x*cells_y),.TRUE.)

            Dx = convert3Dto1D(Dx_3D, cells_x, cells_y, cells_z)
            Dy = convert3Dto1D(Dy_3D, cells_x, cells_y, cells_z)
            Dz = convert3Dto1D(Dz_3D, cells_x, cells_y, cells_z)

            dxe = convert3Dto1D(dxe_3D, cells_x, cells_y, cells_z)
            dxw = convert3Dto1D(dxw_3D, cells_x, cells_y, cells_z)
            dyn = convert3Dto1D(dyn_3D, cells_x, cells_y, cells_z)
            dys = convert3Dto1D(dys_3D, cells_x, cells_y, cells_z)
            dzt = convert3Dto1D(dzt_3D, cells_x, cells_y, cells_z)
            dzb = convert3Dto1D(dzb_3D, cells_x, cells_y, cells_z)

            ! dxe = pack(dxe_3D,.TRUE.)
            ! dxw = pack(dxw_3D,.TRUE.)
            ! dyn = pack(dyn_3D,.TRUE.)
            ! dys = pack(dys_3D,.TRUE.)
            ! dzt = pack(dzt_3D,.TRUE.)
            ! dzb = pack(dzb_3D,.TRUE.)
        else
            Dx = pack(spread(Dx_vec,2,cells_y),.TRUE.)
            Dy = pack(spread(Dy_vec,1,cells_x),.TRUE.)
            dxe = pack(spread(dxe_vec,2,cells_y),.TRUE.)
            dxw = pack(spread(dxw_vec,2,cells_y),.TRUE.)
            dyn = pack(spread(dyn_vec,1,cells_x),.TRUE.)
            dys = pack(spread(dys_vec,1,cells_x),.TRUE.)

            Dz = 1.0
            dzt = 1.0
            dzb = 1.0
        end if

        !------------------------------------------------------------------------------------------
        ! Define global index parameters
        !------------------------------------------------------------------------------------------
        call domain%define
        call domain%define_BC_distances

        ! Define a globar variable with the number of boundaries
        NumberBoundaries = domain%NumberBoundaries

        ! Define top, north, east, west, south and bottom nodes from polar
        nodes_T = domain%nodes_inner + k_boole*x%NumberNodes*y%NumberNodes
        nodes_N = domain%nodes_inner + x%NumberNodes
        nodes_E = domain%nodes_inner + 1
        nodes_P = domain%nodes_inner
        nodes_W = domain%nodes_inner - 1
        nodes_S = domain%nodes_inner - x%NumberNodes
        nodes_B = domain%nodes_inner - k_boole*x%NumberNodes*y%NumberNodes

        ! Associate cell global index (ID) with node global index (glob_index)
        allocate(glob_index(TotalCells))
        ID = 0
        if (flow_2D .eqv. .FALSE.)  then
            do k = z%node_beg, z%node_end
                do j = y%node_beg, y%node_end
                    do i = x%node_beg, x%node_end
                        if (i/=x%node_beg .AND. i/=x%node_end .AND. j/=y%node_beg .AND. j/=y%node_end .AND. k/=z%node_beg .AND. k/=z%node_end) then
                            ID = ID + 1
                            glob_index(ID) = i+(j-1)*x%NumberNodes+(k-1)*x%NumberNodes*y%NumberNodes
                        end if
                    end do
                end do
            end do
        else
            do j = y%node_beg, y%node_end
                do i = x%node_beg, x%node_end
                    if (i/=x%node_beg .AND. i/=x%node_end .AND. j/=y%node_beg .AND. j/=y%node_end) then
                        ID = ID + 1
                        glob_index(ID) = i+(j-1)*x%NumberNodes
                    end if
                end do
            end do
        end if

        !------------------------------------------------------------------------------------------
        ! Create mesh folder and plot
        !------------------------------------------------------------------------------------------
        
        call system("rm -rf ../mesh")
        call system("mkdir -p ../mesh")
        call chdir('../mesh')

        ! Print nodes and faces
        call csvwrite('x_nodes.csv', x%nodes)
        call csvwrite('x_faces.csv', x%faces)
        call csvwrite('y_nodes.csv', y%nodes)
        call csvwrite('y_faces.csv', y%faces)
        if (flow_2D .eqv. .FALSE.)  then
            call csvwrite('z_nodes.csv', z%nodes)
            call csvwrite('z_faces.csv', z%faces)
        end if 
        
        call system('cp ../plot_grid.py plots.py')
        call chdir('../solver')