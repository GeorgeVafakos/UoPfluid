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
        !   File name:          Class_Geometry
        !
        !   Type:               module / class
        !
        !   Description:        In this file the geometry class is created. 
        !
        !

        module Class_Geometry
            use global_variables

            implicit none

            type :: boundary
                integer, allocatable, dimension(:)  :: nodes_bound, cells_bound, cells_inner
                integer, allocatable, dimension(:)  :: nodes_next1, nodes_next2
                real,  allocatable, dimension(:)  :: nodes_dist1, nodes_dist2
            end type

            type :: duct
                integer :: NumberBoundaries
                integer, allocatable, dimension(:) :: NumberNodesBoundary
                type (boundary), allocatable, dimension(:) :: boundary
                integer, allocatable, dimension(:)  :: nodes_inner
            contains
                procedure :: define => define_global_nodes_and_cells
                procedure :: define_BC_distances => define_distances_near_boundaries
                procedure :: allocate_bound_number => allocate_number_of_boundaries
            end type


        contains

            subroutine allocate_number_of_boundaries(domain)
                class (duct) :: domain

                if (flow_2D .eqv. .TRUE.)   then
                    domain%NumberBoundaries = 4
                else 
                    domain%NumberBoundaries = 6
                end if
            end subroutine


            subroutine define_global_nodes_and_cells(domain)
                class (duct) :: domain
                type (duct) :: temp

                if (flow_2D .eqv. .TRUE.)   then
                    
                    allocate(domain%NumberNodesBoundary(domain%NumberBoundaries))
                    allocate(domain%boundary(domain%NumberBoundaries), temp%boundary(domain%NumberBoundaries))

                    do i = 1, domain%NumberBoundaries
                        allocate(temp%boundary(i)%nodes_bound(TotalNodes), temp%boundary(i)%cells_bound(TotalCells), temp%boundary(i)%cells_inner(TotalCells))
                    end do
                    allocate(temp%nodes_inner(TotalNodes))

                    ! Global nodes at the boundaries
                    do i = 1, domain%NumberBoundaries
                        temp%boundary(i)%nodes_bound = 0
                    end do
                    temp%nodes_inner = 0
                    do j = y%node_beg, y%node_end
                        do i = x%node_beg, x%node_end
                            ID = i+(j-1)*x%NumberNodes
                            if (i==x%node_beg .AND. j/=y%node_beg .AND. j/=y%node_end)                          temp%boundary(1)%nodes_bound(ID) = ID
                            if (i==x%node_end .AND. j/=y%node_beg .AND. j/=y%node_end)                          temp%boundary(2)%nodes_bound(ID) = ID
                            if (j==y%node_end .AND. i/=x%node_beg .AND. i/=x%node_end)                          temp%boundary(3)%nodes_bound(ID) = ID
                            if (j==y%node_beg .AND. i/=x%node_beg .AND. i/=x%node_end)                          temp%boundary(4)%nodes_bound(ID) = ID
                            if (i/=x%node_beg .AND. i/=x%node_end .AND. j/=y%node_beg .AND. j/=y%node_end)      temp%nodes_inner(ID) = ID
                        end do
                    end do


                    do i = 1, domain%NumberBoundaries
                        allocate(domain%boundary(i)%nodes_bound(size(pack(temp%boundary(i)%nodes_bound,temp%boundary(i)%nodes_bound/=0))))
                        domain%boundary(i)%nodes_bound = pack(temp%boundary(i)%nodes_bound, temp%boundary(i)%nodes_bound/=0)
                        domain%NumberNodesBoundary(i) = size(domain%boundary(i)%nodes_bound)
                        deallocate(temp%boundary(i)%nodes_bound)
                    end do
                    domain%nodes_inner = pack(temp%nodes_inner, temp%nodes_inner/=0)
                    deallocate(temp%nodes_inner)


                    ! Global cells at the boundary
                    do i = 1, domain%NumberBoundaries
                        temp%boundary(i)%cells_bound = 0
                        temp%boundary(i)%cells_inner = 0
                    end do
                    do j = y%cell_beg, y%cell_end
                        do i = x%cell_beg, x%cell_end

                            ID = i+(j-1)*x%NumberCells
                            ! Left boundary
                            if (i==x%cell_beg) then
                                temp%boundary(1)%cells_bound(ID) = ID
                            else
                                temp%boundary(1)%cells_inner(ID) = ID
                            end if
                            ! Right boundary
                            if (i==x%cell_end) then
                                temp%boundary(2)%cells_bound(ID) = ID
                            else
                                temp%boundary(2)%cells_inner(ID) = ID
                            end if
                            ! Top boundary
                            if (j==y%cell_end) then
                                temp%boundary(3)%cells_bound(ID) = ID
                            else
                                temp%boundary(3)%cells_inner(ID) = ID
                            end if
                            ! Bottom boundary
                            if (j==y%cell_beg) then
                                temp%boundary(4)%cells_bound(ID) = ID
                            else
                                temp%boundary(4)%cells_inner(ID) = ID
                            end if

                        end do
                    end do

                    do i = 1, domain%NumberBoundaries
                        allocate(domain%boundary(i)%cells_bound(size(pack(temp%boundary(i)%cells_bound,temp%boundary(i)%cells_bound/=0))), domain%boundary(i)%cells_inner(size(pack(temp%boundary(i)%cells_inner,temp%boundary(i)%cells_inner/=0))))
                        domain%boundary(i)%cells_bound = pack(temp%boundary(i)%cells_bound, temp%boundary(i)%cells_bound/=0)
                        domain%boundary(i)%cells_inner = pack(temp%boundary(i)%cells_inner, temp%boundary(i)%cells_inner/=0)
                        deallocate(temp%boundary(i)%cells_bound, temp%boundary(i)%cells_inner)
                    end do

                else
                    allocate(domain%NumberNodesBoundary(domain%NumberBoundaries))
                    allocate(domain%boundary(domain%NumberBoundaries), temp%boundary(domain%NumberBoundaries))

                    do i = 1, domain%NumberBoundaries
                        allocate(temp%boundary(i)%nodes_bound(TotalNodes), temp%boundary(i)%cells_bound(TotalCells), temp%boundary(i)%cells_inner(TotalCells))
                    end do
                    allocate(temp%nodes_inner(TotalNodes))

                    ! Global nodes at the boundaries
                    do i = 1, domain%NumberBoundaries
                        temp%boundary(i)%nodes_bound = 0
                    end do
                    temp%nodes_inner = 0
                    do k = z%node_beg, z%node_end
                        do j = y%node_beg, y%node_end
                            do i = x%node_beg, x%node_end
                                ID = i+(j-1)*x%NumberNodes+(k-1)*x%NumberNodes*y%NumberNodes
                                if (k==z%node_beg .AND. i/=x%node_beg .AND. i/=x%node_end .AND. j/=y%node_beg .AND. j/=y%node_end)                          temp%boundary(1)%nodes_bound(ID) = ID
                                if (k==z%node_end .AND. i/=x%node_beg .AND. i/=x%node_end .AND. j/=y%node_beg .AND. j/=y%node_end)                          temp%boundary(2)%nodes_bound(ID) = ID
                                if (j==y%node_end .AND. i/=x%node_beg .AND. i/=x%node_end .AND. k/=z%node_beg .AND. k/=z%node_end)                          temp%boundary(3)%nodes_bound(ID) = ID
                                if (j==y%node_beg .AND. i/=x%node_beg .AND. i/=x%node_end .AND. k/=z%node_beg .AND. k/=z%node_end)                          temp%boundary(4)%nodes_bound(ID) = ID
                                if (i==x%node_beg .AND. j/=y%node_beg .AND. j/=y%node_end .AND. k/=z%node_beg .AND. k/=z%node_end)                          temp%boundary(5)%nodes_bound(ID) = ID
                                if (i==x%node_end .AND. j/=y%node_beg .AND. j/=y%node_end .AND. k/=z%node_beg .AND. k/=z%node_end)                          temp%boundary(6)%nodes_bound(ID) = ID
                                if (i/=x%node_beg .AND. i/=x%node_end .AND. j/=y%node_beg .AND. j/=y%node_end .AND. k/=z%node_beg .AND. k/=z%node_end)      temp%nodes_inner(ID) = ID
                            end do
                        end do
                    end do


                    do i = 1, domain%NumberBoundaries
                        allocate(domain%boundary(i)%nodes_bound(size(pack(temp%boundary(i)%nodes_bound,temp%boundary(i)%nodes_bound/=0))))
                        domain%boundary(i)%nodes_bound = pack(temp%boundary(i)%nodes_bound, temp%boundary(i)%nodes_bound/=0)
                        domain%NumberNodesBoundary(i) = size(domain%boundary(i)%nodes_bound)
                        deallocate(temp%boundary(i)%nodes_bound)
                    end do
                    domain%nodes_inner = pack(temp%nodes_inner, temp%nodes_inner/=0)
                    deallocate(temp%nodes_inner)


                    ! Global cells at the boundary
                    do i = 1, domain%NumberBoundaries
                        temp%boundary(i)%cells_bound = 0
                        temp%boundary(i)%cells_inner = 0
                    end do
                    do k = z%cell_beg, z%cell_end
                        do j = y%cell_beg, y%cell_end
                            do i = x%cell_beg, x%cell_end

                                ID = i+(j-1)*x%NumberCells+(k-1)*x%NumberCells*y%NumberCells
                                ! Left boundary
                                if (k==z%cell_beg) then
                                    temp%boundary(1)%cells_bound(ID) = ID
                                else
                                    temp%boundary(1)%cells_inner(ID) = ID
                                end if
                                ! Right boundary
                                if (k==z%cell_end) then
                                    temp%boundary(2)%cells_bound(ID) = ID
                                else
                                    temp%boundary(2)%cells_inner(ID) = ID
                                end if
                                ! Top boundary
                                if (j==y%cell_end) then
                                    temp%boundary(3)%cells_bound(ID) = ID
                                else
                                    temp%boundary(3)%cells_inner(ID) = ID
                                end if
                                ! Bottom boundary
                                if (j==y%cell_beg) then
                                    temp%boundary(4)%cells_bound(ID) = ID
                                else
                                    temp%boundary(4)%cells_inner(ID) = ID
                                end if
                                ! Front boundary
                                if (i==x%cell_beg) then
                                    temp%boundary(5)%cells_bound(ID) = ID
                                else
                                    temp%boundary(5)%cells_inner(ID) = ID
                                end if
                                ! Back boundary
                                if (i==x%cell_end) then
                                    temp%boundary(6)%cells_bound(ID) = ID
                                else
                                    temp%boundary(6)%cells_inner(ID) = ID
                                end if

                            end do
                        end do
                    end do

                    do i = 1, domain%NumberBoundaries
                        allocate(domain%boundary(i)%cells_bound(size(pack(temp%boundary(i)%cells_bound,temp%boundary(i)%cells_bound/=0))), domain%boundary(i)%cells_inner(size(pack(temp%boundary(i)%cells_inner,temp%boundary(i)%cells_inner/=0))))
                        domain%boundary(i)%cells_bound = pack(temp%boundary(i)%cells_bound, temp%boundary(i)%cells_bound/=0)
                        domain%boundary(i)%cells_inner = pack(temp%boundary(i)%cells_inner, temp%boundary(i)%cells_inner/=0)
                        deallocate(temp%boundary(i)%cells_bound, temp%boundary(i)%cells_inner)
                    end do
                end if

            end subroutine


            subroutine define_distances_near_boundaries(domain)
                class (duct) :: domain

                if (flow_2D .eqv. .TRUE.)   then
                    ! Left boundary
                    domain%boundary(1)%nodes_next1 = domain%boundary(1)%nodes_bound + 1
                    domain%boundary(1)%nodes_next2 = domain%boundary(1)%nodes_bound + 2
                    domain%boundary(1)%nodes_dist1 = dxw(domain%boundary(1)%cells_bound)
                    domain%boundary(1)%nodes_dist2 = dxw(domain%boundary(1)%cells_bound + 1)

                    ! Right boundary
                    domain%boundary(2)%nodes_next1 = domain%boundary(2)%nodes_bound - 1
                    domain%boundary(2)%nodes_next2 = domain%boundary(2)%nodes_bound - 2
                    domain%boundary(2)%nodes_dist1 = dxe(domain%boundary(2)%cells_bound)
                    domain%boundary(2)%nodes_dist2 = dxe(domain%boundary(2)%cells_bound - 1)

                    ! Top boundary
                    domain%boundary(3)%nodes_next1 = domain%boundary(3)%nodes_bound - x%NumberNodes
                    domain%boundary(3)%nodes_next2 = domain%boundary(3)%nodes_bound - 2*x%NumberNodes
                    domain%boundary(3)%nodes_dist1 = dyn(domain%boundary(3)%cells_bound)
                    domain%boundary(3)%nodes_dist2 = dyn(domain%boundary(3)%cells_bound - x%NumberCells)

                    ! Bottom boundary
                    domain%boundary(4)%nodes_next1 = domain%boundary(4)%nodes_bound + x%NumberNodes
                    domain%boundary(4)%nodes_next2 = domain%boundary(4)%nodes_bound + 2*x%NumberNodes
                    domain%boundary(4)%nodes_dist1 = dys(domain%boundary(4)%cells_bound)
                    domain%boundary(4)%nodes_dist2 = dys(domain%boundary(4)%cells_bound + x%NumberCells)
                    
                else
                    ! Left boundary
                    domain%boundary(1)%nodes_next1 = domain%boundary(1)%nodes_bound + x%NumberNodes*y%NumberNodes
                    domain%boundary(1)%nodes_next2 = domain%boundary(1)%nodes_bound + 2*x%NumberNodes*y%NumberNodes
                    domain%boundary(1)%nodes_dist1 = dzb(domain%boundary(1)%cells_bound)
                    domain%boundary(1)%nodes_dist2 = dzb(domain%boundary(1)%cells_bound + x%NumberCells*y%NumberCells)

                    ! Right boundary
                    domain%boundary(2)%nodes_next1 = domain%boundary(2)%nodes_bound - x%NumberNodes*y%NumberNodes
                    domain%boundary(2)%nodes_next2 = domain%boundary(2)%nodes_bound - 2*x%NumberNodes*y%NumberNodes
                    domain%boundary(2)%nodes_dist1 = dzt(domain%boundary(2)%cells_bound)
                    domain%boundary(2)%nodes_dist2 = dzt(domain%boundary(2)%cells_bound - x%NumberCells*y%NumberCells)

                    ! Top boundary
                    domain%boundary(3)%nodes_next1 = domain%boundary(3)%nodes_bound - x%NumberNodes
                    domain%boundary(3)%nodes_next2 = domain%boundary(3)%nodes_bound - 2*x%NumberNodes
                    domain%boundary(3)%nodes_dist1 = dyn(domain%boundary(3)%cells_bound)
                    domain%boundary(3)%nodes_dist2 = dyn(domain%boundary(3)%cells_bound - x%NumberCells)

                    ! Bottom boundary
                    domain%boundary(4)%nodes_next1 = domain%boundary(4)%nodes_bound + x%NumberNodes
                    domain%boundary(4)%nodes_next2 = domain%boundary(4)%nodes_bound + 2*x%NumberNodes
                    domain%boundary(4)%nodes_dist1 = dys(domain%boundary(4)%cells_bound)
                    domain%boundary(4)%nodes_dist2 = dys(domain%boundary(4)%cells_bound + x%NumberCells)

                    ! Front boundary
                    domain%boundary(5)%nodes_next1 = domain%boundary(5)%nodes_bound + 1
                    domain%boundary(5)%nodes_next2 = domain%boundary(5)%nodes_bound + 2
                    domain%boundary(5)%nodes_dist1 = dxw(domain%boundary(5)%cells_bound)
                    domain%boundary(5)%nodes_dist2 = dxw(domain%boundary(5)%cells_bound + 1)

                    ! Back boundary
                    domain%boundary(6)%nodes_next1 = domain%boundary(6)%nodes_bound - 1
                    domain%boundary(6)%nodes_next2 = domain%boundary(6)%nodes_bound - 2
                    domain%boundary(6)%nodes_dist1 = dxe(domain%boundary(6)%cells_bound)
                    domain%boundary(6)%nodes_dist2 = dxe(domain%boundary(6)%cells_bound - 1)

                end if

            end subroutine


        end module



