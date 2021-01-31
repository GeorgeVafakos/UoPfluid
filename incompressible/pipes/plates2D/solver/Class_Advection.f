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
        !   File name:          Class_Advection
        !
        !   Type:               module / class
        !
        !   Description:        In this file the geometry and the grid is generated according to the input file parameters.
        !
        !



        module Class_Advection
            use global_variables
            use Class_Geometry
            use Class_Scalar_Equation
            use Class_Scalar_Variable
            use Class_Vector_Equation
            use Class_Vector_Variable
            implicit none

            type :: Advection
                real*8, allocatable, dimension(:) :: aT
                real*8, allocatable, dimension(:) :: aN
                real*8, allocatable, dimension(:) :: aE
                real*8, allocatable, dimension(:) :: aP
                real*8, allocatable, dimension(:) :: aW
                real*8, allocatable, dimension(:) :: aS
                real*8, allocatable, dimension(:) :: aB
                procedure(CD_scheme), pointer :: create_terms => null()
            contains
                procedure :: allocate => allocate_advection_terms
            end type

            type :: Interpolate_Procedures
            contains
                procedure :: t => interpolate_t
                procedure :: n => interpolate_n
                procedure :: e => interpolate_e
                procedure :: w => interpolate_w
                procedure :: s => interpolate_s
                procedure :: b => interpolate_b
            end type

            ! Declare interpolate precedure class
            type (Interpolate_Procedures) :: Interpolate


        contains

            subroutine allocate_advection_terms(Adv)
                class (Advection) :: Adv

                allocate(Adv%aT(TotalCells), Adv%aN(TotalCells), Adv%aE(TotalCells), Adv%aP(TotalCells), Adv%aW(TotalCells), Adv%aS(TotalCells), Adv%aB(TotalCells))
            end subroutine


            function interpolate_t(Interp, Var)  result(m)
                class (Interpolate_Procedures) :: Interp
                class (Vector_Variable) :: Var
                real*8, dimension(TotalCells) :: m

                m = Var%z(nodes_T)*(0.5*Dz)/dzt + Var%z(nodes_P)*(1-0.5*Dz/dzt)

            end function

            function interpolate_n(Interp, Var)  result(m)
                class (Interpolate_Procedures) :: Interp
                class (Vector_Variable) :: Var
                real*8, dimension(TotalCells) :: m

                m = Var%y(nodes_N)*(0.5*Dy)/dyn + Var%y(nodes_P)*(1-0.5*Dy/dyn)

            end function

            function interpolate_e(Interp, Var)  result(m)
                class (Interpolate_Procedures) :: Interp
                class (Vector_Variable) :: Var
                real*8, dimension(TotalCells) :: m

                m = Var%x(nodes_E)*(0.5*Dx)/dxe + Var%x(nodes_P)*(1-0.5*Dx/dxe)

            end function

            function interpolate_w(Interp, Var)  result(m)
                class (Interpolate_Procedures) :: Interp
                class (Vector_Variable) :: Var
                real*8, dimension(TotalCells) :: m

                m = Var%x(nodes_W)*(0.5*Dx)/dxw + Var%x(nodes_P)*(1-0.5*Dx/dxw)

            end function

            function interpolate_s(Interp, Var)  result(m)
                class (Interpolate_Procedures) :: Interp
                class (Vector_Variable) :: Var
                real*8, dimension(TotalCells) :: m

                m = Var%y(nodes_S)*(0.5*Dy)/dys + Var%y(nodes_P)*(1-0.5*Dy/dys)

            end function

            function interpolate_b(Interp, Var)  result(m)
                class (Interpolate_Procedures) :: Interp
                class (Vector_Variable) :: Var
                real*8, dimension(TotalCells) :: m

                m = Var%z(nodes_B)*(0.5*Dz)/dzb + Var%z(nodes_P)*(1-0.5*Dz/dzb)
            end function


            subroutine CD_scheme(Adv, Var)
                class (Advection) :: Adv
                class (Vector_Variable) :: Var

                Adv%aT = k_boole*interpolate%t(Var)*0.5*Dx*Dy*Dz/dzt
                Adv%aN = interpolate%n(Var)*0.5*Dx*Dy/dyn*(1.0-k_boole+k_boole*Dz)
                Adv%aE = interpolate%e(Var)*0.5*Dx*Dy/dxe*(1.0-k_boole+k_boole*Dz)
                Adv%aP = k_boole*interpolate%t(Var)*(1.0-0.5*Dz/dzt)*Dx*Dy + interpolate%n(Var)*(1.0-0.5*Dy/dyn)*Dx*(1.0-k_boole+k_boole*Dz) + interpolate%e(Var)*(1.0-0.5*Dx/dxe)*Dy*(1.0-k_boole+k_boole*Dz) - interpolate%w(Var)*(1.0-0.5*Dx/dxw)*Dy*(1.0-k_boole+k_boole*Dz) - interpolate%s(Var)*(1.0-0.5*Dy/dys)*Dx*(1.0-k_boole+k_boole*Dz) - k_boole*interpolate%b(Var)*(1.0-0.5*Dz/dzb)*Dx*Dy
                Adv%aW =-interpolate%w(Var)*0.5*Dx*Dy/dxw*(1.0-k_boole+k_boole*Dz)
                Adv%aS =-interpolate%s(Var)*0.5*Dx*Dy/dys*(1.0-k_boole+k_boole*Dz)
                Adv%aB =-k_boole*interpolate%b(Var)*0.5*Dx*Dy*Dz/dzb

            end subroutine


            subroutine upwind_scheme(Adv, Var)
                class (Advection) :: Adv
                class (Vector_Variable) :: Var

                Adv%aT = 0.0
                Adv%aN = 0.0
                Adv%aE = 0.0
                Adv%aP = 0.0
                Adv%aW = 0.0
                Adv%aS = 0.0
                Adv%aB = 0.0

                where (Var%z(nodes_P)<0)  Adv%aT = k_boole*Interpolate%t(Var)*Dx*Dy
                where (Var%y(nodes_P)<0)  Adv%aN = interpolate%n(Var)*Dx*(1.0-k_boole+k_boole*Dz)
                where (Var%x(nodes_P)<0)  Adv%aE = interpolate%e(Var)*Dy*(1.0-k_boole+k_boole*Dz)

                where (Var%z(nodes_P)>=0) Adv%aP = Adv%aP + k_boole*interpolate%t(Var)*Dx*Dy
                where (Var%y(nodes_P)>=0) Adv%aP = Adv%aP + interpolate%n(Var)*Dx*(1.0-k_boole+k_boole*Dz)
                where (Var%x(nodes_P)>=0) Adv%aP = Adv%aP + interpolate%e(Var)*Dy*(1.0-k_boole+k_boole*Dz)
                where (Var%x(nodes_P)<0)  Adv%aP = Adv%aP - interpolate%w(Var)*Dy*(1.0-k_boole+k_boole*Dz)
                where (Var%y(nodes_P)<0)  Adv%aP = Adv%aP - interpolate%s(Var)*Dx*(1.0-k_boole+k_boole*Dz)
                where (Var%z(nodes_P)<0)  Adv%aP = Adv%aP - k_boole*Interpolate%b(Var)*Dx*Dy

                where (Var%x(nodes_P)>=0) Adv%aW =-interpolate%w(Var)*Dy*(1.0-k_boole+k_boole*Dz)
                where (Var%y(nodes_P)>=0) Adv%aS =-interpolate%s(Var)*Dx*(1.0-k_boole+k_boole*Dz)
                where (Var%z(nodes_P)>=0) Adv%aB =-k_boole*Interpolate%b(Var)*Dx*Dy
            end subroutine


            subroutine none(Adv, Var)
                class (Advection) :: Adv
                class (Vector_Variable) :: Var

                Adv%aT = 0.0
                Adv%aN = 0.0
                Adv%aE = 0.0
                Adv%aP = 0.0
                Adv%aW = 0.0
                Adv%aS = 0.0
                Adv%aB = 0.0

            end subroutine

        end module