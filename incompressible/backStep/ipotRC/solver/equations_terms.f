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
        !   File name:          equations_terms
        ! 
        !   Type:               source
        ! 
        !   Description:        This is the file where all finite volume coefficients for the Navier-Stokes and momentum equations terms are being calculated.
        ! 
        ! 
        !------------------------------------------------------------------------------------------
        ! Finite Volume Terms
        !------------------------------------------------------------------------------------------

        ! Time derivative
        aP_time = (Dx*Dy/Dt)*(1.0-k_boole+k_boole*Dz)

        ! Laplacian terms
        aT_lapl = k_boole*Dx*Dy/dzt
        aN_lapl = (Dx/dyn)*(1.0-k_boole+k_boole*Dz)
        aE_lapl = (Dy/dxe)*(1.0-k_boole+k_boole*Dz)
        aP_lapl =-(Dy/dxe + Dy/dxw + Dx/dyn + Dx/dys)*(1.0-k_boole+k_boole*Dz) - k_boole*(Dx*Dy/dzt + Dx*Dy/dzb)
        aW_lapl = (Dy/dxw)*(1.0-k_boole+k_boole*Dz)
        aS_lapl = (Dx/dys)*(1.0-k_boole+k_boole*Dz)
        aB_lapl = k_boole*Dx*Dy/dzb

        ! Gradient terms
        aE_gradx = (0.5*Dx*Dy/dxe)*(1.0-k_boole+k_boole*Dz)
        aP_gradx = (0.5*Dx*Dy*(1.0/dxw - 1.0/dxe))*(1.0-k_boole+k_boole*Dz)
        aW_gradx =-(0.5*Dx*Dy/dxw)*(1.0-k_boole+k_boole*Dz)

        aN_grady = (0.5*Dx*Dy/dyn)*(1.0-k_boole+k_boole*Dz)
        aP_grady = (0.5*Dx*Dy*(1.0/dys - 1.0/dyn))*(1.0-k_boole+k_boole*Dz)
        aS_grady =-(0.5*Dx*Dy/dys)*(1.0-k_boole+k_boole*Dz)

        aT_gradz = k_boole*0.5*Dx*Dy*Dz/dzt
        aP_gradz = k_boole*0.5*Dx*Dy*Dz*(1.0/dzb - 1.0/dzt)
        aB_gradz =-k_boole*0.5*Dx*Dy*Dz/dzb


        ! Divergence Terms
        aE_divx = 0.5*Dx*Dy/dxe*(1.0-k_boole+k_boole*Dz)
        aP_divx = 0.5*Dx*Dy*(1.0/dxw - 1.0/dxe)*(1.0-k_boole+k_boole*Dz)
        aW_divx =-0.5*Dx*Dy/dxw*(1.0-k_boole+k_boole*Dz)

        aN_divy = 0.5*Dx*Dy/dyn*(1.0-k_boole+k_boole*Dz)
        aP_divy = 0.5*Dx*Dy*(1.0/dys - 1.0/dyn)*(1.0-k_boole+k_boole*Dz)
        aS_divy =-0.5*Dx*Dy/dys*(1.0-k_boole+k_boole*Dz)

        aT_divz = 0.5*Dx*Dy*Dz/dzt
        aP_divz = 0.5*Dx*Dy*Dz*(1.0/dzb - 1.0/dzt)
        aB_divz =-0.5*Dx*Dy*Dz/dzb

        !------------------------------------------------------------------------------------------
        ! Navier-Stokes Equation
        !------------------------------------------------------------------------------------------
        ! Advection terms
        call NS_adv_term%create_terms(V)

        ! Navier-Stokes x-axis terms
        NS_eqn%x%aT = -nu*aT_lapl + NS_adv_term%aT
        NS_eqn%x%aN = -nu*aN_lapl + NS_adv_term%aN
        NS_eqn%x%aE = -nu*aE_lapl + NS_adv_term%aE
        NS_eqn%x%aP = -nu*aP_lapl + NS_adv_term%aP + aP_time
        NS_eqn%x%aW = -nu*aW_lapl + NS_adv_term%aW
        NS_eqn%x%aS = -nu*aS_lapl + NS_adv_term%aS
        NS_eqn%x%aB = -nu*aB_lapl + NS_adv_term%aB
        NS_eqn%x%B  =  aP_time*V%x(nodes_P) - (1.0/rho)*gradx(p%field)

        ! Navier-Stokes y-axis terms
        NS_eqn%y%aT = -nu*aT_lapl + NS_adv_term%aT
        NS_eqn%y%aN = -nu*aN_lapl + NS_adv_term%aN
        NS_eqn%y%aE = -nu*aE_lapl + NS_adv_term%aE
        NS_eqn%y%aP = -nu*aP_lapl + NS_adv_term%aP + aP_time
        NS_eqn%y%aW = -nu*aW_lapl + NS_adv_term%aW
        NS_eqn%y%aS = -nu*aS_lapl + NS_adv_term%aS
        NS_eqn%y%aB = -nu*aB_lapl + NS_adv_term%aB
        NS_eqn%y%B  =  aP_time*V%y(nodes_P) - (1.0/rho)*grady(p%field)

        ! Navier-Stokes z-axis terms
        if (flow_2D .eqv. .FALSE.)  then
            NS_eqn%z%aT = -nu*aT_lapl + NS_adv_term%aT
            NS_eqn%z%aN = -nu*aN_lapl + NS_adv_term%aN
            NS_eqn%z%aE = -nu*aE_lapl + NS_adv_term%aE
            NS_eqn%z%aP = -nu*aP_lapl + NS_adv_term%aP + aP_time
            NS_eqn%z%aW = -nu*aW_lapl + NS_adv_term%aW
            NS_eqn%z%aS = -nu*aS_lapl + NS_adv_term%aS
            NS_eqn%z%aB = -nu*aB_lapl + NS_adv_term%aB
            NS_eqn%z%B  =  aP_time*V%z(nodes_P) - (1.0/rho)*gradz(p%field)
        end if


        !------------------------------------------------------------------------------------------
        ! Pressure Equation
        !------------------------------------------------------------------------------------------
        ! Poisson terms
        Pres_eqn%aT = (Dx*Dy*(1.0-k_boole+k_boole*Dz))*aT_lapl
        Pres_eqn%aN = (Dx*Dy*(1.0-k_boole+k_boole*Dz))*aN_lapl
        Pres_eqn%aE = (Dx*Dy*(1.0-k_boole+k_boole*Dz))*aE_lapl
        Pres_eqn%aP = (Dx*Dy*(1.0-k_boole+k_boole*Dz))*aP_lapl
        Pres_eqn%aW = (Dx*Dy*(1.0-k_boole+k_boole*Dz))*aW_lapl
        Pres_eqn%aS = (Dx*Dy*(1.0-k_boole+k_boole*Dz))*aS_lapl
        Pres_eqn%aB = (Dx*Dy*(1.0-k_boole+k_boole*Dz))*aB_lapl
        
        
