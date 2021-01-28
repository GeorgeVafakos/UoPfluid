        !-----------------------------------------------------------------------------------------------------------------------------------------------------------
        !                                                                                     !
        !    ______    ______        _________                                                !   Version: 2.0
        !      ||        ||            ||    \\     ___                             ____      !
        !      ||        ||            ||     ||   // \\ ==||                        ||       !
        !      ||        ||            ||     ||  ||  //   ||                        ||       !
        !      ||        ||    ____    ||____//   ||       ||  ____  ____   ο    ____||       !
        !      ||        ||   //  \\   ||       ==||==     ||   ||    ||  =||   //   ||       !
        !      ||        ||  ||    ||  ||         ||       ||   ||    ||   ||  ||    ||       !
        !       \\______//    \\__// __||__     __||__   __||__  \\__//   _||_  \\___||_      !
        !                                                                                     !
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
        ! call NS_max_tens%create_terms(B)

        ! Navier-Stokes x-axis terms
        NS_eqn%x%aT = -nu*aT_lapl + NS_adv_term%aT
        NS_eqn%x%aN = -nu*aN_lapl + NS_adv_term%aN
        NS_eqn%x%aE = -nu*aE_lapl + NS_adv_term%aE
        NS_eqn%x%aP = -nu*aP_lapl + NS_adv_term%aP + aP_time
        NS_eqn%x%aW = -nu*aW_lapl + NS_adv_term%aW
        NS_eqn%x%aS = -nu*aS_lapl + NS_adv_term%aS
        NS_eqn%x%aB = -nu*aB_lapl + NS_adv_term%aB
        NS_eqn%x%B  =  aP_time*V%x(nodes_P) - (1.0/rho)*gradx(p%field) + (1.0/rho)*crossx(Je,B) - beta*(T%field(nodes_P)-T0)*g_value*e_g(1)

        ! Navier-Stokes y-axis terms
        NS_eqn%y%aT = -nu*aT_lapl + NS_adv_term%aT
        NS_eqn%y%aN = -nu*aN_lapl + NS_adv_term%aN
        NS_eqn%y%aE = -nu*aE_lapl + NS_adv_term%aE
        NS_eqn%y%aP = -nu*aP_lapl + NS_adv_term%aP + aP_time
        NS_eqn%y%aW = -nu*aW_lapl + NS_adv_term%aW
        NS_eqn%y%aS = -nu*aS_lapl + NS_adv_term%aS
        NS_eqn%y%aB = -nu*aB_lapl + NS_adv_term%aB
        NS_eqn%y%B  =  aP_time*V%y(nodes_P) - (1.0/rho)*grady(p%field) + (1.0/rho)*crossy(Je,B) - beta*(T%field(nodes_P)-T0)*g_value*e_g(2)

        ! Navier-Stokes z-axis terms
        if (flow_2D .eqv. .FALSE.)  then
            NS_eqn%z%aT = -nu*aT_lapl + NS_adv_term%aT
            NS_eqn%z%aN = -nu*aN_lapl + NS_adv_term%aN
            NS_eqn%z%aE = -nu*aE_lapl + NS_adv_term%aE
            NS_eqn%z%aP = -nu*aP_lapl + NS_adv_term%aP + aP_time
            NS_eqn%z%aW = -nu*aW_lapl + NS_adv_term%aW
            NS_eqn%z%aS = -nu*aS_lapl + NS_adv_term%aS
            NS_eqn%z%aB = -nu*aB_lapl + NS_adv_term%aB
            NS_eqn%z%B  =  aP_time*V%z(nodes_P) - (1.0/rho)*gradz(p%field) + (1.0/rho)*crossz(Je,B) - beta*(T%field(nodes_P)-T0)*g_value*e_g(3)
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
        
        !------------------------------------------------------------------------------------------
        ! Induction Equation
        !------------------------------------------------------------------------------------------
        ! Advection terms
        call Ind_adv_term%create_terms(V)
        call Ind_Bi_adv%create_terms(Bi)
        call Ind_Bo_adv%create_terms(Bo)

        ! Induction equation x-axis terms
        Ind_eqn%x%aT = -(1.0/(mu*sigma))*aT_lapl + Ind_adv_term%aT
        Ind_eqn%x%aN = -(1.0/(mu*sigma))*aN_lapl + Ind_adv_term%aN
        Ind_eqn%x%aE = -(1.0/(mu*sigma))*aE_lapl + Ind_adv_term%aE
        Ind_eqn%x%aP = -(1.0/(mu*sigma))*aP_lapl + Ind_adv_term%aP + aP_time
        Ind_eqn%x%aW = -(1.0/(mu*sigma))*aW_lapl + Ind_adv_term%aW
        Ind_eqn%x%aS = -(1.0/(mu*sigma))*aS_lapl + Ind_adv_term%aS
        Ind_eqn%x%aB = -(1.0/(mu*sigma))*aB_lapl + Ind_adv_term%aB
        Ind_eqn%x%B  = aP_time*Bi%x(nodes_P) + flux(Ind_Bi_adv,V%x) + flux(Ind_Bo_adv,V%x) - flux(Ind_adv_term,Bo%x)

        ! Induction equation y-axis terms
        Ind_eqn%y%aT = -(1.0/(mu*sigma))*aT_lapl + Ind_adv_term%aT
        Ind_eqn%y%aN = -(1.0/(mu*sigma))*aN_lapl + Ind_adv_term%aN
        Ind_eqn%y%aE = -(1.0/(mu*sigma))*aE_lapl + Ind_adv_term%aE
        Ind_eqn%y%aP = -(1.0/(mu*sigma))*aP_lapl + Ind_adv_term%aP + aP_time
        Ind_eqn%y%aW = -(1.0/(mu*sigma))*aW_lapl + Ind_adv_term%aW
        Ind_eqn%y%aS = -(1.0/(mu*sigma))*aS_lapl + Ind_adv_term%aS
        Ind_eqn%y%aB = -(1.0/(mu*sigma))*aB_lapl + Ind_adv_term%aB
        Ind_eqn%y%B  = aP_time*Bi%y(nodes_P) + flux(Ind_Bi_adv,V%y) + flux(Ind_Bo_adv,V%y) - flux(Ind_adv_term,Bo%y)

        ! Induction equation z-axis terms
        if (flow_2D .eqv. .FALSE.)  then
            Ind_eqn%z%aT = -(1.0/(mu*sigma))*aT_lapl + Ind_adv_term%aT
            Ind_eqn%z%aN = -(1.0/(mu*sigma))*aN_lapl + Ind_adv_term%aN
            Ind_eqn%z%aE = -(1.0/(mu*sigma))*aE_lapl + Ind_adv_term%aE
            Ind_eqn%z%aP = -(1.0/(mu*sigma))*aP_lapl + Ind_adv_term%aP + aP_time
            Ind_eqn%z%aW = -(1.0/(mu*sigma))*aW_lapl + Ind_adv_term%aW
            Ind_eqn%z%aS = -(1.0/(mu*sigma))*aS_lapl + Ind_adv_term%aS
            Ind_eqn%z%aB = -(1.0/(mu*sigma))*aB_lapl + Ind_adv_term%aB
            Ind_eqn%z%B  = aP_time*Bi%z(nodes_P) + flux(Ind_Bi_adv,V%z) + flux(Ind_Bo_adv,V%z) - flux(Ind_adv_term,Bo%z)
        end if

        !------------------------------------------------------------------------------------------
        ! Magnetic Field Correction Equation
        !------------------------------------------------------------------------------------------
        ! Poisson x-axis terms
        dB_Eqn%x%aT = aT_lapl
        dB_Eqn%x%aN = aN_lapl
        dB_Eqn%x%aE = aE_lapl
        dB_Eqn%x%aP = aP_lapl
        dB_Eqn%x%aW = aW_lapl
        dB_Eqn%x%aS = aS_lapl
        dB_Eqn%x%aB = aB_lapl
        dB_Eqn%x%B  = gradx(G%field)

        ! Poisson y-axis terms
        dB_Eqn%y%aT = aT_lapl
        dB_Eqn%y%aN = aN_lapl
        dB_Eqn%y%aE = aE_lapl
        dB_Eqn%y%aP = aP_lapl
        dB_Eqn%y%aW = aW_lapl
        dB_Eqn%y%aS = aS_lapl
        dB_Eqn%y%aB = aB_lapl
        dB_Eqn%y%B  = grady(G%field)

        ! Poisson y-axis terms
        if (flow_2D .eqv. .FALSE.)  then
            dB_Eqn%z%aT = aT_lapl
            dB_Eqn%z%aN = aN_lapl
            dB_Eqn%z%aE = aE_lapl
            dB_Eqn%z%aP = aP_lapl
            dB_Eqn%z%aW = aW_lapl
            dB_Eqn%z%aS = aS_lapl
            dB_Eqn%z%aB = aB_lapl
            dB_Eqn%z%B  = gradz(G%field)
        end if


        !------------------------------------------------------------------------------------------
        ! Heat Transport Equation
        !------------------------------------------------------------------------------------------
        call Heat_conv_term%create_terms(V)

        Heat_eqn%aT = -alpha*aT_lapl + Heat_conv_term%aT
        Heat_eqn%aN = -alpha*aN_lapl + Heat_conv_term%aN
        Heat_eqn%aE = -alpha*aE_lapl + Heat_conv_term%aE
        Heat_eqn%aP = -alpha*aP_lapl + Heat_conv_term%aP + aP_time
        Heat_eqn%aW = -alpha*aW_lapl + Heat_conv_term%aW
        Heat_eqn%aS = -alpha*aS_lapl + Heat_conv_term%aS
        Heat_eqn%aB = -alpha*aB_lapl + Heat_conv_term%aB
        Heat_eqn%B  = aP_time*T%field(nodes_P)