#include <define.h>
#ifdef CROP
module MOD_Irrigation

!   DESCRIPTION:
!       This module has all irrigation related subroutines for irrigated crop at either IGBP/USGS or PFT Land type classification and even in the C and N cycle.
    use MOD_Precision
    USE MOD_TimeManager
    USE MOD_Namelist, only: DEF_simulation_time
    ! ,DEF_IRRIGATION_METHOD
    use MOD_Const_Physical, only: tfrz
    use MOD_Const_PFT, only: irrig_crop
    use MOD_Vars_Global, only: irrig_start_time, irrig_max_depth, irrig_threshold_fraction, irrig_min_cphase, irrig_max_cphase, irrig_time_per_day    
    use MOD_Qsadv, only: qsadv
    use MOD_Vars_TimeInvariants, only: &
#ifdef vanGenuchten_Mualem_SOIL_MODEL
        theta_r, alpha_vgm, n_vgm, L_vgm, sc_vgm, fc_vgm, &
#endif
        porsl, psi0, bsw
    use MOD_Vars_TimeVariables, only : tref, t_soisno, wliq_soisno, zwt, &
        irrig_rate, sum_irrig, sum_irrig_count, n_irrig_steps_left, &
        tairday, usday, vsday, pairday, rnetday, fgrndday, potential_evapotranspiration, &
        waterstorage_supply, uncongwirrig_supply, congwirrig_supply, &
        waterstorage, deficit_irrig, actual_irrig
    use MOD_Vars_PFTimeInvariants, only: pftclass
    use MOD_Vars_PFTimeVariables, only: irrig_method_p
    use MOD_BGC_Vars_PFTimeVariables, only: cphase_p
    use MOD_Vars_1DForcing, only: forc_t, forc_frl, forc_psrf, forc_us, forc_vs
    use MOD_Vars_1DFluxes, only: sabg, sabvsun, sabvsha, olrg, fgrnd
    use MOD_Hydro_SoilFunction, only: soil_vliq_from_psi

    implicit none

    public :: CalIrrigationNeeded
    public :: CalIrrigationApplicationFluxes

    !   local variable
    integer :: irrig_method_drip = 1
    integer :: irrig_method_sprinkler = 2
    integer :: irrig_method_flood = 3
    integer :: irrig_method_paddy = 4

contains
            
    subroutine CalIrrigationNeeded(i,ps,pe,idate,nl_soil,nbedrock,z_soi,dz_soi,zi_soi,deltim,dlon,npcropmin)

        !   DESCRIPTION:
        !   This subroutine is used to calculate how much irrigation needed in each irrigated crop patch
        integer , intent(in) :: i
        integer , intent(in) :: ps, pe
        integer , intent(in) :: idate(3)
        integer , intent(in) :: nl_soil
        integer , intent(in) :: nbedrock
        real(r8), intent(in) :: z_soi(1:nl_soil)
        real(r8), intent(in) :: dz_soi(1:nl_soil)
        real(r8), intent(in) :: zi_soi(1:nl_soil)
        real(r8), intent(in) :: deltim
        real(r8), intent(in) :: dlon
        integer , intent(in) :: npcropmin

        ! local 
        integer :: m
        integer :: irrig_nsteps_per_day
        logical :: check_for_irrig 

        !   initialize irrigation
        deficit_irrig(i) = 0._r8
        actual_irrig(i) = 0._r8

        ! !   calculate last day potential evapotranspiration 
        ! call CalPotentialEvapotranspiration(i,idate,dlon,deltim)

        !   calculate whether irrigation needed
        call PointNeedsCheckForIrrig(i,ps,pe,idate,deltim,dlon,npcropmin,check_for_irrig)

        !   calculate irrigation needed
        if (check_for_irrig) then
            call CalIrrigationPotentialNeeded(i,ps,pe,nl_soil,nbedrock,z_soi,dz_soi)
            call CalIrrigationLimitedSupply(i,nl_soil,nbedrock,zi_soi)
        end if

        !   calculate irrigation rate kg/m2->mm/s
        if ((check_for_irrig) .and. (actual_irrig(i) > 0)) then
            irrig_nsteps_per_day = nint(irrig_time_per_day/deltim)
            irrig_rate(i) = actual_irrig(i)/deltim/irrig_nsteps_per_day
            n_irrig_steps_left(i) = irrig_nsteps_per_day
            sum_irrig(i) = sum_irrig(i) + actual_irrig(i)
            sum_irrig_count(i) = sum_irrig_count(i) + 1._r8
        end if
        
        ! !   zero irrigation at the end of growing season 
        ! do m = ps, pe
        !     if (cphase_p(m) >= 4._r8) then
        !         sum_irrig(i) = 0._r8
        !         sum_irrig_count(i) = 0._r8
        !     end if
        ! end do
    end subroutine CalIrrigationNeeded


    subroutine CalIrrigationPotentialNeeded(i,ps,pe,nl_soil,nbedrock,z_soi,dz_soi)

        !   DESCRIPTION:
        !   This subroutine is used to calculate how much irrigation needed in each irrigated crop patch without water supply restriction
        integer , intent(in) :: i
        integer , intent(in) :: ps, pe
        integer , intent(in) :: nbedrock
        integer , intent(in) :: nl_soil
        real(r8), intent(in) :: z_soi(1:nl_soil)
        real(r8), intent(in) :: dz_soi(1:nl_soil)

        !   local variables
        integer  :: j
        integer  :: m
        logical  :: reached_max_depth
        real(r8) :: h2osoi_liq_tot
        real(r8) :: h2osoi_liq_target_tot
        real(r8) :: h2osoi_liq_wilting_point_tot
        real(r8) :: h2osoi_liq_saturation_capacity_tot
        real(r8) :: h2osoi_liq_wilting_point(1:nl_soil)
        real(r8) :: h2osoi_liq_field_capacity(1:nl_soil)
        real(r8) :: h2osoi_liq_saturation_capacity(1:nl_soil)
        real(r8) :: h2osoi_liq_at_threshold

        real(r8) :: smpswc = -1.5e5
        real(r8) :: smpsfc = -3.3e3  

        !   initialize local variables
        reached_max_depth = .false.
        h2osoi_liq_tot = 0._r8
        h2osoi_liq_target_tot = 0._r8
        h2osoi_liq_wilting_point_tot = 0._r8
        h2osoi_liq_saturation_capacity_tot = 0._r8

        ! !   single site initialization
        ! do m = ps, pe
        !     irrig_method_p(m) = DEF_IRRIGATION_METHOD
        ! enddo

!   calculate wilting point and field capacity
        do j = 1, nl_soil
            if (t_soisno(j,i) > tfrz .and. porsl(j,i) >= 1.e-6) then           
#ifdef Campbell_SOIL_MODEL
                h2osoi_liq_wilting_point(j) = 1000.*dz_soi(j)*porsl(j,i)*((smpswc/psi0(j,i))**(-1/bsw(j,i)))
                h2osoi_liq_field_capacity(j) = 1000.*dz_soi(j)*porsl(j,i)*((smpsfc/psi0(j,i))**(-1/bsw(j,i)))
                h2osoi_liq_saturation_capacity(j) = 1000.*dz_soi(j)*porsl(j,i)
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
                h2osoi_liq_wilting_point(j) = soil_vliq_from_psi(smpswc, porsl(j,i), theta_r(j,i), psi0(j,i), 5, &
                  (/alpha_vgm(j,i), n_vgm(j,i), L_vgm(j,i), sc_vgm(j,i), fc_vgm(j,i)/))
                h2osoi_liq_wilting_point(j) = 1000.*dz_soi(j)*h2osoi_liq_wilting_point(j)
                h2osoi_liq_field_capacity(j) = soil_vliq_from_psi(smpsfc, porsl(j,i), theta_r(j,i), psi0(j,i), 5, &
                  (/alpha_vgm(j,i), n_vgm(j,i), L_vgm(j,i), sc_vgm(j,i), fc_vgm(j,i)/))
                h2osoi_liq_field_capacity(j) = 1000.*dz_soi(j)*h2osoi_liq_field_capacity(j)
                h2osoi_liq_saturation_capacity(j) = 1000.*dz_soi(j)*porsl(j,i)      
#endif
            end if
        end do 

        !   calculate total irrigation needed in all soil layers
        do m = ps, pe
            do j = 1, nl_soil
                if (.not. reached_max_depth) then
                    if (z_soi(j) > irrig_max_depth) then
                        reached_max_depth = .true.
                    else if (j > nbedrock) then
                        reached_max_depth = .true.
                    else if (t_soisno(j,i) <= tfrz) then
                        reached_max_depth = .true.
                    else 
                        h2osoi_liq_tot = h2osoi_liq_tot + wliq_soisno(j,i)
                        h2osoi_liq_wilting_point_tot = h2osoi_liq_wilting_point_tot + h2osoi_liq_wilting_point(j)
                        if (irrig_method_p(m) == irrig_method_drip .or. irrig_method_p(m) == irrig_method_sprinkler) then
                            h2osoi_liq_target_tot = h2osoi_liq_target_tot + h2osoi_liq_field_capacity(j)
                        !   irrigation threshold at field capacity, but irrigation amount at saturation capacity
                        else if (irrig_method_p(m) == irrig_method_flood) then
                            h2osoi_liq_target_tot = h2osoi_liq_target_tot + h2osoi_liq_field_capacity(j)
                            h2osoi_liq_saturation_capacity_tot = h2osoi_liq_saturation_capacity_tot + h2osoi_liq_saturation_capacity(j)
                        else if (irrig_method_p(m) == irrig_method_paddy) then
                            h2osoi_liq_target_tot = h2osoi_liq_target_tot + h2osoi_liq_saturation_capacity(j)
                        else
                            h2osoi_liq_target_tot = h2osoi_liq_target_tot + h2osoi_liq_field_capacity(j)
                        end if
                    end if 
                end if
            end do 
        end do

        !   calculate irrigation threshold
        deficit_irrig(i) = 0._r8
        h2osoi_liq_at_threshold = h2osoi_liq_wilting_point_tot + irrig_threshold_fraction * (h2osoi_liq_target_tot - h2osoi_liq_wilting_point_tot)

        !   calculate total irrigation
        do m = ps, pe
            if (h2osoi_liq_tot < h2osoi_liq_at_threshold) then
                if (irrig_method_p(m) == irrig_method_sprinkler) then 
                    deficit_irrig(i) = h2osoi_liq_target_tot - h2osoi_liq_tot
                    ! deficit_irrig(i) = h2osoi_liq_target_tot - h2osoi_liq_tot + potential_evapotranspiration(i)
                else if (irrig_method_p(m) == irrig_method_flood) then
                    deficit_irrig(i) = h2osoi_liq_saturation_capacity_tot - h2osoi_liq_tot
                else
                    deficit_irrig(i) = h2osoi_liq_at_threshold - h2osoi_liq_tot
                end if
            else
                deficit_irrig(i) = 0
            end if
        end do

    end subroutine CalIrrigationPotentialNeeded

    subroutine CalIrrigationApplicationFluxes(i,ps,pe,deltim,qflx_irrig_drip,qflx_irrig_sprinkler,qflx_irrig_flood,qflx_irrig_paddy)
        !   DESCRIPTION:
        !   This subroutine is used to calculate irrigation application fluxes for each irrigated crop patch
        integer , intent(in) :: i
        integer , intent(in) :: ps, pe
        real(r8), intent(in) :: deltim
        real(r8), intent(out):: qflx_irrig_drip,qflx_irrig_sprinkler,qflx_irrig_flood,qflx_irrig_paddy

        integer :: m 

        qflx_irrig_drip = 0._r8
        qflx_irrig_sprinkler = 0._r8
        qflx_irrig_flood = 0._r8
        qflx_irrig_paddy = 0._r8

        ! !   single site initialization
        ! do m = ps, pe
        !     irrig_method_p(m) = DEF_IRRIGATION_METHOD
        ! enddo

        !   add irrigation fluxes to precipitation or land surface
        do m = ps, pe
            if (n_irrig_steps_left(i) > 0) then
                if (irrig_method_p(m) == irrig_method_drip) then
                    qflx_irrig_drip = irrig_rate(i)
                else if (irrig_method_p(m) == irrig_method_sprinkler) then 
                    qflx_irrig_sprinkler = irrig_rate(i)
                else if (irrig_method_p(m) == irrig_method_flood) then
                    qflx_irrig_flood = irrig_rate(i)
                else if (irrig_method_p(m) == irrig_method_paddy) then
                    qflx_irrig_paddy = irrig_rate(i)
                else
                    qflx_irrig_drip = irrig_rate(i)
                end if
                n_irrig_steps_left(i) = n_irrig_steps_left(i) -1
            else 
                irrig_rate(i) = 0._r8
            end if
        end do
    end subroutine CalIrrigationApplicationFluxes
    
    subroutine PointNeedsCheckForIrrig(i,ps,pe,idate,deltim,dlon,npcropmin,check_for_irrig)
        !   DESCRIPTION:
        !   This subroutine is used to calculate whether irrigation needed in each patch
        integer , intent(in) :: i
        integer , intent(in) :: ps, pe
        integer , intent(in) :: idate(3)
        real(r8), intent(in) :: deltim
        real(r8), intent(in) :: dlon
        integer , intent(in) :: npcropmin
        logical , intent(out):: check_for_irrig

        !   local variable
        integer :: m, ivt
        real(r8):: ldate(3)
        real(r8):: seconds_since_irrig_start_time

        !   adjust flood irrigation in rice to paddy irrigaiton
        do m = ps, pe 
            ivt = pftclass(m)
            if ((ivt == 62) .and. (irrig_method_p(m) == irrig_method_flood))then
                irrig_method_p(m) = irrig_method_paddy
            endif
        enddo

        do m = ps, pe 
            ivt = pftclass(m)
            if ((ivt >= npcropmin) .and. (irrig_crop(ivt)) .and. &
                (cphase_p(m) >= irrig_min_cphase) .and. (cphase_p(m)<irrig_max_cphase)) then
                if (DEF_simulation_time%greenwich) then
                    call gmt2local(idate, dlon, ldate)
                    seconds_since_irrig_start_time = ldate(3) - irrig_start_time + deltim
                else
                    seconds_since_irrig_start_time = idate(3) - irrig_start_time + deltim
                end if
                if ((seconds_since_irrig_start_time >= 0._r8) .and. (seconds_since_irrig_start_time < deltim)) then
                    check_for_irrig = .true.
                else
                    check_for_irrig = .false.
                end if
            else
                check_for_irrig = .false.
            end if
        end do

    end subroutine PointNeedsCheckForIrrig

    ! subroutine CalPotentialEvapotranspiration(i,idate,dlon,deltim)
    !     !   DESCRIPTION:
    !     !   This subroutine is used to calculate daily potential evapotranspiration
    !     integer , intent(in) :: i
    !     integer , intent(in) :: idate(3)
    !     real(r8), intent(in) :: dlon
    !     real(r8), intent(in) :: deltim

    !     !   local variable
    !     real(r8):: ldate(3)
    !     real(r8):: seconds_since_irrig_start_time
    !     real(r8) :: es,esdT,qs,qsdT     ! saturation vapour pressure
    !     real(r8) :: evsat               ! vapour pressure
    !     real(r8) :: ur                  ! wind speed
    !     real(r8) :: delta               ! slope of saturation vapour pressure curve 
    !     real(r8) :: gamma               ! Psychrometric constant

    !     if (DEF_simulation_time%greenwich) then
    !         call gmt2local(idate, dlon, ldate)
    !         seconds_since_irrig_start_time = ldate(3) - irrig_start_time + deltim
    !     else
    !         seconds_since_irrig_start_time = idate(3) - irrig_start_time + deltim
    !     end if

    !     if (((seconds_since_irrig_start_time-deltim) >= 0) .and. ((seconds_since_irrig_start_time-deltim) < deltim)) then
    !         tairday(i) = (forc_t(i)-tfrz)*deltim/86400
    !         usday(i) = forc_us(i)*deltim/86400
    !         vsday(i) = forc_vs(i)*deltim/86400
    !         pairday(i) = forc_psrf(i)*deltim/86400/1000
    !         rnetday(i) = (sabg(i)+sabvsun(i)+sabvsha(i)-olrg(i)+forc_frl(i))*deltim/1000000
    !         fgrndday(i) = fgrnd(i)*deltim/1000000
    !     else
    !         tairday(i) = tairday(i) + (forc_t(i)-tfrz)*deltim/86400
    !         usday(i) = usday(i) + forc_us(i)*deltim/86400
    !         vsday(i) = vsday(i) + forc_vs(i)*deltim/86400
    !         pairday(i) = pairday(i) + forc_psrf(i)*deltim/86400/1000
    !         rnetday(i) = rnetday(i) + (sabg(i)+sabvsun(i)+sabvsha(i)-olrg(i)+forc_frl(i))*deltim/1000000
    !         fgrndday(i) = fgrndday(i) + fgrnd(i)*deltim/1000000
    !     endif

    !     if ((seconds_since_irrig_start_time >= 0) .and. (seconds_since_irrig_start_time < deltim)) then
    !         call qsadv(tairday(i),pairday(i),es,esdT,qs,qsdT)
    !         if (tairday(i) > 0)then
    !             evsat = 0.611*EXP(17.27*tairday(i)/(tairday(i)+237.3))
    !         else
    !             evsat = 0.611*EXP(21.87*tairday(i)/(tairday(i)+265.5))
    !         endif
    !         ur = max(0.1,sqrt(usday(i)*usday(i)+vsday(i)*vsday(i)))
    !         delta = 4098*evsat/((tairday(i)+237.3)*(tairday(i)+237.3))
    !         gamma = 0.665*0.001*pairday(i)
    !         potential_evapotranspiration(i) = (0.408*delta*(rnetday(i)-fgrndday(i))+gamma*(900/(tairday(i)+273))*ur* &
    !             (evsat-es))/(delta+(gamma*(1+0.34*ur)))
    !     end if
    ! end subroutine CalPotentialEvapotranspiration

    subroutine CalIrrigationLimitedSupply(i,nl_soil,nbedrock,zi_soi)
        !   DESCRIPTION:
        !   This subroutine is used to calculate how much irrigation supplied in each irrigated crop patch with water supply restriction
        integer, intent(in) :: i
        integer, intent(in) :: nl_soil
        integer, intent(in) :: nbedrock
        real(r8),intent(in) :: zi_soi(1:nl_soil)
        
        logical :: limitedirrig

        waterstorage_supply(i) = 0._r8
        uncongwirrig_supply(i) = 0._r8
        congwirrig_supply(i) = 0._r8
        limitedirrig = .true.

        !   irrigation withdraw from water storage pool to adjust the unmatched time for different water supply systems
        if (deficit_irrig(i) > 0._r8) then
            waterstorage_supply(i) = min(waterstorage(i),deficit_irrig(i))
            waterstorage(i) = waterstorage(i) - waterstorage_supply(i)
            actual_irrig(i) = actual_irrig(i) + waterstorage_supply(i)
        endif
#ifdef CaMa_Flood
        ! call CalWithdrawReservoirWater() deficit noupdate
        ! call CalWithdrawRiverWater() deficit actual
#endif
        !   irrigation withdraw from unconfined ground water
        if (deficit_irrig(i) > actual_irrig(i)) then
            call CalWithdrawUndergroundWater(i,nl_soil,nbedrock,zi_soi)
            actual_irrig(i) = actual_irrig(i) + uncongwirrig_supply(i)
        end if
        !   irrigation withdraw from confined ground water
        if ((deficit_irrig(i) > actual_irrig(i)) .and. (.not. limitedirrig)) then
            congwirrig_supply(i) = deficit_irrig(i)
            actual_irrig(i) = actual_irrig(i) + congwirrig_supply(i)
        endif
    end subroutine CalIrrigationLimitedSupply


!   地下水抽取部分，需要在这个部分重新call water过程？？？ 还是保留相应的信息，然后在后面计算？？？
!   常规的地下水方案->call subsurfacerunoff；VSF方案->call soilwater_aquifer_exchange
    subroutine CalWithdrawUndergroundWater(i,nl_soil,nbedrock,zi_soi)
        !   DESCRIPTION:
        !   This subroutine is used to calculate how much irrigation supplied in each irrigated crop patch from unconfined ground water
        integer, intent(in) :: i
        integer, intent(in) :: nl_soil
        integer, intent(in) :: nbedrock
        real(r8),intent(in) :: zi_soi(1:nl_soil)

        !   local variable
        integer :: j, jwt
        real(r8):: s_y
        real(r8):: uncongwirrig_demand
        real(r8):: uncongwirrig_supply_layer(1:nl_soil)
            
        ! if (.not. DEF_USE_VARIABLY_SATURATED_FLOW) then
            uncongwirrig_demand = -(deficit_irrig(i) - actual_irrig(i))
            uncongwirrig_supply(i) = 0._r8
            do j = 1, nl_soil 
                uncongwirrig_supply_layer(j) = 0._r8
            end do 
            
            !   The layer index of the first unsaturated layer
            jwt = nl_soil
            !   allow jwt to equal zero when zwt is in top layer
            do j = 1, nl_soil
                if (zwt(i) <= zi_soi(j)) then
                    jwt = j-1
                    exit
                end if
            enddo

            do j = jwt+1, nl_soil
                if (j > nbedrock) then
                    exit
                else
                    ! use analytical expression for specific yield
                    s_y = porsl(j,i) * (1.-(1.-1.e3*zwt(i)/psi0(j,i))**(-1./bsw(j,i)))
                    s_y = max(s_y,0.02)
                    uncongwirrig_supply_layer(j) = max(uncongwirrig_demand,-(s_y*(zi_soi(j) - zwt(i))*1.e3))
                    uncongwirrig_supply_layer(j) = min(uncongwirrig_supply_layer(j),0.)
                    uncongwirrig_demand = uncongwirrig_demand - uncongwirrig_supply_layer(j)

                    if (uncongwirrig_demand >= 0.) then
                        zwt = max(0.,zwt - uncongwirrig_supply_layer(j)/s_y/1000.)
                        exit
                    else
                        zwt(i) = zi_soi(j)
                    endif
                end if
            enddo

            do j = 1, nl_soil
                uncongwirrig_supply(i) = uncongwirrig_supply(i) - uncongwirrig_supply_layer(j)
                wliq_soisno(j,i) = wliq_soisno(j,i) + uncongwirrig_supply_layer(j)
            enddo

        ! else
            ! uncongwirrig_demand = deficit_irrig(i)
            ! do j = 1, nl_soil 
            !     uncongwirrig_supply_layer(j) = 0._r8
            ! end do 
            ! sp_zi(0:nl_soil) = zi_soi(0:nl_soil) * 1000.0   ! from meter to mm

            ! izwt = findloc(zwt >= sp_zi, .true., dim=1, back=.true.)
            ! is_permeable(j) = eff_porosity(j) > max(wimp, theta_r(j))

            ! DO WHILE (uncongwirrig_demand > 0.)
            !     IF (izwt <= nlev) THEN
            !         IF (is_permeable(izwt)) THEN
            !             call get_zwt_from_wa ( &
            !                 porsl(izwt), vl_r(izwt), psi_s(izwt), hksat(izwt), &
            !                 nprm, prms(:,izwt), tol_v, tol_z, -uncongwirrig_demand, zwt, zwtp)
            !             IF (zwtp < sp_zi(izwt)) THEN
            !                 ss_vliq(izwt) = (ss_vliq(izwt)*(zwtmm-sp_zi(izwt-1))  &
            !                     + porsl(izwt)*(zwtp-zwtmm) - uncongwirrig_demand) / (zwtp - sp_zi(izwt-1))
            !                 uncongwirrig_supply_layer(j) = uncongwirrig_demand
            !                 uncongwirrig_demand = 0.
            !                 zwtmm = zwtp
            !             ELSE
            !                 psi  = psi_s(izwt) - (zwtp - 0.5*(sp_zi(izwt) + zwtmm))
            !                 vliq = soil_vliq_from_psi (psi, &
            !                     porsl(izwt), vl_r(izwt), psi_s(izwt), nprm, prms(:,izwt))
            !                 IF (uncongwirrig_demand > (porsl(izwt)-vliq) * (sp_zi(izwt)-zwtmm)) THEN
            !                     ss_vliq(izwt) = (ss_vliq(izwt)*(zwtmm-sp_zi(izwt-1))  &
            !                     + vliq * (sp_zi(izwt)-zwtmm)) / sp_dz(izwt)
            !                     uncongwirrig_supply_layer(j) = (porsl(izwt)-vliq) * (sp_zi(izwt)-zwtmm)
            !                     uncongwirrig_demand = uncongwirrig_demand - (porsl(izwt)-vliq) * (sp_zi(izwt)-zwtmm)
            !                 ELSE
            !                     ss_vliq(izwt) = (ss_vliq(izwt)*(zwtmm-sp_zi(izwt-1))  &
            !                     + porsl(izwt)*(sp_zi(izwt)-zwtmm) - uncongwirrig_demand) / sp_dz(izwt)
            !                     uncongwirrig_supply_layer(j) = uncongwirrig_demand
            !                     uncongwirrig_demand = 0.
            !                 ENDIF
            !                 zwtmm  = sp_zi(izwt)
            !                 izwt = izwt + 1
            !             ENDIF
            !         ELSE
            !             zwtmm  = sp_zi(izwt)
            !             izwt = izwt + 1
            !         ENDIF
            !     ENDIF
            ! ENDDO

            ! ! update the mass of liquid water
            ! DO j = nl_soil, 1, -1
            !     IF (is_permeable(j)) THEN
            !         IF (zwtmm < sp_zi(j)) THEN
            !             IF (zwtmm >= sp_zi(j-1)) THEN
            !                 wliq_soisno(j,i)  = denh2o * ((eff_porosity(j)*(sp_zi(j)-zwtmm))  &
            !                     + vol_liq(j) * (zwtmm - sp_zi(j-1)))/1000.0
            !             ELSE
            !                 wliq_soisno(j,i)  = denh2o * (eff_porosity(j)*(sp_zi(j)-sp_zi(j-1)))/1000.0
            !             ENDIF
            !         ELSE
            !             wliq_soisno(j,i) = denh2o * (vol_liq(j)*(sp_zi(j)-sp_zi(j-1)))/1000.0
            !         ENDIF
            !     ENDIF
            ! ENDDO
            ! zwt = zwtmm/1000.0
        ! end if
        end subroutine CalWithdrawUndergroundWater

end module MOD_Irrigation
#endif
