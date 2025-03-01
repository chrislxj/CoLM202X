#include <define.h>
#ifdef CROP
MODULE MOD_Irrigation

!  DESCRIPTION:
!      This MODULE has all irrigation related subroutines for irrigated crop at either IGBP/USGS or PFT Land type classification and even in the C and N cycle.
   use MOD_Precision
   use MOD_TimeManager
   use MOD_Namelist, only: DEF_simulation_time, DEF_IRRIGATION_ALLOCATION, DEF_IRRIGATION_METHOD, DEF_USE_VariablySaturatedFlow
   use MOD_Const_Physical, only: tfrz, denice, denh2o
   use MOD_Const_PFT, only: irrig_crop
   use MOD_LandPFT, only : patch_pft_s, patch_pft_e
   use MOD_Vars_Global, only: irrig_start_time, irrig_max_depth, irrig_threshold_fraction, irrig_supply_fraction, irrig_min_cphase, irrig_max_cphase, irrig_time_per_day, &
         irrig_method_drip, irrig_method_sprinkler, irrig_method_flood, irrig_method_paddy
   use MOD_Qsadv, only: qsadv
   use MOD_Vars_TimeInvariants, only: pondmx, &
#ifdef vanGenuchten_Mualem_SOIL_MODEL
       theta_r, alpha_vgm, n_vgm, L_vgm, fc_vgm, sc_vgm,&
#endif
       porsl, psi0, bsw
   use MOD_Vars_TimeVariables, only : tref, t_soisno, wliq_soisno, wice_soisno, zwt, wa, &
         irrig_rate, sum_irrig, sum_deficit_irrig, sum_irrig_count, n_irrig_steps_left, &
         tairday, usday, vsday, pairday, rnetday, fgrndday, potential_evapotranspiration,&
         waterstorage_supply, groundwater_demand, groundwater_supply, reservoirriver_demand, reservoirriver_supply, &
         groundwater_supply_rn, groundwater_supply_nrn, reservoir_supply, river_supply, runoff_supply, &
         waterstorage, deficit_irrig, actual_irrig, irrig_gw_alloc, irrig_sw_alloc, zwt_stand, irrig_demand_days, irrig_supply_satisfy_days, irrig_supply_unsatisfy_days
   use MOD_Vars_PFTimeInvariants, only: pftclass
   use MOD_Vars_PFTimeVariables, only: irrig_method_p
   use MOD_BGC_Vars_PFTimeVariables, only: cphase_p
   use MOD_Vars_1DForcing, only: forc_t, forc_frl, forc_psrf, forc_us, forc_vs
   use MOD_Vars_1DFluxes, only: sabg, sabvsun, sabvsha, olrg, fgrnd
   use MOD_Hydro_SoilFunction, only: soil_vliq_from_psi
#ifdef CaMa_Flood
    use MOD_CaMa_Vars, only: withdrawal_cama, withdrawal_dam_cama, withdrawal_riv_cama, withdrawal_rof_cama
#endif
   ! use MOD_SPMD_Task


   implicit none

   public :: CalIrrigationNeeded
   public :: CalIrrigationApplicationFluxes

contains

   subroutine CalIrrigationNeeded(i,idate,nl_soil,nbedrock,zi_soi,dz_soi,deltim,dlon,npcropmin)

      !   DESCRIPTION:
      !   This subroutine is used to calculate how much irrigation needed in each irrigated crop patch
      integer , intent(in) :: i
      integer , intent(in) :: idate(3)
      integer , intent(in) :: nl_soil
      integer , intent(in) :: nbedrock
      real(r8), intent(in) :: zi_soi(1:nl_soil)
      real(r8), intent(in) :: dz_soi(1:nl_soil)
      real(r8), intent(in) :: deltim
      real(r8), intent(in) :: dlon
      integer , intent(in) :: npcropmin

      ! local
      integer :: ps, pe, m
      integer :: irrig_nsteps_per_day
      logical :: check_for_irrig


      ps = patch_pft_s(i)
      pe = patch_pft_e(i)

! #ifdef USEMPI
!             CALL mpi_barrier (p_comm_glb, p_err)
! #endif
      ! if(withdrawal_cama(i).gt.0)then
      !    WRITE(*, '(A, I4, F16.4)'),"LHB debug line68 withdraw error : patch, demand -----> ", i, dirrig_tmp(i)
      !    WRITE(*, '(A, I4, F16.4)'),"LHB debug line68 withdraw error : patch, supply -----> ", i, withdrawal_cama(i)
      ! endif

      !  initialize irrigation
      deficit_irrig(i) = 0._r8
      actual_irrig(i) = 0._r8
      waterstorage_supply(i) = 0._r8
      groundwater_demand(i) = 0._r8
      groundwater_supply(i) = 0._r8
      reservoirriver_demand(i) = 0._r8
      reservoirriver_supply(i) = 0._r8
      reservoir_supply(i) = 0._r8
      river_supply(i) = 0._r8
      runoff_supply(i) = 0._r8
      groundwater_supply_rn(i) = 0._r8
      groundwater_supply_nrn(i) = 0._r8
      
      !   zero irrigation at the begin of the new year
      if (idate(2) == 1 .and. idate(3) == deltim)then
         sum_irrig(i) = 0._r8
         sum_deficit_irrig(i) = 0._r8
         sum_irrig_count(i) = 0._r8
         zwt_stand(i) = zwt(i) + 1._r8
         irrig_demand_days(i) = 0
         irrig_supply_satisfy_days(i) = 0
         irrig_supply_unsatisfy_days(i) = 0
      end if

#ifdef CaMa_Flood
      !   irrigation withdraw from reservoir and river
      reservoirriver_supply(i) = withdrawal_cama(i)
      reservoir_supply(i) = withdrawal_dam_cama(i)
      river_supply(i) = withdrawal_riv_cama(i)
      runoff_supply(i) = withdrawal_rof_cama(i)
      waterstorage(i) = waterstorage(i) + reservoirriver_supply(i)
      withdrawal_cama(i) = 0._r8
      withdrawal_dam_cama(i) = 0._r8
      withdrawal_riv_cama(i) = 0._r8
      withdrawal_rof_cama(i) = 0._r8
#endif 

      ! !   calculate last day potential evapotranspiration
      ! call CalPotentialEvapotranspiration(i,idate,dlon,deltim)

      !   calculate whether irrigation needed
      call PointNeedsCheckForIrrig(i,ps,pe,idate,deltim,dlon,npcropmin,check_for_irrig)

      !   calculate irrigation needed
      if (check_for_irrig) then
         call CalIrrigationPotentialNeeded(i,ps,pe,nl_soil,nbedrock,zi_soi,dz_soi)
         call CalIrrigationLimitedSupply(i,nl_soil,deltim,dz_soi,zi_soi)
      endif

      !   calculate irrigation rate kg/m2->mm/s
      if((check_for_irrig) .and. (deficit_irrig(i) > 0)) then
         sum_deficit_irrig(i) = sum_deficit_irrig(i) + deficit_irrig(i)
      endif

      if ((check_for_irrig) .and. (actual_irrig(i) > 0)) then
         irrig_nsteps_per_day = nint(irrig_time_per_day/deltim)
         irrig_rate(i) = actual_irrig(i)/deltim/irrig_nsteps_per_day
         n_irrig_steps_left(i) = irrig_nsteps_per_day
         sum_irrig(i) = sum_irrig(i) + actual_irrig(i)
         sum_irrig_count(i) = sum_irrig_count(i) + 1._r8
      end if

      ! if(waterstorage(i)<0)then
      !    write(*,*) "LHB debug line121 waterstorage error : patch, waterstorage, deficit_irrig, gw_demand, gw_supply, rr_demand, rr_supply -----> ", i, waterstorage(i), deficit_irrig(i), groundwater_demand(i), groundwater_supply(i), reservoirriver_demand(i), reservoirriver_supply(i)
      ! endif

      if ((check_for_irrig) .and. (deficit_irrig(i) > 0)) then
         irrig_demand_days(i) = irrig_demand_days(i) + 1
         if (abs(actual_irrig(i) - deficit_irrig(i)) <= 1E-6) then 
            irrig_supply_satisfy_days(i) = irrig_supply_satisfy_days(i) + 1
         else
            irrig_supply_unsatisfy_days(i) = irrig_supply_unsatisfy_days(i) + 1
         endif
      end if

   end subroutine CalIrrigationNeeded


   subroutine CalIrrigationPotentialNeeded(i,ps,pe,nl_soil,nbedrock,zi_soi,dz_soi)

   !   DESCRIPTION:
   !   This subroutine is used to calculate how much irrigation needed in each irrigated crop patch without water supply restriction
   integer , intent(in) :: i
   integer , intent(in) :: ps, pe
   integer , intent(in) :: nbedrock
   integer , intent(in) :: nl_soil
   real(r8), intent(in) :: zi_soi(1:nl_soil)
   real(r8), intent(in) :: dz_soi(1:nl_soil)

   !   local variables
   integer  :: j
   integer  :: m
   logical  :: reached_max_depth
   real(r8) :: h2osoi_liq_tot
   real(r8) :: h2osoi_liq_target_tot
   real(r8) :: h2osoi_liq_wilting_point_tot
   real(r8) :: h2osoi_liq_field_capacity_tot
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
      h2osoi_liq_field_capacity_tot = 0._r8
      h2osoi_liq_saturation_capacity_tot = 0._r8

!   calculate wilting point and field capacity
      do j = 1, nl_soil
         if (t_soisno(j,i) > tfrz .and. porsl(j,i) >= 1.e-6) then
#ifdef Campbell_SOIL_MODEL
            h2osoi_liq_wilting_point(j) = denh2o*dz_soi(j)*porsl(j,i)*((smpswc/psi0(j,i))**(-1/bsw(j,i)))
            h2osoi_liq_field_capacity(j) = denh2o*dz_soi(j)*porsl(j,i)*((smpsfc/psi0(j,i))**(-1/bsw(j,i)))
            h2osoi_liq_saturation_capacity(j) = denh2o*dz_soi(j)*porsl(j,i)
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
            h2osoi_liq_wilting_point(j) = soil_vliq_from_psi(smpswc, porsl(j,i), theta_r(j,i), psi0(j,i), 5, &
              (/alpha_vgm(j,i), n_vgm(j,i), L_vgm(j,i), sc_vgm(j,i), fc_vgm(j,i)/))
            h2osoi_liq_wilting_point(j) = denh2o*dz_soi(j)*h2osoi_liq_wilting_point(j)
            h2osoi_liq_field_capacity(j) = soil_vliq_from_psi(smpsfc, porsl(j,i), theta_r(j,i), psi0(j,i), 5, &
              (/alpha_vgm(j,i), n_vgm(j,i), L_vgm(j,i), sc_vgm(j,i), fc_vgm(j,i)/))
            h2osoi_liq_field_capacity(j) = denh2o*dz_soi(j)*h2osoi_liq_field_capacity(j)
            h2osoi_liq_saturation_capacity(j) = denh2o*dz_soi(j)*porsl(j,i)
#endif
         endif
      enddo

      !   calculate total irrigation needed in all soil layers
      do m = ps, pe
         do j = 1, nl_soil
            if (.not. reached_max_depth) then
               if (zi_soi(j) > irrig_max_depth) then
                  reached_max_depth = .true.
               ! if (j > 3) then
               !    reached_max_depth = .true.
               elseif (j > nbedrock) then
                  reached_max_depth = .true.
               elseif (t_soisno(j,i) <= tfrz) then
                  reached_max_depth = .true.
               else
                  h2osoi_liq_tot = h2osoi_liq_tot + wliq_soisno(j,i)
                  h2osoi_liq_wilting_point_tot = h2osoi_liq_wilting_point_tot + h2osoi_liq_wilting_point(j)
                  h2osoi_liq_field_capacity_tot = h2osoi_liq_field_capacity_tot + h2osoi_liq_field_capacity(j)
                  h2osoi_liq_saturation_capacity_tot = h2osoi_liq_saturation_capacity_tot + h2osoi_liq_saturation_capacity(j)
               endif
            endif
         enddo
         if (irrig_method_p(m) == irrig_method_drip .or. irrig_method_p(m) == irrig_method_sprinkler .or. &
               irrig_method_p(m) == irrig_method_flood) then
               !  flood irrigation threshold at field capacity, but irrigation amount at saturation capacity
               h2osoi_liq_target_tot = h2osoi_liq_field_capacity_tot
         elseif (irrig_method_p(m) == irrig_method_paddy) then
               h2osoi_liq_target_tot = h2osoi_liq_saturation_capacity_tot
         else
               !  default irrigation is drip irrigation
               h2osoi_liq_target_tot = h2osoi_liq_field_capacity_tot
         endif
      enddo

      !   calculate irrigation threshold
      deficit_irrig(i) = 0._r8
      h2osoi_liq_at_threshold = h2osoi_liq_wilting_point_tot + irrig_threshold_fraction*(h2osoi_liq_target_tot - h2osoi_liq_wilting_point_tot)

      !   calculate total irrigation
      do m = ps, pe
         if (h2osoi_liq_tot < h2osoi_liq_at_threshold) then
            if (irrig_method_p(m) == irrig_method_sprinkler) then
               deficit_irrig(i) = irrig_supply_fraction*(h2osoi_liq_target_tot - h2osoi_liq_tot)
               ! deficit_irrig(i) = irrig_supply_fraction*h2osoi_liq_target_tot - h2osoi_liq_tot + potential_evapotranspiration(i)
            elseif (irrig_method_p(m) == irrig_method_flood) then
               deficit_irrig(i) = irrig_supply_fraction*(h2osoi_liq_saturation_capacity_tot - h2osoi_liq_tot)
            else
               deficit_irrig(i) = irrig_supply_fraction*(h2osoi_liq_target_tot - h2osoi_liq_tot)
            endif
         else
            deficit_irrig(i) = 0
         endif
      enddo

   end subroutine CalIrrigationPotentialNeeded

   subroutine CalIrrigationApplicationFluxes(i,deltim,qflx_irrig_drip,qflx_irrig_sprinkler,qflx_irrig_flood,qflx_irrig_paddy)
      !   DESCRIPTION:
      !   This subroutine is used to calculate irrigation application fluxes for each irrigated crop patch
      integer , intent(in) :: i
      real(r8), intent(in) :: deltim
      real(r8), intent(out):: qflx_irrig_drip,qflx_irrig_sprinkler,qflx_irrig_flood,qflx_irrig_paddy

      integer :: ps, pe, m


      ps = patch_pft_s(i)
      pe = patch_pft_e(i)

      qflx_irrig_drip = 0._r8
      qflx_irrig_sprinkler = 0._r8
      qflx_irrig_flood = 0._r8
      qflx_irrig_paddy = 0._r8

      !   add irrigation fluxes to precipitation or land surface
      do m = ps, pe
         if (n_irrig_steps_left(i) > 0) then
               n_irrig_steps_left(i) = n_irrig_steps_left(i) -1
               if (waterstorage(i) - irrig_rate(i)*deltim < 0._r8) irrig_rate(i) = waterstorage(i)/deltim
               waterstorage(i) = max(waterstorage(i) - irrig_rate(i)*deltim, 0._r8)
               if (irrig_method_p(m) == irrig_method_drip) then
                  qflx_irrig_drip = irrig_rate(i)
               else if (irrig_method_p(m) == irrig_method_sprinkler) then 
                  qflx_irrig_sprinkler = irrig_rate(i)
               else if (irrig_method_p(m) == irrig_method_flood) then
                  qflx_irrig_flood = irrig_rate(i)
               else if (irrig_method_p(m) == irrig_method_paddy) then
                  qflx_irrig_paddy = irrig_rate(i)
               else
                  qflx_irrig_sprinkler = irrig_rate(i)
               end if
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


      !  adjust flood irrigation in rice to paddy irrigaiton
      do m = ps, pe 
         ivt = pftclass(m)
         if ((ivt == 62) .and. (irrig_method_p(m) == irrig_method_flood))then
            irrig_method_p(m) = irrig_method_paddy
         endif
      enddo
      ! do m = ps, pe 
      !    irrig_method_p(m) = DEF_IRRIGATION_METHOD
      ! enddo
      

      do m = ps, pe
         ivt = pftclass(m)
         if ((ivt >= npcropmin) .and. (irrig_crop(ivt)) .and. &
            (cphase_p(m) >= irrig_min_cphase) .and. (cphase_p(m)<irrig_max_cphase)) then
            if (DEF_simulation_time%greenwich) then
                call gmt2local(idate, dlon, ldate)
                seconds_since_irrig_start_time = ldate(3) - irrig_start_time + deltim
            else
                seconds_since_irrig_start_time = idate(3) - irrig_start_time + deltim
            endif
            if ((seconds_since_irrig_start_time >= 0._r8) .and. (seconds_since_irrig_start_time < deltim)) then
                check_for_irrig = .true.
            else
                check_for_irrig = .false.
            endif
         else
            check_for_irrig = .false.
         endif
      enddo

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
   !     endif

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
   !     endif
   ! end subroutine CalPotentialEvapotranspiration


   subroutine CalIrrigationLimitedSupply(i,nl_soil,deltim,dz_soi,zi_soi)
        !   DESCRIPTION:
        !   This subroutine is used to calculate how much irrigation supplied in each irrigated crop patch with water supply restriction
         integer,  intent(in) :: i
         integer,  intent(in) :: nl_soil      
         real(r8), intent(in) :: deltim       
         real(r8), intent(in) :: dz_soi(1:nl_soil)
         real(r8), intent(in) :: zi_soi(1:nl_soil)
         
         if (deficit_irrig(i) > 0._r8) then
            if (DEF_IRRIGATION_ALLOCATION == 1) then
               actual_irrig(i) = deficit_irrig(i)
               waterstorage(i) = waterstorage(i) + actual_irrig(i)
            elseif (DEF_IRRIGATION_ALLOCATION == 2) then
               waterstorage_supply(i) = min(waterstorage(i), deficit_irrig(i))
               waterstorage_supply(i) = max(waterstorage_supply(i), 0._r8)
               actual_irrig(i) = actual_irrig(i) + waterstorage_supply(i)
#ifdef CaMa_Flood
               if (deficit_irrig(i) > actual_irrig(i)) then
                  reservoirriver_demand(i) = max(deficit_irrig(i) - actual_irrig(i), 0._r8)
               end if
#endif
               !   irrigation withdraw from ground water (unconfined and confined)
               if (deficit_irrig(i) > actual_irrig(i)) then
                  groundwater_demand(i) = max(deficit_irrig(i) - actual_irrig(i), 0._r8)
                  call CalGroudwaterWithdrawal(i,nl_soil,deltim,dz_soi,zi_soi)
                  actual_irrig(i) = actual_irrig(i) + groundwater_supply(i)
                  waterstorage(i) = waterstorage(i) + groundwater_supply(i)
               endif
            elseif (DEF_IRRIGATION_ALLOCATION == 3) then
               waterstorage_supply(i) = min(waterstorage(i), deficit_irrig(i))
               waterstorage_supply(i) = max(waterstorage_supply(i), 0._r8)
               actual_irrig(i) = actual_irrig(i) + waterstorage_supply(i)
               ! write(*,*) "LHB debug line417 irrig_alloc error : patch, irrig_sw_alloc, irrig_gw_alloc -----> ", i, irrig_sw_alloc(i), irrig_gw_alloc(i)
#ifdef CaMa_Flood
               if (deficit_irrig(i) > actual_irrig(i)) then
                  reservoirriver_demand(i) = max((deficit_irrig(i) - actual_irrig(i))*(1.-irrig_gw_alloc(i)), 0._r8)
               end if
#endif
               !   irrigation withdraw from ground water (unconfined and confined)
               if (deficit_irrig(i) > actual_irrig(i)) then
                  groundwater_demand(i) = max((deficit_irrig(i) - actual_irrig(i))*irrig_gw_alloc(i), 0._r8)
                  call CalGroudwaterWithdrawal(i,nl_soil,deltim,dz_soi,zi_soi)
                  actual_irrig(i) = actual_irrig(i) + groundwater_supply(i)
                  waterstorage(i) = waterstorage(i) + groundwater_supply(i)
               endif     
            endif    
         endif
   end subroutine CalIrrigationLimitedSupply


   subroutine CalGroudwaterWithdrawal(i,nl_soil,deltim,dz_soi,zi_soi)
   !   DESCRIPTION:
   !   This subroutine is used to calculate irrigation withdrawals for groudwater
      integer,  intent(in) :: i
      integer,  intent(in) :: nl_soil      
      real(r8), intent(in) :: deltim       
      real(r8), intent(in) :: dz_soi(1:nl_soil)
      real(r8), intent(in) :: zi_soi(1:nl_soil)
      
      if (.not. DEF_USE_VariablySaturatedFlow) then
         call CalWithdrawalWATER(i,nl_soil,deltim,dz_soi,zi_soi)
         ! groundwater_supply_rn(i) = groundwater_supply(i)
         ! groundwater_supply_nrn(i) =  max(groundwater_demand(i) - groundwater_supply(i), 0._r8)
         ! groundwater_supply(i) = groundwater_demand(i)
      else
         groundwater_supply(i) = groundwater_demand(i)
      endif
   end subroutine CalGroudwaterWithdrawal



   subroutine CalWithdrawalWATER(i,nl_soil,deltim,dz_soi,zi_soi)
      !   DESCRIPTION:
      !   This subroutine is used to calculate how much irrigation supplied in each irrigated crop patch with groundwater supply restriction
      IMPLICIT NONE
      integer,  INTENT(in) :: i
      integer,  INTENT(in) :: nl_soil      
      real(r8), INTENT(in) :: deltim              ! land model time step (sec)
      real(r8), INTENT(in) :: dz_soi  (1:nl_soil) ! layer depth (m)
      real(r8), INTENT(in) :: zi_soi  (1:nl_soil) ! interface level below a "z" level (m)

      ! LOCAL ARGUMENTS
      integer  :: j                ! indices
      integer  :: jwt              ! index of the soil layer right above the water table (-)
      real(r8) :: dzmm(1:nl_soil)  ! layer thickness (mm)
      real(r8) :: vol_ice(1:nl_soil)! partitial volume of ice lens in layer
      real(r8) :: eff_porosity(1:nl_soil)! effective porosity = porosity - vol_ice
      real(r8) :: xs               ! water needed to bring soil moisture to watmin (mm)
      real(r8) :: xsi              ! excess soil water above saturation at layer i (mm)
      real(r8) :: xs1              ! excess soil water above saturation at layer 1 (mm)    
      real(r8) :: pump_total     
      real(r8) :: pump_layer   
      real(r8) :: max_groundwater_supply
      real(r8) :: s_y              
      real(r8) :: rous             ! specific yield [-]
   ! -------------------------------------------------------------------------

      do j = 1, nl_soil
         vol_ice(j) = min(porsl(j,i), wice_soisno(j,i)/(dz_soi(j)*denice))
         eff_porosity(j) = max(0.01, porsl(j,i)-vol_ice(j))
      end do
      ! Convert layer thicknesses from m to mm
      do j = 1,nl_soil
         dzmm(j) = dz_soi(j)*1000.
      end do
      ! The layer index of the first unsaturated layer,
      ! i.e., the layer right above the water table
      jwt = nl_soil
      ! allow jwt to equal zero when zwt is in top layer
      do j = 1, nl_soil
         if(zwt(i) <= zi_soi(j)) then
         jwt = j-1
         exit
         end if
      enddo

      rous = porsl(nl_soil,i)*(1.-(1.-1.e3*zwt(i)/psi0(nl_soil,i))**(-1./bsw(nl_soil,i)))
      rous = max(rous,0.02)
   !-- Water table is below the soil column  ----------------------------------------
      if (jwt == nl_soil) then
         max_groundwater_supply = max(1.e3*(zwt_stand(i)-zwt(i))*rous, 0._r8)
         groundwater_supply(i) = min(groundwater_demand(i), max_groundwater_supply)
         wa(i) = wa(i) - groundwater_supply(i)
         zwt(i) = max(0., zwt(i) + groundwater_supply(i)/1000./rous)
         wliq_soisno(nl_soil,i) = wliq_soisno(nl_soil,i) + max(0.,(wa(i)-5000.))
         wa(i) = min(wa(i), 5000.)
      else
   !-- Water table within soil layers 1-9  ------------------------------------------
   !============================== RSUB_TOP =========================================
         !-- Now remove water via pump
         pump_total = - groundwater_demand(i)
         do j = jwt+1, nl_soil
               ! use analytical expression for specific yield
               s_y = porsl(j,i) * ( 1. - (1.-1.e3*zwt(i)/psi0(j,i))**(-1./bsw(j,i)))
               s_y = max(s_y,0.02)
               pump_layer = max(pump_total, -(s_y*(zi_soi(j)-zwt(i))*1.e3))
               pump_layer = min(pump_layer, 0.)
               wliq_soisno(j,i) = wliq_soisno(j,i) + pump_layer
               pump_total = pump_total - pump_layer
               groundwater_supply(i) = groundwater_supply(i) - pump_layer
               if(pump_total >= 0.)then
                  zwt(i) = max(0.,zwt(i) - pump_layer/s_y/1000.)
                  exit
               else
                  zwt(i) = zi_soi(j)
               endif
         enddo
   !-- Remove residual drainage  ------------------------------------------------
         max_groundwater_supply = max(1.e3*(zwt_stand(i)-zwt(i))*rous, 0._r8)
         pump_total = min(-pump_total, max_groundwater_supply)
         pump_total = max(pump_total, 0._r8)
         groundwater_supply(i) = groundwater_supply(i) + pump_total
         zwt(i) = max(0.,zwt(i) + pump_total/1000./rous)
         wa(i) = wa(i) - pump_total

   !-- Recompute jwt  ---------------------------------------------------------------
         ! allow jwt to equal zero when zwt is in top layer
         jwt = nl_soil
         do j = 1, nl_soil
               if(zwt(i) <= zi_soi(j)) then
                  jwt = j-1
                  exit
               end if
         enddo
      end if   ! end of jwt if construct

      zwt(i) = max(0.0,zwt(i))
      zwt(i) = min(80.,zwt(i))

      ! Correction [1]
      ! NON-physically based corection on wliq_soisno
      ! excessive water above saturation added to the above unsaturated layer like a bucket
      ! if column over saturated, excess water goes to runoff
      do j = nl_soil,2,-1
         xsi = max(wliq_soisno(j,i)-eff_porosity(j)*dzmm(j),0.)
         wliq_soisno(j,i) = min(eff_porosity(j)*dzmm(j), wliq_soisno(j,i))
         wliq_soisno(j-1,i) = wliq_soisno(j-1,i) + xsi
      end do
      ! 12/2022, note by yuan: a potential bug below which needs check,
      ! if wice_soisno(1) > pondmx + porsl*dzmm, so xs1>0, in that case,
      ! wliq_soisno(1) will be nagtive, and xs1 is positive.
      xs1 = wliq_soisno(1,i) - (pondmx+porsl(1,i)*dzmm(1)-wice_soisno(1,i))
      if(xs1 > 0.)then
               wliq_soisno(1,i) = pondmx+porsl(1,i)*dzmm(1)-wice_soisno(1,i)
      else
               xs1 = 0.
      endif
      wa(i) = wa(i) + xs1

      ! Correction [2]
      ! NON-physically based corection on wliq_soisno
      ! Limit wliq_soisno to be greater than or equal to watmin.
      ! Get water needed to bring wliq_soisno equal watmin from lower layer.
      ! If insufficient water in soil layers, get from aquifer water
      xs = 0.
      do j = 1, nl_soil
         if (wliq_soisno(j,i) < 0.) then
               xs = xs + wliq_soisno(j,i)
               wliq_soisno(j,i) = 0.
         endif
      enddo
      wa(i) = wa(i) + xs
   end subroutine CalWithdrawalWATER

end MODULE MOD_Irrigation
#endif
