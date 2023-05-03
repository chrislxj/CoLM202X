#include <define.h>

MODULE MOD_2D_BGCFluxes
! ----------------------------------------------------------------------
! perfrom the grid average mapping: average a subgrid input 1d vector 
! of length numpatch to a output 2d array of length [lon_points,lat_points]
!
! Created by Yongjiu Dai, 03/2014
!---------------------------------------------------------------------
#ifdef BGC

   use mod_data_type
   USE GlobalVars

   IMPLICIT NONE
   SAVE

   type(block_data_real8_2d) :: f_leafc              ! leaf carbon display pool  (gC m-2)
   type(block_data_real8_2d) :: f_leafc_storage      ! leaf carbon storage pool  (gC m-2)
   type(block_data_real8_2d) :: f_leafc_xfer         ! leaf carbon transfer pool (gC m-2)
   type(block_data_real8_2d) :: f_frootc             ! fine root carbon display pool  (gC m-2)
   type(block_data_real8_2d) :: f_frootc_storage     ! fine root carbon storage pool  (gC m-2)
   type(block_data_real8_2d) :: f_frootc_xfer        ! fine root carbon transfer pool (gC m-2)
   type(block_data_real8_2d) :: f_livestemc          ! live stem carbon display pool  (gC m-2)
   type(block_data_real8_2d) :: f_livestemc_storage  ! live stem carbon storage pool  (gC m-2)
   type(block_data_real8_2d) :: f_livestemc_xfer     ! live stem carbon transfer pool (gC m-2)
   type(block_data_real8_2d) :: f_deadstemc          ! dead stem carbon display pool  (gC m-2)
   type(block_data_real8_2d) :: f_deadstemc_storage  ! dead stem carbon storage pool  (gC m-2)
   type(block_data_real8_2d) :: f_deadstemc_xfer     ! dead stem carbon transfer pool (gC m-2)
   type(block_data_real8_2d) :: f_livecrootc         ! live coarse root carbon display pool  (gC m-2)
   type(block_data_real8_2d) :: f_livecrootc_storage ! live coarse root carbon storage pool  (gC m-2)
   type(block_data_real8_2d) :: f_livecrootc_xfer    ! live coarse root carbon transfer pool (gC m-2)
   type(block_data_real8_2d) :: f_deadcrootc         ! dead coarse root carbon display pool  (gC m-2)
   type(block_data_real8_2d) :: f_deadcrootc_storage ! dead coarse root carbon storage pool  (gC m-2)
   type(block_data_real8_2d) :: f_deadcrootc_xfer    ! dead coarse root carbon transfer pool (gC m-2)
#ifdef CROP
   type(block_data_real8_2d) :: f_grainc             ! grain carbon display pool (gC m-2)
   type(block_data_real8_2d) :: f_grainc_storage     ! grain carbon storage pool (gC m-2)
   type(block_data_real8_2d) :: f_grainc_xfer        ! grain carbon transfer pool (gC m-2)
#endif
   type(block_data_real8_2d) :: f_leafn              ! leaf nitrogen display pool  (gN m-2)
   type(block_data_real8_2d) :: f_leafn_storage      ! leaf nitrogen storage pool  (gN m-2)
   type(block_data_real8_2d) :: f_leafn_xfer         ! leaf nitrogen transfer pool (gN m-2)
   type(block_data_real8_2d) :: f_frootn             ! fine root nitrogen display pool  (gN m-2)
   type(block_data_real8_2d) :: f_frootn_storage     ! fine root nitrogen storage pool  (gN m-2)
   type(block_data_real8_2d) :: f_frootn_xfer        ! fine root nitrogen transfer pool (gN m-2)
   type(block_data_real8_2d) :: f_livestemn          ! live stem nitrogen display pool  (gN m-2)
   type(block_data_real8_2d) :: f_livestemn_storage  ! live stem nitrogen storage pool  (gN m-2)
   type(block_data_real8_2d) :: f_livestemn_xfer     ! live stem nitrogen transfer pool (gN m-2)
   type(block_data_real8_2d) :: f_deadstemn          ! dead stem nitrogen display pool  (gN m-2)
   type(block_data_real8_2d) :: f_deadstemn_storage  ! dead stem nitrogen storage pool  (gN m-2)
   type(block_data_real8_2d) :: f_deadstemn_xfer     ! dead stem nitrogen transfer pool (gN m-2)
   type(block_data_real8_2d) :: f_livecrootn         ! live coarse root nitrogen display pool  (gN m-2)
   type(block_data_real8_2d) :: f_livecrootn_storage ! live coarse root nitrogen storage pool  (gN m-2)
   type(block_data_real8_2d) :: f_livecrootn_xfer    ! live coarse root nitrogen transfer pool (gN m-2)
   type(block_data_real8_2d) :: f_deadcrootn         ! dead coarse root nitrogen display pool  (gN m-2)
   type(block_data_real8_2d) :: f_deadcrootn_storage ! dead coarse root nitrogen storage pool  (gN m-2)
   type(block_data_real8_2d) :: f_deadcrootn_xfer    ! dead coarse root nitrogen transfer pool (gN m-2)

#ifdef CROP
   type(block_data_real8_2d) :: f_grainn             ! grain nitrogen display pool (gN m-2)
   type(block_data_real8_2d) :: f_grainn_storage     ! grain nitrogen storage pool (gN m-2)
   type(block_data_real8_2d) :: f_grainn_xfer        ! grain nitrogen transfer pool (gN m-2)
#endif
   type(block_data_real8_2d) :: f_retransn           ! retranslocation nitrogen pool (gN m-2)

#ifdef CROP
   type(block_data_real8_2d) :: f_cphase             ! crop phase
   type(block_data_real8_2d) :: f_cropprod1c         ! 1-yr crop production carbon (gC m-2)
   type(block_data_real8_2d) :: f_cropprod1c_loss    ! loss of 1-yr crop production carbon (gC m-2 s-1)
   type(block_data_real8_2d) :: f_cropseedc_deficit  ! crop seed carbon deficit (gC m-2 s-1)
   type(block_data_real8_2d) :: f_grainc_to_cropprodc! grain to crop production carbon (gC m-2 s-1)
   type(block_data_real8_2d) :: f_grainc_to_seed     ! grain to crop seed carbon (gC m-2 s-1)
   type(block_data_real8_2d) :: f_fert_to_sminn      ! fertilization (gN m-2 s-1)
   type(block_data_real8_2d) :: f_plantdate          ! planting date
#endif
   type(block_data_real8_2d) :: f_ndep_to_sminn      ! nitrogen deposition (gN m-2 s-1)

   type(block_data_real8_2d) :: f_gpp                ! net primary production (gC m-2 s-1)
   type(block_data_real8_2d) :: f_downreg            ! gpp downregulation due to N limitation
   type(block_data_real8_2d) :: f_ar                 ! autotrophic respiration (gC m-2 s-1)
   type(block_data_real8_2d) :: f_cwdprod            ! CWD production (gC m-2 s-1)
   type(block_data_real8_2d) :: f_cwddecomp          ! CWD decomposition (gC m-2 s-1)
   type(block_data_real8_2d) :: f_hr                 ! heterotrophic respiration (gC m-2 s-1)
   type(block_data_real8_2d) :: f_fpg                ! fraction of gpp potential
   type(block_data_real8_2d) :: f_fpi                ! fraction of immobalization
   type(block_data_real8_2d) :: f_gpp_enftemp        ! gross primary productivity for needleleaf evergreen temperate tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_enfboreal      ! gross primary productivity for needleleaf evergreen boreal tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_dnfboreal      ! gross primary productivity for needleleaf deciduous boreal tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_ebftrop        ! gross primary productivity for broadleaf evergreen tropical tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_ebftemp        ! gross primary productivity for broadleaf evergreen temperate tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_dbftrop        ! gross primary productivity for broadleaf deciduous tropical tree (gC m-2 s-1) 
   type(block_data_real8_2d) :: f_gpp_dbftemp        ! gross primary productivity for broadleaf deciduous temperate tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_dbfboreal      ! gross primary productivity for broadleaf deciduous boreal tree (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_ebstemp        ! gross primary productivity for broadleaf evergreen temperate shrub (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_dbstemp        ! gross primary productivity for broadleaf deciduous temperate shrub (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_dbsboreal      ! gross primary productivity for broadleaf deciduous boreal shrub (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_c3arcgrass     ! gross primary productivity for c3 arctic grass (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_c3grass        ! gross primary productivity for c3 grass (gC m-2 s-1)
   type(block_data_real8_2d) :: f_gpp_c4grass        ! gross primary productivity for c4 grass (gC m-2 s-1)
   type(block_data_real8_2d) :: f_leafc_enftemp      ! leaf carbon display pool for needleleaf evergreen temperate tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_enfboreal    ! leaf carbon display pool for needleleaf evergreen boreal tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_dnfboreal    ! leaf carbon display pool for needleleaf deciduous boreal tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_ebftrop      ! leaf carbon display pool for broadleaf evergreen tropical tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_ebftemp      ! leaf carbon display pool for broadleaf evergreen temperate tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_dbftrop      ! leaf carbon display pool for broadleaf deciduous tropical tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_dbftemp      ! leaf carbon display pool for broadleaf deciduous temperate tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_dbfboreal    ! leaf carbon display pool for broadleaf deciduous boreal tree (gC m-2)
   type(block_data_real8_2d) :: f_leafc_ebstemp      ! leaf carbon display pool for broadleaf evergreen temperate shrub (gC m-2)
   type(block_data_real8_2d) :: f_leafc_dbstemp      ! leaf carbon display pool for broadleaf deciduous temperate shrub (gC m-2)
   type(block_data_real8_2d) :: f_leafc_dbsboreal    ! leaf carbon display pool for broadleaf deciduous boreal shrub (gC m-2)
   type(block_data_real8_2d) :: f_leafc_c3arcgrass   ! leaf carbon display pool for c3 arctic grass (gC m-2)
   type(block_data_real8_2d) :: f_leafc_c3grass      ! leaf carbon display pool for c3 grass (gC m-2)
   type(block_data_real8_2d) :: f_leafc_c4grass      ! leaf carbon display pool for c4 grass (gC m-2)

   type(block_data_real8_3d) :: f_litr1c_vr          !! litter 1 (metabolic litter) carbon density in soil layers (gC m-2)
   type(block_data_real8_3d) :: f_litr2c_vr          !! litter 2 (cellulosic litte) carbon density in soil layers (gC m-2)
   type(block_data_real8_3d) :: f_litr3c_vr          !! litter 3 (lignin litter) carbon density in soil layers (gC m-2)
   type(block_data_real8_3d) :: f_cwdc_vr            !！coarse woody debris carbon density in soil layers (gC m-2)
   type(block_data_real8_3d) :: f_soil1c_vr          !! soil 1 (active soil organic matter) carbon density in soil layers (gC m-2)
   type(block_data_real8_3d) :: f_soil2c_vr          !! soil 2 (slow soil organic matter) carbon density in soil layers (gC m-2)
   type(block_data_real8_3d) :: f_soil3c_vr          !! soil 3 (passive soil organic matter) carbon density in soil layers (gC m-2)

   type(block_data_real8_3d) :: f_litr1n_vr          !! litter 1 (metabolic litter) carbon density in soil layers (gN m-2)
   type(block_data_real8_3d) :: f_litr2n_vr          !! litter 2 (cellulosic litte) carbon density in soil layers (gN m-2)
   type(block_data_real8_3d) :: f_litr3n_vr          !! litter 3 (lignin litter) carbon density in soil layers (gN m-2)
   type(block_data_real8_3d) :: f_cwdn_vr            !！coarse woody debris carbon density in soil layers (gN m-2)
   type(block_data_real8_3d) :: f_soil1n_vr          !! soil 1 (active soil organic matter) carbon density in soil layers (gN m-2)
   type(block_data_real8_3d) :: f_soil2n_vr          !! soil 2 (slow soil organic matter) carbon density in soil layers (gN m-2)
   type(block_data_real8_3d) :: f_soil3n_vr          !! soil 3 (passive soil organic matter) carbon density in soil layers (gN m-2)
   type(block_data_real8_3d) :: f_sminn_vr           ! mineral nitrogen density in soil layers (gN m-2)

#ifdef NITRIF
   type(block_data_real8_3d) :: f_O2_DECOMP_DEPTH_UNSAT  ! O2 consumption from heterotrophic respiration and autotrophic respiration for non-inundated area (mol m-3 s-1)
   type(block_data_real8_3d) :: f_CONC_O2_UNSAT          ! O2 consumption from heterotrophic respiration and autotrophic respiration for non-inundated area (mol m-3 s-1)
#endif
#ifdef CROP
   type(block_data_real8_2d) :: f_hui                    ! heat unit index
   type(block_data_real8_2d) :: f_vf                     ! vernalization response
   type(block_data_real8_2d) :: f_gddplant               ! gdd since planting (ddays)
   type(block_data_real8_2d) :: f_gddmaturity            ! gdd needed to harvest (ddays)
   type(block_data_real8_2d) :: f_pdcorn                 ! planting date of corn
   type(block_data_real8_2d) :: f_pdswheat               ! planting date of spring wheat
   type(block_data_real8_2d) :: f_pdwwheat               ! planting date of winter wheat
   type(block_data_real8_2d) :: f_pdsoybean              ! planting date of soybean
   type(block_data_real8_2d) :: f_pdcotton               ! planting date of cotton
   type(block_data_real8_2d) :: f_pdrice1                ! planting date of rice1
   type(block_data_real8_2d) :: f_pdrice2                ! planting date of rice2
   type(block_data_real8_2d) :: f_pdsugarcane            ! planting date of sugarcane
   type(block_data_real8_2d) :: f_fertnitro_corn         ! nitrogen fertilizer for corn (gN m-2)
   type(block_data_real8_2d) :: f_fertnitro_swheat       ! nitrogen fertilizer for spring wheat (gN m-2)
   type(block_data_real8_2d) :: f_fertnitro_wwheat       ! nitrogen fertilizer for winter wheat (gN m-2)
   type(block_data_real8_2d) :: f_fertnitro_soybean      ! nitrogen fertilizer for soybean (gN m-2)
   type(block_data_real8_2d) :: f_fertnitro_cotton       ! nitrogen fertilizer for cotton (gN m-2)
   type(block_data_real8_2d) :: f_fertnitro_rice1        ! nitrogen fertilizer for rice1 (gN m-2)
   type(block_data_real8_2d) :: f_fertnitro_rice2        ! nitrogen fertilizer for rice2 (gN m-2)
   type(block_data_real8_2d) :: f_fertnitro_sugarcane    ! nitrogen fertilizer for sugarcane (gN m-2)
#endif
#ifdef Fire
   type(block_data_real8_2d) :: f_abm                    ! peak crop fire month
   type(block_data_real8_2d) :: f_gdp                    !!! global gdp
   type(block_data_real8_2d) :: f_peatf                  !!! global peatf data
   type(block_data_real8_2d) :: f_hdm                    !!! human population density (counts km-2)
   type(block_data_real8_2d) :: f_lnfm                   !!! lightning frequency (counts km-2 hr-1)
#endif

   ! PUBLIC MEMBER FUNCTIONS:
   public :: allocate_2D_BGCFluxes

CONTAINS

   SUBROUTINE allocate_2D_BGCFluxes (grid)
      ! --------------------------------------------------------------------
      ! Allocates memory for CLM 2d [lon_points,lat_points] variables
      ! --------------------------------------------------------------------

      use spmd_task
      use mod_grid
      use mod_data_type
      implicit none

      type(grid_type), intent(in) :: grid

      if (p_is_io) then
         
         call allocate_block_data (grid, f_leafc              ) ! leaf carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_leafc_storage      ) ! leaf carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_leafc_xfer         ) ! leaf carbon transfer pool (gC/m2)
         call allocate_block_data (grid, f_frootc             ) ! fine root carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_frootc_storage     ) ! fine root carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_frootc_xfer        ) ! fine root carbon transfer pool (gC/m2)
         call allocate_block_data (grid, f_livestemc          ) ! live stem carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_livestemc_storage  ) ! live stem carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_livestemc_xfer     ) ! live stem carbon transfer pool (gC/m2)
         call allocate_block_data (grid, f_deadstemc          ) ! dead stem carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_deadstemc_storage  ) ! dead stem carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_deadstemc_xfer     ) ! dead stem carbon transfer pool (gC/m2)
         call allocate_block_data (grid, f_livecrootc         ) ! live coarse root carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_livecrootc_storage ) ! live coarse root carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_livecrootc_xfer    ) ! live coarse root carbon transfer pool (gC/m2)
         call allocate_block_data (grid, f_deadcrootc         ) ! dead coarse root carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_deadcrootc_storage ) ! dead coarse root carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_deadcrootc_xfer    ) ! dead coarse root carbon transfer pool (gC/m2)
#ifdef CROP
         call allocate_block_data (grid, f_grainc             ) ! grain carbon display pool  (gC/m2)
         call allocate_block_data (grid, f_grainc_storage     ) ! grain carbon storage pool  (gC/m2)
         call allocate_block_data (grid, f_grainc_xfer        ) ! grain carbon transfer pool (gC/m2)
#endif
         call allocate_block_data (grid, f_leafn              ) ! leaf nitrogen display pool  (gN/m2)
         call allocate_block_data (grid, f_leafn_storage      ) ! leaf nitrogen storage pool  (gN/m2)
         call allocate_block_data (grid, f_leafn_xfer         ) ! leaf nitrogen transfer pool (gN/m2)
         call allocate_block_data (grid, f_frootn             ) ! fine root nitrogen display pool  (gN/m2)
         call allocate_block_data (grid, f_frootn_storage     ) ! fine root nitrogen storage pool  (gN/m2)
         call allocate_block_data (grid, f_frootn_xfer        ) ! fine root nitrogen transfer pool (gN/m2)
         call allocate_block_data (grid, f_livestemn          ) ! live stem nitrogen display pool  (gN/m2)
         call allocate_block_data (grid, f_livestemn_storage  ) ! live stem nitrogen storage pool  (gN/m2)
         call allocate_block_data (grid, f_livestemn_xfer     ) ! live stem nitrogen transfer pool (gN/m2)
         call allocate_block_data (grid, f_deadstemn          ) ! dead stem nitrogen display pool  (gN/m2)
         call allocate_block_data (grid, f_deadstemn_storage  ) ! dead stem nitrogen storage pool  (gN/m2)
         call allocate_block_data (grid, f_deadstemn_xfer     ) ! dead stem nitrogen transfer pool (gN/m2)
         call allocate_block_data (grid, f_livecrootn         ) ! live coarse root nitrogen display pool  (gN/m2)
         call allocate_block_data (grid, f_livecrootn_storage ) ! live coarse root nitrogen storage pool  (gN/m2)
         call allocate_block_data (grid, f_livecrootn_xfer    ) ! live coarse root nitrogen transfer pool (gN/m2)
         call allocate_block_data (grid, f_deadcrootn         ) ! dead coarse root nitrogen display pool  (gN/m2)
         call allocate_block_data (grid, f_deadcrootn_storage ) ! dead coarse root nitrogen storage pool  (gN/m2)
         call allocate_block_data (grid, f_deadcrootn_xfer    ) ! dead coarse root nitrogen transfer pool (gN/m2)

#ifdef CROP
         call allocate_block_data (grid, f_grainn             ) ! grain nitrogen display pool  (gN/m2)
         call allocate_block_data (grid, f_grainn_storage     ) ! grain nitrogen storage pool  (gN/m2)
         call allocate_block_data (grid, f_grainn_xfer        ) ! grain nitrogen transfer pool (gN/m2)
#endif
         call allocate_block_data (grid, f_retransn           ) ! retranslocation nitrogen pool (gN/m2)

#ifdef CROP
         call allocate_block_data (grid, f_cphase             )  ! crop phase
         call allocate_block_data (grid, f_cropprod1c         )  ! 1-yr crop production carbon
         call allocate_block_data (grid, f_cropprod1c_loss    )  ! loss of 1-yr crop production carbon
         call allocate_block_data (grid, f_cropseedc_deficit  )  ! crop seed carbon deficit
         call allocate_block_data (grid, f_grainc_to_cropprodc ) ! grain to crop production
         call allocate_block_data (grid, f_grainc_to_seed     )  ! grain to crop seed
         call allocate_block_data (grid, f_fert_to_sminn      )  ! grain to crop seed
         call allocate_block_data (grid, f_plantdate          )  ! planting date
#endif
         call allocate_block_data (grid, f_ndep_to_sminn      )  ! grain to crop seed

         call allocate_block_data (grid, f_gpp                ) ! net primary production (gC/m2)
         call allocate_block_data (grid, f_downreg            ) ! gpp downregulation due to N limitation
         call allocate_block_data (grid, f_ar                 )
         call allocate_block_data (grid, f_cwdprod            )
         call allocate_block_data (grid, f_cwddecomp          )
         call allocate_block_data (grid, f_hr                 )
         call allocate_block_data (grid, f_fpg                ) ! fraction of potential gpp
         call allocate_block_data (grid, f_fpi                ) ! fraction of potential immobilization
         call allocate_block_data (grid, f_gpp_enftemp        ) !1
         call allocate_block_data (grid, f_gpp_enfboreal      ) !2
         call allocate_block_data (grid, f_gpp_dnfboreal      ) !3
         call allocate_block_data (grid, f_gpp_ebftrop        ) !4
         call allocate_block_data (grid, f_gpp_ebftemp        ) !5
         call allocate_block_data (grid, f_gpp_dbftrop        ) !6
         call allocate_block_data (grid, f_gpp_dbftemp        ) !7
         call allocate_block_data (grid, f_gpp_dbfboreal      ) !8
         call allocate_block_data (grid, f_gpp_ebstemp        ) !9
         call allocate_block_data (grid, f_gpp_dbstemp        ) !10
         call allocate_block_data (grid, f_gpp_dbsboreal      ) !11
         call allocate_block_data (grid, f_gpp_c3arcgrass     ) !12
         call allocate_block_data (grid, f_gpp_c3grass        ) !13
         call allocate_block_data (grid, f_gpp_c4grass        ) !14
         call allocate_block_data (grid, f_leafc_enftemp      ) !1
         call allocate_block_data (grid, f_leafc_enfboreal    ) !2
         call allocate_block_data (grid, f_leafc_dnfboreal    ) !3
         call allocate_block_data (grid, f_leafc_ebftrop      ) !4
         call allocate_block_data (grid, f_leafc_ebftemp      ) !5
         call allocate_block_data (grid, f_leafc_dbftrop      ) !6
         call allocate_block_data (grid, f_leafc_dbftemp      ) !7
         call allocate_block_data (grid, f_leafc_dbfboreal    ) !8
         call allocate_block_data (grid, f_leafc_ebstemp      ) !9
         call allocate_block_data (grid, f_leafc_dbstemp      ) !10
         call allocate_block_data (grid, f_leafc_dbsboreal    ) !11
         call allocate_block_data (grid, f_leafc_c3arcgrass   ) !12
         call allocate_block_data (grid, f_leafc_c3grass      ) !13
         call allocate_block_data (grid, f_leafc_c4grass      ) !14

         call allocate_block_data (grid, f_litr1c_vr  ,nl_soil)  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_litr2c_vr  ,nl_soil)  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_litr3c_vr  ,nl_soil)  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_cwdc_vr    ,nl_soil)  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_soil1c_vr  ,nl_soil)  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_soil2c_vr  ,nl_soil)  ! soil carbon pool (gC/m2)
         call allocate_block_data (grid, f_soil3c_vr  ,nl_soil)  ! soil carbon pool (gC/m2)

         call allocate_block_data (grid, f_litr1n_vr  ,nl_soil)  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_litr2n_vr  ,nl_soil)  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_litr3n_vr  ,nl_soil)  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_cwdn_vr    ,nl_soil)  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_soil1n_vr  ,nl_soil)  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_soil2n_vr  ,nl_soil)  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_soil3n_vr  ,nl_soil)  ! soil nitrogen pool (gN/m2)
         call allocate_block_data (grid, f_sminn_vr   ,nl_soil)  ! soil mineral nitrogen pool (gN/m2)
#ifdef NITRIF
         call allocate_block_data (grid, f_O2_DECOMP_DEPTH_UNSAT, nl_soil)
         call allocate_block_data (grid, f_CONC_O2_UNSAT        , nl_soil)
#endif
#ifdef CROP
         call allocate_block_data (grid, f_hui                 )
         call allocate_block_data (grid, f_vf                  )
         call allocate_block_data (grid, f_gddmaturity         ) 
         call allocate_block_data (grid, f_gddplant            )   
         call allocate_block_data (grid, f_pdcorn              )
         call allocate_block_data (grid, f_pdswheat            )
         call allocate_block_data (grid, f_pdwwheat            )
         call allocate_block_data (grid, f_pdsoybean           )
         call allocate_block_data (grid, f_pdcotton            )
         call allocate_block_data (grid, f_pdrice1             )
         call allocate_block_data (grid, f_pdrice2             )
         call allocate_block_data (grid, f_pdsugarcane         )
         call allocate_block_data (grid, f_fertnitro_corn      )
         call allocate_block_data (grid, f_fertnitro_swheat    )
         call allocate_block_data (grid, f_fertnitro_wwheat    )
         call allocate_block_data (grid, f_fertnitro_soybean   )
         call allocate_block_data (grid, f_fertnitro_cotton    )
         call allocate_block_data (grid, f_fertnitro_rice1     )
         call allocate_block_data (grid, f_fertnitro_rice2     )
         call allocate_block_data (grid, f_fertnitro_sugarcane )
#endif
#ifdef Fire
         call allocate_block_data (grid, f_abm                 )
         call allocate_block_data (grid, f_gdp                 )
         call allocate_block_data (grid, f_peatf               )
         call allocate_block_data (grid, f_hdm                 )
         call allocate_block_data (grid, f_lnfm                )
#endif

      end if


   END SUBROUTINE allocate_2D_BGCFluxes

#endif
END MODULE MOD_2D_BGCFluxes
