#include <define.h>

#ifdef CROP
MODULE MOD_CropReadin

!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: CROP_readin

CONTAINS

   SUBROUTINE CROP_readin ()
   ! ===========================================================
   ! ! DESCRIPTION:
   ! Read in crop planting date from data, and fertilization from data.
   ! Save these data in patch vector.
   !
   ! Original: Shupeng Zhang, Zhongwang Wei, and Xingjie Lu, 2022
   ! ===========================================================

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_SPMD_Task
   USE MOD_LandPatch
   USE MOD_NetCDFSerial
   USE MOD_NetCDFBlock
   USE MOD_Mapping_Grid2Pset
   USE MOD_Vars_TimeInvariants
   USE MOD_Vars_TimeVariables

   USE MOD_Vars_Global
   USE MOD_LandPFT
   USE MOD_Vars_PFTimeVariables
   USE MOD_RangeCheck
   USE MOD_Block

   IMPLICIT NONE

   character(len=256) :: file_crop
   type(grid_type)    :: grid_crop
   type(block_data_real8_2d)    :: f_xy_crop
   type(mapping_grid2pset_type) :: mg2patch_crop
   type(mapping_grid2pset_type) :: mg2pft_crop
   character(len=256) :: file_irrig
   type(grid_type)    :: grid_irrig
   type(block_data_int32_2d)    :: f_xy_irrig
   type(mapping_grid2pset_type) :: mg2pft_irrig
   
   character(len=256) :: file_irrigalloc
   type(grid_type) :: grid_irrigalloc
   type(block_data_real8_2d) :: f_xy_irrigalloc
   type(mapping_grid2pset_type) :: mg2p_irrigalloc

   real(r8),allocatable :: pdrice2_tmp      (:)
   real(r8),allocatable :: plantdate_tmp    (:)
   real(r8),allocatable :: fertnitro_tmp    (:)
   integer ,allocatable :: irrig_method_tmp (:)
   real(r8),allocatable :: vmax25_opti_tmp  (:)
   real(r8),allocatable :: slatop_opti_tmp  (:)
   real(r8),allocatable :: effcon_opti_tmp  (:)

   character(len=256)            :: file_parameters
   type(block_data_real8_2d)     :: f_xy_para
   type(grid_type)               :: grid_para
   type(mapping_grid2pset_type)  :: mg2pft_para
   
   ! Local variables
   real(r8), allocatable :: lat(:), lon(:)
   real(r8) :: missing_value
   integer  :: cft, npatch, ipft
   character(LEN=2) :: cx
   integer  :: iblkme, iblk, jblk
   integer  :: maxvalue, minvalue

      ! READ in crops

      file_crop = trim(DEF_dir_runtime) // '/crop/plantdt-colm-64cfts-rice2_fillcoast.nc'

      CALL ncio_read_bcast_serial (file_crop, 'lat', lat)
      CALL ncio_read_bcast_serial (file_crop, 'lon', lon)

      CALL grid_crop%define_by_center (lat, lon)

      IF (p_is_io) THEN
         CALL allocate_block_data  (grid_crop, f_xy_crop)
      ENDIF

      ! missing value
      IF (p_is_master) THEN
         CALL ncio_get_attr (file_crop, 'pdrice2', 'missing_value', missing_value)
      ENDIF
#ifdef USEMPI
      CALL mpi_bcast (missing_value, 1, MPI_REAL8, p_root, p_comm_glb, p_err)
#endif
      IF (p_is_io) THEN
         CALL ncio_read_block (file_crop, 'pdrice2', grid_crop, f_xy_crop)
      ENDIF

      CALL mg2patch_crop%build (grid_crop, landpatch, f_xy_crop, missing_value)
      CALL mg2pft_crop%build   (grid_crop, landpft,   f_xy_crop, missing_value)

      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)

      IF (p_is_worker) THEN
         IF (numpatch > 0)  allocate(pdrice2_tmp   (numpatch))
         IF (numpft   > 0)  allocate(plantdate_tmp (numpft))
         IF (numpft   > 0)  allocate(fertnitro_tmp (numpft))
         IF (numpft   > 0)  allocate(irrig_method_tmp (numpft))
         IF (numpft   > 0)  allocate(vmax25_opti_tmp (numpft))
         IF (numpft   > 0)  allocate(slatop_opti_tmp (numpft))
         IF (numpft   > 0)  allocate(effcon_opti_tmp (numpft))
      ENDIF

      ! (1) Read in plant date for rice2.
      file_crop = trim(DEF_dir_runtime) // '/crop/plantdt-colm-64cfts-rice2_fillcoast.nc'
      IF (p_is_io) THEN
         CALL ncio_read_block (file_crop, 'pdrice2', grid_crop, f_xy_crop)
      ENDIF

      CALL mg2patch_crop%map_aweighted (f_xy_crop, pdrice2_tmp)

      IF (p_is_worker) THEN
         DO npatch = 1, numpatch
            IF (pdrice2_tmp(npatch) /= spval) THEN
               pdrice2 (npatch) = int(pdrice2_tmp (npatch))
            ELSE
               pdrice2 (npatch) = 0
            ENDIF
         ENDDO
      ENDIF

#ifdef RangeCheck
      CALL check_vector_data ('plant date value for rice2 ', pdrice2)
#endif

      ! (2) Read in plant date.
      IF (p_is_worker) THEN
         plantdate_p(:) = -99999999._r8
      ENDIF

      file_crop = trim(DEF_dir_runtime) // '/crop/plantdt-colm-64cfts-rice2_fillcoast.nc'
      DO cft = 15, 78
         write(cx, '(i2.2)') cft
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_crop, 'PLANTDATE_CFT_'//trim(cx), grid_crop, 1, f_xy_crop)
         ENDIF

         CALL mg2pft_crop%map_aweighted (f_xy_crop, plantdate_tmp)

         IF (p_is_worker) THEN
            DO ipft = 1, numpft
               IF(landpft%settyp(ipft) .eq. cft)THEN
                  plantdate_p(ipft) = plantdate_tmp(ipft)
                  IF(plantdate_p(ipft) <= 0._r8) THEN
                     plantdate_p(ipft) = -99999999._r8
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO

#ifdef RangeCheck
      CALL check_vector_data ('plantdate_pfts value ', plantdate_p)
#endif

      IF (p_is_worker) THEN
         fertnitro_p(:) = -99999999._r8
      ENDIF

      file_crop = trim(DEF_dir_runtime) // '/crop/fertnitro_fillcoast.nc'
      DO cft = 15, 78
         write(cx, '(i2.2)') cft
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_crop, 'CONST_FERTNITRO_CFT_'//trim(cx), grid_crop, 1, f_xy_crop)
         ENDIF

         CALL mg2pft_crop%map_aweighted (f_xy_crop, fertnitro_tmp)

         IF (p_is_worker) THEN
            DO ipft = 1, numpft
               IF(landpft%settyp(ipft) .eq. cft)THEN
                  fertnitro_p(ipft) = fertnitro_tmp(ipft)
                  IF(fertnitro_p(ipft) <= 0._r8) THEN
                     fertnitro_p(ipft) = 0._r8
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO

#ifdef RangeCheck
      CALL check_vector_data ('fert nitro value ', fertnitro_p)
#endif

      ! (4) Read in irrigation method
!      file_irrig = trim(DEF_dir_runtime) // '/crop/surfdata_irrigation_method.nc'
      file_irrig = trim(DEF_dir_runtime) // '/crop/surfdata_irrigation_method_96x144.nc'

      CALL ncio_read_bcast_serial (file_irrig, 'lat', lat)
      CALL ncio_read_bcast_serial (file_irrig, 'lon', lon)

      CALL grid_irrig%define_by_center (lat, lon)

      IF (p_is_io) THEN
         CALL allocate_block_data  (grid_irrig, f_xy_irrig)
      ENDIF

      CALL mg2pft_irrig%build   (grid_irrig, landpft)

      IF (allocated(lon)) deallocate(lon)
      IF (allocated(lat)) deallocate(lat)

      IF (p_is_worker) THEN
         irrig_method_p(:) = -99999999
      ENDIF

      DO cft = 1, N_CFT
         IF (p_is_io) THEN
            CALL ncio_read_block_time (file_irrig, 'irrigation_method', grid_irrig, cft, f_xy_irrig)
         ENDIF

         CALL mg2pft_irrig%map_max_frequency_2d (f_xy_irrig, irrig_method_tmp)

         IF (p_is_worker) THEN
            DO ipft = 1, numpft

               IF(landpft%settyp(ipft) .eq. cft + 14)THEN
                  irrig_method_p(ipft) = irrig_method_tmp(ipft)
                  IF(irrig_method_p(ipft) < 0) THEN
                     irrig_method_p(ipft) = -99999999
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO

#ifdef RangeCheck
      CALL check_vector_data ('irrigation method ', irrig_method_p)
#endif

      IF (allocated (pdrice2_tmp  ))    deallocate(pdrice2_tmp  )
      IF (allocated (plantdate_tmp))    deallocate(plantdate_tmp)
      IF (allocated (fertnitro_tmp))    deallocate(fertnitro_tmp)
      IF (allocated (irrig_method_tmp)) deallocate(irrig_method_tmp)


      if (DEF_IRRIGATION_ALLOCATION == 3) then 
         ! (5) Read in irrigation allocated to groundwater
         file_irrigalloc = trim(DEF_dir_runtime) // '/crop/surfdata_irrigation_allocation.nc'

         CALL ncio_read_bcast_serial (file_irrigalloc, 'lat', lat)
         CALL ncio_read_bcast_serial (file_irrigalloc, 'lon', lon)

         CALL grid_irrigalloc%define_by_center (lat, lon)

         IF (p_is_io) THEN
            CALL allocate_block_data  (grid_irrigalloc, f_xy_irrigalloc)
         ENDIF

         call mg2p_irrigalloc%build   (grid_irrigalloc, landpatch)

         IF (allocated(lon)) deallocate(lon)
         IF (allocated(lat)) deallocate(lat)

         IF (p_is_worker) THEN
            irrig_gw_alloc(:) = 0._r8
         ENDIF

         IF (p_is_io) THEN
            CALL ncio_read_block (file_irrigalloc, 'irrig_gw_alloc', grid_irrigalloc, f_xy_irrigalloc)
         ENDIF

         call mg2p_irrigalloc%map_aweighted (f_xy_irrigalloc, irrig_gw_alloc)

#ifdef RangeCheck
         CALL check_vector_data ('irrigation goundwater allocation ', irrig_gw_alloc)
#endif

         ! (6) Read in irrigation allocated to surfacewater
         file_irrigalloc = trim(DEF_dir_runtime) // '/crop/surfdata_irrigation_allocation.nc'

         CALL ncio_read_bcast_serial (file_irrigalloc, 'lat', lat)
         CALL ncio_read_bcast_serial (file_irrigalloc, 'lon', lon)

         CALL grid_irrigalloc%define_by_center (lat, lon)

         IF (p_is_io) THEN
            CALL allocate_block_data  (grid_irrigalloc, f_xy_irrigalloc)
         ENDIF

         call mg2p_irrigalloc%build   (grid_irrigalloc, landpatch)

         IF (allocated(lon)) deallocate(lon)
         IF (allocated(lat)) deallocate(lat)

         IF (p_is_worker) THEN
            irrig_sw_alloc(:) = 0._r8
         ENDIF

         IF (p_is_io) THEN
            CALL ncio_read_block (file_irrigalloc, 'irrig_sw_alloc', grid_irrigalloc, f_xy_irrigalloc)
         ENDIF

         call mg2p_irrigalloc%map_aweighted (f_xy_irrigalloc, irrig_sw_alloc)

#ifdef RangeCheck
         CALL check_vector_data ('irrigation surfacewater allocation ', irrig_sw_alloc)
#endif
      endif

      if(DEF_optimization_crop)THEN
#ifndef SinglePoint
         file_parameters = trim(DEF_dir_runtime) // '/crop/optimization_parameters.nc'

         CALL ncio_read_bcast_serial (file_parameters, 'lat', lat)
         CALL ncio_read_bcast_serial (file_parameters, 'lon', lon)

         CALL grid_para%define_by_center (lat, lon)

         IF (p_is_io) THEN
            CALL allocate_block_data  (grid_para, f_xy_para)
         ENDIF

         call mg2pft_para%build   (grid_para, landpft)

         IF (allocated(lon)) deallocate(lon)
         IF (allocated(lat)) deallocate(lat)

         IF (p_is_worker) THEN
            vmax25_opti_p(:) = 0._r8
            slatop_opti_p(:) = 0._r8
            effcon_opti_p(:) = 0._r8
         ENDIF

         DO cft = 1, N_CFT
            IF (p_is_io) THEN
               CALL ncio_read_block_time (file_parameters, 'vmax25', grid_para, cft, f_xy_para)
            ENDIF
            call mg2pft_para%map_aweighted (f_xy_para, vmax25_opti_tmp)

            if(p_is_worker)then
               DO ipft = 1, numpft
                  IF(landpft%settyp(ipft) .eq. cft+14)THEN
                     vmax25_opti_p(ipft) = vmax25_opti_tmp(ipft) * 1.e-6
                     IF(vmax25_opti_p(ipft) < 0) THEN
                        vmax25_opti_p(ipft) = 0._r8
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

         DO cft = 1, N_CFT
            IF (p_is_io) THEN
               CALL ncio_read_block_time (file_parameters, 'slatop', grid_para, cft, f_xy_para)
            ENDIF
            call mg2pft_para%map_aweighted (f_xy_para, slatop_opti_tmp)

            if(p_is_worker)then
               DO ipft = 1, numpft
                  IF(landpft%settyp(ipft) .eq. cft+14)THEN
                     slatop_opti_p(ipft) = slatop_opti_tmp(ipft)
                     IF(slatop_opti_p(ipft) < 0) THEN
                        slatop_opti_p(ipft) = 0._r8
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

         DO cft = 1, N_CFT
            IF (p_is_io) THEN
               CALL ncio_read_block_time (file_parameters, 'effcon', grid_para, cft, f_xy_para)
            ENDIF
            call mg2pft_para%map_aweighted (f_xy_para, effcon_opti_tmp)

            if(p_is_worker)then
               DO ipft = 1, numpft
                  IF(landpft%settyp(ipft) .eq. cft+14)THEN
                     effcon_opti_p(ipft) = effcon_opti_tmp(ipft)
                     IF(effcon_opti_p(ipft) < 0) THEN
                        effcon_opti_p(ipft) = 0._r8
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

         IF (allocated (vmax25_opti_tmp  ))    deallocate(vmax25_opti_tmp  )
         IF (allocated (slatop_opti_tmp  ))    deallocate(slatop_opti_tmp  )
         IF (allocated (effcon_opti_tmp  ))    deallocate(effcon_opti_tmp  )
#else
         file_parameters = trim(DEF_dir_output) //  '/' // 'InputParameters.txt'
         open(10, file=file_parameters, status='old')
         read(10,*) vmax25_opti_p
         read(10,*) slatop_opti_p
         read(10,*) effcon_opti_p
         close(10)
         vmax25_opti_p(:) = vmax25_opti_p(:) * 1.e-6
#endif
      endif
   END SUBROUTINE CROP_readin

END MODULE MOD_CropReadin
#endif
