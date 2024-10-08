MODULE CMF_CALC_DIAG_MOD
!==========================================================
!* PURPOSE: Manage prognostic/diagnostic variables in CaMa-Flood
!
!* CONTAINS:
! -- CMF_PROG_INIT      : Initialize Prognostic variables (include restart data handling)
! -- CMF_DIAG_INIT      : Initialize Diagnostic variables
! -- CMF_DIAG_AVERAGE   : Calculate time-average of Diagnostic Variables
! -- CMF_DIAG_RESET     : Reset Diagnostic Variables (Average & Maximum )
!!
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
! Revised by Zhongwang Wei, add re-infiltation varialbe (WINFILT,D2WINFILTEX) 2023/04/31

! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
   USE PARKIND1,           only: JPIM, JPRM, JPRB
   USE YOS_CMF_INPUT,      only: LOGNAM
   USE YOS_CMF_INPUT,      only: DT, LPTHOUT,  LDAMOUT,  LWEVAP,LWINFILT
   USE YOS_CMF_MAP,        only: NSEQALL,      NPTHOUT,D2GRAREA
   USE YOS_CMF_PROG,       only: D2RIVOUT,     D2FLDOUT,     D1PTHFLW,     D2GDWRTN, &
                              & D2RUNOFF,     D2ROFSUB,     P2DAMINF
   USE YOS_CMF_DIAG,       only: D2OUTFLW,     D2RIVVEL,     D2PTHOUT,     D2PTHINF, &
                              & D2RIVDPH,     D2STORGE,     D2WEVAPEX,   D2WINFILTEX, NADD, &
                              & D2RIVOUT_AVG, D2FLDOUT_AVG, D1PTHFLW_AVG, D2GDWRTN_AVG, D2RUNOFF_AVG, D2ROFSUB_AVG, &
                              & D2OUTFLW_AVG, D2RIVVEL_AVG, D2PTHOUT_AVG, D2DAMINF_AVG, D2WEVAPEX_AVG, D2WINFILTEX_AVG, &
                              & D2OUTFLW_MAX, D2RIVDPH_MAX, D2STORGE_MAX
#ifdef sediment
   USE YOS_CMF_INPUT,      only: LSEDOUT
   USE yos_cmf_sed,        only: d2rivout_sed, d2rivvel_sed, sadd_riv
#endif
   IMPLICIT NONE
CONTAINS 
!####################################################################
! -- CMF_DIAG_AVE_MAX   : Add / Max of diagnostic variables at time step
! -- CMF_DIAG_AVERAGE   : Calculate time-average of Diagnostic Variables
! -- CMF_DIAG_RESET     : Reset Diagnostic Variables (Average & Maximum )
!
!####################################################################
   SUBROUTINE CMF_DIAG_AVEMAX

   IMPLICIT NONE
   integer(KIND=JPIM),SAVE  ::  ISEQ, IPTH
      !====================
      NADD=NADD+DT
      !$OMP PARALLEL DO
      DO ISEQ=1, NSEQALL
         D2RIVOUT_AVG(ISEQ,1)=D2RIVOUT_AVG(ISEQ,1)+D2RIVOUT(ISEQ,1)*DT
         D2FLDOUT_AVG(ISEQ,1)=D2FLDOUT_AVG(ISEQ,1)+D2FLDOUT(ISEQ,1)*DT
         D2RIVVEL_AVG(ISEQ,1)=D2RIVVEL_AVG(ISEQ,1)+D2RIVVEL(ISEQ,1)*DT
         D2OUTFLW_AVG(ISEQ,1)=D2OUTFLW_AVG(ISEQ,1)+D2OUTFLW(ISEQ,1)*DT

         D2PTHOUT_AVG(ISEQ,1)=D2PTHOUT_AVG(ISEQ,1)+D2PTHOUT(ISEQ,1)*DT-D2PTHINF(ISEQ,1)*DT

         D2GDWRTN_AVG(ISEQ,1)=D2GDWRTN_AVG(ISEQ,1)+D2GDWRTN(ISEQ,1)*DT
         D2RUNOFF_AVG(ISEQ,1)=D2RUNOFF_AVG(ISEQ,1)+D2RUNOFF(ISEQ,1)*DT
         D2ROFSUB_AVG(ISEQ,1)=D2ROFSUB_AVG(ISEQ,1)+D2ROFSUB(ISEQ,1)*DT

         D2OUTFLW_MAX(ISEQ,1)=max( D2OUTFLW_MAX(ISEQ,1), abs(D2OUTFLW(ISEQ,1)) )
         D2RIVDPH_MAX(ISEQ,1)=max( D2RIVDPH_MAX(ISEQ,1),     D2RIVDPH(ISEQ,1)  )
         D2STORGE_MAX(ISEQ,1)=max( D2STORGE_MAX(ISEQ,1),     D2STORGE(ISEQ,1)  )
         !recheck here zhongwang@Apr15,2024
         IF( LWEVAP )THEN
            D2WEVAPEX_AVG(ISEQ,1)= (D2WEVAPEX_AVG(ISEQ,1) +D2WEVAPEX(ISEQ,1)*DT) /D2GRAREA(ISEQ,1)
         ENDIF
         IF( LWINFILT )THEN
            D2WINFILTEX_AVG(ISEQ,1)= (D2WINFILTEX_AVG(ISEQ,1) +D2WINFILTEX(ISEQ,1)*DT)/D2GRAREA(ISEQ,1)
         ENDIF
      ENDDO
!$OMP END PARALLEL DO

!! loop for optional variable (separated for computational efficiency)
      IF( LDAMOUT )THEN
!$OMP PARALLEL DO
         DO ISEQ=1, NSEQALL
            D2DAMINF_AVG(ISEQ,1)=D2DAMINF_AVG(ISEQ,1)+P2DAMINF(ISEQ,1)*DT
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      IF( LPTHOUT )THEN
!$OMP PARALLEL DO
         DO IPTH=1, NPTHOUT
            D1PTHFLW_AVG(IPTH,:)=D1PTHFLW_AVG(IPTH,:)+D1PTHFLW(IPTH,:)*DT
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

#ifdef sediment
      !calculate average rivout and rivvel for sediment timestep
      IF( LSEDOUT )THEN
         sadd_riv = sadd_riv + DT
!$OMP PARALLEL DO
         DO ISEQ=1, NSEQALL
            d2rivout_sed(ISEQ) = d2rivout_sed(ISEQ)+D2RIVOUT(ISEQ,1)*DT
            d2rivvel_sed(ISEQ) = d2rivvel_sed(ISEQ)+D2RIVVEL(ISEQ,1)*DT
         ENDDO
!$OMP END PARALLEL DO
      ENDIF
#endif

   END SUBROUTINE CMF_DIAG_AVEMAX
!####################################################################





!####################################################################
   SUBROUTINE CMF_DIAG_AVERAGE
   USE YOS_CMF_TIME,       only: JYYYYMMDD, JHHMM
   IMPLICIT NONE
!================================================
      write(LOGNAM,*) "CMF::DIAG_AVERAGE: time-average", NADD, JYYYYMMDD, JHHMM
      D2RIVOUT_AVG(:,:) = D2RIVOUT_AVG(:,:) / DBLE(NADD)
      D2FLDOUT_AVG(:,:) = D2FLDOUT_AVG(:,:) / DBLE(NADD)
      D2OUTFLW_AVG(:,:) = D2OUTFLW_AVG(:,:) / DBLE(NADD)
      D2RIVVEL_AVG(:,:) = D2RIVVEL_AVG(:,:) / DBLE(NADD)
      D2PTHOUT_AVG(:,:) = D2PTHOUT_AVG(:,:) / DBLE(NADD)
      D2GDWRTN_AVG(:,:) = D2GDWRTN_AVG(:,:) / DBLE(NADD)
      D2RUNOFF_AVG(:,:) = D2RUNOFF_AVG(:,:) / DBLE(NADD)
      D2ROFSUB_AVG(:,:) = D2ROFSUB_AVG(:,:) / DBLE(NADD)
      IF ( LDAMOUT ) THEN
         D2DAMINF_AVG(:,:)  = D2DAMINF_AVG(:,:)  / DBLE(NADD)
      ENDIF
      IF ( LWEVAP ) THEN
         D2WEVAPEX_AVG(:,:) = D2WEVAPEX_AVG(:,:) / DBLE(NADD)
      ENDIF
      D1PTHFLW_AVG(:,:) = D1PTHFLW_AVG(:,:) /DBLE(NADD)
   END SUBROUTINE CMF_DIAG_AVERAGE
!####################################################################





!####################################################################
   SUBROUTINE CMF_DIAG_RESET
   USE YOS_CMF_TIME,       only: JYYYYMMDD, JHHMM
   IMPLICIT NONE
!================================================
      write(LOGNAM,*) "CMF::DIAG_AVERAGE: reset", JYYYYMMDD, JHHMM
      NADD=0
      D2RIVOUT_AVG(:,:) = 0._JPRB
      D2FLDOUT_AVG(:,:) = 0._JPRB
      D2OUTFLW_AVG(:,:) = 0._JPRB
      D2RIVVEL_AVG(:,:) = 0._JPRB
      D2PTHOUT_AVG(:,:) = 0._JPRB
      D2GDWRTN_AVG(:,:) = 0._JPRB
      D2RUNOFF_AVG(:,:) = 0._JPRB
      D2ROFSUB_AVG(:,:) = 0._JPRB
      IF ( LDAMOUT ) THEN
         D2DAMINF_AVG(:,:)  = 0._JPRB
      ENDIF
      IF ( LWEVAP ) THEN
         D2WEVAPEX_AVG(:,:) = 0._JPRB
      ENDIF
      IF ( LWINFILT ) THEN
         D2WINFILTEX_AVG(:,:) = 0._JPRB
      ENDIF      
      

      D1PTHFLW_AVG(:,:) = 0._JPRB 

      D2STORGE_MAX(:,:)=0._JPRB
      D2OUTFLW_MAX(:,:)=0._JPRB
      D2RIVDPH_MAX(:,:)=0._JPRB
   END SUBROUTINE CMF_DIAG_RESET
!####################################################################

END MODULE CMF_CALC_DIAG_MOD
