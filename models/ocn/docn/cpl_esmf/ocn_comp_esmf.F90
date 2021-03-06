!===============================================================================
! SVN: $Id: ocn_comp_esmf.F90 9228 2008-03-14 00:24:55Z tcraig $ 
! SVN: $URL: https://svn-ccsm-models.cgd.ucar.edu/docn7/branch_tags/drva_docn7_070824_tags/drva13_docn7_071129/ocn_comp_esmf.F90 $
!===============================================================================

module ocn_comp_esmf

   use shr_kind_mod, only:  R8=>SHR_KIND_R8, IN=>SHR_KIND_IN, &
                            CS=>SHR_KIND_CS, CL=>SHR_KIND_CL
   use shr_sys_mod   ! shared system calls

   use seq_comm_mct, only: OCNID
   use seq_cdata_mod
   use seq_infodata_mod

   use esmfshr_mod

   use docn_comp_mod
   use esmf_mod
   use perf_mod
   use mct_mod

   implicit none

   public :: ocn_init_esmf
   public :: ocn_run_esmf
   public :: ocn_final_esmf
   public :: ocn_register_esmf

   private ! except

   type(seq_cdata)     :: cdata
   type(seq_infodata_type)  :: infodata
   type(mct_gsMap)     :: gsmap
   type(mct_gGrid)     :: ggrid
   type(mct_aVect)     :: x2d
   type(mct_aVect)     :: d2x

   !----- formats -----
   character(*),parameter :: subName =  "(ocn_comp_esmf) "

   save ! save everything

!
! Author: Fei Liu
! This module is ESMF compliant ocn data component 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================

subroutine ocn_register_esmf(comp, rc)

   implicit none

   type(ESMF_GridComp)  :: comp
   integer, intent(out) :: rc

   rc = ESMF_SUCCESS

   print *, "In ocn register routine"
   ! Register the callback routines.

   call ESMF_GridCompSetEntryPoint(comp, ESMF_SETINIT, ocn_init_esmf, ESMF_SINGLEPHASE, rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   call ESMF_GridCompSetEntryPoint(comp, ESMF_SETRUN, ocn_run_esmf, ESMF_SINGLEPHASE, rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   call ESMF_GridCompSetEntryPoint(comp, ESMF_SETFINAL, ocn_final_esmf, ESMF_SINGLEPHASE, rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

end subroutine

!===============================================================================

subroutine ocn_init_esmf(comp, import_state, export_state, EClock, rc)
   !----------------------------------------------------------
   
   implicit none

   !----- arguments -----
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc
   
   !----- local -----
   character(CL)    :: NLFilename
   type(ESMF_Array) :: Ex2d, Ed2x, Edom
   integer(IN)      :: phase
   
   character(*),parameter :: subName = "(ocn_init_esmf) "
   !----------------------------------------------------------
   
   rc = ESMF_SUCCESS

   NLFilename = 'unused'

   call esmfshr_infodata_state2infodata(export_state,infodata)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   call seq_infodata_GetData(infodata,ocn_phase=phase)

   if (phase == 1) then
      call seq_cdata_init(cdata,OCNID,ggrid,gsmap,infodata,'docn')
   else
      call ESMF_StateGet(import_state, itemName="x2d", array=Ex2d, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call esmf2mct_copy(Ex2d, x2d, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   endif

   call docn_comp_init(EClock, cdata, x2d, d2x, NLFilename)

   call esmfshr_infodata_infodata2state(infodata,export_state,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   if (phase == 1) then
      Edom = mct2esmf_init(ggrid%data,gsmap,name='domain',rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call mct2esmf_copy(ggrid%data,Edom,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      Ed2x = mct2esmf_init(d2x,gsmap,name='d2x',rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      Ex2d = mct2esmf_init(x2d,gsmap,name='x2d',rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call mct2esmf_copy(d2x,Ed2x,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

      call ESMF_StateAdd(export_state,Edom,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call ESMF_StateAdd(export_state,Ed2x,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call ESMF_StateAdd(import_state,Ex2d,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   else
      call ESMF_StateGet(export_state, itemName="d2x", array=Ed2x, rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
      call mct2esmf_copy(d2x,Ed2x,rc=rc)
      if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   endif

   rc = ESMF_SUCCESS

end subroutine ocn_init_esmf

!===============================================================================

subroutine ocn_run_esmf(comp, import_state, export_state, EClock, rc)

   implicit none

   !----- arguments -----
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc

   !----- local -----
   type(ESMF_Array)    :: Ex2d, Ed2x
   
   character(*),parameter :: subName = "(ocn_run_esmf) "
   !----------------------------------------------------------
   
   rc = ESMF_SUCCESS

   call esmfshr_infodata_state2infodata(export_state,infodata)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   call ESMF_StateGet(import_state, itemName="x2d", array=Ex2d, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   call esmf2mct_copy(Ex2d, x2d, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   call docn_comp_run(EClock, cdata, x2d, d2x)

   call esmfshr_infodata_infodata2state(infodata,export_state,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   call ESMF_StateGet(export_state, itemName="d2x", array=Ed2x, rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)
   call mct2esmf_copy(d2x,Ed2x,rc=rc)
   if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, terminationflag=ESMF_ABORT)

   rc = ESMF_SUCCESS

end subroutine ocn_run_esmf

!===============================================================================

subroutine ocn_final_esmf(comp, import_state, export_state, EClock, rc)

   implicit none

   !----- arguments -----
   type(ESMF_GridComp)          :: comp
   type(ESMF_State)             :: import_state
   type(ESMF_State)             :: export_state
   type(ESMF_Clock)             :: EClock
   integer, intent(out)         :: rc

   !----------------------------------------------------------------------------
   ! Finalize routine 
   !----------------------------------------------------------------------------
  
   rc = ESMF_SUCCESS

   call docn_comp_final()
  
end subroutine ocn_final_esmf

!===============================================================================

end module ocn_comp_esmf
