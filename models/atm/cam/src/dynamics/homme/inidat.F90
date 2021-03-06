module inidat
  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: Read initial dataset and spectrally truncate as appropriate.
  !
  ! Method: Initialize one or a few fields at a time, to minimize the 
  !         memory  requirements
  ! 
  ! Author: 
  ! Modified: P. Worley, to implement initialization of subsets
  !           of fields. (8/03)
  ! 
  !-----------------------------------------------------------------------
  use cam_logfile, only : iulog
  use element_mod, only : element_t

  implicit none
  private
  public read_inidat

contains



  subroutine read_inidat( ncid_ini, ncid_topo, dyn_in)
    use dyn_comp,      only: dyn_import_t
    use hybrid_mod, only : hybrid_t
    use parallel_mod,     only: par
    use bndry_mod,     only: bndry_exchangev
    use constituents, only: cnst_name, cnst_read_iv, qmin
    use dimensions_mod,     only: nelemd, nlev, np
    use dof_mod, only           : putUniquePoints
    use edge_mod, only : edgevpack, edgevunpack, InitEdgeBuffer, FreeEdgeBuffer, EdgeBuffer_t
    use ncdio_atm, only : infld
    use shr_vmath_mod, only: shr_vmath_log
    use hycoef,           only: ps0
    use abortutils,     only: endrun
    use pio, only : file_desc_t, io_desc_t, pio_double, pio_get_local_array_size, pio_freedecomp
    use dyn_grid, only : get_horiz_grid_dim_d
    use shr_kind_mod, only: r8 => shr_kind_r8
    use chemistry   , only: chem_implements_cnst, chem_init_cnst
    use aerosol_intr, only: aerosol_implements_cnst, aerosol_init_cnst
    use tracers     , only: tracers_implements_cnst, tracers_init_cnst
    use aoa_tracers , only: aoa_tracers_implements_cnst, aoa_tracers_init_cnst
    use stratiform,   only: stratiform_implements_cnst, stratiform_init_cnst
    use co2_cycle   , only: co2_implements_cnst, co2_init_cnst

    implicit none
    type(file_desc_t),intent(inout) :: ncid_ini, ncid_topo
    type (dyn_import_t), target, intent(out)   :: dyn_in   ! dynamics import

    type(hybrid_t) :: hybrid
    type(element_t), pointer :: elem(:)
    real(r8), allocatable :: tmp(:,:)
    integer, allocatable :: gcols(:)
    integer :: tlncols, ig, ie, start, j, t, k
    character(len=40) :: fieldname
    logical :: found
    integer :: kptr, m_cnst
    type(EdgeBuffer_t) :: edge
    type(io_desc_t) :: iodesc
    integer :: lsize

    integer,parameter :: pcnst = PCNST
    integer, pointer :: gcid(:)


    elem=> dyn_in%elem
    
    call get_dyn_decomp(elem, nlev, pio_double, iodesc)

    lsize = pio_get_local_array_size(iodesc)	

    tlncols = lsize/nlev

    allocate(tmp(tlncols,nlev))    

    if(elem(1)%idxP%NumUniquePts <=0 .or. elem(1)%idxP%NumUniquePts > np*np) then
       write(iulog,*)  elem(1)%idxP%NumUniquePts
       call endrun('inidat')
    end if

    fieldname = 'U'
    call infld(fieldname, ncid_ini, iodesc, tlncols,'lev',  tmp, found)
    if(.not. found) then
       call endrun('Could not find U field on input datafile')
    end if

    start=1
    do ie=1,nelemd
       elem(ie)%state%v=0.0
       call putUniquePoints(elem(ie)%idxP, nlev, tmp(start:,:), &
            elem(ie)%state%v(:,:,1,:,1))
       start=start+elem(ie)%idxP%numUniquePts
    end do

    fieldname = 'V'
    call infld(fieldname, ncid_ini, iodesc, tlncols,'lev',  tmp, found)
    if(.not. found) then
       call endrun('Could not find V field on input datafile')
    end if
    start=1
    do ie=1,nelemd
       call putUniquePoints(elem(ie)%idxP, nlev, tmp(start:,:), &
            elem(ie)%state%v(:,:,2,:,1))
       start=start+elem(ie)%idxP%numUniquePts
    end do

    fieldname = 'T'
    call infld(fieldname, ncid_ini, iodesc, tlncols ,'lev',  tmp, found)
    if(.not. found) then
       call endrun('Could not find T field on input datafile')
    end if
    start=1
    do ie=1,nelemd
       elem(ie)%state%T=0.0
       call putUniquePoints(elem(ie)%idxP, nlev, tmp(start:,:), &
            elem(ie)%state%T(:,:,:,1))
       start=start+elem(ie)%idxP%numUniquePts
    end do

    gcid => get_ldof(elem, 1)

    ! qmin = 1e-12,0,0
    do m_cnst=1,pcnst
       found = .false.
       if(cnst_read_iv(m_cnst)) then
          call infld(cnst_name(m_cnst), ncid_ini, iodesc, tlncols, 'lev', tmp, found)
       end if
       if(.not. found) then

          if(par%masterproc  ) write(iulog,*) 'Field ',cnst_name(m_cnst),' not found on initial dataset'

          if (stratiform_implements_cnst(cnst_name(m_cnst))) then
             call stratiform_init_cnst(cnst_name(m_cnst), tmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "stratiform_init_cnst"'
          else if (chem_implements_cnst(cnst_name(m_cnst))) then
             call chem_init_cnst(cnst_name(m_cnst), tmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "chem_init_cnst"'
          else if (tracers_implements_cnst(cnst_name(m_cnst))) then
             call tracers_init_cnst(cnst_name(m_cnst), tmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "tracers_init_cnst"'
          else if (aoa_tracers_implements_cnst(cnst_name(m_cnst))) then
             call aoa_tracers_init_cnst(cnst_name(m_cnst), tmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "aoa_tracers_init_cnst"'
          else if (aerosol_implements_cnst(cnst_name(m_cnst))) then
             call aerosol_init_cnst(cnst_name(m_cnst), tmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "aerosol_init_cnst"'
          else if (co2_implements_cnst(cnst_name(m_cnst))) then
             call co2_init_cnst(cnst_name(m_cnst), tmp, gcid)
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), &
                   ' initialized by "co2_init_cnst"'
          else
              if(par%masterproc) write(iulog,*) '          ', cnst_name(m_cnst), ' set to 0.'
          end if
       end if
       do k=1,nlev
          do ig=1,tlncols
             tmp(ig,k)=max(qmin(m_cnst),tmp(ig,k))
          end do
       end do
       
       start=1
       do ie=1,nelemd
          elem(ie)%state%Q(:,:,:,m_cnst,:)=0.0
          call putUniquePoints(elem(ie)%idxP, nlev, tmp(start:,:), &
               elem(ie)%state%Q(:,:,:,m_cnst,1))
          start=start+elem(ie)%idxP%numUniquePts
       end do
    end do
    deallocate(gcid)

    call pio_freedecomp(ncid_ini, iodesc)

    call get_dyn_decomp(elem, 1, pio_double, iodesc)


    fieldname = 'PS'
    call infld(fieldname, ncid_ini, iodesc, tmp(:,1), found)
    if(.not. found) then
       call endrun('Could not find PS field on input datafile')
    end if
    start=1
    tmp(:,1)=tmp(:,1)*0.01_r8

    if(minval(tmp(:,1)) < 100) then
       call endrun('Problem reading ps field')
    end if

    do ie=1,nelemd
       elem(ie)%state%ps_v=0.0
       call putUniquePoints(elem(ie)%idxP, tmp(start:start+elem(ie)%idxP%numUniquePts-1,1), &
            elem(ie)%state%ps_v(:,:,1))
       start=start+elem(ie)%idxP%numUniquePts
    end do
    
    fieldname = 'PHIS'
    call infld(fieldname, ncid_topo, iodesc, tmp(:,1), found)
    if(.not. found) then
       call endrun('Could not find PHIS field on input datafile')
    end if
    start=1
    do ie=1,nelemd
       elem(ie)%state%phis=0.0
       call putUniquePoints(elem(ie)%idxP, tmp(start:,1), &
            elem(ie)%state%phis(:,:))
       start=start+elem(ie)%idxP%numUniquePts
    end do

    
    ! once we've read all the fields we do a boundary exchange to 
    ! update the redundent columns in the dynamics

    call initEdgeBuffer(edge, (3+pcnst)*nlev+2)

    do ie=1,nelemd
       kptr=0
       call edgeVpack(edge, elem(ie)%state%ps_v(:,:,1),1,kptr,elem(ie)%desc)
       kptr=kptr+1
       call edgeVpack(edge, elem(ie)%state%phis,1,kptr,elem(ie)%desc)
       kptr=kptr+1
       call edgeVpack(edge, elem(ie)%state%v(:,:,:,:,1),2*nlev,kptr,elem(ie)%desc)
       kptr=kptr+2*nlev
       call edgeVpack(edge, elem(ie)%state%T(:,:,:,1),nlev,kptr,elem(ie)%desc)
       kptr=kptr+nlev
       call edgeVpack(edge, elem(ie)%state%Q(:,:,:,:,1),nlev*pcnst,kptr,elem(ie)%desc)
    end do

    call bndry_exchangeV(par,edge)

    do ie=1,nelemd
       kptr=0
       call edgeVunpack(edge, elem(ie)%state%ps_v(:,:,1),1,kptr,elem(ie)%desc)
       kptr=kptr+1
       call edgeVunpack(edge, elem(ie)%state%phis,1,kptr,elem(ie)%desc)
       kptr=kptr+1
       call edgeVunpack(edge, elem(ie)%state%v(:,:,:,:,1),2*nlev,kptr,elem(ie)%desc)
       kptr=kptr+2*nlev
       call edgeVunpack(edge, elem(ie)%state%T(:,:,:,1),nlev,kptr,elem(ie)%desc)
       kptr=kptr+nlev
       call edgeVunpack(edge, elem(ie)%state%Q(:,:,:,:,1),nlev*pcnst,kptr,elem(ie)%desc)
    end do

!$omp parallel do private(ie, t, m_cnst)
    do ie=1,nelemd
       do t=2,3
          elem(ie)%state%ps_v(:,:,t)=elem(ie)%state%ps_v(:,:,1)
          elem(ie)%state%v(:,:,:,:,t)=elem(ie)%state%v(:,:,:,:,1)
          elem(ie)%state%T(:,:,:,t)=elem(ie)%state%T(:,:,:,1)
          do m_cnst=1,pcnst
             elem(ie)%state%Q(:,:,:,m_cnst,t)=elem(ie)%state%Q(:,:,:,m_cnst,1)
          end do
       end do
       call shr_vmath_log(elem(ie)%state%ps_v,elem(ie)%state%lnps,size(elem(ie)%state%lnps))
    end do

    call pio_freedecomp(ncid_ini, iodesc)
    call FreeEdgeBuffer(edge)
    deallocate(tmp)

  end subroutine read_inidat



  subroutine get_dyn_decomp(elem, nlev, datatype, iodesc)
    use pio, only : io_desc_t, pio_initdecomp
    use cam_pio_utils, only : pio_subsystem
    use dyn_grid, only : get_horiz_grid_dim_d

    type(element_t), pointer :: elem(:)
    integer, intent(in) :: nlev, datatype
    type(io_desc_t), intent(out) :: iodesc
    integer, pointer :: ldof(:)
    integer :: dimlens(2), dimcnt

    dimcnt=1
    call get_horiz_grid_dim_d(dimlens(1)) 
    if(nlev>1) then
       dimlens(2) = nlev
       dimcnt=dimcnt+1
    end if

    ldof => get_ldof(elem, nlev)
    call pio_initdecomp(pio_subsystem, datatype, dimlens(1:dimcnt), ldof, iodesc)
!    call pio_initdecomp(pio_subsystem, datatype, dimlens(1:dimcnt), ldof, iodesc,method=1)
    deallocate(ldof)

  end subroutine get_dyn_decomp


  function get_ldof(elem, nlev) result(ldof)
    use dimensions_mod,     only: nelemd
    use dyn_grid, only : get_horiz_grid_dim_d
    use abortutils,     only: endrun

    type(element_t), pointer :: elem(:)
    integer, intent(in) :: nlev
    integer, pointer :: ldof(:)

    integer :: lcnt, ie, j, k, ig, numpts, offset, hdim


    call get_horiz_grid_dim_d(hdim)

    lcnt = 0
    do ie=1,nelemd
       lcnt = lcnt+nlev*elem(ie)%idxP%NumUniquePts
    end do
    allocate(ldof(lcnt))
    ig=1
    ldof(:) = 0
    do k=1,nlev
       do ie=1,nelemd
          numpts = elem(ie)%idxP%NumUniquePts
          offset = elem(ie)%idxP%UniquePtOffset
          do j=1,numpts
             ldof(ig)=offset+(j-1)+(k-1)*hdim
             ig=ig+1
          end do
       end do
    end do


  end function get_ldof





end module inidat
