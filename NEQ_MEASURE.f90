MODULE NEQ_MEASURE
  USE NEQ_INPUT_VARS
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE ELECTRIC_FIELD
  USE SF_CONSTANTS
  USE SF_MISC, only: assert_shape
  USE SF_IOTOOLS, only:free_unit,str,reg,splot,splot3d
  USE SF_TIMER, only: start_timer,stop_timer,eta
  implicit none

  private

  interface measure_observables
     module procedure :: measure_observables_rank0
  end interface measure_observables
  public :: measure_observables

  interface measure_current
     module procedure :: measure_current_rank0
  end interface measure_current
  public :: measure_current


  !Other stand alone
  public :: measure_dens
  public :: measure_docc
  public :: measure_ekin
  public :: measure_epot
  public :: measure_etot


contains



  !+-------------------------------------------------------------------+
  !PURPOSE: measure some observables and print them
  !+-------------------------------------------------------------------+
  subroutine measure_observables_rank0(g,self)
    type(kb_contour_gf(*,*,*)),intent(in)        :: G
    type(kb_contour_gf(G%N,G%L,G%Lf)),intent(in) :: Self
    integer                                      :: unit,itime
    real(8)                                      :: dens,docc,ekin,epot,etot
    itime = cc_params%Nt
    unit = free_unit()
    open(unit,file="observables.info")
    write(unit,"(8A20)")"time","N","Docc","Ekin","Epot","Etot"
    close(unit)
    open(unit,file="observables.neqipt",position="append")
    dens = measure_dens(g,self)
    docc = measure_docc(g,self)
    ekin = measure_ekin(g,self)
    epot = measure_epot(g,self)
    etot = ekin + epot
    write(unit,"(6F20.12)")cc_params%t(itime),dens,docc,ekin,epot,etot
    close(unit)
  end subroutine measure_observables_rank0

  ! subroutine measure_observables_rank2(params,g,self)
  !   type(kb_contour_gf)     :: g(:,:)
  !   type(kb_contour_gf)     :: self(:,:)
  !   type(kb_contour_params) :: params
  !   integer                 :: unit,itime,io,Nso
  !   real(8)                 :: dens,docc,ekin,epot,etot
  !   Nso = size(G,1)
  !   call assert_shape_kb_contour_gf(G,[Nso,Nso],"measure_observables","G")
  !   call assert_shape_kb_contour_gf(Self,[Nso,Nso],"measure_observables","Self")
  !   itime = params%Nt
  !   unit = free_unit()
  !   open(unit,file="observables.info")
  !   write(unit,"(8A20)")"time","N","Docc","Ekin","Epot","Etot"
  !   close(unit)
  !   do io=1,Nso
  !      open(unit,file="observables_io"//str(io,4)//".neqipt",position="append")
  !      dens = measure_dens(params,g(io,io),self(io,io))
  !      docc = measure_docc(params,g(io,io),self(io,io))
  !      ekin = measure_ekin(params,g(io,io),self(io,io))
  !      epot = measure_epot(params,g(io,io),self(io,io))
  !      etot = ekin + epot
  !      write(unit,"(6F20.12)")params%t(itime),dens,docc,ekin,epot,etot
  !      close(unit)
  !   enddo
  ! end subroutine measure_observables_rank2





  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################





  !+-------------------------------------------------------------------+
  !PURPOSE: measure current
  !+-------------------------------------------------------------------+
  subroutine measure_current_rank0(Gk,Vkt)
    type(kb_contour_gf(*,*,*)),intent(in) :: Gk(:)
    real(8),dimension(:,:,:)              :: Vkt ![Nk][Ntime][Dim]
    integer                               :: unit,itime,Lk,ik,i
    real(8),dimension(size(Vkt,3))        :: Jloc
    real(8)                               :: nkt
    !
    Lk=size(gk)
    itime = cc_params%Nt
    !
    if(size(Vkt,1)/=Lk)stop "neq_measure_current: dim(Vkt,1) != Lk"
    if(size(Vkt,2)<cc_params%Ntime)stop "neq_measure_current: dim(Vkt,2) < Ntime"
    !
    unit = free_unit()
    open(unit,file="current.info")
    write(unit,"(8A20)")"time","Jx","Jy","Jz"
    close(unit)
    !
    open(unit,file="current.neqipt",position="append")
    Jloc=0d0
    do ik=1,Lk
       nkt  = dimag(Gk(ik)%less(itime,itime))
       Jloc = Jloc + Nkt*Vkt(ik,itime,:)/Lk
    enddo
    write(unit,"(4F20.12)")cc_params%t(itime),(Jloc(i),i=1,size(Jloc))
    close(unit)
  end subroutine measure_current_rank0

  ! subroutine measure_current_rank2(params,Gk,Vkt)
  !   type(kb_contour_gf),dimension(:,:,:) :: gk  ![Nso][Nso][Nk]
  !   real(8),dimension(:,:,:,:,:)         :: Vkt ![Nso][Nso][Nk][Ntime][Dim]
  !   type(kb_contour_params)              :: params
  !   integer                              :: unit,itime,Lk,ik,i,Nso,io,dim
  !   real(8),dimension(size(Vkt,5))       :: Jloc
  !   real(8)                              :: nkt
  !   !
  !   Nso = size(Gk,1)
  !   Lk  = size(Gk,3)
  !   dim = size(Vkt,5)
  !   itime = params%Nt
  !   !
  !   call assert_shape_kb_contour_gf(Gk,[Nso,Nso,Lk],"measure_current_rank2","Gk")
  !   call assert_shape(Vkt,[Nso,Nso,Lk,params%Ntime,dim],"measure_current_rank2","Vk(t)")
  !   !
  !   unit = free_unit()
  !   open(unit,file="current.info")
  !   write(unit,"(8A20)")"time","Jx","Jy","Jz"
  !   close(unit)
  !   !
  !   do io=1,Nso
  !      open(unit,file="current_io"//str(io,4)//".neqipt",position="append")
  !      Jloc=0d0
  !      do ik=1,Lk
  !         nkt  = dimag(Gk(io,io,ik)%less(itime,itime))
  !         Jloc = Jloc + Nkt*Vkt(io,io,ik,itime,:)/Lk
  !      enddo
  !      write(unit,"(4F20.12)")params%t(itime),(Jloc(i),i=1,size(Jloc))
  !      close(unit)
  !   enddo
  ! end subroutine measure_current_rank2








  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################




  ! subroutine get_wigner_Aw_rank0(params,G,Aw,it)
  !   type(kb_contour_params)        :: params
  !   type(kb_contour_gf)            :: G
  !   real(8),dimension(:)           :: Aw
  !   integer                        :: it
  !   integer                        :: Nstep,L
  !   complex(8),dimension(size(Aw)) :: Gret
  !   !
  !   Nstep = params%Ntime
  !   L     = params%Nwr
  !   if(it<1 .OR. it>Nstep)stop "ERROR get_wigner_Aw: it !\in [1,Nstep]"
  !   call assert_shape(Aw,[L],"get_wigner_Aw","Aw")
  !   call wigner_transform_time(params,it,G%ret,Gret) ; Gret =Gret+xi*params%dt/2
  !   Aw = -dimag(Gret)/pi
  ! end subroutine get_wigner_Aw_rank0

  ! subroutine get_wigner_Aw_rank2(params,G,Aw,it)
  !   type(kb_contour_params)          :: params
  !   type(kb_contour_gf)              :: G(:,:)
  !   real(8),dimension(:,:,:)         :: Aw
  !   integer                          :: it
  !   integer                          :: Nstep,L
  !   integer                          :: io,jo,Nso
  !   complex(8),dimension(size(Aw,3)) :: Gret
  !   !
  !   Nstep = params%Ntime
  !   L     = params%Nwr
  !   Nso   = size(G,1)
  !   if(it<1 .OR. it>Nstep)stop "ERROR get_wigner_Aw: it !\in [1,Nstep]"
  !   call assert_shape_kb_contour_gf(G,[Nso,Nso],"get_wigner_Aw","G")
  !   call assert_shape(Aw,[Nso,Nso,L],"get_wigner_Aw","Aw")
  !   do io=1,Nso
  !      do jo=1,Nso
  !         call wigner_transform_time(params,it,G(io,jo)%ret,Gret) ; Gret =Gret+xi*params%dt/2
  !         Aw(io,jo,:) = -dimag(Gret(:))/pi
  !      enddo
  !   enddo
  ! end subroutine get_wigner_Aw_rank2


  ! subroutine get_wigner_Nw_rank0(params,G,Nw,it)
  !   type(kb_contour_params)        :: params
  !   type(kb_contour_gf)            :: G
  !   real(8),dimension(:)           :: Nw
  !   integer                        :: it
  !   integer                        :: Nstep,L
  !   complex(8),dimension(size(Nw)) :: Gless
  !   !
  !   Nstep = params%Ntime
  !   L     = params%Nwr
  !   if(it<1 .OR. it>Nstep)stop "ERROR get_wigner_Nw: it !\in [1,Nstep]"
  !   call assert_shape(Nw,[L],"get_wigner_Nw","Nw")
  !   call wigner_transform_time(params,it,G%less,Gless) ; Gless =Gless+xi*params%dt/2
  !   Nw = dimag(Gless)/pi
  ! end subroutine get_wigner_Nw_rank0

  ! subroutine get_wigner_Nw_rank2(params,G,Nw,it)
  !   type(kb_contour_params)          :: params
  !   type(kb_contour_gf)              :: G(:,:)
  !   real(8),dimension(:,:,:)         :: Nw
  !   integer                          :: it
  !   integer                          :: Nstep,L
  !   integer                          :: io,jo,Nso
  !   complex(8),dimension(size(Nw,3)) :: Gless
  !   !
  !   Nstep = params%Ntime
  !   L     = params%Nwr
  !   Nso   = size(G,1)
  !   if(it<1 .OR. it>Nstep)stop "ERROR get_wigner_Nw: it !\in [1,Nstep]"
  !   call assert_shape_kb_contour_gf(G,[Nso,Nso],"get_wigner_Nw","G")
  !   call assert_shape(Nw,[Nso,Nso,L],"get_wigner_Nw","Nw")
  !   do io=1,Nso
  !      do jo=1,Nso
  !         call wigner_transform_time(params,it,G(io,jo)%less,Gless) ; Gless =Gless+xi*params%dt/2
  !         Nw(io,jo,:) = dimag(Gless(:))/pi
  !      enddo
  !   enddo
  ! end subroutine get_wigner_Nw_rank2






  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################




  ! !Evaluate and plot Wigner functions in hybrid representation F(w,t)
  ! ! Spectral function
  ! ! Distribution and Occupation
  ! subroutine get_wigner_functions_main(params,G,suffix,file)
  !   type(kb_contour_params)               :: params
  !   type(kb_contour_gf)                   :: G
  !   character(len=*)                      :: suffix,file
  !   integer                               :: Nstep,L
  !   integer                               :: it
  !   complex(8),dimension(:,:),allocatable :: Gret,Gless
  !   real(8),dimension(:),allocatable      :: Aw,Nw,Fw
  !   !
  !   Nstep = params%Ntime
  !   L     = params%Nwr
  !   !
  !   allocate(Gret(L,Nstep),Gless(L,Nstep),Aw(L),Nw(L),Fw(L))
  !   !
  !   call wigner_transform_fft(params,G%ret,Gret) ;Gret =Gret+xi*params%dt/2
  !   call wigner_transform_fft(params,G%less,Gless);Gless=Gless-xi*params%dt/4
  !   !
  !   !ImG^< = ImG^R*F_fermi(w)
  !   do it=1,Nstep
  !      Aw = -dimag(Gret(:,it))/pi
  !      Nw = dimag(Gless(:,it))/pi
  !      Fw = Nw/Aw
  !      where(abs(params%wr)>params%wmax)Fw=0d0
  !      call splot("WIGNER_"//reg(file)//"/Aw"//reg(suffix)//"_tstep"//str(it,4)//".neqipt",params%wr,Aw)
  !      call splot("WIGNER_"//reg(file)//"/Nw"//reg(suffix)//"_tstep"//str(it,4)//".neqipt",params%wr,Nw)
  !      call splot("WIGNER_"//reg(file)//"/Fw"//reg(suffix)//"_tstep"//str(it,4)//".neqipt",params%wr,Fw)
  !   enddo
  !   deallocate(Aw,Fw,Nw,Gret,Gless)
  ! end subroutine get_wigner_functions_main


  ! subroutine get_wigner_functions_rank0(params,G,file)
  !   type(kb_contour_params) :: params
  !   type(kb_contour_gf)     :: G
  !   character(len=*)        :: file
  !   write(*,"(A)")"Perform Wigner rotation"
  !   call system("mkdir -p WIGNER_"//reg(file))
  !   call get_wigner_functions_main(params,G,"",file)
  ! end subroutine get_wigner_functions_rank0

  ! subroutine get_wigner_functions_rank2(params,G,file,iprint)
  !   type(kb_contour_params) :: params
  !   type(kb_contour_gf)     :: G(:,:)
  !   character(len=*)        :: file
  !   integer,optional        :: iprint
  !   integer                 :: iprint_
  !   integer                 :: Nso,io,jo
  !   iprint_=1;if(present(iprint))iprint_=iprint
  !   Nso = size(G,1)
  !   call assert_shape_kb_contour_gf(G,[Nso,Nso],"get_wigner_functions","G")
  !   write(*,"(A)")"Perform Wigner rotation"
  !   call system("mkdir -p WIGNER_"//reg(file))
  !   select case (iprint_)
  !   case (1)
  !      do io=1,Nso
  !         call get_wigner_functions_main(params,G(io,io),"_io"//str(io),file)
  !      enddo
  !   case default
  !      do io=1,Nso
  !         do jo=1,Nso
  !            call get_wigner_functions_main(params,G(io,io),"_io"//str(io)//"_jo"//str(jo),file)
  !         enddo
  !      enddo
  !   end select
  ! end subroutine get_wigner_functions_rank2






  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################





  ! subroutine get_wigner_Aw_along_Kpath_rank0(params,hkt_model,Self,kpath,Nkpath,Nt)
  !   interface 
  !      function hkt_model(kpoint,time)
  !        real(8),dimension(:)      :: kpoint
  !        real(8)                   :: time
  !        complex(8)                :: hkt_model
  !      end function hkt_model
  !   end interface
  !   type(kb_contour_params)               :: params
  !   type(kb_contour_gf)                   :: Self
  !   complex(8),dimension(:,:),allocatable :: Hkt      ![Nk][Nt]
  !   real(8),dimension(:,:),allocatable    :: kpath
  !   integer                               :: Nt,Nstep,L,Lk,Npts,Ndim,Nkpath
  !   integer                               :: ipts,ik,ic,j,it,itime
  !   type(kb_contour_gf),allocatable       :: Gk(:)
  !   type(kb_contour_dgf),allocatable      :: dGk(:)
  !   type(kb_contour_dgf),allocatable      :: dGk_old(:)
  !   complex(8),dimension(:,:),allocatable :: Gkret
  !   real(8),dimension(:,:),allocatable    :: Akw
  !   real(8),allocatable                   :: kseg(:)
  !   real(8),dimension(size(kpath,2))      :: kstart,kstop,kpoint,kdiff
  !   real(8)                               :: klen
  !   !
  !   Nstep = params%Ntime
  !   L     = params%Nwr
  !   Npts  = size(kpath,1)
  !   Ndim  = size(kpath,2)
  !   Lk    = (Npts-1)*Nkpath
  !   call check_kb_contour_gf(params,Self,"get_wigner_Aw_along_Kpath_main","Self")
  !   !
  !   write(*,"(A)")"Solve Wigner function along the path:"
  !   do ipts=1,Npts
  !      write(*,*)"Point"//str(ipts)//": [",(kpath(ipts,ic),ic=1,size(kpath,2)),"]"
  !   enddo
  !   allocate(kseg(Lk))
  !   allocate(Hkt(Lk,Nstep))
  !   kseg=0d0
  !   ic  =0
  !   do ipts=1,Npts-1
  !      kstart = kpath(ipts,:)
  !      kstop  = kpath(ipts+1,:)
  !      kdiff  = (kstop-kstart)/dble(Nkpath)
  !      do ik=1,Nkpath
  !         ic=ic+1
  !         kseg(ic) = klen
  !         klen     = klen + sqrt(dot_product(kdiff,kdiff))
  !         kpoint   = kstart + (ik-1)*kdiff
  !         do itime=1,Nstep
  !            hkt(ic,itime) = hkt_model(kpoint,params%t(itime))
  !         enddo
  !      enddo
  !   enddo
  !   !
  !   !
  !   allocate(Gk(Lk),dGk(Lk),dGk_old(Lk))
  !   call allocate_kb_contour_gf(params,Gk)
  !   call allocate_kb_contour_gf(params,dGk)
  !   call allocate_kb_contour_gf(params,dGk_old)
  !   write(*,"(A)")"Evolve Gk(t,t') along the path:"
  !   do itime=1,Nstep
  !      params%Nt=itime
  !      dGk_old(:) = dGk(:)
  !      call vide_kb_contour_gf(params,Hkt,Self,Gk,dGk_old,dGk)
  !   enddo
  !   !
  !   call system("mkdir -p WIGNER_KPATH")
  !   !
  !   allocate(Gkret(Lk,L))
  !   allocate(Akw(Lk,L))
  !   call start_timer()
  !   do j=0,Nt
  !      it = Nstep-j
  !      do ik=1,Lk
  !         call wigner_transform_time(params,it,Gk(ik)%ret,Gkret(ik,:)) ; Gkret(ik,:) =Gkret(ik,:)+xi*params%dt/2
  !         Akw(ik,:) = -dimag(Gkret(ik,:))/pi
  !      enddo
  !      call splot3d("WIGNER_KPATH/Akw_ntime"//str(it,4)//".neqipt",kseg,params%wr,Akw)
  !      call eta(j,Nt)
  !   enddo
  !   call stop_timer()
  !   !
  !   call deallocate_kb_contour_gf(Gk)
  !   call deallocate_kb_contour_gf(dGk)
  !   call deallocate_kb_contour_gf(dGk_old)
  !   deallocate(Gk,dGk,dGk_old,Gkret,Akw,Hkt)
  ! end subroutine get_wigner_Aw_along_Kpath_rank0


  ! subroutine get_wigner_Aw_along_Kpath_rank2(params,hkt_model,Self,kpath,Nkpath,Nt)
  !   interface 
  !      function hkt_model(kpoint,time,N)
  !        real(8),dimension(:)      :: kpoint
  !        real(8)                   :: time
  !        integer                   :: N
  !        complex(8),dimension(N,N) :: hkt_model
  !      end function hkt_model
  !   end interface
  !   type(kb_contour_params)                   :: params
  !   type(kb_contour_gf)                       :: Self(:,:) ![Nso,Nso]
  !   complex(8),dimension(:,:,:,:),allocatable :: Hkt      ![Nso,Nso,Nk,Nt]
  !   real(8),dimension(:,:),allocatable        :: kpath
  !   integer                                   :: Nt,Nstep,L,Lk,Npts,Ndim,Nkpath,Nso
  !   integer                                   :: ipts,ik,ic,j,it,itime,io
  !   type(kb_contour_gf),allocatable           :: Gk(:,:,:) ![Nso,Nso,Lk]
  !   type(kb_contour_dgf),allocatable          :: dGk(:,:,:)
  !   type(kb_contour_dgf),allocatable          :: dGk_old(:,:,:)
  !   complex(8),dimension(:,:,:),allocatable   :: Gkret
  !   real(8),dimension(:,:,:),allocatable      :: Akw
  !   real(8),allocatable                       :: kseg(:)
  !   real(8),dimension(size(kpath,2))          :: kstart,kstop,kpoint,kdiff
  !   real(8)                                   :: klen
  !   !
  !   Nstep = params%Ntime
  !   L     = params%Nwr
  !   Nso   = size(Self,1)
  !   Lk    = size(Hkt,3)
  !   Npts  = size(kpath,1)
  !   Ndim  = size(kpath,2)
  !   Lk    = (Npts-1)*Nkpath
  !   call check_kb_contour_gf(params,Self,"get_wigner_Aw_along_Kpath_rank2","Self")
  !   call assert_shape_kb_contour_gf(Self,[Nso,Nso],"get_wigner_Aw_along_Kpath_rank2","Self")
  !   !
  !   write(*,"(A)")"Solve Wigner function along the path:"
  !   do ipts=1,Npts
  !      write(*,*)"Point"//str(ipts)//": [",(kpath(ipts,ic),ic=1,size(kpath,2)),"]"
  !   enddo
  !   allocate(kseg(Lk))
  !   allocate(Hkt(Nso,Nso,Lk,Nstep))
  !   kseg=0d0
  !   ic  =0
  !   do ipts=1,Npts-1
  !      kstart = kpath(ipts,:)
  !      kstop  = kpath(ipts+1,:)
  !      kdiff  = (kstop-kstart)/dble(Nkpath)
  !      do ik=1,Nkpath
  !         ic=ic+1
  !         kseg(ic) = klen
  !         klen     = klen + sqrt(dot_product(kdiff,kdiff))
  !         kpoint   = kstart + (ik-1)*kdiff
  !         do itime=1,Nstep
  !            hkt(:,:,ic,itime) = hkt_model(kpoint,params%t(itime),Nso)
  !         enddo
  !      enddo
  !   enddo
  !   !
  !   !
  !   allocate(Gk(Nso,Nso,Lk),dGk(Nso,Nso,Lk),dGk_old(Nso,Nso,Lk))
  !   call allocate_kb_contour_gf(params,Gk)
  !   call allocate_kb_contour_gf(params,dGk)
  !   call allocate_kb_contour_gf(params,dGk_old)
  !   write(*,"(A)")"Evolve Gk(t,t') along the path:"
  !   do itime=1,params%Ntime
  !      params%Nt=itime
  !      dGk_old = dGk
  !      call vide_kb_contour_gf(params,Hkt,Self,Gk,dGk_old,dGk)
  !   enddo
  !   !
  !   call system("mkdir -p WIGNER_KPATH")
  !   !
  !   allocate(Gkret(Nso,Lk,L))
  !   allocate(Akw(Nso,Lk,L))
  !   call start_timer()
  !   do j=0,Nt
  !      it = Nstep-j
  !      do ik=1,Lk
  !         do io=1,Nso
  !            call wigner_transform_time(params,it,Gk(io,io,ik)%ret,Gkret(io,ik,:)) ; Gkret(io,ik,:) =Gkret(io,ik,:)+xi*params%dt/2
  !            Akw(io,ik,:) = -dimag(Gkret(io,ik,:))/pi/Nso
  !         enddo
  !      enddo
  !      call splot3d("WIGNER_KPATH/Akw_ntime"//str(it,4)//".neqipt",kseg,params%wr,sum(Akw,dim=1))
  !      call eta(j,Nt)
  !   enddo
  !   call stop_timer()
  !   !
  !   call deallocate_kb_contour_gf(Gk)
  !   call deallocate_kb_contour_gf(dGk)
  !   call deallocate_kb_contour_gf(dGk_old)
  !   deallocate(Gk,dGk,dGk_old,Gkret,Akw,Hkt)
  ! end subroutine get_wigner_Aw_along_Kpath_rank2






  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################




  ! !Performs the Wigner transformation of a Keldysh component:
  ! !G^X(w,t) = int_0^T ds G^X(t,t-s)*exp(iws)         
  ! subroutine wigner_transform_fft(params,Gin,Gf)
  !   type(kb_contour_params)   :: params
  !   integer                   :: M,it,iw,is
  !   complex(8),dimension(:,:) :: Gin
  !   complex(8),dimension(:,:) :: Gf
  !   Gf=zero
  !   do it=1,params%Ntime
  !      do iw=1,params%Nwr
  !         Gf(iw,it) = sum(Gin(it,it:1:-1)*exp(xi*params%wr(iw)*params%t(1:it)))
  !      enddo
  !   enddo
  !   Gf=Gf*params%dt
  ! end subroutine wigner_transform_fft





  ! subroutine wigner_transform_time(params,it,Gin,Gf)
  !   type(kb_contour_params)   :: params
  !   integer                   :: it
  !   complex(8),dimension(:,:) :: Gin
  !   complex(8),dimension(:)   :: Gf
  !   integer                   :: M,iw,is
  !   Gf=zero
  !   if(size(Gin,1)/=size(Gin,2))stop "wigner error1"
  !   if(it<1.OR.it>size(Gin,1))stop "Wigner error2"
  !   !G(w,t) = int_0^T ds G(t,t-s)*exp(iws)         
  !   do iw=1,params%Nwr
  !      Gf(iw) = sum(Gin(it,it:1:-1)*exp(xi*params%wr(iw)*params%t(1:it)))
  !   enddo
  !   Gf=Gf*params%dt
  ! end subroutine wigner_transform_time








  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################







  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the density at a given istant of time
  ! n(t)=-xi*G^<(t,t)
  !+-------------------------------------------------------------------+
  function measure_dens(g,self) result(dens)
    type(kb_contour_gf(*,*,*)),intent(in)        :: G
    type(kb_contour_gf(G%N,G%L,G%Lf)),intent(in) :: Self
    real(8)                                      :: dens
    integer                                      :: N
    N = cc_params%Nt
    dens = dimag(G%less(N,N))
  end function measure_dens

  ! function measure_dens_Nso(params,g,self) result(dens)
  !   type(kb_contour_gf)                 :: g(Nspin,Nspin,Norb,Norb)
  !   type(kb_contour_gf)                 :: self(Nspin,Nspin,Norb,Norb)
  !   type(kb_contour_params)             :: params
  !   real(8)                             :: dens(Nspin,Norb)
  !   integer                             :: N,ispin,iorb
  !   N = params%Nt
  !   do ispin=1,Nspin
  !      do iorb=1,Norb
  !         dens(ispin,iorb) = dimag(G(ispin,ispin,iorb,iorb)%less(N,N))
  !      enddo
  !   enddo
  ! end function measure_dens_Nso



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the double occupancy at a given istant of time
  ! d(t)=n_up(t)*n_do(t)-1/U0*[Self^M*G^M]
  !      n_up(t)*n_do(t)-i/U*[Self^R*G^< + Self^<*G^A + Self^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_docc(g,self) result(docc)
    type(kb_contour_gf(*,*,*)),intent(in)        :: G
    type(kb_contour_gf(G%N,G%L,G%Lf)),intent(in) :: Self
    real(8)                                      :: docc
    integer                                      :: i,k,j,N,L
    complex(8),dimension(:),allocatable          :: SxG
    real(8)                                      :: nt
    N = cc_params%Nt
    L = cc_params%Ntau
    !
    nt   = dimag(G%less(N,N))
    allocate(SxG(0:max(N,L)))
    docc = nt**2
    if(N==1)then
       if(ui/=0.d0)then
          do k=0,L
             SxG(k)=Self%mats(L-k)*G%mats(k)
          end do
          docc=docc-1.d0/Ui*cc_params%dtau*kb_trapz(SxG(0:),0,L)
       endif
    else
       if(u/=0.d0)then
          do k=0,L
             SxG(k)=Self%lmix(N,k)*conjg(G%lmix(N,L-k))
          end do
          docc=docc + 1.d0/U*cc_params%dtau*dimag( (-xi)*kb_trapz(SxG(0:),0,L) )
          do k=1,N
             SxG(k)=Self%ret(N,k)*G%less(k,N)
          end do
          docc=docc + 1.d0/U*cc_params%dt*dimag(kb_trapz(SxG(0:),1,N))
          do k=1,N
             SxG(k)=Self%less(N,k)*conjg(G%ret(N,k))
          end do
          docc=docc + 1.d0/U*cc_params%dt*dimag(kb_trapz(SxG(0:),1,N))
       endif
    endif
    deallocate(SxG)
  end function measure_docc


  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! E_k(t)=2*Im[G^R*G^< + G^<*G^A + G^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_ekin(g,self) result(ekin)
    type(kb_contour_gf(*,*,*)),intent(in)        :: G
    type(kb_contour_gf(G%N,G%L,G%Lf)),intent(in) :: Self
    real(8)                                      :: ekin
    integer                                      :: i,k,j,N,L
    complex(8),dimension(:),allocatable          :: Ker
    real(8)                                      :: nt
    N = cc_params%Nt
    L = cc_params%Ntau
    !
    allocate(Ker(0:max(N,L)))
    if(N==1)then
       do k=0,L
          Ker(k)=G%mats(L-k)*G%mats(k)
       end do
       ekin = -2.d0*cc_params%dtau*kb_trapz(Ker(0:),0,L)
    else
       do k=0,L
          Ker(k)=G%lmix(N,k)*conjg(G%lmix(N,L-k))
       end do
       ekin=2.d0*cc_params%dtau*dimag( (-xi)*kb_trapz(Ker(0:),0,L) )
       do k=1,N
          Ker(k)=G%ret(N,k)*G%less(k,N)
       end do
       ekin=ekin + 2.d0*cc_params%dt*dimag(kb_trapz(Ker(0:),1,N))
       do k=1,N
          Ker(k)=G%less(N,k)*conjg(G%ret(N,k))
       end do
       ekin=ekin + 2.d0*cc_params%dt*dimag(kb_trapz(Ker(0:),1,N))
    endif
    deallocate(Ker)
  end function measure_ekin



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! U(t)= U*docc(t) - n(t) + 1/4
  !+-------------------------------------------------------------------+
  function measure_epot(g,self) result(epot)
    type(kb_contour_gf(*,*,*)),intent(in)        :: G
    type(kb_contour_gf(G%N,G%L,G%Lf)),intent(in) :: Self
    real(8)                                      :: epot,docc,nt
    integer                                      :: i,k,j,N,L
    N = cc_params%Nt
    L = cc_params%Ntau
    !
    if(N==1)then
       nt   = measure_dens(g,self)
       docc = measure_docc(g,self)
       epot = Ui*(docc - nt + 0.25d0)
    else
       nt   = measure_dens(g,self)
       docc = measure_docc(g,self)
       epot = U*(docc - nt + 0.25d0)
    endif
  end function measure_epot



  !+-------------------------------------------------------------------+
  !PURPOSE: return the value of the kinetic energy at a given istant of time
  ! E_k(t)=2*Im[G^R*G^< + G^<*G^A + G^\lmix*G^\rmix](t,t)
  !+-------------------------------------------------------------------+
  function measure_etot(g,self) result(etot)
    type(kb_contour_gf(*,*,*)),intent(in)        :: G
    type(kb_contour_gf(G%N,G%L,G%Lf)),intent(in) :: Self
    real(8)                                      :: etot,ekin,epot
    ekin = measure_ekin(g,self)
    epot = measure_epot(g,self)
    etot = ekin + epot
  end function measure_etot




END MODULE NEQ_MEASURE
