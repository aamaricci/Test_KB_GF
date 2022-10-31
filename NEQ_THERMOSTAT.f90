module NEQ_THERMOSTAT
  USE NEQ_INPUT_VARS
  USE NEQ_CONTOUR
  USE NEQ_CONTOUR_GF
  USE SF_ARRAYS
  USE SF_SPECIAL
  USE SF_IOTOOLS
  implicit none
  private

  interface get_thermostat_bath
     module procedure :: get_thermostat_bath_rank0
  end interface get_thermostat_bath
  public :: get_thermostat_bath


contains
  !+-------------------------------------------------------------------+
  !PURPOSE  : Build the Bath part of the system using exclusively time
  !dependent formalism. The bath is not interacting so it is 
  ! time translation invariant.
  !+-------------------------------------------------------------------+
  subroutine get_thermostat_bath_main(S0)
    type(kb_contour_gf(*,*,*)),intent(inout) :: S0
    integer                                  :: N,L
    integer                                  :: iw,i,j,k,ints
    real(8)                                  :: en,w,dw,wmax
    real(8)                                  :: ngtr,nless,arg,Vhopping
    complex(8)                               :: peso,sless,sgtr    
    complex(8),dimension(:,:),allocatable    :: S0gtr
    logical                                  :: check
    real(8),dimension(Lbath)                 :: bath_dens,wfreq
    !
    N = cc_params%Ntime      !<== work with the MAX size of the contour
    L = cc_params%Ntau
    !
    !
    ! Lambda = Vhopping**2/2W
    ! Vhopping**2 = Lambda*2*W
    Vhopping=sqrt(Lambda*2*Wbath)
    write(*,"(A)")"Bath hopping:"//str(Vhopping)    
    wmax  = Wbath
    wfreq = linspace(-wmax,wmax,Lbath,mesh=dw)
    select case(reg(bath_type))
    case("bethe")
       call get_bath_bethe_dos( bath_dens,wfreq)
    case("gaussian")
       call get_bath_gaussian_dos( bath_dens,wfreq)
    case ("flat")
       call get_bath_flat_dos( bath_dens,wfreq)
    case("pgflat")
       call get_bath_pgflat_dos(bath_dens,wfreq)
    case("gapflat")
       wmax  = Wbath+Wgap       !2.d0*Wbath+Wgap
       wfreq = linspace(-wmax,wmax,L,mesh=dw)
       call get_bath_gapflat_dos( bath_dens,wfreq)
    case default
       stop "Bath type: not supported."
    end select
    call splot("DOSbath.neqipt",wfreq,bath_dens)
    !
    !Bath self-energies:
    allocate(S0gtr(N,N))
    S0   = dcmplx(0d0,0d0)
    S0gtr= dcmplx(0d0,0d0)
    if(Lambda==0.d0)return
    do iw=1,Lbath
       en   = wfreq(iw)
       nless= fermi(en,beta)
       ngtr = fermi(en,beta)-1.d0 !it absorbs the minus sign of the greater functions
       do i=1,N
          !Less component:
          do j=1,i
             peso=exp(-xi*(cc_params%t(i)-cc_params%t(j))*en)
             S0%less(i,j)= S0%less(i,j) + xi*Vhopping**2*nless*peso*bath_dens(iw)*dw
             S0gtr(i,j) = S0gtr(i,j)    + xi*Vhopping**2*ngtr*peso*bath_dens(iw)*dw
          enddo
          !Lmix component:
          if(en>=0.d0)then
             do j=0,L
                if(beta*en>20.d0)then
                   peso=exp(-xi*en*cc_params%t(i))*exp(-en*(beta-cc_params%tau(j)))
                else
                   peso=exp(-xi*en*cc_params%t(i))*exp(-en*(beta-cc_params%tau(j)))/(1.d0+exp(-en*beta))
                endif
                S0%lmix(i,j) = S0%lmix(i,j) + xi*Vhopping**2*peso*bath_dens(iw)*dw
             enddo
          else
             do j=0,L
                if(beta*en<-20.d0)then
                   peso=exp(-xi*en*cc_params%t(i))*exp(en*cc_params%tau(j))
                else
                   peso=exp(-xi*en*cc_params%t(i))*exp(-en*(beta-cc_params%tau(j)))/(1.d0+exp(-en*beta))
                endif
                S0%lmix(i,j) = S0%lmix(i,j) + xi*Vhopping**2*peso*bath_dens(iw)*dw
             enddo
          endif
          !
       enddo
    enddo
    !Ret component:
    S0%ret = S0gtr-S0%less
    forall(i=1:cc_params%Ntime,j=1:cc_params%Ntime,i<j)S0%less(i,j)=-conjg(S0%less(j,i))
  end subroutine get_thermostat_bath_main



  subroutine get_thermostat_bath_rank0(S0,iplot)
    type(kb_contour_gf(*,*,*)),intent(inout) :: S0
    logical,optional                         :: iplot
    logical                                  :: iplot_
    iplot_=.true.;if(present(iplot))iplot_=iplot
    S0 = zero
    if(Lambda>0d0)then
       call get_thermostat_bath_main(S0)
       if(iplot_)call plot_kb_contour_gf(S0,"Sbath.neqipt")
    end if
  end subroutine get_thermostat_bath_rank0

  ! subroutine get_thermostat_bath_rank2(params,S0,iplot)
  !   type(kb_contour_gf)     :: S0(:,:)
  !   type(kb_contour_params) :: params
  !   logical,optional        :: iplot
  !   integer                 :: Nso,io
  !   logical                 :: iplot_
  !   iplot_=.true.;if(present(iplot))iplot_=iplot
  !   Nso = size(S0,1)
  !   call assert_shape_kb_contour_gf(S0,[Nso,Nso],"get_thermostat_bath_rank2","S0")
  !   if(Lambda==0d0)then
  !      call del_kb_contour_gf(params,S0)
  !   else
  !      write(*,"(A)")"Generating diagonal Bath:"
  !      do io=1,Nso
  !         call get_thermostat_bath_main(params,S0(io,io))
  !      end do
  !      if(iplot_)call plot_kb_contour_gf(params,S0,"Sbath.neqipt")
  !   end if
  ! end subroutine get_thermostat_bath_rank2




  !+-----------------------------------------------------------------+
  !PURPOSE  : Build constant BATH 
  !+-----------------------------------------------------------------+
  subroutine get_bath_flat_dos(bath_dens,wfreq)
    real(8),dimension(Lbath) :: bath_dens,wfreq
    integer                  :: i
    real(8)                  :: w
    do i=1,Lbath
       w=wfreq(i)
       bath_dens(i)= step(Wbath-abs(w))/(2.d0*Wbath)
    enddo
  end subroutine get_bath_flat_dos

  subroutine get_bath_pgflat_dos( bath_dens,wfreq)
    real(8),dimension(Lbath) ::  bath_dens,wfreq
    integer                  :: i
    real(8)                  :: w,norm
    norm=(Walpha+1.d0)/(2.d0*Wbath**(Walpha+1.d0))
    do i=1,Lbath
       w=wfreq(i)
       if(abs(w)>wbath)then
          bath_dens(i)=0.d0
       else
          bath_dens(i)=norm*abs(w)**Walpha
       end if
    end do
  end subroutine get_bath_pgflat_dos

  subroutine get_bath_gapflat_dos( bath_dens,wfreq)
    real(8),dimension(Lbath) ::  bath_dens,wfreq
    integer                  :: i
    real(8)                  :: w,rho
    rho=1.d0/2.d0/(Wbath)!-Wgap)
    do i=1,Lbath
       w=wfreq(i)
       if(abs(w)<Wbath+Wgap.AND.abs(w)>Wgap)then
          bath_dens(i)=rho
       else
          bath_dens(i)=0.d0         
       end if
    end do
  end subroutine get_bath_gapflat_dos

  subroutine get_bath_gaussian_dos( bath_dens,wfreq)
    real(8),dimension(Lbath) ::  bath_dens,wfreq
    integer                  :: i,ik
    real(8)                  :: w,sig,alpha
    complex(8)               :: gf,zeta
    bath_dens = exp(-0.5d0*(wfreq/Wbath)**2)/(sqrt(pi2)*Wbath) !standard Gaussian
    !bath_dens = exp(-((wfreq)/Wbath)**2)/(sqrt(pi)*Wbath) !Camille's choice
  end subroutine get_bath_gaussian_dos

  subroutine get_bath_bethe_dos( bath_dens,wfreq)
    real(8),dimension(Lbath) ::  bath_dens,wfreq
    integer                  :: i,ik
    real(8)                  :: w,sig,alpha
    complex(8)               :: gf,zeta
    do i=1,Lbath
       w=wfreq(i)
       zeta=cmplx(w,eps,8)
       gf=gfbether(w,zeta,wbath/2.d0)
       bath_dens(i)=-dimag(gf)/pi
    enddo
  end subroutine get_bath_bethe_dos




end module NEQ_THERMOSTAT

