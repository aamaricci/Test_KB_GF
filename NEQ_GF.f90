MODULE NEQ_GF
  USE NEQ_CONTOUR
  USE NEQ_INPUT_VARS
  USE NEQ_GF_COMMON
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_IOTOOLS
  USE SF_LINALG,    only: zeye,inv
  USE SF_MISC, only: assert_shape
  implicit none
  private


  ! KADANOFF-BAYM CONTOUR GREEN'S FUNCTIONS:
  type :: kb_gf
     complex(8),dimension(:,:),allocatable :: less
     complex(8),dimension(:,:),allocatable :: ret
     complex(8),dimension(:,:),allocatable :: lmix
     real(8),dimension(:),allocatable      :: mats
     complex(8),dimension(:),allocatable   :: iw
     logical                               :: status=.false.
   contains
     procedure :: init        => init_kb_gf_params   !constructor
     procedure :: extrapolate => extrapolate_kb_gf   !extrapolate
     procedure :: free        => free_kb_gf          !destructor
     procedure :: check       => check_kb_gf         !check status
     procedure :: del         => del_kb_gf           !set to zero at the boundary
     final :: deallocate_kb_gf                       !destructor final
  end type kb_gf


  interface operator (.plus.)
     module procedure :: kb_gf_plus_gf
  end interface operator (.plus.)

  interface operator (.minus.)
     module procedure :: kb_gf_minus_gf
  end interface operator (.minus.)

  interface operator (+)
     module procedure :: kb_gf_add_gf
  end interface operator (+)

  interface operator (-)
     module procedure :: kb_gf_subtract_gf
  end interface operator (-)


  ! PRODUCT WITH SCALAR (dble,complx):
  interface operator(*)
     module procedure :: kb_gf_left_times_int
     module procedure :: kb_gf_left_times_dble
     module procedure :: kb_gf_left_times_cmplx
  end interface operator(*)


  ! DIVISION BY SCALAR (dble,complx):
  interface operator(/)
     module procedure :: kb_gf_right_division_int
     module procedure :: kb_gf_right_division_dble
     module procedure :: kb_gf_right_division_cmplx
  end interface operator(/)


  !EQUALITY with scalar and function (G=F, G=c)
  interface assignment(=)
     module procedure :: kb_gf_equal_scalar
     module procedure :: kb_gf_equal_gf
  end interface assignment(=)


  intrinsic :: sum
  interface sum
     module procedure :: flat_kb_gf
  end interface sum



  ! CONVOLUTION:
  interface operator(.x.)
     module procedure :: convolute_kb_contour_gf_rank0
  end interface operator(.x.)





  public :: kb_gf
  !
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  public :: operator(.plus.)
  public :: operator(.minus.)
  !
  public :: sum                 !function sum along the boundary
  !
  public :: operator(.x.)       !convolution
  !
  public :: assignment(=)


contains





  ! INIT
  elemental subroutine init_kb_gf_params(G)
    class(kb_gf),intent(inout) :: G
    integer     :: N,L,Lf
    if(allocated(G%less))deallocate(G%less)
    if(allocated(G%ret)) deallocate(G%ret)
    if(allocated(G%lmix))deallocate(G%lmix)
    if(allocated(G%mats))deallocate(G%mats)
    if(allocated(G%iw))deallocate(G%iw)
    N = cc_params%Ntime            !<== allocate at maximum time
    L = cc_params%Ntau
    Lf= cc_params%Niw
    allocate(G%less(N,N))  ; G%less=zero
    allocate(G%ret(N,N))   ; G%ret=zero
    allocate(G%lmix(N,0:L)); G%lmix=zero
    allocate(G%mats(0:L))  ; G%mats=0d0
    allocate(G%iw(Lf))     ; G%iw=zero
    G%status=.true.    
  end subroutine init_kb_gf_params


  !FREE
  elemental subroutine free_kb_gf(G)
    class(kb_gf),intent(inout) :: G
    if(allocated(G%less))deallocate(G%less)
    if(allocated(G%ret)) deallocate(G%ret)
    if(allocated(G%lmix))deallocate(G%lmix)
    if(allocated(G%mats))deallocate(G%mats)
    if(allocated(G%iw))deallocate(G%iw)
    G%status=.false.
  end subroutine free_kb_gf


  elemental subroutine del_kb_gf(G)
    class(kb_gf),intent(inout) :: G
    integer                    :: N
    N   = cc_params%Nt   !<== work with the ACTUAL size of the contour
    if(N==1)then
       G%iw   = zero
       G%mats = 0d0
    endif
    G%ret(N,1:N)   = zero
    G%less(N,1:N)  = zero
    G%lmix(N,0:)   = zero
  end subroutine del_kb_gf



  !DESTRUCTOR
  subroutine deallocate_kb_gf(G)
    type(kb_gf) :: G
    if(allocated(G%less))deallocate(G%less)
    if(allocated(G%ret)) deallocate(G%ret)
    if(allocated(G%lmix))deallocate(G%lmix)
    if(allocated(G%mats))deallocate(G%mats)
    if(allocated(G%iw))deallocate(G%iw)
    G%status=.false.
  end subroutine deallocate_kb_gf





  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################


  !C(t,t')=A(t,t') +/- B(t,t'), with t,t'=0,t_max
  !INTERFACE OPERATORS .PLUS., .MINUS.
  !PLUS
  elemental function kb_gf_plus_gf(a,b) result(c)
    type(kb_gf), intent(in) :: a,b
    type(kb_gf)             :: c
    call c%init()
    c%less = a%less + b%less
    c%ret  = a%ret  + b%ret
    c%lmix = a%lmix + b%lmix
    c%mats = a%mats + b%mats
    c%iw   = a%iw   + b%iw
  end function kb_gf_plus_gf

  !MINUS
  elemental function kb_gf_minus_gf(a,b) result(c)
    type(kb_gf), intent(in) :: a,b
    type(kb_gf)             :: c
    call c%init()
    c%less = a%less - b%less
    c%ret  = a%ret  - b%ret
    c%lmix = a%lmix - b%lmix
    c%mats = a%mats - b%mats
    c%iw   = a%iw   - b%iw
  end function kb_gf_minus_gf




  ! !################################################################################
  ! !################################################################################
  ! !################################################################################
  ! !################################################################################
  ! !################################################################################



  !C(t,t')=A(t,t') +/- B(t,t') along contour, t=t_max && t'=0,t_max
  !INTERFACE OPERATORS +, -
  !ADD:
  elemental function kb_gf_add_gf(a,b) result(c)
    type(kb_gf),intent(in)  :: a,b
    type(kb_gf)             :: c
    integer                 :: N
    N   = cc_params%Nt   !<== work with the ACTUAL size of the contour
    call c%init()
    if(N==1)then
       C%mats(0:) = A%mats(0:) + B%mats(0:)
       C%iw(:)    = A%iw(:)    + B%iw(:)
    endif
    C%ret(N,1:N)   = A%ret(N,1:N)   + B%ret(N,1:N)
    C%less(N,1:N)  = A%less(N,1:N)  + B%less(N,1:N)
    C%lmix(N,0:)   = A%lmix(N,0:)   + B%lmix(N,0:)
    ! !THIS SHOULD NOT BE INVOLVED IN THE CALCULATION:
    ! C%less(1:N-1,N)= -conjg(C%less(N,1:N-1))
  end function kb_gf_add_gf

  !SUBTRACT:
  elemental function kb_gf_subtract_gf(a,b) result(c)
    type(kb_gf), intent(in) :: a,b
    type(kb_gf)             :: c
    integer                 :: N
    N   = cc_params%Nt   !<== work with the ACTUAL size of the contour
    call c%init()
    if(N==1)then
       C%mats(0:) = A%mats(0:) - B%mats(0:)
       C%iw(:)    = A%iw(:)    - B%iw(:)
    endif
    C%ret(N,1:N)   = A%ret(N,1:N)   - B%ret(N,1:N)
    C%less(N,1:N)  = A%less(N,1:N)  - B%less(N,1:N)
    C%lmix(N,0:)   = A%lmix(N,0:)   - B%lmix(N,0:)
    ! !THIS SHOULD NOT BE INVOLVED IN THE CALCULATION:
    ! C%less(1:N-1,N)= -conjg(C%less(N,1:N-1))
  end function kb_gf_subtract_gf





  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################





  !EQUALITY with scalar and function (G=F, G=c)
  elemental subroutine kb_gf_equal_scalar(A,C)
    type(kb_gf),intent(inout) :: A
    complex(8),intent(in)     :: C
    A%less = C
    A%ret  = C
    A%lmix = C
    A%mats = C
    A%iw   = C
  end subroutine kb_gf_equal_scalar

  elemental subroutine kb_gf_equal_gf(A,B)
    type(kb_gf),intent(inout) :: A
    type(kb_gf),intent(in)    :: B
    A%less(:,:)  = B%less(:,:)
    A%ret(:,:)   = B%ret(:,:)
    A%lmix(:,0:) = B%lmix(:,0:)
    A%mats(0:)   = B%mats(0:)
    A%iw(:)      = B%iw(:)
  end subroutine kb_gf_equal_gf



  ! !##################################################################
  ! !##################################################################
  ! !##################################################################
  ! !##################################################################



  ! PRODUCT WITH SCALAR (int,dble,complx):
  elemental function kb_gf_left_times_int(C,Gin) result(Gout)
    integer,intent(in)     :: C
    type(kb_gf),intent(in) :: Gin
    type(kb_gf)            :: Gout
    call Gout%init
    Gout%less(:,:) = C*Gin%less(:,:)
    Gout%ret(:,:)  = C*Gin%ret(:,:)
    Gout%lmix(:,0:)= C*Gin%lmix(:,0:)
    Gout%mats(0:)  = C*Gin%mats(0:)
    Gout%iw(:)     = C*Gin%iw(:)
  end function kb_gf_left_times_int

  elemental function kb_gf_left_times_dble(C,Gin) result(Gout)
    real(8),intent(in)     :: C
    type(kb_gf),intent(in) :: Gin
    type(kb_gf)            :: Gout
    call Gout%init
    Gout%less(:,:) = C*Gin%less(:,:)
    Gout%ret(:,:)  = C*Gin%ret(:,:)
    Gout%lmix(:,0:)= C*Gin%lmix(:,0:)
    Gout%mats(0:)  = C*Gin%mats(0:)
    Gout%iw(:)     = C*Gin%iw(:)
  end function kb_gf_Left_times_dble

  elemental function kb_gf_left_times_cmplx(C,Gin) result(Gout)
    complex(8),intent(in)  :: C
    type(kb_gf),intent(in) :: Gin
    type(kb_gf)            :: Gout
    call Gout%init
    Gout%less(:,:) = C*Gin%less(:,:)
    Gout%ret(:,:)  = C*Gin%ret(:,:)
    Gout%lmix(:,0:)= C*Gin%lmix(:,0:)
    Gout%mats(0:)  = C*Gin%mats(0:)
    Gout%iw(:)     = C*Gin%iw(:)
  end function kb_gf_left_times_cmplx



  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################



  ! DIVISION BY SCALAR (int,dble,complx):
  elemental function kb_gf_right_division_int(Gin,C) result(Gout)
    integer,intent(in)     :: C
    type(kb_gf),intent(in) :: Gin
    type(kb_gf)            :: Gout
    call Gout%init
    Gout%less(:,:) = Gin%less(:,:)/C
    Gout%ret(:,:)  = Gin%ret(:,:)/C
    Gout%lmix(:,0:)= Gin%lmix(:,0:)/C
    Gout%mats(0:)  = Gin%mats(0:)/C
    Gout%iw(:)     = Gin%iw(:)/C
  end function kb_gf_right_division_int

  elemental function kb_gf_right_division_dble(Gin,C) result(Gout)
    real(8),intent(in)     :: C
    type(kb_gf),intent(in) :: Gin
    type(kb_gf)            :: Gout
    call Gout%init
    Gout%less(:,:) = Gin%less(:,:)/C
    Gout%ret(:,:)  = Gin%ret(:,:)/C
    Gout%lmix(:,0:)= Gin%lmix(:,0:)/C
    Gout%mats(0:)  = Gin%mats(0:)/C
    Gout%iw(:)     = Gin%iw(:)/C
  end function kb_gf_right_division_dble

  elemental function kb_gf_right_division_cmplx(Gin,C) result(Gout)
    complex(8),intent(in)  :: C
    type(kb_gf),intent(in) :: Gin
    type(kb_gf)            :: Gout
    call Gout%init
    Gout%less(:,:) = Gin%less(:,:)/C
    Gout%ret(:,:)  = Gin%ret(:,:)/C
    Gout%lmix(:,0:)= Gin%lmix(:,0:)/C
    Gout%mats(0:)  = Gin%mats(0:)/C
    Gout%iw(:)     = Gin%iw(:)/C
  end function kb_gf_right_division_cmplx



  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################
  !################################################################################



  !Flat function
  function flat_kb_gf(A) result(C)
    class(kb_gf),intent(in) :: A(:)
    type(kb_gf)             :: C
    integer                 :: i,N
    !
    call C%init
    N   = cc_params%Nt   !<== work with the ACTUAL size of the contour
    !
    do i=1,size(a)
       if(N==1)then
          C%mats(0:) = C%mats(0:) + A(i)%mats(0:)
          C%iw(:)    = C%iw(:)    + A(i)%iw(:)
       endif
       C%ret(N,1:N)   = C%ret(N,1:N)  + A(i)%ret(N,1:N)
       C%less(N,1:N)  = C%less(N,1:N) + A(i)%less(N,1:N)
       C%lmix(N,0:)   = C%lmix(N,0:)  + A(i)%lmix(N,0:)
    enddo
    ! !THIS SHOULD NOT BE INVOLVED IN THE CALCULATION:
    ! C%less(1:N-1,N)= -conjg(C%less(N,1:N-1))
  end function flat_kb_gf




  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################


  elemental function check_kb_gf(G) result(bool)
    class(kb_gf),intent(in) :: G
    logical                 :: bool
    integer                 :: N,L,Lf
    logical                 :: check(4)
    !
    N=cc_params%Ntime              !<== check size at maximum time
    L=cc_params%Ntau
    !
    check(1) = any(shape(G%less)/=[N,N])
    check(2) = any(shape(G%ret)/=[N,N])
    check(3) = any(shape(G%lmix)/=[N,L+1])
    check(4) = any(shape(G%mats)/=[L+1])
    bool = .not.any(check)
  end function check_kb_gf



  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################







  !EXTRAPOLATE
  elemental subroutine extrapolate_kb_gf(G)
    class(kb_gf),intent(inout) :: G
    integer                    :: i,j,k,N,L
    N   = cc_params%Nt      !<== work with the ACTUAL size of the contour
    L   = cc_params%Ntau
    select case(N)
    case(1)
       return
    case(2)
       !GUESS G AT THE NEXT STEP, GIVEN THE INITIAL CONDITIONS
       do j=1,N
          g%ret(N,j) =g%ret(1,1)
          g%less(N,j)=g%less(1,1)
       end do
       do i=1,N-1
          g%less(i,N)=g%less(1,1)
       end do
       do j=0,L
          g%lmix(N,j)=g%lmix(1,j)
       end do
    case default
       !EXTEND G FROM THE [N-1,N-1] TO THE [N,N] SQUARE TO START DMFT
       !USING QUADRATIC EXTRAPOLATION
       do k=1,N-1
          g%less(N,k)=2.d0*g%less(N-1,k)-g%less(N-2,k)
          g%less(k,N)=2.d0*g%less(k,N-1)-g%less(k,N-2)
       end do
       g%less(N,N)=2.d0*g%less(N-1,N-1)-g%less(N-2,N-2)
       !
       do k=0,L
          g%lmix(N,k)=2.d0*g%lmix(N-1,k)-g%lmix(N-2,k)
       end do
       !
       g%ret(N,N)=-xi
       do k=1,N-2
          g%ret(N,k)=2.d0*g%ret(N-1,k)-g%ret(N-2,k)
       end do
       g%ret(N,N-1)=0.5d0*(g%ret(N,N)+g%ret(N,N-2))
    end select
  end subroutine extrapolate_kb_gf





  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################




  !----------------------------------------------------------------------------
  !C(t,t')=(A*B)(t,t'), with t=t_max && t'=0,t_max
  !----------------------------------------------------------------------------
  function convolute_kb_contour_gf_Rank0(A,B) result(C)
    type(kb_gf), intent(in)             :: a,b
    type(kb_gf)                         :: c
    integer                             :: N,L
    real(8)                             :: dt,dtau
    real(8),dimension(:),allocatable    :: ftau
    complex(8),dimension(:),allocatable :: AxB    
    integer                             :: i,j,k,itau,jtau
    !
    call C%init()
    !
    N   = cc_params%Nt   !<== work with the ACTUAL size of the contour
    L   = cc_params%Ntau
    dt  = cc_params%dt
    dtau= cc_params%dtau
    !
    allocate(AxB(0:max(L,N)))
    !
    if(N==1)then
       !Mats components:
       !C(iw)  =  A(iw)*B(iw) = FT[int_0^beta ds A(tau) B(tau-s)]
       C%iw(:) = A%iw(:)*B%iw(:)
       allocate(ftau(0:Niw))
       call fft_iw2tau(C%iw,ftau(0:),cc_params%beta,notail=.true.)
       call fft_extract_gtau(ftau(0:),C%mats(0:))
       !
       !Ret. component
       !C^R(t,t')=\int_{t'}^t ds A^R(t,s)*B^R(s,t') ==0 
       C%ret(1,1)=zero
       !
       !Lmix. component
       !C^\lmix(t,tau')=\int_0^{beta} ds A^\lmix(t,s)*B^M(s,tau')
       !               = I1
       C%lmix(N,0:L)=zero
       do jtau=0,L
          AxB(0:) = zero
          do k=0,jtau
             AxB(k)=A%lmix(1,k)*B%mats(L+k-jtau)
          end do
          C%lmix(1,jtau)=C%lmix(1,jtau)-dtau*kb_trapz(AxB(0:),0,jtau)
          do k=jtau,L
             AxB(k)=A%lmix(1,k)*B%mats(k-jtau)
          end do
          C%lmix(1,jtau)=C%lmix(1,jtau)+dtau*kb_trapz(AxB(0:),jtau,L)
       enddo
       !
       !Less component
       !C^<(t,t')=-i\int_0^{beta} ds A^\lmix(t,s)*B^\rmix(s,t')
       ! tip (j=1,N-1)
       C%less(1,1)=zero
       AxB(0:) = zero
       do k=0,L
          AxB(k)=A%lmix(1,k)*conjg(B%lmix(1,L-k))
       end do
       C%less(1,1)=C%less(1,1)-xi*dtau*kb_trapz(AxB(0:),0,L)
       !
       return
       !
    endif
    !
    !
    !Ret. component
    !C^R(t,t')=\int_{t'}^t ds A^R(t,s)*B^R(s,t')
    C%ret(N,1:N)=zero
    do j=1,N           !for all t' in 0:t_max
       AxB(0:)  = zero !AxB should be set to zero before integration
       do k=j,N        !store the convolution between t'{=j} and t{=N}
          AxB(k) = A%ret(N,k)*B%ret(k,j)
       enddo
       C%ret(n,j) = C%ret(n,j) + dt*kb_trapz(AxB(0:),j,N)
    enddo
    !
    !Lmix. component
    !C^\lmix(t,tau')=\int_0^{beta} ds A^\lmix(t,s)*B^M(s,tau')
    !                  +\int_0^{t} ds A^R(t,s)*B^\lmix(s,tau')
    !               = I1 + I2
    C%lmix(N,0:L)=zero
    do jtau=0,L
       !I1:
       !break the integral I1 in two parts to take care of the sign of (tau-tau').
       AxB(0:) = zero
       do k=0,jtau
          AxB(k)=A%lmix(N,k)*B%mats(L+k-jtau)
       end do
       C%lmix(N,jtau)=C%lmix(N,jtau)-dtau*kb_trapz(AxB(0:),0,jtau)
       do k=jtau,L
          AxB(k)=A%lmix(N,k)*B%mats(k-jtau)
       end do
       C%lmix(n,jtau)=C%lmix(n,jtau)+dtau*kb_trapz(AxB(0:),jtau,L)
       !
       !I2:
       AxB(0:) = zero
       do k=1,N
          AxB(k) = A%ret(N,k)*B%lmix(k,jtau)
       enddo
       C%lmix(N,jtau) = C%lmix(N,jtau) + dt*kb_trapz(AxB(0:),1,N)
    enddo
    !
    !Less component
    !C^<(t,t')=-i\int_0^{beta} ds A^\lmix(t,s)*B^\rmix(s,t')
    !             +\int_0^{t'} ds A^<(t,s)*B^A(s,t')
    !             +\int_0^{t} ds A^R(t,s)*B^<(s,t')
    ! (t,t')=>(N,j) <==> Vertical side, no tip (j=1,N-1)
    do j=1,N-1
       C%less(N,j)=zero
       !
       AxB(0:) = zero
       do k=0,L
          AxB(k)=A%lmix(N,k)*conjg(B%lmix(j,L-k))
       end do
       C%less(N,j)=C%less(N,j)-xi*dtau*kb_trapz(AxB(0:),0,L)
       !
       AxB(0:) = zero
       do k=1,j
          AxB(k)=A%less(N,k)*conjg(B%ret(j,k))
       end do
       C%less(N,j)=C%less(N,j)+dt*kb_trapz(AxB(0:),1,j)
       !
       AxB(0:) = zero
       do k=1,N
          AxB(k)=A%ret(N,k)*B%less(k,j)
       end do
       C%less(N,j)=C%less(N,j)+dt*kb_trapz(AxB(0:),1,N)
       !
    end do
    !
    ! (t,t')=>(i,N) <==> Horizontal side, w/ tip (i=1,N)
    do i=1,N
       C%less(i,N)=zero
       do k=0,L
          AxB(k)=A%lmix(i,k)*conjg(B%lmix(n,L-k))
       end do
       C%less(i,N)=C%less(i,N)-xi*dtau*kb_trapz(AxB(0:),0,L)
       !
       do k=1,N
          AxB(k)=A%less(i,k)*conjg(B%ret(N,k))
       end do
       C%less(i,N)=C%less(i,N)+dt*kb_trapz(AxB(0:),1,N)
       !
       do k=1,i
          AxB(k)=A%ret(i,k)*B%less(k,N)
       end do
       C%less(i,N)=C%less(i,N)+dt*kb_trapz(AxB(0:),1,i)
    end do
    !      
    deallocate(AxB)
  end function convolute_kb_contour_gf_Rank0

END MODULE NEQ_GF













! !Constructor: this conflicts with assignment(=) overload
! interface kb_gf
!    module procedure :: new_kb_gf_int
!    module procedure :: new_kb_gf_params
! end interface kb_gf

! !CONSTRUCTOR: THIS *CONFLICTS WITH ASSIGNMENT(=) OVERLOAD*
! elemental function new_kb_gf_int(N,L,Lf) result(G)
!   type(kb_gf)        :: G
!   integer,intent(in) :: N,L,Lf
!   if(allocated(G%less))deallocate(G%less)
!   if(allocated(G%ret)) deallocate(G%ret)
!   if(allocated(G%lmix))deallocate(G%lmix)
!   if(allocated(G%mats))deallocate(G%mats)
!   if(allocated(G%iw))deallocate(G%iw)
!   allocate(G%less(N,N))  ; G%less=zero
!   allocate(G%ret(N,N))   ; G%ret=zero
!   allocate(G%lmix(N,0:L)); G%lmix=zero
!   allocate(G%mats(0:L))  ; G%mats=0d0
!   allocate(G%iw(Lf))     ; G%iw=zero
!   G%status=.true.    
! end function new_kb_gf_int

! elemental function new_kb_gf_params() result(G)
!   type(kb_gf) :: G
!   integer     :: N,L,Lf
!   if(allocated(G%less))deallocate(G%less)
!   if(allocated(G%ret)) deallocate(G%ret)
!   if(allocated(G%lmix))deallocate(G%lmix)
!   if(allocated(G%mats))deallocate(G%mats)
!   if(allocated(G%iw))deallocate(G%iw)
!   N = cc_params%Ntime            !<== allocate at maximum time
!   L = cc_params%Ntau
!   Lf= cc_params%Niw
!   allocate(G%less(N,N))  ; G%less=zero
!   allocate(G%ret(N,N))   ; G%ret=zero
!   allocate(G%lmix(N,0:L)); G%lmix=zero
!   allocate(G%mats(0:L))  ; G%mats=0d0
!   allocate(G%iw(Lf))     ; G%iw=zero
!   G%status=.true.    
! end function new_kb_gf_params
