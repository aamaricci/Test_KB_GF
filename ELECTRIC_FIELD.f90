MODULE ELECTRIC_FIELD
  USE NEQ_INPUT_VARS
  USE SF_CONSTANTS
  USE SF_SPECIAL
  USE SF_DERIVATE, only: deriv
  USE SF_IOTOOLS, only: free_unit,txtfy,reg
  implicit none
  private
  public :: Afield
  public :: set_efield_vector


contains

  !+------------------------------------------------------------+
  !PURPOSE: set the normalized electric field versors using given direction
  !+------------------------------------------------------------+
  subroutine set_efield_vector(time)
    real(8),dimension(:) :: time
    real(8)              :: modulo
    integer              :: i
    !Normalize the Electric Field components
    !Keep unaltered the Electric Field Strenght Efield=E0
    modulo=sqrt(dot_product(Evect,Evect))
    if(modulo/=0d0)Evect=Evect/modulo
    modulo=sqrt(dot_product(Evect1,Evect1))
    if(modulo/=0d0)Evect1=Evect1/modulo
    write(0,*)Afield(0d0)
    call print_field(time)
  end subroutine set_efield_vector


  subroutine print_field(t)
    real(8),dimension(:) :: t
    integer              :: i
    real(8),dimension(3) :: A
    real(8),dimension(size(t)) :: Ax,Ay,Ex,Ey
    do i=1,size(t)
       A=Afield(t(i))
       Ax(i)=A(1)
       Ay(i)=A(2)
    enddo
    Ex = deriv(Ax,dt)
    Ey = deriv(Ay,dt)
    open(10,file="Avector_shape.neqipt")
    open(11,file="Efield_shape.neqipt")
    do i=1,size(t)
       write(10,*)t(i),Afield(t(i))
       write(11,*)t(i),Ex(i),Ey(i)
    enddo
    close(10)
    close(11)
    if(field_type=="ac")print*,"Root condition: "//trim(txtfy(bessel_j0(Efield/Omega0)))
  end subroutine print_field



  !+------------------------------------------------------------+
  !PURPOSE : 
  !+------------------------------------------------------------+
  function Afield(t)
    real(8),intent(in)      :: t
    real(8)                 :: ftime,ctime,tau,tau1,time
    real(8),dimension(3)    :: Afield
    complex(8)              :: zp,zm
    select case(field_type)
    case ("dc")                !DC ELECTRIC FIELD:
       ftime=-(step(t-ton)*(t-ton + (toff-t)*step(t-toff) - (toff-ton)*step(ton-toff)))
       Afield=Evect*Efield*ftime

    case("ac")                  !AC ELECTRIC FIELD
       ftime=-sin(Omega0*(t-ton))/Omega0
       Afield=Evect*Efield*ftime 

    case("acdc")                !AC+DC ELECTRIC FIELD (super-bloch)
       !ftime=-(t+sin(Omega0*(t-ton))/Omega0)
       ftime =-sin(Omega0*(t-ton))/Omega0
       ctime =-(step(t-ton)*(t-ton + (toff-t)*step(t-toff) - (toff-ton)*step(ton-toff)))
       Afield=Evect*Efield*ftime + Evect1*E1*ctime

    case("pulse")               !LIGHT PULSE (for Pump&Probe) 
       !Signal function:
       !Cos(\Omega*(t-t0)).Exp(-((t-t0)/tau)**2)
       tau  = tau0/2d0/sqrt(log(2d0))
       time = t-tpulse
       zp   = dcmplx(time/tau, tau*Omega0/2)
       zp   = dcmplx(time/tau,-tau*Omega0/2)
       ftime= -0.5d0*sqrt(pi)*tau*Exp(-(tau*Omega0/2)**2)*(zerf(zp)-zerf(zm))
       Afield=Evect*Efield*ftime

    case("dcpulse")               !LIGHT PULSE + DC field (for Pump&Probe) 
       !Signal function: 
       !Cos(\Omega*(t-t0)).Exp(-((t-t0)/tau)**2) + Theta(t-ton)*Theta(toff-t)
       tau  = tau0/2d0/sqrt(log(2d0))
       time = t-tpulse
       zp   = dcmplx(time/tau, tau*Omega0/2)
       zp   = dcmplx(time/tau,-tau*Omega0/2)
       ftime= -0.5d0*sqrt(pi)*tau*Exp(-(tau*Omega0/2)**2)*(zerf(zp)-zerf(zm))
       ctime= -(step(t-ton)*(t-ton + (toff-t)*step(t-toff) - (toff-ton)*step(ton-toff)))
       Afield=Evect*Efield*ftime + Evect1*E1*ctime

    case ("dcdc")                !DC+DC ELECTRIC FIELD:
       ftime=-(step(t-ton)*(t-ton + (toff-t)*step(t-toff) - (toff-ton)*step(ton-toff)))
       ctime=-(step(t-ton)*(t-ton + (toff-t)*step(t-toff) - (toff-ton)*step(ton-toff)))
       Afield=Evect*Efield*ftime + Evect1*E1*ctime


    case("ramp")                !RAMP TO CONSTANT DC-FIELD:
       ftime=-(24.d0*pi*(t+(t-ton)*step(t-ton)+2.d0*(toff-t)*step(t-ton)*step(t-toff)-&
            2.d0*(ton-toff)*step(t-ton)*step(ton-toff))+                              &
            27.d0*ton*(step(t-ton)-1.d0)*Sin(pi*t/ton) - &
            ton*(step(t-ton)-1.d0)*Sin(3.d0*pi*t/ton))/48.d0/pi
       Afield=Evect*Efield*ftime       !A(t) = E0*F(t)*(e_x + e_y)

       !!add more here:
    case default
       stop "ELECTRIC_FIELD/Afield: wrong field_type. set:dc,ac,acdc,pulse,dcpulse,ramp"
    end select
    !-----------------------------
  end function Afield


end MODULE ELECTRIC_FIELD
