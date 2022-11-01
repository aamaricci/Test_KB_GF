MODULE KB_GF_MAIN
  USE KB_VARS_GLOBAL
  USE KB_CONTOUR
  USE KB_GF_AUX, &
       fft_iw2tau_kb_gf=>fft_iw2tau,  &
       fft_tau2iw_kb_gf=>fft_tau2iw,  &
       fft_extract_tau_kb_gf=> fft_extract_gtau
  USE KB_GF_COMMON
  USE KB_GF_SUM
  USE KB_GF_CONVOLUTE
  USE KB_GF_VIE
  USE KB_GF_DYSON
  USE KB_GF_FREE
  !USE KB_BUBBLES
  private

  public :: kb_gf
  public :: kb_dgf
  !
  public :: operator(+)
  public :: operator(-)
  public :: operator(*)
  public :: operator(/)
  public :: operator(.plus.)
  public :: operator(.minus.)
  public :: operator(.x.)
  public :: assignment(=)
  !
  public :: assert_shape_kb_gf
  public :: reshape_kb_gf
  public :: convolute_kb_gf
  public :: vie_kb_gf
  public :: free_kb_gf
  public :: dyson_kb_gf
  public :: reduce_kb_gf        !reduce as a callable function 
  public :: sum                 !extend intrinsic sum

  !FFT 
  public :: fft_iw2tau_kb_gf
  public :: fft_tau2iw_kb_gf
  public :: fft_extract_tau_kb_gf
END MODULE KB_GF_MAIN
