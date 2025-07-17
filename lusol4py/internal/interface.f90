subroutine lu1fac_c(m, n, nelem, lena, luparm, parmlu, a, indc, indr, &
                    p, q, lenc, lenr, locc, locr, iploc, iqloc, ipinv, &
                    iqinv, w, inform)
    use lusol
    use lusol_precision, only : ip, rp
    implicit none
    integer(ip), intent(in) :: m, n, nelem, lena
    integer(ip), intent(inout) :: luparm(*), indc(*), indr(*), p(*), q(*)
    integer(ip), intent(inout) :: lenc(*), lenr(*), locc(*), locr(*)
    integer(ip), intent(inout) :: iploc(*), iqloc(*), ipinv(*), iqinv(*)
    real(rp), intent(inout) :: parmlu(*), a(*), w(*)
    integer(ip), intent(out) :: inform
  
    call lu1fac(m, n, nelem, lena, luparm, parmlu, a, indc, indr, &
                p, q, lenc, lenr, locc, locr, iploc, iqloc, ipinv, &
                iqinv, w, inform)
  end subroutine lu1fac_c
  