module miev0mod
  use iso_c_binding
  use mathutils
  implicit none
  
  
contains
  
  subroutine miev0easydriver(xx, mre, mim, qext, qsca, gqsc, pmom, mxmdm) bind (c)
    ! dummy parameters
    real(C_DOUBLE), intent(in), value ::  xx, mre, mim
    integer(C_INT), intent(in), value ::  mxmdm
    real(C_DOUBLE), intent(out)       ::  qext, qsca, gqsc
    real(C_DOUBLE), intent(out)       ::  pmom(0:mxmdm,4)
    
    ! actually used parametrs
    LOGICAL  ANYANG, PERFCT, PRNT(2)
    INTEGER  IPOLZN, MOMDIM, NUMANG, NMOM
    REAL     GQSC_S, MIMCUT, QEXT_S, QSCA_S, SPIKE, XMU
    COMPLEX  CREFIN, SFORW, SBACK, S1(1), S2(1), TFORW(2), TBACK(2)
    REAL, ALLOCATABLE ::  PMOM_TMP(:,:)
    complex(kind=dp)  ::  midx
    
    EXTERNAL MIEV0
    INTEGER I, J
    
    
    midx = dcmplx(mre, mim)
    
    MOMDIM = mxmdm
    
    ANYANG = .TRUE.
    PERFCT = .FALSE.
    PRNT   = .FALSE.
    NUMANG = 1
    NMOM   = MOMDIM!3*INT(XX)
    CREFIN = CMPLX(midx)
    MIMCUT = 1.0E-8
    IPOLZN = +1234
    ALLOCATE(PMOM_TMP(0:MOMDIM, 4))
    
    XMU = 0.0
    
    call MIEV0( XX, CREFIN, PERFCT, MIMCUT, ANYANG, NUMANG, XMU,&
                NMOM, IPOLZN, MOMDIM, PRNT, QEXT_S, QSCA_S, GQSC_S,&
                PMOM_TMP, SFORW, SBACK, S1, S2, TFORW, TBACK,&
                SPIKE )
    
    ! COPY VALUES
    qext = DBLE(QEXT_S)
    qsca = DBLE(QSCA_S)
    gqsc = DBLE(GQSC_S)
    
    ! COPY MOMENTS
    DO I=1, 4
      DO J=0, MOMDIM
        pmom(J,I) = DBLE(PMOM_TMP(j,I))
      END DO
    END DO
    
    
    DEALLOCATE(PMOM_TMP)
    
  end subroutine miev0easydriver
  
  subroutine miesdist(xs, ys, N, mre, mim, wl, pmom, momdim, ext, sca, asy, vol, ierr) bind (c)
    real(C_DOUBLE), intent(in)   ::  xs(N), ys(N)
    real(C_DOUBLE), intent(in), value ::  mre, mim, wl
    real(C_DOUBLE), intent(inout)::  pmom(0:momdim, 4)
    real(C_DOUBLE), intent(out)  ::  ext, sca, asy, vol
    integer(C_INT), intent(in), value  ::  momdim, N
    integer(C_INT), intent(out)        ::  ierr
    
    integer     ::  nsize, I,J,L
    complex(kind=dp)  :: midx
    real(kind=dp) ::  k, xx, xsquared, xcubed, qext, qsca, gqsc
    real(kind=dp), allocatable  ::  pmom_tmp(:,:,:), extinction(:),&
                                    scattering(:), asymetry(:), volume(:),&
                                    pmom_i(:,:)
    real(kind=dp) :: norm    
    
    midx = dcmplx(mre, mim)
    
    nsize =N
    
    k = 2.0_dp*pi/wl
    
    allocate(pmom_tmp(nsize, 0:momdim, 4), pmom_i(0:momdim, 4))
    allocate(extinction(nsize), scattering(nsize), asymetry(nsize),&
              volume(nsize))
              
    do I=1, nsize
      xx = k*xs(I)
      
      call miev0easydriver(xx, mre, mim, qext, qsca, gqsc, pmom_i, momdim)
      xsquared = xs(I)*xs(I)*1.0D-12
      xcubed = xsquared*xs(I)*1.0D-6
      
      extinction(I) = qext*xsquared*ys(I)
      scattering(I) = qsca*xsquared*ys(I)
      asymetry(I) = gqsc*scattering(I)
      volume(I) = xcubed*ys(I)
      
      do J=1, 4
        do L=0, momdim
          pmom_tmp(I, L, J) = pmom_i(L,J)*ys(I)
        end do
      end do
    
    end do
    
    norm = trapz(xs, ys)
    
    ext = trapz(xs, extinction) / norm
    
    sca = trapz(xs, scattering) / norm
    
    asy = trapz(xs, asymetry) / norm
    
    vol = trapz(xs, volume) / norm
    
    DO I=1, 4
      DO J=0, momdim
        pmom(J,I) = trapz(xs, pmom_tmp(:, J, I)) / norm
      END DO
    END DO
    
    deallocate(pmom_tmp, pmom_i, extinction, scattering, asymetry, volume)
    ierr = 0
    return
    
  end subroutine miesdist
  
  subroutine wiscombe2evans(wisc, evans)
    real(kind=dp), intent(in) ::  wisc(0:, :)
    real(kind=dp), intent(out)::  evans(0:, :)
    real(kind=dp)             ::  norm, factor
    integer ::  I
    norm = wisc(0,1)+wisc(0,2)
    
    DO I=lbound(wisc,1), ubound(wisc,1)
      factor = dble(2*I+1)
      evans(I,1) = (wisc(I,1)+wisc(I,2))/norm*factor
      evans(I,2) = (wisc(I,2)-wisc(I,1))/norm*factor
      evans(I,3) = 2.0_dp*wisc(I,3)/norm*factor
      evans(I,4) = 2.0_dp*wisc(I,4)/norm*factor
      evans(I,5) = evans(I,1)
      evans(I,6) = evans(I,3)
    END DO 
    
  end subroutine wiscombe2evans
end module miev0mod
