    lbt = lb11
    call rfftwnd_f77_one_complex_to_real(c2r3d,lbt,ignore_me)
    tmp = lvx
    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
    tmp = tmp * lbt * const
    call rfftwnd_f77_one_real_to_complex(r3c3d,tmp,ignore_me)
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      tmp(ii,jj,kk) = eye * kx(ii) * tmp(ii,jj,kk)
    end do
    end do
    end do
    lwb11 = lwb11 - tmp

    tmp = lvy
    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
    tmp = tmp * lbt * const
    call rfftwnd_f77_one_real_to_complex(r3c3d,tmp,ignore_me)
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      tmp(ii,jj,kk) = eye * ky(jj) * tmp(ii,jj,kk)
    end do
    end do
    end do
    lwb11 = lwb11 - tmp

    tmp = lvz
    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
    tmp = tmp * lbt * const
    call rfftwnd_f77_one_real_to_complex(r3c3d,tmp,ignore_me)
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      tmp(ii,jj,kk) = eye * kz(kk) * tmp(ii,jj,kk)
    end do
    end do
    end do
    lwb11 = lwb11 - tmp

    lbt = lb12
    call rfftwnd_f77_one_complex_to_real(c2r3d,lbt,ignore_me)

    tmp = lvx
    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
    tmp = tmp * lbt * const
    call rfftwnd_f77_one_real_to_complex(r3c3d,tmp,ignore_me)
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      tmp(ii,jj,kk) = eye * kx(ii) * tmp(ii,jj,kk)
    end do
    end do
    end do
    lwb12 = lwb12 - tmp

    tmp = lvy
    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
    tmp = tmp * lbt * const
    call rfftwnd_f77_one_real_to_complex(r3c3d,tmp,ignore_me)
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      tmp(ii,jj,kk) = eye * ky(jj) * tmp(ii,jj,kk)
    end do
    end do
    end do
    lwb12 = lwb12 - tmp

    tmp = lvz
    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
    tmp = tmp * lbt * const
    call rfftwnd_f77_one_real_to_complex(r3c3d,tmp,ignore_me)
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      tmp(ii,jj,kk) = eye * kz(kk) * tmp(ii,jj,kk)
    end do
    end do
    end do
    lwb12 = lwb12 - tmp

    lbt = lb13
    call rfftwnd_f77_one_complex_to_real(c2r3d,lbt,ignore_me)

    tmp = lvx
    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
    tmp = tmp * lbt * const
    call rfftwnd_f77_one_real_to_complex(r3c3d,tmp,ignore_me)
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      tmp(ii,jj,kk) = eye * kx(ii) * tmp(ii,jj,kk)
    end do
    end do
    end do
    lwb13 = lwb13 - tmp

    tmp = lvy
    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
    tmp = tmp * lbt * const
    call rfftwnd_f77_one_real_to_complex(r3c3d,tmp,ignore_me)
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      tmp(ii,jj,kk) = eye * ky(jj) * tmp(ii,jj,kk)
    end do
    end do
    end do
    lwb13 = lwb13 - tmp

    tmp = lvz
    call rfftwnd_f77_one_complex_to_real(c2r3d,tmp,ignore_me)
    tmp = tmp * lbt * const
    call rfftwnd_f77_one_real_to_complex(r3c3d,tmp,ignore_me)
    do kk = 1, lz
    do jj = 1, ly
    do ii = 1, lx1
      tmp(ii,jj,kk) = eye * kz(kk) * tmp(ii,jj,kk)
    end do
    end do
    end do
    lwb13 = lwb13 - tmp


