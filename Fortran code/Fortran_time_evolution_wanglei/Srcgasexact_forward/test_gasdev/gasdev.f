      double precision function gasdev()   ! -- gaussian random #
      implicit real*8(a-h,o-z)
      save
* ----------------------------------------------------------------
* Generates a random # distributed according to a gaussian with
* mean 0 and variance 1, i.e., p(x)=1/sqrt[2*pi] exp[-x^2/2]
* ----------------------------------------------------------------
      data iset/0/
      if (iset.eq.0) then
 1              V1=2.d0*rannyu()-1.d0
        V2=2.d0*rannyu()-1.d0
        R=V1**2+V2**2
        if(R.ge.1..or.R.eq.0.)GO TO 1
        fac=dsqrt(-2.d0*dlog(R)/R)
        gset=V1*fac
        gasdev=V2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      end if
      return
      end
