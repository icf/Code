      FUNCTION rannyu(x)
      implicit real*8(a-h,o-z)
      parameter(dismin=1.d-11)
      parameter(tpm12 = 1.d0/4096.d0)
c      common/lrannyu/l1,l2,l3,l4
c      common/nrannyu/n1,n2,n3,n4
      data m1,m2,m3,m4 / 502,1521,4071,2107/
      data l1,l2,l3,l4 /   0,   0,   0, 127/
*      data l1,l2,l3,l4 / 2221, 1813, 1968,  723/
      data n1,n2,n3,n4 /    0,    0, 2896, 1263/
*      data n1,n2,n3,n4 /    0,    0, 2896, 1243/
*      data n1,n2,n3,n4 /    0,    0, 2896, 1237/
*      data n1,n2,n3,n4 /    0,    0, 2896, 1233/
*      data n1,n2,n3,n4 /    0,    0, 2896, 1221/
*      data n1,n2,n3,n4 /    0,    0, 2896, 1197/
      i1=l1*m4+l2*m3+l3*m2+l4*m1 + n1
      i2=l2*m4+l3*m3+l4*m2 + n2
      i3=l3*m4+l4*m3 + n3
      i4=l4*m4 + n4
      l4=and(i4,4095)
      i3=i3+ishft(i4,-12)
      l3=and(i3,4095)
      i2=i2+ishft(i3,-12)
      l2=and(i2,4095)
      l1=and(i1+ishft(i2,-12),4095)
      rannyu=tpm12*(l1+tpm12*(l2+tpm12*(l3+tpm12*l4)))
      if(rannyu.eq.0.d0)rannyu=dismin
      return
      end
