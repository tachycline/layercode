      subroutine fft(a,b,is,n,id)
      real*8  a(0:*),b(0:*)
      integer :: i,j,k,l,is,n,id,lh,ip
      real*8 :: tr,ti,s,c,ur,ui

      j=0
      do i=0,(n-2)*is,is
         if (i.lt.j) then
            tr = a(j)
            a(j) = a(i)
            a(i) = tr
            ti = b(j)
            b(j) = b(i)
            b(i) = ti
         endif
         k=is*n/2
 10      if (k.le.j) then
            j = j-k
            k = k/2
            goto 10
         end if
         j = j+k
      end do
      s = 0.0d0
      c = -1.0d0
      l = is
 30   lh = l
      l = l+l
      ur = 1.0d0
      ui = 0.0d0
      do j=0,lh-is,is
         do i=j,(n-1)*is,l
            ip = i+lh
            tr = a(ip)*ur-b(ip)*ui
            ti = a(ip)*ui+b(ip)*ur
            a(ip) = a(i)-tr
            b(ip) = b(i)-ti
            a(i) = a(i)+tr
            b(i) = b(i)+ti
         enddo
         ti = ur*s+ui*c
         ur = ur*c-ui*s
         ui = ti
      end do
      s = sqrt(0.5d0*(1.0d0-c))*id
      c = sqrt(0.5d0*(1.0d0+c))
      if (l.lt.n*is) goto 30

      return
      end subroutine fft
c  ---------------------------------------------------------------------
