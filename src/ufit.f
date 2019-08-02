C Output from Public domain Ratfor, version 1.0
      subroutine ufit(xk,wk,xmode,x,w,mse,x1,w1,x2,w2,ind,kt,n,goof)
      implicit double precision(a-h,o-z)
      logical goof
      double precision mse
      dimension xk(n), wk(n), x(n), w(n),x1(n), w1(n), x2(n), w2(n), ind
     *(n), kt(n)
      if(xmode .lt. 0)then
      m = n-1
      x0 = 1.5d0
      xmax = -1.d0
      ssemin = 1.d200
      do23002 i = 1,m 
      do23004 j = 1,n 
      x(j) = xk(j)
      w(j) = wk(j)
23004 continue
23005 continue
      call unimode(x,w,x1,w1,x2,w2,ind,kt,x0,n,goof)
      if(goof)then
      return
      endif
      sse = 0.d0
      do23008 j = 1,n 
      sse = sse + (x(j)-xk(j))**2
23008 continue
23009 continue
      if(sse .lt. ssemin)then
      ssemin = sse
      xmax = x0
      endif
      x0 = x0+1.d0
23002 continue
23003 continue
      k1 = int(xmax-0.5d0)
      k2 = int(xmax+0.5d0)
      else
      xmax = xmode
      endif
      do23012 j = 1,n 
      x(j) = xk(j)
      w(j) = wk(j)
23012 continue
23013 continue
      call unimode(x,w,x1,w1,x2,w2,ind,kt,xmax,n,goof)
      if(goof)then
      return
      endif
      if(xmode .lt. 0)then
      mse = ssemin/dble(n)
      if(x(k1).ge.x(k2))then
      xmode=dble(k1)
      else
      xmode=dble(k2)
      endif
      else
      sse = 0.d0
      do23020 j = 1,n 
      sse = sse + (x(j)-xk(j))**2
23020 continue
23021 continue
      mse = sse/dble(n)
      endif
      return
      end
