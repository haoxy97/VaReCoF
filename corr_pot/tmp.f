c frontside only
      if (ipar.eq.1) then
         if ((alpha.gt.120.0d0).or.((phi.gt.90.0d0).and.
     $    (alpha.gt.10.0d0))) c2h3_h=1.0d20
      endif
c backside addition only
      if (ipar.eq.2) then
         if ((alpha.gt.90.0d0).or.(phi.lt.90.0d0).or.
     $    ((phi.gt.90.0d0).and.(alpha.lt.10.0d0))) c2h3_h=1.0d20
      endif
c abstraction only
      if (ipar.eq.3) then
         if ((alpha.lt.90.0d0).or.((alpha.lt.120.0d0).and.
     $    (phi.lt.90.0d0))) c2h3_h=1.0d20
      endif
c     write (6,*) 'alpha,tau',alpha,tau,ipar,c2h3_h

