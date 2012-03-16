program FMA
  !implicit none
  integer, parameter :: n=65
  integer, parameter :: ncycle=2
  integer :: p1, p2, nf, ng, nn, ngrid, nsteps, tot, ns
  real :: d, l, bc, C, dt
  integer ::  i, j, k, jj, jpost, jpre, jcycle !indicies
  real, dimension(2*(n**3)) :: iu, irhs, ires, irho, s
  integer :: t1,t2,clock_rate,clock_max
  nsteps=300
  d = 1.0/(n-1)
  dt = 0.001
  C = 0.001
  p1 = 1
  p2 = 1
  tot=1
  nf=n
  ng=0
  bc = 0 !set constant boundary condition
  !find number of grid levels
  nn=n
  !count the levels of grid resizing
  do while(nn.gt.2)
     ng = ng+1
     nn = (nn+1)/2
  enddo
  !initialize u and the source term
  call initialize(iu(p1),s(p1),nf,d,bc)
  iu = iu+s*dt
  ns=0
  nn=n
  do ns=1, nsteps
     call system_clock(t1, clock_rate, clock_max)
     print*, ns
     irho=iu
     do while (nn.gt.3)
        tot= tot + nn**3 !the running index placement counter
        nn=(nn+1)/2 !the current nc    
        p1=p2 !fine pointer
        p2=tot !coarse pointer  
        call rstrct(irho(p2), irho(p1), nn)
     enddo
     !setup and solve the coarsest grid
     nn=3
     call slvsml(iu(p2), irho(p2), C)
     ngrid=ng     
     do j=2,ngrid
        nn=2*nn-1 !changed to the level j
        p1=tot !coarse, lvl j-1
        tot = tot - nn**3 !changed to the level j
        p2=tot !fine
        call interp(iu(p2), iu(p1), nn)
        !fill rhs at appropriate level
        if (j.ne.ngrid) then
           !copy the j level into irhs
           do i=p2, p1-1
              irhs(i) = irho(i)
           enddo
        else
           do i=p2, p1-1
              !copy the ngrid value in iu
              irhs(i) = iu(i) + s(i)*dt
           enddo
        endif
        !v-cycle at current grid level
        do jcycle=1,ncycle
           nf = nn
           do jj=j,2,-1
              do jpre=1,2
                 call relax(iu(p2), irhs(p2), nf,C)
              enddo
              call resid(ires(p2), iu(p2), irhs(p2), nf, C)
              nf=(nf+1)/2
              call rstrct(irhs(p1), ires(p2), nf)
              do i=p1, (p1+nf**3)-1
                 iu(i)=0.0
              enddo
              p2=p1
              p1=p1+nf**3
           enddo
           call slvsml(iu(p2),irhs(p2),C) !p2 should be the index of the smallest subarray
           nf=3
           do jj=2,j
              nf=2*nf-1
              p1=p2 !p1 points to coarse
              p2=p2-nf**3 !p2 should now point to fine     
              call addint(iu(p2), iu(p1), ires(p2), nf)
              do jpost=1,2
                 call relax(iu(p2),irhs(p2),nf,C)
              enddo
           enddo !jj end   
        enddo !jcycle end
     enddo !j end
     call system_clock(t2,clock_rate,clock_max)
     print*, 'FMA timer:'
     write(*,*) real(t2-t1)/real(clock_rate)
  enddo !nsteps end
end program FMA

subroutine relax(u,rhs,n,C)
  implicit none
  integer :: i, ipass, isw, j, jsw, k, ksw, n
  real :: C
  real, dimension(n,n,n) :: u, uold, rhs
  !jacobi
  uold=u
  do k=2,n-1
     do j=2,n-1
        do i=2,n-1
           u(i,j,k) = C/(6.*C+1.)*(uold(i-1,j,k) + uold(i+1,j,k) + uold(i,j-1,k)&
                + uold(i,j+1,k) + uold(i,j,k-1) + uold(i,j,k+1)) + 1./(6.*C+1.)*rhs(i,j,k)
        enddo
     enddo
  enddo

  !gauss seidel
  !do k=2,n-1
  !   do j=2,n-1
  !      do i=2,n-1
  !         u(i,j,k) = C/(6.*C+1.)*(u(i-1,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i,j+1,k)&
  !              + u(i,j,k-1) + u(i,j,k+1)) + 1./(6.*C+1.)*rhs(i,j,k)
  !      enddo
  !   enddo
  !enddo

  !gauss seidel with red black ordering
  !do ipass=1,2
  !   jsw=ipass
  !   do k=2,n-1
  !      jsw=3-jsw
  !      isw=jsw
  !      do j=2,n-1
  !         isw=3-isw
  !         do i=isw+1,n-1,2
  !            u(i,j,k) = C/(6.*C+1.)*(u(i-1,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i,j+1,k)&
  !                 + u(i,j,k-1) + u(i,j,k+1)) + 1./(6.*C+1.)*(rhs(i,j,k)) 
  !         enddo
  !      enddo
  !   enddo
  !enddo

end subroutine relax

subroutine slvsml(u, rhs, C)
  implicit none
  real, dimension(3,3,3) :: u, rhs
  real :: C
  integer :: i,j,k
  !solve the coarsest level (3x3x3)
  do k=1,3
     do j=1,3
        do i=1,3
           u(i,j,k) = 0.0
        enddo
     enddo
  enddo
  u(2,2,2) = rhs(2,2,2)/(1.+6.*C) 
end subroutine slvsml

subroutine resid(res, u, rhs, n, C)
  implicit none
  integer :: n, i, j, k
  real :: C
  real, dimension(n,n,n) :: res, u, rhs
  !calculate residual and store in res
  do k=2,n-1
     do j=2,n-1
        do i=2,n-1
           res(i,j,k) = -C*(u(i-1,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i,j+1,k)&
                + u(i,j,k-1) + u(i,j,k+1) - 6.*u(i,j,k)) - 1.*u(i,j,k) + 1.
        enddo
     enddo
  enddo
end subroutine resid

subroutine addint(uf, uc, res, nf)
  implicit none
  integer :: nc, nf, i, j, k
  real, dimension(nf,nf,nf) :: uf, res
  real, dimension((nf+1)/2,(nf+1)/2,(nf+1)/2) :: uc
  !make finer grid and add the residual at every point
  call interp(res, uc, nf)
  do k=1,nf
     do j=1,nf
        do i=1,nf
           uf(i,j,k) = uf(i,j,k) + res(i,j,k)
        enddo
     enddo
  enddo
end subroutine addint

subroutine interp(uf, uc, nf)
  implicit none
  integer :: nc, nf, iif, ic, jf, jc, kf, kc
  real, dimension((nf+1)/2,(nf+1)/2,(nf+1)/2) :: uc
  real, dimension(nf,nf,nf) :: uf

  nc=(nf+1)/2
  !elements that are copies
  kf = 1
  do kc=1,nc
     jf=1
     do jc=1,nc
        iif=1
        do ic=1,nc
           uf(iif,jf,kf) = uc(ic,jc,kc)
           iif = iif+2
        enddo
        jf = jf+2
     enddo
     kf=kf+2
  enddo
  !do even rows, every other column, every other sheet
  do kf=1, nf, 2
     do jf=1, nf, 2
        do iif = 2, nf-1, 2
           uf(iif,jf,kf) = (0.5)*(uf(iif+1,jf,kf)+uf(iif-1,jf,kf))
        enddo
     enddo
  enddo
  !do even columns, every row, every other sheet
  do kf=1, nf, 2
     do jf=2, nf-1, 2
        do iif=1, nf
           uf(iif,jf,kf) = (0.5)*(uf(iif,jf+1,kf)+uf(iif,jf-1,kf))
        enddo
     enddo
  enddo
  !do every other sheet
  do kf=2, nf-1, 2
     do jf=1, nf
        do iif=1, nf
           uf(iif,jf,kf) = (0.5)*(uf(iif,jf,kf+1)+uf(iif,jf,kf-1))
        enddo
     enddo
  enddo
  
end subroutine interp
  
subroutine rstrct(uc, uf, nc)
  implicit none
  integer :: nc, nf, iif, ic, jf, jc, kf, kc
  real, dimension((2*nc)-1, (2*nc)-1, (2*nc)-1) :: uf
  real, dimension(nc, nc, nc) :: uc

  nf=(2*nc)-1
  !do the middle points
  kf=3
  do kc=2, nc-1
     jf=3
     do jc=2, nc-1
        iif=3
        do ic=2, nc-1
           uc(ic,jc,kc) = (1./2.)*uf(iif,jf,kf)+(1./12.)*(uf(iif+1,jf,kf)+uf(iif-1,jf,kf)&
                +uf(iif,jf+1,kf)+uf(iif,jf-1,kf)+uf(iif,jf,kf+1)+uf(iif,jf,kf-1))
           iif=iif+2
       enddo
        jf=jf+2
     enddo
     kf=kf+2
  enddo

  !do the boundaries
  !constant i faces
  kf=1
  do kc=1, nc
     jf=1
     do jc=1, nc
        uc(1,jc,kc) = uf(1,jf,kf)
        uc(nc,jc,kc) = uf(nf,jf,kf)
        jf=jf+2        
     enddo
     kf=kf+2
  enddo
  !constant j faces
  kf=1
  do kc=1, nc
     iif=1
     do ic=1, nc
        uc(ic,1,kc) = uf(iif,1,kf)
        uc(ic,nc,kc) = uf(iif,nf,kf)
        iif=iif+2        
     enddo
     kf=kf+2
  enddo
  !constant k faces
  jf=1
  do jc=1, nc
     iif=1
     do ic=1, nc
        uc(ic,jc,1) = uf(iif,jf,1)
        uc(ic,jc,nc) = uf(iif,jf,nf)
        iif=iif+2 
     enddo
     jf=jf+2
  enddo
end subroutine rstrct

subroutine initialize(uf, s, nf, d, bc)
  implicit none
  integer :: nf, i, j, k
  real :: d, xx, yy, zz, bc
  real, dimension(nf,nf,nf) :: uf,s

  !initialize to gaussian
  do k=2, nf-1
     do j=2, nf-1
        do i=2,nf-1
           xx=d*(i-1)
           yy=d*(j-1)
           zz=d*(k-1)
           uf(i,j,k) = exp(-(5*xx - 2.5)**2) * exp(-(5*yy - 2.5)**2) * exp(-(5*zz - 2.5)**2)
           s(i,j,k)=0
        enddo
     enddo
  enddo

  !boundaries
  do k=1, nf
     do j=1, nf
        uf(1,j,k) = bc
        uf(nf,j,k) = bc
        s(1,j,k) = 0.0
        s(nf,j,k) = 0.0
     enddo
  enddo
  do k=1, nf
     do i=1, nf
        uf(i,1,k) = bc
        uf(i,nf,k) = bc
        s(i,1,k) = 0.0
        s(i,nf,k) = 0.0
     enddo
  enddo
  do j=1, nf
     do i=1, nf
        uf(i,j,1) = bc
        uf(i,j,nf) = bc
        s(i,j,1) = 0.0
        s(i,j,nf) = 0.0
     enddo
  enddo

  !add in the source term in the center of the 3D box
  s((nf+1)/2,(nf+1)/2,(nf+1)/2) = 0.0 !should be redundant
end subroutine initialize
