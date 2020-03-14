subroutine electroweakino(M1,M2,mu,tanb)
! 
! Read in electroweakino parameters, M1, M2, mu, tanb and construct the 
! neutralino mass matrix Y following Gounaris, Le Mouel, Porfyriadis (2002).
! For good measure also do the same for the chargino mass matrix, X.
!
! Note that the SLHA paper, Skands et al (2004), (hep-ph/0311123) 
! published in JHEP, says "Different spectrum calculators may disagree 
! on the overall sign of one or more rows in a mixing matrix, 
! owing to different diagonalization algorithms. Such differences 
! correspond to a flip of the sign of the eigenvectors in question and 
! do not lead to inconsistencies. Only the relative sign between entries 
! on the same row is physically significant for processes with 
! interfering amplitudes"
!
! History
!
! 1) I initially implemented this in Fortran as the initial 
!    goal (see basic.f90) was simply to format the neutralino and 
!    chargino mass matrices in formats suitable for subsequent 
!    Mathematica gymnastics. I chose Fortran as it has powerful 
!    formatting abilities (that I am very familiar with). 
!    I found that more modern versions of Fortran beyond f77 had support 
!    for command-line arguments, so decided to use modern Fortran for a change.
!
! 2) In subsequent iterations where I wanted to do the matrix 
!    gymnastics internally, I added Fortran LAPACK/BLAS library features.
!
! 3) In the latest iteration, I decided to use ROOT as the histogramming 
!    package, but kept this underlying Fortran code base for the 
!    calculation. So I have implemented a mixed programming paradigm, 
!    with a C++ wrapper, and am using a struct/common block for 
!    information exchange.
!    See mixedProgramming pdfs for more details.
!
  implicit none

  integer, parameter :: real_8_30 = selected_real_kind(p=8,r=30)
  integer :: i
  integer :: irow,jcol
  integer :: iordered(4)
  real(real_8_30), intent(in) :: M1, M2, mu, tanb
  real(real_8_30) :: mZ, mW, s2W, cW, sW, beta, cb, sb
  real(real_8_30) :: v1(4),v2(4),v3(4),v4(4)               ! Neutralino compositions
  real(real_8_30) :: Y(4,4)
  real(real_8_30) :: X(2,2)
! Chargino mass matrix gymnastics. Should be such that Z = U X V is diagonal.
  real(real_8_30) :: XXT(2,2)   ! X X^T - leads to the U matrix
  real(real_8_30) :: XTX(2,2)   ! X^T X - leads to the V matrix
  real(real_8_30) :: U(2,2)
  real(real_8_30) :: V(2,2)
  real(real_8_30) :: D(2,2)
  real(real_8_30) :: XV(2,2)
  real(real_8_30) :: Identity(2,2)
! LAPACK SSYEV arguments
  character(len=1) :: jobz, uplo
  integer :: info, lda, lwork, n
  real(real_8_30) :: w(4), work(100)
  real(real_8_30) :: eigsu(2),eigsv(2)
! Mass spectrum (ordered)
  real(real_8_30) :: nmass(4)
  real(real_8_30) :: cmass(2)
  real(real_8_30) :: snmass(4)
!  real(real_8_30) :: n11sq,n12sq,n13sq,n14sq
  common/mval/nmass,cmass,v1,v2,v3,v4,U,V,snmass
  
!  logical :: lprint
!  parameter ( lprint = .false. )
  include 'lprint.inc'

! Initialize a 2*2 identity matrix
  Identity(1,1) = 1.0d0
  Identity(1,2) = 0.0d0
  Identity(2,1) = 0.0d0
  Identity(2,2) = 1.0d0

  if(lprint)print *,'M1=',M1,' M2=',M2,' mu=',mu,' tanb=',tanb

! constants
  mZ  = 91.1876d0
  mW  = 80.379d0
  s2W = 0.232d0
  sW = sqrt(s2W)
  cW = sqrt(1.0d0-s2W)
  beta = atan2(tanb,1.0d0)
  cb = cos(beta)
  sb = sin(beta)

! Play with El Kheishen formulae
  call elkheishen(M1,M2,mu,tanb,mZ,beta,s2W)

! Play with Gounaris formulae
  call gounaris(M1,M2,mu,tanb,mZ,beta,s2W)

! The Y matrix is symmetric so it doesn't matter whether 
! we are thinking correctly about column-major or row-major order. 
! Fortran is column-major in terms of storage.
! https://en.wikipedia.org/wiki/Row-_and_column-major_order
! But anyway Y(i,j) means Y(irow, jcol)
!
  Y(1,1) =  M1
  Y(1,2) =  0.0d0
  Y(1,3) = -mZ*sW*cb
  Y(1,4) =  mZ*sW*sb

  Y(2,1) =  0.0d0
  Y(2,2) =  M2
  Y(2,3) =  mZ*cW*cb
  Y(2,4) = -mZ*cW*sb

  Y(3,1) = -mZ*sW*cb
  Y(3,2) =  mZ*cW*cb
  Y(3,3) =  0.0d0
  Y(3,4) = -mu

  Y(4,1) =  mZ*sW*sb
  Y(4,2) = -mZ*cW*sb
  Y(4,3) = -mu
  Y(4,4) =  0.0d0 

! I originally followed the convention for X in equation C9 of 
! Haber and Kane, and eqn A2 of Barnett and Haber, 
! but the X in Steve Martin's SUSY Primer and in eqn 6 of the Fuks et al 
! paper has sin(beta) on the first row. 
! This difference of convention may explain why my U and V matrix seemed 
! different than expected from Martin.

  X(1,1) =  M2
  X(1,2) =  sqrt(2.0d0)*mW*sb
  X(2,1) =  sqrt(2.0d0)*mW*cb
  X(2,2) =  mu

  do irow=1,4,1
     do jcol=1,4,1
        if(lprint)print *,irow,jcol,Y(irow,jcol)
     enddo
  enddo

  if(lprint)then
     do irow=1,4
        if(irow.eq.1)then
           write(6,1001)Y(irow,1),Y(irow,2),Y(irow,3),Y(irow,4)
        elseif(irow.le.3)then
           write(6,1002)Y(irow,1),Y(irow,2),Y(irow,3),Y(irow,4)
        else
           write(6,1003)Y(irow,1),Y(irow,2),Y(irow,3),Y(irow,4)
        endif
     enddo
  endif
 
  if(lprint)then
     do irow=1,2,1
        do jcol=1,2,1
           print *,irow,jcol,X(irow,jcol)
        enddo
     enddo
  endif

  if(lprint)write(6,2001)X(1,1),X(1,2),X(2,1),X(2,2)

! Prepare LAPACK code for the neutralino matrix diagonalization
  jobz = 'V'
  uplo = 'U'
  lda = 4
  n = 4
  lwork = 100
!  call ssyev(jobz,uplo,n,Y,lda,w,work,lwork,info)
  call dsyev('V','U',4,Y,4,w,work,lwork,info)  ! this works too. 
                                               ! less hassle (if done correctly)

! Results
  if(lprint)print *,'info = ',info

  call sortthem(w,iordered)

  if(lprint)then
     do i=1,4
        write(6,1000)i,w(iordered(i))
     enddo
     print *,' '
  endif

! Also the eigenvectors (here sorted by eigenvalue)
  do i=1,4
     jcol = iordered(i)
     nmass(i)=abs(w(jcol))
! Keep track of signs
     snmass(i)=w(jcol)
     if(lprint)then
        write(6,3002)i,abs(w(jcol)),Y(1,jcol),Y(2,jcol),Y(3,jcol),Y(4,jcol)
     endif
     if(i.eq.1)then
        v1(1) = Y(1,jcol)
        v1(2) = Y(2,jcol)
        v1(3) = Y(3,jcol)
        v1(4) = Y(4,jcol)
     elseif(i.eq.2)then
        v2(1) = Y(1,jcol)
        v2(2) = Y(2,jcol)
        v2(3) = Y(3,jcol)
        v2(4) = Y(4,jcol)
     elseif(i.eq.3)then
        v3(1) = Y(1,jcol)
        v3(2) = Y(2,jcol)
        v3(3) = Y(3,jcol)
        v3(4) = Y(4,jcol)
     else
        v4(1) = Y(1,jcol)
        v4(2) = Y(2,jcol)
        v4(3) = Y(3,jcol)
        v4(4) = Y(4,jcol)
     endif
  enddo
! copy LSP components
!  n11sq = v1(1)*v1(1)
!  n12sq = v1(2)*v1(2)
!  n13sq = v1(3)*v1(3)
!  n14sq = v1(4)*v1(4)
 
  call charginomasses(M1,M2,mu,tanb)
  call checkdet(X)

! Now deal with the chargino mass matrices etc.
  call dgemm('N','T',2,2,2,1.0,X,2,X,2,0.0,XXT,2)  ! Form X X^T
  call dgemm('T','N',2,2,2,1.0,X,2,X,2,0.0,XTX,2)  ! Form X^T X
  if(lprint)call printmatrix('  X  ',X)
  if(lprint)call printmatrix('X X^T',XXT)  ! symmetric
  if(lprint)call printmatrix('X^T X',XTX)  ! symmetric
! Now eigenvalues/eigenvectors for each 
! - matrix gets overwritten with the matrix of eigenvectors
  call dsyev('V','U',2,XXT,2,eigsu,work,lwork,info)
  call dsyev('V','U',2,XTX,2,eigsv,work,lwork,info)
  if(lprint)call printvector('  u  ',eigsu)
  if(lprint)call printmatrix('(U)  ',XXT)
  if(lprint)call printvector('  v  ',eigsv)
  if(lprint)call printmatrix('(V)  ',XTX)
  call dgemm('N','N',2,2,2,1.0,XXT,2,Identity,2,0.0,U,2)  ! Make a copy called U 
  call dgemm('N','N',2,2,2,1.0,XTX,2,Identity,2,0.0,V,2)  ! Make a copy called V
  if(lprint)call printmatrix('  U  ',U)
  if(lprint)call printmatrix('  V  ',V)
  call dgemm('N','N',2,2,2,1.0,X,2,V,2,0.0,XV,2)  ! Form XV
  call dgemm('N','N',2,2,2,1.0,U,2,XV,2,0.0,D,2)  ! Form D=UXV
  if(lprint)call printmatrix('  D  ',D)
  cmass(1) = sqrt(eigsu(1))
  cmass(2) = sqrt(eigsu(2))
  if(lprint)then
     write(6,4001)M1,M2,mu,tanb
     write(6,5000)
     write(6,5001)nmass(1),abs(v1(1)),abs(v1(2)),abs(v1(3)),abs(v1(4))
     write(6,5001)nmass(2),abs(v2(1)),abs(v2(2)),abs(v2(3)),abs(v2(4))
     write(6,5001)nmass(3),abs(v3(1)),abs(v3(2)),abs(v3(3)),abs(v3(4))
     write(6,5001)nmass(4),abs(v4(1)),abs(v4(2)),abs(v4(3)),abs(v4(4))
     write(6,5002)
     write(6,5003)cmass(1),abs(U(1,1)),abs(U(2,1)),abs(V(1,1)),abs(V(2,1))  ! could be (2,1) or (1,2) ...
     write(6,5003)cmass(2),abs(U(1,2)),abs(U(2,2)),abs(V(1,2)),abs(V(2,2))  !
     write(6,*)'DMs ',cmass(1)-nmass(1),nmass(2)-cmass(1),nmass(2)-nmass(1) 
  endif
1000 format('Neutralino ',i1,3x,f10.4)
2001 format('X = {{',f10.4,',',f10.4,'} , {',f10.4,',',f10.4,'}}')
1001 format('Y = {{',f10.4,',',f10.4,',',f10.4,',',f10.4,'},')
1002 format('{',f10.4,',',f10.4,',',f10.4,',',f10.4,'},')
1003 format('{',f10.4,',',f10.4,',',f10.4,',',f10.4,'}}')
3001 format('v(',i1,'),(',f10.4,',',f10.4,',',f10.4,',',f10.4,')')
3002 format('v(',i1,')  ',f10.4,' GeV',/,  &
            '      ',f10.4,/,  &
            '      ',f10.4,/,  &
            '      ',f10.4,/,  &
            '      ',f10.4,/)
4001 format(/,' Tree-level Electroweakino Mass Spectrum and Mixings for',&
              ' M1=',f8.2,', M2=',f8.2,', mu=',f8.2,', tanb=',f5.2)
! ,/,   &
!            ' ~chi^0  ',4(1x,f10.3),/' ~chi^+- ',1x,f10.3,1x,f10.3)
5000 format(/,' ~chi^0              m    |~B|  |~W_3|  |~H_1|  |~H_2|')
5001 format(12x,f10.3,4(1x,f7.3))
5002 format(/,' ~chi^+-             m   U:|~W|   |~H|   V:|~W|   |~H|')
5003 format(12x,f10.3,4(1x,f7.3))

end subroutine electroweakino

subroutine elkheishen(M1,M2,mu,tanb,mZ,beta,s2W)
! Evaluate constants in Eqns 24-28 of El Kheishen, Shafik, Aboshousha (PRD45 (1992) 4345)
! This does not yet seem to be doing what it should.
! Suspect Gounaris, LeMouel, Profyriadis is cleaner - so will code that up in 
! gounaris subroutine
implicit none
   
  integer, parameter :: real_8_30 = selected_real_kind(p=8,r=30)
  real(real_8_30), intent(in) :: M1,M2,mu,tanb,mZ,beta,s2W
  real(real_8_30) :: c2,c3,c4
  real(real_8_30) :: msum,mt1,mt2,mt3,c2W
  real(real_8_30) :: D,U,S
  real(real_8_30) :: acubed,a
  real(real_8_30) :: myA, myB, myC, myD
  real(real_8_30) :: mass1, mass2, mass3, mass4

  mt1 = M1*M2 - mZ*mZ - mu*mu
  msum = M1+M2
  c2W = 1.0d0-s2W
  mt2 = M1*c2W + M2*s2W
  mt3 = mu*mZ*mZ*sin(2.0d0*beta)
   
  c2 = mt1 - 0.375d0*msum*msum
  c3 = -0.125d0*msum*msum*msum + 0.5d0*msum*mt1 + msum*mu*mu + mt2*mZ*mZ - mt3
  c4 = mt2*mt3 - M1*M2*mu*mu + 0.25d0*msum*(msum*mu*mu + mt2*mZ*mZ - mt3) 
  c4 = c4 + (msum*msum*mt1/16.0d0) - 3.0d0*(msum**4)/256.0d0 

  U = - ( 4.0d0*c4 + (c2*c2*c2/3.0d0) )
  S = - ( c3*c3 + (2.0d0*c2*c2*c2/27.0d0) - (8.0d0*c2*c4/3.0d0) )
  D = - ( (4.0d0*U*U*U) + (27.0d0*S*S))

  print *,' El Kheishen : c2,c3,c4 ',c2,c3,c4
  print *,' El Kheishen :  U, S, D ',U,S,D

! If D is negative this is OK
  
  acubed = -S - sqrt(-D/27.0d0)
  a = (0.5d0**(1.0d0/3.0d0))*(-acubed)**(1.0d0/3.0d0)
  print *, 'El Kheishen : acubed ',acubed,a

  myD = 0.25d0*msum

  print *,' a, c2, c3, myD ',a,c2,c3,myD

  myA = -sqrt( 0.5d0*a -(c2/6.0d0) )
  myB = -(0.5d0*a + (c2/3.0d0))
  myC = c3/sqrt(8.0d0*a - (8.0d0*c2/3.0d0))

  print *,' A, B, C, D ',myA,myB,myC,myD

! Over-write with correct answers
  myA = -496.1867d0
  myB = 258542.40d0
  myC = -253367.80d0

  print *,' A, B, C, D ',myA,myB,myC,myD

! Why is B so negative ?

  mass1 = myD + myA + sqrt( myB + myC )
  mass2 = myD - myA - sqrt( myB - myC )
  mass3 = myD + myA - sqrt( myB + myC )
  mass4 = myD - myA + sqrt( myB - myC )

  print *,' m1, m2, m3, m4 ',mass1,mass2,mass3,mass4

end subroutine elkheishen

subroutine gounaris(sM1,sM2,smu,stanb,smZ,sbeta,ss2W)
! Follow Gounaris, LeMouel, Profyriadis (Phys.Rev. D65 (2002) 035002) 
! for determining neutralino masses analytically 
implicit none
   
  integer, parameter :: real_8_30 = selected_real_kind(p=8,r=30)
  real(real_8_30), intent(in) :: sM1,sM2,smu,stanb,smZ,sbeta,ss2W
  real(real_8_30) :: M1,M2,mu,tanb,mZ,beta,s2W
  real(real_8_30) :: c2W
  real(real_8_30) :: bigA, bigB, bigC, bigD
  real(real_8_30) :: a,b,Delta
  real(real_8_30) :: c,c3phi,phi,theta
  real(real_8_30) :: E,F
  real(real_8_30) :: alphap,alpham,betap,betam
  real(real_8_30) :: mass1, mass2, mass3, mass4
  real(real_8_30) :: c1, c2, c3, c4
  real(real_8_30) :: d1, d2, d3, d4
  real(real_8_30) :: e1, e2, e3, e4

! Make everything DP from now
! Also make mass units 100 GeV
  M1 = 0.01d0*dble(sM1)
  M2 = 0.01d0*dble(sM2)
  mu = 0.01d0*dble(smu)
  tanb = dble(stanb)
  mZ = 0.01d0*dble(smZ)
  beta = dble(sbeta)
  s2W = dble(ss2W)
  c2W = 1.0d0-s2W

! Eqn 11
  bigA = M1 + M2
  bigB = M1*M2 - mu*mu - mZ*mZ
  bigC = -M1*(mu*mu + mZ*mZ*c2W) -M2*(mu*mu + mZ*mZ*s2W)
  bigD = -M1*M2*mu*mu + mu*mZ*mZ*sin(2.0d0*beta)*(M1*c2W + M2*s2W)

! Eqn 12
  a = -0.25d0*( bigD + (bigB*bigB/12.0d0) - (bigA*bigC/4.0d0) )
  b =  0.25d0*( (-bigA*bigA*bigD/16.0d0) + (bigA*bigB*bigC/48.0d0) - (bigB*bigB*bigB/216.0d0) + &
                (bigB*bigD/6.0d0) - (bigC*bigC/16.0d0) )
! Eqn 13
  Delta = 0.25d0*b*b + (a*a*a/27.0d0)
  
  print *,'Gounaris: A,B,C,D   ',bigA,bigB,bigC,bigD
  print *,'Gounaris: a,b,Delta ',a,b,Delta 

! Eqn 14 (need to check that Delta and a have correct sign ... TODO)
  c = sqrt(3.0d0/abs(a))
  c3phi = -0.5d0*b*c*c*c
  phi = (acos(c3phi))/3.0d0
  theta = 2.0d0*sqrt(abs(a)/3.0d0)*cos(phi)

! Eqn 17
  E = 0.25d0*sqrt(bigA*bigA - (8.0d0*bigB/3.0d0) + 16.0d0*theta)
  F = (0.25d0/E)*(bigC - (bigA*bigB/6.0d0) -2.0d0*bigA*theta)
  print *,'Gounaris: theta,E,F ',theta,E,F

! Eqn 18
  alphap = (0.5d0*bigA + 2.0d0*E)
  alpham = (0.5d0*bigA - 2.0d0*E)
  betap  = -4.0d0*( (bigB/6.0d0) + 2.0d0*theta + F )
  betam  = -4.0d0*( (bigB/6.0d0) + 2.0d0*theta - F )
  print *,'Gounaris: alpha-+ ',alpham,alphap
  print *,'Gounaris:  beta-+ ',betam,betap

  mass1 = 0.5d0*( alpham - sqrt( alpham*alpham + betap ) )
  mass2 = 0.5d0*( alpham + sqrt( alpham*alpham + betap ) )
  mass3 = 0.5d0*( alphap - sqrt( alphap*alphap + betam ) )
  mass4 = 0.5d0*( alphap + sqrt( alphap*alphap + betam ) )

  print *,' m1, m2, m3, m4 (GeV) ',1.0d2*mass1,1.0d2*mass2,1.0d2*mass3,1.0d2*mass4

! Check how well each satisfies the quartic equation (Eqn 10).
  c1 = mass1**4 - bigA*mass1**3 + bigB*mass1**2 - bigC*mass1 + bigD
  c2 = mass2**4 - bigA*mass2**3 + bigB*mass2**2 - bigC*mass2 + bigD
  c3 = mass3**4 - bigA*mass3**3 + bigB*mass3**2 - bigC*mass3 + bigD
  c4 = mass4**4 - bigA*mass4**3 + bigB*mass4**2 - bigC*mass4 + bigD

  print *,' c1, c2, c3, c4       ',c1,c2,c3,c4 

! Increase by 0.01%
  d1 = (mass1*1.0001d0)**4 - bigA*(mass1*1.0001d0)**3 + bigB*(mass1*1.0001d0)**2 - bigC*(mass1*1.0001d0) + bigD
  d2 = (mass2*1.0001d0)**4 - bigA*(mass2*1.0001d0)**3 + bigB*(mass2*1.0001d0)**2 - bigC*(mass2*1.0001d0) + bigD
  d3 = (mass3*1.0001d0)**4 - bigA*(mass3*1.0001d0)**3 + bigB*(mass3*1.0001d0)**2 - bigC*(mass3*1.0001d0) + bigD
  d4 = (mass4*1.0001d0)**4 - bigA*(mass4*1.0001d0)**3 + bigB*(mass4*1.0001d0)**2 - bigC*(mass4*1.0001d0) + bigD
! Decrease by 0.01%
  e1 = (mass1/1.0001d0)**4 - bigA*(mass1/1.0001d0)**3 + bigB*(mass1/1.0001d0)**2 - bigC*(mass1/1.0001d0) + bigD
  e2 = (mass2/1.0001d0)**4 - bigA*(mass2/1.0001d0)**3 + bigB*(mass2/1.0001d0)**2 - bigC*(mass2/1.0001d0) + bigD
  e3 = (mass3/1.0001d0)**4 - bigA*(mass3/1.0001d0)**3 + bigB*(mass3/1.0001d0)**2 - bigC*(mass3/1.0001d0) + bigD
  e4 = (mass4/1.0001d0)**4 - bigA*(mass4/1.0001d0)**3 + bigB*(mass4/1.0001d0)**2 - bigC*(mass4/1.0001d0) + bigD

  print *,' d1, d2, d3, d4       ',d1,d2,d3,d4 
  print *,' e1, e2, e3, e4       ',e1,e2,e3,e4 

end subroutine gounaris

subroutine printvector(vectorname,x)
implicit none
integer, parameter :: real_8_30 = selected_real_kind(p=8,r=30)
logical :: lstandard
parameter (lstandard = .false.)
character(len=5) :: vectorname
real(real_8_30) :: x(2)

print *,' '
print *,'Vector ',vectorname

if(lstandard)then
   write(6,1001)x(1)
   write(6,1001)x(2)
else
   write(6,1002)x(1),sqrt(x(1))
   write(6,1002)x(2),sqrt(x(2))
endif

1001 format(20X,1(4X,f14.4))
1002 format(20X,2(4X,f14.4))

end subroutine printvector

subroutine printmatrix(matrixname,X)
implicit none
integer, parameter :: real_8_30 = selected_real_kind(p=8,r=30)
character(len=5) :: matrixname
real(real_8_30) :: X(2,2)
integer :: i,j

print *,' '
print *,'Matrix ',matrixname
write(6,1001)X(1,1),X(1,2)
write(6,1001)X(2,1),X(2,2)

1001 format(20X,2(4X,f14.4))

end subroutine printmatrix

subroutine checkdet(X)
implicit none
integer, parameter :: real_8_30 = selected_real_kind(p=8,r=30)
real(real_8_30) :: X(2,2)
real(real_8_30) :: a,b,c,d,det
include 'lprint.inc'

a = X(1,1)
b = X(1,2)
c = X(2,1)
d = X(2,2)

det = a*d - b*c

if(lprint)print *,'Determinant is ',det

end subroutine checkdet

subroutine charginomasses(M1,M2,mu,tanb)
implicit none
  integer, parameter :: real_8_30 = selected_real_kind(p=8,r=30)
real(real_8_30) :: M1, M2, mu, tanb
real(real_8_30) :: mZ, mW, s2W, cW, sW, beta, cb, sb, c2b, s2b
real(real_8_30) :: msqsum,msqdif,msqdifsq,mass1,mass2
include 'lprint.inc'

! constants
mZ  = 91.1876
mW  = 80.379
s2W = 0.232
sW = sqrt(s2W)
cW = sqrt(1.0-s2W)
beta = atan2(tanb,1.0d0)
cb = cos(beta)
sb = sin(beta)
c2b = cos(2.0*beta)
s2b = sin(2.0*beta)

! Formula here is from Barnett and Haber, Phys. Rev. D 31, 85 (1985). Eqn A4.
msqsum = M2*M2 + mu*mu + 2.0*mW*mW
msqdifsq = (M2**2-mu**2)**2 + 4.0*mW**4*c2b**2 + 4.0*mW*mW*(M2**2 + mu**2 + 2.0*M2*mu*s2b)
msqdif = sqrt(msqdifsq)

mass1 = sqrt(0.5*(msqsum-msqdif))
mass2 = sqrt(0.5*(msqsum+msqdif))

if(lprint)print *,'Chargino masses ',mass1,mass2
! Checked with Mathematica

end subroutine charginomasses

subroutine sortthem(w,iordered)
! Take the mass eigenvalues and sort them by absolute value
! This is a bit of a clumsy implementation but is likely OK for 4 elements
implicit none
  integer, parameter :: real_8_30 = selected_real_kind(p=8,r=30)
real(real_8_30) :: w(4)
integer :: i
integer :: iordered(4)
integer :: imin,imax
integer :: i2,i3,j2,j3
integer :: checksum
real(real_8_30) :: minmass
real(real_8_30) :: maxmass

do i=1,4
   if(i.eq.1)then
      minmass = abs(w(i))
      imin = i
   else
      if(abs(w(i)).lt.minmass)then
         minmass = abs(w(i))
         imin = i
      endif
   endif
enddo
!print *,' Minimum mass found of ',minmass,' for eigenvalue ',imin

do i=1,4
   if(i.eq.1)then
      maxmass = abs(w(i))
      imax = i
   else
      if(abs(w(i)).gt.maxmass)then
         maxmass = abs(w(i))
         imax = i
      endif
   endif
enddo
!print *,' Maximum mass found of ',maxmass,' for eigenvalue ',imax

i2 = 0
i3 = 0
do i=1,4
   if(i.ne.imin.and.i.ne.imax)then
      if(i2.eq.0)then
         i2 = i
      else
         i3 = i
      endif
   endif
enddo

if( abs(w(i2)).lt.abs(w(i3)) )then
! nominal ordering
else
! swap them
  j2 = i2
  j3 = i3
  i2 = j3
  i3 = j2
endif

checksum = imin + i2 + i3 + imax
if(checksum.ne.10)then
! Issues (1+2+3+4 = 10). So should never occur if above logic is fine.
   print *,'checksum wrong in sortthem'
   print *,'STOPPING!'
   stop
endif

!print *,imin,i2,i3,imax

iordered(1) = imin
iordered(2) = i2
iordered(3) = i3
iordered(4) = imax

!do i=1,4
!   print *,' Mass ',i,abs(w(iordered(i)))
!enddo

end subroutine sortthem
