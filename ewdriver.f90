program ewdriver
implicit none

integer, parameter :: real_8_30 = selected_real_kind(p=8,r=30)
real(real_8_30) :: M1, M2, mu, tanb
real(real_8_30) :: M2sq1, musq1, M2sq2, musq2
real(real_8_30) :: M2sq3, musq3, M2sq4, musq4
real(real_8_30) :: mC1, mC2, mW
real(real_8_30) :: M2mutest, mN

integer :: ichoice

mW  = 80.379d0

M1 =   700.0d0
M2 =   500.0d0
mu =   110.0d0
tanb =  10.0d0

include 'ewdriver_inputs.dat'

!mC1 = 200.0d0
!mC2 = 400.0d0
!mC1 = 145.86422024179075d0
!mC2 = 622.95862345364026d0
mC1 = 200.0d0
mC2 = 400.0d0
tanb=10.0d0

! First check whether M2*mu can be negative and still lead to a solution
M2mutest = mW*mW*sin(2.0d0*atan(tanb)) - mC1*mC2
print *,'M2mutest = ',M2mutest

call calcM2mu(mC1, mC2, tanb, M2sq1, musq1, M2sq2, musq2, M2sq3, musq3, M2sq4, musq4)

!Define solution choice
ichoice = 1

print *,'ichoice = ',ichoice

if(ichoice.eq.1)then
   M2 = sqrt(M2sq1)
   mu = sqrt(musq1)
elseif(ichoice.eq.2)then
   M2 = sqrt(M2sq2)
   mu = sqrt(musq2)
elseif(ichoice.eq.3)then
   M2 = sqrt(M2sq3)
   mu = -sqrt(musq3)
elseif(ichoice.eq.4)then
   M2 = sqrt(M2sq4)
   mu = -sqrt(musq4)
else
   print *,'SCREAM no mu,M2 set based on chargino masses'
endif

print *,'M2, mu set to ',M2,mu

M1=0.5d0*M2

!mN = 129.665d0
!mN = -168.8824d0
mN = -180.0d0
call KMinversion(M2,mu,tanb,mN,M1)

call electroweakino(M1,M2,mu,tanb)

call analyze(M1,M2,mu,tanb)

end program ewdriver

subroutine analyze(M1,M2,mu,tanb)

integer, parameter :: real_8_30 = selected_real_kind(p=8,r=30)
real(real_8_30), intent(in) :: M1,M2,mu,tanb
real(real_8_30) :: a,b,c,d,e,f
real(real_8_30) :: mW,smu,s2W,tW,s2b
parameter (mW = 80.379d0)
parameter (s2W = 0.232d0)

s2b = sin(2.0d0*atan(tanb))
print *,'sin2b = ',s2b

! Sign of mu
smu = 1.0d0
if(mu<0.0d0)smu = -1.0d0

! Not sure where this a,b,c,d,e,f thing came from 
! Looks like it is the three mass differences from 
! Han, Kribs, Martin and Menon, 1401.1235, Phys Rev D89 (2014) 075007
! Case 1: M1 >> M2 > abs(mu)
a = 0.5d0*mW**2*(1.0d0 - smu*sin(2.0d0*atan(tanb)))/(M2 + abs(mu))
b = 0.5d0*mW**2*(1.0d0 + smu*sin(2.0d0*atan(tanb)))/(M2 - abs(mu))
c = mW**2*( mu*sin(2.0d0*atan(tanb)) + M2  )/(M2**2 - mu**2)

! Case 2: M1 >> M2 > mu
tW = tan(asin(sqrt(s2W)))
print *,'tW = ',tW
d = 0.5d0*(mW*tW)**2*(1.0d0 + smu*sin(2.0d0*atan(tanb)))/(M1 - abs(mu))
e = 0.5d0*(mW*tW)**2*(1.0d0 - smu*sin(2.0d0*atan(tanb)))/(M1 + abs(mu))
f = (mW*tW)**2*( mu*sin(2.0d0*atan(tanb)) + M1  )/(M1**2 - mu**2)

print *,'Analytic DMs (Case1) Eqs 4-5',a,b,c
print *,'Analytic DMs (Case2) Eqs 6-7',d,e,f

end subroutine analyze
