program ewdriver
implicit none

real :: M1, M2, mu, tanb
real :: a,b,c,d,e,f
real :: mW
real :: smu
real :: s2W
real :: tW
real :: s2b
parameter (mW = 80.379)
parameter (s2W = 0.232)

!M1 = 150.0
!M2 = 250.0
!mu = 150.0
!tanb = 10.0

!Wagner
!M1 =   63.5
!M2 = -172.0
!mu = -300.0
!tanb = 20.0

!Choi_0202039
!M1 = 100.5
!M2 = 190.8
!mu = 365.1
!tanb = 10.0

M1 =   700.0
M2 =   500.0
mu =   110.0
tanb =  10.0

include 'ewdriver_inputs.dat'

call electroweakino(M1,M2,mu,tanb)

s2b = sin(2.0*atan(tanb))
print *,'sin2b = ',s2b

! Sign of mu
smu = 1.0
if(mu<0.0)smu = -1.0

! Not sure where this a,b,c,d,e,f thing came from 
! Looks like it is the three mass differences from 
! Han, Kribs, Martin and Menon, 1401.1235, Phys Rev D89 (2014) 075007
! Case 1: M1 >> M2 > abs(mu)
a = 0.5*mW**2*(1.0 - smu*sin(2.0*atan(tanb)))/(M2 + abs(mu))
b = 0.5*mW**2*(1.0 + smu*sin(2.0*atan(tanb)))/(M2 - abs(mu))
c = mW**2*( mu*sin(2.0*atan(tanb)) + M2  )/(M2**2 - mu**2)

! Case 2: M1 >> M2 > mu
tW = tan(asin(sqrt(s2w)))
print *,'tW = ',tW
d = 0.5*(mW*tW)**2*(1.0 + smu*sin(2.0*atan(tanb)))/(M1 - abs(mu))
e = 0.5*(mW*tW)**2*(1.0 - smu*sin(2.0*atan(tanb)))/(M1 + abs(mu))
f = (mW*tW)**2*( mu*sin(2.0*atan(tanb)) + M1  )/(M1**2 - mu**2)

print *,'Analytic DMs (Case1) Eqs 4-5',a,b,c
print *,'Analytic DMs (Case2) Eqs 6-7',d,e,f

end program ewdriver
