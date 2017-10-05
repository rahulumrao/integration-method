PROGRAM verlet
!VELOCITY VERLET ALGORITHM FOR THE KINETIC ENERGY
!AND POTENTIAL ENERGY OF SIMPLE HARMONIC OSCILLATOR
INTEGER :: i,n
REAL*8  :: a,v,x
REAL(8) :: k,m,delT,KE,PE
REAL :: Kb,Temp
!F = ma THEN a = F/m & F = -kx
!LET MASS = 1amu ; k = -1.0 ; delT = 0.1
!WRITE(6,*)'time step'
!READ(5,*)delT
OPEN(1,file="verlet_posi_vel")
OPEN(2,file="verlet_energy")
K = 1.0D0 ; n = 10000 
m = 1.0D0 ; delT = 0.01d0 
x = 1.0D0 ; v = 1.0D0
!Kb = 1.3806*10E-23
Kb = 1
!**************************************************************************************************
   DO i = 1,n ! n = no of MD steps 
       A = -k*x/m ! ! Force = -k*x 
       x1 = x + v*delT+(0.5*delT**2)*(A)       !x(t+∆t) = x(t) + v(t)∆t +0.5*a(t)∆t**2
       A1 = -k*x1/m                                     !VELOCITY VERLET ALGORITHM
       v1 = v + (0.5*(A + A1)*delT)           !v(t+∆t) = v(t) + 0.5*(a(t) + a(t+∆t))∆t
       t = delT*i
! 0.50*m*v1*v1 = KE related to temperature 
x = x1
v = v1
!WRITE(6,*)t,v(i),x(i)
!****************************************************************************************************
! THE KINETIC ENERGY , POTENTIAL ENERGY AND TOTAL ENERGY OF
! THE SYSTEM ; TOTAL ENERGY IS CONSTANT
!****************************************************************************************************
       KE = 0.5*(m*(v)**2)
       PE = 0.5*(K*(x)**2)
!****************************************************************************************************
       Temp = (2.0/3.0)*(KE)/Kb
PRINT*, Temp
!****************************************************************************************************
WRITE(1,*)t,x,v
WRITE(2,*)t,KE,PE,KE+PE,Temp
   ENDDO
WRITE(6,*)"Energy of the system are; ",  "KE =",KE,  "PE =",PE,  "Total =", KE + PE
CLOSE(1)
CLOSE(2)
END PROGRAM verlet
