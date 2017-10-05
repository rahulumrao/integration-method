PROGRAM Runga_Kutta_second_order 
!Runga-Kutta second_order method produces more accurate data than Euler`s method by
!calculating the slope at the midpoint of the interval based oïf Euler`s approx-
!imation. This midpoint slope is then used to make a better extrapolation of the
!endpoint of the interval. It is interesting to see how energy is conserved between
!the two methods. Obviously, a simple harmonic oscillator is a conservative sys-
!tem, therefore, we should not get an increase or decrease of energy throughout
!it`s time-development. Energy is conserved in RK2 method.
IMPLICIT NONE
REAL*8  :: x,v,v1,x1,k1
REAL*8  :: delT,h,t,k,KE,PE,m
INTEGER :: i,n
!x(1) = 1.0
         x = 0.1D0 ; v = 0.0D0                   !  initial velocity and positions
         k = 1.0D0 ; m = 1.0D0                   !  spring constant and mass of particle
delT = 0.01                                     !  time step
n = 10000                                        !  number of step
WRITE(6,*)'Enter the number of time steps'
READ(5,*) n
OPEN(1,file='rk_2_posi_vel')
OPEN(2,file='rk_2_energies')
OPEN(3,file='rk_2_trej.xyz')
!*****************************************************************************************
     DO i =1,n
  k1 = DelT*(v + (-k/m)*x*(delT/2))              !k1= dt(v(t) + (-k/m)x(t)dt/2)
   x1 = x + k1                                   !x(t+dt) = x(t) + k2
  k1 = DelT*(-k/m)*(x + v*(delT/2))              !k1= dt(-k/m)(x(t) + v(t)dt/2)
   v1 = v + k1                                   !v(t+dt) = v(t) + k2
v = v1                                           !storing new velocity
x = x1                                           !storing new position
t = delT*i
!*******************************************************************************************
       PE = 0.5*(k*x**2)                         !potential energy
       KE = 0.5*(m*v**2)                         !kinetic energy
      !.... writing trajectory in VMD................
WRITE(3,*) 1       ! 1 indicate that we have one particle
WRITE(3,*)         !here you can write comment or else leave it blank
WRITE(3,*) "m",x,0.d0,0.d0   ! here you can write the position coordinates
!******************************************************************************************
!PRINT*,t,x,v
WRITE(2,*)t,KE,PE,KE+PE
WRITE(1,*)t,x,v
END DO
PRINT*,'THANKS..........'
CLOSE(1)
CLOSE(2)
CLOSE(3)
END PROGRAM Runga_Kutta_second_order

