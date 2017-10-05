PROGRAM rk4
!These methods were developed around 1900 by the German mathematicians C. Runge and M. W. Kutta.
!The most widely known member of the Runge–Kutta family is generally referred to as "RK4",&
!& "classical Runge–Kutta method" or simply as "the Runge–Kutta method".
!The RK4 method is a fourth-order method, meaning that the local truncation error is on the order of O(h^5)
IMPLICIT NONE
REAL*8  :: x,v,v1,x1,v2,x2,v3,x3,v4,x4
REAL*8  :: delT,h,t,k,KE,PE,m
INTEGER :: i,n
!x(1) = 1.0
x = 0.1D0 ; v = 0.0D0                      !  initial velocity and positions
k = 1.0D0 ; m = 1.0D0                      !  spring constant and mass of particle
delT = 0.01D0                             !  time step
n = 10000                                  !  number of step
OPEN(1,file='rk_4_posi_vel')
OPEN(2,file='rk_4_energies')
OPEN(3,file='rk_4_trej.xyz')
!**********************************************************************************
     DO i =1,n
        x1 = v*delT
           v1 = delT*((-k/m)*x)
        x2 = delT*(v+v1/2)
           v2 = delT*((-k/m)*(x+x1/2))
        x3 = delT*(v+v2/2)
           v3 = delT*((-k/m)*(x+x2/2))
        x4 = delT*(v+v3)
           v4 = delT*((-k/m)*(x+x3))

x = x+(x1+2*x2+2*x3+x4)/6
v = v+(v1+2*v2+2*v3+v4)/6 
!***********************************************************************************
       PE = 0.5*(k*x**2)                      !  potential energy
       KE = 0.5*(m*v**2)                      !  kinetic energy
       t = i*delT
!***********************************************************************************
WRITE(1,*)t,x,v,PE    
WRITE(2,*)t,KE,PE,KE+PE     
!WRITE(6,*)t,x,v     
!.... writing trajectory in VMD................
      WRITE(3,*) 1                             ! 1 indicate that we have one particle
      WRITE(3,*)                               !here you can write comment or else leave it blank
      WRITE(3,*) "m",x,0.d0,0.d0              ! here you can write the position coordinates
!***********************************************************************************
     END DO
WRITE(6,*)'THANKS.........'     
WRITE(6,*)'Data Stored in ....     rk_4_posi_vel,      rk_4_energies'
WRITE(6,*)    
END PROGRAM rk4
