program orbit
implicit none

real(8), parameter :: pi =  3.1415926535897932384626433832795028841971
real(8), parameter :: g = 9.80665
real(8), parameter :: mu_0 = pi*4.0d-10 ! kg * km /(s^2 A^2)
real(8), parameter :: Rn = 24765.0 ! km (Radius of Neptune)
real(8) :: t_step, t_final, t	! time
real(8), dimension(6) :: state, temp
real(8), dimension(4,6) :: k ! for runge kutta iterations
real(8), dimension(3) :: accel, b, mag, r_mag
integer :: i, cycles, n
! n is used for runge kutta iterations
real(8) :: mass, grav_const, ke, pe, me, p
! state holds posistion and velocity for 1
! temp holds temp position and velocity for 1
! k# holds velocity and accel for 1 
! p is atmospheric density
! OLD units of AU for ditances, yr for time
! G*M = 4pi**2

! NEW UNITS: USE Length:KM	Time:SEC	Mass:KG
mag = (/0.0d0,0.0d0,0.0d0/) ! magnetic moment
r_mag = (/0.0d0,0.0d0,0.0d0/) ! magnetic moment location
grav_const = 6834483.2d0 ! km^3/s^2
mass = 60.0d0
state = (/231796.0d0,0.0d0,0.0d0,0.0d0,5.43d0,0.0d0/)
!state = (/1.0d0,0.0d0,0.0d0,0.0d0,2.0d0,-1.0d0/)
! state 1-3 is for position, 4-6 is for velocity

!init
temp = 0.0d0
k = 0.0d0
t = 0.0d0
i = 0
b=0.0d0

write(*,*) 'Input length of time step:'
read(*,*) t_step

write(*,*) 'Input time to run for:'
read(*,*) t_final

cycles = ceiling(t_final/t_step)+1

ke = 0.5d0 * mass * dot_product(state(4:6),state(4:6))
pe = (-4.0d0*(pi**2)*mass)/norm2(state(1:3))
me = ke + pe

open(unit=101, file='pos_1.dat')
open(unit=111, file='pos_1_2d.dat')
open(unit=121, file='ke.dat')
open(unit=122, file='pe.dat')
open(unit=123, file='me.dat')

write(101,*) state(1), state(2), state(3)
write(111,*) state(1), state(2)
write(121,*) t, ke
write(122,*) t, pe
write(123,*) t, me



do i = 2,cycles
	t = t + t_step

	temp = state
	do n = 1,4
		! b = (mu_0/4.0d0*pi)*(()-(mag/(temp(1:3)-r_mag)**3))
		p
		accel = (-grav_const*temp(1:3))/(norm2(temp(1:3))**3)
		!THIS LINE HAS THE FORCE ACCEL
		k(n,1:3) = temp(4:6)
		k(n,4:6) = accel
		if (n<=2) then
			temp = state + 0.5d0*t_step*k(n,:)
		else if (n==3) then
			temp = state + t_step*k(n,:)
		else
			continue
		end if
	end do

	state = state + (t_step/6.0d0)*(k(1,:)+(2.0d0*k(2,:))+(2.0d0*k(3,:))+k(4,:))
	
	ke = 0.5d0 * mass * dot_product(state(4:6),state(4:6))
	pe = (-4.0d0*(pi**2)*mass)/norm2(state(1:3))
	me = ke + pe
	
	write(101,*) state(1), state(2), state(3)
	write(111,*) state(1), state(2)
	write(121,*) t, ke
	write(122,*) t, pe
	write(123,*) t, me
	write(*,*) t
end do

close(101)
close(111)
close(121)
close(122)
close(123)
end program orbit

