program main
  !
  ! Control program for the MD update
  !

  implicit none

! This file defines a common block that contains the primary coordinates
! of the particles,
!
!  Nbody	    Number of particles
!  Npair    	Number of particle pairs
!  pos		    Position of the particles
!  r            Distance of particle from central mass
!  velo		    velocity of the particles
!  f		    Forces acting on each particle
!  vis		    viscosity coefficient for each particle
!  mass		    mass of each particle
!  delta_pos	separation vector for each particle pair
!  delta_r		separation for each particle pair
!
      INTEGER Nbody, Npair
      PARAMETER( Nbody=4*1024, Npair=(Nbody*(Nbody-1))/2 )
      INTEGER Xcoord, Ycoord, Zcoord, Ndim
      PARAMETER( Xcoord=1, Ycoord=2, Zcoord=3, Ndim=3 )
      DOUBLE PRECISION r(Nbody)
      DOUBLE PRECISION pos(Nbody,Ndim), velo(Nbody,Ndim)
      DOUBLE PRECISION f(Nbody,Ndim), vis(Nbody), mass(Nbody)
      DOUBLE PRECISION radius(Nbody)
      DOUBLE PRECISION delta_pos(Nbody*Nbody,Ndim)
      DOUBLE PRECISION delta_r(Nbody*Nbody)
      DOUBLE PRECISION wind(Ndim)
      INTEGER collisions

!      COMMON /coord/ f, pos, vis, mass, radius, velo, &
!              delta_pos, delta_r, &
!              r,collisions,wind


!
! Size of central mass.
!
      DOUBLE PRECISION M_central
      PARAMETER( M_central = 1000d0)

!
! Interaction strength for pairwise interactions.
!
      DOUBLE PRECISION G
      PARAMETER( G = 2d0 )
  !INCLUDE 'coord.inc'


  INTEGER i,j
  INTEGER clock_rate, clock_max
  INTEGER tstart, tstop
  INTEGER istart, istop

  !  timestep value
  DOUBLE PRECISION :: dt = 0.02
  !PARAMETER( dt = 0.02 )

  !  number of timesteps to use.
  INTEGER :: Nstep = 100, Nsave = 5
  !    INTEGER Nsave
  !PARAMETER(Nstep=100,Nsave=5)
  CHARACTER(len = 80) :: file_out

  !read the initial data from a file

  wind(XCoord) = 0.9
  wind(YCoord) = 0.4
  wind(ZCoord) = 0.0
  collisions=0

  OPEN(unit=1,file='input.dat',status='OLD')
  DO i=1,Nbody
    READ(1,11) mass(i),vis(i), &
      pos(i,Xcoord),pos(i,Ycoord),pos(i,Zcoord), &
      velo(i,Xcoord),velo(i,Ycoord),velo(i,Zcoord)
    radius(i)=0.5
  END DO
  CLOSE(unit=1)

  !
  ! Run 20 timesteps and time how long it took
  !

  CALL SYSTEM_CLOCK(tstart,clock_rate,clock_max)
  DO j=1,Nsave
    CALL SYSTEM_CLOCK(istart,clock_rate,clock_max)
    CALL evolve(Nstep,dt)
    CALL SYSTEM_CLOCK(istop,clock_rate,clock_max)

    WRITE(*,*) (Nstep), ' timesteps took ', &
      REAL(istop-istart)/REAL(clock_rate)
    write(*,*) collisions, ' collisions'


    ! write result to a file
    WRITE(file_out,12) 'output.dat', j*Nstep
    OPEN(unit=1,file=file_out)
    DO i=1,Nbody
      WRITE(1,11) mass(i),radius(i),vis(i), &
        pos(i,Xcoord),pos(i,Ycoord),pos(i,Zcoord), &
        velo(i,Xcoord),velo(i,Ycoord),velo(i,Zcoord)
    END DO
    CLOSE(unit=1)
  END DO
  CALL SYSTEM_CLOCK(tstop,clock_rate,clock_max)

  WRITE(*,*) (Nsave*Nstep), ' timesteps took ', &
    REAL(tstop-tstart)/REAL(clock_rate)

  11 FORMAT(9E16.8)
  12 FORMAT(A10,I3.3)

contains

  SUBROUTINE evolve( count, dt)
        !include 'coord.inc'
      INTEGER count, step
      DOUBLE PRECISION dt
      DOUBLE PRECISION force
      DOUBLE PRECISION Size
      LOGICAL collided
      INTEGER i,j,k,l

!
! Loop over timesteps.
!
      DO step = 1,count
        write(*,*) 'timestep ',step
        write(*,*) 'collisions ',collisions

! set the viscosity term in the force calculation
        DO j=1,Ndim
          CALL visc_force(Nbody,f(1,j),vis,velo(1,j))
        END DO
! add the wind term in the force calculation
        DO j=1,Ndim
          CALL wind_force(Nbody,f(1,j),vis,wind(j))
        END DO

! calculate distance from central mass
        DO k=1,Nbody
          r(k) = 0.0
        END DO
        DO i=1,Ndim
          call add_norm(Nbody,r,pos(1,i))
        END DO
        DO k=1,Nbody
          r(k) = SQRT(r(k))
        END DO

! calculate central force

        DO i=1,Nbody
          DO l=1,Ndim
                f(i,l) = f(i,l) - &
                   calc_force(G*mass(i)*M_central,pos(i,l),r(i))
          END DO
        END DO

! calculate pairwise separation of particles
        k = 1
        DO i=1,Nbody
          DO j=i+1,Nbody
            DO l=1,Ndim
              delta_pos(k,l) = pos(i,l) - pos(j,l)
            END DO

            k = k + 1
          END DO
        END DO

! calculate norm of seperation vector
        DO k=1,Npair
          delta_r(k) = 0.0
        END DO
        DO i=1,Ndim
          call add_norm(Npair,delta_r,delta_pos(1,i))
        END DO
        DO k=1,Npair
          delta_r(k) = SQRT(delta_r(k))
        END DO

!
! add pairwise forces.
!
        k = 1
        DO i=1,Nbody
          DO j=i+1,Nbody
            collided=.false.
            Size = radius(i) + radius(j)
            DO l=1,Ndim
!  flip force if close in
              IF( delta_r(k) .GE. Size ) THEN
                f(i,l) = f(i,l) - &
                   calc_force(G*mass(i)*mass(j),delta_pos(k,l),delta_r(k))
                f(j,l) = f(j,l) + &
                   calc_force(G*mass(i)*mass(j),delta_pos(k,l),delta_r(k))
              ELSE
                f(i,l) = f(i,l) + &
                    calc_force(G*mass(i)*mass(j),delta_pos(k,l),delta_r(k))
                f(j,l) = f(j,l) - &
                    calc_force(G*mass(i)*mass(j),delta_pos(k,l),delta_r(k))
                collided=.true.
              END IF
            END DO
            IF( collided )THEN
              collisions = collisions + 1
            END IF
            k = k + 1
          END DO
        END DO

! update positions
        DO i=1,Nbody
          DO j=1,Ndim
            pos(i,j) = pos(i,j) + dt * velo(i,j)
          END DO
        END DO

! update velocities
        DO i=1,Nbody
          DO j=1,Ndim
            velo(i,j) = velo(i,j) + dt * (f(i,j)/mass(i))
          END DO
        END DO


      END DO

      END


      SUBROUTINE visc_force(N,f,vis,velo)
      IMPLICIT NONE
      INTEGER N, i
      DOUBLE PRECISION f(N),vis(N),velo(N)
          DO i=1,N
            f(i) = -vis(i) * velo(i)
          END DO
      END

      SUBROUTINE wind_force(N,f,vis,velo)
      IMPLICIT NONE
      INTEGER N, i
      DOUBLE PRECISION f(N),vis(N),velo
          DO i=1,N
            f(i) = f(i) -vis(i) * velo
          END DO
      END

      SUBROUTINE add_norm(N,r,delta)
      IMPLICIT NONE
      INTEGER N, k
      DOUBLE PRECISION r(N),delta(N)
        DO k=1,N
          r(k) = r(k) + (delta(k) * delta(k))
        END DO
      END

      double precision function calc_force(W,delta,r)
      IMPLICIT NONE
      DOUBLE PRECISION :: W,delta,r
        calc_force=W*delta/(r**3)
        RETURN
      END
END
