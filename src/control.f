C
C Control program for the MD update
C

      PROGRAM MD
      IMPLICIT NONE
      INCLUDE 'coord.inc'
      INTEGER i,j
      INTEGER clock_rate, clock_max
      INTEGER tstart, tstop
      INTEGER istart, istop

C  timestep value
      DOUBLE PRECISION dt
      PARAMETER( dt = 0.02 )

C  number of timesteps to use.
      INTEGER Nstep
      INTEGER Nsave
      PARAMETER(Nstep=100,Nsave=5)
      CHARACTER*80 name

C read the initial data from a file

      wind(XCoord) = 0.9
      wind(YCoord) = 0.4
      wind(ZCoord) = 0.0
      collisions=0
      OPEN(unit=1,file='input.dat',status='OLD')
      DO i=1,Nbody
        READ(1,11) mass(i),vis(i),
     $    pos(i,Xcoord),pos(i,Ycoord),pos(i,Zcoord),
     $    velo(i,Xcoord),velo(i,Ycoord),velo(i,Zcoord)
        radius(i)=0.5
      END DO
      CLOSE(unit=1)

C
C Run 20 timesteps and time how long it took
C

      CALL SYSTEM_CLOCK(tstart,clock_rate,clock_max)
      DO j=1,Nsave
        CALL SYSTEM_CLOCK(istart,clock_rate,clock_max)
        CALL evolve(Nstep,dt)
        CALL SYSTEM_CLOCK(istop,clock_rate,clock_max) 

        WRITE(*,*) (Nstep), ' timesteps took ',
     $       REAL(istop-istart)/REAL(clock_rate)
        write(*,*) collisions, ' collisions'


C write result to a file
        WRITE(name,12) 'output.dat', j*Nstep
        OPEN(unit=1,file=name)
        DO i=1,Nbody
          WRITE(1,11) mass(i),radius(i),vis(i),
     $      pos(i,Xcoord),pos(i,Ycoord),pos(i,Zcoord),
     $      velo(i,Xcoord),velo(i,Ycoord),velo(i,Zcoord)
        END DO
        CLOSE(unit=1)
      END DO
      CALL SYSTEM_CLOCK(tstop,clock_rate,clock_max) 

      WRITE(*,*) (Nsave*Nstep), ' timesteps took ',
     $       REAL(tstop-tstart)/REAL(clock_rate)

11    FORMAT(9E16.8)
12    FORMAT(A10,I3.3)

      END

