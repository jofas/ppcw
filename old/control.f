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

C timings to save the benchmark results
      double precision, dimension(Nsave + 1) :: timings

C read the initial data from a file

      wind(XCoord) = 0.9
      wind(YCoord) = 0.4
      wind(ZCoord) = 0.0
      collisions=0
      OPEN(unit=1,file='data/input.dat',status='OLD')
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
C
C Read data back from file, so testing is possible
C
#ifdef test
        IF(j > 1) THEN
          print *, "Reading from last output file"

          WRITE(name,12) 'output.dat', (j-1)*Nstep
          OPEN(unit=1,file=name,status='OLD')
          DO i=1,Nbody
            READ(1,11) mass(i),vis(i),
     $        pos(i,Xcoord),pos(i,Ycoord),pos(i,Zcoord),
     $        velo(i,Xcoord),velo(i,Ycoord),velo(i,Zcoord)
            radius(i)=0.5
          END DO
          CLOSE(unit=1)
        END IF
#endif
        CALL SYSTEM_CLOCK(istart,clock_rate,clock_max)
        CALL evolve(Nstep,dt)
        CALL SYSTEM_CLOCK(istop,clock_rate,clock_max)

        timings(j) = real(istop - istart) / real(clock_rate)

        WRITE(*,*) (Nstep), ' timesteps took ', timings(j)
        write(*,*) collisions, ' collisions'


C write result to a file
        WRITE(name,12) 'output.dat', j*Nstep
        OPEN(unit=1,file=name)
        DO i=1,Nbody
          WRITE(1,11) mass(i),vis(i),
     $      pos(i,Xcoord),pos(i,Ycoord),pos(i,Zcoord),
     $      velo(i,Xcoord),velo(i,Ycoord),velo(i,Zcoord)
        END DO
        CLOSE(unit=1)
      END DO
      CALL SYSTEM_CLOCK(tstop,clock_rate,clock_max)

      timings(Nsave + 1) = real(tstop - tstart) / real(clock_rate)
      WRITE(*,*) (Nsave*Nstep), ' timesteps took ', timings(Nsave + 1)


#ifndef test
      open( unit=1, file="bench_old.out", action="write"
     $    , position="append", status="unknown" )
      write(1, "(6E16.8)") timings
      close(unit=1)
#endif


11    FORMAT(8E16.8)
12    FORMAT(A10,I3.3)

      END

