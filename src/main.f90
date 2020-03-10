program main
  implicit none

  !include "mkl.fi"

  !
  ! Nbody	    number of particles
  ! pos		    position of the particles
  ! r         distance of particle from central mass
  ! velo		  velocity of the particles
  ! f		      forces acting on each particle
  ! vis		    viscosity coefficient for each particle
  ! mass		  mass of each particle
  ! delta_pos	separation vector for each particle pair
  ! delta_r		separation for each particle pair
  !
  integer, parameter :: Nbody = 4 * 1024, Ndim = 3
  integer :: Xcoord = 1, Ycoord = 2, Zcoord = 3

  ! when distance between two particles is smaller than
  ! delta_threshold the force is flipped
  double precision :: delta_threshold = 1.0

  double precision, dimension(Nbody) :: vis, mass
  double precision, dimension(Nbody, Ndim) :: pos, velo, f
  double precision, dimension(Ndim) :: wind
  integer :: collisions

  ! size of central mass
  double precision :: M_central = 1000d0

  ! interaction strength for pairwise interactions
  double precision :: G = 2d0

  integer :: clock_rate, clock_max
  integer :: tstart, tstop
  integer :: istart, istop

  ! timestep value
  double precision :: dt = 0.02

  ! number of timesteps to use
  integer :: Nstep = 100
  integer, parameter :: Nsave = 5

  character(len = *), parameter :: IO_FMT = "(8E16.8)"
  character(len = 80) :: file_out

  double precision, dimension(Nsave + 1) :: timings

  integer :: start = 1
  integer :: i

  f(:,:) = 0.0

  ! create directory for the output files
  call system("mkdir -p out")

  ! remove output from previous run
  call system("rm -f out/*")

  wind(XCoord) = 0.9
  wind(YCoord) = 0.4
  wind(ZCoord) = 0.0

  collisions=0

#ifdef test
  start = Nsave
  call read_data("data/output.dat400")
#else
  call read_data("data/input.dat")
#endif

  call system_clock(tstart, clock_rate, clock_max)

  do i = start, Nsave

#ifdef test_all
    if(i > 1) then
      write(file_out, "(A14,I3.3)") 'out/output.dat', (i - 1) * Nstep
      call read_data(file_out)
    end if
#endif

    call system_clock(istart, clock_rate, clock_max)
    call evolve()
    call system_clock(istop, clock_rate, clock_max)

    timings(i) = real(istop - istart) / real(clock_rate)

    print *, Nstep, " timesteps took ", timings(i)
    print *, collisions, " collisions"

    write(file_out, "(A14,I3.3)") 'out/output.dat', i * Nstep
    call write_data(file_out)
  end do

  call system_clock(tstop, clock_rate, clock_max)

  timings(Nsave + 1) = real(tstop - tstart) / real(clock_rate)

  print *, Nsave * Nstep - (start - 1) * Nstep, ' timesteps took ', &
    timings(Nsave + 1)

#if !defined test && !defined test_all
  open( unit=1, file="bench.out", action="write", position="append" &
      , status="unknown" )
  write(1, "(6E16.8)") timings
  close(unit=1)
#endif


contains


  subroutine evolve() ! {{{
    double precision, dimension(Ndim) :: fi
    integer :: step, i

    do step = 1, Nstep
      print *, 'timestep ', step
      print *, 'collisions ', collisions

      do i = 1, Nbody
        fi(:) = 0.0
        call compute_pairwise_forces(i, fi)
        call update_particles(i, fi)
      end do
    end do
  end ! }}}


  subroutine compute_pairwise_forces(i, fi)
    integer, intent(in) :: i
    double precision, dimension(Ndim), intent(inout) :: fi

    integer :: j
    double precision :: dd_pos
    double precision, dimension(Ndim) :: d_pos, pairwise_force

    !$omp simd private(d_pos, dd_pos, pairwise_force) &
    !$omp      reduction(+:fi) &
    !$omp      reduction(+:collisions)
    do j = i + 1, Nbody

      ! fused and unrolled
      d_pos(:) = pos(i,:) - pos(j,:)
      dd_pos   = sqrt(sum(d_pos(:) ** 2))

      ! unrolled
      pairwise_force(:) = &
        G * mass(i) * mass(j) * d_pos(:) / (dd_pos ** 3)

      if(dd_pos < delta_threshold) then
        ! unrolled
        pairwise_force(:) = -pairwise_force(:)
        collisions = collisions + 1
      end if

      ! fused and unrolled
      fi(:)  = fi(:)  - pairwise_force(:)
      f(j,:) = f(j,:) + pairwise_force(:)
    end do
    !$omp end simd
  end


  subroutine update_particles(i, fi)
    integer, intent(in) :: i
    double precision, dimension(Ndim), intent(inout) :: fi

    integer :: j
    double precision :: r

    ! vectorized with -vec-threshold0
    r = sqrt(sum(pos(i,:) ** 2)) ** 3

    !dir$ vector aligned
    do j = 1, Ndim
      fi(j) = fi(j) + f(i,j)
      fi(j) = fi(j) - vis(i) * (velo(i,j) + wind(j))
      fi(j) = fi(j) - G * mass(i) * M_central * pos(i,j) / r
      fi(j) = fi(j) / mass(i)

      pos(i,j)  = pos(i,j) + dt * velo(i,j)
      velo(i,j) = velo(i,j) + dt * fi(j)
      f(i,j) = 0.0
    end do
  end


  subroutine read_data(filename) ! {{{
    character(len = *), intent(in) :: filename

    integer :: i

    open(unit=1, file=filename, status='OLD')
    do i = 1, Nbody
      read(1, IO_FMT) mass(i), vis(i), &
        pos (i, Xcoord), pos (i, Ycoord), pos (i, Zcoord), &
        velo(i, Xcoord), velo(i, Ycoord), velo(i, Zcoord)
    end do
    close(unit=1)
  end ! }}}


  subroutine write_data(filename) ! {{{
    character(len = *), intent(in) :: filename

    integer :: i

    open(unit=1, file=filename)
    do i = 1, Nbody
      write(1, IO_FMT) mass(i), vis(i), &
        pos (i, Xcoord), pos (i, Ycoord), pos (i, Zcoord), &
        velo(i, Xcoord), velo(i, Ycoord), velo(i, Zcoord)
    end do
    close(unit=1)
  end ! }}}
end
