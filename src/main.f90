program main
  implicit none

  !
  ! Nbody	    number of particles
  ! Npair    	number of particle pairs
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

  double precision, dimension(Nbody) :: vis, mass
  double precision, dimension(Nbody, Ndim) :: f, pos, velo
  !double precision, dimension(Ndim, Nbody) :: f, pos, velo
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

  wind(1) = 0.9
  wind(2) = 0.4
  wind(3) = 0.0

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
    double precision, dimension(Ndim) :: d_pos, force, fi
    double precision :: dd_pos, r
    integer :: step, i, j

    do step = 1, Nstep
      print *, 'timestep ', step
      print *, 'collisions ', collisions

      do i = 1, Nbody
        fi(:) = 0.0

        !$omp simd private(d_pos, dd_pos, force) &
        !$omp   reduction(+:collisions) reduction(+:fi)
        !dir$ vector aligned
        do j = i+1, Nbody
          d_pos(:) = pos(i,:) - pos(j,:)
          !d_pos(:) = pos(:,i) - pos(:,j)
          dd_pos   = sqrt(sum(d_pos(:) ** 2))

          force(:) = -G * mass(i) * mass(j) * d_pos(:) / (dd_pos ** 3)

          if (dd_pos < 1.0) then
            force(:) = -force(:)
            collisions = collisions + 1
          end if

          fi(:) = fi(:) + force(:)
          f(j,:) = f(j,:) - force(:)
          !f(:,j) = f(:,j) - force(:)
        end do

        r = sqrt(sum(pos(i,:) ** 2)) ** 3
        !r = sqrt(sum(pos(:,i) ** 2)) ** 3

        fi(:) = fi(:) + f(i,:) - vis(i) * (velo(i,:) + wind(:)) &
                 - G * mass(i) * M_central * pos(i,:) / r
        !fi(:) = fi(:) + f(:,i) - vis(i) * (velo(:,i) + wind(:)) &
        !         - G * mass(i) * M_central * pos(:,i) / r

        pos(i,:) = pos(i,:) + dt * velo(i,:)
        !pos(:,i) = pos(:,i) + dt * velo(:,i)

        velo(i,:) = velo(i,:) + dt * fi(:) / mass(i)
        !velo(:,i) = velo(:,i) + dt * fi(:) / mass(i)

        f(i,:) = 0.0
        !f(:,i) = 0.0
      end do

    end do
  end ! }}}


  subroutine read_data(filename) ! {{{
    character(len = *), intent(in) :: filename

    integer :: i

    open(unit=1, file=filename, status='OLD')
    do i = 1, Nbody
      read(1, IO_FMT) mass(i), vis(i), &
        pos (i, 1), pos (i, 2), pos (i, 3), &
        velo(i, 1), velo(i, 2), velo(i, 3)
        !pos (1, i), pos (2, i), pos (3, i), &
        !velo(1, i), velo(2, i), velo(3, i)
    end do
    close(unit=1)
  end ! }}}


  subroutine write_data(filename) ! {{{
    character(len = *), intent(in) :: filename

    integer :: i

    open(unit=1, file=filename)
    do i = 1, Nbody
      write(1, IO_FMT) mass(i), vis(i), &
        pos (i, 1), pos (i, 2), pos (i, 3), &
        velo(i, 1), velo(i, 2), velo(i, 3)
        !pos (1, i), pos (2, i), pos (3, i), &
        !velo(1, i), velo(2, i), velo(3, i)
    end do
    close(unit=1)
  end ! }}}


end
