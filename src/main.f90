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
  integer :: Npair = (Nbody * (Nbody - 1)) / 2
  integer :: Xcoord = 1, Ycoord = 2, Zcoord = 3

  double precision, dimension(Nbody) :: r, vis, mass, radius
  double precision, dimension(Nbody, Ndim) :: f, pos, velo
  double precision, dimension(Nbody * Nbody, Ndim) :: delta_pos
  double precision, dimension(Nbody * Nbody) :: delta_r
  double precision, dimension(Ndim) :: wind
  integer :: collisions

  ! size of central mass
  double precision :: M_central = 1000d0

  ! interaction strength for pairwise interactions
  double precision :: G = 2d0


  integer :: i, j
  integer :: clock_rate, clock_max
  integer :: tstart, tstop
  integer :: istart, istop

  ! timestep value
  double precision :: dt = 0.02

  ! number of timesteps to use
  integer :: Nstep = 100, Nsave = 5

  character(len = 80) :: file_out

  wind(XCoord) = 0.9
  wind(YCoord) = 0.4
  wind(ZCoord) = 0.0

  collisions=0

  !read the initial data from a file
  open(unit=1, file='input.dat', status='OLD')
  do i = 1, Nbody

    read(1, 11) mass(i), vis(i), &
      pos (i, Xcoord), pos (i, Ycoord), pos (i, Zcoord), &
      velo(i, Xcoord), velo(i, Ycoord), velo(i, Zcoord)

    ! TODO: radius to constant
    radius(i) = 0.5
  end do
  close(unit=1)


  call system_clock(tstart, clock_rate, clock_max)

  do j = 1, Nsave
    call system_clock(istart, clock_rate, clock_max)
    call evolve()
    call system_clock(istop, clock_rate, clock_max)

    print *, Nstep, " timesteps took ", &
      real(istop - istart) / real(clock_rate)
    print *, collisions, " collisions"

    ! write result to a file
    ! TODO: why format? Try without
    write(file_out, 12) 'output.dat', j * Nstep
    open(unit=1, file=file_out)
    do i = 1, Nbody
      ! TODO: don't write radius to file, unnecessary
      write(1,11) mass(i), radius(i), vis(i), &
        pos (i, Xcoord), pos (i, Ycoord), pos (i, Zcoord), &
        velo(i, Xcoord), velo(i, Ycoord), velo(i, Zcoord)
    end do
    close(unit=1)
  end do

  call system_clock(tstop, clock_rate, clock_max)

  print *, Nsave * Nstep, ' timesteps took ', &
    real(tstop - tstart) / real(clock_rate)

  11 FORMAT(9E16.8)
  12 FORMAT(A10,I3.3)

contains


  subroutine evolve()
    double precision :: force, size
    logical :: collided
    integer :: step, i, j, k, l

    do step = 1, Nstep
      print *, 'timestep ', step
      print *, 'collisions ', collisions

      ! set the viscosity term in the force calculation
      do j = 1, Ndim
        call visc_force(Nbody, f(1,j), vis, velo(1,j))
      end do

      ! add the wind term in the force calculation
      do j = 1, Ndim
        call wind_force(Nbody, f(1,j), vis, wind(j))
      end do

      ! calculate distance from central mass
      do k = 1, Nbody
        r(k) = 0.0
      end do
      do i = 1, Ndim
        call add_norm(Nbody, r, pos(1,i))
      end do
      do k = 1, Nbody
        r(k) = sqrt(r(k))
      end do

      ! calculate central force
      do i = 1, Nbody
        do l = 1, Ndim
          f(i,l) = f(i,l) &
            - calc_force(G * mass(i) * M_central, pos(i,l), r(i))
        end do
      end do

      ! calculate pairwise separation of particles
      k = 1
      do i = 1, Nbody
        do j = i+1, Nbody
          do l = 1, Ndim
            delta_pos(k,l) = pos(i,l) - pos(j,l)
          end do
          k = k + 1
        end do
      end do

      ! calculate norm of seperation vector
      do k = 1, Npair
        delta_r(k) = 0.0
      end do
      do i = 1, Ndim
        call add_norm(Npair, delta_r, delta_pos(1,i))
      end do
      do k = 1, Npair
        delta_r(k) = SQRT(delta_r(k))
      end do

      ! add pairwise forces
      k = 1
      do i = 1, Nbody
        do j = i+1, Nbody
          collided = .false.
          Size = radius(i) + radius(j)
          do l = 1, Ndim

            !  flip force if close in
            if (delta_r(k) .GE. Size) then
              f(i,l) = f(i,l) - calc_force( G * mass(i) * mass(j) &
                                          , delta_pos(k,l) &
                                          , delta_r(k) )
              f(j,l) = f(j,l) + calc_force( G * mass(i) * mass(j) &
                                          , delta_pos(k,l) &
                                          , delta_r(k) )
            else
              f(i,l) = f(i,l) + calc_force( G * mass(i) * mass(j) &
                                          , delta_pos(k,l) &
                                          , delta_r(k) )
              f(j,l) = f(j,l) - calc_force( G * mass(i) * mass(j) &
                                          , delta_pos(k,l) &
                                          , delta_r(k) )
              collided = .true.
            end if

          end do

          if (collided) then
            collisions = collisions + 1
          end if
          k = k + 1
        end do
      end do

      ! update positions
      do i = 1, Nbody
        do j = 1, Ndim
          pos(i,j) = pos(i,j) + dt * velo(i,j)
        end do
      end do

      ! update velocities
      do i = 1, Nbody
        do j = 1, Ndim
          velo(i,j) = velo(i,j) + dt * (f(i,j) / mass(i))
        end do
      end do

    end do
  end


  subroutine visc_force(N, f, vis, velo)
    integer, intent(in) :: N
    double precision, dimension(N), intent(inout) :: f
    double precision, dimension(N), intent(in) :: vis, velo

    integer :: i

    do i = 1, N
      f(i) = -vis(i) * velo(i)
    end do
  end


  subroutine wind_force(N, f, vis, velo)
    integer, intent(in) ::  N
    double precision, dimension(N), intent(inout) :: f
    double precision, dimension(N), intent(in) :: vis
    double precision, intent(in) :: velo

    integer :: i

    do i = 1, N
      f(i) = f(i) - vis(i) * velo
    end do
  end


  subroutine add_norm(N, r, delta)
    integer, intent(in) :: N
    double precision, intent(inout) :: r(N)
    double precision, intent(in) :: delta(N)

    integer :: i

    do i = 1, N
      r(i) = r(i) + delta(i) * delta(i)
    end do
  end


  double precision function calc_force(W, delta, r)
    double precision, intent(in) :: W, delta, r

    calc_force = W * delta / (r**3)
  end
end
