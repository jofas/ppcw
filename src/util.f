
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

      FUNCTION force(W,delta,r)
      IMPLICIT NONE
      DOUBLE PRECISION force,W,delta,r
        force=W*delta/(r**3)
        RETURN
      END

