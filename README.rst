performance programming coursework
==================================

This directory contains a simple molecular dynamics (MD) code

The program implements:

    long range inverse square forces between particles.
    F = G * m1*m2 / r**2

    viscosity term
    F = -u V

If 2 particles approach closer than Size we flip the direction of the
interaction force to approximate a collision.

Coordinates are relative to a large central mass and the entire system
is moving relative to the viscous media.

This program was developed as part of a code optimisation course
and is therefore deliberately inefficient.
