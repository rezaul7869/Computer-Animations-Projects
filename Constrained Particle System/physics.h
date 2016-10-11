/*

  USC/Viterbi/Computer Science
  CONSTRAINED PARTICLE SYSTEM
  Using the -
  "Jello Cube" Assignment 1 starter code

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

void computeAcceleration(struct world * particle, struct point* a);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the particle structure accordingly
void Euler(struct world * particle);
void RK4(struct world * particle);

#endif

