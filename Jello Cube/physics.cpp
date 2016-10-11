/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"


//Calculate elastic spring force between particle at position A and B
void calculateElasticForce(double kElastic, point PositionA, point PositionB, double RestLength, point* ElasticForce)
{
	point PosVector;
	double length;

	pDIFFERENCE(PositionA, PositionB, PosVector);
	pNORMALIZE(PosVector);

	double constTemp = (-kElastic) * (length - RestLength);
	pMULTIPLY(PosVector, constTemp,(*ElasticForce));
}


//Calculate damping spring force between particle at position A and B
void calculateDampingForce(double dElastic, point PositionA, point PositionB, point VelocityA, point VelocityB, point* DampingForce)
{
	point PosVector;
	point RelVelocity;
	double length;

	pDIFFERENCE(PositionA, PositionB, PosVector);
	pDIFFERENCE(VelocityA, VelocityB, RelVelocity);
	pNORMALIZE(PosVector);
	
	double tempVelocityDotPos =  -dElastic * ((RelVelocity.x * PosVector.x) + (RelVelocity.y * PosVector.y) + (RelVelocity.z * PosVector.z));
	pMULTIPLY(PosVector, tempVelocityDotPos, (*DampingForce));
}

// Check and calculate elastic and damping force for the springs
void AddForce(struct world* jello, int i, int j, int k, int di, int dj, int dk, point* totalForce)
{
	int ip = i + di;
    int jp = j + dj;
	int kp = k + dk;
    if(!((ip>7) || (ip<0) || (jp>7) || (jp<0) ||(kp>7) || (kp<0)))
	{
		point ElasticForce = {0, 0, 0};
		point DampingForce = {0, 0, 0};
		
		
		double restLength = sqrt(double(di * di + dj * dj + dk * dk)) / 7;
		
		calculateElasticForce(jello->kElastic, jello->p[i][j][k], jello->p[ip][jp][kp], restLength , &ElasticForce);
		calculateDampingForce(jello->dElastic, jello->p[i][j][k], jello->p[ip][jp][kp], jello->v[i][j][k], jello->v[ip][jp][kp], &DampingForce);

		totalForce->x += (ElasticForce.x + DampingForce.x);
		totalForce->y += (ElasticForce.y + DampingForce.y);
		totalForce->z += (ElasticForce.z + DampingForce.z);
	}
}

// Add interpolated force from the force field to totalForce
void AddInterpolatedForce(struct world* jello, int i, int j, int k, double h, point* totalForce)
{
	point q;
	int qi, qj, qk;
	q.x = qi = floor((jello->p[i][j][k].x + 2.0) / h);
	q.y = qj = floor((jello->p[i][j][k].y + 2.0) / h);
	q.z = qk = floor((jello->p[i][j][k].z + 2.0) / h);

				
	point temp1, temp2;
	point temp3 = {jello->p[i][j][k].x + 2.0, jello->p[i][j][k].y + 2.0, jello->p[i][j][k].z + 2.0};
	pMULTIPLY(q, h, temp1);
	pDIFFERENCE(temp3, temp1, temp2);
	double alpha	= temp2.x / h;
	double beta		= temp2.y / h;
	double gamma	= temp2.z / h;
	if(qi < 0 || qj < 0 || qk <0 || qi > (jello->resolution -1) || qj > (jello->resolution -1) || qk > (jello->resolution -1)) return;

	// Calculate the index in one dimensional array of jello->forcefield
	int index = qi * jello->resolution * jello->resolution + qj * jello->resolution + qk;
	point neighbourForceComp[8] = {0, 0, 0};

	// Perform trilinear interpolation

	pMULTIPLY(jello->forceField[index],((1 - alpha) * (1 - beta) * (1 - gamma)),neighbourForceComp[0]);
	if(qk != (jello->resolution -1)) 
		{pMULTIPLY(jello->forceField[index+1],((alpha) * (1 - beta) * (1 - gamma)),neighbourForceComp[1]);}
	if(qj != (jello->resolution -1)) 
		{pMULTIPLY(jello->forceField[index+jello->resolution],((1 - alpha) * (beta) * (1 - gamma)),neighbourForceComp[2]);}
	if((qk != (jello->resolution -1)) && (qj != (jello->resolution -1))) 
		{pMULTIPLY(jello->forceField[index+jello->resolution+1],((alpha) * (beta) * (1 - gamma)),neighbourForceComp[3]);}
	if(qi != (jello->resolution -1)) 
		{pMULTIPLY(jello->forceField[index+ (jello->resolution * jello->resolution)],((1 - alpha) * (1 - beta) * (gamma)),neighbourForceComp[4]);}
	if((qk != (jello->resolution -1)) && (qi != (jello->resolution -1)))
		{pMULTIPLY(jello->forceField[index+ (jello->resolution * jello->resolution) + 1],((alpha) * (1 - beta) * (gamma)),neighbourForceComp[5]);}
	if((qj != (jello->resolution -1)) && (qi != (jello->resolution -1)))
		{pMULTIPLY(jello->forceField[index+ (jello->resolution * jello->resolution)+ jello->resolution],((1 - alpha) * (beta) * (gamma)),neighbourForceComp[6]);}
	if((qi != (jello->resolution -1)) && (qj != (jello->resolution -1)) && (qk != (jello->resolution -1)))
		{pMULTIPLY(jello->forceField[index+ (jello->resolution * jello->resolution)+jello->resolution+1],((alpha) * (beta) * (gamma)),neighbourForceComp[7]);}

	for(int l = 0; l < 8; l++)
	{
		totalForce->x += neighbourForceComp[l].x;
		totalForce->y += neighbourForceComp[l].y;
		totalForce->z += neighbourForceComp[l].z;
	}
}

// Add collision Force. This function uses the spring elasticity of collision
void AddCollisionForce(struct world* jello, int i, int j, int k, point CollisionContact, point* TotalForce)
{
	point PenaltyForce = {0,0,0}, DampingPenalty = {0,0,0};
	point CollisionVelocity = {0,0,0};
	calculateElasticForce(jello->kCollision, jello->p[i][j][k], CollisionContact, 0, &PenaltyForce);
	calculateDampingForce(jello->dCollision, jello->p[i][j][k], CollisionContact, jello->v[i][j][k], CollisionVelocity, &DampingPenalty); 
	TotalForce->x += (PenaltyForce.x + DampingPenalty.x);
	TotalForce->y += (PenaltyForce.y + DampingPenalty.y);
	TotalForce->z += (PenaltyForce.z + DampingPenalty.z);
}

/* Computes acceleration to every control point of the jello cube, 
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
	double h = 1.0;
	if(jello->resolution > 1)
		h =  (4.0 / (jello->resolution - 1));
	
	//jello->dt = 0.001;
	for(int i = 0; i < 8; i++)
	{
		for(int j = 0; j < 8; j++)
		{
			for(int k = 0; k < 8; k++)
			{
				point TotalForce = { 0, 0, 0};

				double posToIncPlane = 0;
				int ptInsideBox = 1;

				/*Check if inclined plane is present and what is the position of jello with respect to inclined plane*/
				if(jello->incPlanePresent)
				{
					posToIncPlane = jello->a * jello->p[i][j][k].x +
											jello->b * jello->p[i][j][k].y +
											jello->c * jello->p[i][j][k].z +
											jello->d;
				}

				/* Add Spring Forces */
				AddForce(jello, i, j, k, 1, 0, 0, &TotalForce);		// Stuctural Springs
				AddForce(jello, i, j, k, 0, 1, 0, &TotalForce);		// Stuctural Springs 
				AddForce(jello, i, j, k, 0, 0, 1, &TotalForce);		// Stuctural Springs
				
				AddForce(jello, i, j, k, 0, 0, -1, &TotalForce);	// Stuctural Springs
				AddForce(jello, i, j, k, 0, -1, 0, &TotalForce);	// Stuctural Springs
				AddForce(jello, i, j, k, -1, 0, 0, &TotalForce);	// Stuctural Springs

				AddForce(jello, i, j, k, 0, 1, 1, &TotalForce);		// Shear Springs
				AddForce(jello, i, j, k, 0, 1, -1, &TotalForce);	// Shear Springs 
				AddForce(jello, i, j, k, 0, -1, 1, &TotalForce);	// Shear Springs
				AddForce(jello, i, j, k, 0, -1, -1, &TotalForce);	// Shear Springs

				AddForce(jello, i, j, k, 1, 1, 0, &TotalForce);		// Shear Springs
				AddForce(jello, i, j, k, -1, 1, 0, &TotalForce);	// Shear Springs 
				AddForce(jello, i, j, k, 1, -1, 0, &TotalForce);	// Shear Springs
				AddForce(jello, i, j, k, -1, -1, 0, &TotalForce);	// Shear Springs

				AddForce(jello, i, j, k, 1 , 0, 1, &TotalForce);	// Shear Springs
				AddForce(jello, i, j, k, -1, 0, 1, &TotalForce);	// Shear Springs 
				AddForce(jello, i, j, k, 1, 0, -1, &TotalForce);	// Shear Springs
				AddForce(jello, i, j, k, -1, 0, -1, &TotalForce);	// Shear Springs

				AddForce(jello, i, j, k, 1 , 1, 1, &TotalForce);	// Shear Springs
				AddForce(jello, i, j, k, -1, -1, 1, &TotalForce);	// Shear Springs 
				AddForce(jello, i, j, k, 1, -1, 1, &TotalForce);	// Shear Springs
				AddForce(jello, i, j, k, -1, 1, 1, &TotalForce);	// Shear Springs

				AddForce(jello, i, j, k, 1 , 1, -1, &TotalForce);	// Shear Springs
				AddForce(jello, i, j, k, -1, -1, -1, &TotalForce);	// Shear Springs 
				AddForce(jello, i, j, k, 1, -1, -1, &TotalForce);	// Shear Springs
				AddForce(jello, i, j, k, -1, 1, -1, &TotalForce);	// Shear Springs

				AddForce(jello, i, j, k, 0, 0, 2, &TotalForce);		// Bend Springs
				AddForce(jello, i, j, k, 0, 2, 0, &TotalForce);		// Bend Springs 
				AddForce(jello, i, j, k, 2, 0, 0, &TotalForce);		// Bend Springs
				AddForce(jello, i, j, k, 0, 0, -2, &TotalForce);	// Bend Springs
				AddForce(jello, i, j, k, 0, -2, 0, &TotalForce);	// Bend Springs
				AddForce(jello, i, j, k, -2, 0, 0, &TotalForce);	// Bend Springs

				/* Add Collision Forces */
				if(i == 0 || i == 7 || j == 0 || j == 7 || k == 0 || k == 7)
				{
					ptInsideBox = 0;
					if(jello->p[i][j][k].x > 2.0)
					{
						point CollisionContact = {2.0, jello->p[i][j][k].y,jello->p[i][j][k].z};
						AddCollisionForce(jello,i,j,k,CollisionContact, &TotalForce);	
					}
					else if(jello->p[i][j][k].x < -2.0)
					{
						point CollisionContact = {-2.0, jello->p[i][j][k].y,jello->p[i][j][k].z};
						AddCollisionForce(jello,i,j,k,CollisionContact, &TotalForce);							
					}
					else if(jello->p[i][j][k].y > 2.0)
					{
						point CollisionContact = {jello->p[i][j][k].x,2.0,jello->p[i][j][k].z};
						AddCollisionForce(jello,i,j,k,CollisionContact, &TotalForce);
					}
					else if(jello->p[i][j][k].y < -2.0)
					{
						point CollisionContact = {jello->p[i][j][k].x,-2.0,jello->p[i][j][k].z};
						AddCollisionForce(jello,i,j,k,CollisionContact, &TotalForce);				
					}	
					else if(jello->p[i][j][k].z > 2.0)
					{
						point CollisionContact = {jello->p[i][j][k].x,jello->p[i][j][k].y, 2.0};
						AddCollisionForce(jello,i,j,k,CollisionContact, &TotalForce);				
					}
					else if(jello->p[i][j][k].z < -2.0)
					{
						point CollisionContact = {jello->p[i][j][k].x,jello->p[i][j][k].y, -2.0};
						AddCollisionForce(jello,i,j,k,CollisionContact, &TotalForce);				
					}
					//OriginalPos calculated during initialization tells whether the cube has originally positive or negative orientation to the plane
					else if((posToIncPlane * OriginalPos) < 0) 
					{
						// equations taken from http://www.9math.com/book/projection-point-plane

						/* Equation of line's orthogonal intersection to plane */
						double c = posToIncPlane / (jello->a * jello->a + jello->b * jello->b + jello->c * jello->c);
						point CollisionContact = {
													jello->p[i][j][k].x - jello->a * c,
													jello->p[i][j][k].y - jello->b * c,
													jello->p[i][j][k].z - jello->c * c,
													};

						AddCollisionForce(jello,i,j,k,CollisionContact, &TotalForce);
						ptInsideBox = 1;
					}
					else 
							ptInsideBox =  1;

					if(dispGeometry == 1)
					{
						/* Now we check if the particle collision with custom geometry objects added*/
						
						/*Calculation for cylinder*/
						double xdist = jello->p[i][j][k].x - cylinderPos.x;
						double ydist = jello->p[i][j][k].y - cylinderPos.y;
						double zdist = 0;
						double distSQ = (xdist*xdist)+(ydist*ydist)+ (zdist*zdist);

						if(((distSQ) < (cylinderBase*cylinderBase))
							&& ((jello->p[i][j][k].z - cylinderPos.z) < cylinderHeight ))
						{
							//Collision with cylinder
							point ptOfContact = {0,0,0};
							double dist[4], min;
							int minIdx = 0;
							double t = sqrt((cylinderBase*cylinderBase)/ distSQ);
							double x1 = cylinderPos.x + t*xdist;
							double y1 = cylinderPos.y + t*ydist;
							dist[0] = ((x1 - jello->p[i][j][k].x) * (x1 - jello->p[i][j][k].x) + (y1 - jello->p[i][j][k].y)* (y1 - jello->p[i][j][k].y));
							double x2 = cylinderPos.x - t*xdist;
							double y2 = cylinderPos.y - t*ydist;
							dist[1] = ((x2 - jello->p[i][j][k].x) * (x2 - jello->p[i][j][k].x) + (y2 - jello->p[i][j][k].y)* (y2 - jello->p[i][j][k].y));
							dist[2] = (jello->p[i][j][k].z - cylinderPos.z)* (jello->p[i][j][k].z - cylinderPos.z);
							dist[3] = (jello->p[i][j][k].z - cylinderPos.z + cylinderHeight)* (jello->p[i][j][k].z - cylinderPos.z + cylinderHeight);
							min =  dist[minIdx];
							for(int m = 1;m < 4;m++)
							{
								if(dist[m] < min){minIdx = m; min = dist[m];}
							}
							if(minIdx == 0) {ptOfContact.x = x1;ptOfContact.y = y1;ptOfContact.z = jello->p[i][j][k].z;}
							else if(minIdx == 1) {ptOfContact.x = x2;ptOfContact.y = y2;ptOfContact.z = jello->p[i][j][k].z;}
							else if(minIdx == 2) {ptOfContact.x = jello->p[i][j][k].x;ptOfContact.y = jello->p[i][j][k].y;ptOfContact.z = cylinderPos.z;}
							else if(minIdx == 3) {ptOfContact.x = jello->p[i][j][k].x;ptOfContact.y = jello->p[i][j][k].y;ptOfContact.z = cylinderPos.z+cylinderHeight;}

							AddCollisionForce(jello,i,j,k,ptOfContact, &TotalForce);
							ptInsideBox = 0;
						}
						
						/*Calculation for sphere*/
						double xdist1 = jello->p[i][j][k].x - spherePos.x;
						double ydist1 = jello->p[i][j][k].y - spherePos.y;
						double zdist1 = jello->p[i][j][k].z - spherePos.z;
						double distSQ1 = (xdist1*xdist1)+(ydist1*ydist1)+ (zdist1*zdist1);

						if((xdist1*xdist1 + ydist1*ydist1 + zdist1*zdist1) < (sphereRadius*sphereRadius))
						{
							//Collision with a sphere
							double t = sqrt((sphereRadius*sphereRadius)/ distSQ1);
							point ptOfContact1 = {
													spherePos.x + t*xdist1,
													spherePos.y + t*ydist1,
													spherePos.z + t*zdist1,
												};
							point ptOfContact2 = {
													spherePos.x - t*xdist1,
													spherePos.y - t*ydist1,
													spherePos.z - t*zdist1,
												};

							double dist1 = ((ptOfContact1.x - jello->p[i][j][k].x)*(ptOfContact1.x - jello->p[i][j][k].x) +
											(ptOfContact1.y - jello->p[i][j][k].y)*(ptOfContact1.y - jello->p[i][j][k].y) +
											(ptOfContact1.z - jello->p[i][j][k].z)*(ptOfContact1.z - jello->p[i][j][k].z));
							double dist2 = ((ptOfContact2.x - jello->p[i][j][k].x)*(ptOfContact2.x - jello->p[i][j][k].x) +
											(ptOfContact2.y - jello->p[i][j][k].y)*(ptOfContact2.y - jello->p[i][j][k].y) +
											(ptOfContact2.z - jello->p[i][j][k].z)*(ptOfContact2.z - jello->p[i][j][k].z)); 

							if(dist1 < dist2)
								AddCollisionForce(jello,i,j,k,ptOfContact1, &TotalForce);
							else 
								AddCollisionForce(jello,i,j,k,ptOfContact2, &TotalForce);
							//ptInsideBox = 0;
						}
					}
					
				}
				
				/* Add Force Field */
				if(jello->resolution && (ptInsideBox == 1))
				{
					AddInterpolatedForce(jello, i, j, k, h, &TotalForce);
				}

				/* Add External Forces generated by Mouse Drag*/

				TotalForce.x += externalForce.x;
				TotalForce.y += externalForce.y;
				TotalForce.z += externalForce.z;

				/* Compute acceleration from total Force*/
				double massReciprocal = 1 / jello->mass;
				pMULTIPLY(TotalForce,massReciprocal,a[i][j][k]);
			}
		}
	}

}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
