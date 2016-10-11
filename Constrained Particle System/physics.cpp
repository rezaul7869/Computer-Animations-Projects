/*

  USC/Viterbi/Computer Science
  CONSTRAINED PARTICLE SYSTEM
  Using the -
  "Jello Cube" Assignment 1 starter code
*/
//#define DUMP_TO_FILE

#include "particle.h"
#include "physics.h"
#include <gsl/gsl_linalg.h>

double eps=1E-6;

#ifdef DUMP_TO_FILE
	//int frameNo = 0;
#endif

/* Computes acceleration for all the particles, 
   which is in state given by 'particle'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * particle, struct point* a)
{
	double alpha	= particle->alpha;
	double beta		= particle->beta;

	/////////////////MATRIX B////////////////////////////////////////

#ifdef DUMP_TO_FILE
	FILE * file = fopen("C:\\Users\\RezaulAkkdkram\\Desktop\\new.txt", "a+");
#endif

	//this is zero indexed - calculates the index. COnverts 2D to 1D indices
	#define INDEX(x,y)  (((x) * (ncols)) + (y)) 
	//This calculates the index of transpose of deltaC in top right
	#define DELTAC_TRANSPOSE_INDEX(x,y) INDEX((y), ((x) + (2 * particle->nParticles)))
	//This calculates the index of deltaC in left bottom
	#define DELTAC_INDEX(x,y) INDEX(((x) + (2 * particle->nParticles)), (y))

	#define DISPLAY_MATRIX(A, r, c)\
	{\
		for(int _index_r_ = 0;_index_r_ < r;_index_r_++){\
			for(int _index_c_ = 0;_index_c_ < c;_index_c_++){\
				fprintf(file,"%.3f ", A[INDEX((_index_r_),(_index_c_))]);\
			}\
			fprintf(file,"\n");\
		}\
	}\

	#define DISPLAY_MATRIX_B(A, r, c)\
	{\
		for(int _index_r_ = 0;_index_r_ < r;_index_r_++){\
			fprintf(file,"%.3f ", A[(_index_r_)]);\
			fprintf(file,"\n");\
		}\
	}\

#ifdef DUMP_TO_FILE
	fprintf(file,"Ax = B \n MATRIX A\n");
#endif

	int nrows = ((3 * particle->nParticles) + 2);
	int ncols = ((3 * particle->nParticles) + 2);

	double* MatrixA = new double[nrows * ncols];
	memset(MatrixA, 0, sizeof(double) * nrows * ncols);

	//Put M in top - left
	for(int i = 0; i < (2 * particle->nParticles); i++)
	{
		MatrixA[INDEX(i,i)] = particle->mass;
	}

	//Put deltaC_T in top right
	
	MatrixA[DELTAC_TRANSPOSE_INDEX(0,0)] = 1;
	MatrixA[DELTAC_TRANSPOSE_INDEX(1,1)] = 1;

	MatrixA[DELTAC_TRANSPOSE_INDEX((particle->nParticles+1),(2*(particle->nParticles - 1)))] = 2 * particle->p[particle->nParticles -1].x ; // (NO_OF_PARTICLES-1)+2
	MatrixA[DELTAC_TRANSPOSE_INDEX((particle->nParticles+1),((2*(particle->nParticles - 1)))+1)] = 2 * (particle->p[particle->nParticles -1].y + particle->disk.radius);
	
	for(int i = 0; i < (particle->nParticles - 1) ; i++)
	{
		MatrixA[DELTAC_TRANSPOSE_INDEX((i+2),(i*2))] = 2 * (particle->p[i].x - particle->p[i+1].x);
		MatrixA[DELTAC_TRANSPOSE_INDEX((i+2),((i*2)+1))] = 2 * (particle->p[i].y - particle->p[i+1].y);
		MatrixA[DELTAC_TRANSPOSE_INDEX((i+2),((i*2)+2))] = 2 * (particle->p[i+1].x - particle->p[i].x);
		MatrixA[DELTAC_TRANSPOSE_INDEX((i+2),((i*2)+3))] = 2 * (particle->p[i+1].y - particle->p[i].y);
	}


	// This computes the index of deltaC in left bottom
	MatrixA[DELTAC_INDEX(0,0)] = 1;
	MatrixA[DELTAC_INDEX(1,1)] = 1;
	MatrixA[DELTAC_INDEX((particle->nParticles+1),(2*(particle->nParticles - 1)))] = 2 * particle->p[particle->nParticles -1].x ; // (NO_OF_PARTICLES-1)+2
	MatrixA[DELTAC_INDEX((particle->nParticles+1),((2*(particle->nParticles - 1)))+1)] = 2 * (particle->p[particle->nParticles -1].y + particle->disk.radius);
	
	for(int i = 0; i < (particle->nParticles - 1) ; i++)
	{
		MatrixA[DELTAC_INDEX((i+2),(i*2))] = (2 * (particle->p[i].x - particle->p[i+1].x));
		MatrixA[DELTAC_INDEX((i+2),((i*2)+1))] = (2 * (particle->p[i].y - particle->p[i+1].y));
		MatrixA[DELTAC_INDEX((i+2),((i*2)+2))] = (2 * (particle->p[i+1].x - particle->p[i].x));
		MatrixA[DELTAC_INDEX((i+2),((i*2)+3))] = (2 * (particle->p[i+1].y - particle->p[i].y));
	}

#ifdef DUMP_TO_FILE	
	DISPLAY_MATRIX(MatrixA, nrows, ncols);
#endif

	/////////////////MATRIX B////////////////////////////////////////

#ifdef DUMP_TO_FILE
	fprintf(file,"MATRIX B\n");
#endif

	#define INDEX_DC_DQ(x)  (((x) + (2* particle->nParticles))) 
	#define V_X(i) ((particle->v[(i)].x))
	#define V_Y(i) ((particle->v[(i)].y))

	int nrowsB = (3 * particle->nParticles) + 2;
	
	double* MatrixB = new double[nrowsB];
	memset(MatrixB, 0, sizeof(double) * nrowsB);

	//Setup the Fext part
	for(int i = 0; i < particle->nParticles; i++)
	{
		MatrixB[2*i] = particle->force[i].x;
		MatrixB[2*i + 1] = particle->force[i].y;
	}


	//Setup complex part!!!
	MatrixB[INDEX_DC_DQ(0)] = 0;
	MatrixB[INDEX_DC_DQ(1)] = 0;

	for(int i = 0; i < (particle->nParticles - 1) ; i++)
	{
		double value1 = 2 * (V_X((i)) - V_X((i+1)));
		double value2 = 2 * (V_Y((i)) - V_Y((i+1)));
		double value3 = 2 * (V_X((i+1)) - V_X((i)));
		double value4 = 2 * (V_Y((i+1)) - V_Y((i)));
	
		MatrixB[INDEX_DC_DQ(i+2)] -= (value1 * V_X((i)) +
										value2 * V_Y((i)) +
										value3 * V_X((i+1)) +
										value4 * V_Y((i+1)));
	}

	double temp1 = (2 * particle->v[(particle->nParticles - 1)].x);
	double temp2 = (2 * particle->v[(particle->nParticles - 1)].y);

	MatrixB[INDEX_DC_DQ(((particle->nParticles - 1) + 2))] -= (temp1 * particle->v[(particle->nParticles - 1)].x +
															temp2 * particle->v[(particle->nParticles - 1)].y);
#ifdef DUMP_TO_FILE	
	DISPLAY_MATRIX_B(MatrixB, nrows, 1);
#endif

	

	//Add the BAUMGARTE STABILIZATION

	
	MatrixB[INDEX_DC_DQ(0)] -= ((alpha*(particle->v[0].x)) + (beta*(particle->p[0].x)));
	MatrixB[INDEX_DC_DQ(1)] -= ((alpha*(particle->v[0].y)) + (beta*(particle->p[0].y)));;

	double particleDist = (particle->disk.radius * 2) / (particle->nParticles - 1); 

	for(int i = 0; i < (particle->nParticles - 1) ; i++)
	{
		double value1 = 2 * (particle->p[(i+1)].x - particle->p[i].x) *
							(V_X((i+1)) - V_X((i)));
		double value2 = 2 * (particle->p[(i+1)].y - particle->p[i].y) *
							(V_Y((i+1)) - V_Y((i)));
		double value3 =  ((particle->p[(i+1)].x - particle->p[i].x) * (particle->p[(i+1)].x - particle->p[i].x)) + 
							((particle->p[(i+1)].y - particle->p[i].y) * (particle->p[(i+1)].y - particle->p[i].y)) -
							(particleDist * particleDist);

		MatrixB[INDEX_DC_DQ(i+2)] -= ((alpha * (value1 + value2)) +
										(beta * value3));
	}

	temp1 = (2 * (particle->p[(particle->nParticles - 1)].x - particle->disk.center.x) * 
						particle->v[(particle->nParticles - 1)].x);
	temp2 = (2 * (particle->p[(particle->nParticles - 1)].y - particle->disk.center.y) * 
						particle->v[(particle->nParticles - 1)].y);

	double temp3 =  ((particle->p[(particle->nParticles - 1)].x - particle->disk.center.x) * (particle->p[(particle->nParticles - 1)].x - particle->disk.center.x)) + 
							((particle->p[(particle->nParticles - 1)].y - particle->disk.center.y) * (particle->p[(particle->nParticles - 1)].y - particle->disk.center.y)) - 
							(particle->disk.radius * particle->disk.radius);

	MatrixB[INDEX_DC_DQ(((particle->nParticles - 1) + 2))] -= ((alpha *(temp1 + temp2)) + (beta*temp3));
	
#ifdef DUMP_TO_FILE	
	DISPLAY_MATRIX_B(MatrixB, nrows, 1);
#endif
	
	/////////////////////////////////////////////////////////////////////

	//Find the acceleration using gsl

	int s;

	gsl_matrix_view m	= gsl_matrix_view_array(MatrixA, nrows, ncols);
	gsl_vector_view b	= gsl_vector_view_array(MatrixB, nrowsB);
	
	gsl_vector* x = gsl_vector_alloc(nrows);
	gsl_vector* work = gsl_vector_alloc(nrows);
	gsl_vector* S = gsl_vector_alloc(nrows);
	gsl_matrix* V = gsl_matrix_alloc(nrows, nrows);
	
	gsl_linalg_SV_decomp(&m.matrix, V, S, work);

	double maxValue =  gsl_vector_get(S, 0);
	for(int i = 0; i < nrows; i++)
	{
		if(gsl_vector_get(S, i) < (eps * maxValue))
			gsl_vector_set(S, i, 0.0);
	}

	gsl_linalg_SV_solve(&m.matrix, V, S, &b.vector, x);

	for(int i = 0; i < particle->nParticles; i++)
	{
		a[i].x = x->data[i*2];
		a[i].y = x->data[(i*2) + 1];
		a[i].z = 0;
	}

	gsl_vector_free (x);
	gsl_vector_free (S);
	gsl_vector_free (work);
	gsl_matrix_free(V);
	/////////////////////////////////////////////////////////////////////
	delete(MatrixA);
	delete(MatrixB);

#ifdef DUMP_TO_FILE
	fclose(file);
#endif
}

// Compute Acceleration with CRing off


void computeAccelerationWoCRing(struct world * particle, struct point* a)
{
	double alpha	= particle->alpha;
	double beta		= particle->beta;

	/////////////////MATRIX B////////////////////////////////////////

#ifdef DUMP_TO_FILE
	FILE * file = fopen("C:\\Users\\RezaulAkram\\Desktop\\new.txt", "a+");
#endif

	//this is zero indexed - calculates the index. COnverts 2D to 1D indices
	#define INDEX(x,y)  (((x) * (ncols)) + (y)) 
	//This calculates the index of transpose of deltaC in top right
	#define DELTAC_TRANSPOSE_INDEX(x,y) INDEX((y), ((x) + (2 * particle->nParticles)))
	//This calculates the index of deltaC in left bottom
	#define DELTAC_INDEX(x,y) INDEX(((x) + (2 * particle->nParticles)), (y))

	#define DISPLAY_MATRIX(A, r, c)\
	{\
		for(int _index_r_ = 0;_index_r_ < r;_index_r_++){\
			for(int _index_c_ = 0;_index_c_ < c;_index_c_++){\
				fprintf(file,"%.3f ", A[INDEX((_index_r_),(_index_c_))]);\
			}\
			fprintf(file,"\n");\
		}\
	}\

	#define DISPLAY_MATRIX_B(A, r, c)\
	{\
		for(int _index_r_ = 0;_index_r_ < r;_index_r_++){\
			fprintf(file,"%.3f ", A[(_index_r_)]);\
			fprintf(file,"\n");\
		}\
	}\

#ifdef DUMP_TO_FILE
	fprintf(file,"Ax = B \n MATRIX A\n");
#endif

	int nrows = ((3 * particle->nParticles) + 1);
	int ncols = ((3 * particle->nParticles) + 1);

	double* MatrixA = new double[nrows * ncols];
	memset(MatrixA, 0, sizeof(double) * nrows * ncols);

	//Put M in top - left
	for(int i = 0; i < (2 * particle->nParticles); i++)
	{
		MatrixA[INDEX(i,i)] = particle->mass;
	}

	//Put deltaC_T in top right
	
	MatrixA[DELTAC_TRANSPOSE_INDEX(0,0)] = 1;
	MatrixA[DELTAC_TRANSPOSE_INDEX(1,1)] = 1;

	for(int i = 0; i < (particle->nParticles - 1) ; i++)
	{
		MatrixA[DELTAC_TRANSPOSE_INDEX((i+2),(i*2))] = 2 * (particle->p[i].x - particle->p[i+1].x);
		MatrixA[DELTAC_TRANSPOSE_INDEX((i+2),((i*2)+1))] = 2 * (particle->p[i].y - particle->p[i+1].y);
		MatrixA[DELTAC_TRANSPOSE_INDEX((i+2),((i*2)+2))] = 2 * (particle->p[i+1].x - particle->p[i].x);
		MatrixA[DELTAC_TRANSPOSE_INDEX((i+2),((i*2)+3))] = 2 * (particle->p[i+1].y - particle->p[i].y);
	}


	// This computes the index of deltaC in left bottom
	MatrixA[DELTAC_INDEX(0,0)] = 1;
	MatrixA[DELTAC_INDEX(1,1)] = 1;
	
	for(int i = 0; i < (particle->nParticles - 1) ; i++)
	{
		MatrixA[DELTAC_INDEX((i+2),(i*2))] = (2 * (particle->p[i].x - particle->p[i+1].x));
		MatrixA[DELTAC_INDEX((i+2),((i*2)+1))] = (2 * (particle->p[i].y - particle->p[i+1].y));
		MatrixA[DELTAC_INDEX((i+2),((i*2)+2))] = (2 * (particle->p[i+1].x - particle->p[i].x));
		MatrixA[DELTAC_INDEX((i+2),((i*2)+3))] = (2 * (particle->p[i+1].y - particle->p[i].y));
	}

#ifdef DUMP_TO_FILE	
	DISPLAY_MATRIX(MatrixA, nrows, ncols);
#endif

	/////////////////MATRIX B////////////////////////////////////////

#ifdef DUMP_TO_FILE
	fprintf(file,"MATRIX B\n");
#endif

	#define INDEX_DC_DQ(x)  (((x) + (2* particle->nParticles))) 
	#define V_X(i) ((particle->v[(i)].x))
	#define V_Y(i) ((particle->v[(i)].y))

	int nrowsB = (3 * particle->nParticles) + 1;
	
	double* MatrixB = new double[nrowsB];
	memset(MatrixB, 0, sizeof(double) * nrowsB);

	//Setup the Fext part
	for(int i = 0; i < particle->nParticles; i++)
	{
		MatrixB[2*i] = particle->force[i].x;
		MatrixB[2*i + 1] = particle->force[i].y;
	}


	//Setup complex part!!!
	MatrixB[INDEX_DC_DQ(0)] = 0;
	MatrixB[INDEX_DC_DQ(1)] = 0;

	for(int i = 0; i < (particle->nParticles - 1) ; i++)
	{
		double value1 = 2 * (V_X((i)) - V_X((i+1)));
		double value2 = 2 * (V_Y((i)) - V_Y((i+1)));
		double value3 = 2 * (V_X((i+1)) - V_X((i)));
		double value4 = 2 * (V_Y((i+1)) - V_Y((i)));
	
		MatrixB[INDEX_DC_DQ(i+2)] -= (value1 * V_X((i)) +
										value2 * V_Y((i)) +
										value3 * V_X((i+1)) +
										value4 * V_Y((i+1)));
	}

#ifdef DUMP_TO_FILE	
	DISPLAY_MATRIX_B(MatrixB, nrows, 1);
#endif

	

	//Add the BAUMGARTE STABILIZATION

	
	MatrixB[INDEX_DC_DQ(0)] -= ((alpha*(particle->v[0].x)) + (beta*(particle->p[0].x)));
	MatrixB[INDEX_DC_DQ(1)] -= ((alpha*(particle->v[0].y)) + (beta*(particle->p[0].y)));;

	double particleDist = (particle->disk.radius * 2) / (particle->nParticles - 1); 

	for(int i = 0; i < (particle->nParticles - 1) ; i++)
	{
		double value1 = 2 * (particle->p[(i+1)].x - particle->p[i].x) *
							(V_X((i+1)) - V_X((i)));
		double value2 = 2 * (particle->p[(i+1)].y - particle->p[i].y) *
							(V_Y((i+1)) - V_Y((i)));
		double value3 =  ((particle->p[(i+1)].x - particle->p[i].x) * (particle->p[(i+1)].x - particle->p[i].x)) + 
							((particle->p[(i+1)].y - particle->p[i].y) * (particle->p[(i+1)].y - particle->p[i].y)) -
							(particleDist * particleDist);

		MatrixB[INDEX_DC_DQ(i+2)] -= ((alpha * (value1 + value2)) +
										(beta * value3));
	}
	
#ifdef DUMP_TO_FILE	
	DISPLAY_MATRIX_B(MatrixB, nrows, 1);
#endif
	
	/////////////////////////////////////////////////////////////////////

	//Find the acceleration using gsl

	int s;

	gsl_matrix_view m	= gsl_matrix_view_array(MatrixA, nrows, ncols);
	gsl_vector_view b	= gsl_vector_view_array(MatrixB, nrowsB);
	
	gsl_vector* x = gsl_vector_alloc(nrows);
	gsl_vector* work = gsl_vector_alloc(nrows);
	gsl_vector* S = gsl_vector_alloc(nrows);
	gsl_matrix* V = gsl_matrix_alloc(nrows, nrows);
	
	gsl_linalg_SV_decomp(&m.matrix, V, S, work);

	double maxValue =  gsl_vector_get(S, 0);
	for(int i = 0; i < nrows; i++)
	{
		if(gsl_vector_get(S, i) < (eps * maxValue))
			gsl_vector_set(S, i, 0.0);
	}

	gsl_linalg_SV_solve(&m.matrix, V, S, &b.vector, x);

	for(int i = 0; i < particle->nParticles; i++)
	{
		a[i].x = x->data[i*2];
		a[i].y = x->data[(i*2) + 1];
		a[i].z = 0;
	}

	gsl_vector_free (x);
	gsl_vector_free (S);
	gsl_vector_free (work);
	gsl_matrix_free(V);
	/////////////////////////////////////////////////////////////////////
	delete(MatrixA);
	delete(MatrixB);

#ifdef DUMP_TO_FILE
	fclose(file);
#endif
}



/* performs one step of Euler Integration */
/* as a result, updates the particle structure */
void Euler(struct world * particle)
{
  int i;
  point* a = new point[particle->nParticles];
  
#ifdef DUMP_TO_FILE
  FILE * file = fopen("C:\\Users\\RezaulAkram\\Desktop\\new.txt", "a+");
//  fprintf(file,"frame %d\n", frameNo++);
  fclose(file);
#endif
  if(cRing == 1)
	  computeAcceleration(particle, a);
  else
	  computeAccelerationWoCRing(particle, a);
  
  for (i=0; i< particle->nParticles; i++)
  {
        particle->v[i].x += particle->dt * a[i].x;
        particle->v[i].y += particle->dt * a[i].y;
        particle->v[i].z += particle->dt * a[i].z;
		particle->p[i].x += particle->dt * particle->v[i].x;
        particle->p[i].y += particle->dt * particle->v[i].y;
        particle->p[i].z += particle->dt * particle->v[i].z;
  }

	delete(a);
}

/* performs one step of RK4 Integration */
/* as a result, updates the particle structure */
void RK4(struct world * particle)
{
  point* F1p = new point[particle->nParticles];
  point* F1v = new point[particle->nParticles]; 
  point* F2p = new point[particle->nParticles];
  point* F2v = new point[particle->nParticles];
  point* F3p = new point[particle->nParticles];
  point* F3v = new point[particle->nParticles];
  point* F4p = new point[particle->nParticles];
  point* F4v = new point[particle->nParticles];

  point* a = new point[particle->nParticles];


  struct world buffer;
  memcpy(&buffer, particle, sizeof(world));

  buffer.p = new point[particle->nParticles];
  buffer.v = new point[particle->nParticles];
  buffer.force = new point[particle->nParticles];

  int i,j,k;

  //buffer = *particle; // make a copy of particle

  
  memcpy(buffer.p, particle->p, sizeof(point) * particle->nParticles);
  memcpy(buffer.v, particle->v, sizeof(point) * particle->nParticles);
  memcpy(buffer.force, particle->force, sizeof(point) * particle->nParticles);

  computeAcceleration(particle, a);

  for (i=0; i< particle->nParticles; i++)
	{
         pMULTIPLY(particle->v[i],particle->dt,F1p[i]);
         pMULTIPLY(a[i],particle->dt,F1v[i]);
         pMULTIPLY(F1p[i],0.5,buffer.p[i]);
         pMULTIPLY(F1v[i],0.5,buffer.v[i]);
         pSUM(particle->p[i],buffer.p[i],buffer.p[i]);
         pSUM(particle->v[i],buffer.v[i],buffer.v[i]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i< particle->nParticles; i++)
  {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i],particle->dt,F2p[i]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i],particle->dt,F2v[i]);
         pMULTIPLY(F2p[i],0.5,buffer.p[i]);
         pMULTIPLY(F2v[i],0.5,buffer.v[i]);
         pSUM(particle->p[i],buffer.p[i],buffer.p[i]);
         pSUM(particle->v[i],buffer.v[i],buffer.v[i]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<particle->nParticles; i++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i],particle->dt,F3p[i]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i],particle->dt,F3v[i]);
         pMULTIPLY(F3p[i],0.5,buffer.p[i]);
         pMULTIPLY(F3v[i],0.5,buffer.v[i]);
         pSUM(particle->p[i],buffer.p[i],buffer.p[i]);
         pSUM(particle->v[i],buffer.v[i],buffer.v[i]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i< particle->nParticles; i++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i],particle->dt,F4p[i]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i],particle->dt,F4v[i]);

         pMULTIPLY(F2p[i],2,buffer.p[i]);
         pMULTIPLY(F3p[i],2,buffer.v[i]);
         pSUM(buffer.p[i],buffer.v[i],buffer.p[i]);
         pSUM(buffer.p[i],F1p[i],buffer.p[i]);
         pSUM(buffer.p[i],F4p[i],buffer.p[i]);
         pMULTIPLY(buffer.p[i],1.0 / 6,buffer.p[i]);
         pSUM(buffer.p[i],particle->p[i],particle->p[i]);

         pMULTIPLY(F2v[i],2,buffer.p[i]);
         pMULTIPLY(F3v[i],2,buffer.v[i]);
         pSUM(buffer.p[i],buffer.v[i],buffer.p[i]);
         pSUM(buffer.p[i],F1v[i],buffer.p[i]);
         pSUM(buffer.p[i],F4v[i],buffer.p[i]);
         pMULTIPLY(buffer.p[i],1.0 / 6,buffer.p[i]);
         pSUM(buffer.p[i],particle->v[i],particle->v[i]);
      }

  delete(buffer.p);
  delete(buffer.v);
  delete(buffer.force);

  delete(F1p);
  delete(F1v);
  delete(F2p);
  delete(F2v);
  delete(F3p);
  delete(F3v);
  delete(F4p);
  delete(F4v);
  delete(a);

  return;  
}
