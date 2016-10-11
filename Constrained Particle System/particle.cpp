/*

  USC/Viterbi/Computer Science
  CONSTRAINED PARTICLE SYSTEM
  Using the -
  "Jello Cube" Assignment 1 starter code

  Rezaul Akram Barbhuiya

*/

#include "particle.h"
#include "showParticle.h"
#include "input.h"
#include "physics.h"
#include "perlinNoise.h"

// camera parameters
double Theta = pi/2;
double Phi = 0;
double R = 2;

int frameNo;

double PCFreq = 0.0;
__int64 CounterStart = 0;
__int64 InitialTime = 0;

// mouse control
int g_iMenuId;
int g_vMousePos[2];
int g_iLeftMouseButton,g_iMiddleMouseButton,g_iRightMouseButton;

// number of images saved to disk so far
int sprite=0;

// these variables control what is displayed on screen
int pause=0, viewingMode=0, saveScreenToFile=0, singleFrame = 0, No_Of_Particles = NO_OF_PARTICLES, forceMode = 1, rotateCamera = 0, externalForce = 0, cRing = 1;
double underWaterConstant;
//int shear=0, bend=0, structural=1; 

struct world particle;

int windowWidth, windowHeight;

void ShowText(const char *text, int length, int x, int y){
	 glMatrixMode(GL_PROJECTION);
	 double m[16];
	 glGetDoublev(GL_PROJECTION_MATRIX, m);
	 glLoadIdentity();
	 glOrtho(0, 800, 0, 600, -5, 5);
	 glMatrixMode(GL_MODELVIEW);
	 glLoadIdentity();
	 glPushMatrix();
	 glLoadIdentity();
	 glRasterPos2i(x, y);

	 for(int i=0; i<length; i++){
	  glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, (int)text[i]);
	 }

	 glPopMatrix();

	 glMatrixMode(GL_PROJECTION);
	 glLoadMatrixd(m);
	 glMatrixMode(GL_MODELVIEW);
}
 

int DumpError()
{
	double error = 0;
	error += (particle.p[0].x * particle.p[0].x);
	error += (particle.p[0].y * particle.p[0].y);

	double particleDist = (particle.disk.radius * 2) / (particle.nParticles - 1); 
	for(int i = 0;  i < particle.nParticles - 1; i++)
	{
		error += ((particle.p[i+1].x - particle.p[i].x)*(particle.p[i+1].x - particle.p[i].x) +
			(particle.p[i+1].y - particle.p[i].y)*(particle.p[i+1].y - particle.p[i].y) - (particleDist * particleDist)) ;
	}

	if(cRing == 1)
	{
		error += ((particle.p[particle.nParticles - 1].x - particle.disk.center.x)*(particle.p[particle.nParticles - 1].x - particle.disk.center.x) +
			(particle.p[particle.nParticles - 1].y - particle.disk.center.y)*(particle.p[particle.nParticles - 1].y - particle.disk.center.y) - (particle.disk.radius * particle.disk.radius)) ;
	}

	FILE * file = fopen("C:\\Users\\RezaulAkram\\Desktopdd\\new.txt", "a+");
	fprintf(file,"%f\n", error);
	fclose(file);
	

	return 0;
}

void myinit()
{
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(90.0,1.0,0.01,1000.0);

  // set background color to grey
  glClearColor(0.5, 0.5, 0.5, 0.0);

  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_LINE_SMOOTH);

  return; 
}

void reshape(int w, int h) 
{
  // Prevent a divide by zero, when h is zero.
  // You can't make a window of zero height.
  if(h == 0)
    h = 1;

  glViewport(0, 0, w, h);

  // Reset the coordinate system before modifying
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // Set the perspective
  double aspectRatio = 1.0 * w / h;
  gluPerspective(60.0f, aspectRatio, 0.01f, 1000.0f);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity(); 

  windowWidth = w;
  windowHeight = h;

  glutPostRedisplay();
}

void display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // camera parameters are Phi, Theta, R
  gluLookAt(R * cos(Phi) * cos (Theta), R * sin(Phi) * cos (Theta) - 0.5, R * sin (Theta),
	        0.0,-0.5,0.0, 0.0,0.5,0.0);


  /* Lighting */
  /* You are encouraged to change lighting parameters or make improvements/modifications
     to the lighting model . 
     This way, you will personalize your assignment and your assignment will stick out. 
  */

  // global ambient light
  GLfloat aGa[] = { 0.0, 0.0, 0.0, 0.0 };
  
  // light 's ambient, diffuse, specular
  GLfloat lKa1[] = { 0.1, 0.1, 0.1, 1.0 };
  GLfloat lKd1[] = { 0.1, 0.1, 0.1, 1.0 };
  GLfloat lKs1[] = { 0.1, 0.1, 0.1, 1.0 };

  GLfloat lKa0[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat lKd0[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat lKs0[] = { 1.0, 1.0, 1.0, 1.0 };

  GLfloat lKa2[] = { 0.0, 1.0, 1.0, 1.0 };
  GLfloat lKd2[] = { 0.0, 1.0, 1.0, 1.0 };
  GLfloat lKs2[] = { 0.0, 1.0, 1.0, 1.0 };

  GLfloat lKa3[] = { 1.0, 1.0, 0.0, 1.0 };
  GLfloat lKd3[] = { 1.0, 1.0, 0.0, 1.0 };
  GLfloat lKs3[] = { 1.0, 1.0, 0.0, 1.0 };

  GLfloat lKa4[] = { 0.0, 0.0, 0.0, 1.0 };
  GLfloat lKd4[] = { 0.0, 0.0, 1.0, 1.0 };
  GLfloat lKs4[] = { 0.0, 0.0, 1.0, 1.0 };

  GLfloat lKa5[] = { 1.0, 0.0, 1.0, 1.0 };
  GLfloat lKd5[] = { 1.0, 0.0, 1.0, 1.0 };
  GLfloat lKs5[] = { 1.0, 0.0, 1.0, 1.0 };

  GLfloat lKa6[] = { 0.0, 1.0, 1.0, 1.0 };
  GLfloat lKd6[] = { 0.0, 1.0, 1.0, 1.0 };
  GLfloat lKs6[] = { 0.0, 1.0, 1.0, 1.0 };

  GLfloat lKa7[] = { 0.0, 0.0, 1.0, 1.0 };
  GLfloat lKd7[] = { 0.0, 0.0, 1.0, 1.0 };
  GLfloat lKs7[] = { 0.0, 0.0, 1.0, 1.0 };

  // light positions and directions
  GLfloat lP0[] = { -8.0,-8.0, 3, 1.0 };
  GLfloat lP1[] = { -8.0,-8.0, -3, 1.0 };
  GLfloat lP2[] = { -8.0,8.0, 3, 1.0 };
  GLfloat lP3[] = { -8.0,8.0, -3, 1.0 };
  GLfloat lP4[] = { 8.0,-8.0, 3, 1.0 };
  GLfloat lP5[] = { 8.0,-8.0, -3, 1.0 };
  GLfloat lP6[] = { 8.0,8.0, 3, 1.0 };
  GLfloat lP7[] = { 8.0,8.0, -3, 1.0 };
  
  // jelly material color

  GLfloat mKa[] = { 0.2, 0.2, 0.2, 1.0 };
  GLfloat mKd[] = { 0.4, 0.4, 0.4, 1.0 };
  GLfloat mKs[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mKe[] = { 0.1, 0.1, 0.1, 1.0 };

  /* set up lighting */
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, aGa);
  glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

  // set up cube color
  glMaterialfv(GL_FRONT, GL_AMBIENT, mKa);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, mKd);
  glMaterialfv(GL_FRONT, GL_SPECULAR, mKs);
  glMaterialfv(GL_FRONT, GL_EMISSION, mKe);
  glMaterialf(GL_FRONT, GL_SHININESS, 50);
    
  // macro to set up light i
  #define LIGHTSETUP(i)\
  glLightfv(GL_LIGHT##i, GL_POSITION, lP##i);\
  glLightfv(GL_LIGHT##i, GL_AMBIENT, lKa##i);\
  glLightfv(GL_LIGHT##i, GL_DIFFUSE, lKd##i);\
  glLightfv(GL_LIGHT##i, GL_SPECULAR, lKs##i);\
  glEnable(GL_LIGHT##i)
  
  LIGHTSETUP (0);
  LIGHTSETUP (1);
  LIGHTSETUP (2);
  LIGHTSETUP (3);
  LIGHTSETUP (4);
  LIGHTSETUP (5);
  LIGHTSETUP (6);
  LIGHTSETUP (7);

  // enable lighting
  glEnable(GL_LIGHTING);    
  glEnable(GL_DEPTH_TEST);
  glShadeModel(GL_SMOOTH);
  glEnable(GL_COLOR_MATERIAL);

  if(particle.nParticles != No_Of_Particles)
	  ReinitializeWorld(&particle);

  // show the cube
  showParticles(&particle);
  //showDisc(particle.disk);
  showRoom();
  initTorus(particle.disk);
  glCallList(theTorus);
  char text[80];
  ShowText("Keys:",5,50,220);
  sprintf(text, "e : Reset Camera");
  ShowText(text,strlen(text),50,200);
  sprintf(text, "a : Add Left Force");
  ShowText(text,strlen(text),50,180);
  sprintf(text, "d : Add Right Force");
  ShowText(text,strlen(text),50,160);
  sprintf(text, "t : Toggle Integrater Euler/RK4");
  ShowText(text,strlen(text),50,140);
  sprintf(text, "v : Toggle Camera Rotate");
  ShowText(text,strlen(text),50,120);
  sprintf(text, "p : Pause");
  ShowText(text,strlen(text),50,100);
  sprintf(text, "j : Single Frame");
  ShowText(text,strlen(text),50,80);
  sprintf(text, "z/x : Zoom In/Out");
  ShowText(text,strlen(text),50,60);
  sprintf(text, "c : Toggle CRing Constraint");
  ShowText(text,strlen(text),50,40);
  sprintf(text, "b/n/m- Force = Constant/Perlin/UnderWater");
  ShowText(text,strlen(text),50,20);

  LARGE_INTEGER li;
  QueryPerformanceCounter(&li);
  double currentTime = (li.QuadPart-CounterStart)/PCFreq;
  CounterStart = li.QuadPart;
  sprintf(text, "Frame Rate : %f",((1000/currentTime) * particle.n));
  ShowText(text,strlen(text),550,200);
  sprintf(text, "Time Elapsed : %.3f",(li.QuadPart-InitialTime)/PCFreq);
  ShowText(text,strlen(text),550,180);
  sprintf(text, "Frames Elapsed : %d",frameNo);
  ShowText(text,strlen(text),550,160);
  sprintf(text, "Integrator : %s",(particle.integrator == 0) ? "Euler":"RK4");
  ShowText(text,strlen(text),550,140);

  glFlush();

  glDisable(GL_LIGHTING);

  glutSwapBuffers();
}

void doIdle()
{

  char s[20]="picxxxx.ppm";
  int i;
  frameNo++;
  // save screen to file
  s[3] = 48 + (sprite / 1000);
  s[4] = 48 + (sprite % 1000) / 100;
  s[5] = 48 + (sprite % 100 ) / 10;
  s[6] = 48 + sprite % 10;

  if (saveScreenToFile==1)
  {
    saveScreenshot(windowWidth, windowHeight, s);
    saveScreenToFile=0; // save only once, change this if you want continuos image generation (i.e. animation)
    sprite++;
  }

  if (sprite >= 300) // allow only 300 snapshots
  {
    exit(0);	
  }

  if(rotateCamera == 1)
	  Theta += 0.02;

  if (pause == 0 || singleFrame == 1)
  {
	  if(forceMode == 1)
	  {
			for(int i = 0; i < particle.nParticles; i++)
			{
			particle.force[i].y = particle.mass * (-1);
			particle.force[i].z = particle.force[i].x = 0;
			}
			double randomForce = PerlinNoise::PerlinNoise1D((rand()*1.0/RAND_MAX), 5, 0.8);
			int q =   (floor(randomForce));
			particle.force[particle.nParticles - 1].x = ((randomForce - q * 1) - 0.5) * 2 * particle.mass;
	  }
	  else if(forceMode == 2)
	  {
			for(int i = 0; i < particle.nParticles; i++)
			{
				particle.force[i].x -= (underWaterConstant * particle.v[i].x);
				particle.force[i].y -= (underWaterConstant * particle.v[i].y);
			}
	  }

	  particle.force[particle.nParticles - 1].x += (externalForce * particle.mass);
	  externalForce = 0;
	  // insert code which appropriately performs one step of the cube simulation:
	  for(int i = 0; i < particle.n; i++)
	  {
		  if(particle.integrator == 0)
			Euler(&particle);
		  else if(particle.integrator == 1)
			RK4(&particle);
		  else {}
		if(frameNo == 100)
			cRing = 0;
		else if(frameNo == 200)
			cRing = 1;
		if(frameNo >=  200 && frameNo <  1000)
			DumpError();
	  }


	  singleFrame = 0;
  }

  glutPostRedisplay();
}

int InitializeDisk(struct Disk* disk)
{
	frameNo = 0;
	double step = (2 * pi)/ disk->no_of_segments;
	if(disk->p) free(disk->p);
	disk->p = (point*) malloc(sizeof(point)*(disk->no_of_segments+1));
	if(!disk->p)
		return -1;
	memset(disk->p, 0, sizeof(point)*(disk->no_of_segments+1));

	for(int i = 0; i <= disk->no_of_segments; i++)
	{
		double r = step *i;
		disk->p[i].x = disk->center.x + disk->radius * cos(r);
		disk->p[i].y = disk->center.y + disk->radius * sin(r);
	}
	return 0;
}

void Deinitializedisk(struct Disk* disk)
{
	free(disk->p);
}

int ReinitializeWorld(struct world* particle)
{
	point* tempPoints =  new point[particle->nParticles];
	point* tempVelocity =  new point[particle->nParticles];
	point* tempForce =  new point[particle->nParticles];
	memcpy(tempPoints, particle->p, sizeof(point) * particle->nParticles);
	memcpy(tempVelocity, particle->v, sizeof(point) * particle->nParticles);
	memcpy(tempForce, particle->force, sizeof(point) * particle->nParticles);
	int prevNParticles =  particle->nParticles;

	particle->nParticles = No_Of_Particles;
	delete(particle->p);
	delete(particle->v);
	delete(particle->force);

	particle->p = new point[particle->nParticles];
	particle->v = new point[particle->nParticles];
	particle->force = new point[particle->nParticles];
	particle->alpha = 0.2;
	particle->beta = 0.5;

	//particle->alpha = 0;
	//particle->beta = 0;

	particle->mass = 1.0 / particle->nParticles;
	
	double particleDist = (particle->disk.radius * 2) / (particle->nParticles - 1); 
	

	double ratio = ((prevNParticles-1) * 1.0)/((particle->nParticles - 1) * 1.0);
	for(int i = 0; i < particle->nParticles; i++)
	{
		double position =  i * ratio;
		int beadLoc1 = floor(position);
		int beadLoc2 = beadLoc1 + 1;
		double decPosition =  position - beadLoc1;

		
		if(decPosition == 0)
		{
			particle->p[i] = tempPoints[beadLoc1];
			particle->v[i] = tempVelocity[beadLoc1];
			particle->force[i] = tempForce[beadLoc1];
		}
		else
		{
			particle->p[i].x =  tempPoints[beadLoc1].x * (1 - decPosition) + tempPoints[beadLoc2].x * decPosition;
			particle->p[i].y =  tempPoints[beadLoc1].y * (1 - decPosition) + tempPoints[beadLoc2].y * decPosition;
			particle->p[i].z =  tempPoints[beadLoc1].z * (1 - decPosition) + tempPoints[beadLoc2].z * decPosition;

			particle->v[i].x =  tempVelocity[beadLoc1].x * (1 - decPosition) + tempVelocity[beadLoc2].x * decPosition;
			particle->v[i].y =  tempVelocity[beadLoc1].y * (1 - decPosition) + tempVelocity[beadLoc2].y * decPosition;
			particle->v[i].z =  tempVelocity[beadLoc1].z * (1 - decPosition) + tempVelocity[beadLoc2].z * decPosition;

			/*
			particle->force[i].x =  tempForce[beadLoc1].x * (1 - decPosition) + tempForce[beadLoc2].x * decPosition;
			particle->force[i].y =  tempForce[beadLoc1].y * (1 - decPosition) + tempForce[beadLoc2].y * decPosition;
			particle->force[i].z =  tempForce[beadLoc1].z * (1 - decPosition) + tempForce[beadLoc2].z * decPosition;
			*/
			particle->force[i].z = 0;
			particle->force[i].y = particle->mass * (-1);
		}

		if(forceMode == 1)
		{
			for(int i = 0; i < particle->nParticles; i++)
			{
				particle->force[i].y = particle->mass * (-1);
				particle->force[i].z = particle->force[i].x = 0;
			}
				double randomForce = PerlinNoise::PerlinNoise1D(10.0, 5, 0.8);
				int q =   (floor(randomForce));
				particle->force[particle->nParticles - 1].x = ((randomForce - q * 1) - 0.5) * 2 * particle->mass;
		}
		else if(forceMode == 0)
		{
			particle->force[particle->nParticles - 1].x = particle->mass * 1.5;
		}
	}

	delete(tempPoints);
	delete(tempVelocity);
	delete(tempForce);

	return 0;
}


int InitializeWorld(struct world* particle)
{
	particle->nParticles = NO_OF_PARTICLES;

	particle->p = new point[particle->nParticles];
	particle->v = new point[particle->nParticles];
	particle->force = new point[particle->nParticles];
	underWaterConstant =  0.05;
	particle->dt = 0.01;
	particle->integrator = 0;
	particle->mass = 1.0 / particle->nParticles;
	particle->n = 5;
	particle->alpha = 0.2;
	particle->beta = 0.5;

	//particle->alpha = 0;
	//particle->beta = 0;

	particle->disk.center.x = 0;
	particle->disk.center.y = -0.5;
	particle->disk.center.z = 0;

	particle->disk.radius = DISK_RADIUS;
	particle->disk.no_of_segments = 200;
	particle->disk.p = NULL;

	InitializeDisk(&particle->disk);

	double particleDist = (particle->disk.radius * 2) / (particle->nParticles - 1); 
	

	for(int i = 0; i < particle->nParticles; i++)
	{
		particle->p[i].x = 0;
		particle->p[i].y = -(i * particleDist);
		particle->p[i].z = 0;
		particle->v[i].x = particle->v[i].y = particle->v[i].z = 0;

		particle->force[i].x =  particle->force[i].z = 0;
		particle->force[i].y = particle->mass * (-1);
	}


	if(forceMode == 1)
	{
		for(int i = 0; i < particle->nParticles; i++)
		{
			particle->force[i].y = particle->mass * (-1);
			particle->force[i].z = particle->force[i].x = 0;
		}
		double randomForce = PerlinNoise::PerlinNoise1D(10.0, 5, 0.8);
		int q =   (floor(randomForce));
		particle->force[particle->nParticles - 1].x = ((randomForce - q * 1) - 0.5) * 2 * particle->mass;
	}
	return 0;
}

void DeinitializeWorld(struct world* particle)
{
	Deinitializedisk(&particle->disk);
	free(particle->p);
	free(particle->v);
	free(particle->force);
}

int main (int argc, char ** argv)
{

	LARGE_INTEGER li;
    if(!QueryPerformanceFrequency(&li))
	printf("QueryPerformanceFrequency failed!\n");
	PCFreq = double(li.QuadPart)/1000.0;
	QueryPerformanceCounter(&li);
    CounterStart = li.QuadPart;
	InitialTime = CounterStart;

  InitializeWorld(&particle);

  glutInit(&argc,argv);
  
  /* double buffered window, use depth testing, 640x480 */
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  
  windowWidth = 640;
  windowHeight = 480;
  glutInitWindowSize (windowWidth, windowHeight);
  glutInitWindowPosition (0,0);
  glutCreateWindow ("Particle Constraint System");

  /* tells glut to use a particular display function to redraw */
  glutDisplayFunc(display);

  /* replace with any animate code */
  glutIdleFunc(doIdle);

  /* callback for mouse drags */
  glutMotionFunc(mouseMotionDrag);

  /* callback for window size changes */
  glutReshapeFunc(reshape);

  /* callback for mouse movement */
  glutPassiveMotionFunc(mouseMotion);

  /* callback for mouse button changes */
  glutMouseFunc(mouseButton);

  /* register for keyboard events */
  glutKeyboardFunc(keyboardFunc);

  /* do initialization */
  myinit();

  /* forever sink in the black hole */
  glutMainLoop();

  DeinitializeWorld(&particle);
  return(0);
}

