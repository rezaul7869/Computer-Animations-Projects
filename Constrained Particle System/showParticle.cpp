/*

  USC/Viterbi/Computer Science
  CONSTRAINED PARTICLE SYSTEM
  Using the -
  "Jello Cube" Assignment 1 starter code

*/

#include "particle.h"
#include "showParticle.h"
#include "texture.h"

GLuint theTorus;


void drawSphere(point center, double r, int lat_div, int long_div) {
	int i, j;
	for(i = 0; i <= lat_div; i++) 
	{
		double lat0 = 2 * pi * (-0.5 + (double) (i - 1) / lat_div);
		
		double lat1 = 2 * pi * (-0.5 + (double) i / lat_div);
		
		glBegin(GL_QUAD_STRIP);
		for(j = 0; j <= long_div; j++) 
		{
			double lng = 2 * pi * (double) (j - 1) / long_div;
			
			glNormal3f(r * (cos(lng) * cos(lat0))  + center.x, r * (sin(lng) * cos(lat0))  + center.y, r * sin(lat0) + center.z);
			glVertex3f(r * (cos(lng) * cos(lat0))  + center.x, r * (sin(lng) * cos(lat0))  + center.y, r * sin(lat0) + center.z);
			glNormal3f(r * (cos(lng) * cos(lat1))  + center.x, r * (sin(lng) * cos(lat1)) +  + center.y, r * sin(lat1) + center.z);
			glVertex3f(r * (cos(lng) * cos(lat1))  + center.x, r * (sin(lng) * cos(lat1)) +  + center.y, r * sin(lat1) + center.z);
		}
		glEnd();
    }
}

void drawCylinder(point a, point b, double R, int resolution)
{
	double length =  (a.x - b.x)*(a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.x - b.x) * (a.z - b.z);
	double two_pi = 2 * pi;

	glPushMatrix();
	glDisable(GL_CULL_FACE);
	glEnable(GL_COLOR_MATERIAL);
	glTranslatef(-a.x, -a.y, -a.z);
	glRotatef(-90, 1, 0, 0);
	glColor4f(0.7,0.0,0.0,0.9);
	glBegin(GL_QUAD_STRIP);
	for(int i = 0; i < resolution; i++)
	{
		glVertex3f(R * cos( i * two_pi/resolution),R * sin(i * two_pi/resolution), 0.0);
		glVertex3f(R * cos(i * two_pi/resolution),R * sin(i* two_pi/resolution), sqrt(length));
	}
	glEnd();
	glDisable(GL_CULL_FACE);
	glDisable(GL_COLOR_MATERIAL);
	glPopMatrix();
}

// create a quad strip for torus with R = radius from center of torus to center of tube
// r = radius of tube
void torus(point center, double R, double r,int numQuads, int numComponents)
{
   int i, j, k;
   double s, t, x, y, z, two_pi;
   glEnable(GL_COLOR_MATERIAL);
   two_pi = 2 * (double)pi;
   for (i = 0; i < numQuads; i++) {
      glBegin(GL_QUAD_STRIP);
      for (j = 0; j <= numComponents; j++) {
		  if(((j) / 10) % 2   == 0)
			  glColor4f(0.7,0,0,0.9);
		  else
			  glColor4f(1,1,1,0.9);
         for (k = 1; k >= 0; k--) {
            s = (i + k) % numQuads + 0.5;
            t = j % numComponents;

			x = (R + r * cos(s * two_pi / numQuads)) * cos(t * two_pi / numComponents) + center.x;
			y = (R + r * cos(s * two_pi/numQuads)) * sin(t * two_pi / numComponents) + center.y;
            z = r * sin(s * two_pi / numQuads)  + center.z;
            glVertex3f(x, y, z);
         }
      }
      glEnd();
   }
   glDisable(GL_COLOR_MATERIAL);
}


void initTorus(struct Disk disk)
{
   theTorus = glGenLists (1);
   glNewList(theTorus, GL_COMPILE);
   torus(disk.center, (disk.radius + 0.05), 0.05, 30, 200);
   glEndList();

   point p1 = {disk.center.x, disk.center.y + disk.radius - 0.025, disk.center.z};
   point p2 = {-8,-8, disk.center.z};
   drawCylinder(p1,p2, 0.025, 50);
}


void showDisc(struct Disk disk)
{
	glLineWidth(12);
	glPointSize(2);
	glBegin(GL_LINES);
	for(int i =0; i < disk.no_of_segments; i++)
	{
		glVertex3f(disk.p[i].x,disk.p[i].y,disk.p[i].z);\
		glVertex3f(disk.p[i+1].x,disk.p[i+1].y,disk.p[i+1].z);\
	}
	glEnd();

}

void showParticles(struct world * particle)
{
	double particleRadius =  0.03;
	int particleDivs = 12;


	#define PROCESS_NEIGHBOUR(i, di) \
    ip=i+(di);\
    if\
	((i >= 0) || (ip< particle->nParticles))\
    {\
		glBegin(GL_LINES);\
		glVertex3f(particle->p[i].x,particle->p[i].y,particle->p[i].z);\
      glVertex3f(particle->p[ip].x,particle->p[ip].y,particle->p[ip].z);\
	  glEnd();\
    }\

	glLineWidth(2);
    glPointSize(5);
    //glDisable(GL_LIGHTING);

	GLfloat white[] = {0.8f, 0.8f, 0.8f, 1.0f};
	GLfloat cyan[] = {0.f, .8f, .8f, 1.f};
	glMaterialfv(GL_FRONT, GL_DIFFUSE, cyan);
	glMaterialfv(GL_FRONT, GL_SPECULAR, white);
	GLfloat shininess[] = {50};
	glMaterialfv(GL_FRONT, GL_SHININESS, shininess);

	glColor4f(0,1.0,0,0);  
		for(int i = 0; i < particle->nParticles; i++)
		{
			drawSphere(particle->p[i], particleRadius, particleDivs, particleDivs);
		}

	glEnd();

	glDisable(GL_LIGHTING);
	glColor4f(0,0,1.0,0);
	for(int i = 0; i < particle->nParticles - 1; i++)
	{
		int ip;
		PROCESS_NEIGHBOUR(i, 1);       
	}
	glEnable(GL_LIGHTING);
}

void showRoom()
{
	glDisable(GL_CULL_FACE);
	glPushMatrix();

#define DRAW_VERTEX(x,y,z, u,v, r,g,b)\
	if(textured == 1)\
		glTexCoord2f (u,v);\
	else\
		glColor3f(r,g,b);\
	glVertex3f(x,y,z);

	int textured = 1;
	int texId = -1;
	if(texId = loadBMP("tex.bmp") == -1)
		textured = 0;
	glEnable(GL_LIGHTING);
	glColor3f(1.0,1.0,1.0);
	glEnable(GL_TEXTURE_2D);							// Enable Texture Mapping ( NEW )
	glBegin(GL_QUADS);
	
	
	DRAW_VERTEX(-8,-8,-15, 0.0, 0.0, 0.0, 0.5, 0.2 );
	DRAW_VERTEX(8,-8,-15, 1.0, 0.0, 0.0, 0.5, 0.2 );
	DRAW_VERTEX(8,-8,15, 1.0, 1.0, 0.0, 0.5, 0.2 );
	DRAW_VERTEX(-8,-8,15, 0.0, 1.0, 0.0, 0.5, 0.2 );
	
	/* Top */

	DRAW_VERTEX(-8,8,-15, 0.0, 0.0, 0.6, 0.1, 0.2 );
	DRAW_VERTEX(8,8,-15, 1.0, 0.0, 0.6, 0.1,0.2 );
	DRAW_VERTEX(8,8,15, 1.0, 1.0, 0.6, 0.1, 0.2 );
	DRAW_VERTEX(-8,8,15, 0.0, 1.0, 0.6, 0.1, 0.2 );

	
	/* Side Walls */

	DRAW_VERTEX(-8,-8,-15, 0.0, 0.0, 0.0, 1.0, 0.1 );
	DRAW_VERTEX(8,-8,-15, 1.0, 0.0, 0.0, 1.0, 0.1 );
	DRAW_VERTEX(8,8,-15, 1.0, 1.0, 0.0, 1.0, 0.1 );
	DRAW_VERTEX(-8,8,-15, 0.0, 1.0, 0.0, 1.0, 0.1 );

	
	DRAW_VERTEX(8,8,15, 0.0, 0.0, 0.0, 0.9, 0.4 );
	DRAW_VERTEX(8,-8,15, 1.0, 0.0, 0.0, 0.9, 0.4 );
	DRAW_VERTEX(8,-8,-15, 1.0, 1.0, 0.0, 0.9, 0.4 );
	DRAW_VERTEX(8,8,-15, 0.0, 1.0, 0.0, 0.9, 0.4 );
	
	DRAW_VERTEX(-8,8,15, 0.0, 0.0, 0.1, 1.0, 0.2 );
	DRAW_VERTEX(-8,-8,15, 1.0, 0.0, 0.1, 1.0, 0.2 );
	DRAW_VERTEX(-8,-8,-15, 1.0, 1.0, 0.1, 1.0, 0.2 );
	DRAW_VERTEX(-8,8,-15, 0.0, 1.0, 0.1, 1.0, 0.2 );

	DRAW_VERTEX(-8,-8,15, 0.0, 0.0, 0.7, 0.0, 0.2 );
	DRAW_VERTEX(8,-8,15, 1.0, 0.0, 0.7, 0.0, 0.2 );
	DRAW_VERTEX(8,8,15, 1.0, 1.0, 0.7, 0.0, 0.2 );
	DRAW_VERTEX(-8,8,15, 0.0, 1.0, 0.7, 0.0, 0.2 );

	glEnd();

	glPopMatrix();
	glEnable(GL_CULL_FACE);
}