#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "motion.h"
#include "interpolator.h"
#include "types.h"


#define PRECISION 10000000

Interpolator::Interpolator()
{
  //Set default interpolation type
  m_InterpolationType = LINEAR;

  //set default angle representation to use for interpolation
  m_AngleRepresentation = EULER;
}

Interpolator::~Interpolator()
{
}

//Create interpolated motion
void Interpolator::Interpolate(Motion * pInputMotion, Motion ** pOutputMotion, int N) 
{
  //Allocate new motion
  *pOutputMotion = new Motion(pInputMotion->GetNumFrames(), pInputMotion->GetSkeleton()); 

  //Perform the interpolation
  if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == EULER))
    LinearInterpolationEuler(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == LINEAR) && (m_AngleRepresentation == QUATERNION))
    LinearInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == EULER))
    BezierInterpolationEuler(pInputMotion, *pOutputMotion, N);
  else if ((m_InterpolationType == BEZIER) && (m_AngleRepresentation == QUATERNION))
    BezierInterpolationQuaternion(pInputMotion, *pOutputMotion, N);
  else
  {
    printf("Error: unknown interpolation / angle representation type.\n");
    exit(1);
  }
}

void Interpolator::LinearInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  int startKeyframe = 0;

  LARGE_INTEGER li;
  if(!QueryPerformanceFrequency(&li))
		printf("QueryPerformanceFrequency failed!\n");
  double PCFreq = double(li.QuadPart)/1000.0;
  double totalTime = 0;
  int ncases = 0;

  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
	  QueryPerformanceCounter(&li);
	  double CounterStart = li.QuadPart;

      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);

      // interpolate root position
      interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
        interpolatedPosture.bone_rotation[bone] = startPosture->bone_rotation[bone] * (1-t) + endPosture->bone_rotation[bone] * t;

      pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);

	  LARGE_INTEGER li;
	  QueryPerformanceCounter(&li);
      double timeTaken = (li.QuadPart-CounterStart)/PCFreq;
	  totalTime += timeTaken;
	  ncases++;
    }

    startKeyframe = endKeyframe;
  }

  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));

  double avgTime = totalTime / ncases;
  printf("Average computation time per frame is %f milliseconds\n", avgTime);
}

void Interpolator::Rotation2Euler(double R[9], double angles[3])
{
  double cy = sqrt(R[0]*R[0] + R[3]*R[3]);

  if (cy > 16*DBL_EPSILON) 
  {
    angles[0] = atan2(R[7], R[8]);
    angles[1] = atan2(-R[6], cy);
    angles[2] = atan2(R[3], R[0]);
  } 
  else 
  {
    angles[0] = atan2(-R[5], R[4]);
    angles[1] = atan2(-R[6], cy);
    angles[2] = 0;
  }

  for(int i=0; i<3; i++)
    angles[i] *= 180 / M_PI;
}

void Interpolator::Euler2Rotation(double angles[3], double R[9])
{
  // Change the Euler angles to rotation matrix using standard conversion rules
	double A[3];

	// Convert degree to radians
	for(int i=0; i<3; i++)
		A[i] = angles[i] * M_PI / 180;

	R[0] = cos(A[1]) * cos(A[2]);
	R[1] = sin(A[0]) * sin(A[1]) * cos(A[2]) - cos(A[0]) * sin(A[2]);
	R[2] = cos(A[0]) * sin(A[1]) * cos(A[2]) + sin(A[0]) * sin(A[2]);

	R[3] = cos(A[1]) * sin(A[2]);
	R[4] = sin(A[0]) * sin(A[1]) * sin(A[2]) + cos(A[0]) * cos(A[2]);
	R[5] = cos(A[0]) * sin(A[1]) * sin(A[2]) - sin(A[0]) * cos(A[2]);

	R[6] = -sin(A[1]);
	R[7] = sin(A[0]) * cos(A[1]);
	R[8] = cos(A[0]) * cos(A[1]);
}

void Interpolator::BezierInterpolationEuler(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
  // Interpolate Euler angles by Bezier interpolations
	int inputLength = pInputMotion->GetNumFrames();

  int startKeyframe = 0;

  /* Initialize variables for time checking*/
  LARGE_INTEGER li;
  if(!QueryPerformanceFrequency(&li))
		printf("QueryPerformanceFrequency failed!\n");
  double PCFreq = double(li.QuadPart)/1000.0;
  double totalTime = 0;
  int ncases = 0;

  while (startKeyframe + N + 1 < inputLength)
  {
	  //Considering the interpolation is between q_n and q_n+1
    int endKeyframe = startKeyframe + N + 1;	// q_n
	int nextKeyframe = endKeyframe + N + 1;		// q_n+2
	int prevKeyframe = startKeyframe - N - 1;	// q_n-1


    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

	Posture *prevPosture, *nextPosture;
	if(nextKeyframe < inputLength)
	{
		nextPosture = pInputMotion->GetPosture(nextKeyframe);
	}

	if(prevKeyframe >= 0)
	{
		prevPosture = pInputMotion->GetPosture(prevKeyframe);
	}

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);
    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
		QueryPerformanceCounter(&li);
		double CounterStart = li.QuadPart;

		Posture interpolatedPosture;
		double t = 1.0 * frame / (N+1);

		// interpolate root position
		vector p1 = startPosture->root_pos;
		vector p2 = endPosture->root_pos;
		vector m1, n2; 

		/*For the last part of interpolation calculating bn*/
		if(nextKeyframe >= inputLength)
		{
			if(prevKeyframe < 0)
				printf("Error : No previous frame and no next frame");
			vector p0 = prevPosture->root_pos;
			vector m_ = p1 * 2 - p0;
			double u = (1.0/3);
			n2 = p2 * (1-u) + m_ * u;
		}
		else
		{
			vector p3 = nextPosture->root_pos;
			vector m__ = p2 * 2 - p1;
			vector m_= p3 * 0.5 + m__ * 0.5;
			double u = (-1.0/3);
			n2 = p2 * (1-u) + m_ * u;
		}

		/*For the first part of interpolation calculating a0*/
		if(prevKeyframe < 0)
		{
			if(nextKeyframe >= inputLength)
				printf("Error : No previous frame and no next frame");
			vector p3 = nextPosture->root_pos;
			vector m_ = p2 * 2 - p3;
			double u = (-1.0/3);
			m1 = p1 * (1-u) + m_ * u;
		}
		else
		{
			vector p0 = prevPosture->root_pos;
			vector m__ = p1 * 2 - p0;
			vector m_= p2 * 0.5 + m__ * 0.5;
			double u = (1.0/3);
			m1 = p1 * (1-u) + m_ * u;
		}
		interpolatedPosture.root_pos = DeCasteljauEuler(t, p1, m1, n2, p2);

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
	  {
		    
		    vector q1 =  startPosture->bone_rotation[bone];
		    vector q2 =  endPosture->bone_rotation[bone];
			vector a1, b2;

			/*For the last part of interpolation calculating bn*/
			if(nextKeyframe >= inputLength)
			{
				if(prevKeyframe < 0)
					printf("Error : No previous frame and no next frame");
				vector q0 = prevPosture->bone_rotation[bone];
				vector a_ = q1 * 2 - q0;
				double u = (1.0/3);
				b2 = q2 * (1-u) + a_ * u;
			}
			else
			{
				vector q3 = nextPosture->bone_rotation[bone];
				vector a__ = q2 * 2 - q1;
				vector a_= q3 * 0.5 + a__ * 0.5;
				double u = (-1.0/3);
				b2 = q2 * (1-u) + a_ * u;
			}
			/*For the first part of interpolation calculating a0*/
			if(prevKeyframe < 0)
			{
				if(nextKeyframe >= inputLength)
					printf("Error : No previous frame and no next frame");
				vector q3 = nextPosture->bone_rotation[bone];
				vector a_ = q3 * 2 - q2;
				double u = (-1.0/3);
				a1 = q1 * (1-u) + a_ * u;
			}
			else
			{
				vector q0 = prevPosture->bone_rotation[bone];
				vector a__ = q1 * 2 - q0;
				vector a_= q2 * 0.5 + a__ * 0.5;
				double u = (1.0/3);
				a1 = q1 * (1-u) + a_ * u;
			}
		
			interpolatedPosture.bone_rotation[bone] = DeCasteljauEuler(t, q1, a1, b2, q2);

	  }

	  pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);

	  /*Compute average time taken*/
      LARGE_INTEGER li;
	  QueryPerformanceCounter(&li);
      double timeTaken = (li.QuadPart-CounterStart)/PCFreq;
	  totalTime += timeTaken;
	  ncases++;
    }
    startKeyframe = endKeyframe;
  }

  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));

  double avgTime = totalTime / ncases;
  printf("Average computation time per frame is %f milliseconds\n", avgTime);
}

void Interpolator::LinearInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
	// Interpolate quaternion by SLERP interpolations
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  /* Initialize variables for time checking*/
  LARGE_INTEGER li;
  if(!QueryPerformanceFrequency(&li))
		printf("QueryPerformanceFrequency failed!\n");
  double PCFreq = double(li.QuadPart)/1000.0;
  double totalTime = 0;
  int ncases = 0;


  int startKeyframe = 0;
  while (startKeyframe + N + 1 < inputLength)
  {
    int endKeyframe = startKeyframe + N + 1;

    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
		/*Calculate time for each frame*/
	  QueryPerformanceCounter(&li);
	  double CounterStart = li.QuadPart;
      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);

      // interpolate root position
      interpolatedPosture.root_pos = startPosture->root_pos * (1-t) + endPosture->root_pos * t;

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
	  {
		double rotation1[3] = {0}, rotation2[3]= {0};
		double temprotation1[3] = {0}, temprotation2[3]= {0};
		startPosture->bone_rotation[bone].getValue(rotation1);
		endPosture->bone_rotation[bone].getValue(rotation2);
		Quaternion<double> p,q;
		double tAngles[3] = {0};

		/*Get quaternion for each Euler rotation*/
		Euler2Quaternion(rotation1, p);
		Euler2Quaternion(rotation2, q);
		
		/*Do SLERP*/
		Quaternion<double> interpolatedP = Slerp(t, p, q);

		/*Get Euler rotation for the interpolated quaternion*/
		Quaternion2Euler(interpolatedP, tAngles);

		interpolatedPosture.bone_rotation[bone].setValue(tAngles);
	  }
      
	  pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);

	  /*Compute average time taken*/
	  LARGE_INTEGER li;
	  QueryPerformanceCounter(&li);
      double timeTaken = (li.QuadPart-CounterStart)/PCFreq;
	  totalTime += timeTaken;
	  ncases++;
    }

    startKeyframe = endKeyframe;
  }

  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));

  double avgTime = totalTime / ncases;
  printf("Average computation time per frame is %f milliseconds\n", avgTime);
}

void Interpolator::BezierInterpolationQuaternion(Motion * pInputMotion, Motion * pOutputMotion, int N)
{
	// Interpolate Quaternion by Bezier SLERP interpolations
  int inputLength = pInputMotion->GetNumFrames(); // frames are indexed 0, ..., inputLength-1

  int startKeyframe = 0;
  int prevKeyframe = 0;
  int nextKeyframe = 0;

  /* Initialize variables for time checking*/
  LARGE_INTEGER li;
  if(!QueryPerformanceFrequency(&li))
		printf("QueryPerformanceFrequency failed!\n");
  double PCFreq = double(li.QuadPart)/1000.0;
  double totalTime = 0;
  int ncases = 0;

  while (startKeyframe + N + 1 < inputLength)
  {
	  //Considering the interpolation is between q_n and q_n+1
    int endKeyframe = startKeyframe + N + 1; // q_n
	int nextKeyframe = endKeyframe + N + 1;		//q_n+2
	int prevKeyframe = startKeyframe - N - 1;	//q_n-1


    Posture * startPosture = pInputMotion->GetPosture(startKeyframe);
    Posture * endPosture = pInputMotion->GetPosture(endKeyframe);

	Posture *prevPosture, *nextPosture;
	if(nextKeyframe < inputLength)
	{
		nextPosture = pInputMotion->GetPosture(nextKeyframe);
	}

	if(prevKeyframe >= 0)
	{
		prevPosture = pInputMotion->GetPosture(prevKeyframe);
	}

    // copy start and end keyframe
    pOutputMotion->SetPosture(startKeyframe, *startPosture);
    pOutputMotion->SetPosture(endKeyframe, *endPosture);

    // interpolate in between
    for(int frame=1; frame<=N; frame++)
    {
		/*Calculate time for each frame*/
	  QueryPerformanceCounter(&li);
	  double CounterStart = li.QuadPart;

      Posture interpolatedPosture;
      double t = 1.0 * frame / (N+1);

      // interpolate root position
        vector p1 = startPosture->root_pos;
		vector p2 = endPosture->root_pos;
		vector m1, n2; 

		/*For the last part of interpolation calculating bn*/
		if(nextKeyframe >= inputLength)
		{
			if(prevKeyframe < 0)
				printf("Error : No previous frame and no next frame");
			vector p0 = prevPosture->root_pos;
			vector m_ = p1 * 2 - p0;
			double u = (1.0/3);
			n2 = p2 * (1-u) + m_ * u;
		}
		else
		{
			vector p3 = nextPosture->root_pos;
			vector m__ = p2 * 2 - p1;
			vector m_= p3 * 0.5 + m__ * 0.5;
			double u = (-1.0/3);
			n2 = p2 * (1-u) + m_ * u;
		}

		/*For the first part of interpolation calculating a0*/
		if(prevKeyframe < 0)
		{
			if(nextKeyframe >= inputLength)
				printf("Error : No previous frame and no next frame");
			vector p3 = nextPosture->root_pos;
			vector m_ = p2 * 2 - p3;
			double u = (-1.0/3);
			m1 = p1 * (1-u) + m_ * u;
		}
		else
		{
			vector p0 = prevPosture->root_pos;
			vector m__ = p1 * 2 - p0;
			vector m_= p2 * 0.5 + m__ * 0.5;
			double u = (1.0/3);
			m1 = p1 * (1-u) + m_ * u;
		}
		interpolatedPosture.root_pos = DeCasteljauEuler(t, p1, m1, n2, p2);

      // interpolate bone rotations
      for (int bone = 0; bone < MAX_BONES_IN_ASF_FILE; bone++)
	  {
		double rotation1[3] = {0}, rotation2[3]= {0};
		double rotation0[3] = {0}, rotation3[3] ={0};

		startPosture->bone_rotation[bone].getValue(rotation1);
		endPosture->bone_rotation[bone].getValue(rotation2);
		Quaternion<double> q1, a1, b2, q2;
		double tAngles[3] = {0};
		Euler2Quaternion(rotation1, q1);
		Euler2Quaternion(rotation2, q2);

		/*For the last part of interpolation calculating bn*/
		if(nextKeyframe >= inputLength)
		{
			if(prevKeyframe < 0)
				printf("Error : No previous frame and no next frame");
			prevPosture->bone_rotation[bone].getValue(rotation0);
			Quaternion<double> q0;
			Euler2Quaternion(rotation0, q0);
			b2 = Slerp((1.0/3), q2, Double(q0, q1));
		}
		else
		{
			nextPosture->bone_rotation[bone].getValue(rotation3);
			Quaternion<double> q3;
			Euler2Quaternion(rotation3, q3);
			Quaternion<double> a_=Slerp(0.5, Double(q1, q2), q3);
			b2 = Slerp((-1.0/3),q2,a_);
		}

		/*For the first part of interpolation calculating a0*/
		if(prevKeyframe < 0)
		{
			if(nextKeyframe >= inputLength)
				printf("Error : No previous frame and no next frame");
			nextPosture->bone_rotation[bone].getValue(rotation3);
			Quaternion<double> q3;
			Euler2Quaternion(rotation3, q3);
			a1 = Slerp((1.0/3),q1,Double(q2, q3));
		}
		else
		{
			prevPosture->bone_rotation[bone].getValue(rotation0);
			Quaternion<double> q0;
			Euler2Quaternion(rotation0, q0);
			Quaternion<double> a_=Slerp(0.5, Double(q0, q1), q2);
			a1 = Slerp((1.0/3),q1,a_);
		}
		
		Quaternion<double> interpolatedP = DeCasteljauQuaternion(t, q1, a1, b2, q2);
		
		Quaternion2Euler(interpolatedP, tAngles);

		interpolatedPosture.bone_rotation[bone].setValue(tAngles);
	  }

	  pOutputMotion->SetPosture(startKeyframe + frame, interpolatedPosture);
      
	  /*Compute average time taken*/
	  LARGE_INTEGER li;
	  QueryPerformanceCounter(&li);
      double timeTaken = (li.QuadPart-CounterStart)/PCFreq;
	  totalTime += timeTaken;
	  ncases++;
    }

    startKeyframe = endKeyframe;
  }

  for(int frame=startKeyframe+1; frame<inputLength; frame++)
    pOutputMotion->SetPosture(frame, *(pInputMotion->GetPosture(frame)));

  double avgTime = totalTime / ncases;
  printf("Average computation time per frame is %f milliseconds\n", avgTime);
}

void Interpolator::Euler2Quaternion(double angles[3], Quaternion<double> & q) 
{
  // DONE : students should implement this

	double A[3];

	for(int i=0; i<3; i++)
		A[i] = angles[i] * M_PI / 180;

	double cos0 = cos(A[0]/2);
	double cos1 = cos(A[1]/2);
	double cos2 = cos(A[2]/2);
	double sin0 = sin(A[0]/2);
	double sin1 = sin(A[1]/2);
	double sin2 = sin(A[2]/2);

	double q0 = cos0 * cos1 * cos2 + sin0 * sin1 * sin2;
	double q1 = sin0 * cos1 * cos2 - cos0 * sin1 * sin2;
	double q2 = cos0 * sin1 * cos2 + sin0 * cos1 * sin2;
	double q3 = cos0 * cos1 * sin2 - sin0 * sin1 * cos2;
	
	q.Set(q0, q1, q2, q3);
	

	/*
	double R[9];
	Euler2Rotation(angles, R);
	q = Quaternion<double>::Matrix2Quaternion(R);
	*/
}

void Interpolator::Quaternion2Euler(Quaternion<double> & q, double angles[3]) 
{
  // DONE : students should implement this
	double q0 = q.Gets();
	double q1 = q.Getx();
	double q2 = q.Gety();
	double q3 = q.Getz();
	angles[0] =  atan2( (2*(q0 * q1 + q2 * q3)) , (1- (2 * (q1 * q1 + q2 * q2))) );
	angles[1] =  asin(2*(q0 * q2 - q3 * q1));
	angles[2] =  atan2( (2*(q0 * q3 + q1 * q2)) , (1- (2 * (q2 * q2 + q3 * q3))) );

	for(int i=0; i<3; i++)
		angles[i] *= 180 / M_PI;
	
	//double R[9];
	//q.Quaternion2Matrix(R);
	//Rotation2Euler(R,angles);
}

Quaternion<double> Interpolator::Slerp(double t, Quaternion<double> & qStart, Quaternion<double> & qEnd_)
{
  // DONE : students should implement this
  Quaternion<double> result;

  double dotProd = qStart.Gets() * qEnd_.Gets() +
						qStart.Getx() * qEnd_.Getx() +
						qStart.Gety() * qEnd_.Gety() +
						qStart.Getz() * qEnd_.Getz();

  double angle =  acos((floor(dotProd*PRECISION))/PRECISION);
  
  if(angle == 0)
  {
	result = ((1-t) * qStart) + (t * qEnd_);
  }
  else
	result = ((sin((1-t) * angle) / sin(angle)) * qStart) + ((sin((t) * angle) / sin(angle)) * qEnd_);

  result.Normalize();
  
  return result;
}

Quaternion<double> Interpolator::Double(Quaternion<double> p, Quaternion<double> q)
{
  // DONE : students should implement this
  Quaternion<double> result;

   double dotProd =  (p.Gets() * q.Gets() +
						p.Getx() * q.Getx() +
						p.Gety() * q.Gety() +
						p.Getz() * q.Getz());

  result =  (2 * dotProd) * q - p;
  result.Normalize();
  return result;
}

vector Interpolator::DeCasteljauEuler(double t, vector p0, vector p1, vector p2, vector p3)
{
  // students should implement this
  vector result;
  
  vector q0 = p0 * (1-t) + p1 * t;
  vector q1 = p1 * (1-t) + p2 * t;
  vector q2 = p2 * (1-t) + p3 * t;

  vector r0 = q0 * (1-t) + q1 * t;
  vector r1 = q1 * (1-t) + q2 * t;

  result = r0 * (1-t) + r1 * t;
  return result;
}

Quaternion<double> Interpolator::DeCasteljauQuaternion(double t, Quaternion<double> p0, Quaternion<double> p1, Quaternion<double> p2, Quaternion<double> p3)
{
  // students should implement this
  Quaternion<double> result;

  Quaternion<double>  q0 = Slerp(t, p0, p1);
  Quaternion<double>  q1 = Slerp(t, p1, p2);
  Quaternion<double>  q2 = Slerp(t, p2, p3);

  Quaternion<double>  r0 = Slerp(t, q0, q1);
  Quaternion<double>  r1 = Slerp(t, q1, q2);

  result = Slerp(t, r0, r1);
  return result;
}

