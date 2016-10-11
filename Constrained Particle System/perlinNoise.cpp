#include <math.h>
#include "perlinNoise.h"
#include "particle.h"

double PerlinNoise::Noise(int x)
{
  x = (x<<13) ^ x;
  return ( 1.0 - ( (x * (x * x * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0);    
}

double PerlinNoise::SmoothedNoise(int x)
{
  return (Noise(x)/2  +  Noise(x-1)/4  +  Noise(x+1)/4);
}

double PerlinNoise::InterpolatedNoise(double x)
{
  int integer_X = int(x);
  double fractional_X = x - integer_X;

  double v1 = SmoothedNoise(integer_X);
  double v2 = SmoothedNoise(integer_X + 1);

  return CosineInterpolate(v1 , v2 , fractional_X);
}

double PerlinNoise::CosineInterpolate(double a, double b, double x)
{
  double ft = x * pi;
  double f = 0.5 * (1 - cos(ft));
  return  a*(1-f) + b*f;
}

double PerlinNoise::PerlinNoise1D(double x, int numFrequencies, double persistence)
{
  double noise = 0.0;
  double frequency = 1.0;
  double amplitude = persistence;
  for(int i=0; i<numFrequencies - 1; i++)
  {
    noise += InterpolatedNoise(x * frequency) * amplitude;
    frequency *= 2;
    amplitude *= persistence;
  }

  return noise;
}

