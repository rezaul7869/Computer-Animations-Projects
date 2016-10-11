#ifndef _PERLINNOISE_H_
#define _PERLINNOISE_H_

/*

  Jernej Barbic, MIT, 2007-2009

  Perlin Noise

  This C++ implementation has been created using the material available at the following webpage:
  http://freespace.virgin.net/hugo.elias/models/m_perlin.htm

*/

class PerlinNoise
{
public:
  // "numFrequencies" determines how much spatial variability there will be in the noise
  // "persistence" should be on the interval (0,1]. It measures how quickly high frequencies decay.
  static double PerlinNoise1D(double x, int numFrequencies, double persistence);

protected:
  static double Noise(int x);
  static double SmoothedNoise(int x);
  static double InterpolatedNoise(double x);
  static double CosineInterpolate(double a, double b, double x);
};

#endif

