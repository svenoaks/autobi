/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.cuny.qc.speech.AuToBI.util;

import edu.cuny.qc.speech.AuToBI.core.Pair;
import edu.cuny.qc.speech.AuToBI.core.Triple;
import edu.cuny.qc.speech.AuToBI.core.WavData;

/**
 *
 * @author vitisoto
 */
public class SoundUtils {
    
    
    private final static int Vector_VALUE_INTERPOLATION_NEAREST = 0;
    private final static int Vector_VALUE_INTERPOLATION_LINEAR = 1;
    private final static int Vector_VALUE_INTERPOLATION_CUBIC = 2;
    private final static int Vector_VALUE_INTERPOLATION_SINC70 = 3;
    private final static int Vector_VALUE_INTERPOLATION_SINC700 = 4;
    
 
    
    public static double getAbsoluteExtremum (WavData sound,double xmin, double xmax, int interpolation) {
        double minimum = Math.abs (getMinimum (sound, xmin, xmax, interpolation));
        double maximum = Math.abs (getMaximum (sound,xmin, xmax, interpolation));
        return minimum > maximum ? minimum : maximum;
    }
  
    public static double getMinimum (WavData sound,double xmin, double xmax, int interpolation) {
	Triple<Double, Double, Double> minimum = getMinimumAndXAndChannel (sound,xmin, xmax, interpolation);
        return minimum.first;
    }
  
    public static double getMaximum (WavData sound,double xmin, double xmax, int interpolation) {
	Triple<Double, Double, Double> maximum = getMaximumAndXAndChannel (sound,xmin, xmax, interpolation);
        return maximum.first;
    }
  
    public static Triple<Double,Double,Double> getMinimumAndXAndChannel (WavData sound,double xmin, double xmax, int interpolation){
	double minimum, xOfMinimum;
	long channelOfMinimum = 0;
	Pair<Double,Double> pair = getMinimumAndX (sound, xmin, xmax, 0, interpolation);
        minimum = pair.first;
        xOfMinimum = pair.second;
        
	for (int channel = 1; channel < sound.numberOfChannels ; channel ++) {
		double minimumOfChannel, xOfMinimumOfChannel;
		pair = getMinimumAndX (sound, xmin, xmax, channel, interpolation);
                minimumOfChannel = pair.first;
                xOfMinimumOfChannel = pair.second;
		if (minimumOfChannel < minimum) {
			minimum = minimumOfChannel;
			xOfMinimum = xOfMinimumOfChannel;
			channelOfMinimum = channel;
		}
	}
        
        Triple<Double,Double,Double> tri = new Triple(minimum, xOfMinimum, channelOfMinimum);
	
        return tri;
    }
  
    public static Triple<Double,Double,Double> getMaximumAndXAndChannel (WavData sound,double xmin, double xmax, int interpolation){
	double maximum, xOfMaximum;
	long channelOfMaximum = 0;
	Pair<Double,Double> pair = getMaximumAndX (sound, xmin, xmax, 0, interpolation);
        maximum = pair.first;
        xOfMaximum = pair.second;
        
	for (int channel = 1; channel < sound.numberOfChannels; channel ++) {
		double maximumOfChannel, xOfMaximumOfChannel;
		pair = getMaximumAndX (sound, xmin, xmax, channel, interpolation);
                maximumOfChannel = pair.first;
                xOfMaximumOfChannel = pair.second;
		if (maximumOfChannel > maximum) {
			maximum = maximumOfChannel;
			xOfMaximum = xOfMaximumOfChannel;
			channelOfMaximum = channel;
		}
	}
        
        Triple<Double,Double,Double> tri = new Triple(maximum, xOfMaximum, channelOfMaximum);
	
        return tri;
  }
  
  public static Pair<Double,Double> getMinimumAndX (WavData sound, double xmin, double xmax, int channel, int interpolation){
      int imin, imax, n = sound.getNumSamples();
      //Melder_assert (channel >= 1 && channel <= my ny);
      double[] y = sound.getSamples(channel);
      //double *y = my z [channel];
      double minimum, x;
      if (xmax <= xmin) { 
          System.err.println("xmin and xmax incorrect in SoundUtils.getMininmumAndX()");
      }
      
      Triple<Integer, Integer, Integer> tri = SampledUtils.getWindowSamples (sound.t0,sound.getFrameSize(),sound.getNumSamples(), xmin, xmax);
      imin = tri.first;imax = tri.second;
      if(tri.third==0.0) {
          /*
           * No samples between xmin and xmax.
           * Try to return the lesser of the values at these two points.
           */
          double yleft = NUMUtils.Vector_getValueAtX (sound, xmin, channel,interpolation > Vector_VALUE_INTERPOLATION_NEAREST ? Vector_VALUE_INTERPOLATION_LINEAR : Vector_VALUE_INTERPOLATION_NEAREST);
          double yright = NUMUtils.Vector_getValueAtX (sound, xmax, channel, interpolation > Vector_VALUE_INTERPOLATION_NEAREST ? Vector_VALUE_INTERPOLATION_LINEAR : Vector_VALUE_INTERPOLATION_NEAREST);
          minimum = yleft < yright ? yleft : yright;
          x = yleft == yright ? (xmin + xmax) / 2 : yleft < yright ? xmin : xmax;
                
      } else {
          minimum = y [imin]; x = imin;
          if (y [imax] < minimum){ minimum = y [imax]; x = imax;}
          if (imin == 0) imin ++;
          if (imax == y.length-1) imax --;
          for (int i = imin; i <= imax; i ++) {
              if (y [i] < y [i - 1] && y [i] <= y [i + 1]) {
                  
                  Pair<Double, Double> pair = NUMUtils.improveMinimum (y, n, i, interpolation);
                  double localMinimum = pair.first;
                  double i_real = pair.second;
                  
                  if (localMinimum < minimum){minimum = localMinimum; x = i_real;}
              }
          }
          x = sound.t0 + x*sound.getFrameSize();   /* Convert sample to x. */
          if (x < xmin) x = xmin; else if (x > xmax) x = xmax;
      }
      Pair<Double, Double> pair = new Pair(minimum,x);
      return pair;
  }
  
  public static Pair<Double,Double> getMaximumAndX (WavData sound, double xmin, double xmax, int channel, int interpolation){
      int imin, imax, n = sound.getNumSamples();
      //Melder_assert (channel >= 1 && channel <= my ny);
      double[] y = sound.getSamples(channel);
      //double *y = my z [channel];
      double maximum, x;
      //if (xmax <= xmin) { xmin = my xmin; xmax = my xmax; }
      Triple<Integer, Integer, Integer> tri = SampledUtils.getWindowSamples (sound.t0,sound.getFrameSize(),sound.getNumSamples(), xmin, xmax);
      imin = tri.first;imax = tri.second;
      if(tri.third==0.0) {
          /*
           * No samples between xmin and xmax.
           * Try to return the lesser of the values at these two points.
           */
          double yleft = NUMUtils.Vector_getValueAtX (sound, xmin, channel,interpolation > Vector_VALUE_INTERPOLATION_NEAREST ? Vector_VALUE_INTERPOLATION_LINEAR : Vector_VALUE_INTERPOLATION_NEAREST);
          double yright = NUMUtils.Vector_getValueAtX (sound, xmax, channel, interpolation > Vector_VALUE_INTERPOLATION_NEAREST ? Vector_VALUE_INTERPOLATION_LINEAR : Vector_VALUE_INTERPOLATION_NEAREST);
          maximum = yleft > yright ? yleft : yright;
          x = yleft == yright ? (xmin + xmax) / 2 : yleft > yright ? xmin : xmax;
                
      } else {
          maximum = y [imin]; x = imin;
          if (y [imax] > maximum){maximum = y [imax]; x = imax;}
          if (imin == 0) imin ++;
          if (imax == y.length-1) imax --;
          for (int i = imin; i <= imax; i ++) {
              if (y [i] > y [i - 1] && y [i] >= y [i + 1]) {
                  
                  Pair<Double, Double> pair = NUMUtils.improveMaximum (y, n, i, interpolation);
                  double localMaximum = pair.first;
                  double i_real = pair.second;
                  
                  if (localMaximum > maximum){maximum = localMaximum; x = i_real;}
              }
          }
          x = sound.t0 + x *sound.getFrameSize();   /* Convert sample to x. */
          if (x < xmin) x = xmin; else if (x > xmax) x = xmax;
          
      }
      Pair<Double, Double> pair = new Pair(maximum,x);
      return pair;
  }
  
  

    
}
