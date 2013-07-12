/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.cuny.qc.speech.AuToBI.util;

import edu.cuny.qc.speech.AuToBI.core.Contour;
import edu.cuny.qc.speech.AuToBI.core.Pair;
import edu.cuny.qc.speech.AuToBI.core.Triple;
import edu.cuny.qc.speech.AuToBI.core.WavData;

/**
 *
 * @author vitisoto
 */
public class NUMUtils {
    
    
    private final static int NUM_VALUE_INTERPOLATE_NEAREST=0;
    private final static int NUM_VALUE_INTERPOLATE_LINEAR =1;
    private final static int NUM_VALUE_INTERPOLATE_CUBIC=2;
    // Higher values than 2 yield a true sinc interpolation. Here are some examples:
    private final static int NUM_VALUE_INTERPOLATE_SINC70=70;
    private final static int NUM_VALUE_INTERPOLATE_SINC700=700;
    
    
    public final static int NUM_PEAK_INTERPOLATE_NONE=0;
    public final static int NUM_PEAK_INTERPOLATE_PARABOLIC=1;
    public final static int NUM_PEAK_INTERPOLATE_CUBIC=2;
    public final static int NUM_PEAK_INTERPOLATE_SINC70=3;
    public final static int NUM_PEAK_INTERPOLATE_SINC700=4;
    
    private final static double NUM_goldenSection=0.6180339887498948482045868343656381177203;
    
    private final static int Vector_CHANNEL_AVERAGE = 0;
    
    private final static int Vector_VALUE_INTERPOLATION_NEAREST = 0;
    private final static int Vector_VALUE_INTERPOLATION_LINEAR = 1;
    private final static int Vector_VALUE_INTERPOLATION_CUBIC = 2;
    private final static int Vector_VALUE_INTERPOLATION_SINC70 = 3;
    private final static int Vector_VALUE_INTERPOLATION_SINC700 = 4;
 
    
    public final static int Pitch_STRENGTH_UNIT_min = 0;
    public final static int Pitch_STRENGTH_UNIT_AUTOCORRELATION = 0;
    public final static int Pitch_STRENGTH_UNIT_NOISE_HARMONICS_RATIO = 1;
    public final static int Pitch_STRENGTH_UNIT_HARMONICS_NOISE_DB = 2;
    public final static int Pitch_STRENGTH_UNIT_max = 2;
    
    public final static double NUMundefined = Double.NaN;
    

    
    
    private static double calculateMachineEpsilonDouble() {
        double machEps = 1.0f;
 
        do {
           machEps /= 2.0f;
        }
        while ((double)(1.0 + (machEps/2.0)) != 1.0);
 
        return machEps;
    }    
    
    public static double interpolate_sinc (double y [], int nx, double x, int maxDepth) {
	int left, right;
        int ix, midleft = (int)Math.floor(x), midright = midleft + 1;
	double result = 0.0, a, halfsina, aa, daa;
	//NUM_interpolate_simple_cases
        if (nx < 1) return NUMundefined;
	if (x > nx) return y[nx];
	if (x < 1) return y[1];
	if (x == midleft) return y [midleft];
	/* 1 < x < nx && x not integer: interpolate. */
	if (maxDepth > midright - 1) maxDepth = midright - 1;
	if (maxDepth > nx - midleft) maxDepth = nx - midleft;
	if (maxDepth <= NUM_VALUE_INTERPOLATE_NEAREST) return y [(int) Math.floor (x + 0.5)];
	if (maxDepth == NUM_VALUE_INTERPOLATE_LINEAR) return y [midleft] + (x - midleft) * (y [midright] - y [midleft]); 
	if (maxDepth == NUM_VALUE_INTERPOLATE_CUBIC) { 
		double yl = y [midleft], yr = y [midright]; 
		double dyl = 0.5 * (yr - y [midleft - 1]), dyr = 0.5 * (y [midright + 1] - yl); 
		double fil = x - midleft, fir = midright - x; 
		return yl * fir + yr * fil - fil * fir * (0.5 * (dyr - dyl) + (fil - 0.5) * (dyl + dyr - 2 * (yr - yl))); 
	}
        
        
	left = midright - maxDepth;
        right = midleft + maxDepth;
	a = Math.PI * (x - midleft);
	halfsina = 0.5 * Math.sin(a);
	aa = a / (x - left + 1);
	daa = Math.PI / (x - left + 1);
	for (ix = midleft; ix >= left; ix --) {
		double d = halfsina / a * (1.0 + Math.cos(aa));
		result += y [ix] * d;
		a += Math.PI;
		aa += daa;
		halfsina = - halfsina;
	}
	a = Math.PI * (midright - x);
	halfsina = 0.5 * Math.sin(a);
	aa = a / (right - x + 1);
	daa = Math.PI / (right - x + 1);
	for (ix = midright; ix <= right; ix ++) {
		double d = halfsina / a * (1.0 + Math.cos (aa));
		result += y [ix] * d;
		a += Math.PI;
		aa += daa;
		halfsina = - halfsina;
	}
	return result;
    }
    
    private static double improve_evaluate (double x,double y[], int ixmax, int depth, boolean isMaximum) {
	
	double newy = interpolate_sinc (y, ixmax, x, depth);
	return isMaximum ? - newy : newy;
    }
    
    private static Pair<Double,Double> minimize_brent (double a, double b, double[] y, int depth, int nx, boolean isMaximum, double tol) {
	double x, v, fv, w, fw;
	double golden = 1 - NUM_goldenSection;
	double sqrt_epsilon = Math.sqrt(calculateMachineEpsilonDouble());
	long itermax = 60;

	if(!(tol > 0 && a < b))
            System.err.println("WTF?!?!?!");

	/* First step - golden section */
	v = a + golden * (b - a);
	fv = improve_evaluate(v, y, nx, depth, isMaximum);
	x = v;
        w = v;
	double fx = fv;
        fw = fv;

	for (long iter = 1; iter <= itermax; iter++) {
            double range = b - a;
            double middle_range = (a + b) / 2;
            double tol_act = sqrt_epsilon * Math.abs (x) + tol / 3;
            double new_step; /* Step at this iteration */
            
            if (Math.abs (x - middle_range) + range / 2 <= 2 * tol_act) {
                return new Pair(x,fx);
            }
            
            /* Obtain the golden section step */
            new_step = golden * (x < middle_range ? b - x : a - x);

            /* Decide if the parabolic interpolation can be tried	*/
            if (Math.abs (x - w) >= tol_act) {
                /*Interpolation step is calculated as p/q;
                 * division operation is delayed until last moment.
                */
                double p, q, t;
                t = (x - w) * (fx - fv);
                q = (x - v) * (fx - fw);
                p = (x - v) * q - (x - w) * t;
                q = 2 * (q - t);
                if (q > 0) {
                    p = -p;
                } else {
                    q = -q;
                }
                
                /*
                 * If x+p/q falls in [a,b], not too close to a and b,
                 * and isn't too large, it is accepted.
                 * If p/q is too large then the golden section procedure can
                 * reduce [a,b] range.
                 */
                if (Math.abs (p) < Math.abs (new_step * q) &&
                        p > q * (a - x + 2 * tol_act) &&
                        p < q * (b - x - 2 * tol_act)) {
                    new_step = p / q;
                }
            }
            /* Adjust the step to be not less than tolerance. */
            if (Math.abs (new_step) < tol_act) {
                new_step = new_step > 0 ? tol_act : - tol_act;
            }
            /* Obtain the next approximation to min	and reduce the enveloping range */
            {
                double t = x + new_step;	/* Tentative point for the min	*/
                double ft = improve_evaluate(t, y, nx, depth, isMaximum);
                /*
                 * If t is a better approximation, reduce the range so that
                 * t would fall within it. If x remains the best, reduce the range
                 * so that x falls within it.
                 */
                if (ft <= fx) {
                    if (t < x) {
                        b = x;
                    } else {
                        a = x;
                    }
                    v = w;  w = x;  x = t;
                    fv = fw;  fw = fx;  fx = ft;
                } else {
                    if (t < x) {
                        a = t;
                    } else {
                        b = t;
                    }
                    if (ft <= fw || w == x) {
                        v = w; w = t;
                        fv = fw; fw = ft;
                    } else if (ft <= fv || v == x || v == w) {
                        v = t;
                        fv = ft;
                    }
                }
            }
        }
        return new Pair(x,fx);
    }
    
  public static Pair<Double,Double> improveMinimum (double[] y, int nx, int ixmid, int interpolation){ 
      return improveExtremum (y, nx, ixmid, interpolation, false);
  }
  
  public static Pair<Double,Double> improveMaximum (double[] y, int nx, int ixmid, int interpolation){
      return improveExtremum (y, nx, ixmid, interpolation , true);
  }
        
  public static Pair<Double, Double> improveExtremum (double[] y, int nx, int ixmid, int interpolation, boolean isMaximum) {
	//struct improve_params params;
	double result;
        double ixmid_real;
	if (ixmid <=0) {return new Pair(y[0],0.0);}
	if (ixmid >=(nx-1)) { return new Pair(y[nx-1],(nx-1)*1.0); }
	if (interpolation <= NUM_PEAK_INTERPOLATE_NONE){ return new Pair(y[ixmid],1.0*ixmid);}
	if (interpolation == NUM_PEAK_INTERPOLATE_PARABOLIC) {
		double dy = 0.5 * (y [ixmid + 1] - y [ixmid - 1]);
		double d2y = 2 * y [ixmid] - y [ixmid - 1] - y [ixmid + 1];
		ixmid_real = ixmid + dy / d2y;
		return new Pair(y [ixmid] + 0.5 * dy * dy / d2y, ixmid_real);
	}
	/* Sinc interpolation. */
	
	int depth = interpolation == NUM_PEAK_INTERPOLATE_SINC70 ? 70 : 700;
	int ixmax = nx;
        Pair<Double,Double> pair = minimize_brent (ixmid - 1, ixmid + 1,y, depth, ixmax, isMaximum, 1e-10);
        ixmid_real = pair.first;
        result = pair.second;
	return new Pair(isMaximum ? - result : result, ixmid_real);
}
  
  public static double Vector_getValueAtX (WavData sound, double x, int ilevel, int interpolation) {
      
        double tN = sound.t0 + sound.getFrameSize()*(sound.getNumSamples()-1);
	double leftEdge = sound.t0 - 0.5 * sound.getFrameSize(), rightEdge = leftEdge + sound.sampleSize * sound.getFrameSize();
	if (x <  leftEdge || x > rightEdge) return NUMundefined;
	if (ilevel >= Vector_CHANNEL_AVERAGE) {
            if(ilevel >= sound.numberOfChannels)
                System.err.println("WTF>!>!>");
            return interpolate_sinc (sound.getSamples(ilevel), sound.getNumSamples(), SampledUtils.xToIndex (sound.t0,sound.getFrameSize(),x),
			interpolation == Vector_VALUE_INTERPOLATION_SINC70 ? NUM_VALUE_INTERPOLATE_SINC70 :
			interpolation == Vector_VALUE_INTERPOLATION_SINC700 ? NUM_VALUE_INTERPOLATE_SINC700 :
			interpolation);
	}
	double sum = 0.0;
	for (int channel = 0; channel < sound.numberOfChannels; channel ++) {
            sum += interpolate_sinc (sound.getSamples(channel), sound.getNumSamples(), SampledUtils.xToIndex (sound.t0,sound.getFrameSize(),x),
			interpolation == Vector_VALUE_INTERPOLATION_SINC70 ? NUM_VALUE_INTERPOLATE_SINC70 :
			interpolation == Vector_VALUE_INTERPOLATION_SINC700 ? NUM_VALUE_INTERPOLATE_SINC700 :
			interpolation);
	}
	return sum / sound.numberOfChannels;
  }
 
  //checked
  public static double hertzToMel (double hertz) {
      return hertz < 0 ? NUMundefined : 550.0 * Math.log (1.0 + hertz / 550.0);
  }
  
  //checked
  public static double hertzToErb (double hertz) {
     return hertz < 0 ? NUMundefined : 11.17 *  Math.log ((hertz + 312.0) / (hertz + 14680.0)) + 43.0;
  }
  
  
  //chcked
  public static boolean defined(double num){
      return (!Double.isNaN(num));
  }
  
  //checked
  public static Triple<Boolean,Double, Double> Function_intersectRangeWithDomain (Contour contour,double x1, double x2) {
	if (x1 == x2) return new Triple(false, x1, x2);
	if (x1 < x2) {
		if (x1 < contour.getStart()) x1 = contour.getStart();   // intersect requested range with logical domain
		if (x2 > contour.getEnd()) x2 = contour.getEnd();
		if (x2 <= x1) return new Triple(false,x1,x2);   // requested range and logical domain do not intersect
	} else {
		if (x2 < contour.getStart()) x2 = contour.getStart();   // intersect requested range with logical domain
		if (x1 > contour.getEnd()) x1 = contour.getEnd();
		if (x1 <= x2) return new Triple(false,x1,x2);   // requested range and logical domain do not intersect
	}
	return new Triple(true,x1,x2);
}


    
}
