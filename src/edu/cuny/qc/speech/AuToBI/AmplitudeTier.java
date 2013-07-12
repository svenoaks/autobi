/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.cuny.qc.speech.AuToBI;

import edu.cuny.qc.speech.AuToBI.core.Pair;
import edu.cuny.qc.speech.AuToBI.core.Triple;
import edu.cuny.qc.speech.AuToBI.core.WavData;
import edu.cuny.qc.speech.AuToBI.util.NUMUtils;
import java.util.ArrayList;

/**
 *
 * @author vitisoto
 */
public class AmplitudeTier {
    double xmin;
    double xmax;
    
    ArrayList<Pair<Double, Double>> points;
    
    public AmplitudeTier(double tmin, double tmax){
        xmin = tmin;
        xmax = tmax;
    }
    
    public AmplitudeTier (PointProcess me, WavData thee, double tmin, double tmax, double pmin, double pmax, double maximumPeriodFactor){
        int imin, imax, numberOfPeaks;
        if (tmax <= tmin){tmin = me.xmin; tmax = me.xmax;}
        
        Triple<Integer,Integer,Integer> tri = me.getWindowPoints (tmin, tmax);
        numberOfPeaks = tri.first;
        imin = tri.second;
        imax = tri.third;

        //if (numberOfPeaks < 3) Melder_throw ("Too few pulses between ", tmin, " and ", tmax, " seconds.");
        this.xmin = tmin;
        this.xmax = tmax;
        this.points  = new ArrayList<Pair<Double,Double>>();

        for (int i = imin + 1; i < imax; i ++) {
            double p1 = me.getPulse(i) - me.getPulse(i - 1);
            double p2 = me.getPulse(i + 1) - me.getPulse(i);
            double intervalFactor = p1 > p2 ? p1 / p2 : p2 / p1;
            if (pmin == pmax || (p1 >= pmin && p1 <= pmax && p2 >= pmin && p2 <= pmax && intervalFactor <= maximumPeriodFactor)) {
                double peak = thee.getHannWindowedRms (me.getPulse(i), 0.2 * p1, 0.2 * p2);
                if (NUMUtils.defined (peak) && peak > 0.0)
                    points.add(new Pair(me.getPulse(i), peak));
            }
        }
    }
    
    public double getShimmerLocal (double pmin, double pmax, double maximumAmplitudeFactor) {
	long numberOfPeaks = 0;
	double numerator = 0.0, denominator = 0.0;
        //Pair<Double>
	//RealPoint *points = (RealPoint *) my points -> item;
	for (int i = 1; i < points.size(); i ++) {
            double p = points.get(i).first - points.get(i - 1).first;
            if (pmin == pmax || (p >= pmin && p <= pmax)) {
                double a1 = points.get(i - 1).second, a2 = points.get(i).second;
                double amplitudeFactor = a1 > a2 ? a1 / a2 : a2 / a1;
                if (amplitudeFactor <= maximumAmplitudeFactor) {
                    numerator += Math.abs (a1 - a2);
                    numberOfPeaks ++;
                }
            }
        }
	if (numberOfPeaks < 1) return Pitch.NUMundefined;
	numerator /= numberOfPeaks;
	numberOfPeaks = 0;
	for (int i = 0; i < points.size() - 1; i ++) {
		denominator += points.get(i).second;
		numberOfPeaks ++;
	}
	denominator /= numberOfPeaks;
	if (denominator == 0.0) return Pitch.NUMundefined;
	return numerator / denominator;
    }
    
    public double getShimmerLocalDB (double pmin, double pmax, double maximumAmplitudeFactor) {
	long numberOfPeaks = 0;
	double result = 0.0;
	
	for (int i = 1; i < points.size(); i ++) {
            double p = points.get(i).first - points.get(i - 1).first;
            if (pmin == pmax || (p >= pmin && p <= pmax)) {
                double a1 = points.get(i - 1).second, a2 = points.get(i).second;
                double amplitudeFactor = a1 > a2 ? a1 / a2 : a2 / a1;
                if (amplitudeFactor <= maximumAmplitudeFactor) {
                    result += Math.abs (Math.log10 (a1 / a2));
                    numberOfPeaks ++;
                }
            }
        }
	if (numberOfPeaks < 1) return Pitch.NUMundefined;
	result /= numberOfPeaks;
	return 20.0 * result;
    }
    
    public double getShimmerAPQ3(double pmin, double pmax, double maximumAmplitudeFactor) {
	long numberOfPeaks = 0;
	double numerator = 0.0, denominator = 0.0;
	
	for (int i = 1; i <= points.size() - 2; i ++) {
            double p1 = points.get(i).first - points.get(i - 1).first;
            double p2 = points.get(i + 1).first - points.get(i).first;
            if (pmin == pmax || (p1 >= pmin && p1 <= pmax && p2 >= pmin && p2 <= pmax)) {
                double a1 = points.get(i - 1).second, a2 = points.get(i).second, a3 = points.get(i + 1).second;
                double f1 = a1 > a2 ? a1 / a2 : a2 / a1, f2 = a2 > a3 ? a2 / a3 : a3 / a2;
                if (f1 <= maximumAmplitudeFactor && f2 <= maximumAmplitudeFactor) {
                    double threePointAverage = (a1 + a2 + a3) / 3.0;
                    numerator += Math.abs (a2 - threePointAverage);
                    numberOfPeaks ++;
                }
            }
        }
	if (numberOfPeaks < 1) return Pitch.NUMundefined;
	numerator /= numberOfPeaks;
	numberOfPeaks = 0;
	for (int i = 0; i < points.size()-1; i ++) {
            denominator += points.get(i).second;
            numberOfPeaks ++;
	}
	denominator /= numberOfPeaks;
	if (denominator == 0.0) return Pitch.NUMundefined;
	return numerator / denominator;
    }
    
    
    public double getShimmerAPQ5 (double pmin, double pmax, double maximumAmplitudeFactor) {
	long numberOfPeaks = 0;
	double numerator = 0.0, denominator = 0.0;
	
	for (int i = 2; i <= points.size()-3; i ++) {
            double p1 = points.get(i - 1).first - points.get(i - 2).first;
            double p2 = points.get(i).first - points.get(i - 1).first;
            double p3 = points.get(i + 1).first - points.get(i).first;
            double p4 = points.get(i + 2).first - points.get(i + 1).first;
            if (pmin == pmax || (p1 >= pmin && p1 <= pmax && p2 >= pmin && p2 <= pmax && p3 >= pmin && p3 <= pmax && p4 >= pmin && p4 <= pmax))
            {
                double a1 = points.get(i - 2).second, a2 = points.get(i - 1).second, a3 = points.get(i).second,
				a4 = points.get(i + 1).second, a5 = points.get(i + 2).second;
			double f1 = a1 > a2 ? a1 / a2 : a2 / a1, f2 = a2 > a3 ? a2 / a3 : a3 / a2,
				f3 = a3 > a4 ? a3 / a4 : a4 / a3, f4 = a4 > a5 ? a4 / a5 : a5 / a4;
			if (f1 <= maximumAmplitudeFactor && f2 <= maximumAmplitudeFactor &&
			    f3 <= maximumAmplitudeFactor && f4 <= maximumAmplitudeFactor)
			{
				double fivePointAverage = (a1 + a2 + a3 + a4 + a5) / 5.0;
				numerator += Math.abs (a3 - fivePointAverage);
				numberOfPeaks ++;
			}
		}
	}
	if (numberOfPeaks < 1) return Pitch.NUMundefined;
	numerator /= numberOfPeaks;
	numberOfPeaks = 0;
	for (int i = 0; i < points.size() - 1; i ++) {
            denominator += points.get(i).second;
            numberOfPeaks ++;
	}
	denominator /= numberOfPeaks;
	if (denominator == 0.0) return Pitch.NUMundefined;
	return numerator / denominator;
    }
    
    public double getShimmerAPQ11 (double pmin, double pmax, double maximumAmplitudeFactor) {
	long numberOfPeaks = 0;
	double numerator = 0.0, denominator = 0.0;
	
	for (int i = 5; i <= points.size() - 6; i ++) {
            double
                    p1 = points.get(i - 4).first - points.get(i - 5).first,
                    p2 = points.get(i - 3).first - points.get(i - 4).first,
                    p3 = points.get(i - 2).first - points.get(i - 3).first,
                    p4 = points.get(i - 1).first - points.get(i - 2).first,
                    p5 = points.get(i).first - points.get(i - 1).first,
                    p6 = points.get(i + 1).first - points.get(i).first,
                    p7 = points.get(i + 2).first - points.get(i + 1).first,
                    p8 = points.get(i + 3).first - points.get(i + 2).first,
                    p9 = points.get(i + 4).first - points.get(i + 3).first,
                    p10 = points.get(i + 5).first - points.get(i + 4).first;
            if (pmin == pmax || (p1 >= pmin && p1 <= pmax && p2 >= pmin && p2 <= pmax
                    && p3 >= pmin && p3 <= pmax && p4 >= pmin && p4 <= pmax && p5 >= pmin && p5 <= pmax
                    && p6 >= pmin && p6 <= pmax && p7 >= pmin && p7 <= pmax && p8 >= pmin && p8 <= pmax
                    && p9 >= pmin && p9 <= pmax && p10 >= pmin && p10 <= pmax))
            {
                double a1 = points.get(i - 5).second, a2 = points.get(i - 4).second, a3 = points.get(i - 3).second,
                        a4 = points.get(i - 2).second, a5 = points.get(i - 1).second, a6 = points.get(i).second,
                        a7 = points.get(i + 1).second, a8 = points.get(i + 2).second, a9 = points.get(i + 3).second,
                        a10 = points.get(i + 4).second, a11 = points.get(i + 5).second;
                double f1 = a1 > a2 ? a1 / a2 : a2 / a1, f2 = a2 > a3 ? a2 / a3 : a3 / a2,
                        f3 = a3 > a4 ? a3 / a4 : a4 / a3, f4 = a4 > a5 ? a4 / a5 : a5 / a4,
                        f5 = a5 > a6 ? a5 / a6 : a6 / a5, f6 = a6 > a7 ? a6 / a7 : a7 / a6,
                        f7 = a7 > a8 ? a7 / a8 : a8 / a7, f8 = a8 > a9 ? a8 / a9 : a9 / a8,
                        f9 = a9 > a10 ? a9 / a10 : a10 / a9, f10 = a10 > a11 ? a10 / a11 : a11 / a10;
                if (f1 <= maximumAmplitudeFactor && f2 <= maximumAmplitudeFactor &&
                        f3 <= maximumAmplitudeFactor && f4 <= maximumAmplitudeFactor &&
                        f5 <= maximumAmplitudeFactor && f6 <= maximumAmplitudeFactor &&
                        f7 <= maximumAmplitudeFactor && f8 <= maximumAmplitudeFactor &&
                        f9 <= maximumAmplitudeFactor && f10 <= maximumAmplitudeFactor)
                {
                    double elevenPointAverage = (a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9 + a10 + a11) / 11.0;
                    numerator += Math.abs (a6 - elevenPointAverage);
                    numberOfPeaks ++;
                }
            }
        }
	if (numberOfPeaks < 1) return Pitch.NUMundefined;
	numerator /= numberOfPeaks;
	numberOfPeaks = 0;
	for (int i = 0; i < points.size()-1; i ++) {
            denominator += points.get(i).second;
            numberOfPeaks ++;
	}
	denominator /= numberOfPeaks;
	if (denominator == 0.0) return Pitch.NUMundefined;
	return numerator / denominator;
    }

    
}
