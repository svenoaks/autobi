/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.cuny.qc.speech.AuToBI;


import edu.cuny.qc.speech.AuToBI.core.Triple;
import edu.cuny.qc.speech.AuToBI.core.WavData;
import edu.cuny.qc.speech.AuToBI.util.NUMUtils;
import edu.cuny.qc.speech.AuToBI.util.SampledUtils;
import edu.cuny.qc.speech.AuToBI.util.SoundUtils;

/**
 *
 * @author vitisoto
 */
public class PointProcess {
    
    public double xmin;
    public double xmax;
    private int maxnt;
    private int nt;
    private double[] t;
    
    public PointProcess(double tmin, double tmax, int initialMaxnt){
        if (initialMaxnt < 1) initialMaxnt = 1;
        
        this.xmin = tmin;
        this.xmax = tmax;
	this.maxnt = initialMaxnt;
	this.nt = 0;
	this.t = new double[maxnt];
    }
    
    public PointProcess(WavData sound, Pitch pitch) {
        this(sound.t0, sound.t0+sound.getDuration(), 10);
        
            
        double t = pitch.getStart();
        double addedRight = -1e300;
        double globalPeak = SoundUtils.getAbsoluteExtremum(sound, sound.t0, sound.t0+sound.getDuration(), 0);
        double peak;
        
        int numAdd = 0;
        /*
         * Cycle over all voiced intervals.
         */
        //autoMelderProgress progress (L"Sound & Pitch: To PointProcess...");
        for (;;) {
            
            Triple<Boolean, Double, Double> triVoiced = pitch.getVoicedIntervalAfter (t);
            if (!triVoiced.first) break;
            double tleft = triVoiced.second;
            double tright =triVoiced.third;

            //Go to the middle of the voice stretch.
            double tmiddle = (tleft + tright) / 2;
            //Melder_progress ((tmiddle - sound -> xmin) / (sound -> xmax - sound -> xmin), L"Sound & Pitch to PointProcess");
            double f0middle = pitch.getValueAtTime (tmiddle,Pitch.kPitch_unit_HERTZ, Pitch.Pitch_LINEAR);
            // Our first point is near this middle.

            if (Double.isNaN(f0middle)) {
                System.err.println("Sound_Pitch_to_PointProcess_cc: tleft "+tleft+", tright "+tright+", f0middle "+f0middle);
            }

            double tmax = findExtremum (sound, tmiddle - 0.5 / f0middle, tmiddle + 0.5 / f0middle, true, true);
            
            if(!NUMUtils.defined(tmax))
                System.err.println("tmax undefined");
            
            addPoint (tmax);
            numAdd++;
            double tsave = tmax;
            
            for (;;) {
                double f0 = pitch.getValueAtTime (tmax, Pitch.kPitch_unit_HERTZ, Pitch.Pitch_LINEAR);
                double correlation;
                
                if (Double.isNaN(f0)) break;
                
                Triple<Double, Double, Double> triCorr = findMaximumCorrelation (sound, tmax, 1.0 / f0, tmax - 1.25 / f0, tmax - 0.8 / f0);

                correlation = triCorr.first;
                tmax = triCorr.second;
                peak = triCorr.third;

                if (correlation == -1) /*break*/ tmax -= 1.0 / f0;   /* This one period will drop out. */
                if (tmax < tleft) {
                    if (correlation > 0.7 && peak > 0.023333 * globalPeak && tmax - addedRight > 0.8 / f0) {
                        addPoint (tmax);
                        numAdd++;
                    }
                    break;
                }

                if (correlation > 0.3 && (peak == 0.0 || peak > 0.01 * globalPeak)) {
                    if (tmax - addedRight > 0.8 / f0) {   // do not fill in a short originally unvoiced interval twice
                        addPoint (tmax);
                        numAdd++;
                    }
                }
            }

            tmax = tsave;
            
            
            for (;;) {
                double f0 = pitch.getValueAtTime (tmax, Pitch.kPitch_unit_HERTZ, Pitch.Pitch_LINEAR);
                double correlation;
                if (Double.isNaN(f0)) break;
                Triple<Double, Double, Double> triCorr = findMaximumCorrelation (sound, tmax, 1.0 / f0, tmax + 0.8 / f0, tmax + 1.25 / f0);
                correlation = triCorr.first;
                tmax = triCorr.second;
                peak = triCorr.third;

                if (correlation == -1) /*break*/ tmax += 1.0 / f0;
                if (tmax > tright) {
                    if (correlation > 0.7 && peak > 0.023333 * globalPeak) {
                        addPoint (tmax);
                        numAdd++;
                        addedRight = tmax;
                    }

                    break;
                }

                if (correlation > 0.3 && (peak == 0.0 || peak > 0.01 * globalPeak)) {
                    addPoint (tmax);
                    numAdd++;
                    addedRight = tmax;
                }
            }

            t = tright;
        }   
       
    }
    
    public int size(){
        return this.nt;
    }
    
    public double get(int i){
        return t[i];
    }
    
    private static double findExtremum_3 (double[] channel1_base, double[] channel2_base, int d, int n, boolean includeMaxima, boolean includeMinima) {
	
        
	boolean includeAll = includeMaxima == includeMinima;
	int imin = 0, imax = 0, i, iextr;
	double minimum, maximum;
	if (n < 3) {
		if (n <= 0) return -1.0;   /* Outside. */
		else if (n == 1) return 0.0;
		else {   /* n == 2 */
			double x1 = channel2_base!=null ? 0.5 * (channel1_base [d] + channel2_base [d]) : channel1_base [d];
			double x2 = channel2_base!=null ? 0.5 * (channel1_base [d+1] + channel2_base [d+1]) : channel1_base [d+1];
			double xleft = includeAll ? Math.abs (x1) : includeMaxima ? x1 : - x1;
			double xright = includeAll ? Math.abs (x2) : includeMaxima ? x2 : - x2;
			if (xleft > xright) return 0.0;
			else if (xleft < xright) return 1.0;
			else return 0.5;
		}
	}
	minimum = maximum = channel2_base!=null ? 0.5 * (channel1_base [d] + channel2_base [d]) : channel1_base[d];
	for (i = 1; i < n; i ++) {
		double value = channel2_base!=null ? 0.5 * (channel1_base[d+i] + channel2_base [d+i]) : channel1_base[d+i];
		if (value < minimum) { minimum = value; imin = i; }
		if (value > maximum) { maximum = value; imax = i; }
	}
	if (minimum == maximum) {
		return 0.5 * (n + 1.0)-1;   /* All equal. */
	}
	iextr = includeAll ? ( Math.abs (minimum) > Math.abs (maximum) ? imin : imax ) : includeMaxima ? imax : imin;
	if (iextr == 0) return 0.0;
	if (iextr == n-1) return (double) n-1;
	/* Parabolic interpolation. */
	/* We do NOT need fabs here: we look for a genuine extremum. */
	double valueMid = channel2_base!=null ? 0.5 * (channel1_base [d+iextr] + channel2_base [d+iextr]) : channel1_base [d+iextr];
	double valueLeft = channel2_base!=null ? 0.5 * (channel1_base [d+iextr - 1] + channel2_base [d+iextr - 1]) : channel1_base [d+iextr - 1];
	double valueRight = channel2_base!=null ? 0.5 * (channel1_base [d+iextr + 1] + channel2_base [d+iextr + 1]) : channel1_base [d+iextr + 1];
	return iextr + 0.5 * (valueRight - valueLeft) / (2 * valueMid - valueLeft - valueRight);
    }

    private static double findExtremum (WavData me, double tmin, double tmax, boolean includeMaxima, boolean includeMinima) {
	int imin = SampledUtils.xToLowIndex (me.t0,me.getFrameSize(), tmin), imax = SampledUtils.xToHighIndex (me.t0,me.getFrameSize(),tmax);
	//Melder_assert (NUMdefined (tmin));
	//Melder_assert (NUMdefined (tmax));
	if (imin < 0) imin = 0;
	if (imax >=  me.getNumSamples()) imax = me.getNumSamples()-1;
	double iextremum = findExtremum_3 (me.getSamples(0), me.numberOfChannels > 1 ? me.getSamples(1) : null, imin, imax - imin + 1, includeMaxima, includeMinima);
	if (iextremum >= 0.0)
		return me.t0 + (imin + iextremum) * me.getFrameSize();
	else
		return (tmin + tmax) / 2;
    }

    private void addPoint(double t) {
    
        if (Double.isNaN(t))
            return;
        
        if (this.nt >= this.maxnt) {
            /*
             * Create without change.
             */
            double[] dum = new double[2 * this.maxnt];
            
            for(int i=0;i<this.nt;i++)
                dum[i]=this.t[i];
            
            
            this.t = dum;
            this.maxnt *= 2;
        }
        
        if (this.nt == 0 || t >= this.t [this.nt-1]) {   // special case that often occurs in practice
            this.t [this.nt] = t;
            this.nt++;
        } else {
            int left = getLowIndex (t);
            if (left == -1 || this.t [left] != t) {
                for (int i = this.nt-1; i > left; i --) this.t [i + 1] = this.t [i];
                this.nt ++;
                this.t [left+1] = t;
            }
        }
    }
    
    private int getLowIndex (double t) {
	if (this.nt == 0 || t < this.t[0])
		return -1;
	if (t >= this.t [this.nt-1])   /* Special case that often occurs in practice. */
		return this.nt-1;
	//Melder_assert (this.nt != 1);   /* May fail if t or my t [1] is NaN. */
	/* Start binary search. */
	int left = 0, right = this.nt-1;
	while (left < right - 1) {
		int mid = (left + right) / 2;
		if (t >= this.t [mid]) left = mid; else right = mid;
	}
	//Melder_assert (right == left + 1);
	return left;
    }
    
    private Triple<Double,Double,Double> findMaximumCorrelation (WavData me, double t1, double windowLength, double tmin2, double tmax2) {
	double maximumCorrelation = -1.0, r1 = 0.0, r2 = 0.0, r3 = 0.0, r1_best=0.0, r3_best=0.0, ir=0.0;
	double halfWindowLength = 0.5 * windowLength;
	int ileft1 = SampledUtils.xToNearestIndex (me.t0,me.getFrameSize(),t1 - halfWindowLength);
	int iright1 = SampledUtils.xToNearestIndex (me.t0,me.getFrameSize(), t1 + halfWindowLength);
	int ileft2min = SampledUtils.xToLowIndex (me.t0,me.getFrameSize(), tmin2 - halfWindowLength);
	int ileft2max = SampledUtils.xToHighIndex (me.t0,me.getFrameSize(), tmax2 - halfWindowLength);
        
	double peak = 0.0;   /* Default. */
	for (int ileft2 = ileft2min; ileft2 <= ileft2max; ileft2 ++) {
		double norm1 = 0.0, norm2 = 0.0, product = 0.0, localPeak = 0.0;
		if (me.numberOfChannels == 1) {
			for (int i1 = ileft1, i2 = ileft2; i1 <= iright1; i1 ++, i2 ++) {
				if (i1 < 0 || i1 >= me.getNumSamples() || i2 < 0 || i2 >= me.getNumSamples()) continue;
				double amp1 = me.getSample(0,i1);
                                double amp2 = me.getSample(0,i2);
				norm1 += amp1 * amp1;
				norm2 += amp2 * amp2;
				product += amp1 * amp2;
				if (Math.abs (amp2) > localPeak)
					localPeak = Math.abs (amp2);
			}
		} else {
			for (int i1 = ileft1, i2 = ileft2; i1 <= iright1; i1 ++, i2 ++) {
				if (i1 < 0 || i1 >= me.getNumSamples() || i2 < 0 || i2 >= me.getNumSamples()) continue;
				double amp1 = 0.5 * (me.getSample(0,i1) + me.getSample(1,i1));
                                double amp2 = 0.5 * (me.getSample(0,i2) + me.getSample(1,i2));
				norm1 += amp1 * amp1;
				norm2 += amp2 * amp2;
				product += amp1 * amp2;
				if (Math.abs (amp2) > localPeak)
					localPeak = Math.abs(amp2);
			}
		}
		r1 = r2;
		r2 = r3;
		r3 = product!=0 ? product / (Math.sqrt (norm1 * norm2)) : 0.0;
		if (r2 > maximumCorrelation && r2 >= r1 && r2 >= r3) {
			r1_best = r1;
			maximumCorrelation = r2;
			r3_best = r3;
			ir = ileft2 - 1;
			peak = localPeak;  
		}
	}
	/*
	 * Improve the result by means of parabolic interpolation.
	 */
        double tout=0.0;
	if (maximumCorrelation > -1.0) {
		double d2r = 2 * maximumCorrelation - r1_best - r3_best;
		if (d2r != 0.0) {
			double dr = 0.5 * (r3_best - r1_best);
			maximumCorrelation += 0.5 * dr * dr / d2r;
			ir += dr / d2r;
		}
		tout = t1 + (ir - ileft1) * me.getFrameSize();
	}
        return new Triple(maximumCorrelation,tout,peak);
    }
    
    public Triple<Integer,Integer,Integer> getWindowPoints (double tmin, double tmax) {
	int imin = getHighIndex (tmin);
	int imax = getLowIndex (tmax);
        return new Triple(imax-imin+1,imin,imax);
	
    }
    
    public int getNumberOfPeriods (double tmin, double tmax,
            double minimumPeriod, double maximumPeriod, double maximumPeriodFactor)
    {
        int imin, imax, numberOfPeriods, i;
	if (tmax <= tmin){tmin = this.xmin; tmax =this.xmax;}   /* Autowindowing. */
	Triple<Integer,Integer,Integer> tri = getWindowPoints (tmin, tmax);
        numberOfPeriods = tri.first-1;
        imin = tri.second;
        imax = tri.third;
        
	if (numberOfPeriods < 1) return 0;
	for (i = imin; i < imax; i ++) {
		if (isPeriod (i, minimumPeriod, maximumPeriod, maximumPeriodFactor)) {
			//(void) 0;   /* This interval counts as a period. */
		} else {
			numberOfPeriods --;   /* This interval does not count as a period. */
		}
	}
	return numberOfPeriods;
    }
    
    
    private int getHighIndex (double t) {
	if (this.nt == 0)
		return -1;
	if (t <= this.t [0])
		return 0;
	if (t > this.t [this.nt-1])
		return this.nt;
	/* Start binary search. */
	int left = 0, right = this.nt-1;
	while (left < right - 1) {
		int mid = (left + right) / 2;
		if (t > this.t [mid]) left = mid; else right = mid;
	}
	//Melder_assert (right == left + 1);
	return right;
    }
    
    private boolean isPeriod (int ileft, double minimumPeriod, double maximumPeriod, double maximumPeriodFactor) {
	/*
	 * This function answers the question: is the interval from point 'ileft' to point 'ileft+1' a period?
	 */
	int iright = ileft + 1;
	/*
	 * Period condition 1: both 'ileft' and 'iright' have to be within the point process.
	 */
	if (ileft < 0 || iright > this.nt-1) {
		return false;
	} else {
		/*
		 * Period condition 2: the interval has to be within the boundaries, if specified.
		 */
		if (minimumPeriod == maximumPeriod) {
			return true;   /* All intervals count as periods, irrespective of absolute size and relative size. */
		} else {
			double interval = this.t [iright] - this.t [ileft];
			if (interval <= 0.0 || interval < minimumPeriod || interval > maximumPeriod) {
				return false;
			} else if (!NUMUtils.defined (maximumPeriodFactor) || maximumPeriodFactor < 1.0) {
				return true;
			} else {
				/*
				 * Period condition 3: the interval cannot be too different from both of its neigbours, if any.
				 */
				double previousInterval = ileft <= 0 ? Pitch.NUMundefined : this.t [ileft] - this.t [ileft - 1];
				double nextInterval = iright >= this.nt-1 ? Pitch.NUMundefined : this.t [iright + 1] - this.t [iright];
				double previousIntervalFactor = NUMUtils.defined (previousInterval) && previousInterval > 0.0 ? interval / previousInterval : Pitch.NUMundefined;
				double nextIntervalFactor = NUMUtils.defined (nextInterval) && nextInterval > 0.0 ? interval / nextInterval : Pitch.NUMundefined;
				if (! NUMUtils.defined (previousIntervalFactor) && ! NUMUtils.defined (nextIntervalFactor)) {
					return true;   /* No neighbours: this is a period. */
				}
				if (NUMUtils.defined (previousIntervalFactor) && previousIntervalFactor > 0.0 && previousIntervalFactor < 1.0) {
					previousIntervalFactor = 1.0 / previousIntervalFactor;
				}
				if (NUMUtils.defined (nextIntervalFactor) && nextIntervalFactor > 0.0 && nextIntervalFactor < 1.0) {
					nextIntervalFactor = 1.0 / nextIntervalFactor;
				}
				if (NUMUtils.defined (previousIntervalFactor) && previousIntervalFactor > maximumPeriodFactor &&
					NUMUtils.defined (nextIntervalFactor) && nextIntervalFactor > maximumPeriodFactor)
				{
					return false;
				}
			}
		}
	}
	return true;
    }
    
    
    public double getMeanPeriod (double tmin, double tmax, double minimumPeriod, double maximumPeriod, double maximumPeriodFactor){
        
	if (tmax <= tmin){tmin = this.xmin; tmax = this.xmax;}   /* Autowindowing. */
	
	Triple<Integer,Integer,Integer> tri = getWindowPoints (tmin, tmax);
        int numberOfPeriods = tri.first - 1;
        int imin = tri.second;
        int imax = tri.third;
        
	if (numberOfPeriods < 1) return Pitch.NUMundefined;
	double sum = 0.0;
	for (int i = imin; i < imax; i ++) {
		if (isPeriod (i, minimumPeriod, maximumPeriod, maximumPeriodFactor)) {
			sum += this.t [i + 1] - this.t [i];   /* This interval counts as a period. */
		} else {
			numberOfPeriods --;   /* This interval does not count as a period. */
		}
	}
	return numberOfPeriods > 0 ? sum / numberOfPeriods : Pitch.NUMundefined;
    }
    
    public double getStdevPeriod (double tmin, double tmax, double minimumPeriod, double maximumPeriod, double maximumPeriodFactor){
        if (tmax <= tmin){tmin = this.xmin; tmax = this.xmax;}   /* Autowindowing. */
	
        Triple<Integer,Integer,Integer> tri = getWindowPoints (tmin, tmax);
        int numberOfPeriods = tri.first - 1;
        int imin = tri.second;
        int imax = tri.third;
        
	if (numberOfPeriods < 2) return Pitch.NUMundefined;
	
	double sum = 0.0;
	for (int i = imin; i < imax; i ++) {
		if (isPeriod (i, minimumPeriod, maximumPeriod, maximumPeriodFactor)) {
			sum += this.t [i + 1] - this.t [i];   /* This interval counts as a period. */
		} else {
			numberOfPeriods --;   /* This interval does not count as a period. */
		}
	}
        
	if (numberOfPeriods < 2) return Pitch.NUMundefined;
        
	double mean = sum / numberOfPeriods;
	
	double sum2 = 0.0;
	for (int i = imin; i < imax; i ++) {
		if (isPeriod (i, minimumPeriod, maximumPeriod, maximumPeriodFactor)) {
			double dperiod = this.t [i + 1] - this.t [i] - mean;
			sum2 += dperiod * dperiod;
		}
	}
	
	return Math.sqrt (sum2 / (numberOfPeriods - 1));
    }
    
    public double getPulse(int i){
        return t[i];
    }
    
    public double getJitterLocal (double tmin, double tmax,double pmin, double pmax, double maximumPeriodFactor)
    {
        double jitterLocalAbsolute = getJitterLocalAbsolute(tmin,tmax,pmin,pmax,maximumPeriodFactor);
        return NUMUtils.defined(jitterLocalAbsolute)? jitterLocalAbsolute/getMeanPeriod(tmin, tmax, pmin, pmax, maximumPeriodFactor) : Pitch.NUMundefined;
    }
    
    public double getJitterLocalAbsolute (double tmin, double tmax, double pmin, double pmax, double maximumPeriodFactor)
    {
        double sum = 0.0;
        if (tmax <= tmin){ tmin = this.xmin; tmax = this.xmax;}   /* Autowindowing. */
        
        Triple<Integer,Integer,Integer> tri = getWindowPoints(tmin,tmax);
        
        int numberOfPeriods = tri.first - 1;
        int imin = tri.second;
        int imax = tri.third;
        
        if (numberOfPeriods < 2) return Pitch.NUMundefined;
        for (int i = imin + 1; i < imax; i ++) {
            double p1 = this.t [i] - this.t [i - 1], p2 = this.t [i + 1] - this.t [i];
            double intervalFactor = p1 > p2 ? p1 / p2 : p2 / p1;
            if (pmin == pmax || (p1 >= pmin && p1 <= pmax && p2 >= pmin && p2 <= pmax && intervalFactor <= maximumPeriodFactor)) {
                sum += Math.abs (p1 - p2);
            } else {
                numberOfPeriods --;
            }
        }
        if (numberOfPeriods < 2) return Pitch.NUMundefined;
        return sum / (numberOfPeriods - 1);
    }
    
    public double getJitterRap (double tmin, double tmax, double pmin, double pmax, double maximumPeriodFactor){
        
	if (tmax <= tmin){ tmin = this.xmin; tmax = this.xmax; }  /* Autowindowing. */
        
        Triple<Integer,Integer,Integer> tri = getWindowPoints(tmin,tmax);
        
        int numberOfPeriods = tri.first - 1;
        int imin = tri.second;
        int imax = tri.third;
        
	if (numberOfPeriods < 3) return Pitch.NUMundefined;
	double sum = 0.0;
	for (int i = imin + 2; i < imax; i ++) {
            double p1 = this.t [i - 1] - this.t [i - 2], p2 = this.t [i] - this.t [i - 1], p3 = this.t [i + 1] - this.t [i];
            double intervalFactor1 = p1 > p2 ? p1 / p2 : p2 / p1, intervalFactor2 = p2 > p3 ? p2 / p3 : p3 / p2;
            if (pmin == pmax || (p1 >= pmin && p1 <= pmax && p2 >= pmin && p2 <= pmax && p3 >= pmin && p3 <= pmax
		    && intervalFactor1 <= maximumPeriodFactor && intervalFactor2 <= maximumPeriodFactor))
		{
			sum += Math.abs (p2 - (p1 + p2 + p3) / 3.0);
		} else {
			numberOfPeriods --;
		}
	}
	if (numberOfPeriods < 3) return Pitch.NUMundefined;
	return sum / (numberOfPeriods - 2) / getMeanPeriod (tmin, tmax, pmin, pmax, maximumPeriodFactor);
    }
    
    public double getJitterPPQ5 (double tmin, double tmax, double pmin, double pmax, double maximumPeriodFactor)
{
	if (tmax <= tmin){tmin = this.xmin; tmax = this.xmax;}   /* Autowindowing. */
	Triple<Integer,Integer,Integer> tri = getWindowPoints(tmin,tmax);
        
        int numberOfPeriods = tri.first - 1;
        int imin = tri.second;
        int imax = tri.third;
	if (numberOfPeriods < 5) return Pitch.NUMundefined;
        
	double sum = 0.0;
	for (int i = imin + 5; i <= imax; i ++) {
            double
			p1 = this.t [i - 4] - this.t [i - 5],
			p2 = this.t [i - 3] - this.t [i - 4],
			p3 = this.t [i - 2] - this.t [i - 3],
			p4 = this.t [i - 1] - this.t [i - 2],
			p5 = this.t [i] - this.t [i - 1];
		double
			f1 = p1 > p2 ? p1 / p2 : p2 / p1,
			f2 = p2 > p3 ? p2 / p3 : p3 / p2,
			f3 = p3 > p4 ? p3 / p4 : p4 / p3,
			f4 = p4 > p5 ? p4 / p5 : p5 / p4;
		if (pmin == pmax || (p1 >= pmin && p1 <= pmax && p2 >= pmin && p2 <= pmax && p3 >= pmin && p3 <= pmax &&
			p4 >= pmin && p4 <= pmax && p5 >= pmin && p5 <= pmax &&
			f1 <= maximumPeriodFactor && f2 <= maximumPeriodFactor && f3 <= maximumPeriodFactor && f4 <= maximumPeriodFactor))
		{
			sum += Math.abs (p3 - (p1 + p2 + p3 + p4 + p5) / 5.0);
		} else {
			numberOfPeriods --;
		}
	}
	if (numberOfPeriods < 5) return Pitch.NUMundefined;
	return sum / (numberOfPeriods - 4) / getMeanPeriod (tmin, tmax, pmin, pmax, maximumPeriodFactor);
}


}
