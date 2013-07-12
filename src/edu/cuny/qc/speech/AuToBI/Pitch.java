/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.cuny.qc.speech.AuToBI;

import edu.cuny.qc.speech.AuToBI.core.Contour;
import edu.cuny.qc.speech.AuToBI.core.Pair;
import edu.cuny.qc.speech.AuToBI.core.PitchFrame;
import edu.cuny.qc.speech.AuToBI.core.Triple;
import edu.cuny.qc.speech.AuToBI.core.WavData;
import edu.cuny.qc.speech.AuToBI.util.NUMUtils;
import edu.cuny.qc.speech.AuToBI.util.SampledUtils;
import java.util.ArrayList;
import java.util.Arrays;

/**
 *
 * @author vitisoto
 */
public class Pitch {

    public final static int Pitch_LEVEL_FREQUENCY = 0;
    public final static int Pitch_LEVEL_STRENGTH = 1;
    public final static int kPitch_unit_HERTZ = 0;
    public final static int kPitch_unit_HERTZ_LOGARITHMIC = 1;
    public final static int kPitch_unit_MEL = 2;
    public final static int kPitch_unit_LOG_HERTZ = 3;
    public final static int kPitch_unit_SEMITONES_1 = 4;
    public final static int kPitch_unit_SEMITONES_100 = 5;
    public final static int kPitch_unit_SEMITONES_200 = 6;
    public final static int kPitch_unit_SEMITONES_440 = 7;
    public final static int kPitch_unit_ERB = 8;
    public final static double NUMundefined = Double.NaN;
    public final static int Pitch_STRENGTH_UNIT_min = 0;
    public final static int Pitch_STRENGTH_UNIT_AUTOCORRELATION = 0;
    public final static int Pitch_STRENGTH_UNIT_NOISE_HARMONICS_RATIO = 1;
    public final static int Pitch_STRENGTH_UNIT_HARMONICS_NOISE_DB = 2;
    public final static int Pitch_STRENGTH_UNIT_max = 2;
    public final static int Pitch_NEAREST = 0;
    public final static int Pitch_LINEAR = 1;
    public final static double NUMln2 = Math.log(2);
    private ArrayList<PitchFrame> pitchFrames;
    private Contour contour;
    private double ceiling;
    private double x1;        
    private double xmin;      
    private double xmax;      
    private double dx;        
    
    private int n;

    //Need candidates!!!!!!!!
    //checked
    public Pitch(WavData sound,ArrayList<PitchFrame> pitchFrames, Contour contour) {
        this.ceiling = contour.getCeiling();
        this.pitchFrames = pitchFrames;
        this.xmin = sound.t0;
        this.xmax = sound.t0+sound.getDuration();
        this.x1 = contour.getStart();
        this.dx = contour.getStep();
        this.n = contour.size();
        this.contour = contour;
    }

    public void normalizeZScore(double ta, double tb) {
        
        
        double meanPitchUtt = this.getMean(ta,tb,Pitch.kPitch_unit_HERTZ);
        double stdPitchUtt = this.getStandardDeviation(ta,tb,Pitch.kPitch_unit_HERTZ);
        
        Triple<Integer, Integer, Integer> tri = SampledUtils.getWindowSamples(this.x1, this.dx, this.n, ta, tb);
        int imin = tri.first;
        int imax = tri.second;
        

        if(tri.third==0)
            System.err.println("HOLY SHIT!!!");
        
        
        for (int isamp = imin; isamp <= imax; isamp++) {
            double pitch = pitchFrames.get(isamp).getCandidate(0).frequency;
            if (NUMUtils.defined(pitch)) {
                pitchFrames.get(isamp).getCandidate(0).frequency= (pitch-meanPitchUtt)/stdPitchUtt; 
            }
        }
        
    }
    
    //Checked
    public double getValueAtSample(int iframe, int ilevel, int unit) {
        if (iframe < 0 || iframe >= pitchFrames.size()) {
            return NUMundefined;
        }
        double f = pitchFrames.get(iframe).getCandidate(0).frequency;
        if (f <= 0.0 || f >= this.ceiling) {
            return NUMundefined;   // frequency out of range (or NUMundefined)? Voiceless
        }
        return convertStandardToSpecialUnit(ilevel == Pitch_LEVEL_FREQUENCY ? f : pitchFrames.get(iframe).getCandidate(0).strength, ilevel, unit);
    }

    //checked
    private static double convertStandardToSpecialUnit(double value, long ilevel, int unit) {
        if (ilevel == Pitch_LEVEL_FREQUENCY) {
            return unit == kPitch_unit_HERTZ ? value
                    : unit == kPitch_unit_HERTZ_LOGARITHMIC ? value <= 0.0 ? NUMundefined : Math.log10(value)
                    : unit == kPitch_unit_MEL ? NUMUtils.hertzToMel(value)
                    : unit == kPitch_unit_LOG_HERTZ ? value <= 0.0 ? NUMundefined : Math.log10(value)
                    : unit == kPitch_unit_SEMITONES_1 ? value <= 0.0 ? NUMundefined : 12.0 * Math.log(value / 1.0) / NUMln2
                    : unit == kPitch_unit_SEMITONES_100 ? value <= 0.0 ? NUMundefined : 12.0 * Math.log(value / 100.0) / NUMln2
                    : unit == kPitch_unit_SEMITONES_200 ? value <= 0.0 ? NUMundefined : 12.0 * Math.log(value / 200.0) / NUMln2
                    : unit == kPitch_unit_SEMITONES_440 ? value <= 0.0 ? NUMundefined : 12.0 * Math.log(value / 440.0) / NUMln2
                    : unit == kPitch_unit_ERB ? NUMUtils.hertzToErb(value)
                    : NUMundefined;
        } else {
            return unit == Pitch_STRENGTH_UNIT_AUTOCORRELATION ? value
                    : unit == Pitch_STRENGTH_UNIT_NOISE_HARMONICS_RATIO
                    ? value <= 1e-15 ? 1e15 : value > 1.0 - 1e-15 ? 1e-15 : (1.0 - value) / value
                    : /* Before losing precision. */ unit == Pitch_STRENGTH_UNIT_HARMONICS_NOISE_DB
                    ? value <= 1e-15 ? -150.0 : value > 1.0 - 1e-15 ? 150.0 : 10 * Math.log10(value / (1.0 - value))
                    : /* Before losing precision. */ NUMundefined;
        }
    }

    public boolean isVoiced_i(int iframe) {
        return NUMUtils.defined(getValueAtSample(iframe, Pitch_LEVEL_FREQUENCY, kPitch_unit_HERTZ));
    }

    public Triple<Boolean, Double, Double> getVoicedIntervalAfter(double after) {
        double tleft = 0.0, tright = 0.0;

        int ileft = SampledUtils.xToHighIndex(this.x1, this.dx, after);
        int iright;

        if (ileft >= this.n) {
            return new Triple(false, tleft, tright);   /* Offright. */
        }
        if (ileft < 0) {
            ileft = 0;   /* Offleft. */
        }

        /* Search for first voiced frame. */
        for (; ileft < this.n; ileft++) {
            if (isVoiced_i(ileft)) {
                break;
            }
        }

        if (ileft >= this.n) {
            return new Triple(false, tleft, tright);   /* Offright. */
        }

        /* Search for last voiced frame. */
        for (iright = ileft; iright < this.n; iright++) {
            if (!isVoiced_i(iright)) {
                break;
            }
        }
        iright--;

        tleft = SampledUtils.indexToX(this.x1, this.dx, ileft) - 0.5 * this.dx;   /* The whole frame is considered voiced. */
        tright = SampledUtils.indexToX(this.x1, this.dx, iright) + 0.5 * this.dx;


        if (tleft >= this.xmax - 0.5 * this.dx) {
            return new Triple(false, tleft, tright);
        }
        if (tleft < this.xmin) {
            tleft = this.xmin;
        }
        if (tright > this.xmax) {
            tright = this.xmax;
        }

        return new Triple(true, tleft, tright);
    }

    public double getValueAtTime(double time, int unit, int interpolate) {
        return getValueAtX(time, Pitch_LEVEL_FREQUENCY, unit, interpolate);
    }

    public double getValueAtX(double x, int ilevel, int unit, int interpolate) {
        if (x < this.xmin || x > this.xmax) {
            return NUMundefined;
        }
        if (interpolate > 0) {
            double ireal = SampledUtils.xToIndex(this.x1, this.dx, x);
            int ileft = (int) Math.floor(ireal);
            int inear, ifar;
            double phase = ireal - ileft;
            if (phase < 0.5) {
                inear = ileft;
                ifar = ileft + 1;
            } else {
                ifar = ileft;
                inear = ileft + 1;
                phase = 1.0 - phase;
            }

            if (inear < 0 || inear >= this.n) {
                return NUMundefined;   // x out of range?
            }
            double fnear = getValueAtSample(inear, ilevel, unit);
            if (!NUMUtils.defined(fnear)) {
                return NUMundefined;   // function value not defined?
            }
            if (ifar < 0 || ifar >= this.n) {
                return fnear;   // at edge? Extrapolate
            }
            double ffar = getValueAtSample(ifar, ilevel, unit);
            if (!NUMUtils.defined(ffar)) {
                return fnear;   // neighbour undefined? Extrapolate
            }
            return fnear + phase * (ffar - fnear);   // interpolate
        }
        return getValueAtSample(SampledUtils.xToNearestIndex(this.x1, this.dx, x), ilevel, unit);
    }

    public Pair<Double, Double> getSumAndDefinitionRange(double xmin, double xmax, int ilevel, int unit, boolean interpolate) {
        /*
        This function computes the area under the linearly interpolated curve between xmin and xmax.
        Outside [x1-dx/2, xN+dx/2], the curve is undefined and neither times nor values are counted.
        In [x1-dx/2,x1] and [xN,xN+dx/2], the curve is linearly extrapolated.
         */
        int imin, imax, isamp;
        double sum = 0.0, definitionRange = 0.0;

        if (xmin >= xmax) {
            xmin = this.xmin;
            xmax = this.xmax;
        }

        Triple<Boolean, Double, Double> triInter = NUMUtils.Function_intersectRangeWithDomain(contour, xmin, xmax);
        xmin = triInter.second;
        xmax = triInter.third;
        
        if (triInter.first) {
            if (interpolate) {
                Triple<Integer, Integer, Integer> tri = SampledUtils.getWindowSamples(this.x1, this.dx, this.n, xmin, xmax);
                
                imin = tri.first;
                imax = tri.second;
                
                if (tri.third != 0) {
                    double leftEdge = this.x1 - 0.5 * this.dx;
                    double rightEdge = leftEdge + this.n * this.dx;
                    
                    
                    for (isamp = imin; isamp <= imax; isamp++) {
                        double value = getValueAtSample(isamp, ilevel, unit);   /* A fast way to integrate a linearly interpolated curve; works everywhere except at the edges. */
                        if (NUMUtils.defined(value)) {
                            definitionRange += 1.0;
                            sum += value;
                        }           
                    }
                    
                    /*
                     * Corrections within the first and last sampling intervals.
                     */
                    if (xmin > leftEdge) {   /* Otherwise, constant extrapolation over 0.5 sample is OK. */
                        double phase = (this.x1 + imin * this.dx - xmin) / this.dx;   /* This fraction of sampling interval is still to be determined. */
                        double rightValue = getValueAtSample(imin, ilevel, unit);
                        double leftValue = getValueAtSample(imin - 1, ilevel, unit);
                        if (NUMUtils.defined(rightValue)) {
                            definitionRange -= 0.5;   /* Delete constant extrapolation over 0.5 sample. */
                            sum -= 0.5 * rightValue;
                            if (NUMUtils.defined(leftValue)) {
                                definitionRange += phase;   /* Add current fraction. */
                                sum += phase * (rightValue + 0.5 * phase * (leftValue - rightValue));   /* Interpolate to outside sample. */
                            } else {
                                if (phase > 0.5) {
                                    phase = 0.5;
                                }
                                definitionRange += phase;   /* Add current fraction, but never more than 0.5. */
                                sum += phase * rightValue;
                            }
                        } else if (NUMUtils.defined(leftValue) && phase > 0.5) {
                            definitionRange += phase - 0.5;
                            sum += (phase - 0.5) * leftValue;
                        }
                    }
                    if (xmax < rightEdge) {   /* Otherwise, constant extrapolation is OK. */
                        double phase = (xmax - (this.x1 + imax * this.dx)) / this.dx;   /* This fraction of sampling interval is still to be determined. */
                        double leftValue = getValueAtSample(imax, ilevel, unit);
                        double rightValue = getValueAtSample(imax + 1, ilevel, unit);
                        if (NUMUtils.defined(leftValue)) {
                            definitionRange -= 0.5;   /* Delete constant extrapolation over 0.5 sample. */
                            sum -= 0.5 * leftValue;
                            if (NUMUtils.defined(rightValue)) {
                                definitionRange += phase;   /* Add current fraction. */
                                sum += phase * (leftValue + 0.5 * phase * (rightValue - leftValue));   /* Interpolate to outside sample. */
                            } else {
                                if (phase > 0.5) {
                                    phase = 0.5;
                                }
                                definitionRange += phase;   /* Add current fraction, but never more than 0.5. */
                                sum += phase * leftValue;
                            }
                        } else if (NUMUtils.defined(rightValue) && phase > 0.5) {
                            definitionRange += phase - 0.5;
                            sum += (phase - 0.5) * rightValue;
                        }
                    }
                } else {
                    /* No sample centres between xmin and xmax. */
                    /*
                     * Try to return the mean of the interpolated values at these two points.
                     * Thus, a small (xmin, xmax) range gives the same value as the (xmin+xmax)/2 point.
                     */
                    double leftValue = getValueAtSample(imax, ilevel, unit);
                    double rightValue = getValueAtSample(imin, ilevel, unit);
                    double phase1 = (xmin - (this.x1 + imax * this.dx)) / this.dx;
                    double phase2 = (xmax - (this.x1 + imax * this.dx)) / this.dx;
                    if (imin == imax + 1) {   /* Not too far from sample definition region. */
                        if (NUMUtils.defined(leftValue)) {
                            if (NUMUtils.defined(rightValue)) {
                                definitionRange += phase2 - phase1;
                                sum += (phase2 - phase1) * (leftValue + 0.5 * (phase1 + phase2) * (rightValue - leftValue));
                            } else if (phase1 < 0.5) {
                                if (phase2 > 0.5) {
                                    phase2 = 0.5;
                                }
                                definitionRange += phase2 - phase1;
                                sum += (phase2 - phase1) * leftValue;
                            }
                        } else if (NUMUtils.defined(rightValue) && phase2 > 0.5) {
                            if (phase1 < 0.5) {
                                phase1 = 0.5;
                            }
                            definitionRange += phase2 - phase1;
                            sum += (phase2 - phase1) * rightValue;
                        }
                    }
                }
            } else {   /* No interpolation. */
                double rimin = SampledUtils.xToIndex(this.x1, this.dx, xmin), rimax = SampledUtils.xToIndex(this.x1, this.dx, xmax);
                if (rimax >= -0.5 && rimin < this.n - 0.5) {
                    imin = rimin < -0.5 ? -1 : (int) Math.floor(rimin + 0.5);
                    imax = rimax >= this.n - 0.5 ? this.n : (int) Math.floor(rimax + 0.5);
                    for (isamp = imin + 1; isamp < imax; isamp++) {
                        double value = getValueAtSample(isamp, ilevel, unit);
                        if (NUMUtils.defined(value)) {
                            definitionRange += 1.0;
                            sum += value;
                        }
                    }
                    if (imin == imax) {
                        double value = getValueAtSample(imin, ilevel, unit);
                        if (NUMUtils.defined(value)) {
                            double phase = rimax - rimin;
                            definitionRange += phase;
                            sum += phase * value;
                        }
                    } else {
                        if (imin >= 0) {
                            double value = getValueAtSample(imin, ilevel, unit);
                            if (NUMUtils.defined(value)) {
                                double phase = imin - rimin + 0.5;
                                definitionRange += phase;
                                sum += phase * value;
                            }
                        }
                        if (imax < contour.size()) {
                            double value = getValueAtSample(imax, ilevel, unit);
                            if (NUMUtils.defined(value)) {
                                double phase = rimax - imax + 0.5;
                                definitionRange += phase;
                                sum += phase * value;
                            }
                        }
                    }
                }
            }
        }
        return new Pair(sum, definitionRange);
    }

    public double getMeanStrength(double tmin, double tmax, int unit) {
        return getMean(tmin, tmax, Pitch_LEVEL_STRENGTH, unit, true);
    }
    
    public double getStdStrength(double tmin, double tmax, int unit) {
        return getStandardDeviation(tmin, tmax, Pitch_LEVEL_STRENGTH, unit, true);
    }

    public double getMean(double tmin, double tmax, int unit) {
        return getMean(tmin, tmax, Pitch_LEVEL_FREQUENCY, unit, true);
    }

    public double getMean(double xmin, double xmax, int ilevel, int unit, boolean interpolate) {
        double sum, definitionRange;
        Pair<Double, Double> pair = getSumAndDefinitionRange(xmin, xmax, ilevel, unit, interpolate);
        sum = pair.first;
        definitionRange = pair.second;
        
        return definitionRange <= 0.0 ? NUMundefined : sum / definitionRange;
    }

    //Checked
    public double getStandardDeviation(double xmin, double xmax, int ilevel, int unit, boolean interpolate) {
        double sum, definitionRange;
        Pair<Double, Double> pair = getSumAndDefinitionRange(xmin, xmax, ilevel, unit, interpolate);
        sum = pair.first;
        definitionRange = pair.second;
        if (definitionRange < 1.0) {
            //System.out.println("xmin: "+xmin+" xmax: "+xmax+" sum:"+sum+" definitionRange:"+definitionRange);
            return Double.NaN;
        }else if(definitionRange < 2.0){
            return 0.0;
        }
        Pair<Double, Double> pair2 = getSum2AndDefinitionRange(xmin, xmax, ilevel, unit, sum / definitionRange, interpolate);
        return Math.sqrt(pair2.first / (pair2.second - 1.0));
    }

    public double getStandardDeviation(double tmin, double tmax, int unit) {
        return getStandardDeviation(tmin, tmax, Pitch_LEVEL_FREQUENCY, unit, true);
    }

    public Pair<Double, Double> getSum2AndDefinitionRange(double xmin, double xmax, int ilevel, int unit, double mean, boolean interpolate) {
        /*
         * This function computes the area under the linearly interpolated squared difference curve between xmin and xmax.
         * Outside [x1-dx/2, xN+dx/2], the curve is undefined and neither times nor values are counted.
         * In [x1-dx/2,x1] and [xN,xN+dx/2], the curve is linearly extrapolated.
         */
        int imin, imax, isamp;
        double sum2 = 0.0, definitionRange = 0.0;
        if (xmin >= xmax) {
            xmin = this.xmin;
            xmax = this.xmax;
        }

        Triple<Boolean, Double, Double> triInter = NUMUtils.Function_intersectRangeWithDomain(contour, xmin, xmax);
        xmin = triInter.second;
        xmax = triInter.third;
        if (triInter.first) {
            if (interpolate) {
                Triple<Integer, Integer, Integer> tri = SampledUtils.getWindowSamples(this.x1, this.dx, this.n, xmin, xmax);
                imin = tri.first;
                imax = tri.second;

                if (tri.third != 0) {
                    double leftEdge = this.x1 - 0.5 * this.dx, rightEdge = leftEdge + this.n * this.dx;
                    for (isamp = imin; isamp <= imax; isamp++) {
                        double value = getValueAtSample(isamp, ilevel, unit);   // a fast way to integrate a linearly interpolated curve; works everywhere except at the edges
                        if (NUMUtils.defined(value)) {
                            value -= mean;
                            value *= value;
                            definitionRange += 1.0;
                            sum2 += value;
                        }
                    }
                    /*
                     * Corrections within the first and last sampling intervals.
                     */
                    if (xmin > leftEdge) {   // otherwise, constant extrapolation over 0.5 sample is OK
                        double phase = (this.x1 + imin* this.dx - xmin) / this.dx;   // this fraction of sampling interval is still to be determined
                        double rightValue = getValueAtSample(imin, ilevel, unit);
                        double leftValue = getValueAtSample(imin - 1, ilevel, unit);
                        if (NUMUtils.defined(rightValue)) {
                            rightValue -= mean;
                            rightValue *= rightValue;
                            definitionRange -= 0.5;   // delete constant extrapolation over 0.5 sample
                            sum2 -= 0.5 * rightValue;
                            if (NUMUtils.defined(leftValue)) {
                                leftValue -= mean;
                                leftValue *= leftValue;
                                definitionRange += phase;   // add current fraction
                                sum2 += phase * (rightValue + 0.5 * phase * (leftValue - rightValue));   // interpolate to outside sample
                            } else {
                                if (phase > 0.5) {
                                    phase = 0.5;
                                }
                                definitionRange += phase;   // add current fraction, but never more than 0.5
                                sum2 += phase * rightValue;
                            }
                        } else if (NUMUtils.defined(leftValue) && phase > 0.5) {
                            leftValue -= mean;
                            leftValue *= leftValue;
                            definitionRange += phase - 0.5;
                            sum2 += (phase - 0.5) * leftValue;
                        }
                    }
                    if (xmax < rightEdge) {   // otherwise, constant extrapolation is OK
                        double phase = (xmax - (this.x1 + imax * this.dx)) / this.dx;   // this fraction of sampling interval is still to be determined
                        double leftValue = getValueAtSample(imax, ilevel, unit);
                        double rightValue = getValueAtSample(imax + 1, ilevel, unit);
                        if (NUMUtils.defined(leftValue)) {
                            leftValue -= mean;
                            leftValue *= leftValue;
                            definitionRange -= 0.5;   // delete constant extrapolation over 0.5 sample
                            sum2 -= 0.5 * leftValue;
                            if (NUMUtils.defined(rightValue)) {
                                rightValue -= mean;
                                rightValue *= rightValue;
                                definitionRange += phase;   // add current fraction
                                sum2 += phase * (leftValue + 0.5 * phase * (rightValue - leftValue));   // interpolate to outside sample
                            } else {
                                if (phase > 0.5) {
                                    phase = 0.5;
                                }
                                definitionRange += phase;   // add current fraction, but never more than 0.5
                                sum2 += phase * leftValue;
                            }
                        } else if (NUMUtils.defined(rightValue) && phase > 0.5) {
                            rightValue -= mean;
                            rightValue *= rightValue;
                            definitionRange += phase - 0.5;
                            sum2 += (phase - 0.5) * rightValue;
                        }
                    }
                } else {   // no sample centres between xmin and xmax
				/*
                     * Try to return the mean of the interpolated values at these two points.
                     * Thus, a small (xmin, xmax) range gives the same value as the (xmin+xmax)/2 point.
                     */
                    double leftValue = getValueAtSample(imax, ilevel, unit);
                    double rightValue = getValueAtSample(imin, ilevel, unit);
                    double phase1 = (xmin - (this.x1 + imax * this.dx)) / this.dx;
                    double phase2 = (xmax - (this.x1 + imax * this.dx)) / this.dx;
                    if (imin == imax + 1) {   // not too far from sample definition region
                        if (NUMUtils.defined(leftValue)) {
                            leftValue -= mean;
                            leftValue *= leftValue;
                            if (NUMUtils.defined(rightValue)) {
                                rightValue -= mean;
                                rightValue *= rightValue;
                                definitionRange += phase2 - phase1;
                                sum2 += (phase2 - phase1) * (leftValue + 0.5 * (phase1 + phase2) * (rightValue - leftValue));
                            } else if (phase1 < 0.5) {
                                if (phase2 > 0.5) {
                                    phase2 = 0.5;
                                }
                                definitionRange += phase2 - phase1;
                                sum2 += (phase2 - phase1) * leftValue;
                            }
                        } else if (NUMUtils.defined(rightValue) && phase2 > 0.5) {
                            rightValue -= mean;
                            rightValue *= rightValue;
                            if (phase1 < 0.5) {
                                phase1 = 0.5;
                            }
                            definitionRange += phase2 - phase1;
                            sum2 += (phase2 - phase1) * rightValue;
                        }
                    }
                }
            } else {   // no interpolation
                double rimin = SampledUtils.xToIndex(this.x1, this.dx, xmin), rimax = SampledUtils.xToIndex(this.x1, this.dx, xmax);
                if (rimax >= -0.5 && rimin < this.n - 0.5) {
                    imin = rimin < -0.5 ? -1 : (int) Math.floor(rimin + 0.5);
                    imax = rimax >= this.n - 0.5 ? this.n : (int) Math.floor(rimax + 0.5);
                    for (isamp = imin + 1; isamp < imax; isamp++) {
                        double value = getValueAtSample(isamp, ilevel, unit);
                        if (NUMUtils.defined(value)) {
                            value -= mean;
                            value *= value;
                            definitionRange += 1.0;
                            sum2 += value;
                        }
                    }
                    if (imin == imax) {
                        double value = getValueAtSample(imin, ilevel, unit);
                        if (NUMUtils.defined(value)) {
                            double phase = rimax - rimin;
                            value -= mean;
                            value *= value;
                            definitionRange += phase;
                            sum2 += phase * value;
                        }
                    } else {
                        if (imin >= 0) {
                            double value = getValueAtSample(imin, ilevel, unit);
                            if (NUMUtils.defined(value)) {
                                double phase = imin - rimin + 0.5;
                                value -= mean;
                                value *= value;
                                definitionRange += phase;
                                sum2 += phase * value;
                            }
                        }
                        if (imax < contour.size()) {
                            double value = getValueAtSample(imax, ilevel, unit);
                            if (NUMUtils.defined(value)) {
                                double phase = rimax - imax + 0.5;
                                value -= mean;
                                value *= value;
                                definitionRange += phase;
                                sum2 += phase * value;
                            }
                        }
                    }
                }
            }
        }
        return new Pair(sum2, definitionRange);

    }
    
    public double[] extractVoicedFrames(){
        
        double[] voicedFrames = new double[countVoicedFrames()];
        
        int numberOfDefinedSamples = 0;
        for (int isamp = 0; isamp < this.n; isamp++) {
            double value = getValueAtSample(isamp,Pitch_LEVEL_FREQUENCY, kPitch_unit_HERTZ);
            if (NUMUtils.defined(value)) {
                voicedFrames[numberOfDefinedSamples++] = value;
            }
        }
        
        return voicedFrames;
    }

    public int countVoicedFrames() {
        return countDefinedSamples(Pitch_LEVEL_FREQUENCY, kPitch_unit_HERTZ);

    }
    
    

    public int countDefinedSamples(int ilevel, int unit) {

        int numberOfDefinedSamples = 0;
        for (int isamp = 0; isamp < this.n; isamp++) {
            double value = getValueAtSample(isamp, ilevel, unit);
            if (!NUMUtils.defined(value)) {
                continue;
            }
            numberOfDefinedSamples += 1;
        }
        return numberOfDefinedSamples;
    }

    public double getStart() {
        return x1;
    }

    //public double getEnd() {
    //    return xN;
    //}

    public double getStep() {
        return dx;
    }

    public double getCeiling() {
        return ceiling;
    }
    
    public int size(){
        return n;
    }

    public ArrayList<PitchFrame> getFrames() {
        return pitchFrames;
    }

    public Double getQuantile(double quantile) {
        return getQuantile(this.xmin, this.xmax, quantile, Pitch_LEVEL_FREQUENCY, kPitch_unit_HERTZ);
    }

    public Double getQuantile(double xmin, double xmax, double quantile, int ilevel, int unit) {

        double[] values = new double[this.n];
        if (xmin >= xmax) {
            xmin = this.xmin;
            xmax = this.xmax;
        }

        if (xmin == xmax) {
            return Pitch.NUMundefined;
        }

        if (xmin < xmax) {
            if (xmin < this.xmin) {
                xmin = this.xmin;   // intersect requested range with logical domain
            }
            if (xmax > this.xmax) {
                xmax = this.xmax;
            }
            if (xmax <= xmin) {
                return Pitch.NUMundefined;   // requested range and logical domain do not intersect
            }
        } else {
            if (xmax < this.xmin) {
                xmax = this.xmin;   // intersect requested range with logical domain
            }
            if (xmin > this.xmax) {
                xmin = this.xmax;
            }
            if (xmin <= xmax) {
                return Pitch.NUMundefined;   // requested range and logical domain do not intersect
            }
        }

        int imin, imax, numberOfDefinedSamples = 0;

        double rixmin = Math.ceil((xmin - this.x1) / this.dx);
        double rixmax = Math.floor((xmax - this.x1) / this.dx);
        imin = rixmin < 0 ? 0 : (int) rixmin;
        imax = rixmax >= (double) this.n ? this.n - 1 : (int) rixmax;
        
        for (int i = imin; i <= imax; i++) {
            double value = getValueAtSample(i,ilevel,unit);
            if (!Double.isNaN(value)) {
                values[numberOfDefinedSamples++] = value;
            }
        }

        double[] endValues = Arrays.copyOfRange(values, 0, numberOfDefinedSamples);
        
        
        double result = Pitch.NUMundefined;

        if (numberOfDefinedSamples >= 1) {
            Arrays.sort(endValues);
            

            double place = quantile * endValues.length - 0.5;
            int left = (int) Math.round(place);

            if (numberOfDefinedSamples == 1) {
                return endValues[0];
            }
            if (left < 0) {
                left = 0;
            }
            if (left >= numberOfDefinedSamples - 1) {
                left = numberOfDefinedSamples - 2;
            }
            if (endValues[left + 1] == endValues[left]) {
                return endValues[left];
            }
            return endValues[left] + (place - left) * (endValues[left + 1] - endValues[left]);
        }

        return result;
    }

    


    public double getMaximum(double xmin, double xmax, int ilevel, int unit, int interpolate) {
        int imin, imax, i;
        
        double maximum = -1e301;
       

        if (!NUMUtils.defined(xmin) || !NUMUtils.defined(xmax)) {
            maximum =  NUMundefined;
            return maximum;
        }
        
        //NUMUtils.Function_unidirectionalAutowindow(me,  & xmin,  & xmax);
        
        Triple<Boolean, Double, Double> tri = NUMUtils.Function_intersectRangeWithDomain(contour,xmin,xmax);
        xmin = tri.second;
        xmax = tri.third;
        
        if (!tri.first) {
            
            maximum = NUMundefined;   // requested range and logical domain do not intersect
            return maximum;
        }
        
        Triple <Integer, Integer, Integer> tri2 = SampledUtils.getWindowSamples(this.x1,this.dx,this.n,xmin,xmax);
        imin = tri2.first;
        imax = tri2.second;
        if (tri.third==0) {
            /*
             * No sample centres between tmin and tmax.
             * Try to return the greater of the values at these two points.
             */
            double fleft = getValueAtX(xmin, ilevel, unit, interpolate);
            double fright = getValueAtX( xmax, ilevel, unit, interpolate);
            if (NUMUtils.defined(fleft) && fleft > maximum) 
                maximum = fleft;
            if (NUMUtils.defined(fright) && fright > maximum)
                maximum = fright;
            
        } else {
            for (i = imin; i <= imax; i++) {
                double fmid = getValueAtSample(i, ilevel, unit);
                if (!NUMUtils.defined(fmid)) {
                    continue;
                }
                if (interpolate == 0) {
                    if (fmid > maximum)
                        maximum = fmid;
                    
                } else {
                    /*
                     * Try an interpolation, possibly even taking into account a sample just outside the selection.
                     */
                    double fleft = i <= 0 ? NUMundefined : getValueAtSample(i - 1, ilevel, unit);
                    double fright = i >= this.n-1 ? NUMundefined : getValueAtSample(i + 1, ilevel, unit);
                    if (!NUMUtils.defined(fleft) || !NUMUtils.defined(fright)) {
                        if (fmid > maximum)
                            maximum = fmid;
                        
                    } else if (fmid > fleft && fmid >= fright) {
                        double[] y = new double[4];
                        double localMaximum;
                        y [0] = fleft; y [1] = fmid; y [2] = fright;
                        Pair<Double, Double> pair = NUMUtils.improveMaximum (y, 3, 2, NUMUtils.NUM_PEAK_INTERPOLATE_PARABOLIC);
                        localMaximum = pair.first;
                        
                        if (localMaximum > maximum){
                            maximum = localMaximum;
            
                        }
                        
                    }
                }
            }
            
            /* Check boundary values. */
            if (interpolate!=0) {
                double fleft = getValueAtX(xmin, ilevel, unit, interpolate);
                double fright = getValueAtX( xmax, ilevel, unit, interpolate);
                if (NUMUtils.defined(fleft) && fleft > maximum)
                    maximum = fleft;
                
                
                if (NUMUtils.defined(fright) && fright > maximum) {
                    maximum = fright;
                }
            }
        }
        
        if (maximum == -1e301) {
            maximum = NUMundefined;
        }
        
        return maximum;
    }
    
    
    public double getMinimum(double xmin, double xmax, int ilevel, int unit, int interpolate) {
        int imin, imax, i;
        
        double minimum = 1e301;
        if (Double.isNaN(xmin) || Double.isNaN(xmax)) {
            minimum  = NUMundefined;
            
            return minimum;
        }
        
        //NUMUtils.Function_unidirectionalAutowindow(me,  & xmin,  & xmax);
        
        Triple<Boolean, Double, Double> tri = NUMUtils.Function_intersectRangeWithDomain(contour,xmin,xmax);
        xmin = tri.second;
        xmax = tri.third;
        
        if (!tri.first) {
            minimum =  NUMundefined;   // requested range and logical domain do not intersect
            return minimum;
        }
        Triple <Integer, Integer, Integer> tri2 = SampledUtils.getWindowSamples(this.x1,this.dx,this.n,xmin,xmax);
        imin = tri2.first;
        imax = tri2.second;
        if (tri.third==0) {
            /*
             * No sample centres between tmin and tmax.
             * Try to return the greater of the values at these two points.
             */
            double fleft = getValueAtX(xmin, ilevel, unit, interpolate);
            double fright = getValueAtX( xmax, ilevel, unit, interpolate);
            if (NUMUtils.defined(fleft) && fleft < minimum) 
                minimum = fleft;
            
            if (NUMUtils.defined(fright) && fright < minimum)
                minimum = fright;
            
        } else {
            for (i = imin; i <= imax; i++) {
                double fmid = getValueAtSample(i, ilevel, unit);
                if (!NUMUtils.defined(fmid)) {
                    continue;
                }
                if (interpolate == 0) {
                    if (fmid < minimum)
                        minimum = fmid;
                    
                } else {
                    /*
                     * Try an interpolation, possibly even taking into account a sample just outside the selection.
                     */
                    double fleft = i <= 0 ? NUMundefined : getValueAtSample(i - 1, ilevel, unit);
                    double fright = i >= this.n-1 ? NUMundefined : getValueAtSample(i + 1, ilevel, unit);
                    if (!NUMUtils.defined(fleft) || !NUMUtils.defined(fright)) {
                        if (fmid < minimum)
                            minimum = fmid;
                        
                    }else if (fmid < fleft && fmid <= fright) {
                        double[] y = new double[4];
                        double localMinimum;
                        y [0] = fleft; y [1] = fmid; y [2] = fright;
                        Pair<Double, Double> pair = NUMUtils.improveMinimum (y, 3, 2, NUMUtils.NUM_PEAK_INTERPOLATE_PARABOLIC);
                        localMinimum = pair.first;
                        
                        if (localMinimum < minimum){
                            minimum = localMinimum;
                        }
                    }
                }
            }
            /* Check boundary values. */
            if (interpolate!=0) {
                double fleft = getValueAtX(xmin, ilevel, unit, interpolate);
                double fright = getValueAtX( xmax, ilevel, unit, interpolate);
                if (NUMUtils.defined(fleft) && fleft < minimum)
                    minimum = fleft;
                
                
                if (NUMUtils.defined(fright) && fright < minimum) {
                    minimum = fright;
                }
            }
        }
        if (minimum == 1e301) {
            minimum = NUMundefined;
        }
        
        return minimum;
    }
}
