/*  IntensityExtractor.java

    Copyright (c) 2009-2010 Andrew Rosenberg
    
    Based on Sound_to_Intensity.c distributed as part of the Praat package Copyright (C) 1992-2008 Paul Boersma

    This file is part of the AuToBI prosodic analysis package.

    AuToBI is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AuToBI is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AuToBI.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.cuny.qc.speech.AuToBI;

import edu.cuny.qc.speech.AuToBI.core.*;
import edu.cuny.qc.speech.AuToBI.io.WavReader;

import javax.sound.sampled.UnsupportedAudioFileException;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import java.io.IOException;
import java.io.File;
import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.util.ArrayList;

/**
 * Intensity extractor an intensity contour from a wave file.
 * <p/>
 * Implements Paul Boersma's Sound_To_Intensity algorithm which is included in the Praat tool.
 */
public class IntensityExtractor extends SampledDataAnalyzer {
  public IntensityExtractor(WavData wav) {
    this.wav = wav;
  }


  /**
   * Generate an intensity contour from a wav file using default parameters.
   * <p/>
   * The default settings are a minimum intensity of 75dB, and a timestep of 0.1 secs (100ms).
   *
   * @return a list of intensity points.
   */
  public Contour soundToIntensity() {
    return soundToIntensity(75, 0.01, true);
  }

  /**
   * Generate an intensity contour from a wav file using default parameters.
   * <p/>
   * Implementation of Paul Boersma's Sound_to_Intensity function
   * <p/>
   * <p/>
   * Uses the RMS value amplitude of the signal weighted by a bessel-function defined window.
   * <p/>
   * The minimum pitch is used to determine the analysis window size, and a default time step if one isn't provided.
   *
   * @param min_pitch              The minimum pitch.
   * @param time_step              the time_step for the analysis
   * @param subtract_mean_pressure Wether or not to subtract the mean pressure from the intensity contour.
   * @return a list of intensity points.
   */
  public Contour soundToIntensity(double min_pitch, double time_step, boolean subtract_mean_pressure) {
    if (time_step <= 0.0) time_step = 0.8 / min_pitch;   /* Default: four times oversampling Hanning-wise. */

    int i, iframe, numberOfFrames;
    double windowDuration = 6.4 / min_pitch, t0;
    double halfWindowDuration = 0.5 * windowDuration;
    int halfWindowSamples = (int) (halfWindowDuration / wav.getFrameSize());
    int windowSamples = halfWindowSamples * 2 + 1;
    double amplitude[] = new double[windowSamples];
    double window[] = new double[windowSamples];

    for (i = 0; i < windowSamples; i++) {
      double x = i * wav.getFrameSize() / halfWindowDuration, root = 1 - x * x;
      window[i] = root <= 0.0 ? 0.0 : besselI0((2 * Math.PI * Math.PI + 0.5) * Math.sqrt(root));
    }

    // Identify the number of frames and initial time.
    Pair<Integer, Double> pair = shortTermAnalysis(time_step, windowDuration);
    numberOfFrames = pair.first;
    t0 = pair.second;

    Contour contour = new Contour(t0, time_step, numberOfFrames);
    for (iframe = 0; iframe < numberOfFrames; iframe++) {
      double midTime = indexToX(t0, time_step, iframe);
      int midSample = xToNearestIndex(0, wav.getFrameSize(), midTime);
      int leftSample = Math.max(0, midSample - halfWindowSamples);
      int rightSample = Math.min(wav.getNumSamples(), midSample + halfWindowSamples + 1);
      double sumxw = 0.0, sumw = 0.0, intensity;
      if (leftSample < 1) leftSample = 0;
      if (rightSample > wav.getNumSamples()) rightSample = wav.getNumSamples();

      for (int channel = 0; channel < wav.numberOfChannels; channel++) {
        for (i = leftSample; i < rightSample; i++) {
          amplitude[i - leftSample] = wav.getSample(channel, i);
        }
        if (subtract_mean_pressure) {
          double sum = 0.0;
          for (i = leftSample; i < rightSample; i++) {
            sum += amplitude[i - leftSample];
          }
          double mean = sum / (rightSample - leftSample + 1);
          for (i = leftSample; i < rightSample; i++) {
            amplitude[i - leftSample] -= mean;
          }
        }
        for (i = leftSample; i < rightSample; i++) {
          sumxw += amplitude[i - leftSample] * amplitude[i - leftSample] * window[i - leftSample];
          sumw += window[i - leftSample];
        }
      }
      intensity = sumxw / sumw;
      if (intensity != 0.0) intensity /= 4e-10;
      contour.set(midTime, intensity < 1e-30 ? -300 : 10 * Math.log10(intensity));
    }

    return contour;
  }

  /**
   * Modified Bessel function I0. Abramowicz & Stegun, p. 378.
   * <p/>
   * In this class, this function is used to define and weight an analysis window.
   *
   * @param x The value to evaluate
   * @return The bessel function evaluated at x
   */
  private double besselI0(double x) {
    double t;
    if (x < 0.0) return besselI0(-x);
    if (x < 3.75) {
      /* Formula 9.8.1. Accuracy 1.6e-7. */
      t = x / 3.75;
      t *= t;
      return 1.0 + t * (3.5156229 + t * (3.0899424 + t * (1.2067492
          + t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))));
    }
    /*
      otherwise: x >= 3.75
    */
    /* Formula 9.8.2. Accuracy of the polynomial factor 1.9e-7. */
    t = 3.75 / x;   /* <= 1.0 */
    return Math.exp(x) / Math.sqrt(x) * (0.39894228 + t * (0.01328592
        + t * (0.00225319 + t * (-0.00157565 + t * (0.00916281
        + t * (-0.02057706 + t * (0.02635537 + t * (-0.01647633
        + t * 0.00392377))))))));
  }

  public static void main(String[] args) {
    File file = new File(args[0]);
    AudioInputStream soundIn = null;
    try {
      soundIn = AudioSystem.getAudioInputStream(new BufferedInputStream(new FileInputStream(file)));
    } catch (UnsupportedAudioFileException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }

    WavReader reader = new WavReader();
    WavData wav;
    try {

      if (args.length > 1)
        wav = reader.read(soundIn, Double.parseDouble(args[1]), Double.parseDouble(args[2]));

      else
        wav = reader.read(soundIn);
      System.out.println(wav.sampleRate);
      System.out.println(wav.sampleSize);
      System.out.println(wav.getFrameSize());
      System.out.println(wav.getDuration());

      IntensityExtractor intensityFactory = new IntensityExtractor(wav);
      Contour intensity = intensityFactory.soundToIntensity();

      System.out.println("intensity points:" + intensity.size());

      for (int i = 0; i < intensity.size(); ++i) {
        System.out
            .println(
                "intensity point[" + i + "]: " + intensity.get(i) + " -- " + intensity.timeFromIndex(i));
      }
    } catch (AuToBIException e) {
      e.printStackTrace();
    }

  }
  
  public ArrayList<Pair<Double,Double>> getListSoundsSilences(double minPitch, double timeStep, double silenceThresholdDB, double minSilenceDuration, double minSoundingDuration){
  
        Contour contour = soundToIntensity(minPitch,timeStep,true);
      
	double duration = this.wav.getDuration(), time;
        
        if (silenceThresholdDB >= 0) {
            System.err.println("The silence threshold w.r.t. the maximum intensity should be a negative number.");
	}
        
        if (minSilenceDuration > duration) {
            return null;
	}

	double intensity_max_db, intensity_min_db;
        intensity_max_db = contour.getMaximum();
        intensity_min_db = contour.getMinimum();
        
	double intensity_dbRange = intensity_max_db - intensity_min_db;

	if (intensity_dbRange < 10)
               System.err.println("The loudest and softest part in your sound only differ by "+intensity_dbRange+" dB.");

	double intensityThreshold = intensity_max_db - Math.abs(silenceThresholdDB);

	if (minSilenceDuration > duration || intensityThreshold < intensity_min_db) {
			return null;
	}

	boolean inSilenceInterval = contour.get(1) < intensityThreshold;
	long iinterval = 1;
        double label =-1;
        int numBoundaries = 0;
        ArrayList<Pair<Double, Double>> intervals = new ArrayList<Pair<Double,Double>>();
        
        for (long i = 2; i <= contour.size(); i++) {
		boolean addBoundary = false;
		if (contour.get(i) < intensityThreshold) {
			if (!inSilenceInterval) { // Start of silence
				addBoundary = true;
				inSilenceInterval = true;
				label = 1;
			}
                        
		} else {
			if (inSilenceInterval) { // End of silence
				addBoundary = true;
				inSilenceInterval = false;
				label = 0;
			}
		}

		if (addBoundary) {
                        numBoundaries++;
			time = contour.getStart() + (i - 1) * contour.getStep();
			//IntervalTier_addBoundaryUnsorted (it, iinterval, time, label);
                        intervals.add(new Pair<Double,Double>(time,label));
                        iinterval++;
		}
	}
        
	// (re)label last interval */
        label = inSilenceInterval ? 0 : 1;
        intervals.get(numBoundaries).second=label;
        
        
        // First remove short non-silence intervals in-between silence intervals and
	// then remove the remaining short silence intervals.
	// This works much better than first removing short silence intervals and
	// then short non-silence intervals.
        
        intervals = cutIntervalsMinDuration(intervals, 1, minSoundingDuration);
        intervals = cutIntervalsLabelMatch(intervals, 0);
        intervals = cutIntervalsMinDuration(intervals, 0, minSilenceDuration);
        intervals = cutIntervalsLabelMatch(intervals, 1);
        
	return intervals;
	
    }
  
  
    private ArrayList<Pair<Double,Double>> cutIntervalsMinDuration(ArrayList<Pair<Double,Double>> it, double label, double minDuration){
        
        double previous_time=0.0;
        ArrayList<Pair<Double,Double>> newIntervals = new ArrayList<Pair<Double,Double>>();
        
        for(int i=0;i<it.size();i++){
            if(it.get(i).second == label){
                double d = it.get(i).first-previous_time;
                if(d >= minDuration){
                    newIntervals.add(it.get(i));
                }
            }else{
                newIntervals.add(it.get(i));
            }
            previous_time = it.get(i).first;
        }
        
        return newIntervals;
    }
    
    private ArrayList<Pair<Double,Double>> cutIntervalsLabelMatch(ArrayList<Pair<Double,Double>> it, double label){
        
        int i=1;
        for(;;){
            if(it.get(i).second == it.get(i-1).second){
                it.remove(i-1);
            }else{
                i++;
            }
            if(i==it.size()){
                break;
            }
            
        }
        return it;
    }
  

}
