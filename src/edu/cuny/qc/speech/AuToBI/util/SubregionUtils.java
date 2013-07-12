/*  SubregionUtils.java

    Copyright 2009-2010 Andrew Rosenberg

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
package edu.cuny.qc.speech.AuToBI.util;

import edu.cuny.qc.speech.AuToBI.core.AuToBIException;
import edu.cuny.qc.speech.AuToBI.core.Region;
import edu.cuny.qc.speech.AuToBI.core.WavData;
import edu.cuny.qc.speech.AuToBI.core.Word;
import edu.cuny.qc.speech.AuToBI.featureextractor.FeatureExtractorException;

import java.util.List;

/**
 * SubregionUtils is a utility class to house static functions for the processing of subregions.
 */
public class SubregionUtils {

  // Utility classes cannot be constructed.
  private void SubregionUtils() {
    throw new AssertionError();
  }

  /**
   * Aligns the subregion that is identified within a word to be the representative pseudosyllable.
   * <p/>
   * We say that the subregion that covers the longest duration of a word to be representative.
   * <p/>
   * This is an approximation of the accent bearing syllable, and may introduce errors.
   * <p/>
   * This is used for pitch accent type classification.
   * <p/>
   * NOTE: if pseudosyllables are very long relative to the length of words, the same pseudosyllable can be assigned to
   * the same word.
   *
   * @param words        the words
   * @param subregions   the pseudosyllables
   * @param feature_name the feature name for storing the selected pseudosyllable region
   */
  public static void alignLongestSubregionsToWords(List<Word> words, List<Region> subregions,
                                                   String feature_name) {
    int i = 0;
    for (Word w : words) {

      // Always try the last subregionregion.
      while (i >= subregions.size() && i > 0) {
        --i;
      }
      Region current_subregion = subregions.get(i);
      while (current_subregion.getStart() > w.getStart() && i > 0) {
        current_subregion = subregions.get(--i);
      }

      while (current_subregion != null && current_subregion.getEnd() <= w.getStart()) {
        ++i;
        if (i == subregions.size()) {
          // There are no subregions that are within the current word.
          current_subregion = null;
        } else {
          current_subregion = subregions.get(i);
        }
      }

      if (current_subregion != null) {
        Region best_subregion = current_subregion;
        double max_overlap = -Double.MAX_VALUE;
        while (current_subregion.getStart() < w.getEnd() && i < subregions.size()) {
          // Calculate the amount of overlapping material
          double overlap = Math.min(w.getEnd(), current_subregion.getEnd()) -
              Math.max(w.getStart(), current_subregion.getStart());
          if (overlap > max_overlap) {
            max_overlap = overlap;
            best_subregion = current_subregion;
          }
          ++i;
          if (i < subregions.size())
            current_subregion = subregions.get(i);
        }
        w.setAttribute(feature_name, best_subregion);
      }
    }
  }

  /**
   * Given a string describing the subregion length, return the desired number of seconds.
   * <p/>
   * Currently only parses seconds and miliseconds.
   * <p/>
   * 400ms 1s
   *
   * @param subregion_name the name of the subregion to parse.
   * @return the number of seconds described by the string
   * @throws FeatureExtractorException if the subregion label is unparseable
   */
  public static Double parseSubregionName(String subregion_name) throws FeatureExtractorException {

    if (!subregion_name.matches("^\\d+(ms|s)$")) {
      throw new FeatureExtractorException("Cannot parse the subregion: " + subregion_name);
    }

    boolean milliseconds = false;
    if (subregion_name.matches(".*ms$")) {
      milliseconds = true;
    }

    String number;
    if (milliseconds) {
      number = subregion_name.substring(0, subregion_name.indexOf('m'));
      return Double.parseDouble(number) / 1000;
    } else {
      number = subregion_name.substring(0, subregion_name.indexOf('s'));
      return Double.parseDouble(number);
    }
  }

  /**
   * Copies a feature from a region to the specified subregion.
   *
   * @param regions             the regions to copy attributes from
   * @param subregion_attribute the name of the subregion attribute that is the destination of the attribute
   * @param attribute_name      the name of the copied attribute
   */
  public static void assignFeatureToSubregions(List regions, String subregion_attribute,
                                               String attribute_name) {
    for (Region r : (List<Region>) regions) {
      if (r.hasAttribute(subregion_attribute) && r.hasAttribute(attribute_name)) {
        Region subregion = (Region) r.getAttribute(subregion_attribute);
        subregion.setAttribute(attribute_name, r.getAttribute(attribute_name));
      }
    }
  }

  /**
   * Copies a feature from a region to the specified list of subregions.
   *
   * @param regions             the regions to copy attributes from
   * @param subregion_attribute the name of the subregion attribute that is the destination of the attribute
   * @param attribute_name      the name of the copied attribute
   */
  public static void assignFeatureToAllSubregions(List regions, String subregion_attribute,
                                                  String attribute_name) {
    for (Region r : (List<Region>) regions) {
      if (r.hasAttribute(subregion_attribute) && r.hasAttribute(attribute_name)) {
        for (Region subregion : (List<Region>) r.getAttribute(subregion_attribute)) {
          subregion.setAttribute(attribute_name, r.getAttribute(attribute_name));
        }
      }
    }
  }

  /**
   * Extracts a subregion from a wave file.  Creating a new WavData object containing a subset of the original wav
   * file.
   * <p/>
   * The frame rate and sample size are unchanged from the original.
   *
   * @param wav_data the input WavData object
   * @param start    the start time
   * @param end      the end time
   * @return a new WavData object
   * @throws edu.cuny.qc.speech.AuToBI.core.AuToBIException
   *          if the sizes are inappropriate
   */
  public static WavData getSlice(WavData wav_data, double start, double end) throws AuToBIException {
    if (start >= end) {
      throw new AuToBIException("Negative Sized Slice requested.");
    }
    if (end < wav_data.t0 || start > wav_data.t0 + wav_data.getNumSamples() / wav_data.sampleRate) {
      return null;
    }
    WavData sub_wav = new WavData();
    sub_wav.numberOfChannels = wav_data.numberOfChannels;
    sub_wav.sampleRate = wav_data.sampleRate;
    sub_wav.sampleSize = wav_data.sampleSize;

    int start_idx = Math.max(0, (int) Math.ceil((start - wav_data.t0) * sub_wav.sampleRate));
    sub_wav.t0 = wav_data.t0 + start_idx / sub_wav.sampleRate;
    int end_idx =
        Math.min(wav_data.samples[0].length - 1, (int) Math.floor((end - wav_data.t0) * sub_wav.sampleRate));

    int num_frames = end_idx - start_idx + 1;
    sub_wav.samples = new double[wav_data.numberOfChannels][num_frames];

    for (int channel = 0; channel < wav_data.numberOfChannels; ++channel) {
      double[] norm = wav_data.getSamples(channel);
      System.arraycopy(norm, start_idx, sub_wav.samples[channel], start_idx - start_idx, num_frames);
    }

    return sub_wav;
  }
  
  /**
   * Extracts a subregion from a wave file.  Creating a new WavData object containing a subset of the original wav
   * file.
   * <p/>
   * The frame rate and sample size are unchanged from the original.
   *
   * @param wav_data the input WavData object
   * @param start    the start time
   * @param end      the end time
   * @return a new WavData object
   * @throws edu.cuny.qc.speech.AuToBI.core.AuToBIException
   *          if the sizes are inappropriate
   */
  public static WavData concatSlices(WavData[] slices) throws AuToBIException {
    WavData wd = new WavData();
    
    int totalFrames = 0;
    for(int i=0;i<slices.length;i++)
        totalFrames+=slices[i].samples[0].length;
    
    
    
    
    wd.numberOfChannels = slices[0].numberOfChannels;
    wd.sampleRate = slices[0].sampleRate;
    wd.sampleSize = slices[0].sampleSize;
    wd.t0 = 0;
    wd.samples = new double[wd.numberOfChannels][totalFrames];
    
    int acum_idx = 0;
    for(int n = 0; n < slices.length;n++){
        
        
        int num_frames = slices[n].samples[0].length;
        
        for (int channel = 0; channel < wd.numberOfChannels; ++channel) {
        
            double[] norm = slices[n].getSamples(channel);
            System.arraycopy(norm, 0, wd.samples[channel], acum_idx, num_frames);
        }
        acum_idx+=num_frames;
    }
    

    return wd;
  }
  
}
