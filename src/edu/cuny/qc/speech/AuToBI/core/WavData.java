/*  WavData.java

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

package edu.cuny.qc.speech.AuToBI.core;

import edu.cuny.qc.speech.AuToBI.Pitch;
import edu.cuny.qc.speech.AuToBI.util.SampledUtils;

/**
 * WavData is used to store Wav file data.
 */
public class WavData {
  public double[][] samples;   // Normalized Audio Data
  public int numberOfChannels; // Number of stored channels
  public int sampleSize;       // Size of each sample in bits
  public float sampleRate;     // Number of raw_samples per second.
  public double t0;            // The time of the first sample.
  private String filename;     // The filename containing the this audio data.

  /**
   * Constructs a new WavData object with no data.
   */
  public WavData() {
    this.t0 = 0.0;
  }

  /**
   * Gets the duration of the wav file in seconds.
   *
   * @return the duration of the file
   */
  public double getDuration() {
    return samples[0].length / sampleRate;
  }

  /**
   * Gets the frame size in seconds.
   *
   * @return the frame size
   */
  public float getFrameSize() {
    return 1 / sampleRate;
  }

  /**
   * Gets the total number of samples in the file.
   *
   * @return the number of samples
   */
  public int getNumSamples() {
    return samples[0].length;
  }

  /**
   * Retries a specific sample indexed by channel and index number.
   *
   * @param channel the desired channel
   * @param index   the desired index
   * @return the sample stored at the specified channel and index
   */
  public double getSample(int channel, int index) {
    return samples[channel][index];
  }

  /**
   * Retrieves the full list of samples stored in a give channel.
   *
   * @param channel the desired channel
   * @return the wav samples in the specified channel
   */
  public double[] getSamples(int channel) {
    return samples[channel];
  }

  /**
   * Get the filename where the data was found.
   *
   * @return the filename
   */
  public String getFilename() {
    return filename;
  }

  /**
   * Sets the filename.
   *
   * @param filename the filename
   */
  public void setFilename(String filename) {
    this.filename = filename;
  }
  
  public double getHannWindowedRms (double tmid, double widthLeft, double widthRight) {
      double sumOfSquares = 0.0, windowSumOfSquares = 0.0;
      
      Triple<Integer,Integer,Integer> tri = SampledUtils.getWindowSamples(t0,1.0/this.sampleRate,this.samples[0].length,tmid - widthLeft, tmid + widthRight);
      int imin = tri.first;
      int imax = tri.second;
      
      if (tri.third < 3) return Pitch.NUMundefined;
      for (int i = imin; i <= imax; i ++) {
          double t = this.t0 + (i - 1) / this.sampleRate;
          double width = t < tmid ? widthLeft : widthRight;
          double windowPhase = (t - tmid) / width;   /* in [-1 .. 1] */
          double window = 0.5 + 0.5 * Math.cos (Math.PI * windowPhase);   /* Hann */
          double windowedValue = ( this.numberOfChannels == 1 ? this.samples[0] [i] : 0.5 * (this.samples[0][i] + this.samples[1][i])) * window;
          sumOfSquares += windowedValue * windowedValue;
          windowSumOfSquares += window * window;
      }
      return Math.sqrt(sumOfSquares / windowSumOfSquares);
  }

  
}
