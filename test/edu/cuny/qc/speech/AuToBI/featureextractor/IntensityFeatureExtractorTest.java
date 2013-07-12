/*  IntensityFeatureExtractorTest.java

    Copyright 2012 Andrew Rosenberg

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
package edu.cuny.qc.speech.AuToBI.featureextractor;

import edu.cuny.qc.speech.AuToBI.core.*;
import edu.cuny.qc.speech.AuToBI.io.WavReader;
import org.junit.Test;

import javax.sound.sampled.UnsupportedAudioFileException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import static junit.framework.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Test class for IntensityFeatureExtractor
 *
 * @see IntensityFeatureExtractor
 */
public class IntensityFeatureExtractorTest {

  private static final String TEST_DIR = "/Users/andrew/code/AuToBI/release/test_data/";

  @Test
  public void testConstructorSetsExtractedFeaturesCorrectly() {
    IntensityFeatureExtractor fe =
        new IntensityFeatureExtractor("I");

    assertEquals(1, fe.getExtractedFeatures().size());
    assertTrue(fe.getExtractedFeatures().contains("I"));
  }

  @Test
  public void testConstructorSetsRequiredFeaturesCorrectly() {
    IntensityFeatureExtractor fe =
        new IntensityFeatureExtractor("I");

    assertEquals(1, fe.getRequiredFeatures().size());
    assertTrue(fe.getRequiredFeatures().contains("wav"));
  }

  @Test
  public void testExtractFeaturesExtractsFeatures() {
    IntensityFeatureExtractor fe =
        new IntensityFeatureExtractor("I");

    List<Region> regions = new ArrayList<Region>();
    Word w = new Word(0.0, 1.0, "test");
    regions.add(w);

    WavReader reader = new WavReader();

    WavData wav = null;
    try {
      wav = reader.read(TEST_DIR + "/test.wav");
    } catch (UnsupportedAudioFileException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    } catch (AuToBIException e) {
      e.printStackTrace();
    }
    w.setAttribute("wav", wav);

    try {
      fe.extractFeatures(regions);
      assertTrue(w.hasAttribute("I"));
    } catch (FeatureExtractorException e) {
      e.printStackTrace();
    }
  }

  @Test
  public void testExtractFeaturesExtractsFeaturesCorrectly() {

    IntensityFeatureExtractor fe =
        new IntensityFeatureExtractor("I");

    List<Region> regions = new ArrayList<Region>();
    Word w = new Word(0.0, 1.0, "test");
    regions.add(w);
    WavReader reader = new WavReader();

    WavData wav = null;
    try {
      wav = reader.read(TEST_DIR + "/test.wav");
    } catch (UnsupportedAudioFileException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    } catch (AuToBIException e) {
      e.printStackTrace();
    }
    w.setAttribute("wav", wav);

    try {
      fe.extractFeatures(regions);
      Contour c = (Contour) w.getAttribute("I");
      assertEquals(100, c.size());
      assertEquals(0.01, c.getStep(), 0.0001);
      assertEquals(0.005, c.getStart(), 0.0001);
      // Assume that the intensity extraction algorithm is tested in IntensityExtractor.
      // Here we'll make sure that the generated contour passes some sanity checks.
    } catch (FeatureExtractorException e) {
      e.printStackTrace();
    }
  }
}
