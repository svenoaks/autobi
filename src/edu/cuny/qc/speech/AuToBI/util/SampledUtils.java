/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.cuny.qc.speech.AuToBI.util;

import edu.cuny.qc.speech.AuToBI.Pitch;
import edu.cuny.qc.speech.AuToBI.core.Contour;
import edu.cuny.qc.speech.AuToBI.core.Pair;
import edu.cuny.qc.speech.AuToBI.core.Triple;

/**
 *
 * @author vitisoto
 */
public class SampledUtils {

    //Checked
    public static Triple<Integer, Integer, Integer> getWindowSamples(double start, double step, int numSamples, double xmin, double xmax) {

        double rixmin = Math.ceil((xmin - start) / step);
        double rixmax = Math.floor((xmax - start) / step);
        int ixmin = rixmin < 0.0 ? 0 : (int) rixmin;
        int ixmax = rixmax >= numSamples ? (numSamples - 1) : (int) rixmax;
        Triple<Integer, Integer, Integer> tri = new Triple(ixmin, ixmax, ixmin <= ixmax ? ixmax - ixmin + 1 : 0);
        return tri;
    }

    //Checked
    public static int xToHighIndex(double start, double step, double x) {
        return (int) Math.ceil((x - start) / step);
    }

    //Checked
    public static int xToLowIndex(double start, double step, double x) {
        return (int) Math.floor((x - start) / step);
    }

    //Checked
    public static double xToIndex(double start, double step, double x) {
        return (x - start) / step;
    }

    //Checked
    public static double indexToX(double start, double step, long i) {
        return start + i * step;
    }

    //Checked
    public static int xToNearestIndex(double start, double step, double x) {
        return (int) Math.floor((x - start) / step + 1.5) - 1;
    }

    
}
