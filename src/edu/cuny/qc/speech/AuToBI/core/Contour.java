/*  Contour.java

    Copyright 2009-2011 Andrew Rosenberg

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
import java.util.Arrays;

/**
 * Contour is used to store time series data.
 * <p/>
 * Data is stored at an evenly sampled interval.
 * <p/>
 * Data is retrieved by requesting a time or index into the contour.  This can be done manually or using an iterator to
 * iterate over the whole contour.  The iterator returns a pair of time and value information.
 * <p/>
 * Empty entries can be placed in the contour.  When iterating, these entries are skipped.
 */
public class Contour implements Iterable<Pair<Double, Double>> {
  private double x0;        // The starting time
  private double dx;        // The time step
  private double xN;        // The time step
  private int n;            // The number of steps in the contour
  
  private double ceiling;
  
  
  private double[] values;  // The values of the contour

  private int num_empty;
  private boolean[] empty_values;

  /**
   * Constructs a new contour, with a specified starting value, time step, and initial values.
   *
   * @param x0     the starting time
   * @param dx     the time step
   * @param values the initial values
   */
  public Contour(double x0, double dx, double[] values) {
    this.x0 = x0;
    this.dx = dx;
    this.n = values.length;
    this.xN = x0 + dx*this.n;
    this.values = values;
    this.empty_values = new boolean[values.length];
    for (int i = 0; i < empty_values.length; ++i) {
      empty_values[i] = false;
    }
    this.num_empty = 0;
  }

  /**
   * Constructs a new empty contour with n values, a starting time x0 and a time step dx.
   *
   * @param x0 the starting time
   * @param dx the time step
   * @param n  the number of values
   */
  public Contour(double x0, double dx, int n) {
    this.x0 = x0;
    this.dx = dx;
    this.n = n;
    this.xN = x0 + dx*(this.n-1);
    this.values = new double[n];

    this.empty_values = new boolean[n];
    for (int i = 0; i < n; ++i) {
      empty_values[i] = true;
    }
    this.num_empty = n;
  }
  

  /**
   * Retrieves the starting time of the contour.
   *
   * @return the starting time
   */
  public double getStart() {
    return x0;
  }
  
  public double getEnd() {
    return xN;
  }

  /**
   * Retrieves the time step of the contour
   *
   * @return the time step
   */
  public double getStep() {
    return dx;
  }

  /**
   * Retrieves the size of the contour.
   *
   * @return the size of the contour
   */
  public int size() {
    return n;
  }

  /**
   * Calculates the closest index for a given time.
   *
   * @param time the time
   * @return the closest integral index
   */
  public int indexFromTime(double time) {
    return (int) (Math.round((time - x0) / dx));
  }

  /**
   * Calculates the closest index greater than time.
   *
   * @param time the time
   * @return the associated index.
   */
  public int indexFromTimeCeil(double time) {
    return (int) (Math.ceil((time - x0) / dx));
  }

  /**
   * Calculates the closest index lower than time.
   *
   * @param time the time
   * @return the associated index
   */
  public int indexFromTimeFloor(double time) {
    return (int) (Math.floor((time - x0) / dx));
  }

  /**
   * Calculates the time associated with an index.
   * <p/>
   * The associated time is the lower bound of the bin associated with the index.  That is, with a time step of 1s, the
   * time associated with index 0 is 0.0s, not 0.5s.
   *
   * @param index the index.
   * @return the time
   */
  public double timeFromIndex(int index) {
    return x0 + index * dx;
  }

  /**
   * Gets the value associated with a time.
   * <p/>
   * Returns NaN if the time has no associated value.
   *
   * @param time the time
   * @return the value stored at that time.
   */
  public double get(double time) {
    return get(indexFromTime(time));
  }

  /**
   * Gets the value associated with an index.
   * <p/>
   * Returns NaN if the index has no associated value
   *
   * @param index the index
   * @return the value stored at the index
   */
  public double get(int index) {
    if (index < 0 || index >= n || isEmpty(index))
      return Double.NaN;
    return values[index];
  }
  

  /**
   * Sets the value for a given time.
   * <p/>
   * This will overwrite the value currently associated with that time.
   * <p/>
   * Be careful, two distinct times can be mapped to the same index in the contour.  This can lead to accidentally
   * overwriting a contour value if the time steps are not consistent.
   *
   * @param time  the time
   * @param value the value
   */
  public void set(double time, double value) {
    set(indexFromTime(time), value);
  }

  /**
   * Sets the value for a given index.
   * <p/>
   * Will resize the contour array if necessary.
   *
   * @param index the index
   * @param value the value
   */
  public void set(int index, double value) {
    if (index > values.length) {
      // Resize the array
      int origlen = values.length;
      values = Arrays.copyOf(values, index + 1);
      empty_values = Arrays.copyOf(empty_values, index + 1);
      n = values.length;
      for (int i = origlen; i < n; ++i) {
        setEmpty(i);
      }
    }
    values[index] = value;
    if (isEmpty(index)) {
      setEmpty(index, false);
    }
  }

  /**
   * Sets an index to the desired empty state.
   * <p/>
   * Also maintains the num_empty variable.
   *
   * @param index the index to assign
   * @param b     the boolean value to set
   */
  private void setEmpty(int index, boolean b) {
    if (empty_values[index] && !b)
      --num_empty;
    if (!empty_values[index] && b)
      ++num_empty;
    empty_values[index] = b;
  }

  /**
   * Sets the value associated with a time to be empty.
   * <p/>
   * Be careful, two distinct times can be mapped to the same index in the contour.  This can lead to accidentally
   * overwriting a contour value if the time steps are not consistent.
   *
   * @param time the time
   */
  public void setEmpty(double time) {
    setEmpty(indexFromTime(time), true);
  }

  /**
   * Sets the value associated win an index to be empty.
   *
   * @param index the index
   */
  public void setEmpty(int index) {
    setEmpty(index, true);
  }

  /**
   * Retrieces a ContourIterator to traverse the contour.
   *
   * @return the iterator
   */
  public ContourIterator iterator() {
    return new ContourIterator(this);
  }

  /**
   * Determines if an index is empty.
   * <p/>
   * Will return true if i is out of bounds.
   *
   * @param i the index
   * @return true if empty, false otherwise
   */
  public boolean isEmpty(int i) {
    return i < 0 || i >= n || empty_values[i];
  }

  /**
   * Retrieves a time value pair associated with an index.
   * <p/>
   * Returns null if the time is empty.
   *
   * @param i the index
   * @return a pair of time and value
   */
  public Pair<Double, Double> getPair(int i) {
    if (i < 0 || i >= n || isEmpty(i))
      return null;
    return new Pair<Double, Double>(timeFromIndex(i), values[i]);
  }


  /**
   * Returns the number of non-empty values in the Contour.
   *
   * @return the number of elements in the contour.
   */
  public int contentSize() {
    return n - num_empty;
  }
  
  public void setCeiling(double ceiling){
      this.ceiling = ceiling;
  }
  
  public double getCeiling(){
      return this.ceiling;
  }
  
  public double getMaximum(){
      double maximum = values[0];
      for(int i=1;i<values.length;i++){
          if(maximum < values[i]){
              maximum = values[i];
          }
      }
      
      return maximum;
  }
  
  public double getMinimum(){
      double minimum = values[0];
      for(int i=1;i<values.length;i++){
          if(minimum > values[i]){
              minimum = values[i];
          }
      }
      
      return minimum;
  }
  
  public double[] getValues(){
      return values;
  }
  

  
}
