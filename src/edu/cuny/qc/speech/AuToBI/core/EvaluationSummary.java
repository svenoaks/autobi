/*  EvaluationSummary.java

    Copyright (c) 2011 Andrew Rosenberg

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

import edu.cuny.qc.speech.AuToBI.util.AuToBIUtils;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * A class for storing Evaluation data from a set of classification experiments. This class serves as a wrapper around a
 * set of annotated contingency matrices stored as EvaluationResults.
 */
public class EvaluationSummary {
  ArrayList<EvaluationResults> results;

  /**
   * Constructs an empty EvaluationSummary.
   */
  public EvaluationSummary() {
    results = new ArrayList<EvaluationResults>();
  }

  /**
   * Constructs a new EvaluationSummary with a single result in it.
   *
   * @param eval the evaluation result
   */
  public EvaluationSummary(EvaluationResults eval) {
    results = new ArrayList<EvaluationResults>();
    try {
      add(eval);
    } catch (AuToBIException e) {
      AuToBIUtils.error("FATAL ERROR: EvaluationSummary Construction failed.");
    }
  }

  /**
   * Add the results of an evaluation to this collection
   *
   * @param eval the results object to add
   * @throws AuToBIException when the new eval object is inconsistent with existing ones
   */
  public void add(EvaluationResults eval) throws AuToBIException {
    if (results.size() == 0) {
      results.add(eval);
    } else {
      // check that the new results are consistent with the existing results
      if (consistentClasses(eval, results.get(0))) {
        results.add(eval);
      } else {
        throw new AuToBIException(
            "new evaluation results are inconsistent with existing results\nnew number of classes: " +
                eval.getNumClasses() +
                "\nexisting number of classes: " + results.get(0).getNumClasses());
      }
    }
  }

  /**
   * Returns true if the class labels of two EvaluationResults are consistent, false otherwise.
   * <p/>
   * Two arrays of class labels are consistent if they contain the same elements.
   *
   * @param eval_a first EvaluationResult
   * @param eval_b second EvaluationResult
   * @return true if consistent, false otherwise.
   */
  private boolean consistentClasses(EvaluationResults eval_a, EvaluationResults eval_b) {
    if (eval_a.getNumClasses() == eval_b.getNumClasses()) {
      for (String s : eval_a.getClassNames()) {
        if (!Arrays.asList(eval_b.getClassNames()).contains(s)) {
          return false;
        }
      }
      return true;
    } else {
      return false;
    }
  }

  /**
   * Collapse the contingency tables and extract the accuracy over all stored evaluations
   *
   * @return the overall accuracy
   */
  public double getAccuracy() {
    double correct = 0;
    double n = 0;
    for (EvaluationResults eval : results) {
      correct += eval.getNumCorrect();
      n += eval.getNumInstances();
    }

    return correct / n;
  }

  /**
   * Calculate the standard deviation of accuracy from all stored evaluations
   *
   * @return standard deviation of accuracy
   */
  public double getStdevAccuracy() {
    Aggregation agg = new Aggregation();
    for (EvaluationResults eval : results) {
      agg.insert(eval.getPctCorrect());
    }
    return agg.getStdev();
  }


  /**
   * Collapse the contingency tables and extract the f-measure for a class over all evaluations
   *
   * @param class_name The f-measure class
   * @return the overall f-measure
   * @throws AuToBIException if the class does not exist
   */
  public double getFMeasure(String class_name) throws AuToBIException {
    double correct = 0;
    double classified_as = 0;
    double class_instances = 0;

    for (EvaluationResults eval : results) {
      correct += eval.getNumCorrectClassInstances(class_name);
      classified_as += eval.getNumClassifiedAs(class_name);
      class_instances += eval.getNumClassInstances(class_name);
    }

    double precision = correct / classified_as;
    double recall = correct / class_instances;
    return (2 * precision * recall) / (precision + recall);
  }

  /**
   * Return the precision of the class
   *
   * @param class_name the class to query
   * @return preision of class
   * @throws AuToBIException if the class does not exist
   */
  public double getPrecision(String class_name) throws AuToBIException {
    double correct = 0;
    double classified_as = 0;

    for (EvaluationResults eval : results) {
      correct += eval.getNumCorrectClassInstances(class_name);
      classified_as += eval.getNumClassifiedAs(class_name);
    }

    return correct / classified_as;
  }

  /**
   * Return the recall of the class
   *
   * @param class_name the class to query
   * @return the recall of the class
   * @throws AuToBIException if the class does not exist
   */
  public double getRecall(String class_name) throws AuToBIException {
    double correct = 0;
    double class_instances = 0;

    for (EvaluationResults eval : results) {
      correct += eval.getNumCorrectClassInstances(class_name);
      class_instances += eval.getNumClassInstances(class_name);
    }
    return correct / class_instances;
  }

  /**
   * Generate a string representation of the collapsed contingency matrix
   *
   * @return string representation of the contingency matrix
   */
  private String printContingencyMatrix() {
    String s = "";
    if (results.size() > 0) {
      s += String.format(" %1$10s", "");
      for (String c : results.get(0).getClassNames()) {
        s += String.format("%1$6s ", c.substring(0, Math.min(5, c.length())));
      }
      s += "\n";
      for (String c1 : results.get(0).getClassNames()) {
        s += "-" + String.format("%1$10s ", c1);
        for (String c2 : results.get(0).getClassNames()) {
          int count = 0;
          for (EvaluationResults result : results) {
            try {
              count += result.getInstances(c1, c2);
            } catch (AuToBIException e) {
              AuToBIUtils.error("FATAL EROR: Inconsistent EvaluationSummary.");
              return s;
            }
          }
          s += String.format("%1$6d ", count);
        }
        s += "\n";
      }
    }
    return s;
  }

  /**
   * Construct a result string with no class value to report f-measure for
   *
   * @return a description of the aggregated evaluation results
   */
  public String toString() {
    return toString(null);
  }

  /**
   * Construct a result string with f-measure based on the given class values.
   * <p/>
   * Class values should be a comma separated string of class names.
   *
   * @param class_values the class values to report on
   * @return result string
   */
  public String toString(String class_values) {
    String s = "Accuracy -- Mean: " + getAccuracy() + "\n";
    s += "         -- Stdev: " + getStdevAccuracy() + "\n";
    s += "         -- Sterr: " + getStdevAccuracy() / Math.sqrt(results.size()) + "\n";
    s += "         -- Conf : " + getStdevAccuracy() * 1.64 / Math.sqrt(results.size()) + "\n";

    if (class_values != null)
      for (String class_value : class_values.split(",")) {
        if (class_value != null && class_value.length() > 0) {
          try {
            s += "FMeasure (" + class_value + ") -- Mean: " + getFMeasure(class_value) + "\n";
            s += "Precision -- " + getPrecision(class_value) + "\n";
            s += "Recall -- " + getRecall(class_value) + "\n";
          } catch (AuToBIException e) {
            s += "FMeasure -- " + e.getMessage();
          }
        }
      }
    s += "Contingency Matrix\n";
    s += printContingencyMatrix() + "\n";

    if (results.size() > 0) {
      for (String class_name : results.get(0).getClassNames()) {
        try {
          s += class_name + " - FMeasure: " + getFMeasure(class_name) + "\n";
        } catch (AuToBIException e) {
          e.printStackTrace();
        }
      }
    }
    s += "Mutual Information: " + getMutualInformation() + "\n";
    return s;
  }

  /**
   * Calculates the mean Mutual Information between class values and hypotheses across all results.
   *
   * @return the mean MI.
   */
  public Double getMutualInformation() {
    Double mi = 0.0;
    for (EvaluationResults result : results) {
      mi += result.getMutualInformation();
    }
    mi /= results.size();
    return mi;
  }

  /**
   * Get the number of points with true value equal to class name across all contained evaluation results
   *
   * @param class_name the class to aggregate
   * @return the number of points with class equal to class_name
   * @throws AuToBIException if the class is not found in results.
   */
  public Integer getNumClassInstances(String class_name) throws AuToBIException {
    Integer n = 0;

    for (EvaluationResults r : results) {
      n += r.getNumClassInstances(class_name);
    }
    return n;
  }
}
