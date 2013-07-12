/*  AuToBI.java

    Copyright (c) 2009-2012 Andrew Rosenberg

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

import edu.cuny.qc.speech.AuToBI.classifier.AuToBIClassifier;
import edu.cuny.qc.speech.AuToBI.core.*;
import edu.cuny.qc.speech.AuToBI.featureextractor.*;
import edu.cuny.qc.speech.AuToBI.featureset.*;
import edu.cuny.qc.speech.AuToBI.io.*;
import edu.cuny.qc.speech.AuToBI.util.AuToBIUtils;
import edu.cuny.qc.speech.AuToBI.util.ClassifierUtils;
import edu.cuny.qc.speech.AuToBI.util.WordReaderUtils;
import org.apache.log4j.PropertyConfigurator;
import org.apache.log4j.BasicConfigurator;

import javax.sound.sampled.UnsupportedAudioFileException;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

import static java.util.concurrent.Executors.newFixedThreadPool;

/**
 * This is the main class for the AuToBI system.
 * <p/>
 * AuToBI generates hypothesized ToBI tones based on manual segmentation of words and an audio file.
 * <p/>
 * The prediction is divided into six tasks. 1) detection of pitch accents 2) classification of pitch accent types 3)
 * detection of intonational phrase boundaries 4) classification of intonational phrase ending tones (phrase accent and
 * boundary tone pairs) 5) detection of intermediate phrase boundaries 6) classification of intermediate phrase ending
 * tones (phrase accents)
 * <p/>
 * Each of these tasks require distinct features to be extracted from the words that are being analyzed. To perform the
 * feature extraction, each task has an associated FeatureSet which describes the features required for classification.
 * AuToBI maintains a registry of FeatureExtractors describing the Features that they will calculate when processing a
 * word.  When extracting features for a FeatureSet, AuToBI only calls those FeatureExtractors necessary to generate the
 * required features.
 * <p/>
 * This class manages command line parameters, the execution of feature extractors and the generation of hypothesized
 * ToBI tones.
 *
 * @see FeatureSet
 * @see FeatureExtractor
 */
public class AuToBIDisfluencies {
  private AuToBIParameters params;  // Command line parameters

  // A map from feature names to the FeatureExtractors that generate them
  private Map<String, FeatureExtractor> feature_registry;

  // A set of FeatureExtractors that have already executed
  protected Set<FeatureExtractor> executed_feature_extractors;

  // A map from input filenames to serialized speaker normalization parameter files.
  private Map<String, String> speaker_norm_file_mapping;

  // A map of the number of times each feature is needed.
  private HashMap<String, Integer> reference_count;

  // A set of features to delete on the next garbage collection call
  private Set<String> dead_features;

  // A list of AuToBITasks to be executed.
  protected HashMap<String, AuToBITask> tasks;

  /**
   * Constructs a new AuToBI object.
   */
  public AuToBIDisfluencies() {
    params = new AuToBIParameters();
    feature_registry = new HashMap<String, FeatureExtractor>();
    executed_feature_extractors = new HashSet<FeatureExtractor>();
    speaker_norm_file_mapping = new HashMap<String, String>();
    tasks = new HashMap<String, AuToBITask>();
  }

  /**
   * Parses command line parameters and sets up log4j logging.
   *
   * @param args command line arguments.
   */
  public void init(String[] args) {
    params.readParameters(args);
    if (params.hasParameter("log4j_config_file")) {
      try {
        PropertyConfigurator.configure(params.getParameter("log4j_config_file"));
      } catch (Exception e) {
        BasicConfigurator.configure();
        AuToBIUtils.log("Couldn't read -log4j_config_file. BasicConfigurator logging to console");
      }
    } else {
      BasicConfigurator.configure();
    }
  }

  /**
   * Gets the requested parameter.
   *
   * @param parameter the parameter name
   * @return the parameter value
   * @throws AuToBIException if the parameter was not set
   */
  public String getParameter(String parameter) throws AuToBIException {
    return params.getParameter(parameter);
  }

  /**
   * Gets the requested parameter.
   *
   * @param parameter the parameter name
   * @return the parameter value or null if the parameter was not set
   */
  public String getOptionalParameter(String parameter) {
    return params.getOptionalParameter(parameter);
  }

  /**
   * Gets the requested boolean parameter with default value if the parameter is not set.
   *
   * @param parameter     the parameter name
   * @param default_value the default boolean value if not explicitly set.
   * @return the parameter value or null if the parameter was not set
   */
  public Boolean getBooleanParameter(String parameter, boolean default_value) {
    return params.booleanParameter(parameter, default_value);
  }

  /**
   * Gets the requested parameter with a default value if the parameter has not been set.
   *
   * @param parameter     the parameter name
   * @param default_value a default value
   * @return The parameter value or default_value if the parameter has not been set.
   */
  public String getOptionalParameter(String parameter, String default_value) {
    return params.getOptionalParameter(parameter, default_value);
  }

  /**
   * Determines if a parameter has been set or not.
   *
   * @param parameter the parameter name
   * @return true if the parameter has been set, false otherwise
   */
  public boolean hasParameter(String parameter) {
    return params.hasParameter(parameter);
  }

  /**
   * Evaluates the performance of a particular task.
   *
   * @param task a task identifier
   * @param fs   the feature set
   * @return a string representation of the evaluation
   * @throws AuToBIException if there is no task corresponding to the identifier
   */
  protected String evaluateTaskPerformance(String task, FeatureSet fs) throws AuToBIException {
    if (tasks.containsKey(task)) {
      AuToBITask autobi_task = tasks.get(task);
      return ClassifierUtils.evaluateClassification(autobi_task.getHypFeature(), autobi_task.getTrueFeature(), fs);
    } else throw new AuToBIException(
        "Task, " + task + ", is unavailable. Either no classifier has been set, or there is some other problem");
  }
  
  protected String evaluateDetectionPerformance(String task, FeatureSet fs) throws AuToBIException {
      
      
      String result = "";
      if (tasks.containsKey(task)) {
          AuToBITask autobi_task = tasks.get(task);
          result = ClassifierUtils.evaluateClassification(autobi_task.getHypFeature(), autobi_task.getTrueFeature(), fs);
          
      }
      return result;
  }

  /**
   * Generates predictions corresponding to a particular task
   *
   * @param task A task identifier
   * @param fs   The feature set
   * @throws AuToBIException if there is no task associated with the identifier.
   */
  public void generatePredictions(String task, FeatureSet fs) throws AuToBIException {
    if (tasks.containsKey(task)) {
      AuToBITask autobi_task = tasks.get(task);
      ClassifierUtils
          .generatePredictions(autobi_task.getClassifier(),
              autobi_task.getHypFeature(),
              autobi_task.getDefaultValue(),
              fs);
    } else throw new AuToBIException(
        "Task, " + task + ", is unavailable. Either no classifier has been set, or there is some other problem");
  }

  /**
   * Generates predictions with confidence scores corresponding to a particular task
   *
   * @param task A task identifier
   * @param fs   The feature set
   * @throws AuToBIException if there is no task associated with the identifier.
   */
  public void generatePredictionsWithConfidenceScores(String task, FeatureSet fs) throws AuToBIException {
    if (tasks.containsKey(task)) {
      AuToBITask autobi_task = tasks.get(task);
      ClassifierUtils
          .generatePredictionsWithConfidenceScores(autobi_task.getClassifier(),
              autobi_task.getHypFeature(),
              autobi_task.getConfFeature(),
              autobi_task.getDefaultValue(),
              fs);
    } else throw new AuToBIException(
        "Task, " + task + ", is unavailable. Either no classifier has been set, or there is some other problem");
  }

  /**
   * Retrieves an empty feature set for the given task.
   *
   * @param task a task identifier.
   * @return a corresponding FeatureSet object
   * @throws AuToBIException If there is no FeatureSet defined for the task identifier
   */
  public FeatureSet getTaskFeatureSet(String task) throws AuToBIException {
    if (tasks.containsKey(task)) {
      AuToBITask autobi_task = tasks.get(task);
      if (autobi_task.getFeatureSet() == null) {
        throw new AuToBIException("Task, " + task + ", does not have an associated feature set");
      }
      return autobi_task.getFeatureSet().newInstance();
    } else throw new AuToBIException(
        "Task, " + task + ", is unavailable. Either no classifier has been set, or there is some other problem");
  }

  /**
   * Retrieves a default hypothesized feature name set for the given task.
   *
   * @param task a task identifier.
   * @return a string for the hypothesized name
   * @throws AuToBIException If there is no FeatureSet defined for the task identifier
   */
  public String getHypotheizedFeature(String task) throws AuToBIException {
    if (tasks.containsKey(task)) {
      AuToBITask autobi_task = tasks.get(task);
      return autobi_task.getHypFeature();
    } else throw new AuToBIException(
        "Task, " + task + ", is unavailable. Either no classifier has been set, or there is some other problem");
  }

  /**
   * Retrieves a default feature distribution name set for the given task.
   *
   * @param task a task identifier.
   * @return a string for the hypothesized name
   * @throws AuToBIException If there is no FeatureSet defined for the task identifier
   */
  public String getDistributionFeature(String task) throws AuToBIException {
    if (tasks.containsKey(task)) {
      AuToBITask autobi_task = tasks.get(task);
      return autobi_task.getDistFeature();
    } else throw new AuToBIException(
        "Task, " + task + ", is unavailable. Either no classifier has been set, or there is some other problem");
  }

  /**
   * Retrieves a previously loaded AuToBIClassifier for the given task.
   *
   * @param task a task identifier.
   * @return a corresponding FeatureSet object
   * @throws AuToBIException If there is no FeatureSet defined for the task identifier
   */
  public AuToBIClassifier getTaskClassifier(String task) throws AuToBIException {
    if (tasks.containsKey(task)) {
      AuToBITask autobi_task = tasks.get(task);
      return autobi_task.getClassifier();
    } else throw new AuToBIException(
        "Task, " + task + ", is unavailable. Either no classifier has been set, or there is some other problem");
  }

  
  /**
   * Initializes the Autobi task list.  This is driven by loading classifiers from serialized objects.
   * <p/>
   * When a classifier is correctly loaded a corresponding task object is created to handle the appropriate bookkeeping
   * <p/>
   * Only those classifiers which have been specified using the following command line parameters are loaded:
   * -pitch_accent_detector -pitch_accent_classifier -IP_detector -ip_detector -phrase_accent_classifier
   * -boundary_tone_classifier
   */
  public void initializeAuToBITasks() {
    try {
      String pad_filename = getParameter("pitch_accent_detector");
      AuToBIClassifier pitch_accent_detector = ClassifierUtils.readAuToBIClassifier(pad_filename);
      AuToBITask task = new AuToBITask();
      task.setClassifier(pitch_accent_detector);
      task.setTrueFeature("nominal_PitchAccent");
      task.setHypFeature("hyp_pitch_accent");
      task.setConfFeature("hyp_pitch_accent_conf");
      task.setDistFeature("hyp_pitch_accent_dist");
      task.setFeatureSet(new PitchAccentDetectionFeatureSet());
      tasks.put("pitch_accent_detection", task);
    } catch (AuToBIException ignored) {
    }

    try {
      String pac_filename = getParameter("pitch_accent_classifier");
      AuToBIClassifier pitch_accent_classifier = ClassifierUtils.readAuToBIClassifier(pac_filename);
      AuToBITask task = new AuToBITask();
      task.setClassifier(pitch_accent_classifier);
      task.setTrueFeature("nominal_PitchAccentType");
      task.setHypFeature("hyp_pitch_accent_type");
      task.setConfFeature("hyp_pitch_accent_type_conf");
      task.setDistFeature("hyp_pitch_accent_type_dist");
      task.setFeatureSet(new PitchAccentClassificationFeatureSet());
      tasks.put("pitch_accent_classification", task);
    } catch (AuToBIException ignored) {
    }

    try {
      String intonational_phrase_detector_filename = getParameter("IP_detector");
      AuToBIClassifier intonational_phrase_boundary_detector =
          ClassifierUtils.readAuToBIClassifier(intonational_phrase_detector_filename);
      AuToBITask task = new AuToBITask();
      task.setClassifier(intonational_phrase_boundary_detector);
      task.setTrueFeature("nominal_IntonationalPhraseBoundary");
      task.setHypFeature("hyp_IP_location");
      task.setConfFeature("hyp_IP_location_conf");
      task.setDistFeature("hyp_IP_location_dist");
      task.setFeatureSet(new IntonationalPhraseBoundaryDetectionFeatureSet());
      tasks.put("intonational_phrase_boundary_detection", task);
    } catch (AuToBIException ignored) {
    }

    try {
      String intermediate_phrase_detector_filename = getParameter("ip_detector");
      AuToBIClassifier intermediate_phrase_boundary_detector =
          ClassifierUtils.readAuToBIClassifier(intermediate_phrase_detector_filename);
      AuToBITask task = new AuToBITask();
      task.setClassifier(intermediate_phrase_boundary_detector);
      task.setTrueFeature("nominal_IntermediatePhraseBoundary");
      task.setHypFeature("hyp_ip_location");
      task.setConfFeature("hyp_ip_location_conf");
      task.setDistFeature("hyp_ip_location_dist");
      task.setFeatureSet(new IntermediatePhraseBoundaryDetectionFeatureSet());
      tasks.put("intermediate_phrase_boundary_detection", task);
    } catch (AuToBIException ignored) {
    }

    try {
      String phrase_accent_classifier_name = getParameter("phrase_accent_classifier");
      AuToBIClassifier phrase_accent_classifier =
          ClassifierUtils.readAuToBIClassifier(phrase_accent_classifier_name);
      AuToBITask task = new AuToBITask();
      task.setClassifier(phrase_accent_classifier);
      task.setTrueFeature("nominal_PhraseAccent");
      task.setHypFeature("hyp_phrase_accent");
      task.setConfFeature("hyp_phrase_accent_conf");
      task.setDistFeature("hyp_phrase_accent_dist");
      task.setFeatureSet(new PhraseAccentClassificationFeatureSet());
      tasks.put("phrase_accent_classification", task);
    } catch (AuToBIException ignored) {
    }

    try {
      String boundary_tone_classifier_name = getParameter("boundary_tone_classifier");
      AuToBIClassifier boundary_tone_classifier =
          ClassifierUtils.readAuToBIClassifier(boundary_tone_classifier_name);
      AuToBITask task = new AuToBITask();
      task.setClassifier(boundary_tone_classifier);
      task.setTrueFeature("nominal_PhraseAccentBoundaryTone");
      task.setHypFeature("hyp_pabt");
      task.setConfFeature("hyp_pabt_conf");
      task.setDistFeature("hyp_pabt_dist");
      task.setFeatureSet(new PhraseAccentBoundaryToneClassificationFeatureSet());
      tasks.put("boundary_tone_classification", task);
    } catch (AuToBIException ignored) {
    }
  }

  /**
   * Retrieves a list of classification task identifiers corresponding to the tasks to be performed.
   * <p/>
   * Only those tasks corresponding to loaded classifiers are executed.
   *
   * @return a list of task identifiers.
   */
  public Set<String> getClassificationTasks() {
    return tasks.keySet();
  }

  /**
   * Writes a TextGrid file containing words and hypothesized ToBI labels.
   *
   * @param words    The words
   * @param out_file The destination file
   * @throws IOException If there is a problem writing to the destination file.
   */
  public void writeTextGrid(List<Word> words, String out_file) throws IOException {
    String text_grid = generateTextGridString(words);

    AuToBIFileWriter writer = new AuToBIFileWriter(out_file);
    writer.write(text_grid);
    writer.close();
  }

  /**
   * Generates a TextGrid representation of the words and hypothesized ToBI labels.
   *
   * @param words the words to output
   * @return a string representing the textgrid contents of the words.
   */
  public String generateTextGridString(List<Word> words) {
    String text_grid = "File type = \"ooTextFile\"\n";
    text_grid += "Object class = \"TextGrid\"\n";
    text_grid += "xmin = 0\n";
    text_grid += "xmax = " + words.get(words.size() - 1).getEnd() + "\n";
    text_grid += "tiers? <exists>\n";
    text_grid += "size = 3\n";
    text_grid += "item []:\n";
    text_grid += "item [1]:\n";
    text_grid += "class = \"IntervalTier\"\n";
    text_grid += "name = \"words\"\n";
    text_grid += "xmin = 0\n";
    text_grid += "xmax = " + words.get(words.size() - 1).getEnd() + "\n";
    text_grid += "intervals: size = " + words.size() + "\n";
    for (int i = 0; i < words.size(); ++i) {
      Word w = words.get(i);
      text_grid += "intervals [" + (i + 1) + "]:\n";
      text_grid += "xmin = " + w.getStart() + "\n";
      text_grid += "xmax = " + w.getEnd() + "\n";
      text_grid += "text = \"" + w.getLabel() + "\"\n";
    }


    text_grid += "item [2]:\n";
    text_grid += "class = \"IntervalTier\"\n";
    text_grid += "name = \"pitch_accent_hypothesis\"\n";
    text_grid += "xmin = 0\n";
    text_grid += "xmax = " + words.get(words.size() - 1).getEnd() + "\n";
    text_grid += "intervals: size = " + words.size() + "\n";
    for (int i = 0; i < words.size(); ++i) {
      Word w = words.get(i);
      String text = "";
      if (getBooleanParameter("distributions", false)) {
        String det_dist_feature = tasks.get("pitch_accent_detection").getDistFeature();
        String class_dist_feature = tasks.get("pitch_accent_classification").getDistFeature();
        if (w.hasAttribute(det_dist_feature)) {
          text = w.getAttribute(det_dist_feature).toString();
        }
        if (w.hasAttribute(class_dist_feature)) {
          text += w.getAttribute(class_dist_feature).toString();
        }
      } else {
        if (tasks.containsKey("pitch_accent_detection") && w.hasAttribute(tasks.get("pitch_accent_detection").getHypFeature())) {
          text = w.getAttribute(tasks.get("pitch_accent_detection").getHypFeature()).toString();
        }
      }

      text_grid += "intervals [" + (i + 1) + "]:\n";
      text_grid += "xmin = " + w.getStart() + "\n";
      text_grid += "xmax = " + w.getEnd() + "\n";
      text_grid += "text = \"" + text + "\"\n";
    }

    text_grid += "item [3]:\n";
    text_grid += "class = \"IntervalTier\"\n";
    text_grid += "name = \"phrase_hypothesis\"\n";
    text_grid += "xmin = 0\n";
    text_grid += "xmax = " + words.get(words.size() - 1).getEnd() + "\n";
    text_grid += "intervals: size = " + words.size() + "\n";
    for (int i = 0; i < words.size(); ++i) {
      Word w = words.get(i);

      String text = "";
      if (getBooleanParameter("distributions", false)) {
        if (w.hasAttribute(tasks.get("intonational_phrase_boundary_detection").getDistFeature())) {
          text = w.getAttribute(tasks.get("intonational_phrase_boundary_detection").getDistFeature()).toString();
        }
        if (w.hasAttribute(tasks.get("intermediate_phrase_boundary_detection").getDistFeature())) {
          text = w.getAttribute(tasks.get("intermediate_phrase_boundary_detection").getDistFeature()).toString();
        }
        if (w.hasAttribute(tasks.get("boundary_tone_classification").getDistFeature())) {
          text += w.getAttribute(tasks.get("boundary_tone_classification").getDistFeature()).toString();
        }
        if (w.hasAttribute(tasks.get("phrase_accent_classification").getDistFeature())) {
          text += w.getAttribute(tasks.get("phrase_accent_classification").getDistFeature()).toString();
        }
      } else {
        if (w.hasAttribute("hyp_phrase_boundary")) {
          text = w.getAttribute("hyp_phrase_boundary").toString();
        }
      }

      text_grid += "intervals [" + (i + 1) + "]:\n";
      text_grid += "xmin = " + w.getStart() + "\n";
      text_grid += "xmax = " + w.getEnd() + "\n";
      text_grid += "text = \"" + text + "\"\n";
    }
    return text_grid;
  }

  public static void main(String[] args) {
    AuToBIDisfluencies autobi = new AuToBIDisfluencies();
    autobi.init(args);

    autobi.run();
  }

  public void run() {
    try {


      String filename = getOptionalParameter("input_file");
      AuToBIWordReader word_reader;
      FormattedFile file;

      if (hasParameter("input_file")) {
        // Let the FormattedFile constructor determine the file based on the extension or other file name conventions
        file = new FormattedFile(getOptionalParameter("input_file"));
        word_reader = WordReaderUtils.getAppropriateReader(file, getParameters());
      } else {
        AuToBIUtils.info(
            "No -input_file filename specified.");
        return;
      }
      AuToBIUtils.log("Reading words from: " + filename);

      if (word_reader == null) {
        AuToBIUtils.error("Unable to create wordreader for file: " + filename + "\n\tCheck the file extension.");
        return;
      }

      if (hasParameter("silence_regex")) {
        word_reader.setSilenceRegex(getParameter("silence_regex"));
      }

      
      List<Word> words = word_reader.readWords();
      double lastTime = 0.0;
      String lastBreak = null;
      int num3 = 0,sil3=0, sil=0;
      int num4 = 0,sil4=0,dis=0;
      for(Word w:words){
          //There is a silence
          if(w.getStart()!=lastTime){
              
              if(w.getBreakBefore()!=null && w.getBreakBefore().equals("3"))
                sil3++;
              else if(w.getBreakBefore()!=null && w.getBreakBefore().equals("4"))
                sil4++;
              else dis++;
              sil++;
          }
          lastTime = w.getEnd();
          lastBreak = w.getBreakAfter();
          
          if(lastBreak.equals("3"))
              num3++;
          if(lastBreak.equals("4"))
              num4++;
          
      }

      if (hasParameter("out_file")) {
        AuToBIFileWriter writer = new AuToBIFileWriter(getOptionalParameter("out_file"));
        writer.write("Number of Intermediate Boundaries: "+num3+"\n");
        writer.write("Number of Intonational Boundaries: "+num4+"\n");
        writer.write("Number of silences: "+sil+"\n");
        writer.write("Number of intermediate silences: "+sil3+" "+(sil3*1.0/sil)+"\n");
        writer.write("Number of intonational silences: "+sil4+" "+(sil4*1.0/sil)+"\n");
        writer.write("Number of disfluencies: "+dis+" "+(dis*1.0/sil)+"\n");
        writer.close();
      }
    } catch (AuToBIException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    } 
  }

  /**
   * Registers a null feature extractor with the registry.
   * <p/>
   * This implies that the feature will be manually set by the user outside the typical feature extraction process.
   * <p/>
   * This is used to satisfy feature requirements for feature extractors without writing a feature extractor for the
   * requirement.  This is often used for features that are assigned by a file reader.
   *
   * @param s the feature name
   */

  public void registerNullFeatureExtractor(String s) {
    this.feature_registry.put(s, null);
  }

  /**
   * Gets command line parameters.
   *
   * @return an AuToBIParameters object
   */
  public AuToBIParameters getParameters() {
    return params;
  }

  /**
   * Sets AuToBIParameters.
   *
   * @param params the parameters.
   */
  public void setParameters(AuToBIParameters params) {
    this.params = params;
  }
}
