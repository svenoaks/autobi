/*  AuToBIParameters.java

    Copyright (c) 2009-2010 Andrew Rosenberg

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

import edu.cuny.qc.speech.AuToBI.core.AuToBIException;
import edu.cuny.qc.speech.AuToBI.util.AuToBIUtils;

import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.util.HashMap;

/**
 * AuToBIParameters stores and provides access to command line arguments sent to AuToBI.
 */
public class AuToBIParameters {
  // A map of parameters and values
  private HashMap<String, String> parameters = new HashMap<String, String>();

  /**
   * Parses command line arguments into parameters.
   * <p/>
   * parameters are of the form -parameter_name=parameter_value
   * <p/>
   * or for boolean parameters -boolean_parameter_name
   *
   * @param args command line parameters
   */
  public void readParameters(String[] args) {
    Pattern boolean_pattern = Pattern.compile("-(\\w*)");
    Pattern parameter_pattern = Pattern.compile("-(\\w*)=(.*)");
    for (String param : args) {
      Matcher boolean_matcher = boolean_pattern.matcher(param);
      Matcher parameter_matcher = parameter_pattern.matcher(param);
      if (boolean_matcher.matches()) {
        parameters.put(boolean_matcher.group(1), "true");
      }
      if (parameter_matcher.matches()) {
        parameters.put(parameter_matcher.group(1), parameter_matcher.group(2));
      }
    }
  }

  /**
   * Sets a parameter value manually.
   *
   * @param parameter the parameter name
   * @param value     the parameter value
   */
  public void setParameter(String parameter, String value) {
    parameters.put(parameter, value);
  }

  /**
   * Retrieves a value of a parameter
   *
   * @param parameter_name the parameter name
   * @return the parameter value
   * @throws AuToBIException if the parameter was not set
   */
  public String getParameter(String parameter_name) throws AuToBIException {
    if (parameters.containsKey(parameter_name)) return parameters.get(parameter_name);
    else {
      throw new AuToBIException("No parameter, " + parameter_name + ", defined.");
    }
  }

  /**
   * Retrieves a value of a parameter.
   *
   * @param parameter_name the parameter naem
   * @return the parameter value or null if the parameter was not set.
   */
  public String getOptionalParameter(String parameter_name) {
    if (parameters.containsKey(parameter_name)) {
      return parameters.get(parameter_name);
    } else {
      return null;
    }
  }

  /**
   * Retrieves a value of a parameter or a default value if the parameter hasn't been set.
   *
   * @param parameter_name the parameter naem
   * @param default_value  the default value.
   * @return the parameter value or null if the parameter was not set.
   */
  public String getOptionalParameter(String parameter_name, String default_value) {
    if (hasParameter(parameter_name)) {
      return parameters.get(parameter_name);
    }
    return default_value;
  }

  /**
   * Determines if the parameter has been set.
   *
   * @param parameter_name the parameter name
   * @return true if the parameter has been set, false otherwise.
   */
  public Boolean hasParameter(String parameter_name) {
    return parameters.containsKey(parameter_name);
  }

  /**
   * Get the value of a boolean parameter
   * <p/>
   * This will raise a warning if a parameter that has not been set as boolean is interpreted as such.
   *
   * @param parameter_name the parameter name
   * @param default_value  a default value if the parameter was not set.
   * @return the value of the boolean parameter
   */
  public Boolean booleanParameter(String parameter_name, Boolean default_value) {
    try {
      if (getParameter(parameter_name).equalsIgnoreCase("true")) return true;
      if (getParameter(parameter_name).equalsIgnoreCase("false")) return false;
    } catch (AuToBIException e) {
      return default_value;
    }
    return default_value;
  }
}
