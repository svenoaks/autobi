/*  IntonationalPhraseBoundaryDetectionFeatureSet.java

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
package edu.cuny.qc.speech.AuToBI.featureset;

import edu.cuny.qc.speech.AuToBI.core.FeatureSet;

/**
 * IntonationalPhraseBoundaryDetectionFeatureSet describes the required features and class attribute for the
 * Intonational Phrase Boundary detection task.
 */
public class IntonationalPhraseBoundaryDetectionFeatureSet extends FeatureSet {

    /**
     * Constructs a new IntonationalPhraseBoundaryDetectionFeatureSet
     */
    public IntonationalPhraseBoundaryDetectionFeatureSet() {
        super();

        insertRequiredFeature("nominal_precedesSilence");

        insertRequiredFeature("duration__duration");
        insertRequiredFeature("duration__follPause");
        insertRequiredFeature("diff_duration__duration");
        
        for (String acoustic : new String[]{"f0", "I"}) {
        
            for (String diff : new String[]{"", "diff_"}) {
                for (String norm : new String[]{"", "norm_"}) {
                    for (String slope : new String[]{"", "delta_"}) {
                        for (String agg : new String[]{"max", "mean", "stdev", "zMax"}) {
                            insertRequiredFeature(diff + slope + norm + acoustic + "__" + agg);
                        }
                    }
                }
            }
        }

        class_attribute = "nominal_IntonationalPhraseBoundary";
    }
}
