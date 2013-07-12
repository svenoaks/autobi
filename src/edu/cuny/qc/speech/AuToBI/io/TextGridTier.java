/*  TextGridTier.java

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
package edu.cuny.qc.speech.AuToBI.io;

import edu.cuny.qc.speech.AuToBI.core.Region;
import edu.cuny.qc.speech.AuToBI.util.AuToBIReaderUtils;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;


/**
 * TextGridTier is a Tier object tailored to the format of TextGrid tiers.
 *
 * @see Tier
 */
public class TextGridTier extends Tier {

    double lastTime = 0.0;
  /**
   * Reads the information from the reader into the Tier regions, and sets the name of the tier.
   * <p/>
   * Returns true if there are additional items in the reader.
   *
   * @param reader the reader to read from
   * @return true if there are more items.
   * @throws IOException                  if there is an input output problem
   * @throws TextGridSyntaxErrorException If there is a formatting problem
   */
  public boolean readTier(AuToBIFileReader reader) throws TextGridSyntaxErrorException, IOException {
    String line;
    while ((line = reader.readLine()) != null) {
      line = AuToBIReaderUtils.removeTabsAndTrim(line);
      if (line.matches("item\\s\\[[\\d]+\\]:")) {
        return true;
      }
      if (line.contains("class = \"IntervalTier\"")) {
        is_point_tier = false;
      }
      if (line.contains("class = \"TextTier\"")) {
        is_point_tier = true;
      } else if (line.contains("name = ")) {
        name = line.replaceFirst("name = \"(.*?)\"", "$1");
      } else if (line.matches("points\\s\\[[\\d]+\\]:") || line.matches("accents\\s\\[[\\d]+\\]:") || line.matches("tones\\s\\[[\\d]+\\]:")) {
        if (!is_point_tier) {
          throw new TextGridSyntaxErrorException("Found point label in interval tier.");
        }
        addPoint(reader);
      } else if (line.matches("intervals\\s\\[[\\d]+\\]:")) {
        if (is_point_tier) {
          throw new TextGridSyntaxErrorException("Found interval label in point tier.");
        }
        addInterval(reader);
      } else if (line.matches("words\\s\\[[\\d]+\\]:")) {
        if (!is_point_tier) {
          throw new TextGridSyntaxErrorException("Found interval label in point tier.");
        }
        addDirndlInterval(reader);
      }
      
    }
    return false;
  }

  /**
   * Reads a point region from the AuToBIFileReader.
   *
   * @param reader the active AuToBIFileReader
   * @throws IOException                  If there is a problem reading the file
   * @throws TextGridSyntaxErrorException If there is a formatting problem
   */
  protected void addPoint(AuToBIFileReader reader)
      throws IOException, TextGridSyntaxErrorException {
    
     
    
    String time = AuToBIReaderUtils.removeTabsAndTrim(reader.readLine());
    if (time == null || !(time.contains("time =") || time.contains("number =")))
      throw new TextGridSyntaxErrorException("missing point at line: " + reader.getLineNumber());

    time = time.replace("time = ", "").trim();
    time = time.replace("number = ", "").trim();
    
    Region region = new Region(Double.valueOf(time));
    region.setFile(reader.getFilename());

    String mark = AuToBIReaderUtils.removeTabsAndTrim(reader.readLine());
    if (mark == null || !mark.contains("mark ="))
      throw new TextGridSyntaxErrorException("missing mark at line: " + reader.readLine());

    mark = mark.replaceAll("mark = \"(.*?)\"", "$1").trim();
    if(!mark.equals("-")){
        region.setLabel(mark);
        regions.add(region);
    }
  }

  /**
   * Reads an interval region from the AuToBIFileReader.
   *
   * @param reader the active AuToBIFileReader
   * @throws IOException                  If there is a problem reading the file
   * @throws TextGridSyntaxErrorException If there is a formatting problem
   */
  protected void addInterval(AuToBIFileReader reader)
      throws IOException, TextGridSyntaxErrorException {
    String xmin = AuToBIReaderUtils.removeTabsAndTrim(reader.readLine());
    if (xmin == null || !xmin.contains("xmin ="))
      throw new TextGridSyntaxErrorException("missing xmin at line: " + reader.getLineNumber());

    xmin = xmin.replace("xmin = ", "").trim();

    String xmax = AuToBIReaderUtils.removeTabsAndTrim(reader.readLine());
    if (xmax == null || !xmax.contains("xmax ="))
      throw new TextGridSyntaxErrorException("missing xmax at line: " + reader.getLineNumber());

    xmax = xmax.replace("xmax = ", "").trim();

    Region region = new Region(Double.valueOf(xmin), Double.valueOf(xmax));
    region.setFile(reader.getFilename());

    String text = AuToBIReaderUtils.removeTabsAndTrim(reader.readLine());
    if (text == null || !text.contains("text ="))
      throw new TextGridSyntaxErrorException("missing mark at line: " + reader.getLineNumber());

    text = text.replaceAll("text = \"(.*?)\"", "$1").trim();
    
    if(text.equals("@BG"))
        text="sil";
    region.setLabel(text);

    regions.add(region);
  }
  
  protected void addDirndlInterval(AuToBIFileReader reader)
      throws IOException, TextGridSyntaxErrorException {
    String time = AuToBIReaderUtils.removeTabsAndTrim(reader.readLine());
    
    if (time == null || !time.contains("time ="))
      throw new TextGridSyntaxErrorException("missing time at line: " + reader.getLineNumber());

    time = time.replace("time = ", "").trim();

    String mark = AuToBIReaderUtils.removeTabsAndTrim(reader.readLine());
    
    if (mark == null || !mark.contains("mark ="))
      throw new TextGridSyntaxErrorException("missing mark at line: " + reader.getLineNumber());

    mark = mark.replaceAll("mark = \"(.*?)\"", "$1").trim();
    
   
    if(!mark.matches("\\[[\\w@]+\\]")){
    
        Region region = new Region(lastTime, Double.valueOf(time));
        region.setFile(reader.getFilename());
        
        if(mark.equals("<P>"))
            mark = "sil";
        region.setLabel(mark);
        regions.add(region);
    
    }
    lastTime = Double.valueOf(time);

  }
  
  public void wordsToBreaksForDirndl(Tier wordTier){
      is_point_tier = true;
      name = "breaks";
      
      for(Region r: wordTier.getRegions()){
          
          if(r.getLabel().equals("sil")){
              if(regions.size()>0)
                regions.get(regions.size()-1).setLabel("4");
              
          }else{
              
              Region region = new Region(r.getEnd());
              region.setFile(r.getFile());
              region.setLabel("1");
              regions.add(region);
          }
          
      }
      
      
  }
  
  public void mergeWith(Tier t){
      is_point_tier = true;
      name = "breaks";
      
      List<Region> oldRegions = regions;
      regions =new LinkedList<Region>();
      
      int i=0;
      int j=0;
      
      while(j<oldRegions.size() && i<t.getRegions().size()){
          
          if(oldRegions.get(j).getStart()>t.getRegions().get(i).getStart()){
              regions.add(t.getRegions().get(i));
              i++;
          }
          else{
              regions.add(oldRegions.get(j));
              j++;
          }
      }
      
      
      if(j==oldRegions.size()){
          while(i<t.getRegions().size())
              regions.add(t.getRegions().get(i++));
      }else if(i==t.getRegions().size()){
          while(j<oldRegions.size())
              regions.add(oldRegions.get(j++));
      }
      
     
      
  }
  
  
}
