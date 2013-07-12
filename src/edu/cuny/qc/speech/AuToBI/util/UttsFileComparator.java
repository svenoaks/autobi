/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.cuny.qc.speech.AuToBI.util;

import java.util.Comparator;

/**
 *
 * @author vs2411
 */
public class UttsFileComparator implements Comparator<String>{
    
    
    @Override
    public int compare(String a,String b){
                
        double tsa = Double.parseDouble(a.substring(a.lastIndexOf("Line") + 5, a.lastIndexOf("_")));
        double tsb = Double.parseDouble(b.substring(b.lastIndexOf("Line") + 5, b.lastIndexOf("_")));
        
        if(tsa >= tsb){
            return 1;
        }else{
            return -1;
        }
    }
    
}
