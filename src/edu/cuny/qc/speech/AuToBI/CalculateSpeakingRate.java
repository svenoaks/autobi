///*
// * To change this template, choose Tools | Templates
// * and open the template in the editor.
// */
//package edu.cuny.qc.speech.AuToBI;
//
//import edu.cuny.qc.speech.AuToBI.core.Contour;
//import edu.cuny.qc.speech.AuToBI.core.Pair;
//import edu.cuny.qc.speech.AuToBI.core.WavData;
//import edu.cuny.qc.speech.AuToBI.io.WavReader;
//import java.util.ArrayList;
//
///**
// *
// * @author vitisoto
// */
//public class CalculateSpeakingRate {
//    
//    double silenceThresholdDB = -25;
//    double minimumDipBetweenPeaks = 1;
//    double minimumPauseDurationS = 1;
//    boolean keepSoundFilesAndTextGrids = false;
//    String fileName = null;
//    
//    WavReader wr;
//    
//    public void CalculateSpeakingRate(){
//        wr=new WavReader();
//    }
//    
//    public void init(String[] args) {
//       
//    
//    }
//    
//    public static void main(String[] args) {
//        CalculateSpeakingRate csr = new CalculateSpeakingRate();
//        csr.init(args);
//        csr.run();
//    }
//    
//    public void run(){
//        
//        try{
//            //read the .wav file
//            WavData wd = wr.read(fileName);
//            //extract the pitch contour
//            PitchExtractor pe = new PitchExtractor(wd);
//            //In Praat script: 0.02 30 4 no 0.03 0.25 0.01 0.35 0.25 450
//            Contour pitchContour = pe.soundToPitchAc(0.02,30,3,4,0.03,0.25,0.01,0.35,0.25,450);
//            
//            double duration = wd.getDuration();
//            double startTime = 0.0;//wd.getgetDuration();
//            
//            IntensityExtractor ie = new IntensityExtractor(wd);
//            Contour intensityContour = ie.soundToIntensity(50, 0, true);
//            
//            double start = intensityContour.getStart();
//            int nframes = intensityContour.size();
//            double end = intensityContour.get(nframes-1);
//            
//            double minint = intensityContour.getPercentile(0);
//            double maxint = intensityContour.getPercentile(1);
//            double max99int = intensityContour.getPercentile(0.99);
//            
//            double threshold = max99int + silenceThresholdDB;
//            double threshold2 = maxint - max99int;
//            double threshold3 = silenceThresholdDB - threshold2;
//            if(threshold < minint)
//                threshold = minint;
//                
//            
//            //            # get pauses (silences) and speakingtime
////   select soundid
////   To TextGrid (silences)... 100 0 threshold3 minpause 0.1 silent sounding
////   textgridid = selected("TextGrid")
////   silencetierid = Extract tier... 1
////   silencetableid = Down to TableOfReal... sounding
////   nsounding = Get number of rows
////   npauses = 'nsounding'
////   speakingtot = 0
////   for ipause from 1 to npauses
////      beginsound = Get value... 'ipause' 1
////      endsound = Get value... 'ipause' 2
////      speakingdur = 'endsound' - 'beginsound'
////      speakingtot = 'speakingdur' + 'speakingtot'
////   endfor
//            ArrayList<Pair<Double,Double>> intervals =ie.getListSoundsSilences(100, 0 , threshold3, minimumPauseDurationS,0.1);
//            double speakingtot=0.0;
//            double lastTime=wd.t0;
//            
//            for(int i=0;i<intervals.size();i++){
//                
//                if(intervals.get(i).second==1)
//                    speakingtot +=(intervals.get(i).second-lastTime);
//                
//                lastTime = intervals.get(i).first;
//            }
//            
//
////   select 'intid'
////   Down to Matrix
////   matid = selected("Matrix")
////   # Convert intensity to sound
////   To Sound (slice)... 1
////   sndintid = selected("Sound")
//
//            
//            
////   # use total duration, not end time, to find out duration of intdur
////   # in order to allow nonzero starting times.
////   intdur = Get total duration
////   intmax = Get maximum... 0 0 Parabolic
//     double intdur;       
//     double intmax;
//     
////
////   # estimate peak positions (all peaks)
////   To PointProcess (extrema)... Left yes no Sinc70
////   ppid = selected("PointProcess")
////
////   numpeaks = Get number of points
////
////   # fill array with time points
////   for i from 1 to numpeaks
////       t'i' = Get time from index... 'i'
////   endfor
////
////
////   # fill array with intensity values
////   select 'sndintid'
////   peakcount = 0
////   for i from 1 to numpeaks
////       value = Get value at time... t'i' Cubic
////       if value > threshold
////             peakcount += 1
////             int'peakcount' = value
////             timepeaks'peakcount' = t'i'
////       endif
////   endfor
////
////
////   # fill array with valid peaks: only intensity values if preceding
////   # dip in intensity is greater than mindip
////   select 'intid'
////   validpeakcount = 0
////   currenttime = timepeaks1
////   currentint = int1
////
////   for p to peakcount-1
////      following = p + 1
////      followingtime = timepeaks'following'
////      dip = Get minimum... 'currenttime' 'followingtime' None
////      diffint = abs(currentint - dip)
////
////      if diffint > mindip
////         validpeakcount += 1
////         validtime'validpeakcount' = timepeaks'p'
////      endif
////         currenttime = timepeaks'following'
////         currentint = Get value at time... timepeaks'following' Cubic
////   endfor
////
////
////   # Look for only voiced parts
////   select 'soundid'
//// #  To Pitch (ac)... 0.02 30 4 no 0.03 0.25 0.01 0.35 0.25 450
////   # keep track of id of Pitch
//// #  pitchid = selected("Pitch")
////
////   voicedcount = 0
////   for i from 1 to validpeakcount
////      querytime = validtime'i'
////
////      select 'textgridid'
////      whichinterval = Get interval at time... 1 'querytime'
////      whichlabel$ = Get label of interval... 1 'whichinterval'
////
////      select 'pitchid'
////      value = Get value at time... 'querytime' Hertz Linear
////
////      if value <> undefined
////         if whichlabel$ = "sounding"
////             voicedcount = voicedcount + 1
////             voicedpeak'voicedcount' = validtime'i'
////         endif
////      endif
////   endfor
////
////  
////   # calculate time correction due to shift in time for Sound object versus
////   # intensity object
////   timecorrection = originaldur/intdur
////
////   # Insert voiced peaks in TextGrid
////   if showtext > 0
////      select 'textgridid'
////      Insert point tier... 1 syllables
////     
////      for i from 1 to voicedcount
////          position = voicedpeak'i' * timecorrection
////          Insert point... 1 position 'i'
////      endfor
////   endif
////
////   # clean up before next sound file is opened
////    select 'intid'
////    plus 'matid'
////    plus 'sndintid'
////    plus 'ppid'
////    plus 'silencetierid'
////    plus 'silencetableid'
////
////    Remove
////    if showtext < 1
////       select 'textgridid'
////       Remove
////    endif
////
////# summarize results in Info window
////   speakingrate = 'voicedcount'/'originaldur'
////   articulationrate = 'voicedcount'/'speakingtot'
////   npause = 'npauses'-1
////   asd = 'speakingtot'/'voicedcount'
////
////   speakingrate$ = "'speakingrate'"
////   voicedcount$ = "'voicedcount'"
////   articulationrate$ = "'articulationrate'"
////   speakingtot$ = "'speakingtot'"
////
////# voicedcount$ > 'directory$'/'name$'.nsyll
////# speakingtot$ > 'directory$'/'name$'.st
////   resultline$ = "='articulationrate' sr= 'speakingrate' vc = 'voicedcount' od='originaldur'ar='articulationrate' st='speakingtot'"
////   resultline$ > articulationrate.inf
////
////printline 'filename$','articulationrate'
////
////select Pitch sound
////Remove
////select Sound sound
////Remove
//
//            
//            
//            
//            
//            
//        }catch(Exception e){
//            
//        }
//        
//    }
//    
//}
