/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.cuny.qc.speech.AuToBI;

import edu.cuny.qc.speech.AuToBI.core.Contour;
import edu.cuny.qc.speech.AuToBI.core.WavData;
import edu.cuny.qc.speech.AuToBI.io.AuToBIFileWriter;
import edu.cuny.qc.speech.AuToBI.io.WavReader;
import edu.cuny.qc.speech.AuToBI.util.SubregionUtils;
import java.sql.Timestamp;
import java.util.Arrays;


/**
 *
 * @author vitisoto
 */
public class ExtractAcousticsSS {
    
    
    private String UNDEF = "undef";
    private String EXSTAC_PRAAT_SCRIPT = "/proj/speech/tools/speechlab/praat_scripts/extractStandardAcoustics2.praat";
    private String PRAAT = "/proj/speech/tools/praat/linux/praat";
    
    private static String WAV_FILE;
    private static double MIN_PITCH;
    private static double MAX_PITCH;
    private static double TRIM_PERCENTAGE;
    private static int N;
    private static double[] FROM;
    private static double[] TO;
    
    private static int OK=0;
    private static int ERR=-1;

    public ExtractAcousticsSS(){
        
    }
    
    public int init(String[] args) {        
        WAV_FILE = args[0];
        MIN_PITCH = Double.valueOf(args[1]);
        MAX_PITCH = Double.valueOf(args[2]);
        TRIM_PERCENTAGE = Double.valueOf(args[3]);
        
        if((args.length-4) % 2 !=0)
            return ERR;
        
        N = (args.length-4)/2;
        
        if(N<1)
            return ERR;
        
        
        FROM = new double[N];
        TO = new double[N];
        for(int i=0;i<N;i++){
            FROM[i]=Double.valueOf(args[2*i+4]);
            System.err.println((i+1)+" "+FROM[i]);
            TO[i]=Double.valueOf(args[2*i+1+4]);
            System.err.println((i+1)+" "+TO[i]);
            
            if(FROM[i]>=TO[i])
                return ERR;
            
            if((i>0) && FROM[i]<TO[i-1])
                return ERR;
        }
        System.err.println("YES");
        return OK;
        
    }
    
    private static double[] trimContourIntensity(Contour contour, Double alpha){
        
        int n = (int)(alpha*contour.size())/100;
        double[] values=getSortedArray(contour.getValues());
        
        //for(int k=0;k<values.length;k++){
        //    System.err.print(" "+values[k]);
        //}
        //System.err.println(" "+n);
        
        return Arrays.copyOfRange(values,n,values.length-n);
        
    }
    
    public static double[] getSortedArray(double[] values){
      double[] sortedArray = Arrays.copyOf(values, values.length);
      Arrays.sort(sortedArray);
      return sortedArray;
      
    }
    
    private static double[] trimContourPitch(Pitch pitch, Double alpha){
        
        //int voicedFrames = pitch.countVoicedFrames();
        double[] voicedFrames = pitch.extractVoicedFrames(); 
        
        int n = (int)(alpha*voicedFrames.length)/100;
        
        
        double[] trimmedValues=getSortedArray(voicedFrames);
        
        //for(int k=0;k<trimmedValues.length;k++){
        //    System.err.print(" "+trimmedValues[k]);
        //}
        //System.err.println(" "+n);
        
        return Arrays.copyOfRange(trimmedValues,n,trimmedValues.length-n);
        
    }
    
    private static double computeMean(double[] values){
        double mean =0;
        for(int i=0;i<values.length;i++){
            mean+=values[i];
        }
        mean=mean/values.length;
        return mean;
    
    }
    
    private static double computeStd(double[] values, double mean){
        double std = 0;
        for(int i=0;i<values.length;i++){
            std+=Math.pow(values[i]-mean,2);
        }
        std = Math.sqrt(std/(values.length-1));
        return std;
    }
    
    public static void main(String[] args){
        
 //       die "Usage: extract_acoustics.pl WavFile MinPitch MaxPitch Trim% From_1 To_1 From_2 To_2...\n"
        ExtractAcousticsSS ea = new ExtractAcousticsSS();
        if(ea.init(args)==ERR){
            System.out.println("This script receives several parameters: "
                + "WavFile   full path of the wav file to process"
                + "MinPitch  speaker's estimated minimum pitch (eg, 50 for male, 75 for female)"
                + "MaxPitch  speaker's estimated maximum pitch (eg, 300 for male, 500 for female)"
                + "Trim%     percentage of the data to trim on each side (eg, 5 will trim the upper 5% and the lower 5% --or 10% in total)."
                + "From_1    start point of interval 1"
                + "To_1      end point of interval 1"
                + "From_2    start point of interval 2"
                + "To_2      end point of interval 2"
                + "and computes a set of standard acoustic features over the given intervals only, ignoring the rest of the wav file."
                + "The intervals must not be empty (i.e. From_i < To_i), and they must"
                + "be sorted incrementally (i.e. To_i <= From_i+1)."
                + "There may be an arbitrary number of intervals, but at least one."
                + "Example:"
                + "./extract_acoustics.pl /proj/speech/projects/games/data/session_01/s01.objects.0.A.wav 50 300 5 69.41 70.66 75.47 76.03");
            
            return;
        }

        //unless ($wav_file and $min_pitch and $max_pitch and $TRIM_PERCENTAGE ne '');
        //die "File not found: '$wav_file'\n" unless -e $wav_file;
        //die "File path must be absolute: '$wav_file'\n" unless substr($wav_file, 0, 1) eq "/";
        
        //Now prepare the praat script to chunk all the intervals
        //and concatenate them together.
        
        //stem for temporary files
        java.util.Date date= new java.util.Date();
	Timestamp ts = new Timestamp(date.getTime());
        String tmp = "/tmp/acoustics-"+ts;
        
        try{
            //open FILE, ">$tmp.praat";
            //print FILE "Open long sound file... $wav_file\n";
            //print FILE "Rename... main\n";
            WavReader wr=new WavReader();
            WavData main=wr.read(WAV_FILE);
            
            
            //for my $i (0..$num_intervals-1) {
            //	print FILE "select LongSound main\n";
            //	print FILE "Extract part... $from[$i] $to[$i] yes\n";
            //	print FILE "Rename... part_$i\n";
            //}
            
            WavData[] part = new WavData[N];
            for(int i=0; i<N;i++){
                part[i] = SubregionUtils.getSlice(main, FROM[i],TO[i]);
            }
            
            //print FILE "select Sound part_0\n";
            //for my $i (1..$num_intervals-1) {
            //	print FILE "plus Sound part_$i\n";
            //}
            //
            //print FILE qq^
            //Concatenate
            //Rename... sound
            WavData sound = SubregionUtils.concatSlices(part);
            
            
            //select Sound sound
            //endTime = Get end time
            //To Pitch... 0 $min_pitch $max_pitch
            //vcd2tot_frames = undefined
            //vcd_frames = Count voiced frames
            //tot_frames = Get number of frames
            //if tot_frames > 0
            //	vcd2tot_frames = vcd_frames / tot_frames
            //endif
            //Down to PitchTier
            //Write to text file... $tmp.PitchTier
            
            //double endTime = sound.getDuration();
            //System.err.println("NUMBER OF CHANNELS: "+sound.numberOfChannels);
            //System.err.println("SAMPLE RATE: "+sound.sampleRate);
            //System.err.println("SAMPLE SIZE: "+sound.sampleSize);
            //System.err.println("T0: "+sound.t0);
            //System.err.println("END TIME "+endTime);
            
            PitchExtractor pe = new PitchExtractor(sound);
            Contour contourPitch = pe.soundToPitch(0, MIN_PITCH, MAX_PITCH);
            Pitch pitch = new Pitch(sound,pe.getPitchFrames(),contourPitch);
            
            double vcd2tot_frames = 0.0;
            int vcd_frames = pitch.countVoicedFrames();
            //System.err.println("EMPTY VALUES: "+contourPitch.size()+" "+contourPitch.contentSize());
            System.err.println("VCD FRAMES "+vcd_frames);
            int tot_frames  = contourPitch.size();
            //System.err.println("TOTAL FRAMES "+tot_frames);
            if(tot_frames>0)
                vcd2tot_frames = vcd_frames*1.0/tot_frames;
            
            //System.err.println("VCD2TOT FRAMES "+vcd2tot_frames);
            //select Sound sound
            //To Intensity... $min_pitch 0 yes
            //Down to IntensityTier
            //Write to text file... $tmp.IntensityTier
            //
            //printline 'vcd2tot_frames:3'
            //^;
            //close FILE;
            IntensityExtractor ie = new IntensityExtractor(sound);
            Contour contourIntensity = ie.soundToIntensity(MIN_PITCH, 0, true);
            
            
            //# Run the script we just wrote, to create the PraatTiers of all
            //# the chunks concatenated, and get the voiced2total frames ratio
            //chomp $vcd2tot_frames;


            //# Read the praat tiers
            //my %PitchTier     = read_PraatTier("$tmp.PitchTier");
            //my %IntensityTier = read_PraatTier("$tmp.IntensityTier");

            //# Remove the top and bottom $TRIM_PERCENTAGE of the data
            //# from the Pitch and Intensity tiers.
            //trim_PraatTier(\%PitchTier,     $TRIM_PERCENTAGE);
            //trim_PraatTier(\%IntensityTier, $TRIM_PERCENTAGE);
            double[] pitchValues = trimContourPitch(pitch, TRIM_PERCENTAGE);
            double[] intensityValues = trimContourIntensity(contourIntensity, TRIM_PERCENTAGE);

            double minPitch = pitchValues[0];
            double maxPitch = pitchValues[pitchValues.length-1];
            double meanPitch = computeMean(pitchValues);
            double stdPitch = computeStd(pitchValues, meanPitch);
            
            
            double minIntensity = intensityValues[0];
            double maxIntensity = intensityValues[intensityValues.length-1];
            double meanIntensity = computeMean(intensityValues);
            double stdIntensity = computeStd(intensityValues, meanIntensity);
            
            String resultsFile = WAV_FILE.substring(0, WAV_FILE.length()-4)+".txt";
            AuToBIFileWriter fw = new AuToBIFileWriter(resultsFile);
            fw.write("duration:--undefined--\n");
            fw.write("F0_MIN:"+minPitch+"\n");
            fw.write("F0_MAX:"+maxPitch+"\n");
            fw.write("F0_MEAN:"+meanPitch+"\n");
            fw.write("F0_MEDIAN:--undefined--\n");
            fw.write("F0_STDV:"+stdPitch+"\n");
            fw.write("F0_MAS:--undefined--\n");
            fw.write("ENG_MAX:"+maxIntensity+"\n");
            fw.write("ENG_MIN:"+minIntensity+"\n");
            fw.write("ENG_MEAN:"+meanIntensity+"\n");
            fw.write("ENG_STDV:"+stdIntensity+"\n");
            fw.write("VCD2TOT_FRAMES:"+vcd2tot_frames+"\n");
            fw.close();

        }catch(Exception e){
            
        }
        


        
        
    }
    
}
