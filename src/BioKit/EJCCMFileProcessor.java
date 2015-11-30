package BioKit;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.TreeSet;
import java.util.Vector;
import java.util.zip.GZIPInputStream;

public class EJCCMFileProcessor 
{
	private String ejccmFile;
	private ReadMappingEvaluatorLight readMappingEvaluator;
	private MappedElementProcessingEJCCM mappedElementProcessing;
	
	public EJCCMFileProcessor(String ejccmFile, ReadMappingEvaluatorLight evaluator,
			MappedElementProcessingEJCCM mappedElementProcessing)
	{
		this.ejccmFile = ejccmFile;
		this.readMappingEvaluator = evaluator;
		this.mappedElementProcessing = mappedElementProcessing;
	}
	
	public void run() throws IOException
	{
		BufferedReader reader = null;
		if (ejccmFile.endsWith(".gz")) 
		{
			GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(ejccmFile));
			reader = new BufferedReader(new InputStreamReader(gzip));
			gzip.close();
		}
		else reader = new BufferedReader(new FileReader(ejccmFile));
		
		String curLine = null;
		int lineCounter = 0;				
		while ((curLine = reader.readLine()) != null)
		{
			if(lineCounter % 100000 == 0)
			{
				System.out.println("Processed " + lineCounter + " lines...");
			}			
			lineCounter++;

			// the lines starting with a hash symbol
			// contain general informatin that can be ignore here.
			if (curLine.startsWith("#")) continue;
			
			//DEBUG
			if (curLine.startsWith("*")) continue;
			//END DEBUG
			
			String[] curLineContent = curLine.split("\\s+");
			String combinationID = curLineContent[SAMBAMCompression.COL_EJCCM_COMB_ID - 1];
			int counts = Integer.parseInt(curLineContent[SAMBAMCompression.COL_EJCCM_COUNTS - 1]);
			
			Vector<String> mappedElementsCurCombID = getTranscriptsForCombinationID(combinationID);
			mappedElementProcessing.processElements(mappedElementsCurCombID, combinationID, counts);
		}
		
		mappedElementProcessing.finishProcessing();
	}
		
	private Vector<String> getTranscriptsForCombinationID(String combinationID)
	{
		Vector<String> transcriptsCompatibleWithCombID = new Vector<String>();
				
		TreeSet<String> compatibleTranscripts = getTranscriptsCompatibleWithCombination(combinationID, 
				this.readMappingEvaluator, false);
		
		transcriptsCompatibleWithCombID.addAll(compatibleTranscripts);	
		return transcriptsCompatibleWithCombID;
	}
	
	public TreeSet<String> getTranscriptsCompatibleWithCombination(String combinationID, ReadMappingEvaluatorLight evaluator,
			boolean transcriptsMustContainExonMatchingExonicRegionExactly)
	{
		String[] features = combinationID.split(SAMBAMCompression.ELEMENT_ID_SEPARATOR);
		
		HashMap<String, Integer> transcriptsToNumMappedFeatures = new HashMap<String, Integer>();
		
		//DEBUG
		boolean debug = false;
		String debugString = "";
		//END DEBUG
		
		for (String featureID : features)
		{
			String[] featureIDParts = featureID.split(",");
			// If the featureID does not contain any ',' then
			// the featureID describes an exon.
			// If the feature ID describes a junction then the 
			// junction start position and stop position are
			// separated by a comma s.t. the array has length 2.
			if (featureIDParts.length == 1)
			{
				featureIDParts = featureID.split(":");
				String seqName = featureIDParts[0];
				String exonicRegionDescription = null;
				exonicRegionDescription = featureIDParts[1].replaceAll("\\(", "").replaceAll("\\)", ""); 
				String[] exonStartStop = exonicRegionDescription.split("-"); 
				int exonicRegionStart = Integer.parseInt(exonStartStop[0]);
				int exonicRegionStop = Integer.parseInt(exonStartStop[1]);

				//DEBUG
				debugString += "Exonic region: " + seqName + ":" + exonicRegionStart + "-" + exonicRegionStop + "\n";
				debugString += "Mapped transcripts:\n";
				//END DEBUG
				
				TreeSet<String> mappedTranscripts = evaluator.mappingToCompatibleTranscripts(exonicRegionStart, exonicRegionStop, seqName,
						transcriptsMustContainExonMatchingExonicRegionExactly);					
				if (mappedTranscripts == null) continue;
				
				for (String transcript : mappedTranscripts)
				{
					//DEBUG
					debugString += "\t" + transcript + "\n";
					//END DEBUG
					if (!transcriptsToNumMappedFeatures.containsKey(transcript)) transcriptsToNumMappedFeatures.put(transcript, 0);
					transcriptsToNumMappedFeatures.put(transcript, transcriptsToNumMappedFeatures.get(transcript) + 1);
				}
				
			}
			else if (featureIDParts.length == 2)
			{
				String junctionSeqName = null;
				String junction = null;
				Integer junctionStart = null;
				Integer junctionStop = null;
				for (String junctionInfo : featureIDParts)
				{
					String[] junctionInfoParts = junctionInfo.split(":");
					if (junctionSeqName == null) junctionSeqName = junctionInfoParts[0];
					if (junction == null)
					{
						junction = junctionInfoParts[1];
						junctionStart = Integer.parseInt(junctionInfoParts[1]);
					}
					else
					{
						junction += "-" + junctionInfoParts[1];
						junctionStop = Integer.parseInt(junctionInfoParts[1]);
					}
				}
				//DEBUG
				debugString += "Junction: " + junctionSeqName + ":" + junctionStart + "-" + junctionStop + "\n";
				debugString += "Mapped transcripts:\n";
				//END DEBUG
				
				TreeSet<String> mappedTranscripts = evaluator.mappingToCompatibleTranscripts(junction, junctionSeqName);
				if (mappedTranscripts == null) continue;
				
				for (String transcript : mappedTranscripts)
				{
					//DEBUG
					debugString += "\t" + transcript + "\n";
					//END DEBUG
					if (!transcriptsToNumMappedFeatures.containsKey(transcript)) transcriptsToNumMappedFeatures.put(transcript, 0);
					transcriptsToNumMappedFeatures.put(transcript, transcriptsToNumMappedFeatures.get(transcript) + 1);
				}
			}	
			else
			{
				System.err.println("Cannot read feature ID : " + featureID);
				System.exit(0);
			}
		}
		
		/* Write transcripts mapped to combination to file.
		 * Originally the hits-file for MMSeq contains for each
		 * read the transcripts that could be mapped to that read.
		 * 
		 * The compressed sam/bam file basically counts how many reads
		 * mapped to each combination. In order to translate that 
		 * information to the hits-file format the combination with
		 * the mapped transcripts is written as often to the hits-file
		 * as the counts specify. The combination ID is then extended
		 * with a number indicating the repeat. 
		 */
		//DEBUG
		debugString += "Final mapped transcripts:\n";
		//END DEBUG
		TreeSet<String> transcriptsMappedToCombination = new TreeSet<String>();
		for (String transcript : transcriptsToNumMappedFeatures.keySet())
		{
			if (transcriptsToNumMappedFeatures.get(transcript) == features.length)
			{
				transcriptsMappedToCombination.add(transcript);
				//DEBUG
				debugString += "\t" + transcript + "\n";
				//END DEBUG
			}
		}
		
		//DEBUG
		if (debug) System.out.println(debugString);
		//END DEBUG
		
		return transcriptsMappedToCombination;
	}
}
