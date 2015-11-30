package BioKit;

import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.Vector;

public class MappedElementProcessingDebug extends MappedElementProcessingSAMBAMCompression
{
	private ReadMappingEvaluatorLight readMappingEvaluator;
	private EJCCMFileProcessor ejccmFileProcessor;
	private int caseCounter;
	private HashMap<String, Integer[]> transcriptToBAMAndEJCCMReadCounts;
	private final int indexBAMReadCounts = 0;
	private final int indexEJCCMReadCounts = 1;
	private int readsWithSeveralCombinationIDs;
	private String transcriptReadCountOutputFile;
	//DEBUG
	public String debuggingMessage;
	//END DEBUG
	
	public MappedElementProcessingDebug(ReadMappingEvaluatorLight evaluator,
			String transcriptReadCountOutputFile) 
	{
		this.readMappingEvaluator = evaluator;
		this.ejccmFileProcessor = new EJCCMFileProcessor(null, evaluator, null);
		this.caseCounter = 0;
		this.transcriptToBAMAndEJCCMReadCounts = new HashMap<String, Integer[]>();
		this.transcriptReadCountOutputFile = transcriptReadCountOutputFile;
	}
	
	@Override
	public void storeProcessingResultsPairedEnd(
			HashSet<String> genesMappedByBothReads,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsFirstMate,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsSecondMate,
			String referenceFirstMate, String referenceSecondMate,
			boolean ignoreReadMappingsToDifferentGenes) 
	{
		// Get transcripts directly
		HashSet<String> transcriptsFoundDirectly = getTranscriptsDirectlyFromPairedEndRead(genesMappedByBothReads, 
				genesToCoveredElementsFirstMate, genesToCoveredElementsSecondMate, referenceFirstMate, 
				referenceSecondMate, ignoreReadMappingsToDifferentGenes);
		// Get transcripts via EJCCM
		HashSet<String> transcriptsFoundViaEJCCM = getTranscriptsViaEJCCMFromPairedEndRead(genesMappedByBothReads, 
				genesToCoveredElementsFirstMate, genesToCoveredElementsSecondMate, referenceFirstMate, 
				referenceSecondMate, ignoreReadMappingsToDifferentGenes);
		// Check if the same
		compareTranscriptsFromDirectAndEJCCM(transcriptsFoundDirectly, transcriptsFoundViaEJCCM);
		// if not alarm!!
		//DEBUG
		this.debuggingMessage = "";
		//END DEBUG

	}
	
	public HashSet<String> getTranscriptsDirectlyFromPairedEndRead(HashSet<String> genesMappedByBothReads,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsFirstMate,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsSecondMate,
			String referenceFirstMate, String referenceSecondMate,
			boolean ignoreReadMappingsToDifferentGenes)
	{
		HashSet<String> transcripts = new HashSet<String>();
		for (String gene : genesMappedByBothReads)
		{
			//DEBUG
			this.debuggingMessage += "\nCASE--------" + this.caseCounter + "--------\n";
			this.debuggingMessage += "Directly mapped transcripts: \n";
			//END DEBUG
			HashMap<String, HashSet<String>> transcriptToMappableElements = new HashMap<String, HashSet<String>>();
			HashSet<String> mappableElements = new HashSet<String>();
			//DEBUG
			this.debuggingMessage += "\tgene: " + gene + "\n";
			//END DEBUG
			for (MappableElement m : genesToCoveredElementsFirstMate.get(gene))
			{
				mappableElements.add(m.getDescription());
				//DEBUG
				this.debuggingMessage += "\tFirst Read: \n";
				this.debuggingMessage += "\tMappableElement: " + m.getDescription() + "\n";
				//END DEBUG
				
				for (String transcript : m.getMappedTranscripts())
				{
					//DEBUG
					this.debuggingMessage += "\t\t" + transcript + "\n";
					//END DEBUG
					if (!transcriptToMappableElements.containsKey(transcript))
					{
						transcriptToMappableElements.put(transcript, new HashSet<String>());
					}
					transcriptToMappableElements.get(transcript).add(m.getDescription());
				}
			}
			
			for (MappableElement m : genesToCoveredElementsSecondMate.get(gene))
			{
				mappableElements.add(m.getDescription());
				//DEBUG
				this.debuggingMessage += "\tSecond Read: \n";
				this.debuggingMessage += "\tMappableElement: " + m.getDescription() + "\n";
				//END DEBUG
				for (String transcript : m.getMappedTranscripts())
				{
					//DEBUG
					this.debuggingMessage += "\t\t" + transcript + "\n";
					//END DEBUG
					if (!transcriptToMappableElements.containsKey(transcript))
					{
						transcriptToMappableElements.put(transcript, new HashSet<String>());
					}
					transcriptToMappableElements.get(transcript).add(m.getDescription());
				}
			}
			
			//DEBUG
			this.debuggingMessage += "Final mapped transcripts:\n";
			//END DEBUG
			for (String transcript : transcriptToMappableElements.keySet())
			{
				if (transcriptToMappableElements.get(transcript).size() == mappableElements.size())
				{
					transcripts.add(transcript);
					//DEBUG
					this.debuggingMessage += "\t" + transcript + "\n";
					//END DEBUG
				}
			}
		}
		
		//DEBUG
		this.debuggingMessage += "Totally mapped transcripts\n";
		for (String transcript : transcripts)
		{
			this.debuggingMessage += "\t" + transcript + "\n";
		}
		//END DEBUG
		return transcripts;
	}
	
	public HashSet<String> getTranscriptsViaEJCCMFromPairedEndRead(HashSet<String> genesMappedByBothReads,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsFirstMate,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsSecondMate,
			String referenceFirstMate, String referenceSecondMate,
			boolean ignoreReadMappingsToDifferentGenes)
	{
		HashSet<String> transcripts = new HashSet<String>();
				
		HashSet<String> combIDs = generateCombinationIDsPairedEnd(genesMappedByBothReads, 
				genesToCoveredElementsFirstMate, genesToCoveredElementsSecondMate, 
				referenceFirstMate, referenceSecondMate, ignoreReadMappingsToDifferentGenes);
		
		if (combIDs.size() > 1) readsWithSeveralCombinationIDs++;
		
		for (String combID : combIDs)
		{
			transcripts.addAll(getTranscriptsForCombinationID(combID));
			//DEBUG
			this.debuggingMessage += "#######\n";
			this.debuggingMessage += "Combination ID: " + combID + "\n";
			this.debuggingMessage += "Transcripts for combination ID:\n";
			for (String transcript : getTranscriptsForCombinationID(combID))
			{
				this.debuggingMessage += "\t" + transcript + "\n";
			}
			this.debuggingMessage += "#######\n";
			//END DEBUG
		}
		
		return transcripts;
	}	

	private HashSet<String> generateCombinationIDsPairedEnd(HashSet<String> genesMappedByBothReads,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsFirstMate,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsSecondMate,
			String referenceFirstMate, String referenceSecondMate,
			boolean ignoreReadMappingsToDifferentGenes)
	{
		HashSet<String> combinationIDs = new HashSet<String>();
		/* The mates of a paired-end read can map to different genes or chromosomes.
		 * If the set of genes the mates map to overlaps then the mappings to
		 * the genes in the intersection are taken for combination ID generation.
		 * In this case a combination ID per gene is generated. 
		 * If there exists no gene the mates commonly map to then combination IDs
		 * for the combination of the genes are generated, since this can represent
		 * gene fusion events.
		 */
		if (referenceFirstMate.equals("*"))
		{
			referenceFirstMate = referenceSecondMate;
		}
		if (genesMappedByBothReads.isEmpty() && !ignoreReadMappingsToDifferentGenes
				&& !genesToCoveredElementsFirstMate.isEmpty() && !genesToCoveredElementsSecondMate.isEmpty())
		{
			for (String geneFirstMate : genesToCoveredElementsFirstMate.keySet())
			{
				HashMap<Integer, Integer> minStartMaxStop = new HashMap<Integer, Integer>();
				String combinationID = generateCombinationID(genesToCoveredElementsFirstMate.get(geneFirstMate), 
						minStartMaxStop);
				
				for (String geneSecondMate : genesToCoveredElementsSecondMate.keySet())
				{
					combinationID += SAMBAMCompression.ELEMENT_ID_SEPARATOR + 
							generateCombinationID(genesToCoveredElementsSecondMate.get(geneSecondMate), 
							minStartMaxStop);
					
					combinationIDs.add(combinationID);
				}
			}		
		}
		else if (genesToCoveredElementsFirstMate.isEmpty())
		{			
			for (String geneSecondMate : genesToCoveredElementsSecondMate.keySet())
			{
				HashMap<Integer, Integer> minStartMaxStop = new HashMap<Integer, Integer>();
				String combinationID = generateCombinationID(genesToCoveredElementsSecondMate.get(geneSecondMate), 
						minStartMaxStop);

				combinationIDs.add(combinationID);
			}
		}
		else if (genesToCoveredElementsSecondMate.isEmpty())
		{			
			for (String geneFirstMate : genesToCoveredElementsFirstMate.keySet())
			{
				HashMap<Integer, Integer> minStartMaxStop = new HashMap<Integer, Integer>();
				String combinationID = generateCombinationID(genesToCoveredElementsFirstMate.get(geneFirstMate), 
						minStartMaxStop);

				combinationIDs.add(combinationID);
			}
		}
		else
		{
			for (String gene : genesMappedByBothReads)
			{					
				HashMap<Integer, Integer> minStartMaxStop = new HashMap<Integer, Integer>();
				// It is made sure that each covered element is stored only once.
				HashSet<String> coveredElementsAlreadySeen = new HashSet<String>();
				HashSet<MappableElement> coveredElementsUnique = new HashSet<MappableElement>();
				for (MappableElement m : genesToCoveredElementsFirstMate.get(gene))
				{
					coveredElementsAlreadySeen.add(m.getDescription());
					coveredElementsUnique.add(m);
				}
				for (MappableElement m : genesToCoveredElementsSecondMate.get(gene))
				{
					if (!coveredElementsAlreadySeen.contains(m.getDescription()))
					{
						coveredElementsAlreadySeen.add(m.getDescription());
						coveredElementsUnique.add(m);
					}
				}				
				String combinationID = generateCombinationID(coveredElementsUnique, minStartMaxStop);
	
				combinationIDs.add(combinationID);						
			}
		}
		return combinationIDs;
	}
	
	private Vector<String> getTranscriptsForCombinationID(String combinationID)
	{
		Vector<String> transcriptsCompatibleWithCombID = new Vector<String>();
				
		TreeSet<String> compatibleTranscripts = this.ejccmFileProcessor.getTranscriptsCompatibleWithCombination(combinationID, 
				this.readMappingEvaluator, false);
		
		for (String transcript : compatibleTranscripts)
		{
			transcriptsCompatibleWithCombID.add(transcript);
		}
			
		return transcriptsCompatibleWithCombID;
	}

	@Override
	public void storeProcessingResultsSingleEnd(
			HashMap<String, HashSet<MappableElement>> genesToCoveredElements,
			String reference) 
	{
		// Get transcripts directly
		HashSet<String> transcriptsFoundDirectly = getTranscriptsDirectlyFromSingleEndRead(genesToCoveredElements, reference);
		// Get transcripts via EJCCM
		HashSet<String> transcriptsFoundViaEJCCM = getTranscriptsViaEJCCMFromSingleEndRead(genesToCoveredElements, reference);
		
		// Check if the same
		compareTranscriptsFromDirectAndEJCCM(transcriptsFoundDirectly, transcriptsFoundViaEJCCM);
		//DEBUG
		debuggingMessage = "";
		//END DEBUG
	}
	
	private HashSet<String> generateCombinationIDsSingleEnd(HashMap<String, HashSet<MappableElement>> genesToCoveredElements,
			String reference)
	{
		HashSet<String> combinationIDs = new HashSet<String>();
		for (String gene : genesToCoveredElements.keySet())
		{
			HashMap<Integer, Integer> minStartMaxStop = new HashMap<Integer, Integer>();
			String combinationID = generateCombinationID(genesToCoveredElements.get(gene), 
					minStartMaxStop);

			combinationIDs.add(combinationID);
		}
		return combinationIDs;
	}
	
	public HashSet<String> getTranscriptsDirectlyFromSingleEndRead(HashMap<String, HashSet<MappableElement>> genesToCoveredElements,
			String reference)
	{
		HashSet<String> transcripts = new HashSet<String>();
		for (String gene : genesToCoveredElements.keySet())
		{			
			//DEBUG
			debuggingMessage += "Directly mapped transcripts: \n";
			//END DEBUG
			/* In case of junction read there are several mappable
			 * elements returned. The exonic regions before the
			 * junction, the exonic after the junction and the
			 * junction itsself. Thus, also here is has to be 
			 * counted to how many mappable elements the transcripts
			 * belong. 
			 */
			HashMap<String, HashSet<String>> transcriptToMappableElements = new HashMap<String, HashSet<String>>();
			HashSet<String> mappableElements = new HashSet<String>();
			//DEBUG
			debuggingMessage += "\tgene: " + gene + "\n";
			//END DEBUG
			for (MappableElement m : genesToCoveredElements.get(gene))
			{
				mappableElements.add(m.getDescription());
				//DEBUG
				debuggingMessage += "\tMappableElement: " + m.getDescription() + "\n";
				//END DEBUG
				
				for (String transcript : m.getMappedTranscripts())
				{
					//DEBUG
					debuggingMessage += "\t\t" + transcript + "\n";
					//END DEBUG
					if (!transcriptToMappableElements.containsKey(transcript))
					{
						transcriptToMappableElements.put(transcript, new HashSet<String>());
					}
					transcriptToMappableElements.get(transcript).add(m.getDescription());
				}
			}
			
			//DEBUG
			debuggingMessage += "Final mapped transcripts:\n";
			//END DEBUG
			for (String transcript : transcriptToMappableElements.keySet())
			{
				if (transcriptToMappableElements.get(transcript).size() == mappableElements.size())
				{
					transcripts.add(transcript);
					//DEBUG
					debuggingMessage += "\t" + transcript + "\n";
					//END DEBUG
				}
			}
		}
		return transcripts;
	}
	
	public HashSet<String> getTranscriptsViaEJCCMFromSingleEndRead(HashMap<String, HashSet<MappableElement>> genesToCoveredElements,
			String reference)
	{
		HashSet<String> transcripts = new HashSet<String>();
		HashSet<String> combIDs = generateCombinationIDsSingleEnd(genesToCoveredElements, reference);
		for (String combID : combIDs)
		{
			transcripts.addAll(getTranscriptsForCombinationID(combID));
		}
		return transcripts;
	}
	
	private void compareTranscriptsFromDirectAndEJCCM(HashSet<String> transcriptsDirect,
			HashSet<String> transcriptsViaEJCCM)
	{
		this.caseCounter++;
		
		if (this.caseCounter % 100000 == 0)
		{
			System.out.println("Compared " + this.caseCounter + " cases...");
		}
		
		for (String transcript : transcriptsDirect)
		{
			if (!transcriptsViaEJCCM.contains(transcript))
			{
				printErrorMessage(transcriptsDirect, transcriptsViaEJCCM);
//				System.exit(0);
			}
			
			if (!transcriptToBAMAndEJCCMReadCounts.containsKey(transcript))
			{
				transcriptToBAMAndEJCCMReadCounts.put(transcript, new Integer[2]);
				transcriptToBAMAndEJCCMReadCounts.get(transcript)[indexBAMReadCounts] = 0;
				transcriptToBAMAndEJCCMReadCounts.get(transcript)[indexEJCCMReadCounts] = 0;
			}
			transcriptToBAMAndEJCCMReadCounts.get(transcript)[indexBAMReadCounts]++;
		}
		
		for (String transcript : transcriptsViaEJCCM)
		{
			if (!transcriptsDirect.contains(transcript))
			{
				printErrorMessage(transcriptsDirect, transcriptsViaEJCCM);
//				System.exit(0);
			}
			
			if (!transcriptToBAMAndEJCCMReadCounts.containsKey(transcript))
			{
				transcriptToBAMAndEJCCMReadCounts.put(transcript, new Integer[2]);
				transcriptToBAMAndEJCCMReadCounts.get(transcript)[indexBAMReadCounts] = 0;
				transcriptToBAMAndEJCCMReadCounts.get(transcript)[indexEJCCMReadCounts] = 0;
			}
			transcriptToBAMAndEJCCMReadCounts.get(transcript)[indexEJCCMReadCounts]++;
		}
		
//		if (caseCounter == 5183713)
//		{
//			System.out.println(this.debuggingMessage);
//			System.out.println();
//			System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!NOOOOOOOOOOOOOOOOOOOOOOOOO DIFFERENCE FOUND!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
//			System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!Case " + this.caseCounter + "!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
//			System.out.println("Transcripts Direct");
//			for (String transcript: transcriptsDirect)
//			{
//				System.out.println("\t" + transcript);
//			}
//			System.out.println("Transcripts via EJCCM");
//			for (String transcript: transcriptsViaEJCCM)
//			{
//				System.out.println("\t" + transcript);
//			}
//		}
	}
	
	private void printErrorMessage(HashSet<String> transcriptsDirect,
			HashSet<String> transcriptsViaEJCCM)
	{
		System.out.println(this.debuggingMessage);
		System.out.println();
		System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!DIFFERENCE FOUND!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!Case " + this.caseCounter + "!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		System.out.println("Transcripts Direct");
		for (String transcript: transcriptsDirect)
		{
			System.out.println("\t" + transcript);
		}
		System.out.println("Transcripts via EJCCM");
		for (String transcript: transcriptsViaEJCCM)
		{
			System.out.println("\t" + transcript);
		}
	}

	@Override
	public void finishProcessing() 
	{
		System.out.println("Number of reads with several combination IDs: " + readsWithSeveralCombinationIDs);
		
		System.out.println("Writing transcript read count output file " + this.transcriptReadCountOutputFile);
		try
		{
			PrintWriter writer = new PrintWriter(this.transcriptReadCountOutputFile);
			writer.println("TRANSCRIPT\tBAM_COUNT\tEJCCM_COUNT");
			for (String transcript : this.transcriptToBAMAndEJCCMReadCounts.keySet())
			{
				Integer[] transcriptCounts = this.transcriptToBAMAndEJCCMReadCounts.get(transcript);
				writer.println(transcript + "\t" + transcriptCounts[indexBAMReadCounts] + "\t" + transcriptCounts[indexEJCCMReadCounts]);
			}
			writer.flush();
			writer.close();
		}
		catch (Exception e)
		{
			e.printStackTrace();
			System.exit(0);
		}		
	}

	@Override
	public void setReadID(String readID) {
		// TODO Auto-generated method stub
		
	}
}
