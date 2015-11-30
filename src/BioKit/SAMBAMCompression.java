package BioKit;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;

public class SAMBAMCompression 
{
	public static final int COL_EJCCM_SEQNAME = 1;
	public static final int COL_EJCCM_START = 2;
	public static final int COL_EJCCM_STOP = 3;
	public static final int COL_EJCCM_COMB_ID = 4;
	public static final int COL_EJCCM_COUNTS = 5;
	public static final String ELEMENT_ID_SEPARATOR = ";";
	
	private String samBamSortedFile;
	private String outputFile;
	private String gtfSortedFile;
	private boolean processMatesInAscendingStartPositionOrder;
	private boolean ignoreReadMappingsToDifferentChromosomes;
	protected boolean ignoreReadMappingsToDifferentGenes;
	protected boolean onlyUsePairedEndReadsWithBothMatesMapped;
	protected boolean generateStrandSpecificCombIDs;
	
	public static final int TRANSCRIPT_MMSEQ_SCORES = 0;
	public static final int TRANSCRIPT_ISOFORM_LEVELS = 1;
	public static final int GENE_GENE_SCORES = 2;
	
	public SAMBAMCompression(String samBamSortedFile, String gtfSortedFile, String outputFile, 
			boolean processMatesInAscendingStartPositionOrder, 
			boolean ignoreReadMappingsToDifferentChromosomes,
			boolean ignoreReadMappingsToDifferentGenes,
			boolean onlyUsePairedEndReadsWithBothMatesMapped,
			boolean generateStrandSpecificCombIDs) throws Exception
	{
		this.samBamSortedFile = samBamSortedFile;
		this.outputFile = outputFile;
		this.processMatesInAscendingStartPositionOrder = processMatesInAscendingStartPositionOrder;
		this.ignoreReadMappingsToDifferentChromosomes = ignoreReadMappingsToDifferentChromosomes;
		this.ignoreReadMappingsToDifferentGenes = ignoreReadMappingsToDifferentGenes;
		this.onlyUsePairedEndReadsWithBothMatesMapped = onlyUsePairedEndReadsWithBothMatesMapped;
		this.generateStrandSpecificCombIDs = generateStrandSpecificCombIDs;
		this.gtfSortedFile = gtfSortedFile;
	}
	
	public void runEJCCMCompression() throws IOException
	{
		System.out.println("Running EJCCM compression");
		MappedElementProcessingSAMBAMEJCCM mappingStoringMethod = new MappedElementProcessingSAMBAMEJCCM();
		runCompression(mappingStoringMethod);
	}
	
	public void runTCGACompression()
	{
		System.out.println("Running TCGA compression");
		MappedElementProcessingSAMBAMTCGA mappingStoringMethod = new MappedElementProcessingSAMBAMTCGA();
		runCompression(mappingStoringMethod);
	}
	
	private void runCompression(MappedElementProcessingSAMBAMCompression mappingStoringMethod)
	{
		GTFMapperSorted gtfMapperSorted = new GTFMapperSorted(gtfSortedFile);
		gtfMapperSorted.setSetCorrespondingExonsForExonicRegionForGenes(false);
		gtfMapperSorted.setSetCorrespondingTranscriptsForExonicRegionForGenes(false);
		SortedSAMBAMReadProcessor readProcessor = new SortedSAMBAMReadProcessor(samBamSortedFile, gtfMapperSorted, 
				mappingStoringMethod, processMatesInAscendingStartPositionOrder, ignoreReadMappingsToDifferentChromosomes, 
				ignoreReadMappingsToDifferentGenes, onlyUsePairedEndReadsWithBothMatesMapped, generateStrandSpecificCombIDs);
		readProcessor.setSetCorrespondingExonsForJunctions(false);
		readProcessor.setSetCorrespondingTranscriptsForJunctions(false);
		HashMap<String, Integer[]> featureCombinationToStartStopAndOccurrence = mappingStoringMethod.getCombinationIDsStartTopAndOccurrenceData();
		try 
		{
			readProcessor.run();
		}
		catch (IOException e) 
		{
			e.printStackTrace();
			System.exit(0);
		}
		HashMap<String, String> featureCombinationToSeqName = mappingStoringMethod.getCombinationIDSeqNames();
		printToFile(this.outputFile, featureCombinationToStartStopAndOccurrence,
				featureCombinationToSeqName, readProcessor.isPairedEndData());
	}
	
	/**This method print the results tab-delimited to a text file and subsequently
	 * prepares this file for tabix.
	 * 
	 * @param pathToOutputFile
	 * @param featureToStartStopAndOccurrence
	 * @param featureToSeqName
	 */
	private void printToFile(String pathToOutputFile,
			HashMap<String, Integer[]> featureCombinationToStartStopAndOccurrence,
			HashMap<String, String> featureCombinationToSeqName, 
			boolean pairedEndData)
	{
		System.out.println("Printing results to file: " + pathToOutputFile);
		//columns in compressed file:
		//SEQ_NAME	START	STOP	MAPPED_FEATURE_IDs	OCCURRENCE
		//If a read maps to the several features the IDs of these features are concatenated
		//and thereby separated by ";".
		try
		{
			PrintWriter writer = new PrintWriter(pathToOutputFile);
			writer.println("# Original SAM/BAM file:\t" + this.samBamSortedFile);
			
			if (pairedEndData) writer.println("# Read type:\tpaired-end");
			else writer.println("# Read type:\tsingle-end");
			
			writer.println("# Reference transcriptome:\t" + this.gtfSortedFile);
			for (String combID : featureCombinationToStartStopAndOccurrence.keySet())
			{
				Integer[] combInfo = featureCombinationToStartStopAndOccurrence.get(combID);	
				writer.println(featureCombinationToSeqName.get(combID) + "\t" 
						+ combInfo[MappedElementProcessingSAMBAMCompression.INDEX_START] + "\t" 
						+ combInfo[MappedElementProcessingSAMBAMCompression.INDEX_STOP] 
						+ "\t" + combID + "\t" + combInfo[MappedElementProcessingSAMBAMCompression.INDEX_OCCURRENCE]);
			}
			writer.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.exit(0);
		}
	}
}

