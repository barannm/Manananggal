package BioKit;

import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;

public class MappedElementProcessingSAMBAMMMSeq implements MappedElementProcessingSAMBAM 
{
	private PrintWriter hitsFileWriter;
	private String curReadID;
	private int numOfStoredMappings;
	
	public MappedElementProcessingSAMBAMMMSeq(String hitsFile, String unsortedGTFFile)
	{
		initPrintWriter(hitsFile, unsortedGTFFile);
		this.curReadID = null;
		this.numOfStoredMappings = 0;
	}
	
	private void initPrintWriter(String hitsFile, String unsortedGTFFile)
	{
		try 
		{
			System.out.println("Init hits file");
			hitsFileWriter = new PrintWriter(hitsFile);
			// writer mmseq header information
			GTFParser gtfParser = new GTFParser(unsortedGTFFile);
			
			String geneIsoformsInfo = "";
			String transcriptMetaData = "";
			
			System.out.println("Parse GTF file at " + unsortedGTFFile + " to obtain transcript structures");
			
			GTFGene gtfGene = null;
			
			while((gtfGene = gtfParser.nextGene()) != null)
			{
				Gene gene = gtfGene.createGene();
				
				geneIsoformsInfo += "@GeneIsoforms\t" + gene.getGeneID();
				
				for(String transcript : gene.getArrayOfGeneProductNames())
				{
					geneIsoformsInfo += "\t" + transcript;
					transcriptMetaData += "@TranscriptMetaData\t" + transcript + "\t" + gene.getTranscriptCDNALength(transcript) + "\n";
				}
				geneIsoformsInfo += "\n";
			}
			gtfParser.closeReader();
			
			hitsFileWriter.print(transcriptMetaData);
			hitsFileWriter.print(geneIsoformsInfo);
		} 
		catch (Exception e) 
		{
			e.printStackTrace();
			System.exit(0);
		} 
	}
		
	public void storeProcessingResultsPairedEnd(
			HashSet<String> genesMappedByBothReads,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsFirstMate,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsSecondMate,
			String referenceFirstMate, String referenceSecondMate,
			boolean ignoreReadMappingsToDifferentGenes) 
	{
		String mappedTranscripts = "";
		
		for (String gene : genesMappedByBothReads)
		{
			HashMap<String, HashSet<String>> transcriptToMappableElements = new HashMap<String, HashSet<String>>();
			HashSet<String> mappableElements = new HashSet<String>();
			
			for (MappableElement m : genesToCoveredElementsFirstMate.get(gene))
			{
				mappableElements.add(m.getDescription());
								
				for (String transcript : m.getMappedTranscripts())
				{
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
				
				for (String transcript : m.getMappedTranscripts())
				{					
					if (!transcriptToMappableElements.containsKey(transcript))
					{
						transcriptToMappableElements.put(transcript, new HashSet<String>());
					}
					transcriptToMappableElements.get(transcript).add(m.getDescription());
				}
			}
			
			for (String transcript : transcriptToMappableElements.keySet())
			{
				if (transcriptToMappableElements.get(transcript).size() == mappableElements.size())
				{
					mappedTranscripts += transcript + "\n";
				}
			}
		}
		
		if (!mappedTranscripts.equals(""))
		{
			hitsFileWriter.println(">" + this.curReadID);
			hitsFileWriter.print(mappedTranscripts);
			numOfStoredMappings++;
		}		
	}

	public void storeProcessingResultsSingleEnd(
			HashMap<String, HashSet<MappableElement>> genesToCoveredElements,
			String reference) 
	{
		String mappedTranscripts = "";
		for (String gene : genesToCoveredElements.keySet())
		{		
			/* In case of junction read there are several mappable
			 * elements returned. The exonic regions before the
			 * junction, the exonic after the junction and the
			 * junction itsself. Thus, also here is has to be 
			 * counted to how many mappable elements the transcripts
			 * belong. 
			 */
			HashMap<String, HashSet<String>> transcriptToMappableElements = new HashMap<String, HashSet<String>>();
			HashSet<String> mappableElements = new HashSet<String>();
			
			for (MappableElement m : genesToCoveredElements.get(gene))
			{
				mappableElements.add(m.getDescription());
				
				for (String transcript : m.getMappedTranscripts())
				{
					if (!transcriptToMappableElements.containsKey(transcript))
					{
						transcriptToMappableElements.put(transcript, new HashSet<String>());
					}
					transcriptToMappableElements.get(transcript).add(m.getDescription());
				}
			}
			
			for (String transcript : transcriptToMappableElements.keySet())
			{
				if (transcriptToMappableElements.get(transcript).size() == mappableElements.size())
				{
					mappedTranscripts += transcript + "\n";
				}
			}
		}
		
		if (!mappedTranscripts.equals(""))
		{
			hitsFileWriter.println(">" + this.curReadID);
			hitsFileWriter.print(mappedTranscripts);
			numOfStoredMappings++;
		}		
	}
	
	@Override
	public int getNumOfStoredElements() 
	{
		return numOfStoredMappings;
	}

	@Override
	public void setReadID(String readID) 
	{
		this.curReadID = readID;
	}
	
	public void finishProcessing()
	{
		this.hitsFileWriter.flush();
		this.hitsFileWriter.close();
	}
}
