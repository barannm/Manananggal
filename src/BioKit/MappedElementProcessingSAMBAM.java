package BioKit;

import java.util.HashMap;
import java.util.HashSet;

public interface MappedElementProcessingSAMBAM 
{
	public void storeProcessingResultsPairedEnd(HashSet<String> genesMappedByBothReads,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsFirstMate,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsSecondMate,
			String referenceFirstMate, String referenceSecondMate,
			boolean ignoreReadMappingsToDifferentGenes);
	
	public void storeProcessingResultsSingleEnd(HashMap<String, HashSet<MappableElement>> genesToCoveredElements,
			String reference);
	
	public int getNumOfStoredElements();
	
	public void setReadID(String readID);
	
	public void finishProcessing();
}
