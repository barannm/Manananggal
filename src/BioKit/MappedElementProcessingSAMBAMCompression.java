package BioKit;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public abstract class MappedElementProcessingSAMBAMCompression implements MappedElementProcessingSAMBAM 
{
	private HashMap<String, Integer[]> featureCombinationToStartStopAndOccurrence;
	private HashMap<String, String> featureCombinationToSeqName;
	public static final int INDEX_START = 0;
	public static final int INDEX_STOP = 1;
	public static final int INDEX_OCCURRENCE = 2;
	
	public MappedElementProcessingSAMBAMCompression()
	{
		this.featureCombinationToStartStopAndOccurrence = new HashMap<String, Integer[]>();
		this.featureCombinationToSeqName = new HashMap<String, String>();
	}
	
	public HashMap<String, Integer[]> getCombinationIDsStartTopAndOccurrenceData()
	{
		return this.featureCombinationToStartStopAndOccurrence;
	}
	
	public HashMap<String, String> getCombinationIDSeqNames()
	{
		return this.featureCombinationToSeqName;
	}
	
	protected String generateCombinationID(HashSet<MappableElement> coveredElements, 
			HashMap<Integer, Integer> minStartMaxStop)
	{
		String combinationID = null;		
		
		HashMap<String, MappableElement> descriptionToMappableElement = new HashMap<String, MappableElement>();
		for (MappableElement curElement : coveredElements)
		{
			descriptionToMappableElement.put(curElement.getDescription(), curElement);
		}
		
		String[] descriptions = new String[descriptionToMappableElement.keySet().size()];
		descriptions = descriptionToMappableElement.keySet().toArray(descriptions);
		Arrays.sort(descriptions);
		
		for (String description : descriptions)
		{
			MappableElement curElement = descriptionToMappableElement.get(description);
			String curElementID = curElement.getDescription();
			if (combinationID == null) combinationID = curElementID;
			else combinationID += SAMBAMCompression.ELEMENT_ID_SEPARATOR + curElementID;
			if (minStartMaxStop.get(MappedElementProcessingSAMBAMCompression.INDEX_START) == null 
					|| curElement.getCodingStart() < minStartMaxStop.get(MappedElementProcessingSAMBAMCompression.INDEX_START))
				minStartMaxStop.put(MappedElementProcessingSAMBAMCompression.INDEX_START, curElement.getCodingStart());
			if (minStartMaxStop.get(MappedElementProcessingSAMBAMCompression.INDEX_STOP) == null
					|| curElement.getCodingStop() > minStartMaxStop.get(MappedElementProcessingSAMBAMCompression.INDEX_STOP))
				minStartMaxStop.put(MappedElementProcessingSAMBAMCompression.INDEX_STOP, curElement.getCodingStop());
		}
				
		return combinationID;
	}
	
	protected void storeFeatureInfo(String featureID, String reference, 
			int start, int stop)
	{
		if (this.featureCombinationToStartStopAndOccurrence.containsKey(featureID))
		{
			this.featureCombinationToStartStopAndOccurrence.get(featureID)[MappedElementProcessingSAMBAMCompression.INDEX_OCCURRENCE]++;
		}
		else
		{
			Integer[] combInfo = new Integer[3];
			combInfo[MappedElementProcessingSAMBAMCompression.INDEX_START] = start;
			combInfo[MappedElementProcessingSAMBAMCompression.INDEX_STOP] = stop;
			combInfo[MappedElementProcessingSAMBAMCompression.INDEX_OCCURRENCE] = 1;
			this.featureCombinationToStartStopAndOccurrence.put(featureID, combInfo);
			this.featureCombinationToSeqName.put(featureID, reference);
		}
	}
	
	@Override
	public int getNumOfStoredElements() 
	{
		return this.featureCombinationToStartStopAndOccurrence.keySet().size();
	}
}
