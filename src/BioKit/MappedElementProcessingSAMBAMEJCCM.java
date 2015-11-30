package BioKit;

import java.util.HashMap;
import java.util.HashSet;

public class MappedElementProcessingSAMBAMEJCCM extends MappedElementProcessingSAMBAMCompression 
{
	@Override
	public void storeProcessingResultsPairedEnd(
			HashSet<String> genesMappedByBothReads,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsFirstMate,
			HashMap<String, HashSet<MappableElement>> genesToCoveredElementsSecondMate,
			String referenceFirstMate, String referenceSecondMate,
			boolean ignoreReadMappingsToDifferentGenes) 
	{		
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
					
					storeFeatureInfo(combinationID, referenceFirstMate,
							minStartMaxStop.get(MappedElementProcessingSAMBAMCompression.INDEX_START), 
							minStartMaxStop.get(MappedElementProcessingSAMBAMCompression.INDEX_STOP));
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

				storeFeatureInfo(combinationID, referenceSecondMate,
						minStartMaxStop.get(MappedElementProcessingSAMBAMCompression.INDEX_START), 
						minStartMaxStop.get(MappedElementProcessingSAMBAMCompression.INDEX_STOP));
			}
		}
		else if (genesToCoveredElementsSecondMate.isEmpty())
		{			
			for (String geneFirstMate : genesToCoveredElementsFirstMate.keySet())
			{
				HashMap<Integer, Integer> minStartMaxStop = new HashMap<Integer, Integer>();
				String combinationID = generateCombinationID(genesToCoveredElementsFirstMate.get(geneFirstMate), 
						minStartMaxStop);

				storeFeatureInfo(combinationID, referenceFirstMate,
						minStartMaxStop.get(MappedElementProcessingSAMBAMCompression.INDEX_START), 
						minStartMaxStop.get(MappedElementProcessingSAMBAMCompression.INDEX_STOP));
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
	
				storeFeatureInfo(combinationID, referenceSecondMate, 
						minStartMaxStop.get(MappedElementProcessingSAMBAMCompression.INDEX_START), 
						minStartMaxStop.get(MappedElementProcessingSAMBAMCompression.INDEX_STOP));							
			}
		}
	}

	@Override
	public void storeProcessingResultsSingleEnd(
			HashMap<String, HashSet<MappableElement>> genesToCoveredElements,
			String reference) 
	{
		for (String gene : genesToCoveredElements.keySet())
		{
			HashMap<Integer, Integer> minStartMaxStop = new HashMap<Integer, Integer>();
			String combinationID = generateCombinationID(genesToCoveredElements.get(gene), 
					minStartMaxStop);

			storeFeatureInfo(combinationID, reference,
					minStartMaxStop.get(MappedElementProcessingSAMBAMCompression.INDEX_START), 
					minStartMaxStop.get(MappedElementProcessingSAMBAMCompression.INDEX_STOP));
		}
	}

	@Override
	public void finishProcessing() {
		// TODO Auto-generated method stub
	}

	@Override
	public void setReadID(String readID) {
		// TODO Auto-generated method stub
	}
}
