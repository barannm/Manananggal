package BioKit;

import java.util.HashMap;
import java.util.HashSet;

public class MappedElementProcessingSAMBAMTCGA extends MappedElementProcessingSAMBAMCompression 
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
		if (genesMappedByBothReads.isEmpty() && !ignoreReadMappingsToDifferentGenes)
		{
			for (String gene : genesToCoveredElementsFirstMate.keySet())
			{
				for (MappableElement coveredElement : genesToCoveredElementsFirstMate.get(gene))
				{
					storeCoveredElement(coveredElement);
				}
			}
			
			for (String gene : genesToCoveredElementsSecondMate.keySet())
			{
				for (MappableElement coveredElement : genesToCoveredElementsSecondMate.get(gene))
				{
					storeCoveredElement(coveredElement);
				}
			}							
		}
		else
		{
			for (String gene : genesMappedByBothReads)
			{
				for (MappableElement coveredElement : genesToCoveredElementsFirstMate.get(gene))
				{
					storeCoveredElement(coveredElement);
				} 
				
				for (MappableElement coveredElement : genesToCoveredElementsSecondMate.get(gene))
				{
					storeCoveredElement(coveredElement);
				}
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
			for (MappableElement coveredElement : genesToCoveredElements.get(gene))
			{
				storeCoveredElement(coveredElement);
			}		
		}
	}
	
	private void storeCoveredElement(MappableElement coveredElement)
	{
		for (Exon e : coveredElement.getMappedExons())
		{
			storeFeatureInfo(e.getDescription(), e.getReference(), 
				e.getCodingStart(), 
				e.getCodingStop());
		}
		
		if (coveredElement instanceof Junction)
		{
			storeFeatureInfo(coveredElement.getDescription(), coveredElement.getReference(), 
					coveredElement.getCodingStart(), 
					coveredElement.getCodingStop());
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
