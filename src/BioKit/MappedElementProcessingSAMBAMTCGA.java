/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author Fabian Birzele
 */
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
