package BioKit;

import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.Vector;

public class SlimExonicRegionGene implements Comparable<SlimExonicRegionGene>, Serializable 
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	// start position of the gene
	private int start;
	
	// end position of the gene
	private int stop;
	
	// chromosome of the gene
	private String chromosome;
	
	private String name;
	
	// strand
	private boolean plusStrand;
	
	private Integer geneID;
	
	private boolean setCorrespondingExonsForExonicRegion;
	private boolean setCorrespondingTranscriptsForExonicRegion;
	
	private Vector<Exon> exonsSorted;
	private GeneticRegion[] exonicRegions; 
	private HashMap<String, Exon> exonIDToExon;
		
	public SlimExonicRegionGene(String name, Integer geneID, String chromosome,
			boolean plusStrand)
	{
		this.exonsSorted = new Vector<Exon>();
		this.exonIDToExon = new HashMap<String, Exon>();
		this.exonicRegions = null;
		this.chromosome = chromosome;
		this.plusStrand = plusStrand;
		this.start = -1;
		this.stop = -1;
		this.name = name;
		this.geneID = geneID;
		this.setCorrespondingExonsForExonicRegion = true;
		this.setCorrespondingTranscriptsForExonicRegion = true;
	}

	public void addExon(Exon e)
	{
		/* Exons come in Sorted therefore the
		 * beginning of the first exon is the
		 * beginning of the genes coding region.
		 */
		if (exonsSorted.isEmpty())
		{
			this.start = e.getCodingStart();
			this.stop = e.getCodingStop();
		}
		else if (this.stop < e.getCodingStop())
		{
			this.stop = e.getCodingStop();
		}
		
		String exonID = e.getCodingStart() + "-" + e.getCodingStop();
		
		this.exonsSorted.add(e);
		this.exonIDToExon.put(exonID, e);
	}
		
	public GeneticRegion[] getExonicRegions()
	{
		if (this.exonicRegions == null) 
		{
			ExonGroup[] exonGroups = computeOverlappingExonGroups();
			int exonicRegionCounter = 0;
			Vector<GeneticRegion> exonicRegions = new Vector<GeneticRegion>();
			for(ExonGroup g : exonGroups)
			{
				HashMap<Integer, HashSet<String>> positionToExonID = new HashMap<Integer, HashSet<String>>();
				HashSet<Integer> positions = new HashSet<Integer>();
				// There exist transcripts annotated with exons containing only one
				// base. These exons disappear when adding the positions to the
				// upper hash set and cause errors. In order to circumvent that
				// these positions are memorized in the hash below an are considered
				// twice later.
				HashSet<Integer> positionsToConsiderTwice = new HashSet<Integer>();
				Vector<Exon> exonsOfGroup = g.getExons();
				
				for (Exon exon : exonsOfGroup)
				{
					positions.add(exon.getGenomicStart());
					positions.add(exon.getGenomicStop());
					// Dealing with exon containing only one base.
					if (exon.getGenomicStart() == exon.getGenomicStop())
						positionsToConsiderTwice.add(exon.getGenomicStart());
					
					for (int pos = exon.getGenomicStart(); pos <= exon.getGenomicStop(); pos++)
					{
						if (!positionToExonID.containsKey(pos))
						{
							positionToExonID.put(pos, new HashSet<String>());
						}
						positionToExonID.get(pos).add(exon.getGenomicStart() + "-" + exon.getGenomicStop());
					}
				}			
				
				Integer[] sortedPositions = new Integer[positions.size() + positionsToConsiderTwice.size()];
				int i = 0;
				for (Integer position : positions)	sortedPositions[i++] = position;
				for (Integer position : positionsToConsiderTwice)	sortedPositions[i++] = position;
				Arrays.sort(sortedPositions);

				//DEBUG
//				if (this.getName().equals("ENSG00000187634_1"))
//				{
//					System.out.println("ENSG00000187634_1\tSorted positions: " + sortedPositions);
//					System.out.println("ENSG00000187634_1\tExonic regions: ");
//				}
				//DEBUG
				
				int startCurrentGeneticRegion = sortedPositions[0];
				for (i = 1; i < sortedPositions.length; i++)
				{
					GeneticRegion exonicRegion = new GeneticRegion(startCurrentGeneticRegion, sortedPositions[i]);
					exonicRegion.setElementID(exonicRegionCounter);
					exonicRegion.setReference(getChromosome());
					exonicRegion.setID(exonicRegionCounter);
					exonicRegions.add(exonicRegion);
					startCurrentGeneticRegion = sortedPositions[i];
					//DEBUG
//					if (this.getName().equals("ENSG00000187634_1"))
//					{
//						System.out.println("\nENSG00000187634_1\t" + exonicRegion.getCodingStart() + "-" + exonicRegion.getCodingStop());
//					}
//					System.out.println("\nDEBUG:\tCurrent exonic region: " + exonicRegion.getDescription());
					//END DEBUG
					
					/* Add mapped exons and mapped transcripts to exonic region.
					 * Therefore it is checked which exons overlap with the start and
					 * stop of the current exonic region and only those which overlap
					 * both are taken.
					 */
					if (this.setCorrespondingExonsForExonicRegion || this.setCorrespondingTranscriptsForExonicRegion)
					{
						//DEBUG
//						System.out.println("DEBUG:\tCheck exons: ");
						//END DEBUG
						/* Find out overlapping exons with the current exonic region.
						 * Since the start and stop positions overlap by one position
						 * they can be covered by two exons simultanously whereas
						 * for one exon it is the last position. Therefore the exonIDs
						 * stored for a position within the exonic region is taken
						 * for determining the covering exons.
						 */
//						if (this.getName().equals("ENSG00000187634_1"))
//						{
//							System.out.println("ENSG00000187634_1\tAdding exons to exonic region " + exonicRegion.getDescription());
//						}
						int inbetweenPos = exonicRegion.getCodingStart() + ((exonicRegion.getCodingStop() - exonicRegion.getCodingStart())/2);
						for (String exonID : positionToExonID.get(inbetweenPos))
						{
							//DEBUG

//							System.out.println("DEBUG:\tExon overlapping start " + exonicRegion.getCodingStart() + ": " + exonID);
//							if (!positionToExonID.get(exonicRegion.getCodingStop()).contains(exonID))
//							{
//								System.out.println("DEBUG:\tExon does not overlap stop " + exonicRegion.getCodingStop() + ": " + exonID);
//							}
//							else
//							{
//								System.out.println("DEBUG:\tExon also overlaps stop " + exonicRegion.getCodingStop() + ": " + exonID);
//							}
//							//END DEBUG
							if (!positionToExonID.get(exonicRegion.getCodingStart()).contains(exonID)) continue;
							
							Exon currentExon = this.exonIDToExon.get(exonID);
							if (this.setCorrespondingExonsForExonicRegion)
							{
								exonicRegion.addMappedExon(currentExon);
//								if (this.getName().equals("ENSG00000187634_1"))
//								{
//									System.out.println("ENSG00000187634_1\tCurrent exon: " + exonID);
//								}
							}
//							//DEBUG
//							System.out.print("DEBUG:\tExons of region now:");
//							for (Exon exon : exonicRegion.getMappedExons())
//							{
//								System.out.print(" " + exon.getDescription());
//							}
//							System.out.println();
//							//END DEBUG
							if (this.setCorrespondingTranscriptsForExonicRegion)
							{
								//DEBUG
//								System.out.println("DEBUG:\tTRANSCRIPTS:");
//								if (this.getName().equals("ENSG00000187634_1"))
//								{
//									System.out.println("ENSG00000187634_1\tAdding transcripts to exonic region from exon " + exonID);
//								}
								//END DEBUG
								for (String transcript : currentExon.getMappedTranscripts())
								{
//									if (this.getName().equals("ENSG00000187634_1"))
//									{
//										System.out.println("ENSG00000187634_1\t" + transcript);
//									}

									exonicRegion.addMappedTranscript(transcript);
									//DEBUG
//									System.out.println("DEBUG:\t" + transcript);
									//END DEBUG
								}
								//DEBUG
//								System.out.print("DEBUG:\tTranscripts of region now:");
//								for (String transcript : exonicRegion.getMappedTranscripts())
//								{
//									System.out.print(" " + transcript);
//								}
//								System.out.println();
								//END DEBUG
							}
						}		
					}
					exonicRegionCounter++;
				}
			}
			this.exonicRegions = exonicRegions.toArray(new GeneticRegion[exonicRegions.size()]); 
		}
		return this.exonicRegions;
	}
	
	private ExonGroup[] computeOverlappingExonGroups()
	{
		TreeSet<ExonGroup> exonGroups = new TreeSet<ExonGroup>();
		
		Exon[] sortedExons = new Exon[exonsSorted.size()];
		sortedExons = this.exonsSorted.toArray(sortedExons);
		
		TreeSet<NumericSorterContainer<Exon>> exons = new TreeSet<NumericSorterContainer<Exon>>();
		
		for(Exon e : sortedExons)
		{
			exons.add(new NumericSorterContainer<Exon>(e.getCodingStop(), e, true));
		}
		
		int groupCounter = 0;
		
		int counter = 0;
		
		for(NumericSorterContainer<Exon> c : exons.descendingSet())
		{
			Exon exon = c.getContent();
			//System.out.println("Process Exon " + exon.getGenomicLength() + " " + exon.getGenomicStart() + "-" + exon.getGenomicStop());
			
			boolean added = false;
			
			// check if exon intersects with known other exon
			for(ExonGroup group : exonGroups)
			{
				if(group.exonIntersectsWithExonsInGroup(exon))
				{
					// if exon intersects
					
					if(!added)
					{
						group.addExon(counter, exon.getExonID(),  exon);
						
						//System.out.println("Added to group: " + group.toShortString());
						
						added = true;
					}
					else if(added)
					{
						//System.out.println("Exon is already added and fits to another group! " + group.toShortString());
					}
				}
			}
			
			// create new group if not added so far
			if(!added)
			{
				ExonGroup group = new ExonGroup(groupCounter, isPlusStrand());
				group.addExon(counter, exon.getExonID(), exon);
				
				//System.out.println("New group: " + group.toShortString());
				
				exonGroups.add(group);
				
				groupCounter++;
			}
			
			counter++;
		}
		
		return exonGroups.toArray(new ExonGroup[exonGroups.size()]);
	}
	
	public void setSetCorrespondingExonsForExonicRegion(
			boolean setCorrespondingExonsForExonicRegion) 
	{
		this.setCorrespondingExonsForExonicRegion = setCorrespondingExonsForExonicRegion;
	}

	public void setSetCorrespondingTranscriptsForExonicRegion(
			boolean setCorrespondingTranscriptsForExonicRegion) 
	{
		this.setCorrespondingTranscriptsForExonicRegion = setCorrespondingTranscriptsForExonicRegion;
	}

	public int getStart() 
	{
		return start;
	}

	public int getStop() 
	{
		return stop;
	}

	public String getChromosome() 
	{
		return chromosome;
	}

	public Vector<Exon> getExons() 
	{
		return this.exonsSorted;
	}

	public boolean isPlusStrand()
	{
		return plusStrand;
	}
	
	public Integer getGeneID()
	{
		return this.geneID;
	}
	
	public String getName()
	{
		return this.name;
	}
	
	@Override
	public int compareTo(SlimExonicRegionGene otherGene) 
	{
		return this.geneID.compareTo(otherGene.getGeneID());
	}
}
