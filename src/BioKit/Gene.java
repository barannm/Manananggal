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

import java.io.Serializable;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

/**
 * This class represents a gene in a genome. A gene has a unique location on the genome specified
 * by its position, its strand and its chromosome. It consists of a set of exons. Exons can be combined
 * in various ways leading to a set of possible gene products, i.e. proteins. 
 *  
 * @author Fabian Birzele
 *
 */
public class Gene implements Comparable<Gene>, Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	/**
	 * Name for the Exon property
	 */
	public static final String PROPERTY_EXONS = "PROPERTY_EXONS";
	
	// start position of the gene
	private int start;
	
	// end position of the gene
	private int stop;
	
	// chromosome of the gene
	private String chromosome;
	
	// strand
	private boolean plusStrand;
	
	// unique gene id
	private String geneID;
	
	// gene symbol
	private String geneName;
	
	// gene type, i.e. miRNA, protein_coding, etc.
	private String geneType;
	
	// set of exons comprising the gene
	private TreeSet<Exon> exons;
	private TreeSet<Intron> introns;
	
	// set of gene products that are known for this gene
	private HashMap<String, Protein> geneProducts;
	
	// allows to add properties to the gene
	private HashMap<String, Object> properties;
	
	private TreeMap<String, Integer> mapCDSstartToIsoforms;
	private TreeMap<String, Integer> mapCDSendToIsoforms;
	
	public Gene(String geneID, int start, int stop, String chromosome, boolean plusStrand)
	{
		this.geneID = geneID;
		this.geneName = null;
		this.geneType = null;
		this.start = start;
		this.stop = stop;
		this.chromosome = chromosome;
		this.plusStrand = plusStrand;
		
		geneProducts = new HashMap<String, Protein>();
		exons = new TreeSet<Exon>();
		introns = new TreeSet<Intron>();
	}
	
	public Gene(String geneID, String geneName, String geneType, int start, int stop, String chromosome, boolean plusStrand)
	{
		this.geneID = geneID; 
		this.geneName = geneName;
		this.geneType = geneType;
		this.start = start;
		this.stop = stop;
		this.chromosome = chromosome;
		this.plusStrand = plusStrand;		
		
		geneProducts = new HashMap<String, Protein>();
		exons = new TreeSet<Exon>();
		introns = new TreeSet<Intron>();
	}
	
	public boolean ParseFromRefFlat(String[] pColumns)
	{
		if (pColumns.length != 16)
		{
			System.out.println("Error: invalid refFlat line -> requires 16 columns.");
			return false;
		}
		
		this.geneID = pColumns[1];
		this.geneName = pColumns[12];
		
		this.chromosome = pColumns[2];
		this.start = (Integer.parseInt(pColumns[4]) + 1);
		this.stop = Integer.parseInt(pColumns[5]);
		
		if (pColumns[3].equals("+"))
		  this.plusStrand = true;
		else
		  this.plusStrand = false;

		int nExons = Integer.parseInt(pColumns[8]);
		String[] pExonStarts = pColumns[9].split(",");
		String[] pExonEnds = pColumns[10].split(",");
	    
		for (int i = 0; i < nExons; i++)
	    {
	      Exon exon = new Exon(Integer.parseInt(pExonStarts[i]) + 1, Integer.parseInt(pExonEnds[i]));
	      addExon(exon);
	    }
	    return true;
	}
	
	public void SetCDSPositions(TreeMap<String, Integer> mapCDSstart, TreeMap<String, Integer> mapCDSend)
	{
		mapCDSstartToIsoforms = mapCDSstart;
		mapCDSendToIsoforms	  = mapCDSend;
	}
	
	public int GetCodingStartForIsoform(String strIsoform)
	{
		if(mapCDSstartToIsoforms.containsKey(strIsoform))
			return mapCDSstartToIsoforms.get(strIsoform);
		else
		{
			if(isPlusStrand())
				return Integer.MIN_VALUE;
			else
				return Integer.MAX_VALUE;
		}			
	}
	
	public int GetCodingEndForIsoform(String strIsoform)
	{
		if(mapCDSendToIsoforms.containsKey(strIsoform))
			return mapCDSendToIsoforms.get(strIsoform);
		else
		{
			if(isPlusStrand())
				return Integer.MAX_VALUE;
			else
				return Integer.MIN_VALUE;
		}	
	}
	
	/**
	 * Method allows to add properties to the gene
	 * 
	 * @param id
	 * @param property
	 */
	public void addProperty(String id, Object property)
	{
		if(properties == null)
			properties = new HashMap<String, Object>();
		
		properties.put(id, property);
	}
	
	/**
	 * Method returns the object stored for a specified property name or null if the property does not exist
	 * 
	 * @param id
	 * @return
	 */
	public Object getProperty(String id)
	{
		if(properties == null)
			return null;
		
		return properties.get(id);
	}
	
	/**
	 * Method adds an exon to the gene. The exon will automatically be sorted into the row of existing
	 * exons with respect to its compareTo method. Start and stop positions of the gene are updated
	 * acordingly (if necessary)
	 * 
	 * @param exon
	 */
	public void addExon(Exon exon)
	{
		exons.add(exon);
	
		if(exon.getGenomicStart() < start)
			start = exon.getGenomicStart();
		
		if(exon.getGenomicStop() > stop)
			stop = exon.getGenomicStop();
	}
	
	/**
	 * Removes the specified exon from the gene
	 * 
	 * @param exon
	 */
	public void removeExon(Exon exon)
	{
		exons.remove(exon);
		
		// update start and end positions
		if(exons.size() > 0)
		{
			start = exons.first().getGenomicStart();
			stop = exons.last().getGenomicStop();
		}
	}
	
	/**
	 * Method returns an array containing all exons as sorted by their compareTo method. Usually this
	 * should be their position in the gene, represented by the exons chromosomal position in the plus
	 * strand.
	 * 
	 * @return
	 */
	public Exon[] getSortedExons()
	{
		return (Exon[])exons.toArray(new Exon[exons.size()]);
	}
	
	/**
	 * Method adds an intron to the gene. The intron will automatically be sorted into the row of existing
	 * introns with respect to its compareTo method. 
	 * 
	 * @param exon
	 */
	public void addIntron(Intron intron)
	{
		introns.add(intron);
	}
	
	/**
	 * Removes the specified intron from the gene
	 * 
	 * @param exon
	 */
	public void removeIntron(Intron intron)
	{
		introns.remove(intron);
	}
	
	/**
	 * Method returns an array containing all introns as sorted by their compareTo method. Usually this
	 * should be their position in the gene, represented by the Introns chromosomal position in the plus
	 * strand.
	 * 
	 * @return
	 */
	public Intron[] getSortedIntrons()
	{
		return (Intron[])introns.toArray(new Intron[introns.size()]);
	}
	
	/**
	 * Method returns a sorted array containing all introns and exons of the gene sorted by their compartTo method.
	 * 
	 * @return
	 */
	public GeneticElement[] getSortedGeneticElements()
	{
		TreeSet<GeneticElement> elements = new TreeSet<GeneticElement>();
		elements.addAll(exons);
		elements.addAll(introns);
		
		return (GeneticElement[])elements.toArray(new GeneticElement[elements.size()]);
	}
	
	/**
	 * Method returns an array of exon groups where each group represents one genomic locus
	 * which may result in several, overlapping cassette exons.
	 * 
	 * @param exons
	 * @return
	 */
	public ExonGroup[] computeOverlappingExonGroups()
	{
		TreeSet<ExonGroup> exonGroups = new TreeSet<ExonGroup>();
		
		Exon[] sortedExons = getSortedExons();
		
		TreeSet<NumericSorterContainer<Exon>> exons = new TreeSet<NumericSorterContainer<Exon>>();
		
		for(Exon e : sortedExons)
		{
//			exons.add(new NumericSorterContainer<Exon>(e.getGenomicLength(), e, true));
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
						break; // just break the loop here, the added=true case is not executed anyway
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
	
	/**
	 * Similar to above, the method returns an array of exon groups where each group represents one genomic locus
	 * which may result in several, overlapping cassette exons.
	 * 
	 * HOWEVER: This function does not reverse the exon order.
	 * 
	 * @param exons
	 * @return
	 */
	public ExonGroup[] computeOverlappingExonGroupsNormalOrder()
	{
		TreeSet<ExonGroup> exonGroups = new TreeSet<ExonGroup>();
		
		Exon[] sortedExons = getSortedExons();
		
		int groupCounter = 0;
		
		int counter = 0;
		
		for(Exon exon : sortedExons)
		{
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
						break; // just break the loop here, the added=true case is not executed anyway
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
	
	/**
	 * Method computes the composite gene length, so the length of all exonic positions of the gene in any of the transcript models.
	 * @return
	 */
	public int computeCompositeGeneLength()
	{
		ExonGroup[] groups = computeOverlappingExonGroups();
		
		int sum = 0;
		
		for(ExonGroup g : groups)
			sum += g.getExonGroupLengthInBp();
		
		return sum;
	}
	
	/**
	 * This method returns a TreeSet which contains all splice junctions found in this gene.
	 * Each junction is represented as from-to string.
	 *
	 * @return
	 */
	public TreeSet<String> getSpliceJunctionInformation()
	{
		TreeSet<String> junctions = new TreeSet<String>();
		
		for(String transcript : geneProducts.keySet())
		{
			// get sorted exons, reverse order if minus strand gene
			Exon[] sorted = getSortedExonsForGeneProduct(transcript);
			
			if(!isPlusStrand())
			{
				Exon[] temp = new Exon[sorted.length];
				
				int counter = 0;
				
				for(int k=sorted.length-1; k>=0; k--)
				{
					temp[counter] = sorted[k];
					counter++;
				}
				
				sorted = temp;
			}
			
			Exon lastExon = null;
			
			for(Exon currentExon : sorted)
			{
				// draw junction
				if(lastExon != null)
				{
					String junctionID = lastExon.getGenomicStop() + "-" +currentExon.getGenomicStart();
					
					if(!isPlusStrand())
						junctionID = currentExon.getGenomicStop() + "-" + lastExon.getGenomicStart();
						
					junctions.add(junctionID);
				}
				
				lastExon = currentExon;
			}
		}
		
		return junctions;
	}
	
	public TreeSet<Junction> getSpliceJunctionInformationForGeneProductAsJunctions(String strGeneProductID)
	{
		TreeSet<Junction> junctions = new TreeSet<Junction>();
		
		// get sorted exons, reverse order if minus strand gene
		Exon[] sorted = getSortedExonsForGeneProduct(strGeneProductID);
		
		if(!isPlusStrand())
		{
			Exon[] temp = new Exon[sorted.length];
			
			int counter = 0;
			
			for(int k=sorted.length-1; k>=0; k--)
			{
				temp[counter] = sorted[k];
				counter++;
			}
			
			sorted = temp;
		}
		
		Exon lastExon = null;
		
		for(Exon currentExon : sorted)
		{
			// draw junction
			if(lastExon != null)
			{
				if(isPlusStrand())
				{
					Junction jun = new Junction(chromosome, lastExon.getGenomicStop(), currentExon.getGenomicStart());
					junctions.add(jun);
				}
				else
				{
					Junction jun = new Junction(chromosome, currentExon.getGenomicStop(), lastExon.getGenomicStart());
					junctions.add(jun);
				}
			}
			
			lastExon = currentExon;
		}
		
		return junctions;
	}
	
	public TreeSet<String> getSpliceJunctionInformationForGeneProduct(String strGeneProductID)
	{
		TreeSet<String> junctions = new TreeSet<String>();
		
		// get sorted exons, reverse order if minus strand gene
		Exon[] sorted = getSortedExonsForGeneProduct(strGeneProductID);
		
		if(!isPlusStrand())
		{
			Exon[] temp = new Exon[sorted.length];
			
			int counter = 0;
			
			for(int k=sorted.length-1; k>=0; k--)
			{
				temp[counter] = sorted[k];
				counter++;
			}
			
			sorted = temp;
		}
		
		Exon lastExon = null;
		
		for(Exon currentExon : sorted)
		{
			// draw junction
			if(lastExon != null)
			{
				String junctionID = lastExon.getGenomicStop() + "-" +currentExon.getGenomicStart();
				
				if(!isPlusStrand())
					junctionID = currentExon.getGenomicStop() + "-" + lastExon.getGenomicStart();
					
				junctions.add(junctionID);
			}
			
			lastExon = currentExon;
		}
		
		return junctions;
	}
	
	/**
	 * get unique junctions for gene (i.e. junctions that occur in one isoform only).
	 * return null if no unique junctions were identified. You may restrict the search
	 * for unique exons by providing a list of isoform (gene product) names.
	 */
	public TreeSet<String> getUniqueJunctions(TreeSet<String> vcIsoformsToConsider)
	{
		TreeSet<String> res = null;
		
		// count how many times each junction occurs in all isoforms
		TreeMap<String, Integer> mapJunctionCounts = new TreeMap<String, Integer>();	// key = junction position, value = count
		
		for(String strIsoform : geneProducts.keySet())
		{
			// if a list of isoforms was provided, use only isoforms contained in that list
			if(vcIsoformsToConsider != null && !vcIsoformsToConsider.contains(strIsoform))
				continue;
			
			TreeSet<String> vcJunctions = getSpliceJunctionInformationForGeneProduct(strIsoform);
			
			for(String strJunction : vcJunctions)
			{
				if(mapJunctionCounts.containsKey(strJunction))
				{
					int nVal = mapJunctionCounts.get(strJunction);
					mapJunctionCounts.put(strJunction, nVal+1);
				}
				else
				{
					mapJunctionCounts.put(strJunction, 1);
				}
			}
		}
		
		for(Map.Entry<String, Integer> e : mapJunctionCounts.entrySet())
		{
			if(e.getValue() == 1)
			{
				if(res == null) res =  new TreeSet<String>();
				
				res.add(e.getKey());
			}
		}
		
		return res;
	}
	
	/**
	* get unique junctions of the specified gene product, return null if there are no unique features.
	* You may restrict the search for unique exons by providing a list of isoform (gene product) names.
	*/
	public TreeSet<String> getUniqueJunctionsForGeneProduct(String strGeneProductID, TreeSet<String> vcIsoformsToConsider)
	{
		TreeSet<String> vcIsoformSpecificJunctions = null;
		
		// if a list of isoforms is considered, the current isoform must be part of it
		if(vcIsoformsToConsider != null && !vcIsoformsToConsider.contains(strGeneProductID))
			return null;
		
		// get all unique junctions
		TreeSet<String> vcUniqueJunctions = getUniqueJunctions(vcIsoformsToConsider);
		
		if(vcUniqueJunctions == null)
			return null;
		
		// get current isoform junctions
		TreeSet<String> vcJunctions = getSpliceJunctionInformationForGeneProduct(strGeneProductID);
		
		for(String strJunction : vcJunctions)
		{
			if(vcUniqueJunctions.contains(strJunction))
			{
				if(vcIsoformSpecificJunctions == null)
				{
					vcIsoformSpecificJunctions = new TreeSet<String>();
					vcIsoformSpecificJunctions.add(strJunction);
				}
			}
		}

		return vcIsoformSpecificJunctions;
	}
	
	/**
	* get unique exons for the gene (i.e. exons that occur in one isoform only).
	* return null if no unique exons were identified.  You may restrict the search
	* for unique exons by providing a list of isoform (gene product) names.
	*/
	public TreeSet<Exon> getUniqueExons(TreeSet<String> vcIsoformsToConsider)
	{
		TreeSet<Exon> vcUniqueExons = null;
		
		TreeMap<Exon, Integer> mapExonCounts = new TreeMap<Exon, Integer>(); // key = exon, value = count 
		
		for(String strGeneProductID : this.geneProducts.keySet())
		{
			// if a list of isoforms was provided, use only isoforms contained in that list
			if(vcIsoformsToConsider != null && !vcIsoformsToConsider.contains(strGeneProductID))
				continue;
			
			Exon pExons[] = getSortedExonsForGeneProduct(strGeneProductID);
			
			for(Exon ex : pExons)
			{
				if(mapExonCounts.containsKey(ex))
				{
					int nVal = mapExonCounts.get(ex);
					mapExonCounts.put(ex, nVal+1);
				}
				else
				{
					mapExonCounts.put(ex, 1);
				}
			}
		}
		
		for(Map.Entry<Exon, Integer> e : mapExonCounts.entrySet())
		{
			if(e.getValue() == 1)
			{
				if(vcUniqueExons == null)
					vcUniqueExons = new TreeSet<Exon>();
				
				vcUniqueExons.add(e.getKey());
			}
		}
		
		return vcUniqueExons;
	}
	
	/**
	* get unique exons of the specified gene product, return null if there are no unique exons.
	 */
	public TreeSet<Exon> getUniqueExonsForGeneProduct(String strGeneProductID, TreeSet<String> vcIsoformsToConsider)
	{
		// if a list of isoforms is considered, the current isoform must be part of it
		if(vcIsoformsToConsider != null && !vcIsoformsToConsider.contains(strGeneProductID))
			return null;
				
		TreeSet<Exon> vcIsoformSpecificExons = null;
		
		// get unique exons
		TreeSet<Exon> vcUniqueExons = getUniqueExons(vcIsoformsToConsider);
		
		if(vcUniqueExons == null)
			return null;
		
		// get current isoform exons
		Exon pExons[] = getSortedExonsForGeneProduct(strGeneProductID);
		
		for(Exon ex : pExons)
		{
			if(vcUniqueExons.contains(ex))
			{
				if(vcIsoformSpecificExons == null)
					vcIsoformSpecificExons = new TreeSet<Exon>();
				
				vcIsoformSpecificExons.add(ex);
			}
		}
		
		return vcIsoformSpecificExons;
	}
	
	/**
	 * Method returns an array of intron groups where each group represents one genomic locus
	 * which may result in several, overlapping cassette exons.
	 * 
	 * @param exons
	 * @return
	 */
	public IntronGroup[] computeOverlappingIntronGroups()
	{
		TreeSet<IntronGroup> intronGroups = new TreeSet<IntronGroup>();
		
		Intron[] sortedIntrons = getSortedIntrons();
		
		TreeSet<NumericSorterContainer<Intron>> introns = new TreeSet<NumericSorterContainer<Intron>>();
		
		for(Intron i : sortedIntrons)
		{
			introns.add(new NumericSorterContainer<Intron>(i.getGenomicLength(), i, true));
		}
		
		int groupCounter = 0;
		
		int counter = 0;
		
		for(NumericSorterContainer<Intron> c : introns.descendingSet())
		{
			Intron intron = c.getContent();
			
			//System.out.println("Process Exon " + exon.getGenomicLength() + " " + exon.getGenomicStart() + "-" + exon.getGenomicStop());
			
			boolean added = false;
			
			// check if exon intersects with known other exon
			for(IntronGroup group : intronGroups)
			{
				if(group.intronIntersectsWithIntronsInGroup(intron))
				{
					// if exon intersects
					if(!added)
					{
						group.addIntron(counter, intron.getID(), intron);
						
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
				IntronGroup group = new IntronGroup(groupCounter, isPlusStrand());
				group.addIntron(counter, intron.getID(), intron);
				
				//System.out.println("New group: " + group.toShortString());
				
				intronGroups.add(group);
				
				groupCounter++;
			}
			
			counter++;
		}
		
		return intronGroups.toArray(new IntronGroup[intronGroups.size()]);
	}
	
	/**
	 * This method maps exons to their corresponding groups given both arrays 
	 * 
	 * @param groups
	 * @param exons
	 * @return
	 */
	public HashMap<Integer, Integer> mapExonsToGroups(ExonGroup[] groups, Exon[] exons)
	{
		HashMap<Integer, Integer> map = new HashMap<Integer, Integer>();
		
		for(int i=0; i<exons.length; i++)
		{
			for(ExonGroup group : groups)
			{
				if(group.groupContainsExon(i))
				{
					map.put(i, group.getGroupID());
					break;
				}
			}
		}
		
		return map;
	}
	
	/**
	 * Method returns a hashmap mapping each exon onto its position in the gene
	 * 
	 * @return
	 */
	public HashMap<Exon, Integer> getExonPositionsInSortedGeneExonArray()
	{
		Exon[] sortedExons = (Exon[])exons.toArray(new Exon[exons.size()]);
		
		HashMap<Exon, Integer> toReturn = new HashMap<Exon, Integer>();
		
		for(int i=0; i<sortedExons.length; i++)
		{
			toReturn.put(sortedExons[i], i);
		}
		
		return toReturn;
	}
	
	/**
	 * Method returns an iterator instance to iterate over all exons
	 * 
	 * @return
	 */
	public Iterator<Exon> getExonIterator()
	{
		return exons.iterator();
	}
	
	/**
	* Method returns an iterator instance to iterate over all introns
	*/
	public Iterator<Intron> getIntronIterator()
	{
		return introns.iterator();
	}
	
	/**
	 * Method adds a product of the gene to the gene. All exons required for the gene product are added to the 
	 * genes exon list. Every exon is only added once for the gene. Gene boundaries are updated accordingly
	 * with each protein exon, if new start or stop position is found.
	 * 
	 * @param ensemblProteinID
	 * @param proteinSequence
	 * @param proteinExons
	 */
	public void addGeneProduct(String ensemblProteinID, String proteinSequence, Exon[] proteinExons)
	{
		// create new protein
		Protein newProtein = new Protein(ensemblProteinID, null, proteinSequence);
		
		// store references to exons
		newProtein.addProperty(PROPERTY_EXONS, proteinExons);
		
		// add exons to gene
		for(int i=0; i<proteinExons.length; i++)
		{
			addExon(proteinExons[i]);
		}
		
		// add protein to gene products
		geneProducts.put(ensemblProteinID, newProtein);
	}
	
	/**
	 * Method returns an iterator instance for the gene products of the gene
	 * 
	 * @return
	 */
	public Iterator<String> getGeneProductNames()
	{
		return geneProducts.keySet().iterator();
	}
	
	/**
	 * Method returns an array of gene product names
	 * @return
	 */
	public String[] getArrayOfGeneProductNames()
	{
		return (String[])geneProducts.keySet().toArray(new String[0]);
	}
	
	/**
	 * Method returns an array of gene products
	 * 
	 * @return
	 */
	public Protein[] getGeneProducts()
	{
		return geneProducts.values().toArray(new Protein[geneProducts.size()]);
	}
	
	/**
	 * Method returns the gene product for a specific name or null, if the gene
	 * product does not exist.
	 * 
	 * @param name
	 * @return
	 */
	public Protein getGeneProductForName(String name)
	{
		return geneProducts.get(name);
	}
	
	/**
	 * Method returns the number of gene products of the gene
	 * 
	 * @return
	 */
	public int getNumberOfGeneProducts()
	{
		return geneProducts.size();
	}
	
	/**
	 * Method returns the number of exons for the gene
	 * @return
	 */
	public int getNumberOfExons()
	{
		return exons.size();
	}
	
	/**
	 * Method returns true if the gene is a plus strand gene, false otherwise
	 * 
	 * @return
	 */
	public boolean isPlusStrand()
	{
		return plusStrand;
	}

	/**
	 * Method returns the chromosome id of the gene
	 * 
	 * @return
	 */
	public String getChromosome()
	{
		return chromosome;
	}

	/**
	 * Method returns the start position of the gene
	 * 
	 * @return
	 */
	public int getStart()
	{
		return start;
	}

	/**
	 * Method returns the stop position of the gene
	 * 
	 * @return
	 */
	public int getStop()
	{
		return stop;
	}

	/**
	 * Method returns the ID of the gene
	 * 
	 * @return
	 */
	public String getGeneID()
	{
		return geneID;
	}
	
	/**
	 * Method returns gene symbol of the gene
	 * 
	 * @return
	 */
	public String getGeneName()
	{
		return geneName;
	}
	
	public void setGeneName(String strGeneName)
	{
		geneName = strGeneName;
	}
	
	public String getGeneType()
	{
		return geneType;
	}
	
	/**
	 * Method returns the strand identifier of the gene (+ or -)
	 * 
	 * @return
	 */
	public char getStrandStringID()
	{
		if(isPlusStrand())
			return '+';
		else 
			return '-';
	}
	
	/**
	 * Method returns exon with specified id or null if exon does not exist.
	 * 
	 * @param exonID
	 * @return
	 */
	public Exon getExonForExonID(int exonID)
	{
		for(Exon exon : exons)
		{
			if(exon.getExonID() == exonID)
				return exon;
		}
		
		return null;
	}
	
	/**
	 * This method returns an array containing the stop position of each exon in the sequence at the
	 * specifiedposition ofthe exon in the protein sequence. I.e. if exon 4 stops at position 200 in the
	 * protein sequence, position 3 in the array (start index from 0) has value 200. For minus strand
	 * genes the values have to be reversed!
	 * 
	 * @param geneProductID
	 * @return
	 */
	public int[] getExonStopPositionsInSequenceForGeneProduct(String geneProductID)
	{
		Protein protein = geneProducts.get(geneProductID);
	
		int[] exonStartPositionsInSequence = getExonStartPositionsInSequenceForGeneProduct(geneProductID);
		int[] exonStopPositionsInSequence = new int[exonStartPositionsInSequence.length];
		
		if(isPlusStrand())
		{
			for(int i=0; i<exonStartPositionsInSequence.length-1; i++)
			{
				exonStopPositionsInSequence[i] = exonStartPositionsInSequence[i+1] - 1;
			}
			
			exonStopPositionsInSequence[exonStopPositionsInSequence.length-1] = protein.getSequenceLength()-1;
		}
		else
		{
			for(int i=exonStartPositionsInSequence.length-1; i>0; i--)
			{
				exonStopPositionsInSequence[i] = exonStartPositionsInSequence[i-1] - 1;
			}
			
			exonStopPositionsInSequence[0] = protein.getSequenceLength()-1;
		}
		
		return exonStopPositionsInSequence;
	}
	
	/**
	 * Method returns the mRNA sequence of the gene product
	 * 
	 * @param geneProductID
	 * @return
	 */
	public String getmRNASequenceForGeneProduct(String geneProductID)
	{
		Protein protein = geneProducts.get(geneProductID);
		StringBuffer buffer = new StringBuffer();
		
		Exon[] exons = (Exon[])protein.getProperty(PROPERTY_EXONS);
		
		// sort exons to make sure that they are sorted!
		TreeSet<Exon> sortedExons = new TreeSet<Exon>();
		
		for(Exon e : exons)
		{
			sortedExons.add(e);
		}
		
		exons = sortedExons.toArray(new Exon[sortedExons.size()]);
		
		if(isPlusStrand())
		{
			for(int e=0; e<exons.length; e++)
			{
				buffer.append(exons[e].getCodingNucleotideSequence());
			}
		}
		else
		{
			for(int e=exons.length-1; e>=0; e--)
			{
				buffer.append(exons[e].getCodingNucleotideSequence());
			}
		}
		
		// check translation
		/*if(!protein.getSequence().equals(GeneticCodeTranslator.translateSequenceToAminoAcid(buffer.toString()).replaceAll("\\*", "")))
		{
			System.out.println("ERROR: Translation and stored sequence disagree for transcript " + protein.getID() + " and gene: " + getGeneID());
			System.out.println("Protein: " + protein.getSequence());
			System.out.println("Transla: " + GeneticCodeTranslator.translateSequenceToAminoAcid(buffer.toString()));
			
		}*/
		
		return buffer.toString();
			
	}
	
	public Exon[] getSortedExonsForGeneProduct(String productID)
	{
		Protein protein = geneProducts.get(productID);
		
		Exon[] exons = (Exon[])protein.getProperty(PROPERTY_EXONS);
		
		// sort exons to make sure that they are sorted!
		TreeSet<Exon> sortedExons = new TreeSet<Exon>();
		
		for(Exon e : exons)
		{
			sortedExons.add(e);
		}
		
		exons = sortedExons.toArray(new Exon[sortedExons.size()]);
		
		return exons;
	}
	
	public Exon[] getUnsortedExonsForGeneProduct(String productID)
	{
		Protein protein = geneProducts.get(productID);
		Exon[] exons = (Exon[])protein.getProperty(PROPERTY_EXONS);
		
		return exons;
	}
	
	public Exon[] getSortedExonsForGeneProduct(String productID, boolean bSafeMode)
	{
		Protein protein = (Protein)this.geneProducts.get(productID);
		Exon[] exons = (Exon[])protein.getProperty("PROPERTY_EXONS");
    
		TreeSet<Exon> sortedExons = new TreeSet<Exon>();
		Exon[] arrayOfExon1;
		int j = (arrayOfExon1 = exons).length;
		for (int i = 0; i < j; i++)
		{
			Exon e = arrayOfExon1[i];
      
			sortedExons.add(e);
		}
		exons = (Exon[])sortedExons.toArray(new Exon[sortedExons.size()]);
		if (bSafeMode)
		{
			Exon[] pAllExons = getSortedExons();
			Exon[] arrayOfExon2;
			int k = (arrayOfExon2 = exons).length;
			for (j = 0; j < k; j++)
			{
				Exon ex = arrayOfExon2[j];
				Exon[] arrayOfExon3;
				int n = (arrayOfExon3 = pAllExons).length;
				for (int m = 0; m < n; m++)
				{
					Exon ex2 = arrayOfExon3[m];
					if ((ex.getCodingStart() == ex2.getCodingStart()) && (ex.getCodingStop() == ex2.getCodingStop()))
					{
						ex.setID(ex2.getExonID());
						break;
					}
				}
			}
		}
		return exons;
  }
	
	/**
	 * Method returns an array of integer values which correspond to the exon which contributes
	 * the base at the respective position
	 * 
	 * @param geneProductID
	 * @return
	 */
	public int[] getExonPositionsForGeneProduct(String geneProductID)
	{
		Protein protein = geneProducts.get(geneProductID);
		
		Exon[] exons = (Exon[])protein.getProperty(PROPERTY_EXONS);
		
		int cDNALength = 0;
		
		for(int i=0; i<exons.length; i++)
			cDNALength += exons[i].getGenomicLength();
		
		int[] exonPositions = new int[cDNALength];
		
		// sort exons to make sure that they are sorted!
		TreeSet<Exon> sortedExons = new TreeSet<Exon>();
		
		for(Exon e : exons)
		{
			sortedExons.add(e);
		}
		
		exons = sortedExons.toArray(new Exon[sortedExons.size()]);
		
		int positionInArray = 0;
		
		if(isPlusStrand())
		{
			for(int e=0; e<exons.length; e++)
			{
				int exonID = exons[e].getExonID();
				
				for(int i=0; i<exons[e].getGenomicLength(); i++)
				{
					exonPositions[positionInArray] = exonID;
					positionInArray++;
				}
			}
		}
		else
		{
			for(int e=exons.length-1; e>=0; e--)
			{
				int exonID = exons[e].getExonID();
			
				for(int i=0; i<exons[e].getGenomicLength(); i++)
				{
					exonPositions[positionInArray] = exonID;
					positionInArray++;
				}				
			}
		}
		
		return exonPositions;
			
	}
	
	public int getGenomicCodingStartForGeneProduct(String geneProductID)
	{
		Protein protein = geneProducts.get(geneProductID);
		
		Exon[] exons = (Exon[])protein.getProperty(PROPERTY_EXONS);
		
		// sort exons to make sure that they are sorted!
		TreeSet<Exon> sortedExons = new TreeSet<Exon>();
		
		for(Exon e : exons)
		{
			sortedExons.add(e);
		}
		
		exons = sortedExons.toArray(new Exon[sortedExons.size()]);
		
		if(isPlusStrand())
		{
			for(Exon e : exons)
			{
				if(e.getCodingStart() != -1)
					return e.getCodingStart();
			}
		}
		else
		{
			for(int i=exons.length-1;i>=0;i--)
			{
				if(exons[i].getCodingStop() != -1)
					return exons[i].getCodingStop();
			}
		}
		
		return -1;
	}
	
	/**
	 * Method returns an array of integer values which correspond to the genomic position of
	 * the base at the respective position
	 * 
	 * @param geneProductID
	 * @return
	 */
	public int[] getGenomicPositionsForGeneProduct(String geneProductID)
	{
		Protein protein = geneProducts.get(geneProductID);
		
		Exon[] exons = (Exon[])protein.getProperty(PROPERTY_EXONS);
		
		int cDNALength = 0;
		
		for(int i=0; i<exons.length; i++)
			cDNALength += exons[i].getGenomicLength();
		
		int[] genomicPositions = new int[cDNALength];
		
		// sort exons to make sure that they are sorted!
		TreeSet<Exon> sortedExons = new TreeSet<Exon>();
		
		for(Exon e : exons)
		{
			sortedExons.add(e);
		}
		
		exons = sortedExons.toArray(new Exon[sortedExons.size()]);
		
		int positionInArray = 0;
		
		if(isPlusStrand())
		{
			for(int e=0; e<exons.length; e++)
			{
				for(int i=0; i<exons[e].getGenomicLength(); i++)
				{
					genomicPositions[positionInArray] = exons[e].getGenomicStart()+i;
					positionInArray++;
				}
			}
		}
		else
		{
			for(int e=exons.length-1; e>=0; e--)
			{
				for(int i=0; i<exons[e].getGenomicLength(); i++)
				{
					genomicPositions[positionInArray] = exons[e].getGenomicStop()-i;
					positionInArray++;
				}				
			}
		}
		
		return genomicPositions;
			
	}
	
	/**
	 * Method computes the length of the cDNA of the corresponding transcript
	 * @param geneProductID
	 * @return
	 */
	public int getTranscriptCDNALength(String geneProductID)
	{
		Protein protein = geneProducts.get(geneProductID);
		Exon[] exons = (Exon[])protein.getProperty(PROPERTY_EXONS);
		
		int length=0;
		
		for(Exon e : exons)
		{
			length += e.getGenomicLength();
		}
		
		return length;
	}
	
	/**
	 * Method returns the mRNA sequence of the gene product
	 * 
	 * @param geneProductID
	 * @return
	 */
	public String getcDNASequenceForGeneProduct(String geneProductID)
	{
		Protein protein = geneProducts.get(geneProductID);
		StringBuffer buffer = new StringBuffer();
		
		Exon[] exons = (Exon[])protein.getProperty(PROPERTY_EXONS);
		
		// sort exons to make sure that they are sorted!
		TreeSet<Exon> sortedExons = new TreeSet<Exon>();
		
		for(Exon e : exons)
		{
			sortedExons.add(e);
		}
		
		exons = sortedExons.toArray(new Exon[sortedExons.size()]);
		
		if(isPlusStrand())
		{
			for(int e=0; e<exons.length; e++)
			{
				buffer.append(exons[e].getGenomicNucleotideSequence());
			}
		}
		else
		{
			for(int e=exons.length-1; e>=0; e--)
			{
				buffer.append(exons[e].getGenomicNucleotideSequence());
			}
		}
		
		return buffer.toString();
			
	}
	
	/**
	 * This method returns an array containing the start position of each exon in the sequence at the specified
	 * position of the exon in the protein sequence i.e. if exon 4 starts at position 120 in the protein sequence,
	 * position 3 in the array (start index from 0) has value 120. Please notice that for minus strand genes, the
	 * values have to be reversed, starting with 0 at position length-1.
	 * 
	 * @param geneProductID
	 * @return
	 */
	public int[] getExonStartPositionsInSequenceForGeneProduct(String geneProductID)
	{
		Protein protein = geneProducts.get(geneProductID);
		
		Exon[] exons = (Exon[])protein.getProperty(PROPERTY_EXONS);
		
		// sort exons to make sure that they are sorted!
		TreeSet<Exon> sortedExons = new TreeSet<Exon>();
		
		for(Exon e : exons)
		{
			sortedExons.add(e);
		}
		
		exons = sortedExons.toArray(new Exon[sortedExons.size()]);
		
		int[] exonStartPositionsInSequence = new int[exons.length];
		
		int startInProteinSequence = 0;
		int lastExonLeftOver = 0;
		
		if(isPlusStrand())
		{
			for(int e=0; e<exons.length; e++)
			{
				// add last exon left over to it such that no nucleotides are lost
				int exonLength = (exons[e].getCodingStop() - exons[e].getCodingStart())+lastExonLeftOver+1;
				int completeCodons = exonLength / 3;
				lastExonLeftOver = exonLength % 3;

				exonStartPositionsInSequence[e] = startInProteinSequence;
				
				startInProteinSequence = startInProteinSequence + completeCodons;
			}
		}
		else
		{
			for(int e=exons.length-1; e>=0; e--)
			{
				// add last exon left over to it such that no nucleotides are lost
				int exonLength = (exons[e].getCodingStop() - exons[e].getCodingStart())+lastExonLeftOver+1;
				int completeCodons = exonLength / 3;
				lastExonLeftOver = exonLength % 3;
				
				exonStartPositionsInSequence[e] = startInProteinSequence;
			
				startInProteinSequence = startInProteinSequence + completeCodons;
			}
		}
		
		return exonStartPositionsInSequence;
	}
	
	/**
	 * Method returns all sequences of exons where the exons are the keys and their sequences in the
	 * different gene products are the values.
	 * 
	 * @return
	 */
	public TreeMap<Exon, TreeSet<String>> getSequencesOfAllExonsExonToSequencesMap()
	{
		TreeMap<Exon, TreeSet<String>> exonToSequences = new TreeMap<Exon, TreeSet<String>>();
		
		Iterator<String> iter = geneProducts.keySet().iterator();
		
		while(iter.hasNext())
		{
			String geneProductID = iter.next();
		
			TreeMap<Exon, String> exonSequences = getExonSequencesForGeneProduct(geneProductID);
			
			Iterator<Exon> exons = exonSequences.keySet().iterator();
			
			while(exons.hasNext())
			{
				Exon current = exons.next();
				String sequence = exonSequences.get(current);
				
				if(exonToSequences.containsKey(current))
				{
					exonToSequences.get(current).add(sequence);
				}
				else
				{
					TreeSet<String> set = new TreeSet<String>();
					set.add(sequence);
					exonToSequences.put(current, set);
				}
			}
		}
		
		return exonToSequences;
	}
	
	/**
	 * Method returns all sequences of exons as strings (keys) and their corresponding exon as value
	 * 
	 * @return
	 */
	public TreeMap<String, Exon> getSequencesOfAllExonsSequenceToExonMap()
	{
		TreeMap<String, Exon> exonSequenceToExon = new TreeMap<String, Exon>();
		
		Iterator<String> iter = geneProducts.keySet().iterator();
		
		while(iter.hasNext())
		{
			String geneProductID = iter.next();
		
			TreeMap<Exon, String> exonSequences = getExonSequencesForGeneProduct(geneProductID);
			
			Iterator<Exon> exons = exonSequences.keySet().iterator();
			
			while(exons.hasNext())
			{
				Exon current = exons.next();
				exonSequenceToExon.put(exonSequences.get(current), current);
			}
		}
		
		return exonSequenceToExon;
	}
	
	/**
	 * Method takes the ID of a gene product of the gene and slices the protein sequence in pieces such that 
	 * each exon is assigned to its protein sequence. Residues that are not fully created by one exon are 
	 * assigned to the exon that fills them up such that they are created. 
	 * 
	 * @param geneProductID
	 * @return
	 */
	public TreeMap<Exon, String> getExonSequencesForGeneProduct(String geneProductID)
	{
		Protein protein = geneProducts.get(geneProductID);
		
		TreeMap<Exon, String> exonSequences = new TreeMap<Exon, String>();
		
		this.getExonsForGeneProduct(geneProductID);
		
		Exon[] exons = (Exon[])protein.getProperty(PROPERTY_EXONS);
		
		System.out.println("Gene Product: " + geneProductID + " Exons: " + exons.length);
		
		// sort exons to make sure that they are sorted!
		TreeSet<Exon> sortedExons = new TreeSet<Exon>();
		
		for(Exon e : exons)
		{
			sortedExons.add(e);
		}
		
		exons = sortedExons.toArray(new Exon[sortedExons.size()]);
		
		int startInProteinSequence = 0;
		int lastExonLeftOver = 0;
		
		if(isPlusStrand())
		{
			for(int e=0; e<exons.length; e++)
			{
				// add last exon left over to it such that no nucleotides are lost
				int exonLength = (exons[e].getCodingStop() - exons[e].getCodingStart())+lastExonLeftOver+1;
				int completeCodons = exonLength / 3;
				lastExonLeftOver = exonLength % 3;
				
				System.out.println("Start in Protein sequence: " + startInProteinSequence + " exonLength: " + exonLength + " protein length: " + protein.getAminoAcids().length);
				
				String exonAminoAcidSequence = protein.getSequence().substring(startInProteinSequence, Math.min(protein.getSequence().length(), startInProteinSequence + completeCodons));
				
				exonSequences.put(exons[e], exonAminoAcidSequence);
				
				startInProteinSequence = startInProteinSequence + completeCodons;
			}
		}
		else
		{
			for(int e=exons.length-1; e>=0; e--)
			//for(int e=0; e<exons.length; e++)
			{
				// add last exon left over to it such that no nucleotides are lost
				int exonLength = (exons[e].getCodingStop() - exons[e].getCodingStart())+lastExonLeftOver+1;
				int completeCodons = exonLength / 3;
				lastExonLeftOver = exonLength % 3;
				
				String exonAminoAcidSequence = protein.getSequence().substring(startInProteinSequence, 
					Math.min(protein.getSequence().length(), startInProteinSequence + completeCodons));
				
				exonSequences.put(exons[e], exonAminoAcidSequence);
				
				startInProteinSequence = startInProteinSequence + completeCodons;
			}
		}
		
		return exonSequences;
	}
	
	/**
	 * Method takes the ID of a gene product of the gene and slices the protein sequence in pieces such that 
	 * each exon is assigned to its protein sequence. Residues that are not fully created by one exon are 
	 * assigned to the exon that fills them up such that they are created. 
	 * 
	 * Please note that the position the exon sequence is stored to corresponds to the exon position in the gene sequence!!
	 * Therefore, minus strand genes have their exons reverted and the last exon contributes the first protein sequence part!!
	 * 
	 * @param geneProductID
	 * @return
	 */
	public String[] getExonSequencesForGeneProductInArray(String geneProductID)
	{
		Protein protein = geneProducts.get(geneProductID);
		
		String[] exonSequences = new String[this.getExonsForGeneProduct(geneProductID).length];
		
		Exon[] exons = (Exon[])protein.getProperty(PROPERTY_EXONS);
		
		// sort exons to make sure that they are sorted!
		TreeSet<Exon> sortedExons = new TreeSet<Exon>();
		
		for(Exon e : exons)
		{
			sortedExons.add(e);
		}
		
		exons = sortedExons.toArray(new Exon[sortedExons.size()]);
		
		int startInProteinSequence = 0;
		int lastExonLeftOver = 0;
		
		if(isPlusStrand())
		{
			for(int e=0; e<exons.length; e++)
			{
				// add last exon left over to it such that no nucleotides are lost
				int exonLength = (exons[e].getCodingStop() - exons[e].getCodingStart())+lastExonLeftOver+1;
				int completeCodons = exonLength / 3;
				lastExonLeftOver = exonLength % 3;
				
				String exonAminoAcidSequence = protein.getSequence().substring(startInProteinSequence, Math.min(protein.getSequence().length(), startInProteinSequence + completeCodons));
				
				exonSequences[e] = exonAminoAcidSequence;
				
				startInProteinSequence = startInProteinSequence + completeCodons;
			}
		}
		else
		{
			for(int e=exons.length-1; e>=0; e--)
			//for(int e=0; e<exons.length; e++)
			{
				// add last exon left over to it such that no nucleotides are lost
				int exonLength = (exons[e].getCodingStop() - exons[e].getCodingStart())+lastExonLeftOver+1;
				int completeCodons = exonLength / 3;
				lastExonLeftOver = exonLength % 3;
				
				String exonAminoAcidSequence = protein.getSequence().substring(startInProteinSequence, Math.min(protein.getSequence().length(), startInProteinSequence + completeCodons));

				exonSequences[e] = exonAminoAcidSequence;
				
				startInProteinSequence = startInProteinSequence + completeCodons;
			}
		}
		
		return exonSequences;
	}
	
	/**
	 * This method returns an array of integer values which correspond to the positions of the exon
	 * of a specific gene product in the gene product sequence. So if a gene consists of two exons,
	 * both being 4 residues (12 Nucleotides) long, the array will look like this 00001111
	 * 
	 * @param geneProduct
	 * @return
	 */
	public int[] exonNumbersPerSequencePosition(String geneProduct)
	{	
		int[] startPositions = this.getExonStartPositionsInSequenceForGeneProduct(geneProduct);
		int[] stopPositions = this.getExonStopPositionsInSequenceForGeneProduct(geneProduct);
		
		Protein protein = geneProducts.get(geneProduct);
		
		Exon[] exons = (Exon[])protein.getProperty(Gene.PROPERTY_EXONS);
		
		// sort exons to make sure that they are sorted!
		TreeSet<Exon> sortedExons = new TreeSet<Exon>();
		
		for(Exon e : exons)
		{
			sortedExons.add(e);
		}
		
		exons = sortedExons.toArray(new Exon[sortedExons.size()]);
		
		HashMap<Exon, Integer> exonPositions = this.getExonPositionsInSortedGeneExonArray();
		
		int[] exonNumbersPerSequencePosition = new int[protein.getAminoAcids().length];
		
		if(this.isPlusStrand())
		{
			for(int i=0; i<exons.length; i++)
			{
				Exon e = exons[i];
				int exonPositionInSortedArray = exonPositions.get(e);
				int start = startPositions[i];
				int stop = stopPositions[i];
				
				for(int j=start; j<=stop; j++)
				{
					try
					{
						exonNumbersPerSequencePosition[j] = exonPositionInSortedArray;
					}
					catch(Exception ex)
					{
						// do nothing
					}
				}
			}
		}
		else
		{
			for(int i=exons.length-1; i>=0; i--)
			{
				Exon e = exons[i];
				int exonPositionInSortedArray = exonPositions.get(e);
				int start = startPositions[i];
				int stop = stopPositions[i];
				
				for(int j=start; j<=stop; j++)
				{
					try
					{
						exonNumbersPerSequencePosition[j] = exonPositionInSortedArray;
					}
					catch(Exception ex)
					{
						// do nothing
					}
				}
			}
		}
		
		return exonNumbersPerSequencePosition;
	}
	
	/**
	 * Method returns true if gene product contains the specified exon
	 * @param geneProductName
	 * @param exonID
	 * @return
	 */
	public boolean geneProductContainsExon(String geneProductName, int exonID)
	{
		Exon[] proteinExons = (Exon[])(geneProducts.get(geneProductName).getProperty(PROPERTY_EXONS));
		
		for(Exon exon : proteinExons)
		{
			if(exon.getExonID() == exonID)
				return true;
		}
		
		return false;
	}
	
	/**
	 * Method returns an array containing the positions of the exons that build the gene product.
	 * The positions are relative to the total number of exons of the gene at the time of calling
	 * this method.
	 * 
	 * @param geneProductName
	 * @return
	 */
	public int[] getExonsForGeneProduct(String geneProductName)
	{
		Protein geneProduct = geneProducts.get(geneProductName);
		
		Exon[] proteinExons = (Exon[])geneProduct.getProperty(PROPERTY_EXONS);
		
		// sort exons to make sure that they are sorted!
		TreeSet<Exon> sortedProteinExons = new TreeSet<Exon>();
		
		for(Exon e : proteinExons)
		{
			sortedProteinExons.add(e);
		}
		
		proteinExons = sortedProteinExons.toArray(new Exon[sortedProteinExons.size()]);
		
		Exon[] geneExons = getSortedExons();
		
		int[] exonPositionsForGeneProduct = new int[proteinExons.length];
		
		NEXTEXON: for(int j=0; j<proteinExons.length; j++)
		{
			Exon currentExon = proteinExons[j];
			
			for(int i=0; i<geneExons.length; i++)
			{
				if(currentExon.equals(geneExons[i]))
				{
					exonPositionsForGeneProduct[j] = i;
					continue NEXTEXON;
				}
			}
		}
		
		return exonPositionsForGeneProduct;
	}
	
	/**
	 * Method returns the distance of the specified location from the 3' end of the gene
	 * @param location
	 * @return
	 */
	public int computeDistanceToThreePrimeEndOfGene(int location)
	{
		return Math.abs(stop-location);
	}
	
	/**
	 * Method returns the genomic loacion of the three prime end of a transcript
	 * @param productName
	 * @return
	 */
	public int getThreePrimeEndOfGeneProduct(String productName)
	{
		int[] exonsOfProduct = getExonsForGeneProduct(productName);
		Exon[] sortedExons = getSortedExons();
		
		if(isPlusStrand())
		{
			int threePrimeStop = Integer.MIN_VALUE;
			
			for(int i : exonsOfProduct)
			{
				if(sortedExons[i].getGenomicStop() > threePrimeStop)
				{
					threePrimeStop = sortedExons[i].getGenomicStop();
				}
			}
			
			return threePrimeStop;
		}
		else
		{
			int threePrimeStop = Integer.MAX_VALUE;
			
			for(int i : exonsOfProduct)
			{
				if(sortedExons[i].getGenomicStart() < threePrimeStop)
				{
					threePrimeStop = sortedExons[i].getGenomicStart();
				}
			}
			
			return threePrimeStop;
		}
	}
	
	/**
	 * Method checks if the specified exon is a known three prime exon in some transcript
	 * and returns the transcript name. If multiple transcripts share the same 3' end, 
	 * the FIRST HIT is returned. If the exon is not a known 3' end, null is returned. 
	 *
	 * @param exonID
	 * @return
	 */
	public String isExonKnownThreePrimeEndOfTranscript(int exonID)
	{
		String[] names = getArrayOfGeneProductNames();
		
		for(String s : names)
		{
			if(exonID == getThreePrimeExonOfGeneProduct(s))
				return s;
		}
		
		return null;
	}
	
	
	/**
	 * Method returns the id of the exon which makes up the 3' end of the specified transcript
	 * 
	 * @param geneProduct
	 * @return
	 */
	public int getThreePrimeExonOfGeneProduct(String geneProduct)
	{
		int[] exonsOfProduct = getExonsForGeneProduct(geneProduct);
		Exon[] sortedExons = getSortedExons();
		
		int exonID = -1;
		
		if(isPlusStrand())
		{
			int threePrimeStop = Integer.MIN_VALUE;
			
			for(int i : exonsOfProduct)
			{
				if(sortedExons[i].getGenomicStop() > threePrimeStop)
				{
					threePrimeStop = sortedExons[i].getGenomicStop();
					exonID = sortedExons[i].getExonID();
				}
			}
			
			return exonID;
		}
		else
		{
			int threePrimeStop = Integer.MAX_VALUE;
			
			for(int i : exonsOfProduct)
			{
				if(sortedExons[i].getGenomicStart() < threePrimeStop)
				{
					threePrimeStop = sortedExons[i].getGenomicStart();
					exonID = sortedExons[i].getExonID();
				}
			}
			
			return exonID;
		}
	}
	
	/**
	 * This method computes the minimal distance to a known three prime end of a known transcript and
	 * returns the transcript with the minimal distance to the calling method
	 * 
	 * @param location
	 * @return
	 */
	public String getTranscriptWithMinimalThreePrimeDistanceToLocation(int location)
	{
		String[] names = this.getArrayOfGeneProductNames();
		
		String transcript = null;
		int minDistance = Integer.MAX_VALUE;
		
		for(String product : names)
		{
			int threePrimeStop = getThreePrimeEndOfGeneProduct(product);
				
			int distance = Math.abs(threePrimeStop-location);
				
			if(distance < minDistance)
			{
				minDistance = distance;
				transcript = product;
			}
		}
		
		return transcript;
	}
	
	/**
	 * Method retuns a short String representation of the gene
	 * @return
	 */
	public String toShortString()
	{
		StringBuffer buffer = new StringBuffer();
		
		buffer.append("> ");
		buffer.append(geneID);
		buffer.append("\t");
		buffer.append(start);
		buffer.append("-");
		buffer.append(stop);
		buffer.append("\t");
		buffer.append(chromosome);
		buffer.append("(");
		
		if(plusStrand)
			buffer.append("+");
		else
			buffer.append("-");
		
		buffer.append(")");
		
		buffer.append(" Exons: ");
		buffer.append(exons.size());
		buffer.append(" Gene Products: ");
		buffer.append(geneProducts.size());
		
		return buffer.toString();
	}
	
	/**
	 * Method returns a string representation of the object
	 * 
	 * @return String
	 */
	public String toString()
	{
		StringBuffer buffer = new StringBuffer();
		
		buffer.append("> ");
		buffer.append(geneID);
		buffer.append("\t");
		buffer.append(start);
		buffer.append("-");
		buffer.append(stop);
		buffer.append("\t");
		buffer.append(chromosome);
		buffer.append("(");
		
		if(plusStrand)
			buffer.append("+");
		else
			buffer.append("-");
		
		buffer.append(")");
		
		Iterator<Exon> iter = exons.iterator();
		
		int counter = 0;
		
		while(iter.hasNext())
		{
			Exon exon = iter.next();
			buffer.append("\n\t");
			buffer.append(counter);
			buffer.append(": ");
			buffer.append(exon.toString());
			counter++;
		}
		
		Iterator<String> proteinIter = geneProducts.keySet().iterator();
		
		while(proteinIter.hasNext())
		{
			buffer.append("\n#");
			
			String proteinID = proteinIter.next();
			Protein geneProduct = geneProducts.get(proteinID);
			
			buffer.append(proteinID);
			buffer.append("\t"); 
			buffer.append(geneProduct.getSequence());
			buffer.append("\t"); 
			
			Exon[] proteinExons = (Exon[])geneProduct.getProperty(PROPERTY_EXONS);
			
			Exon[] geneExons = getSortedExons();
			
			NEXTEXON: for(int j=0; j<proteinExons.length; j++)
			{
				Exon currentExon = proteinExons[j];
				
				for(int i=0; i<geneExons.length; i++)
				{
					if(currentExon.equals(geneExons[i]))
					{
						buffer.append(i +" ");
						continue NEXTEXON;
					}
				}
			}
		}
		
		return buffer.toString();
	}

	/**
	 * Method compares two gene instances based on their chromosomal positions
	 * 
	 * @param o
	 * @return
	 */
	public int compareTo(Gene o)
	{
		return geneID.compareTo(o.getGeneID());
	}
	
	
}
