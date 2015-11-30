package BioKit;

import java.util.Vector;

/**
 * Class representing a group of exons, for example a group of overlapping exons
 * or a group of exons showing identical behaviour in their usage in different 
 * transcripts.
 * 
 * @author Fabian Birzele
 */
public class ExonGroup implements Comparable<ExonGroup>
{
	private int groupID;
	private Vector<Integer> exonsIDs;
	private Vector<Integer> exonsDBIDs;
	private Vector<Exon> exons;
	private int minStart;
	private int maxStop;
	private boolean isPlusStrand;
	
	private String genomicSequenceOfGroup;
	
	public ExonGroup(int groupID, boolean isPlusStrand)
	{
		this.groupID = groupID;
		exonsIDs = new Vector<Integer>();
		exons = new Vector<Exon>();
		exonsDBIDs = new Vector<Integer>();
		
		this.isPlusStrand = isPlusStrand;
		
		minStart = Integer.MAX_VALUE;
		maxStop = Integer.MIN_VALUE;
	}
	
	public void addExon(int exonID, int exonDBID, Exon exon)
	{
		exonsIDs.add(exonID);
		exons.add(exon);
		exonsDBIDs.add(exonDBID);
		
		if(minStart > exon.getGenomicStart())
			minStart = exon.getGenomicStart();
		
		if(maxStop < exon.getGenomicStop())
			maxStop = exon.getGenomicStop();
	}
	
	/**
	 * Method returns true if at least one exon in the group is coding. 
	 * Otherwise it returns false.
	 * 
	 * @return
	 */
	public boolean isCodingGroup()
	{
		for(Exon e : exons)
		{
			if(e.isCoding())
				return true;
		}
		
		return false;
	}
	
	/**
	 * Method returns the genomic sequence of the group
	 * 
	 * @return
	 */
	public String getGenomicSequenceOfGroup()
	{
		if(genomicSequenceOfGroup != null)
			return genomicSequenceOfGroup;
		
		if(exons.size() == 1)
			return exons.firstElement().getGenomicNucleotideSequence();
		
		char[] sequence = new char[getExonGroupLengthInBp()];
		
		for(int i=0;i<sequence.length; i++)
			sequence[i] = 'N';
		
		for(Exon e : exons)
		{
			int offset = -1;
			
			if(isPlusStrand)
				offset = e.getGenomicStart() - minStart;
			else
				offset = maxStop - e.getGenomicStop();
			
			char[] es = e.getGenomicNucleotideSequence().toCharArray();
			
			for(int i=0; i<es.length; i++)
			{
				if(sequence[offset+i] == 'N')
				{
					sequence[offset+i] = es[i];
				}
				// if there is a discrepancy between two exons covering the same position, return empty sequence and print warning
				else if(sequence[offset+i] != es[i])
				{
					System.err.println("Error when computing genomic sequence of group " + getExonDBIDString() + "\t" + isPlusStrand +"\t" + minStart + "\t" + maxStop + "\t" + getExonGroupLengthInBp());
					
					return "";
				}
			}
		}
		
		StringBuffer b = new StringBuffer();
		
		for(char c : sequence)
			b.append(c);
		
		genomicSequenceOfGroup = b.toString();
		
		return b.toString();
	}
	
	public int[] getGenomicPositionsOfGroup()
	{
		int[] positions = new int[getGenomicSequenceOfGroup().length()];
		
		if(isPlusStrand)
		{
			int counter = 0;
			
			for(int i=minStart; i<=maxStop; i++)
			{
				positions[counter] = i;
				counter++;
			}
		}
		else
		{
			int counter = 0;
			
			for(int i=maxStop; i>=minStart; i--)
			{
				positions[counter] = i;
				counter++;
			}
		}
		
		return positions;
	}
	
	/**
	 * Method returns the coding sequence of the group. This means the method returns an String of the length of the exon group.
	 * All non-coding parts of the sequence are marked by N.
	 * 
	 * @return
	 */
	public String getCodingSequenceOfGroup()
	{
		char[] sequence = new char[getExonGroupLengthInBp()];
		
		for(int i=0;i<sequence.length; i++)
			sequence[i] = 'N';
		
		for(Exon e : exons)
		{
			int offset = -1;
			
			if(isPlusStrand)
				offset = e.getCodingStart() - minStart;
			else
				offset = maxStop - e.getCodingStop();
			
			char[] es = e.getCodingNucleotideSequence().toCharArray();
			
			for(int i=0; i<es.length; i++)
			{
				if(sequence[offset+i] == 'N')
				{
					sequence[offset+i] = es[i];
				}
				// if there is a discrepancy between two exons covering the same position, return empty sequence and print warning
				else if(sequence[offset+i] != es[i])
				{
					System.err.println("Error when computing genomic sequence of group " + getExonDBIDString() + "\t" + isPlusStrand +"\t" + minStart + "\t" + maxStop + "\t" + getExonGroupLengthInBp());
					
					return "";
				}
			}
		}
		
		StringBuffer b = new StringBuffer();
		
		for(char c : sequence)
			b.append(c);
		
		return b.toString();
	}
	
	/**
	 * Method returns the number of coding positions in the group.
	 * 
	 * @return
	 */
	public int getNumberOfCodingPositionsInGroup()
	{
		boolean[] sequence = new boolean[getExonGroupLengthInBp()];
		
		for(Exon e : exons)
		{
			int offset = -1;
			
			if(isPlusStrand)
				offset = e.getCodingStart() - minStart;
			else
				offset = maxStop - e.getCodingStop();
			
			int codingLength = e.getCodingStop()-e.getCodingStart();
			
			for(int i=0; i<codingLength; i++)
			{
				sequence[offset+i] = true;
			}
		}
		
		int counter = 0;
		
		for(boolean f : sequence)
		{
			if(f)
				counter++;
		}
		
		return counter;
	}
	
	/**
	 * Method returns the minimal genomic start of any exon in this group
	 * 
	 * @return
	 */
	public int getGenomicStartOfGroup()
	{
		return minStart;
	}
	
	/**
	 * Method returns the maximal genomic stop of any exon in this group
	 * 
	 * @return
	 */
	public int getGenomicStopOfGroup()
	{
		return maxStop;
	}
	
	/**
	 * Method returns the maximal length of the exons in this group, i.e. the genomic
	 * location covered by this group
	 * 
	 * @return
	 */
	public int getExonGroupLengthInBp()
	{
		return Math.abs(maxStop - minStart)+1;
	}
	
	/**
	 * Method checks if the groups overlaps with the specified exon
	 * 
	 * @param exon
	 * @return
	 */
	public boolean exonIntersectsWithExonsInGroup(Exon exon)
	{
		for(Exon groupExon : exons)
		{
			if(groupExon.overlaps(exon))
				return true;
		}
		
		return false;
	}
	
	/**
	 * Method returns true if the group contains only small exons below a specified size measured in
	 * amino acids
	 * 
	 * @param size
	 * @return
	 */
	public boolean containsOnlySmallExons(int size)
	{
		for(Exon e : exons)
		{
			if(e.getLength()/3 > size)
				return false;
		}
		
		return true;
	}
	
	/**
	 * Method returns the ids of the exons in this group
	 * 
	 * @return
	 */
	public Vector<Integer> getExonIDs()
	{
		return exonsIDs;
	}
	
	/**
	 * Method returns the db ids of the exons contained in this group
	 * @return
	 */
	public Vector<Integer> getExonsDBIDs() 
	{
		return exonsDBIDs;
	}

	/**
	 * Method returns the exons of this group
	 * 
	 * @return
	 */
	public Vector<Exon> getExons()
	{
		return exons;
	}
	
	/**
	 * Method returns a string containing all exon ids separated by ,
	 * 
	 * @return
	 */
	public String getExonIDString()
	{
		StringBuffer b = new StringBuffer();
		
		for(int i : exonsIDs)
			b.append(i + ",");
		
		return b.toString();
	}
	
	/**
	 * Method returns a string containing all exon database ids separated by ,
	 * @return
	 */
	public String getExonDBIDString()
	{
		StringBuffer b = new StringBuffer();
		
		for(int i : exonsDBIDs)
			b.append(i + ",");
		
		return b.toString();
	}
	
	/**
	 * Method returns true if the specified exon is contained in this group, false otherwise.
	 * 
	 * @param exon
	 * @return
	 */
	public boolean groupContainsExon(int exon)
	{
		return exonsIDs.contains(exon);
	}
	
	/**
	 * Method returns true if the exon group contains the specified exon (identified by its DB id)
	 * 
	 * @param exonDBID
	 * @return
	 */
	public boolean groupContainsExonWithDBID(int exonDBID)
	{
		return exonsDBIDs.contains(exonDBID);
	}
	
	public int getGroupID()
	{
		return groupID;
	}
	
	public void setGroupID(int nGroupID)
	{
		groupID = nGroupID;
	}
	
	public int size()
	{
		return exonsIDs.size();
	}
	
	public String toShortString()
	{

		StringBuffer buffer = new StringBuffer();
		
		buffer.append(groupID + " (");
		
		for(int i=0; i<exonsIDs.size(); i++)
		{
			buffer.append(exonsIDs.get(i));
			
			if(i < exonsIDs.size()-1)
				buffer.append(" ");
		}
		
		buffer.append(") ");
		
		buffer.append(minStart + " - " + maxStop);
		
		return buffer.toString();
	}
	
	public String toString()
	{
		StringBuffer buffer = new StringBuffer();
		
		buffer.append("Group " + groupID + "\n");
		
		for(int i=0; i<exons.size(); i++)
		{
			buffer.append(exonsIDs.get(i) + " | " + exons.get(i) + "\n");
		}
		
		buffer.append("\n");
		
		return buffer.toString();
	}
	
	/**
	 * Compare exon groups based on their position in the genome.
	 * 
	 * @param o
	 * @return
	 */
	public int compareTo(ExonGroup other)
	{
		if(minStart < other.minStart)
	           return -1;
	    else if(minStart == other.minStart)
	    {
	        // genomic regions are identical
	    	if(maxStop == other.maxStop)
	    			return 0;
	    	else if(maxStop < other.maxStop)
	    		return -1;
	    	else
	    		return 1;
	    }
	    else
	    	return 1;
	}
}