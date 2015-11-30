package BioKit;

import java.util.Vector;

/**
 * 
 * @author Fabian Birzele
 */
public class IntronGroup implements Comparable<IntronGroup>
{
	private int groupID;
	private Vector<Integer> intronIDs;
	private Vector<Integer> intronDBIDs;
	private Vector<Intron> introns;
	private int minStart;
	private int maxStop;
	private boolean isPlusStrand;
	
	private String genomicSequenceOfGroup;
	
	public IntronGroup(int groupID, boolean isPlusStrand)
	{
		this.groupID = groupID;
		intronIDs = new Vector<Integer>();
		introns = new Vector<Intron>();
		intronDBIDs = new Vector<Integer>();
		
		this.isPlusStrand = isPlusStrand;
		
		minStart = Integer.MAX_VALUE;
		maxStop = Integer.MIN_VALUE;
	}
	
	public void addIntron(int intronID, int intronDBID, Intron intron)
	{
		intronIDs.add(intronID);
		introns.add(intron);
		intronDBIDs.add(intronDBID);
		
		if(minStart > intron.getGenomicStart())
			minStart = intron.getGenomicStart();
		
		if(maxStop < intron.getGenomicStop())
			maxStop = intron.getGenomicStop();
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
		
		if(introns.size() == 1)
			return introns.firstElement().getNucleotideSequence();
		
		char[] sequence = new char[getIntronGroupLengthInBp()];
		
		for(int i=0;i<sequence.length; i++)
			sequence[i] = 'N';
		
		for(Intron e : introns)
		{
			int offset = -1;
			
			if(isPlusStrand)
				offset = e.getGenomicStart() - minStart;
			else
				offset = maxStop - e.getGenomicStop();
			
			char[] es = e.getNucleotideSequence().toCharArray();
			
			for(int i=0; i<es.length; i++)
			{
				if(sequence[offset+i] == 'N')
				{
					sequence[offset+i] = es[i];
				}
				// if there is a discrepancy between two introns covering the same position, return empty sequence and print warning
				else if(sequence[offset+i] != es[i])
				{
					System.err.println("Error when computing genomic sequence of group " + groupID + "\t" + isPlusStrand +"\t" + minStart + "\t" + maxStop + "\t" + getIntronGroupLengthInBp());
					
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
	

	/**
	 * Method returns the genomic sequence of the group
	 * 
	 * @return
	 */
	public String getShortestIntronSequence()
	{
		int shortestLength = Integer.MAX_VALUE;
		Intron shortestLengthIntron = null;
		
		for(Intron i : introns)
		{
			if(shortestLength > i.getGenomicLength())
			{
				shortestLength = i.getGenomicLength();
				shortestLengthIntron = i;
			}
		}
		
		return shortestLengthIntron.getNucleotideSequence();
	}
	
	/**
	 * Method returns the minimal genomic start of any intron in this group
	 * 
	 * @return
	 */
	public int getGenomicStartOfGroup()
	{
		return minStart;
	}
	
	/**
	 * Method returns the maximal genomic stop of any intron in this group
	 * 
	 * @return
	 */
	public int getGenomicStopOfGroup()
	{
		return maxStop;
	}
	
	/**
	 * Method returns the maximal length of the introns in this group, i.e. the genomic
	 * location covered by this group
	 * 
	 * @return
	 */
	public int getIntronGroupLengthInBp()
	{
		return Math.abs(maxStop - minStart)+1;
	}
	
	/**
	 * Method checks if the groups overlaps with the specified intron
	 * 
	 * @param exon
	 * @return
	 */
	public boolean intronIntersectsWithIntronsInGroup(Intron intron)
	{
		for(Intron groupIntron : introns)
		{
			if(groupIntron.overlaps(intron))
				return true;
		}
		
		return false;
	}
	
	/**
	 * Method returns the ids of the exons in this group
	 * 
	 * @return
	 */
	public Vector<Integer> getIntronIDs()
	{
		return intronIDs;
	}
	
	/**
	 * Method returns the db ids of the exons contained in this group
	 * @return
	 */
	public Vector<Integer> getIntronDBIDs() 
	{
		return intronDBIDs;
	}

	/**
	 * Method returns the exons of this group
	 * 
	 * @return
	 */
	public Vector<Intron> getIntrons()
	{
		return introns;
	}
	
	/**
	 * Method returns true if the specified exon is contained in this group, false otherwise.
	 * 
	 * @param exon
	 * @return
	 */
	public boolean groupContainsIntron(int intron)
	{
		return intronIDs.contains(intron);
	}
	
	/**
	 * Method returns true if the exon group contains the specified exon (identified by its DB id)
	 * 
	 * @param exonDBID
	 * @return
	 */
	public boolean groupContainsIntronWithDBID(int intronDBID)
	{
		return intronDBIDs.contains(intronDBID);
	}
	
	public int getGroupID()
	{
		return groupID;
	}
	
	public int size()
	{
		return intronIDs.size();
	}
	
	public String toShortString()
	{

		StringBuffer buffer = new StringBuffer();
		
		buffer.append(groupID + " (");
		
		for(int i=0; i<intronIDs.size(); i++)
		{
			buffer.append(intronIDs.get(i));
			
			if(i < intronIDs.size()-1)
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
		
		for(int i=0; i<introns.size(); i++)
		{
			buffer.append(intronIDs.get(i) + " | " + introns.get(i) + "\n");
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
	public int compareTo(IntronGroup other)
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