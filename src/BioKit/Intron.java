package BioKit;

import java.io.IOException;
import java.io.Serializable;

/**
 * This class represents an intron wrapping start, stop and sequence information.
 * 
 * @author Fabian Birzele
 *
 */
public class Intron extends GeneticElement implements Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private int id;
	private CompressedNucleotideSequence nucleotideSequence;
	
	/**
	 * Standard constructor, initializes the exon
	 * 
	 * @param start
	 * @param stop
	 * @param nucleotideSequence
	 * @throws Exception
	 */
	public Intron(int start, int stop, String nucleotideSequence) throws IOException
	{
		super(start, stop, start, stop, false);
		this.nucleotideSequence = new CompressedNucleotideSequence(nucleotideSequence);
	}
	
	public Intron(int start, int stop)
	{
		super(start, stop, start, stop, false);
	}
	
	/**
	 * Method sets id of exon;
	 * 
	 * @param id
	 */
	public void setID(int id)
	{
		this.id = id;
	}
	
	/**
	 * Method returns exon id
	 * @return
	 */
	public int getID()
	{
		return id;
	}
	
	/**
	 * Method returns the nucleotide sequence of the exon
	 * 
	 * @return
	 */
	public String getNucleotideSequence()
	{
		if(nucleotideSequence != null)
			return nucleotideSequence.toString();
		else
			return "";
	}
	
	/**
	 * set sequence
	 * 
	 * @param nucleotideSequence
	 */
	public void setNucleotideSequence(String nucleotideSequence) throws Exception
	{
		this.nucleotideSequence = new CompressedNucleotideSequence(nucleotideSequence);
	}
	
	/**
	 * Method returns string representation of exon
	 */
	public String toString()
	{
		StringBuffer buffer = new StringBuffer();
		
		buffer.append(id);
		buffer.append("\t");
		buffer.append(getGenomicStart());
		buffer.append("\t");
		buffer.append(getGenomicStop());
		buffer.append("\t");
		buffer.append(getNucleotideSequence());
		
		return buffer.toString();
	}
	
	public int getGenomicLength()
	{
		return Math.abs(getGenomicStart() - getGenomicStop())+1;
	}
	
	/**
	 * Method returns true if one of the intron is fully contained in the other intron
	 * @param other
	 * @return
	 */
	public boolean intersects(Intron other)
	{
		if((getGenomicStart() >= other.getGenomicStart() && getGenomicStop() <= other.getGenomicStop()) || (getGenomicStart() <= other.getGenomicStart() && getGenomicStop() >= other.getGenomicStop()))
			return true;
		else 
			return false;
	}
	
	public boolean overlaps(Intron other)
	{
		// test for fully contained
		if(intersects(other))
			return true;
		// test for overlap, start of one exon lies within borders of the other
		if((getGenomicStart() <= other.getGenomicStart() && getGenomicStop() >= other.getGenomicStart()) || (getGenomicStart() >= other.getGenomicStart() && getGenomicStart() <= other.getGenomicStop()))
			return true;
		else 
			return false;
	}
}
