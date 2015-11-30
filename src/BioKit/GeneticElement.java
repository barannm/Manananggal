package BioKit;

import java.io.Serializable;

/**
 * Class acts as superclass of all genetic elements like introns and exons. To sort elements of a gene a comparable method is provided
 * 
 * @author Fabian Birzele
 *
 */
public abstract class GeneticElement implements Comparable<GeneticElement>, Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private int genomicStart;
	private int genomicStop;
	private int codingStart;
	private int codingStop;
	private boolean coding;
	private int id;
	
	/**
	 * Constructor
	 * 
	 * @param genomicStart
	 * @param genomicStop
	 * @param codingStart
	 * @param codingStop
	 * @param coding
	 */
	public GeneticElement(int genomicStart, int genomicStop, int codingStart, int codingStop, boolean coding)
	{
		this.genomicStart = genomicStart;
		this.genomicStop = genomicStop;
		this.codingStart = codingStart;
		this.codingStop = codingStop;
		this.coding = coding;
	}
	
	public GeneticElement()
	{
	}
	
	/**
	 * Method sets the id of this genetic element
	 * @param id
	 */
	public void setElementID(int id)
	{
		this.id = id;
	}
	
	/**
	 * Method returns the id of this genetic element
	 * 
	 * @return
	 */
	public int getElementID()
	{
		return id;
	}

	/**
	 * Method returns the genomic start position of the element
	 * 
	 * @return
	 */
	public int getGenomicStart()
	{
		return genomicStart;
	}

	/**
	 * Method returns the genomic stop position of the element
	 * 
	 * @return
	 */
	public int getGenomicStop()
	{
		return genomicStop;
	}

	/**
	 * Method returns the coding start position of the element
	 * 
	 * @return
	 */
	public int getCodingStart()
	{
		return codingStart;
	}

	/**
	 * Method returns the coding stop position of the element
	 * 
	 * @return
	 */
	public int getCodingStop()
	{
		return codingStop;
	}

	/**
	 * Method returns true if the element is coding, false otherwise
	 * 
	 * @return
	 */
	public boolean isCoding()
	{
		return coding;
	}
	
	/**
	 * Standard compare to, compares genetic elements according to their start and end positions. Elements are sorted
	 * in ascending order with respect to their start positions. Additionally, if ids are set for the corresponding 
	 * elements, the id is checked additionally.
	 * 
	 * @param o
	 * @return
	 */
	public int compareTo(GeneticElement other)
	{
		if(genomicStart < other.genomicStart)
	           return -1;
	    else if(genomicStart == other.genomicStart)
	    {
	        // genomic regions are identical
	    	if(genomicStop == other.genomicStop)
	        {
	    		if(codingStart == other.codingStart && codingStop == other.codingStop)
	            {
	    			if(id != -1 && id == other.id)
	    				return 0;
	    			else if(id != -1 && id < other.id)
	    				return -1;
	    			else if(id != -1 && id > other.id)
	    				return 1;
	    			
	    			return 0;
	            }
	            else if(codingStart < other.codingStart)
	            	return -1;
	            else if(codingStop < other.codingStop)
	            	return -1;
	            else
	            	return 1;
	        }
	    	else if(genomicStop < other.genomicStop)
	    		return -1;
	    	else
	    		return 1;
	    }
	    else
	    	return 1;
	}
	
	/**
	 * Standard equals, compares genetic elements according to their start and end positions. 
	 * Additionally, if ids are set for the corresponding elements, the id is checked additionally.
	 * 
	 * @param o
	 * @return
	 */
	public boolean equals(GeneticElement other)
	{
		if(genomicStart != other.genomicStart) 
			return false;
	    else if(genomicStop != other.genomicStop)
	    	return false;
	    else if (codingStart != other.codingStart)
	    	return false;
	    else if (codingStop != other.codingStop)
	    	return false;
	    else if (id != -1 && id != other.id)
	    	return false;
	    else return true;
	}
}
