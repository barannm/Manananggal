package BioKit;

/**
 * Interface representing a sequence. This interface is required to represent any sequence to some algorithms
 * of Biokit. 
 * 
 * @author Fabian Birzele
 */
public interface Sequence
{
	/**
	 * Method returns the element of the alphabet stored at position x in the sequence.
	 * Sequences start with index 0.
	 * 
	 * @param position
	 * @return
	 */
	public String getSequenceElement(int position);
	
	/**
	 * Method returns the object of the sequence stored at position x in the sequence.
	 * Sequences start with index 0;
	 * 
	 * @param position
	 * @return
	 */
	public Object getSequenceObjectAt(int position);
	
	
	/**
	 * Method returns the length of the sequence. Since the sequence starts at position
	 * 0 the last element can be retrieved at position sequenceLength-1.
	 * 
	 * @return
	 */
	public int getSequenceLength();
	
	/**
	 * Method returns the ID of the sequence
	 * 
	 * @return
	 */
	public String getID();
}
