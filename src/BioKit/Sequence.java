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
