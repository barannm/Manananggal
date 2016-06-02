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

public class DNA implements Sequence
{
	private String id;
	private String sequence;
	
	public DNA(String id, String sequence)
	{
		this.id = id;
		this.sequence = sequence;
	}
	
	public String getSequence()
	{
		return sequence;
	}
	
	public String getID() 
	{
		return id;
	}

	public String getSequenceElement(int position) 
	{
		return sequence.charAt(position)+"";
	}

	public int getSequenceLength() 
	{
		return sequence.length();
	}

	public Object getSequenceObjectAt(int position) 
	{
		return new Nucleotide(sequence.charAt(position));
	}
	
	public Nucleotide[] getNucleotideSequence()
	{
		Nucleotide[] nucleotides = new Nucleotide[sequence.length()];
		
		for(int i=0; i<nucleotides.length; i++)
		{
			nucleotides[i] = new Nucleotide(sequence.charAt(i));
		}
		
		return nucleotides;
	}
	
	public DNA getReverseCompDNASequence()
	{
		StringBuffer rev = new StringBuffer();
	
		char[] seqArr = sequence.toCharArray();
		
		for(int i=seqArr.length-1; i>=0; i--)
		{
			if(seqArr[i]=='A')
				rev.append('T'); 
			else if(seqArr[i]=='T')
				rev.append('A'); 
			else if(seqArr[i]=='G')
				rev.append('C'); 
			else if(seqArr[i]=='C')
				rev.append('G'); 
			else if(seqArr[i]=='N')
				rev.append('N');
		}
		
		return new DNA(id + "_revComp", rev.toString());
	}
	

}
