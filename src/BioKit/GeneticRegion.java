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

import java.io.IOException;
import java.io.Serializable;
import java.util.HashSet;

/**
 * This class represents a exonic region wrapping start, stop and sequence information.
 * 
 * @author Lucia Puchbauer
 *
 */
public class GeneticRegion extends GeneticElement implements MappableElement, Sequence, Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private int id;
	private String reference;
	private CompressedNucleotideSequence codingNucleotideSequence;
	private CompressedNucleotideSequence completeNucleotideSequence;
	private Nucleotide[] nucleotides;
	private HashSet<Exon> coveredExons;
	private HashSet<String> coveredTranscripts;
	
	/**
	 * Standard constructor, initializes the exon
	 * 
	 * @param codingStart
	 * @param codingStop
	 * @param nucleotideSequence
	 * @throws Exception
	 */
	public GeneticRegion(int codingStart, int codingStop, String codingNucleotideSequence) throws IOException
	{
		super(codingStart, codingStop, codingStart, codingStop, true);
		this.codingNucleotideSequence = new CompressedNucleotideSequence(codingNucleotideSequence);
		this.coveredExons = new HashSet<Exon>();
		this.coveredTranscripts = new HashSet<String>();
	}
	
	public GeneticRegion(int codingStart, int codingStop)
	{
		super(codingStart, codingStop, codingStart, codingStop, true);
		this.coveredExons = new HashSet<Exon>();
		this.coveredTranscripts = new HashSet<String>();
	}
	
	public GeneticRegion(int codingStart, int codingStop, int genomicStart, int genomicStop)
	{
		super(genomicStart, genomicStop, codingStart, codingStop, true);
		this.coveredExons = new HashSet<Exon>();
		this.coveredTranscripts = new HashSet<String>();
	}
	
	public GeneticRegion(int codingStart, int codingStop, int genomicStart, int genomicStop, String codingNucleotideSequence) throws IOException, Exception
	{
		super(genomicStart, genomicStop, codingStart, codingStop, true);
		
		if(codingStart == -1 && codingStop == -1)
		{
			// do nothing
		}
		else if(codingStop-codingStart+1 != codingNucleotideSequence.length())
		{
			throw new Exception("Coding sequence positions and sequence length do not agree... position length: " + (codingStop-codingStart+1) 
				+ ", sequence length: " + codingNucleotideSequence.length());
		}
			
		this.codingNucleotideSequence = new CompressedNucleotideSequence(codingNucleotideSequence);
		this.coveredExons = new HashSet<Exon>();
		this.coveredTranscripts = new HashSet<String>();
	}
	
	/**
	 * Method sets id of exon;
	 * 
	 * @param id
	 */
	public void setID(int id)
	{
		this.id = id;
		super.setElementID(id);
	}
	
	/**
	 * Method returns length of the exons in nucleotides
	 * 
	 * @return
	 */
	public int getLength()
	{
		return Math.abs(getCodingStart() - getCodingStop())+1;
	}
	
	public int getGenomicLength()
	{
		return Math.abs(getGenomicStart() - getGenomicStop())+1;
	}
	
	/**
	 * Method returns the phase of the exon
	 * 
	 * @return
	 */
	public int getPhase()
	{
		return ((Math.abs(getCodingStart() - getCodingStop())+1)%3);
	}
	
	/**
	 * Method returns the nucleotide sequence of the exon
	 * 
	 * @return
	 */
	public String getCodingNucleotideSequence()
	{
		if(codingNucleotideSequence != null)
			return codingNucleotideSequence.toString();
		else
			return "";
	}
	
	/**
	 * Method returns true if the exon is coding (coding Sequence != null)
	 * and false otherwise.
	 * 
	 * @return
	 */
	public boolean isCoding()
	{
		if(codingNucleotideSequence != null && codingNucleotideSequence.toString().length() > 0)
			return true;
		else 
			return false;
	}
	
	/**
	 * Method returns the complete nucleotide sequence of the exon (including non-coding parts)
	 * 
	 * @return
	 */
	public String getGenomicNucleotideSequence()
	{
		if(completeNucleotideSequence != null)
			return completeNucleotideSequence.toString();
		else
			return "";
	}
	
	/**
	 * set coding sequence of the exon
	 * 
	 * @param nucleotideSequence
	 */
	public void setCodingNucleotideSequence(String nucleotideSequence) throws Exception
	{
		this.codingNucleotideSequence = new CompressedNucleotideSequence(nucleotideSequence);
	}
	
	/**
	 * set complete sequence of the exon (including non-coding parts)
	 * 
	 * @param nucleotideSequence
	 */
	public void setGenomicNucleotideSequence(String nucleotideSequence) throws Exception
	{
		this.completeNucleotideSequence = new CompressedNucleotideSequence(nucleotideSequence);
	}
	
	/**
	 * Method returns the array of nucleotide instances for the exon
	 * 
	 * @return
	 */
	public Nucleotide[] getNucleotides()
	{
		if(nucleotides != null)
			return nucleotides;
		
		nucleotides = new Nucleotide[this.getLength()];
		
		char[] seq = getCodingNucleotideSequence().toCharArray();
		
		for(int i=0; i<seq.length; i++)
		{
			nucleotides[i] = new Nucleotide(seq[i]);
		}
		
		return nucleotides;
	}
	
	/**
	 * Method returns string representation of exon
	 */
	public String toString()
	{
		StringBuffer buffer = new StringBuffer();
		
		buffer.append(id);
		buffer.append("\t");
		buffer.append(getGenomicNucleotideSequence().length());
		buffer.append("\t");
		buffer.append(getCodingStart());
		buffer.append("\t");
		buffer.append(getCodingStop());
		buffer.append("\t");
		buffer.append(getGenomicStart());
		buffer.append("\t");
		buffer.append(getGenomicStop());
		buffer.append("\t");
		buffer.append(getCodingNucleotideSequence());
		
		return buffer.toString();
	}
	
	/**
	 * Method returns true if one of the exons is fully contained in the other exon
	 * @param other
	 * @return
	 */
	public boolean intersects(Exon other)
	{
		if((getGenomicStart() >= other.getGenomicStart() && getGenomicStop() <= other.getGenomicStop()) || (getGenomicStart() <= other.getGenomicStart() && getGenomicStop() >= other.getGenomicStop()))
			return true;
		else 
			return false;
	}
	
	public boolean overlaps(Exon other)
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
	
	/**
	 * Method returns the nucleotide char at the specified position
	 */
	public String getSequenceElement(int position)
	{
		return "" + getNucleotides()[position].getNucleotide();
	}

	/**
	 * Method returns the length of the exon
	 */
	public int getSequenceLength()
	{
		return getLength();
	}

	/**
	 * Method returns the object at the specified position
	 */
	public Object getSequenceObjectAt(int position)
	{
		return getNucleotides()[position];
	}

	/**
	 * Method returns the string representation of the id
	 */
	public String getID()
	{
		return getID()+"";
	}
		
	public int getIntegerID()
	{
		return this.id;
	}

	@Override
	public String getReference() 
	{
		return this.reference;
	}

	@Override
	public void setReference(String reference) 
	{
		this.reference = reference;
	}
	
	@Override
	/* This method is needed for the SAMBAMCompression method.
	 * (non-Javadoc)
	 * @see biokit.ngs.rnaseq.compression.mapper.MappableElement#getDescription()
	 */
	public String getDescription() 
	{
		return getReference() + ":(" + getCodingStart() + "-" + getCodingStop() + ")";
	}

	@Override
	/* Exons that OVERLAP with this region
	 * are added.
	 */
	public void addMappedExon(Exon e) 
	{
		this.coveredExons.add(e);
	}

	@Override
	/* Transcripts that OVERLAP with this region
	 * are added.
	 */
	public void addMappedTranscript(String transcript) 
	{
		this.coveredTranscripts.add(transcript);
	}

	@Override
	/* Returns exons that OVERLAP with this region.
	 */
	public HashSet<Exon> getMappedExons() 
	{
		return this.coveredExons;
	}

	@Override
	/* Returns transcripts that OVERLAP with this region.
	 */
	public HashSet<String> getMappedTranscripts() 
	{
		return this.coveredTranscripts;
	}
}
