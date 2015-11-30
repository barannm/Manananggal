package BioKit;

import java.io.IOException;
import java.io.Serializable;
import java.util.HashSet;

import BioKit.Sequence;
import BioKit.MappableElement;


/**
 * This class represents a exon wrapping start, stop and sequence information.
 * 
 * @author Fabian Birzele
 *
 */
public class Exon extends GeneticElement implements MappableElement, Sequence, Serializable
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private String ensemblExonID;
	
	private int id;
	private String reference;
	private CompressedNucleotideSequence codingNucleotideSequence;
	private CompressedNucleotideSequence completeNucleotideSequence;
	private Nucleotide[] nucleotides;
	private HashSet<String> coveredTranscripts;
	
	/**
	 * Standard constructor, initializes the exon
	 * 
	 * @param codingStart
	 * @param codingStop
	 * @param nucleotideSequence
	 * @throws Exception
	 */
	public Exon(int codingStart, int codingStop, String codingNucleotideSequence) throws IOException
	{
		super(codingStart, codingStop, codingStart, codingStop, true);
		this.codingNucleotideSequence = new CompressedNucleotideSequence(codingNucleotideSequence);
		this.coveredTranscripts = new HashSet<String>();
	}
	
	public Exon(int codingStart, int codingStop)
	{
		super(codingStart, codingStop, codingStart, codingStop, true);
		this.coveredTranscripts = new HashSet<String>();
	}
	
	public Exon(int codingStart, int codingStop, int genomicStart, int genomicStop)
	{
		super(genomicStart, genomicStop, codingStart, codingStop, true);
		this.coveredTranscripts = new HashSet<String>();
	}
	
	public Exon(int codingStart, int codingStop, int genomicStart, int genomicStop, String codingNucleotideSequence) throws IOException, Exception
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
	 * Method sets the ensembl exon id
	 * 
	 * @param exonID
	 */
	public void setEnsemblExonID(String exonID)
	{
		this.ensemblExonID = exonID;
	}
	
	/**
	 * Method returns the ensembl exon id
	 * 
	 * @return
	 */
	public String getEnsemblExonID()
	{
		return ensemblExonID;
	}
	
	/**
	 * Method returns exon id
	 * @return
	 */
	public int getExonID()
	{
		return id;
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
		buffer.append(ensemblExonID);
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
		
		int nStart	= getGenomicStart();
		int nEnd	= getGenomicStop();
		int nOtherStart = other.getGenomicStart();
		int nOtherEnd	= other.getGenomicStop();
		
		// check orientation
		boolean bThisIsForward	= false;
		boolean bOtherIsForward	= false;
		
		if(getGenomicStart() < getGenomicStop())
			bThisIsForward = true;
		
		if(other.getGenomicStart() < other.getGenomicStop())
			bOtherIsForward = true;
		
		if(bThisIsForward && bOtherIsForward)
		{
			return(nStart <= nOtherEnd && nEnd >= nOtherStart);
		}
		else if(bThisIsForward && !bOtherIsForward)
		{
			return(nStart <= nOtherStart && nEnd >= nOtherEnd);
		}
		else if(!bThisIsForward && bOtherIsForward)
		{
			return(nEnd <= nOtherEnd && nStart >= nOtherStart);
		}
		else // both on second strand
		{
			return(nEnd <= nOtherStart && nStart >= nOtherEnd);
		}
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
		return getExonID()+"";
	}
	
	@Override
    public boolean equals(Object object)
    {
	    if(object instanceof Exon)
	    {
	        Exon other = (Exon)object;
	        
	        if(this.getElementID() == other.getElementID() && this.getExonID() == other.getExonID())
	        	return true;
	    }
	    return false;
    }
	
	@Override
	public int getIntegerID() 
	{
		return this.getElementID();
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
	public void addMappedExon(Exon e) 
	{
	}

	@Override
	/* Transcripts that CONTAIN this exon are added.
	 */
	public void addMappedTranscript(String transcript) 
	{
		this.coveredTranscripts.add(transcript);
	}

	@Override
	public HashSet<Exon> getMappedExons() 
	{
		HashSet<Exon> result = new HashSet<Exon>();
		result.add(this);
		return result;
	}

	@Override
	/* Returns transcripts that CONTAIN this exon.
	 */
	public HashSet<String> getMappedTranscripts() 
	{
		return this.coveredTranscripts;
	}
}
