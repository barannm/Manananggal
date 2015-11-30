package BioKit;

import java.util.HashSet;

public class Junction extends GeneticElement implements MappableElement 
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private int junctionID;
	private String reference;
	private HashSet<Exon> surroundingExons;
	private HashSet<String> transcriptsContainingJunction;
	
	public Junction(String reference, int junctionStart, int junctionStop)
	{
		super(junctionStart, junctionStop, junctionStart, junctionStop, true);
		this.reference = reference;
		this.surroundingExons = new HashSet<Exon>();
		this.transcriptsContainingJunction = new HashSet<String>();
	}
		
	@Override
	public String getID() 
	{
		return getDescription();
	}

	@Override
	public int getIntegerID() 
	{
		return junctionID;
	}

	@Override
	public int getCodingStart() 
	{
		return super.getCodingStart();
	}

	@Override
	public int getCodingStop() 
	{
		return super.getCodingStop();
	}

	@Override
	public String getReference() 
	{
		return reference;
	}

	@Override
	public void setReference(String reference) 
	{
		this.reference = reference;
	}
	
	public void setID(int id)
	{
		this.junctionID = id;
	}
		
	@Override
	/* This method is needed for the SAMBAMCompression method.
	 */
	public String getDescription() 
	{
		return reference + ":" + getCodingStart() + "," + reference + ":" + getCodingStop();
	}

	@Override
	public HashSet<Exon> getMappedExons() 
	{

		return this.surroundingExons;
	}

	@Override
	public HashSet<String> getMappedTranscripts() 
	{
		return this.transcriptsContainingJunction;
	}

	@Override
	public void addMappedExon(Exon e) 
	{
		this.surroundingExons.add(e);	
	}

	@Override
	public void addMappedTranscript(String transcript) 
	{
		this.transcriptsContainingJunction.add(transcript);
	}
}

