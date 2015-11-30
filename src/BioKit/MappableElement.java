package BioKit;

import java.util.HashSet;

public interface MappableElement 
{
	public String getID();
	public int getIntegerID();
	public int getCodingStart();
	public int getCodingStop();
	public String getReference();
	public void setReference(String reference);
	public String getDescription();
	public void addMappedExon(Exon e);
	public void addMappedTranscript(String transcript);
	public HashSet<Exon> getMappedExons();
	public HashSet<String> getMappedTranscripts();
}
