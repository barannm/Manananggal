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
