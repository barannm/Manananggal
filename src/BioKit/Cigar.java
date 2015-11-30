package BioKit;

import java.util.Vector;

/**
 * Biokit representation of a Cigar string. Provides some additional functionality not available e.g. in Picard Cigar representation
 * 
 * @author Fabian Birzele
 *
 */
public class Cigar
{
	public static final char ANY_MATCH = 'M';
	public static final char DELETION_IN_REFERENCE = 'I';
	public static final char INSERTION_IN_REFERENCE = 'D';
	public static final char MATCH = '=';
	public static final char MISMATCH = 'X';
	public static final char SOFT_CLIP = 'S';
	public static final char HARD_CLIP = 'H';
	public static final char PADDING = 'P';
	public static final char SKIP_IN_REFERNCE = 'N';
	
	private String originalCigar;
	
	private Vector<Integer> lengths;
	private Vector<Character> symbols;
	
	public Cigar(String cigar)
	{
		this.originalCigar = cigar;
		
		parseCigar(cigar);
	}
	
	public String getOriginalCigar()
	{
		return originalCigar;
	}
	
	public int getNumberOfGroups()
	{
		return lengths.size();
	}
	
	public int getLengthOfGroup(int group)
	{
		return lengths.get(group);
	}
	
	public char getSymbolOfGroup(int group)
	{
		return symbols.get(group);
	}
	
	public Cigar reverseSoftClip()
	{
		StringBuffer b = new StringBuffer();
		
		boolean startsWithSoftclip = false;
		boolean endsWithSoftclip = false;
		String startSoftClip = null;
		String stopSoftClip = null;
		
		String value = null;
		
		for(int i=0; i<symbols.size(); i++)
		{
			value = lengths.get(i) + "" + symbols.get(i);
			
			if(i == 0 && symbols.get(i).equals(SOFT_CLIP))
			{
				startSoftClip = value;
				startsWithSoftclip = true;
			}
			else if(i == symbols.size()-1 && symbols.get(i).equals(SOFT_CLIP))
			{
				stopSoftClip = value;
				endsWithSoftclip = true;
			}
			else 
			{
				b.append(value);
			}
		}
		
		value = b.toString();
		
		if(startsWithSoftclip)
			value = value + startSoftClip;
		
		if(endsWithSoftclip)
			value = stopSoftClip + value;
		
		return new Cigar(value);
	}
	
	public int getCigarGroupLength()
	{
		int length = 0;
		
        for(int i=0;i<symbols.size();i++) 
        {
            length += lengths.get(i);
        }
        
        return length;
	}
	
	public int getReferencePositionsCovered()
	{
		int length = 0;
        for(int i=0;i<symbols.size();i++) 
        {
            switch ((char)symbols.get(i)) 
            {
                case ANY_MATCH:
                case INSERTION_IN_REFERENCE:
                case SKIP_IN_REFERNCE:
                case MATCH:
                case MISMATCH:
                    length += lengths.get(i);
            }
        }
        return length;

	}
	
	public int getReadPositionsCovered()
	{
		int length = 0;
        for(int i=0;i<symbols.size();i++) 
        {
            switch ((char)symbols.get(i)) 
            {
                case ANY_MATCH:
                case DELETION_IN_REFERENCE:
                case MATCH:
                case MISMATCH:
                case SOFT_CLIP:
                    length += lengths.get(i);
            }
        }
        
        return length;
	}
	
	public String cigarSubstring(int referencePositionStart, int referencePositionStop)
	{
		StringBuffer b = new StringBuffer();
		
		char[] cigar = getCigarInformationPerPosition();
		
		int startInCigar = 0;
		int coveredReferencePositionCounter = 0;
		
		for(int i=0; i<cigar.length; i++)
		{
			startInCigar = i;

			if(coveredReferencePositionCounter == referencePositionStart)
				break;
			
			switch((char)cigar[i]) 
            {
                case ANY_MATCH:
                case INSERTION_IN_REFERENCE:
                case SKIP_IN_REFERNCE:
                case MATCH:
                case MISMATCH:
                    coveredReferencePositionCounter++;
            }
		}
		
		//System.out.println("Request: " + referencePositionStart + " - " + referencePositionStop);
		//System.out.println("Start in Cigar: " + startInCigar);
		
		int positionInCigar = startInCigar-1;
		int positionInReference = referencePositionStart;
		
		char currentGroup = '-';
		int currentLength = 0;
		
		//StringBuffer positionsVisited = new StringBuffer();
		
		while(positionInReference <= referencePositionStop)
		{
			positionInCigar++;
			
			if(currentGroup != cigar[positionInCigar] && currentGroup != '-')
			{
				//System.out.println("\tAdded group " + currentLength + "" + currentGroup);
				b.append(currentLength + "" + currentGroup);
				currentGroup = '-';
				currentLength = 0;
			}
			
			//positionsVisited.append(cigar[positionInCigar]);
			
			currentGroup = cigar[positionInCigar];
			currentLength++;
			
			if(cigar[positionInCigar] == DELETION_IN_REFERENCE)
			{
				//System.out.println("Deletion in reference at position " + positionInCigar);
				continue;
			}
			else if(cigar[positionInCigar] == HARD_CLIP)
			{
				continue;
			}
			else if(cigar[positionInCigar] == SOFT_CLIP)
			{
				continue;
			}
			
			positionInReference++;
		}
		
		b.append(currentLength + "" + currentGroup);
		//System.out.println("\tAdded group " + currentLength + "" + currentGroup);
		//System.out.println("\tPositions visited " + positionsVisited.length() + " " + positionsVisited.toString());
		
		return b.toString();
	}
	
	public char[] getCigarInformationPerPosition()
	{
		char[] array = new char[getCigarGroupLength()];
		
		//StringBuffer positionsVisited = new StringBuffer();
		
		int counter = 0;
		
		for(int i=0;i<symbols.size();i++) 
        {
			for(int j=0; j<lengths.get(i); j++)
			{
				array[counter] = symbols.get(i);
				//positionsVisited.append(symbols.get(i));
				counter++;
			}
        }
		
		//System.out.println("Cigar: " + getOriginalCigar());
		//System.out.println("Position array: " + positionsVisited);
		
		return array;
	}
	
	public int getNumberOfReadPositionsAlignedWithReference()
	{
		int length = 0;
        for(int i=0;i<symbols.size();i++) 
        {
            switch ((char)symbols.get(i)) 
            {
                case ANY_MATCH:
                case MATCH:
                case MISMATCH:
                    length += lengths.get(i);
            }
        }
        return length;
	}
	
	private void parseCigar(String cigar)
	{
		lengths = new Vector<Integer>();
		symbols = new Vector<Character>();
		
		StringBuffer currentValue = new StringBuffer();
		boolean lastWasDigit = false;
		
		for(char c : cigar.toCharArray())
		{
			if(Character.isDigit(c))
			{
				if(currentValue.length() != 0 && !lastWasDigit)
				{
					symbols.add(currentValue.toString().charAt(0));
					currentValue = new StringBuffer();
				}
				
				currentValue.append(c);
				lastWasDigit = true;
			}
			else
			{
				if(currentValue.length() != 0 && lastWasDigit)
				{
					lengths.add(Integer.parseInt(currentValue.toString()));
					currentValue = new StringBuffer();
				}
			
				currentValue.append(c);
				lastWasDigit = false;
			}
		}
		
		if(currentValue.length() != 0 && !lastWasDigit)
		{
			symbols.add(currentValue.toString().charAt(0));
		}
	}
	
	public static void main(String[] args)
	{
		Cigar cigar = new Cigar("10M5I5M10D20M");
		
		System.out.println("Cigar: " + cigar.getOriginalCigar());
		System.out.println("Ref 0 - 14: " + cigar.cigarSubstring(0, 15));
		System.out.println("Ref 5 - 19: " + cigar.cigarSubstring(5, 20));
		
		Cigar c = new Cigar("1S73M3N4M2S");
		
		System.out.println(c.reverseSoftClip().getOriginalCigar());
	}
}
