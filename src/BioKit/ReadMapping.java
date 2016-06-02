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

import java.util.Vector;

/**
 * Biokit representation of a SAM Read mapping
 * 
 * @author Fabian Birzele
 *
 */
public class ReadMapping 
{
	public static String[] ZERO_MAPPING_BIT_MEANING = new String[]{"SINGLE_END", "SEGMENTS_NOT_ALL_PROPERLY_ALIGNED", 
		"SEGMENT_MAPPED", "NEXT_SEGMENT_MAPPED", "ALIGNED_AS_IS", "NEXT_SEGMENT_ALIGNED_AS_IS", "NOT_FIRST_SEGMENT_IN_ALIGNMENT", 
		"NOT_LAST_SEGMENT_IN_ALIGNMENT", "PRIMARY_ALIGNMENT", "PASSING_QC", "NO_DUPLICATE"};
	
	public static String[] ONE_MAPPING_BIT_MEANING = new String[]{"PAIRED_END", "SEGMENTS_ALL_PROPERLY_ALIGNED", 
		"SEGMENT_NOT_MAPPED", "NEXT_SEGMENT_NOT_MAPPED", "SEGMENT_ALIGNED_REVERSE_COMPLEMENT", "NEXT_SEGMENT_ALIGNED_REVERSE_COMPLEMENT", 
		"FIRST_SEGMENT_IN_ALIGNMENT", "LAST_SEGMENT_IN_ALIGNMENT", "SECONDARY_ALIGNMENT", "NOT_PASSING_QC", "PCR_DUPLICATE"};
	
	private String readID;
	private char strand;
	private int mappingBit;
	private int startPosition;
	private int stopPosition;
	private String reference;
	private String referenceOfPair;
	private String readSequence;
	private String qualitySequence;
	private String mappingCIGAR;
	private float readMappingQuality;
	private int startPositionOfPairedEnd;
	private String readIDOfPairedEnd;
	private String featureString;
	private char[] bits;
	
	public ReadMapping(String readID, int mappingBit, char strand, int startPosition, int stopPosition, String reference, String referenceOfPair, String readSequence,
		String qualitySequence, String mappingCIGAR, float readMappingQuality) 
	{
		this.readID = readID;
		this.mappingBit = mappingBit;
		this.strand = strand;
		this.startPosition = startPosition;
		this.stopPosition = stopPosition;
		this.reference = reference;
		this.readSequence = readSequence;
		this.qualitySequence = qualitySequence;
		this.mappingCIGAR = mappingCIGAR;
		this.readMappingQuality = readMappingQuality;
		this.readIDOfPairedEnd = "*";
		this.startPositionOfPairedEnd = 0;
		this.featureString = "";
		this.referenceOfPair = referenceOfPair;
		
		bits = Integer.toBinaryString(mappingBit).toCharArray();
	}
	
	public void setMappingBit(boolean pairedEnd, boolean allSegmentsProperlyAligned, boolean segmentNotMapped, boolean nextSegmentNotMapped, boolean alignedReverseComplement, boolean nextSegmentAlignedReverseComplement, 
		boolean firstSegmentInAlignment, boolean lastSegmentInAlignment, boolean secondaryAlignment, boolean notPassingQC, boolean pcrDuplicate)
	{
		StringBuffer b = new StringBuffer();
		
		if(pcrDuplicate)
			b.append(1);
		else
			b.append(0);
		
		if(notPassingQC)
			b.append(1);
		else
			b.append(0);
		
		if(secondaryAlignment)
			b.append(1);
		else
			b.append(0);
		
		if(lastSegmentInAlignment)
			b.append(1);
		else
			b.append(0);
		
		if(firstSegmentInAlignment)
			b.append(1);
		else
			b.append(0);
		
		if(nextSegmentAlignedReverseComplement)
			b.append(1);
		else
			b.append(0);
		
		if(alignedReverseComplement)
			b.append(1);
		else
			b.append(0);
		
		if(nextSegmentNotMapped)
			b.append(1);
		else
			b.append(0);
		
		if(segmentNotMapped)
			b.append(1);
		else
			b.append(0);
		
		if(allSegmentsProperlyAligned)
			b.append(1);
		else
			b.append(0);
		
		if(pairedEnd)
			b.append(1);
		else
			b.append(0);
		
		mappingBit = Integer.parseInt(b.toString(), 2);
		bits = Integer.toBinaryString(mappingBit).toCharArray();
	}
	
	public String explainMappingBit()
	{
		StringBuffer b = new StringBuffer();
		
		int position = 0;
		
		for(int i=bits.length-1; i>=0; i--)
		{
			if(bits[i] == '1')
				b.append(ONE_MAPPING_BIT_MEANING[position] + " ");
			else
				b.append(ZERO_MAPPING_BIT_MEANING[position] + " ");
			
			position++;
		}
		
		return b.toString();
	}
	
	public String getReferenceOfPair()
	{
		return referenceOfPair;
	}
	
	public void setReferenceOfPair(String pairReference)
	{
		referenceOfPair = pairReference;
	}
	
	public boolean pairedEndRead()
	{
		try
		{
			if(bits[bits.length-1] == '1')
				return true;
		
			return false;
		}
		catch(Exception ex)
		{
			return false;
		}
	}
	
	public boolean allReadsProperlyAligned()
	{
		try
		{
			if(bits[bits.length-2] == '1')
				return true;
			
			return false;
		}
		catch(Exception ex)
		{
			return false;
		}
	}
	
	public boolean isUnmapped()
	{
		try
		{
			if(bits[bits.length-3] == '1')
				return true;
		}
		catch(Exception ex)
		{
			return false;
		}
		
		return false;
	}
	
	public boolean pairUnmapped()
	{
		try
		{
			if(bits[bits.length-4] == '1')
				return true;
		}
		catch(Exception ex)
		{
			return false;
		}
		
		return false;
	}
	
	public boolean readAlignedAsReverseComplement()
	{
		try
		{
			if(bits[bits.length-5] == '1')
				return true;
		}
		catch(Exception ex)
		{
			return false;
		}
		
		return false;
	}
	
	public boolean pairAlignedAsReverseComplement()
	{
		try
		{
			if(bits[bits.length-6] == '1')
				return true;
		}
		catch(Exception ex)
		{
			return false;
		}
		
		return false;
	}
	
	public boolean firstSegmentInAlignment()
	{
		try
		{
			if(bits[bits.length-7] == '1')
				return true;
		}
		catch(Exception ex)
		{
			return false;
		}
		
		return false;
	}
	
	public boolean lastSegmentInAlignment()
	{
		try
		{
			if(bits[bits.length-8] == '1')
				return true;
		}
		catch(Exception ex)
		{
			return false;
		}
		
		return false;
	}
	
	public boolean secondaryAlignment()
	{
		try
		{
			if(bits[bits.length-9] == '1')
				return true;
		}
		catch(Exception ex)
		{
			return false;
		}
		
		return false;
	}
	
	public boolean alignmentNotPassingQC()
	{
		try
		{
			if(bits[bits.length-10] == '1')
				return true;
		}
		catch(Exception ex)
		{
			return false;
		}
		
		return false;
	}
	
	public boolean pcrDuplicate()
	{
		try
		{
			if(bits[bits.length-11] == '1')
				return true;
		}
		catch(Exception ex)
		{
			return false;
		}
		
		return false;
	}
	
	public void setReadID(String readID)
	{
		this.readID = readID;
	}
	
	public int getStartPositionOfPairedEnd() 
	{
		return startPositionOfPairedEnd;
	}

	public void setStartPositionOfPairedEnd(int startPositionOfPairedEnd)
	{
		this.startPositionOfPairedEnd = startPositionOfPairedEnd;
	}

	public String getReadIDOfPairedEnd() 
	{
		return readIDOfPairedEnd;
	}

	public void setReadIDOfPairedEnd(String readIDOfPairedEnd)
	{
		this.readIDOfPairedEnd = readIDOfPairedEnd;
	}

	public void setFeatureString(String feature)
	{
		this.featureString = feature;
	}
	
	public String getFeatureString()
	{
		return featureString;
	}
	
	public void addFeature(String feature)
	{
		if(featureString.equals(""))
			featureString += feature;
		else
			featureString += "\t" + feature;
	}
	
	public String toString()
	{
		if(featureString.length() != 0)
			return readID + "\t" + mappingBit + "\t" + reference + "\t" + startPosition + "\t" + (int)readMappingQuality +"\t" + mappingCIGAR + "\t"+referenceOfPair+"\t"+startPositionOfPairedEnd+"\t0\t" + readSequence + "\t" + qualitySequence + "\t" + featureString;
		else
			return readID + "\t" + mappingBit + "\t" + reference + "\t" + startPosition + "\t" + (int)readMappingQuality +"\t" + mappingCIGAR + "\t"+referenceOfPair+"\t"+startPositionOfPairedEnd+"\t0\t" + readSequence + "\t" + qualitySequence;
	}
	
	public String getReadID() 
	{
		return readID;
	}
	
	public char getStrand() 
	{
		return strand;
	}
	
	public int getStartPosition() 
	{
		return startPosition;
	}
	
	public int getStopPosition()
	{
		return stopPosition;
	}
	
	public boolean isSplicedRead()
	{
		if(mappingCIGAR.contains("N"))
			return true;
		else
			return false;
	}
	
	/**
	 * Method returns all matched reference segments that are matched by the read.
	 * 
	 * They are returned as a vector of strings that represent: start1-stop1, start2-stop2, ..., startN-stopN
	 * 
	 * @return
	 */
	public Vector<String> getMatchedReferenceSegments()
	{
		Vector<String> groups = new Vector<String>();
		
		Cigar cigar = new Cigar(mappingCIGAR);
		
		int currentPosition = startPosition;
		
		for(int i=0; i<cigar.getNumberOfGroups(); i++)
		{
			if(cigar.getSymbolOfGroup(i) == Cigar.ANY_MATCH)
			{
				groups.add(currentPosition + "-"  + (currentPosition+cigar.getLengthOfGroup(i)-1));
				
				currentPosition += cigar.getLengthOfGroup(i);
			}
			else if(cigar.getSymbolOfGroup(i) == Cigar.SKIP_IN_REFERNCE)
			{
				currentPosition += cigar.getLengthOfGroup(i);
			}
		}
		
		return groups;
	}
	
	public String getReference() 
	{
		return reference;
	}
	
	public String getReadSequence()
	{
		return readSequence;
	}
	
	public String getQualitySequence()
	{
		return qualitySequence;
	}
	
	public String getMappingCIGARString()
	{
		return mappingCIGAR;
	}
	
	public Cigar getMappingCIGAR()
	{
		return new Cigar(mappingCIGAR);
	}
	
	public float getReadMappingQuality()
	{
		return readMappingQuality;
	}
	
	public int getMappingBit()
	{
		return mappingBit;
	}
	
	public String getFeatureValue(String featureTag)
	{
		String[] split = featureString.split("\\s+");
		
		for(String f : split)
		{
			if(f.startsWith(featureTag))
				return f.split(":")[2];
		}
		
		return null;
	}
	
	public static void main(String[] args)
	{
		ReadMapping m = new ReadMapping("", 153, '+', 0, 0, "bla", "bla", "ABC", "ABC", "3M", 255);
		System.out.println("153: " + Integer.toBinaryString(153) + " " + m.explainMappingBit());
		m = new ReadMapping("", 101, '+', 0, 0, "bla", "bla", "ABC", "ABC", "3M", 255);
		System.out.println("101: " + Integer.toBinaryString(101) + " " + m.explainMappingBit());
	}
}
