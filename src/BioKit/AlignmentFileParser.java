package BioKit;

import java.util.Iterator;

import net.sf.samtools.SAMRecord;

/**
 * Parser for SAM alignment format and to convert Picard SAM records to Biokit ones.
 * 
 * @author Fabian Birzele
 */
public class AlignmentFileParser 
{
	public static ReadMapping convertPicardSAMRecord(SAMRecord samRecord)
	{
		char strand = '+';
		
		if(samRecord.getReadNegativeStrandFlag())
			strand = '-';	
		
		ReadMapping mapping = new ReadMapping(samRecord.getReadName(), samRecord.getFlags(), strand, 
				samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), samRecord.getReferenceName(), samRecord.getMateReferenceName(),
				samRecord.getReadString(), samRecord.getBaseQualityString(), samRecord.getCigarString(), 
				samRecord.getMappingQuality());
		
		mapping.setStartPositionOfPairedEnd(samRecord.getMateAlignmentStart());
		
		Iterator<SAMRecord.SAMTagAndValue> iter = samRecord.getAttributes().iterator();
		
		while(iter.hasNext())
		{
			SAMRecord.SAMTagAndValue value = iter.next();
			
			if(value.tag.equals("MD"))
				mapping.addFeature(value.tag + ":Z:" + value.value);
			else if(value.tag.equals("XA"))
				mapping.addFeature(value.tag + ":i:" + value.value);
			else if(value.tag.equals("NM"))
				mapping.addFeature(value.tag + ":i:" + value.value);
			else
				mapping.addFeature(value.tag + ":Z:" + value.value);
		}
		
		
		
		return mapping;
	}
	
	public static ReadMapping parseSAMFormat(String line, boolean returnNullForUnmapped)
	{
		if(line.startsWith("@"))
			return null;
		
		String[] split = line.split("\t");
		
		try
		{
			String readID = split[0];
				
			int mappingFlag = Integer.parseInt(split[1]);
				
			// Bit 004 set indicates an unmapped read
			if(returnNullForUnmapped)
			{
				if(mappingFlag == 4)
					return null;
				else if(mappingFlag == 77)
					return null;
				else if(mappingFlag == 141)
					return null;
			}
				
			char strand;
			
			// Bit 010 indicates a reverse complemented mapping
			if(mappingFlag == 16)
				strand = '-';
			else if(mappingFlag == 147)
				strand = '-';
			else if(mappingFlag == 83)
				strand = '-';
			else
				strand = '+';
				
			// only use the first part of a header as id. Split at any space, tab... 
			String referenceName = split[2].trim().split("\\s+")[0];
			
			if(referenceName.startsWith("chr"))
				referenceName = referenceName.replace("chr", "");
			else if(referenceName.startsWith("NC_"))
				referenceName = translateHumanNCBIChromosomeID(referenceName);
			
			String mateReferenceName = split[6];
			
			if(mateReferenceName.startsWith("chr"))
				mateReferenceName = mateReferenceName.replace("chr", "");
			else if(mateReferenceName.startsWith("NC_"))
				mateReferenceName = translateHumanNCBIChromosomeID(mateReferenceName);
			
			String sequence = split[9];
			String quality = split[10];
			Cigar cigar = new Cigar(split[5]);
			
			int startPosition = Integer.parseInt(split[3]);
			int stopPosition = startPosition + cigar.getReferencePositionsCovered() - 1;
			int mapQ = Integer.parseInt(split[4]);
			int startPositionOfPairedEnd = Integer.parseInt(split[7]);
			
			ReadMapping mapping = new ReadMapping(readID, mappingFlag, strand, startPosition, stopPosition, referenceName, 
				mateReferenceName, sequence, quality, cigar.getOriginalCigar(), mapQ);
			mapping.setStartPositionOfPairedEnd(startPositionOfPairedEnd);
			
			for(int i=11; i<split.length; i++)
			{
				if(split[i].startsWith("MD"))
					mapping.addFeature(split[i]);
				else if(split[i].startsWith("NM"))
					mapping.addFeature(split[i]);
			}
			
			return mapping;
		}
		catch(Exception e)
		{
			if(split.length > 0)
				System.out.println("Read " + split[0] + " invalid line: " + line);
			else
				System.out.println("invalid line: " + line);
		}
		
		return null;
	}
	
	private static String translateHumanNCBIChromosomeID(String matchID)
	{
		String id = matchID.replaceAll(".+NC_0+", "");
		
		int value = Integer.parseInt(id.split("\\.")[0]);
		
		if(value <= 22)
			return value + "";
		else if(value == 23)
			return "X";
		else if(value == 24)
			return "Y";
		else if(value == 1807)
			return "MT";
		
		return "unknown";
	}

}
