package BioKit;

import java.util.HashMap;
import java.util.Vector;

public class GenomeReadMappingReport 
{
	private String chromosome;
	private int regionStart;
	private int regionStop;
	private int numberOfReads;
	private int readLength;
	
	private HashMap<Integer, Integer> genomePositionToCoverageCount;
	private HashMap<String, Integer> spliceJunctionToReadCount;
	private HashMap<String, Integer> connectedExonsToReadCount;
	
	public GenomeReadMappingReport(String chromosome, int regionStart, int regionStop)
	{
		this.chromosome = chromosome;
		this.regionStart = regionStart;
		this.regionStop = regionStop;
		
		numberOfReads = 0;
		
		genomePositionToCoverageCount 	= new HashMap<Integer, Integer>();
		spliceJunctionToReadCount 		= new HashMap<String, Integer>();
		connectedExonsToReadCount		= new HashMap<String, Integer>();
	}
	
	public void setReadLength(int readLength)
	{
		this.readLength = readLength;
	}
	
	public int getReadLength()
	{
		return readLength;
	}
	
	public void increaseNumberOfReads()
	{
		numberOfReads++;
	}
	
	public int getNumberOfReads()
	{
		return numberOfReads;
	}
	
	public void increaseCountAtPosition(int position)
	{
		if(genomePositionToCoverageCount.containsKey(position))
			genomePositionToCoverageCount.put(position, genomePositionToCoverageCount.get(position) + 1);
		else
			genomePositionToCoverageCount.put(position, 1);
	}
	
	public void increaseCountsInRange(int startPosition, int stopPosition)
	{
		for(int i=startPosition; i<=stopPosition; i++)
		{
			if(genomePositionToCoverageCount.containsKey(i))
				genomePositionToCoverageCount.put(i, genomePositionToCoverageCount.get(i) + 1);
			else
				genomePositionToCoverageCount.put(i, 1);
		}
		
	}
	
	public void increaseCountsInRangeByValue(int startPosition, int stopPosition, int value)
	{
		for(int i=startPosition; i<=stopPosition; i++)
		{
			if(genomePositionToCoverageCount.containsKey(i))
				genomePositionToCoverageCount.put(i, genomePositionToCoverageCount.get(i) + value);
			else
				genomePositionToCoverageCount.put(i, value);
		}
		
	}
	
	public void addSpliceJunction(int lastPositionFirstExon, int firstPositionSecondExon)
	{
		// make sure that id works always from the left to the right
		if(lastPositionFirstExon > firstPositionSecondExon)
		{
			int temp = lastPositionFirstExon;
			lastPositionFirstExon = firstPositionSecondExon;
			firstPositionSecondExon = temp;
		}
		
		String id = lastPositionFirstExon + "-" + firstPositionSecondExon;
		
		if(spliceJunctionToReadCount.containsKey(id))
			spliceJunctionToReadCount.put(id, spliceJunctionToReadCount.get(id)+1);
		else
			spliceJunctionToReadCount.put(id, 1);
	}
	
	public void addConnectedExons(int nExStart1, int nExEnd1, int nExStart2, int nExEnd2)
	{
		String id = null;
		if(nExStart1 < nExStart2)
			id = nExStart1 + "_" + nExEnd1 + "_" + nExStart2 + "_" + nExEnd2;
		else
			id = nExStart2 + "_" + nExEnd2 + "_" + nExStart1 + "_" + nExEnd1;
		
		if(connectedExonsToReadCount.containsKey(id))
			connectedExonsToReadCount.put(id, spliceJunctionToReadCount.get(id)+1);
		else
			connectedExonsToReadCount.put(id, 1);
	}
	
	int getReadCountForConnectedExons(String name)
	{
		if(connectedExonsToReadCount.containsKey(name))
			return connectedExonsToReadCount.get(name);
		else
			return 0;
	}
	
	public HashMap<String, Integer> getAllReadCountsForConnectedExons()
	{
		return connectedExonsToReadCount;
	}

	public String getChromosome() 
	{
		return chromosome;
	}

	public int getRegionStart()
	{
		return regionStart;
	}

	public int getRegionStop() 
	{
		return regionStop;
	}
	
	public int getReadCountAtPosition(int position)
	{
		if(genomePositionToCoverageCount.containsKey(position))
			return genomePositionToCoverageCount.get(position);
		else
			return 0;
	}
	
	public Vector<String> getSpliceJunctionNames()
	{
		Vector<String> v = new Vector<String>();
		v.addAll(spliceJunctionToReadCount.keySet());
		
		return v;
	}
	
	public int getReadCountForSpliceJunction(String name)
	{
		if(spliceJunctionToReadCount.containsKey(name))
			return spliceJunctionToReadCount.get(name);
		else
			return 0;
	}
	
	public int getNumberOfSpliceJunctionsWithCoverage()
	{
		return spliceJunctionToReadCount.size();
	}
	
	public int getNumberOfCoveredPositions()
	{
		return genomePositionToCoverageCount.size();
	}
	
	public int getReadCountForSpliceJunction(int lastPositionFirstExon, int firstPositionSecondExon)
	{
		// make sure that id works always from the left to the right
		if(lastPositionFirstExon > firstPositionSecondExon)
		{
			int temp = lastPositionFirstExon;
			lastPositionFirstExon = firstPositionSecondExon;
			firstPositionSecondExon = temp;
		}
		
		String id = lastPositionFirstExon + "-" + firstPositionSecondExon;
		
		if(spliceJunctionToReadCount.containsKey(id))
			return spliceJunctionToReadCount.get(id);
		else
			return 0;
		
	}
}
