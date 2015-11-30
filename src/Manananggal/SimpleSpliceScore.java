package Manananggal;

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.TreeSet;

import BioKit.Exon;
import BioKit.Gene;

public class SimpleSpliceScore implements Comparable<SimpleSpliceScore>
{
	String			m_strGeneID;
	String			m_strID;
	int				m_nRating;
	CountElement	m_JunctionInclusion;
	CountElement	m_JunctionExclusion;
	Exon			m_Exon;
	double			m_fPValue;
	double			m_fInclusionChange;
	String			m_strConditionA;
	String			m_strConditionB;
	TreeSet<String> m_vcValidIsoforms;
	boolean 		m_bIsNovel;
	
	SimpleSpliceScore()
	{
		m_strGeneID			= "?";
		m_nRating			= 0;
		m_JunctionInclusion = null;
		m_JunctionExclusion = null;
		m_Exon				= null;
		m_fPValue			= 0.0;
		m_fInclusionChange	= 0.0;
		m_strConditionA	 	= "?";
		m_strConditionB		= "?";
		m_strID				= "NA";
		m_bIsNovel			= false;
		
		m_vcValidIsoforms	= new TreeSet<String>();
	}
	
	SimpleSpliceScore(String strGeneID, CountElement junctionInclusion, CountElement junctionExclusion, Exon exon, double fPValue, double fInclusionChange, String strConditionA, String ConditionB)
	{
		m_strGeneID			= strGeneID;
		m_nRating			= 0;
		m_JunctionInclusion = junctionInclusion;
		m_JunctionExclusion = junctionExclusion;
		m_Exon				= exon;
		m_fPValue			= fPValue;
		m_fInclusionChange	= fInclusionChange;
		m_strConditionA	 	= strConditionA;
		m_strConditionB		= ConditionB;
		m_bIsNovel			= false;
		
		m_vcValidIsoforms	= new TreeSet<String>();
		
		m_strID = m_nRating + "_" + m_Exon.getCodingStart() + "_" + m_Exon.getCodingStop() + junctionInclusion.m_nStart + "_" + junctionInclusion.m_nEnd + "_" + junctionExclusion.m_nStart + "_" + junctionExclusion.m_nEnd + "_" + m_strConditionA + "_" + m_strConditionB + "_" + m_bIsNovel;
	}

	@Override
	public int compareTo(SimpleSpliceScore other)
	{
		int nRes = m_JunctionInclusion.compareTo(other.m_JunctionInclusion);
		
		if(nRes == 0)
			nRes = m_strConditionA.compareTo(other.m_strConditionA);
		
		if(nRes == 0)
			nRes = m_strConditionB.compareTo(other.m_strConditionB);
		
		if(nRes == 0) 
			nRes = m_JunctionExclusion.compareTo(other.m_JunctionExclusion);
		
		if(nRes == 0)
			nRes = m_Exon.compareTo(other.m_Exon);
		
		if(nRes == 0)
		{
			if(m_fPValue < other.m_fPValue)
				return -1;
			else if(m_fPValue > other.m_fPValue)
				return 1;
			else
				return 0;
		}
		
		return nRes;
	}
	
	public String toString()
	{
		String strOut = m_strConditionA + " vs. " + m_strConditionB + " variable_exon=" + m_Exon.getCodingStart() + "-" + m_Exon.getCodingStop();
		strOut += " incl_jun=" + m_JunctionInclusion.m_nStart + "-" + m_JunctionInclusion.m_nEnd;
		strOut += " excl_jun=" + m_JunctionExclusion.m_nStart + "-" + m_JunctionExclusion.m_nEnd;
		strOut += " pValue=" + m_fPValue;
		strOut += " incl_change=" + m_fInclusionChange + " novel=" + m_bIsNovel + "\n";
		
		return strOut;
	}
	
	public void GetValidIsoforms(Gene gene)
	{
		boolean bInclusionKnown = false;
		boolean bExclusionKnown = false;
		for(String strIsoform : gene.getArrayOfGeneProductNames())
		{
			TreeSet<String> strJunctions = gene.getSpliceJunctionInformationForGeneProduct(strIsoform);

			boolean bIncludesInclusionJunction = false;
			boolean bIncludesExclusionJunction = false;
			
			if(strJunctions.contains(m_JunctionInclusion.m_nStart + "-" + m_JunctionInclusion.m_nEnd))
			{
				bIncludesInclusionJunction 	= true;
				bInclusionKnown				= true;
			}
			
			if(strJunctions.contains(m_JunctionExclusion.m_nStart + "-" + m_JunctionExclusion.m_nEnd))
			{
				bIncludesExclusionJunction	= true;
				bExclusionKnown				= true;
			}
			
			if(bIncludesInclusionJunction || bIncludesExclusionJunction)
			{
				m_vcValidIsoforms.add(strIsoform.split("\\.")[0]);
			}
		}
		
		if(!bInclusionKnown || !bExclusionKnown)
			m_bIsNovel = true;
	}
	
	public void WriteToFile(RandomAccessFile pOut)
	{
		try
		{
			pOut.writeUTF(m_strGeneID);
			pOut.writeUTF(m_strID);
			pOut.writeBoolean(m_bIsNovel);
			
			pOut.writeInt(m_nRating);
			
			pOut.writeInt(m_JunctionInclusion.m_nStart);
			pOut.writeInt(m_JunctionInclusion.m_nEnd);
			
			pOut.writeInt(m_JunctionExclusion.m_nStart);
			pOut.writeInt(m_JunctionExclusion.m_nEnd);
			
			pOut.writeInt(m_Exon.getIntegerID());
			pOut.writeInt(m_Exon.getCodingStart());
			pOut.writeInt(m_Exon.getCodingStop());			
			
			pOut.writeDouble(m_fPValue);
			pOut.writeDouble(m_fInclusionChange);
			
			pOut.writeUTF(m_strConditionA);
			pOut.writeUTF(m_strConditionB);
			
			pOut.writeInt(m_vcValidIsoforms.size());
			
			for(String strIsoform : m_vcValidIsoforms)
			{
				pOut.writeUTF(strIsoform);
			}
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}
	
	public boolean ReadFromFile(RandomAccessFile pIn)
	{	
		try
		{
			m_strGeneID 		= pIn.readUTF();
			m_strID				= pIn.readUTF();
			m_bIsNovel			= pIn.readBoolean();
			
			m_nRating			= pIn.readInt();
			
			m_JunctionInclusion = new CountElement();
			m_JunctionInclusion.m_nStart = pIn.readInt();
			m_JunctionInclusion.m_nEnd	 = pIn.readInt();
			
			m_JunctionExclusion = new CountElement();
			m_JunctionExclusion.m_nStart = pIn.readInt();
			m_JunctionExclusion.m_nEnd	 = pIn.readInt();
			
			int nExID			= pIn.readInt();
			int nExStart 		= pIn.readInt();
			int nExEnd			= pIn.readInt();
			m_Exon				= new Exon(nExStart, nExEnd);
			m_Exon.setID(nExID);
			
			m_fPValue			= pIn.readDouble();
			m_fInclusionChange	= pIn.readDouble();
			m_strConditionA	 	= pIn.readUTF();
			m_strConditionB		= pIn.readUTF();
			
			int nIsoforms		= pIn.readInt();
			m_vcValidIsoforms	= new TreeSet<String>();
			for(int i=0; i<nIsoforms; i++)
			{
				m_vcValidIsoforms.add(pIn.readUTF());
			}
			
			m_strID = m_nRating + "_" + m_JunctionInclusion.m_nStart + "_" + m_JunctionInclusion.m_nEnd + "_" + m_JunctionExclusion.m_nStart + "_" + m_JunctionExclusion.m_nEnd + "_" + m_strConditionA + "_" + m_strConditionB + "_" + m_bIsNovel;
		}
		catch (IOException e)
		{
			e.printStackTrace();
			return false;
		}
		return true;
	}
}
