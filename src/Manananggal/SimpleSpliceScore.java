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
 * @author Matthias Barann
 */
package Manananggal;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.TreeSet;

import BioKit.Exon;
import BioKit.Gene;

/**
 *    This class stores all values that are associated with
 *    a PSI score result (e.g. the inclusion and exclusion
 *    junction and p-value).
 */
public class SimpleSpliceScore implements Comparable<SimpleSpliceScore>
{
	String			m_strGeneID;
	String			m_strID;
	CountElement	m_JunctionInclusion;
	CountElement	m_JunctionExclusion;
	double			m_fPValue;
	boolean			m_bSignificant;
	double			m_fInclusionChange;
	String			m_strConditionA;
	String			m_strConditionB;
	boolean 		m_bIsNovel;
	int				m_nType;
	
	SimpleSpliceScore()
	{
		m_strGeneID			= "?";
		m_JunctionInclusion = null;
		m_JunctionExclusion = null;
		m_fPValue			= 0.0;
		m_bSignificant		= false;
		m_fInclusionChange	= 0.0;
		m_strConditionA	 	= "?";
		m_strConditionB		= "?";
		m_strID				= "NA";
		m_bIsNovel			= false;
		m_nType				= -1;
	}
	
	SimpleSpliceScore(String strGeneID, CountElement junctionInclusion, CountElement junctionExclusion, double fPValue, double fInclusionChange, String strConditionA, String ConditionB, int nType)
	{
		m_strGeneID			= strGeneID;
		m_JunctionInclusion = junctionInclusion;
		m_JunctionExclusion = junctionExclusion;
		m_fPValue			= fPValue;
		m_bSignificant		= false;
		m_fInclusionChange	= fInclusionChange;
		m_strConditionA	 	= strConditionA;
		m_strConditionB		= ConditionB;
		m_bIsNovel			= false;
		m_nType				= nType;
		
		m_strID = junctionInclusion.m_nStart + "_" + junctionInclusion.m_nEnd + "_" + junctionExclusion.m_nStart + "_" + junctionExclusion.m_nEnd + "_" + m_strConditionA + "_" + m_strConditionB + "_" + m_bIsNovel + "_" + m_nType;
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
	
	public int GetType()
	{
		return m_nType;
	}
	
	public String GetTypeAsString()
	{
		String strType = "?";
		switch(m_nType)
		{
			case SplicingWebApp.AS_TYPE_EXON_SKIPPING: 			strType = "exn_skipping"; 			break;
			case SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN: 	strType = "alt_start_unique_jun"; 	break;
			case SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN: 	strType = "alt_end_unique_jun"; 	break;
			case SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN: 	strType = "alt_start_shared_jun"; 	break;
			case SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN: 	strType = "alt_end_shared_jun"; 	break;
			case SplicingWebApp.AS_TYPE_RETAINED_INTRON:		strType = "retained_intron";		break;
			case SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN_DOUBLE: 	strType = "alt_start_unique_jun_double"; 	break;
			case SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN_DOUBLE: 		strType = "alt_end_unique_jun_double"; 		break;
			case SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN_DOUBLE: 	strType = "alt_start_shared_jun_double"; 	break;
			case SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN_DOUBLE: 		strType = "alt_end_shared_jun_double"; 		break;
			case SplicingWebApp.AS_TYPE_ALT_3_PRIME_EXON_END:			strType = "alt_3_prime_exon_end";			break;
			case SplicingWebApp.AS_TYPE_ALT_5_PRIME_EXON_END:			strType = "alt_5_prime_exon_end";			break;
			default: strType = "" + m_nType; 			break;
		}
		
		return strType;
	}
	
	public String toString()
	{
		String strOut = m_strConditionA + " vs. " + m_strConditionB + " type=" + m_nType;
		strOut += " incl_jun=" + m_JunctionInclusion.m_nStart + "-" + m_JunctionInclusion.m_nEnd;
		strOut += " excl_jun=" + m_JunctionExclusion.m_nStart + "-" + m_JunctionExclusion.m_nEnd;
		strOut += " pValue=" + m_fPValue + " significant=" + m_bSignificant;
		strOut += " incl_change=" + m_fInclusionChange + " novel=" + m_bIsNovel + "\n";
		
		return strOut;
	}
	
	/** 
	 *    Checks whether the result includes any unannotated junctions and,
	 *    thus, is referring to a novel splicing event
	 */
	public void CheckIfNovel(Gene gene)
	{
		boolean bInclusionKnown = false;
		boolean bExclusionKnown = false;
		for(String strIsoform : gene.getArrayOfGeneProductNames())
		{
			TreeSet<String> strJunctions = gene.getSpliceJunctionInformationForGeneProduct(strIsoform);
			
			if(strJunctions.contains(m_JunctionInclusion.m_nStart + "-" + m_JunctionInclusion.m_nEnd))
			{
				bInclusionKnown				= true;
			}
			
			if(strJunctions.contains(m_JunctionExclusion.m_nStart + "-" + m_JunctionExclusion.m_nEnd))
			{
				bExclusionKnown				= true;
			}
		}
		
		// if either junction is novel, the result is novel
		if(!bInclusionKnown || !bExclusionKnown)
			m_bIsNovel = true;
	}
	
	/** Writes the score to a binary output file */
	public void WriteToFile(FileOutputStream pOut) throws IOException
	{
		SplicingWebApp.WriteStringToFileOutputStream(m_strID, pOut);
		SplicingWebApp.WriteStringToFileOutputStream(m_strGeneID, pOut);
		SplicingWebApp.WriteStringToFileOutputStream(m_strConditionA, pOut);
		SplicingWebApp.WriteStringToFileOutputStream(m_strConditionB, pOut);
		
		// convert boolean to int
		int nIsNovel = (m_bIsNovel)? 1 : 0;
		int nIsSignificant = (m_bSignificant)? 1 : 0;
		
		ByteBuffer bb = ByteBuffer.allocate(Integer.BYTES * 5 + Double.BYTES * 5 + Character.BYTES * 2);
		bb.putInt(nIsNovel);
		bb.putInt(nIsSignificant);
		bb.putInt(m_nType);
		bb.putChar(m_JunctionInclusion.m_chStrand);
		bb.putInt(m_JunctionInclusion.m_nStart);
		bb.putInt(m_JunctionInclusion.m_nEnd);
		bb.putChar(m_JunctionExclusion.m_chStrand);
		bb.putInt(m_JunctionExclusion.m_nStart);
		bb.putInt(m_JunctionExclusion.m_nEnd);
		bb.putDouble(m_fPValue);
		bb.putDouble(m_fInclusionChange);
		pOut.write(bb.array());
	}
	
	/** Reads the score from a binary input file */
	public void ReadFromFile(FileInputStream pIn) throws IOException
	{
		m_strID 		= SplicingWebApp.ReadStringFromFileInputStream(pIn);
		m_strGeneID 	= SplicingWebApp.ReadStringFromFileInputStream(pIn);
		m_strConditionA	= SplicingWebApp.ReadStringFromFileInputStream(pIn);
		m_strConditionB	= SplicingWebApp.ReadStringFromFileInputStream(pIn);
		
		byte pBytes[] = new byte[Integer.BYTES * 5 + Double.BYTES * 5 + Character.BYTES * 2];
		if(pIn.read(pBytes) == -1) return;
		ByteBuffer bb = ByteBuffer.wrap(pBytes);
		
		m_JunctionInclusion				= new CountElement();
		m_JunctionExclusion				= new CountElement();
		
		int nIsNovel 					= bb.getInt();
		int nIsSignificant 				= bb.getInt();
		m_nType							= bb.getInt();
		m_JunctionInclusion.m_chStrand	= bb.getChar();
		m_JunctionInclusion.m_nStart	= bb.getInt();
		m_JunctionInclusion.m_nEnd		= bb.getInt();
		m_JunctionExclusion.m_chStrand	= bb.getChar();
		m_JunctionExclusion.m_nStart	= bb.getInt();
		m_JunctionExclusion.m_nEnd		= bb.getInt();
		m_fPValue						= bb.getDouble();
		m_fInclusionChange				= bb.getDouble();

		if(nIsNovel == 1)
			m_bIsNovel = true;
		if(nIsSignificant == 1)
			m_bSignificant = true; 
	}
	
	/** Returns whether a given exon is linked to the splicing event */
	public boolean AffectsExon(Exon ex)
	{
		if((m_JunctionInclusion.m_nEnd == ex.getCodingStart() || m_JunctionInclusion.m_nStart == ex.getCodingStop()) && m_JunctionExclusion.m_nStart < ex.getCodingStart() && m_JunctionExclusion.m_nEnd > ex.getCodingStop() )
			return true;
			
		return false;
	}
	
	/** Returns whether a given exon (defined by start and end position) is linked to the splicing event */
	public boolean AffectsExon(int nExStart, int nExEnd)
	{
		if((m_JunctionInclusion.m_nEnd == nExStart || m_JunctionInclusion.m_nStart == nExEnd) && m_JunctionExclusion.m_nStart < nExStart && m_JunctionExclusion.m_nEnd > nExEnd)
			return true;
			
		return false;
	}
}
