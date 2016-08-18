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
import java.util.Locale;
import java.util.TreeSet;

import BioKit.Exon;

/**
 *    This class holds an alternative splicing result, which includes various
 *    information on an alternatively spliced exon.
 *    
 *    1. The settings used to generate this result are stored
 *    2. The isoforms that were selected
 *    3. The type of the alternative splicing event, the alternatively
 *    spliced exon as AlternativeSplicingExon object and the PSI result
 *    as SimpleSpliceScore object.
 *    
 *    A second AlternativeSplicingExon object might be filled if the event
 *    refers to terminal exons (e.g. two start or end exons show a difference)
 *    
 *    The AlternativeSplicingExon objects refer to coverage ratio changes and
 *    the SimpleSpliceScore to a difference in the PSI score.
 *    
 *    Functions to read or write the data from/to a binary file are implemented.
*/
public class AnalysisResult implements Comparable<AnalysisResult>
{
	private int				m_nRating;
	private String			m_strRef;
	private String			m_strGeneID;
	private String			m_strGeneSymbol;
	private int				m_nMinJunctionReads;
	private int				m_nMinCovPerBase;
	private double			m_fMinCoveredBases;
	private double			m_fVariableExonThreshold;
	private TreeSet<String> m_vcSelectedIsoforms;
	private String 			m_strFileGTF;
	private String 			m_strComment;

	private int							m_nType;
	private int 						m_nID;
	private AlternativeSplicingExon 	m_altExonRatioChangeA;
	private AlternativeSplicingExon 	m_altExonRatioChangeB;	// only used for alternative terminal exons
	private SimpleSpliceScore 			m_spliceScore;
	
	static final int RESULT_EVIDENCE_TYPE_COMBINED			= 1;
	static final int RESULT_EVIDENCE_TYPE_RATIO_ONLY		= 2;
	static final int RESULT_EVIDENCE_TYPE_SPLIT_READ_ONLY	= 3;
	
	public AnalysisResult()
	{
		m_nRating 				= 0;
		m_nID 					= -1;
		m_altExonRatioChangeA	= null;
		m_altExonRatioChangeB	= null;
		m_spliceScore		 	= null;
		
		m_nType					= -1;
				
		m_strRef				= "";
		m_strGeneID				= "";
		m_strGeneSymbol			= "";		
		m_vcSelectedIsoforms	= null;
		
		m_strFileGTF			= "";
		m_strComment			= "";
		
		m_vcSelectedIsoforms	= new TreeSet<String>();
	}
	
	public AnalysisResult(AnalysisResult res)
	{
		m_nRating 				= res.m_nRating;
		m_nID 					= res.m_nID;
		m_altExonRatioChangeA	= res.m_altExonRatioChangeA;
		m_altExonRatioChangeB	= res.m_altExonRatioChangeB;
		m_spliceScore		 	= res.m_spliceScore;
		
		m_nType					= res.m_nType;
				
		m_strRef				= res.m_strRef;
		m_strGeneID				= res.m_strGeneID;
		m_strGeneSymbol			= res.m_strGeneSymbol;
		m_vcSelectedIsoforms	= res.m_vcSelectedIsoforms;
		
		m_strFileGTF			= res.m_strFileGTF;
		m_strComment			= res.m_strComment;
		
		m_nMinJunctionReads		= res.m_nMinJunctionReads;
		m_nMinCovPerBase		= res.m_nMinCovPerBase;
		m_fMinCoveredBases		= res.m_fMinCoveredBases;
		m_fVariableExonThreshold= res.m_fVariableExonThreshold;
		
		m_vcSelectedIsoforms	= new TreeSet<String>();
		m_vcSelectedIsoforms.addAll(res.m_vcSelectedIsoforms);
	}

	public AnalysisResult(SplicingWebApp app, int nID, String strRef, String strGeneID, String strGeneSymbol, AlternativeSplicingExon altExonRatioChangeA, AlternativeSplicingExon altExonRatioChangeB, SimpleSpliceScore spliceScore)
	{
		m_nID 					= nID;
		m_strRef				= strRef;
		m_strGeneID				= strGeneID;
		m_strGeneSymbol			= strGeneSymbol;
		m_altExonRatioChangeA	= altExonRatioChangeA;
		m_altExonRatioChangeB	= altExonRatioChangeB;
		m_spliceScore		 	= spliceScore;
		m_strFileGTF			= app.GetGTFFile();
		
		m_nMinJunctionReads		= app.GetMinimumJunctionReads();
		m_nMinCovPerBase		= app.GetMinimumCoveragePerBase();
		m_fMinCoveredBases		= app.GetMinimumCoveredBases();
		m_fVariableExonThreshold = app.GetVariableExonThreshold();
		
		m_vcSelectedIsoforms 	= new TreeSet<String>();
		for(String strIsoform : app.GetValidIsoforms())
			m_vcSelectedIsoforms.add(strIsoform);
		
		//###############################
		//         define type
		//###############################
		
		if(altExonRatioChangeA != null)
			m_nType = altExonRatioChangeA.GetType();
		else
			m_nType = spliceScore.GetType();
		
		m_strComment = "";
	}
	
	public String toString()
	{
		String strOut = "ID_" + m_nID;
		
		if(m_altExonRatioChangeA != null && m_spliceScore != null)
		{			
			strOut += "\t" + m_strGeneID + "\t" + m_strGeneSymbol;
			strOut += "\t" + GetResultTypeAsString();
			strOut += String.format(Locale.ENGLISH, "\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e", m_altExonRatioChangeA.GetTypeAsString(), m_altExonRatioChangeA.m_strExonGroup, m_altExonRatioChangeA.m_strConditionA, m_altExonRatioChangeA.m_strConditionB, m_altExonRatioChangeA.m_fAltExonCovPerBaseA, m_altExonRatioChangeA.m_fAltExonCovPerBaseB, Math.abs(m_altExonRatioChangeA.m_fAltExonCovPerBaseA-m_altExonRatioChangeA.m_fAltExonCovPerBaseB), m_altExonRatioChangeA.m_fFractionChangeAbsolute, m_altExonRatioChangeA.m_fFractionChangeRelative, m_altExonRatioChangeA.m_fPValue);
			strOut += String.format(Locale.ENGLISH, "\t%d\t%d\t%d\t%d\t%.2f\t%.3e\t%s\t%b", m_spliceScore.m_JunctionInclusion.m_nStart+1, m_spliceScore.m_JunctionInclusion.m_nEnd-1, m_spliceScore.m_JunctionExclusion.m_nStart+1, m_spliceScore.m_JunctionExclusion.m_nEnd-1, m_spliceScore.m_fInclusionChange, m_spliceScore.m_fPValue, m_spliceScore.m_bSignificant, m_spliceScore.m_bIsNovel);
		}
		else if(m_altExonRatioChangeA != null)
		{
			strOut += "\t" + m_strGeneID + "\t" + m_strGeneSymbol;
			strOut += "\t" + GetResultTypeAsString();			
			strOut += String.format(Locale.ENGLISH, "\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e", m_altExonRatioChangeA.GetTypeAsString(), m_altExonRatioChangeA.m_strExonGroup, m_altExonRatioChangeA.m_strConditionA, m_altExonRatioChangeA.m_strConditionB, m_altExonRatioChangeA.m_fAltExonCovPerBaseA, m_altExonRatioChangeA.m_fAltExonCovPerBaseB, Math.abs(m_altExonRatioChangeA.m_fAltExonCovPerBaseA-m_altExonRatioChangeA.m_fAltExonCovPerBaseB), m_altExonRatioChangeA.m_fFractionChangeAbsolute, m_altExonRatioChangeA.m_fFractionChangeRelative, m_altExonRatioChangeA.m_fPValue);
			strOut += String.format(Locale.ENGLISH, "\t-\t-\t-\t-\t-\t-\t-\t-");
		}
		else
		{
			strOut += "\t" + m_strGeneID + "\t" + m_strGeneSymbol;
			strOut += "\t" + GetResultTypeAsString();
			strOut += String.format(Locale.ENGLISH, "\t%s\t-\t-\t-\t-\t-\t-\t-\t-\t-", m_spliceScore.GetTypeAsString());
			strOut += String.format(Locale.ENGLISH, "\t%d\t%d\t%d\t%d\t%.2f\t%.3e\t%s\t%b", m_spliceScore.m_JunctionInclusion.m_nStart+1, m_spliceScore.m_JunctionInclusion.m_nEnd-1, m_spliceScore.m_JunctionExclusion.m_nStart+1, m_spliceScore.m_JunctionExclusion.m_nEnd-1, m_spliceScore.m_fInclusionChange, m_spliceScore.m_fPValue, m_spliceScore.m_bSignificant, m_spliceScore.m_bIsNovel); 
		}
		return strOut;
	}
	
	public String toID()
	{
		String strOut = "ID_" + m_nID;
		
		if(m_altExonRatioChangeA != null && m_spliceScore != null)
		{
			strOut += "\tcombined";
			strOut += "\t" + m_nRating + "\t" + m_strGeneID + "\t" + m_strGeneSymbol + "\t" + m_strFileGTF + "\t" + m_nMinJunctionReads + "\t" + m_nMinCovPerBase + "\t" + m_fMinCoveredBases + "\t" + m_fVariableExonThreshold + "\t" + m_strComment;
			strOut += String.format(Locale.ENGLISH, "\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e", m_altExonRatioChangeA.GetTypeAsString(), m_altExonRatioChangeA.m_strExonGroup, m_altExonRatioChangeA.m_strConditionA, m_altExonRatioChangeA.m_strConditionB, m_altExonRatioChangeA.m_fAltExonCovPerBaseA, m_altExonRatioChangeA.m_fAltExonCovPerBaseB, Math.abs(m_altExonRatioChangeA.m_fAltExonCovPerBaseA-m_altExonRatioChangeA.m_fAltExonCovPerBaseB), m_altExonRatioChangeA.m_fFractionChangeAbsolute, m_altExonRatioChangeA.m_fFractionChangeRelative, m_altExonRatioChangeA.m_fPValue);
			strOut += String.format(Locale.ENGLISH, "\t%d\t%d\t%d\t%d\t%.2f\t%.3e\t%s\t%b", m_spliceScore.m_JunctionInclusion.m_nStart+1, m_spliceScore.m_JunctionInclusion.m_nEnd-1, m_spliceScore.m_JunctionExclusion.m_nStart+1, m_spliceScore.m_JunctionExclusion.m_nEnd-1, m_spliceScore.m_fInclusionChange, m_spliceScore.m_fPValue, m_spliceScore.m_bSignificant, m_spliceScore.m_bIsNovel);
		}
		else if(m_altExonRatioChangeA != null)
		{
			strOut += "\t" + m_nRating + "\t" + m_strGeneID + "\t" + m_strGeneSymbol + "\t" + m_strFileGTF + "\t" + m_nMinJunctionReads + "\t" + m_nMinCovPerBase + "\t" + m_fMinCoveredBases + "\t" + m_fVariableExonThreshold + "\t" + m_strComment + "\n";
			strOut += "\tratio_only";			
			strOut += String.format(Locale.ENGLISH, "\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e", m_altExonRatioChangeA.GetTypeAsString(), m_altExonRatioChangeA.m_strExonGroup, m_altExonRatioChangeA.m_strConditionA, m_altExonRatioChangeA.m_strConditionB, m_altExonRatioChangeA.m_fAltExonCovPerBaseA, m_altExonRatioChangeA.m_fAltExonCovPerBaseB, Math.abs(m_altExonRatioChangeA.m_fAltExonCovPerBaseA-m_altExonRatioChangeA.m_fAltExonCovPerBaseB), m_altExonRatioChangeA.m_fFractionChangeAbsolute, m_altExonRatioChangeA.m_fFractionChangeRelative, m_altExonRatioChangeA.m_fPValue);
			strOut += String.format(Locale.ENGLISH, "\t-\t-\t-\t-\t-\t-\t-\t-");
		}
		else
		{
			strOut += "\t" + m_nRating + "\t" + m_strGeneID + "\t" + m_strGeneSymbol + "\t" + m_strFileGTF + "\t" + m_nMinJunctionReads + "\t" + m_nMinCovPerBase + "\t" + m_fMinCoveredBases + "\t" + m_fVariableExonThreshold + "\t" + m_strComment + "\n";
			strOut += "\tsplit_read_only";
			strOut += String.format(Locale.ENGLISH, "\t%s\t-\t-\t-\t-\t-\t-\t-\t-\t-", m_spliceScore.GetTypeAsString());
			strOut += String.format(Locale.ENGLISH, "\t%d\t%d\t%d\t%d\t%.2f\t%.3e\t%s\t%b", m_spliceScore.m_JunctionInclusion.m_nStart+1, m_spliceScore.m_JunctionInclusion.m_nEnd-1, m_spliceScore.m_JunctionExclusion.m_nStart+1, m_spliceScore.m_JunctionExclusion.m_nEnd-1, m_spliceScore.m_fInclusionChange, m_spliceScore.m_fPValue, m_spliceScore.m_bSignificant, m_spliceScore.m_bIsNovel); 
		}
		return strOut;
	}
	
	public int GetType()
	{
		return m_nType;
	}

	/**
	 *    This function tests whether the alternative splicing event overlaps as region,
	 *    specified by reference name, start and end position.
	 *    Additionally, the type of the alternative splicing event must be specified
	*/
	public boolean OverlapsRegion(String strRef, int nStart, int nEnd, int nType)
	{
		if(m_nType == nType)
		{
			if(m_altExonRatioChangeA != null)
			{
				if(m_strRef.equals(strRef) && nStart <= m_altExonRatioChangeA.GetStart() && nEnd >= m_altExonRatioChangeA.GetEnd())
					return true;
			}
			
			if(m_altExonRatioChangeB != null)
			{
				if(m_strRef.equals(strRef) && nStart <= m_altExonRatioChangeB.GetStart() && nEnd >= m_altExonRatioChangeB.GetEnd())
					return true;
			}
			
			if(m_spliceScore != null)
			{
				if(m_strRef.equals(strRef) && nStart == m_spliceScore.m_JunctionInclusion.m_nEnd && nEnd == m_spliceScore.m_JunctionInclusion.m_nStart)
					return true;
			}
		}
		
		return false;
	}
	
	/**
	 *    This function checks whether any of the result components (coverage ratio change in exon A,
	 *    coverage ratio change in exon B or change in PSI score) refers to a given exon.
	 *    It requires three parameters:
	 *    - The reference name of the query exon (because the BioKit Exon class does not include it)
	 *    - The query exon
	 *    - the type of alternative splicing event
	 */
	public boolean IncludesASExon(String strRef, Exon ex, int nType)
	{
		int nStart = ex.getCodingStart();
		int nEnd   = ex.getCodingStop();
		
		if(m_nType == nType)
		{			
			if(m_altExonRatioChangeA != null)
			{
				if(m_strRef.equals(strRef) && nStart <= m_altExonRatioChangeA.GetStart() && nEnd >= m_altExonRatioChangeA.GetEnd())
					return true;
			}
			
			if(m_altExonRatioChangeB != null)
			{
				if(m_strRef.equals(strRef) && nStart <= m_altExonRatioChangeB.GetStart() && nEnd >= m_altExonRatioChangeB.GetEnd())
					return true;
			}
			
			if(m_spliceScore != null)
			{
				if(m_strRef.equals(strRef) && (nStart == m_spliceScore.m_JunctionInclusion.m_nEnd || nEnd == m_spliceScore.m_JunctionInclusion.m_nStart) && nStart > m_spliceScore.m_JunctionExclusion.m_nStart && nEnd < m_spliceScore.m_JunctionExclusion.m_nEnd)
					return true;
			}
		}
		
		return false;
	}
	
	public boolean UsesNovelJunction()
	{
		if(m_spliceScore != null)
			return m_spliceScore.m_bIsNovel;
		
		return false;
	}
	
	public boolean HasSignificantPSIScore()
	{
		if(m_spliceScore != null)
			return m_spliceScore.m_bSignificant;
		
		return false;
	}
	
	public boolean HasAltExonA()
	{
		if(m_altExonRatioChangeA != null)
			return true;
		
		return false;
	}
	
	public boolean HasAltExonB()
	{
		if(m_altExonRatioChangeB != null)
			return true;
		
		return false;
	}
	
	public boolean HasPSIScore()
	{
		if(m_spliceScore != null)
			return true;
		
		return false;
	}
	
	public AlternativeSplicingExon GetAltExonA()
	{
		return m_altExonRatioChangeA;
	}
	
	public AlternativeSplicingExon GetAltExonB()
	{
		return m_altExonRatioChangeB;
	}

	public SimpleSpliceScore GetPSIScore()
	{
		return m_spliceScore;
	}
	
	public double GetAbsoluteChangeA()
	{
		if(m_altExonRatioChangeA != null)
			return m_altExonRatioChangeA.m_fFractionChangeAbsolute;
		else
			return 0;
	}
		
	public double GetAbsoluteChangeB()
	{
		if(m_altExonRatioChangeB != null)
			return m_altExonRatioChangeB.m_fFractionChangeAbsolute;
		else
			return 0;
	}
	
	public double GetRelativeChangeA()
	{
		if(m_altExonRatioChangeA != null)
			return m_altExonRatioChangeA.m_fFractionChangeRelative;
		else
			return 0.0;
	}
	
	public double GetRelativeChangeB()
	{
		if(m_altExonRatioChangeB != null)
			return m_altExonRatioChangeB.m_fFractionChangeRelative;
		else
			return 0.0;
	}
	
	public double GetInclusionChange()
	{
		if(m_spliceScore != null)
			return m_spliceScore.m_fInclusionChange;
		else
			return 0.0;
	}
	
	public String GetConditionA()
	{
		if(m_altExonRatioChangeA != null)
		{
			return m_altExonRatioChangeA.m_strConditionA;
		}
		
		if(m_spliceScore != null)
			return m_spliceScore.m_strConditionA;
		
		return "";
	}
	
	public String GetConditionB()
	{
		if(m_altExonRatioChangeA != null)
		{
			return m_altExonRatioChangeA.m_strConditionB;
		}
		
		if(m_spliceScore != null)
			return m_spliceScore.m_strConditionB;
		
		return "";
	}
	
	public double GetPValueA()
	{
		if(m_altExonRatioChangeA != null)
			return m_altExonRatioChangeA.GetPValue();
		else
			return 99.9;
	}
	
	public double GetPValueB()
	{
		if(m_altExonRatioChangeB != null)
			return m_altExonRatioChangeB.GetPValue();
		else
			return 99.9;
	}
	
	public double GetPValuePSI()
	{
		if(m_spliceScore != null)
			return m_spliceScore.m_fPValue;
		else
			return 99.9;
	}
	
	public int GetStartA()
	{
		return m_altExonRatioChangeA.GetStart();
	}
	
	public int GetStartB()
	{
		return m_altExonRatioChangeB.GetStart();
	}
	
	public int GetEndA()
	{
		return m_altExonRatioChangeA.GetEnd();
	}
	
	public int GetEndB()
	{
		return m_altExonRatioChangeB.GetEnd();
	}
	
	public int GetInclusionJunctionStart()
	{
		if(m_spliceScore != null)
			return m_spliceScore.m_JunctionInclusion.m_nStart;
		
		return -1;
	}
	
	public int GetInclusionJunctionEnd()
	{
		if(m_spliceScore != null)
			return m_spliceScore.m_JunctionInclusion.m_nEnd;
		
		return -1;
	}
	
	public int GetExclusionJunctionStart()
	{
		if(m_spliceScore != null)
			return m_spliceScore.m_JunctionExclusion.m_nStart;
		
		return -1;
	}
	
	public int GetExclusionJunctionEnd()
	{
		if(m_spliceScore != null)
			return m_spliceScore.m_JunctionExclusion.m_nEnd;
		
		return -1;
	}
	
	public CountElement GetInclusionJunction()
	{
		if(m_spliceScore != null)
			return m_spliceScore.m_JunctionInclusion;
		
		return null;
	}
	
	public CountElement GetExclusionJunction()
	{
		if(m_spliceScore != null)
			return m_spliceScore.m_JunctionExclusion;
		
		return null;
	}
	
	public String GetReferenceName()
	{
		return m_strRef;
	}
	
	public String GetGeneID()
	{
		return m_strGeneID;
	}
	
	public String GetGeneSymbol()
	{
		return m_strGeneSymbol;
	}
	
	public String GetComment()
	{
		return m_strComment;
	}
	
	public int GetRating()
	{
		return m_nRating;
	}
	
	public int GetMinCovPerBase()
	{
		return m_nMinCovPerBase;
	}
	
	public double GetMinCoveredBases()
	{
		return m_fMinCoveredBases;
	}
	
	public int GetMinJunctionReads()
	{
		return m_nMinJunctionReads;
	}
	
	public double GetVariableExonThreshold()
	{
		return m_fVariableExonThreshold;
	}
	
	public String GetGTFFile()
	{
		return m_strFileGTF;
	}
	
	public TreeSet<String> GetSelectedIsoforms()
	{
		return m_vcSelectedIsoforms;
	}
	
	public TreeSet<String> GetValidIsoforms(DataSupplier data)
	{
		TreeSet<String> vcResult = new TreeSet<String>();
		
		for(String strIsoform : data.GetIsoformNames())
		{			
			// ignore unselected Isoforms
			if(!m_vcSelectedIsoforms.contains(strIsoform))
				continue;
			
			// select isoforms that contain the correct junctions (if available)
			if(m_spliceScore != null)
			{
				boolean bIncludesInclusionJunction = false;
				boolean bIncludesExclusionJunction = false;
				
				TreeSet<CountElement> vcJunctions = data.GetJunctionsForIsoform(strIsoform);

				for(CountElement jun : vcJunctions)
				{
					if(jun.m_nStart == m_spliceScore.m_JunctionInclusion.m_nStart && jun.m_nEnd == m_spliceScore.m_JunctionInclusion.m_nEnd)
						bIncludesInclusionJunction = true;
					
					if(jun.m_nStart == m_spliceScore.m_JunctionExclusion.m_nStart && jun.m_nEnd == m_spliceScore.m_JunctionExclusion.m_nEnd)
						bIncludesExclusionJunction = true;
				}
				
				if(bIncludesInclusionJunction || bIncludesExclusionJunction)
				{
					vcResult.add(strIsoform);
				}
			}
			else
				vcResult.add(strIsoform);
		}
		
		return vcResult;
	}
	
	public String GetTypeAsString()
	{
		String strType = "?";
		switch(m_nType)
		{
			case SplicingWebApp.AS_TYPE_EXON_SKIPPING: 					strType = "exn_skipping"; 					break;
			case SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN: 			strType = "alt_start_unique_jun"; 			break;
			case SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN: 			strType = "alt_end_unique_jun"; 			break;
			case SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN: 			strType = "alt_start_shared_jun"; 			break;
			case SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN: 			strType = "alt_end_shared_jun"; 			break;
			case SplicingWebApp.AS_TYPE_RETAINED_INTRON:				strType = "retained_intron";				break;
			case SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN_DOUBLE: 	strType = "alt_start_unique_jun_double"; 	break;
			case SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN_DOUBLE: 		strType = "alt_end_unique_jun_double"; 		break;
			case SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN_DOUBLE: 	strType = "alt_start_shared_jun_double"; 	break;
			case SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN_DOUBLE: 		strType = "alt_end_shared_jun_double"; 		break;
			case SplicingWebApp.AS_TYPE_ALT_3_PRIME_EXON_END:			strType = "alt_3_prime_exon_end";			break;
			case SplicingWebApp.AS_TYPE_ALT_5_PRIME_EXON_END:			strType = "alt_5_prime_exon_end";			break;
			default:													strType = "" + m_nType; 					break;
		}
		
		return strType;
	}
	
	public String GetResultTypeAsString()
	{
		String strOut = "?";
		
		if(m_altExonRatioChangeA != null && m_spliceScore != null)
		{
			strOut = "combined";
		}
		else if(m_altExonRatioChangeA != null)
		{
			strOut = "ratio_only";
		}
		else
		{
			strOut = "split_read_only";
		}
		
		return strOut;
	}
	
	public int GetResultType()
	{
		if(m_altExonRatioChangeA != null && m_spliceScore != null)
			return RESULT_EVIDENCE_TYPE_COMBINED;

		else if(m_altExonRatioChangeA != null)
			return RESULT_EVIDENCE_TYPE_RATIO_ONLY;

		else if(m_spliceScore != null)
			return RESULT_EVIDENCE_TYPE_SPLIT_READ_ONLY;
		
		return 0;
	}
	
	public void SetID(int nID)
	{
		m_nID = nID;
	}
	
	public void SetRating(int nRating)
	{
		m_nRating = nRating;
	}
	
	public void SetComment(String strComment)
	{
		m_strComment = strComment;
	}

	public int GetID()
	{
		return m_nID;
	}

	@Override
	public int compareTo(AnalysisResult other)
	{
		int nComp = m_strRef.compareTo(other.m_strRef);
		if(nComp != 0) return nComp;
		
		nComp = m_strGeneID.compareTo(other.m_strGeneID);
		if(nComp != 0) return nComp;		
		
		if(m_vcSelectedIsoforms.size() < other.m_vcSelectedIsoforms.size())
			return -1;
		else if(m_vcSelectedIsoforms.size() > other.m_vcSelectedIsoforms.size())
			return 1;
		
		if(m_nType < other.m_nType)
			return -1;
		else if(m_nType > other.m_nType)
			return 1;
		
		if(!HasAltExonA() && other.HasAltExonA())
			return -1;
		
		if(HasAltExonA() && !other.HasAltExonA())
			return 1;
		
		if(HasAltExonA() && other.HasAltExonA())
		{
			nComp = m_altExonRatioChangeA.compareTo(other.m_altExonRatioChangeA);
			if(nComp != 0) return nComp;
		}
		
		if(!HasAltExonB() && other.HasAltExonB())
			return -1;
		
		if(HasAltExonB() && !other.HasAltExonB())
			return 1;
		
		if(HasAltExonB() && other.HasAltExonB())
		{
			nComp = m_altExonRatioChangeB.compareTo(other.m_altExonRatioChangeB);
			if(nComp != 0) return nComp;
		}
		
		if(!HasPSIScore() && other.HasPSIScore())
			return -1;
		
		if(HasPSIScore() && !other.HasPSIScore())
			return 1;
		
		if(HasPSIScore() && other.HasPSIScore())
		{
			nComp = m_spliceScore.compareTo(other.m_spliceScore);
			if(nComp != 0) return nComp;
		}

		return 0;
	}

	public void WriteToFile(FileOutputStream pOut) throws IOException
	{
		SplicingWebApp.WriteStringToFileOutputStream(m_strRef, pOut);
		SplicingWebApp.WriteStringToFileOutputStream(m_strGeneID, pOut);
		SplicingWebApp.WriteStringToFileOutputStream(m_strGeneSymbol, pOut);
		
		int nHasAltExonA = HasAltExonA() ? 1 : 0;
		int nHasAltExonB = HasAltExonB() ? 1 : 0;
		int nHasPSIScore = HasPSIScore() ? 1 : 0;
		
		ByteBuffer bb = ByteBuffer.allocate(Integer.BYTES * 9 + Double.BYTES * 2);
		bb.putInt(m_nRating);
		bb.putInt(m_nMinJunctionReads);
		bb.putInt(m_nMinCovPerBase);
		bb.putInt(m_nType);
		bb.putInt(m_nID);
		bb.putDouble(m_fMinCoveredBases);
		bb.putDouble(m_fVariableExonThreshold);
		bb.putInt(m_vcSelectedIsoforms.size());
		bb.putInt(nHasAltExonA);
		bb.putInt(nHasAltExonB);
		bb.putInt(nHasPSIScore);
		pOut.write(bb.array());

		for(String strIsoform : m_vcSelectedIsoforms)
		{
			SplicingWebApp.WriteStringToFileOutputStream(strIsoform, pOut);
		}
		
		SplicingWebApp.WriteStringToFileOutputStream(m_strFileGTF, pOut);
		SplicingWebApp.WriteStringToFileOutputStream(m_strComment, pOut);

		if(HasAltExonA())
			m_altExonRatioChangeA.WriteToFile(pOut);
		if(HasAltExonB())
			m_altExonRatioChangeB.WriteToFile(pOut);
		if(HasPSIScore())
			m_spliceScore.WriteToFile(pOut);
	}
	
	public void ReadFromFile(FileInputStream pIn) throws IOException
	{
		m_strRef		= SplicingWebApp.ReadStringFromFileInputStream(pIn);
		m_strGeneID		= SplicingWebApp.ReadStringFromFileInputStream(pIn);
		m_strGeneSymbol	= SplicingWebApp.ReadStringFromFileInputStream(pIn);

		byte pBytes[] = new byte[Integer.BYTES * 9 + Double.BYTES * 2];
		if(pIn.read(pBytes) == -1) return;
		ByteBuffer bb = ByteBuffer.wrap(pBytes);
		m_nRating 			= bb.getInt();
		m_nMinJunctionReads = bb.getInt();
		m_nMinCovPerBase 	= bb.getInt();
		m_nType 			= bb.getInt();
		m_nID 				= bb.getInt();
		m_fMinCoveredBases			= bb.getDouble();
		m_fVariableExonThreshold	= bb.getDouble();
		int nIsoforms 		= bb.getInt();
		int nHasAltExonA 	= bb.getInt();
		int nHasAltExonB 	= bb.getInt();
		int nHasPSIScore	= bb.getInt();		

		m_vcSelectedIsoforms = new TreeSet<String>();
		for(int i=0; i<nIsoforms; i++)
		{
			String strIsoform = SplicingWebApp.ReadStringFromFileInputStream(pIn);	
			m_vcSelectedIsoforms.add(strIsoform);
		}

		m_strFileGTF = SplicingWebApp.ReadStringFromFileInputStream(pIn);
		m_strComment = SplicingWebApp.ReadStringFromFileInputStream(pIn);
		
		if(nHasAltExonA == 1)
		{	
			m_altExonRatioChangeA = new AlternativeSplicingExon();
			m_altExonRatioChangeA.ReadFromFile(pIn);
		}
		else
		{
			m_altExonRatioChangeA = null;
		}
		
		if(nHasAltExonB == 1)
		{
			m_altExonRatioChangeB = new AlternativeSplicingExon();
			m_altExonRatioChangeB.ReadFromFile(pIn);
		}
		else
		{
			m_altExonRatioChangeB = null;
		}
		
		if(nHasPSIScore == 1)
		{
			m_spliceScore = new SimpleSpliceScore();
			m_spliceScore.ReadFromFile(pIn);
		}
		else
		{
			m_spliceScore = null;
		}
	}
}
