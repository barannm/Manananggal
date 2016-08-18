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

import java.util.Locale;
import java.util.TreeSet;

/**
 *    Most functionality was replaced by "AnalysisResult".
 *    
 *    This class is used during the detection of exon skipping event and holds
 *    a TreeSet of potentially alternatively spliced exons and information
 *    on the settings used for the detection.
 */
public class AlternativeSplicingHit implements Comparable<AlternativeSplicingHit>
{
	private String			m_strGeneID;
	private String			m_strGeneSymbol;
	private int				m_nMinJunctionReads;
	private int				m_nMinCovPerBase;
	private double			m_fMinCoveredBases;
	private double			m_fVariableExonThreshold;
	private TreeSet<String> m_vcSelectedIsoforms;
	private String 			m_strFileGTF;
	private String 			m_strComment;
	
	TreeSet<AlternativeSplicingExon> m_vcAlternativeSplicedExons;

	public AlternativeSplicingHit()
	{
		m_strGeneID					= "?";
		m_strGeneSymbol				= "?";
		m_nMinJunctionReads			= 0;
		m_nMinCovPerBase			= 0;
		m_fMinCoveredBases			= 0.0f;
		m_fVariableExonThreshold 	= 0.0f;
		m_vcSelectedIsoforms	 	= new TreeSet<String>();
		m_strFileGTF				= "?";
		m_strComment				= "?";

		
		m_vcAlternativeSplicedExons = new TreeSet<AlternativeSplicingExon>();
	}
	
	public AlternativeSplicingHit(GeneIdentifier gid, String strRef, int nMinJunctionReads, int nMinCovPerBase,
										double fMinCoveredBases, double fVariableExonThreshold,
										TreeSet<String> vcSelectedIsoforms, String strFileGTF, String strComment)
	{
		m_strGeneID			= gid.m_strEnsemblGeneID;
		m_strGeneSymbol		= gid.m_strApprovedGeneSymbol;
		m_nMinJunctionReads	= nMinJunctionReads;
		m_nMinCovPerBase	= nMinCovPerBase;
		m_fMinCoveredBases	= fMinCoveredBases;
		m_fVariableExonThreshold = fVariableExonThreshold;
		m_vcSelectedIsoforms 	 = vcSelectedIsoforms;
		m_strFileGTF		= strFileGTF;
		m_strComment		= strComment;
		
		m_vcAlternativeSplicedExons = new TreeSet<AlternativeSplicingExon>();
	}
	
	// adds an alternatively spliced exon to this result
	public void AddASExons(int nType, String strExonGroup, int nStart, int nEnd, String strConditionA, String strConditionB, double fAltExonCovPerBaseA, double fAltExonCovPerBaseB, double fFractionChangeAbsolute, double fFractionChangeRelative, double fPValue, double[] pFractionTestedExon, double[] pFractionOtherExons)
	{		
		AlternativeSplicingExon hit = new AlternativeSplicingExon(nType, strExonGroup, nStart, nEnd, strConditionA, strConditionB, fAltExonCovPerBaseA, fAltExonCovPerBaseB, fFractionChangeAbsolute, fFractionChangeRelative, fPValue, pFractionTestedExon, pFractionOtherExons);
		m_vcAlternativeSplicedExons.add(hit);
	}
	
	// adds an alternatively spliced exon to this result
	public void AddASExon(AlternativeSplicingExon hit)
	{
		m_vcAlternativeSplicedExons.add(hit);
	}

	@Override
	public int compareTo(AlternativeSplicingHit other)
	{
		return m_strGeneID.compareTo(other.m_strGeneID);
	}

	public TreeSet<String> GetIsoforms()
	{
		return m_vcSelectedIsoforms;
	}

	public int GetNumHits()
	{
		return m_vcAlternativeSplicedExons.size();
	}
	
	public AlternativeSplicingExon GetResult(int nIdx)
	{
		int i=0;
		for(AlternativeSplicingExon hit : m_vcAlternativeSplicedExons)
		{
			if(i == nIdx)
				return hit;
			i++;
		}
		
		return null;
	}

	@Override
	public String toString()
	{	
		String strRes = m_strGeneID + "\t" + m_strGeneSymbol + "\t" + m_strFileGTF + "\t" + m_nMinJunctionReads + "\t" + m_nMinCovPerBase + "\t" + m_fMinCoveredBases + "\t" + m_fVariableExonThreshold + "\t" + m_strComment + "\n";
		strRes += "Isoforms: ";
		for(String strIsoform : m_vcSelectedIsoforms)
			strRes += "\t" + strIsoform;
		strRes += "\n";
		
		for(AlternativeSplicingExon hit : m_vcAlternativeSplicedExons)
		{
			String strType = "?";
			switch(hit.m_nType)
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
				default: strType = "" + hit.m_nType; 			break;
			}
			
			strRes += String.format(Locale.ENGLISH, "\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e\n", strType, hit.m_strExonGroup, hit.m_strConditionA, hit.m_strConditionB, hit.m_fAltExonCovPerBaseA, hit.m_fAltExonCovPerBaseB, hit.m_fFractionChangeAbsolute, hit.m_fFractionChangeRelative, hit.m_fPValue);
		}
		
		return strRes;
	}

}
