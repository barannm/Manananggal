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

import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.Locale;
import java.util.TreeSet;

import BioKit.Exon;
import BioKit.ExonGroup;

public class AlternativeSplicingHit implements Comparable<AlternativeSplicingHit>
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
	
	TreeSet<AlternativeSplicingExon> m_vcAlternativeSplicedExons;

	public AlternativeSplicingHit()
	{
		m_nRating					= 0;
		m_strRef					= "?";
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
	
	public AlternativeSplicingHit(int nRating, GeneIdentifier gid, String strRef, int nMinJunctionReads, int nMinCovPerBase,
										double fMinCoveredBases, double fVariableExonThreshold,
										TreeSet<String> vcSelectedIsoforms, String strFileGTF, String strComment)
	{
		m_nRating			= nRating;
		m_strRef			= strRef;
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
	
	public void AddASExons(int nType, String strExonGroup, int nStart, int nEnd, String strConditionA, String strConditionB, double fAltExonCovPerBaseA, double fAltExonCovPerBaseB, double fFractionChangeAbsolute, double fFractionChangeRelative, double fPValue, double[] pFractionTestedExon, double[] pFractionOtherExons)
	{		
		AlternativeSplicingExon hit = new AlternativeSplicingExon(nType, strExonGroup, nStart, nEnd, strConditionA, strConditionB, fAltExonCovPerBaseA, fAltExonCovPerBaseB, fFractionChangeAbsolute, fFractionChangeRelative, fPValue, pFractionTestedExon, pFractionOtherExons);
		m_vcAlternativeSplicedExons.add(hit);
	}
	
	public void AddASExon(AlternativeSplicingExon hit)
	{
		m_vcAlternativeSplicedExons.add(hit);
	}
	
	public void ParseString(String strLine)
	{
		String pSplit[] = strLine.split("\t");

		m_nRating			= Integer.parseInt(pSplit[0]);
		m_strGeneID			= pSplit[1];
		m_strGeneSymbol		= pSplit[2];
		m_strFileGTF		= pSplit[3];
		m_nMinJunctionReads = Integer.parseInt(pSplit[4]);
		m_nMinCovPerBase	= Integer.parseInt(pSplit[5]);
		m_fMinCoveredBases	= Double.parseDouble(pSplit[6]);
		m_fVariableExonThreshold = Double.parseDouble(pSplit[7]);
		m_strComment		= pSplit[8];
		
		if(m_vcSelectedIsoforms == null)
			m_vcSelectedIsoforms = new TreeSet<String>();
		
		for(int i=9; i<pSplit.length; i++)
		{
			m_vcSelectedIsoforms.add(pSplit[i]);
		}
	}

	@Override
	public int compareTo(AlternativeSplicingHit other)
	{
		return m_strGeneID.compareTo(other.m_strGeneID);
	}
	
	public void WriteToFile(RandomAccessFile pFile) throws IOException
	{
		if(pFile == null)
			return;

		pFile.writeUTF(m_strRef);
		pFile.writeUTF(m_strGeneID);
		pFile.writeUTF(m_strGeneSymbol);
		
		pFile.writeInt(m_nRating);
		
		pFile.writeInt(m_nMinJunctionReads);
		pFile.writeInt(m_nMinCovPerBase);
		pFile.writeDouble(m_fMinCoveredBases);
		pFile.writeDouble(m_fVariableExonThreshold);
		
		pFile.writeInt(m_vcSelectedIsoforms.size());
		for(String strIsoform : m_vcSelectedIsoforms)
			pFile.writeUTF(strIsoform);

		pFile.writeUTF(m_strFileGTF);
		pFile.writeUTF(m_strComment);
		
		pFile.writeInt(m_vcAlternativeSplicedExons.size());
		
		for(AlternativeSplicingExon hit : m_vcAlternativeSplicedExons)
		{
			pFile.writeInt(hit.m_nType);
			
			pFile.writeUTF(hit.m_strExonGroup);
			
			pFile.writeInt(hit.m_nStart);
			pFile.writeInt(hit.m_nEnd);
			
			pFile.writeUTF(hit.m_strConditionA);
			pFile.writeUTF(hit.m_strConditionB);
			
			pFile.writeDouble(hit.m_fAltExonCovPerBaseA);
			pFile.writeDouble(hit.m_fAltExonCovPerBaseB);
			
			pFile.writeDouble(hit.m_fFractionChangeAbsolute);
			pFile.writeDouble(hit.m_fFractionChangeRelative);
			pFile.writeDouble(hit.m_fPValue);
			
			pFile.writeInt(hit.m_pFractionTestedExon.length);
			for(double fVal : hit.m_pFractionTestedExon)
				pFile.writeDouble(fVal);
			
			pFile.writeInt(hit.m_pFractionOtherExons.length);
			for(double fVal : hit.m_pFractionOtherExons)
				pFile.writeDouble(fVal);
		}
	}
	
	public boolean ReadFromFile(RandomAccessFile pFile)
	{
		if(pFile == null)
			return false;

		try
		{
			m_strRef		= pFile.readUTF();
			m_strGeneID		= pFile.readUTF();
			m_strGeneSymbol = pFile.readUTF();
			
			m_nRating		= pFile.readInt();
			
			m_nMinJunctionReads = pFile.readInt();
			m_nMinCovPerBase	= pFile.readInt();
			m_fMinCoveredBases	= pFile.readDouble();
			m_fVariableExonThreshold = pFile.readDouble();
			
			m_vcSelectedIsoforms = new TreeSet<String>();
			int nCount = pFile.readInt();
			for(int i=0; i<nCount; i++)
				m_vcSelectedIsoforms.add(pFile.readUTF());
	
			m_strFileGTF = pFile.readUTF();
			m_strComment = pFile.readUTF();
			
			int nHits = pFile.readInt();
			
			for(int i=0; i<nHits; i++)
			{
				AlternativeSplicingExon hit = new AlternativeSplicingExon();
				
				hit.m_nType			= pFile.readInt();
				
				hit.m_strExonGroup	= pFile.readUTF();
				hit.m_nStart		= pFile.readInt();
				hit.m_nEnd			= pFile.readInt();
				
				hit.m_strConditionA	= pFile.readUTF();
				hit.m_strConditionB	= pFile.readUTF();
		
				hit.m_fFractionChangeAbsolute = pFile.readDouble();
				hit.m_fFractionChangeRelative = pFile.readDouble();
				hit.m_fPValue 				  = pFile.readDouble();
				
				hit.m_fAltExonCovPerBaseA	  = pFile.readDouble();
				hit.m_fAltExonCovPerBaseB	  = pFile.readDouble();
				
				int nFractions = pFile.readInt();
				hit.m_pFractionTestedExon = new double[nFractions];
				for(int j=0; j<nFractions; j++)
					hit.m_pFractionTestedExon[j] = pFile.readDouble();
				
				nFractions = pFile.readInt();
				hit.m_pFractionOtherExons = new double[nFractions];
				for(int j=0; j<nFractions; j++)
					hit.m_pFractionOtherExons[j] = pFile.readDouble();
				
				hit.m_strID = hit.m_nType + "_" + hit.m_strExonGroup + "_" + hit.m_strConditionA + "_" + hit.m_strConditionB;
			}
		}
		catch(IOException ex)
		{
			System.out.println(this);
			
			System.out.println(ex.toString());
			System.out.println(ex.getStackTrace());
			System.out.println(ex.getMessage());
			return false;
		}
		
		return true;
	}

	public boolean ContainsExonGroup(String strQuery)
	{
		for(AlternativeSplicingExon hit : m_vcAlternativeSplicedExons)
		{
			if(hit.m_strExonGroup.equals(strQuery))
				return true;
		}
		
		return false;
	}
	
	public boolean ContainsRegion(String strRef, int nStart, int nEnd)
	{
		for(AlternativeSplicingExon hit : m_vcAlternativeSplicedExons)
		{
			String strCurRef 	= m_strRef;
			int nCurStart 		= hit.m_nStart;
			int nCurEnd			= hit.m_nEnd;
			
			if(strRef.equals(strCurRef))
			{
				if(nStart >= nCurStart && nEnd <= nCurEnd)
					return true;
			}
		}
		
		return false;
	}
	
	public boolean IsContainedInRegion(String strRef, int nStart, int nEnd)
	{
		for(AlternativeSplicingExon hit : m_vcAlternativeSplicedExons)
		{
			String strCurRef 	= m_strRef;
			int nCurStart 		= hit.m_nStart;
			int nCurEnd			= hit.m_nEnd;
			
			if(strRef.equals(strCurRef))
			{
				if(nStart <= nCurStart && nEnd >= nCurEnd)
					return true;
			}
		}
		
		return false;
	}
	
	public boolean ContainsExon(Exon ex)
	{
		return ContainsRegion(ex.getReference(), ex.getCodingStart(), ex.getCodingStop());
	}
	
	public String GetGeneID()
	{
		return m_strGeneID;
	}
	
	public String GetGeneSymbol()
	{
		return m_strGeneSymbol;
	}
	
	public String GetGTFFile()
	{
		return m_strFileGTF;
	}
	
	public int GetRating()
	{
		return m_nRating;
	}
	
	public String GetComment()
	{
		return m_strComment;
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
	
	public TreeSet<String> GetIsoforms()
	{
		return m_vcSelectedIsoforms ;
	}

	public int GetNumHits()
	{
		return m_vcAlternativeSplicedExons.size();
	}
	
	public String GetExonGroup(ExonGroup grp)
	{
		for(AlternativeSplicingExon hit : m_vcAlternativeSplicedExons)
		{
			if(grp.getGenomicStartOfGroup() == hit.m_nStart && grp.getGenomicStopOfGroup() == hit.m_nEnd)
				return hit.m_strExonGroup;
		}
		
		return null;
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
	
	public TreeSet<AlternativeSplicingExon> GetResultsForRegion(int nStart, int nEnd)
	{
		TreeSet<AlternativeSplicingExon> vcRes = new TreeSet<AlternativeSplicingExon>();
		
		for(AlternativeSplicingExon hit : m_vcAlternativeSplicedExons)
		{
			if(nStart > hit.m_nEnd || nEnd < hit.m_nStart)
				continue;

			vcRes.add(hit);
		}
		
		return vcRes;
	}
	
	public String GetDataForExonGroup(String strID)
	{
		String strRes = String.format(Locale.ENGLISH, "%s\t%s", m_strGeneID, m_strGeneSymbol);
		
		for(AlternativeSplicingExon hit : m_vcAlternativeSplicedExons)
		{
			if(!hit.m_strID.equals(strID))
				continue;

			strRes += String.format(Locale.ENGLISH, "\t%s\t%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4e", hit.GetTypeAsString(), hit.m_strExonGroup, hit.m_strConditionA, hit.m_strConditionB, hit.m_fAltExonCovPerBaseA, hit.m_fAltExonCovPerBaseB, Math.abs(hit.m_fAltExonCovPerBaseA-hit.m_fAltExonCovPerBaseB), hit.m_fFractionChangeAbsolute, hit.m_fFractionChangeRelative, hit.m_fPValue);
			break;
		}

		return strRes;
	}
	
	public AlternativeSplicingExon GetExonGroup(String strID)
	{
		for(AlternativeSplicingExon hit : m_vcAlternativeSplicedExons)
		{
			if(!hit.m_strID.equals(strID))
				continue;
			
			return hit;
		}
		
		return null;
	}
	
	@Override
	public String toString()
	{	
		String strRes = m_nRating + "\t" + m_strGeneID + "\t" + m_strGeneSymbol + "\t" + m_strFileGTF + "\t" + m_nMinJunctionReads + "\t" + m_nMinCovPerBase + "\t" + m_fMinCoveredBases + "\t" + m_fVariableExonThreshold + "\t" + m_strComment + "\n";
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
	
	public void SetRating(int nRating)
	{
		m_nRating = nRating;
	}
	
	public void SetIsoforms(TreeSet<String> vcIsoforms)
	{
		m_vcSelectedIsoforms = vcIsoforms;
	}
	
	public void SetComment(String strComment)
	{
		m_strComment = strComment;
	}

}
