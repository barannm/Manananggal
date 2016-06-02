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
import java.io.File;
import java.util.Scanner;
import java.util.TreeSet;
import java.util.Vector;

import BioKit.Gene;


public class GeneIdentifierHandler
{
	Vector<GeneIdentifier>		m_vcGeneIdentifiers;
	Vector<GeneIdentifierRange> m_vcRangesEnsemblIDs;
	Vector<GeneIdentifierRange> m_vcRangesSymbols;
	
	private class GeneIdentifierRange
	{
		String m_strFirstString;
		String m_strLastString;
		Vector<Integer> m_vcIndices;
		
		private GeneIdentifierRange()
		{
			m_strFirstString	= new String();
			m_strLastString		= new String();
			m_vcIndices			= new Vector<Integer>();
		}
	}
	
	public GeneIdentifierHandler()
	{
		m_vcGeneIdentifiers 	= null; 
		m_vcRangesEnsemblIDs 	= null;
		m_vcRangesSymbols		= null;
	}
	
	public void Init(String strFile)
	{
		m_vcGeneIdentifiers		= new Vector<GeneIdentifier>(); 
		m_vcRangesEnsemblIDs	= new Vector<GeneIdentifierRange>();
		m_vcRangesSymbols		= new Vector<GeneIdentifierRange>();
		
		File pFile = new File(strFile);
		
		if(!pFile.exists())
		{
			System.out.println("ERROR: File does not exist: " + strFile);
			return;
		}
		
		Scanner pScanner = null;
		try
		{
			pScanner = new Scanner(pFile);
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}

		// read all gene identifiers, use a set for sorting...
		TreeSet<GeneIdentifier> vcGIDs = new TreeSet<GeneIdentifier>(); 
		while(pScanner.hasNextLine())
		{
			String strLine = pScanner.nextLine();

			GeneIdentifier gid = GeneIdentifier.ParseFromLineBioMart(strLine);
			if(gid != null && !gid.m_strApprovedGeneSymbol.contains("withdrawn"))
			{
				vcGIDs.add(gid);
			}
		}
		pScanner.close();
		
		//... then convert to vector
		m_vcGeneIdentifiers = new Vector<GeneIdentifier>(vcGIDs);

		// index gene identifiers, group 1000 together (ensemblID) or split by first letter (gene symbols)
		GeneIdentifierRange currentRange = null;
		
		int nIdx = 0;
		for(GeneIdentifier gid : m_vcGeneIdentifiers)
		{
			// group size = 1000 gene identifiers
			if(nIdx % 1000 == 0)
			{
				currentRange = new GeneIdentifierRange();
				currentRange.m_strFirstString	= gid.m_strEnsemblGeneID;
				m_vcRangesEnsemblIDs.add(currentRange);
			}
			
			currentRange.m_strLastString = gid.m_strEnsemblGeneID;
			currentRange.m_vcIndices.add(nIdx);
			
			// also index by first letter of the gene symbol
			if(!gid.m_strApprovedGeneSymbol.isEmpty())
			{
				boolean bAdded = false;
				for(GeneIdentifierRange range : m_vcRangesSymbols)
				{
					if(range.m_strFirstString.charAt(0) == gid.m_strApprovedGeneSymbol.charAt(0))
					{
						range.m_vcIndices.add(nIdx);
						bAdded = true;
					}
				}
				
				if(!bAdded)
				{
					GeneIdentifierRange newRange = new GeneIdentifierRange();
					newRange.m_strFirstString = "" + gid.m_strApprovedGeneSymbol.charAt(0);
					newRange.m_strLastString  = "" + gid.m_strApprovedGeneSymbol.charAt(0);
					newRange.m_vcIndices.add(nIdx);
					m_vcRangesSymbols.add(newRange);
				}
			}
			
			nIdx++;
		}
		
		System.out.println("read " + nIdx + " gene identifiers from cross reference file");
	}
	
	// search for gene identifiers (strID may be the gene ID or gene symbol)
	// - 'gene' is an optional parameter and may be null
	public GeneIdentifier GetGeneIdentifierForGene(String strID, Gene gene)
	{		
		// search by ID
		for(GeneIdentifierRange range : m_vcRangesEnsemblIDs)
		{
			int nTestA = range.m_strFirstString.compareTo(strID.split("\\.")[0]);
			int nTestB = range.m_strLastString.compareTo(strID.split("\\.")[0]);
			
			if(nTestA <= 0 && nTestB >= 0)
			{
				// test by exact gene ID
				for(int nIdx : range.m_vcIndices)
				{
					GeneIdentifier gid =  m_vcGeneIdentifiers.get(nIdx);
					
					if(gid.EqualsGene(strID) && gid.m_bIsValid)
					{
						return gid;
					}
				}
				
				// test again, but remove any version numbers
				for(int nIdx : range.m_vcIndices)
				{
					GeneIdentifier gid =  m_vcGeneIdentifiers.get(nIdx);					
					
					if(gid.EqualsGene(strID.split("\\.")[0]) && gid.m_bIsValid)
					{
						return gid;
					}
				}
			}
		}
		
		// if no identifier was returned yet, try by gene symbols
		for(GeneIdentifierRange range : m_vcRangesSymbols)
		{
			if(range.m_strFirstString.charAt(0) == strID.charAt(0))
			{
				for(int nIdx : range.m_vcIndices)
				{
					GeneIdentifier gid =  m_vcGeneIdentifiers.get(nIdx);
					
					if(gid.EqualsGene(strID) && gid.m_bIsValid)
					{				
						return gid;
					}
				}
			}
		}
		
		if(gene != null)
		{
			GeneIdentifier gid = new GeneIdentifier();
			gid.m_strEntrezGeneID = "?";
			gid.m_strEnsemblTranscriptID 	= gene.getGeneID().split("\\.")[0];
			gid.m_strApprovedGeneName		= gene.getGeneID().split("\\.")[0];
			gid.m_strApprovedGeneSymbol 	= gene.getGeneName();
			return gid;
		}

		return null;
	}
	
	// this function 'requires' the gene parameter
	public GeneIdentifier GetGeneIdentifierForTranscript(String strIsoformID, Gene gene)
	{
		// search by ID
		for(GeneIdentifierRange range : m_vcRangesEnsemblIDs)
		{
			int nTestA = range.m_strFirstString.compareTo(gene.getGeneID().split("\\.")[0]);
			int nTestB = range.m_strLastString.compareTo(gene.getGeneID().split("\\.")[0]);
			
			if(nTestA <= 0 && nTestB >= 0)
			{
				// test by exact gene ID
				for(int nIdx : range.m_vcIndices)
				{
					GeneIdentifier gid =  m_vcGeneIdentifiers.get(nIdx);
					
					if(gid.EqualsGene(gene.getGeneID()) && gid.m_bIsValid)
					{
						if(gid.m_strEnsemblTranscriptID.equals(strIsoformID.split("\\.")[0]))
							return gid;
					}
				}
				
				// test again, but remove any version numbers
				for(int nIdx : range.m_vcIndices)
				{
					GeneIdentifier gid =  m_vcGeneIdentifiers.get(nIdx);					
					
					if(gid.EqualsGene(gene.getGeneID().split("\\.")[0]) && gid.m_bIsValid)
					{
						if(gid.m_strEnsemblTranscriptID.equals(strIsoformID.split("\\.")[0]))
						{
							return gid;
						}
					}
				}
			}
		}
		
		// if no identifier was returned yet, try by gene symbols
		for(GeneIdentifierRange range : m_vcRangesSymbols)
		{
			if(range.m_strFirstString.charAt(0) == gene.getGeneName().charAt(0))
			{
				for(int nIdx : range.m_vcIndices)
				{
					GeneIdentifier gid =  m_vcGeneIdentifiers.get(nIdx);
					
					if(gid.EqualsGene(gene.getGeneName()) && gid.m_bIsValid)
					{
						if(gid.m_strApprovedGeneSymbol.equals(strIsoformID))
							return gid;
					}
				}
			}
		}
		
		
		GeneIdentifier gid = new GeneIdentifier();
		gid.m_strEntrezGeneID = "?";
		gid.m_strEnsemblTranscriptID 	= strIsoformID.split("\\.")[0];
		gid.m_strApprovedGeneName		= strIsoformID.split("\\.")[0];
		
		if(gene != null)
		{
			gid.m_strApprovedGeneSymbol 	= gene.getGeneName();
		}
		return gid;
	}
	
	public Vector<GeneIdentifier> GetAllGeneIdentifiers()
	{
		return m_vcGeneIdentifiers;
	}

	public Vector<GeneIdentifier> GetAllGeneIdentifiersStartingWith(String strQuery)
	{
		Vector<GeneIdentifier> vcResult = new Vector<GeneIdentifier>();
		
		TreeSet<String> vcIDsAdded = new TreeSet<String>();
		
		// search by ID
		for(GeneIdentifierRange range : m_vcRangesEnsemblIDs)
		{
			int nTestA = range.m_strFirstString.compareTo(strQuery.split("\\.")[0]);
			int nTestB = range.m_strLastString.compareTo(strQuery.split("\\.")[0]);
			
			if(nTestA <= 0 && nTestB >= 0)
			{
				// test by exact gene ID
				for(int nIdx : range.m_vcIndices)
				{
					GeneIdentifier gid =  m_vcGeneIdentifiers.get(nIdx);
					
					if(vcIDsAdded.contains(gid.m_strEnsemblGeneID))
						continue;
					
					if(gid.EqualsGene(strQuery) && gid.m_bIsValid)
					{
						vcIDsAdded.add(gid.m_strEnsemblGeneID);
						vcResult.add(gid);
					}
				}
				
				// test again, but remove any version numbers
				for(int nIdx : range.m_vcIndices)
				{
					GeneIdentifier gid =  m_vcGeneIdentifiers.get(nIdx);
					
					if(vcIDsAdded.contains(gid.m_strEnsemblGeneID))
						continue;
					
					if(gid.EqualsGene(strQuery.split("\\.")[0]) && gid.m_bIsValid)
					{
						vcIDsAdded.add(gid.m_strEnsemblGeneID);
						vcResult.add(gid);
					}
				}
			}
		}
		
		for(GeneIdentifierRange range : m_vcRangesSymbols)
		{
			if(range.m_strFirstString.charAt(0) == strQuery.charAt(0))
			{
				for(int nIdx : range.m_vcIndices)
				{
					GeneIdentifier gid =  m_vcGeneIdentifiers.get(nIdx);
					
					if(vcIDsAdded.contains(gid.m_strEnsemblGeneID))
						continue;
					
					if(gid.m_bIsValid && gid.m_strApprovedGeneSymbol.startsWith(strQuery))
					{
						vcIDsAdded.add(gid.m_strEnsemblGeneID);
						vcResult.add(gid);
					}
				}
			}
		}

		return vcResult;
	}
	
	public void EnableAll()
	{
		for(GeneIdentifier gid : m_vcGeneIdentifiers)
			gid.m_bIsValid = true;
	}
}
