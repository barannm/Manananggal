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

/**
 *    The GeneIdentifier class stores various gene IDs for
 *    a gene that are used for cross referencing. Valid
 *    gene identifiers are included in the cross reference
 *    file AND the GTF file.
 */
public class GeneIdentifier implements Comparable<GeneIdentifier>
{
	String m_strEntrezGeneID;
	String m_strEnsemblGeneID;
	String m_strEnsemblTranscriptID;
	String m_strRefGeneID;	
	String m_strApprovedGeneSymbol;
	String m_strApprovedGeneName;
	String m_strSynonyms;
	boolean m_bIsValid;
	
	GeneIdentifier()
	{
		m_strEntrezGeneID 		= "?";
		m_strEnsemblTranscriptID= "?";
		m_strEnsemblGeneID 		= "?";
		m_strRefGeneID 			= "?";	
		m_strApprovedGeneSymbol	= "?";
		m_strApprovedGeneName 	= "?";
		m_strSynonyms 			= "?";
		m_bIsValid				= false;
	}
	
	public static GeneIdentifier ParseFromLine(String strLine)
	{
		GeneIdentifier res = new GeneIdentifier();
		String split[] = strLine.split("\t", -1);
		
		if(split.length != 11)
		{
			System.out.println("warning: invalid line: " + split.length + " " + strLine);
			return null;
		}
		
		res.m_strApprovedGeneSymbol	= split[1];
		res.m_strApprovedGeneName 	= split[2];
		res.m_strSynonyms			= split[5];
		res.m_strEntrezGeneID		= split[8];
		res.m_strEnsemblGeneID		= split[9];
		res.m_strRefGeneID			= split[10];		
		
		return res;
	}
	
	public static GeneIdentifier ParseFromLineBioMart(String strLine)
	{
		GeneIdentifier res = new GeneIdentifier();
		String split[] = strLine.split("\t", -1);
		
		if(split.length < 6)
		{
			System.out.println("warning: invalid line: " + split.length + " " + strLine);
			return null;
		}
				
		res.m_strEnsemblGeneID		= split[0];
		res.m_strEnsemblTranscriptID= split[1];
		res.m_strRefGeneID			= split[2];
		res.m_strEntrezGeneID		= split[3];
		res.m_strApprovedGeneSymbol	= split[4];		
		res.m_strApprovedGeneName	= split[5];
		
		if(split.length == 7)
			res.m_strSynonyms		= split[6];
		
		return res;
	}
	
	public boolean EqualsGene(String strID)
	{
		strID = strID.trim().toUpperCase();
		
		if(m_strApprovedGeneSymbol.toUpperCase().equals(strID))
			return true;
		
		if(m_strEnsemblGeneID.toUpperCase().equals(strID))
			return true;
		
		if(m_strEntrezGeneID.toUpperCase().equals(strID))
			return true;
		
		if(m_strRefGeneID.toUpperCase().equals(strID))
			return true;
		
		if(m_strEnsemblTranscriptID.toUpperCase().equals(strID))
			return true;
		
		/*
		 * Don't search for synonyms, this causes problems! e.g. for HP
		String split[] = m_strSynonyms.split(",");
		for(String strSyn : split)
		{
			if(strSyn.trim().equals(strID))
				return true;
		}
		*/
			
		return false;
	}
	
	public String toString()
	{
		return m_strApprovedGeneSymbol + "\t" + m_strEntrezGeneID + "\t" + m_strEnsemblGeneID + "\t" + m_strEnsemblTranscriptID + "\t" + m_strRefGeneID + "\t" + m_strSynonyms + "\t" + m_strApprovedGeneName;
	}

	@Override
	public int compareTo(GeneIdentifier other)
	{
		int nComparison = m_strEnsemblGeneID.compareTo(other.m_strEnsemblGeneID);
		
		if(nComparison != 0)
			return nComparison;
		
		nComparison = m_strEnsemblTranscriptID.compareTo(other.m_strEnsemblTranscriptID);
		
		return nComparison;
	}
}
