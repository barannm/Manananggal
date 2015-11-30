package Manananggal;

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
		return m_strApprovedGeneSymbol.compareTo(other.m_strApprovedGeneSymbol);
	}
}
