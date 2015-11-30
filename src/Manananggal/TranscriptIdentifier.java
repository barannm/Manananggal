package Manananggal;

import java.util.Vector;

public class TranscriptIdentifier
{
	String m_strEnsemblGeneID;
	String m_strRefGeneID;
	
	TranscriptIdentifier()
	{
		m_strEnsemblGeneID 	= "?";
		m_strRefGeneID 		= "?";
	}
	
	public static Vector<TranscriptIdentifier> ParseFromLine(String strLine)
	{
		Vector<TranscriptIdentifier> vcResult = new Vector<TranscriptIdentifier>();
		
		String pSplit[] = strLine.split("\\s+");
		
		if(pSplit.length > 1)
		{
			String pSplit2[] = pSplit[1].split(",");		
			for(String strID : pSplit2)
			{
				TranscriptIdentifier tid = new TranscriptIdentifier();
				tid.m_strEnsemblGeneID = pSplit[0];
				tid.m_strRefGeneID = strID;
				
				vcResult.add(tid);
			}
		}
		else
		{
			TranscriptIdentifier tid = new TranscriptIdentifier();
			tid.m_strEnsemblGeneID 	 = pSplit[0];
			tid.m_strRefGeneID 		 = "";
			vcResult.add(tid);
		}
		
		return vcResult;
	}
	
	public boolean Equals(String strID)
	{
		if(m_strEnsemblGeneID.equals(strID))
			return true;
		
		if(m_strRefGeneID.equals(strID))
			return true;
		
		return false;
	}
	
	public String toString()
	{
		return m_strEnsemblGeneID + " " + m_strRefGeneID;
	}
}
