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
