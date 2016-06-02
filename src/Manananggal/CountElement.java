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

public class CountElement implements Comparable<CountElement>
{
	String 	m_strRef;
	int 	m_nStart;
	int 	m_nEnd;
	String 	m_strGeneID;
	boolean m_bKnown;
	boolean m_bIsExon;
	char 	m_chStrand;
	
	Vector<Integer> m_vcCounts;
	
	CountElement()
	{
		m_strRef 	= "?";
		m_nStart 	= -1;
		m_nEnd		= -1;
		m_strGeneID	= "?";
		m_bKnown	= false;
		m_chStrand	= '?';
		m_bIsExon	= false;
		
		m_vcCounts = new Vector<Integer>();
	}
	
	@Override
	public boolean equals(Object o)
	{
		if(o instanceof CountElement)
		{			
			CountElement other = (CountElement) o;			
			
			if(!m_strRef.equals(other.m_strRef)) return false;
			if(m_chStrand != other.m_chStrand) return false;
			if(m_nStart != other.m_nStart || m_nEnd != other.m_nEnd) return false;
			
			return true;
		}
		else
			return false;
	}
	
	// use this function when reading from single counts files
	void ParseLine(String strLine)
	{
		String vcSplit[] = strLine.split("\t");
		if(vcSplit.length < 9)
		{
			System.out.println("Could not parse splice junction: invalid formated line: " + strLine);
			return;
		}
		
		m_strGeneID	= vcSplit[1].split("\\.")[0];
		
		if(vcSplit[2].equals("true"))
			m_chStrand = '+';
		else if(vcSplit[2].equals("false"))
			m_chStrand = '-';
		else
			m_chStrand = '?';
		
		m_strRef	= vcSplit[4];
		m_nStart	= Integer.parseInt(vcSplit[5]);
		m_nEnd		= Integer.parseInt(vcSplit[6]);			

		if(vcSplit[0].equals("Novel_Junction"))
			m_bKnown = false;
		else
			m_bKnown = true;
		
		m_vcCounts.addElement(Integer.parseInt(vcSplit[8]));
	}
	
	int GetStart() {return m_nStart;}
	int GetEnd() {return m_nEnd;}
	int GetLength() {return m_nEnd - m_nStart +1;}
	String GetRef() {return m_strRef;}
	String GetGeneID() {return m_strGeneID;}

	public int compareTo(CountElement other)
	{
		if(m_strRef.compareTo(other.m_strRef) != 0)
			return m_strRef.compareTo(other.m_strRef);
		
		int nDiff = m_nStart - other.GetStart();
		if(nDiff != 0)
			return nDiff;
		
		nDiff = m_nEnd - other.GetEnd();
		if(nDiff != 0)
			return nDiff;
			
		if(m_chStrand == other.m_chStrand)
			return 0;
		else if(m_chStrand == '+')
			return 1;
		else
			return -1;
	}
	
	public String toString()
	{
		return(m_strRef + "\t" + m_nStart + "\t" + m_nEnd + "\t" + m_chStrand + "\t" + m_bKnown + "\t" + m_bIsExon + "\t" + m_vcCounts);
	}
	
	public String toSimpleString()
	{
		return(m_strRef + "\t" + m_nStart + "\t" + m_nEnd + "\t" + m_chStrand + "\t" + m_bKnown);
	}
	
	public String toOutputString()
	{
		String strOut = m_strGeneID + "\t" + m_strRef + "\t" + m_nStart + "\t" + m_nEnd + "\t" + m_chStrand + "\t" + m_bKnown;
		
		return strOut;
	}
}