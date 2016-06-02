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

import java.util.TreeSet;

public class TextFilter
{
	private TreeSet<String> m_vcStrings;
	private TreeSet<String> m_vcSelectedStrings;
	
	TextFilter(TreeSet<String> vcStrings)
	{			
		m_vcStrings			= vcStrings;
		m_vcSelectedStrings = new TreeSet<String>();

		for(String strString : vcStrings)
		{
			m_vcStrings.add(strString);
			m_vcSelectedStrings.add(strString);
		}
	}
	
	public void SelectString(String strString)
	{
		m_vcSelectedStrings.add(strString);
	}
	
	public void UnselectString(String strString)
	{
		m_vcSelectedStrings.remove(strString);
	}
	
	public void SelectAll()
	{
		for(String strString : m_vcStrings)
			m_vcSelectedStrings.add(strString);
	}
	
	public void UnselectAll()
	{
		m_vcSelectedStrings.clear();
	}
	
	public boolean IsStringSelected(String strString)
	{
		if(m_vcSelectedStrings.contains(strString))
			return true;
		
		return false;
	}

	public TreeSet<String> GetStrings()
	{
		return m_vcStrings;
	}
}
