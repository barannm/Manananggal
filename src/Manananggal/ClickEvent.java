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

import BioKit.Exon;

public class ClickEvent
{
	String m_strClickedA;	// name of first clicked element
	String m_strClickedB;	// name of second clicked element
	
	TreeSet<Exon> m_vcSelectedExons;
	TreeSet<String> m_vcTargetIsoforms;
	
	int m_nExonGroupStartA;
	int m_nExonGroupStartB;
	int m_nExonGroupEndA;
	int m_nExonGroupEndB;
	
	ClickEvent()
	{
		m_strClickedA = null;
		m_strClickedB = null;

		m_nExonGroupStartA	= -1;
		m_nExonGroupStartB	= -1;
		m_nExonGroupEndA	= -1;
		m_nExonGroupEndB	= -1;
		
		m_vcSelectedExons = new TreeSet<Exon>();
		m_vcTargetIsoforms = new TreeSet<String>();
	}
	
	void clear()
	{
		m_strClickedA = null;
		m_strClickedB = null;

		m_nExonGroupStartA	= -1;
		m_nExonGroupStartB	= -1;
		m_nExonGroupEndA	= -1;
		m_nExonGroupEndB	= -1;
		
		m_vcSelectedExons.clear();
		m_vcTargetIsoforms.clear();
	}
}
