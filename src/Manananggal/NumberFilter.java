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

public class NumberFilter
{
	private float m_fMinValue;
	private float m_fMaxValue;
	
	public NumberFilter()
	{
		m_fMinValue = Float.NEGATIVE_INFINITY;
		m_fMaxValue	= Float.POSITIVE_INFINITY;
	}
	
	public void SetMinValue(float fValue)
	{
		m_fMinValue = fValue;
	}
	
	public void SetMaxValue(float fValue)
	{
		m_fMaxValue = fValue;
	}
	
	public double GetMinValue()
	{
		return m_fMinValue;
	}
	
	public double GetMaxValue()
	{
		return m_fMaxValue;
	}
	
	public void Disable()
	{
		m_fMinValue = Float.NEGATIVE_INFINITY;
		m_fMaxValue	= Float.POSITIVE_INFINITY;
	}
	
	public boolean IsInRange(double fVal)
	{		
		if(fVal >= m_fMinValue && fVal <= m_fMaxValue)
			return true;

		return false;
	}
	
	public boolean IsInRange(int nVal)
	{
		if(nVal >= m_fMinValue && nVal <= m_fMaxValue)
			return true;
		
		return false;
	}
}
