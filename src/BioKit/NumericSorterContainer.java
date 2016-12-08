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
 * @author Fabian Birzele
 */
package BioKit;

public class NumericSorterContainer<E> implements Comparable<NumericSorterContainer<E>>
{
	private double value;
	private E content;
	boolean ignoreEqual;
	
	public NumericSorterContainer(double value, E content, boolean ignoreEqual)
	{
		this.value = value;
		this.content = content;
		this.ignoreEqual = ignoreEqual;
	}
	
	public E getContent()
	{
		return content;
	}
	
	public void setContent(E content)
	{
		this.content = content;
	}
	
	public double getValue()
	{
		return value;
	}

	public int compareTo(NumericSorterContainer<E> o)
	{
		if(value > o.value)
			return 1;
		else if(value == o.value && !ignoreEqual)
			return 0;
		else
			return -1;
	}	
}
