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

import java.util.HashSet;

public interface MappableElement 
{
	public String getID();
	public int getIntegerID();
	public int getCodingStart();
	public int getCodingStop();
	public String getReference();
	public void setReference(String reference);
	public String getDescription();
	public void addMappedExon(Exon e);
	public void addMappedTranscript(String transcript);
	public HashSet<Exon> getMappedExons();
	public HashSet<String> getMappedTranscripts();
}
