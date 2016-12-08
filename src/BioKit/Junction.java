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

public class Junction extends GeneticElement implements MappableElement 
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private int junctionID;
	private String reference;
	private HashSet<Exon> surroundingExons;
	private HashSet<String> transcriptsContainingJunction;
	
	public Junction(String reference, int junctionStart, int junctionStop)
	{
		super(junctionStart, junctionStop, junctionStart, junctionStop, true);
		this.reference = reference;
		this.surroundingExons = new HashSet<Exon>();
		this.transcriptsContainingJunction = new HashSet<String>();
	}
		
	@Override
	public String getID() 
	{
		return getDescription();
	}

	@Override
	public int getIntegerID() 
	{
		return junctionID;
	}

	@Override
	public int getCodingStart() 
	{
		return super.getCodingStart();
	}

	@Override
	public int getCodingStop() 
	{
		return super.getCodingStop();
	}

	@Override
	public String getReference() 
	{
		return reference;
	}

	@Override
	public void setReference(String reference) 
	{
		this.reference = reference;
	}
	
	public void setID(int id)
	{
		this.junctionID = id;
	}
		
	@Override
	/* This method is needed for the SAMBAMCompression method.
	 */
	public String getDescription() 
	{
		return reference + ":" + getCodingStart() + "," + reference + ":" + getCodingStop();
	}

	@Override
	public HashSet<Exon> getMappedExons() 
	{

		return this.surroundingExons;
	}

	@Override
	public HashSet<String> getMappedTranscripts() 
	{
		return this.transcriptsContainingJunction;
	}

	@Override
	public void addMappedExon(Exon e) 
	{
		this.surroundingExons.add(e);	
	}

	@Override
	public void addMappedTranscript(String transcript) 
	{
		this.transcriptsContainingJunction.add(transcript);
	}
}

