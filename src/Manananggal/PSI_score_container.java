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
import java.util.Vector;

import BioKit.Gene;

/**
 *    The PSI_score_container includes all PSI score results
 *    that were generated during the analysis of the current gene.
 */
public class PSI_score_container
{
	TreeSet<PSI_score_result> 	m_vcResults;
	String						m_strGeneID;
	String 						m_strGeneSymbol;
	
	public PSI_score_container()
	{
		m_vcResults = new TreeSet<PSI_score_result>();
	}
	
	public void AddResult(TreeSet<CountElement> vcPathA, TreeSet<CountElement> vcPathB, double fPValue, Gene gene, Vector<String> vcConditionsInOrder, Vector<double[]> vcScoresPerCondition)
	{
		PSI_score_result newRes = new PSI_score_result(vcPathA, vcPathB, fPValue, gene, vcConditionsInOrder, vcScoresPerCondition);
		
		m_strGeneID 	= gene.getGeneID().split("\\.")[0];
		m_strGeneSymbol = gene.getGeneName();
		
		boolean bContained = false;
		for(PSI_score_result res : m_vcResults)
		{
			if(res.equals(newRes))
			{
				bContained = true;
				break;
			}
		}

		if(!bContained)
		{
			m_vcResults.add(newRes);
		}
	}
	
	public void Report()
	{
		for(PSI_score_result r : m_vcResults)
		{
			r.Print();
		}
	}
	
	public void Clear()
	{
		m_vcResults.clear();
	}
	
	public void Print()
	{
		for(PSI_score_result res : m_vcResults)
			res.Print();
	}
}
