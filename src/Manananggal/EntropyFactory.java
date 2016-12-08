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

import java.util.Collections;
import java.util.Map;
import java.util.TreeMap;
import java.util.Vector;

import org.apache.commons.math3.stat.StatUtils;

/**
 *    The Entropy factory calculates entropy values for tissue specific data.
 *    By default, it will use the Gini index to identify genes or exonic parts
 *    that are tissue specific, but other entropy functions for the Gini
 *    Simpson index, Shannon index and Theil index were also implemented.
 */
public class EntropyFactory
{
	public TreeMap<String, Double> CalculateGiniIndexForExons(TreeMap<String, TreeMap<String, Vector<Double>>> mapCountsPerExonAndTissue)
	{
		TreeMap<String, Double> mapEntropyToExonicParts = new TreeMap<String, Double>();
		
		for(String strPos : mapCountsPerExonAndTissue.keySet())
		{
			Vector<Double> vcMedianPerTissue = new Vector<Double>();
			for(Map.Entry<String, Vector<Double>> e : mapCountsPerExonAndTissue.get(strPos).entrySet())
			{
				double pValues[] = new double[e.getValue().size()];
				int i=0;
				for(double fVal : e.getValue())
				{
					pValues[i] = fVal;
					i++;
				}
				
				double fMedian = StatUtils.percentile(pValues, 50.0);
				vcMedianPerTissue.add(fMedian);
			}
			
			// sort by ascending order
			Collections.sort(vcMedianPerTissue);
			
			// calculate gini index
			int nTissues = vcMedianPerTissue.size();
			double fCountSum = 0.0;
			double fGiniSum  = 0.0;
			for(int i=0; i<nTissues; i++)
			{
				double fVal = vcMedianPerTissue.get(i);
				fGiniSum  += (nTissues+1-(i+1)) * fVal;
				fCountSum += fVal;
			}
			
			if(fGiniSum == 0.0)
			{
				mapEntropyToExonicParts.put(strPos, 0.0);
			}
			else
			{
				double fGini = (nTissues+1)/nTissues - (2*fGiniSum) / (nTissues*fCountSum);
				mapEntropyToExonicParts.put(strPos, fGini);
			}
		}
		
		return mapEntropyToExonicParts;
	}
	
	public TreeMap<String, Double> CalculateGiniSimpsonIndexForExons(TreeMap<String, TreeMap<String, Vector<Double>>> mapCountsPerExonAndTissue)
	{
		TreeMap<String, Double> mapEntropyToExonicParts = new TreeMap<String, Double>();
		
		for(String strPos : mapCountsPerExonAndTissue.keySet())
		{
			int nTotalCount = 0;
			TreeMap<Double, String> sumPerTissue = new TreeMap<Double, String>();
			for(Map.Entry<String, Vector<Double>> e : mapCountsPerExonAndTissue.get(strPos).entrySet())
			{
				double pValues[] = new double[e.getValue().size()];
				int i=0;
				for(double fVal : e.getValue())
				{
					pValues[i] = fVal;
					i++;
					
					nTotalCount += fVal;
				}
				
				double fSum = StatUtils.sum(pValues);
				sumPerTissue.put(fSum, e.getKey());
			}
			
			double fDiversity = 0.0;
			for(double fSum : sumPerTissue.keySet())
			{
				fDiversity += (fSum/nTotalCount)*(fSum/nTotalCount);
			}
			
			mapEntropyToExonicParts.put(strPos, fDiversity);
		}
		
		return mapEntropyToExonicParts;
	}

	public TreeMap<String, Double> CalculateShannonIndexForExons(TreeMap<String, TreeMap<String, Vector<Double>>> mapCountsPerExonAndTissue)
	{
		TreeMap<String, Double> mapEntropyToExonicParts = new TreeMap<String, Double>();
		
		for(String strPos : mapCountsPerExonAndTissue.keySet())
		{
			TreeMap<Double, String> medianPerTissue = new TreeMap<Double, String>();
			for(Map.Entry<String, Vector<Double>> e : mapCountsPerExonAndTissue.get(strPos).entrySet())
			{
				double pValues[] = new double[e.getValue().size()];
				int i=0;
				for(double fVal : e.getValue())
				{
					pValues[i] = fVal;
					i++;
				}
				
				double fMedian = StatUtils.percentile(pValues, 50.0);
				medianPerTissue.put(fMedian, e.getKey());
			}
			
			// calculate shannon index
			// get sum for all tissues
			double fSum = 0.0;
			for(Double fVal : medianPerTissue.keySet())
			{
				fSum += fVal;
			}

			// get a sum of fraction * ln fraction
			double fSum2 = 0.0;
			for(Double fVal : medianPerTissue.keySet())
			{
				double fPI = (fVal/fSum);
				if(fPI != 0)
					fSum2 += fPI * Math.log(fPI);
			}
			
			int nTissues = medianPerTissue.size();

			// convert to eveneness
			fSum2 = Math.exp(-fSum2) / nTissues;

			mapEntropyToExonicParts.put(strPos, 1-fSum2);
		}
		
		return mapEntropyToExonicParts;
	}

	public TreeMap<String, Double> CalculateTheilIndexForExons(TreeMap<String, TreeMap<String, Vector<Double>>> mapCountsPerExonAndTissue)
	{
		TreeMap<String, Double> mapEntropyToExonicParts = new TreeMap<String, Double>();
		
		for(String strPos : mapCountsPerExonAndTissue.keySet())
		{
			double fEntropy = 0;
			
			int nTotalSamples = 0;
			double fSumAllSamples = 0;
			for(String strTissue : mapCountsPerExonAndTissue.get(strPos).keySet())
			{
				Vector<Double> vcValues = mapCountsPerExonAndTissue.get(strPos).get(strTissue);
				nTotalSamples += vcValues.size();
				
				for(double fVal : vcValues)
					fSumAllSamples += fVal;
			}
			double fMeanAllSamples = fSumAllSamples / nTotalSamples;
				
			for(Map.Entry<String, Vector<Double>> e : mapCountsPerExonAndTissue.get(strPos).entrySet())
			{
				double pValues[] = new double[e.getValue().size()];
				int i=0;
				for(double fVal : e.getValue())
				{
					pValues[i] = fVal;
					i++;
				}
				
				int nSamples = pValues.length;
				double fMean = StatUtils.mean(pValues);

				double fA = (double)nSamples / (double)nTotalSamples;
				double fB = fMean / fMeanAllSamples;
				double fC = 0.0;
				if(fB != 0.0)
					fC = Math.log(fB);
				
				fEntropy += fA * fB * fC;
			}
			
			// normalize theils index
			fEntropy = 1 - Math.exp(-fEntropy);

			mapEntropyToExonicParts.put(strPos, fEntropy);
		}
		
		return mapEntropyToExonicParts;
	}
}
