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

import java.io.IOException;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.TTest;

import BioKit.Exon;
import BioKit.ExonGroup;

/**
 *    This class is used to identify alternatively spliced exon extensions
 *    The detection is solely based on changes in PSI scores.
 */
public class AnalyzerExonExtensions
{
	public TreeSet<SimpleSpliceScore> IdentifyExonExtensionEvents(SplicingWebApp app, boolean bDebug, TreeSet<SimpleSpliceScore> results) throws IOException
	{	
		//###################################################################
		//   retrieve data for the current selections and gene information
		//###################################################################
		TreeSet<String> vcValidIsoforms 	= app.GetValidIsoforms();
		String strSelectedConditionType 	= app.GetSelectedConditionType();
		TreeSet<String> vcSelectedSamples 	= app.GetSelectedSamples();		
		ProjectModel projectModel			= app.GetProjectModel();
		DataSupplier data					= app.GetDataSupplier();
		//###################################################################
		
		// create new result container if necessary
		if(results == null)
			results = new TreeSet<SimpleSpliceScore>();
		
		// get samples per Condition
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = projectModel.GetSelectedSamplesPerCondition(strSelectedConditionType, vcSelectedSamples);
		
		// get a list of all junctions that are included in the junction count file
		TreeSet<CountElement> vcJunctionsWithCoverage = data.GetJunctions();

		// create new exon groups, discarding retained introns
		// add all novel junctions to the valid junctions
		TreeSet<CountElement> vcValidJunctions = new TreeSet<CountElement>();
		for(CountElement e : vcJunctionsWithCoverage)
		{
			if(!e.m_bKnown)
				vcValidJunctions.add(e);
		}
		TreeSet<Exon> vcExons = new TreeSet<Exon>();
		for(String strIsoform : vcValidIsoforms)
		{
			Exon pIsoformExons[] = data.GetExonsForIsoform(strIsoform);
			
			for(int i=0; i<pIsoformExons.length; i++)
			{
				Exon ex = pIsoformExons[i];
				
				// skip retained introns
				if(data.IsRetainedIntron(ex))
					continue;

				vcExons.add(ex);
			}
			
			// get splice junctions
			TreeSet<CountElement> vcJunctions =	data.GetJunctionsForIsoform(strIsoform);
			for(CountElement e : vcJunctions)
			{
				if(vcJunctionsWithCoverage.contains(e))
					vcValidJunctions.add(e);
			}
		}
		ExonGroup pGroups[] = data.RecalculateExonGroups(vcExons);

		//################################
		//   group junctions by exons
		//################################
		TreeMap<ExonGroup, Vector<CountElement>> map5PrimeJunctionsToExons = new TreeMap<ExonGroup, Vector<CountElement>>();
		TreeMap<ExonGroup, Vector<CountElement>> map3PrimeJunctionsToExons = new TreeMap<ExonGroup, Vector<CountElement>>();
				
		for(CountElement e : vcValidJunctions)
		{
			// check whether the junction is valid (passes the threshold filter)
			if(data.IsInvalidJunction(e))
				continue;
			
			for(ExonGroup ex : pGroups)
			{
				if(e.m_nStart >= ex.getGenomicStartOfGroup() && e.m_nStart <= ex.getGenomicStopOfGroup())
				{
					if(data.GetStrand() == '+' && e.m_nEnd > ex.getGenomicStopOfGroup())
					{						
						if(map3PrimeJunctionsToExons.containsKey(ex))
						{
							Vector<CountElement> vcTmp = map3PrimeJunctionsToExons.get(ex);
							vcTmp.add(e);
							
							map3PrimeJunctionsToExons.put(ex, vcTmp);
						}
						else
						{
							Vector<CountElement> vcTmp = new Vector<CountElement>();
							vcTmp.add(e);
							
							map3PrimeJunctionsToExons.put(ex, vcTmp);
						}
					}
					else if(data.GetStrand() == '-' && e.m_nEnd > ex.getGenomicStopOfGroup())
					{
						if(map5PrimeJunctionsToExons.containsKey(ex))
						{
							Vector<CountElement> vcTmp = map5PrimeJunctionsToExons.get(ex);
							vcTmp.add(e);
							
							map5PrimeJunctionsToExons.put(ex, vcTmp);
						}
						else
						{
							Vector<CountElement> vcTmp = new Vector<CountElement>();
							vcTmp.add(e);
							
							map5PrimeJunctionsToExons.put(ex, vcTmp);
						}
					}
				}
				
				if(e.m_nEnd >= ex.getGenomicStartOfGroup() && e.m_nEnd <= ex.getGenomicStopOfGroup())
				{
					if(data.GetStrand() == '+' && e.m_nStart < ex.getGenomicStartOfGroup())
					{						
						if(map5PrimeJunctionsToExons.containsKey(ex))
						{
							Vector<CountElement> vcTmp = map5PrimeJunctionsToExons.get(ex);
							vcTmp.add(e);
							
							map5PrimeJunctionsToExons.put(ex, vcTmp);
						}
						else
						{
							Vector<CountElement> vcTmp = new Vector<CountElement>();
							vcTmp.add(e);
							
							map5PrimeJunctionsToExons.put(ex, vcTmp);
						}
					}
					else if(data.GetStrand() == '-' && e.m_nStart < ex.getGenomicStartOfGroup())
					{
						if(map3PrimeJunctionsToExons.containsKey(ex))
						{
							Vector<CountElement> vcTmp = map3PrimeJunctionsToExons.get(ex);
							vcTmp.add(e);
							
							map3PrimeJunctionsToExons.put(ex, vcTmp);
						}
						else
						{
							Vector<CountElement> vcTmp = new Vector<CountElement>();
							vcTmp.add(e);
							
							map3PrimeJunctionsToExons.put(ex, vcTmp);
						}
					}
				}
			}
		}
		
		TreeSet<ExonGroup> vcInvalidGroup = new TreeSet<ExonGroup>();
		for(ExonGroup grp : map3PrimeJunctionsToExons.keySet())
		{
			if(map3PrimeJunctionsToExons.get(grp).size() < 2)
				vcInvalidGroup.add(grp);
		}
		for(ExonGroup grp : vcInvalidGroup)
			map3PrimeJunctionsToExons.remove(grp);
		
		vcInvalidGroup.clear();
		for(ExonGroup grp : map5PrimeJunctionsToExons.keySet())
		{
			if(map5PrimeJunctionsToExons.get(grp).size() < 2)
				vcInvalidGroup.add(grp);
		}
		for(ExonGroup grp : vcInvalidGroup)
			map5PrimeJunctionsToExons.remove(grp);
/*	
		System.out.println("junction to exons map");
		System.out.println("3'");
		for(ExonGroup grp : map3PrimeJunctionsToExons.keySet())
		{
			System.out.println(grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup());
			System.out.println(map3PrimeJunctionsToExons.get(grp));
		}
		System.out.println("5'");
		for(ExonGroup grp : map5PrimeJunctionsToExons.keySet())
		{
			System.out.println(grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup());
			System.out.println(map5PrimeJunctionsToExons.get(grp));
		}
*/
		// get number of conditions
		int nConditions = mapSamplesToConditions.keySet().size();
		
		String pConditions[] = new String[nConditions];
		mapSamplesToConditions.keySet().toArray(pConditions);
		
		// list of all p-values
		Vector<Double> vcPValues = new Vector<Double>();

		// proceed through all exons
		for(int k=0; k<2; k++)
		{
			TreeMap<ExonGroup, Vector<CountElement>> mapJunctionsToExons;
			
			if(k == 0)
			{
				mapJunctionsToExons = map5PrimeJunctionsToExons;
			}
			else
			{
				mapJunctionsToExons = map3PrimeJunctionsToExons;
			}
			
			for(ExonGroup ex : mapJunctionsToExons.keySet())
			{
				Vector<CountElement> vcTmp = mapJunctionsToExons.get(ex);
				int nJunctions = vcTmp.size();
				
				for(int i=0; i<nJunctions-1; i++)
				{
					for(int j=i+1; j<nJunctions; j++)
					{
						CountElement jun1 = vcTmp.get(i);
						CountElement jun2 = vcTmp.get(j);
						
//						if(k == 1)
//							System.out.println("testing: " + jun1 + " vs. " + jun2);

						// make sure both ends are linked to the same exons while this ignores some special
						// cases of exon extensions, it will reduce the number of false positives that would
						// be detected otherwise (e.g. actual exon skipping events)
						// Also, the event must not link to the same position in one of the shared exon groups
						
						ExonGroup grpCommonLeft		= null;
						ExonGroup grpCommonRight	= null;
						
						for(ExonGroup ex2 : mapJunctionsToExons.keySet())
						{
							if(jun1.m_nStart <= ex2.getGenomicStopOfGroup() && jun1.m_nStart >= ex2.getGenomicStartOfGroup())
							{
								if(jun2.m_nStart <= ex2.getGenomicStopOfGroup() && jun2.m_nStart >= ex2.getGenomicStartOfGroup())
								{
									if(jun1.m_nStart != jun2.m_nStart)
										grpCommonLeft = ex2;
								}
							}

							if(jun1.m_nEnd <= ex2.getGenomicStopOfGroup() && jun1.m_nEnd >= ex2.getGenomicStartOfGroup())
							{
								if(jun2.m_nEnd <= ex2.getGenomicStopOfGroup() && jun2.m_nEnd >= ex2.getGenomicStartOfGroup())
								{
									if(jun1.m_nEnd != jun2.m_nEnd)
										grpCommonRight = ex2;
								}
							}
						}
						
						// skip events that do not have any common group (in which the junctions map to different positions)
						if(grpCommonRight == null && grpCommonLeft == null)
							continue;
						
						// test all condition combinations
						for(int nCondA=0; nCondA<nConditions; nCondA++)
						{
							for(int nCondB=nCondA+1; nCondB<nConditions; nCondB++)
							{
								String strConditionA = pConditions[nCondA];
								String strConditionB = pConditions[nCondB];
								
								double pRatios[][] = data.RetrieveJunctionRatios(app, strConditionA, strConditionB, jun1, jun2);
								double pRatiosA[] = pRatios[0];
								double pRatiosB[] = pRatios[1];

								// perform t-test (only if there are at least 3 samples in each group)
								if(pRatiosA.length >= 3 && pRatiosB.length >= 3)
								{									
									double fPValue = Double.NaN;
									TTest test = new TTest();
									
									if(projectModel.ConditionsHavePairedData(strSelectedConditionType, strConditionA, strConditionB))
									{
										// sort data by individual id
										fPValue = test.pairedTTest(pRatiosA, pRatiosB);
									}
									else
									{
										fPValue = test.tTest(pRatiosA, pRatiosB);
									}
									
									
									/*
									int pIndices[] = new int[pRatiosA.length + pRatiosB.length];
									double pValues[] = new double[pRatiosA.length + pRatiosB.length];
									
									int nIdx = 0;
									for(int z=0; z<pRatiosA.length; z++)
									{
										pValues[nIdx] = pRatiosA[z];
										pIndices[nIdx] = 1;
										nIdx++;
									}
									for(int z=0; z<pRatiosB.length; z++)
									{
										pValues[nIdx] = pRatiosB[z];
										pIndices[nIdx] = 2;
										nIdx++;
									}
									
									double pRes[] = DistributionTest.kruskal_wallis_test(pValues, pIndices);
									double fPValue = pRes[1];
									*/
									
									if(!Double.isNaN(fPValue))
									{
										vcPValues.add(fPValue);
									
										if(fPValue <= 0.05)
										{
											double fIncLevelChange = Math.abs(StatUtils.mean(pRatiosA) - StatUtils.mean(pRatiosB));
		
											if(k == 0)
											{
												SimpleSpliceScore score = new SimpleSpliceScore(data.GetGeneID(), jun1, jun2, fPValue, fIncLevelChange, strConditionA, strConditionB, SplicingWebApp.AS_TYPE_ALT_5_PRIME_EXON_END);
												score.CheckIfNovel(data.GetGene());
												results.add(score);
											}
											else
											{
												SimpleSpliceScore score = new SimpleSpliceScore(data.GetGeneID(), jun1, jun2, fPValue, fIncLevelChange, strConditionA, strConditionB, SplicingWebApp.AS_TYPE_ALT_3_PRIME_EXON_END);
												score.CheckIfNovel(data.GetGene());
												results.add(score);
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		
		results = app.AdjustPValuesPsiScores(results, vcPValues);

		return results;
	}
}
