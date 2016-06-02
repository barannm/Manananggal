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

public class AnalyzerAltStartOrEndEvents
{
	public AlternativeSplicingHit IdentifyVariableStartOrEndExons(SplicingWebApp app, AlternativeSplicingHit result) throws IOException
	{	
		//###################################################################
		//   retrieve data for the current selections and gene information
		//###################################################################
		TreeSet<String> vcValidIsoforms 	= app.GetValidIsoforms();		
		String strSelectedConditionType 	= app.GetSelectedConditionType();
		TreeSet<String> vcSelectedSamples 	= app.GetSelectedSamples();
		double fVariableExonThreshold 		= app.GetVariableExonThreshold();		
		ProjectModel projectModel			= app.GetProjectModel();
		DataSupplier data					= app.GetDataSupplier();		
		//###################################################################
		
		//###################################################################
		//             retrieve data from the data supplier
		//###################################################################
		TreeSet<String> vcInvalidConditions = data.GetInvalidConditions();
		//###################################################################
		
		// get samples per condition
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = projectModel.GetSelectedSamplesPerCondition(strSelectedConditionType, vcSelectedSamples);
		int nConditions = mapSamplesToConditions.size();
		
		// get conditions
		String pConditions[] = new String[nConditions];
		mapSamplesToConditions.keySet().toArray(pConditions);
		
		TreeSet<Exon> vcCurrentFirstExons = new TreeSet<Exon>();
		TreeSet<Exon> vcCurrentLastExons = new TreeSet<Exon>();
		
		for(String strIsoform : data.GetIsoformNames())
		{
			if(vcValidIsoforms.contains(strIsoform))
			{
				Exon[] pIsoformExons = data.GetExonsForIsoform(strIsoform);
				
				if(data.GetStrand() == '+')
				{
					vcCurrentFirstExons.add(pIsoformExons[0]);
					vcCurrentLastExons.add(pIsoformExons[pIsoformExons.length-1]);
				}
				else
				{
					vcCurrentFirstExons.add(pIsoformExons[pIsoformExons.length-1]);
					vcCurrentLastExons.add(pIsoformExons[0]);
				}
			}
		}

		// determine which exons share their junction with a 'middle' exon
		TreeSet<Exon> vcAltStartEndWithSharedJunction = new TreeSet<Exon>();
		for(String strIsoform : data.GetIsoformNames())
		{
			if(vcValidIsoforms.contains(strIsoform))
			{
				Exon[] pIsoformExons = data.GetExonsForIsoform(strIsoform);
				
				// proceed through list of first exons
				for(Exon ex1 : vcCurrentFirstExons)
				{
					// check for any middle exon that matches the "first exon" end position
					for(int i=1; i<pIsoformExons.length-1; i++)
					{
						Exon ex2 = pIsoformExons[i];
						
						if(data.GetStrand() == '+')
						{
							if(ex1.getCodingStop() == ex2.getCodingStop())
								vcAltStartEndWithSharedJunction.add(ex1);
						}
						else
						{
							if(ex1.getCodingStart() == ex2.getCodingStart())
								vcAltStartEndWithSharedJunction.add(ex1);
						}
					}
				}
				
				// proceed through list of last exons
				for(Exon ex1 : vcCurrentLastExons)
				{
					// check for any middle exon that matches the "last exon" start position
					for(int i=1; i<pIsoformExons.length-1; i++)
					{
						Exon ex2 = pIsoformExons[i];
						
						if(data.GetStrand() == '+')
						{
							if(ex1.getCodingStart() == ex2.getCodingStart())
								vcAltStartEndWithSharedJunction.add(ex1);
						}
						else
						{
							if(ex1.getCodingStop() == ex2.getCodingStop())
								vcAltStartEndWithSharedJunction.add(ex1);
						}
					}
				}
			}
		}

		//###########################################
		//    test for all condition combinations
		//###########################################
		for(int i=0; i<nConditions; i++)
		{
			if(vcInvalidConditions.contains(pConditions[i]))
				continue;
			
			for(int j=i+1; j<nConditions; j++)
			{
				if(vcInvalidConditions.contains(pConditions[j]))
					continue;

				String strConditionA = pConditions[i];
				String strConditionB = pConditions[j];
				
				// get samples for each condition
				TreeSet<String> vcSamplesA = mapSamplesToConditions.get(strConditionA);
				TreeSet<String> vcSamplesB = mapSamplesToConditions.get(strConditionB);
		
				// prepare fraction map for exon groups and exons
				TreeMap<ExonGroup, double[]> 	mapCovFractionsPerExonGroup = new TreeMap<ExonGroup, double[]>();
				
				// calculate fractions for the selected conditions
				for(ExonGroup grp : data.GetMiddleExonGroups())
				{
					Vector<Double> vcFractions = new Vector<Double>();

					// calculate for exon groups
					for(int x=0; x<grp.getExonGroupLengthInBp(); x++)
					{
						//##########################################################
						//   get mean coverage for this exon group in condition A
						//##########################################################
						double pCov[] = new double[vcSamplesA.size()];
						
						int nSample = 0;
						for(String strSample : vcSamplesA)
						{
							pCov[nSample] = data.GetCoverageArrayForMiddleExonGroupAndSample(grp, strSample)[x];
							nSample++;
						}
						double fMeanCoverageA = StatUtils.mean(pCov);

						//##########################################################
						//   get mean coverage for this exon group in condition B
						//##########################################################
						pCov = new double[vcSamplesB.size()];
						
						nSample = 0;
						for(String strSample : vcSamplesB)
						{
							pCov[nSample] = data.GetCoverageArrayForMiddleExonGroupAndSample(grp, strSample)[x];
							nSample++;
						}
						
						double fMeanCoverageB = StatUtils.mean(pCov);

						//############################################
						//   get total coverage and calculate ratio
						//############################################
						double fTotalCoverage   = fMeanCoverageA + fMeanCoverageB;

						// calculate ratio, skip uncovered positions
						if(fTotalCoverage != 0)
							vcFractions.add(fMeanCoverageA / fTotalCoverage);
					}
					
					double[] pFractions = new double[vcFractions.size()];
					for(int k=0; k<vcFractions.size(); k++)
						pFractions[k] = vcFractions.get(k);
					
					mapCovFractionsPerExonGroup.put(grp, pFractions);
				}
				
				// get average coverage ratio for all 'middle' exon groups
				Vector<Double> vcAllGroups = new Vector<Double>();
				for(ExonGroup grp2 : mapCovFractionsPerExonGroup.keySet())
				{					
					double pFractionsGroup2[] = mapCovFractionsPerExonGroup.get(grp2);
					vcAllGroups.add(StatUtils.mean(pFractionsGroup2));
				}
				
				double pAllGroups[] = new double[vcAllGroups.size()];
				for(int l=0; l<vcAllGroups.size(); l++)
					pAllGroups[l] = vcAllGroups.get(l);
				
				double fAvgOtherExonGroups = StatUtils.mean(pAllGroups);

				for(int k=0; k<2; k++)
				{					
					TreeMap<Exon, double[]> mapCovFractionsPerExon = new TreeMap<Exon, double[]>();
					
					TreeSet<Exon> vcExons = null;
					if(k==0)
						vcExons = vcCurrentFirstExons;
					else
						vcExons = vcCurrentLastExons;

					TreeSet<AlternativeSplicingExon> vcHits = new TreeSet<AlternativeSplicingExon>();
					
					// calculate for exons
					for(Exon ex : vcExons)
					{
						double[] pFractions = new double[ex.getGenomicLength()];
						
						for(int x=0; x<ex.getGenomicLength(); x++)
						{
							// get mean coverage for this exon in condition A
							double pCov[] = new double[vcSamplesA.size()];
							
							int nSample = 0;
							for(String strSample : vcSamplesA)
							{
								pCov[nSample] = data.GetCoverageArrayForExonAndSample(ex, strSample)[x];
								nSample++;
							}

							double fMeanCoverageA = StatUtils.mean(pCov);
							
							// get mean coverage for this exon in condition B						
							pCov = new double[vcSamplesB.size()];
							
							nSample = 0;
							for(String strSample : vcSamplesB)
							{ 
								pCov[nSample] = data.GetCoverageArrayForExonAndSample(ex, strSample)[x];
								nSample++;
							}
							double fMeanCoverageB = StatUtils.mean(pCov);

							// get total coverage
							double fTotalCoverage   = fMeanCoverageA + fMeanCoverageB;

							if(fTotalCoverage == 0)
								continue;
							
							// calculate ratio
							pFractions[x] = fMeanCoverageA / fTotalCoverage;
						}

						mapCovFractionsPerExon.put(ex, pFractions);
					}
					
					// perform test for each exon group to identify alternatively spliced exons
					for(Exon ex : mapCovFractionsPerExon.keySet())
					{
						int nType = -1;

						// determine whether it is a alternative first or last exon						
						nType = SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN;
						
						// check whether the junction is unique to the alternative start or end exon
						// (otherweise it is an alternative start/end that shares its junction with a 'middle' exon)
						if(vcAltStartEndWithSharedJunction.contains(ex))
							nType = SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN;
						
						if(k == 1)
						{
							nType = SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN;
							
							if(vcAltStartEndWithSharedJunction.contains(ex))
								nType = SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN;
						}

						double[] pFractionsGroup1 = mapCovFractionsPerExon.get(ex);
						double fAvgGroup1 = StatUtils.mean(pFractionsGroup1);
						
						double fPValue = 1.0;
						if(pFractionsGroup1.length > 1 && pAllGroups.length > 1) // ignore 1bp exons and 2-exon transcripts
						{
							TTest test = new TTest();
							fPValue = test.tTest(pFractionsGroup1, pAllGroups);
						}
	
						double fRelDifference = Math.abs(1.0 - fAvgGroup1 / fAvgOtherExonGroups);
						double fAbsDifference = fAvgGroup1 - fAvgOtherExonGroups;
						
						if(fPValue < 0.05 && Math.abs(fAbsDifference) >= fVariableExonThreshold)
						{
							double fMeanCoveragePerBaseA = 0.0;
							double fMeanCoveragePerBaseB = 0.0;
									
							// get mean coverage per base for condition A
							double pMeanCovA[] = new double[vcSamplesA.size()];
							int nSample=0;
							for(String strSample : vcSamplesA)
							{
								double pCov[] = data.GetCoverageArrayForExonAndSample(ex, strSample);
								pMeanCovA[nSample] = StatUtils.mean(pCov);
								nSample++;
							}
							fMeanCoveragePerBaseA = StatUtils.mean(pMeanCovA);
							
							// get coverage array for samples in condition B
							double pMeanCovB[] = new double[vcSamplesB.size()];
							nSample=0;
							for(String strSample : vcSamplesB)
							{
								double pCov[] = data.GetCoverageArrayForExonAndSample(ex, strSample);
								pMeanCovB[nSample] = StatUtils.mean(pCov);
								nSample++;
							}
							fMeanCoveragePerBaseB = StatUtils.mean(pMeanCovB);
							
							AlternativeSplicingExon hit = new AlternativeSplicingExon(nType, data.GetReferenceName() + ":" + ex.getCodingStart() + "-" + ex.getCodingStop(), ex.getCodingStart(), ex.getCodingStop(),
									strConditionA, strConditionB, fMeanCoveragePerBaseA, fMeanCoveragePerBaseB, fAbsDifference, fRelDifference, fPValue, pFractionsGroup1, pAllGroups);
							
							vcHits.add(hit);
							/*
							result.AddASExons(nType, data.GetReferenceName() + ":" + ex.getCodingStart() + "-" + ex.getCodingStop(), ex.getCodingStart(), ex.getCodingStop(),
									strConditionA, strConditionB, fMeanCoveragePerBaseA, fMeanCoveragePerBaseB, fAbsDifference, fRelDifference, fPValue, pFractionsGroup1, pAllGroups);
							*/
						}
					}
					
					// map hits to positions
					TreeMap<Integer, TreeSet<AlternativeSplicingExon>> mapHitsToSplicePosition = new TreeMap<Integer, TreeSet<AlternativeSplicingExon>>();
					
					if(k == 0)
					{
						// k == 0 -> alternative start exons -> splice donors
						for(AlternativeSplicingExon hit : vcHits)
						{
							int nPosition = 0;
							
							if(data.GetStrand() == '+')
								nPosition = hit.GetEnd();
							else
								nPosition = hit.GetStart();
							
							if(mapHitsToSplicePosition.containsKey(nPosition))
							{
								mapHitsToSplicePosition.get(nPosition).add(hit);
							}
							else
							{
								TreeSet<AlternativeSplicingExon> vcTmp = new TreeSet<AlternativeSplicingExon>();
								vcTmp.add(hit);
								mapHitsToSplicePosition.put(nPosition, vcTmp);
							}
						}
					}
					else
					{
						// k == 0 -> alternative end exons -> splice acceptors						
						for(AlternativeSplicingExon hit : vcHits)
						{
							int nPosition = 0;
							
							if(data.GetStrand() == '+')
								nPosition = hit.GetStart();
							else
								nPosition = hit.GetEnd();
							
							if(mapHitsToSplicePosition.containsKey(nPosition))
							{
								mapHitsToSplicePosition.get(nPosition).add(hit);
							}
							else
							{
								TreeSet<AlternativeSplicingExon> vcTmp = new TreeSet<AlternativeSplicingExon>();
								vcTmp.add(hit);
								mapHitsToSplicePosition.put(nPosition, vcTmp);
							}
						}
					}
					
					// keep the best result for each alternative splice acceptor or donor site
					for(int nPosition : mapHitsToSplicePosition.keySet())
					{
						AlternativeSplicingExon bestHit = null;
						
						for(AlternativeSplicingExon hit : mapHitsToSplicePosition.get(nPosition))
						{
							if(bestHit == null || hit.m_fPValue < bestHit.m_fPValue)
							{
								bestHit = hit;
							}
						}
						
						result.AddASExon(bestHit);
					}
				}
			}
		}
		
		return result;
	}

	public TreeSet<SimpleSpliceScore> CalculateJunctionInclusionScoresForFirstAndLastExons(SplicingWebApp app, boolean bDebug, TreeSet<SimpleSpliceScore> results) throws IOException
	{	
		//###################################################################
		//   retrieve data for the current selections and gene information
		//###################################################################
		TreeSet<String> vcValidIsoforms 	= app.GetValidIsoforms();
		String strSelectedConditionType 	= app.GetSelectedConditionType();
		int nMinJunctionReads 				= app.GetMinimumJunctionReads();
		TreeSet<String> vcSelectedSamples 	= app.GetSelectedSamples();		
		ProjectModel projectModel			= app.GetProjectModel();
		DataSupplier data					= app.GetDataSupplier();
		Exon pExons[]						= data.GetExons();
		//###################################################################
				
		// create new result container if necessary
		if(results == null)
			results = new TreeSet<SimpleSpliceScore>();
		
		// get samples per Condition
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = projectModel.GetSelectedSamplesPerCondition(strSelectedConditionType, vcSelectedSamples);
		
		// get number of conditions
		int nConditions = mapSamplesToConditions.keySet().size();
		
		//################################
		//    get all first/last exons
		//################################
		TreeSet<Exon> vcFirstExons = new TreeSet<Exon>();
		TreeSet<Exon> vcLastExons  = new TreeSet<Exon>();
		
		for(String strIsoform : vcValidIsoforms)
		{
			if(vcValidIsoforms.contains(strIsoform))
			{
				Exon pIsoformExons[] = data.GetExonsForIsoform(strIsoform);

				if(data.GetStrand() == '+')
				{
					vcFirstExons.add(pIsoformExons[0]);
					vcLastExons.add(pIsoformExons[pIsoformExons.length-1]);
				}
				else
				{
					vcFirstExons.add(pExons[pIsoformExons.length-1]);
					vcLastExons.add(pExons[0]);
				}
			}
		}
		
		//##################################################################
		//    filter junctions that are connected to first or last exons
		//##################################################################											= new TreeSet<CountElement>();
		TreeMap<CountElement, TreeMap<String, Integer>> mapJunctionCountsFirstExons = new TreeMap<CountElement, TreeMap<String, Integer>>();
		TreeMap<CountElement, TreeMap<String, Integer>> mapJunctionCountsLastExons  = new TreeMap<CountElement, TreeMap<String, Integer>>();
		
		for(String strIsoform : vcValidIsoforms)
		{
			TreeSet<CountElement> vcIsoformJunctions = data.GetJunctionsForIsoform(strIsoform);
			
			// skip single exon genes and genes that do not have any junction counts
			if(vcIsoformJunctions == null || vcIsoformJunctions.isEmpty())
				continue;
			
			CountElement firstJunction	= null;
			CountElement lastJunction	= null;
			
			if(data.GetStrand() == '+')
			{
				firstJunction = vcIsoformJunctions.first();
				lastJunction  = vcIsoformJunctions.last();
			}
			else
			{
				firstJunction = vcIsoformJunctions.last();
				lastJunction  = vcIsoformJunctions.first();
			}
			
			// check whether the junction passes the threshold filter
			boolean bFirstJunIsValid = false;
			boolean bLastJunIsValid = false;
			for(String strCondition : mapSamplesToConditions.keySet())
			{
				double fSumFirst = 0.0f;
				double fSumLast = 0.0f;
				int nSamples = 0;
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					TreeMap<String, Integer> mapCountsToSamples = data.GetCountsForJunction(firstJunction);					
					
					if(mapCountsToSamples != null && mapCountsToSamples.containsKey(strSample))
						fSumFirst += mapCountsToSamples.get(strSample);
					
					mapCountsToSamples = data.GetCountsForJunction(lastJunction);
					
					if(mapCountsToSamples != null && mapCountsToSamples.containsKey(strSample))
						fSumLast += mapCountsToSamples.get(strSample);
					
					nSamples++;
				}
				
				double fMeanFirst = fSumFirst/nSamples;
				double fMeanLast = fSumLast/nSamples;
				
				if(fMeanFirst >= nMinJunctionReads)
					bFirstJunIsValid = true;
				
				if(fMeanLast >= nMinJunctionReads)
					bLastJunIsValid = true;
			}
			
			if(data.GetCountsForJunction(firstJunction) != null && bFirstJunIsValid)
				mapJunctionCountsFirstExons.put(firstJunction, data.GetCountsForJunction(firstJunction));
			
			if(data.GetCountsForJunction(lastJunction) != null && bLastJunIsValid)
				mapJunctionCountsLastExons.put(lastJunction,   data.GetCountsForJunction(lastJunction));
		}

		String pConditions[] = new String[nConditions];
		mapSamplesToConditions.keySet().toArray(pConditions);
		
		// prepare list of p-values
		Vector<Double> vcPValues = new Vector<Double>();
		
		CountElement[] pJunctionsFirst = new CountElement[mapJunctionCountsFirstExons.keySet().size()];
		mapJunctionCountsFirstExons.keySet().toArray(pJunctionsFirst);
		
		// compare all combinations of first-exon junctions
		for(int i=0; i<mapJunctionCountsFirstExons.size()-1; i++)
		{
			for(int j=i+1; j<mapJunctionCountsFirstExons.size(); j++)
			{
				CountElement jun1 = pJunctionsFirst[i];
				CountElement jun2 = pJunctionsFirst[j];

				if(data.GetStrand() == '+')
				{
					if(jun1.m_nStart == jun2.m_nStart)
						continue;
				}
				else
				{
					if(jun1.m_nEnd == jun2.m_nEnd)
						continue;
				}
				
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
							
							if(!Double.isNaN(fPValue))
								vcPValues.add(fPValue);
							
							if(fPValue <= 0.05)
							{
								double fIncLevelChange = Math.abs(StatUtils.mean(pRatiosA) - StatUtils.mean(pRatiosB));

								if(data.GetStrand() == '+')
								{
									SimpleSpliceScore score = new SimpleSpliceScore(data.GetGeneID(), jun1, jun2, fPValue, fIncLevelChange, strConditionA, strConditionB, SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN);
									score.GetValidIsoforms(data.GetGene());
									results.add(score);
								}
								else
								{
									SimpleSpliceScore score = new SimpleSpliceScore(data.GetGeneID(), jun1, jun2, fPValue, fIncLevelChange, strConditionA, strConditionB, SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN);
									score.GetValidIsoforms(data.GetGene());
									results.add(score);
								}
							}
						}
					}
				}
			}
		}
		
		CountElement[] pJunctionsLast = new CountElement[mapJunctionCountsLastExons.keySet().size()];
		mapJunctionCountsLastExons.keySet().toArray(pJunctionsLast);
		
		// compare all combinations of last-exon junctions
		for(int i=0; i<mapJunctionCountsLastExons.size()-1; i++)
		{
			for(int j=i+1; j<mapJunctionCountsLastExons.size(); j++)
			{
				CountElement jun1 = pJunctionsLast[i];
				CountElement jun2 = pJunctionsLast[j];
	
				if(data.GetStrand() == '+')
				{
					if(jun1.m_nEnd == jun2.m_nEnd)
						continue;
				}
				else
				{
					if(jun1.m_nStart == jun2.m_nStart)
						continue;
				}

				for(int nCondA=0; nCondA<nConditions; nCondA++)
				{
					for(int nCondB=nCondA+1; nCondB<nConditions; nCondB++)
					{						
						String strConditionA = pConditions[nCondA];
						String strConditionB = pConditions[nCondB];
						
						double pRatios[][] = data.RetrieveJunctionRatios(app, strConditionA, strConditionB, jun1, jun2);
						double pRatiosA[] = pRatios[0];
						double pRatiosB[] = pRatios[1];

						// there must be at least two values for each group
						if(pRatiosA.length > 3 && pRatiosB.length > 3)
						{
							TTest test = new TTest();
							
							double fPValue = Double.NaN;
							if(projectModel.ConditionsHavePairedData(strSelectedConditionType, strConditionA, strConditionB))
							{
								fPValue = test.pairedTTest(pRatiosA, pRatiosB);
							}
							else
							{
								fPValue = test.tTest(pRatiosA, pRatiosB);
							}

							if(!Double.isNaN(fPValue))
							{
								vcPValues.add(fPValue);

								if(fPValue <= 0.05)
								{
									double fIncLevelChange = Math.abs(StatUtils.mean(pRatiosA) - StatUtils.mean(pRatiosB));
	
									if(data.GetStrand() == '+')
									{
										SimpleSpliceScore score = new SimpleSpliceScore(data.GetGeneID(), jun1, jun2, fPValue, fIncLevelChange, strConditionA, strConditionB, SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN);
										score.GetValidIsoforms(data.GetGene());
										results.add(score);
									}
									else
									{
										SimpleSpliceScore score = new SimpleSpliceScore(data.GetGeneID(), jun1, jun2, fPValue, fIncLevelChange, strConditionA, strConditionB, SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN);
										score.GetValidIsoforms(data.GetGene());
										results.add(score);
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

