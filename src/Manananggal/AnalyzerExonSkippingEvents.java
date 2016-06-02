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

public class AnalyzerExonSkippingEvents
{
	public AlternativeSplicingHit IdentifyVariableExons(SplicingWebApp app, boolean bIgnoreFirstAndLastExons) throws IOException
	{
		//###################################################################
		//   retrieve data for the current selections and gene information
		//###################################################################
		TreeSet<String> vcValidIsoforms 	= app.GetValidIsoforms();
		String strSelectedConditionType 	= app.GetSelectedConditionType();
		TreeSet<String> vcSelectedSamples 	= app.GetSelectedSamples();
		int nMinCovPerBase 					= app.GetMinimumCoveragePerBase();
		int nMinJunctionReads 				= app.GetMinimumJunctionReads();
		double fMinCoveredBases				= app.GetMinimumCoveredBases();
		double fVariableExonThreshold 		= app.GetVariableExonThreshold();
		String strFileGTF					= app.GetGTFFile();
		ProjectModel projectModel			= app.GetProjectModel();
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = projectModel.GetSelectedSamplesPerCondition(strSelectedConditionType, vcSelectedSamples);
//		TreeMap<String, TreeSet<String>> mapSamplesToIndividuals = projectModel.GetSamplesPerIndividual();
		DataSupplier data					= app.GetDataSupplier();
		//###################################################################
		
		//###################################################################
		//             retrieve data from the data supplier
		//###################################################################
		TreeSet<String> vcInvalidConditions = data.GetInvalidConditions();
		//###################################################################

		// number of conditions
		int nConditions = mapSamplesToConditions.size();

		//###################################################################
		//                   get gene identifier
		//###################################################################
		GeneIdentifier gid = data.GetGeneIdentifier();
		
		//###################################################################
		//                    prepare result container
		//###################################################################
		AlternativeSplicingHit result = null;
		result = new AlternativeSplicingHit(0, gid, data.GetReferenceName(), nMinJunctionReads, nMinCovPerBase,
				fMinCoveredBases, fVariableExonThreshold, vcValidIsoforms, strFileGTF, "");

		//###########################################
		//    test for all condition combinations
		//###########################################
		String pConditions[] = new String[nConditions];
		mapSamplesToConditions.keySet().toArray(pConditions);
		
		ExonGroup pExonGroups[] = null;
		
		if(bIgnoreFirstAndLastExons)
		{
			pExonGroups = data.GetMiddleExonGroups();
		}
		else
		{
			pExonGroups = data.GetExonGroups();
		}
		
		for(int i=0; i<nConditions; i++)
		{
			// skip invalid conditions
			if(vcInvalidConditions.contains(pConditions[i]))
				continue;

			for(int j=i+1; j<nConditions; j++)
			{
				double fMeanExonCoverageConditionA = 0.0;
				double fMeanExonCoverageConditionB = 0.0;
				
				// skip invalid conditions
				if(vcInvalidConditions.contains(pConditions[j]))
					continue;
				
				String strConditionA = pConditions[i];
				String strConditionB = pConditions[j];
				
				// get samples for each condition
				TreeSet<String> vcSamplesA = mapSamplesToConditions.get(strConditionA);
				TreeSet<String> vcSamplesB = mapSamplesToConditions.get(strConditionB);
		
				// prepare fraction map for exon groups and exons
				TreeMap<ExonGroup, double[]> mapCovFractionsPerExonGroup = new TreeMap<ExonGroup, double[]>();
				TreeMap<Exon, double[]> 	 mapCovFractionsPerExon		 = new TreeMap<Exon, double[]>();
				
				// calculate fractions for the selected conditions
				for(ExonGroup grp : pExonGroups)
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
							if(bIgnoreFirstAndLastExons)
							{
								pCov[nSample] = data.GetCoverageArrayForMiddleExonGroupAndSample(grp, strSample)[x];
							}
							else
							{
								pCov[nSample] = data.GetCoverageArrayForExonGroupAndSample(app, grp, strSample)[x];
							}
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
							if(bIgnoreFirstAndLastExons)
							{
								pCov[nSample] = data.GetCoverageArrayForMiddleExonGroupAndSample(grp, strSample)[x];
							}
							else
							{
								pCov[nSample] = data.GetCoverageArrayForExonGroupAndSample(app, grp, strSample)[x];
							}
							
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
					
					// calculate for exons
					for(Exon ex : grp.getExons())
					{
						Vector<Double> vcExonFractions = new Vector<Double>();
						
						// skip exons that are not 'middle' exons
						if(!data.GetMiddleExons().contains(ex))
							continue;

						// get mean coverage per base for the whole exon
						double pMeanCovA[] = new double[vcSamplesA.size()];
						double pMeanCovB[] = new double[vcSamplesB.size()];
						int nSample = 0;
						for(String strSample : vcSamplesA)
						{
							pMeanCovA[nSample] = StatUtils.mean(data.GetCoverageArrayForExonAndSample(ex, strSample));
							nSample++;
						}
						nSample = 0;
						for(String strSample : vcSamplesB)
						{
							pMeanCovB[nSample] = StatUtils.mean(data.GetCoverageArrayForExonAndSample(ex, strSample));
							nSample++;
						}

						fMeanExonCoverageConditionA = StatUtils.mean(pMeanCovA);
						fMeanExonCoverageConditionB = StatUtils.mean(pMeanCovB);
						
						// ignore result that have insufficient coverage
						if(fMeanExonCoverageConditionA < nMinCovPerBase && fMeanExonCoverageConditionB < nMinCovPerBase)
							continue;
						
						for(int x=0; x<ex.getGenomicLength(); x++)
						{
							//##########################################################
							//   get mean coverage for this exon group in condition A
							//##########################################################
							double pCov[] = new double[vcSamplesA.size()];							
							
							nSample = 0;
							for(String strSample : vcSamplesA)
							{
								pCov[nSample] = data.GetCoverageArrayForExonAndSample(ex, strSample)[x];								
								nSample++;
							}
							
							// get mean coverage at current position
							double fMeanCoverageA = StatUtils.mean(pCov);

							//##########################################################
							//   get mean coverage for this exon group in condition B
							//##########################################################
							pCov = new double[vcSamplesB.size()];
							
							nSample = 0;
							for(String strSample : vcSamplesB)
							{
								pCov[nSample] = data.GetCoverageArrayForExonAndSample(ex, strSample)[x];
								nSample++;
							}
							
							// get mean coverage at current position
							double fMeanCoverageB = StatUtils.mean(pCov);
							
							//############################################
							//   get total coverage and calculate ratio
							//############################################
							double fTotalCoverage   = fMeanCoverageA + fMeanCoverageB;

							// calculate ratio, skip uncovered positions
							if(fTotalCoverage != 0)
								vcExonFractions.add(fMeanCoverageA / fTotalCoverage);
						}
						
						double[] pExonFractions = new double[vcExonFractions.size()];
						for(int k=0; k<vcExonFractions.size(); k++)
							pExonFractions[k] = vcExonFractions.get(k);
					
						mapCovFractionsPerExon.put(ex, pExonFractions);
					}
				}

				// perform test for each exon to identify alternatively spliced exons
				for(Exon ex : mapCovFractionsPerExon.keySet())
				{
					double[] pFractionsGroup1 = mapCovFractionsPerExon.get(ex);
					double fAvgGroup1 = StatUtils.mean(pFractionsGroup1);
	
					// get average coverage ratio for all other exons
					Vector<Double> vcAllGroups = new Vector<Double>();
					for(ExonGroup grp2 : mapCovFractionsPerExonGroup.keySet())
					{					
						if(!grp2.groupContainsExon(ex.getExonID()))
						{
							double pFractionsGroup2[] = mapCovFractionsPerExonGroup.get(grp2);
							vcAllGroups.add(StatUtils.mean(pFractionsGroup2));
						}
					}
					
					double pAllGroups[] = new double[vcAllGroups.size()];
					for(int k=0; k<vcAllGroups.size(); k++)
						pAllGroups[k] = vcAllGroups.get(k);
					
					double fAvgOtherExonGroups = StatUtils.mean(pAllGroups);
					
					double fPValue = 1.0;
					if(pFractionsGroup1.length > 1 && pAllGroups.length > 1) // ignore 1bp exons and 2-exon transcripts
					{
						TTest test = new TTest();
						fPValue = test.tTest(pFractionsGroup1, pAllGroups);
					}

					double fRelDifference = Math.abs(1.0 - fAvgGroup1 / fAvgOtherExonGroups);
					double fAbsDifference = fAvgGroup1 - fAvgOtherExonGroups;
					/*
					if(ex.getCodingStart() == 203702351 && ex.getCodingStop() == 203702528 && strConditionA.equals("GBM"))
					{
						System.out.println("testing: " + strConditionA + " vs " + strConditionB);
						System.out.println(Arrays.toString(pFractionsGroup1));
						System.out.println(Arrays.toString(pAllGroups));						
						System.out.println("p-value:" + fPValue);
						System.out.println("abs. change: " + fAbsDifference);
						System.out.println("threshold: " + m_fVariableExonThreshold);
					}
					*/
					
					if(fPValue < 0.05 && Math.abs(fAbsDifference) >= fVariableExonThreshold)
					{
						/*
						if(ex.getCodingStart() == 203702351 && ex.getCodingStop() == 203702528)
						{
							System.out.println("adding result for: " + strConditionA + " vs " + strConditionB);
						}
						*/
						
						result.AddASExons(SplicingWebApp.AS_TYPE_EXON_SKIPPING, data.GetReferenceName() + ":" + ex.getCodingStart() + "-" + ex.getCodingStop(), ex.getCodingStart(), ex.getCodingStop(),
								strConditionA, strConditionB, fMeanExonCoverageConditionA, fMeanExonCoverageConditionB, fAbsDifference, fRelDifference, fPValue, pFractionsGroup1, pAllGroups);
					}
				}
			}
		}
		
		return result;
	}
	
	public TreeSet<SimpleSpliceScore> CalculateJunctionPSIScores(SplicingWebApp app, boolean bDebug, TreeSet<SimpleSpliceScore> results) throws IOException
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
		
		// get number of conditions
		int nConditions = mapSamplesToConditions.keySet().size();
		
		// list of all p-values
		Vector<Double> vcPValues = new Vector<Double>();
		
		// get a list of all junctions that are included in the junction count file
		TreeSet<CountElement> vcJunctionsWithCoverage = data.GetJunctions();
		
		//#####################################################################
		//        get all inclusion/exclusion junction combinations
		//#####################################################################
		// key = inclusion junction, value = exclusion junctions
		TreeMap<CountElement, TreeSet<CountElement>> mapCombinations = new TreeMap<CountElement, TreeSet<CountElement>>();
		
		// add all novel junctions to the valid junctions
		TreeSet<CountElement> vcValidJunctions = new TreeSet<CountElement>();
		for(CountElement e : vcJunctionsWithCoverage)
		{
			if(!e.m_bKnown)
				vcValidJunctions.add(e);
		}
		
		TreeSet<Exon> vcValidExons = new TreeSet<Exon>();
		for(String strIsoform : vcValidIsoforms)
		{
			Exon pIsoformExons[] = data.GetExonsForIsoform(strIsoform);
			for(Exon ex : pIsoformExons)
				vcValidExons.add(ex);
			
			// get splice junctions
			TreeSet<CountElement> vcJunctions =	data.GetJunctionsForIsoform(strIsoform);
			for(CountElement e : vcJunctions)
			{
				if(vcJunctionsWithCoverage.contains(e))
				{
					vcValidJunctions.add(e);
				}
			}
		}
		
		Exon pExons[] = new Exon[vcValidExons.size()];
		vcValidExons.toArray(pExons);

		for(Exon ex : pExons)
		{
			TreeSet<CountElement> vcInclusionJunctions = new TreeSet<CountElement>();
			TreeSet<CountElement> vcExclusionJunctions = new TreeSet<CountElement>();
			
			for(CountElement jun : vcValidJunctions)
			{
				// check whether the junction is valid (passes the threshold filter)
				if(data.IsInvalidJunction(jun))
					continue;

				// inclusion junction
				if(jun.m_nStart == ex.getCodingStop() || jun.m_nEnd == ex.getCodingStart())
				{
					vcInclusionJunctions.add(jun);
				}
				// exclusion junction
				else if(jun.m_nStart < ex.getCodingStart() && jun.m_nEnd > ex.getCodingStop())
				{
					vcExclusionJunctions.add(jun);
				}
			}
			
			// skip exons that only have one inclusion or exclusion junction
			if(vcInclusionJunctions.size() < 1 || vcExclusionJunctions.size() < 1)
				continue;
			
			for(CountElement inclJun : vcInclusionJunctions)
			{
				if(mapCombinations.containsKey(inclJun))
				{					
					for(CountElement exclJun : vcExclusionJunctions)
					{
						mapCombinations.get(inclJun).add(exclJun);
					}
				}
				else
				{
					TreeSet<CountElement> vcTmp = new TreeSet<CountElement>();
					for(CountElement exclJun : vcExclusionJunctions)
					{
						vcTmp.add(exclJun);
					}
					mapCombinations.put(inclJun, vcTmp);
				}
			}			
		}
		
		// Get splice scores
		TreeSet<SimpleSpliceScore> vcSpliceScores = new TreeSet<SimpleSpliceScore>();

		for(CountElement junIncl : mapCombinations.keySet())
		{
			for(CountElement junExcl : mapCombinations.get(junIncl))
			{
				// test all combinations
				String pConditions[] = new String[nConditions];
				mapSamplesToConditions.keySet().toArray(pConditions);
		
				for(int i=0; i<nConditions; i++)
				{
					for(int j=i+1; j<nConditions; j++)
					{					
						String strConditionA = pConditions[i];
						String strConditionB = pConditions[j];

						double pRatios[][] = data.RetrieveJunctionRatios(app, strConditionA, strConditionB, junIncl, junExcl);
						double pRatiosA[] = pRatios[0];
						double pRatiosB[] = pRatios[1];

						// there must be at least two values for each group
						if(pRatiosA.length >= 3 && pRatiosB.length >= 3)
						{
							TTest test = new TTest();
							double fPValue = Double.NaN;
							if(projectModel.ConditionsHavePairedData(strSelectedConditionType, strConditionA, strConditionB))
							{
								// sort data by individual id
								fPValue = test.pairedTTest(pRatiosA, pRatiosB);
							}
							else
							{
								fPValue = test.tTest(pRatiosA, pRatiosB);
							}
							
							double fIncLevel = Math.abs(StatUtils.mean(pRatiosA) - StatUtils.mean(pRatiosB));

							if(!Double.isNaN(fPValue))
							{
								vcPValues.add(fPValue);
							
								if(fPValue <= 0.05)
								{
									SimpleSpliceScore res = new SimpleSpliceScore(data.GetGeneID(), junIncl, junExcl, fPValue, fIncLevel, strConditionA, strConditionB, SplicingWebApp.AS_TYPE_EXON_SKIPPING);
									res.GetValidIsoforms(data.GetGene());
									vcSpliceScores.add(res);
								}
							}
						}
					}
				}
			}
		}
		
		// adjust p-values
		TreeSet<SimpleSpliceScore> vcAdjustedScores = new TreeSet<SimpleSpliceScore>();
		vcAdjustedScores = app.AdjustPValuesPsiScores(vcSpliceScores, vcPValues);
		
		// map results to inclusion junctions
		TreeMap<CountElement, TreeSet<SimpleSpliceScore>> mapResToInclJun = new TreeMap<CountElement, TreeSet<SimpleSpliceScore>>();
		
		for(SimpleSpliceScore score : vcAdjustedScores)
		{
			if(mapResToInclJun.containsKey(score.m_JunctionInclusion))
			{
				mapResToInclJun.get(score.m_JunctionInclusion).add(score);
			}
			else
			{
				TreeSet<SimpleSpliceScore> vcRes = new TreeSet<SimpleSpliceScore>();
				vcRes.add(score);
				mapResToInclJun.put(score.m_JunctionInclusion, vcRes);
			}
		}
		
		// keep only two PSI scores per inclusion junction
		// one for 'known' junctions and another for 'novel' junctions
		// only keep the best p-value result for each
		for(CountElement junIncl : mapResToInclJun.keySet())
		{
			// find best p-value hits
			double fPValueKnown = Double.MAX_VALUE;
			double fPValueNovel = Double.MAX_VALUE; 
			SimpleSpliceScore bestScoreNovel = null;
			SimpleSpliceScore bestScoreKnown = null;
			
			for(SimpleSpliceScore score : mapResToInclJun.get(junIncl))
			{
				if(!score.m_bIsNovel)
				{
					if(score.m_fPValue < fPValueKnown)
					{
						fPValueKnown = score.m_fPValue;
						bestScoreKnown = score;
					}
				}
				else
				{
					if(score.m_fPValue < fPValueNovel)
					{
						fPValueNovel = score.m_fPValue;
						bestScoreNovel = score;
					}
				}
			}
			
			if(bestScoreKnown != null)
				results.add(bestScoreKnown);
			
			if(bestScoreNovel != null && (bestScoreKnown == null || (bestScoreKnown != null && bestScoreKnown.m_fPValue > bestScoreNovel.m_fPValue)))
				results.add(bestScoreNovel);
		}
	
		return results;
	}
}
