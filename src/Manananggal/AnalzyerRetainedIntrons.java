package Manananggal;

import java.io.IOException;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.TTest;

import BioKit.Exon;
import BioKit.ExonGroup;

public class AnalzyerRetainedIntrons {

	public AlternativeSplicingHit IdentifyIntronRetentionEvents(SplicingWebApp app, AlternativeSplicingHit result, boolean bIgnoreFirstAndLastExons) throws IOException
	{
		//###################################################################
		//   retrieve data for the current selections and gene information
		//###################################################################
		String strSelectedConditionType 	= app.GetSelectedConditionType();
		TreeSet<String> vcSelectedSamples 	= app.GetSelectedSamples();
		int nMinCovPerBase 					= app.GetMinimumCoveragePerBase();
		double fMinCoveredBases				= app.GetMinimumCoveredBases();
		double fVariableExonThreshold 		= app.GetVariableExonThreshold();
		ProjectModel projectModel 			= app.GetProjectModel();
		DataSupplier data					= app.GetDataSupplier();
		//###################################################################
		
		//###################################################################
		//             retrieve data from the data supplier
		//###################################################################
		TreeSet<String> vcInvalidConditions = data.GetInvalidConditions();
		//###################################################################
		
		// define the length of the intron fragments that will be investigated, 3 fragments will be analyzed in total: 3'-flank, 5'-flank and middle of the intron
//		int nMaxIntronSectionSize = 200;
		
		// get samples per condition
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = projectModel.GetSelectedSamplesPerCondition(strSelectedConditionType, vcSelectedSamples);
		int nConditions = mapSamplesToConditions.size();
		
		// get conditions
		String pConditions[] = new String[nConditions];
		mapSamplesToConditions.keySet().toArray(pConditions);
		
		//###############################################################
		//        estimate background coverage in introns
		//###############################################################
		TreeMap<String, Double> mapNoiseToCondition = new TreeMap<String, Double>();
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			double pMeanBackgroundPerIntron[] = new double[data.GetIntrons().size()];
			int nIntronIdx = 0;
			
			for(Exon intron : data.GetIntrons())
			{
				double fBackgroundIntronCoverageSum = 0.0;
				int nSamples = 0;
				
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					double pCov[] = data.GetCoverageArrayForIntronAndSample(intron, strSample);
					double fMean = StatUtils.mean(pCov);
					fBackgroundIntronCoverageSum += fMean;
					nSamples++;
				}

				pMeanBackgroundPerIntron[nIntronIdx] = fBackgroundIntronCoverageSum / (double)nSamples;
				
				nIntronIdx++;
			}
			
			mapNoiseToCondition.put(strCondition, StatUtils.mean(pMeanBackgroundPerIntron));
		}
			
		//###########################################
		//    test for all condition combinations
		//###########################################
		for(int i=0; i<nConditions; i++)
		{
			// skip invalid conditions
			if(vcInvalidConditions.contains(pConditions[i]))
				continue;
			
			for(int j=i+1; j<nConditions; j++)
			{
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

				ExonGroup pExonGroups[] = null;
				if(bIgnoreFirstAndLastExons)
					pExonGroups = data.GetMiddleExonGroups();
				else
					pExonGroups = data.GetExonGroups();
				
				// calculate fractions for the selected conditions
				for(ExonGroup grp : pExonGroups)
				{
					Vector<Double> vcFractions = new Vector<Double>();
					
					//########################################################
					//    now calculate the coverage ratio for exon groups
					//########################################################
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
								pCov[nSample] = data.GetCoverageArrayForMiddleExonGroupAndSample(grp, strSample)[x];
							else
								pCov[nSample] = data.GetCoverageArrayForExonGroupAndSample(app, grp, strSample)[x];
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
								pCov[nSample] = data.GetCoverageArrayForMiddleExonGroupAndSample(grp, strSample)[x];
							else
								pCov[nSample] = data.GetCoverageArrayForExonGroupAndSample(app, grp, strSample)[x];
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

				//########################################################
				//    now calculate the coverage ratio for the introns
				//########################################################
				TreeMap<Exon, double[]> mapCovFractionsPerIntron = new TreeMap<Exon, double[]>();
				for(Exon intron : data.GetIntrons())
				{
					//##############################################################
					//   check whether the intron is valid (sufficiently covered)
					//##############################################################

					double pCovPerBaseValues[] = new double[vcSamplesA.size()];
					double pUncoveredBases[]   = new double[vcSamplesA.size()];
					
					int nSample = 0;
					for(String strSample : vcSamplesA)
					{
						double pCov[] = data.GetCoverageArrayForIntronAndSample(intron, strSample);
						
						int nUncoveredBases = 0;
						for(double fVal : pCov)
						{
							if(fVal == 0)
								nUncoveredBases++;
						}
						
						pCovPerBaseValues[nSample] = StatUtils.mean(pCov);
						pUncoveredBases[nSample] = nUncoveredBases;
						
						nSample++;
					}
					double fMeanCovPerBaseA = StatUtils.mean(pCovPerBaseValues);
					double fMeanUncoveredBasesA = StatUtils.mean(pUncoveredBases);
					
					pCovPerBaseValues = new double[vcSamplesB.size()];
					pUncoveredBases   = new double[vcSamplesB.size()];
					nSample = 0;
					for(String strSample : vcSamplesB)
					{
						double pCov[] = data.GetCoverageArrayForIntronAndSample(intron, strSample);
						
						int nUncoveredBases = 0;
						for(double fVal : pCov)
						{
							if(fVal == 0)
								nUncoveredBases++;
						}
						
						pCovPerBaseValues[nSample] = StatUtils.mean(pCov);
						pUncoveredBases[nSample] = nUncoveredBases;
						
						nSample++;
					}
					double fMeanCovPerBaseB = StatUtils.mean(pCovPerBaseValues);
					double fMeanUncoveredBasesB = StatUtils.mean(pUncoveredBases);
					
					if( (fMeanCovPerBaseA < nMinCovPerBase || fMeanCovPerBaseA < mapNoiseToCondition.get(strConditionA)*2) &&
						(fMeanCovPerBaseB < nMinCovPerBase || fMeanCovPerBaseB < mapNoiseToCondition.get(strConditionB)*2))
						continue;
					
					if(1 - fMeanUncoveredBasesA / intron.getGenomicLength() < fMinCoveredBases && 1 - fMeanUncoveredBasesB / intron.getGenomicLength() < fMinCoveredBases)
						continue;

					//######################################
					//           calculate ratios
					//######################################
					Vector<Double> vcIntronFractions = new Vector<Double>();
					
					if(intron.getGenomicLength() < 2000)
					{
						for(int x=0; x<intron.getGenomicLength(); x++)
						{
							//##########################################################
							//   get mean coverage for this intron in condition A
							//##########################################################
							double pCov[] = new double[vcSamplesA.size()];
							nSample = 0;
							for(String strSample : vcSamplesA)
							{
								pCov[nSample] = data.GetCoverageArrayForIntronAndSample(intron, strSample)[x];
								nSample++;
							}
							double fMeanCoverageA = StatUtils.mean(pCov);
	
							//##########################################################
							//   get mean coverage for this intron in condition B
							//##########################################################
							pCov = new double[vcSamplesB.size()];

							nSample = 0;
							for(String strSample : vcSamplesB)
							{
								pCov[nSample] = data.GetCoverageArrayForIntronAndSample(intron, strSample)[x];
								nSample++;
							}
							
							double fMeanCoverageB = StatUtils.mean(pCov);
	
							//############################################
							//   get total coverage and calculate ratio
							//############################################
							double fTotalCoverage   = fMeanCoverageA + fMeanCoverageB;
	
							// calculate ratio, skip uncovered positions
							if(fTotalCoverage != 0)
								vcIntronFractions.add(fMeanCoverageA / fTotalCoverage);
						}
					}
					else
					{
						// use an approximation for large introns
						//##########################################################
						//   get mean coverage for this intron in condition A
						//##########################################################
						double pCov[] = new double[vcSamplesA.size()];
						nSample = 0;
						for(String strSample : vcSamplesA)
						{
							pCov[nSample] = StatUtils.mean(data.GetCoverageArrayForIntronAndSample(intron, strSample));
							nSample++;
						}
						double fMeanCoverageA = StatUtils.mean(pCov);

						//##########################################################
						//   get mean coverage for this intron in condition B
						//##########################################################
						pCov = new double[vcSamplesB.size()];
						
						nSample = 0;
						for(String strSample : vcSamplesB)
						{
							pCov[nSample] = StatUtils.mean(data.GetCoverageArrayForIntronAndSample(intron, strSample));
							nSample++;
						}
						
						double fMeanCoverageB = StatUtils.mean(pCov);

						//############################################
						//   get total coverage and calculate ratio
						//############################################
						double fTotalCoverage   = fMeanCoverageA + fMeanCoverageB;

						// calculate ratio, skip uncovered positions
						if(fTotalCoverage != 0)
							vcIntronFractions.add(fMeanCoverageA / fTotalCoverage);
					}
					
					double[] pIntronFractions = new double[vcIntronFractions.size()];
					for(int k=0; k<vcIntronFractions.size(); k++)
						pIntronFractions[k] = vcIntronFractions.get(k);
					
					mapCovFractionsPerIntron.put(intron, pIntronFractions);
				}
				
				// perform test for each exon group to identify alternatively spliced exons
				for(Exon intron : mapCovFractionsPerIntron.keySet())
				{
					double[] pFractionsGroup1 = mapCovFractionsPerIntron.get(intron);
					double fAvgGroup1 = StatUtils.mean(pFractionsGroup1);
	
					// get average coverage ratio for all other exons
					Vector<Double> vcAllGroups = new Vector<Double>();
					for(ExonGroup grp : mapCovFractionsPerExonGroup.keySet())
					{					
						double pFractionsGroup[] = mapCovFractionsPerExonGroup.get(grp);
						vcAllGroups.add(StatUtils.mean(pFractionsGroup));
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
					
					if(fPValue < 0.05 && Math.abs(fAbsDifference) >= fVariableExonThreshold)
					{
						double fMeanCoveragePerBaseA = 0.0;
						double fMeanCoveragePerBaseB = 0.0;
						
						// get mean coverage per base for condition A
						double pMeanCovA[] = new double[vcSamplesA.size()];
						int nSample=0;
						for(String strSample : vcSamplesA)
						{
							double pCov[] = data.GetCoverageArrayForIntronAndSample(intron, strSample);
							pMeanCovA[nSample] = StatUtils.mean(pCov);
							nSample++;
						}
						fMeanCoveragePerBaseA = StatUtils.mean(pMeanCovA);
						
						// get coverage array for samples in condition B
						double pMeanCovB[] = new double[vcSamplesB.size()];
						nSample=0;
						for(String strSample : vcSamplesB)
						{
							double pCov[] = data.GetCoverageArrayForIntronAndSample(intron, strSample);
							pMeanCovB[nSample] = StatUtils.mean(pCov);
							nSample++;
						}
						fMeanCoveragePerBaseB = StatUtils.mean(pMeanCovB);
						
						/*
						if(ex.getCodingStart() == 203702351 && ex.getCodingStop() == 203702528)
						{
							System.out.println("adding result for: " + strConditionA + " vs " + strConditionB);
						}
						*/
						
						result.AddASExons(SplicingWebApp.AS_TYPE_RETAINED_INTRON, data.GetReferenceName() + ":" + intron.getCodingStart() + "-" + intron.getCodingStop(), intron.getCodingStart(), intron.getCodingStop(),
								strConditionA, strConditionB, fMeanCoveragePerBaseA, fMeanCoveragePerBaseB, fAbsDifference, fRelDifference, fPValue, pFractionsGroup1, pAllGroups);
					}
				}
			}
		}

		return result;
	}
}
