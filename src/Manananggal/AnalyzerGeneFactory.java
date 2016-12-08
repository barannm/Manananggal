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

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Calendar;
import java.util.TreeSet;

import BioKit.GTFGene;
import BioKit.GTFParser;
import BioKit.Gene;
import BioKit.RandomAccessGFFReader;

/**
 *    The analyzer gene factory (I should rename it some time) includes the console application
 *    functions that are necessary to search for AS events in all genes. 
 */
public class AnalyzerGeneFactory
{
	/**
	 *    RunCompleteAnalysis starts the analysis for all genes included in the GTF file (specified as parameter) for a given condition type.
	 *    Optionally, if strFirstGene is specified, the analysis will start at the specified gene entry, thus skipping all previous genes.
	 *    If debug is specified, only a single gene (must be specified via strFirstGene) will be analyzed.
	 */
	public void RunCompleteAnalysis(SplicingWebApp app, String strFileGTF, String strFileProject, String strConditionType, String strFirstGene, boolean bSkipFirstAndLastExon, boolean bDebug, int nThreads) throws IOException
	{
		//###################################################################
		//   retrieve data for the current selections and gene information
		//###################################################################
		TreeSet<String> vcSelectedSamples 	= app.GetSelectedSamples();
		int nMinCovPerBase 					= app.GetMinimumCoveragePerBase();
		int nMinJunctionReads 				= app.GetMinimumJunctionReads();
		double fMinCoveredBases				= app.GetMinimumCoveredBases();
		double fVariableExonThreshold 		= app.GetVariableExonThreshold();		
		//###################################################################
		
		//TRUE
		boolean bShowProcessingTime = false;
		
		// for debugging
		long nTimeA = 0;
		long nTimeB = 0;
		
		boolean bOkay 	= true;
		if(strFirstGene != null)
			bOkay = false;
		
		System.out.println("GTF file: " + strFileGTF);
		System.out.println("project: " + strFileProject);
		System.out.println("Condition type: " + strConditionType);
		System.out.println("skipping first and last exons: " + bSkipFirstAndLastExon);
		System.out.println("debug mode: " + bDebug);
		System.out.println("starting at: " + strFirstGene);
		
		app.SetConditionType(strConditionType);
		
		// init project		
		app.InitProjectModel(strFileProject);
		ProjectModel projectModel = app.GetProjectModel();
		
		app.SetNumberOfThreads(nThreads);
		
		for(String strSample : projectModel.GetSamples())
			vcSelectedSamples.add(strSample);
		
		PrintWriter pWriterSE 	= null;
		PrintWriter pWriterATE 	= null;
		PrintWriter pWriterRE 	= null;
		PrintWriter pWriterEE 	= null;
		
		for(int i=0; i<4; i++)
		{
			String strTimeStamp = new SimpleDateFormat("yyyy_MM_dd").format(Calendar.getInstance().getTime());
			String strFileName = strTimeStamp + "_" + projectModel.GetProjectName() + "_" + nMinJunctionReads +"_" + nMinCovPerBase + "_" + fMinCoveredBases + "_" + fVariableExonThreshold + "_" + bSkipFirstAndLastExon;
			
			PrintWriter pOut = null;			
			switch(i)
			{
				case 0:
				{
					strFileName += "_SE_output.txt";
					
					if(!bOkay)
						pWriterSE = new PrintWriter(new FileWriter(strFileName, true)); 
					else
						pWriterSE = new PrintWriter(new FileWriter(strFileName, false));
					
					pOut = pWriterSE;
					break;
				}
				case 1:
				{
					strFileName += "_ATE_output.txt";
					
					if(!bOkay)
						pWriterATE = new PrintWriter(new FileWriter(strFileName, true)); 
					else
						pWriterATE = new PrintWriter(new FileWriter(strFileName, false));
					
					pOut = pWriterATE;
					break;
				}
				case 2:
				{
					strFileName += "_RE_output.txt";
					
					if(!bOkay)
						pWriterRE = new PrintWriter(new FileWriter(strFileName, true)); 
					else
						pWriterRE = new PrintWriter(new FileWriter(strFileName, false));
					
					pOut = pWriterRE;
					break;
				}
				case 3:
				{
					strFileName += "_EE_output.txt";
					
					if(!bOkay)
						pWriterEE = new PrintWriter(new FileWriter(strFileName, true)); 
					else
						pWriterEE = new PrintWriter(new FileWriter(strFileName, false));
					
					pOut = pWriterEE;
					break;
				}
			}
			
			strTimeStamp = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(Calendar.getInstance().getTime());
			pOut.println("## " + strTimeStamp);
			pOut.println("## minimum number of junction reads: " 			+ nMinJunctionReads);
			pOut.println("## minimum coverage per base per exon: " 			+ nMinCovPerBase);
			pOut.println("## minimum fraction of covered bases per exon: " 	+ fMinCoveredBases);
			pOut.println("## coverage ratio threshold: " 					+ fVariableExonThreshold + "\n");
			pOut.println("#ID\tENSEMBL_ID\tgene_symbol\ttype\tAS_type\tgrp_position\tconditionA\tconditionB\talt_exon_cov_per_base_A\talt_exon_cov_per_base_B\talt_exon_cov_per_base_diff\tabs_ratio_changes\trel_ratio_changes\tp.values_ratio\tincl_jun_start\tincl_jun_end\texcl_jun_start\texcl_jun_end\teffect_PSI\tp.value_PSI\tis_significant\tuses_novel_junction");
		}
		
		if(strFileGTF.toLowerCase().contains("gff") || strFileGTF.toLowerCase().contains("gtf"))
		{
			GTFParser parser = null;
			
			int nID = 0;
		
			if(strFirstGene != null)
			{
				GeneIdentifier gid = app.GetGeneIdentifierHandler().GetGeneIdentifierForGene(strFirstGene, null);
				
				// check if the index file exists
				File pFileIndex = new File(strFileGTF + ".idx");
		
				try
				{
					RandomAccessGFFReader gffReader = new RandomAccessGFFReader(new File(strFileGTF), pFileIndex);
					long nOffset = gffReader.GetIndexForGene(gid.m_strEnsemblGeneID);
					
					if(nOffset != -1)
					{
						RandomAccessFile pFile = new RandomAccessFile(strFileGTF, "r");
						parser = new GTFParser(pFile, nOffset);
					}
				}
				catch(Exception e)
				{
					e.printStackTrace();
				}
			}
			
			if(parser == null)
				parser = new GTFParser(strFileGTF);
			
			// analyze all genes in the GTF file			
			while(true)
			{
				if(bShowProcessingTime)
				{
					nTimeA = System.currentTimeMillis();
				}
				
				app.ClearCurrentGeneData();
				
				GTFGene gtf_gene = parser.nextGene();
				if(gtf_gene == null)
					break;
				
				if(bDebug || !bOkay)
				{					
					if(gtf_gene.getGeneName() != null)
					{
						if(!gtf_gene.getGeneName().equals(strFirstGene) && !gtf_gene.getId().equals(strFirstGene))				
							continue;
					}
					else if(!gtf_gene.getId().equals(strFirstGene))
					{
						continue;
					}
					
					bOkay = true;
				}
				
				if(!bOkay)
					continue;
				
				Gene gene = gtf_gene.createGene();

				System.out.print("current gene: " + gene.getGeneName() + " (" + gene.getGeneID() + ")                    \r");

				if(bShowProcessingTime)
				{
					nTimeB = System.currentTimeMillis();
					
					System.out.println("time for gene data preparation: " + (nTimeB - nTimeA));
					
					nTimeA = nTimeB;
				}
				
				// select all isoforms
				app.InitValidIsoforms(gene);

				app.PrepareDataSupplier(gene);
				TreeSet<AnalysisResult> vcResults = AnalyzeCurrentGene(app, bSkipFirstAndLastExon, true, nID, bDebug);
				
				if(bShowProcessingTime)
				{
					nTimeB = System.currentTimeMillis();
					
					System.out.println("time for gene processing: " + (nTimeB - nTimeA));
					
					nTimeA = nTimeB;
				}
				
				nID += vcResults.size();
				
				for(AnalysisResult res : vcResults)
				{
					switch(res.GetType())
					{
						case SplicingWebApp.AS_TYPE_EXON_SKIPPING:
						{
							pWriterSE.println(res.toString());
							break;
						}
						
						case SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN:
						case SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN:
						case SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN_DOUBLE:
						case SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN_DOUBLE:
						case SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN:
						case SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN:
						case SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN_DOUBLE:
						case SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN_DOUBLE:
						{
							pWriterATE.println(res.toString());
							break;
						}
						
						case SplicingWebApp.AS_TYPE_RETAINED_INTRON:
						{
							pWriterRE.println(res.toString());
							break;
						}
						
						case SplicingWebApp.AS_TYPE_ALT_3_PRIME_EXON_END:
						case SplicingWebApp.AS_TYPE_ALT_5_PRIME_EXON_END:
						{
							pWriterEE.println(res.toString());
							break;
						}
					}
				}
					 

				pWriterSE.flush();
				pWriterATE.flush();
				pWriterRE.flush();
				pWriterEE.flush();
				
				if(bDebug)
					break;
			}
		}
		else
		{
			// NOT SUPPORTED
		}
		
		String strTimeStamp = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(Calendar.getInstance().getTime());
		
		pWriterSE.println("finished: " + strTimeStamp);
		pWriterATE.println("finished: " + strTimeStamp);
		pWriterRE.println("finished: " + strTimeStamp);
		pWriterEE.println("finished: " + strTimeStamp);
		
		pWriterSE.close();
		pWriterATE.close();
		pWriterRE.close();
		pWriterEE.close();
	}
	
	/** Invoked from "RunCompleteAnalysis". This functions analyzes the current gene. */
	public TreeSet<AnalysisResult> AnalyzeCurrentGene(SplicingWebApp app, boolean bSkipFirstAndLastExon, boolean bIdentifyValidIsoforms, int nID, boolean bDebug) throws IOException
	{
		//TODO
		// for debugging
		boolean bShowProcessingTime = false;
		long nTimeA = 0;
		long nTimeB = 0;
		
		if(bShowProcessingTime)
		{
			nTimeA = System.currentTimeMillis();
		}
		
		TreeSet<AnalysisResult> vcResults = new TreeSet<AnalysisResult>();
		
		Gene gene = app.GetDataSupplier().GetGene();
		
		if(bShowProcessingTime)
		{
			nTimeB = System.currentTimeMillis();
			System.out.println(" - time for first data prep step: " + (nTimeB - nTimeA));
			nTimeA = nTimeB;
		}
		
		TreeSet<String> vcValidIsoforms = app.GetValidIsoforms();
		if(bDebug)
		{
			System.out.println("\n######### ISOFORMS #########");
			for(String strIsoform : vcValidIsoforms)
				System.out.println(strIsoform);
		}
		
		// identify valid isoforms (used only in the console application, the GUI uses the user selected isoforms)
		if(bIdentifyValidIsoforms)
		{			
			// remove irrelevant isoforms
			if(bDebug) System.out.println("removing irrelevant isoforms");
			if(!app.HideIrrelevantIsoforms(bSkipFirstAndLastExon, bDebug))
			{
				System.out.println("failed to hide irrelevant isoforms");
				return null;
			}
			
			if(bShowProcessingTime)
			{
				nTimeB = System.currentTimeMillis();
				System.out.println(" - time for new isoform detection: " + (nTimeB - nTimeA));
				nTimeA = nTimeB;
			}

			vcValidIsoforms = app.GetValidIsoforms();
			if(bDebug)
			{
				System.out.println("remaining isoforms:");
				for(String strIsoform : vcValidIsoforms)
					System.out.println(strIsoform);
				System.out.print("\n");
			}
	
			// recalculate exon groups based on valid isoforms
			if(bDebug) System.out.println("recalculating exon groups");
			app.GetDataSupplier().RetrieveCoverageData(app);
		}
		
		if(bShowProcessingTime)
		{
			nTimeB = System.currentTimeMillis();
			System.out.println(" - time for new exon group recalculation: " + (nTimeB - nTimeA));
			nTimeA = nTimeB;
		}
		
		// recalculate exon group values
		app.PrepareDataSupplier(gene);
		
		if(bDebug)
		{
			System.out.println("\n######### EXON GROUPS ######### ");
			System.out.println(Arrays.toString(app.GetDataSupplier().GetExonGroups()) + "\n");
		}
		
		if(bShowProcessingTime)
		{
			nTimeB = System.currentTimeMillis();
			System.out.println(" - time for remaining data prep: " + (nTimeB - nTimeA));
			nTimeA = nTimeB;
		}

		//###############################################################
		//                 calculate variable exons
		//###############################################################		

		// get skipped exons
		if(bDebug) System.out.println("calculating variable exons");
		AnalyzerExonSkippingEvents analyzer = new AnalyzerExonSkippingEvents();
		AlternativeSplicingHit res = analyzer.IdentifyVariableExons(app, bSkipFirstAndLastExon);
		
		// calculate and add number of PSI score results and the best p-value				
		TreeSet<SimpleSpliceScore> psi_scores = analyzer.CalculateJunctionPSIScores(app, bDebug, null);
		
		if(bShowProcessingTime)
		{
			nTimeB = System.currentTimeMillis();
			System.out.println(" - time for exon skipping events: " + (nTimeB - nTimeA));
			nTimeA = nTimeB;
		}

		// get alternative transcript start or end exons
		if(bDebug) System.out.println("identifying alternative terminal exons");
		AnalyzerAltStartOrEndEvents analyzer2 = new AnalyzerAltStartOrEndEvents();
		analyzer2.IdentifyVariableStartOrEndExons(app, res);
		TreeSet<SimpleSpliceScore> psi_scores_terminal_exons = analyzer2.CalculateJunctionInclusionScoresForFirstAndLastExons(app, bDebug, null);

		if(bShowProcessingTime)
		{
			nTimeB = System.currentTimeMillis();
			System.out.println(" - time for alternative terminal exons: " + (nTimeB - nTimeA));
			nTimeA = nTimeB;
		}
		
		// get retained introns
		if(app.DetectIntronRetentionEvents())
		{
			if(bDebug) System.out.println("calculating retained introns");
			AnalzyerRetainedIntrons analyzer3 = new AnalzyerRetainedIntrons();
			analyzer3.IdentifyIntronRetentionEvents(app, res, bSkipFirstAndLastExon);
			
			if(bShowProcessingTime)
			{
				nTimeB = System.currentTimeMillis();
				System.out.println(" - time for retained introns: " + (nTimeB - nTimeA));
				nTimeA = nTimeB;
			}
		}
		
		if(bDebug)
		{
			if(res != null)
			{
				System.out.println("########### variable exon coverage ###########");
				System.out.println(res);
			}
		}
		
		// get exon extensions
		if(bDebug) System.out.println("identifying exon extensions");
		AnalyzerExonExtensions analyzer4 = new AnalyzerExonExtensions();
		TreeSet<SimpleSpliceScore> psi_scores_exon_extensions = analyzer4.IdentifyExonExtensionEvents(app, bDebug, null);
		
		if(bShowProcessingTime)
		{
			nTimeB = System.currentTimeMillis();
			System.out.println(" - time for exon extensions: " + (nTimeB - nTimeA));
			nTimeA = nTimeB;
		}

		if(bDebug)
		{
			System.out.println("########### psi_scores ###########");
			System.out.println(psi_scores);
			
			System.out.println("########### psi_scores first + last exons ###########");
			System.out.println(psi_scores_terminal_exons);
			
			System.out.println("########### psi exon extensions ###########");
			System.out.println(psi_scores_exon_extensions);
		}

		// add to output
		if(bDebug) System.out.println("writing output");			
		TreeSet<SimpleSpliceScore> vcUsedSpliceEvent = new TreeSet<SimpleSpliceScore>();
		TreeSet<SimpleSpliceScore> vcUsedSpliceEventTerminalExons = new TreeSet<SimpleSpliceScore>();

		if(res != null)
		{
			int nResults = res.GetNumHits();
			
			for(int k=0; k<nResults; k++)
			{
				AlternativeSplicingExon hit = res.GetResult(k);
				int nStart = hit.GetStart();
				int nEnd   = hit.GetEnd();

				if(bDebug)
				{
					System.out.println("hit: " + nStart + " - " + nEnd);
				}
			
				boolean bOverlapWithSpliceScore = false;
				SimpleSpliceScore matchScore 	= null;
				
				boolean bDoubleSupport = false;
				
				AlternativeSplicingExon hit2	= null;
				AlternativeSplicingExon hitTmp	= null;

				for(SimpleSpliceScore r : psi_scores_terminal_exons)
				{
					if((hit.GetType() == SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN || hit.GetType() == SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN) && r.GetType() == SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN)
					{
						int nValidType = 0;
						
						if(gene.getStrandStringID() == '+')
						{
							if(r.m_JunctionExclusion.m_nStart == nEnd)
								nValidType = 1;
							
							if(r.m_JunctionInclusion.m_nStart == nEnd)
								nValidType = 2;
						}
						else
						{
							if(r.m_JunctionExclusion.m_nEnd == nStart)
								nValidType = 1;
							
							if(r.m_JunctionInclusion.m_nEnd == nStart)
								nValidType = 2;
						}
						
						boolean bHasDoubleSupport = false;
						
						if(	(nValidType != 0) &&
								( (r.m_strConditionA.equals(hit.GetCondition(true)) && r.m_strConditionB.equals(hit.GetCondition(false)) ||
								  (r.m_strConditionB.equals(hit.GetCondition(true)) && r.m_strConditionA.equals(hit.GetCondition(false))))))
						{
							/*
							if(bDebug)
							{
								System.out.println("matching junction hit: " + r + " valid_type = " + nValidType);
							}
							*/
							
							// the other junction should also be connected to an exon that has a coverage ratio change									
							for(int l=0; l<nResults; l++)
							{
								hitTmp = res.GetResult(l);
								
								if(hit.equals(hitTmp))
									continue;
								
								if(hitTmp.GetType() != SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN && hitTmp.GetType() != SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN)
									continue;
								
								if(!(r.m_strConditionA.equals(hitTmp.GetCondition(true)) && r.m_strConditionB.equals(hitTmp.GetCondition(false)) ) ||
									(r.m_strConditionB.equals(hitTmp.GetCondition(true)) && r.m_strConditionA.equals(hitTmp.GetCondition(false)) ))
									continue;
								
								// check for double junction support
								if(gene.getStrandStringID() == '+')
								{
									if(nValidType == 1 && r.m_JunctionInclusion.m_nStart == hitTmp.GetEnd())
									{
										bHasDoubleSupport = true;
										break;
									}
									
									if(nValidType == 2 && r.m_JunctionExclusion.m_nStart == hitTmp.GetEnd())
									{
										bHasDoubleSupport = true;
										break;
									}
								}
								else
								{
									if(nValidType == 1 && r.m_JunctionInclusion.m_nEnd == hitTmp.GetStart())
									{
										bHasDoubleSupport = true;
										break;
									}
									
									if(nValidType == 2 && r.m_JunctionExclusion.m_nEnd == hitTmp.GetStart())
									{
										bHasDoubleSupport = true;
										break;
									}
								}
							}
							
							if(matchScore == null || r.m_fPValue < matchScore.m_fPValue || (!bDoubleSupport && bHasDoubleSupport))
							{
								// if there is no support of both regions, just use the new psi_score result, because it has the better p-value
								// if a previous result had suporrt for both regions, only overwrite it if the current result also has support for both regions 
								if(bHasDoubleSupport)
								{
									if(hit2 == null || hit2.GetType() != SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN || hitTmp.GetType() == SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN)
									{
										bDoubleSupport = true;
										matchScore = r;
										
										hit2 = hitTmp;
										/*
										if(bDebug)
										{
											System.out.println("HAS DOUBLE SUPPORT");
											System.out.println(hitTmp.GetStart() + "-" + hitTmp.GetEnd());
										}
										*/
									}
								}
								else if(!bDoubleSupport)
									matchScore = r;
							}
							
							vcUsedSpliceEventTerminalExons.add(r);
						}
					}
					else if((hit.GetType() == SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN || hit.GetType() == SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN) && r.GetType() == SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN)
					{
						int nValidType = 0;
							
						if(gene.getStrandStringID() == '-')
						{
							if(r.m_JunctionExclusion.m_nStart == nEnd)
								nValidType = 1;
							
							if(r.m_JunctionInclusion.m_nStart == nEnd)
								nValidType = 2;
						}
						else
						{
							if(r.m_JunctionExclusion.m_nEnd == nStart)
								nValidType = 1;
							
							if(r.m_JunctionInclusion.m_nEnd == nStart)
								nValidType = 2;
						}
														
						boolean bHasDoubleSupport = false;

						if(	(nValidType != 0) &&
								( (r.m_strConditionA.equals(hit.GetCondition(true)) && r.m_strConditionB.equals(hit.GetCondition(false)) ||
								  (r.m_strConditionB.equals(hit.GetCondition(true)) && r.m_strConditionA.equals(hit.GetCondition(false))))))
						{
							/*
							if(bDebug)
							{
								System.out.println("matching junction hit: " + r + " valid_type = " + nValidType);
							}
							*/
							
							// the other junction should also be connected to an exon that has a coverage ratio change
							for(int l=0; l<nResults; l++)
							{
								hitTmp = res.GetResult(l);
								
								if(hit.equals(hitTmp))
									continue;
								
								if(hitTmp.GetType() != SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN && hitTmp.GetType() != SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN)
									continue;
								
								if(!(r.m_strConditionA.equals(hitTmp.GetCondition(true)) && r.m_strConditionB.equals(hitTmp.GetCondition(false)) ) ||
									(r.m_strConditionB.equals(hitTmp.GetCondition(true)) && r.m_strConditionA.equals(hitTmp.GetCondition(false)) ))
									continue;
								
								if(gene.getStrandStringID() == '-')
								{
									if(nValidType == 1 && r.m_JunctionInclusion.m_nStart == hitTmp.GetEnd())
									{
										bHasDoubleSupport = true;
										break;
									}
									
									if(nValidType == 2 && r.m_JunctionExclusion.m_nStart == hitTmp.GetEnd())
									{
										bHasDoubleSupport = true;
										break;
									}
								}
								else
								{
									if(nValidType == 1 && r.m_JunctionInclusion.m_nEnd == hitTmp.GetStart())
									{
										bHasDoubleSupport = true;
										break;
									}
									
									if(nValidType == 2 && r.m_JunctionExclusion.m_nEnd == hitTmp.GetStart())
									{
										bHasDoubleSupport = true;
										break;
									}
								}
							}

							if(matchScore == null || r.m_fPValue < matchScore.m_fPValue || (!bDoubleSupport && bHasDoubleSupport))
							{
								// if there is no support of both regions, just use the new psi_score result, because it has the better p-value
								// if a previous result had suporrt for both regions, only overwrite it if the current result also has support for both regions 
								if(bHasDoubleSupport)
								{
									if(hit2 == null || hit2.GetType() != SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN || hitTmp.GetType() == SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN)
									{
										bDoubleSupport = true;
										matchScore = r;
										
										hit2 = hitTmp;
										
										/*
										if(bDebug)
										{
											System.out.println("HAS DOUBLE SUPPORT");
											System.out.println(hitTmp.GetStart() + "-" + hitTmp.GetEnd());
										}
										*/
									}
								}
								else if(!bDoubleSupport)
									matchScore = r;
							}
							
							vcUsedSpliceEventTerminalExons.add(r);
						}
					}
				}
					
				if(matchScore != null)
				{
					if(bDoubleSupport)
					{
//						System.out.println("hit1: " + hit.GetStart() + "-" + hit.GetEnd() + " " + hit.GetType());
//						System.out.println("hit2: " + hit2.GetStart() + "-" + hit2.GetEnd() + " " + hit2.GetType());
//						System.out.println("changing type from: " + hit.GetType());
						
						switch(hit.GetType())
						{
							case SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN:
							{
								if(hit2.GetType() == SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN)
									hit.SetType(SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN_DOUBLE);
								else
									hit.SetType(SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN_DOUBLE);
								break;
							}
							case SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN: { hit.SetType(SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN_DOUBLE); break; }
							case SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN:
							{
								if(hit2.GetType() == SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN)
									hit.SetType(SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN_DOUBLE);
								else
									hit.SetType(SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN_DOUBLE);
								break;
							}
							case SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN:   { hit.SetType(SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN_DOUBLE); break; }
						}
						
//						System.out.println("to: " + hit.GetType());
					}
					
					bOverlapWithSpliceScore = true;
					
					AnalysisResult result = new AnalysisResult(app, nID, gene.getChromosome(), gene.getGeneID(), gene.getGeneName(), hit, hit2, matchScore);
					vcResults.add(result);
					nID++;
				}
				
				for(SimpleSpliceScore r : psi_scores)
				{
					if(hit.GetType() != SplicingWebApp.AS_TYPE_EXON_SKIPPING || r.GetType() != SplicingWebApp.AS_TYPE_EXON_SKIPPING)
						continue;
					
					// Find the PSI scores with best p-values, one using only known junctions, one including novel junctions
					if(r.AffectsExon(nStart, nEnd) &&
							( (r.m_strConditionA.equals(hit.GetCondition(true)) && r.m_strConditionB.equals(hit.GetCondition(false)) ||
							  (r.m_strConditionB.equals(hit.GetCondition(true)) && r.m_strConditionA.equals(hit.GetCondition(false)) ))))
					{
						bOverlapWithSpliceScore = true;
						
						AnalysisResult result = new AnalysisResult(app, nID, gene.getChromosome(), gene.getGeneID(), gene.getGeneName(), hit, null, r);						
						vcResults.add(result);
						vcUsedSpliceEvent.add(r);
						nID++;
					}
				}
				
				if(!bOverlapWithSpliceScore)
				{
					AnalysisResult result = new AnalysisResult(app, nID, gene.getChromosome(), gene.getGeneID(), gene.getGeneName(), hit, null, null);
					vcResults.add(result);
					nID++;
				}
			}
		}
		
		for(SimpleSpliceScore r : psi_scores)
		{
			// the result should have an change of at least 1% 
			if(!vcUsedSpliceEvent.contains(r) && r.m_fInclusionChange >= 0.01)
			{
				AnalysisResult result = new AnalysisResult(app, nID, gene.getChromosome(), gene.getGeneID(), gene.getGeneName(), null, null, r);
				vcResults.add(result);
				nID++;
			}
		}
		
		for(SimpleSpliceScore r : psi_scores_terminal_exons)
		{
			// the result should have an change of at least 1%
			if(!vcUsedSpliceEvent.contains(r) && r.m_fInclusionChange >= 0.01)
			{
				AnalysisResult result = new AnalysisResult(app, nID, gene.getChromosome(), gene.getGeneID(), gene.getGeneName(), null, null, r);
				vcResults.add(result);
				nID++;
			}
		}
		
		for(SimpleSpliceScore r : psi_scores_exon_extensions)
		{
			// the result should have an change of at least 1%
			if(!vcUsedSpliceEvent.contains(r) && r.m_fInclusionChange >= 0.01)
			{
				AnalysisResult result = new AnalysisResult(app, nID, gene.getChromosome(), gene.getGeneID(), gene.getGeneName(), null, null, r);
				vcResults.add(result);
				nID++;
			}
		}
		
		if(bShowProcessingTime)
		{
			nTimeB = System.currentTimeMillis();
			System.out.println(" - time for result generation: " + (nTimeB - nTimeA));
			nTimeA = nTimeB;
		}
		
		return vcResults;
	}
}
