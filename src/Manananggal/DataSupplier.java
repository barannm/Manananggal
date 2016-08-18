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
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.math3.stat.StatUtils;
import org.zkoss.zul.Messagebox;

import BioKit.Exon;
import BioKit.ExonGroup;
import BioKit.Gene;

/**
 *    The DataSupplier is a crucial class that provides coverage and junction
 *    count data to other functions. It is first filled with data on gene
 *    selection. Changing the isoform selection will also update the data
 *    stored in the data supplier (e.g. exon groups must be re-calculated if
 *    the selected isoforms change).
 */
public class DataSupplier
{
	private Gene												m_Gene;
	private Exon[]												m_pExons;			// all exons
	private ExonGroup[]											m_pExonGroups;		// all exon groups
	private GeneIdentifier										m_GID;				// gene identifier (multiple annotations)
	
	private TreeSet<Exon>										m_vcIntrons; 		// set of introns
	private TreeSet<Exon> 										m_vcMiddleExons; 	// set of all middle exons 
	private TreeSet<Exon> 										m_vcFirstExons;		// set of 5' terminal exons
	private TreeSet<Exon> 										m_vcLastExons;		// set of 3' terminal exons
	private ExonGroup[]											m_pMiddleExonGroups; // array of exon groups, excluding all terminal exons

	// coverage container for each exon group, exon and intron
	// 1st key = exon/group, 2nd key = sample, value = coverage (adjusted by size factors) 
	private TreeMap<String, double[]>							m_mapGeneCoveragePerSample;
	private TreeMap<ExonGroup, TreeMap<String, double[]>> 		m_mapCoverageToExonGroupsPerSample;
	private TreeMap<ExonGroup, TreeMap<String, double[]>> 		m_mapCoverageToMiddleExonGroupsPerSample;
	private TreeMap<Exon, TreeMap<String, double[]>> 			m_mapCoverageToExonsPerSample;
	private TreeMap<Exon, TreeMap<String, double[]>> 			m_mapCoverageToIntronsPerSample;
	
	private TreeMap<Exon, TreeMap<String, Double>>				m_mapRawExonCoveragePerSample;		// raw values, based on wig coverage
	private TreeMap<Exon, TreeMap<String, Double>>				m_mapAdjustedExonCoveragePerSample;	// adjusted coverage values (taking into account ambigious exons)
	private TreeMap<Exon, TreeMap<String, Double>>				m_mapExonCoveredBasesPerSample;
	private TreeMap<Exon, TreeMap<String, TreeSet<Exon>>>		m_mapIndistinguishableExonPerSample;
	
	private TreeMap<CountElement, TreeMap<String, Integer>> 	m_mapJunctionCounts;
	private TreeMap<String, TreeSet<CountElement>>				m_mapJunctionsToIsoforms;
	private TreeSet<CountElement>								m_vcNovelJunctions;
	private TreeSet<CountElement>								m_vcInvalidJunctions;	// invalid due to too low coverage
	private TreeSet<Exon>										m_vcRetainedIntrons;
	
	private TreeMap<String, Exon[]>								m_mapExonsToIsoforms;
	
	private TreeSet<String> 									m_vcInvalidConditions; // set of all conditions with insufficient coverage for the current gene
	
	public DataSupplier()
	{
		m_Gene				= null;
		m_pExons 			= null;
		m_pExonGroups		= null;
		m_GID				= null;
		
		// prepare a set of introns
		m_vcIntrons = new TreeSet<Exon>();
		
		// prepare a set of all exons split by 5' (first) and 3' (last) and middle exons 
		m_vcMiddleExons = new TreeSet<Exon>();
		m_vcFirstExons = new TreeSet<Exon>();
		m_vcLastExons  = new TreeSet<Exon>();
		
		// prepare array of exon groups, excluding all terminal exons
		m_pMiddleExonGroups = null;
		
		// prepare coverage container
		m_mapGeneCoveragePerSample					= new TreeMap<String, double[]>();
		m_mapCoverageToExonGroupsPerSample 			= new TreeMap<ExonGroup, TreeMap<String, double[]>>();
		m_mapCoverageToMiddleExonGroupsPerSample 	= new TreeMap<ExonGroup, TreeMap<String, double[]>>();
		m_mapCoverageToExonsPerSample	 			= new TreeMap<Exon, TreeMap<String, double[]>>();
		m_mapCoverageToIntronsPerSample 			= new TreeMap<Exon, TreeMap<String, double[]>>();
		
		m_mapRawExonCoveragePerSample			= new TreeMap<Exon, TreeMap<String, Double>>();
		m_mapAdjustedExonCoveragePerSample		= new TreeMap<Exon, TreeMap<String, Double>>();
		m_mapExonCoveredBasesPerSample			= new TreeMap<Exon, TreeMap<String, Double>>();
		m_mapIndistinguishableExonPerSample		= new TreeMap<Exon, TreeMap<String, TreeSet<Exon>>>();
		
		m_mapJunctionCounts					= new TreeMap<CountElement, TreeMap<String, Integer>>();
		m_mapJunctionsToIsoforms			= new TreeMap<String, TreeSet<CountElement>>();
		m_vcNovelJunctions					= new TreeSet<CountElement>();
		m_vcInvalidJunctions				= new TreeSet<CountElement>();
		
		m_mapExonsToIsoforms				= new TreeMap<String, Exon[]>();
		
		m_vcInvalidConditions = new TreeSet<String>();
	}
	
	public DataSupplier(SplicingWebApp app, Gene gene)
	{
		m_Gene				= gene;
		m_pExons 			= null;
		m_pExonGroups		= null;
		
		m_GID 				= app.GetGeneIdentifierHandler().GetGeneIdentifierForGene(m_Gene.getGeneID(), m_Gene);
		
		// prepare a set of introns
		m_vcIntrons = new TreeSet<Exon>();
		
		// prepare a set of all exons split by 5' (first) and 3' (last) and middle exons 
		m_vcMiddleExons = new TreeSet<Exon>();
		m_vcFirstExons = new TreeSet<Exon>();
		m_vcLastExons  = new TreeSet<Exon>();
		
		// prepare array of exon groups, excluding all terminal exons
		m_pMiddleExonGroups = null;
		
		// prepare coverage container
		m_mapGeneCoveragePerSample			= new TreeMap<String, double[]>();
		m_mapCoverageToExonGroupsPerSample 	= new TreeMap<ExonGroup, TreeMap<String, double[]>>();
		m_mapCoverageToMiddleExonGroupsPerSample 	= new TreeMap<ExonGroup, TreeMap<String, double[]>>();
		m_mapCoverageToExonsPerSample 		= new TreeMap<Exon, TreeMap<String, double[]>>();
		m_mapCoverageToIntronsPerSample 	= new TreeMap<Exon, TreeMap<String, double[]>>();
		
		m_mapRawExonCoveragePerSample		= new TreeMap<Exon, TreeMap<String, Double>>();
		m_mapAdjustedExonCoveragePerSample	= new TreeMap<Exon, TreeMap<String, Double>>();
		m_mapExonCoveredBasesPerSample		= new TreeMap<Exon, TreeMap<String, Double>>();
		m_mapIndistinguishableExonPerSample	= new TreeMap<Exon, TreeMap<String, TreeSet<Exon>>>();
		
		m_mapJunctionCounts					= new TreeMap<CountElement, TreeMap<String, Integer>>();
		m_mapJunctionsToIsoforms			= new TreeMap<String, TreeSet<CountElement>>();
		m_vcNovelJunctions					= new TreeSet<CountElement>();
		m_vcInvalidJunctions				= new TreeSet<CountElement>();
		
		m_mapExonsToIsoforms				= new TreeMap<String, Exon[]>();
		
		m_vcInvalidConditions 				= new TreeSet<String>();
		
		// process gene information (get introns, terminal exons, etc.)
		ProcessGeneStructure();
	}
	
	/** clears all data stored in the data supplier */
	public void Clear()
	{
		m_Gene				= null;
		m_pExons 			= null;
		m_pExonGroups		= null;
		
		// prepare a set of introns
		m_vcIntrons = new TreeSet<Exon>();
		
		// prepare a set of all exons split by 5' (first) and 3' (last) and middle exons 
		m_vcMiddleExons = new TreeSet<Exon>();
		m_vcFirstExons = new TreeSet<Exon>();
		m_vcLastExons  = new TreeSet<Exon>();
		
		// prepare array of exon groups, excluding all terminal exons
		m_pMiddleExonGroups = null;
		
		// prepare coverage container
		m_mapGeneCoveragePerSample			= new TreeMap<String, double[]>();
		m_mapCoverageToExonGroupsPerSample 	= new TreeMap<ExonGroup, TreeMap<String, double[]>>();
		m_mapCoverageToMiddleExonGroupsPerSample 	= new TreeMap<ExonGroup, TreeMap<String, double[]>>();
		m_mapCoverageToExonsPerSample 		= new TreeMap<Exon, TreeMap<String, double[]>>();
		m_mapCoverageToIntronsPerSample 	= new TreeMap<Exon, TreeMap<String, double[]>>();
		
		m_mapJunctionCounts					= new TreeMap<CountElement, TreeMap<String, Integer>>();
		m_mapJunctionsToIsoforms			= new TreeMap<String, TreeSet<CountElement>>();
		
		m_mapExonsToIsoforms				= new TreeMap<String, Exon[]>();
		m_vcNovelJunctions					= new TreeSet<CountElement>();
		
		m_vcInvalidConditions = new TreeSet<String>();
	}
	
	/**
	 *    Maps junctions to isoform names.
	 *    Used to retrieve all junctions associated
	 *    with a certain isoform.
	 */
	private void MapJunctionsToIsoforms()
	{
		for(String strIsoform : m_Gene.getArrayOfGeneProductNames())
		{
			TreeSet<String> vcJunctionStrings = m_Gene.getSpliceJunctionInformationForGeneProduct(strIsoform);
			TreeSet<CountElement> vcJunctions = new TreeSet<CountElement>();
			
			for(String strJunction : vcJunctionStrings)
			{				
				String pSplit[] = strJunction.split("-");
				
				CountElement jun = new CountElement();
				jun.m_nStart = Integer.parseInt(pSplit[0]);
				jun.m_nEnd	 = Integer.parseInt(pSplit[1]);
				jun.m_chStrand = m_Gene.getStrandStringID();
				jun.m_bKnown	= true;
				
				vcJunctions.add(jun);
			}
			
			m_mapJunctionsToIsoforms.put(strIsoform, vcJunctions);
		}
	}

	/**
	 *    Maps exons to isoforms.
	 *    Used to retrieve all exons that belong
	 *    to a specific isoform.
	 */
	private void MapExonsToIsoforms()
	{
		for(String strIsoform : m_Gene.getArrayOfGeneProductNames())
		{
			Exon pExons[] = m_Gene.getSortedExonsForGeneProduct(strIsoform, true);
			m_mapExonsToIsoforms.put(strIsoform, pExons);
		}
	}
	
	/**
	 *    Identifies annotated exons that are retained introns.
	 *    Exons that represent retained introns are those exons
	 *    that overlap with two unjoined exons or that contain
	 *    a complete splice junction (from start to end).
	 */
	private void IdentifyRetainedIntrons()
	{
		m_vcRetainedIntrons = new TreeSet<Exon>();
		
		Exon pExons[] 					  = GetExons();
		TreeSet<CountElement> vcJunctions = GetJunctions();
		
		// there are two possibilities for an exon to be a retained intron:
		// 1. The exon has an (annotated) internal junction
		// 2. The exon includes two other exons that are not overlapping
		
		// identify retained introns by 'known' internal junctions
		for(Exon ex : pExons)
		{			
			for(CountElement jun : vcJunctions)
			{
				if(IsNovelJunction(jun))
					continue;
				
				if(jun.m_nStart > ex.getGenomicStart() && jun.m_nEnd < ex.getGenomicStop())
					m_vcRetainedIntrons.add(ex);
			}
		}
		
		// identify retained introns by overlap of multiple non-overlapping exons
		for(Exon ex1 : pExons)
		{
			// skip previously identified retained introns
			if(m_vcRetainedIntrons.contains(ex1))
				continue;
			
			Exon overlappingExon = null;
			
			for(Exon ex2 : pExons)
			{
				if(m_vcRetainedIntrons.contains(ex2))
					continue;
				
				if(!ex1.equals(ex2) && ex1.overlaps(ex2) && (overlappingExon == null || !ex2.overlaps(overlappingExon)))
				{
					if(overlappingExon != null)
					{
						m_vcRetainedIntrons.add(ex1);
						break;
					}
					else
						overlappingExon = ex2;
				}
			}
		}
	}
	
	/**
	 *    Processes a gene, which includes retrieving a sorted array
	 *    of exons, determination of exon groups, mapping of junction
	 *    counts and exon coverages to isoforms, identification of
	 *    exons that are retained introns and sub-grouping of exons into
	 *    terminal exons, middle exons and introns.
	 */
	private void ProcessGeneStructure()
	{
		//#################################
		//    get exons and exon groups
		//#################################
		m_pExons = m_Gene.getSortedExons();
		
		for(int i=0; i<m_pExons.length; i++)
			m_pExons[i].setID(i);
		
		// get exon groups
		m_pExonGroups = m_Gene.computeOverlappingExonGroupsNormalOrder();
		
		MapJunctionsToIsoforms();
		MapExonsToIsoforms();
		IdentifyRetainedIntrons();
				
		//###########################################################
		//    prepare sets of arrays, introns and exon groups
		//###########################################################
		for(String strIsoform : m_Gene.getArrayOfGeneProductNames())
		{
			Exon[] pIsoformExons 		= m_mapExonsToIsoforms.get(strIsoform);
			
			// collect middle exons and introns
			for(int i=0; i<pIsoformExons.length; i++)
			{				
				if(pIsoformExons.length > 2 && i != 0 && i != pIsoformExons.length-1)
				{
					// add middle exon
					m_vcMiddleExons.add(pIsoformExons[i]);
						
					// left exon
					Exon exon1 = pIsoformExons[i];
						
					// right exon
					Exon exon2 = pIsoformExons[i+1];
						
					// add the intron
					Exon intron = new Exon(exon1.getGenomicStop()+1, exon2.getGenomicStart()-1);
					m_vcIntrons.add(intron);
				}
			}
			
			// collect terminal exons
			if(m_Gene.isPlusStrand())
			{
				m_vcFirstExons.add(pIsoformExons[0]);
				m_vcLastExons.add(pIsoformExons[pIsoformExons.length-1]);
			}
			else
			{
				m_vcFirstExons.add(pIsoformExons[pIsoformExons.length-1]);
				m_vcLastExons.add(pIsoformExons[0]);
			}
		}
	}

	/** Returns whether the specified junction is included in the gene annotation */
	public boolean IsNovelJunction(CountElement jun)
	{
		if(m_vcNovelJunctions.contains(jun))
			return true;
		else
			return false;
	}
	
	/** Returns whether the junction is a valid junctions (i.e. has enough coverage)*/
	public boolean IsInvalidJunction(CountElement jun)
	{
		if(m_vcInvalidJunctions.contains(jun))
			return true;
		else
			return false;
	}
	
	/** Returns whether the exon is a retained intron */
	public boolean IsRetainedIntron(Exon ex)
	{
		if(m_vcRetainedIntrons.contains(ex))
			return true;
		else
			return false;
	}
	
	/** Updates the exon groups after isoform selection change */
	private void UpdateExonGroups(SplicingWebApp app)
	{
		//###################################################################
		//   retrieve data for the current selections and gene information
		//###################################################################
		TreeSet<String> vcValidIsoforms = app.GetValidIsoforms();
		//###################################################################
		
		//###########################################################
		//    prepare sets of arrays, introns and exon groups
		//###########################################################
		
		TreeSet<Exon> vcMiddleExons = new TreeSet<Exon>();
		TreeSet<Exon> vcAllExons	= new TreeSet<Exon>();
		
		for(String strIsoform : vcValidIsoforms)
		{
			Exon[] pIsoformExons 		= m_mapExonsToIsoforms.get(strIsoform);

			for(int i=0; i<pIsoformExons.length; i++)
			{				
				if(pIsoformExons.length > 2 && i != 0 && i != pIsoformExons.length-1)
				{
					// add middle exon
					vcMiddleExons.add(pIsoformExons[i]);
				}
				
				vcAllExons.add(pIsoformExons[i]);
			}
		}
		
		m_pMiddleExonGroups = RecalculateExonGroups(vcMiddleExons);
		m_pExonGroups		= RecalculateExonGroups(vcAllExons);
	}
	
	/**
	 *    Class that is used as thread to fill the coverage arrays for the current gene.
	 *    If the gene is opened for the first time or if samples have been changed, the
	 *    data supplier retrieves a coverage array for the whole gene. The array is stored
	 *    if retained introns should be identified (because then intron coverage is needed).
	 *    The coverage is used to retrieve the coverage for each exon group and, if a gene
	 *    is loaded for the first time, also to retrieve the coverage of exons and (if required)
	 *    introns.
	 */
	private class ThreadFillCoverageArrays implements Runnable
	{
		SplicingWebApp m_app;
		String m_strSample;
		boolean m_bDetectIntronDetection;
			
		ThreadFillCoverageArrays(SplicingWebApp app, String strSample)
		{
			m_app						= app;
			m_strSample					= strSample;
			m_bDetectIntronDetection	= app.DetectIntronRetentionEvents();
		}
			
		@Override
		public void run()
		{
			boolean bFirstInitialization = false;
			
			// retrieve coverage for the whole gene
			double pCoverage[] = null;

			synchronized(m_mapGeneCoveragePerSample)
			{
				if(!m_mapGeneCoveragePerSample.containsKey(m_strSample))
				{
					bFirstInitialization = true;
					try
					{
						pCoverage = m_app.GetCoverageForRegion(m_Gene.getChromosome(), m_Gene.getStart(), m_Gene.getStop(), m_strSample, true);
					}
					catch (IOException e)
					{
						System.out.println("failed to retrieve coverage from bigwig file for sample: " + m_strSample);
						System.out.println("region: " + m_Gene.getChromosome() + ":" + m_Gene.getStart() + "-" + m_Gene.getStop());
						e.printStackTrace();
						return;
					}
					
					if(m_bDetectIntronDetection)
					{
						m_mapGeneCoveragePerSample.put(m_strSample, pCoverage);
					}
				}
				else
					pCoverage = m_mapGeneCoveragePerSample.get(m_strSample);
			}
			
			// obtain coverage for exon groups (must be redone if exon groups change)
			for(ExonGroup grp : m_pExonGroups)
			{
				int nLength = grp.getGenomicStopOfGroup() - grp.getGenomicStartOfGroup() + 1;
				double pCov[] = new double[nLength];
				
				for(int i=0; i<nLength; i++)
				{
					int nOffset = grp.getGenomicStartOfGroup() - m_Gene.getStart();
					pCov[i] = pCoverage[i + nOffset];
				}
				
				synchronized(m_mapCoverageToExonGroupsPerSample)
				{
					m_mapCoverageToExonGroupsPerSample.get(grp).put(m_strSample, pCov);
				}
			}
			
			// obtain coverage for middle exon groups (they might be slightly different to the normal exon groups)
			for(ExonGroup grp : m_pMiddleExonGroups)
			{
				int nLength = grp.getGenomicStopOfGroup() - grp.getGenomicStartOfGroup() + 1;
				double pCov[] = new double[nLength];
				
				for(int i=0; i<nLength; i++)
				{
					int nOffset = grp.getGenomicStartOfGroup() - m_Gene.getStart();
					pCov[i] = pCoverage[i + nOffset];
				}
				
				synchronized(m_mapCoverageToMiddleExonGroupsPerSample)
				{
					m_mapCoverageToMiddleExonGroupsPerSample.get(grp).put(m_strSample, pCov);
				}
			}
			
			// this only has to be done once
			if(bFirstInitialization)
			{
				// obtain coverage for first exons
				for(Exon ex : m_vcFirstExons)
				{
					int nLength = ex.getGenomicStop() - ex.getGenomicStart() +1;
					double pCov[] = new double[nLength];
					
					for(int i=0; i<nLength; i++)
					{
						int nOffset = ex.getGenomicStart() - m_Gene.getStart();
						pCov[i] = pCoverage[i + nOffset];
					}
					
					synchronized(m_mapCoverageToExonsPerSample)
					{
						m_mapCoverageToExonsPerSample.get(ex).put(m_strSample, pCov);
					}
					
					// add coverage mean for the current exon
					synchronized(m_mapRawExonCoveragePerSample)
					{
						m_mapRawExonCoveragePerSample.get(ex).put(m_strSample, StatUtils.mean(pCov));
					}
					
					// calculate coverage fraction for current exon
					double nBasesBelowThreshold = 0.0;
					for(double fVal : pCov)
					{
						if(fVal < m_app.GetMinimumCoveragePerBase())
							nBasesBelowThreshold += 1;
					}
					double fCoveredFraction = 1.0-(double)(nBasesBelowThreshold) / pCov.length;
					
					// add coverage fraction for current exon
					synchronized(m_mapExonCoveredBasesPerSample)
					{
						m_mapExonCoveredBasesPerSample.get(ex).put(m_strSample, fCoveredFraction);
					}
				}
				
				// obtain coverage for last exons
				for(Exon ex : m_vcLastExons)
				{
					int nLength = ex.getGenomicStop() - ex.getGenomicStart() +1;
					double pCov[] = new double[nLength];
					
					for(int i=0; i<nLength; i++)
					{
						int nOffset = ex.getGenomicStart() - m_Gene.getStart();
						pCov[i] = pCoverage[i + nOffset];
					}
					
					synchronized(m_mapCoverageToExonsPerSample)
					{
						m_mapCoverageToExonsPerSample.get(ex).put(m_strSample, pCov);
					}
						
				 	// add coverage mean for the current exon
					synchronized(m_mapRawExonCoveragePerSample)
					{
						m_mapRawExonCoveragePerSample.get(ex).put(m_strSample, StatUtils.mean(pCov));
					}
					
					// calculate coverage fraction for current exon
					double nBasesBelowThreshold = 0.0;
					for(double fVal : pCov)
					{
						if(fVal < m_app.GetMinimumCoveragePerBase())
							nBasesBelowThreshold += 1;
					}
					double fCoveredFraction = 1.0-(double)(nBasesBelowThreshold) / pCov.length;
					
					// add coverage fraction for current exon
					synchronized(m_mapExonCoveredBasesPerSample)
					{
						m_mapExonCoveredBasesPerSample.get(ex).put(m_strSample, fCoveredFraction);
					}
				}
				
				// obtain coverage for middle exons
				for(Exon ex : m_vcMiddleExons)
				{
					int nLength = ex.getGenomicStop() - ex.getGenomicStart() +1;
					double pCov[] = new double[nLength];
					
					for(int i=0; i<nLength; i++)
					{
						int nOffset = ex.getGenomicStart() - m_Gene.getStart();
						pCov[i] = pCoverage[i + nOffset];
					}
					
					// get coverage array
					synchronized(m_mapCoverageToExonsPerSample)
					{
						m_mapCoverageToExonsPerSample.get(ex).put(m_strSample, pCov);
					}
					
					// add coverage mean for the current exon
					synchronized(m_mapRawExonCoveragePerSample)
					{
						m_mapRawExonCoveragePerSample.get(ex).put(m_strSample, StatUtils.mean(pCov));
					}
					
					// calculate coverage fraction for current exon
					double nBasesBelowThreshold = 0.0;
					for(double fVal : pCov)
					{
						if(fVal < m_app.GetMinimumCoveragePerBase())
							nBasesBelowThreshold += 1;
					}
					double fCoveredFraction = 1.0-(double)(nBasesBelowThreshold) / pCov.length;
					
					// add coverage fraction for current exon
					synchronized(m_mapExonCoveredBasesPerSample)
					{
						m_mapExonCoveredBasesPerSample.get(ex).put(m_strSample, fCoveredFraction);
					}
				}
				
				if(m_bDetectIntronDetection)
				{
					// obtain coverage for introns
					for(Exon ex : m_vcIntrons)
					{
						int nLength = ex.getGenomicStop() - ex.getGenomicStart() +1;
						double pCov[] = new double[nLength];
						
						for(int i=0; i<nLength; i++)
						{
							int nOffset = ex.getGenomicStart() - m_Gene.getStart();
							pCov[i] = pCoverage[i + nOffset];
						}
						
						synchronized(m_mapCoverageToIntronsPerSample)
						{
							m_mapCoverageToIntronsPerSample.get(ex).put(m_strSample, pCov);
						}
					}
				}
			}

			// calculate 'real' exon coverage
			for(ExonGroup grp : m_pExonGroups)
			{
				CalculateNormalExonCoverage(m_app, grp, m_strSample);
			}
			
			pCoverage = null;
		}
	}
	
	/**
	 *    Invoked from the ThreadFillCoverageArrays class. This function calculate the 
	 *    'normal' coverage for each exon of a given exon group for a given sample.
	 *    
	 *    First 'raw' coverage values are mapped to each exon. Then exons are divided
	 *    into exons that are indisinguishable from other exons (because they differ only
	 *    by a few bases, which is insufficient to assign good coverage estimates to both
	 *    of them) and exons that differ by a sufficient number of bases. The threshold is
	 *    set to 50 bp.
	 *    
	 *    Afterwards coverage is distributed across all exons, with the exception that
	 *    indistinguishable exons receive the same coverage vale. For the distribution,
	 *    each position of the exon group is assigned an exon. Each position may only be
	 *    assigned a single exon and the assignment starts with the highest covered exon.
	 *    
	 *    Exon coverages are then recalculated based on the coverage at all positions that
	 *    were assigned to the exon (representing the unique coverage for each exon).
	 *    
	 *    The unique coverage is then used to calculate the 'final'(adjusted) exon coverage
	 *    that is calculated by subsequently removing the coverage of the exon with the
	 *    lowest unique coverage from the exon group coverage array.
	 */
	private void CalculateNormalExonCoverage(SplicingWebApp app, ExonGroup exonGroup, String strSample)
	{
		// get all exons of the target group & get coverage of all of them for the selected sample
		TreeMap<Exon, Double> mapCoverageToExons = new TreeMap<Exon, Double>();
		Vector<Exon> vcExons = new Vector<Exon>();
		
		synchronized(m_mapRawExonCoveragePerSample)
		{
			for(Exon ex : exonGroup.getExons())
			{
				vcExons.add(ex);
				mapCoverageToExons.put(ex, m_mapRawExonCoveragePerSample.get(ex).get(strSample));
			}
		}

		//###################################################
		//            group too similar exons
		//    use a minimum of 50 bp to distinguish exons
		//###################################################
		TreeMap<Exon, TreeSet<Exon>> mapGroupedExons = GroupSimilarExons(vcExons, 50);
		
		//###############################################################################
		//    Also merge exons with similar coverage that are completely overlapping.
		//###############################################################################
		MergeSimilarExpressedContainedExons(mapGroupedExons, mapCoverageToExons);

		//##############################################
		//    create list of indistinguishable exons
		//##############################################
		synchronized(m_mapIndistinguishableExonPerSample)
		{
			for(Exon ex : mapGroupedExons.keySet())
			{
				if(m_mapIndistinguishableExonPerSample.containsKey(ex))
				{
					TreeMap<String, TreeSet<Exon>> vcTmp = m_mapIndistinguishableExonPerSample.get(ex);
					vcTmp.put(strSample, mapGroupedExons.get(ex));
				}
				else
				{
					TreeMap<String, TreeSet<Exon>> vcTmp = new TreeMap<String, TreeSet<Exon>>();
					vcTmp.put(strSample, mapGroupedExons.get(ex));
					m_mapIndistinguishableExonPerSample.put(ex, vcTmp);
				}
			}
		}

		//###########################################################################
		//    Track which base positions are claimed by which exon, starting with
		//    the highest expressed exon.
		//###########################################################################
		TreeMap<Integer, TreeSet<Exon>> mapCoveragePositions = DefineExonSpecificCoverageRegions(mapCoverageToExons, exonGroup.getGenomicStartOfGroup());
		
		double pCoverage[] = null;
		synchronized(m_mapCoverageToExonGroupsPerSample)
		{
			pCoverage = m_mapCoverageToExonGroupsPerSample.get(exonGroup).get(strSample);
		}
		
		if(pCoverage == null)
		{
			if(!m_mapCoverageToExonGroupsPerSample.containsKey(exonGroup))
			{
				System.out.println(m_Gene);
				System.out.println("ERROR: missing coverage array for group: " + exonGroup);
			}
			else
			{
				if(m_mapCoverageToExonGroupsPerSample.get(exonGroup).keySet().contains(strSample) == false)
				{
					System.out.println(m_Gene);
					System.out.println("ERROR: missing coverage array for sample " + strSample + " in group: " + exonGroup);
				
					for(String strSample2 : m_mapCoverageToExonGroupsPerSample.get(exonGroup).keySet())
					{
						System.out.println(strSample2);
						System.out.println(strSample.equals(strSample2));
						
						if(strSample.equals(strSample2))
							System.out.println(m_mapCoverageToExonGroupsPerSample.get(exonGroup).keySet().contains(strSample));
					}
				}
				else
				{
					if(m_mapCoverageToExonGroupsPerSample.get(exonGroup).get(strSample) == null)
						System.out.println("coverage array is empty! " + strSample + " "  + exonGroup);
				}
			}
		}

		//###############################################################################
		//    Recalculate exon coverages. Already used exon positions will be ignored
		//###############################################################################
		TreeMap<Exon, Double> mapUniqueCoverage = CalculateUniqueCoveragePerExon(mapCoveragePositions, mapGroupedExons, exonGroup.getGenomicStartOfGroup(), pCoverage);
		
		//###############################################################################
		//    Now that we have defined 'unique' coverage regions per exon, we can use
		//    these regions to calculate the final exon coverages.
		//###############################################################################
		TreeMap<Exon, Double> mapFinalExonCoverage = CalculateFinalCoveragePerExon(mapUniqueCoverage, exonGroup.getGenomicStartOfGroup(), pCoverage);
		
		synchronized(m_mapAdjustedExonCoveragePerSample)
		{
			for(Exon ex : mapFinalExonCoverage.keySet())
			{
				if(m_mapAdjustedExonCoveragePerSample.containsKey(ex))
				{
					TreeMap<String, Double> mapCoveragePerSample = m_mapAdjustedExonCoveragePerSample.get(ex);
					mapCoveragePerSample.put(strSample, mapFinalExonCoverage.get(ex));
				}
				else
				{
					TreeMap<String, Double> mapCoveragePerSample = new TreeMap<String, Double>();
					mapCoveragePerSample.put(strSample, mapFinalExonCoverage.get(ex));
					m_mapAdjustedExonCoveragePerSample.put(ex, mapCoveragePerSample);
				}
			}
		}
	}

	/** Recalculates the exon groups based on a given set of exon */
	public ExonGroup[] RecalculateExonGroups(TreeSet<Exon> vcExons)
	{
		int nGrpID = 0;
		TreeSet<ExonGroup> vcExonGroups = new TreeSet<ExonGroup>();
		for(Exon exon : vcExons)
		{
			boolean bAdded = false;
			for(ExonGroup grp : vcExonGroups)
			{
				if(grp.groupContainsExon(exon.getExonID()))
				{
					bAdded = true;
					break;
				}
				
				if(grp.exonIntersectsWithExonsInGroup(exon))
				{					
					grp.addExon(exon.getExonID(), exon.getExonID(), exon);
					bAdded = true;
					break;
				}
			}
			
			if(!bAdded)
			{
				ExonGroup grp = new ExonGroup(nGrpID, m_Gene.isPlusStrand());
				grp.addExon(exon.getExonID(), exon.getExonID(), exon);
				vcExonGroups.add(grp);
				nGrpID += 1;
			}
		}
		
		// the previously associated exon group IDs are incorrect, thus modify them
		int nGroupID = 0;
		for(ExonGroup grp : vcExonGroups)
		{
			grp.setGroupID(nGroupID);
			nGroupID += 1;
		}
		
		ExonGroup pExonGroups[] = new ExonGroup[vcExonGroups.size()];
		
		int i=0;
		for(ExonGroup grp : vcExonGroups)
		{
			pExonGroups[i] = grp;
			i++;
		}

		return pExonGroups;
	}
	
	/**
	 *    Prepares coverage containers and retrieves coverage arrays for
	 *    all samples from big wig files. Each sample is processed in its
	 *    own thread.
	 */
	private boolean FillCoverageArraysThreaded(SplicingWebApp app)
	{
		ProjectModel projectModel = app.GetProjectModel();
		
		//###################################################################
		//   retrieve data for the current selections and gene information
		//###################################################################
		String strSelectedConditionType 	= app.GetSelectedConditionType();
		TreeSet<String> vcSelectedSamples 	= app.GetSelectedSamples();
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = projectModel.GetSelectedSamplesPerCondition(strSelectedConditionType, vcSelectedSamples);
		//###################################################################
		
		TreeMap<String, String> mapBigWigFiles 	= projectModel.GetBigWigFilesForSamples();
		
		if(mapBigWigFiles == null)
		{
			Messagebox.show("ERROR: no valid bigwig or bam files detected");
			return false;
		}
		
		// only get the raw coverage once for each exon
		boolean bFirstInitialization = false;
		
		if(m_mapGeneCoveragePerSample.isEmpty())
			bFirstInitialization = true;
		
		//###################################################################
		//  prepare coverage container for each exon group, exon and intron
		//###################################################################
		for(ExonGroup grp : m_pExonGroups)
		{
			TreeMap<String, double[]> mapTmp = new TreeMap<String, double[]>();
			m_mapCoverageToExonGroupsPerSample.put(grp, mapTmp);
		}
		
		for(ExonGroup grp : m_pMiddleExonGroups)
		{
			TreeMap<String, double[]> mapTmp = new TreeMap<String, double[]>();
			m_mapCoverageToMiddleExonGroupsPerSample.put(grp, mapTmp);
		}
		
		for(Exon ex : m_vcFirstExons)
		{
			if(bFirstInitialization)
			{
				TreeMap<String, double[]> mapTmp = new TreeMap<String, double[]>();
				m_mapCoverageToExonsPerSample.put(ex, mapTmp);
				
				TreeMap<String, Double> mapTmp2 = new TreeMap<String, Double>();
				m_mapRawExonCoveragePerSample.put(ex, mapTmp2);
				
				mapTmp2 = new TreeMap<String, Double>();
				m_mapExonCoveredBasesPerSample.put(ex, mapTmp2);
			}
			
			TreeMap<String, Double> mapTmp2 = new TreeMap<String, Double>();
			m_mapAdjustedExonCoveragePerSample.put(ex, mapTmp2);
		}
		
		for(Exon ex : m_vcLastExons)
		{
			if(bFirstInitialization)
			{
				TreeMap<String, double[]> mapTmp = new TreeMap<String, double[]>();
				m_mapCoverageToExonsPerSample.put(ex, mapTmp);
				
				TreeMap<String, Double> mapTmp2 = new TreeMap<String, Double>();
				m_mapRawExonCoveragePerSample.put(ex, mapTmp2);
				
				mapTmp2 = new TreeMap<String, Double>();
				m_mapExonCoveredBasesPerSample.put(ex, mapTmp2);
			}
			
			TreeMap<String, Double> mapTmp2 = new TreeMap<String, Double>();
			m_mapAdjustedExonCoveragePerSample.put(ex, mapTmp2);
		}
		
		for(Exon ex : m_vcMiddleExons)
		{
			if(bFirstInitialization)
			{
				TreeMap<String, double[]> mapTmp = new TreeMap<String, double[]>();
				m_mapCoverageToExonsPerSample.put(ex, mapTmp);
			
				TreeMap<String, Double> mapTmp2 = new TreeMap<String, Double>();
				m_mapRawExonCoveragePerSample.put(ex, mapTmp2);
				
				mapTmp2 = new TreeMap<String, Double>();
				m_mapExonCoveredBasesPerSample.put(ex, mapTmp2);
			}
			
			TreeMap<String, Double> mapTmp2 = new TreeMap<String, Double>();
			m_mapAdjustedExonCoveragePerSample.put(ex, mapTmp2);
		}
		
		if(bFirstInitialization)
		{
			for(Exon ex : m_vcIntrons)
			{
				TreeMap<String, double[]> mapTmp = new TreeMap<String, double[]>();
				m_mapCoverageToIntronsPerSample.put(ex, mapTmp);
			}
		}
		
		m_mapIndistinguishableExonPerSample = new TreeMap<Exon, TreeMap<String, TreeSet<Exon>>>();
		
		//###################################################################
		//    load coverage data from bigwig files and prepare coverage
		//###################################################################
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			ExecutorService executor = Executors.newFixedThreadPool(app.GetMaximumThreads());
			
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				// skip samples that are currently not selected
				if(!vcSelectedSamples.contains(strSample))
					continue;
				
				Runnable r = new ThreadFillCoverageArrays(app, strSample);
				executor.execute(r);
			}
			
			executor.shutdown();
			while(!executor.isTerminated()) {}
		}

		return true;
	}
	
	/**
	 * This function identifies conditions that have insufficient coverage
	 * for all exons of the current gene. These conditions won't be tested
 	 * for alternative splicing (regarding the current gene). 
	 */
	private void IdentifyInvalidConditions(SplicingWebApp app)
	{
		ProjectModel projectModel = app.GetProjectModel();
		
		//###################################################################
		//   retrieve data for the current selections and gene information
		//###################################################################
		String strSelectedConditionType 	= app.GetSelectedConditionType();
		TreeSet<String> vcSelectedSamples 	= app.GetSelectedSamples();
		int nMinCovPerBase 					= app.GetMinimumCoveragePerBase();
		TreeMap<String, TreeSet<String>> vcSamplesAndConditions = projectModel.GetSamplesPerCondition(strSelectedConditionType);
		//###################################################################
		
		// identify conditions with insufficient coverage for all of its exons
		for(String strCondition : vcSamplesAndConditions.keySet())
		{
			double fMaxCoverage = 0.0;
			
			for(ExonGroup grp : m_mapCoverageToExonGroupsPerSample.keySet())
			{
				Vector<Double> vcValues = new Vector<Double>();
				
				for(String strSample : vcSamplesAndConditions.get(strCondition))
				{
					if(!vcSelectedSamples.contains(strSample))
						continue;

					if(m_mapCoverageToExonGroupsPerSample.get(grp).keySet().size() > 0)
					{
						if(m_mapCoverageToExonGroupsPerSample.get(grp).containsKey(strSample))
							vcValues.add(StatUtils.mean(m_mapCoverageToExonGroupsPerSample.get(grp).get(strSample)));
						else
							vcValues.add(0.0);
					}
					else
					{
						vcValues.add(0.0);
					}
				}
				
				double[] pValues =  new double[vcValues.size()];
				int nIdx = 0;
				for(double val : vcValues)
				{
					pValues[nIdx] = val;
					nIdx++;
				}
				
				fMaxCoverage = Math.max(fMaxCoverage, StatUtils.percentile(pValues, 50.0));
			}
			
			if(fMaxCoverage < nMinCovPerBase)
			{
				m_vcInvalidConditions.add(strCondition);
			}
		}
	}
	
	/**
	 *    Similar to above, but generates a list of conditions that don't
	 *    have any usable samples. Requires that IdentifyInvalidConditions()
	 *    has been invoked first.
	 */
	private void IdentifyConditionsWithNoValidSamples(SplicingWebApp app)
	{
		ProjectModel projectModel = app.GetProjectModel();
		
		//###################################################################
		//   retrieve data for the current selections and gene information
		//###################################################################
		String strSelectedConditionType 	= app.GetSelectedConditionType();
		TreeSet<String> vcSelectedSamples 	= app.GetSelectedSamples();
		TreeMap<String, TreeSet<String>> vcSamplesAndConditions = projectModel.GetSamplesPerCondition(strSelectedConditionType);
		//###################################################################
				
		//##############################################################
		//   get set of conditions for which no samples were selected
		//##############################################################
		for(String strCondition : vcSamplesAndConditions.keySet())
		{	
			int nValidSamples = 0;
			for(String strSample : vcSamplesAndConditions.get(strCondition))
			{
				if(vcSelectedSamples.contains(strSample))
					nValidSamples += 1;
			}

			if(nValidSamples == 0)
			{
				m_vcInvalidConditions.add(strCondition);
			}
		}
	}
	
	/**
	 *    Wrapper function with calls to recalculate the exon groups,
	 *    retrieve coverage arrays and identify invalid conditions.
	 */
	public boolean RetrieveCoverageData(SplicingWebApp app)
	{
		// update exon groups for isoform selection
		UpdateExonGroups(app);
		
		// identify conditions for which no samples were selected
		IdentifyConditionsWithNoValidSamples(app);
		
		// fill coverage arrays
		FillCoverageArraysThreaded(app);
		
		// identify conditions that should not be tested due to insufficient gene coverage
		IdentifyInvalidConditions(app);
		
		// success
		return true;
	}

	/**
	 *    Retrieves all junction counts overlapping the current gene locus.
	 *    Also collects and flags novel junctions, and identifies retained introns.
	 */
	public boolean RetrieveJunctionCounts(SplicingWebApp app)
	{
		ProjectModel projectModel = app.GetProjectModel();
		
		try
		{
			m_mapJunctionCounts = projectModel.GetJunctionsForRangeAsCountElements(m_Gene.getChromosome(), m_Gene.getStart(), m_Gene.getStop());
		}
		catch (IOException e)
		{
			System.out.println("failed to obtain junction counts for gene: " + m_Gene);
			e.printStackTrace();
			return false;
		}
		
		// identify novel junctions
		for(CountElement jun : m_mapJunctionCounts.keySet())
		{
			boolean bNovel = true;
			
			for(String strIsoform : GetIsoformNames())
			{
				TreeSet<CountElement> vcJunctions = GetJunctionsForIsoform(strIsoform);
				if(vcJunctions.contains(jun))
				{
					bNovel = false;
					break;
				}
			}
			
			if(bNovel) m_vcNovelJunctions.add(jun);
		}
		
		// flag known junctions
		for(CountElement jun : m_mapJunctionCounts.keySet())
		{
			if(!m_vcNovelJunctions.contains(jun))
				jun.m_bKnown = true;
		}
		
		// identify retained introns now that we have the junction counts
		IdentifyRetainedIntrons();

		return true;
	}
	
	/** Identifies lowly covered (=invalid) junctions */
	public void IdentifyInvalidJunctions(SplicingWebApp app, boolean bDebug)
	{
		ProjectModel projectModel 			= app.GetProjectModel();
		int nMinJunctionReads 				= app.GetMinimumJunctionReads();
		String strSelectedConditionType 	= app.GetSelectedConditionType();
		
		// clear old invalid junctions
		m_vcInvalidJunctions.clear();
		
		// get samples per Condition
		TreeMap<String, TreeSet<String>> mapSamplesToConditions =  projectModel.GetSamplesPerCondition(strSelectedConditionType);
		
		// identify junctions with insufficient coverage
		for(CountElement jun : m_mapJunctionCounts.keySet())
		{
			TreeMap<String, double[]> mapCoverageToConditions = new TreeMap<String, double[]>();
			for(String strCondition : mapSamplesToConditions.keySet())
			{
				int nSamples = mapSamplesToConditions.get(strCondition).size();
				mapCoverageToConditions.put(strCondition, new double[nSamples]);
			}
		
			int nInvalid = 0;
			for(String strCondition : mapSamplesToConditions.keySet())
			{							
				int nSampleIdx = 0;
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					if(m_mapJunctionCounts.get(jun).containsKey(strSample))
					{
						mapCoverageToConditions.get(strCondition)[nSampleIdx] = m_mapJunctionCounts.get(jun).get(strSample);									
					}
					else
					{
						mapCoverageToConditions.get(strCondition)[nSampleIdx] = 0;
					}
					nSampleIdx++;
				}
				double fMean = StatUtils.mean(mapCoverageToConditions.get(strCondition));
				
				if(fMean < nMinJunctionReads)
				{
					nInvalid++;
				}
			}
			
			if(nInvalid == mapSamplesToConditions.keySet().size())
			{
				m_vcInvalidJunctions.add(jun);
			}
		}
		
		if(bDebug)
		{
			System.out.println("identified " + m_vcInvalidJunctions.size() + " junctions with insufficient coverage");
			for(CountElement jun : m_vcInvalidJunctions)
			{
				System.out.println("invalid junction: " + jun.toSimpleString());
			}
		}
	}

	/** Returns the raw coverage for a given exon and sample */
	double GetRawCoverageForExonAndSample(Exon ex, String strSample)
	{
		if(m_mapRawExonCoveragePerSample.containsKey(ex))
		{
			if(m_mapRawExonCoveragePerSample.get(ex).containsKey(strSample))
			{
				return m_mapRawExonCoveragePerSample.get(ex).get(strSample);
			}
		}
			
		return 0.0;
	}
	
	/** Retrieves the coverage fraction (number of bases covered of an exon in %) for a given exon and sample */
	double GetCoverageFractionForExonAndSample(Exon ex, String strSample)
	{
		if(m_mapExonCoveredBasesPerSample.containsKey(ex))
		{
			if(m_mapExonCoveredBasesPerSample.get(ex).containsKey(strSample))
			{
				return m_mapExonCoveredBasesPerSample.get(ex).get(strSample);
			}
		}
			
		return 0.0;
	}

	/** Returns the list of indistinguishable exons */
	TreeMap<Exon, TreeMap<String, TreeSet<Exon>>> GetIndistinguishableExons()
	{
		return m_mapIndistinguishableExonPerSample;
	}
	
	/** Groups exons of an exon group if they differ by less than nMinUniqueBases */
	public TreeMap<Exon, TreeSet<Exon>> GroupSimilarExons(Vector<Exon> vcExons, int nMinUniqueBases)
	{
		TreeMap<Exon, TreeSet<Exon>> mapGroupedExons = new TreeMap<Exon, TreeSet<Exon>>(); 
		for(Exon ex : vcExons)
		{
			boolean bNewExon = true;
			// check whether the current exon has sufficient unique bases
			for(Exon exIncluded : mapGroupedExons.keySet())
			{
				int nExonStart	= ex.getGenomicStart();
				int nExonEnd	= ex.getGenomicStop();
				
				int nMetaExonStart 	= exIncluded.getGenomicStart();
				int nMetaExonEnd 	= exIncluded.getGenomicStop();
				
				// calculate number of differing bases
				int nDiff = Math.abs(nMetaExonStart - nExonStart) + Math.abs(nExonEnd - nMetaExonEnd);
				
				if(nDiff < nMinUniqueBases)
				{
					mapGroupedExons.get(exIncluded).add(ex);
					bNewExon = false;
					break;
				}
			}
			
			if(bNewExon)
			{
				TreeSet<Exon> exGrp = new TreeSet<Exon>();
				exGrp.add(ex);
				mapGroupedExons.put(ex, exGrp);
			}
		}
		
		return mapGroupedExons;
	}

	/** Merges all exons with similar expression where one exon is completely contained in the other */
	public void MergeSimilarExpressedContainedExons(TreeMap<Exon, TreeSet<Exon>> mapGroupedExons, TreeMap<Exon, Double> mapCoveragePerExon)
	{
		boolean bMerged = true;
		while(bMerged)
		{
			bMerged = false;
	
			outer_loop:
			for(Exon ex1 : mapGroupedExons.keySet())
			{
				for(Exon ex2 : mapGroupedExons.keySet())
				{
					Double fCov1 = 0.0;
					Double fCov2 = 0.0;
					
					for(Map.Entry<Exon, Double> e : mapCoveragePerExon.entrySet())
					{
						if(e.getValue().equals(ex1))
							fCov1 = e.getValue();
						
						if(e.getValue().equals(ex2))
							fCov2 = e.getValue();
					}

					// different exons with similar coverage where one exon contains the other
					if(!ex1.equals(ex2) && ex1.intersects(ex2) && Math.abs(fCov1 / (fCov1 + fCov2) -0.5) < 0.03)
					{
						TreeSet<Exon> vcExons1 = mapGroupedExons.get(ex1);
						TreeSet<Exon> vcExons2 = mapGroupedExons.get(ex2);
						
						vcExons1.addAll(vcExons2);
						
						mapGroupedExons.remove(ex1);
						mapGroupedExons.remove(ex2);
						
						if(ex1.getGenomicLength() > ex2.getGenomicLength())
							mapGroupedExons.put(ex1, vcExons1);
						else
							mapGroupedExons.put(ex2, vcExons1);
						
						bMerged = true;
						
						break outer_loop;
					}
				}
			}
		}
	}

	/** Assigns exons to exon positions  */
	public TreeMap<Integer, TreeSet<Exon>> DefineExonSpecificCoverageRegions(TreeMap<Exon, Double> mapCoverageToExon, int nExonGrpStart)
	{
		// track if a base position has been used by a previous exon already
		TreeMap<Integer, TreeSet<Exon>> mapCoveragePositions = new TreeMap<Integer, TreeSet<Exon>>();

		// define exon specific regions in the coverage array
		while(mapCoverageToExon.size() > 0)
		{
			double fMaxCoverage = 0;
			
			// get highest coverage
			for(Map.Entry<Exon, Double> e : mapCoverageToExon.entrySet())
			{
				if(fMaxCoverage < e.getValue())
				{
					fMaxCoverage 	= e.getValue();
				}
			}
			
			// get all exons with that coverage
			TreeSet<Exon> vcHighestCovExons = new TreeSet<Exon>();
			for(Map.Entry<Exon, Double> e : mapCoverageToExon.entrySet())
			{
				if(e.getValue() == fMaxCoverage)
				{
					vcHighestCovExons.add(e.getKey());
				}
			}
			
			// assign coverage positions to the highest covered exon(s)
			// it should be sufficient to use the first exon position
			Exon firstExon = vcHighestCovExons.first();
			
			int nExStart = firstExon.getGenomicStart();
			int nExStop  = firstExon.getGenomicStop();
			
			for(int i = nExStart; i<=nExStop; i++)
			{
				int nPos = i - nExonGrpStart;
				
				if(mapCoveragePositions.containsKey(nPos))
				{
					// do nothing, position claimed by higher covered exon
				}
				else
				{
					TreeSet<Exon> vcTmp  = new TreeSet<Exon>();
					// add all exons to the position
					for(Exon ex : vcHighestCovExons)
						vcTmp.add(ex);
					mapCoveragePositions.put(nPos, vcTmp);
				} 
			}

			// remove exons from list
			for(Exon ex : vcHighestCovExons)
			{
				mapCoverageToExon.remove(ex);
			}			
		}
		
		return mapCoveragePositions;
	}

	/**
	 *    Calculates the 'unique' coverage per exon. The first parameter is a map
	 *    where each position of the coverage array (last parameter) is assigned
	 *    to exactly one exon. For each exon unambigious exon (key of the grouped
	 *    exon map) the unique coverage is calculated.
	 *    
	 *    PLEASE BE AWARE that, in this context, 'unique' coverage may also
	 *    originate from regions that are shared among multiple exons but were
	 *    assigned to a specifc exon due to it's superior overall coverage.
	*/
	public TreeMap<Exon, Double> CalculateUniqueCoveragePerExon(TreeMap<Integer, TreeSet<Exon>> mapCoveragePositions, TreeMap<Exon, TreeSet<Exon>> mapGroupedExons, int nExonGrpStart, double[] pCoverage)
	{
		// for each exon, calculate the unique coverage
		TreeMap<Exon, Double> mapUniqueCoverage = new TreeMap<Exon, Double>();
		
		for(Exon ex : mapGroupedExons.keySet())
		{		
			int nExStart = ex.getGenomicStart();
			int nExEnd   = ex.getGenomicStop();
			for(Exon ex2 : mapGroupedExons.get(ex))
			{
				if(ex2.getGenomicStart() < nExStart)
					nExStart = ex2.getGenomicStart();
				
				if(ex2.getGenomicStop() > nExEnd)
					nExEnd = ex2.getGenomicStop();
			}

			double fCoverage = 0.0;
			int nUniqueBases = 0;
									
			for(int i = nExStart; i<=nExEnd; i++)
			{
				int nPos = i - nExonGrpStart;								

				// only use positions previously assigned to the exon 
				if(mapCoveragePositions.containsKey(nPos))
				{
					boolean bValid = false;
					for(Exon ex2 : mapGroupedExons.get(ex))
					{
						if(mapCoveragePositions.get(nPos).contains(ex2))
						{
							bValid = true;
							break;
						}
					}
					
					if(bValid)
					{
						//TODO
						try
						{
							fCoverage += pCoverage[nPos];
						}
						catch(Exception e)
						{
							System.out.println("pos: " + nPos);
							System.out.println(Arrays.toString(pCoverage));
						}
						nUniqueBases += 1;
					}
				}
				else
				{
					// no coverage available
				}
			}

			if(nUniqueBases > 0)
			{
				fCoverage /= nUniqueBases;
				
				// add exon
				mapUniqueCoverage.put(ex, fCoverage);
			}
			else
			{
				mapUniqueCoverage.put(ex, 0.0);
			}
		}
		
		return mapUniqueCoverage;
	}
	
	/**
	 *    Calculates the final coverage per exon based on the 'unique' exon coverage.
	 *    At first, all exons that have 'unique' coverage are assigned to the respective
	 *    base positions in the coverage array (see above for the definition of 'unique'
	 *    coverage). The coverage array starts with the first base of the exon group.
	 */
	public TreeMap<Exon, Double> CalculateFinalCoveragePerExon(TreeMap<Exon, Double> mapUniqueCoverage, int nExonGroupStart, double pCoverage[])
	{
		TreeMap<Exon, Double> vcFinalExonCoverage = new TreeMap<Exon, Double>();
		TreeMap<Exon, Double> vcTmpCoverage = new TreeMap<Exon, Double>();
		
		// assign exons to positions in the coverage array
		TreeMap<Integer, TreeSet<Exon>> mapExonsToPositions = new TreeMap<Integer, TreeSet<Exon>>(); 
		for(Exon ex : mapUniqueCoverage.keySet())
		{
			int nExStart = ex.getGenomicStart();
			int nExEnd   = ex.getGenomicStop();
			
			for(int i = nExStart; i<=nExEnd; i++)
			{
				int nPos = i - nExonGroupStart;
				
				if(mapExonsToPositions.containsKey(nPos))
				{
					mapExonsToPositions.get(nPos).add(ex);
				}
				else
				{
					TreeSet<Exon> vcTmp = new TreeSet<Exon>();
					vcTmp.add(ex);
					mapExonsToPositions.put(nPos, vcTmp);
				}
			}							
		}
		
		// calculate coverage
		for(int nPos : mapExonsToPositions.keySet())
		{
			Vector<Double> vcCoverages = new Vector<Double>();
			Vector<Exon>   vcExons = new Vector<Exon>();
			
			Vector<Double> vcCovTmp = new Vector<Double>();
			Vector<Exon>   vcExnTmp = new Vector<Exon>();
			
			for(Exon ex : mapExonsToPositions.get(nPos))
			{
				vcExnTmp.add(ex);
				vcCovTmp.add(mapUniqueCoverage.get(ex));
			}
			
			// add the lowest coverage first
			while(vcExnTmp.size() > 0)
			{
				int nLowestIdx = -1;
				Double nLowestCov = Double.MAX_VALUE;
				
				for(int i=0; i<vcCovTmp.size(); i++)
				{
					if(vcCovTmp.get(i) < nLowestCov)
					{
						nLowestIdx = i;
						nLowestCov = vcCovTmp.get(i); 
					}
				}
				
				vcCoverages.add(vcCovTmp.get(nLowestIdx));
				vcExons.add(vcExnTmp.get(nLowestIdx));
				
				vcExnTmp.remove(nLowestIdx);
				vcCovTmp.remove(nLowestIdx);
			}

			// get maximum coverage 
			double fMaxCoverage = 0;
			for(double fCov : vcCoverages)
			{
				fMaxCoverage = Math.max(fMaxCoverage, fCov);
			}
			
			for(int nExon=0; nExon<vcExons.size(); nExon++)
			{
				double fCoverageFraction = 0;
				
				Exon exCur = vcExons.get(nExon);
				
				if(fMaxCoverage == 0)
				{
					fCoverageFraction = 0.0;
				}
				else
				{
					if(nExon == 0)
					{									
						fCoverageFraction = pCoverage[nPos]*(vcCoverages.get(nExon)/fMaxCoverage);
					}
					else
					{
						fCoverageFraction = pCoverage[nPos]*(vcCoverages.get(nExon)-vcCoverages.get(nExon-1))/fMaxCoverage;
					}
				}

				if(vcTmpCoverage.containsKey(exCur))
				{
					double fValue = vcTmpCoverage.get(exCur);
					fValue += fCoverageFraction;
					vcTmpCoverage.put(exCur, fValue);
				}
				else
				{
					double fValue = fCoverageFraction;
					vcTmpCoverage.put(exCur, fValue);
				}
			}
		}

		for(Exon ex : vcTmpCoverage.keySet())
		{
			int nExLength = ex.getLength();
			double fCoverage = vcTmpCoverage.get(ex)/nExLength;
			
			vcFinalExonCoverage.put(ex, fCoverage);
		}
		
		return vcFinalExonCoverage;
	}
	
	/**
	 *    Different sets of ambigious (indistinguishable) exons may be created for each sample.
	 *    This function retrieves the set of ambigious exons that is most frequent across all
	 *    samples of a condition.
	 */
	public TreeMap<String, TreeSet<Exon>> GetMostFrequentSetOfAmbigiousExons(SplicingWebApp app, Exon exon)
	{
		//###################################################################
		//   retrieve data for the current selections and gene information
		//###################################################################
		String strSelectedConditionType 	= app.GetSelectedConditionType();
		TreeSet<String> vcSelectedSamples 	= app.GetSelectedSamples();
		ProjectModel projectModel 			= app.GetProjectModel();
		//###################################################################
				
		TreeMap<String, TreeSet<Exon>> mapAmbigiousExonsPerCondition = new TreeMap<String, TreeSet<Exon>>();
		
		TreeMap<String, TreeSet<String>> vcSamplesPerCondition = projectModel.GetSamplesPerCondition(strSelectedConditionType);
		
		for(String strCondition : vcSamplesPerCondition.keySet())
		{
			HashMap<TreeSet<Exon>, Integer> mapExonOccurrences = new HashMap<TreeSet<Exon>, Integer>(); 
			for(String strSample : vcSamplesPerCondition.get(strCondition))
			{
				// skip unselected samples
				if(!vcSelectedSamples.contains(strSample))
					continue;
				
				// proceed through all sets of indistinguishable exons
				for(Exon tarEx : m_mapIndistinguishableExonPerSample.keySet())
				{
					if(m_mapIndistinguishableExonPerSample.get(tarEx).containsKey(strSample))
					{
						// check whether the query exon is contained in the current set of indistinguishable exons
						if(m_mapIndistinguishableExonPerSample.get(tarEx).get(strSample).contains(exon))
						{
							TreeSet<Exon> vcExons = m_mapIndistinguishableExonPerSample.get(tarEx).get(strSample);
							
							// check if the set of exons has been added already
							boolean bSameSet = false;
							loop_check_occurence:
							for(TreeSet<Exon> vcCurExons : mapExonOccurrences.keySet())
							{
								bSameSet = true;
								
								// exon sets must have the same size
								if(vcExons.size() != vcCurExons.size())
								{
									bSameSet = false;
									continue loop_check_occurence;
								}
		
								// all exons must be the same
								for(Exon ex : vcExons)
								{
									if(!vcCurExons.contains(ex))
									{
										bSameSet = false;
										continue loop_check_occurence;
									}
								}
								
								// it's the same set, thus just increase the counter
								if(bSameSet)
								{
									mapExonOccurrences.put(vcCurExons, mapExonOccurrences.get(vcCurExons)+1);
									break;
								}
							}
		
							// it's a new set
							if(!bSameSet)
							{
								mapExonOccurrences.put(vcExons, 1);
							}
						}
					}
				}
			}
			
			// now identify the set which occurs most often
			int nMax = -1;
			TreeSet<Exon> vcTarget = new TreeSet<Exon>();
			for(Map.Entry<TreeSet<Exon>, Integer> e : mapExonOccurrences.entrySet())
			{
				if(nMax < e.getValue())
				{
					vcTarget = e.getKey();
					nMax = e.getValue();
				}
			}
			
			mapAmbigiousExonsPerCondition.put(strCondition, vcTarget);
		}
		
		return mapAmbigiousExonsPerCondition;
	}

	/**
	 *    Returns a exon coverage values for each sample fir a specified exon.
	 *    The function first calls GetMostFrequentSetOfAmbigiousExons() to identify
	 *    which exons it must query and then returns the respective 'final' (adjusted)
	 *    coverage value.
	 */
	public TreeMap<String, Double> GetAbsoluteCoveragePerSample(SplicingWebApp app, Exon exon)
	{
		//###################################################################
		//   retrieve data for the current selections and gene information
		//###################################################################
		String strSelectedConditionType 	= app.GetSelectedConditionType();
		TreeSet<String> vcSelectedSamples 	= app.GetSelectedSamples();
		ProjectModel projectModel 			= app.GetProjectModel();
		//###################################################################
		
		TreeMap<String, TreeSet<String>> vcSamplesPerCondition = projectModel.GetSamplesPerCondition(strSelectedConditionType);
		TreeMap<String, Double> mapAbsoluteExonCoveragePerSample = new TreeMap<String, Double>();
		
		// Get a majority vote on how to calculate the expression of this exon
		// The most frequent exon separation will be used
		TreeMap<String, Exon> mapTargetExonsPerSample = new TreeMap<String, Exon>();
		
		TreeMap<String, TreeSet<Exon>> vcAmbigiousExonSetsPerCondition = GetMostFrequentSetOfAmbigiousExons(app, exon);

		for(String strCondition : vcSamplesPerCondition.keySet())
		{		
			TreeSet<Exon> vcTarget = vcAmbigiousExonSetsPerCondition.get(strCondition);
			
			// get the target exon per sample
			for(String strSample : vcSamplesPerCondition.get(strCondition))
			{
				// use only selected samples
				if(!vcSelectedSamples.contains(strSample))
					continue;
				
				// get indistinguishable exons of the current exon
				for(Exon tarEx : m_mapIndistinguishableExonPerSample.keySet())
				{
					if(m_mapIndistinguishableExonPerSample.get(tarEx).containsKey(strSample))
					{
						TreeSet<Exon> vcExons = m_mapIndistinguishableExonPerSample.get(tarEx).get(strSample);
						
						boolean bIsSame = true;
						if(vcTarget.size() != vcExons.size())
							bIsSame = false;
						
						for(Exon ex : vcExons)
						{
							if(!vcTarget.contains(ex))
								bIsSame = false;
						}
						
						if(bIsSame)
							mapTargetExonsPerSample.put(strSample, tarEx);
					}
				}
				
			}
		}

		for(String strCondition : vcSamplesPerCondition.keySet())
		{
			for(String strSample : vcSamplesPerCondition.get(strCondition))
			{
				if(!vcSelectedSamples.contains(strSample))
					continue;
				
				if(!mapTargetExonsPerSample.containsKey(strSample))
					continue;
				
				Exon targetExon = mapTargetExonsPerSample.get(strSample);

				if(vcSamplesPerCondition.get(strCondition).contains(strSample))
				{
					if(m_mapAdjustedExonCoveragePerSample.containsKey(targetExon))
					{
						mapAbsoluteExonCoveragePerSample.put(strSample, m_mapAdjustedExonCoveragePerSample.get(targetExon).get(strSample));
					}
				}
			}
		}

		return mapAbsoluteExonCoveragePerSample;
	}

	/** Returns an array with all exons of the gene */
	public Exon[] GetExons()
	{
		return m_pExons;
	}
	
	/** Returns an array with all exon groups of the gene */
	public ExonGroup[] GetExonGroups()
	{
		return m_pExonGroups;
	}
	
	/** Returns an array of all 'middle' (non-terminal) exon groups */
	public ExonGroup[] GetMiddleExonGroups()
	{
		return m_pMiddleExonGroups;
	}
	
	/** Returns a list of all splice junctions of the gene */
	public TreeSet<CountElement> GetJunctions()
	{
		TreeSet<CountElement> vcJunctions = new TreeSet<CountElement>();
		for(CountElement jun : m_mapJunctionCounts.keySet())
			vcJunctions.add(jun);
		
		return vcJunctions;
	}
	
	/** Returns a list of all junctions for a specific isoform */
	public TreeSet<CountElement> GetJunctionsForIsoform(String strIsoform)
	{
		return m_mapJunctionsToIsoforms.get(strIsoform);
	}
	
	/** Returns the coverage counts for a specific junction */
	public TreeMap<String, Integer> GetCountsForJunction(CountElement jun)
	{
		if(m_mapJunctionCounts.containsKey(jun))
			return m_mapJunctionCounts.get(jun);
		
		return null;
	}
	
	/** 
	 *    Returns the coverage counts for a junction that is specified
	 *    by start and end position
	 */
	public TreeMap<String, Integer> GetCountsForJunction(int nStart, int nEnd)
	{
		for(CountElement jun : m_mapJunctionCounts.keySet())
		{
			if(jun.m_nStart == nStart && jun.m_nEnd == nEnd)
			{				
				return m_mapJunctionCounts.get(jun);
			}
		}
		
		return null;
	}

	/** Returns an array of exons for the specified isoform */
	public Exon[] GetExonsForIsoform(String strIsoform)
	{
		return m_mapExonsToIsoforms.get(strIsoform);
	}
	 
	/** Returns an array with all isoform names */
	public String[] GetIsoformNames()
	{
		return m_Gene.getArrayOfGeneProductNames();
	}
	
	/** Returns a list with all first (start) exons */
	public TreeSet<Exon> GetAllFirstExons()
	{
		return m_vcFirstExons;
	}
	
	/** Returns a list with all last (end) exons */
	public TreeSet<Exon> GetAllLastExons()
	{
		return m_vcLastExons;
	}
	
	/** Returns a list of all middle (non-terminal) exons */
	public TreeSet<Exon> GetMiddleExons()
	{
		return m_vcMiddleExons;
	}
	
	/** Returns a list of all introns */
	public TreeSet<Exon> GetIntrons()
	{
		return m_vcIntrons;
	}
	
	/** Returns the gene's strand */
	public char GetStrand()
	{
		return m_Gene.getStrandStringID();
	}
	
	/** Returns the name of the reference on which the gene is located*/
	public String GetReferenceName()
	{
		return m_Gene.getChromosome();
	}
	
	/** Returns the gene ID*/
	public String GetGeneID()
	{
		return m_Gene.getGeneID();
	}
	
	/** Returns the gene name (symbol) */
	public String GetGeneName()
	{
		return m_Gene.getGeneName();
	}
	
	/** Returns a list of all unique exons */
	public TreeSet<Exon> GetUniqueExons(TreeSet<String> vcIsoformsToConsider)
	{
		return m_Gene.getUniqueExons(vcIsoformsToConsider);
	}
	
	/** Returns a list of all unique junctions */
	public TreeSet<String> GetUniqueJunctions(TreeSet<String> vcIsoformsToConsider)
	{
		return m_Gene.getUniqueJunctions(vcIsoformsToConsider);
	}
	
	/** Returns the gene */
	public Gene GetGene()
	{
		return m_Gene;
	}
	
	/** Retrieves the coverage array for the whole gene for the specified sample */
	public double[] GetGeneCoverageArrayForSample(String strSample)
	{
		if(m_mapGeneCoveragePerSample.containsKey(strSample))
			return m_mapGeneCoveragePerSample.get(strSample);
		else
			return null;
	}
	
	/** Retrieves the coverage array for a specific exon group for the specified sample */
	public double[] GetCoverageArrayForExonGroupAndSample(SplicingWebApp app, ExonGroup grp, String strSample)
	{		
		if(m_mapCoverageToExonGroupsPerSample.containsKey(grp))
		{
			if(m_mapCoverageToExonGroupsPerSample.get(grp).keySet().contains(strSample))
			{
				return m_mapCoverageToExonGroupsPerSample.get(grp).get(strSample);
			}
		}
		
		// if the above did not find coverage for the specified exon group it might not be included in the list of prepared exon groups
		int nRegionStart = m_Gene.getStart();
		
		// groups must be with gene locus
		if(grp.getGenomicStopOfGroup() > m_Gene.getStop() || grp.getGenomicStartOfGroup() < m_Gene.getStart())
			return null;
		
		int nGrpStart  = grp.getGenomicStartOfGroup();
		int nGrpEnd	   = grp.getGenomicStopOfGroup();
		int nGrpLength = nGrpEnd - nGrpStart + 1;
		
		// group length must be positive
		if(nGrpLength < 1)
			return null;
		
		double pCoverage[] = GetGeneCoverageArrayForSample(strSample);
		if(pCoverage == null)
		{
			try
			{
				pCoverage = app.GetCoverageForRegion(GetReferenceName(), grp.getGenomicStartOfGroup(), grp.getGenomicStopOfGroup(), strSample, true);
				nRegionStart = grp.getGenomicStartOfGroup();
			}
			catch(Exception e)
			{
				System.out.println("ERROR: failed to retrieve coverage for region " + GetReferenceName() + ":" + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup() + " for sample " + strSample);
				e.printStackTrace();
				return null;
			}
		}
		double[] pGrpCoverage = new double[nGrpLength];
		
		int nArrayPos = 0;
		for(int i=nGrpStart; i<=nGrpEnd; i++)
		{
			pGrpCoverage[nArrayPos] = pCoverage[i - nRegionStart];
			nArrayPos++;
		}
		
		return pGrpCoverage;
	}
	
	/** Retrieves the coverage array for a specific 'middle' exon group for the specified sample */
	public double[] GetCoverageArrayForMiddleExonGroupAndSample(ExonGroup grp, String strSample)
	{		
		if(m_mapCoverageToMiddleExonGroupsPerSample.containsKey(grp))
		{
			if(m_mapCoverageToMiddleExonGroupsPerSample.get(grp).keySet().contains(strSample))
			{				
				return m_mapCoverageToMiddleExonGroupsPerSample.get(grp).get(strSample);
			}
		}
		
		return null;
	}
	
	/** Retrieves the coverage array for a specific exon for the specified sample */
	public double[] GetCoverageArrayForExonAndSample(Exon ex, String strSample)
	{
		if(m_mapCoverageToExonsPerSample.containsKey(ex))
		{
			if(m_mapCoverageToExonsPerSample.get(ex).keySet().contains(strSample))
			{				
				return m_mapCoverageToExonsPerSample.get(ex).get(strSample);
			}
		}
		
		return null;
	}
	
	/** Retrieves the coverage array for a specific intron for the specified sample */
	public double[] GetCoverageArrayForIntronAndSample(Exon intron, String strSample)
	{
		if(m_mapCoverageToIntronsPerSample.containsKey(intron))
		{
			if(m_mapCoverageToIntronsPerSample.get(intron).keySet().contains(strSample))
			{			
				return m_mapCoverageToIntronsPerSample.get(intron).get(strSample);
			}
		}
		
		return null;
	}
	
	/** Returns the gene start position */
	public int GetGeneStart()
	{
		return m_Gene.getStart();
	}

	/** Returns the gene end position */
	public int GetGeneEnd()
	{
		return m_Gene.getStop();
	}

	/** Returns the gene identifier associated with the gene */
	public GeneIdentifier GetGeneIdentifier()
	{
		return m_GID;
	}
	
	/** Returns the list of conditions with insufficient coverage */
	public TreeSet<String> GetInvalidConditions()
	{
		return m_vcInvalidConditions;
	}

	/**
	 *    Calculates and returns junction coverage ratios for two conditions.
	 *    Ratios are calculated for each sample and stored in a two-dimensional
	 *    result array, where each dimension corresponds to a different condition.
	 *     
	 *    If the project includes paired-samples the function makes sure to
	 *    add the ratios for both conditions in the same order. This is required
	 *    for the paired t-test that follows.
	 */
	public double[][] RetrieveJunctionRatios(SplicingWebApp app, String strConditionA, String strConditionB, CountElement jun1, CountElement jun2)
	{
		String strConditionType = app.GetSelectedConditionType();
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = app.GetProjectModel().GetSamplesPerCondition(strConditionType);
		TreeMap<String, TreeSet<String>> mapSamplesToIndividuals = app.GetProjectModel().GetSamplesPerIndividual();
		boolean bDataIsPaired = app.GetProjectModel().ConditionsHavePairedData(strConditionType, strConditionA, strConditionB);

		// test the inclusion junction reads against the exclusion junction reads
		Vector<Double> vcRatiosA = new Vector<Double>();
		Vector<Double> vcRatiosB = new Vector<Double>();
		
		TreeMap<String, Integer> mapReadsToSamplesForJun1 = GetCountsForJunction(jun1);
		TreeMap<String, Integer> mapReadsToSamplesForJun2 = GetCountsForJunction(jun2);
		
		// if paired, make sure to add the samples in the right order
		if(bDataIsPaired)
		{
			for(String strIndividual : mapSamplesToIndividuals.keySet())
			{
				double fValueSampleA = 0.0;
				double fValueSampleB = 0.0;
				
				for(String strSample : mapSamplesToIndividuals.get(strIndividual))
				{					
					int nReadsJun1 = 0;
					int nReadsJun2 = 0;
					
					// get reads for junction A
					if(mapReadsToSamplesForJun1 != null && mapReadsToSamplesForJun1.containsKey(strSample))
					{
						nReadsJun1 = mapReadsToSamplesForJun1.get(strSample);
					}
					else
					{
						nReadsJun1 = 0;
					}
					
					// get reads for junction B
					if(mapReadsToSamplesForJun2 != null && mapReadsToSamplesForJun2.containsKey(strSample))
					{
						nReadsJun2 = mapReadsToSamplesForJun2.get(strSample);
					}
					else
					{
						nReadsJun2 = 0;
					}
					
					double fSum = nReadsJun1 + nReadsJun2;
					
					// skip low coverage samples
					if(fSum < app.GetMinimumJunctionReads())
						continue;
					
					if(mapSamplesToConditions.get(strConditionA).contains(strSample))
					{
						if(fSum > 0)
							fValueSampleA = (double) nReadsJun1 / fSum;
						else
							fValueSampleA = Double.NaN;
					}
					else if(mapSamplesToConditions.get(strConditionB).contains(strSample))
					{
						if(fSum > 0)
							fValueSampleB = (double) nReadsJun1 / fSum;
						else
							fValueSampleB = Double.NaN;
					}
				}
				
				if(fValueSampleA != Double.NaN && fValueSampleB != Double.NaN)
				{
					vcRatiosA.add(fValueSampleA);
					vcRatiosB.add(fValueSampleB);
				}
			}
		}
		else
		{
			for(String strSample : mapSamplesToConditions.get(strConditionA))
			{
				int nInclusionReads = 0;
				int nExclusionReads = 0;

				if(mapReadsToSamplesForJun1.containsKey(strSample))
					nInclusionReads = mapReadsToSamplesForJun1.get(strSample);
				
				if(mapReadsToSamplesForJun2.containsKey(strSample))
					nExclusionReads = mapReadsToSamplesForJun2.get(strSample);
											
				double fSum = nInclusionReads + nExclusionReads;
				
				if(fSum == 0.0)
				{
					vcRatiosA.add(Double.NaN);
				}
				else
					vcRatiosA.add((double)nInclusionReads / fSum);
			}
			
			for(String strSample : mapSamplesToConditions.get(strConditionB))
			{
				int nInclusionReads = 0;
				int nExclusionReads = 0;
				
				if(mapReadsToSamplesForJun1.containsKey(strSample))
					nInclusionReads = mapReadsToSamplesForJun1.get(strSample);
				
				if(mapReadsToSamplesForJun2.containsKey(strSample))
					nExclusionReads = mapReadsToSamplesForJun2.get(strSample);
				
				double fSum = nInclusionReads + nExclusionReads;

				if(fSum == 0.0)
				{
					vcRatiosB.add(Double.NaN);
				}
				else
					vcRatiosB.add(nInclusionReads / fSum);
			}
		}
		
		double pRes[][] = new double[2][];
		pRes[0] = new double[vcRatiosA.size()];
		pRes[1] = new double[vcRatiosB.size()];
		
		for(int i=0; i<vcRatiosA.size(); i++)
			pRes[0][i] = vcRatiosA.get(i);
		
		for(int i=0; i<vcRatiosB.size(); i++)
			pRes[1][i] = vcRatiosB.get(i);
		
		return pRes;
	}

	/** Returns the start position of the coding region for a specified isoform */
	public int GetCodingStartForIsoform(String strIsoform)
	{
		return m_Gene.GetCodingStartForIsoform(strIsoform);
	}
	
	/** Returns the end position of the coding region for a specified isoform */
	public int GetCodingEndForIsoform(String strIsoform)
	{
		return m_Gene.GetCodingEndForIsoform(strIsoform); 
	}
}
