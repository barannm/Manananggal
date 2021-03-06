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

import java.awt.Color;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.LinearGradientPaint;
import java.awt.MultipleGradientPaint;
import java.awt.image.BufferedImage;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.nio.ByteBuffer;
import java.nio.charset.Charset;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.Locale;
import java.util.Map;
import java.util.Properties;
import java.util.Scanner;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.imageio.ImageIO;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.log4j.PropertyConfigurator;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;
import org.zkoss.image.AImage;
import org.zkoss.io.Files;
import org.zkoss.zhtml.Filedownload;
import org.zkoss.zk.ui.Component;
import org.zkoss.zk.ui.Executions;
import org.zkoss.zk.ui.event.ClientInfoEvent;
import org.zkoss.zk.ui.event.Event;
import org.zkoss.zk.ui.event.EventListener;
import org.zkoss.zk.ui.event.Events;
import org.zkoss.zk.ui.event.MouseEvent;
import org.zkoss.zk.ui.ext.AfterCompose;
import org.zkoss.zk.ui.util.Clients;
import org.zkoss.zul.Area;
import org.zkoss.zul.Bandbox;
import org.zkoss.zul.Bandpopup;
import org.zkoss.zul.Button;
import org.zkoss.zul.Checkbox;
import org.zkoss.zul.Combobox;
import org.zkoss.zul.Comboitem;
import org.zkoss.zul.Groupbox;
import org.zkoss.zul.Hbox;
import org.zkoss.zul.Hlayout;
import org.zkoss.zul.Html;
import org.zkoss.zul.Image;
import org.zkoss.zul.Imagemap;
import org.zkoss.zul.Label;
import org.zkoss.zul.Listbox;
import org.zkoss.zul.Listcell;
import org.zkoss.zul.Listhead;
import org.zkoss.zul.Listheader;
import org.zkoss.zul.Listitem;
import org.zkoss.zul.Messagebox;
import org.zkoss.zul.Paging;
import org.zkoss.zul.Radio;
import org.zkoss.zul.Radiogroup;
import org.zkoss.zul.Separator;
import org.zkoss.zul.Tabbox;
import org.zkoss.zul.Tabpanels;
import org.zkoss.zul.Tabs;
import org.zkoss.zul.Textbox;
import org.zkoss.zul.Tree;
import org.zkoss.zul.Treecell;
import org.zkoss.zul.Treechildren;
import org.zkoss.zul.Treecol;
import org.zkoss.zul.Treecols;
import org.zkoss.zul.Treeitem;
import org.zkoss.zul.Treerow;
import org.zkoss.zul.Vlayout;
import org.zkoss.zul.Window;
import org.zkoss.zul.event.PagingEvent;

import BioKit.Exon;
import BioKit.ExonGroup;
import BioKit.GTFGene;
import BioKit.GTFParser;
import BioKit.Gene;
import BioKit.RandomAccessGFFReader;
import BioKit.Utils;

/**
 *    The SplicingWebApp class is the core of the GUI application. It creates all GUI elements and processes GUI interactions.
 *    Currently, it also provides some functions for data processing (e.g. the selection of valid isoforms) that might be moved to
 *    the DataSupplier in later versions of the program.
 */
public class SplicingWebApp extends Window implements AfterCompose
{
	private static final long serialVersionUID = 1L;
	
	static final int		PROJECT_PAGING_SIZE		= 10;
	
	static final int		GTF_REFERENCE_FILE 		= 1;
	static final int		REFFLAT_REFERENCE_FILE 	= 2;
	
	static final int		AS_TYPE_EXON_SKIPPING				= 1;
	static final int		AS_TYPE_ALT_START_UNIQUE_JUN		= 2;	// alternative start exon with unique junction
	static final int		AS_TYPE_ALT_END_UNIQUE_JUN			= 3;	// alternative end exon with unique junction
	static final int		AS_TYPE_ALT_START_SHARED_JUN		= 4;	// alternative start exon that shares the junction with a 'middle' exon
	static final int		AS_TYPE_ALT_END_SHARED_JUN			= 5;	// alternative start exon that shares the junction with a 'middle' exon
	static final int		AS_TYPE_ALT_START_UNIQUE_JUN_DOUBLE	= 6;	// this type is only set if there are at least two start exons that have coverage ratio changes
	static final int		AS_TYPE_ALT_END_UNIQUE_JUN_DOUBLE	= 7;	// this type is only set if there are at least two end exons that have coverage ratio changes
	static final int		AS_TYPE_ALT_START_SHARED_JUN_DOUBLE	= 8;	// this type is only set if there are at least two start exons that have coverage ratio changes
	static final int		AS_TYPE_ALT_END_SHARED_JUN_DOUBLE	= 9;	// this type is only set if there are at least two end exons that have coverage ratio changes
	static final int		AS_TYPE_RETAINED_INTRON				= 10;
	static final int		AS_TYPE_ALT_5_PRIME_EXON_END		= 11;
	static final int		AS_TYPE_ALT_3_PRIME_EXON_END		= 12;
	
	static final int		IMPORTED_DATA_TYPE_MANA   			= 1;
	static final int 		IMPORTED_DATA_TYPE_RMATS  			= 2;
	static final int 		IMPORTED_DATA_TYPE_DEXSEQ 			= 3;
	static final int 		IMPORTED_DATA_TYPE_JSPLICE 			= 4;
	static final int 		IMPORTED_DATA_TYPE_CUFFDIFF			= 5;
	
	private int				m_nThreads;
	private String			m_strPathReferences;
	private String			m_strPathInput;
	private String			m_strPathHitLists;
	private String			m_strPathMMSeqResults;
	private String			m_strScreenshotPath;
	private String			m_strMMSeqExecutable;
	private String			m_strFileGTF;
	private String			m_strTmpFolder;
	private int				m_nReferenceType;
	RandomAccessGFFReader 	m_gffReader;
	
	private GeneIdentifierHandler 	m_geneIdentifierHandler;

	private int						m_nMinCovPerBase;			// Used as threshold to exclude exons with low coverage
	private double					m_fMinCoveredBases;			// Used as threshold to exclude exons that are only partially covered
	private int						m_nMinJunctionReads; 		// minimum coverage for junctions
	private double					m_fVariableExonThreshold;	// threshold at which exons are designated as differentially expressed (e.g. 0.3 if there must be a 30% difference in relative coverage)
	
	private Textbox					m_textboxMinCovPerBase;
	private Textbox					m_textboxMinCoveredBases;
	private Textbox					m_textboxMinJunctionReads;
	private Textbox					m_textboxVariableExonThreshold;
	
	private Textbox					m_textboxMinExonThreshold;	// threshold for the removal of small isoforms (viewer only)
	
	// handle to application and popup window
	SplicingWebApp 					m_hWindow;
	Window							m_windowPopup;
	
	private Vlayout					m_layoutCoveragePlot;
	private Vlayout					m_layoutIsoformPlot;
	
	// plot factories
	PlotFactory						m_plotFactory;
	PlotFactoryGTEX					m_plotFactoryGTEX;
	
	ResultListHandler				m_resultListHandler;
	
	// stores information when clicking an element in the isoform plot
	Parameters					m_parameters;
	
	// grants access to the SQL junction and exon count data
	private ProjectModel 		m_projectModel;
	
	// available options
	private boolean				m_bSkipFirstAndLastExon;
	private boolean				m_bDetectIntronRetentionEvents;
	
	private Textbox				m_txtboxURL;
	private Bandbox				m_bandboxProjects;
	private Bandpopup 			m_bandpopupProjects;
	private Listbox				m_listboxProjects;
	private Paging				m_pagingProjects;
	private String				m_strSelectedCondition;
	private Combobox			m_comboboxSelectedCondition;
	private Combobox			m_comboboxSelectedConditionType;
	private Combobox			m_comboboxSelectedGeneAnnotation;
	private String				m_strSelectedConditionType;
	private Bandbox				m_bandboxSelectedGene;
	private Listbox				m_listboxSelectedGene;
	private Tree				m_treeSelectedIsoforms;	
	private Textbox				m_textboxWindowWidth;
	
	private Checkbox 			m_checkboxUseMedian;
	private Checkbox			m_checkboxUseGeometricMean;
	private Checkbox 			m_checkboxUseMean;
	private Checkbox 			m_checkboxShowRelativeCoverage;
	private Checkbox 			m_checkboxUseLog2;
	private Checkbox			m_checkboxShowQuartiles;
	private Checkbox			m_checkboxColorExonsAndJunctionsByCoverage;
	private	Checkbox			m_checkboxHighlightUniqueExons;
	private Checkbox			m_checkboxCoverageGrid;
	private Checkbox			m_checkboxSwitchStrand;

	private TreeSet<String>		m_vcSelectedSamples;	
	private boolean				m_bUseReducedDataSet; 
	private Checkbox			m_checkboxUseReducedDataSet;
	private Checkbox			m_checkboxSkipFirstAndLastExon;
	private Checkbox			m_checkboxSkipIntronRetentionDetection;
	
	private boolean				m_bMMSeqAvailable;
	private boolean				m_bGTEXAvailableForGenes;
	private boolean				m_bGTEXAvailableForExons;
	
	private boolean				m_bIsoformSelectionChanged;
	private TreeSet<String>		m_vcValidIsoforms;
	private TreeSet<ProjectModel> m_vcProjectInfos;
	
	DataSupplier 				m_dataSupplier;
	
	private Combobox			m_comboboxSelectedEntropyIndex;
	private boolean				m_bShowEntropyData;	// Gini index, Gini-Simpson index, Theil index, Atkinson, Generalized entropy index,
	private boolean				m_bWindowSizeChanged;
	
	private Imagemap 			m_imgMapCoverage;
	private Imagemap 			m_imgMapIsoforms;
	private Imagemap 			m_imgMapColors;
	
	private boolean				m_bSampleSelectionChanged;
	private Tree				m_treeSelectedSamples;
	private Textbox				m_textboxSampleSelectionList;
	private Groupbox 			m_ColorSelectionGroupBox;
	
	// values important for the popup windows	
	private TreeSet<String> 				m_vcPreviouslySelectedSamples;	
	private TreeSet<String>					m_vcPreviouslySelectedIsoforms;
	private AnalysisResultHandler			m_vcASResults;
	private AnalysisResult					m_selectedResult;	
	
	private String 												m_strPathToExonData;
	private String												m_strFileDEXSeqIDs;
	private TreeMap<String, TreeMap<String, Vector<Double>>> 	m_mapCountsPerExonAndTissue;
	private TreeMap<String, Double> 							m_mapEntropyToExonicPart;
	
	/** Helper class used to store URL parameters */
	private class Parameters
	{
		ProjectModel 	m_project;
		long			m_nSelectedIsoforms;
		int				m_nProjectID;
		GeneIdentifier 	m_gid;
		boolean			m_bShowScreenshot;	// show just a screenshot instead of the web interface (or create it if necessary)
		TreeSet<Integer> m_vcSelectedSamples;
		String			m_strSelectedSamples;
		int				m_nReferenceID;
		
		public Parameters()
		{
			m_project 			= null;
			m_nSelectedIsoforms = 0;
			m_nProjectID		= 0;
			m_gid 				= null;
			m_bShowScreenshot	= false;
			m_vcSelectedSamples	= new TreeSet<Integer>();
			m_strSelectedSamples = null;
			m_nReferenceID		= -1;
		}
	};
	
	/**
	 *     Helper class used as thread to identify invalid isoforms based on the coverage of split reads
	 *     - run once for each isoform
	 */
	public class ThreadRemoveIrrelevantIsoformsBasedOnSplitReads implements Runnable
	{
		String 										m_strIsoform;
		String										m_strCondition;
		TreeSet<String> 							m_vcIrrelevantIsoforms;
		TreeMap<String, TreeSet<String>> 			m_vcSamplesPerCondition;
//		TreeMap<String, TreeMap<String, Integer>> 	m_mapJunctionReadCounts;
		boolean										m_bDebug;
		boolean										m_bSkipFirstAndLastExon;
		
		ThreadRemoveIrrelevantIsoformsBasedOnSplitReads(String strIsoform, String strCondition, TreeSet<String> vcIrrelevantIsoforms, TreeMap<String, TreeSet<String>> vcSamplesPerCondition, boolean bSkipFirstAndLastExon, boolean bDebug)
		{
			m_strIsoform 			= strIsoform;
			m_strCondition			= strCondition;
			m_vcIrrelevantIsoforms	= vcIrrelevantIsoforms;
			m_vcSamplesPerCondition = vcSamplesPerCondition;
//			m_mapJunctionReadCounts = mapJunctionReadCounts;
			m_bDebug				= bDebug;
			m_bSkipFirstAndLastExon = bSkipFirstAndLastExon;
		}
		
		@Override
		public void run()
		{
			//###################################
			//    get number of valid samples
			//###################################
			int nSamples = 0;
			for(String strSample : m_vcSamplesPerCondition.get(m_strCondition))
			{
				if(m_vcSelectedSamples.contains(strSample))
					nSamples++;
			}
			
			//##########################################################################################
			//          check whether all junctions of the isoform have sufficient coverage
			//##########################################################################################
			TreeSet<CountElement> vcJunctions = m_dataSupplier.GetJunctionsForIsoform(m_strIsoform);
			int nJunIdx = 0;
			
			for(CountElement jun : vcJunctions)
			{
				if(m_bSkipFirstAndLastExon)
				{
					// skip first and last junction
					if(nJunIdx == 0 || nJunIdx == vcJunctions.size()-1)
					{
						if(nJunIdx == 0)
							nJunIdx = 1;
						continue;
					}
				}				

				// get counts for the junction
				TreeMap<String, Integer> mapCountsToJunction = m_dataSupplier.GetCountsForJunction(jun);

				double pCounts[] = new double[nSamples];
				int nIdx = 0;
				for(String strSample : m_vcSamplesPerCondition.get(m_strCondition))
				{	
					if(m_vcSelectedSamples.contains(strSample))
					{
						if(mapCountsToJunction == null)
							pCounts[nIdx] = 0;
						else
						{
								
							if(mapCountsToJunction.containsKey(strSample))
							{
								pCounts[nIdx] = mapCountsToJunction.get(strSample);
							}
							else
							{
								pCounts[nIdx] = 0;
							}
						}
						nIdx++;
					}
				}
				
				//################################################
				//    check whether the coverage is sufficient
				//################################################
				double fMedian = StatUtils.percentile(pCounts, 50.0);
		
				if(fMedian < m_nMinJunctionReads)
				{
					synchronized(m_vcIrrelevantIsoforms)
					{
						if(m_bDebug)
						{
							System.out.println("[min junctions reads] -> removing " + m_strIsoform);
						}
						m_vcIrrelevantIsoforms.add(m_strIsoform);
						return;
					}
				}
				
				nJunIdx++;
			}
		}
	}

	/** 
	 *    Helper class used as thread to identify invalid isoforms based on insufficient [exon] coverage
     *    - run once for each isoform
	 */
	public class ThreadRemoveIrrelevantIsoformsBasedOnCoverage implements Runnable
	{
		String 										m_strIsoform;
		String										m_strCondition;
		TreeSet<String> 							m_vcIrrelevantIsoforms;
		TreeMap<String, TreeSet<String>> 			m_vcSamplesPerCondition;
		TreeMap<String, TreeMap<String, Integer>> 	m_mapJunctionReadCounts;
		boolean										m_bDebug;
		boolean										m_bSkipFirstAndLastExon;
		
		ThreadRemoveIrrelevantIsoformsBasedOnCoverage(String strIsoform, String strCondition, TreeSet<String> vcIrrelevantIsoforms, TreeMap<String, TreeSet<String>> vcSamplesPerCondition, boolean bSkipFirstAndLastExon, boolean bDebug)
		{
			m_strIsoform 			= strIsoform;
			m_strCondition			= strCondition;
			m_vcIrrelevantIsoforms	= vcIrrelevantIsoforms;
			m_vcSamplesPerCondition = vcSamplesPerCondition;
			m_bDebug 				= bDebug;
			m_bSkipFirstAndLastExon	= bSkipFirstAndLastExon;
		}
		
		@Override
		public void run()
		{
			// skip already irrelevant isoforms
			if(m_vcIrrelevantIsoforms.contains(m_strIsoform))
				return;
			
			int nSamples = 0;
			for(String strSample : m_vcSamplesPerCondition.get(m_strCondition))
			{
				if(m_vcSelectedSamples.contains(strSample))
					nSamples++;
			}
			
			//#######################################################################
			//    remove all isoforms where the coverage of some exon groups is much
			//    lower than for the other exon groups of the isoform
			//    IGNORE first and last exons
			//#######################################################################
			Exon pExons[] = m_dataSupplier.GetExonsForIsoform(m_strIsoform);
		
			for(Exon ex : pExons)
			{
				boolean bIsFirstExon = false;
				boolean bIsLastExon  = false;
				if(ex.equals(pExons[0]))
					bIsFirstExon = true;
				
				if(ex.equals(pExons[pExons.length-1]))
					bIsLastExon = true;
				
				// skip first and last exons?
				if(m_bSkipFirstAndLastExon && (bIsFirstExon || bIsLastExon))
					continue;
				
				// get average expression of the exon for the given condition
				double pCoverages[] 	= new double[nSamples];
				double pCoveredBases[] 	= new double[nSamples];
				
				int nSample = 0;
				for(String strSample : m_vcSamplesPerCondition.get(m_strCondition))
				{
					if(m_vcSelectedSamples.contains(strSample))
					{					
						// combine coverage of all exons within the exon group for first and last exons
						if(bIsFirstExon || bIsLastExon)
						{
							ExonGroup[] grps = m_dataSupplier.GetExonGroups();
							for(ExonGroup grp : grps)
							{
								if(grp.groupContainsExon(ex.getExonID()))
								{
									double fCoverage = 0.0f;
									double fCoverageLength = 0.0;
									
									for(Exon exon : grp.getExons())
									{
										// Only merge first with other first exons and last with other last exons
										String vcIsoforms[] = m_dataSupplier.GetIsoformNames();
										boolean bOkay = false;
										
										for(String strIsoform : vcIsoforms)
										{
											if(!m_vcIrrelevantIsoforms.contains(strIsoform))
											{
												Exon pIsoformExons[] = m_dataSupplier.GetExonsForIsoform(strIsoform);
												
												if((pIsoformExons[0] == exon && bIsFirstExon) || (pIsoformExons[pIsoformExons.length-1] == exon && bIsLastExon))
												{
													// exons must also overlap
													if(exon.getCodingStart() <= ex.getCodingStop() && exon.getCodingStop() >= ex.getCodingStart())
													{
														bOkay=true;
													}
												}
												// also keep it if it the last/first exon is completely contained in another exon
												else if(exon.getCodingStart() <= ex.getCodingStart() && exon.getCodingStop() >= ex.getCodingStop())
												{
													bOkay=true;
												}
											}
										}
										
										if(bOkay)
										{
											fCoverage		 += m_dataSupplier.GetRawCoverageForExonAndSample(exon, strSample);
											fCoverageLength  += m_dataSupplier.GetCoverageFractionForExonAndSample(exon, strSample);
										}
									}

									pCoverages[nSample] 	= fCoverage;
									pCoveredBases[nSample]	= fCoverageLength;
									
									break;
								}
							}
						}
						else
						{
							pCoverages[nSample] 	= m_dataSupplier.GetRawCoverageForExonAndSample(ex, strSample);
							pCoveredBases[nSample]	= m_dataSupplier.GetCoverageFractionForExonAndSample(ex, strSample);
						}
						nSample += 1;
					}
				}
				
				double fCoverageMedian		 = StatUtils.percentile(pCoverages, 50.0);
				double fCoveredBasesMedian 	 = StatUtils.percentile(pCoveredBases, 50.0);
				
				// do not apply number-of-covered-bases-test to first/last exons
				if((fCoveredBasesMedian < m_fMinCoveredBases && !bIsFirstExon && !bIsLastExon) || fCoverageMedian < m_nMinCovPerBase)
				{
					synchronized(m_vcIrrelevantIsoforms)
					{
						m_vcIrrelevantIsoforms.add(m_strIsoform);
					}						

					if(m_bDebug)
					{
						System.out.println("(" + m_strCondition + ") isoform irrelevant: " + m_strIsoform);
						if(fCoveredBasesMedian < m_fMinCoveredBases)
						{
							System.out.println("first: " + bIsFirstExon);
							System.out.println("last: " + bIsFirstExon);
							System.out.println("(" + m_strCondition + ") Reason: insufficiently covered exon: " + ex.getCodingStart() + "-" + ex.getCodingStop() + " " + fCoveredBasesMedian + " < " + m_fMinCoveredBases);
//							System.out.println(Arrays.toString(pCoveredBases));
						}
						else
						{
							System.out.println("first: " + bIsFirstExon);
							System.out.println("last: " + bIsFirstExon);
							System.out.println("(" + m_strCondition + ") Reason: coverage too low for exon: " + ex.getCodingStart() + "-" + ex.getCodingStop() + " " + fCoverageMedian + " < " + m_nMinCovPerBase);
						}
					}
					return;
				}
			}
		}
	}

	/** 
	 *    Helper class used as thread to read the coverage for a gene from bigwig files
	 *    - run once for each sample
	 */
	public class ThreadReadBigWigCoverage implements Runnable
	{
		Gene 					m_Gene;
		String 					m_strSample;
		TreeMap<String, int[]>	m_mapBigWigCoverageToSamples;
		String					m_strFile;
		
		ThreadReadBigWigCoverage(Gene gene, String strSample, String strFile, TreeMap<String, int[]> mapBigWigCoverageToSamples)
		{
			m_strSample 					= strSample;
			m_Gene							= gene;
			m_mapBigWigCoverageToSamples 	= mapBigWigCoverageToSamples;
			m_strFile						= strFile;
		}
		
		@Override
		public void run()
		{
			int nStart 	= m_Gene.getStart();
			int nEnd	= m_Gene.getStop();
			
			BigWigReader reader = null;
			
			try
			{
				reader = new BigWigReader(m_strFile);
			}
			catch(IOException ex)
			{
				System.out.println(ex.getMessage());
				return;
			}

			BigWigIterator it = reader.getBigWigIterator(m_Gene.getChromosome(), nStart-1, m_Gene.getChromosome(), nEnd, false);

			int pValues[] = new int[nEnd-nStart+1];
			
			while(it.hasNext())
			{
				WigItem item = it.next();
				int nIdx = item.getEndBase() - nStart;						
				int nValue = (int)item.getWigValue();
				
				pValues[nIdx] = nValue;
			}
			
			try
			{
				reader.close();
			}
			catch(IOException ex)
			{
				System.out.println(ex.getMessage());
			}
			
			synchronized(m_mapBigWigCoverageToSamples)
			{
				m_mapBigWigCoverageToSamples.put(m_strSample, pValues);
			}
		}
	}

	/** Inits the application without the GUI */
	public SplicingWebApp(boolean bEmpty) throws IOException
	{
		Configure4JLogger();
//		ReadAppPaths();
		
		m_hWindow 			= this;
		m_gffReader			= null;
		
		m_dataSupplier 		= new DataSupplier();

		m_vcValidIsoforms	= new TreeSet<String>();

		m_projectModel 		= new ProjectModel();
		m_vcProjectInfos 	= new TreeSet<ProjectModel>();
		
		m_vcSelectedSamples		= new TreeSet<String>();
		m_treeSelectedIsoforms	= null;
		m_bandboxSelectedGene	= null;
		
		m_imgMapIsoforms	= null;
		m_imgMapCoverage	= null;
		
		m_bMMSeqAvailable		 = false;
		m_bGTEXAvailableForGenes = false;
		m_bGTEXAvailableForExons = false;

		m_vcASResults 		= new AnalysisResultHandler();
		
		m_bWindowSizeChanged = false;
	}

	/** inits the application with GUI */
	public SplicingWebApp() throws Exception
	{
		Configure4JLogger();
		
		m_hWindow			= this;
		m_plotFactory	 	= new PlotFactory(m_hWindow);
		m_plotFactoryGTEX 	= new PlotFactoryGTEX(m_hWindow);
		m_gffReader			= null;

		ReadAppPaths();	
		
		m_dataSupplier 		= new DataSupplier();

		m_imgMapIsoforms	= new Imagemap();
		m_imgMapCoverage	= new Imagemap();
		m_imgMapColors		= null;
			
		m_bShowEntropyData		 = false;
		
		m_vcValidIsoforms					= new TreeSet<String>();
		m_bIsoformSelectionChanged			= false;
		m_bSkipFirstAndLastExon				= true;
		m_bDetectIntronRetentionEvents		= true;
		m_bWindowSizeChanged				= false;

		m_nMinJunctionReads					= 3;
		m_nMinCovPerBase					= 5;
		m_fMinCoveredBases					= 0.7;
		m_fVariableExonThreshold			= 0.05;		
		
		m_projectModel 						= new ProjectModel();
		
		m_vcPreviouslySelectedSamples		= new TreeSet<String>();
		m_vcPreviouslySelectedIsoforms		= new TreeSet<String>();
		
		m_bUseReducedDataSet		= true;
		m_bSampleSelectionChanged	= false;
		m_treeSelectedSamples		= null;
		m_vcSelectedSamples			= new TreeSet<String>();
		
		m_strSelectedCondition		= null;
		m_strSelectedConditionType 	= null;
		
		m_vcASResults 				= new AnalysisResultHandler();				
		m_geneIdentifierHandler 	= new GeneIdentifierHandler();
		m_vcProjectInfos 			= new TreeSet<ProjectModel>();
		
		// read human cross reference file by default
		String strFileEntrezGene = m_strPathReferences + "/biomart_cross_ref_human.txt";
		m_geneIdentifierHandler.Init(strFileEntrezGene);		
		
		addEventListener(Events.ON_CLIENT_INFO, new EventListener<ClientInfoEvent>()
		{
			public void onEvent(ClientInfoEvent event) throws Exception
			{
				Clients.showBusy("Please wait");
				Events.echoEvent("onLater", m_hWindow, null);
				
				int nClientWindowWidth  = event.getDesktopWidth();
				int nClientWindowHeight = event.getDesktopHeight();

				m_plotFactory.SetClientSize(nClientWindowWidth, nClientWindowHeight);
				
				m_textboxWindowWidth.setText("" + nClientWindowWidth);
				
				// some dirty code to enforce a maximum width on the project bandbox
				m_bandpopupProjects.setHflex("min");
				String strWidth = "width: " + (m_plotFactory.GetClientSize()[0] * 0.75) + "px;";
				String strMaxWidth = "max-" + strWidth;
				m_bandpopupProjects.setStyle(strMaxWidth);
				m_bandpopupProjects.setZclass(strMaxWidth);

				Clients.resize(m_bandpopupProjects);

				m_pagingProjects.setHflex("min");
				Clients.resize(m_pagingProjects);
			}
		});
		
		addEventListener("onLater", new EventListener<Event>()
		{
			public void onEvent(Event evt) throws Exception
			{
				// process parameters
				if(m_parameters.m_project == null)
				{
					Clients.clearBusy();
					return;
				}
				
				// if a screenshot is requested, check if it exists
				String strFileScreenshot = null;
				if(m_parameters.m_bShowScreenshot)
				{
					strFileScreenshot = m_strScreenshotPath + "/" + m_parameters.m_nProjectID + "_" + m_parameters.m_nReferenceID + "_" + m_parameters.m_gid.m_strEnsemblGeneID;
					if(m_parameters.m_nSelectedIsoforms != 0)
						strFileScreenshot += "_" + m_parameters.m_nSelectedIsoforms;
					
					if(m_parameters.m_strSelectedSamples != null)
						strFileScreenshot += "_samples_" + m_parameters.m_strSelectedSamples;
					
					strFileScreenshot += ".png";

					File pScreenshot = new File(strFileScreenshot);
					
					if(pScreenshot.exists())
					{
						m_windowPopup = new Window();
						m_windowPopup.setTitle("coverage");
						m_windowPopup.setSizable(true);
						m_windowPopup.setMinimizable(true);
						m_windowPopup.setMaximizable(true);
						m_windowPopup.setClosable(true);
						m_windowPopup.setBorder(true);
						m_windowPopup.doPopup();
						m_windowPopup.setTopmost();
						m_windowPopup.setParent(m_hWindow);
						
						int pClientSize[] = m_plotFactory.GetClientSize();
						
						Vlayout layout = new Vlayout();
						layout.setParent(m_windowPopup);
						layout.setStyle("overflow: auto; width: 100%; height: 100%; max-width: " + pClientSize[0] + "px; max-height: " + pClientSize[1] + "px;");

						AImage img = new AImage(strFileScreenshot);

						Image imgSS = new Image();
						imgSS.setContent(img);
						imgSS.setParent(layout);
						
						m_windowPopup.setVflex("0");
						m_windowPopup.setHflex("0");

						m_windowPopup.setWidth("100%");
						m_windowPopup.setHeight("100%");
						m_windowPopup.setPosition("left,top");
						m_windowPopup.setVisible(true);

						Clients.clearBusy();
						return;
					}
				}
				
				// open the project
				try
				{
					if(m_parameters.m_vcSelectedSamples.size() != 0)
						OnProjectChange(m_parameters.m_project.GetFullPathOfProjectFile(), true);
					else
						OnProjectChange(m_parameters.m_project.GetFullPathOfProjectFile(), false);
				}
				catch(Exception e)
				{
					ErrorMessage.ShowError("ERROR: could not open project file: " + m_parameters.m_project.GetFullPathOfProjectFile());
					e.printStackTrace();
					Clients.clearBusy();
					return;
				}
				
				m_bandboxProjects.setText(m_parameters.m_project.GetProjectName());
				
				if(m_parameters.m_gid == null)
				{
					Clients.clearBusy();
					return;
				}
				
				// select gene annotation file
				if(m_parameters.m_nReferenceID != -1)
				{
					System.out.println("selected reference ID: " + m_parameters.m_nReferenceID);
					for(Comboitem item : m_comboboxSelectedGeneAnnotation.getItems())
					{
						if(item.getValue().hashCode() == m_parameters.m_nReferenceID)
						{
							m_comboboxSelectedGeneAnnotation.setSelectedItem(item);
							m_strFileGTF = item.getValue();
							LoadGeneAnnotation(m_strFileGTF);
							break;
						}
					}
				}
				
				m_bandboxSelectedGene.setText(m_parameters.m_gid.m_strEnsemblGeneID);
				m_bandboxSelectedGene.setValue(m_parameters.m_gid.m_strEnsemblGeneID);
				try
				{
					OnGeneChange(false, false);
				}
				catch(Exception e)
				{
					ErrorMessage.ShowError("ERROR: could not open gene: " + m_parameters.m_gid.m_strApprovedGeneSymbol + " (" + m_parameters.m_gid.m_strEnsemblGeneID + ")");
					e.printStackTrace();
					Clients.clearBusy();
					return;
				}
				
				boolean bRequiresRedraw = false;
				// select samples
				if(m_parameters.m_vcSelectedSamples.size() != 0)
				{
					m_treeSelectedSamples.clearSelection();
					m_vcPreviouslySelectedSamples.clear();
					m_vcSelectedSamples.clear();
					
					int i=0;
					for(Treeitem item : m_treeSelectedSamples.getItems())
					{
						if(m_parameters.m_vcSelectedSamples.contains(i))
						{
							m_treeSelectedSamples.addItemToSelection(item);
							m_vcPreviouslySelectedSamples.add(item.getLabel());
							m_vcSelectedSamples.add(item.getLabel());
						}
						i++;
					}
					bRequiresRedraw = true;
					UpdateComboboxForSampleSelection();
				}
				
				// LONG supports up to 64 isoforms!
				if(m_parameters.m_nSelectedIsoforms != 0)
				{
					long nSelectedIsoforms = m_parameters.m_nSelectedIsoforms;
					m_vcValidIsoforms.clear();
					
					// check if specific isoforms were selected
					try
					{
						m_treeSelectedIsoforms.clearSelection();
						
						Treeitem pItems[] = new Treeitem[m_treeSelectedIsoforms.getItemCount()];
						m_treeSelectedIsoforms.getItems().toArray(pItems);
	
						for(int i=pItems.length; i>=0; i--)
						{
							long nFlag = (long) Math.pow(2, i);
							
							if( (nSelectedIsoforms & nFlag) == nFlag)
							{
								nSelectedIsoforms -= nFlag;
								m_treeSelectedIsoforms.addItemToSelection(pItems[i]);
								m_vcValidIsoforms.add((String)pItems[i].getValue());
							}
						}
						bRequiresRedraw = true;
					}
					catch(Exception e)
					{
						ErrorMessage.ShowError("ERROR: could not select isoforms");
						e.printStackTrace();
						Clients.clearBusy();
						return;
					}
				}
				
				if(bRequiresRedraw)
				{
					// now get new results for the selected isoforms and samples
					OnIsoformChange(false);
					m_plotFactory.RequestCoverageRedraw();
					m_plotFactory.DrawPlots();
				}
				
				// 
				if(m_parameters.m_bShowScreenshot)
				{
					// generate
					BufferedImage imgCoverage = ImageIO.read(m_imgMapCoverage.getContent().getStreamData());
					BufferedImage imgIsoforms = ImageIO.read(m_imgMapIsoforms.getContent().getStreamData());
					
					int nWidth = Math.max(imgCoverage.getWidth(), imgIsoforms.getWidth());
					int nHeight = imgCoverage.getHeight() + imgIsoforms.getHeight();
					BufferedImage combined = new BufferedImage(nWidth, nHeight, BufferedImage.TYPE_INT_ARGB);
					
					Graphics g = combined.getGraphics();
					g.drawImage(imgCoverage, 0, 0, null);
					g.drawImage(imgIsoforms, 0, imgCoverage.getHeight(), null);

					// Save as new image
					ByteArrayOutputStream os = new ByteArrayOutputStream();
					ImageIO.write(combined, "png", os);
					InputStream is = new ByteArrayInputStream(os.toByteArray());

					File pOut = new File(strFileScreenshot);
					Files.copy(pOut, is);
					
					// show screenshot
					m_windowPopup = new Window();
					m_windowPopup.setTitle("coverage");
					m_windowPopup.setSizable(true);
					m_windowPopup.setMinimizable(true);
					m_windowPopup.setMaximizable(true);
					m_windowPopup.setClosable(true);
					m_windowPopup.setBorder(true);
					m_windowPopup.doPopup();
					m_windowPopup.setTopmost();
					m_windowPopup.setParent(m_hWindow);
					
					int pClientSize[] = m_plotFactory.GetClientSize();
					
					Vlayout layout = new Vlayout();
					layout.setParent(m_windowPopup);
					layout.setStyle("overflow: auto; width: 100%; height: 100%; max-width: " + pClientSize[0] + "px; max-height: " + pClientSize[1] + "px;");

					AImage img = new AImage(strFileScreenshot);

					Image imgSS = new Image();
					imgSS.setContent(img);
					imgSS.setParent(layout);
					
					m_windowPopup.setVflex("0");
					m_windowPopup.setHflex("0");

					m_windowPopup.setWidth("100%");
					m_windowPopup.setHeight("100%");
					m_windowPopup.setPosition("left,top");
					m_windowPopup.setVisible(true);
				}
				
				Clients.clearBusy();
			}
		});
		
		//######################################
		//           add gui elements
		//######################################
		Vlayout mainLayout = new Vlayout();
		mainLayout.setParent(this);

		//#################################
		//      add buttons/options
		//#################################
		Hlayout layoutH = new Hlayout();
		layoutH.setParent(mainLayout);
		layoutH.setHeight("570px");
		layoutH.setStyle("overflow:auto; margin-top: 5px; margin-right: 10px; margin-bottom: 0px;");
		
		Vlayout layoutV = new Vlayout();		
		layoutV.setVflex("min");
		layoutV.setParent(layoutH);
		
		//####################################
		//    Add div and layout for coverage plot
		//####################################
		Groupbox grpBox = new Groupbox();
		grpBox.setTitle("Coverage tracks");
		grpBox.setParent(mainLayout);
		grpBox.setMold("3d");
		grpBox.setClosable(false);
		m_layoutCoveragePlot = new Vlayout();
		m_layoutCoveragePlot.setId("CoveragePlot");
		m_layoutCoveragePlot.setParent(grpBox);
		m_layoutCoveragePlot.setAttribute("org.zkoss.zul.nativebar", "true");
		m_layoutCoveragePlot.setStyle("overflow: auto;");
		m_layoutCoveragePlot.setHeight("240px");

		//####################################
		//    Add div and layout for isoform plot
		//####################################
		grpBox = new Groupbox();
		grpBox.setTitle("Isoform View");
		grpBox.setParent(mainLayout);
		grpBox.setMold("3d");
		grpBox.setClosable(false);
		m_layoutIsoformPlot = new Vlayout();
		m_layoutIsoformPlot.setId("IsoformPlot");
		m_layoutIsoformPlot.setParent(grpBox);		
//		m_layoutIsoformPlot.setAttribute("org.zkoss.zul.nativebar", "true");
		m_layoutIsoformPlot.setStyle("overflow: auto;");
		m_layoutIsoformPlot.setHeight("600px");

		//########################################
		//    add buttons for setting selection 
		//########################################
		AddSettingsMenu(layoutV);

		//########################################
		//    add buttons for option selection 
		//########################################
		AddOptionsMenu(layoutV);
		
		//########################################
		//   add combobox for sample selection 
		//########################################
		layoutV = new Vlayout();		
		layoutV.setVflex("min");
		layoutV.setParent(layoutH);
		AddComboboxForSamples(layoutV);
		
		//########################################
		//   add combobox for color selection 
		//########################################
		layoutV = new Vlayout();		
		layoutV.setVflex("min");
		layoutV.setParent(layoutH);
		AddColorSelectionMenu(layoutV);

		//########################################
		//   add combobox for isoform selection 
		//########################################
		layoutV = new Vlayout();		
		layoutV.setVflex("min");
		layoutV.setParent(layoutH);
		AddComboboxForIsoformSelection(layoutV);

		//########################################
		//        add advanced options 
		//########################################
		layoutV = new Vlayout();		
		layoutV.setVflex("min");
		layoutV.setParent(layoutH);
		AddAdvancedOptionsMenu(layoutV);
		
		//########################################
		//        add list of AS hits
		//########################################
		Tabbox tabbox = new Tabbox();
		tabbox.setVflex("min");
		tabbox.setParent(layoutH);
		
		Tabs tabs = new Tabs();
		tabs.setParent(tabbox);
		
		Tabpanels panels = new Tabpanels();
		panels.setParent(tabbox);
		
		tabs.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Tabs tabs = (Tabs)event.getTarget();
				tabs.invalidate();
			}
		});
		
		m_resultListHandler = new ResultListHandler(tabs, panels, m_hWindow);
		
		//############################
		//    prepare popup window
		//############################
		m_windowPopup = new Window();
		m_windowPopup.setTitle("coverage");
		m_windowPopup.setSizable(true);
		m_windowPopup.setClosable(true);
		m_windowPopup.setMinimizable(true);
		m_windowPopup.setMaximizable(true);
		m_windowPopup.setBorder(true);
		m_windowPopup.doPopup();
		m_windowPopup.setTopmost();
		
		// everything loaded, thus check if any parameters were set
		ProcessParameters();
	}
	
	/** Retrieves and stores the URL parameters */
	public void ProcessParameters()
	{
		m_parameters = new Parameters();
				
		// check if the project exists
		String strParameter = Executions.getCurrent().getParameter("ID");
		if(strParameter == null)
			return;
		
		int nID = Integer.parseInt(strParameter);
		for(ProjectModel project : m_vcProjectInfos)
		{			
			if(project.GetProjectName().hashCode() == nID)
			{
				m_parameters.m_nProjectID = nID;
				m_parameters.m_project = project;
				break;
			}
		}
				
		if(m_parameters.m_project == null)
		{
			ErrorMessage.ShowError("Invalid project identifier: " + nID);
			return;
		}
		
		// check if gene annotation file exists
		strParameter = Executions.getCurrent().getParameter("reference");
		if(strParameter != null)
		{
			m_parameters.m_nReferenceID = Integer.parseInt(strParameter);
		}
				
		// check if the gene is valid
		strParameter = Executions.getCurrent().getParameter("gene");
		if(strParameter == null)
			return;
		
		GeneIdentifier gid = m_geneIdentifierHandler.GetGeneIdentifierForGene(strParameter, null);
		if(gid == null)
		{
			ErrorMessage.ShowError("unknown gene identifier: " + strParameter);
			return;
		}
		m_parameters.m_gid = gid;

		strParameter = Executions.getCurrent().getParameter("isoforms");
		if(strParameter != null)
			m_parameters.m_nSelectedIsoforms = Long.parseLong(strParameter);

		// 64 samples can  be stored per long, thus we will need multiple
		strParameter = Executions.getCurrent().getParameter("samples");
		
		if(strParameter != null)
		{
			m_parameters.m_strSelectedSamples = strParameter.replace(",", "_");
			
			String pSamples[] = strParameter.split(",");
			for(int i=0; i<pSamples.length; i++)
			{
				long nSelectedSamples = Long.parseLong(pSamples[i]);
				int nOffset = i * 63;
				
				for(int j=62; j>=0; j--)
				{
					long nFlag = (long) Math.pow(2, j);
					
					if( (nSelectedSamples & nFlag) == nFlag)
					{
						nSelectedSamples -= nFlag;
						m_parameters.m_vcSelectedSamples.add(nOffset + j);
					}
				}
			}
		}
		
		strParameter = Executions.getCurrent().getParameter("screenshot");
		if(strParameter != null)
		{
			m_parameters.m_bShowScreenshot = true;
		}
	}
	
	/** Returns the GTF/GFF reader */
	public RandomAccessGFFReader GetGFFReader()
	{
		return m_gffReader;
	}
	
	/** Returns the region (part of the GUI) that shows the coverage plots */
	public Vlayout GetCoveragePlotRegion()
	{
		return m_layoutCoveragePlot;
	}

	/** Returns the region (part of the GUI) that shows the isoform plots */
	public Vlayout GetIsoformPlotRegion()
	{
		return m_layoutIsoformPlot;
	}
	
	/** Returns the image map associated with the coverage plot */
	public Imagemap GetImageMapForCoveragePlot()
	{
		return m_imgMapCoverage;
	}
	
	/** Returns the image map associated with the isoform plot */
	public Imagemap GetImageMapForIsoforms()
	{
		return m_imgMapIsoforms;
	}
	
	/** Returns whether intron retention event detection is active */
	public boolean DetectIntronRetentionEvents()
	{
		return m_bDetectIntronRetentionEvents;
	}
	
	/** Returns the list of valid isoforms (=isoforms that pass the automatic selection) */
	public TreeSet<String> GetValidIsoforms()
	{
		return m_vcValidIsoforms;
	}
	
	/** Returns the maximum allowed threads */
	public int GetMaximumThreads()
	{
		return m_nThreads;
	}
	
	/** Returns the selected condition type */
	public String GetSelectedConditionType()
	{
		return m_strSelectedConditionType;
	}
	
	/** Returns the selected condition */
	public String GetSelectedCondition()
	{
		return m_strSelectedCondition;
	}
	
	/** Returns the path to the MMSeq results */
	public String GetPathToMMSeqResults()
	{
		return m_strPathMMSeqResults;
	}
		
	/** Returns a list of all selected samples */
	public TreeSet<String> GetSelectedSamples()
	{
		return m_vcSelectedSamples;
	}
	
	/** Returns the minimum coverage per base used to identify valid isoforms */
	public int GetMinimumCoveragePerBase()
	{
		return m_nMinCovPerBase;
	}
	
	/** Returns the minimum number of junction spanning reads used to identify valid isoforms */
	public int GetMinimumJunctionReads()
	{
		return m_nMinJunctionReads;
	}
	
	/** Returns the minimum fraction (e.g. 0.7) of bases of a gene that must be covered to identify valid isoforms */
	public double GetMinimumCoveredBases()
	{
		return m_fMinCoveredBases;
	}
		
	/** Returns the coverage ratio threshold used to identify "variable" (=potentially alternatively spliced) exon */
	public double GetVariableExonThreshold()
	{
		return m_fVariableExonThreshold;
	}
	
	/** Returns the currently selected result from the result list */
	public AnalysisResult GetSelectedResult()
	{
		return m_selectedResult;
	}
	
	/** Returns the GTF file */
	public String GetGTFFile()
	{
		return m_strFileGTF;
	}
	
	/** Returns the project model that stores all information associated with the project */
	public ProjectModel GetProjectModel()
	{
		return m_projectModel;
	}
	
	/** Returns the data supplier. The data supplier stores coverage and read data for the currently selected gene */
	public DataSupplier GetDataSupplier()
	{
		return m_dataSupplier;
	}
	
	/** Returns the result handler that stores the potential AS events for the current view */
	public AnalysisResultHandler GetResultHandler()
	{
		return m_vcASResults;
	}
	
	/** Returns the result list handler that manages the result grid */
	public ResultListHandler GetResultListHandler()
	{
		return m_resultListHandler;
	}
	
	/** Returns the gene identifier handler that stores gene symbols and various gene IDs per gene */ 
	public GeneIdentifierHandler GetGeneIdentifierHandler()
	{
		return m_geneIdentifierHandler;
	}
	
	/** 
	 *    Sets the number of threads available to the application
	 *    - only used by the console application (the GUI reads this value from the app_config file)
	 */
	public void SetNumberOfThreads(int nThreads)
	{
		m_nThreads = nThreads;
		System.out.println("set max. number of threads to " + nThreads);
	}
	
	/** 
	 *    Sets the currently selected condition type
	 *    - only used by the console application
	 */
	public void SetConditionType(String strConditionType)
	{
		m_strSelectedConditionType = strConditionType;
	}
	
	/**
	 *    Removes data associated with the current gene
	 *    - only used by the console application
	 */
	public void ClearCurrentGeneData()
	{
		m_dataSupplier.Clear();
		m_vcValidIsoforms.clear();
	}
	
	/** Returns whether entropy values should be plotted */
	public boolean IsEntropyEnabled()
	{
		return m_bShowEntropyData;
	}
	
	/**
	 *    Returns whether entropy values should be plotted per exonic part
	 *    - redundant to the above, should be removed in a future version
	 */
	public TreeMap<String, Double> GetEntropyPerExonicPart()
	{
		return m_mapEntropyToExonicPart;
	}
	
	/** Returns the read counts of a specified exonic part for all tissues */	
	public TreeMap<String, Vector<Double>> GetCountsForGTEXExonicPart(String strExonicPartID)
	{
		return m_mapCountsPerExonAndTissue.get(strExonicPartID);
	}

	/** Init the Log4j Logger to get reid of the warning message... */
	public void Configure4JLogger()
	{
		// just disable it, we are not using it anyway
		Properties properties = new Properties();
		
		properties.setProperty("log4j.debug","FALSE");
	    properties.setProperty("log4j.rootLogger","OFF");

	    /*
	    properties.setProperty("log4j.rootLogger","TRACE,stdout");
	    properties.setProperty("log4j.rootCategory","TRACE");
	      
	    properties.setProperty("log4j.appender.stdout",     "org.apache.log4j.ConsoleAppender");
	    properties.setProperty("log4j.appender.stdout.layout",  "org.apache.log4j.PatternLayout");
	    properties.setProperty("log4j.appender.stdout.layout.ConversionPattern","%d{yyyy/MM/dd HH:mm:ss.SSS} [%5p] %t (%F) - %m%n");
	    */

	    PropertyConfigurator.configure(properties);
	}

	/** 
	 *    Sets the application paths, e.g. tells the application where to find
	 *    third party tools, reference sequences and project data
	 */
	public void ReadAppPaths() throws IOException
	{
		String strFileName = null;
		File pIn = null;
				
		// search for an environmental variable first
		Map<String, String> env = System.getenv();
		for(String strVariable : env.keySet())
        {
        	if(strVariable.equals("MANANANGGAL_CONFIG"))
        	{
        		strFileName = env.get(strVariable);
        		pIn = new File(strFileName);        		
        	}
        }

		if(pIn == null || !pIn.exists())
		{
			strFileName = System.getProperty("user.home") + "/app_paths.config";
			pIn = new File(strFileName);
			System.out.println("... searching for app paths in home directory.");
		}
		else
		{
			System.out.println("using environmental variable for app paths.");
		}
		
		if(!pIn.exists())
		{
			System.out.println("FATAL ERROR: could not open app_paths.config.");
			return;
		}
		
		// use 1 thread by default
		m_nThreads = 1;
		String strPathGeneGTEX	= null;
		m_strPathMMSeqResults	= null;
		m_strMMSeqExecutable	= null;
		m_strScreenshotPath		= null;

		if(pIn.exists())
		{
			Scanner pScanner = new Scanner(pIn);

			while(pScanner.hasNextLine())
			{
				String strLine = pScanner.nextLine();
				
				String pSplit[] = strLine.split("=");
				
				switch(pSplit[0])
				{
					case "path_references": 		m_strPathReferences = pSplit[1].trim().replace("\"", ""); break;
					case "path_projects": 			m_strPathInput		= pSplit[1].trim().replace("\"", ""); break;
					case "path_GTEX_exon_IDs":		m_strFileDEXSeqIDs	= pSplit[1].trim().replace("\"", ""); break;
					case "path_GTEX_gene_counts": 	strPathGeneGTEX 	= pSplit[1].trim().replace("\"", ""); break;
					case "path_GTEX_exon_counts": 	m_strPathToExonData = pSplit[1].trim().replace("\"", ""); break;
					case "path_tmp": 				m_strTmpFolder 		= pSplit[1].trim().replace("\"", ""); break;
					case "path_screenshots":		m_strScreenshotPath	= pSplit[1].trim().replace("\"", ""); break;
					case "path_MMSEQ_results": 		m_strPathMMSeqResults = pSplit[1].trim().replace("\"", ""); break;
					case "num_threads":				m_nThreads 			= Integer.parseInt(pSplit[1].trim());
					case "path_hits":				m_strPathHitLists	= pSplit[1].trim().replace("\"", ""); break;
					case "path_MMSEQ_executable":	m_strMMSeqExecutable = pSplit[1].trim().replace("\"", ""); break;
				}
			}
			
			m_bMMSeqAvailable		 = false;
			m_bGTEXAvailableForGenes = false;
			m_bGTEXAvailableForExons = false;
			
			// check if MMseq is available
			File pTestFile = new File(m_strMMSeqExecutable);
			if(pTestFile.exists())
				m_bMMSeqAvailable = true;
			else
				System.out.println("could not detect MMSeq executable at specified path: " + m_strMMSeqExecutable);
			
			// check if GTEX gene counts are available
			File pTestFolder = new File(strPathGeneGTEX);
			if(pTestFolder.exists())
			{
				m_bGTEXAvailableForGenes = true;
			}
			else
				System.out.println("could not detect GTEX gene data at specified path: " + strPathGeneGTEX);
			
			// check if GTEX exon counts are available
			pTestFolder = new File(m_strPathToExonData);
			if(pTestFolder.exists())
			{
				m_bGTEXAvailableForExons = true;
			}
			else
				System.out.println("could not detect GTEX exon data at specified path: " + strPathGeneGTEX);
			
			System.out.println("input directory: " + m_strPathInput);
			
			pScanner.close();
		}
		else
		{
			System.out.println("failed to read: " + strFileName);
		}
		
		m_plotFactoryGTEX.SetPaths(strPathGeneGTEX);
	}
	
	/** Reads a gene (specified by gene identifier) from a GTF file	 */
	public Gene ReadGeneFromFile(GeneIdentifier gid, String strFileGTF) throws Exception
	{
		if(strFileGTF == null)
		{
			Messagebox.show("No gene annotation database selected!");
			return null;
		}
		
		if(m_gffReader == null || !m_gffReader.GetFileName().equals(strFileGTF))
		{
			LoadGeneAnnotation(strFileGTF);
		}
		
		if(m_gffReader == null)
			return null;

		Gene gene = null;
		
		String strGene = gid.m_strEnsemblGeneID;
		gene = m_gffReader.ReadGene(strGene);
		
		// try again with m_strEntrezGeneID
		if(gene == null)
		{
			strGene = gid.m_strEntrezGeneID;
			gene = m_gffReader.ReadGene(strGene);
		}
		
		if(gene == null)
		{
			ErrorMessage.ShowError("Could not find " + gid.m_strApprovedGeneSymbol + " (" + gid.m_strEnsemblGeneID  + ") in GTF file.");
			return null;
		}
		
		return gene;
	}
	
	/** Resets the list of valid isoforms (selects all available isoforms) */
	public void InitValidIsoforms(Gene gene)
	{
		m_vcValidIsoforms.clear();
		for(String strIsoform : gene.getArrayOfGeneProductNames())
		{
			if(m_vcPreviouslySelectedIsoforms != null)
				m_vcPreviouslySelectedIsoforms.add(strIsoform);
			m_vcValidIsoforms.add(strIsoform);
		}
	}
	
	/** 
	 *    Calls ReadGeneFromFile (unless the 'gene' parameter is not null) and then uses the gene to fill the data supplier with data.
	 *    Also resets the valid isoforms.
	 *    the gene parameter is optional and may be null
	 */
	public void ProcessGeneInformation(GeneIdentifier gid, String strFileGTF, Gene gene) throws Exception
	{
		if(gene == null)
			gene = ReadGeneFromFile(gid, strFileGTF);
		
		if(gene == null)
			return;
		
		// reset valid isoforms
		InitValidIsoforms(gene);
		
		m_bIsoformSelectionChanged = true;
	
		// get coverage from bigwig files and junction counts, also detects invalid samples and conditions
		PrepareDataSupplier(gene);

		// update isoform selection window
		OnIsoformChange(true);
	}
	
	/**
	 *    Fills the data supplier with data such as:
	 *   - exon coverage from bigwig files,
	 *   - junction counts from count files
	 *   Also identifies junctions with insufficient coverage.
	 */
	public void PrepareDataSupplier(Gene gene)
	{
		// prepare exon and exon group arrays
		m_dataSupplier = new DataSupplier(m_hWindow, gene);
		
		m_dataSupplier.RetrieveCoverageData(m_hWindow);
		m_dataSupplier.RetrieveJunctionCounts(m_hWindow);
		m_dataSupplier.IdentifyInvalidJunctions(m_hWindow, false);
	}

	//TODO should be moved to the data supplier
	/** Retrieves the coverage of a genomic region specified by reference name, start and end position for a given sample */
	public double[] GetCoverageForRegion(String strRef, int nStart, int nEnd, String strSample, boolean bAdjustForSizeFactor) throws IOException
	{
		// get size factor
		TreeMap<String, Double> mapSizeFactors = m_projectModel.GetSizeFactors();
		double fSizeFactor = mapSizeFactors.get(strSample);
		
		// init empty array
		int nLength = nEnd - nStart + 1;
		double pCoverage[] = new double[nLength];
		Arrays.fill(pCoverage, 0.0);
	
		double pBigWigCoverage[] = null;
		
		pBigWigCoverage = m_dataSupplier.GetGeneCoverageArrayForSample(strSample);
		int nArrayStart = m_dataSupplier.GetGeneStart();
		
		if(pBigWigCoverage == null || nStart < m_dataSupplier.GetGeneStart() || nEnd > m_dataSupplier.GetGeneEnd())
		{
			nArrayStart = nStart;
			
			if(m_dataSupplier.GetGene() == null)
			{
				ErrorMessage.ShowError("GetBigWigCoverageForGene - No gene selected");
				return null;
			}
			
			TreeMap<String, String> mapBigWigFiles 	= m_projectModel.GetBigWigFilesForSamples();
			
			BigWigReader reader = null;
			try
			{
				reader = new BigWigReader(mapBigWigFiles.get(strSample));
			}
			catch(IOException ex)
			{
				System.out.println(ex.getMessage());
				return null;
			}
	
			BigWigIterator it = reader.getBigWigIterator(m_dataSupplier.GetReferenceName(), nStart-1, m_dataSupplier.GetReferenceName(), nEnd, false);

			double pValues[] = new double[nEnd-nStart+1];
			
			while(it.hasNext())
			{
				WigItem item = it.next();
				int nIdx = item.getEndBase() - nStart;
				int nValue = (int)item.getWigValue();				
				pValues[nIdx] = nValue;
			}
			
			try
			{
				reader.close();
			}
			catch(IOException ex)
			{
				System.out.println(ex.getMessage());
			}
			
			pBigWigCoverage = pValues;
		}
		
		for(int i=nStart; i<=nEnd; i++)
		{
			int nPos = i - nArrayStart;
			int nIdx = i - nStart;
			
			pCoverage[nIdx] = pBigWigCoverage[nPos];
			
			if(bAdjustForSizeFactor)
				pCoverage[nIdx] /= fSizeFactor;
		}

		return pCoverage;
	}

	/**
	 *    This function is triggered when a new project is selected. It populates the project model with new data,
	 *    clears all values associated with the previous project/gene and updates several GUI elements 
	 */
	public void OnProjectChange(String strProjectFile, boolean bHideMessage) throws Exception
	{
		m_strSelectedCondition 			= null;
		m_strSelectedConditionType		= null;
		
		m_projectModel.clear();
		if(!m_projectModel.Init(strProjectFile, -1, -1, -1.0, -1, true, false))
		{
			ErrorMessage.ShowError("Failed to open project: " + strProjectFile);
			return;
		}
		
		if(m_projectModel.GetSamples().size() > 50)
		{
			if(!bHideMessage)
				Messagebox.show("The project has more than 50 samples. To speed up things, detection of retained introns was disabled. You can enable retained intron detection in the advanced options tab again.");

			m_bDetectIntronRetentionEvents = false;
			m_checkboxSkipIntronRetentionDetection.setChecked(m_bDetectIntronRetentionEvents);
		}
		
		// read last used gene annotation file hashcode from file
		String strFileIn = m_projectModel.GetFullPathOfProjectFile() + ".last_reference";
		File pIn = new File(strFileIn);
		if(pIn.exists())
		{
			Scanner pRefIn = new Scanner(pIn);

			String strRef = pRefIn.nextLine();
			int nReferenceID = Integer.parseInt(strRef);
			pRefIn.close();
			
			// select gene annotation file
			for(Comboitem item : m_comboboxSelectedGeneAnnotation.getItems())
			{
				if(item.getValue().hashCode() == nReferenceID)
				{
					m_comboboxSelectedGeneAnnotation.setSelectedItem(item);
					m_strFileGTF = item.getValue();
					LoadGeneAnnotation(m_strFileGTF);
					break;
				}
			}
		}
		
		// select first available condition type and condition
		UpdateComboboxesForCondition();
		if(m_comboboxSelectedConditionType.getItemCount() > 0)
		{
			m_strSelectedConditionType = m_projectModel.GetFirstConditionType();
//					m_comboboxSelectedConditionType.setSelectedIndex(0);
//					m_strSelectedConditionType = m_comboboxSelectedConditionType.getSelectedItem().getLabel();
			
			for(Comboitem item : m_comboboxSelectedConditionType.getItems())
			{
				if(item.getLabel().equals(m_strSelectedConditionType))
				{
					m_comboboxSelectedConditionType.setSelectedItem(item);
					break;
				}
			}
			
			OnConditionTypeChange();
			
			// select first available condition
			m_comboboxSelectedCondition.setSelectedIndex(0);
			m_strSelectedCondition = m_comboboxSelectedCondition.getSelectedItem().getLabel();
		}
		
		// update colors
		UpdateColorSelection();

		// clear the gene selection
		m_bandboxSelectedGene.setText("Select Gene");
		
		// clear previous AS results and data in the data supplier
		m_vcASResults.Clear();
		m_dataSupplier.Clear();		
		
		//###########################################################################
		//    load previously calculated alternative splicing hits, if available
		//###########################################################################
		String strResultFile = m_projectModel.GetProjectName() + "_splicing_results.dat";
		pIn = new File(m_strPathHitLists + "/" + strResultFile);
		if(pIn.exists())
		{
			m_vcASResults.LoadHitList(m_strPathHitLists, strResultFile);			
			m_resultListHandler.UpdatePermanentResultList(m_vcASResults);
		}
		else
		{
			m_resultListHandler.Clear();
		}
		
		OnGeneChange(true, true);
	}
	
	/**
	 *    This function is triggered when a new condition type is selected.
	 *    It updates the available conditions in the condition combo box. 
	 */	
	public void OnConditionTypeChange()
	{
		//############################################
		//       update condition selection box
		//############################################
		
		//remove old conditions
		m_comboboxSelectedCondition.getItems().clear();
		
		// add conditions
		TreeSet<String> vcConditions = m_projectModel.GetConditions(m_strSelectedConditionType);
		for(String strCondition : vcConditions)
		{
			Comboitem item = new Comboitem(strCondition);
			item.setValue(strCondition);
			m_comboboxSelectedCondition.appendChild(item);
		}
		m_plotFactory.RequestCoverageRedraw();
				
		m_plotFactory.GetColorsPerCondition();
		AddComboboxForSamples(null);
	}

	/**
	 *    This function is triggered when a new gene is selected. It clears data associated with the previous gene,
	 *    processes the new gene information and updates some GUI elements.
	 */
	public void OnGeneChange(boolean bProjectChanged, boolean bIdentifyASExons) throws Exception
	{
		if(!m_projectModel.IsReady())
		{
			Messagebox.show("Did you forget to select a project?");
			return;
		}
		
		if(m_strSelectedCondition == null || m_strSelectedCondition.equals("Select Condition"))
		{
			if(!bProjectChanged)
				Messagebox.show("Did you forget to select a condition?");
			return;
		}
		
		// clear old gene data
		if(m_vcValidIsoforms == null)
			m_vcValidIsoforms = new TreeSet<String>();
		else
			m_vcValidIsoforms.clear();
		
		// clear AS results
		m_vcASResults.ClearTemporaryResults();
		
		// clear highlighted results
		m_selectedResult = null;
		
		if(m_mapCountsPerExonAndTissue != null)
			m_mapCountsPerExonAndTissue.clear();
		else
			m_mapCountsPerExonAndTissue = new TreeMap<String, TreeMap<String, Vector<Double>>>();
		
		if(m_mapEntropyToExonicPart != null)
			m_mapEntropyToExonicPart.clear();
		else
			m_mapEntropyToExonicPart = new TreeMap<String, Double>();
		
		m_plotFactory.RequestCoverageRedraw();
		
		if(m_bandboxSelectedGene.getText().equals("Select Gene"))
			return;

		GeneIdentifier gid = m_geneIdentifierHandler.GetGeneIdentifierForGene(m_bandboxSelectedGene.getText(), null);
		if(gid == null)
		{
			Messagebox.show("Gene not found: " + m_bandboxSelectedGene.getText());
			return;
		}

		ProcessGeneInformation(gid, m_strFileGTF, null);
		
		if(m_bShowEntropyData)
		{
			Events.postEvent("onSelect", m_comboboxSelectedEntropyIndex, null);
		}
		
		m_plotFactory.DrawPlots();

		InitComboboxForIsoformSelection();
	}

	/**
	 *    This function is triggered when the isoforms change. If a new gene was selected it updates the result list
	 *    of alternative splicing events, or else, it invokes an update of the data stored in the data supplier to reflect
	 *    the data associated with the remaining isoforms.
	 */
	public boolean OnIsoformChange(boolean bNewGene) throws IOException, ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		m_plotFactory.ClearUnknownJunctions();
		
		//###########################################################################
		//                  identify alternatively spliced exons
		//###########################################################################
		if(bNewGene)
		{
			AnalyzerGeneFactory analyzer = new AnalyzerGeneFactory();
			TreeSet<AnalysisResult> vcResults = analyzer.AnalyzeCurrentGene(m_hWindow, m_bSkipFirstAndLastExon, false, 0, false);
			m_vcASResults.ClearTemporaryResults();
			m_vcASResults.AddTemporaryResultsForCurrentlySelectedIsoforms(vcResults);

			m_resultListHandler.UpdateTemporaryResultList(m_vcASResults);
		}
		else
		{
			// isoforms have been removed, thus recalculate the exon groups and group coverage values
			m_dataSupplier.RetrieveCoverageData(m_hWindow);
		}

		return true;
	}
	
	/** Helper function to remove a single isoform from the isoform selection. */
	public void RemoveIsoformFromSelection(String strIsoform)
	{		
		for(Treeitem item : m_treeSelectedIsoforms.getItems())
		{
			if(item.getLabel().equals(strIsoform.split("\\.")[0]))
			{
				m_treeSelectedIsoforms.removeItemFromSelection(item);
			}
		}
		
		m_bIsoformSelectionChanged = true;
		m_vcValidIsoforms.remove(strIsoform);
		
		// redraw isoform plot
		try
		{
			OnIsoformChange(false);
			m_plotFactory.RequestCoverageRedraw();
			m_plotFactory.DrawPlots();
		}
		catch(Exception e)
		{
			System.out.println("failed to draw isoform plot");
			e.printStackTrace();
		}
	}
	
	/** Prepares the option menu */
	public void AddSettingsMenu(Vlayout parentLayout) throws Exception
	{
		Groupbox grpBox = new Groupbox();
		grpBox.setTitle("Settings");
		grpBox.setMold("3d");
		grpBox.setParent(parentLayout);
		grpBox.setWidth("530px");
		grpBox.setStyle("margin-left: 10px;");
		grpBox.setClosable(false);
		
		Hlayout layoutH = new Hlayout();
		layoutH.setParent(grpBox);
		
		Vlayout layoutV = new Vlayout();
		layoutV.setParent(layoutH);

		// add label for window size edit box
		Label lab = new Label("Window width in pixel:");
		lab.setStyle("margin-left: 10px");
		lab.setParent(layoutV);
		
		Hlayout layoutH2 = new Hlayout();
		layoutH2.setParent(layoutV);
		
		// add box for window size specification
		m_textboxWindowWidth = new Textbox("Enter value");
		m_textboxWindowWidth.setParent(layoutH2);
		m_textboxWindowWidth.setWidth("177px");
		m_textboxWindowWidth.setStyle("margin-left: 10px;");
		
		m_textboxWindowWidth.addEventListener(Events.ON_CHANGING, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bWindowSizeChanged = true;
			}
		});
		
		m_textboxWindowWidth.addEventListener(Events.ON_BLUR, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				if(!m_bWindowSizeChanged)
					return;
				
				m_bWindowSizeChanged = false;
				
				Textbox box = (Textbox)event.getTarget();
				String strValue = box.getValue();
				if(strValue != null)
				{
					try
					{
						m_plotFactory.SetMaxWidth(Integer.parseInt(strValue));
						m_plotFactory.RequestCoverageRedraw();
						
						if(m_dataSupplier.GetGene() != null)
							m_plotFactory.DrawPlots();
					}
					catch(NumberFormatException e)
					{
						if(m_dataSupplier.GetGene() == null)
						{
							Messagebox.show("Invalid input: " + strValue + ". The window width must be a valid number.");
							return;
						}
					}
				}
			}
		});
		
		m_textboxWindowWidth.addEventListener(Events.ON_OK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				if(!m_bWindowSizeChanged)
					return;
				
				m_bWindowSizeChanged = false;
				
				Textbox box = (Textbox)event.getTarget();
				String strValue = box.getValue();
				if(strValue != null)
				{
					try
					{
						m_plotFactory.SetMaxWidth(Integer.parseInt(strValue));
						m_plotFactory.RequestCoverageRedraw();
						
						if(m_dataSupplier.GetGene() != null)
							m_plotFactory.DrawPlots();
					}
					catch(NumberFormatException e)
					{
						if(m_dataSupplier.GetGene() == null)
						{
							Messagebox.show("Invalid input: " + strValue + ". The window width must be a valid number.");
							return;
						}
					}
				}
			}
		});

		Image btnResetWindowWidth = new Image("/img/red_cross.png");
		btnResetWindowWidth.setWidth("18px");
		btnResetWindowWidth.setHeight("18px");
		btnResetWindowWidth.setStyle("margin-top: 4px; cursor:pointer;");
		btnResetWindowWidth.setParent(layoutH2);
		
		btnResetWindowWidth.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{	
				int pSize[] = m_plotFactory.GetClientSize();
				m_plotFactory.SetMaxWidth(pSize[0]);
				m_plotFactory.RequestCoverageRedraw();
						
				if(m_dataSupplier.GetGene() != null)
					m_plotFactory.DrawPlots();

				m_textboxWindowWidth.setText("" + pSize[0]);
			}
		});
				
		layoutV = new Vlayout();
		layoutV.setStyle("margin-left: 40px");
		layoutV.setParent(layoutH);
		
		// add label for project selection
		lab = new Label("Selected Project:");
		lab.setParent(layoutV);
		
		// add bandbox for project selection
		m_bandboxProjects = new Bandbox("Select Project");
		m_bandboxProjects.setParent(layoutV);
		m_bandboxProjects.setWidth("200px");
		
		// add bandpopup
		m_bandpopupProjects = new Bandpopup();
		m_bandpopupProjects.setParent(m_bandboxProjects);
		m_bandpopupProjects.setWidth("75%");
		m_bandpopupProjects.setHeight("75%");
		
		// create tree set of project information
		File pFolder = new File(m_strPathInput);
		for(File pFile : pFolder.listFiles())
		{
			String strFile = pFile.getName();
			
			if(strFile.endsWith("project"))// || strFile.endsWith("project2"))
			{
				ProjectModel project = new ProjectModel();
				project.Init(m_strPathInput + "/" + strFile, 0, 0, 0.0, 0, false, false);
				m_vcProjectInfos.add(project);
			}
		}
		
		m_listboxProjects = new Listbox();
		m_listboxProjects.setParent(m_bandpopupProjects);		
		m_listboxProjects.setHflex("min");
		m_listboxProjects.setStyle("autoWidth:true");

		Listhead listheadProject = new Listhead();
		listheadProject.setSizable(true);
		listheadProject.setParent(m_listboxProjects);
		
		Listheader listheaderID = new Listheader("Project ID");
		listheaderID.setWidth("200px");
		listheaderID.setHflex("min");		
		listheaderID.setParent(listheadProject);
		
		Listheader listheaderProject = new Listheader("Project name");
		listheaderProject.setWidth("200px");
		listheaderProject.setHflex("min");		
		listheaderProject.setParent(listheadProject);
		
		listheaderProject = new Listheader("Creation date");
		listheaderProject.setWidth("200px");
		listheaderProject.setHflex("min");
		listheaderProject.setParent(listheadProject);
		
		listheaderProject = new Listheader("#samples");
		listheaderProject.setWidth("200px");
		listheaderProject.setHflex("min");
		listheaderProject.setParent(listheadProject);
		
		listheaderProject = new Listheader("paired data");
		listheaderProject.setWidth("60px");
		listheaderProject.setHflex("min");
		listheaderProject.setParent(listheadProject);
		
		listheaderProject = new Listheader("Conditions");
		listheaderProject.setWidth("200px");
		listheaderProject.setHflex("min");
		listheaderProject.setParent(listheadProject);
	
		m_listboxProjects.addEventListener(Events.ON_SELECT, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Listbox box = (Listbox)event.getTarget();

				if(box.getSelectedItem().getValue() != null)
				{
					ProjectModel model = box.getSelectedItem().getValue();
					OnProjectChange(model.GetFullPathOfProjectFile(), false);
					Bandpopup popup = (Bandpopup)box.getParent();
					Bandbox bandbox = (Bandbox)popup.getParent();					
					bandbox.setText(model.GetProjectName());
//					bandbox.setOpen(false);
					bandbox.close();
				}
			}
		});
		
		// fill the first page with project data
		int nCurProject = 0;
		for(ProjectModel project : m_vcProjectInfos)
		{
			if(nCurProject == PROJECT_PAGING_SIZE)
				break;
			
			Listitem item = new Listitem();
			item.setValue(project);
			
			Listcell cell = new Listcell(""+project.GetProjectName().hashCode()); cell.setParent(item);
			cell = new Listcell(project.GetProjectName()); cell.setParent(item);
			cell = new Listcell(project.GetCreationDate()); cell.setParent(item);
			cell = new Listcell(""+project.GetSamples().size()); cell.setParent(item);
			
			if(project.ProjectHasPairedData())
			{
				cell = new Listcell("yes"); cell.setParent(item);
			}
			else
			{
				cell = new Listcell("no"); cell.setParent(item);
			}
			
			int nIdx = 0;
			cell = new Listcell();			
			for(String strConditionType : project.GetConditionsToConditionTypes().keySet())
			{
				if(nIdx != 0)
				{
					Separator sep = new Separator();
					sep.setParent(cell);
				}
				String strLabel = strConditionType + ": " + project.GetConditionsToConditionTypes().get(strConditionType);
				lab = new Label(strLabel); lab.setParent(cell);
				nIdx++;
			}

			cell.setParent(item);
//			""+project.GetConditionsToConditionTypes()); cell.setParent(item);
			
			m_listboxProjects.appendChild(item);
			
			nCurProject++;
		}
		
		// Add paging component
		m_pagingProjects = new Paging();
		m_pagingProjects.setParent(m_bandpopupProjects);
		m_pagingProjects.setDetailed(true);

		m_pagingProjects.addEventListener("onPaging", new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				PagingEvent pe = (PagingEvent) event;
				int nCurrentPage = pe.getActivePage();

				int nStartItem = nCurrentPage * PROJECT_PAGING_SIZE;
				m_listboxProjects.getItems().clear();
				
				// fill the current page with project data
				int nCurProject = 0;
				for(ProjectModel project : m_vcProjectInfos)
				{
					if(nCurProject < nStartItem)
					{
						nCurProject++;
						continue;
					}
					
					if(nCurProject == nStartItem+PROJECT_PAGING_SIZE)
						break;
					
					Listitem item = new Listitem();
					item.setValue(project);
					
					Listcell cell = new Listcell(""+project.GetProjectName().hashCode()); cell.setParent(item);
					cell = new Listcell(project.GetProjectName()); cell.setParent(item);
					cell = new Listcell(project.GetCreationDate()); cell.setParent(item);
					cell = new Listcell(""+project.GetSamples().size()); cell.setParent(item);
					
					if(project.ProjectHasPairedData())
					{
						cell = new Listcell("yes"); cell.setParent(item);
					}
					else
					{
						cell = new Listcell("no"); cell.setParent(item);
					}
					
					int nIdx = 0;
					cell = new Listcell();			
					for(String strConditionType : project.GetConditionsToConditionTypes().keySet())
					{
						if(nIdx != 0)
						{
							Separator sep = new Separator();
							sep.setParent(cell);
						}
						String strLabel = strConditionType + ": " + project.GetConditionsToConditionTypes().get(strConditionType);
						Label lab = new Label(strLabel); lab.setParent(cell);
						nIdx++;
					}

					cell.setParent(item);
					
					m_listboxProjects.appendChild(item);
					
					nCurProject++;
				}
				
				// some dirty code to enforce a maximum width on the bandbox
				m_bandpopupProjects.setHflex("min");
				String strMaxWidth = "max-width: " + (m_plotFactory.GetClientSize()[0] * 0.75) + "px;";
				m_bandpopupProjects.setStyle(strMaxWidth);
				m_bandpopupProjects.setZclass(strMaxWidth);

				Clients.resize(m_bandpopupProjects);
				
				// also adjust the paging object
//				m_pagingProjects.setHflex("min");
				Clients.resize(m_pagingProjects);
			}
		});
		
		m_bandboxProjects.addEventListener(Events.ON_OK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Bandbox box = (Bandbox)event.getTarget();
				if(box.getValue() != null && !box.getValue().isEmpty())
				{
					OnProjectChange(box.getValue(), false);
				}
			}
		});

		//#####################################
		//       condition selection
		//#####################################
		layoutH = new Hlayout();
		layoutH.setStyle("margin-left: 10px");
		layoutH.setParent(grpBox);
		
		AddComboboxesForCondition(layoutH);
		
		//#####################################
		//     gene identifier selection
		//#####################################
		layoutH = new Hlayout();
		layoutH.setStyle("margin-top: 10px");
		layoutH.setParent(grpBox);
		
		layoutV = new Vlayout();
		layoutV.setStyle("margin-left: 10px");
		layoutV.setParent(layoutH);
		
		// add label for gene identifier selection
		lab = new Label("Selected Gene:");
		lab.setParent(layoutV);
		
		// add text box for gene identifier selection
		m_bandboxSelectedGene = new Bandbox();
		m_bandboxSelectedGene.setParent(layoutV);
		m_bandboxSelectedGene.setMold("rounded");
//		m_bandboxSelectedGene.setStyle("padding-left: 10px;");
		m_bandboxSelectedGene.setAutodrop(true);
		m_bandboxSelectedGene.setVflex("true");
		m_bandboxSelectedGene.setWidth("200px");
		
		Bandpopup popup = new Bandpopup();
		popup.setParent(m_bandboxSelectedGene);
		
		m_listboxSelectedGene = new Listbox();
		m_listboxSelectedGene.setParent(popup);
		m_listboxSelectedGene.setHeight("300px");
		m_listboxSelectedGene.setMold("paging");
		m_listboxSelectedGene.setVflex(true);
		m_listboxSelectedGene.setAutopaging(true);
		m_listboxSelectedGene.setSizedByContent(true);
		
		Listhead listhead = new Listhead();
		Listheader listheader = new Listheader("gene symbol"); listheader.setParent(listhead);
		listheader = new Listheader("Ensembl ID"); listheader.setParent(listhead);
		listheader = new Listheader("RefSeq ID"); listheader.setParent(listhead);
		listheader = new Listheader("Entrez Gene"); listheader.setParent(listhead);
		listheader = new Listheader("Synonyms"); listheader.setParent(listhead);
		listheader = new Listheader("Description"); listheader.setParent(listhead);
		
		listhead.setParent(m_listboxSelectedGene);

		layoutV = new Vlayout();
		layoutV.setStyle("margin-left: 40px");
		layoutV.setParent(layoutH);
		
		// add label for project selection
		lab = new Label("Selected Gene Annotation:");
		lab.setParent(layoutV);
		
		m_comboboxSelectedGeneAnnotation = new Combobox("Select Gene Annotation");
		m_comboboxSelectedGeneAnnotation.setParent(layoutV);
		m_comboboxSelectedGeneAnnotation.setWidth("200px");
		pFolder = new File(m_strPathReferences);
		int nIdx = -1;
		int i = 0;
		for(File pFile : pFolder.listFiles())
		{
			String strFile = pFile.getName();
			
			if(strFile.toLowerCase().endsWith("gtf") || strFile.toLowerCase().endsWith("gff") || strFile.toLowerCase().endsWith("refflat"))
			{
				Comboitem item = new Comboitem(strFile);
				item.setValue(pFile.getAbsolutePath());
				m_comboboxSelectedGeneAnnotation.appendChild(item);
				
				if(nIdx == -1 && strFile.toLowerCase().contains("gencode"))
					nIdx=i;
				i++;
			}
		}
		// select gencode gene annotation if present, otherwise select first available gene annotation
		if(nIdx != -1)
			m_comboboxSelectedGeneAnnotation.setSelectedIndex(nIdx);
		else
			m_comboboxSelectedGeneAnnotation.setSelectedIndex(0);
				
		m_strFileGTF = m_comboboxSelectedGeneAnnotation.getSelectedItem().getValue();
		LoadGeneAnnotation(m_strFileGTF);
		
		m_comboboxSelectedGeneAnnotation.addEventListener(Events.ON_SELECT, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Combobox box = (Combobox)event.getTarget();
				if(box.getSelectedItem() != null)
				{
					String strText = (String)box.getSelectedItem().getValue();
					m_strFileGTF = strText;
					LoadGeneAnnotation(m_strFileGTF);
					
					// save last used gene annotation for project
					if(m_projectModel == null || !m_projectModel.IsReady())
						return;
					
					// save reference ID to file. The program will try to load the last reference
					// used on changing the project.
					String strOut = m_projectModel.GetFullPathOfProjectFile() + ".last_reference";
					PrintWriter pOut = new PrintWriter(new File(strOut));
					pOut.println(m_strFileGTF.hashCode());
					pOut.close();
				}
			}
		});
		
		m_bandboxSelectedGene.addEventListener(Events.ON_OK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				OnGeneChange(false, true);
			}
		});
		
		m_bandboxSelectedGene.addEventListener(Events.ON_CHANGING, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
//				Bandbox box = (Bandbox)event.getTarget();
				org.zkoss.zk.ui.event.InputEvent inputEvent=(org.zkoss.zk.ui.event.InputEvent)event;

				String strInput = inputEvent.getValue();
				strInput.trim();

				// clear old items
				m_listboxSelectedGene.getItems().clear();

				// get first letter
				if(strInput != null && strInput.length() > 0)
				{
					TreeSet<GeneIdentifier> vcGeneIdentifiers = new TreeSet<GeneIdentifier>();

					Vector<GeneIdentifier> vcIdentifiers = m_geneIdentifierHandler.GetAllUniqueGeneIdentifiers();

					// get list of valid gene identifiers
					for(GeneIdentifier gid : vcIdentifiers)
					{
						if(gid.m_bIsValid)
						{
							for(int i=0; i<4; i++)
							{
								String strID = "?";
								
								switch(i)
								{
									case 0: strID = gid.m_strApprovedGeneSymbol; break;
									case 1: strID = gid.m_strEnsemblGeneID; break;
									case 2: strID = gid.m_strRefGeneID; break;
									case 3: strID = gid.m_strEntrezGeneID; break;
								}
								
								if(strID.startsWith(strInput))
								{
									vcGeneIdentifiers.add(gid);
								}
							}
						}
						
						if(gid.m_strSynonyms.contains(strInput))
						{
							vcGeneIdentifiers.add(gid);
						}
					}
					
					// sort gene identifiers by:
					// #1 complete matches
					// #2 same start string
					// #3 synonyms
					
					TreeSet<GeneIdentifier> vcFullMatch 	= new TreeSet<GeneIdentifier>();
					TreeSet<GeneIdentifier> vcPartialMatch  = new TreeSet<GeneIdentifier>();
					TreeSet<GeneIdentifier> vcSynonyms		= new TreeSet<GeneIdentifier>();
					
					for(GeneIdentifier gid : vcGeneIdentifiers)
					{
						if(gid.m_strApprovedGeneSymbol.startsWith(strInput))
						{
							if(gid.m_strApprovedGeneSymbol.equals(strInput))
								vcFullMatch.add(gid);
							else
								vcPartialMatch.add(gid);
							continue;
						}
						
						if(gid.m_strSynonyms.contains(strInput))
						{
							boolean bFullMatch = false;
							String pSplit[] = gid.m_strSynonyms.split(", ");
							for(String strName : pSplit)
							{
								// separate matching synonyms from the others
								if(strName.equals(strInput))
								{
									vcFullMatch.add(gid);
									bFullMatch = true;
								}
							}
							
							if(!bFullMatch)
								vcSynonyms.add(gid);
						}
					}

					// add 'real' matches first, then use synonyms
					for(int i=0; i<3; i++)
					{
						switch(i)
						{
							case 0: vcGeneIdentifiers = vcFullMatch; break;
							case 1: vcGeneIdentifiers = vcPartialMatch; break;
							case 2: vcGeneIdentifiers = vcSynonyms; break;
						}
						
						// now add the sorted list of gene identifiers to the selection
						for(GeneIdentifier gid : vcGeneIdentifiers)
						{
							// valid items must have the correct gene symbol or synonym
							boolean bIsValid = false;
							
							Listitem item = new Listitem();
							Listcell cell = new Listcell(gid.m_strApprovedGeneSymbol); cell.setParent(item);
							if(gid.m_strApprovedGeneSymbol.startsWith(strInput))
							{
								cell.setStyle("font-weight: bold; color: blue");
								bIsValid = true;
							}
							item.setValue(gid.m_strEnsemblGeneID);
							cell = new Listcell(gid.m_strEnsemblGeneID); cell.setParent(item);						
							cell = new Listcell(gid.m_strRefGeneID); cell.setParent(item);
							cell = new Listcell(gid.m_strEntrezGeneID); cell.setParent(item);
							if(gid.m_strSynonyms.contains(strInput))
							{
								cell = new Listcell();
								cell.setParent(item);
								
								Vector<String> vcHit = new Vector<String>();
								Vector<String> vcOther = new Vector<String>();
								String pSplit[] = gid.m_strSynonyms.split(", ");
								
								for(String strName : pSplit)
								{
									// separate matching synonyms from the others
									if(strName.startsWith(strInput))
									{
										vcHit.add(strName);
										bIsValid = true;
									}
									else
									{
										vcOther.add(strName);
									}
								}
								
								String strText = "<b><span style=\"color: blue\">";
								boolean bFirst = true;
								for(String strName : vcHit)
								{
									if(!bFirst)
										strText += ", ";
	
									strText += strName;
								}
								
								 strText += "</span></b>";
								
								for(String strName : vcOther)
								{
									strText += ", " + strName;
								}
	
								Html html = new Html(strText);
								html.setParent(cell);
							}
							else
							{
								cell = new Listcell(gid.m_strSynonyms); cell.setParent(item);
							}
							cell = new Listcell(gid.m_strApprovedGeneName); cell.setParent(item);
							
							if(bIsValid)
								m_listboxSelectedGene.appendChild(item);
						}
					}
				}
			}
		});
			
		m_listboxSelectedGene.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Listbox box = (Listbox)event.getTarget();
				
				// skip the rest if the selection has not changed
				if(box.getSelectedItem() == null)
					return;
				
				String strText = box.getSelectedItem().getLabel();
				String strValue = (String)box.getSelectedItem().getValue();
				
				m_bandboxSelectedGene.close();				
				m_bandboxSelectedGene.setText(strText);
				m_bandboxSelectedGene.setValue(strValue);

				OnGeneChange(false, true);
			}
		});

		//##################################################
		
		layoutH = new Hlayout();
		layoutH.setStyle("margin-top: 10px");
		layoutH.setParent(grpBox);
		
		Radiogroup rgroup = new Radiogroup();
		
		Radio r = new Radio("human");
		r.setParent(rgroup);
		r.setChecked(true);
		
		r = new Radio("mouse");
		r.setParent(rgroup);

		rgroup.addEventListener(Events.ON_CHECK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Radio btn = (Radio)event.getTarget();
				
				if(btn.getLabel().equals("human"))
				{
					String strFileCrossReference = m_strPathReferences + "/biomart_cross_ref_human.txt";
					m_geneIdentifierHandler.Init(strFileCrossReference);
					LoadGeneAnnotation(m_strFileGTF);
				}
				else if(btn.getLabel().equals("mouse"))
				{
					String strFileCrossReference = m_strPathReferences + "/biomart_cross_ref_mouse.txt";
					m_geneIdentifierHandler.Init(strFileCrossReference);
					LoadGeneAnnotation(m_strFileGTF);
				}
			}
		});

		rgroup.setParent(layoutH);
	}

	/** Prepares the options menu */
	public void AddOptionsMenu(Vlayout parentLayout)
	{
		Groupbox grpBox = new Groupbox();
		grpBox.setTitle("Options");
		grpBox.setMold("3d");
		grpBox.setParent(parentLayout);
		grpBox.setWidth("530px");
		grpBox.setStyle("margin-left: 10px;");
		grpBox.setVflex("min");
		grpBox.setClosable(false);
	
		Hlayout layoutH = new Hlayout();
		layoutH.setVflex("min");
		layoutH.setParent(grpBox);
		
		Vlayout layoutV = new Vlayout();
		layoutV.setHflex("min");
		layoutV.setParent(layoutH);
		
		Label label = new Label("Isoform view options");
		label.setStyle("text-decoration:underline");
		label.setParent(layoutV);
		
		m_checkboxHighlightUniqueExons = new Checkbox("Show unique isoform elements");
		m_checkboxHighlightUniqueExons.setParent(layoutV);
		m_checkboxHighlightUniqueExons.setSclass("checkbox");
		m_checkboxHighlightUniqueExons.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				if(m_dataSupplier.GetGene() == null)
				{
					Messagebox.show("No gene selected!");
					m_checkboxHighlightUniqueExons.setChecked(false);
					return;
				}
				m_plotFactory.ShowUniqueFeatures(m_checkboxHighlightUniqueExons.isChecked());			
				m_plotFactory.DrawIsoforms();
			}
		});
		
		Checkbox checkBox = new Checkbox("Show tissue specificity index");
		checkBox.setParent(layoutV);
		checkBox.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bShowEntropyData = !m_bShowEntropyData;
				
				if(!m_bShowEntropyData)
				{
					m_plotFactory.DrawIsoforms();
					return;
				}
				
				m_mapCountsPerExonAndTissue = GetGTEXDataForExonicParts();
				EntropyFactory entropyFactory = new EntropyFactory();
				
				switch(m_comboboxSelectedEntropyIndex.getSelectedItem().getLabel())
				{
					case "None":
						m_bShowEntropyData = false;
						break;
						
					case "Gini":
						m_mapEntropyToExonicPart = entropyFactory.CalculateGiniIndexForExons(m_mapCountsPerExonAndTissue);
						m_bShowEntropyData = true;
						break;
						
					case "Simpson":
						m_mapEntropyToExonicPart = entropyFactory.CalculateGiniSimpsonIndexForExons(m_mapCountsPerExonAndTissue);
						m_bShowEntropyData = true;
						break;
						
					case "Shannon":
						m_mapEntropyToExonicPart = entropyFactory.CalculateShannonIndexForExons(m_mapCountsPerExonAndTissue);
						m_bShowEntropyData = true;
						break;
						
					case "Theil":
						m_mapEntropyToExonicPart = entropyFactory.CalculateTheilIndexForExons(m_mapCountsPerExonAndTissue);
						m_bShowEntropyData = true;
						break;
				}
				
				// force isoform redraw
				m_plotFactory.DrawIsoforms();
			}
		});
		if(!m_bGTEXAvailableForExons)
		{
			checkBox.setDisabled(true);
		}
		m_checkboxColorExonsAndJunctionsByCoverage = new Checkbox("Color exons and junctions by expression level");
		m_checkboxColorExonsAndJunctionsByCoverage.setParent(layoutV);
		m_checkboxColorExonsAndJunctionsByCoverage.setChecked(true);
		m_checkboxColorExonsAndJunctionsByCoverage.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_plotFactory.ShowColoredExonsAndJunctions(m_checkboxColorExonsAndJunctionsByCoverage.isChecked());
				m_plotFactory.DrawPlots();
			}
		});
		
		label = new Label("Coverage view options");
		label.setStyle("margin-top: 10px; display:inline-block; text-decoration:underline");
		label.setParent(layoutV);
		
		m_checkboxShowRelativeCoverage = new Checkbox("Show relative coverage");
		m_checkboxShowRelativeCoverage.setWidth("100%");
		m_checkboxShowRelativeCoverage.setParent(layoutV);
		m_checkboxShowRelativeCoverage.setChecked(m_plotFactory.IsRelativeCoverageEnabled());
		m_checkboxShowRelativeCoverage.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_plotFactory.ShowRelativeCoverage(m_checkboxShowRelativeCoverage.isChecked());
			}
		});	
		
		m_checkboxUseLog2 = new Checkbox("Show log2 transformed");
		m_checkboxUseLog2.setWidth("100%");
		m_checkboxUseLog2.setParent(layoutV);
		m_checkboxUseLog2.setChecked(m_plotFactory.IsLog2Enabled());
		m_checkboxUseLog2.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_plotFactory.ShowAsLog2(m_checkboxUseLog2.isChecked());
			}
		});
		
		label = new Label("Isoform direction");
		label.setStyle("margin-top: 10px; display:inline-block; text-decoration:underline");
		label.setParent(layoutV);
		
		m_checkboxSwitchStrand = new Checkbox("swap strand automatically");
		m_checkboxSwitchStrand.setWidth("100%");
		m_checkboxSwitchStrand.setParent(layoutV);
		m_checkboxSwitchStrand.setChecked(m_plotFactory.IsRelativeCoverageEnabled());
		m_checkboxSwitchStrand.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_plotFactory.ToggleSwitchStrand(m_checkboxSwitchStrand.isChecked());
				m_plotFactory.RequestCoverageRedraw();
				m_plotFactory.DrawPlots();
			}
		});	
		
		Hbox layoutH2 = new Hbox();
		layoutH2.setParent(layoutV);
		layoutH2.setStyle("width: 100%; margin-top: 10px;");
		
		Image img = new Image();
		img.setSrc("./img/help.png");
		img.setParent(layoutH2);
		img.setWidth("26px");
		img.setStyle("display: inline-block; cursor:pointer;");
		img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Executions.getCurrent().sendRedirect("./Manual/index.html", "_blank");
			}
		});
		
		Label labHelp = new Label("Manual");
		labHelp.setParent(layoutH2);
		labHelp.setStyle("display: inline-block; font-size: 16px; font-weight: 600; color: hsl(211, 50%, 44%); margin-top: 11px; cursor: pointer;");
		labHelp.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Executions.getCurrent().sendRedirect("./Manual/index.html", "_blank");
			}
		});
		
		layoutV = new Vlayout();
//		layoutV.setHflex("true");
		layoutV.setHeight("262px");
		layoutV.setParent(layoutH);
		
		//################################################
		//    Add button to remove irrelevant isoforms
		//################################################
		Button btn = new Button("Show relevant isoforms only");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.setSclass("button blue");
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
			
				if(m_dataSupplier.GetGene() == null)
				{
					Messagebox.show("No gene specified. Please confirm the gene you entered by pressing 'Enter' or select one gene from the list of genes offered while typing.");
					return;
				}
				
				HideIrrelevantIsoforms(m_bSkipFirstAndLastExon, false);
				m_plotFactory.RequestCoverageRedraw();

				OnIsoformChange(false);
				m_plotFactory.DrawPlots();
			}
		});
		
		//############################################
		//     Add button to view GTEX expression
		//############################################		
		btn = new Button("Show GTEX data for gene");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.setSclass("button blue");
		if(!m_bGTEXAvailableForGenes)
			btn.setDisabled(true);
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				if(m_dataSupplier.GetGene() == null)
				{
					Messagebox.show("No gene specified. Please confirm the gene you entered by pressing 'Enter' or select one gene from the list of genes offered while typing.");
					return;
				}
				
				m_plotFactoryGTEX.ShowGTEXDataForGene(m_dataSupplier.GetGeneID());
			}
		});

		//############################################
		//     Add button to show MMSeq estimates
		//############################################
		btn = new Button("Calculate MMSeq estimates");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.setSclass("button blue");
		if(!m_bMMSeqAvailable)
			btn.setDisabled(true);
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				if(m_dataSupplier.GetGene() == null)
				{
					Messagebox.show("No gene specified");
					return;
				}
				
				GenerateMMSeqEstimates();
				m_plotFactory.DrawPlots();
			}
		});

		//############################################
		//     Add button to save a screenshot
		//############################################
		btn = new Button("Save screenshot");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.setSclass("button blue");
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				GeneIdentifier gid = m_dataSupplier.GetGeneIdentifier();

				BufferedImage imgCoverage = ImageIO.read(m_imgMapCoverage.getContent().getStreamData());
				BufferedImage imgIsoforms = ImageIO.read(m_imgMapIsoforms.getContent().getStreamData());
				
				int nWidth = Math.max(imgCoverage.getWidth(), imgIsoforms.getWidth());
				int nHeight = imgCoverage.getHeight() + imgIsoforms.getHeight();
				BufferedImage combined = new BufferedImage(nWidth, nHeight, BufferedImage.TYPE_INT_ARGB);
				
				Graphics g = combined.getGraphics();
				g.drawImage(imgCoverage, 0, 0, null);
				g.drawImage(imgIsoforms, 0, imgCoverage.getHeight(), null);

				// Save as new image
				ByteArrayOutputStream os = new ByteArrayOutputStream();
				ImageIO.write(combined, "png", os);
				InputStream is = new ByteArrayInputStream(os.toByteArray());

				Filedownload.save(is, "image/png", gid.m_strApprovedGeneSymbol + ".png");
			}
		});
		
		//############################################
		//     Add button to calculate PSI scores
		//############################################
		btn = new Button("Show splicing heat map");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.setSclass("button blue");
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event arg0) throws Exception
			{
				m_plotFactory.DrawJunctionHeatmap();
			}
		});

		//#############################################################
		//    Add button to reanalyze using currently selected data
		//#############################################################
		btn = new Button("Reanalyze gene");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.setSclass("button orange");
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				if(m_dataSupplier.GetGene() == null)
				{
					Messagebox.show("No gene specified");
					return;
				}

				OnIsoformChange(true);
				m_plotFactory.DrawPlots();
			}
		});
		
		//############################################
		//     Add button to redraw isoforms
		//############################################
		btn = new Button("Redraw plots");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.setSclass("button orange");
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				if(m_dataSupplier.GetGene() == null)
				{
					Messagebox.show("No gene specified");
					return;
				}
								
				if(m_bIsoformSelectionChanged)
				{
					OnIsoformChange(false);
				}
				
				if(m_bSampleSelectionChanged)
				{
					m_dataSupplier.RetrieveCoverageData(m_hWindow);
				}
				
				m_selectedResult = null;
				m_plotFactory.RequestCoverageRedraw();
				m_plotFactory.DrawPlots();
			}
		});

	}
	
	/** Prepares the advanced options menu */
	public void AddAdvancedOptionsMenu(Vlayout parentLayout)
	{
		Groupbox grpBox = new Groupbox();
		grpBox.setTitle("Advanced Options");
		grpBox.setMold("3d");
		grpBox.setParent(parentLayout);
		grpBox.setWidth("440px");
		grpBox.setHeight("490px");
		
		//#############################################################
		Hlayout layoutH = new Hlayout();
		layoutH.setParent(grpBox);
		
		Vlayout layoutV = new Vlayout();
		layoutV.setWidth("220px");
		layoutV.setParent(layoutH);
		
		m_checkboxCoverageGrid = new Checkbox("Show coverage plot grid");
		m_checkboxCoverageGrid.setParent(layoutV);
		m_checkboxCoverageGrid.setChecked(true);
		m_checkboxCoverageGrid.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_plotFactory.ShowCoverageGrid(m_checkboxCoverageGrid.isChecked());
				m_plotFactory.RequestCoverageRedraw();
				m_plotFactory.DrawPlots();
			}
		});
		
		m_checkboxUseReducedDataSet = new Checkbox("Reduce sample size");
		m_checkboxUseReducedDataSet.setParent(layoutV);
		m_checkboxUseReducedDataSet.setChecked(m_bUseReducedDataSet);
		m_checkboxUseReducedDataSet.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bUseReducedDataSet = m_checkboxUseReducedDataSet.isChecked();				
				
				if(m_bUseReducedDataSet)
				{
					ReduceSampleSize();
				}
				else
				{					
					// get coverage for all samples
//					GetBigWigCoverageForGene();
				}
			}
		});
		
		m_checkboxSkipFirstAndLastExon = new Checkbox("auto. removal: skip first/last exons");
		m_checkboxSkipFirstAndLastExon.setParent(layoutV);
		m_checkboxSkipFirstAndLastExon.setChecked(m_bSkipFirstAndLastExon);
		m_checkboxSkipFirstAndLastExon.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bSkipFirstAndLastExon = m_checkboxSkipFirstAndLastExon.isChecked();
			}
		});

		m_checkboxSkipIntronRetentionDetection = new Checkbox("Detect intron retention events");
		m_checkboxSkipIntronRetentionDetection.setParent(layoutV);
		m_checkboxSkipIntronRetentionDetection.setChecked(m_bSkipFirstAndLastExon);
		m_checkboxSkipIntronRetentionDetection.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bDetectIntronRetentionEvents = m_checkboxSkipIntronRetentionDetection.isChecked();
			}
		});
		
		//####################################################
		//           add coverage plot options
		//####################################################
		// add buttons to option layout
		layoutV = new Vlayout();
		layoutV.setParent(layoutH);
		layoutV.setWidth("200px");
		
		m_checkboxUseMedian = new Checkbox("Show median per condition");
		m_checkboxUseMedian.setWidth("100%");
		m_checkboxUseMedian.setParent(layoutV);
		m_checkboxUseMedian.setChecked(m_plotFactory.IsMedianEnabled());
		m_checkboxUseMedian.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_plotFactory.ShowMedianCoverage(m_checkboxUseMedian.isChecked());
				
				m_checkboxUseMean.setChecked(false);
				m_checkboxUseGeometricMean.setChecked(false);
			}
		});
		
		m_checkboxUseGeometricMean = new Checkbox("Show geometric mean per condition");
		m_checkboxUseGeometricMean.setWidth("100%");
		m_checkboxUseGeometricMean.setParent(layoutV);
		m_checkboxUseGeometricMean.setChecked(m_plotFactory.IsGeometricMeanEnabled());
		m_checkboxUseGeometricMean.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_plotFactory.ShowGeometricMeanCoverage(m_checkboxUseGeometricMean.isChecked());

				m_checkboxUseMean.setChecked(false);
				m_checkboxUseMedian.setChecked(false);
			}
		});
		
		m_checkboxUseMean = new Checkbox("Show mean per condition");
		m_checkboxUseMean.setWidth("100%");
		m_checkboxUseMean.setParent(layoutV);
		m_checkboxUseMean.setChecked(m_plotFactory.IsMeanEnabled());
		m_checkboxUseMean.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_plotFactory.ShowMeanCoverage(m_checkboxUseMean.isChecked());

				m_checkboxUseMedian.setChecked(false);
				m_checkboxUseGeometricMean.setChecked(false);
			}
		});
		
		m_checkboxShowQuartiles = new Checkbox("Show quartiles");
		m_checkboxShowQuartiles.setWidth("100%");
		m_checkboxShowQuartiles.setParent(layoutV);
		m_checkboxShowQuartiles.setChecked(m_plotFactory.IsQuartilesEnabled());
		m_checkboxShowQuartiles.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_plotFactory.ShowQuartiles(m_checkboxShowQuartiles.isChecked());
			}
		});
		
		//####################################################
		//            add valid isoform options
		//####################################################
		layoutH = new Hlayout();
		layoutH.setStyle("margin-top: 20px");
		layoutH.setParent(grpBox);
		
		layoutV = new Vlayout();
		layoutV.setWidth("220px");
		layoutV.setStyle("margin-left: 5px;");
		layoutV.setParent(layoutH);
		
		Label lab = new Label("Minimum junction reads:");
		lab.setParent(layoutV);
		
		m_textboxMinJunctionReads = new Textbox("" + m_nMinJunctionReads);
		m_textboxMinJunctionReads.setParent(layoutV);
		m_textboxMinJunctionReads.addEventListener(Events.ON_CHANGE, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Textbox box = (Textbox)event.getTarget();
				try
				{
					m_nMinJunctionReads = Integer.parseInt(box.getText());
					m_dataSupplier.IdentifyInvalidJunctions(m_hWindow, false);
				}
				catch(NumberFormatException e)
				{
					Messagebox.show("The minimum junction read must be a non-negative number");
				}
			}
		});
		
		lab = new Label("Minimum coverage per base:");
		lab.setParent(layoutV);
		
		m_textboxMinCovPerBase = new Textbox("" + m_nMinCovPerBase);
		m_textboxMinCovPerBase.setParent(layoutV);
		m_textboxMinCovPerBase.addEventListener(Events.ON_CHANGE, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Textbox box = (Textbox)event.getTarget();
				try
				{
					m_nMinCovPerBase = Integer.parseInt(box.getText());
				}
				catch(NumberFormatException e)
				{
					Messagebox.show("The minimum coverage per base must be a non-negative number");
				}
			}
		});
		
		layoutV = new Vlayout();
		layoutV.setWidth("220px");
		layoutV.setParent(layoutH);
		
		lab = new Label("Minimum fraction of covered bases:");
		lab.setParent(layoutV);
		
		m_textboxMinCoveredBases = new Textbox(String.format(Locale.ENGLISH, "%.2f", m_fMinCoveredBases));
		m_textboxMinCoveredBases.setParent(layoutV);
		m_textboxMinCoveredBases.addEventListener(Events.ON_CHANGE, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Textbox box = (Textbox)event.getTarget();
				try
				{
					m_fMinCoveredBases = Double.parseDouble(box.getText());
					if(m_fMinCoveredBases < 0.0 || m_fMinCoveredBases > 1.0)
					{
						Messagebox.show("The minimum fraction of covered bases must be a number between 0 and 1");
						m_fMinCoveredBases = 0.7;
					}
					
					System.out.println("new covered bases: " + m_fMinCoveredBases);
				}
				catch(NumberFormatException e)
				{
					Messagebox.show("The minimum fraction of covered bases must be a number between 0 and 1");
				}
			}
		});
		
		lab = new Label("Alternative splicing threshold:");
		lab.setParent(layoutV);
		
		m_textboxVariableExonThreshold = new Textbox(String.format(Locale.ENGLISH, "%.2f", m_fVariableExonThreshold));
		m_textboxVariableExonThreshold.setParent(layoutV);
		m_textboxVariableExonThreshold.addEventListener(Events.ON_CHANGE, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Textbox box = (Textbox)event.getTarget();
				try
				{
					m_fVariableExonThreshold = Double.parseDouble(box.getText());
					if(m_fVariableExonThreshold <= 0.0)
					{
						Messagebox.show("The minimum difference to report alternative spliced exons must be a number > 0");
						m_fVariableExonThreshold = 0.4;
					}
				}
				catch(NumberFormatException e)
				{
					Messagebox.show("The minimum fraction of covered bases must be a number between 0 and 1");
				}
			}
		});

		//#############################################################	
		layoutH = new Hlayout();
		layoutH.setStyle("margin-top: 20px;");
		layoutH.setParent(grpBox);
		
		layoutV = new Vlayout();
		layoutV.setStyle("margin-left: 5px;");
		layoutV.setParent(layoutH);
		
		lab = new Label("Selected entropy index:");
		lab.setParent(layoutV);
		
		// add combo box for different entropy/index measurements
		m_comboboxSelectedEntropyIndex = new Combobox("Gini");
		m_comboboxSelectedEntropyIndex.setWidth("160px");
		m_comboboxSelectedEntropyIndex.setParent(layoutV);
		for(int i=0; i<5; i++)
		{
			String strText = "";
			switch(i)
			{
				case 0: strText = "None"; break;
				case 1: strText = "Gini"; break;
				case 2: strText = "Simpson"; break;
				case 3: strText = "Theil"; break;
				case 4: strText = "Shannon"; break;
			}
			Comboitem item = new Comboitem(strText);
			item.setValue(strText);
			
			m_comboboxSelectedEntropyIndex.appendChild(item);
		}
		m_comboboxSelectedEntropyIndex.addEventListener(Events.ON_SELECT, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Combobox box = (Combobox)event.getTarget();
				if(box.getSelectedItem() != null)
				{
					if(m_mapCountsPerExonAndTissue != null)
						m_mapCountsPerExonAndTissue.clear();
					else
						m_mapCountsPerExonAndTissue = new TreeMap<String, TreeMap<String, Vector<Double>>>();
					
					if(m_mapEntropyToExonicPart != null)
						m_mapEntropyToExonicPart.clear();
					else
						m_mapEntropyToExonicPart = new TreeMap<String, Double>();
					
					m_mapCountsPerExonAndTissue = GetGTEXDataForExonicParts();

					String strSelection = (String)box.getSelectedItem().getValue();
					
					EntropyFactory entropyFactory = new EntropyFactory();
					
					switch(strSelection)
					{
						case "None":
							m_bShowEntropyData = false;							
							break;
							
						case "Gini":
							m_mapEntropyToExonicPart = entropyFactory.CalculateGiniIndexForExons(m_mapCountsPerExonAndTissue);
							m_bShowEntropyData = true;
							break;
							
						case "Simpson":
							m_mapEntropyToExonicPart = entropyFactory.CalculateGiniSimpsonIndexForExons(m_mapCountsPerExonAndTissue);
							m_bShowEntropyData = true;
							break;
							
						case "Shannon":
							m_mapEntropyToExonicPart = entropyFactory.CalculateShannonIndexForExons(m_mapCountsPerExonAndTissue);
							m_bShowEntropyData = true;
							break;
							
						case "Theil":
							m_mapEntropyToExonicPart = entropyFactory.CalculateTheilIndexForExons(m_mapCountsPerExonAndTissue);
							m_bShowEntropyData = true;
							break;
					}
					
					// force isoform redraw
					m_plotFactory.DrawIsoforms();
				}
			}
		});
	
		Button btnGetURL = new Button("Generate HTML link for current view");
		btnGetURL.setParent(layoutV);
		btnGetURL.setClass("button blue");
		btnGetURL.setStyle("margin-top: 80px; margin-bottom: 10px;");
		
		m_txtboxURL = new Textbox();
		m_txtboxURL.setParent(layoutV);
		m_txtboxURL.setStyle("width: 420px;");
		
		btnGetURL.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				// get current path
				String strServer 		= Executions.getCurrent().getServerName();
				int nServerPort 		= Executions.getCurrent().getServerPort();
				String strContentsPath 	= Executions.getCurrent().getContextPath();
				String strRequestPath	= Executions.getCurrent().getDesktop().getRequestPath();
				
				// get project ID
				if(m_projectModel == null || m_projectModel.GetProjectName() == null || m_projectModel.GetProjectName().length() == 0)
				{
					Messagebox.show("You need to select a project first.");
					return;
				}

				int nProjectID = m_projectModel.GetProjectName().hashCode();
				
				// get selected gene
				if(m_dataSupplier.GetGene() == null)
				{
					Messagebox.show("You need to select a gene first.");
					return;
				}
				
				String strSelectedGene = m_dataSupplier.GetGeneID();
				
				// get selected isoforms
				Set<Treeitem> vcSelectedItems = m_treeSelectedIsoforms.getSelectedItems();
				long nSelectedIsoforms = 0;
				int i=0;
				for(Treeitem item : m_treeSelectedIsoforms.getItems())
				{
					if(vcSelectedItems.contains(item))
						nSelectedIsoforms += (long)Math.pow(2, i);
					
					i++;
				}
				
				// get selected samples
				Collection<Treeitem> vcAllSamples = m_treeSelectedSamples.getItems();
				Set<Treeitem> vcSelectedSamples = m_treeSelectedSamples.getSelectedItems();
				String strSelectedSamples	= "";
				long nSelectedSamples		= 0;
				i=0;
				
				for(Treeitem item : vcAllSamples)
				{
					if(vcSelectedSamples.contains(item))
					{
						nSelectedSamples += (long)Math.pow(2, i);
					}
					
					i++;
					
					if(i==63)
					{
						strSelectedSamples += nSelectedSamples+",";
						nSelectedSamples = 0;
						i=0;
					}
				}
				
				strSelectedSamples += nSelectedSamples;
				
				String strURL = "?";
				if(nServerPort != 443)
					strURL = strServer + ":" + nServerPort + strContentsPath + strRequestPath + "?ID=" + nProjectID + "&reference=" + m_strFileGTF.hashCode() + "&gene=" + strSelectedGene + "&isoforms=" + nSelectedIsoforms + "&samples=" + strSelectedSamples;
				else
					strURL = strServer + strContentsPath + strRequestPath + "?ID=" + nProjectID + "&gene=" + strSelectedGene + "&isoforms=" + nSelectedIsoforms + "&samples=" + strSelectedSamples;
//				Messagebox.show(strURL);
				m_txtboxURL.setText(strURL);
				
				/*
				// copy to clipboard -> ONLY SERVER SIDED
				Clipboard clpbrd = Toolkit.getDefaultToolkit().getSystemClipboard();
				clpbrd.setContents(new StringSelection(strURL), null);
				*/
			}
		});
	}

	/** Prepares the color selection menu */
	public void AddColorSelectionMenu(Vlayout parentLayout)
	{
		m_ColorSelectionGroupBox = new Groupbox();
		m_ColorSelectionGroupBox.setTitle("Colors");
		m_ColorSelectionGroupBox.setParent(parentLayout);		
		m_ColorSelectionGroupBox.setWidth("200px");
		m_ColorSelectionGroupBox.setHeight("490px");
		m_ColorSelectionGroupBox.setContentStyle("overflow:auto;");
		m_ColorSelectionGroupBox.setMold("3d");
		
		Button btn_SaveColorScheme = new Button("Save Color Scheme");
		btn_SaveColorScheme.setParent(parentLayout);
		btn_SaveColorScheme.setClass("button green");
		btn_SaveColorScheme.setStyle("display:inline-block; margin-left: 40px; margin-bottom: 10px;");
		
		btn_SaveColorScheme.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				if(m_projectModel == null || !m_projectModel.IsReady())
					return;
				
				String strOut = m_projectModel.GetFullPathOfProjectFile() + ".color_scheme";
				PrintWriter pOut = new PrintWriter(new File(strOut));
				
				TreeMap<String, TreeSet<String>> mapConditionsPerConditionType = m_projectModel.GetConditionsToConditionTypes();
				for(String strConditionType : mapConditionsPerConditionType.keySet())
				{
					for(String strCondition : mapConditionsPerConditionType.get(strConditionType))
					{
						Color clr = m_plotFactory.GetColorForCondition(strCondition);
						pOut.println(strCondition + "\t" + clr.getRed() + "\t" + clr.getGreen() + "\t" + clr.getBlue());
					}
				}
				
				pOut.close();
				
				Messagebox.show("Saved color scheme successfully");
			}
		});
	}
	
	/** Updates the color selection menu*/
	public void UpdateColorSelection()
	{
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = m_projectModel.GetSelectedSamplesPerCondition(m_strSelectedConditionType, m_vcSelectedSamples);
		
		if(mapSamplesToConditions.keySet().size() == 0)
			return;
		
		int nHeight = mapSamplesToConditions.keySet().size()*30;
		BufferedImage img = new BufferedImage(200, nHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D graph = img.createGraphics();
		
		m_ColorSelectionGroupBox.getChildren().clear();

		m_imgMapColors = new Imagemap();
		m_imgMapColors.setWidth(200+"px");
		m_imgMapColors.setHeight(nHeight+"px");
		m_imgMapColors.setParent(m_ColorSelectionGroupBox);

		// define selected conditions
		TreeSet<String> vcSelectedConditions = new TreeSet<String>();
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			if(mapSamplesToConditions.get(strCondition).size() > 0)
				vcSelectedConditions.add(strCondition);
		}
		
		// fill background
		graph.setColor(Color.WHITE);
		graph.fillRect(0, 0, 200, nHeight);
		
		int nOffset = 0;
		for(String strCondition : vcSelectedConditions)
		{
			Color clr = m_plotFactory.GetColorForCondition(strCondition);
			graph.setColor(clr);
			graph.fillRect(3, nOffset, 20, 20);
			graph.drawString(strCondition, 28, nOffset+15);
			
			graph.setColor(Color.BLACK);
			graph.drawRect(3, nOffset, 20, 20);
			
			// add clickable area
			Area area = new Area();
			area.setId("clr�" + strCondition);
			area.setShape("rectangle");
			area.setCoords(3 + ", " + nOffset + ", " + (3+20) + ", " + (nOffset+20));				
			area.setTooltiptext("Color for: " + strCondition);
			area.setParent(m_imgMapColors);

			nOffset+=30;
		}

		m_imgMapColors.setContent(img);
		
		m_imgMapColors.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				MouseEvent evnt = (MouseEvent) event;

				String strCondition = evnt.getArea().split("�")[1];

				m_windowPopup.setParent(m_hWindow);
				m_windowPopup.doPopup();
				m_windowPopup.setTitle("Color Selection");
				m_windowPopup.setStyle("relative;left:500px;top:100px;");
				m_windowPopup.getChildren().clear();
				
				int nWidth = 240;
				int nHeight = 260;
				
				m_windowPopup.setWidth(nWidth + "px");
				m_windowPopup.setHeight(nHeight +"px");
				
				Imagemap imgMap = new Imagemap();
				imgMap.setWidth(nWidth + "px");
				imgMap.setHeight(nHeight + "px");
				imgMap.setParent(m_windowPopup);
	
				UpdateColorSelectionPopup(imgMap, nWidth, nHeight, m_plotFactory.GetColorForCondition(strCondition), strCondition);
			}
		});
	}
	
	/** Updates the color selection popup window */
	private void UpdateColorSelectionPopup(Imagemap imgMap, int nWidth, int nHeight, Color clr, String strCondition)
	{
		for(Component c : imgMap.getChildren())
		{
			imgMap.removeChild(c);
		}
			
		BufferedImage img = new BufferedImage(nWidth, nHeight, BufferedImage.TYPE_INT_ARGB);
		Graphics2D graph = img.createGraphics();
		
		// fill background
		graph.setColor(Color.WHITE);
		graph.fillRect(0, 0, 800, 500);
		
		// add saturation gradient
		GradientPaint  clrGradient1 = new GradientPaint(0.0f, 0.0f, Color.WHITE, 200.0f, 0.0f, clr);
		GradientPaint  clrGradient2 = new GradientPaint(0.0f, 0.0f, new Color(0,0,0,0), 0.0f, 200.0f, new Color(0,0,0,255));
		graph.setPaint(clrGradient1);
		graph.fillRect(5, 5, 200, 200);
		graph.setPaint(clrGradient2);
		graph.fillRect(5, 5, 200, 200);
		
		// add color gradient
		float pFractions[] = new float[] {0.00f, 0.167f, 0.333f, 0.50f, 0.666f, 0.833f, 1.0f};
		Color pColors[] = new Color[] {new Color(255, 0, 0), new Color(255, 255, 0), new Color(0,255,0), new Color(0, 255, 255), new Color(0,0,255), new Color(255, 0, 255), new Color(255, 0 , 0)};
		
		LinearGradientPaint clrGradient = new LinearGradientPaint(0, 0, 0, 200, pFractions, pColors, MultipleGradientPaint.CycleMethod.NO_CYCLE);
		
		graph.setPaint(clrGradient);
		graph.fillRect(210, 5, 20, 200);
		
		graph.setColor(Color.BLACK);
		graph.drawRect(210, 5, 20, 200);
		
		imgMap.setId(clr.getRed() + "$" + clr.getGreen() + "$" + clr.getBlue() + "$" + strCondition);
		imgMap.setContent(img);
		
		imgMap.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				MouseEvent evnt = (MouseEvent) event;
				
				int nMouseX = evnt.getX(); 
				int nMouseY = evnt.getY();
				
				// was the color hit?
				if(nMouseX >= 210 && nMouseX <= 230 && nMouseY >= 5 && nMouseY <= 205)
				{
					Imagemap imgMap = (Imagemap)evnt.getTarget();
					String pSplit[] = imgMap.getId().split("\\$");
					int r = Integer.parseInt(pSplit[0]);
					int g = Integer.parseInt(pSplit[1]);
					int b = Integer.parseInt(pSplit[2]);
					String strCondition = pSplit[3];
					
					BufferedImage img = new BufferedImage(200, 200, BufferedImage.TYPE_INT_ARGB);
					Graphics2D graph = img.createGraphics();
					
					float pFractions[] = new float[] {0.00f, 0.167f, 0.333f, 0.50f, 0.666f, 0.833f, 1.0f};
					Color pColors[] = new Color[] {new Color(255, 0, 0), new Color(255, 255, 0), new Color(0,255,0), new Color(0, 255, 255), new Color(0,0,255), new Color(255, 0, 255), new Color(255, 0 , 0)};
					
					LinearGradientPaint clrGradient = new LinearGradientPaint(0, 0, 0, 200, pFractions, pColors, MultipleGradientPaint.CycleMethod.NO_CYCLE);
					
					graph.setPaint(clrGradient);
					graph.fillRect(0, 0, 200, 200);
					
					// remove offsets
					nMouseX -= 210;
					nMouseY -= 5;
					
					int rgb = img.getRGB(nMouseX, nMouseY);
					r = (rgb & 0x00ff0000) >> 16;
					g = (rgb & 0x0000ff00) >> 8;
					b = rgb & 0x000000ff;
					
					Color clr = new Color(r, g, b);
					
					evnt.stopPropagation();

					imgMap.setId(r + "$" + g + "$" + b + "$" + strCondition);
					UpdateColorSelectionPopup(imgMap, 240, 260, clr, strCondition);
				}
				// was the saturation hit?
				else if(nMouseX >= 5 && nMouseX < 205 && nMouseY >= 5 && nMouseY < 205)
				{
					Imagemap imgMap = (Imagemap)evnt.getTarget();
					String pSplit[] = imgMap.getId().split("\\$");
					int r = Integer.parseInt(pSplit[0]);
					int g = Integer.parseInt(pSplit[1]);
					int b = Integer.parseInt(pSplit[2]);
					String strCondition = pSplit[3];
					
					Color clr = new Color(r,g,b);

					BufferedImage img = new BufferedImage(200, 200, BufferedImage.TYPE_INT_ARGB);
					Graphics2D graph = img.createGraphics();
					
					// remove offsets
					nMouseX -= 5;
					nMouseY -= 5;
					
					// add saturation gradient
					GradientPaint  clrGradient1 = new GradientPaint(0.0f, 0.0f, Color.WHITE, 200.0f, 0.0f, clr);
					GradientPaint  clrGradient2 = new GradientPaint(0.0f, 0.0f, new Color(0,0,0,0), 0.0f, 200.0f, new Color(0,0,0,255));
					graph.setPaint(clrGradient1);
					graph.fillRect(0, 0, 200, 200);
					graph.setPaint(clrGradient2);
					graph.fillRect(0, 0, 200, 200);
					
					int rgb = img.getRGB(nMouseX, nMouseY);
					r = (rgb & 0x00ff0000) >> 16;
					g = (rgb & 0x0000ff00) >> 8;
					b = rgb & 0x000000ff;
					
					clr = new Color(img.getRGB(nMouseX, nMouseY));

					m_plotFactory.UpdateColorSelection(strCondition, clr);
					UpdateColorSelection();
				}
			}
		});
	}

	/** Processes a stored alternative splicing result and retrieves all data necessary for its visualization */
	public void PrepareHitForVisualization(AnalysisResult res) throws Exception
	{
		String strOldGeneID = "";
		if(m_dataSupplier.GetGene() != null)
			strOldGeneID = m_dataSupplier.GetGeneID();
		
		m_selectedResult 		 = res;
		
		m_nMinCovPerBase 		 = m_selectedResult.GetMinCovPerBase();
		m_fMinCoveredBases 		 = m_selectedResult.GetMinCoveredBases();
		m_nMinJunctionReads 	 = m_selectedResult.GetMinJunctionReads();
		m_fVariableExonThreshold = m_selectedResult.GetVariableExonThreshold();
		
		m_textboxMinCovPerBase.setText(""+m_nMinCovPerBase);
		m_textboxMinCoveredBases.setText(""+m_fMinCoveredBases);
		m_textboxMinJunctionReads.setText(""+m_nMinJunctionReads);
		m_textboxVariableExonThreshold.setText(""+m_fVariableExonThreshold);
		
		GeneIdentifier gid = null;
		if(m_selectedResult.GetGeneID().isEmpty())
		{
			gid = m_geneIdentifierHandler.GetGeneIdentifierForGene(m_selectedResult.GetGeneSymbol(), null);
			
			if(gid == null)
			{
				ErrorMessage.ShowError("Can't find gene ID for gene: " + m_selectedResult.GetGeneSymbol());
				return;
			}
			
			m_selectedResult.SetGeneID(gid.m_strEnsemblGeneID);
		}
		else
			gid = m_geneIdentifierHandler.GetGeneIdentifierForGene(m_selectedResult.GetGeneID().split("\\.")[0], null);
		
		if(m_selectedResult.GetGeneSymbol().isEmpty() || !gid.m_strApprovedGeneSymbol.equals(m_selectedResult.GetGeneSymbol()))
			m_bandboxSelectedGene.setValue(gid.m_strEnsemblGeneID);
		else
			m_bandboxSelectedGene.setValue(m_selectedResult.GetGeneSymbol());
		
		
		// switch gene annotation if the result specified a different one
		if(!m_selectedResult.GetGTFFile().isEmpty() && !m_strFileGTF.equals(m_selectedResult.GetGTFFile()))
		{
			m_strFileGTF = m_selectedResult.GetGTFFile();
			LoadGeneAnnotation(m_strFileGTF);
			
			for(Comboitem item : m_comboboxSelectedGeneAnnotation.getItems())
			{
				File pFile = new File(m_strFileGTF);
				if(item.getLabel().equals(pFile.getName()))
					m_comboboxSelectedGeneAnnotation.setSelectedItem(item);
			}
		}
		
		// check whether the selected gene is a new gene
		boolean bSameGene = m_selectedResult.GetGeneID().equals(strOldGeneID);
		if(!bSameGene)
		{
			OnGeneChange(false, false);
			m_selectedResult = res;
		}
		
		m_vcValidIsoforms = m_selectedResult.GetValidIsoforms(m_dataSupplier);

		if(m_vcValidIsoforms.size() > 0)
		{
			m_vcPreviouslySelectedIsoforms.clear();
			m_treeSelectedIsoforms.clearSelection();
			
			for(Treeitem item : m_treeSelectedIsoforms.getItems())
			{
				String strID = item.getLabel();
				
				for(String strIsoform : m_vcValidIsoforms)
				{
					m_vcPreviouslySelectedIsoforms.add(strIsoform);
					
					if(strID.equals(strIsoform.split("\\.")[0]))
					{
						m_treeSelectedIsoforms.addItemToSelection(item);
					}
				}
			}			
		}
		else // select find isoforms matching to the variable exon in the result
		{
			boolean bIsoformSelectionFinished = false;
			Vector<String> vcValidIsoforms = new Vector<String>();
			
			// step 1: Try to find matching isoforms based on junction reads
			if(res.HasPSIScore())
			{
				boolean bFoundIsoformForInclusion = false;
				boolean bFoundIsoformForExclusion = false;
				
				CountElement junIncl = res.GetInclusionJunction();
				CountElement junExcl = res.GetExclusionJunction();
				
				for(String strIsoform : m_dataSupplier.GetIsoformNames())
				{
					TreeSet<CountElement> vcJunctions = m_dataSupplier.GetJunctionsForIsoform(strIsoform);
					if(vcJunctions == null)
						continue;
					
					if(vcJunctions.contains(junIncl))
					{
						bFoundIsoformForInclusion = true;
						vcValidIsoforms.add(strIsoform);
						break;
					}
					
					if(vcJunctions.contains(junExcl))
					{
						bFoundIsoformForExclusion = true;
						vcValidIsoforms.add(strIsoform);
						break;
					}
				}
				
				// add isoform for unknown junctions
				if(bFoundIsoformForExclusion && !bFoundIsoformForInclusion)
				{
					m_plotFactory.AddUnknownJunction(junIncl);
					bIsoformSelectionFinished = true;
				}
				
				if(bFoundIsoformForInclusion && !bFoundIsoformForExclusion)
				{
					m_plotFactory.AddUnknownJunction(junExcl);
					bIsoformSelectionFinished = false;
				}
			}
			
			if(!bIsoformSelectionFinished && res.HasAltExonA())
			{
				Vector<String> vcFullMatchInclusion = new Vector<String>();
				Vector<String> vcOverlappingInclusion = new Vector<String>();
				Vector<String> vcSkipping = new Vector<String>();
				
				int nOldDist = 0; // defines how similar the closes matching exon is
				Exon exonClosestMatch = null;
				
				// step 2: Find matching isoforms based on the overlap to a meta exon
				for(String strIsoform : m_dataSupplier.GetIsoformNames())
				{
					Exon pExons[] = m_dataSupplier.GetExonsForIsoform(strIsoform);					
					
					// find overlap to isoform exons
					for(Exon ex : pExons)
					{
						if(res.GetStartA() == ex.getCodingStart() && res.GetEndA() == ex.getCodingStop())
						{
							vcFullMatchInclusion.add(strIsoform);
							break;
						}
						else if(res.GetStartA() >= ex.getCodingStart() && res.GetEndA() <= ex.getCodingStop())
						{
							if(exonClosestMatch == null)
							{
								exonClosestMatch = ex;
								nOldDist = Math.abs(ex.getCodingStart() - res.GetStartA()) + Math.abs(ex.getCodingStop() - res.GetEndA());
								vcOverlappingInclusion.add(strIsoform);
							}
							else
							{
								int nNewDist = Math.abs(ex.getCodingStart() - res.GetStartA()) + Math.abs(ex.getCodingStop() - res.GetEndA());
								
								if(nNewDist < nOldDist)
								{
									exonClosestMatch = ex;
									nOldDist = nNewDist;
									vcOverlappingInclusion.clear();
									vcOverlappingInclusion.add(strIsoform);
								}
								else if(nNewDist == nOldDist)
								{
									vcOverlappingInclusion.add(strIsoform);
								}
							}
							
							break;
						}
					}

					// find junctions that support skipping of the exon
					TreeSet<CountElement> vcJunctions = m_dataSupplier.GetJunctionsForIsoform(strIsoform);
					for(CountElement jun : vcJunctions)
					{
						if(vcJunctions == null)
							continue;
						
						if(res.GetStartA() > jun.m_nStart && res.GetEndA() < jun.m_nEnd)
						{
							vcSkipping.add(strIsoform);
							break;
						}
					}
				}
				
				if(vcFullMatchInclusion.size() > 0)
				{
					vcValidIsoforms.addAll(vcFullMatchInclusion);
					bIsoformSelectionFinished = true;
				}
				else if(vcOverlappingInclusion.size() > 0)
				{
					vcValidIsoforms.addAll(vcOverlappingInclusion);
					bIsoformSelectionFinished = true;
					
					res.GetAltExonA().m_nStart = exonClosestMatch.getCodingStart();
					res.GetAltExonA().m_nEnd   = exonClosestMatch.getCodingStop();
				}
				
				if(bIsoformSelectionFinished)
					vcValidIsoforms.addAll(vcSkipping);
			}
			
			if(!bIsoformSelectionFinished)
			{
				for(String strIsoform : m_dataSupplier.GetIsoformNames())
					vcValidIsoforms.add(strIsoform);
			}
			
			m_vcValidIsoforms.clear();
			m_vcPreviouslySelectedIsoforms.clear();
			m_treeSelectedIsoforms.clearSelection();
		
			for(String strIsoform : vcValidIsoforms)
			{
				m_vcValidIsoforms.add(strIsoform);
				m_vcPreviouslySelectedIsoforms.add(strIsoform);
			}

			for(Treeitem item : m_treeSelectedIsoforms.getItems())
			{
				if(m_vcValidIsoforms.contains(item.getValue()))
					m_treeSelectedIsoforms.addItemToSelection(item);
			}
		}
		OnIsoformChange(!bSameGene);
		
		// check if any of the junctions is not annotated (not available for imported results)
		if(res.HasPSIScore() && res.GetResultSource() == -1)
		{
			if(m_dataSupplier.IsNovelJunction(res.GetInclusionJunction()))
			{
				m_plotFactory.AddUnknownJunction(res.GetInclusionJunction());
			}
			
			if(m_dataSupplier.IsNovelJunction(res.GetExclusionJunction()))
			{
				m_plotFactory.AddUnknownJunction(res.GetExclusionJunction());
			}
		}
	}
	
	/** Prepares the sample selection menu */
	public void AddComboboxForSamples(Vlayout parentLayout)
	{		
		if(m_treeSelectedSamples == null && parentLayout != null)
		{
			//####################################################
			//    add tree for sample selection to options box
			//####################################################
			m_treeSelectedSamples = new Tree();
			m_treeSelectedSamples.setId("sample_tree");
			m_treeSelectedSamples.setParent(parentLayout);
			m_treeSelectedSamples.setSizedByContent(true);
			m_treeSelectedSamples.setHeight("490px");
			m_treeSelectedSamples.setWidth("280px");
			m_treeSelectedSamples.setMultiple(true);
			m_treeSelectedSamples.setCheckmark(true);
			
			m_treeSelectedSamples.addEventListener(Events.ON_SELECT, new EventListener<Event>()
			{
				@Override
				public void onEvent(Event event) throws Exception
				{
					if(m_bUseReducedDataSet)
					{
						Messagebox.show("You are about to change the sample selection although the \"reduced sample size\" option is selected.\n This will deselect the option and will cause the application to process the bigwig files of all samples. This might take quite some time for large data sets.",
								null,
								Messagebox.OK | Messagebox.CANCEL,
								Messagebox.EXCLAMATION,
								new EventListener<Event>()
								{
									@Override
									public void onEvent(Event event)
									{
										switch (((Integer) event.getData()).intValue())
										{
											case Messagebox.OK:
											{
												UpdateComboboxForSampleSelection();
												m_bUseReducedDataSet = false;
												m_checkboxUseReducedDataSet.setChecked(false);
												break;
											}
											
											case Messagebox.CANCEL:
											{
												// re-select the item that was clicked
												for(Treeitem item : m_treeSelectedSamples.getItems())
												{
													if(m_vcPreviouslySelectedSamples.contains(item.getLabel()))
														m_treeSelectedSamples.addItemToSelection(item);
												}
												break;
											}
										}
									}
								});
					}
					else
					{
						UpdateComboboxForSampleSelection();
					}
				}
			});
		}
		else
		{
			// clear old data
			m_treeSelectedSamples.getChildren().clear();
		}

		// add header
		Treecols cols = new Treecols();
		Treecol col = new Treecol("Samples"); col.setParent(cols);
		cols.setParent(m_treeSelectedSamples);
		
		//############################################
		//         populate tree with samples
		//############################################
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		Treechildren children = new Treechildren();
		children.setParent(m_treeSelectedSamples);
		
		m_treeSelectedSamples.invalidate();
		m_vcPreviouslySelectedSamples.clear();
		m_vcSelectedSamples.clear();

		for(String strCondition : mapSamplesToConditions.keySet())
		{
			Treeitem conditionItem = new Treeitem();
			conditionItem.setParent(children);
			conditionItem.setLabel(strCondition);
			
			Treechildren conditionChildren = new Treechildren();
			conditionChildren.setParent(conditionItem);
			
			// add condition to previous selection
			m_vcPreviouslySelectedSamples.add(strCondition);
			
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				Treeitem sampleItem = new Treeitem();
				sampleItem.setParent(conditionChildren);

				m_treeSelectedSamples.addItemToSelection(sampleItem);
				m_vcPreviouslySelectedSamples.add(strSample);
				m_vcSelectedSamples.add(strSample);
	
				Treerow row = new Treerow();
				row.setParent(sampleItem);
				
				Treecell cell = new Treecell(strSample); cell.setParent(row);
			}
			
			m_vcPreviouslySelectedSamples.add(strCondition);
		}
		
		m_treeSelectedSamples.setSizedByContent(true);
		
		Hlayout layout = new Hlayout();
		layout.setHeight("40px");
		layout.setParent(parentLayout);
			
		// add button to unselect all samples
		Button btnUnselectAll = new Button("Unselect all samples");
		btnUnselectAll.setSclass("button green");
		btnUnselectAll.setParent(layout);		

		btnUnselectAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				m_treeSelectedSamples.clearSelection();				
				m_vcPreviouslySelectedSamples.clear();
				m_vcSelectedSamples.clear();
				
				// uncheck the "use reduced data set" checkbox
				m_checkboxUseReducedDataSet.setChecked(false);
				m_bUseReducedDataSet = false;
				
				UpdateComboboxForSampleSelection();
			}
		});
		
		// use a reduced data set if there are more than 10 samples per condition. Scales up to 50 samples in total.
		if(m_bUseReducedDataSet)
		{
			ReduceSampleSize();
		}
		else
		{
			UpdateComboboxForSampleSelection();
		}
		
		// add button to open a popup to select a list of samples
		Button btnSelectListOfSamples = new Button("Select by sample list");
		btnSelectListOfSamples.setSclass("button green");
		btnSelectListOfSamples.setParent(layout);

		btnSelectListOfSamples.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				//############################
				//    prepare popup window
				//############################
				Window popup = new Window();
				popup.setTitle("Sample selection");
				popup.setPosition("center");
				popup.setSizable(true);
				popup.setClosable(true);
				popup.setMinimizable(true);
				popup.setMaximizable(true);
				popup.setBorder(true);
				popup.doPopup();
				popup.setTopmost();
				popup.setParent(m_hWindow);				
				
				Vlayout layout = new Vlayout();
				layout.setParent(popup);
				
				m_textboxSampleSelectionList = new Textbox("Paste sample names here");
				m_textboxSampleSelectionList.setMultiline(true);
				m_textboxSampleSelectionList.setWidth("200px");
				m_textboxSampleSelectionList.setHeight("400px");
				m_textboxSampleSelectionList.setParent(layout);
				
				Button btnSelect = new Button("Select");
				btnSelect.setParent(layout);
				
				btnSelect.addEventListener(Events.ON_CLICK, new EventListener<Event>()
				{
					@Override
					public void onEvent(Event event) throws Exception
					{
						m_treeSelectedSamples.clearSelection();				
						m_vcPreviouslySelectedSamples.clear();
						m_vcSelectedSamples.clear();
						
						TreeSet<String> vcAllSamples = m_projectModel.GetSamples();
						
						String pSamples[] = m_textboxSampleSelectionList.getText().split("[\\r?\\n?\u0000]+");
						for(String strSample : pSamples)
						{
							strSample = strSample.trim();

							if(!vcAllSamples.contains(strSample))
								continue;
	
							m_vcSelectedSamples.add(strSample);

							for(Treeitem item : m_treeSelectedSamples.getItems())
							{
								if(item.getLabel().equals(strSample))
								{
									m_treeSelectedSamples.addItemToSelection(item);
								}
							}
						}
						UpdateComboboxForSampleSelection();
					}
				});
			}
		});
		
		// use a reduced data set if there are more than 10 samples per condition. Scales up to 50 samples in total.
		if(m_bUseReducedDataSet)
		{
			ReduceSampleSize();
		}
		else
		{
			UpdateComboboxForSampleSelection();
		}
	}

	/** Updates the sample selection menu */
	private void UpdateComboboxForSampleSelection()
	{
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		// sample selection changed, this is important to retrieve new coverage data on "redraw" or "reanalysis"
		m_bSampleSelectionChanged = true;

		// check which items were newly selected/unselected
		for(Treeitem item : m_treeSelectedSamples.getSelectedItems())
		{
			if(!m_vcPreviouslySelectedSamples.contains(item.getLabel()))
			{
				// check whether the item is a condition item
				boolean bIsConditionItem = false;
				
				for(String strCondition : mapSamplesToConditions.keySet())
				{
					if(strCondition.equals(item.getLabel()))
					{						
						bIsConditionItem = true;
						// if so, select all samples of this condition type
						for(String strSample : mapSamplesToConditions.get(strCondition))
						{
							m_vcSelectedSamples.add(strSample);
						}
						m_vcSelectedSamples.add(item.getLabel());
						break;
					}
				}

				if(!bIsConditionItem)
				{
					m_vcSelectedSamples.add(item.getLabel());
				}
			}
		}

		// get newly unselected items
		for(String strItem : m_vcPreviouslySelectedSamples)
		{
			boolean bIncluded = false;
			for(Treeitem item : m_treeSelectedSamples.getSelectedItems())
			{
				if(item.getLabel().equals(strItem))
					bIncluded = true;
			}

			// it was previously selected but now it's not?
			if(!bIncluded)
			{				
				for(Treeitem item : m_treeSelectedSamples.getItems())
				{
					if(item.getLabel().equals(strItem))
					{
						// check whether the item is a condition item
						boolean bIsConditionItem = false;
						
						for(String strCondition : mapSamplesToConditions.keySet())
						{
							if(strCondition.equals(strItem))
							{								
								bIsConditionItem = true;
								// if so, select all samples of this condition type
								for(String strSample : mapSamplesToConditions.get(strCondition))
								{
									m_vcSelectedSamples.remove(strSample);
								}
								m_vcSelectedSamples.remove(item.getLabel());
								break;
							}
						}
						
						if(!bIsConditionItem)
						{
							m_vcSelectedSamples.remove(strItem);
						}
					}
				}
			}
		}
		
		// update selected categories
		TreeSet<String> vcFullySelectedCategories 		= new TreeSet<String>();
		TreeSet<String> vcPartiallySelectedCategories 	= new TreeSet<String>();
		
		// also track which conditions have any selected samples
		TreeSet<String> vcConditionsWithSelectedSamples = new TreeSet<String>();
		for(Treeitem item : m_treeSelectedSamples.getItems())
		{
			for(String strCondition : mapSamplesToConditions.keySet())
			{
				if(strCondition.equals(item.getLabel()))
				{					
					int nValidSamples = 0;
					for(String strSample : mapSamplesToConditions.get(strCondition))
					{
						if(m_vcSelectedSamples.contains(strSample))
						{
							nValidSamples++;
						}
					}
					
					if(nValidSamples > 0)
					{
						vcConditionsWithSelectedSamples.add(strCondition);
					
						if(nValidSamples == mapSamplesToConditions.get(strCondition).size())
						{
							vcFullySelectedCategories.add(strCondition);
						}
						else
						{
							vcPartiallySelectedCategories.add(strCondition);
						}
					}
					break;
				}
			}
		}
		
		// check all selected samples
		m_vcPreviouslySelectedSamples.clear();
		m_treeSelectedSamples.clearSelection();
		for(Treeitem item : m_treeSelectedSamples.getItems())
		{
			if(vcFullySelectedCategories.contains(item.getLabel()))
			{
				m_treeSelectedSamples.addItemToSelection(item);
				item.setOpen(false);
				m_vcPreviouslySelectedSamples.add(item.getLabel());
			}
			else if(vcPartiallySelectedCategories.contains(item.getLabel()))
			{
				m_treeSelectedSamples.removeItemFromSelection(item);
				item.setOpen(true);
			}
			else if(mapSamplesToConditions.keySet().contains(item.getLabel()))
			{
				m_treeSelectedSamples.removeItemFromSelection(item);
				item.setOpen(false);
			}
			else if(m_vcSelectedSamples.contains(item.getLabel()))
			{
				m_treeSelectedSamples.addItemToSelection(item);
				m_vcPreviouslySelectedSamples.add(item.getLabel());
			}
		}

		if(vcConditionsWithSelectedSamples.size() > 0)
		{
			if(m_strSelectedCondition == null || !vcConditionsWithSelectedSamples.contains(m_strSelectedCondition))
			{
				m_strSelectedCondition = vcConditionsWithSelectedSamples.first();
				for(Comboitem item : m_comboboxSelectedCondition.getItems())
				{
					if(item.getLabel().equals(m_strSelectedCondition))
					{
						m_comboboxSelectedCondition.setSelectedItem(item);
						break;
					}
				}
			}
		}

		UpdateColorSelection();
	}
	
	/** Prepares the condition selection combo box */
	public void AddComboboxesForCondition(Hlayout layoutH)
	{
		Vlayout layoutV = new Vlayout();
		// add label for condition type selection
		layoutV.setStyle("margin-top: 10px");
		layoutV.setParent(layoutH);
		
		// add label for condition type selection
		Label lab = new Label("Selected condition type:");		
		lab.setParent(layoutV);
		
		if(m_comboboxSelectedCondition == null)
		{
			m_comboboxSelectedConditionType = new Combobox("Select Condition Type");
			m_comboboxSelectedConditionType.setWidth("200px");
			
			m_comboboxSelectedConditionType.setParent(layoutV);
			m_comboboxSelectedConditionType.addEventListener(Events.ON_SELECT, new EventListener<Event>()
			{
				public void onEvent(Event event) throws Exception
				{
					Combobox box = (Combobox)event.getTarget();
					if(box.getSelectedItem() != null)
					{
						// get condition type
						m_strSelectedConditionType = box.getSelectedItem().getValue();
						OnConditionTypeChange();
					}
				}
			});
		}
		
		layoutV = new Vlayout();
		layoutV.setStyle("margin-top: 10px; margin-left: 40px");
		layoutV.setParent(layoutH);
		
		// add label for condition selection
		lab = new Label("Selected condition:");
		lab.setParent(layoutV);
		
		if(m_comboboxSelectedCondition == null)
		{
			m_comboboxSelectedCondition = new Combobox("Select Condition");
			m_comboboxSelectedCondition.setWidth("200px");
			m_comboboxSelectedCondition.setParent(layoutV);
			m_comboboxSelectedCondition.addEventListener(Events.ON_SELECT, new EventListener<Event>()
			{
				public void onEvent(Event event) throws Exception
				{
					Combobox box = (Combobox)event.getTarget();
					if(box.getSelectedItem() != null)
					{
						m_strSelectedCondition = box.getSelectedItem().getValue();
						
						if(m_dataSupplier.GetGene() != null)
							m_plotFactory.DrawPlots();
					}
				}
			});
		}
	}

	/** Updates the condition selection combo box */
	public void UpdateComboboxesForCondition()
	{
		m_comboboxSelectedCondition.getChildren().clear();
		m_comboboxSelectedConditionType.getChildren().clear();
		
		TreeMap<String, TreeSet<String>> mapConditionsToConditonTypes = m_projectModel.GetConditionsToConditionTypes();
		
		for(String strConditionType : mapConditionsToConditonTypes.keySet())
		{
			Comboitem item = new Comboitem(strConditionType);
			item.setValue(strConditionType);
			
			m_comboboxSelectedConditionType.appendChild(item);
		}
	}
	
	/** Prepares the isoform selection menu */
	public void AddComboboxForIsoformSelection(Vlayout parentLayout) throws IOException
	{
		if(m_treeSelectedIsoforms == null)
		{
			m_treeSelectedIsoforms = new Tree();
			m_treeSelectedIsoforms.setId("isoform_tree");
			m_treeSelectedIsoforms.setParent(parentLayout);
			m_treeSelectedIsoforms.setWidth("400px");
			m_treeSelectedIsoforms.setHeight("490px");
//			m_treeSelectedIsoforms.setAutopaging(true);
			m_treeSelectedIsoforms.setMultiple(true);
			m_treeSelectedIsoforms.setCheckmark(true);
			
			Treecols cols = new Treecols();
			Treecol col = new Treecol("Ensembl ID"); col.setParent(cols);
			col = new Treecol("RefSeq ID"); col.setParent(cols);
			col = new Treecol("Length [bp]"); col.setParent(cols);
			col = new Treecol("Exons [N]"); col.setParent(cols);
			cols.setParent(m_treeSelectedIsoforms);
			
			m_treeSelectedIsoforms.addEventListener(Events.ON_SELECT, new EventListener<Event>()
			{
				@Override
				public void onEvent(Event event) throws Exception
				{					
					// get newly selected items
					for(Treeitem item : m_treeSelectedIsoforms.getSelectedItems())
					{
						if(!m_vcPreviouslySelectedIsoforms.contains(item.getValue().toString()))
						{
							m_treeSelectedIsoforms.addItemToSelection(item);
							m_vcValidIsoforms.add(item.getValue().toString());
						}
					}
					
					// get newly unselected items
					for(String strItem : m_vcPreviouslySelectedIsoforms)
					{
						boolean bIncluded = false;
						for(Treeitem item : m_treeSelectedIsoforms.getSelectedItems())
						{
							if(item.getValue().toString().equals(strItem))
								bIncluded = true;
						}
						
						if(!bIncluded)
						{
							for(Treeitem item : m_treeSelectedIsoforms.getItems())
							{
								if(item.getValue().toString().equals(strItem))
								{
									m_treeSelectedIsoforms.removeItemFromSelection(item);
									m_vcValidIsoforms.remove(strItem);
									break;
								}
							}
						}
					}
					
					m_vcPreviouslySelectedIsoforms.clear();
					
					for(String strIsoform : m_vcValidIsoforms)
						m_vcPreviouslySelectedIsoforms.add(strIsoform);
					
					m_bIsoformSelectionChanged = true;
				}
			});
			
			// add buttons
			Hlayout layoutH = new Hlayout();
			layoutH.setHeight("40px");
			layoutH.setStyle("text-align: center; valign: middle;");			
			layoutH.setParent(parentLayout);
			
			// add button to unselect all isoforms
			Button btnUnselectAll = new Button("Unselect all isoforms");
			btnUnselectAll.setSclass("button green");
			btnUnselectAll.setStyle("margin-left: 10px;");
			btnUnselectAll.setParent(layoutH);

			btnUnselectAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
			{
				@Override
				public void onEvent(Event event) throws Exception
				{
					m_vcValidIsoforms.clear();
					m_vcPreviouslySelectedIsoforms.clear();

					// remove all items from selection
					for(Treeitem item : m_treeSelectedIsoforms.getItems())
					{
						m_treeSelectedIsoforms.removeItemFromSelection(item);
					}
					
					m_bIsoformSelectionChanged = true;
				}
			});
			
			// add button to unselect short isoforms
			Button btnUnselectSmallIsoforms = new Button("Unselect all isoforms <=");
			btnUnselectSmallIsoforms.setSclass("button green");
			btnUnselectSmallIsoforms.setStyle("margin-left: 5px");
			btnUnselectSmallIsoforms.setParent(layoutH);

			btnUnselectSmallIsoforms.addEventListener(Events.ON_CLICK, new EventListener<Event>()
			{
				@Override
				public void onEvent(Event event) throws Exception
				{
					int nThreshold = Integer.parseInt(m_textboxMinExonThreshold.getText());
					TreeSet<String> strUnSelectedIsoforms = new TreeSet<String>();
					
					String[] pIsoforms = m_dataSupplier.GetIsoformNames();

					// find isoforms that should be unselected
					for(String strIsoform : m_vcValidIsoforms)
					{
						for(String strGeneIsoforms : pIsoforms)
						{
							if(strGeneIsoforms.equals(strIsoform))
							{								
								if(m_dataSupplier.GetExonsForIsoform(strGeneIsoforms).length <= nThreshold)
									strUnSelectedIsoforms.add(strIsoform);
							}
						}
					}
					
					// unselect isoforms
					for(String strIsoform : strUnSelectedIsoforms)
					{
						m_vcValidIsoforms.remove(strIsoform);
						m_vcPreviouslySelectedIsoforms.remove(strIsoform);
					}
					
					for(Treeitem item : m_treeSelectedIsoforms.getItems())
					{
						if(m_vcValidIsoforms.contains(item.getValue().toString()))
							m_treeSelectedIsoforms.addItemToSelection(item);
						else
							m_treeSelectedIsoforms.removeItemFromSelection(item);
					}
					
					m_bIsoformSelectionChanged = true;
				}
			});
		
			// add textbox to control which isoforms are too small
			m_textboxMinExonThreshold = new Textbox("3");
			m_textboxMinExonThreshold.setConstraint("/^[0-9]+$/");
			m_textboxMinExonThreshold.setParent(layoutH);
			m_textboxMinExonThreshold.setWidth("20px");
			m_textboxMinExonThreshold.setStyle("margin-top: 4px");
			
			Label lblText = new Label("exons");
			lblText.setStyle("display:inline-block; margin-top: 10px;");
			lblText.setParent(layoutH);
		}
	}

	/** Initializes the isoform combobox with data after a gene change */
	public void InitComboboxForIsoformSelection()
	{
		m_treeSelectedIsoforms.invalidate();
		
		// clear old list
		m_treeSelectedIsoforms.getChildren().clear();
		
		Treecols cols = new Treecols();
		Treecol col = new Treecol("Ensembl ID"); col.setParent(cols);
		col = new Treecol("RefSeq ID"); col.setParent(cols);
		col = new Treecol("Length [bp]"); col.setParent(cols);
		col = new Treecol("Exons [N]"); col.setParent(cols);
		cols.setParent(m_treeSelectedIsoforms);
		
		Treechildren children = new Treechildren();
		children.setParent(m_treeSelectedIsoforms);
		
		TreeSet<String> vcIdentifierUsed = new TreeSet<String>();
		
		String pIsoforms[] = m_dataSupplier.GetIsoformNames();
		Arrays.sort(pIsoforms);
		
		for(String strIsoform : pIsoforms)
		{
			GeneIdentifier gid = m_geneIdentifierHandler.GetGeneIdentifierForTranscript(strIsoform, m_dataSupplier.GetGene());
					
			String strIdentifier = gid.m_strEnsemblTranscriptID + "@" + gid.m_strRefGeneID + "@" + gid.m_strEntrezGeneID;

			if(!vcIdentifierUsed.contains(strIdentifier))
			{
				Treeitem item = new Treeitem();
				item.setParent(children);
				item.setCheckable(true);
				item.setValue(strIsoform);
				
				Treerow row = new Treerow();
				row.setParent(item);
				
				vcIdentifierUsed.add(strIdentifier);
				row.setId(gid.m_strEnsemblTranscriptID + "@" + gid.m_strRefGeneID + "@" + gid.m_strEntrezGeneID);
				Treecell cell = new Treecell(gid.m_strEnsemblTranscriptID); cell.setParent(row);
				cell = new Treecell(gid.m_strRefGeneID); cell.setParent(row);
				
				Exon pExons[] = m_dataSupplier.GetExonsForIsoform(strIsoform);

				int nExons = pExons.length;
				int nLength = 0;
				for(Exon ex : pExons)
				{
					nLength += ex.getLength();
				}

				cell = new Treecell(String.format(Locale.ENGLISH, "%d", nLength)); cell.setParent(row);
				cell = new Treecell(String.format(Locale.ENGLISH, "%d", nExons)); cell.setParent(row);
				
				// add to selection
				m_treeSelectedIsoforms.addItemToSelection(item);
			}
		}
		
		m_treeSelectedIsoforms.setSizedByContent(true);
	}
	
	/** Helper function that is called by HideIrrelevantIsoforms (defined below) to use multi-threading to identify invalid isoforms*/
	public TreeSet<String> RemoveIrrelevantIsoforms(String strCondition, int nMode, TreeSet<String> vcIrrelevantIsoforms, boolean bSkipFirstAndLastExon, boolean bDebug) throws IOException
	{
		// make a copy of the irrelevant isoform vector, because we do not want to change it directly
		TreeSet<String> vcResult = new TreeSet<String>();
		for(String strIsoform : vcIrrelevantIsoforms)
			vcResult.add(strIsoform);
		
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);

		//#####################################################################
		//    remove all isoforms that contain exons for which there are no
		//    splice junction reads
		//#####################################################################		
		
		ExecutorService executor = Executors.newFixedThreadPool(m_nThreads);
		for(String strIsoform : m_dataSupplier.GetIsoformNames())
		{
			if(!m_vcValidIsoforms.contains(strIsoform))
				continue;
			
			// start a thread for each isoform
			if(nMode == 0)
			{
				Runnable r = new ThreadRemoveIrrelevantIsoformsBasedOnSplitReads(strIsoform, strCondition, vcResult, mapSamplesToConditions, bSkipFirstAndLastExon, bDebug);
				executor.execute(r);
			}
			else
			{
				Runnable r = new ThreadRemoveIrrelevantIsoformsBasedOnCoverage(strIsoform, strCondition, vcResult, mapSamplesToConditions, bSkipFirstAndLastExon, bDebug);
				executor.execute(r);
			}
		}
		executor.shutdown();
		while(!executor.isTerminated()) {}
		
		return vcResult;
	}

	/**
	 *    Invokes MMSeq to calculates the isoform expression value per sample 
	 *    Returns a map of the form <isoform_id, <sample, coverage_value>>
	 */
	public void GenerateMMSeqEstimates() throws Exception
	{
		// check whether the executable is valid
		File pExecutable = new File(m_strMMSeqExecutable);
		if(!pExecutable.exists())
		{
			return;
		}
		
		DateFormat dateFormat = new SimpleDateFormat("yyyy_MM_dd_HH_mm_ss");
		Calendar cal = Calendar.getInstance();
		String strTimeStamp = dateFormat.format(cal.getTime());
		
		TreeMap<String, TreeSet<String>> vcSamplesPerCondition = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		TreeSet<String> vcIsoforms 	= new TreeSet<String>();
		String[] pIsoforms 			= m_dataSupplier.GetIsoformNames();

		/*
		// get only the selected isoforms
		for(Treeitem item : m_treeSelectedIsoforms.getSelectedItems())
		{
			String strIsoform = item.getLabel();
			
			for(String strTmp : pIsoforms)
			{
				if(strTmp.contains(strIsoform))
					vcIsoforms.add(strTmp);
			}
		}
		*/
		
		// get all isoforms
		for(String strIsoform : pIsoforms)
			vcIsoforms.add(strIsoform);

		// calculate which exon belongs to which position
		TreeMap<String, Integer> mapLengthToIsoforms = new TreeMap<String, Integer>();
		TreeMap<Integer, TreeSet<String>> mapIsoformsToPositions = new TreeMap<Integer, TreeSet<String>>();
		for(String strIsoform : pIsoforms)
		{
			Exon pExons[] = m_dataSupplier.GetExonsForIsoform(strIsoform);

			int nGeneStart = m_dataSupplier.GetGeneStart();
			int nLength = 0;
			
			for(Exon ex : pExons)
			{
				nLength += ex.getLength();
				for(int i = ex.getCodingStart(); i<ex.getCodingStop(); i++)
				{
					int nPos = i - nGeneStart;
					
					if(mapIsoformsToPositions.containsKey(nPos))
					{
						mapIsoformsToPositions.get(nPos).add(strIsoform);
					}
					else
					{
						TreeSet<String> vcTmp = new TreeSet<String>();
						vcTmp.add(strIsoform);
						mapIsoformsToPositions.put(nPos, vcTmp);
					}
				}
			}
			
			mapLengthToIsoforms.put(strIsoform, nLength);
		}

		//#######################################
		//           write hits data
		//#######################################		
		// <Isoform, <Sample, Value>>
		TreeMap<String, TreeMap<String, Double>> mapRelativeExpressionToIsoform = new TreeMap<String, TreeMap<String, Double>>();
		
		for(String strCondition : vcSamplesPerCondition.keySet())
		{
			for(String strSample : vcSamplesPerCondition.get(strCondition))
			{
				// get coverage for the current sample
				double pCoverage[] = GetCoverageForRegion(m_dataSupplier.GetReferenceName(), m_dataSupplier.GetGeneStart(), m_dataSupplier.GetGeneEnd(), strSample, true);
				
				if(pCoverage == null)
				{
					System.out.println("skipping sample: " + strSample);
					continue;
				}
				
				//######################################
				//          prepare hits file
				//######################################
				String strFileHits = m_strTmpFolder + "/" + strTimeStamp + "_" + strSample + "_" + m_dataSupplier.GetGeneID() + "_hits";
				strFileHits = strFileHits.replaceAll("\\s+", "_");
				
				//######################################
				//    write mmseq header information
				//######################################		
				PrintWriter pWriter = new PrintWriter(new FileWriter(strFileHits));
				
				for(String strIsoform : vcIsoforms)
				{
					pWriter.print("@TranscriptMetaData\t" + strIsoform + "\t" + mapLengthToIsoforms.get(strIsoform) + "\t" + mapLengthToIsoforms.get(strIsoform) + "\n");
				}
				
				pWriter.print("@GeneIsoforms\t");
				pWriter.print(m_dataSupplier.GetGeneID());
					
				for(String strIsoform : vcIsoforms)
					pWriter.print("\t" + strIsoform);
				pWriter.print("\n");
				
				
				for(int nPos : mapIsoformsToPositions.keySet())
				{
					long nCov = Math.round(pCoverage[nPos]);
					
					// add one entry per coverage
					for(int nReadID = 1; nReadID<nCov; nReadID++)
					{
						pWriter.print(">" + nPos + "_read_" + nReadID + "\n");
						
						for(String strIsoform : mapIsoformsToPositions.get(nPos))
							pWriter.print(strIsoform + "\n");
					}
				}
				pWriter.flush();
				pWriter.close();
				
				// execute MMSeq
				String strTmpOut = m_strTmpFolder + "/" + strTimeStamp + "_tmp_out";
				String strCommand = m_strMMSeqExecutable + " " + strFileHits + " " + strTmpOut;
				Utils.runProcess(strCommand, true, false, false);
				
				// obtain isoform expression value
				Scanner pScanner = new Scanner(new File(strTmpOut + ".mmseq"));
				while(pScanner.hasNextLine())
				{
					String strLine = pScanner.nextLine();
					
					// skip header
					if(strLine.contains("feature_id") || strLine.startsWith("#"))
						continue;
							
					String pSplit[] = strLine.split("\\s+");

					String strIsoform = pSplit[0];
					double fValue = Double.parseDouble(pSplit[8]);
					
					if(mapRelativeExpressionToIsoform.containsKey(strIsoform))
					{
						TreeMap<String, Double> mapExpressionToSample = mapRelativeExpressionToIsoform.get(strIsoform);
						mapExpressionToSample.put(strSample, fValue);
					}
					else
					{
						TreeMap<String, Double> mapExpressionToSample = new TreeMap<String, Double>();
						mapExpressionToSample.put(strSample, fValue);
						mapRelativeExpressionToIsoform.put(strIsoform, mapExpressionToSample);
					}
				}
				pScanner.close();
				
				// delete hits and result file
				Utils.removeFile(strFileHits);
				Utils.removeFile(strTmpOut + ".k");
				Utils.removeFile(strTmpOut + ".M");
				Utils.removeFile(strTmpOut + ".mmseq");
				Utils.removeFile(strTmpOut + ".prop.trace_gibbs.gz");
				Utils.removeFile(strTmpOut + ".trace_gibbs.gz");
				Utils.removeFile(strTmpOut + ".identical.trace_gibbs.gz");
				Utils.removeFile(strTmpOut + ".identical.mmseq");
				Utils.removeFile(strTmpOut + ".gene.trace_gibbs.gz");
				Utils.removeFile(strTmpOut + ".gene.mmseq");
			}
		}

		m_projectModel.AddIsoformExpressionToDataBase(mapRelativeExpressionToIsoform, m_dataSupplier.GetGeneID(), m_strPathMMSeqResults);
	}

	/** Calculates the sample size factors */
	public boolean CalculateSizeFactors(String strFileGTF, String strFileProject) throws IOException
	{
		// init project
		try
		{
			if(!m_projectModel.Init(strFileProject, -1, -1, -1.0, -1, true, true))
				return false;
		}
		catch(Exception e)
		{
			ErrorMessage.ShowError("failed to process project file: " + strFileProject);
			System.out.println(ExceptionUtils.getStackTrace(e));
			return false;
		}
		
		TreeMap<String, String> mapBigWigFiles 	= m_projectModel.GetBigWigFilesForSamples();
		TreeSet<String> vcSamples = m_projectModel.GetSamples(); 
		int nSamples = vcSamples.size();

		if(this.m_nReferenceType == REFFLAT_REFERENCE_FILE)
		{
			ErrorMessage.ShowError("can't calculate size factors with refflat file -> GTF format required");
			return false;
		}
		
		// get a list of all genes
		TreeSet<Gene> vcGenes = new TreeSet<Gene>();
		GTFParser parser = new GTFParser(strFileGTF);
		while(true)
		{
			GTFGene gtf_gene = parser.nextGene();
			if(gtf_gene == null)
				break;
			
			vcGenes.add(gtf_gene.createGene());
		}
		
		Vector<Vector<Double>> vcQuotient = new Vector<Vector<Double>>();
		for(int i=0; i<nSamples; i++)
			vcQuotient.add(new Vector<Double>());

		// <gene, <sample, coverage>>
		TreeMap<Gene, double[]> mapCoverageToGenePerSample = new TreeMap<Gene, double[]>();
		
		int nSampleIdx = 0;
		for(String strSample : vcSamples)
		{
			System.out.println("calculating average gene coverage for: " + strSample);
			BigWigReader reader = null;
			reader = new BigWigReader(mapBigWigFiles.get(strSample));

			for(Gene gene : vcGenes)
			{		
				// get coverage for the whole gene
				int nStart 	= gene.getStart();
				int nEnd	= gene.getStop();
				BigWigIterator it = reader.getBigWigIterator(gene.getChromosome(), nStart-1, gene.getChromosome(), nEnd, false);
				
				double pValues[] = new double[nEnd-nStart+1];
				while(it.hasNext())
				{
					WigItem item = it.next();
					int nIdx = item.getEndBase() - nStart;
					pValues[nIdx] = item.getWigValue();				
				}
				
				// get coverage values for all exon positions
				Vector<Double> vcCoverage = new Vector<Double>();
				
				ExonGroup grps[] = gene.computeOverlappingExonGroups();
				
				for(ExonGroup grp : grps)
				{
					int nExStart= grp.getGenomicStartOfGroup();
					int nExEnd = grp.getGenomicStopOfGroup();
					
					for(int i=nExStart; i<nExEnd; i++)
					{
						int nPos = i - nStart;
						
						vcCoverage.add(pValues[nPos]);
					}
				}
				
				double pCoverage[] = new double[vcCoverage.size()];
				int nIdx = 0;
				for(Double fVal : vcCoverage)
				{
					pCoverage[nIdx] = fVal;
					nIdx++;
				}

				if(mapCoverageToGenePerSample.containsKey(gene))
				{
					double vcGeneReadCounts[] = mapCoverageToGenePerSample.get(gene);
					vcGeneReadCounts[nSampleIdx] = StatUtils.mean(pCoverage);
				}
				else
				{
					double vcGeneReadCounts[] = new double[nSamples];
					vcGeneReadCounts[nSampleIdx] = StatUtils.mean(pCoverage);
					mapCoverageToGenePerSample.put(gene, vcGeneReadCounts);
				}				
			}
			
			nSampleIdx += 1;
			reader.close();
		}

		int nGenes = 0;
		for(Gene gene : mapCoverageToGenePerSample.keySet())
		{			
			double vcGeneReadCounts[] = mapCoverageToGenePerSample.get(gene);
			
			double fMean = StatUtils.geometricMean(vcGeneReadCounts);
			
			Vector<Double> vcCurrentQuotients = new Vector<Double>();
			
			// calculate per sample difference to geometric mean
			for(int i=0; i<nSamples; i++)
			{
				if(fMean != 0)
				{
					vcCurrentQuotients.add(vcGeneReadCounts[i] / fMean);
				}
			}
			
			// test fMean to skip unexpressed genes
			if(fMean != 0 && vcCurrentQuotients.size() == nSamples)
			{
				for(int i=0; i<nSamples; i++)
				{
					vcQuotient.get(i).add(vcCurrentQuotients.get(i));
				}
				
				nGenes += 1;
			}
		}
		
		for(int i=0; i<nSamples; i++)
		{
			double[] pVals = new double[nGenes];
			for(int j=0; j<nGenes; j++)
				pVals[j] = vcQuotient.get(i).get(j);
			
			double fSizeFactor = StatUtils.percentile(pVals, 50.0);
			
			String strSampleName = null;
			int nIdx = 0;
			
			Iterator<String> it = vcSamples.iterator();
			int k = -1;
			while(it.hasNext() && k < i)
			{
				strSampleName = it.next();
				k++;
			}
			
			if(strSampleName == null)
				System.out.println("ERROR: invalid sample " + i + " " +  nIdx);

			m_projectModel.SetSizeFactor(strSampleName, fSizeFactor);
		}
		
		m_projectModel.WriteSizeFactors();
		
		return true;
	}

	/**
	 *    Identifies isoforms that do not pass the filter criteria and flags them as irrelevant.
	 *     
 	 *    Hiding irrelevant isoforms is done in two steps:
	 *
	 *    First, isoforms containing exons that are not supported by junction reads are removed.
	 *
	 *    Second, the coverage for the exons of the remaining isoforms are calculate and
	 *    isoforms with too low coverage per base or coverage length are discarded.
	 */
	public boolean HideIrrelevantIsoforms(boolean bSkipFirstAndLastExon, boolean bDebug) throws IOException
	{
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		int nConditions = mapSamplesToConditions.size();
		
		//##################################
		//    remove irrelevant isoforms
		//##################################
		TreeMap<String, TreeSet<String>> mapIrrelevantIsoforms = new TreeMap<String, TreeSet<String>>();
		
		// repeat twice, first remove isoforms based on junction reads, then recalculate exon coverages
		// then remove isoforms based on coverage values
		TreeSet<String> vcIrrelevantIsoforms = new TreeSet<String>();
		 
		for(int i=0; i<2; i++)
		{	
			for(String strCondition : mapSamplesToConditions.keySet())
			{
				TreeSet<String> vcIrrelevantIsoformsTmp = RemoveIrrelevantIsoforms(strCondition, i, vcIrrelevantIsoforms, bSkipFirstAndLastExon, bDebug);
				if(vcIrrelevantIsoformsTmp == null)
					return false;

				for(String strIsoform : vcIrrelevantIsoformsTmp)
				{
					if(mapIrrelevantIsoforms.containsKey(strIsoform))
					{
						mapIrrelevantIsoforms.get(strIsoform).add(strCondition);
					}
					else
					{
						TreeSet<String> vcInvalidInConditions = new TreeSet<String>();
						vcInvalidInConditions.add(strCondition);
						mapIrrelevantIsoforms.put(strIsoform, vcInvalidInConditions);
					}
				}
			}
			
			for(Map.Entry<String, TreeSet<String>> e : mapIrrelevantIsoforms.entrySet())
			{
				if(e.getValue().size() == nConditions)
				{
					m_vcValidIsoforms.remove(e.getKey());
					vcIrrelevantIsoforms.add(e.getKey());
				}
			}

			// recalculate exon expression based on newly selected exon groups			
			m_dataSupplier.RetrieveCoverageData(m_hWindow);
		}
		
		if(bDebug)
		{
			System.out.println("############# irrelevant isoforms ##############");
			System.out.println(mapIrrelevantIsoforms);
		}

		//######################################################################################
		//    while the above hides transcripts based on raw coverages (i.e. they are taken
		//    directly form the bigwig file), the code below removes exons that do not have
		//    enough coverage after exon specific coverages have been assigned.
		//    In principle, the code above ignores the fact that different exons might be
		//    overlapping, while the code below adjusts coverages for overlapping exons and
		//    then removes exons based on those values.
		//######################################################################################
		
		String pIsoforms[] = m_dataSupplier.GetIsoformNames();

		// remove isoforms that include exons with too low coverage
		for(String strIsoform : pIsoforms)
		{
			Exon pExons[] = m_dataSupplier.GetExonsForIsoform(strIsoform);
			
			int nExonIdx = 0;
			for(Exon exon : pExons)
			{
				boolean bIsFirstExon = false;
				boolean bIsLastExon	 = false;
				
				if(nExonIdx == 0)
					bIsFirstExon = true;

				if(nExonIdx == pExons.length-1)
					bIsLastExon = true;
				
				if(bSkipFirstAndLastExon && (bIsFirstExon || bIsLastExon))
					continue;

				int nInvalidInConditions = 0;
				
				TreeMap<String, Double> mapExonCoverageToSample = new TreeMap<String, Double>();
				if(bIsFirstExon || bIsLastExon)
				{					
					ExonGroup[] grps = m_dataSupplier.GetExonGroups();
					for(ExonGroup grp : grps)
					{
						if(grp.groupContainsExon(exon.getExonID()))
						{
							for(Exon curExon : grp.getExons())
							{
								// Only merge first with other first exons and last with other last exons
								String vcIsoforms[] = m_dataSupplier.GetIsoformNames();
								boolean bOkay = false;
								for(String strCurrentIsoform : vcIsoforms)
								{
									if(!vcIrrelevantIsoforms.contains(strCurrentIsoform))
									{
										Exon pIsoformExons[] = m_dataSupplier.GetExonsForIsoform(strCurrentIsoform);
										if((pIsoformExons[0] == curExon && bIsFirstExon) || (pIsoformExons[pIsoformExons.length-1] ==curExon && bIsLastExon))
										{
											// exons must also overlap
											if(curExon.getCodingStart() <= exon.getCodingStop() && curExon.getCodingStop() >= exon.getCodingStart())
											{
												bOkay=true;
											}
										}
										// or the first/last exon must be fully contained in another exon
										else if(curExon.getCodingStart() <= exon.getCodingStart() && curExon.getCodingStop() >= exon.getCodingStop())
										{
											bOkay=true;
										}
										
									}
								}
								
								if(bOkay)
								{
									TreeMap<String, Double> mapExonCoverageToSampleTmp = m_dataSupplier.GetAbsoluteCoveragePerSample(m_hWindow, curExon);

									// sum up coverages
									for(String strSample : mapExonCoverageToSampleTmp.keySet())
									{
										if(mapExonCoverageToSample.containsKey(strSample))
											mapExonCoverageToSample.put(strSample, mapExonCoverageToSample.get(strSample) + mapExonCoverageToSampleTmp.get(strSample));
										else
											mapExonCoverageToSample.put(strSample, mapExonCoverageToSampleTmp.get(strSample));
									}
								}
							}
						}
					}
					
					nExonIdx++;
				}
				else
				{
					mapExonCoverageToSample = m_dataSupplier.GetAbsoluteCoveragePerSample(m_hWindow, exon);
				}
				
				// get exon coverage
				double[] pCoverage = null;
				
				// test each condition
				for(String strCondition : mapSamplesToConditions.keySet())
				{
					int nSampleIdx = 0;
					// get number of valid samples for this condition
					for(String strSample : mapSamplesToConditions.get(strCondition))
					{
						if(m_vcSelectedSamples.contains(strSample))
							nSampleIdx++;
					}
					
					pCoverage = new double[nSampleIdx];
					
					nSampleIdx = 0;
					for(String strSample : mapSamplesToConditions.get(strCondition))
					{
						if(mapExonCoverageToSample.containsKey(strSample))
						{
							pCoverage[nSampleIdx] = mapExonCoverageToSample.get(strSample);
							nSampleIdx++;
						}
					}
					
					// calculate how many samples are invalid				
					double fMedianCov = StatUtils.percentile(pCoverage, 50.0);

					if(fMedianCov < m_nMinCovPerBase)
					{
						nInvalidInConditions += 1;
					}
				}
				
				// remove if invalid in all conditions
				if(nInvalidInConditions == mapSamplesToConditions.size())
				{
					if(bDebug)
					{
						System.out.println("removing isoform due to insufficient exon coverage: " + strIsoform + " (" + exon.getCodingStart() + "-" + exon.getCodingStop()  + ")");
					}
					m_vcValidIsoforms.remove(strIsoform);
					break;
				}
			}
		}

		// only adjust the isoform selection if the GUI was initialized
		if(m_treeSelectedIsoforms != null)
		{
			m_vcPreviouslySelectedIsoforms.clear();
			for(Treeitem item : m_treeSelectedIsoforms.getItems())
			{
				if(m_vcValidIsoforms.contains(item.getValue().toString()))
				{
					m_vcPreviouslySelectedIsoforms.add(item.getValue().toString());
					m_treeSelectedIsoforms.addItemToSelection(item);
				}
				else
					m_treeSelectedIsoforms.removeItemFromSelection(item);
			}
		}

		return true;
	}

	/**
	 *     Adjust PSI score p-values by Benjamini-Hochberg multiple-testing correction
	 *     - This function would also fit else where better
	 */
	protected TreeSet<SimpleSpliceScore> AdjustPValuesPsiScores(TreeSet<SimpleSpliceScore> vcScores, Vector<Double> vcPValues)
	{
		int nTests = vcPValues.size();
		
		// sort p-values ascendingly
		Collections.sort(vcPValues);
		
		// find the last significant hit with highest p-value
		int nLastSignificantHit = -1;
		for(int i=0; i<nTests; i++)
		{
			double fCriticalValue = ((double)i / (double)nTests) * 0.05;
			double fPValue = vcPValues.get(i);
	
			if(fPValue < fCriticalValue)
				nLastSignificantHit = i;
		}
		
		// create a set of significant p-values
		TreeSet<Double> vcSignificantPValues = new TreeSet<Double>();
		int i=0;
		for(double fPValue : vcPValues)
		{			
			if(i > nLastSignificantHit)
				break;
			i++;

			vcSignificantPValues.add(fPValue);
		}
		
		// add significance flag
		for(SimpleSpliceScore score : vcScores)
		{
			if(vcSignificantPValues.contains(score.m_fPValue))
			{
				score.m_bSignificant = true;
			}
		}
		
		return vcScores;
	}
	
	/** Sets the application parameters, which is only necessary for the console application. */
	public void SetParameters(int nMinJunctionReads, int nMinCovPerBase, double fMinFracCoveredBases, double fMinRatio, int nThreads, String strFileXref, String strFileGTF, boolean bDetectIntronRetention)
	{
		m_nMinJunctionReads					= nMinJunctionReads;
		m_nMinCovPerBase					= nMinCovPerBase;
		m_fMinCoveredBases					= fMinFracCoveredBases;
		m_fVariableExonThreshold			= fMinRatio;
		m_nThreads							= nThreads;
		
		m_strFileGTF						= strFileGTF;
		m_nReferenceType 					= GTF_REFERENCE_FILE;
		
		m_geneIdentifierHandler				= new GeneIdentifierHandler();
		
		m_geneIdentifierHandler.Init(strFileXref);
		m_geneIdentifierHandler.EnableAll();
		
		m_bDetectIntronRetentionEvents	= bDetectIntronRetention;
		
		if(bDetectIntronRetention)
			System.out.println(" - intron retention events will not be discovered.");
	}

	/** Loads the gene annotation file */
	public void LoadGeneAnnotation(String strFile)
	{
		if(strFile.toLowerCase().endsWith(".gtf") || strFile.toLowerCase().endsWith(".gff"))
		{
			m_strFileGTF = strFile;
			String strFileIdx = m_strFileGTF + ".idx";
			File pFile = new File(strFileIdx);
			
			m_nReferenceType = GTF_REFERENCE_FILE;

			// check whether the file is accessible
			if(m_gffReader == null || !m_gffReader.GetFileName().equals(m_strFileGTF))
			{
				try
				{
					// set the GTF/GFF reader member variable
					m_gffReader = new RandomAccessGFFReader(new File(m_strFileGTF), pFile);
				}
				catch(Exception e)
				{
					ErrorMessage.ShowError("failed to open GFF file: " + m_strFileGTF);
					e.printStackTrace();
					return;
				}
			}

			// Now open the index file for the GTF file and collect all gene IDs
			// Missing GTF index files are automatically generated by the RandomAccessGFFReader
			Scanner pIn = null;
			try
			{
				pIn = new Scanner(pFile);
			}
			catch(Exception e)
			{
				ErrorMessage.ShowError("failed to open file: " + strFileIdx);
				e.printStackTrace();
				return;
			}
			
			TreeSet<String> vcGeneIDs = new TreeSet<String>();
			while(pIn.hasNextLine())
			{
				String strLine = pIn.nextLine();
				vcGeneIDs.add(strLine.split("\\s+")[0]);
			}
			pIn.close();
			
			// Now check which genes in the GTF are also included in the cross reference file that was used to generate the list of gene identifiers
			for(GeneIdentifier gid : m_geneIdentifierHandler.GetAllGeneIdentifiers())
			{
				// valid gene identifiers require an ensembl or entrez gene ID
				if(vcGeneIDs.contains(gid.m_strEnsemblGeneID) || vcGeneIDs.contains(gid.m_strEntrezGeneID))
					gid.m_bIsValid = true;
				else
					gid.m_bIsValid = false;
			}
		}
		else
		{
			ErrorMessage.ShowError("ERROR: Only GTF files are supported.");
		}
	}

	/** Retrieves the coverage for an exonic part from tissue specific data */
	public TreeMap<String, TreeMap<String, Vector<Double>>> GetGTEXDataForExonicParts()
	{
		TreeMap<String, TreeMap<String, Vector<Double>>> mapCountsPerExonicPartAndTissue = new TreeMap<String, TreeMap<String, Vector<Double>>>();
	
		GeneIdentifier id = m_dataSupplier.GetGeneIdentifier();
		
		if(id == null)
			return null;
		
		String strEntrezID = id.m_strEntrezGeneID;
		
		String strGTEXID = "?";

		// open file with DEXSeq exon IDs
		File pFileDEXSeqIDs = new File(m_strFileDEXSeqIDs);
		
		Scanner pScanner = null;
		try
		{
			pScanner = new Scanner(pFileDEXSeqIDs);
		}
		catch(Exception e)
		{
			ErrorMessage.ShowError("failed to topen file: " + m_strFileDEXSeqIDs);
			e.printStackTrace();
			return null;
		}
		
		// map exons to DEXSeq exon IDs
		TreeMap<String, String> mapIDsToPositions = new TreeMap<String, String>(); 
		while(pScanner.hasNextLine())
		{
			String strLine = pScanner.nextLine();
			
			String pSplit[] = strLine.split("\t");
			
			if(pSplit.length < 3)
			{
				System.out.println("Warning: invalid line: " + strLine);
				continue;
			}
			
			if(pSplit[2].equals("aggregate_gene"))
			{
				String strIDs = pSplit[8].split("\"")[1];
				String pGeneIDs[] = strIDs.split("\\+");
				
				boolean bOkay = false;
				for(String strID : pGeneIDs)
				{
					if(strID.equals(strEntrezID))
					{
						strGTEXID = strIDs;
						bOkay=true;
						break;
					}
				}

				if(bOkay)
				{
					while(true)
					{
						if(!pScanner.hasNextLine())
							break;

						strLine = pScanner.nextLine();
						String pSplit2[] = strLine.split("\t");
						
						if(!pSplit2[2].equals("exonic_part"))
							break;
						
						int nStart = Integer.parseInt(pSplit2[3]);
						int nEnd   = Integer.parseInt(pSplit2[4]);
						
						String pSplit3[] = pSplit2[8].split(";");
						String strID = pSplit3[1].split("\"")[1]; 

						mapIDsToPositions.put(nStart + "-" + nEnd, strID);
					}
				}
			}
		}
		pScanner.close();
		
		//#######################################
		//    get read counts per exonic part
		//#######################################
		File pFolder = new File(m_strPathToExonData);
		if(!pFolder.exists())
		{
			ErrorMessage.ShowError("ERROR: Invalid GTEX folder specified: " + m_strPathToExonData);
			return null;
		}
		
		for(String strPos : mapIDsToPositions.keySet())
		{
			String strFile = m_strPathToExonData + "//GTEX_" + strGTEXID + "_" + mapIDsToPositions.get(strPos) + ".tsv";
			File pFile = new File(strFile);
			
			if(!pFile.exists())
			{
				ErrorMessage.ShowError("ERROR: No GTEX data available for exonic part: " + strGTEXID + " " + strPos + " " + mapIDsToPositions.get(strPos) + " -> " + strFile);
				continue;
			}
			
			try
			{
				pScanner = new Scanner(pFile);
			}
			catch(Exception e)
			{
				ErrorMessage.ShowError("failed to open file: " + strFile);
				e.printStackTrace();
			}
			
			// skip first line (header)
			if(pScanner.hasNextLine())
				pScanner.nextLine();
			
			TreeMap<String, Vector<Double>> mapCountsPerTissue = new TreeMap<String, Vector<Double>>();
			// read data
			while(pScanner.hasNextLine())
			{
				String strLine = pScanner.nextLine();
				
				String split[] = strLine.split("\t");
				
				String strTissue = split[1];
				double nCount 	 = Double.parseDouble(split[2]);
				
				if(mapCountsPerTissue.containsKey(strTissue))
				{
					Vector<Double> vcCounts = mapCountsPerTissue.get(strTissue);
					vcCounts.add(nCount);
				}
				else
				{
					Vector<Double> vcCounts = new Vector<Double>();
					vcCounts.add(nCount);
					mapCountsPerTissue.put(strTissue, vcCounts);
				}
			}
			pScanner.close();
			
			mapCountsPerExonicPartAndTissue.put(strPos, mapCountsPerTissue);
		}
		
		return mapCountsPerExonicPartAndTissue;
	}	
	
	/** Limits the number of initially selected samples for large projects (>50 samples) */
	public void ReduceSampleSize()
	{
		int nMaxTotalSamples = 50;
		
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);	
		m_treeSelectedSamples.clearSelection();				
		m_vcPreviouslySelectedSamples.clear();
		m_vcSelectedSamples.clear();
		
		int nTotalSamples = m_projectModel.GetSamples().size();
		
		TreeSet<String> vcSelectedItems = new TreeSet<String>();
		
		if(nTotalSamples > nMaxTotalSamples)
		{
			int nMaxSamplesPerCondition = nMaxTotalSamples / mapSamplesToConditions.keySet().size();

			for(String strCondition : mapSamplesToConditions.keySet())
			{				
				TreeSet<String> vcSamples = mapSamplesToConditions.get(strCondition);
				if(vcSamples.size() < nMaxSamplesPerCondition)
				{
					for(String strSample : vcSamples)
					{
						vcSelectedItems.add(strSample);
					}
				}
				else
				{
					// index all samples
					TreeSet<Integer> vcIndices = new TreeSet<Integer>();
					for(int i=0; i<vcSamples.size(); i++)
						vcIndices.add(i);
					
					// remove random indices until the number of samples is okay
					while(vcIndices.size() > nMaxSamplesPerCondition)
					{
						int nIdx = (int)(Math.random() * vcIndices.size());
						
						int nCurrent = 0;
						for(int n : vcIndices)
						{
							if(nCurrent == nIdx)
							{
								vcIndices.remove(n);
								break;
							}
							
							nCurrent++;
						}
					}
					
					// get an array of the treeset
					String[] pSamples = new String[vcSamples.size()];
					int nIdx = 0;
					for(String strSample : vcSamples)
					{
						pSamples[nIdx] = strSample;
						nIdx++;
					}
					
					// select remaining samples
					for(int idx : vcIndices)
					{
						vcSelectedItems.add(pSamples[idx]);
					}
				}
			}
		}
		else
		{
			for(String strSample : m_projectModel.GetSamples())
			{
				vcSelectedItems.add(strSample);
			}
			
			for(String strCondition : mapSamplesToConditions.keySet())
			{
				vcSelectedItems.add(strCondition);
			}
		}
		
		// update selection
		for(String strItem : vcSelectedItems)
		{			
			m_vcSelectedSamples.add(strItem);

			for(Treeitem item : m_treeSelectedSamples.getItems())
			{
				if(item.getLabel().equals(strItem))
				{
					m_treeSelectedSamples.addItemToSelection(item);
				}
			}
		}

		UpdateComboboxForSampleSelection();
	}

	/**
	 *    Helper function that opens a project file.
	 *    - only invoked by the console application (AnalyzerGeneFactory)
	 */
	public void InitProjectModel(String strFileProject)
	{
		try
		{
			m_projectModel.Init(strFileProject, -1, -1, -1.0, -1, true, false);
		}
		catch(Exception e)
		{
			ErrorMessage.ShowError("Failed to open project file: " + strFileProject);
			return;
		}
	}

	/**
	 *    Saves the list of 'permanent' alterantive splicing events to a file.
	 *    - Called by pressing the "Save changes" button located below the result list
	 */
	public void SaveSplicingResults()
	{
		String strResultFile = m_projectModel.GetProjectName() + "_splicing_results.dat";
		m_vcASResults.SaveToFile(m_strPathHitLists, strResultFile);
	}
	
	/**
	 *    Writes the size in bytes to the output file and then adds the binary data
	 *    or, if there is no data to be written, adds a byte size of 0 to the output file
	 */
	public static final void WriteStringToFileOutputStream(String strOutputString, FileOutputStream pOut) throws IOException
	{
		if(!strOutputString.isEmpty())
		{
			byte pBytes[] = strOutputString.getBytes("UTF-8");
			
			// write string length
			ByteBuffer bb = ByteBuffer.allocate(Integer.BYTES);
			bb.putInt(pBytes.length);
			pOut.write(bb.array());
			
			// write string
			pOut.write(pBytes);
		}
		else
		{
			// write string length
			ByteBuffer bb = ByteBuffer.allocate(Integer.BYTES);
			bb.putInt(0);
			pOut.write(bb.array());
		}
	}
	
	/**
	 *    Reads data from a binary file. First, the length of the data is read to reserve sufficient memory, then the data is read to a buffer
	 *    and converted to a UTF-8 string. If the data has a length of 0 bytes an empty string is returned.
	 */
	public static final String ReadStringFromFileInputStream(FileInputStream pIn) throws IOException
	{
		byte pIntBuffer[] = new byte[Integer.BYTES];
		if(pIn.read(pIntBuffer) == -1) return "";
		ByteBuffer bb = ByteBuffer.wrap(pIntBuffer);
		int nLength = bb.getInt();
		
		String strRes = "";
		if(nLength > 0)
		{
			byte pBytes[] = new byte[nLength];
			if(pIn.read(pBytes) == -1) return "";
			strRes = new String(pBytes, Charset.forName("UTF-8"));
		}
		
		return strRes;
	}

	@Override
	public void afterCompose()
	{
		m_pagingProjects.setPageSize(PROJECT_PAGING_SIZE);
		m_pagingProjects.setTotalSize(m_vcProjectInfos.size());		
	}

}
