package Manananggal;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.LinearGradientPaint;
import java.awt.MultipleGradientPaint;
import java.awt.MultipleGradientPaint.CycleMethod;
import java.awt.RadialGradientPaint;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.sql.SQLException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.imageio.ImageIO;

import org.apache.commons.lang3.exception.ExceptionUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.OneWayAnova;
import org.apache.commons.math3.stat.inference.TTest;
import org.apache.log4j.BasicConfigurator;
import org.broad.igv.bbfile.BigWigIterator;
import org.broad.igv.bbfile.WigItem;
import org.zkoss.zhtml.Filedownload;
import org.zkoss.zk.ui.Component;
import org.zkoss.zk.ui.Components;
import org.zkoss.zk.ui.Executions;
import org.zkoss.zk.ui.event.ClientInfoEvent;
import org.zkoss.zk.ui.event.Event;
import org.zkoss.zk.ui.event.EventListener;
import org.zkoss.zk.ui.event.Events;
import org.zkoss.zk.ui.event.MouseEvent;
import org.zkoss.zul.Area;
import org.zkoss.zul.Bandbox;
import org.zkoss.zul.Bandpopup;
import org.zkoss.zul.Button;
import org.zkoss.zul.Cell;
import org.zkoss.zul.Checkbox;
import org.zkoss.zul.Column;
import org.zkoss.zul.Columns;
import org.zkoss.zul.Combobox;
import org.zkoss.zul.Comboitem;
import org.zkoss.zul.East;
import org.zkoss.zul.Grid;
import org.zkoss.zul.Groupbox;
import org.zkoss.zul.Hlayout;
import org.zkoss.zul.Html;
import org.zkoss.zul.Image;
import org.zkoss.zul.Imagemap;
import org.zkoss.zul.Borderlayout;
import org.zkoss.zul.Label;
import org.zkoss.zul.Listbox;
import org.zkoss.zul.Listcell;
import org.zkoss.zul.Listhead;
import org.zkoss.zul.Listheader;
import org.zkoss.zul.Listitem;
import org.zkoss.zul.Messagebox;
import org.zkoss.zul.North;
import org.zkoss.zul.Center;
import org.zkoss.zul.Popup;
import org.zkoss.zul.Radio;
import org.zkoss.zul.Radiogroup;
import org.zkoss.zul.Row;
import org.zkoss.zul.Rows;
import org.zkoss.zul.South;
import org.zkoss.zul.Tab;
import org.zkoss.zul.Tabbox;
import org.zkoss.zul.Tabpanel;
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
import org.zkoss.zul.West;
import org.zkoss.zul.Window;

import BioKit.Exon;
import BioKit.ExonGroup;
import BioKit.GTFGene;
import BioKit.GTFParser;
import BioKit.Gene;
import BioKit.RandomAccessGFFReader;
import BioKit.Utils;
import Manananggal.AlternativeSplicingHit.AlternativeSplicingExon;

// <?link rel="stylesheet" type="text/css" href="/bootstrap/v3/bootstrap/css/bootstrap.css" ?>

public class SplicingWebApp extends Window
{
	private int				m_nThreads;
	private String			m_strPathReferences;
	private String			m_strPathInput;
	private String			m_strPathGeneGTEX;
	private String			m_strPathExonGTEX;
	private String			m_strFileDEXSeqIDs;
	private String			m_strPathHitLists;
	private String			m_strPathSplicingHitList;
	private String			m_strPathMMSeqResults;
	private String			m_strMMSeqExecutable;
	private String			m_strFileGTF;
	private String			m_strTmpFolder;
	private int				m_nReferenceType;
	
	static final int		GTF_REFERENCE_FILE 		= 1;
	static final int		REFFLAT_REFERENCE_FILE 	= 2;
	
	private Gene 			m_Gene;
	private Exon[] 			m_pExons;
	private ExonGroup[] 	m_pExonGroups;
	private Vector<GeneIdentifier> m_vcGeneIdentifier;
	
	// this layout defines regions for buttons/options (north), the splice graph (west) and main plots (center)
	private Borderlayout 			m_layout;
	private North					m_layoutNorth;
	private West					m_layoutWest;
	private Center					m_layoutCenter;
	
	private int						m_nClientWindowWidth;
	private int						m_nClientWindowHeight;
	
	private int 					m_nMaxWidth;				// standard window size
	private int						m_nMinCovPerBase;			// Used as threshold to exclude exons with low coverage
	private double					m_fMinCoveredBases;			// Used as threshold to exclude exons that are only partially covered
	private int						m_nMinJunctionReads; 		// minimum coverage for junctions
	private double					m_fVariableExonThreshold;	// threshold at which exons are designated as differentially expressed (e.g. 0.3 if there must be a 30% difference in relative coverage)
	
	private Textbox					m_textboxMinCovPerBase;
	private Textbox					m_textboxMinCoveredBases;
	private Textbox					m_textboxMinJunctionReads;
	private Textbox					m_textboxVariableExonThreshold;
	
	private Textbox					m_textboxMinExonThreshold;	// threshold for the removal of small isoforms (viewer only)

	// this layout is for the main plots and will be used for the isoform, coverage and new junction plots
	private Borderlayout			m_layoutMainPlots;
	private North					m_layoutMainTop;
	private South					m_layoutMainSouth;
	private Vlayout					m_layoutMainBottom;
	
	// handle to application
	SplicingWebApp 					m_hWindow;
	
	// popup window
	private Window 				m_windowPopup;
	
	// stores information when clicking an element in the isoform plot
	private ClickEvent 			m_ClickEvent;
	
	// grants access to the SQL junction and exon count data
	private ProjectModel 		m_projectModel;
	
	// available options
	private boolean				m_bSkipFirstAndLastExon;
	private boolean				m_bShowUniqueFeatures;
	private boolean				m_bShowAmbigiousUniqueExons;
	private boolean				m_bCoverageRequiresRedraw;
	private boolean				m_bColorExonsByRelativeCoverage;
	private boolean				m_bColorJunctionPaths;
	private boolean				m_bShowRelativeCoverageNumbersJunctions;
	private boolean				m_bShowASEvents;
	private boolean				m_bShowSecondCoveragePlot;
	private boolean				m_bShowRelativeCoverageNumbersExons;
	private boolean				m_bUseRelativeExonCoverageMode1;
	private boolean				m_bUseRelativeExonCoverageMode2;
	private boolean				m_bColorExonsAndJunctionsByCoverage;
	
	@SuppressWarnings("unused")
	private int					m_nOrganismType;	// 0 = human, 1 = mouse
	
	private String				m_strSelectedCondition;
	private Combobox			m_comboboxSelectedCondition;
	private Combobox			m_comboboxSelectedConditionType;
	private Combobox			m_comboboxSelectedGeneAnnotation;
	private String				m_strSelectedConditionType;
	private Bandbox				m_bandboxSelectedGene;
	private Listbox				m_listboxSelectedGene;
	private Tree				m_treeSelectedIsoforms;	
	
	private Checkbox 			m_checkboxUseMedian;
	private Checkbox			m_checkboxUseGeometricMean;
	private Checkbox 			m_checkboxUseMean;
	private Checkbox 			m_checkboxShowRelativeCoverage;
	private Checkbox 			m_checkboxUseLog2;

	private TreeSet<String>		m_vcSelectedSamples;	
	private boolean				m_bUseReducedDataSet; 
	private Checkbox			m_checkboxUseReducedDataSet;
	private Checkbox			m_checkboxSkipFirstAndLastExon;
	
	private boolean 				m_bHideIrrelevantIsoforms;
	private boolean					m_bIsoformSelectionChanged;
	private TreeSet<String>			m_vcValidIsoforms;
	private TreeMap<String, int[]> 	m_mapExonsToIsoforms;
	
	private Combobox			m_comboboxSelectedEntropyIndex;
	private boolean				m_bShowEntropyData;	// Gini index, Gini-Simpson index, Theil index, Atkinson, Generalized entropy index, 
	
	private Imagemap 			m_imgMapCoverage;
	private Imagemap 			m_imgMapIsoforms;
	private Imagemap 			m_imgMapColors;
	
	private boolean				m_bSampleSelectionChanged;
	private Tree				m_treeSelectedSamples;
	private Groupbox 			m_ColorSelectionGroupBox;
	private Grid				m_gridHitList;
	private Grid				m_gridPSIHitList;
	
	private Checkbox			m_checkboxRelativeCoverage1;
	private Checkbox			m_checkboxRelativeCoverage2;	
	
	// only filled when experiment specific data was processed
	private TreeMap<Exon, TreeMap<String, Double>> 		  m_mapExonCoveragePerSample;
	private TreeMap<Exon, TreeMap<String, Double>>   	  m_mapExonCoveragePerSampleRawFromBigWig;
	private TreeMap<Exon, TreeMap<String, Double>>   	  m_mapExonCoverageLengthPerSampleRawFromBigWig;
	private TreeMap<Exon, TreeMap<String, TreeSet<Exon>>> m_mapIndistinguishableExonPerSample;
	private TreeMap<String, int[]>						  m_mapBigWigCoverageToSamples;

	// This map stores a set of junctions for each isoform representing the path between two exons.
	// Each junction set is stored in a TreeMap using a unique path ID
	// junctions are stored in the format (start-end).
	// <isoform, set_of_junctions>
	private TreeMap<String, TreeMap<Integer, TreeSet<String>>> 	m_mapAnalysedJunctionPaths;
	
	// The relative coverage for each path and sample is stored in this map: <sample, <path_id, ratio>>
	private TreeMap<String, TreeMap<Integer, Double>> 			m_mapCoveragePerJunctionPathAndSample;

	private TreeMap<String, Double> 							m_mapEntropyToExonicPart;
	private TreeMap<String, TreeMap<String, Vector<Double>>> 	m_mapCountsPerExonAndTissue;	
	
	// values important for the popup windows	
	private boolean							m_bCoveragePlotUseLog2;
	private boolean							m_bCoveragePlotShowMean;
	private boolean							m_bCoveragePlotShowGeometricMean;
	private boolean							m_bCoveragePlotShowMedian;
	private boolean							m_bCoveragePlotShowRelativeCoverage;
	private boolean							m_bCoveragePlotShowQuartiles;
	private TreeSet<String> 				m_vcPreviouslySelectedTreeItems;
	private TreeSet<String>					m_vcPreviouslySelectedGTEXTissues;
	private TreeSet<AlternativeSplicingHit> m_vcASHits;
	private AlternativeSplicingHit			m_vcASCurrentHit;
	private TreeSet<SimpleSpliceScore>		m_vcSplicingScores;
	private TreeSet<CountElement>			m_vcHighlightedJunctions;
	private TreeSet<Exon>					m_vcHighlightedExons;
	
	private TreeMap<String, Color> m_mapColorsToConditions;
	
	@SuppressWarnings("unused")
	private class ThreadGetExonGroupCoverage implements Runnable
	{
		double					 						m_fSizeFactor;
		String 											m_strSample;
		int[]											m_pCoverage;
		TreeMap<String, TreeMap<ExonGroup, double[]>> 	m_mapCoveragePerExonGroupToSamples;
		
		ThreadGetExonGroupCoverage(int[] pCoverage, TreeMap<String, TreeMap<ExonGroup, double[]>> mapCoveragePerExonGroupToSamples, String strSample, double fSizeFactor, TreeMap<String, TreeMap<String, double[]>> mapCoverageToSamplesAndExons)
		{
			m_fSizeFactor 						= fSizeFactor;
			m_strSample							= strSample;
			m_mapCoveragePerExonGroupToSamples 	= mapCoveragePerExonGroupToSamples;
			m_pCoverage							= pCoverage;
		}
		
		@Override
		public void run()
		{
			TreeMap<ExonGroup, double[]> mapCoveragePerExonGroup = new TreeMap<ExonGroup, double[]>();

			for(ExonGroup grp : m_pExonGroups)
			{
				int nStart = grp.getGenomicStartOfGroup();
				int nEnd   = grp.getGenomicStopOfGroup();
				
				double pValues[] = new double[nEnd-nStart+1];
				
				for(int i=nStart; i<nEnd; i++)
				{
					int nPos = i - m_Gene.getStart();
					pValues[i-nStart] = m_pCoverage[nPos] * m_fSizeFactor;
				}
									
				mapCoveragePerExonGroup.put(grp, pValues);
			}
			
			synchronized(m_mapCoveragePerExonGroupToSamples)
			{
				m_mapCoveragePerExonGroupToSamples.put(m_strSample, mapCoveragePerExonGroup);
			}
		}
	}

	private class ThreadGetRawExonCoverage implements Runnable
	{
		Exon	m_Exon;
		String	m_strSample;
		boolean m_bAdjustForSizeFactor;
		TreeMap<Exon, TreeMap<String, Double>> m_mapExonCoveragePerSampleRawFromBigWig;
		TreeMap<Exon, TreeMap<String, Double>> m_mapExonCoverageLengthPerSampleRawFromBigWig;
		
		ThreadGetRawExonCoverage(Exon ex, String strSample, boolean bAdjustForSizeFactor,
				TreeMap<Exon, TreeMap<String, Double>> mapExonCoveragePerSampleRawFromBigWig,
				TreeMap<Exon, TreeMap<String, Double>> mapExonCoverageLengthPerSampleRawFromBigWig)
		{
			m_Exon = ex;
			m_strSample = strSample;
			m_bAdjustForSizeFactor = bAdjustForSizeFactor;
			m_mapExonCoveragePerSampleRawFromBigWig = mapExonCoveragePerSampleRawFromBigWig;
			m_mapExonCoverageLengthPerSampleRawFromBigWig = mapExonCoverageLengthPerSampleRawFromBigWig;
		}
		
		@Override
		public void run()
		{
			try
			{
				// get raw coverage for exon
				double pCov[] = GetCoverageForExon(m_Exon, m_strSample, m_bAdjustForSizeFactor);
				
				int nBasesBelowThreshold = 0;
				
				double fSum = 0.0;
				for(double fVal : pCov)
				{
					fSum += fVal;
					if(fVal < m_nMinCovPerBase)
						nBasesBelowThreshold += 1;
				}
				
				synchronized(m_mapExonCoveragePerSampleRawFromBigWig)
				{
					if(m_mapExonCoveragePerSampleRawFromBigWig.containsKey(m_Exon))
					{
						m_mapExonCoveragePerSampleRawFromBigWig.get(m_Exon).put(m_strSample, fSum / pCov.length);
					}
					else
					{
						TreeMap<String, Double> mapTmp =  new TreeMap<String, Double>();
						mapTmp.put(m_strSample, fSum / pCov.length);
						m_mapExonCoveragePerSampleRawFromBigWig.put(m_Exon, mapTmp);
					}
					
					if(m_mapExonCoverageLengthPerSampleRawFromBigWig.containsKey(m_Exon))
					{
						m_mapExonCoverageLengthPerSampleRawFromBigWig.get(m_Exon).put(m_strSample, 1.0-(double)(nBasesBelowThreshold) / pCov.length);
					}
					else
					{
						TreeMap<String, Double> mapTmp =  new TreeMap<String, Double>();
						mapTmp.put(m_strSample, 1.0-(double)(nBasesBelowThreshold) / pCov.length);
						m_mapExonCoverageLengthPerSampleRawFromBigWig.put(m_Exon, mapTmp);
					}
				}
			}
			catch (IOException e)
			{
				e.printStackTrace();
			}
		}
	}
	
	private class ThreadCalculateExonCoverage implements Runnable
	{
		ExonGroup m_ExonGroup;
		String m_strSample;
		
		TreeMap<Exon, TreeMap<String, Double>>			m_mapExonCoveragePerSample;
		TreeMap<Exon, TreeMap<String, TreeSet<Exon>>> 	m_mapIndistinguishableExonPerSample;

		ThreadCalculateExonCoverage(String strSample, ExonGroup exonGroup, TreeMap<Exon, TreeMap<String, Double>> mapExonCoveragePerSample, TreeMap<Exon, TreeMap<String, TreeSet<Exon>>> mapIndistinguishableExonPerSample)
		{
			m_strSample 						= strSample;
			m_ExonGroup 						= exonGroup;
			m_mapExonCoveragePerSample 			= mapExonCoveragePerSample;
			m_mapIndistinguishableExonPerSample	= mapIndistinguishableExonPerSample;
		}
		
		@Override
		public void run()
		{	
			// get a copy of the exons
			Vector<Exon> vcExons = new Vector<Exon>();
			for(Exon ex : m_ExonGroup.getExons())
				vcExons.add(ex);

			//###################################################################
			//    use only exons that have support by junction spanning reads
			//###################################################################
//			TreeSet<Exon> vcUnsupportedExons = RemoveUnsupportedExons(vcExons, m_mapJunctionReadCounts, m_strSample);
	
			//###################################
			//    get coverage for exon group
			//###################################
			double pCoverage[];
			try
			{
				pCoverage = GetCoverageForExonGroup(m_ExonGroup, m_strSample, true);
				
			}
			catch (IOException e)
			{
				e.printStackTrace();
				return;
			}
			if(pCoverage == null)
			{
				System.out.println("Skipping sample: " + m_strSample);
				return;
			}
/*
			synchronized(m_mapExonCoveragePerSample)
			{
				// add 0 coverage for unsupported exons
				for(Exon ex : vcUnsupportedExons)
				{
					if(m_mapExonCoveragePerSample.containsKey(ex))
					{
						m_mapExonCoveragePerSample.get(ex).put(m_strSample, 0.0);
					}
					else
					{
						TreeMap<String, Double> mapTmp = new TreeMap<String, Double>();
						mapTmp.put(m_strSample, 0.0);
						m_mapExonCoveragePerSample.put(ex, mapTmp);
					}
					
					// add as 'unique' exon to the mapIndistinguishableExonPerSample list
					TreeSet<Exon> vcTmpExons = new TreeSet<Exon>();
					vcTmpExons.add(ex);
					
					if(m_mapIndistinguishableExonPerSample.containsKey(ex))
					{
						TreeMap<String, TreeSet<Exon>> vcTmp = m_mapIndistinguishableExonPerSample.get(ex);							
						vcTmp.put(m_strSample, vcTmpExons);
					}
					else
					{
						TreeMap<String, TreeSet<Exon>> vcTmp = new TreeMap<String, TreeSet<Exon>>();
						vcTmp.put(m_strSample, vcTmpExons);
						m_mapIndistinguishableExonPerSample.put(ex, vcTmp);
					}
				}
			}
*/			
			//###############################
			//    group too similar exons
			//###############################
			TreeMap<Exon, TreeSet<Exon>> mapGroupedExons = GroupSimilarExons(vcExons, 50);
			
			//###############################################################
			//    Calculate the coverage per base for each exon and store
			//    them in a hashmap (key = coverage, value = exon)
			//###############################################################
			TreeSet<Exon> vcTmpExons = new TreeSet<Exon>();
			vcTmpExons.addAll(mapGroupedExons.keySet());
			HashMap<Exon, Double> mapCoverageToExon = CalculateExonCoverageGroupedByCoverage(vcTmpExons, m_ExonGroup.getGenomicStartOfGroup(), pCoverage);
			
			//###############################################################################
			//    Also merge exons with similar coverage that are completely overlapping.
			//###############################################################################
			MergeSimilarExpressedContainedExons(mapGroupedExons, mapCoverageToExon);

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
						vcTmp.put(m_strSample, mapGroupedExons.get(ex));
					}
					else
					{
						TreeMap<String, TreeSet<Exon>> vcTmp = new TreeMap<String, TreeSet<Exon>>();
						vcTmp.put(m_strSample, mapGroupedExons.get(ex));
						m_mapIndistinguishableExonPerSample.put(ex, vcTmp);
					}
				}
			}

			//###########################################################################
			//    Track which base positions are claimed by which exon, starting with
			//    the highest expressed exon.
			//###########################################################################
			TreeMap<Integer, TreeSet<Exon>> mapCoveragePositions = DefineExonSpecificCoverageRegions(mapCoverageToExon, m_ExonGroup.getGenomicStartOfGroup());

			//###############################################################################
			//    Recalculate exon coverages. Already used exon positions will be ignored
			//###############################################################################
			TreeMap<Exon, Double> mapUniqueCoverage = CalculateUniqueCoveragePerExon(mapCoveragePositions, mapGroupedExons, m_ExonGroup.getGenomicStartOfGroup(), pCoverage);
			
			//###############################################################################
			//    Now that we have defined 'unique' coverage regions per exon, we can use
			//    these regions to calculate the relative exon coverages.
			//###############################################################################
			TreeMap<Exon, Double> mapFinalExonCoverage = CalculateFinalCoveragePerExon(mapUniqueCoverage, m_ExonGroup.getGenomicStartOfGroup(), pCoverage);
			
			synchronized(m_mapExonCoveragePerSample)
			{
				for(Exon ex : mapFinalExonCoverage.keySet())
				{
					if(m_mapExonCoveragePerSample.containsKey(ex))
					{
						TreeMap<String, Double> mapCoveragePerSample = m_mapExonCoveragePerSample.get(ex);
						mapCoveragePerSample.put(m_strSample, mapFinalExonCoverage.get(ex));
					}
					else
					{
						TreeMap<String, Double> mapCoveragePerSample = new TreeMap<String, Double>();
						mapCoveragePerSample.put(m_strSample, mapFinalExonCoverage.get(ex));
						m_mapExonCoveragePerSample.put(ex, mapCoveragePerSample);
					}
				}
			}
		}
	}
	
	public class ThreadRemoveIrrelevantIsoformsBasedOnSplitReads implements Runnable
	{
		String 										m_strIsoform;
		String										m_strCondition;
		TreeSet<String> 							m_vcIrrelevantIsoforms;
		TreeMap<String, TreeSet<String>> 			m_vcSamplesPerCondition;
		TreeMap<String, TreeMap<String, Integer>> 	m_mapJunctionReadCounts;
		boolean										m_bDebug;
		boolean										m_bSkipFirstAndLastExon;
		
		ThreadRemoveIrrelevantIsoformsBasedOnSplitReads(String strIsoform, String strCondition, TreeSet<String> vcIrrelevantIsoforms, TreeMap<String, TreeSet<String>> vcSamplesPerCondition, TreeMap<String, TreeMap<String, Integer>> mapJunctionReadCounts, boolean bSkipFirstAndLastExon, boolean bDebug)
		{
			m_strIsoform 			= strIsoform;
			m_strCondition			= strCondition;
			m_vcIrrelevantIsoforms	= vcIrrelevantIsoforms;
			m_vcSamplesPerCondition = vcSamplesPerCondition;
			m_mapJunctionReadCounts = mapJunctionReadCounts;
			m_bDebug				= bDebug;
			m_bSkipFirstAndLastExon = bSkipFirstAndLastExon; 
		}
		
		@Override
		public void run()
		{			
			int nSamples = 0;
			for(String strSample : m_vcSamplesPerCondition.get(m_strCondition))
			{
				if(m_vcSelectedSamples.contains(strSample))
					nSamples++;
			}
			
			//##########################################################################################
			//          check whether all junctions of the isoform have sufficient coverage
			//##########################################################################################
			TreeSet<String> vcJunctions = m_Gene.getSpliceJunctionInformationForGeneProduct(m_strIsoform);
			int nJunIdx = 0;
			for(String strJunction : vcJunctions)
			{
				if(m_bSkipFirstAndLastExon)
				{
					// skip first and last junction
					if(nJunIdx == 0 || nJunIdx == vcJunctions.size()-1)
						continue;
				}
				
				String pSplit[] = strJunction.split("-");
				String strJunNew = pSplit[0] + "_" + pSplit[1];
				double pCounts[] = new double[nSamples];
				
				int nIdx = 0;
				for(String strSample : m_vcSamplesPerCondition.get(m_strCondition))
				{
					if(m_vcSelectedSamples.contains(strSample))
					{
						if(m_mapJunctionReadCounts.containsKey(strJunNew))
						{
							if(!m_mapJunctionReadCounts.get(strJunNew).containsKey(strSample))
							{
								System.out.println("no junction counts found for sample: " + strSample + " (" + strJunNew + ")");
								continue;
							}
							pCounts[nIdx] = m_mapJunctionReadCounts.get(strJunNew).get(strSample);
						}
						else
						{
							System.out.println("could not detect any counts for the following junction: " + strJunNew);
						}
						nIdx++;
					}
				}
				
				double fMedian = StatUtils.percentile(pCounts, 50.0);
		
				if(fMedian < m_nMinJunctionReads)
				{
					synchronized(m_vcIrrelevantIsoforms)
					{
						if(m_bDebug)
						{
							System.out.println("[min junctions reads] -> removing " + m_strIsoform);
						}
						m_vcIrrelevantIsoforms.add(m_strIsoform.split("\\.")[0]);
						return;
					}
				}
				
				nJunIdx++;
			}
		}
	}

	public class ThreadRemoveIrrelevantIsoformsBasedOnCoverage implements Runnable
	{
		String 										m_strIsoform;
		String										m_strCondition;
		TreeSet<String> 							m_vcIrrelevantIsoforms;
		TreeMap<String, TreeSet<String>> 			m_vcSamplesPerCondition;
		TreeMap<String, TreeMap<String, Integer>> 	m_mapJunctionReadCounts;
		boolean										m_bDebug;
		boolean										m_bSkipFirstAndLastExon;
		
		ThreadRemoveIrrelevantIsoformsBasedOnCoverage(String strIsoform, String strCondition, TreeSet<String> vcIrrelevantIsoforms, TreeMap<String, TreeSet<String>> vcSamplesPerCondition, TreeMap<String, TreeMap<String, Integer>> mapJunctionReadCounts, boolean bSkipFirstAndLastExon, boolean bDebug)
		{
			m_strIsoform 			= strIsoform;
			m_strCondition			= strCondition;
			m_vcIrrelevantIsoforms	= vcIrrelevantIsoforms;
			m_vcSamplesPerCondition = vcSamplesPerCondition;
			m_mapJunctionReadCounts = mapJunctionReadCounts;
			m_bDebug 				= bDebug;
			m_bSkipFirstAndLastExon	= bSkipFirstAndLastExon;
		}
		
		@Override
		public void run()
		{
			// skip already irrelevant isoforms
			if(m_vcIrrelevantIsoforms.contains(m_strIsoform.split("\\.")[0]))
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
			Exon pExons[] = m_Gene.getSortedExonsForGeneProduct(m_strIsoform, true);
		
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
						if(!m_mapExonCoveragePerSampleRawFromBigWig.containsKey(ex))
						{
							System.out.println("failed to retrieve raw exon coverage for exon [" + ex + "] and sample: " + strSample);
							System.out.println(m_mapExonCoveragePerSampleRawFromBigWig);
							continue;
						}
						
						// combine coverage of all exons within the exon group for first and last exons
						if(bIsFirstExon || bIsLastExon)
						{
							ExonGroup[] grps = m_Gene.computeOverlappingExonGroups();
							for(ExonGroup grp : grps)
							{
								if(grp.groupContainsExon(ex.getExonID()))
								{
									double fCoverage = 0.0f;
									double fCoverageLength = 0.0;
									
									for(Exon exon : grp.getExons())
									{
										// Only merge first with other first exons and last with other last exons
										Iterator<String> it = m_Gene.getGeneProductNames();
										boolean bOkay = false;
										while(it.hasNext())
										{
											String strIsoform = it.next();
											if(!m_vcIrrelevantIsoforms.contains(strIsoform.split("\\.")[0]))
											{
												Exon pIsoformExons[] = m_Gene.getSortedExonsForGeneProduct(strIsoform);
												if((pIsoformExons[0] == exon && bIsFirstExon) || (pIsoformExons[pIsoformExons.length-1] == exon && bIsLastExon))
												{
													// exons must also overlap
													if(exon.getCodingStart() <= ex.getCodingStop() && exon.getCodingStop() >= ex.getCodingStart())
													{
														/*//T_ODO
														if(m_strIsoform.contains("525604"))
															System.out.println("exon valid: " + strIsoform.split("\\.")[0]);
														*/
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
											/*
											//T_ODO
											if(m_strIsoform.contains("525604"))
												System.out.println("adding cov for: " + exon + m_mapExonCoveragePerSampleRawFromBigWig.get(exon).get(strSample));
											 */
			
											if(m_mapExonCoveragePerSampleRawFromBigWig == null)
												System.out.println("ERROR: m_mapExonCoveragePerSampleRawFromBigWig is: " + m_mapExonCoveragePerSampleRawFromBigWig);
											
											if(!m_mapExonCoveragePerSampleRawFromBigWig.containsKey(exon))
											{
												System.out.println("ERROR: m_mapExonCoveragePerSampleRawFromBigWig is missing: " + exon);
												System.out.println(m_mapExonCoveragePerSampleRawFromBigWig);
											}
											
											if(!m_mapExonCoveragePerSampleRawFromBigWig.get(exon).containsKey(strSample))
											{
												System.out.println("ERROR: m_mapExonCoveragePerSampleRawFromBigWig is missing: " + strSample);
												System.out.println(m_mapExonCoveragePerSampleRawFromBigWig.get(exon));
											}
											
											fCoverage		+= m_mapExonCoveragePerSampleRawFromBigWig.get(exon).get(strSample);
											fCoverageLength += m_mapExonCoverageLengthPerSampleRawFromBigWig.get(exon).get(strSample);
										}
									}

									/*
									//T_ODO
									if(m_strIsoform.contains("525604"))
									{
										System.out.println(m_strIsoform);
										System.out.println("coverage: " + fCoverage);
										System.out.println("length: " + fCoverageLength);
										System.out.println(ex);
										System.out.println(grp);
										System.out.println("-------------");
									}
									*/

									pCoverages[nSample] 	= fCoverage;
									pCoveredBases[nSample]	= fCoverageLength;
									
									break;
								}
							}
						}
						else
						{
							pCoverages[nSample] 	= m_mapExonCoveragePerSampleRawFromBigWig.get(ex).get(strSample);
							pCoveredBases[nSample]	= m_mapExonCoverageLengthPerSampleRawFromBigWig.get(ex).get(strSample);
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
						m_vcIrrelevantIsoforms.add(m_strIsoform.split("\\.")[0]);
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

	private class DrawingOffsets
	{
		int m_nBorderOffset;
		int m_nIsoformNameOffset;
		int m_nMaxExonLength;
		int m_nMMSeqTextLength;
		int m_nIntronLength;
		int m_nTotalWidth;
	};
	
	//##########################################################################
	// https://eleanormaclure.files.wordpress.com/2011/03/colour-coding.pdf
	// Kelly's colors of maximum contrast
	//##########################################################################
	private Color m_pColors[] = {
	    new Color(0xFFFFB300, true), // Vivid Yellow								1
	    new Color(0xFF803E75, true), // Strong Purple								2
	    new Color(0xFFFF6800, true), // Vivid Orange								3
	    new Color(0xFFA6BDD7, true), // Very Light Blue								4
	    new Color(0xFFC10020, true), // Vivid Red									5
	    new Color(0xFFCEA262, true), // Grayish Yellow								6
	    new Color(0xFF817066, true), // Medium Gray									7

	    // The following will not be good for people with defective color vision	
	    new Color(0xFF007D34, true), // Vivid Green									8
	    new Color(0xFFF6768E, true), // Strong Purplish Pink						9
	    new Color(0xFF00538A, true), // Strong Blue									10
	    new Color(0xFFFF7A5C, true), // Strong Yellowish Pink						11
	    new Color(0xFF53377A, true), // Strong Violet								12
	    new Color(0xFFFF8E00, true), // Vivid Orange Yellow							13
	    new Color(0xFFB32851, true), // Strong Purplish Red							14
	    new Color(0xFFF4C800, true), // Vivid Greenish Yellow						15
	    new Color(0xFF7F180D, true), // Strong Reddish Brown						16
	    new Color(0xFF93AA00, true), // Vivid Yellowish Green						17
	    new Color(0xFF593315, true), // Deep Yellowish Brown						18
	    new Color(0xFFF13A13, true), // Vivid Reddish Orange						19
	    new Color(0xFF232C16, true), // Dark Olive Green							20
	    
	    // And some more added by myself (Matthias Barann)
	    new Color(  0, 255, 153, 127), // Teal-Teal-Cyan							21
	    new Color(  0, 255,  51, 127),  // Green-Green-Teal							22
	    new Color(  0, 255,  51, 127),  // Green-Green-Teal							23
	    new Color(  0,   0, 153, 127),  // Dark Faded Blue							24
	    new Color(  0, 204, 255, 127),  // Cyan-Cyan-Azure							25
	    new Color(255,   0, 153, 127),  // Pink-Pink-Magenta						26
	    new Color(204, 204,  51, 127),  // Medium Faded Yellow						27
	    new Color(204,  51,  51, 127),  // Medium Faded Red							28
	    new Color(255, 153,   0, 127),  // Orange-Orange-Yellow						29
	    new Color(204,   0, 255, 127)   // Magenta-Magenta-Violet					30
	};
	
	private static final long serialVersionUID = 1L;
	
	private class ClickEvent
	{
		String m_strClickedA;	// name of first clicked element
		String m_strClickedB;	// name of second clicked element
		
		TreeSet<Exon> m_vcSelectedExons;
		TreeSet<String> m_vcTargetIsoforms;
		
		int m_nExonGroupStartA;
		int m_nExonGroupStartB;
		int m_nExonGroupEndA;
		int m_nExonGroupEndB;
		
		ClickEvent()
		{
			m_strClickedA = null;
			m_strClickedB = null;

			m_nExonGroupStartA	= -1;
			m_nExonGroupStartB	= -1;
			m_nExonGroupEndA	= -1;
			m_nExonGroupEndB	= -1;
			
			m_vcSelectedExons = new TreeSet<Exon>();
			m_vcTargetIsoforms = new TreeSet<String>();
		}
		
		void clear()
		{
			m_strClickedA = null;
			m_strClickedB = null;

			m_nExonGroupStartA	= -1;
			m_nExonGroupStartB	= -1;
			m_nExonGroupEndA	= -1;
			m_nExonGroupEndB	= -1;
			
			m_vcSelectedExons.clear();
			m_vcTargetIsoforms.clear();
		}
	}

	public SplicingWebApp(boolean bEmpty) throws IOException
	{	
		ReadAppPaths();
				
		m_nMaxWidth			= 1800;
		
		m_Gene				= null;
		m_pExons 			= null;
		m_pExonGroups		= null;
		
		m_mapBigWigCoverageToSamples = null;

		m_vcValidIsoforms							= new TreeSet<String>();
		m_mapExonsToIsoforms 						= new TreeMap<String, int[]>();
		
		m_mapExonCoveragePerSample					= null;
		m_mapExonCoveragePerSampleRawFromBigWig		= null;		
		m_mapIndistinguishableExonPerSample 		= null;
		m_mapExonCoverageLengthPerSampleRawFromBigWig = null;
		
		m_mapAnalysedJunctionPaths	 = new TreeMap<String, TreeMap<Integer, TreeSet<String>>>();
		m_mapCoveragePerJunctionPathAndSample = new TreeMap<String, TreeMap<Integer, Double>>();
		m_mapEntropyToExonicPart	= null;
		m_mapCountsPerExonAndTissue	= null;
		
		m_nOrganismType		= 0;
		
		m_ClickEvent		= new ClickEvent();
		m_projectModel 		= new ProjectModel();
		
		m_vcSelectedSamples = new TreeSet<String>();
		m_treeSelectedIsoforms = new Tree();
		m_bandboxSelectedGene = new Bandbox();
		
		m_imgMapIsoforms	= null;
		m_imgMapCoverage	= null;
		
		m_layoutNorth		= new North();
		m_layoutWest		= new West();
		m_layoutCenter		= new Center();
				
		m_layoutMainTop	 	= new North();
		m_layoutMainSouth	= new South();	
		m_layoutMainPlots 	= new Borderlayout();	
		m_layoutMainBottom 	= new Vlayout();
		
		m_vcSplicingScores = new TreeSet<SimpleSpliceScore>();
	}

	public SplicingWebApp() throws Exception
	{
		BasicConfigurator.configure();
		ReadAppPaths();
		
		m_hWindow = this;
		
		m_nOrganismType		= 0;
		
		m_Gene				= null;
		m_pExons 			= null;
		m_pExonGroups		= null;
		
		m_mapBigWigCoverageToSamples = null;

		m_imgMapIsoforms	= null;
		m_imgMapCoverage	= null;
		m_imgMapColors		= null;
		
		m_layoutNorth		= new North();
		m_layoutWest		= new West();
		m_layoutCenter		= new Center();
				
		m_layoutMainTop	 	= new North();
		m_layoutMainSouth	= new South();		

		m_bShowAmbigiousUniqueExons			= false;
		m_bShowUniqueFeatures 				= false;
		m_bCoverageRequiresRedraw 			= true;
		m_bColorExonsByRelativeCoverage 	= false;
		m_bShowRelativeCoverageNumbersExons	= false;
		m_bUseRelativeExonCoverageMode1 	= false;
		m_bUseRelativeExonCoverageMode2		= false;
		m_bColorJunctionPaths				= false;
		m_bShowRelativeCoverageNumbersJunctions = false;
		m_bShowEntropyData					= false;
		m_bShowASEvents						= true;
		m_bShowSecondCoveragePlot			= false;
		m_bColorExonsAndJunctionsByCoverage	= true;
		
		m_vcValidIsoforms					= new TreeSet<String>();
		m_mapExonsToIsoforms 				= new TreeMap<String, int[]>();
		m_bIsoformSelectionChanged			= false;
		m_bHideIrrelevantIsoforms			= false;
		m_bSkipFirstAndLastExon				= true;

		m_nMinJunctionReads					= 1;
		m_nMinCovPerBase					= 3;
		m_fMinCoveredBases					= 0.7;
		m_fVariableExonThreshold			= 0.15;		
		
		m_windowPopup							= null;
		m_bCoveragePlotUseLog2					= false;
		m_bCoveragePlotShowRelativeCoverage 	= false;
		m_bCoveragePlotShowGeometricMean		= false;
		m_bCoveragePlotShowMedian				= true;
		m_bCoveragePlotShowMean					= false;
		m_bCoveragePlotShowQuartiles			= true;
		
		m_mapExonCoveragePerSample					= null;
		m_mapExonCoveragePerSampleRawFromBigWig		= null;		
		m_mapIndistinguishableExonPerSample 		= null;
		m_mapExonCoverageLengthPerSampleRawFromBigWig = null;
		
		m_mapAnalysedJunctionPaths	 = new TreeMap<String, TreeMap<Integer, TreeSet<String>>>();
		m_mapCoveragePerJunctionPathAndSample = new TreeMap<String, TreeMap<Integer, Double>>();
		m_mapEntropyToExonicPart	= null;
		m_mapCountsPerExonAndTissue	= null;
		
		m_ClickEvent		= new ClickEvent();
		
		m_projectModel 		= new ProjectModel();
		
		m_vcPreviouslySelectedGTEXTissues	= null;
		m_vcPreviouslySelectedTreeItems		= new TreeSet<String>();
		
		m_bUseReducedDataSet		= true;
		m_bSampleSelectionChanged	= false;
		m_treeSelectedSamples		= null;
		m_vcSelectedSamples			= new TreeSet<String>();
		
		m_strSelectedCondition		= null;
		m_strSelectedConditionType 	= null;
		
		m_vcASHits					= new TreeSet<AlternativeSplicingHit>();
		m_vcASCurrentHit			= null; 
		m_vcSplicingScores			= new TreeSet<SimpleSpliceScore>();
		m_vcHighlightedJunctions	= new TreeSet<CountElement>();
		m_vcHighlightedExons		= new TreeSet<Exon>();
		
		m_vcGeneIdentifier = new Vector<GeneIdentifier>();
		
		// read human cross reference file by default
		String strFileEntrezGene = m_strPathReferences + "/biomart_cross_ref_human.txt";
		ReadGeneIdentifier(strFileEntrezGene);
		
		addEventListener(Events.ON_CLIENT_INFO, new EventListener<ClientInfoEvent>()
		{
			public void onEvent(ClientInfoEvent event) throws Exception
			{
				m_nMaxWidth = event.getDesktopWidth();
				
				m_nClientWindowWidth  = event.getDesktopWidth();
				m_nClientWindowHeight = event.getDesktopHeight();
			}
		});
		
		//######################################
		//           add gui elements
		//######################################
		m_layout = new Borderlayout();
		m_layout.setVflex("min");
		m_layout.setHeight("100%");
		m_layout.setWidth("100%");
		m_layout.setParent(this);

		m_layoutNorth.setHflex("min");
		m_layoutNorth.setParent(m_layout);
		m_layoutNorth.setMargins("0,5,5,5");		
	
		//#################################################################################################
		//     add main plot layouts (top=coverage track, center = isoforms, bottom = novel junctions)
		//#################################################################################################
		m_layoutCenter.setVflex("min");
		m_layoutCenter.setParent(m_layout);

		m_layoutMainPlots 		 = new Borderlayout();
		m_layoutMainPlots.setVflex("min");
		m_layoutMainPlots.setParent(m_layoutCenter);		

		m_layoutMainTop.setTitle("Coverage Track");
		m_layoutMainTop.setHeight("280px");
		m_layoutMainTop.setAutoscroll(true);
		m_layoutMainTop.setAttribute("org.zkoss.zul.nativebar", "true");
		m_layoutMainTop.setParent(m_layoutMainPlots);

		m_layoutMainSouth.setTitle("Isoform View");
		m_layoutMainSouth.setHeight("600px");
		m_layoutMainSouth.setParent(m_layoutMainPlots);
		
		m_layoutMainBottom 	= new Vlayout();
		m_layoutMainBottom.setParent(m_layoutMainSouth);
		
		m_layoutMainBottom.setAttribute("org.zkoss.zul.nativebar", "true");
		m_layoutMainBottom.setStyle("overflow: auto;");

		//#################################
		//      add buttons/options
		//#################################
		Hlayout layoutH = new Hlayout();
		layoutH.setParent(m_layoutNorth);
		layoutH.setStyle("overflow:auto; margin-top: 10px; margin-right: 10px; margin-bottom: 10px;");
		
		Vlayout layoutV = new Vlayout();
		layoutV.setHeight("200px");
		layoutV.setVflex("min");
		layoutV.setParent(layoutH);

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
		tabbox.setHeight("470px");
		tabbox.setParent(layoutH);
		
		Tabs tabs = new Tabs();
		tabs.setParent(tabbox);
		
		Tabpanels panels = new Tabpanels();
		panels.setParent(tabbox);
		
		AddHitList(tabs, panels);
		AddPSIHitList(tabs, panels);
		
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
	}

	public void ReadAppPaths() throws IOException
	{
		String strFileName = System.getProperty("user.home") + "/app_paths.config"; 
		File pIn = new File(strFileName);
		
		// use 1 thread by default
		m_nThreads = 1;
		m_strPathGeneGTEX		= null;
		m_strPathExonGTEX		= null;
		m_strPathMMSeqResults	= null;
		m_strMMSeqExecutable	= null;

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
					case "path_GTEX_gene_counts": 	m_strPathGeneGTEX 	= pSplit[1].trim().replace("\"", ""); break;
					case "path_GTEX_exon_counts": 	m_strPathExonGTEX 	= pSplit[1].trim().replace("\"", ""); break;
					case "path_tmp": 				m_strTmpFolder 		= pSplit[1].trim().replace("\"", ""); break;
					case "path_MMSEQ_results": 		m_strPathMMSeqResults = pSplit[1].trim().replace("\"", ""); break;
					case "num_threads":				m_nThreads 			= Integer.parseInt(pSplit[1].trim());
					case "path_hits":				m_strPathHitLists	= pSplit[1].trim().replace("\"", ""); break;
					case "path_PSI_hits":			m_strPathSplicingHitList	= pSplit[1].trim().replace("\"", ""); break;
					
					case "path_MMSEQ_executable":	m_strMMSeqExecutable = "/home/softadm/ssdadm/cellodbtools/tools/mmseq/mmseq-1.0.8-Linux";
				}
				
				System.out.println("input directory: " + m_strPathInput);
			}
			
			pScanner.close();
		}
		else
		{
			System.out.println("failed to read: " + strFileName);
		}
	}
	
	public void ProcessGeneInformation(GeneIdentifier gid, String strFileGTF) throws Exception
	{
		m_Gene	 = null;
		m_pExons = null;
		
		if(strFileGTF == null)
		{
			Messagebox.show("No gene annotation database selected!");
			return;
		}
	
		// check whether it is a GTF/GFF file or refFlat file
		File pFile = new File(strFileGTF);
		if(m_nReferenceType == GTF_REFERENCE_FILE)
		{
			// check if the index file exists
			File pFileIndex = new File(strFileGTF + ".idx");
	
			String strGene = gid.m_strEnsemblGeneID;
			try
			{
				RandomAccessGFFReader gffReader = new RandomAccessGFFReader(new File(strFileGTF), pFileIndex);
				m_Gene = gffReader.ReadGene(strGene);
			}
			catch(Exception e)
			{
				System.out.println("failed to open GTF file: " + strFileGTF);
				System.out.println(e.getMessage());
			}
			
			// try again with m_strEntrezGeneID
			if(m_Gene == null)
			{
				strGene = gid.m_strEntrezGeneID;
				try
				{
					RandomAccessGFFReader gffReader = new RandomAccessGFFReader(new File(strFileGTF), pFileIndex);
					m_Gene = gffReader.ReadGene(strGene);
				}
				catch(Exception e)
				{
					System.out.println("failed to open GTF file: " + strFileGTF);
					System.out.println(e.getMessage());
				}
			}
		}
		else if(m_nReferenceType == REFFLAT_REFERENCE_FILE)
		{
			Scanner pIn = new Scanner(pFile);
			
			String strGene = gid.m_strRefGeneID;
			
			while(pIn.hasNextLine())
			{
				String strLine = pIn.nextLine();
				
				if(strLine.startsWith("#"))
					continue;
				
				String pSplit[] = strLine.split("\\s+");
				
				String strGeneID = pSplit[1];
				String strRef = pSplit[2];
				int nStart  = Integer.parseInt(pSplit[4])+1; // must be made 1-based
				int nEnd	= Integer.parseInt(pSplit[5]);  
				String strGeneSymbol = pSplit[12];
				
				if(strGeneID.equals(strGene) || (m_Gene != null && m_Gene.getGeneName().equals(strGeneSymbol) && m_Gene.getStart() <= nEnd && m_Gene.getStop() >= nStart && m_Gene.getChromosome().equals(strRef)))
				{
					boolean bFirstStrand = false;
					if(pSplit[3].equals("+"))
						bFirstStrand = true;
					
					Gene trans = new Gene(strGeneID, strGeneSymbol, "?", nStart, nEnd, strRef, bFirstStrand);
					trans.ParseFromRefFlat(pSplit);

					if(m_Gene == null)
						m_Gene = trans;
					
					// add new isoform to current gene
					m_Gene.addGeneProduct(strGeneID,  "", trans.getSortedExons());
				}
			}
			
			pIn.close();
		}
		else
		{
			Messagebox.show("Invalid reference file detected");
		}
		
		//############################
		//    get gene information
		//############################
		if(m_Gene == null)
		{
			Messagebox.show("Could not find " + gid.m_strApprovedGeneSymbol + " (" + gid.m_strEnsemblGeneID  + ") in GTF file.");
			return;
		}

		Exon pExonsTmp[] = m_Gene.getSortedExons();
		
		// reverse exon order
		m_pExons = new Exon[pExonsTmp.length];
		for(int i=0; i<pExonsTmp.length; i++)
		{
			m_pExons[i] = pExonsTmp[pExonsTmp.length-i-1];
		}
		
		// get exon groups
		m_pExonGroups = null;
		m_pExonGroups = m_Gene.computeOverlappingExonGroupsNormalOrder();
		
		// get bigwig coverages
		GetBigWigCoverageForGene();
		
		// clear old gene data
		if(m_vcValidIsoforms == null)
			m_vcValidIsoforms = new TreeSet<String>();
		else
			m_vcValidIsoforms.clear();
		
		if(m_mapExonsToIsoforms == null)
			m_mapExonsToIsoforms = new TreeMap<String, int[]>();
		else
			m_mapExonsToIsoforms.clear();
		
		//#####################################
		//    get all isoforms for the gene
		//#####################################
		Iterator<String> it = m_Gene.getGeneProductNames();
		while(it.hasNext())
		{
			String strGeneProductID = it.next();
			int[] pExonIDs = m_Gene.getExonsForGeneProduct(strGeneProductID);
			
			m_mapExonsToIsoforms.put(strGeneProductID, pExonIDs);
			
			// add all isoforms on gene change
			m_vcValidIsoforms.add(strGeneProductID.split("\\.")[0]);
		}
		
		UpdateComboboxForIsoformSelection(true);
	}

	ExonGroup[] RecalculateExonGroups(TreeMap<String, int[]> mapValidIsoforms)
	{
		TreeSet<Exon> vcAllExons = new TreeSet<Exon>();
		for(Map.Entry<String, int[]> e : mapValidIsoforms.entrySet())
		{
			for(int nExonID : e.getValue())
			{
				// convert exon ID, because the exons are sorted in the reverse order
				int nID = Math.abs(nExonID - (m_pExons.length-1));
				m_pExons[nID].setID(nID);
				vcAllExons.add(m_pExons[nID]);
			}
		}
		
		int nGrpID = 0;
		TreeSet<ExonGroup> vcExonGroups = new TreeSet<ExonGroup>();
		for(Exon exon : vcAllExons)
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

	private DrawingOffsets ComputeDrawingOffsets()
	{
		DrawingOffsets res = new DrawingOffsets();
		
		TreeMap<String, TreeSet<String>> vcSamplesAndConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		// general 'fixed' offsets
		int nOffsetX = 20;
		int nIntronLength = 20;
		
		int nMaxWidth	= m_nMaxWidth;
		// the height is not important right now
		int nMaxHeight  = 200;
		
		// prepare graph
		BufferedImage img = new BufferedImage(nMaxWidth, nMaxHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D graph = img.createGraphics();
		
		// set font size
		Font fontBoldText	= new Font("Courier", Font.BOLD, 13);
		graph.setFont(fontBoldText);
		
		// define selected conditions
		TreeSet<String> vcSelectedConditions = new TreeSet<String>();
		for(String strCondition : vcSamplesAndConditions.keySet())
		{
			for(String strSample : m_vcSelectedSamples)
			{
				if(vcSamplesAndConditions.get(strCondition).contains(strSample))
					vcSelectedConditions.add(strCondition);
			}
		}
		
		// get longest condition name and adjust isoform offset if necessary
		int nIsoformNameOffset 	= 100; // save 100 px for the isoform and condition names (or more if condition names are longer)
		int nMMSeqTextLength 	= 0; // this size is added to the graph for MMseq results
		for(String strCondition : vcSelectedConditions)
		{
			FontMetrics fontMetrics = graph.getFontMetrics(fontBoldText);
			int nStringWidth = fontMetrics.stringWidth(strCondition);
			
			nIsoformNameOffset = Math.max(nIsoformNameOffset, nStringWidth+20);
			nMMSeqTextLength += nStringWidth + 20;
		}
		
		// calculate exon group positions (maximum screen size - transcript lengths - 2xborder - isoform name length)
		int nMaxExonLength = nMaxWidth - (nIntronLength*m_pExonGroups.length) - (nOffsetX*2) - nIsoformNameOffset;
		
		res.m_nBorderOffset 		= nOffsetX;
		res.m_nIntronLength			= nIntronLength;
		res.m_nIsoformNameOffset	= nIsoformNameOffset;
		res.m_nMaxExonLength		= nMaxExonLength;
		res.m_nMMSeqTextLength		= nMMSeqTextLength;
		res.m_nTotalWidth			= nMaxWidth + nMMSeqTextLength;
		
		return res;
	}
	
	private void DrawExtendedCoveragePlot() throws IOException
	{
		TreeMap<String, Double> mapSizeFactors	= m_projectModel.GetSizeFactors();
		TreeMap<String, String> mapBigWigFiles 	= m_projectModel.GetBigWigFilesForSamples();
		TreeSet<String> vcSamples 				= m_projectModel.GetSamples();
		
		if(mapBigWigFiles == null)
		{
			System.out.println("ERROR: no valid bigwig or bam files detected");
			return;
		}
		
		// get reference name
		String strRef = m_Gene.getChromosome();
		
		// get samples per condition
		final TreeMap<String, TreeSet<String>> vcSamplesAndConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		//########################################
		//      prepare coverage containers
		//########################################
		// key 1 = exon_name, key 2 = sample, value = coverage array
		TreeMap<String, TreeMap<String, double[]>> mapCoverageToSamplesAndExons = new TreeMap<String, TreeMap<String, double[]>>();
		for(ExonGroup grp : m_pExonGroups)
		{			
			int nStart	= grp.getGenomicStartOfGroup();
			int nEnd	= grp.getGenomicStopOfGroup();
			int nLength = nEnd-nStart+1;
			
			String strName = strRef + ":" + nStart + "-" + nEnd;
			
			mapCoverageToSamplesAndExons.put(strName, new TreeMap<String, double[]>());
			
			for(String strSample: vcSamples)
			{
				// use only selected samples				
				if(m_vcSelectedSamples.contains(strSample))
				{
					TreeMap<String, double[]> tmp = mapCoverageToSamplesAndExons.get(strName);
					tmp.put(strSample, new double[nLength]);
				}
			}
		}
		
		//################################################
		//         define selected conditions
		//################################################
		TreeSet<String> vcSelectedConditions = new TreeSet<String>();
		for(String strCondition : vcSamplesAndConditions.keySet())
		{
			for(String strSample : m_vcSelectedSamples)
			{
				if(vcSamplesAndConditions.get(strCondition).contains(strSample))
					vcSelectedConditions.add(strCondition);
			}
		}
		
		//######################################################
		//                 prepare graph
		//######################################################
		DrawingOffsets offsets = ComputeDrawingOffsets();
		
		int nMaxHeight  = 220;	// 200 for the track, 20 as border
		
		// prepare graph
		BufferedImage img = new BufferedImage(offsets.m_nTotalWidth, nMaxHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D graph = img.createGraphics();
		
		// set font size
		Font fontBoldText	= new Font("Courier", Font.BOLD, 13);
		Font fontNormalText = new Font("Courier", Font.PLAIN, 11);
		
		graph.setFont(fontBoldText);
		
		// fill background
		graph.setColor(Color.white);
		graph.fillRect(0, 0, offsets.m_nTotalWidth, nMaxHeight);
		
		// reserve some space for the isoform name
		int nOffsetX = offsets.m_nBorderOffset + offsets.m_nIsoformNameOffset;
		
		int nCombinedExonLength = 0;
		for(ExonGroup grp : m_pExonGroups)
		{
			nCombinedExonLength += grp.getExonGroupLengthInBp();
		}
		
		// draw condition name in associated color		
		int nOffset = 0;
		for(String strCondition : vcSamplesAndConditions.keySet())
		{
			// check that condition is selected
			if(!vcSelectedConditions.contains(strCondition))
				continue;
			
			graph.setColor(m_mapColorsToConditions.get(strCondition));
			graph.drawString(strCondition, 5, 12 + nOffset);
			nOffset += 13;
		}

		// calculate shrinkage factor for window width
		double fXShrinkageFactor = (double)offsets.m_nMaxExonLength / (double)nCombinedExonLength;
				
		//####################################################################
		//     obtain coverage for all exons and estimate maximum coverage
		//####################################################################
		double fMaxCoverage = 1.0;
		
		for(String strCondition : vcSamplesAndConditions.keySet())
		{
			for(String strSample : vcSamplesAndConditions.get(strCondition))
			{
				if(!m_vcSelectedSamples.contains(strSample))
					continue;

				// get size factor for current gene
				if(mapSizeFactors == null || !mapSizeFactors.containsKey(strSample))
				{
					System.out.println("warning: no size factor available for sample: " + strSample);
					continue;
				}
	
				for(ExonGroup grp : m_pExonGroups)
				{
					int nStart	= grp.getGenomicStartOfGroup();
					int nEnd	= grp.getGenomicStopOfGroup();
					String strName = strRef + ":" + nStart + "-" + nEnd;
					
					double[] pCoverage = GetCoverageForExonGroup(grp, strSample, true);
					
					if(m_bCoveragePlotUseLog2)
					{
						for(int i=0; i<pCoverage.length; i++)
						{
							pCoverage[i] = Math.log(pCoverage[i]+1) / Math.log(2);
						}
					}
					
					mapCoverageToSamplesAndExons.get(strName).put(strSample, pCoverage);
					
					if(!m_bCoveragePlotShowMedian && !m_bCoveragePlotShowGeometricMean && !m_bCoveragePlotShowMean && m_vcSelectedSamples.contains(strSample))
						fMaxCoverage = Math.max(fMaxCoverage, StatUtils.max(pCoverage)); 
				}
			}
		}
				
		if(m_bCoveragePlotShowMedian || m_bCoveragePlotShowGeometricMean || m_bCoveragePlotShowMean)
		{
			for(ExonGroup grp : m_pExonGroups)
			{
				// get correct group coverage data
				String strName = strRef + ":" + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup();
				TreeMap<String, double[]> mapData = mapCoverageToSamplesAndExons.get(strName);
	
				// combine coverage per condition
				for(String strCondition : vcSamplesAndConditions.keySet())
				{
					graph.setColor(this.m_mapColorsToConditions.get(strCondition));
					
					int nSamples = 0;
					for(String strSample : vcSamplesAndConditions.get(strCondition))
					{
						if(m_vcSelectedSamples.contains(strSample))
							nSamples += 1;
					}
	
					double[][] pValues =  new double[grp.getExonGroupLengthInBp()+1][nSamples];
					
					int nCurrentSample = 0;
					for(String strSample : vcSamplesAndConditions.get(strCondition))
					{
						if(!m_vcSelectedSamples.contains(strSample))
							continue;
						
						double pData[] =  mapData.get(strSample);
						
						for(int x=0; x<pData.length; x++)
						{
							pValues[x][nCurrentSample] = pData[x];
						}
						
						nCurrentSample++;
					}
					
					if(nCurrentSample != 0)
					{
						// get maximum median coverage
						for(int x=0; x<pValues.length; x++)
						{
							if(m_bCoveragePlotShowQuartiles)
								fMaxCoverage = Math.max(fMaxCoverage, StatUtils.percentile(pValues[x], 75));
							else if(m_bCoveragePlotShowMedian)
								fMaxCoverage = Math.max(fMaxCoverage, StatUtils.percentile(pValues[x], 50));
							else
								fMaxCoverage = Math.max(fMaxCoverage, StatUtils.geometricMean(pValues[x]));
						}
					}
				}
			}
		}
		
		double fYShrinkageFactor = (nMaxHeight-20) / fMaxCoverage;	// subtract border size (20)

		//#######################################################
		//                  draw coverage
		//#######################################################
		if(m_bCoveragePlotShowRelativeCoverage)
		{
			//######################################################
			//                   draw axis
			//######################################################
			graph.setFont(fontNormalText);
			graph.setColor(Color.BLACK);
			
			// y-axis
			int nYAxisPos = nOffsetX-10;
			graph.drawLine(nYAxisPos,  10, nYAxisPos,      210);

			for(int i=0; i<11; i++)
			{
				int nYPos = (int) Math.ceil(210.0 - (i*0.1 * 200.0));
				
				graph.drawLine(nYAxisPos-2, nYPos, nYAxisPos, nYPos);			
				
				String strValue = String.format(Locale.ENGLISH, "%1f", 0.1*i);
				FontMetrics fontMetrics = graph.getFontMetrics(fontNormalText);
				int nStringWidth = fontMetrics.stringWidth(strValue);
				int nStringHeight = fontMetrics.getHeight();
				
				graph.drawString(strValue, (int)(nYAxisPos-3-nStringWidth), (int)(nYPos + nStringHeight*0.25));
			}
			
			for(ExonGroup grp : m_pExonGroups)
			{
				int nGrpStart = nOffsetX;
				
				// get correct group coverage data
				String strName = strRef + ":" + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup();
				TreeMap<String, double[]> mapData = mapCoverageToSamplesAndExons.get(strName);
				
				int nWidth = (int)(grp.getExonGroupLengthInBp() * fXShrinkageFactor);
				
				// key = condition, value = coverage per position for each sample
				TreeMap<String, double[][]> mapCoveragePerSamplePerCondition = new TreeMap<String, double[][]>();

				for(String strCondition : vcSamplesAndConditions.keySet())
				{
					int nValidSamples = 0;
					for(String strSample : vcSamplesAndConditions.get(strCondition))
					{
						if(m_vcSelectedSamples.contains(strSample))
							nValidSamples += 1;
					}
	
					if(nValidSamples > 0)
					{
						double[][] pValues =  new double[grp.getExonGroupLengthInBp()+1][nValidSamples];
						
						int nCurrentSample = 0;
						for(String strSample : vcSamplesAndConditions.get(strCondition))
						{
							if(!m_vcSelectedSamples.contains(strSample))
								continue;
							
							double pData[] =  mapData.get(strSample);						
							for(int x=0; x<pData.length; x++)
							{
								double fValue = pData[x];
								pValues[x][nCurrentSample] = fValue;
							}
							
							nCurrentSample += 1;
						}
						
						if(nCurrentSample != 0)
							mapCoveragePerSamplePerCondition.put(strCondition, pValues);
					}
				}
				
				for(int x=1; x<grp.getExonGroupLengthInBp(); x++)
				{			
					// get median coverage for all conditions
					TreeMap<String, Double> mapMedianCoverageToCondition = new TreeMap<String, Double>();
					
					// get coverage sum
					double fTotalCoverage = 0;

					for(String strCondition : mapCoveragePerSamplePerCondition.keySet())
					{
						double pValues[][] = mapCoveragePerSamplePerCondition.get(strCondition);					
						double pCov[] = new double[pValues[0].length];
						
						for(int nSample = 0; nSample<pValues[0].length; nSample++)
						{
							pCov[nSample] = pValues[x][nSample]; 
						}

//						double fMedian = StatUtils.percentile(pCov, 50.0);
						double fMedian = StatUtils.mean(pCov);
						mapMedianCoverageToCondition.put(strCondition, fMedian);
						fTotalCoverage += fMedian;
					}
					
					if(fTotalCoverage == 0)
						continue;

					// calculate contribution of each condition to the coverage
					double fBottom = 210.0;
					for(String strCondition : mapMedianCoverageToCondition.keySet())
					{
						graph.setColor(this.m_mapColorsToConditions.get(strCondition));						
						double fValue = (mapMedianCoverageToCondition.get(strCondition) / fTotalCoverage)*200.0;	// 200 is the plot height
						
						Shape rect = new Rectangle2D.Double(nGrpStart+(x*fXShrinkageFactor), fBottom-fValue, 1, fValue);
						graph.draw(rect);
						fBottom -= fValue;
					}
				}				

				nOffsetX += nWidth + offsets.m_nIntronLength;
			}
		}
		else
		{
			//######################################################
			//                   draw axis
			//######################################################
			graph.setFont(fontNormalText);
			graph.setColor(Color.BLACK);
			
			// y-axis
			int nYAxisPos = nOffsetX-10;
			graph.drawLine(nYAxisPos,  10, nYAxisPos,      210);
			int nStep = (int)Math.ceil(fMaxCoverage / 10.0);	// results in ~10 ticks
			
			if(nStep > 9999) nStep = nStep - (nStep%1000);
			else if(nStep > 999) nStep = nStep - (nStep%100);
			else if(nStep > 10) nStep = nStep - (nStep%10);

			int nVal = 0;		
			while(nVal < fMaxCoverage)
			{
				int nYPos = (int) Math.ceil(210.0 - (nVal / fMaxCoverage * 200.0));
				graph.drawLine(nYAxisPos-2, nYPos, nYAxisPos, nYPos);			
				
				String strValue = ""+nVal;
				FontMetrics fontMetrics = graph.getFontMetrics(fontNormalText);
				int nStringWidth = fontMetrics.stringWidth(strValue);
				int nStringHeight = fontMetrics.getHeight();
				
				graph.drawString(strValue, (int)(nYAxisPos-3-nStringWidth), (int)(nYPos + nStringHeight*0.25));
				
				nVal += nStep;
			}
			
			for(ExonGroup grp : m_pExonGroups)
			{
				int nGrpStart = nOffsetX;
				
				// get correct group coverage data
				String strName = strRef + ":" + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup();
				TreeMap<String, double[]> mapData = mapCoverageToSamplesAndExons.get(strName);
				
				int nWidth = (int)(grp.getExonGroupLengthInBp() * fXShrinkageFactor);
				
				TreeMap<String, double[][]> mapCoverageValuesToConditions = new TreeMap<String, double[][]>();
				
				//##########################################
				//##   combine coverage per condition     ##
				//##   and draw coverage for all samples  ##
				//##   if required.                      ##
				//##########################################
				for(String strCondition : vcSamplesAndConditions.keySet())
				{
					graph.setColor(m_mapColorsToConditions.get(strCondition));
					
					int nSamples = 0;
					for(String strSample : vcSamplesAndConditions.get(strCondition))
					{
						if(m_vcSelectedSamples.contains(strSample))
							nSamples += 1;
					}
	
					if(nSamples > 0)
					{
						double[][] pValues =  new double[grp.getExonGroupLengthInBp()+1][nSamples];
						
						int nCurrentSample = 0;
	
						for(String strSample : vcSamplesAndConditions.get(strCondition))
						{
							if(!m_vcSelectedSamples.contains(strSample))
								continue;
							
							double pData[] =  mapData.get(strSample);
												
							GeneralPath polyline = new GeneralPath(GeneralPath.WIND_EVEN_ODD, pValues.length);
							
							for(int x=0; x<pData.length; x++)
							{
								double fValue = pData[x];
								
								if(m_bCoveragePlotShowMedian || m_bCoveragePlotShowGeometricMean || m_bCoveragePlotShowMean)
								{
									pValues[x][nCurrentSample] = fValue;
								}
								else
								{									
									if(x == 0)
									{
										polyline.moveTo(nGrpStart+(x*fXShrinkageFactor), 210);
									}
									polyline.lineTo(nGrpStart+(x*fXShrinkageFactor), 210-fValue*fYShrinkageFactor);
								}
								
							}
							
							if(!m_bCoveragePlotShowMedian && !m_bCoveragePlotShowGeometricMean && !m_bCoveragePlotShowMean)
							{							
								// move back to 0
								polyline.lineTo(nGrpStart+(pData.length*fXShrinkageFactor), 210);
								graph.draw(polyline);
							}
							
							nCurrentSample += 1;
						}
						mapCoverageValuesToConditions.put(strCondition, pValues);
					}
				}

				if(m_bCoveragePlotShowMedian || m_bCoveragePlotShowGeometricMean || m_bCoveragePlotShowMean)
				{
					//####################################
					//##         draw quantiles         ##
					//####################################
					if(m_bCoveragePlotShowQuartiles)
					{
						for(String strCondition : vcSamplesAndConditions.keySet())
						{
							if(mapCoverageValuesToConditions.containsKey(strCondition))
							{
								double[][] pValues = mapCoverageValuesToConditions.get(strCondition);
							
								Color clrAlpha = new Color(m_mapColorsToConditions.get(strCondition).getRed(), m_mapColorsToConditions.get(strCondition).getGreen(), m_mapColorsToConditions.get(strCondition).getBlue(), 40);
								graph.setColor(clrAlpha);
								
								GeneralPath polyline = new GeneralPath(GeneralPath.WIND_EVEN_ODD, pValues.length);

								for(int x=0; x<pValues.length; x++)
								{
									// get median coverage
									double fQ75 = StatUtils.percentile(pValues[x], 75);
									
									if(x == 0)
									{
										polyline.moveTo(nGrpStart+(x*fXShrinkageFactor), 210);
									}
									
									polyline.lineTo(nGrpStart+(x*fXShrinkageFactor), 210-fQ75*fYShrinkageFactor);
								}
								
								for(int x=(pValues.length-1); x>=0; x--)
								{
									// get median coverage
									double fQ25 = StatUtils.percentile(pValues[x], 25);

									if(x == 0)
									{
										polyline.moveTo(nGrpStart+(x*fXShrinkageFactor), 210);
									}
									
									polyline.lineTo(nGrpStart+(x*fXShrinkageFactor), 210-fQ25*fYShrinkageFactor);
								}
								
								graph.fill(polyline);
							}
						}
					}

					//####################################
					//##           draw lines           ##
					//####################################
					for(String strCondition : vcSamplesAndConditions.keySet())
					{
						if(mapCoverageValuesToConditions.containsKey(strCondition))
						{
							double[][] pValues = mapCoverageValuesToConditions.get(strCondition);
			
							GeneralPath polyline = new GeneralPath(GeneralPath.WIND_EVEN_ODD, pValues.length);
							
							for(int x=0; x<pValues.length; x++)
							{
								// get median coverage
								double fQ50 = 0.0;
								
								if(m_bCoveragePlotShowMedian)
								{
									fQ50 = StatUtils.percentile(pValues[x], 50);
								}
								else if(m_bCoveragePlotShowGeometricMean)
									fQ50 = StatUtils.geometricMean(pValues[x]);
								else
									fQ50 = StatUtils.mean(pValues[x]);
		
								if(x == 0)
								{
									polyline.moveTo(nGrpStart+(x*fXShrinkageFactor), 210);
								}
								
								polyline.lineTo(nGrpStart+(x*fXShrinkageFactor), 210-fQ50*fYShrinkageFactor);
							}
							
							graph.setColor(this.m_mapColorsToConditions.get(strCondition));
							polyline.lineTo(nGrpStart+(pValues.length*fXShrinkageFactor), 210);
							graph.draw(polyline);
						}
					}
				}

				nOffsetX += nWidth + offsets.m_nIntronLength;
			}
		}

		// Show second coverage plot attached to the isoform plot
		if(m_bShowSecondCoveragePlot)
		{
			Imagemap imgMap = new Imagemap();
			imgMap.setWidth(offsets.m_nTotalWidth+"px");
			imgMap.setHeight(nMaxHeight+"px");
			imgMap.setContent(img);
			imgMap.setParent(m_layoutMainBottom);
		}
		
		m_imgMapCoverage = new Imagemap();
		m_imgMapCoverage.setWidth(offsets.m_nTotalWidth+"px");
		m_imgMapCoverage.setHeight(nMaxHeight+"px");
		m_imgMapCoverage.setContent(img);
		m_imgMapCoverage.setParent(m_layoutMainTop);
		
		m_bCoverageRequiresRedraw = false;
	}
	
	private void DrawIsoforms() throws IOException, SQLException
	{
		if(m_pExonGroups == null || m_pExonGroups.length == 0)
			return;
		
		TreeMap<String, TreeSet<String>> mapSamplesToConditions =  m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		TreeMap<String, TreeMap<String, Integer>> mapJunCountsToJunctions = null;
		double fMaxJunCoverage = 0.0;
		if(m_bColorExonsAndJunctionsByCoverage)
		{
			// get junction counts
			mapJunCountsToJunctions = m_projectModel.GetJunctionCountsForGene(m_Gene);
			if(mapJunCountsToJunctions == null)
			{
				System.out.println("failed to get junction read counts for gene: " + m_Gene.getGeneID() + " (" + m_Gene.getGeneName() + ")");
				return;
			}
			
			// get maximum junction coverage
			for(String strJunction : m_Gene.getSpliceJunctionInformation())
			{
				// get median for the current condition
				int nSamples = 0;
				
				for(String strSample : mapSamplesToConditions.get(m_strSelectedCondition))
				{
					if(m_vcSelectedSamples.contains(strSample))
						nSamples++;
				}
				double pCoverages[] = new double[nSamples];
				
				String split[] = strJunction.split("-");
				int nStart = Integer.parseInt(split[0].trim());
				int nEnd   = Integer.parseInt(split[1].trim());
				
				double fCoverage = 0.0;

				if(!mapJunCountsToJunctions.containsKey(nStart + "_" + nEnd))
				{
					//System.out.println("No junction counts available for junction: " + strJunction + ". Are you sure, you have used the correct gene annotation database for the stored junction counts?");
				}	
				else
				{
					TreeMap<String, Integer> mapJunCountsToSamples = mapJunCountsToJunctions.get(nStart + "_" + nEnd);
					
					int nIdx = 0;
					for(String strSample : mapSamplesToConditions.get(m_strSelectedCondition))
					{
						if(m_vcSelectedSamples.contains(strSample))
						{
							if(mapJunCountsToSamples.containsKey(strSample))
								pCoverages[nIdx] = mapJunCountsToSamples.get(strSample);
							else
								pCoverages[nIdx] = 0.0;
							nIdx++;
						}
					}
					
					fCoverage = StatUtils.percentile(pCoverages, 50.0);
					fMaxJunCoverage = Math.max(fMaxJunCoverage, fCoverage);
				}
			}
		}
		
		// get first and second selected exon position
		ExonGroup SelectedExonGrpA = null;
		ExonGroup SelectedExonGrpB = null;

		for(ExonGroup grp : m_pExonGroups)
		{
			if(grp.getGenomicStartOfGroup() >= m_ClickEvent.m_nExonGroupStartA && grp.getGenomicStopOfGroup() <= m_ClickEvent.m_nExonGroupEndA)
			{
				SelectedExonGrpA = grp;
			}
			
			if(grp.getGenomicStartOfGroup() >= m_ClickEvent.m_nExonGroupStartB && grp.getGenomicStopOfGroup() <= m_ClickEvent.m_nExonGroupEndB)
			{
				SelectedExonGrpB = grp;
			}
		}
		
		// get mmseq results (if available)
		TreeMap<String, TreeMap<String, Double>> mapIsoformExpressionValues = m_projectModel.GetIsoformExpressions(m_Gene.getGeneID(), m_strPathMMSeqResults);
		
		// clear all target isoforms
		m_ClickEvent.m_vcTargetIsoforms.clear();

		// get all currently selected/visible isoforms
		TreeMap<String, int[]> mapValidIsoforms = new TreeMap<String, int[]>();
		for(String strIsoform : m_vcValidIsoforms)
		{
			for(String strIncludedIsoform : m_mapExonsToIsoforms.keySet())
			{
				if(strIsoform.equals(strIncludedIsoform) || strIsoform.equals(strIncludedIsoform.split("\\.")[0]))
				{
					mapValidIsoforms.put(strIsoform, m_mapExonsToIsoforms.get(strIncludedIsoform));
					break;
				}
			}
//				mapValidIsoforms.put(strIsoform, m_mapExonsToIsoforms.get(strIsoform));
			m_ClickEvent.m_vcTargetIsoforms.add(strIsoform);
		}

		int nVerticalSpaceBetweenIsoforms = 30;
		
		// if numbers are added to the graph, we need more space between the isoforms
		if(m_bShowRelativeCoverageNumbersExons || m_bColorExonsAndJunctionsByCoverage)
			nVerticalSpaceBetweenIsoforms += 20;
		
		// define selected conditions
		TreeSet<String> vcSelectedConditions = new TreeSet<String>();
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			for(String strSample : m_vcSelectedSamples)
			{
				if(mapSamplesToConditions.get(strCondition).contains(strSample))
					vcSelectedConditions.add(strCondition);
			}
		}
		
		DrawingOffsets offsets = ComputeDrawingOffsets();
		// an additional 30 pixel as margin, and three more isoforms for the meta transcript and overlapping sense + antisense transcripts
		int nMaxHeight  = (mapValidIsoforms.size()+3) * nVerticalSpaceBetweenIsoforms + 20 + 30;
		
		// prepare graph
		BufferedImage img	= new BufferedImage(offsets.m_nTotalWidth, nMaxHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D graph	= img.createGraphics();
		graph.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		graph.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		// set font size
		int nFontSize    = 11;
		Font fontBoldText	= new Font("Courier", Font.BOLD, 13);
		graph.setFont(new Font("Courier", Font.BOLD, nFontSize));
			
		m_layout.setHeight((nMaxHeight+440) + "px");
		m_layout.setVflex("min");
		m_layoutMainBottom.setHeight("100%");
		
		// fill background
		graph.setColor(Color.white);
		graph.fillRect(0, 0, offsets.m_nTotalWidth, nMaxHeight);
		
		// reserve some space for the isoform name
		int nOffsetX = offsets.m_nIsoformNameOffset + offsets.m_nBorderOffset;

		int nCombinedExonLength = 0;
		for(ExonGroup grp : m_pExonGroups)
		{
			nCombinedExonLength += grp.getExonGroupLengthInBp();
		}
		
		m_imgMapIsoforms = new Imagemap();

		// calculate shrinkage factor
		double fShrinkageFactor = (double)offsets.m_nMaxExonLength / (double)nCombinedExonLength;
		
		// calculate start position (=value) of each exon (=key)
		TreeMap<Integer, Integer> mapExonStarts = new TreeMap<Integer, Integer>(); 	// key = exon_id, value = exon_start_pos
		
		int nX = nOffsetX;
		for(ExonGroup grp : m_pExonGroups)
		{
			int nGrpStart = nX;
			Vector<Exon> vcExons = grp.getExons();
			
			for(int i=0; i<vcExons.size(); i++)
			{
				Exon exon = vcExons.get(i);
				int nExonID = grp.getExonIDs().get(i);
				int nOffset = (int) ((exon.getCodingStart() - grp.getGenomicStartOfGroup()) * fShrinkageFactor);
				mapExonStarts.put(nExonID, nGrpStart + nOffset);
			}

			nX += (grp.getExonGroupLengthInBp() * fShrinkageFactor) + offsets.m_nIntronLength;
		}

		Color clrGradient1 = new Color(255,150,150);	// RED
		Color clrGradient2 = new Color(200,0,0);		// RED
		Color clrGradient3 = new Color(217,217,255);	// BLUE
		Color clrGradient4 = new Color(120,120,230);	// BLUE
		Color clrGradient5 = new Color(255,255,255);	// WHITE
		Color clrGradient6 = new Color(150,150,150);	// GRAY
		
		Color clrGradient7 = new Color(255,240,200);	// LIGHT ORANGE
		Color clrGradient8 = new Color(230,180,0);		// ORANGE
		Color clrGradient9 = new Color(200, 80, 0);		// DAKR ORANGE
		
		//######################################
		//    check for alternative splicing
		//######################################
		if(m_bShowASEvents)
		{
			m_vcASCurrentHit = CalculateVariableExons();
		}
		
		//###################################
		//              prepare font
		//###################################
		Font fontNumbers = new Font("Courier", Font.BOLD, 11);
		graph.setFont(fontNumbers);
		FontMetrics fontMetrics = graph.getFontMetrics(fontNumbers);
		
		//#######################################################
		//    get list of unique features of the current gene
		//#######################################################
		TreeSet<Exon> vcUniqueExons = m_Gene.getUniqueExons(null);
		TreeSet<String> vcUniqueJunctions = m_Gene.getUniqueJunctions(null);
		
		//###################################
		//          specify offsets
		//###################################
		int nY = 20; 	// top margin
		nX = nOffsetX;	// left margin
		
		//##############################################
		//         add overlapping transcripts
		//##############################################
		nY = DrawOverlappingTranscripts(nX, nY, offsets.m_nIntronLength, fShrinkageFactor, nVerticalSpaceBetweenIsoforms, graph, fontMetrics);
		
		//#######################################################
		//           add description to meta transcript
		//#######################################################
		graph.setColor(Color.BLACK);
		String strText = "meta exons";
		int nTextWidth = fontMetrics.stringWidth(strText);
		int nTextHeight = fontMetrics.getHeight();
		graph.drawString(strText, nX - nTextWidth - 5, nY + 10 + (int)(nTextHeight*0.25));
		
		//##############################################
		//            draw meta transcript
		//##############################################
		for(ExonGroup grp : m_pExonGroups)
		{
			int nExonGrpID = grp.getGroupID()+1;

			int nWidth = (int)(grp.getExonGroupLengthInBp() * fShrinkageFactor);
			
			String strQuery = m_Gene.getChromosome() + ":" + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup();
			
			Color pClrs[] = new Color[3];
			pClrs[0] = new Color(120, 120, 120);
			pClrs[1] = new Color(200, 200, 200);
			pClrs[2] = new Color(120, 120, 125);
			
			float pFractions[] = new float[3];
			pFractions[0] = 0.0f;
			pFractions[1] = 0.5f;
			pFractions[2] = 1.0f;
			
			boolean bUseRadiantGradient = false;
			if(SelectedExonGrpA != null && SelectedExonGrpA.getGroupID() == grp.getGroupID() || SelectedExonGrpB != null && SelectedExonGrpB.getGroupID() == grp.getGroupID())
			{
				pClrs = new Color[3];
				pClrs[0] = new Color(255, 242, 0);
				pClrs[1] = new Color(200, 190, 0);
				pClrs[2] = new Color(255, 242, 0);
				
				pFractions[0] = 0.8f;
				pFractions[1] = 0.9f;
				pFractions[2] = 1.0f;
				
				Point2D.Double center = new Point2D.Double((float)(nX + nWidth * 0.5), (double)(nY+10));
				float fRadius = nWidth*0.48f;
				
				RadialGradientPaint clrGradient = new RadialGradientPaint(center, fRadius, pFractions, pClrs, CycleMethod.NO_CYCLE);
				
				// draw exon
				graph.setPaint(clrGradient);
				graph.fillRect(nX, nY, nWidth, 20);
				
				bUseRadiantGradient = true;
			}

			// show alternatively spliced exon groups
			boolean bIsASHit = false;
			if(m_bShowASEvents && m_vcASCurrentHit != null)
			{
//				if(m_vcASCurrentHit.ContainsExonGroup(strQuery))
				if(m_vcASCurrentHit.IsContainedInRegion(m_Gene.getChromosome(), grp.getGenomicStartOfGroup(), grp.getGenomicStopOfGroup()))
				{
					pClrs[0] = new Color(255, 180, 180);
					pClrs[1] = new Color(255,  60,  60);
					pClrs[2] = new Color(255, 180, 180);

					pFractions[0] = 0.0f;
					pFractions[1] = 0.5f;
					pFractions[2] = 1.0f;
					bIsASHit = true;
				}
			}
			
			if(!bUseRadiantGradient)
			{
				LinearGradientPaint clrGradient = new LinearGradientPaint(nX, nY, nX, nY+20, pFractions, pClrs);
				
				// draw exon
				graph.setPaint(clrGradient);
				graph.fillRect(nX, nY, nWidth, 20);
			}

			// add entropy levels if available
			if(m_bShowEntropyData && m_mapEntropyToExonicPart != null)
			{
				// clear gradient
				graph.setPaint(Color.WHITE);
				graph.fillRect(nX, nY, nWidth, 20);
				
				for(Map.Entry<String, Double> e : m_mapEntropyToExonicPart.entrySet())
				{
					String pSplit[] = e.getKey().split("-");
					int nStart = Integer.parseInt(pSplit[0]);
					int nEnd   = Integer.parseInt(pSplit[1]);
					
					if(nStart >= grp.getGenomicStartOfGroup() && nEnd <= grp.getGenomicStopOfGroup())
					{
						int nOffset = (int) ((nStart - grp.getGenomicStartOfGroup()) * fShrinkageFactor);
						int nWidthPart = (int)Math.ceil((nEnd - nStart+1)*fShrinkageFactor);
						
						double fVal = e.getValue();
						
						int nColorRed = (int)Math.max(0, Math.min(255.0, 255.0*fVal));
						int nColorBlue = (int)Math.max(0, Math.min(255.0, 128.0-128.0*fVal));
						
						float pWeights[] = {0.0f,0.2f,0.8f,1.0f};
						Color clrA = new Color(nColorRed, 0, nColorBlue);
						Color clrB = new Color((int)(clrA.getRed()*0.8), (int)(clrA.getGreen()*0.8), (int)(clrA.getBlue()*0.8));
						Color pColors[] = {clrB, clrA, clrA, clrB};
						
						LinearGradientPaint clrGradientEntropy = new LinearGradientPaint((int)(nX+nWidth*0.5), nY, (int)(nX+nWidth*0.5), nY+20, pWeights, pColors);
						
						graph.setPaint(clrGradientEntropy);
						graph.fillRect(nX+nOffset, nY, nWidthPart, 20);
						
						// add text
						strText = String.format(Locale.ENGLISH, "%.2f", fVal);
						nTextWidth = fontMetrics.stringWidth(strText);
						nTextHeight = fontMetrics.getHeight();
						graph.drawString(strText, nX+nOffset + (int)(nWidthPart*0.5) - (int)(nTextWidth*0.5), (int)(nY-nTextHeight*0.5));
						
						// add clickable area
						Area area = new Area();
						area.setId("exonic_part " + e.getKey());
						area.setShape("rectangle");
						area.setCoords(nX+nOffset + ", " + nY + ", " + (nX+nOffset+nWidthPart) + ", " + (nY+20));				
						area.setTooltiptext("Show exonic part coverage " + e.getKey() + " entropy: " + e.getValue());
						area.setParent(m_imgMapIsoforms);
					}
				}
			}
			
			// draw exon border
			graph.setColor(Color.BLACK);
			graph.drawRect(nX, nY, nWidth, 20);
			
			// add orientation
			if(nWidth >= 20) // only for exons with a minimum size
			{
				strText = "";
				
				int nTextLength = (int)Math.floor((nWidth-12)/7);
				
				for(int i=0; i<nTextLength; i++)
				{
					if(m_Gene.isPlusStrand())
						strText += ">";
					else
						strText += "<";
				}
				
				nTextWidth = fontMetrics.stringWidth(strText);
				nTextHeight = fontMetrics.getHeight();
				
				int nXOffset = (int)((nWidth - nTextWidth) * 0.5);
				
				graph.drawString(strText, nX + nXOffset, (int)(nY+6+nTextHeight*0.5));
			}
			
			// add clickable area
			Area area = new Area();
			if(bIsASHit)
				area.setId("highlighted_exon_group " + nExonGrpID + " " + grp.getGenomicStartOfGroup() + " " + grp.getGenomicStopOfGroup());
			else
				area.setId("exon_group " + nExonGrpID + " " + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup());
			area.setShape("rectangle");
			area.setCoords(nX + ", " + nY + ", " + (nX+nWidth) + ", " + (nY+20));
			
			String strToolTip = "exon group " + nExonGrpID + ", " + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup();
			if(m_vcASCurrentHit != null && m_vcASCurrentHit.ContainsExonGroup(strQuery))
				strToolTip += " Rating: " + m_vcASCurrentHit.GetRating();
			area.setTooltiptext(strToolTip);
			area.setParent(m_imgMapIsoforms);
			
			nX += nWidth + offsets.m_nIntronLength;
		}
		
		nY += nVerticalSpaceBetweenIsoforms;
		//##############################################
		//               draw isoforms
		//##############################################
		m_ClickEvent.m_vcSelectedExons.clear();
		for(Map.Entry<String, int[]> e : mapValidIsoforms.entrySet())
		{
			int nPreviousExonEnd   = -1;
			int nPreviousExonID    = -1;
			int nPreviousExonGrpID = -1;
			boolean bIsFirstExon = true;
			
			//#######################################################
			//                add transcript name
			//#######################################################
			graph.setColor(Color.BLACK);
			strText = e.getKey();
			nTextWidth = fontMetrics.stringWidth(strText);
			nTextHeight = fontMetrics.getHeight();
			graph.drawString(strText, nOffsetX - nTextWidth - 5, nY + 8 + (int)(nTextHeight*0.5));
			
			// get maximum exon value
			double fMaxAbsoluteCoverage = 0.0;
			if(m_bColorExonsAndJunctionsByCoverage)
			{
				for(int nExon : e.getValue())
				{
					int nExonID = Math.abs(nExon - (m_pExons.length-1));
					
					// get exon group ID
					ExonGroup exonGroup = null;
					for(ExonGroup grp : m_pExonGroups)
					{
						if(grp.groupContainsExon(nExonID))
						{
							exonGroup = grp;
						}
					}
					
					TreeMap<String, Double> mapAbsoluteCoveragesPerSample = GetAbsoluteCoveragePerSample(m_pExons[nExonID], exonGroup);
					
					int nSamples = 0;
					for(String strSample : mapSamplesToConditions.get(m_strSelectedCondition))
					{
						if(m_vcSelectedSamples.contains(strSample))
							nSamples++;
					}
					double pCoverages[] = new double[nSamples];
					
					int nIdx = 0;
					for(String strSample : mapSamplesToConditions.get(m_strSelectedCondition))
					{
						if(m_vcSelectedSamples.contains(strSample))
						{
							if(mapAbsoluteCoveragesPerSample.containsKey(strSample))
								pCoverages[nIdx] = mapAbsoluteCoveragesPerSample.get(strSample);
							else
								pCoverages[nIdx] = 0.0;
							nIdx++;
						}
					}
					
					double fCoverage = StatUtils.percentile(pCoverages, 50.0);
					
					fMaxAbsoluteCoverage = Math.max(fMaxAbsoluteCoverage, fCoverage);
				}
			}
			
			for(int nExon : e.getValue())
			{
				int nExonID = Math.abs(nExon - (m_pExons.length-1));

				// get ambigious exons within an exon group per condition
				TreeMap<String, TreeSet<Exon>> mapAmbigiousExonsPerCondition = null;
				if(m_mapIndistinguishableExonPerSample != null)
				{
					mapAmbigiousExonsPerCondition = GetMostFrequentSetOfAmbigiousExons(m_pExons[nExonID]);
				}
				
				// get exon group ID
				ExonGroup exonGroup = null;
				int nExonGrpID = -1;
				for(ExonGroup grp : m_pExonGroups)
				{
					if(grp.groupContainsExon(nExonID))
					{
						nExonGrpID = grp.getGroupID()+1;
						exonGroup = grp;
					}
				}
				
				boolean bIsHighlighted = false;
				if(m_vcHighlightedExons.size() > 0 && m_vcHighlightedExons.contains(m_pExons[nExonID]))
				{
					bIsHighlighted = true;
				}

				nX = mapExonStarts.get(nExonID);
				int nWidth = (int)(m_pExons[nExonID].getLength() * fShrinkageFactor);

				// exon fill color depends on strand orientation
				GradientPaint clrGradient = null;

				if(SelectedExonGrpA != null && SelectedExonGrpA.getGroupID() == exonGroup.getGroupID())
				{
					m_ClickEvent.m_vcSelectedExons.add(m_pExons[nExonID]);
				}
				else if(SelectedExonGrpA != null && SelectedExonGrpA.getGroupID() == exonGroup.getGroupID())
				{
					m_ClickEvent.m_vcSelectedExons.add(m_pExons[nExonID]);
				}
				clrGradient = new GradientPaint((int)(nX+nWidth*0.5), nY, clrGradient3, (int)(nX+nWidth*0.5), nY+20, clrGradient4);
				
				if(m_bShowUniqueFeatures)
				{
					if(vcUniqueExons.contains(m_pExons[nExonID]))
					{
						if(mapAmbigiousExonsPerCondition != null && m_bShowAmbigiousUniqueExons)
						{
							TreeSet<Exon> vcAmbigiousExons = mapAmbigiousExonsPerCondition.get(m_strSelectedCondition);
						
							if(vcAmbigiousExons.size() == 1)
							{
								clrGradient = new GradientPaint((int)(nX+nWidth*0.5), nY, clrGradient7, (int)(nX+nWidth*0.5), nY+20, clrGradient8);
							}
							else
							{
								for(Exon ex : vcAmbigiousExons)
								{
									if(ex.equals(m_pExons[nExonID]))
									{
										clrGradient = new GradientPaint((int)(nX+nWidth*0.5), nY, clrGradient8, (int)(nX+nWidth*0.5), nY+20, clrGradient9);
										break;
									}
								}
							}
						}
						else
						{
							clrGradient = new GradientPaint((int)(nX+nWidth*0.5), nY, clrGradient7, (int)(nX+nWidth*0.5), nY+20, clrGradient8);
						}
					}
				}
				
				graph.setPaint(clrGradient);
				//########################################################
				//          add additional exon information
				//########################################################
				if(m_mapExonCoveragePerSample != null && m_mapIndistinguishableExonPerSample != null)
				{
					TreeMap<String, Double> mapRelativeCoveragesPerSample 		= CalculateRelativeExonExpressionPerSample(m_pExons[nExonID], exonGroup);
					TreeMap<String, Double> mapRelativeCoveragesPerCondition 	= CalculateRelativeExonExpressionPerCondition(mapRelativeCoveragesPerSample);					
					
					if(mapRelativeCoveragesPerCondition != null)
					{
						if(m_bShowRelativeCoverageNumbersExons)
						{
							graph.setColor(Color.BLACK);
							
							strText = "";
							String strCondition = m_strSelectedCondition;
							{
								double fCoverage = mapRelativeCoveragesPerCondition.get(strCondition);
								if(fCoverage == 100.0)
									strText += String.format(Locale.ENGLISH, "100");
								else
									strText += String.format(Locale.ENGLISH, "%.2f", fCoverage);
							}

							// get text width to center the text above the exon
							nTextWidth = fontMetrics.stringWidth(strText);
							graph.drawString(strText, nX + (int)(nWidth*0.5) - (int)(nTextWidth*0.5), nY-2);
						}
						else if(m_bColorExonsAndJunctionsByCoverage)
						{
							graph.setColor(Color.BLACK);
							
							TreeMap<String, Double> mapAbsoluteCoveragesPerSample = GetAbsoluteCoveragePerSample(m_pExons[nExonID], exonGroup);
																					
//							System.out.println(mapAbsoluteCoveragesPerSample);
							
							//##################################
							//            add text
							//##################################
							strText = "";
							
							int nSamples = 0;
							for(String strSample : mapSamplesToConditions.get(m_strSelectedCondition))
							{
								if(m_vcSelectedSamples.contains(strSample))
									nSamples++;
							}
							
							double pCoverages[] = new double[nSamples];
							
							int nIdx = 0;
							for(String strSample : mapSamplesToConditions.get(m_strSelectedCondition))
							{
								if(m_vcSelectedSamples.contains(strSample))
								{
									if(mapAbsoluteCoveragesPerSample.containsKey(strSample))
										pCoverages[nIdx] = mapAbsoluteCoveragesPerSample.get(strSample);
									else
										pCoverages[nIdx] = 0.0;
									nIdx++;
								}
							}
							
							double fCoverage = StatUtils.percentile(pCoverages, 50.0);

							strText += String.format(Locale.ENGLISH, "%.1f", fCoverage);

							// get text width to center the text above the exon
							nTextWidth = fontMetrics.stringWidth(strText);
							graph.setColor(Color.BLACK);
							graph.drawString(strText, nX + (int)(nWidth*0.5) - (int)(nTextWidth*0.5), nY-2);
							
							//##################################
							//             add color
							//##################################
							double fVal = fCoverage / fMaxAbsoluteCoverage;
								
							if(fVal < 0.01 && !bIsHighlighted)
							{
								graph.setPaint(Color.WHITE);
							}
							else
							{
								int nColorValue = 255-Math.min(255, (int)(fVal*255.0));
								
								if(bIsHighlighted)
								{
									Color clrExpressionLight = new Color(nColorValue, 0, 0);	// LIGHT color
									nColorValue *= 0.8;
									Color clrExpressionDark  = new Color(nColorValue, 0, 0);	// DARK color 
									
									clrGradient = new GradientPaint((int)(nX+nWidth*0.5), nY, clrExpressionLight, (int)(nX+nWidth*0.5), nY+20, clrExpressionDark);
									graph.setPaint(clrGradient);
								}
								else
								{
									Color clrExpressionLight = new Color(nColorValue, nColorValue, nColorValue);	// LIGHT color
									nColorValue *= 0.8;
									Color clrExpressionDark  = new Color(nColorValue, nColorValue, nColorValue);	// DARK color 
									
									clrGradient = new GradientPaint((int)(nX+nWidth*0.5), nY, clrExpressionLight, (int)(nX+nWidth*0.5), nY+20, clrExpressionDark);
									graph.setPaint(clrGradient);
								}
							}
						}

						if(m_bColorExonsByRelativeCoverage)
						{
							graph.setColor(Color.BLACK);
							
							double fVal = mapRelativeCoveragesPerCondition.get(m_strSelectedCondition);
							
							if(fVal == 0.0)
							{
								graph.setPaint(Color.WHITE);
							}
							else
							{
								int nColorValue = 255-Math.min(255, (int)(fVal/100.0*255.0));
								
								if(bIsHighlighted)
								{
									Color clrExpressionLight = new Color(nColorValue, 0, 0);	// LIGHT color
									nColorValue *= 0.8;
									Color clrExpressionDark  = new Color(nColorValue, 0, 0);	// DARK color 
									
									clrGradient = new GradientPaint((int)(nX+nWidth*0.5), nY, clrExpressionLight, (int)(nX+nWidth*0.5), nY+20, clrExpressionDark);
									graph.setPaint(clrGradient);
								}
								else
								{
									Color clrExpressionLight = new Color(nColorValue, nColorValue, nColorValue);	// LIGHT color
									nColorValue *= 0.8;
									Color clrExpressionDark  = new Color(nColorValue, nColorValue, nColorValue);	// DARK color 
									
									clrGradient = new GradientPaint((int)(nX+nWidth*0.5), nY, clrExpressionLight, (int)(nX+nWidth*0.5), nY+20, clrExpressionDark);
									graph.setPaint(clrGradient);
								}
							}
						}
					}
				}
				
				// do not color exons that are not highlighted, if there are any highlighted exons
				if(!bIsHighlighted && m_vcHighlightedExons.size() > 0)
				{
					graph.setPaint(Color.WHITE);
				}
				
				//###################################
				//             draw exon
				//###################################				
				graph.fillRect(nX, nY, nWidth, 20);

				graph.setColor(Color.BLACK);
				graph.drawRect(nX, nY, nWidth, 20);
				
				Area area = new Area();
				area.setId("exon_" + nExonGrpID + " " + e.getKey() + " " + nExonID);
				area.setShape("rectangle");
				area.setCoords(nX + ", " + nY + ", " + (nX+nWidth) + ", " + (nY+20));				
				area.setTooltiptext(e.getKey() + " exon_" + nExonGrpID + ", " + m_pExons[nExonID].getCodingStart() + "-" + m_pExons[nExonID].getCodingStop() + ", (" + nExonID + ")");
				area.setParent(m_imgMapIsoforms);

				//#############################################
				//               draw introns 
				//#############################################
				if(bIsFirstExon)
					bIsFirstExon = false;
				else
				{
					int nIntronStart = nPreviousExonEnd;
					int nIntrnLength = nX-nPreviousExonEnd;	// = length - 1
					
					String strJunctionName = (m_pExons[nPreviousExonID].getCodingStop()+1) + " - " + (m_pExons[nExonID].getCodingStart()-1);
					String strJunctionName2 = (m_pExons[nPreviousExonID].getCodingStop()) + "-" + (m_pExons[nExonID].getCodingStart());
					
					if(m_ClickEvent.m_strClickedA != null && m_ClickEvent.m_strClickedA.equals("edge_"+nPreviousExonGrpID+"_"+nExonGrpID))
					{
						clrGradient = new GradientPaint((int)(nIntronStart+(nIntrnLength)*0.5), nY+7, clrGradient1, (int)(nIntronStart+(nIntrnLength)*0.5), nY+13, clrGradient2);
					}
					else if(m_ClickEvent.m_strClickedA != null && m_ClickEvent.m_strClickedA.equals("edge_"+nPreviousExonGrpID+"_"+nExonGrpID))
					{
						clrGradient = new GradientPaint((int)(nIntronStart+(nIntrnLength)*0.5), nY+7, clrGradient2, (int)(nIntronStart+(nIntrnLength)*0.5), nY+13, clrGradient1);
					}
					else
					{
						clrGradient = new GradientPaint((int)(nIntronStart+(nIntrnLength)*0.5), nY+7, clrGradient5, (int)(nIntronStart+(nIntrnLength)*0.5), nY+13, clrGradient6);
					}
					
					if(m_bShowUniqueFeatures)
					{
						if(vcUniqueJunctions.contains(strJunctionName2))
							clrGradient = new GradientPaint((int)(nIntronStart+(nIntrnLength)*0.5), nY+7, clrGradient7, (int)(nIntronStart+(nIntrnLength)*0.5), nY+13, clrGradient8);
					}
					
					bIsHighlighted = false;
					if(m_vcHighlightedJunctions.size() > 0)
					{
						for(CountElement jun : m_vcHighlightedJunctions)
						{
							String pSplit[] = strJunctionName.split("-");
							int nJunctionStart = Integer.parseInt(pSplit[0].trim())-1;
							int nJunctionEnd = Integer.parseInt(pSplit[1].trim())+1;
							
							if(jun.m_nStart == nJunctionStart && jun.m_nEnd == nJunctionEnd)
							{
								bIsHighlighted = true;
								break;
							}
						}
					}

					if(m_bColorJunctionPaths || m_bShowRelativeCoverageNumbersJunctions)
					{
						String strIsoform = e.getKey();
						if(m_mapAnalysedJunctionPaths.containsKey(strIsoform))
						{
							TreeMap<Integer, TreeSet<String>> mapJunctionPath = m_mapAnalysedJunctionPaths.get(strIsoform); 
							int nPathID = mapJunctionPath.firstKey();

							String split[] = strJunctionName.split("-");
							int nStart = Integer.parseInt(split[0].trim())-1;
							int nEnd   = Integer.parseInt(split[1].trim())+1;							
							if(mapJunctionPath.get(nPathID).contains(nStart + "_" + nEnd))
							{
								double fVal = CalculateRelativeJunctionExpressionPerCondition(nPathID).get(m_strSelectedCondition)*100.0;
								
								if(m_bShowRelativeCoverageNumbersJunctions)
								{
									graph.setPaint(Color.BLUE);
									strText = String.format(Locale.ENGLISH, "%.1f", fVal);
									if(fVal == 100.0)
										strText = "100";

									nTextWidth = fontMetrics.stringWidth(strText);
									nTextHeight = fontMetrics.getHeight();
									graph.drawString(strText, nIntronStart + (int)(nIntrnLength*0.5) - (int)(nTextWidth*0.5), nY+14+nTextHeight);
									graph.setPaint(Color.BLACK);
								}

								if(m_bColorJunctionPaths)
								{
									if(fVal == 0.0)
									{
										graph.setPaint(Color.WHITE);
									}
									else
									{
										int nColorValue = 255-(int)(fVal/100*255);
										
										Color clrExpressionLight = new Color(nColorValue, nColorValue, nColorValue);	// LIGHT color
										nColorValue *= 0.8;
										Color clrExpressionDark  = new Color(nColorValue, nColorValue, nColorValue);	// DARK color 
										
										clrGradient = new GradientPaint((int)(nX+nWidth*0.5), nY, clrExpressionLight, (int)(nX+nWidth*0.5), nY+20, clrExpressionDark);
										graph.setPaint(clrGradient);
									}
	
									graph.fillRect(nIntronStart, nY+7, nIntrnLength, 6);
								}
							}
							else if(m_bColorJunctionPaths)
							{
								graph.setPaint(Color.BLUE);
								graph.fillRect(nIntronStart, nY+7, nIntrnLength, 6);
							}
						}
						else
						{
							if(m_bColorJunctionPaths)
							{
								graph.setPaint(Color.WHITE);
								graph.fillRect(nIntronStart, nY+7, nIntrnLength, 6);
							}
							
							if(m_bShowRelativeCoverageNumbersJunctions)
							{
								graph.setPaint(Color.BLUE);
								strText = "n.e.";
								nTextWidth = fontMetrics.stringWidth(strText);
								nTextHeight = fontMetrics.getHeight();
								graph.drawString(strText, nIntronStart + (int)(nIntrnLength*0.5) - (int)(nTextWidth*0.5), nY+14+nTextHeight);
								graph.setPaint(Color.BLACK);
							}
						}
					}
					else if(m_bColorExonsAndJunctionsByCoverage)
					{
						//##################################
						//             add text
						//##################################
						String split[] = strJunctionName.split("-");
						int nStart = Integer.parseInt(split[0].trim())-1;
						int nEnd   = Integer.parseInt(split[1].trim())+1;
						
						double fCoverage = 0.0;
						
						TreeMap<String, Integer> mapJunCountsToSamples = mapJunCountsToJunctions.get(nStart + "_" + nEnd);
						
						if(!mapJunCountsToJunctions.containsKey(nStart + "_" + nEnd))
						{
							//System.out.println("No junction counts available for junction: " + strJunctionName + ". Are you sure, you have used the correct gene annotation database for the stored junction counts?");
						}
						else
						{
							// get median for the current condition
							int nSamples = 0;
							for(String strSample : mapSamplesToConditions.get(m_strSelectedCondition))
							{
								if(m_vcSelectedSamples.contains(strSample))
									nSamples++;
							}
	
							double pCoverages[] = new double[nSamples];
							
							int nIdx = 0;
							for(String strSample : mapSamplesToConditions.get(m_strSelectedCondition))
							{
								if(m_vcSelectedSamples.contains(strSample))
								{
									if(mapJunCountsToSamples.containsKey(strSample))
										pCoverages[nIdx] = mapJunCountsToSamples.get(strSample);
									else
										pCoverages[nIdx] = 0.0;
									nIdx++;
								}
							}
							
							fCoverage = StatUtils.percentile(pCoverages, 50.0);
						}
						
						graph.setPaint(Color.BLUE);
						strText = String.format(Locale.ENGLISH, "%.0f", fCoverage);

						nTextWidth = fontMetrics.stringWidth(strText);
						nTextHeight = fontMetrics.getHeight();
						graph.drawString(strText, nIntronStart + (int)(nIntrnLength*0.5) - (int)(nTextWidth*0.5), nY+18+nTextHeight);
						graph.setPaint(Color.BLACK);
						
						//##################################
						//             add color
						//##################################
						double fRelCoverage = fCoverage/fMaxJunCoverage;
						
						if(fRelCoverage == 0)
						{
							graph.setPaint(Color.WHITE);
						}
						else
						{
							int nColorValue = 255-(int)(fRelCoverage*255);

							if(bIsHighlighted)
							{								
								Color clrExpressionLight = new Color(nColorValue, 0, 0);	// LIGHT color
								nColorValue *= 0.8;
								Color clrExpressionDark  = new Color(nColorValue, 0, 0);	// DARK color 
								
								clrGradient = new GradientPaint((int)(nX+nWidth*0.5), nY, clrExpressionLight, (int)(nX+nWidth*0.5), nY+20, clrExpressionDark);
								graph.setPaint(clrGradient);
							}
							else
							{
								Color clrExpressionLight = new Color(nColorValue, nColorValue, nColorValue);	// LIGHT color
								nColorValue *= 0.8;
								Color clrExpressionDark  = new Color(nColorValue, nColorValue, nColorValue);	// DARK color 
								
								clrGradient = new GradientPaint((int)(nX+nWidth*0.5), nY, clrExpressionLight, (int)(nX+nWidth*0.5), nY+20, clrExpressionDark);
								graph.setPaint(clrGradient);
							}
						}
					}
					else
					{
						graph.setPaint(clrGradient);
					}
					
					// do not color junctions that were not highlighted (if any where highlighted)
					if(!bIsHighlighted && m_vcHighlightedJunctions.size() > 0)
					{
						graph.setPaint(Color.WHITE);
					}
					
					// fill the junction
					graph.fillRect(nIntronStart, nY+7, nIntrnLength, 6);
					
					// intron border
					graph.setColor(Color.BLACK);
					graph.drawRect(nIntronStart, nY+7, nIntrnLength, 6);
					
					area = new Area();
					area.setId("edge_" + nPreviousExonGrpID + "_" + nExonGrpID + " " + e.getKey() + " " + nPreviousExonID + " " + nExonID);
					area.setShape("rectangle");
					area.setCoords(nIntronStart + ", " + (nY+7) + ", " + (nIntronStart+nIntrnLength) + ", " + (nY+13));				
					area.setTooltiptext(e.getKey() + " edge_" + nPreviousExonGrpID + "_" + nExonGrpID + ", " + strJunctionName);
					area.setParent(m_imgMapIsoforms);
				}
				
				nPreviousExonGrpID 	= nExonGrpID;
				nPreviousExonID		= nExonID;
				nPreviousExonEnd 	= nX + nWidth;
			}

			//#######################################
			//    add mmseq estimate if available
			//#######################################
			if(mapIsoformExpressionValues != null)
			{				
				TreeMap<String, Double> mapIsoformExpressionToSamples = null;
				if(mapIsoformExpressionValues.containsKey(e.getKey()))
				{
					mapIsoformExpressionToSamples = mapIsoformExpressionValues.get(e.getKey());
				}
				else
				{
					System.out.println("failed to get mmseq results for isoform: " + e.getKey());
					System.out.println(mapIsoformExpressionValues);
				}
				
				if(mapIsoformExpressionToSamples != null)
				{
					strText = "(";
					boolean bFirst = true;
					boolean bHasLargerPercentage = false;
					
					for(String strCondition : mapSamplesToConditions.keySet())
					{
						double fSum = 0.0;
						int nSamples = 0;
						for(String strSample : mapSamplesToConditions.get(strCondition))
						{
							if(!mapIsoformExpressionToSamples.containsKey(strSample))
							{
								System.out.println("failed to get mmseq results for sample: " + strSample);
								System.out.println(mapIsoformExpressionToSamples);
								continue;
							}
							fSum += mapIsoformExpressionToSamples.get(strSample);
							nSamples += 1;
						}
						
						double fAvg = fSum / nSamples;
						if(fAvg >= 0.1)
							bHasLargerPercentage = true;
						
						if(bFirst)
						{
							strText += strCondition + ": " + String.format(Locale.ENGLISH, "%.2f%%", fAvg*100.0);						
							bFirst = false;
						}
						else
						{
							strText += "; " + strCondition + ": " + String.format(Locale.ENGLISH, "%.2f%%", fAvg*100.0);						
						}					
					}
					
					strText += ")";
					
					if(!bHasLargerPercentage)
					{
						fontNumbers = new Font("Courier", Font.PLAIN, 11);
						graph.setFont(fontNumbers);
					}
					
					nTextWidth = fontMetrics.stringWidth(strText);
					nTextHeight = fontMetrics.getHeight();
					graph.drawString(strText, nPreviousExonEnd+10, nY + 9 + (int)(nTextHeight*0.5));
					
					if(!bHasLargerPercentage)
					{
						fontNumbers = new Font("Courier", Font.BOLD, 11);  
						graph.setFont(fontNumbers);
					}
				}
			}
			
			nY += nVerticalSpaceBetweenIsoforms;
		}
				
		m_imgMapIsoforms.setWidth(offsets.m_nTotalWidth+"px");
		m_imgMapIsoforms.setHeight(nMaxHeight+"px");
		m_imgMapIsoforms.setContent(img);
		m_imgMapIsoforms.setParent(m_layoutMainBottom);
		
		Area areaBackground = new Area();
		areaBackground.setCoords("0, 0, " + offsets.m_nTotalWidth + ", " + nMaxHeight);
		areaBackground.setId("background_isoforms");
		areaBackground.setParent(m_imgMapIsoforms);

		m_imgMapIsoforms.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				MouseEvent evnt = (MouseEvent) event;

				// find selected exon groups
				// get first and second selected exon position
				ExonGroup SelectedExonGrpA = null;
				ExonGroup SelectedExonGrpB = null;

				for(ExonGroup grp : m_pExonGroups)
				{
					if(grp.getGenomicStartOfGroup() >= m_ClickEvent.m_nExonGroupStartA && grp.getGenomicStopOfGroup() <= m_ClickEvent.m_nExonGroupEndA)
					{
						SelectedExonGrpA = grp;
					}
					
					if(grp.getGenomicStartOfGroup() >= m_ClickEvent.m_nExonGroupStartB && grp.getGenomicStopOfGroup() <= m_ClickEvent.m_nExonGroupEndB)
					{
						SelectedExonGrpB = grp;
					}
				}
				
				if(evnt.getArea().startsWith("edge"))
				{
					m_ClickEvent.m_strClickedA = evnt.getArea().split("\\s")[0];

					m_ClickEvent.m_nExonGroupStartA = -1;
					m_ClickEvent.m_nExonGroupEndA 	= -1;
					
					m_ClickEvent.m_strClickedB 		= null;					
					m_ClickEvent.m_nExonGroupStartB = -1;
					m_ClickEvent.m_nExonGroupEndB 	= -1;
				}
				else if(evnt.getArea().startsWith("background_isoforms"))
				{
					m_ClickEvent.m_strClickedA 		= null;
					m_ClickEvent.m_nExonGroupStartA = -1;
					m_ClickEvent.m_nExonGroupEndA 	= -1;
					
					m_ClickEvent.m_strClickedB 		= null;					
					m_ClickEvent.m_nExonGroupStartB = -1;
					m_ClickEvent.m_nExonGroupEndB 	= -1;
					
					m_vcHighlightedExons.clear();
					m_vcHighlightedJunctions.clear();
				}
				else if(evnt.getArea().startsWith("exonic_part"))
				{
					ShowGTEXDataForExon(evnt.getArea());
					return;
				}
				else if(evnt.getArea().startsWith("exon_group"))
				{
					String strExonGroupPosition = evnt.getArea().split("\\s")[2];
					int nClickedGrpStart = Integer.parseInt(strExonGroupPosition.split("-")[0]);
					int nClickedGrpEnd   = Integer.parseInt(strExonGroupPosition.split("-")[1]);
					
					if(m_ClickEvent.m_strClickedA != null && evnt.getArea().contains("edge")) m_ClickEvent.m_strClickedA = null;
					if(m_ClickEvent.m_strClickedB != null && evnt.getArea().contains("edge")) m_ClickEvent.m_strClickedB = null;
					
					// unselect previously selected exon groups on reselect
					boolean bReSelected = false;
					String strClicked = evnt.getArea().split("\\s")[0] + " " + evnt.getArea().split("\\s")[1];
					if(m_ClickEvent.m_strClickedA != null && nClickedGrpStart == SelectedExonGrpA.getGenomicStartOfGroup() && nClickedGrpEnd == SelectedExonGrpA.getGenomicStopOfGroup())
					{
						if(m_ClickEvent.m_strClickedB != null)
						{
							m_ClickEvent.m_strClickedA 	= m_ClickEvent.m_strClickedB;
							
							m_ClickEvent.m_nExonGroupStartA = m_ClickEvent.m_nExonGroupStartB;
							m_ClickEvent.m_nExonGroupEndA	= m_ClickEvent.m_nExonGroupEndB;
							
							m_ClickEvent.m_strClickedB 		= null;
							m_ClickEvent.m_nExonGroupStartB = -1;
							m_ClickEvent.m_nExonGroupEndB 	= -1;
						}
						else
						{
							m_ClickEvent.m_strClickedA 		= null;
							m_ClickEvent.m_nExonGroupStartA = -1;
							m_ClickEvent.m_nExonGroupEndA 	= -1;
						}
						bReSelected = true;
					}
					
					if(m_ClickEvent.m_strClickedB != null && nClickedGrpStart == SelectedExonGrpB.getGenomicStartOfGroup() && nClickedGrpEnd == SelectedExonGrpB.getGenomicStopOfGroup())
					{
						m_ClickEvent.m_strClickedB 		= null;
						m_ClickEvent.m_nExonGroupStartB = -1;
						m_ClickEvent.m_nExonGroupEndB 	= -1;
						bReSelected = true;
					}
					
					// select new exon
					if(!bReSelected)
					{
						m_ClickEvent.m_strClickedB = m_ClickEvent.m_strClickedA;
						m_ClickEvent.m_strClickedA = strClicked;

						m_ClickEvent.m_nExonGroupStartB = m_ClickEvent.m_nExonGroupStartA;
						m_ClickEvent.m_nExonGroupEndB	= m_ClickEvent.m_nExonGroupEndA;
						
						m_ClickEvent.m_nExonGroupStartA = nClickedGrpStart;
						m_ClickEvent.m_nExonGroupEndA 	= nClickedGrpEnd;
					}
					
					if(m_Gene != null)
					{
						UpdateComboboxForIsoformSelection(false);
						DrawPlots();
					}
				}
				else if(evnt.getArea().startsWith("highlighted_exon_group"))
				{
					String pSplit[] = evnt.getArea().split("\\s");
					
					int nGroupStart  = Integer.parseInt(pSplit[2]);
					int nGroupEnd	 = Integer.parseInt(pSplit[3]);					
					
					Popup popup = new Popup();
					popup.setParent(m_hWindow);
					popup.open(m_hWindow, "after_pointer");
					popup.setHflex("min");
					popup.setVflex("min");
					
					Grid grid = new Grid();
					grid.setMold("default");
					grid.setWidth("800px");
					grid.setVflex("min");
					grid.setParent(popup);
					
					Columns cols = new Columns();
					
					Column col = new Column("Condition A");
					col.setHflex("min");
					col.setParent(cols);

					col = new Column("Condition B");
					col.setHflex("min");
					col.setParent(cols);
					
					col = new Column("Absolute Change");
					col.setHflex("min");
					col.setParent(cols);
					
					col = new Column("Relative Change");
					col.setHflex("min");
					col.setParent(cols);
					
					col = new Column("P-Value");
					col.setHflex("min");
					col.setParent(cols);
				
					cols.setParent(grid);
					
					Rows rows = new Rows();
					rows.setParent(grid);

					TreeSet<AlternativeSplicingExon> vcHits = m_vcASCurrentHit.GetResultsForRegion(nGroupStart, nGroupEnd);
					
					for(AlternativeSplicingExon ex : vcHits)
					{
						Row row = new Row();
						row.setParent(rows);
						
						Label label = new Label(ex.GetCondition(true));
						label.setParent(row);
						
						label = new Label(ex.GetCondition(false));
						label.setParent(row);
						
						label = new Label(String.format(Locale.ENGLISH, "%.3f", ex.GetAbsoluteChange()));
						label.setParent(row);
						
						label = new Label(String.format(Locale.ENGLISH, "%.3f", ex.GetRelativeChange()));
						label.setParent(row);
						
						label = new Label(String.format(Locale.ENGLISH, "%.3e", ex.GetPValue()));
						label.setParent(row);
					}
				}
			}
		});
	}

	private void DrawJunctionHeatmap() throws IOException
	{
		CreateCoverageMatrix();
		
		//###################################
		//       prepare popup window
		//###################################
		m_windowPopup.setParent(this);
		m_windowPopup.setWidth((m_nClientWindowWidth*0.5) +"px");
		m_windowPopup.setHeight((m_nClientWindowHeight*0.5) +"px");
		m_windowPopup.setMaximizable(true);
		m_windowPopup.setMinimizable(true);
		m_windowPopup.setPosition("center,center");
		m_windowPopup.setVisible(true);
		
		// clear old popup
		if(m_windowPopup != null)
		{
			Components.removeAllChildren(m_windowPopup);
		}

		//####################################################
		//           add options to options box
		//####################################################	
		Hlayout layout = new Hlayout();
		
		layout.setHflex("100%");
		layout.setHeight("100%");
		layout.setStyle("overflow:auto;");

		TreeMap<String, Double> mapSizeFactors	= m_projectModel.GetSizeFactors();
		TreeMap<String, String> mapBigWigFiles 	= m_projectModel.GetBigWigFilesForSamples();
		
		if(mapBigWigFiles == null)
		{
			System.out.println("ERROR: no valid bigwig or bam files detected");
			return;
		}
		
		// get samples per Condition
		TreeMap<String, TreeSet<String>> mapSamplesToConditions =  m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);

		// get reference name
		String strRef = m_Gene.getChromosome();
		
		//##########################################################
		//                 get junction counts
		//##########################################################
		TreeMap<String, TreeMap<String, Integer>> mapJunCountsToJunctions = m_projectModel.GetJunctionCountsForGene(m_Gene);
		
		// get maximum junction count value
		int nSamples = m_projectModel.GetSamples().size();
		double fMaxCount = 0.0;
		for(String strJunction : mapJunCountsToJunctions.keySet())
		{
			for(String strSample : mapJunCountsToJunctions.get(strJunction).keySet())
			{
				double fValue = mapJunCountsToJunctions.get(strJunction).get(strSample) * mapSizeFactors.get(strSample);
				fMaxCount = Math.max(fMaxCount, fValue);
			}
		}
		
		//##########################################################
		//          remove junctions with too low coverage
		//##########################################################
		boolean bModified = true;
		while(bModified)
		{
			bModified = false;
			
			for(String strJunction : mapJunCountsToJunctions.keySet())
			{
				double fRowMaxCount = 0;
				int nValidSamples 	= 0;
				for(String strSample : mapJunCountsToJunctions.get(strJunction).keySet())
				{
					double fValue = mapJunCountsToJunctions.get(strJunction).get(strSample) * mapSizeFactors.get(strSample);
					fRowMaxCount = Math.max(fRowMaxCount, fValue);
					
					if(fValue >= m_nMinJunctionReads)
						nValidSamples++;
				}
				
				if(fRowMaxCount == 0 || nValidSamples < 3)
				{
					bModified = true;
					mapJunCountsToJunctions.remove(strJunction);
					break;
				}
			}
		}
		
		//######################################################
		//                 prepare graph
		//######################################################
		int nNameOffset = 180; // subtract 100 px for the isoform names
		int nMargin = 20;
		int nSymbolWidth = 20;
		int nMaxWidth	= nMargin+30 + nSamples*20 + nNameOffset + nSymbolWidth;
		int nMaxHeight  = (mapJunCountsToJunctions.size()*40) + 30;
		
		// prepare graph
		BufferedImage img = new BufferedImage(nMaxWidth, nMaxHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D graph = img.createGraphics();
		
		if(m_nClientWindowWidth*0.5 > nMaxWidth)
			m_windowPopup.setWidth((nMaxWidth+20)+"px");
		
		if(m_nClientWindowHeight*0.5 > nMaxHeight)
			m_windowPopup.setHeight((nMaxHeight+60)+"px");
		
		// set font size
		Font fontNormalText = new Font("Courier", Font.PLAIN, 11);
		graph.setFont(fontNormalText);
		
		FontMetrics fontMetrics = graph.getFontMetrics(fontNormalText);
		
		// fill background
		graph.setColor(Color.white);
		graph.fillRect(0, 0, nMaxWidth, nMaxHeight);
		
		// reserve some space for the isoform name
		int nOffsetX = nMargin + nNameOffset;
		
		//###################################################
		// draw heatmap
		//###################################################
		int nY = 0;
		
		for(String strJunction : mapJunCountsToJunctions.keySet())
		{
			double fRowMaxCount = 0;
			for(String strSample : mapJunCountsToJunctions.get(strJunction).keySet())
			{
				double fValue = mapJunCountsToJunctions.get(strJunction).get(strSample) * mapSizeFactors.get(strSample);
				fRowMaxCount = Math.max(fRowMaxCount, fValue);
			}

			int nStart = Integer.parseInt(strJunction.split("_")[0])+1;
			int nEnd   = Integer.parseInt(strJunction.split("_")[1])-1;
			
			String strText = strRef + ":" + nStart + "-" + nEnd;
			graph.setColor(Color.BLACK);
 			graph.drawString(strText, nMargin, nY+fontMetrics.getHeight());
			
			int nX = nOffsetX;
			for(String strSample : mapJunCountsToJunctions.get(strJunction).keySet())
			{
				// determine color
				double fValue = mapJunCountsToJunctions.get(strJunction).get(strSample) * mapSizeFactors.get(strSample);
				
				Shape rect = new Rectangle2D.Double(nX, nY, 20, 20);
				graph.setColor(new Color((int)(fValue/fMaxCount*255.0), 0, 0));
				graph.fill(rect);
				
				nX += 20;				
			}
			nY += 20;
		}
		
		nY += 20;
		for(String strJunction : mapJunCountsToJunctions.keySet())
		{
			double fRowMaxCount = 0;
			for(String strSample : mapJunCountsToJunctions.get(strJunction).keySet())
			{
				double fValue = mapJunCountsToJunctions.get(strJunction).get(strSample) * mapSizeFactors.get(strSample);
				fRowMaxCount = Math.max(fRowMaxCount, fValue);
			}
			
			int nStart = Integer.parseInt(strJunction.split("_")[0])+1;
			int nEnd   = Integer.parseInt(strJunction.split("_")[1])-1;
			
			String strText = strRef + ":" + nStart + "-" + nEnd;
			graph.setColor(Color.BLACK);
			graph.drawString(strText, nMargin, nY+fontMetrics.getHeight());
			
			int nX = nOffsetX;
			Vector<double[]> vcCounts = new Vector<double[]>();
			
			for(String strCondition : mapSamplesToConditions.keySet())
			{
				int nCurSamples = mapSamplesToConditions.get(strCondition).size();
				double pCounts[] = new double[nCurSamples];
				
				int nIdx = 0;
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					// determine color
					double fValue = mapJunCountsToJunctions.get(strJunction).get(strSample) * mapSizeFactors.get(strSample);
					
					pCounts[nIdx] = fValue;
					
					float fRes = (float)(fValue / fRowMaxCount);
					
					Shape rect = new Rectangle2D.Double(nX, nY, 20, 20);
					graph.setColor(new Color((int)(fRes*255.0), 0, 0));
	//				graph.setColor(getHeatMapColor(fRes));
					
					graph.fill(rect);
					
					nX += 20;
					
					nIdx += 1;
				}
				
				// gap for each condition
				nX += 2;
				
				vcCounts.add(pCounts);
			}
			
			// calculate anova		
			OneWayAnova anova = new OneWayAnova();
			double fPValue = anova.anovaPValue(vcCounts);
			
			if(fPValue < 0.05)
			{
				String strPath = Executions.getCurrent().getDesktop().getWebApp().getRealPath("/");
				BufferedImage imgExclamation = ImageIO.read(new File(strPath + "/img/exclamation_sign.png"));
				BufferedImage resizedImg = new BufferedImage(20, 18, BufferedImage.TYPE_INT_ARGB);
				
				Graphics2D g = resizedImg.createGraphics();
				g.drawImage(imgExclamation, 0, 0, 20, 18, null);
				g.dispose();
				
				graph.drawImage(resizedImg, nX+5, nY, null);
			}
//			System.out.println(strJunction + ": " + fPValue);
			
			nY += 20;
		}
		
		layout.setParent(m_windowPopup);
		
		Imagemap imgMap = new Imagemap();
		imgMap.setWidth(nMaxWidth+"px");
		imgMap.setHeight((nMaxHeight+60)+"px");
		imgMap.setContent(img);
		imgMap.setParent(layout);
	}
	
	// calculates the exon coverage just using their posititon and the coverage in the bigwig files
	public void CalculateRawExonCoverage() throws IOException
	{		
		if(m_mapExonCoveragePerSampleRawFromBigWig == null)
			m_mapExonCoveragePerSampleRawFromBigWig = new TreeMap<Exon, TreeMap<String, Double>>();
		
		if(m_mapExonCoverageLengthPerSampleRawFromBigWig == null)
			m_mapExonCoverageLengthPerSampleRawFromBigWig = new TreeMap<Exon, TreeMap<String, Double>>();
		
		TreeMap<String, TreeSet<String>> mapSamples = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		for(String strCondition : mapSamples.keySet())
		{	
			// and each sample with the condition
			for(String strSample : mapSamples.get(strCondition))
			{
				if(m_vcSelectedSamples.contains(strSample))
				{
					// proceed through all exon groups
					Exon pExons[] = m_Gene.getSortedExons();
					
					ExecutorService executor = Executors.newFixedThreadPool(m_nThreads);
					
					// start a thread for each exon
					for(int i=0; i<pExons.length; i++)
					{
						Runnable r = new ThreadGetRawExonCoverage(pExons[i], strSample, true, m_mapExonCoveragePerSampleRawFromBigWig, m_mapExonCoverageLengthPerSampleRawFromBigWig);
						executor.execute(r);
					}
					
					executor.shutdown();
					
					while(!executor.isTerminated()) {}
				}
			}
		}
	}

	public boolean CalculateExonCoverage() throws IOException
	{
		TreeMap<String, TreeSet<String>> mapSamples = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);

		// check which exons are actually expressed
		TreeMap<Exon, TreeMap<String, Double>> mapExonCoveragePerSample   = new TreeMap<Exon, TreeMap<String, Double>>();
		TreeMap<Exon, TreeMap<String, TreeSet<Exon>>> mapIndistinguishableExonPerSample = new TreeMap<Exon, TreeMap<String, TreeSet<Exon>>>();
		
		// get junction read counts
		for(String strCondition : mapSamples.keySet())
		{	
			// and each sample with the condition
			for(String strSample : mapSamples.get(strCondition))
			{
				// only use selected samples
				if(m_vcSelectedSamples.contains(strSample))
				{
					ExecutorService executor = Executors.newFixedThreadPool(m_nThreads);
	
					// proceed through all exon groups
					for(ExonGroup grp : m_pExonGroups)
					{
						//TODO exclude invalid exons (based on intron coverage)
						Runnable r = new ThreadCalculateExonCoverage(strSample, grp, mapExonCoveragePerSample, mapIndistinguishableExonPerSample);
						executor.execute(r);
					}
					
					executor.shutdown();
					while(!executor.isTerminated()) {}
				}
			}
		}

		// fill missing values
		for(Exon ex : mapExonCoveragePerSample.keySet())
		{
			for(String strCondition : mapSamples.keySet())
			{	
				// and each sample with the condition
				for(String strSample : mapSamples.get(strCondition))
				{
					if(!mapExonCoveragePerSample.get(ex).containsKey(strSample))
					{
						mapExonCoveragePerSample.get(ex).put(strSample, 0.0);
					}
				}
			}
		}

		m_mapIndistinguishableExonPerSample = mapIndistinguishableExonPerSample;
		m_mapExonCoveragePerSample = mapExonCoveragePerSample;
		
		return true;
	}

	// Given an exon, get the most frequent set of ambigious exons (=values) per condition (=key) is returned
	public TreeMap<String, TreeSet<Exon>> GetMostFrequentSetOfAmbigiousExons(Exon exon)
	{
		TreeMap<String, TreeSet<Exon>> mapAmbigiousExonsPerCondition = new TreeMap<String, TreeSet<Exon>>();
		
		if(m_strSelectedConditionType == null)
		{
			Messagebox.show("No condition type selected");
			return null;
		}
		
		TreeMap<String, TreeSet<String>> vcSamplesPerCondition = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		for(String strCondition : vcSamplesPerCondition.keySet())
		{
			HashMap<TreeSet<Exon>, Integer> mapExonOccurrences = new HashMap<TreeSet<Exon>, Integer>(); 
			for(String strSample : vcSamplesPerCondition.get(strCondition))
			{
				// skip unselected samples
				if(!m_vcSelectedSamples.contains(strSample))
					continue;
				
				// proceed through all sets of indistinguishable exons
				for(Exon tarEx : m_mapIndistinguishableExonPerSample.keySet())
				{
//					if(m_mapIndistinguishableExonPerSample.get(tarEx).containsKey(strSample))
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
	
	//################################################################################################
	//    Calculates the expression of an exon relative to all other exons in the same exon group.
	//    Returns an map with the median relative coverage (value) per condition (key).
	//################################################################################################
	public TreeMap<String, Double> CalculateRelativeExonExpressionPerCondition(TreeMap<String, Double> mapRelativeExonCoveragePerSample)
	{
		TreeMap<String, TreeSet<String>> vcSamplesPerCondition = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		TreeMap<String, Double> mapRelativeExonCoveragePerCondition = new TreeMap<String, Double>();

		for(String strCondition : vcSamplesPerCondition.keySet())
		{
			Vector<Double> vcCoveragesPerSample = new Vector<Double>();

			for(String strSample : vcSamplesPerCondition.get(strCondition))
			{
				if(mapRelativeExonCoveragePerSample.containsKey(strSample))
				{
					vcCoveragesPerSample.add(mapRelativeExonCoveragePerSample.get(strSample));
				}
			}
			
			if(vcCoveragesPerSample.size() == 0)
			{
				mapRelativeExonCoveragePerCondition.put(strCondition, 0.0);
			}
			else
			{
				double[] pCoveragesPerSample = new double[vcCoveragesPerSample.size()];
				for(int i=0; i<vcCoveragesPerSample.size(); i++)
					pCoveragesPerSample[i] = vcCoveragesPerSample.get(i);
				
				double fMean = StatUtils.mean(pCoveragesPerSample);
				mapRelativeExonCoveragePerCondition.put(strCondition, fMean);
				
//				double fMedian = StatUtils.geometricMean(pCoveragesPerSample);
//				mapRelativeExonCoveragePerCondition.put(strCondition, fMedian);
			}
		}
		
		return mapRelativeExonCoveragePerCondition;
	}
	
	public TreeMap<String, Double> GetAbsoluteCoveragePerSample(Exon exon, ExonGroup exonGroup)
	{
		TreeMap<String, TreeSet<String>> vcSamplesPerCondition = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		TreeMap<String, Double> mapAbsoluteExonCoveragePerSample = new TreeMap<String, Double>();
		
		// Get a majority vote on how to calculate the (relative) expression of this exon
		// The most frequent exon separation will be used
		TreeMap<String, Exon> mapTargetExonsPerSample = new TreeMap<String, Exon>();
		
		TreeMap<String, TreeSet<Exon>> vcAmbigiousExonSetsPerCondition = GetMostFrequentSetOfAmbigiousExons(exon);

		for(String strCondition : vcSamplesPerCondition.keySet())
		{		
			TreeSet<Exon> vcTarget = vcAmbigiousExonSetsPerCondition.get(strCondition);
			
			// get the target exon per sample
			for(String strSample : vcSamplesPerCondition.get(strCondition))
			{
				// use only selected samples
				if(!m_vcSelectedSamples.contains(strSample))
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
				if(!m_vcSelectedSamples.contains(strSample))
					continue;
				
				if(!mapTargetExonsPerSample.containsKey(strSample))
					continue;
				
				Exon targetExon = mapTargetExonsPerSample.get(strSample);

				if(vcSamplesPerCondition.get(strCondition).contains(strSample))
				{
					if(m_mapExonCoveragePerSample.containsKey(targetExon))
					{
						mapAbsoluteExonCoveragePerSample.put(strSample, m_mapExonCoveragePerSample.get(targetExon).get(strSample));
					}
				}
			}
		}
		
		return mapAbsoluteExonCoveragePerSample;
	}
	
	public TreeMap<String, Double> CalculateRelativeExonExpressionPerSample(Exon exon, ExonGroup exonGroup)
	{	
		/*
		int nLogStart = 35236399;
		int nLogStop = 35236408;
		*/
		
		TreeMap<String, TreeSet<String>> vcSamplesPerCondition = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		TreeMap<String, Double> mapRelativeExonCoveragePerSample = new TreeMap<String, Double>();
		
		// Get a majority vote on how to calculate the relative expression of this exon
		// The most frequent exon separation will be used
		TreeMap<String, Exon> mapTargetExonsPerSample = new TreeMap<String, Exon>();
		
		TreeMap<String, TreeSet<Exon>> vcAmbigiousExonSetsPerCondition = GetMostFrequentSetOfAmbigiousExons(exon);
		
		for(String strCondition : vcSamplesPerCondition.keySet())
		{		
			TreeSet<Exon> vcTarget = vcAmbigiousExonSetsPerCondition.get(strCondition);
			
			// get the target exon per sample
			for(String strSample : vcSamplesPerCondition.get(strCondition))
			{				
				for(Exon tarEx : m_mapIndistinguishableExonPerSample.keySet())
				{
					/*
					if(exon.getCodingStart() == 6644468)
						System.out.println(m_mapIndistinguishableExonPerSample.get(tarEx));
					*/
					
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
				
//			mapTargetExonsPerCondition.put(strCondition, vcTarget);
		}

		/*
		if(exon.getCodingStart() == nLogStart && exon.getCodingStop() == nLogStop)
		{
			System.out.println("most frequent sets: ");
			for(String strCondition : mapTargetExonsPerSample.keySet())
			{
				System.out.println(strCondition + "\t" + mapTargetExonsPerSample.get(strCondition));
			}
		}
		*/
	
		for(String strCondition : vcSamplesPerCondition.keySet())
		{
			for(String strSample : vcSamplesPerCondition.get(strCondition))
			{
				if(!mapTargetExonsPerSample.containsKey(strSample))
					continue;
				
				Exon targetExon = mapTargetExonsPerSample.get(strSample);

				/*
				if(exon.getCodingStart() == nLogStart && exon.getCodingStop() == nLogStop)
				{
					System.out.println(strCondition + " target exon: " + targetExon.getCodingStart() + "-" + targetExon.getCodingStop());
				}
				*/
				
				if(vcSamplesPerCondition.get(strCondition).contains(strSample))
				{
					if(m_mapExonCoveragePerSample.containsKey(targetExon))
					{
						double fMaxCoverage = 0.0;
						double fTotalCoverage = 0.0;
						
						// get total coverage for the current exon group
						for(Exon ex : exonGroup.getExons())
						{
							if(m_mapExonCoveragePerSample.containsKey(ex))
							{
								fTotalCoverage += m_mapExonCoveragePerSample.get(ex).get(strSample);
								fMaxCoverage = Math.max(fMaxCoverage, m_mapExonCoveragePerSample.get(ex).get(strSample));
							}
						}
		
						if(fTotalCoverage > 0)
						{
							double fCoverage = 0.0;
							if(m_bUseRelativeExonCoverageMode1)
							{
								fCoverage = m_mapExonCoveragePerSample.get(targetExon).get(strSample) / fTotalCoverage * 100.0;
							}
							else if(m_bUseRelativeExonCoverageMode2)
							{
								fCoverage = m_mapExonCoveragePerSample.get(targetExon).get(strSample) / fMaxCoverage * 100.0;
							}
							mapRelativeExonCoveragePerSample.put(strSample, fCoverage);
							
							/*
							if(exon.getCodingStart() == nLogStart && exon.getCodingStop() == nLogStop)
								System.out.println(targetExon.getCodingStart() + "-" + targetExon.getCodingStop() + "\t" + exon.getCodingStart() + "-" + exon.getCodingStop() + "\t" + strSample + " " + fCoverage);
							*/
						}
						else
						{
							mapRelativeExonCoveragePerSample.put(strSample, 0.0);
						}
					}
				}
			}
		}

		/*
		if(exon.getCodingStart() == nLogStart && exon.getCodingStop() == nLogStop)
		{
			for(String strSample : mapRelativeExonCoveragePerSample.keySet())
				System.out.println(strSample + ": " + mapRelativeExonCoveragePerSample.get(strSample));
		}
		*/

		return mapRelativeExonCoveragePerSample;
	}

	public double[] GetCoverageForRegion(String strRef, int nStart, int nEnd, String strSample, boolean bAdjustForSizeFactor) throws IOException
	{		
		// check whether bigwig and bam files are available for the samples
		TreeMap<String, Double> mapSizeFactors		= m_projectModel.GetSizeFactors();
		
		int nLength = nEnd - nStart + 1;
		double pCoverage[] = new double[nLength];
		for(int i=0; i<pCoverage.length; i++)
			pCoverage[i] = 0.0;
		
		// get size factor
		double fSizeFactor = mapSizeFactors.get(strSample);
		
		int pBigWigCoverage[] = m_mapBigWigCoverageToSamples.get(strSample);
		
		for(int i=nStart; i<=nEnd; i++)
		{
			int nPos = i - m_Gene.getStart();
			int nIdx = i - nStart;
			
			pCoverage[nIdx] = pBigWigCoverage[nPos];
			
			if(bAdjustForSizeFactor)
				pCoverage[nIdx] *= fSizeFactor;
		}

		return pCoverage;
	}
	
	public double[] GetCoverageForExonGroup(ExonGroup grp, String strSample, boolean bAdjustForSizeFactor) throws IOException
	{
		//###################################
		//    get coverage for exon group
		//###################################
		String strRef = m_Gene.getChromosome();
		int nGrpStart	= grp.getGenomicStartOfGroup();
		int nGrpEnd		= grp.getGenomicStopOfGroup();
		
		return(GetCoverageForRegion(strRef, nGrpStart, nGrpEnd, strSample, bAdjustForSizeFactor));
	}
	
	public double[] GetCoverageForExon(Exon ex, String strSample, boolean bAdjustForSizeFactor) throws IOException
	{
		//###################################
		//    get coverage for single exon
		//###################################
		String strRef 	= m_Gene.getChromosome();
		int nStart		= ex.getCodingStart();
		int nEnd		= ex.getCodingStop();
		
		return(GetCoverageForRegion(strRef, nStart, nEnd, strSample, bAdjustForSizeFactor));
	}
	
	//#################################################################################
	//    Groups exons of an exon group if they differ by less than nMinUniqueBases
	//#################################################################################
	public TreeMap<Exon, TreeSet<Exon>> GroupSimilarExons(Vector<Exon> vcExons, int nMinUniqueBases)
	{
		TreeMap<Exon, TreeSet<Exon>> mapGroupedExons = new TreeMap<Exon, TreeSet<Exon>>(); 
		for(Exon ex : vcExons)
		{
			boolean bNewExon = true;
			// check whether the current exon has sufficient unique bases
			for(Exon exIncluded : mapGroupedExons.keySet())
			{
				int nExonStart	= ex.getCodingStart();
				int nExonEnd	= ex.getCodingStop();
				
				int nMetaExonStart 	= exIncluded.getCodingStart();
				int nMetaExonEnd 	= exIncluded.getCodingStop();
				
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
	
	//################################################################################
	//    Takes a set of exons and calculates the coverage per base for each exon.
	//    Requires a coverage array where the first value corresponds to the first
	//    base of the corresponding exon group to which the exon belongs.
	//################################################################################
	public HashMap<Exon, Double> CalculateExonCoverageGroupedByCoverage(TreeSet<Exon> vcExons, int nExonGroupStart, double pCoverage[])
	{
		HashMap<Exon, Double> mapCoveragePerExon = new HashMap<Exon, Double>();
		for(Exon ex : vcExons)
		{
			int nExStart = ex.getCodingStart();
			int nExStop  = ex.getCodingStop();
			
			double fCoverage = 0.0;
			
			for(int i = nExStart; i<=nExStop; i++)
			{
				int nPos = i - nExonGroupStart;
				
				fCoverage += pCoverage[nPos];
			}
			
			fCoverage /= ex.getLength();
			
			mapCoveragePerExon.put(ex, fCoverage);
		}
		
		return mapCoveragePerExon;
	}
	
	//####################################################################################################
	//    Merges all exons with similar expression where one exon is completely contained in the other
	//####################################################################################################
	public void MergeSimilarExpressedContainedExons(TreeMap<Exon, TreeSet<Exon>> mapGroupedExons, HashMap<Exon, Double> mapCoveragePerExon)
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

	//##############################################################################
	//    Takes a map of coverage values and exons to calculate which positions
	//    are blocked by which exon. To be more precise, beginning with the
	//    highest covered exon, all base positions covered by the exon are assigned
	//	  the respective exon. The next highest covered exon may then only claim
	//    those base positions that have not been claimed before. This results
	//    in a map where each position is assigned to a specific exon.
	//##############################################################################
	public TreeMap<Integer, TreeSet<Exon>> DefineExonSpecificCoverageRegions(HashMap<Exon, Double> mapCoverageToExon, int nExonGrpStart)
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
			
			int nExStart = firstExon.getCodingStart();
			int nExStop  = firstExon.getCodingStop();
			
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
	
	//###############################################################################
	//    Calculates the 'unique' coverage per exon. The first parameter is a map
	//    where each position of the coverage array (last parameter) is assigned
	//    to exactly one exon. For each exon unambigious exon (key of the grouped
	//    exon map) the unique coverage is calculated.
	//
	//    PLEASE BE AWARE that, in this context, 'unique' coverage may also
	//    originate from regions that are shared among multiple exons but were
	//    assigned to a specifc exon due to it's superior overall coverage.
	//###############################################################################
	public TreeMap<Exon, Double> CalculateUniqueCoveragePerExon(TreeMap<Integer, TreeSet<Exon>> mapCoveragePositions, TreeMap<Exon, TreeSet<Exon>> mapGroupedExons, int nExonGrpStart, double pCoverage[])
	{
		// for each exon, calculate the unique coverage
		TreeMap<Exon, Double> mapUniqueCoverage = new TreeMap<Exon, Double>();
		
		for(Exon ex : mapGroupedExons.keySet())
		{
			int nExStart = ex.getCodingStart();
			int nExEnd  = ex.getCodingStop();
			
			double fCoverage = 0.0;
			int nUniqueBases = 0;
									
			for(int i = nExStart; i<=nExEnd; i++)
			{
				int nPos = i - nExonGrpStart;								
				
				// only use this position if it hasn't been blocked by another transcript 
				if(mapCoveragePositions.containsKey(nPos))
				{
					if(mapCoveragePositions.get(nPos).contains(ex))
					{					
						fCoverage += pCoverage[nPos];
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
	
	//###################################################################################
	//    Calculates the final coverage per exon based on the 'unique' exon coverage.
	//    At first, all exons that have 'unique' coverage are assigned to the respective
	//    base positions in the coverage array (see above for the definition of 'unique'
	//    coverage). The coverage array starts with the first base of the exon group.
	//###################################################################################
	public TreeMap<Exon, Double> CalculateFinalCoveragePerExon(TreeMap<Exon, Double> mapUniqueCoverage, int nExonGroupStart, double pCoverage[])
	{
		TreeMap<Exon, Double> vcFinalExonCoverage = new TreeMap<Exon, Double>();
		TreeMap<Exon, Double> vcTmpCoverage = new TreeMap<Exon, Double>();
		
		// assign exons to positions in the coverage array
		TreeMap<Integer, TreeSet<Exon>> mapExonsToPositions = new TreeMap<Integer, TreeSet<Exon>>(); 
		for(Exon ex : mapUniqueCoverage.keySet())
		{
			int nExStart = ex.getCodingStart();
			int nExEnd   = ex.getCodingStop();
			
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
	
	public void CalculateJunctionCoverage() throws IOException
	{
		// clear old values
		m_mapAnalysedJunctionPaths.clear();
		m_mapCoveragePerJunctionPathAndSample.clear();
		
		TreeMap<String, TreeMap<String, Integer>> mapJunctionCounts = m_projectModel.GetJunctionCountsForGene(m_Gene);
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		TreeMap<String, int[]> mapSelectedIsoforms 	= GetSelectedIsoforms();
		
		// get first and second selected exon position
		ExonGroup grpA = null;
		ExonGroup grpB = null;

		for(ExonGroup grp : m_pExonGroups)
		{
			if(grp.getGenomicStartOfGroup() >= m_ClickEvent.m_nExonGroupStartA && grp.getGenomicStopOfGroup() <= m_ClickEvent.m_nExonGroupEndA)
			{
				grpA = grp;
			}
			
			if(grp.getGenomicStartOfGroup() >= m_ClickEvent.m_nExonGroupStartB && grp.getGenomicStopOfGroup() <= m_ClickEvent.m_nExonGroupEndB)
			{
				grpB = grp;
			}
		}
		
		if(grpA == null || grpB == null)
			return;
				
		// get all possible paths between the two selected exon groups: <path_id, included_junctions> 
		TreeMap<Integer, TreeSet<String>> mapPaths = new TreeMap<Integer, TreeSet<String>>();
		for(String strIsoform : mapSelectedIsoforms.keySet())
		{
			// get exons for the isoform
			int[] pExonIDs = mapSelectedIsoforms.get(strIsoform);
			
			int nRegionStart = Math.min(grpA.getGenomicStartOfGroup(), grpB.getGenomicStartOfGroup());
			int nRegionEnd	 = Math.max(grpA.getGenomicStopOfGroup(), grpB.getGenomicStopOfGroup());
			
			// get a list of all exons within the specified region
			TreeSet<Exon> vcExons = new TreeSet<Exon>();
			for(int nID : pExonIDs)
			{
				int nExonID = Math.abs(nID - (m_pExons.length-1));
				Exon ex = m_pExons[nExonID];
			
				if(ex.getCodingStart() >= nRegionStart && ex.getGenomicStop() <= nRegionEnd)
				{
					vcExons.add(ex);
				}
			}
			
			// define junctions
			TreeSet<String> vcJunctions = new TreeSet<String>();
			Exon exPrevious = null;
			boolean bValidPath = true;
			for(Exon ex : vcExons)
			{
				if(exPrevious != null)
				{
					String strJunction = exPrevious.getCodingStop() + "_" + ex.getCodingStart(); 
//					System.out.println(strIsoform + ": " + strJunction);
					vcJunctions.add(strJunction);
					
					// invalid path if one of the junctions is not expressed
					if(!mapJunctionCounts.containsKey(strJunction))
					{
						bValidPath = false;
						System.out.println("junction not detected: " + strJunction + " (" + m_Gene.getGeneName() + ")");
					}
				}
				exPrevious = ex;
			}
			
			if(bValidPath)
			{
				String strID = "path";
				for(String strJun : vcJunctions)
					strID += "_" + strJun;			

				TreeMap<Integer, TreeSet<String>> mapIsoformPath = new TreeMap<Integer, TreeSet<String>>();
				mapIsoformPath.put(strID.hashCode(), vcJunctions);
//				System.out.println(strIsoform + ": " + strID + " " + mapIsoformPath);
				m_mapAnalysedJunctionPaths.put(strIsoform, mapIsoformPath);

				mapPaths.put(strID.hashCode(), vcJunctions);
			}
		}
		
		// calculate coverage ratios of the paths
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			// <sample, <junction, ratio>>
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				// <path_hash_id, junction>
				HashMap<Integer, String> mapLowestCoveredJunction = new HashMap<Integer, String>();
				for(Map.Entry<Integer, TreeSet<String>> e : mapPaths.entrySet())
				{
					// get lowest covered junction per path
					int nMin			= Integer.MAX_VALUE;
					String strLowest 	= "";
					
					for(String strJun : e.getValue())
					{						
						int nCount = mapJunctionCounts.get(strJun).get(strSample);
						if(nCount < nMin)
						{
							nMin = nCount;
							strLowest = strJun;
						}
					}
					mapLowestCoveredJunction.put(e.getKey(), strLowest);
				}
				
				// count how many times the identified junction was the lowest covered junction
				TreeMap<String, Integer> mapCountLowest = new TreeMap<String, Integer>();
				for(Map.Entry<Integer, String> e : mapLowestCoveredJunction.entrySet())
				{
					if(mapCountLowest.containsKey(e.getValue()))
					{
						mapCountLowest.put(e.getValue(), mapCountLowest.get(e.getValue())+1);
					}
					else
					{
						mapCountLowest.put(e.getValue(), 1);
					}
				}
				
				// get lowest junction coverage per path, adjusted for the number of times the junction occurs
				double fReadCountSum = 0;
				TreeMap<Integer, Double> mapRecalculatedCounts = new TreeMap<Integer, Double>();
				for(Map.Entry<Integer, TreeSet<String>> e : mapPaths.entrySet())
				{
					// get lowest covered junction of the path
					String strLowestJunction = mapLowestCoveredJunction.get(e.getKey());

					// how many times was the lowest covered junction used?
					int nCount = mapCountLowest.get(strLowestJunction);
					
					// get read count of lowest covered junction
					int nReadCount = mapJunctionCounts.get(strLowestJunction).get(strSample);
					
					// recalculate count
					double fValue = (double)nReadCount / (double)nCount;
					mapRecalculatedCounts.put(e.getKey(), fValue);
					
//					System.out.println(strSample + " lowest jun: " + strLowestJunction + "  | " + nCount + " | " + nReadCount + " | " + fValue);
					
					fReadCountSum += fValue;
				}
				
				// calculate ratios
				for(Map.Entry<Integer, TreeSet<String>> e : mapPaths.entrySet())
				{
					double fRatio = 0.0;
					if(fReadCountSum != 0)
						fRatio = mapRecalculatedCounts.get(e.getKey()) / fReadCountSum;
					
//					System.out.println(strSample + " " + e.getKey() + " " + fRatio);
					
					if(m_mapCoveragePerJunctionPathAndSample.containsKey(strSample))
					{
						TreeMap<Integer, Double> mapValue = m_mapCoveragePerJunctionPathAndSample.get(strSample); 
						mapValue.put(e.getKey(), fRatio);
						m_mapCoveragePerJunctionPathAndSample.put(strSample, mapValue);
					}
					else
					{
						TreeMap<Integer, Double> mapValue = new TreeMap<Integer, Double>(); 
						mapValue.put(e.getKey(), fRatio);
						m_mapCoveragePerJunctionPathAndSample.put(strSample, mapValue);
					}
				}
			}
		}
	}
	
	public TreeMap<String, Double> CalculateRelativeJunctionExpressionPerCondition(int nJunctionPathID)
	{
		TreeMap<String, Double> mapRelativeJunctionCoveragePerCondition = new TreeMap<String, Double>();
		TreeMap<String, TreeSet<String>> mapSamplesPerCondition = m_projectModel.GetSamplesPerCondition(this.m_strSelectedConditionType);
		
		for(String strCondition : mapSamplesPerCondition.keySet())
		{
			double pValues[] = new double[mapSamplesPerCondition.get(strCondition).size()];
			
			int nSample = 0;
			for(String strSample : mapSamplesPerCondition.get(strCondition))
			{
				pValues[nSample] = m_mapCoveragePerJunctionPathAndSample.get(strSample).get(nJunctionPathID);
				nSample += 1;
			}
			
			double fMedian = StatUtils.percentile(pValues, 50);
			mapRelativeJunctionCoveragePerCondition.put(strCondition, fMedian);
		}
		
		return mapRelativeJunctionCoveragePerCondition;
	}
	
	public TreeMap<String, int[]> GetSelectedIsoforms()
	{
		// <isoform_name, exon_ids>
		TreeMap<String, int[]> mapValidIsoforms = new TreeMap<String, int[]>();
		
		// get first and second selected exon position
		ExonGroup grpA = null;
		ExonGroup grpB = null;

		if(m_pExonGroups != null)
		{
			for(ExonGroup grp : m_pExonGroups)
			{
				if(grp.getGenomicStartOfGroup() >= m_ClickEvent.m_nExonGroupStartA && grp.getGenomicStopOfGroup() <= m_ClickEvent.m_nExonGroupEndA)
				{
					grpA = grp;
				}
				
				if(grp.getGenomicStartOfGroup() >= m_ClickEvent.m_nExonGroupStartB && grp.getGenomicStopOfGroup() <= m_ClickEvent.m_nExonGroupEndB)
				{
					grpB = grp;
				}
			}
		}
		
		boolean bRequiresBothExons = false;
		if(grpB != null)
			bRequiresBothExons = true;

		if(grpA != null && (grpB != null || !bRequiresBothExons))
		{
			HashMap<String, int[]> mapAllIsoforms   = new HashMap<String, int[]>();
			
			Iterator<String> it = m_Gene.getGeneProductNames();
			while(it.hasNext())
			{
				String strGeneProductID = it.next();
				int[] pExonIDs = m_Gene.getExonsForGeneProduct(strGeneProductID);
				
				mapAllIsoforms.put(strGeneProductID, pExonIDs);
			}			

			it = m_Gene.getGeneProductNames();
			while(it.hasNext())
			{
				String strIsoform = it.next();		
				Exon pExons[] = m_Gene.getSortedExonsForGeneProduct(strIsoform, true);
										
				// get exon group for each exon
				boolean bIncludesExonA = false;
				boolean bIncludesExonB = false;
				
				for(Exon ex : pExons)
				{
					if(ex.getGenomicStart() >= grpA.getGenomicStartOfGroup() && ex.getGenomicStop() <= grpA.getGenomicStopOfGroup())
						bIncludesExonA = true;
					
					if(grpB != null && ex.getGenomicStart() >= grpB.getGenomicStartOfGroup() && ex.getGenomicStop() <= grpB.getGenomicStopOfGroup())
						bIncludesExonB = true;

					if(bIncludesExonA && (bIncludesExonB || !bRequiresBothExons))
						break;
				}
				
				if(bIncludesExonA && (bIncludesExonB || !bRequiresBothExons))
				{
					mapValidIsoforms.put(strIsoform, m_Gene.getExonsForGeneProduct(strIsoform));
				}
			}
		}
		return mapValidIsoforms;
	}
	
	public boolean ProcessCountData() throws IOException
	{		
		// exons must be recalculated afterwards
		if(!CalculateExonCoverage()) return false;
		
		CalculateJunctionCoverage();
		
		if(m_ClickEvent.m_nExonGroupStartA != -1 && m_ClickEvent.m_nExonGroupStartB != -1)
		{
			CalculateJunctionCoverage();
		}
		
		return true;
	}

	public boolean DrawPlots() throws IOException, SQLException
	{
		if(m_Gene == null)
		{
			Messagebox.show("No gene specified");
			return false;
		}
		
		// clear old plots
		m_layoutWest.getChildren().clear();
		m_layoutMainBottom.getChildren().clear();		
		
		if(m_bCoverageRequiresRedraw)
		{
			m_layoutMainTop.getChildren().clear();
			DrawExtendedCoveragePlot();
		}
		
//		DrawSplicingGraph();
		DrawIsoforms();
	
		return true;
	}

	public boolean ReadGeneIdentifier(String strFile) throws FileNotFoundException
	{
		m_vcGeneIdentifier.clear();
		
		File pFile = new File(strFile);
		
		if(!pFile.exists())
		{
			System.out.println("ERROR: File does not exist: " + strFile);
			return false;
		}
		
		Scanner pScanner = new Scanner(pFile);
		while(pScanner.hasNextLine())
		{
			String strLine = pScanner.nextLine();

			GeneIdentifier gid = GeneIdentifier.ParseFromLineBioMart(strLine);
			if(gid != null && !gid.m_strApprovedGeneSymbol.contains("withdrawn"))
			{
				m_vcGeneIdentifier.add(gid);
			}
		}
		
		pScanner.close();
		
		return true;
	}

	public GeneIdentifier GetGeneIdentifierForGene(String strGene)
	{
		if(strGene == null)
			return null;
		
		for(GeneIdentifier gid : m_vcGeneIdentifier)
		{
			if(gid.EqualsGene(strGene) && gid.m_bIsValid)
			{				
				return gid;
			}
		}
		
		// nothing found yet? try to remove and dots (.) from the gene identifier
		for(GeneIdentifier gid : m_vcGeneIdentifier)
		{
			if(gid.EqualsGene(strGene.split("\\.")[0]) && gid.m_bIsValid)
			{
				return gid;
			}
		}
		
		System.out.println("failed to get gene identifier for: " + strGene);
		return null;
	}

	public TreeMap<String, Vector<Double>> GetGTEXDataForGene(String strGene) throws FileNotFoundException
	{
		strGene = strGene.split("\\.")[0];
		GeneIdentifier id = GetGeneIdentifierForGene(strGene);
		if(id == null)
		{
			Messagebox.show("failed to get gene identifier for gene: " + strGene);
			return null;
		}
		String strEntrezID = id.m_strEntrezGeneID;

		File pFolder = new File(m_strPathGeneGTEX);
		if(!pFolder.exists())
		{
			Messagebox.show("ERROR: Invalid GTEX folder specified");
			return null;
		}

		String strFile = m_strPathGeneGTEX + "/GTEX_" + strEntrezID + ".tsv";
		File pFile = new File(strFile);
		
		if(!pFile.exists())
		{
			Messagebox.show("ERROR: No GTEX data available for gene: " + id);
			return null;
		}
		
		// store counts per condition <tissue, counts> 
		TreeMap<String, Vector<Double>> mapCountsToTissues = new TreeMap<String, Vector<Double>>(); 
		
		Scanner pScanner = new Scanner(pFile);
		
		// skip first line (header)
		if(pScanner.hasNextLine())
			pScanner.nextLine();
		
		while(pScanner.hasNextLine())
		{
			String strLine = pScanner.nextLine();
			
			String split[] = strLine.split("\t");
			
			String strTissue = split[1];
			double nCount 	 = Double.parseDouble(split[2]);
			
			if(mapCountsToTissues.containsKey(strTissue))
			{
				Vector<Double> vcCounts = mapCountsToTissues.get(strTissue);
				vcCounts.add(nCount);
			}
			else
			{
				Vector<Double> vcCounts = new Vector<Double>();
				vcCounts.add(nCount);
				mapCountsToTissues.put(strTissue, vcCounts);
			}
		}
		
		pScanner.close();
		
		return mapCountsToTissues;
	}
	
	public void DrawGTEXBoxPlot(final TreeMap<String, Vector<Double>> mapCountsPerTissue)
	{
		int nLeftMargin = 80;
		int nMaxWidth	= 1000+nLeftMargin;
		int nMaxHeight  = 460;	// 400 for the track, 20 as border
		int nYAxisPos	= 60;
		
		//######################################
		//    get number of different organs
		//######################################
		TreeMap<String, TreeSet<String>> mapTissuesToOrgans = new TreeMap<String, TreeSet<String>>();
		for(String strTissue : mapCountsPerTissue.keySet())
		{
			String split[] = strTissue.split(" - ");
			String strOrgan = split[0].trim();
			
			if(mapTissuesToOrgans.containsKey(strOrgan))
			{
				TreeSet<String> vcTissues = mapTissuesToOrgans.get(strOrgan);
				if(split.length > 1)
					vcTissues.add(split[1]);
				else
					vcTissues.add(split[0]);
			}
			else
			{
				TreeSet<String> vcTissues = new TreeSet<String>();
				if(split.length > 1)
					vcTissues.add(split[1]);
				
				mapTissuesToOrgans.put(strOrgan, vcTissues);
			}
		}
		
		//#################################
		//    prepare tissue categories
		//#################################
		final TreeMap<String, String> mapCategoryToTissue = new TreeMap<String, String>();
		int nCategory = 1;
		for(String strOrgan : mapTissuesToOrgans.keySet())
		{
			int nNum = mapTissuesToOrgans.get(strOrgan).size();
			String strCategory = ""+nCategory;
			
			if(nNum > 0)
			{
				strCategory += "-" + (nCategory+nNum-1);
				mapCategoryToTissue.put(strCategory, strOrgan);
				
				for(String strTissue : mapTissuesToOrgans.get(strOrgan))
				{
					mapCategoryToTissue.put(""+nCategory, strOrgan + " - " + strTissue);
					nCategory +=1;
				}
			}
			else
			{
				mapCategoryToTissue.put(strCategory, strOrgan);
				nCategory +=1;
			}
		}
		
		// adjust selected checkboxes on first run
		if(m_vcPreviouslySelectedGTEXTissues == null)
		{
			m_vcPreviouslySelectedGTEXTissues = new TreeSet<String>();

			for(Map.Entry<String, String> e : mapCategoryToTissue.entrySet())
			{
				m_vcPreviouslySelectedGTEXTissues.add(e.getKey());
			}
		}
		
		//####################################
		//    get list of selected tissues
		//####################################
		TreeSet<String> vcSelectedTissues = new TreeSet<String>();
		TreeSet<String> vcSelectedOrgans  = new TreeSet<String>();

		for(String strOrgan : mapTissuesToOrgans.keySet())
		{
			// get organ category
			String strCategory = "?";
			for(Map.Entry<String, String> e : mapCategoryToTissue.entrySet())
			{
				if(e.getValue().equals(strOrgan))
				{
					strCategory = e.getKey();
					break;
				}
			}
			
			if(m_vcPreviouslySelectedGTEXTissues.contains(strCategory))
				vcSelectedOrgans.add(strOrgan);
			
			if(mapTissuesToOrgans.get(strOrgan).size() == 0)
			{
				if(m_vcPreviouslySelectedGTEXTissues.contains(strCategory))
					vcSelectedTissues.add(strOrgan);
			}
			else
			{
				for(String strTissue : mapTissuesToOrgans.get(strOrgan))
				{
					strCategory = "?";
					for(Map.Entry<String, String> e : mapCategoryToTissue.entrySet())
					{
						if(e.getValue().equals(strOrgan + " - " + strTissue))
						{
							strCategory = e.getKey();
							break;
						}
					}
	
					if(m_vcPreviouslySelectedGTEXTissues.contains(strCategory))
					{
						vcSelectedTissues.add(strOrgan + " - " + strTissue);
					}
				}
			}
		}
		
		//###################################
		//          prepare graph
		//###################################
		BufferedImage img = new BufferedImage(nMaxWidth, nMaxHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D graph = img.createGraphics();
		
		graph.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		graph.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		// set font size
		int nFontSize    = 10;
		Font font = new Font("Arial", Font.BOLD, nFontSize);
		graph.setFont(font);
		
		// fill background
		graph.setColor(Color.white);
		graph.fillRect(0, 0, nMaxWidth, nMaxHeight);
		
		int nTissues = vcSelectedTissues.size();//mapCountsPerTissue.keySet().size();
		int nWidthPerTissue = ((nMaxWidth-nLeftMargin) - 5*nTissues) / nTissues;	// use 20 pixel as spacer between tissues and 40 pixels for the border/axis of the graph 

		// get a unique color per organ
		TreeMap<String, Color> mapOrganColors = new TreeMap<String, Color>();
		int nOrgan = 0;
		for(String strOrgan : mapTissuesToOrgans.keySet())
		{
			if(nOrgan < 30)
				mapOrganColors.put(strOrgan, m_pColors[nOrgan]);
			else
				mapOrganColors.put(strOrgan, new Color((int)(Math.random()%255), (int)(Math.random()%255), (int)(Math.random()%255)));
			nOrgan+=1;
		}
		
		// get maximum value
		double fMax = 0.0;
		for(String strTissue : vcSelectedTissues)
		{
			if(!mapCountsPerTissue.containsKey(strTissue))
			{
				System.out.println(vcSelectedTissues);
				System.out.println(strTissue);
				System.out.println(mapCountsPerTissue);
			}
			Vector<Double> vcValues = mapCountsPerTissue.get(strTissue);
			for(double fVal : vcValues)
			{
				if(fVal > fMax)
				{
					fMax = fVal;
				}
			}
		}
		
		graph.setColor(Color.BLACK);
		
		// draw axis
		int nX = nYAxisPos;
		graph.drawLine(nYAxisPos,  20, nYAxisPos,      420);
		graph.drawLine(nYAxisPos, 420, nYAxisPos+1000, 420);
		
		// add tick marks
		nX = nYAxisPos + 12;
		
		BasicStroke strokeBasic = new BasicStroke(1);
		BasicStroke strokeBold 	= new BasicStroke(2);
		
		// x-axis
		graph.setStroke(strokeBold);
		for(String strOrgan : vcSelectedOrgans)
		{
			int nNum = mapTissuesToOrgans.get(strOrgan).size();
			if(nNum == 0)
				nNum = 1;
			int nWidth  = (int)(nNum * (5+nWidthPerTissue)-5);
			
			graph.drawLine(nX+2, 425, nX+nWidth-2, 425);
			nX += nWidth+5;
		}
		
		graph.setStroke(strokeBasic);
		nX = nYAxisPos+12;
		for(String strTissue : vcSelectedTissues)
		{
			graph.drawLine((int)(nX+nWidthPerTissue*0.5), 425, (int)(nX+nWidthPerTissue*0.5), 427);
			
			String strCategory = "?";
			for(Map.Entry<String, String> e : mapCategoryToTissue.entrySet())
			{
				if(e.getValue().equals(strTissue))
				{
					strCategory = e.getKey();
					break;
				}
			}
			
			FontMetrics fontMetrics = graph.getFontMetrics(font);
			int nStringWidth = fontMetrics.stringWidth(strCategory);
			int nStringHeight = fontMetrics.getHeight();
			
			graph.drawString(strCategory, (int)(nX+nWidthPerTissue*0.5-nStringWidth*0.5), 426+nStringHeight);
			
			nX += 5 + nWidthPerTissue;
		}
		
		// y-axis
		int nStep = (int)Math.ceil(fMax / 20.0);	// results in ~20 ticks
		
		if(nStep > 9999) nStep = nStep - (nStep%1000);
		else if(nStep > 999) nStep = nStep - (nStep%100);
		else if(nStep > 10) nStep = nStep - (nStep%10);

		int nVal = 0;
		while(nVal < fMax)
		{
			int nYPos = (int) Math.ceil(420.0 - (nVal / fMax * 400.0));
			graph.drawLine(nYAxisPos-2, nYPos, nYAxisPos, nYPos);			
			
			String strValue = ""+nVal;
			FontMetrics fontMetrics = graph.getFontMetrics(font);
			int nStringWidth = fontMetrics.stringWidth(strValue);
			int nStringHeight = fontMetrics.getHeight();
			
			graph.drawString(strValue, (int)(nYAxisPos-3-nStringWidth), (int)(nYPos + nStringHeight*0.25));
			
			nVal += nStep;
		}

		nX = nYAxisPos+12;
		for(String strTissue : vcSelectedTissues)
		{
			Vector<Double> vcValues = mapCountsPerTissue.get(strTissue); 
			double pValues[] = new double[vcValues.size()];
			
			int nIdx = 0;
			for(double fVal : vcValues)
			{
				pValues[nIdx] = fVal;
				nIdx++;
			}
			
			double fQ25 = StatUtils.percentile(pValues, 25.0f);
			double fQ50 = StatUtils.percentile(pValues, 50.0f);
			double fQ75 = StatUtils.percentile(pValues, 75.0f);
			
			double fIQR 	= fQ75 - fQ25;
			double fMinQR 	= fQ25 -1.5 * fIQR;
			double fMaxQR 	= fQ75 +1.5 * fIQR;
			
			double fMinWhisker = Double.MAX_VALUE;
			double fMaxWhisker = Double.MIN_VALUE;
			
			Vector<Double> vcOutliers = new Vector<Double>();
			
			for(double fVal : pValues)
			{
				if(fVal >= fMinQR && fVal < fMinWhisker)
					fMinWhisker = fVal;
				
				if(fVal <= fMaxQR && fVal > fMaxWhisker)
					fMaxWhisker = fVal;
				
				// add outlier
				if(fVal > fMaxQR) vcOutliers.add(fVal);
				if(fVal < fMinQR) vcOutliers.add(fVal);
			}
			
			// adjust for maximum value
			fQ25 = Math.ceil(420.0 - (fQ25 / fMax * 400.0));
			fQ50 = Math.ceil(420.0 - (fQ50 / fMax * 400.0));
			fQ75 = Math.ceil(420.0 - (fQ75 / fMax * 400.0));
			fMinWhisker = Math.ceil(420.0 - (fMinWhisker / fMax * 400.0));
			fMaxWhisker = Math.ceil(420.0 - (fMaxWhisker / fMax * 400.0));
			
			// draw whishers
			graph.setColor(Color.BLACK);
			graph.drawLine((int)(nX+nWidthPerTissue*0.5), (int)fMinWhisker, (int)(nX+nWidthPerTissue*0.5), (int)fQ25);
			graph.drawLine((int)(nX+nWidthPerTissue*0.5), (int)fQ75, (int)(nX+nWidthPerTissue*0.5), (int)fMaxWhisker);
			
			String split[] = strTissue.split("-");
			
			// draw box
			if(fQ25 - fQ75 == 0.0)
			{
				graph.fillOval((int)(nX+nWidthPerTissue*0.5-2), 418, 4, 4);
				graph.drawOval((int)(nX+nWidthPerTissue*0.5-2), 418, 4, 4);				
			}
			else
			{				
				graph.setColor(mapOrganColors.get(split[0].trim()));
				graph.fillRect(nX, (int)fQ75, nWidthPerTissue, (int)(fQ50-fQ75));
				graph.fillRect(nX, (int)fQ50, nWidthPerTissue, (int)(fQ25-fQ50));
				
				graph.setColor(Color.BLACK);
				graph.drawRect(nX, (int)fQ75, nWidthPerTissue, (int)(fQ50-fQ75));
				graph.drawRect(nX, (int)fQ50, nWidthPerTissue, (int)(fQ25-fQ50));
			}
						
			// draw outlier
			for(double fVal : vcOutliers)
			{
				fVal = 420.0 - (fVal / fMax * 400.0);
				
				graph.setColor(mapOrganColors.get(split[0].trim()));
				graph.fillOval((int)(nX+nWidthPerTissue*0.5-2), (int)fVal-2, 4, 4);
				
				graph.setColor(Color.BLACK);
				graph.drawOval((int)(nX+nWidthPerTissue*0.5-2), (int)fVal-2, 4, 4);
			}
			
			nX += 5 + nWidthPerTissue;
		}
		
		//###################################
		//       prepare popup window
		//###################################
		m_windowPopup.setParent(this);
		m_windowPopup.setVflex("0");
		m_windowPopup.setHflex("0");

		m_windowPopup.setPosition("left,top");
		m_windowPopup.setVisible(true);
		
		// clear old popup
		if(m_windowPopup != null)
		{
			Components.removeAllChildren(m_windowPopup);
		}

		//####################################################
		//           add options to options box
		//####################################################	
		Borderlayout layout = new Borderlayout();
		West regionPlot 	= new West();
		East regionOptions 	= new East();
		Vlayout	layoutOptions = new Vlayout();
		
		// add buttons to option layout
		layoutOptions.setParent(regionOptions);
		layoutOptions.setWidth("200px");
		layoutOptions.setHflex("true");
		
		Button btnRedraw = new Button("Redraw");
		btnRedraw.setWidth("100%");
		btnRedraw.setParent(layoutOptions);
		btnRedraw.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				DrawGTEXBoxPlot(mapCountsPerTissue);
			}
		});
		
		//####################################################
		//    add tree for sample selection to options box
		//####################################################
		final Tree tree = new Tree();
		tree.setParent(layoutOptions);
		tree.setMultiple(true);
		tree.setCheckmark(true);
		tree.setSizedByContent(true);
		tree.setAutopaging(true);
		
		layout.setVflex("true");
		tree.setVflex("true");
		tree.setHflex("true");
		tree.setHeight("480px");
		
		regionOptions.setHflex("true");
		regionOptions.setVflex("true");
		
		regionPlot.setHflex("true");
		regionPlot.setVflex("true");
		
		// tree columns
		Treecols columns = new Treecols();
		columns.setParent(tree);
		Treecol column = new Treecol("Category");
		column.setParent(columns);
		column = new Treecol("Organ/Tissue");
		column.setParent(columns);
		
		// add organ nodes
		Treechildren childrenOrgan = new Treechildren();
		childrenOrgan.setParent(tree);

		//#################################
		//    add tissue selection tree
		//#################################
		
		for(String strOrgan : mapTissuesToOrgans.keySet())
		{
			int nNum = mapTissuesToOrgans.get(strOrgan).size();
			
			Treeitem itemOrgan = new Treeitem();
			itemOrgan.setParent(childrenOrgan);

			// get category
			String strCategory = "?";
			for(Map.Entry<String, String> e : mapCategoryToTissue.entrySet())
			{
				if(e.getValue().equals(strOrgan))
					strCategory = e.getKey();
			}
			
			itemOrgan.setOpen(false);
			if(m_vcPreviouslySelectedGTEXTissues.contains(strCategory))
			{
				tree.addItemToSelection(itemOrgan);
			}
			
			// add new organ row
			Treerow rowOrgan = new Treerow(strCategory);
			rowOrgan.setParent(itemOrgan);

			Treecell cellOrgan = new Treecell(strOrgan);
			cellOrgan.setParent(rowOrgan);
			
			cellOrgan = new Treecell(strOrgan);
			cellOrgan.setParent(rowOrgan);

			// add tissues
			if(nNum > 0)
			{
				Treechildren nodeTissue = new Treechildren();
				nodeTissue.setParent(itemOrgan);

				for(String strTissue : mapTissuesToOrgans.get(strOrgan))
				{					
					Treeitem itemTissue = new Treeitem();
					itemTissue.setParent(nodeTissue);
					
					// add row for tissue
					Treerow rowTissue = new Treerow();
					rowTissue.setParent(itemTissue);
					
					strCategory = "?";
					for(Map.Entry<String, String> e : mapCategoryToTissue.entrySet())
					{
						if(e.getValue().equals(strOrgan + " - " + strTissue))
							strCategory = e.getKey();
					}
					
					Treecell cellTissue = new Treecell(strCategory);
					cellTissue.setParent(rowTissue);
					cellTissue = new Treecell(strTissue);
					cellTissue.setParent(rowTissue);
					
					if(m_vcPreviouslySelectedGTEXTissues.contains(strCategory))
					{
						tree.addItemToSelection(itemTissue);
					}
				}
			}
		}

		// add event listener
		tree.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Tree clickedTree = (Tree)event.getTarget();

				TreeSet<String> vcNewSelectedTreeItems = new TreeSet<String>();
				for(Treeitem item : clickedTree.getSelectedItems())
					vcNewSelectedTreeItems.add(item.getLabel());
				
				String strToggledItem = null;
				boolean bSelected = false;	// selected or unselected

				//##################################################
				//    get the item that was clicked in this call
				//##################################################
				if(vcNewSelectedTreeItems.size() == 0)
				{
					// do nothing
					return;
				}
				
				if(vcNewSelectedTreeItems.size() > m_vcPreviouslySelectedGTEXTissues.size())
				{
					if(!m_vcPreviouslySelectedGTEXTissues.isEmpty())
					{
						vcNewSelectedTreeItems.removeAll(m_vcPreviouslySelectedGTEXTissues); 
					}
					bSelected = true;
					
					strToggledItem = vcNewSelectedTreeItems.first();
				}
				else
				{
					m_vcPreviouslySelectedGTEXTissues.removeAll(vcNewSelectedTreeItems);
					
					if(m_vcPreviouslySelectedGTEXTissues.size() > 0)
						strToggledItem = m_vcPreviouslySelectedGTEXTissues.first();
					else
						return;
				}

				// update previously selected items
				m_vcPreviouslySelectedGTEXTissues.clear();
				for(Treeitem item : clickedTree.getSelectedItems())
					m_vcPreviouslySelectedGTEXTissues.add(item.getLabel());

				// check whether a condition was checked
				if(mapCategoryToTissue.containsKey(strToggledItem))
				{					
					String split[] = strToggledItem.split("-");
					int nFirst = Integer.parseInt(split[0]);
					int nLast = -1;
					
					if(split.length > 1)
					{
						nLast = Integer.parseInt(split[1]); 
					}
					
					// select/unselect all samples for the current condition
					for(Treeitem item : tree.getItems())
					{					
						if(!item.getLabel().contains("-"))
						{
							String strIdx = item.getLabel();
							int nIdx = Integer.parseInt(strIdx);
	
							if(nIdx >= nFirst &&  nIdx <= nLast || nIdx == nFirst && nLast == -1)
							{
								if(bSelected)
								{
									tree.addItemToSelection(item);
									m_vcPreviouslySelectedGTEXTissues.add(item.getLabel());
								}
								else
								{
									tree.removeItemFromSelection(item);
									m_vcPreviouslySelectedGTEXTissues.remove(item.getLabel());
								}
							}
						}
					}
				}
			}
		});

		layout.setParent(m_windowPopup);
		m_windowPopup.setHeight((nMaxHeight + 100) + "px");
		m_windowPopup.setWidth ((nMaxWidth  + 220) + "px");
		m_windowPopup.setSizable(false);

		regionPlot.setParent(layout);
		regionOptions.setParent(layout);
		
		Imagemap imgMap = new Imagemap();
		imgMap.setWidth(nMaxWidth+"px");
		imgMap.setHeight((nMaxHeight+60)+"px");
		imgMap.setContent(img);
		imgMap.setParent(regionPlot);
	}
	
	public void ShowGTEXDataForGene(String strGene) throws FileNotFoundException
	{
		TreeMap<String, Vector<Double>> mapCountsPerTissue = GetGTEXDataForGene(strGene);
		if(mapCountsPerTissue != null)
		{
			DrawGTEXBoxPlot(mapCountsPerTissue);
		}
	}
	
	public void ShowGTEXDataForExon(String strExonicPart) throws FileNotFoundException
	{
		String pSplit[] = strExonicPart.split(" ");
		TreeMap<String, Vector<Double>> mapCountsPerTissue = m_mapCountsPerExonAndTissue.get(pSplit[1]);
		if(mapCountsPerTissue != null)
		{
			DrawGTEXBoxPlot(mapCountsPerTissue);
		}
	}
	
	public TreeMap<String, TreeMap<String, Vector<Double>>> GetGTEXDataForExonicParts() throws FileNotFoundException
	{
		TreeMap<String, TreeMap<String, Vector<Double>>> mapCountsPerExonicPartAndTissue = new TreeMap<String, TreeMap<String, Vector<Double>>>();
	
		GeneIdentifier id = GetGeneIdentifierForGene(m_Gene.getGeneID().split("\\.")[0]);
		if(id == null)
		{
			Messagebox.show("failed to get gene identifier for gene: " + m_Gene.getGeneID());
			return null;
		}
		
		String strEntrezID = id.m_strEntrezGeneID;
		
		String strGTEXID = "?";

		// open file with DEXSeq exon IDs
		File pFileDEXSeqIDs = new File(m_strFileDEXSeqIDs);
		Scanner pScanner = new Scanner(pFileDEXSeqIDs);
		
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
				
//				if(pSplit[8].equals("gene_id \"" + strEntrezID +"\""))
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
		File pFolder = new File(m_strPathExonGTEX);
		if(!pFolder.exists())
		{
			Messagebox.show("ERROR: Invalid GTEX folder specified: " + m_strPathExonGTEX);
			return null;
		}
		
		for(String strPos : mapIDsToPositions.keySet())
		{
			String strFile = m_strPathExonGTEX + "//GTEX_" + strGTEXID + "_" + mapIDsToPositions.get(strPos) + ".tsv";
			File pFile = new File(strFile);
			
			if(!pFile.exists())
			{
				Messagebox.show("ERROR: No GTEX data available for exonic part: " + strGTEXID + " " + strPos + " " + mapIDsToPositions.get(strPos) + " -> " + strFile);
				continue;
			}
			
			pScanner = new Scanner(pFile);
			
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

	public void OnProjectChange(String strProjectFile) throws Exception
	{
		m_strSelectedCondition 			= null;
		m_strSelectedConditionType		= null;
		
		m_projectModel.clear();
		m_projectModel.Init(strProjectFile, -1, -1, -1.0, -1, true);		

		m_Gene							= null;
		m_bandboxSelectedGene.setText("Select Gene");
		
		UpdateComboboxesForCondition();
		
		// select first available condition type
		m_comboboxSelectedConditionType.setSelectedIndex(0);
		m_strSelectedConditionType = m_comboboxSelectedConditionType.getSelectedItem().getLabel();
		OnConditionTypeChange();
		
		// select first available condition
		m_comboboxSelectedCondition.setSelectedIndex(0);
		m_strSelectedCondition = m_comboboxSelectedCondition.getSelectedItem().getLabel();
		
		OnGeneChange(true);
		
		m_vcASHits = m_projectModel.LoadHitList(m_strPathHitLists);
		UpdateHitsList();
		
		// use a reduced data set if there are more than 10 samples per condition. Scales up to 50 samples in total.
		if(m_bUseReducedDataSet)
		{
			ReduceSampleSize();
		}
	}
	
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
		m_bCoverageRequiresRedraw = true;
		
		// intially select all samples		
		for(String strSample : m_projectModel.GetSamples())
		{
			m_vcSelectedSamples.add(strSample);
			m_vcPreviouslySelectedTreeItems.add(strSample);
		}
		
		for(String strCondition : m_projectModel.GetConditions(m_strSelectedConditionType))
			m_vcPreviouslySelectedTreeItems.add(strCondition);
		
		GetColorsPerCondition();
		UpdateComboboxForSampleSelection(true);
		UpdateColorSelection();
	}
	
	public void OnGeneChange(boolean bProjectChanged) throws Exception
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
		
		// clear highlighted exon
		m_vcHighlightedExons.clear();
		m_vcHighlightedJunctions.clear();
		
		m_pExons		= null;
		m_pExonGroups 	= null;
		
		m_mapAnalysedJunctionPaths	 			= new TreeMap<String, TreeMap<Integer, TreeSet<String>>>();
		m_mapCoveragePerJunctionPathAndSample 	= new TreeMap<String, TreeMap<Integer, Double>>();
		
		m_mapExonCoveragePerSample			= null;
		m_mapIndistinguishableExonPerSample = null;
				
		m_mapEntropyToExonicPart 			= null;
		m_mapCountsPerExonAndTissue			= null;
		
		m_bCoverageRequiresRedraw = true;
		
		if(m_bandboxSelectedGene.getText().equals("Select Gene"))
			return;
		
		m_ClickEvent.clear();

		GeneIdentifier gid = GetGeneIdentifierForGene(m_bandboxSelectedGene.getValue());
		if(gid == null)
		{
			Messagebox.show("Gene not found: " + m_bandboxSelectedGene.getValue() + "(" + m_bandboxSelectedGene.getText() + ")");
			return;
		}

		ProcessGeneInformation(gid, m_strFileGTF);
		
		// add splicing hits if available
		m_vcSplicingScores = m_projectModel.GetSplicingHitForGene(m_strPathSplicingHitList, m_Gene);
		
		if(m_Gene != null)
		{
			OnIsoformChange();
			if(m_bShowEntropyData)
			{
				Events.postEvent("onSelect", m_comboboxSelectedEntropyIndex, null);
			}
			
			DrawPlots();
		}
		
		UpdatePSIHitList();
	}

	public boolean OnIsoformChange() throws IOException
	{
		TreeMap<String, int[]> mapExonsToValidIsoforms = new TreeMap<String, int[]>();

		// get exons for valid isoform 
		for(String strIsoform : m_vcValidIsoforms)
		{
			for(String strIncludedIsoform : m_mapExonsToIsoforms.keySet())
			{
				if(strIsoform.equals(strIncludedIsoform) || strIsoform.equals(strIncludedIsoform.split("\\.")[0]))
				{
					mapExonsToValidIsoforms.put(strIsoform, m_mapExonsToIsoforms.get(strIncludedIsoform));
					break;
				}
			}			
		}

		// isoforms have been removed, thus recalculate the exon groups
		m_pExonGroups = RecalculateExonGroups(mapExonsToValidIsoforms);

		// recalculate expression values
		if(!ProcessCountData())
			return false;
		
		return true;
	}
	
	public void AddSettingsMenu(Vlayout parentLayout) throws Exception
	{
		Groupbox grpBox = new Groupbox();
		grpBox.setTitle("Settings");
		grpBox.setMold("3d");
		grpBox.setParent(parentLayout);
		grpBox.setWidth("500px");
		grpBox.setStyle("margin-left: 10px;");
		
		Hlayout layoutH = new Hlayout();
		layoutH.setParent(grpBox);
		
		Vlayout layoutV = new Vlayout();
		layoutV.setParent(layoutH);
		
		// add label for window size edit box
		Label lab = new Label("Window width in pixel:");
		lab.setStyle("margin-left: 10px");
		lab.setParent(layoutV);
		
		// add box for window size specification
		Textbox m_textboxWindowWidth 	= new Textbox("Enter value");
		m_textboxWindowWidth.setParent(layoutV);
		m_textboxWindowWidth.setWidth("180px");
		m_textboxWindowWidth.setStyle("margin-left: 10px;");
		m_textboxWindowWidth.addEventListener(Events.ON_OK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Textbox box = (Textbox)event.getTarget();
				String strValue = box.getValue();
				if(strValue != null)
				{
					try
					{
						m_nMaxWidth = Integer.parseInt(strValue);
						m_bCoverageRequiresRedraw = true;
						
						if(m_Gene != null)
							DrawPlots();
					}
					catch(NumberFormatException e)
					{
						if(m_Gene == null)
						{
							Messagebox.show("Invalid input: " + strValue + ". The window width must be a valid number.");
							return;
						}
					}
				}
			}
		});

		layoutV = new Vlayout();
		layoutV.setStyle("margin-left: 40px");
		layoutV.setParent(layoutH);
		
		// add label for project selection
		lab = new Label("Selected Project:");
		lab.setParent(layoutV);
		
		// add combobox for project selection
		Combobox combobox = new Combobox("Select Project");
		combobox.setParent(layoutV);
		combobox.setWidth("180px");
		
		// create tree set of project names to get a sorted list
		TreeMap<String, String> vcProjects = new TreeMap<String, String>();
		File pFolder = new File(m_strPathInput);
		for(File pFile : pFolder.listFiles())
		{
			String strFile = pFile.getName();
			
			vcProjects.put(strFile, pFile.getAbsolutePath());
		}
		
		for(String strFile : vcProjects.keySet())
		{
			if(strFile.endsWith("project"))
			{
				Comboitem item = new Comboitem(strFile);
				item.setValue(vcProjects.get(strFile));
				combobox.appendChild(item);
			}
		}
		combobox.addEventListener(Events.ON_SELECT, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Combobox box = (Combobox)event.getTarget();
				if(box.getSelectedItem() != null)
				{
					OnProjectChange((String)box.getSelectedItem().getValue());
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
		layoutV.setParent(layoutH);
		
		// add label for gene identifier selection
		lab = new Label("Selected Gene:");
		lab.setStyle("margin-left: 10px");
		lab.setParent(layoutV);
		
		// add text box for gene identifier selection
		m_bandboxSelectedGene = new Bandbox();
		m_bandboxSelectedGene.setParent(layoutV);
		m_bandboxSelectedGene.setMold("rounded");
		m_bandboxSelectedGene.setStyle("padding-left: 10px;");
		m_bandboxSelectedGene.setAutodrop(true);
		m_bandboxSelectedGene.setVflex("true");
//		m_bandboxSelectedGene.setWidth("170px");
		
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
		m_comboboxSelectedGeneAnnotation.setWidth("180px");
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
				}
			}
		});
		
		m_bandboxSelectedGene.addEventListener(Events.ON_OK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				OnGeneChange(false);
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
					TreeSet<GeneIdentifier> vcGeneIdentifier = new TreeSet<GeneIdentifier>();

					// get list of valid gene identifiers
					for(GeneIdentifier gid : m_vcGeneIdentifier)
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
									vcGeneIdentifier.add(gid);
								}
							}
						}
					}
					
					// now add the sorted list of gene identifiers to the selection
					for(GeneIdentifier gid : vcGeneIdentifier)
					{
						Listitem item = new Listitem();
						Listcell cell = new Listcell(gid.m_strApprovedGeneSymbol); cell.setParent(item);
						item.setValue(gid.m_strEnsemblGeneID);
						cell = new Listcell(gid.m_strEnsemblGeneID); cell.setParent(item);
						cell = new Listcell(gid.m_strRefGeneID); cell.setParent(item);
						cell = new Listcell(gid.m_strEntrezGeneID); cell.setParent(item);
						cell = new Listcell(gid.m_strSynonyms); cell.setParent(item);
						cell = new Listcell(gid.m_strApprovedGeneName); cell.setParent(item);
						m_listboxSelectedGene.appendChild(item);
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

				OnGeneChange(false);
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
					m_nOrganismType = 0;
					String strFileCrossReference = m_strPathReferences + "/biomart_cross_ref_human.txt";
					ReadGeneIdentifier(strFileCrossReference);
					LoadGeneAnnotation(m_strFileGTF);
				}
				else if(btn.getLabel().equals("mouse"))
				{
					m_nOrganismType = 1;
					String strFileCrossReference = m_strPathReferences + "/biomart_cross_ref_mouse.txt";
					ReadGeneIdentifier(strFileCrossReference);
					LoadGeneAnnotation(m_strFileGTF);
				}
			}
		});

		
		rgroup.setParent(layoutH);
	}

	public void AddOptionsMenu(Vlayout parentLayout)
	{
		Groupbox grpBox = new Groupbox();
		grpBox.setTitle("Options");
		grpBox.setMold("3d");
		grpBox.setParent(parentLayout);
		grpBox.setWidth("500px");
		grpBox.setStyle("margin-left: 10px;");
		grpBox.setHeight("270px");
		
		Hlayout layoutH = new Hlayout();
		layoutH.setParent(grpBox);
		
		Vlayout layoutV = new Vlayout();
		layoutV.setWidth("220px");
		layoutV.setParent(layoutH);
		
		Label label = new Label("Isoform view options");
		label.setStyle("text-decoration:underline");
		label.setParent(layoutV);
		
		Checkbox checkBox = new Checkbox("Show unique isoform elements");
		checkBox.setParent(layoutV);
		checkBox.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bShowUniqueFeatures = !m_bShowUniqueFeatures;
				m_bShowAmbigiousUniqueExons = !m_bShowAmbigiousUniqueExons;
				DrawPlots();
			}
		});
		
		checkBox = new Checkbox("Show tissue specificity index");
		checkBox.setParent(layoutV);
		checkBox.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bShowEntropyData = !m_bShowEntropyData;
				
				if(m_bShowEntropyData)
				{
					m_mapCountsPerExonAndTissue = GetGTEXDataForExonicParts();
					
					if(m_mapCountsPerExonAndTissue == null)
						return;
					
					m_mapEntropyToExonicPart = CalculateTheilIndexForExons(m_mapCountsPerExonAndTissue);
				}
				DrawPlots();
			}
		});
		
		checkBox = new Checkbox("Color exons and junctions by expression level");
		checkBox.setParent(layoutV);
		checkBox.setChecked(true);
		checkBox.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bColorExonsAndJunctionsByCoverage = !m_bColorExonsAndJunctionsByCoverage;
				m_bColorExonsByRelativeCoverage		= false;
				m_bUseRelativeExonCoverageMode1		= false;
				m_bUseRelativeExonCoverageMode2		= false;
				
				m_checkboxRelativeCoverage1.setChecked(false);
				m_checkboxRelativeCoverage2.setChecked(false);
				DrawPlots();
			}
		});
		
		label = new Label("Coverage view options");
		label.setStyle("margin-top: 10px; display:inline-block; text-decoration:underline");
		label.setParent(layoutV);
		
		m_checkboxShowRelativeCoverage = new Checkbox("Show relative coverage");
		m_checkboxShowRelativeCoverage.setWidth("100%");
		m_checkboxShowRelativeCoverage.setParent(layoutV);
		m_checkboxShowRelativeCoverage.setChecked(m_bCoveragePlotShowRelativeCoverage);
		m_checkboxShowRelativeCoverage.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bCoveragePlotShowRelativeCoverage = m_checkboxShowRelativeCoverage.isChecked();
			}
		});	
		
		m_checkboxUseLog2 = new Checkbox("Show log2 tranformed");
		m_checkboxUseLog2.setWidth("100%");
		m_checkboxUseLog2.setParent(layoutV);
		m_checkboxUseLog2.setChecked(m_bCoveragePlotUseLog2);
		m_checkboxUseLog2.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bCoveragePlotUseLog2 = m_checkboxUseLog2.isChecked();
			}
		});	
		
		layoutV = new Vlayout();
		layoutV.setWidth("230px");
		layoutV.setParent(layoutH);
		
		//################################################
		//    Add button to remove irrelevant isoforms
		//################################################
		Button btn = new Button("Show relevant isoforms only");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bHideIrrelevantIsoforms = true;				
				UpdateComboboxForIsoformSelection(false);
				OnIsoformChange();
				DrawPlots();
			}
		});
		
		//############################################
		//     Add button to view GTEX expression
		//############################################
		btn = new Button("Show GTEX data for gene");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				if(m_Gene == null)
				{
					Messagebox.show("No gene specified");
					return;
				}
				
				ShowGTEXDataForGene(m_Gene.getGeneID());
			}
		});

		//############################################
		//     Add button to show MMSeq estimates
		//############################################
		btn = new Button("Calculate MMSeq estimates");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				if(m_Gene == null)
				{
					Messagebox.show("No gene specified");
					return;
				}
				
				GenerateMMSeqEstimates();
				DrawPlots();
			}
		});
		
		//############################################
		//     Add button to redraw isoforms
		//############################################
		btn = new Button("Redraw plots");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				if(m_Gene == null)
				{
					Messagebox.show("No gene specified");
					return;
				}
								
				if(m_bIsoformSelectionChanged)
				{
					UpdateComboboxForIsoformSelection(false);
					OnIsoformChange();
				}
				
				if(m_bSampleSelectionChanged)
				{
					UpdateComboboxForSampleSelection(true);
				}
				
				m_bCoverageRequiresRedraw = true;
				DrawPlots();
			}
		});
		
		//############################################
		//     Add button to redraw isoforms
		//############################################
		btn = new Button("Save screenshot");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				GeneIdentifier gid = GetGeneIdentifierForGene(m_Gene.getGeneID());

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
		btn = new Button("Calculate splicing scores");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event arg0) throws Exception
			{
				if(m_vcSplicingScores != null)
					m_vcSplicingScores.clear();
				
				m_vcSplicingScores = CalculateNovelJunctionPSIScores(false);
				
				if(m_vcSplicingScores == null || m_vcSplicingScores.size() == 0)
				{
					Messagebox.show("No alternative splicing events detected\n");
				}
				else
				{
					m_projectModel.SaveSplicingHitsToFile(m_strPathSplicingHitList, m_vcSplicingScores);
				}
				
				UpdatePSIHitList();
			}
		});
		
		//############################################
		//     Add button to calculate PSI scores
		//############################################
		btn = new Button("Show splicing heat map");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event arg0) throws Exception
			{
				DrawJunctionHeatmap();
			}
		});
			
		//############################################
		//    Add button to save hits to hit list
		//############################################
		btn = new Button("Add to list of alternative splicing events");
		btn.setParent(layoutV);
		btn.setWidth("230px");
		btn.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				if(m_Gene == null)
				{
					Messagebox.show("No gene specified");
					return;
				}

				// get all currently selected/visible isoforms
				TreeSet<String> vcValidIsoforms = new TreeSet<String>();
				for(String strIsoform : m_vcValidIsoforms)
				{
					if(m_vcValidIsoforms.contains(strIsoform))
					{
						for(String strIncludedIsoform : m_mapExonsToIsoforms.keySet())
						{
							if(strIsoform.equals(strIncludedIsoform) || strIsoform.equals(strIncludedIsoform.split("\\.")[0]))
							{
								vcValidIsoforms.add(strIsoform);
								break;
							}
						}
					}
				}
				
				if(m_vcASCurrentHit == null)
				{
					GeneIdentifier gid = GetGeneIdentifierForGene(m_Gene.getGeneID());
					if(gid == null)
					{
						gid = GetGeneIdentifierForGene(m_Gene.getGeneName());
						
						if(gid == null)
						{
							gid = new GeneIdentifier();
							gid.m_strEntrezGeneID = m_Gene.getGeneID();
							gid.m_strApprovedGeneName = m_Gene.getGeneID();
							gid.m_strApprovedGeneSymbol = m_Gene.getGeneName();
						}
					}
					
					//TODO
//					m_vcASCurrentHit = new AlternativeSplicingHit(0, gid, m_nMinJunctionReads, m_nMinCovPerBase,
//							m_fMinCoveredBases, m_fVariableExonThreshold, m_vcValidIsoforms, m_strFileGTF, "", m_projectModel.GetConditions(m_strSelectedConditionType));
				}
				
				m_vcASCurrentHit.SetIsoforms(vcValidIsoforms);
				
				if(m_vcASHits == null)
				{
					m_vcASHits = new TreeSet<AlternativeSplicingHit>();
				}
				
				if(m_vcASHits.contains(m_vcASCurrentHit))
					m_vcASHits.remove(m_vcASCurrentHit);
				m_vcASHits.add(m_vcASCurrentHit);

				UpdateHitsList();
			}
		});
		
	}
	
	public void AddAdvancedOptionsMenu(Vlayout parentLayout)
	{
		Groupbox grpBox = new Groupbox();
		grpBox.setTitle("Advanced Options");
		grpBox.setMold("3d");
		grpBox.setParent(parentLayout);
		grpBox.setWidth("500px");
		grpBox.setHeight("470px");

		Hlayout layoutH = new Hlayout();
		layoutH.setParent(grpBox);
		
		Vlayout layoutV = new Vlayout();
		layoutV.setWidth("220px");
		layoutV.setParent(layoutH);

		//#############################################################
		layoutH = new Hlayout();
		layoutH.setParent(grpBox);
		
		layoutV = new Vlayout();
		layoutV.setWidth("200px");
		layoutV.setParent(layoutH);
		
		m_checkboxRelativeCoverage1 = new Checkbox("Show exon coverage relative to total coverage");
		m_checkboxRelativeCoverage1.setParent(layoutV);
		m_checkboxRelativeCoverage1.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{				
				m_bUseRelativeExonCoverageMode1 = !m_bUseRelativeExonCoverageMode1;
				m_bColorExonsByRelativeCoverage = m_bUseRelativeExonCoverageMode1;
				
				if(m_bUseRelativeExonCoverageMode1)
				{
					m_bShowRelativeCoverageNumbersExons = true;
				}
				
				m_bUseRelativeExonCoverageMode2 = false;
				m_checkboxRelativeCoverage2.setChecked(false);
				DrawPlots();
			}
		});
		
		m_checkboxRelativeCoverage2 = new Checkbox("Show exon coverage relative to maximum coverage");
		m_checkboxRelativeCoverage2.setParent(layoutV);
		m_checkboxRelativeCoverage2.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bUseRelativeExonCoverageMode2 = !m_bUseRelativeExonCoverageMode2;
				m_bColorExonsByRelativeCoverage = m_bUseRelativeExonCoverageMode2;
				
				if(m_bUseRelativeExonCoverageMode2)
				{
					m_bShowRelativeCoverageNumbersExons = true;
				}
				
				m_bUseRelativeExonCoverageMode1 = false;
				m_checkboxRelativeCoverage1.setChecked(false);
				DrawPlots();
			}
		});

		
		//#############################################################
		layoutV = new Vlayout();
		layoutV.setStyle("margin-left: 20px");
		layoutV.setWidth("220px");
		layoutV.setParent(layoutH);
		
		Checkbox checkBox = new Checkbox("Color junctions by relative coverage");
		checkBox.setParent(layoutV);
		checkBox.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bColorJunctionPaths = !m_bColorJunctionPaths;
				DrawPlots();
			}
		});
		
		checkBox = new Checkbox("Show relative coverage for junctions");
		checkBox.setParent(layoutV);
		checkBox.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bShowRelativeCoverageNumbersJunctions = !m_bShowRelativeCoverageNumbersJunctions;
				DrawPlots();
			}
		});
		
		//#############################################################
		layoutH = new Hlayout();
		layoutH.setStyle("margin-top: 20px");
		layoutH.setParent(grpBox);
		
		layoutV = new Vlayout();
		layoutV.setWidth("220px");
		layoutV.setParent(layoutH);
		
		checkBox = new Checkbox("Color alternative spliced exons");
		checkBox.setParent(layoutV);
		checkBox.setChecked(m_bShowASEvents);
		checkBox.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bShowASEvents = !m_bShowASEvents;
				DrawPlots();
			}
		});

		checkBox = new Checkbox("Attach coverage plot to isoforms");
		checkBox.setParent(layoutV);
		checkBox.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bShowSecondCoveragePlot = !m_bShowSecondCoveragePlot;
				m_bCoverageRequiresRedraw = true;
				DrawPlots();
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
					GetBigWigCoverageForGene();
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
		m_checkboxUseMedian.setChecked(m_bCoveragePlotShowMedian);
		m_checkboxUseMedian.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bCoveragePlotShowMedian = m_checkboxUseMedian.isChecked();
				
				m_bCoveragePlotShowMean				= false;
				m_bCoveragePlotShowGeometricMean	= false;
				m_checkboxUseMean.setChecked(false);
				m_checkboxUseGeometricMean.setChecked(false);
			}
		});
		
		m_checkboxUseGeometricMean = new Checkbox("Show geometric mean per condition");
		m_checkboxUseGeometricMean.setWidth("100%");
		m_checkboxUseGeometricMean.setParent(layoutV);
		m_checkboxUseGeometricMean.setChecked(m_bCoveragePlotShowGeometricMean);
		m_checkboxUseGeometricMean.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bCoveragePlotShowGeometricMean = m_checkboxUseGeometricMean.isChecked();
				
				m_bCoveragePlotShowMean		= false;
				m_bCoveragePlotShowMedian	= false;
				m_checkboxUseMean.setChecked(false);
				m_checkboxUseMedian.setChecked(false);
			}
		});
		
		m_checkboxUseMean = new Checkbox("Show mean per condition");
		m_checkboxUseMean.setWidth("100%");
		m_checkboxUseMean.setParent(layoutV);
		m_checkboxUseMean.setChecked(m_bCoveragePlotShowMean);
		m_checkboxUseMean.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bCoveragePlotShowMean = m_checkboxUseMean.isChecked();
				
				m_bCoveragePlotShowMedian			= false;
				m_bCoveragePlotShowGeometricMean	= false;
				m_checkboxUseMedian.setChecked(false);
				m_checkboxUseGeometricMean.setChecked(false);
			}
		});
		
		Checkbox checkboxShowQuartiles = new Checkbox("Show quartiles");
		checkboxShowQuartiles.setWidth("100%");
		checkboxShowQuartiles.setParent(layoutV);
		checkboxShowQuartiles.setChecked(m_bCoveragePlotShowQuartiles);
		checkboxShowQuartiles.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_bCoveragePlotShowQuartiles = !m_bCoveragePlotShowQuartiles;
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
		
		m_textboxMinJunctionReads = new Textbox("1");
		m_textboxMinJunctionReads.setParent(layoutV);
		m_textboxMinJunctionReads.addEventListener(Events.ON_CHANGE, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Textbox box = (Textbox)event.getTarget();
				try
				{
					m_nMinJunctionReads = Integer.parseInt(box.getText());
				}
				catch(NumberFormatException e)
				{
					Messagebox.show("The minimum junction read must be a non-negative number");
				}
			}
		});
		
		lab = new Label("Minimum coverage per base:");
		lab.setParent(layoutV);
		
		m_textboxMinCovPerBase = new Textbox("3");
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
		
		m_textboxMinCoveredBases = new Textbox("0.7");
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
		
		m_textboxVariableExonThreshold = new Textbox("0.15");
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
		layoutH.setStyle("margin-top: 20px");
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
					m_mapCountsPerExonAndTissue = GetGTEXDataForExonicParts();
					
					if(m_mapCountsPerExonAndTissue == null)
						return;
					
					if(m_mapEntropyToExonicPart != null)
						m_mapEntropyToExonicPart.clear();
					else
					{
						m_mapEntropyToExonicPart = new TreeMap<String, Double>();
					}
	
					switch((String)box.getSelectedItem().getValue())
					{
						case "None":
							m_bShowEntropyData = false;
							break;
							
						case "Gini":
							m_mapEntropyToExonicPart = CalculateGiniIndexForExons(m_mapCountsPerExonAndTissue);
							m_bShowEntropyData = true;
							break;
							
						case "Simpson":
							m_mapEntropyToExonicPart = CalculateGiniSimpsonIndexForExons(m_mapCountsPerExonAndTissue);
							m_bShowEntropyData = true;
							break;
							
						case "Shannon":
							m_mapEntropyToExonicPart = CalculateShannonIndexForExons(m_mapCountsPerExonAndTissue);
							m_bShowEntropyData = true;
							break;
							
						case "Theil":
							m_mapEntropyToExonicPart = CalculateTheilIndexForExons(m_mapCountsPerExonAndTissue);
							m_bShowEntropyData = true;
							break;
	
						default:
							m_bShowEntropyData = true;
							break;
					}
					
					// requires redraw
					DrawPlots();
				}
			}
		});
	}

	public void AddColorSelectionMenu(Vlayout parentLayout)
	{
		m_ColorSelectionGroupBox = new Groupbox();
		m_ColorSelectionGroupBox.setTitle("Colors");
		m_ColorSelectionGroupBox.setParent(parentLayout);		
		m_ColorSelectionGroupBox.setWidth("200px");
		m_ColorSelectionGroupBox.setHeight("470px");
		m_ColorSelectionGroupBox.setContentStyle("overflow:auto;");
		m_ColorSelectionGroupBox.setMold("3d");
	}
	
	public void UpdateColorSelection()
	{				
		TreeMap<String, TreeSet<String>> vcSamplesAndConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		int nHeight = vcSamplesAndConditions.keySet().size()*30;
		BufferedImage img = new BufferedImage(200, nHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D graph = img.createGraphics();
		
		m_ColorSelectionGroupBox.getChildren().clear();

		m_imgMapColors = new Imagemap();
		m_imgMapColors.setWidth(200+"px");
		m_imgMapColors.setHeight(nHeight+"px");
		m_imgMapColors.setParent(m_ColorSelectionGroupBox);

		// define selected conditions
		TreeSet<String> vcSelectedConditions = new TreeSet<String>();
		for(String strCondition : vcSamplesAndConditions.keySet())
		{
			for(String strSample : m_vcSelectedSamples)
			{
				if(vcSamplesAndConditions.get(strCondition).contains(strSample))
					vcSelectedConditions.add(strCondition);
			}
		}
		
		// fill background
		graph.setColor(Color.WHITE);
		graph.fillRect(0, 0, 200, nHeight);
		
		int nOffset = 0;
		for(String strCondition : vcSelectedConditions)
		{
			Color clr = m_mapColorsToConditions.get(strCondition);
			graph.setColor(clr);
			graph.fillRect(3, nOffset, 20, 20);
			graph.drawString(strCondition, 28, nOffset+15);
			
			graph.setColor(Color.BLACK);
			graph.drawRect(3, nOffset, 20, 20);
			
			// add clickable area
			Area area = new Area();
			area.setId("clr§" + strCondition);
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

				String strCondition = evnt.getArea().split("§")[1];
				
				m_windowPopup = new Window();
				m_windowPopup.setParent(m_hWindow);
				m_windowPopup.doPopup();
				m_windowPopup.setStyle("relative;left:500px;top:100px;");
				
				int nWidth = 240;
				int nHeight = 220;
				
				m_windowPopup.setWidth(nWidth + "px");
				m_windowPopup.setHeight(nHeight +"px");
				
				Imagemap imgMap = new Imagemap();
				imgMap.setWidth(nWidth + "px");
				imgMap.setHeight(nHeight + "px");
				
				UpdateColorSelectionPopup(imgMap, nWidth, nHeight, m_mapColorsToConditions.get(strCondition), strCondition);
			}
		});
	}
	
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
		imgMap.setParent(m_windowPopup);
		
		imgMap.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				MouseEvent evnt = (MouseEvent) event;
				
				int nMouseX = evnt.getX(); 
				int nMouseY = evnt.getY();
				
				// was the saturation hit?
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
					UpdateColorSelectionPopup(imgMap, 240, 220, clr, strCondition);
				}
				// was the saturation hit?
				else if(nMouseX >= 5 && nMouseX <= 205 && nMouseY >= 5 && nMouseY <= 205)
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
					graph.fillRect(5, 5, 200, 200);
					graph.setPaint(clrGradient2);
					graph.fillRect(5, 5, 200, 200);
					
					int rgb = img.getRGB(nMouseX, nMouseY);
					r = (rgb & 0x00ff0000) >> 16;
					g = (rgb & 0x0000ff00) >> 8;
					b = rgb & 0x000000ff;
					
					clr = new Color(img.getRGB(nMouseX, nMouseY));
					
					m_mapColorsToConditions.put(strCondition, clr);
					
					UpdateColorSelection();
				}
			}
		});
	}
	
	public void AddHitList(Tabs tabs, Tabpanels panels)
	{
		Tab tab = new Tab("Hits");
		tab.setParent(tabs);
		
		Tabpanel panel = new Tabpanel();
		panel.setParent(panels);
		
		Vlayout layoutV = new Vlayout();
		layoutV.setParent(panel);
		
		m_gridHitList = new Grid();
		m_gridHitList.setId("hitlist");
		m_gridHitList.setParent(layoutV);
		m_gridHitList.setWidth("830px");
		m_gridHitList.setHeight("380px");
		m_gridHitList.setVflex("min");
		m_gridHitList.setMold("paging");
		m_gridHitList.setPageSize(10);
		
		Button btnSaveHitList = new Button("Save changes");
		btnSaveHitList.setParent(layoutV);
		btnSaveHitList.setStyle("margin-left: 10px;");
		btnSaveHitList.setWidth("200px");
		btnSaveHitList.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				m_projectModel.SaveHitListToFile(m_strPathHitLists, m_vcASHits);
				Messagebox.show("Changes saved.");
			}
		});

	}
		
	public void UpdateHitsList() throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		if(m_vcASHits == null)
			return;
		
		m_gridHitList.setMold("default");
		// clear old rows
		m_gridHitList.getChildren().clear();
		m_gridHitList.setMold("paging");
		m_gridHitList.setPageSize(10);

		//##################################
		//           add header
		//##################################
		Columns cols = new Columns();

		Column col = new Column("");
		col.sort(true);
		col.setWidth("40px");
		col.setParent(cols);
		
		col = new Column("Rating");
//		col.sort(true);
//		col.setSort("auto");
		col.setSortAscending(getAscRatingComparator());
		col.setSortDescending(getDescRatingComparator());
		col.setWidth("110px");
		col.setParent(cols);

		col = new Column("Gene ID");
		col.setSort("auto");
		col.setWidth("140px");
		col.setParent(cols);
		
		col = new Column("Gene Symbol");
		col.setSort("auto");
		col.setWidth("100px");
		col.setParent(cols);
		
		col = new Column("Comment");
		col.setSort("auto");
		col.setWidth("410px");
		col.setParent(cols);
		
		col = new Column("");
		col.setWidth("25px");
		col.setParent(cols);
	
		cols.setParent(m_gridHitList);
		
		Rows rows = new Rows();
		rows.setParent(m_gridHitList);
		
		//##################################
		//             add rows
		//##################################	
		for(AlternativeSplicingHit hit : m_vcASHits)
		{
			Row row = new Row();
			row.setId("rating_" + hit.GetGeneID() + "_" + hit.GetRating());
			row.setParent(rows);
			
			Image img = new Image("/img/view.png");
			img.setHeight("12px");
			img.setParent(row);
			img.setId("view_" + hit.GetGeneID());
			img.setStyle("display: block; margin-left: auto; margin-right: auto; margin-top: 4px; margin-bottom: auto; cursor:hand; cursor:pointer;");
			
			img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
			{
				public void onEvent(Event event) throws Exception
				{
					Image img = (Image)event.getTarget();
					String strText = img.getId().split("_")[1];
					
					ShowHit(strText);						
				}
			});
			
			Hlayout layoutH = new Hlayout();
			layoutH.setParent(row);
			
			for(int i=0; i<5; i++)
			{
				img = null;
				
				if(i<hit.GetRating())
				{
					img = new Image("/img/good_rating.png");
				}
				else
				{
					img = new Image("/img/neutral_rating.png");
				}
				
				img.setHeight("15px");
				img.setParent(layoutH);
				img.setId("rating_" + (i+1) + "_" + hit.GetGeneID());
				img.setStyle("cursor:hand; cursor:pointer;");
				
				img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
				{
					public void onEvent(Event event) throws Exception
					{
						Image img = (Image)event.getTarget();
						
						String strID = img.getId();
						
						String pSplit[] = strID.split("_");
						int nNewRating = Integer.parseInt(pSplit[1]);
						String strGene = pSplit[2];
						
						int i=0;
						for(Component c : img.getParent().getChildren())
						{
							Image imgRating = (Image) c;
							imgRating.invalidate();
							
							if(i < nNewRating)
							{
								imgRating.setSrc("/img/good_rating.png");
							}
							else
							{
								imgRating.setSrc("/img/neutral_rating.png");
							}							
							i++;
						}
						
						for(AlternativeSplicingHit hit : m_vcASHits)
						{
							if(hit.GetGeneID().equals(strGene))
							{
								hit.SetRating(nNewRating);
								img.getParent().getParent().setId("rating_" + hit.GetGeneID() + "_" + hit.GetRating());
								break;
							}
						}
					}
				});
			}
			
			Label label = new Label(hit.GetGeneID());			
			label.setParent(row);
			
			label = new Label(hit.GetGeneSymbol());
			label.setParent(row);
			
			Textbox boxComment = new Textbox();
			boxComment.setParent(row);
			boxComment.setInplace(true);
			boxComment.setWidth("99%");
			boxComment.setId("commentField_" + hit.GetGeneID());
			boxComment.setText(hit.GetComment());
			boxComment.setTooltiptext(hit.GetComment());
			
			boxComment.addEventListener(Events.ON_CHANGE, new EventListener<Event>()
			{
				public void onEvent(Event event) throws Exception
				{
					Textbox box = (Textbox)event.getTarget();
					
					String strID = box.getId();
					
					String pSplit[] = strID.split("_");
					String strGene = pSplit[1];
					
					for(AlternativeSplicingHit hit : m_vcASHits)
					{
						if(hit.GetGeneID().equals(strGene))
						{
							hit.SetComment(box.getText());
						}
					}				
				}
			});
			
			img = new Image("/img/red_cross.png");
			img.setHeight("15px");
			img.setParent(row);
			img.setId("delete_" + hit.GetGeneID());
			img.setStyle("display: block; margin-left: auto; margin-right: auto; margin-top: auto; margin-bottom: auto; cursor:hand; cursor:pointer;");
			
			img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
			{
				public void onEvent(Event event) throws Exception
				{
					Image img = (Image)event.getTarget();
					
					String strID = img.getId();
					
					String pSplit[] = strID.split("_");
					String strGene = pSplit[1];
					
					for(AlternativeSplicingHit hit : m_vcASHits)
					{
						if(hit.GetGeneID().equals(strGene))
						{
							m_vcASHits.remove(hit);
							break;
						}
					}
				
					UpdateHitsList();					
				}
			});
		}
	}
	
	public void AddPSIHitList(Tabs tabs, Tabpanels panels)
	{
		Tab tab = new Tab("Splicing scores");
		tab.setParent(tabs);
		
		Tabpanel panel = new Tabpanel();
		panel.setParent(panels);
		
		Vlayout layoutV = new Vlayout();
		layoutV.setParent(panel);
		
		m_gridPSIHitList = new Grid();
		m_gridPSIHitList.setId("psi_hitlist");
		m_gridPSIHitList.setParent(layoutV);
		m_gridPSIHitList.setWidth("830px");
		m_gridPSIHitList.setHeight("380px");
		m_gridPSIHitList.setVflex("min");
		m_gridPSIHitList.setMold("paging");
		m_gridPSIHitList.setPageSize(10);

		Button btnSaveHitList = new Button("Save changes");
		btnSaveHitList.setParent(layoutV);
		btnSaveHitList.setStyle("margin-left: 10px;");
		btnSaveHitList.setWidth("200px");
		btnSaveHitList.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				m_projectModel.SaveSplicingHitsToFile(m_strPathSplicingHitList, m_vcSplicingScores);
				Messagebox.show("Changes saved.");
			}
		});
	}
	
	public void UpdatePSIHitList() throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		// clear old rows
		m_gridPSIHitList.setMold("default");
		m_gridPSIHitList.getChildren().clear();
		
		m_gridPSIHitList.setMold("paging");
		m_gridPSIHitList.setPageSize(10);
	
		if(m_vcSplicingScores == null)
			return;	
		
//		m_SplicingScores.Report();

		//##################################
		//           add header
		//##################################
		Columns cols = new Columns();
		cols.setMenupopup("auto");

		Column col = new Column("");
		col.sort(true);
		col.setWidth("40px");
		col.setParent(cols);
		
		col = new Column("Rating");
		col.setWidth("110px");
		col.setParent(cols);

		col = new Column("Info");
		col.setWidth("30px");
		col.setParent(cols);
		
		col = new Column("Novel");
		col.setSort("auto");
		col.setWidth("60px");
		col.setParent(cols);

		col = new Column("Gene Symbol");
		col.setSort("auto");
		col.setWidth("100px");
		col.setParent(cols);
		
		col = new Column("Condition A");
		col.setSort("auto");
		col.setWidth("80px");
		col.setParent(cols);
		
		col = new Column("Condition B");
		col.setSort("auto");
		col.setWidth("80px");
		col.setParent(cols);
		
		col = new Column("Exon");
		col.setSort("auto");
		col.setWidth("180px");
		col.setParent(cols);
		
		col = new Column("effect");
		col.setSort("client(number)");
		col.setWidth("65px");
		col.setParent(cols);
		
		col = new Column("P-Value");
		col.setSort("client(number)");
		col.setWidth("82px");
		col.setParent(cols);
	
		cols.setParent(m_gridPSIHitList);
		
		Rows rows = new Rows();
		rows.setParent(m_gridPSIHitList);
		
		// get a tree set of p-values for sorting
		TreeSet<Double> vcPValues = new TreeSet<Double>();
		for(SimpleSpliceScore res : m_vcSplicingScores)
			vcPValues.add(res.m_fPValue);
		
		for(SimpleSpliceScore res : m_vcSplicingScores)
		{
			/*
			SimpleSpliceScore res = null;
			for(SimpleSpliceScore r : m_vcSplicingScores)
			{
				if(r.m_fPValue == fVal)
				{
					res = r;
					break;
				}
			}
			*/
				
			Row row = new Row();
			row.setId("rating_" + res.m_strID);
			row.setParent(rows);
			
			Image img = new Image("/img/view.png");
			img.setHeight("12px");
			img.setParent(row);
			img.setId("view_splicing" + res.m_strID);
			img.setStyle("display: block; margin-left: auto; margin-right: auto; margin-top: 4px; margin-bottom: auto; cursor:hand; cursor:pointer;");
			
			img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
			{
				public void onEvent(Event event) throws Exception
				{
					Image img = (Image)event.getTarget();
					
					String strID = img.getId();
					String strEventID = strID.substring(13, strID.length());
					
					for(SimpleSpliceScore res : m_vcSplicingScores)
					{
						if(res.m_strID.equals(strEventID))
						{
							m_vcHighlightedJunctions.clear();
							m_vcHighlightedExons.clear();
							
							m_vcHighlightedJunctions.add(res.m_JunctionInclusion);
							m_vcHighlightedJunctions.add(res.m_JunctionExclusion);
							
							for(String strIsoform : m_Gene.getArrayOfGeneProductNames())
							{
								Exon pExons[] = m_Gene.getSortedExonsForGeneProduct(strIsoform);
								
								for(Exon ex : pExons)
								{
									if(ex.getCodingStart() == res.m_Exon.getCodingStart() || ex.getCodingStop() == res.m_Exon.getCodingStop())
										m_vcHighlightedExons.add(ex);
								}
							}

							m_treeSelectedIsoforms.clearSelection();
							
							for(Treeitem item : m_treeSelectedIsoforms.getItems())
							{
								component_loop:
								for(Component c : item.getChildren())
								{		
									Treerow row = (Treerow) c;
									
									for(String strIsoform : res.m_vcValidIsoforms)
									{
										String pSplit[] = row.getId().split("@");
										String strRefSeqID  = "?";
										String strEnsemblID = pSplit[0];
										
										if(pSplit.length > 1)
										{
											strRefSeqID = pSplit[1];
										}
										
										if(strEnsemblID.equals(strIsoform) || strRefSeqID.equals(strIsoform))
										{
											m_treeSelectedIsoforms.addItemToSelection(item);
											break component_loop;
										}
									}
								}
							}
							
							UpdateComboboxForIsoformSelection(false);
							DrawPlots();
						}
					}
				}
			});

			
			Hlayout layoutH = new Hlayout();
			layoutH.setParent(row);
			
			for(int i=0; i<5; i++)
			{
				img = null;
				
				if(i<res.m_nRating)
				{
					img = new Image("/img/good_rating.png");
				}
				else
				{
					img = new Image("/img/neutral_rating.png");
				}
				
				img.setHeight("15px");
				img.setParent(layoutH);
				img.setId("rating_" + (i+1) + "_" + res.m_strID);
				img.setStyle("cursor:hand; cursor:pointer;");
				
				img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
				{
					public void onEvent(Event event) throws Exception
					{
						Image img = (Image)event.getTarget();
						
						String strID = img.getId();
						
						String pSplit[] = strID.split("_");
						int nNewRating = Integer.parseInt(pSplit[1]);
						String strEventID = strID.substring(9, strID.length());
						
						int i=0;
						for(Component c : img.getParent().getChildren())
						{
							Image imgRating = (Image) c;
							imgRating.invalidate();
							
							if(i < nNewRating)
							{
								imgRating.setSrc("/img/good_rating.png");
							}
							else
							{
								imgRating.setSrc("/img/neutral_rating.png");
							}							
							i++;
						}
						
						for(SimpleSpliceScore res : m_vcSplicingScores)
						{
							if(res.m_strID.equals(strEventID))
							{
								res.m_nRating = nNewRating;
								img.getParent().getParent().setId("rating_" + (i+1) + "_" + res.m_strID);
								break;
							}
						}
					}
				});
			}	
			
			layoutH = new Hlayout();
			layoutH.setParent(row);
			
			img = new Image("/img/info.png");
			img.setHeight("20px");
			img.setId("info_" + res.m_strID);
			img.setStyle("display: block; margin-left: 4px; margin-right: auto; margin-top: 4px; margin-bottom: auto; cursor:hand; cursor:pointer;");
			img.setParent(layoutH);
			
			img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
			{
				public void onEvent(Event event) throws Exception
				{
					Image img = (Image)event.getTarget();
					
					String strID = img.getId();
					String strEventID = strID.substring(5, strID.length());

					for(SimpleSpliceScore res : m_vcSplicingScores)
					{
						if(res.m_strID.equals(strEventID))
						{
							Popup popup = new Popup();
							popup.setParent(img.getParent().getParent());
							popup.open(img.getParent().getParent(), "after_start");
							popup.setHflex("min");
							popup.setVflex("min");
							
							Grid grid = new Grid();
							grid.setMold("default");
							grid.setWidth("400px");
							grid.setVflex("min");
							grid.setParent(popup);
							
							Columns cols = new Columns();
							
							Column col = new Column("Attribute");
							col.setHflex("min");
							col.setParent(cols);

							col = new Column("Value");
							col.setHflex("min");
							col.setParent(cols);
						
							cols.setParent(grid);
							
							Rows rows = new Rows();
							rows.setParent(grid);
							
							//---------------------------------------------------------
							Row row = new Row();
							row.setParent(rows);
							
							Label label = new Label("Gene ID");
							label.setParent(row);
							
							label = new Label(m_Gene.getGeneID());
							label.setParent(row);
							
							//---------------------------------------------------------
							row = new Row();
							row.setParent(rows);
							
							label = new Label("Gene Symbol");
							label.setParent(row);
							
							label = new Label(m_Gene.getGeneName());
							label.setParent(row);
							
							//---------------------------------------------------------
							row = new Row();
							row.setParent(rows);
							
							label = new Label("Inclusion Junction");
							label.setParent(row);
							
							label = new Label(m_Gene.getChromosome() + ":" + res.m_JunctionInclusion.m_nStart + "-" + res.m_JunctionInclusion.m_nEnd);
							label.setParent(row);
							
							//---------------------------------------------------------
							row = new Row();
							row.setParent(rows);
							
							label = new Label("Exclusion Junction");
							label.setParent(row);
							
							label = new Label(m_Gene.getChromosome() + ":" + res.m_JunctionExclusion.m_nStart + "-" + res.m_JunctionExclusion.m_nEnd);
							label.setParent(row);
							
							//---------------------------------------------------------
							row = new Row();
							row.setParent(rows);
							
							label = new Label("Exon");
							label.setParent(row);
							
							label = new Label(m_Gene.getChromosome() + ":" + res.m_Exon.getCodingStart() + "-" + res.m_Exon.getCodingStop());
							label.setParent(row);
							
							//---------------------------------------------------------
							row = new Row();
							row.setParent(rows);
			
							label = new Label("Inclusion Change");
							label.setParent(row);
							
							row = new Row();
							row.setParent(rows);
							
							label = new Label("" + res.m_fInclusionChange);
							label.setParent(row);
							break;
						}
					}
				}
			});
			
			if(res.m_bIsNovel)
			{
				Cell cell = new Cell();
				cell.setParent(row);
				
				Html html = new Html("<span style=\"color:#F01414\"><b>yes</b>");
				html.setParent(cell);
			}
			else
			{
				Label label = new Label("no");
				label.setParent(row);
			}
			
			Label label = new Label(m_Gene.getGeneName());
			label.setParent(row);
			
			label = new Label(res.m_strConditionA);
			label.setParent(row);
			
			label = new Label(res.m_strConditionB);
			label.setParent(row);
			
			label = new Label(m_Gene.getChromosome() + ":" + res.m_Exon.getCodingStart() + "-" + res.m_Exon.getCodingStop());
			label.setParent(row);
			
			label = new Label(String.format(Locale.ENGLISH, "%.3f", res.m_fInclusionChange));
			label.setParent(row);
			
			label = new Label(String.format(Locale.ENGLISH, "%.3e", res.m_fPValue));
			label.setParent(row);
		}
		
		col.sort(true);
	}
	
	public void ShowHit(String strID) throws Exception
	{
		for(AlternativeSplicingHit hit : m_vcASHits)
		{
			if(hit.GetGeneID().equals(strID))
			{
				m_nMinCovPerBase 		 = hit.GetMinCovPerBase();
				m_fMinCoveredBases 		 = hit.GetMinCoveredBases();
				m_nMinJunctionReads 	 = hit.GetMinJunctionReads();
				m_fVariableExonThreshold = hit.GetVariableExonThreshold();
				
				m_textboxMinCovPerBase.setText(""+m_nMinCovPerBase);
				m_textboxMinCoveredBases.setText(""+m_fMinCoveredBases);
				m_textboxMinJunctionReads.setText(""+m_nMinJunctionReads);
				m_textboxVariableExonThreshold.setText(""+m_fVariableExonThreshold);
				
				m_bandboxSelectedGene.setValue(strID);
				
				if(!m_strFileGTF.equals(hit.GetGTFFile()))
				{
					m_strFileGTF = hit.GetGTFFile();
					LoadGeneAnnotation(m_strFileGTF);
					
					for(Comboitem item : m_comboboxSelectedGeneAnnotation.getItems())
					{
						File pFile = new File(m_strFileGTF);
						if(item.getLabel().equals(pFile.getName()))
							m_comboboxSelectedGeneAnnotation.setSelectedItem(item);
					}
				}
				
				OnGeneChange(false);
				
				m_treeSelectedIsoforms.clearSelection();
				
				for(Treeitem item : m_treeSelectedIsoforms.getItems())
				{
					component_loop:
					for(Component c : item.getChildren())
					{		
						Treerow row = (Treerow) c;
						
						for(String strIsoform : hit.GetIsoforms())
						{
							String pSplit[] = row.getId().split("@");
							String strRefSeqID  = "?";
							String strEnsemblID = pSplit[0];
							
							if(pSplit.length > 1)
							{
								strRefSeqID = pSplit[1];
							}
							
							if(strEnsemblID.equals(strIsoform) || strRefSeqID.equals(strIsoform))
							{
								m_treeSelectedIsoforms.addItemToSelection(item);
								break component_loop;
							}
						}
					}
				}
				
				UpdateComboboxForIsoformSelection(false);
				DrawPlots();
			}
		}
	}
	
	public void AddComboboxForSamples(Vlayout parentLayout)
	{
		if(m_treeSelectedSamples == null)
		{
			//####################################################
			//    add tree for sample selection to options box
			//####################################################
			m_treeSelectedSamples = new Tree();
			m_treeSelectedSamples.setId("sample_tree");
			m_treeSelectedSamples.setParent(parentLayout);
			m_treeSelectedSamples.setSizedByContent(true);
			m_treeSelectedSamples.setHeight("470px");
			m_treeSelectedSamples.setWidth("200px");
			m_treeSelectedSamples.setMultiple(true);
			m_treeSelectedSamples.setCheckmark(true);
	
			m_treeSelectedSamples.addEventListener(Events.ON_SELECT, new EventListener<Event>()
			{
				@Override
				public void onEvent(Event event) throws Exception
				{
					// changing the sample selection will disable the reduced data set option
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
												// check which items were newly selected/unselected
												TreeSet<String> vcSelectedItems = new TreeSet<String>();
												TreeSet<String> vcUnselectedItems = new TreeSet<String>();
												
												// get newly selected items
												for(Treeitem item : m_treeSelectedSamples.getSelectedItems())
												{
													if(!m_vcPreviouslySelectedTreeItems.contains(item.getLabel()))
														vcSelectedItems.add(item.getLabel());
												}
												
												// get newly unselected items
												for(String strItem : m_vcPreviouslySelectedTreeItems)
												{
													boolean bIncluded = false;
													for(Treeitem item : m_treeSelectedSamples.getSelectedItems())
													{
														if(item.getLabel().equals(strItem))
															bIncluded = true;
													}
													if(!bIncluded)
														vcUnselectedItems.add(strItem);
												}
												
												m_bUseReducedDataSet = false;
												m_checkboxUseReducedDataSet.setChecked(false);
												ProcessNewSampleSelection(vcSelectedItems, vcUnselectedItems);
												break;
											}
											
											case Messagebox.CANCEL:
											{
												ProcessNewSampleSelection(m_vcPreviouslySelectedTreeItems, null);
												break;
											}
										}
									}
								});
					}
					else
					{
						// check which items were newly selected/unselected
						final TreeSet<String> vcSelectedItems = new TreeSet<String>();
						final TreeSet<String> vcUnselectedItems = new TreeSet<String>();
						
						// get newly selected items
						for(Treeitem item : m_treeSelectedSamples.getSelectedItems())
						{
							if(!m_vcPreviouslySelectedTreeItems.contains(item.getLabel()))
								vcSelectedItems.add(item.getLabel());
						}
						
						// get newly unselected items
						for(String strItem : m_vcPreviouslySelectedTreeItems)
						{
							boolean bIncluded = false;
							for(Treeitem item : m_treeSelectedSamples.getSelectedItems())
							{
								if(item.getLabel().equals(strItem))
									bIncluded = true;
							}
							if(!bIncluded)
								vcUnselectedItems.add(strItem);
						}
						
						ProcessNewSampleSelection(vcSelectedItems, vcUnselectedItems);
					}
				}
			});
			
			// add button to unselect all samples
			Button btnUnselectAll = new Button("Unselect all samples");
			btnUnselectAll.setStyle("margin-left: 30px;");
			btnUnselectAll.setParent(parentLayout);

			btnUnselectAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
			{
				@Override
				public void onEvent(Event event) throws Exception
				{
					TreeSet<String> vcSelectedItems = new TreeSet<String>();
					TreeSet<String> vcUnselectedItems = new TreeSet<String>();
					
					for(String strCondition : m_projectModel.GetConditions(m_strSelectedConditionType))
						vcUnselectedItems.add(strCondition);

					ProcessNewSampleSelection(vcSelectedItems, vcUnselectedItems);
					UpdateComboboxForSampleSelection(true);
					
					m_vcPreviouslySelectedTreeItems.clear();
					
					// uncheck the "use reduced data set" checkbox
					m_checkboxUseReducedDataSet.setChecked(false);
					m_bUseReducedDataSet = false;
				}
			});
		}
	}

	public void ProcessNewSampleSelection(TreeSet<String> vcSelectedItems, TreeSet<String> vcUnselectedItems)
	{
		if(vcSelectedItems == null)
			return;
		
		m_vcSelectedSamples.clear();

		TreeMap<String, TreeSet<String>> mapSamplesToConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			if(vcSelectedItems.contains(strCondition))
			{
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					for(Treeitem item : m_treeSelectedSamples.getItems())
					{
						if(item.getLabel().equals(strSample))
							m_treeSelectedSamples.addItemToSelection(item);
					}
				}
			}
			else if(vcUnselectedItems.contains(strCondition))
			{
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					for(Treeitem item : m_treeSelectedSamples.getItems())
					{
						if(item.getLabel().equals(strSample))
							m_treeSelectedSamples.removeItemFromSelection(item);
					}
				}
			}
		}

		// update selected samples and previously selected items (items include conditions!)
		m_vcPreviouslySelectedTreeItems.clear();
		for(Treeitem item : m_treeSelectedSamples.getSelectedItems())
		{
			m_vcPreviouslySelectedTreeItems.add(item.getLabel());

			if(m_projectModel.GetSamples().contains(item.getLabel()))
				m_vcSelectedSamples.add(item.getLabel());
		}
		
		m_bCoverageRequiresRedraw = true;
		m_bSampleSelectionChanged = true;
		
		UpdateColorSelection();
	}
	
	public void UpdateComboboxForSampleSelection(boolean bSamplesChanged)
	{
		if(m_treeSelectedSamples == null)
			return;
		
		m_bCoverageRequiresRedraw = true;
		
		// clear old samples
		m_treeSelectedSamples.getChildren().clear();
		
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
		
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			Treeitem conditionItem = new Treeitem();
			conditionItem.setParent(children);
			conditionItem.setLabel(strCondition);
			
			Treechildren conditionChildren = new Treechildren();
			conditionChildren.setParent(conditionItem);
			
			int nSamplesSelected = 0;
			
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				Treeitem sampleItem = new Treeitem();
				sampleItem.setParent(conditionChildren);
				
				if(m_vcSelectedSamples.contains(strSample))
				{
					m_treeSelectedSamples.addItemToSelection(sampleItem);
					nSamplesSelected++;
				}
				
				Treerow row = new Treerow();
				row.setParent(sampleItem);
				
				Treecell cell = new Treecell(strSample); cell.setParent(row);
			}
			
			if(nSamplesSelected == 0)
			{
				conditionItem.setOpen(false);
				m_treeSelectedSamples.removeItemFromSelection(conditionItem);
			}
			else if(mapSamplesToConditions.get(strCondition).size() == nSamplesSelected)
			{
				conditionItem.setOpen(false);
//				conditionItem.setSelected(true);
				m_vcPreviouslySelectedTreeItems.add(conditionItem.getLabel());
				m_treeSelectedSamples.addItemToSelection(conditionItem);
			}
			else
			{
				if(m_vcPreviouslySelectedTreeItems.contains(strCondition))
				{
					m_vcPreviouslySelectedTreeItems.remove(strCondition);
				}
			}
		}
		
		if(bSamplesChanged && m_Gene != null)
		{
			System.out.println("recalculating exon values");
			// get coverage for all samples
			GetBigWigCoverageForGene();
			
			// recalculate expression values based on selected samples
			try
			{
				ProcessCountData();
			} catch (IOException e)
			{
				Messagebox.show("failed to recalculate exon coverage");
				e.printStackTrace();
			}
		}
		
		m_treeSelectedSamples.setSizedByContent(true);
	}

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
			m_comboboxSelectedConditionType.setWidth("180px");
			
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
			m_comboboxSelectedCondition.setWidth("180px");
			m_comboboxSelectedCondition.setParent(layoutV);
			m_comboboxSelectedCondition.addEventListener(Events.ON_SELECT, new EventListener<Event>()
			{
				public void onEvent(Event event) throws Exception
				{
					Combobox box = (Combobox)event.getTarget();
					if(box.getSelectedItem() != null)
					{
						m_strSelectedCondition = box.getSelectedItem().getValue();
						
						if(m_Gene != null)
							DrawPlots();
					}
				}
			});
		}
	}
	
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
	
	public void AddComboboxForIsoformSelection(Vlayout parentLayout) throws IOException
	{
		if(m_treeSelectedIsoforms == null)
		{
			m_treeSelectedIsoforms = new Tree();
			m_treeSelectedIsoforms.setId("isoform_tree");
			m_treeSelectedIsoforms.setParent(parentLayout);
			m_treeSelectedIsoforms.setWidth("400px");
			m_treeSelectedIsoforms.setHeight("470px");
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
					m_bIsoformSelectionChanged = true;
				}
			});
			
			// add buttons
			Hlayout layoutH = new Hlayout();
			layoutH.setHflex("min");
			layoutH.setStyle("text-align: center; valign: middle;");			
			layoutH.setParent(parentLayout);
			
			// add button to unselect all isoforms
			Button btnUnselectAll = new Button("Unselect all isoforms");
			btnUnselectAll.setStyle("margin-left: 10px;");
			btnUnselectAll.setParent(layoutH);

			btnUnselectAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
			{
				@Override
				public void onEvent(Event event) throws Exception
				{
					m_vcValidIsoforms.clear();
					UpdateComboboxForIsoformSelection(true);
					
					m_bIsoformSelectionChanged = true;
				}
			});
			
			// add button to unselect short isoforms
			Button btnUnselectSmallIsoforms = new Button("Unselect all isoforms <=");
			btnUnselectSmallIsoforms.setStyle("margin-left: 5px;");
			btnUnselectSmallIsoforms.setParent(layoutH);

			btnUnselectSmallIsoforms.addEventListener(Events.ON_CLICK, new EventListener<Event>()
			{
				@Override
				public void onEvent(Event event) throws Exception
				{
					int nThreshold = Integer.parseInt(m_textboxMinExonThreshold.getText());
					TreeSet<String> strUnSelectedIsoforms = new TreeSet<String>();
					
					String[] pIsoforms = m_Gene.getArrayOfGeneProductNames();

					// find isoforms that should be unselected
					for(String strIsoform : m_vcValidIsoforms)
					{
						for(String strGeneIsoforms : pIsoforms)
						{
							if(strGeneIsoforms.split("\\.")[0].equals(strIsoform))
							{
								if(m_Gene.getSortedExonsForGeneProduct(strGeneIsoforms).length <= nThreshold)
									strUnSelectedIsoforms.add(strIsoform);
							}
						}
					}
					
					// unselect isoforms
					for(String strIsoform : strUnSelectedIsoforms)
						m_vcValidIsoforms.remove(strIsoform);
						
					UpdateComboboxForIsoformSelection(true);
					
					m_bIsoformSelectionChanged = true;
				}
			});
		
			// add textbox to control which isoforms are too small
			m_textboxMinExonThreshold = new Textbox("3");
			m_textboxMinExonThreshold.setConstraint("/^[0-9]+$/");
			m_textboxMinExonThreshold.setParent(layoutH);
			m_textboxMinExonThreshold.setWidth("20px");
			
			Label lblText = new Label("exons");
			lblText.setStyle("display:inline-block; margin-top: 6px;");
			lblText.setParent(layoutH);
		}
	}
	
	public void UpdateComboboxForIsoformSelection(boolean bNewGene) throws IOException
	{
		if(m_treeSelectedIsoforms == null)
			return;
		
		if(!bNewGene)
		{
			GetValidIsoforms();
			
			if(m_bHideIrrelevantIsoforms)
			{
				HideIrrelevantIsoforms(m_bSkipFirstAndLastExon, false);
				m_bHideIrrelevantIsoforms = false;
			}
		}

		m_bCoverageRequiresRedraw = true;
		
		m_treeSelectedIsoforms.getChildren().clear();

		if(m_Gene != null)
		{
			Treecols cols = new Treecols();
			Treecol col = new Treecol("Ensembl ID"); col.setParent(cols);
			col = new Treecol("RefSeq ID"); col.setParent(cols);
			col = new Treecol("Length [bp]"); col.setParent(cols);
			col = new Treecol("Exons [N]"); col.setParent(cols);
			cols.setParent(m_treeSelectedIsoforms);
			
			Treechildren children = new Treechildren();
			children.setParent(m_treeSelectedIsoforms);
			
			m_treeSelectedIsoforms.invalidate();
			
			TreeSet<String> vcIdentifierUsed = new TreeSet<String>();
			
			for(String strIsoform : m_mapExonsToIsoforms.keySet())
			{
				String strIsoformID = strIsoform.split("\\.")[0];

				GeneIdentifier gid = GetGeneIdentifierForGene(strIsoform);
				if(gid != null)
				{					
					String strIdentifier = gid.m_strEnsemblTranscriptID + "@" + gid.m_strRefGeneID + "@" + gid.m_strEntrezGeneID;

					if(!vcIdentifierUsed.contains(strIdentifier))
					{
						Treeitem item = new Treeitem();
						item.setParent(children);
						item.setCheckable(true);
						item.setValue(strIsoformID);
						
						Treerow row = new Treerow();
						row.setParent(item);
						
						vcIdentifierUsed.add(strIdentifier);
						row.setId(gid.m_strEnsemblTranscriptID + "@" + gid.m_strRefGeneID + "@" + gid.m_strEntrezGeneID);
						Treecell cell = new Treecell(gid.m_strEnsemblTranscriptID); cell.setParent(row);
						cell = new Treecell(gid.m_strRefGeneID); cell.setParent(row);
						
						Exon pExons[] = m_Gene.getSortedExonsForGeneProduct(strIsoform, true);
	
						int nExons = pExons.length;
						int nLength = 0;
						for(Exon ex : pExons)
						{
							nLength += ex.getLength();
						}
	
						cell = new Treecell(String.format(Locale.ENGLISH, "%d", nLength)); cell.setParent(row);
						cell = new Treecell(String.format(Locale.ENGLISH, "%d", nExons)); cell.setParent(row);
						
						// add to selection
						if(m_vcValidIsoforms.contains(strIsoformID))
							m_treeSelectedIsoforms.addItemToSelection(item);
					}
				}
				else
				{
					Treeitem item = new Treeitem();
					item.setParent(children);
					item.setCheckable(true);
					item.setValue(strIsoformID);
					
					Treerow row = new Treerow();
					row.setParent(item);
					
					Treecell cell = new Treecell(strIsoformID); cell.setParent(row);
					cell = new Treecell("?"); cell.setParent(row);
					cell = new Treecell("?"); cell.setParent(row);
					cell = new Treecell("?"); cell.setParent(row);
					
					row.setId(strIsoformID);
					
					// add to selection
					if(m_vcValidIsoforms.contains(item.getLabel()))
						m_treeSelectedIsoforms.addItemToSelection(item);
				}
			}
			
			m_treeSelectedIsoforms.setSizedByContent(true);
		}
	}
	
	public AlternativeSplicingHit CalculateVariableExons() throws IOException
	{
		AlternativeSplicingHit result = null;
		
		TreeMap<String, String> mapBigWigFiles 	= m_projectModel.GetBigWigFilesForSamples();
		
		if(mapBigWigFiles == null)
		{
			Messagebox.show("ERROR: no valid bigwig or bam files detected");
			return null;
		}
		
		TreeMap<String, Double> mapSizeFactorsToSampleNames = m_projectModel.GetSizeFactors();
		TreeMap<Integer, Double> mapSizeFactorsToSampleIDs 	= new TreeMap<Integer, Double>();
		
		// create new set of exon groups excluding all first and last exons
		//TODO alternative first + last exons
		TreeMap<String, int[]> mapExonsToIsoforms = new TreeMap<String, int[]>();
		
		for(String strIsoform : m_mapExonsToIsoforms.keySet())
		{
			if(m_vcValidIsoforms.contains(strIsoform.split("\\.")[0]))
			{
				int[] pIsoformExons = m_mapExonsToIsoforms.get(strIsoform);
				
				if(pIsoformExons.length > 2)
				{
					int[] pExonIDs = new int[pIsoformExons.length-2];
					
					for(int i=1; i<pIsoformExons.length-1; i++)
					{
						pExonIDs[i-1] = pIsoformExons[i];
					}
					
					mapExonsToIsoforms.put(strIsoform, pExonIDs);
				}
			}
		}
		
		ExonGroup[] pExonGroups = RecalculateExonGroups(mapExonsToIsoforms);
		
		// get reference name
		String strRef = m_Gene.getChromosome();
		
		// get samples per condition
		final TreeMap<String, TreeSet<String>> vcSamplesAndConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		int nConditions = vcSamplesAndConditions.size();
		
		//########################################
		//      prepare coverage containers
		//########################################
		// key 1 = exon_name, key 2 = sample, value = coverage array
		TreeMap<String, TreeMap<String, double[]>> mapCoverageToSamplesAndExons = new TreeMap<String, TreeMap<String, double[]>>();
		for(ExonGroup grp : pExonGroups)
		{			
			int nStart	= grp.getGenomicStartOfGroup();
			int nEnd	= grp.getGenomicStopOfGroup();
			
			String strName = strRef + ":" + nStart + "-" + nEnd;
			
			mapCoverageToSamplesAndExons.put(strName, new TreeMap<String, double[]>());
		}
		
		//###########################################
		//    obtain coverage for all exon groups
		//###########################################
		for(String strCondition : vcSamplesAndConditions.keySet())
		{
			int nSampleID = 0;
			for(String strSample : vcSamplesAndConditions.get(strCondition))
			{
				if(!m_vcSelectedSamples.contains(strSample))
					continue;
	
				for(ExonGroup grp : pExonGroups)
				{
					int nStart	= grp.getGenomicStartOfGroup();
					int nEnd	= grp.getGenomicStopOfGroup();
					String strGrpName = strRef + ":" + nStart + "-" + nEnd;
					
					double[] pCoverage = GetCoverageForExonGroup(grp, strSample, true);
					
					mapCoverageToSamplesAndExons.get(strGrpName).put(strSample, pCoverage);
				}
				
				mapSizeFactorsToSampleIDs.put(nSampleID, mapSizeFactorsToSampleNames.get(strSample));
				nSampleID++;
			}
		}
		
		// identify conditions with insufficient coverage for all of its exons
		TreeSet<String> vcInvalidConditions = new TreeSet<String>();
		for(String strCondition : vcSamplesAndConditions.keySet())
		{
			double fMaxCoverage = 0.0;
			
			for(String strGrp : mapCoverageToSamplesAndExons.keySet())
			{
				Vector<Double> vcValues = new Vector<Double>();
				
				for(String strSample : vcSamplesAndConditions.get(strCondition))
				{
					if(!m_vcSelectedSamples.contains(strSample))
						continue;
					
					// get the mean coverage for the current exon group and sample
					vcValues.add(StatUtils.mean(mapCoverageToSamplesAndExons.get(strGrp).get(strSample)));
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
			
			if(fMaxCoverage < m_nMinCovPerBase)
			{
				vcInvalidConditions.add(strCondition);
			}
		}
		
		// get gene identifier
		GeneIdentifier gid = GetGeneIdentifierForGene(m_Gene.getGeneID());
		if(gid == null)
		{
			gid = GetGeneIdentifierForGene(m_Gene.getGeneName());
			
			if(gid == null)
			{
				gid = new GeneIdentifier();
				gid.m_strEntrezGeneID		= m_Gene.getGeneID();
				gid.m_strApprovedGeneName	= m_Gene.getGeneID();
				gid.m_strApprovedGeneSymbol = m_Gene.getGeneName();
			}
		}
		
		// prepare result
		result = new AlternativeSplicingHit(0, gid, m_Gene.getChromosome(), m_nMinJunctionReads, m_nMinCovPerBase,
				m_fMinCoveredBases, m_fVariableExonThreshold, m_vcValidIsoforms, m_strFileGTF, "");
		
		// prepare coverage arrays for all conditions
		
		// key1 = exon group, key2 = condition, value = coverage per position for each sample
		TreeMap<ExonGroup, TreeMap<String, double[][]>> mapCoveragePerSamplePerCondition = new TreeMap<ExonGroup, TreeMap<String, double[][]>>();

		for(ExonGroup grp : pExonGroups)
		{
			String strName = strRef + ":" + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup();
			TreeMap<String, double[]> mapData = mapCoverageToSamplesAndExons.get(strName);
	
			for(String strCondition : vcSamplesAndConditions.keySet())
			{	
				int nValidSamples = 0;
				for(String strSample : vcSamplesAndConditions.get(strCondition))
				{
					if(m_vcSelectedSamples.contains(strSample))
						nValidSamples += 1;
				}

				if(nValidSamples > 0)
				{
					double[][] pValues =  new double[grp.getExonGroupLengthInBp()][nValidSamples];
					
					int nCurrentSample = 0;
					for(String strSample : vcSamplesAndConditions.get(strCondition))
					{
						if(!m_vcSelectedSamples.contains(strSample))
							continue;
						
						double pData[] =  mapData.get(strSample);						
						for(int x=0; x<pData.length; x++)
						{
							double fValue = pData[x];
							pValues[x][nCurrentSample] = fValue;
						}
						
						nCurrentSample += 1;
					}
					
					if(nCurrentSample != 0)
					{
						if(mapCoveragePerSamplePerCondition.containsKey(grp))
						{
							TreeMap<String, double[][]> mapTmp = mapCoveragePerSamplePerCondition.get(grp);
							mapTmp.put(strCondition, pValues);
						}
						else
						{
							TreeMap<String, double[][]> mapTmp = new TreeMap<String, double[][]>();
							mapTmp.put(strCondition, pValues);
							mapCoveragePerSamplePerCondition.put(grp, mapTmp);
						}
					}
				}
				else
				{
					vcInvalidConditions.add(strCondition);
				}
			}
		}
		
		// test for all condition combinations		
		String pConditions[] = new String[nConditions];
		vcSamplesAndConditions.keySet().toArray(pConditions);
		
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
		
				// prepare fraction map for exon groups and exons
				TreeMap<ExonGroup, double[]> mapCovFractionsPerExonGroup = new TreeMap<ExonGroup, double[]>();
				TreeMap<Exon, double[]> mapCovFractionsPerExon		 	 = new TreeMap<Exon, double[]>();
				
				// calculate fractions for the selected conditions
				for(ExonGroup grp : pExonGroups)
				{
					Vector<Double> vcFractions = new Vector<Double>();
					
					double pValuesA[][] = mapCoveragePerSamplePerCondition.get(grp).get(strConditionA);
					double pValuesB[][] = mapCoveragePerSamplePerCondition.get(grp).get(strConditionB);

					// calculate for exon groups
					for(int x=0; x<grp.getExonGroupLengthInBp(); x++)
					{
						//#################################
						//   first, calculate the sum
						//#################################
						// get mean coverage for this exon group in condition A
						
						double pCov[] = new double[pValuesA[0].length];
						
						for(int nSample = 0; nSample<pValuesA[0].length; nSample++)
							pCov[nSample] = pValuesA[x][nSample]; 
		
						double fMeanCoverageA = StatUtils.mean(pCov);
						
						// get median coverage for this exon group in condition B						
						pCov = new double[pValuesB[0].length];
						
						for(int nSample = 0; nSample<pValuesB[0].length; nSample++)
							pCov[nSample] = pValuesB[x][nSample]; 

						double fMeanCoverageB = StatUtils.mean(pCov);

						// get total coverage
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
						pFractions = new double[ex.getGenomicLength()];
						
						int nStart = ex.getCodingStart() - grp.getGenomicStartOfGroup();
						int nEnd   = nStart + ex.getGenomicLength();
						
						for(int x=nStart; x<nEnd; x++)
						{
							//#################################
							//   first, calculate the sum
							//#################################
							// get mean coverage for this exon group in condition A
							
							double pCov[] = new double[pValuesA[0].length];
							
							for(int nSample = 0; nSample<pValuesA[0].length; nSample++)
							{
								pCov[nSample] = pValuesA[x][nSample]; 
							}

							double fMeanCoverageA = StatUtils.mean(pCov);
							
							// get median coverage for this exon group in condition B						
							pCov = new double[pValuesB[0].length];
							
							for(int nSample = 0; nSample<pValuesB[0].length; nSample++)
							{
								pCov[nSample] = pValuesB[x][nSample]; 
							}

							double fMeanCoverageB = StatUtils.mean(pCov);

							// get total coverage
							double fTotalCoverage   = fMeanCoverageA + fMeanCoverageB;

							if(fTotalCoverage == 0)
								continue;
							
							// calculate ratio
							pFractions[x-nStart] = fMeanCoverageA / fTotalCoverage;
						}
						
						mapCovFractionsPerExon.put(ex, pFractions);
					}
				}

				// perform test for each exon group to identify alternatively spliced exons
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
					
					if(ex.getCodingStart() == 203702351 && ex.getCodingStop() == 203702528 && strConditionA.equals("GBM"))
					{
						System.out.println("testing: " + strConditionA + " vs " + strConditionB);
						System.out.println(Arrays.toString(pFractionsGroup1));
						System.out.println(Arrays.toString(pAllGroups));						
						System.out.println("p-value:" + fPValue);
						System.out.println("abs. change: " + fAbsDifference);
						System.out.println("threshold: " + m_fVariableExonThreshold);
					}
					
					if(fPValue < 0.05 && Math.abs(fAbsDifference) >= m_fVariableExonThreshold)
					{
						if(ex.getCodingStart() == 203702351 && ex.getCodingStop() == 203702528 && strConditionA.equals("GBM"))
							System.out.println("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$");
						double fCoveragePerBaseA = 0.0;
						double fCoveragePerBaseB = 0.0;
						
						for(ExonGroup grp : pExonGroups)
						{
							if(grp.groupContainsExon(ex.getExonID()))
							{
								double pValuesA[][] = mapCoveragePerSamplePerCondition.get(grp).get(strConditionA);
								double pValuesB[][] = mapCoveragePerSamplePerCondition.get(grp).get(strConditionB);
								
								int nStart = ex.getCodingStart() - grp.getGenomicStartOfGroup();
								int nEnd   = nStart + ex.getGenomicLength();
								
								for(int x=nStart; x<nEnd; x++)
								{
									double pCovA[] = new double[pValuesA[0].length];
									
									for(int nSample = 0; nSample<pValuesA[0].length; nSample++)
										pCovA[nSample] = pValuesA[x][nSample] * mapSizeFactorsToSampleIDs.get(nSample);
									
									fCoveragePerBaseA += StatUtils.mean(pCovA);
								}
								
								for(int x=nStart; x<nEnd; x++)
								{
									double pCovB[] = new double[pValuesB[0].length];
									
									for(int nSample = 0; nSample<pValuesB[0].length; nSample++)
										pCovB[nSample] = pValuesB[x][nSample] * mapSizeFactorsToSampleIDs.get(nSample);
									
									fCoveragePerBaseB += StatUtils.mean(pCovB);
								}
								
								fCoveragePerBaseA /= (nEnd-nStart+1);
								fCoveragePerBaseB /= (nEnd-nStart+1);
								break;
							}
						}
						
						if(ex.getCodingStart() == 203702351 && ex.getCodingStop() == 203702528)
						{
							System.out.println("adding result for: " + strConditionA + " vs " + strConditionB);
						}
						
						result.AddASExons(m_Gene.getChromosome() + ":" + ex.getCodingStart() + "-" + ex.getCodingStop(), ex.getCodingStart(), ex.getCodingStop(),
								strConditionA, strConditionB, fCoveragePerBaseA, fCoveragePerBaseB, fAbsDifference, fRelDifference, fPValue, pFractionsGroup1, pAllGroups);
					}
				}
			}
		}
		
//		System.out.println(result);
		
		return result;
	}
	
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
		TreeMap<String, TreeMap<String, Integer>> mapJunctionReadCounts =  m_projectModel.GetJunctionCountsForGene(m_Gene);

		if(mapJunctionReadCounts == null)
		{
			System.out.println("failed to get junction read counts for gene: " + m_Gene.getGeneID() + " (" + m_Gene.getGeneName() + ")");
			return null;
		}
		
		ExecutorService executor = Executors.newFixedThreadPool(m_nThreads);
		for(String strIsoform : m_Gene.getArrayOfGeneProductNames())
		{
			if(!m_vcValidIsoforms.contains(strIsoform.split("\\.")[0]))
				continue;
			
			// start a thread for each isoform
			if(nMode == 0)
			{
				Runnable r = new ThreadRemoveIrrelevantIsoformsBasedOnSplitReads(strIsoform, strCondition, vcResult, mapSamplesToConditions, mapJunctionReadCounts, bSkipFirstAndLastExon, bDebug);
				executor.execute(r);
			}
			else
			{
				Runnable r = new ThreadRemoveIrrelevantIsoformsBasedOnCoverage(strIsoform, strCondition, vcResult, mapSamplesToConditions, mapJunctionReadCounts, bSkipFirstAndLastExon, bDebug);
				executor.execute(r);
			}
		}
		executor.shutdown();
		while(!executor.isTerminated()) {}
		
		return vcResult;
	}

	//########################################################################
	//    Calculates the isoform expression value per sample
	//    Returns a map of the form <isoform_id, <sample, coverage_value>>
	//########################################################################
	public void GenerateMMSeqEstimates() throws Exception
	{	
		DateFormat dateFormat = new SimpleDateFormat("yyyy_MM_dd_HH_mm_ss");
		Calendar cal = Calendar.getInstance();
		String strTimeStamp = dateFormat.format(cal.getTime());
		
		TreeMap<String, TreeSet<String>> vcSamplesPerCondition = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		TreeSet<String> vcIsoforms 	= new TreeSet<String>();
		String[] pIsoforms 			= m_Gene.getArrayOfGeneProductNames();

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
			Exon pExons[] = m_Gene.getSortedExonsForGeneProduct(strIsoform, true);
			int nGeneStart = m_Gene.getStart();
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
//		System.out.println(mapIsoformsToPositions);
		
		// <Isoform, <Sample, Value>>
		TreeMap<String, TreeMap<String, Double>> mapRelativeExpressionToIsoform = new TreeMap<String, TreeMap<String, Double>>();
		
		for(String strCondition : vcSamplesPerCondition.keySet())
		{
			for(String strSample : vcSamplesPerCondition.get(strCondition))
			{
				// get coverage for the current sample
				double pCoverage[] = GetCoverageForRegion(m_Gene.getChromosome(), m_Gene.getStart(), m_Gene.getStop(), strSample, true);
				
				//######################################
				//          prepare hits file
				//######################################
				String strFileHits = m_strTmpFolder + "/" + strTimeStamp + "_" + strSample + "_" + m_Gene.getGeneID() + "_hits";
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
				pWriter.print(m_Gene.getGeneID());
					
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
					
					if(mapRelativeExpressionToIsoform.containsKey(strIsoform.split("\\.")[0]))
					{
						TreeMap<String, Double> mapExpressionToSample = mapRelativeExpressionToIsoform.get(strIsoform.split("\\.")[0]);
						mapExpressionToSample.put(strSample, fValue);
					}
					else
					{
						TreeMap<String, Double> mapExpressionToSample = new TreeMap<String, Double>();
						mapExpressionToSample.put(strSample, fValue);
						mapRelativeExpressionToIsoform.put(strIsoform.split("\\.")[0], mapExpressionToSample);
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

		m_projectModel.AddIsoformExpressionToDataBase(mapRelativeExpressionToIsoform, m_Gene.getGeneID(), m_strPathMMSeqResults);
	}

	//########################################################################
	//                Calculates the sample size factors
	//########################################################################
	public boolean CalculateSizeFactors(String strFileGTF, String strFileProject) throws IOException
	{
		// init project
		try
		{
			if(!m_projectModel.Init(strFileProject, -1, -1, -1.0, -1, true))
				return false;
		}
		catch(Exception e)
		{
			System.out.println("failed to process project file: " + strFileProject);
			System.out.println(ExceptionUtils.getStackTrace(e));
			return false;
		}
		
		TreeMap<String, String> mapBigWigFiles 	= m_projectModel.GetBigWigFilesForSamples();
		TreeSet<String> vcSamples = m_projectModel.GetSamples(); 
		int nSamples = vcSamples.size();

		if(this.m_nReferenceType == REFFLAT_REFERENCE_FILE)
		{
			System.out.println("can't calculate size factors with refflat file -> GTF format required");
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
	
	//##############################################################################################
	//    selected isoforms are those isoforms that are actually drawn and used for calculations
	//##############################################################################################
	public void GetValidIsoforms() throws IOException
	{
		//######################################################################
		//    Remove isoforms that have are not compatible with the selected
		//    exons from the list of available isoforms.
		//    All isoforms must contain exon group A if it has been selected.
		//######################################################################
		if(m_ClickEvent.m_strClickedA != null)
		{
			if(m_ClickEvent.m_strClickedA.contains("exon_group"))
			{
				m_vcValidIsoforms.clear();
				TreeMap<String, int[]> mapExonsToIsoforms = GetSelectedIsoforms();
				
				for(String strIsoform : mapExonsToIsoforms.keySet())
					m_vcValidIsoforms.add(strIsoform);
			}
			else if(m_ClickEvent.m_strClickedA.contains("edge"))
			{
				String pSplitA[] = m_ClickEvent.m_strClickedA.split("_");
			
				HashMap<Integer, Integer> mapExonsToExonGroups = m_Gene.mapExonsToGroups(m_pExonGroups, m_pExons);
				
				int nClickedEdgeStart = Integer.parseInt(pSplitA[1])-1;
				int nClickedEdgeEnd   = Integer.parseInt(pSplitA[2])-1;
				
				Iterator<String> it = m_Gene.getGeneProductNames();
				while(it.hasNext())
				{						
					String strGeneProductID = it.next();
					int[] pExonIDs = m_Gene.getExonsForGeneProduct(strGeneProductID);
					
					Vector<Integer> vcExonGroupIDs = new Vector<Integer>();
					// get exon group for each exon
					for(int nID : pExonIDs)
					{
						if(mapExonsToExonGroups.containsKey(nID))
						{
							int nExonGrpID = mapExonsToExonGroups.get(nID);
							vcExonGroupIDs.add(nExonGrpID);
						}
					}
					
					for(int i=0; i<vcExonGroupIDs.size()-1; i++)
					{
						if(vcExonGroupIDs.get(i) == nClickedEdgeStart && vcExonGroupIDs.get(i+1) == nClickedEdgeEnd)
							m_vcValidIsoforms.add(strGeneProductID);	
					}
				}
			}
		}
		else
		{
			//##############################################
			//    get all selected isoforms for the gene
			//##############################################
			m_vcValidIsoforms.clear();
			for(Treeitem item : m_treeSelectedIsoforms.getSelectedItems())
			{
				m_vcValidIsoforms.add((String)item.getValue());
			}
		}

		OnIsoformChange();
	}

	public void GetColorsPerCondition()
	{	
		if(m_mapColorsToConditions == null)
			m_mapColorsToConditions = new TreeMap<String, Color>();
		else
			m_mapColorsToConditions.clear();

		TreeMap<String, TreeSet<String>> vcSamplesAndConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		// get a unique color per condition		
		int nCondition = 0;
		for(String strCondition : vcSamplesAndConditions.keySet())
		{
			if(nCondition < 30)
				m_mapColorsToConditions.put(strCondition, m_pColors[nCondition]);
			else
			{
				Random rand = new Random();
				int r = rand.nextInt(255);
				int g = rand.nextInt(255);
				int b = rand.nextInt(255);
				m_mapColorsToConditions.put(strCondition, new Color(r, g, b));
			}
			nCondition+=1;
		}
	}
	
	public void GetBigWigCoverageForGene()
	{		
		if(m_Gene == null)
		{
			Messagebox.show("GetBigWigCoverageForGene - No gene selected");
			return;
		}
		
		if(m_mapBigWigCoverageToSamples == null)
			m_mapBigWigCoverageToSamples = new TreeMap<String, int[]>();
		else
			m_mapBigWigCoverageToSamples.clear();
		
		TreeMap<String, String> mapBigWigFiles 	= m_projectModel.GetBigWigFilesForSamples();
		TreeSet<String> vcSamples 				= m_vcSelectedSamples; //m_projectModel.GetSamples();
		
		/*
		ExecutorService executor = Executors.newFixedThreadPool(m_nThreads);
		
		// start a thread for each sample
		for(String strSample : vcSamples)
		{
			Runnable r = new ThreadReadBigWigCoverage(m_Gene, strSample, mapBigWigFiles.get(strSample), m_mapBigWigCoverageToSamples);
			executor.execute(r);
		}
		
		executor.shutdown();
		
		while(!executor.isTerminated()) {}
		*/
		
		int nStart 	= m_Gene.getStart();
		int nEnd	= m_Gene.getStop();
	
		for(String strSample : vcSamples)
		{
			BigWigReader reader = null;
			try
			{
				reader = new BigWigReader(mapBigWigFiles.get(strSample));
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
			
			m_mapBigWigCoverageToSamples.put(strSample, pValues);
		}
	}
		
	// Some exons share the same coverage value because they are too similar to be discriminated.
	// This function returns the coverage value for a given exon (or the ambigious counter part) for all samples
	public double[] GetFinalCoverageForExon(Exon exon, TreeMap<String, Double> mapExonCoverageToSample)
	{
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		int nSamples = 0;
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				if(m_vcSelectedSamples.contains(strSample))
					nSamples++;
			}
		}

		//#####################################################################################################
		// Identify which exon has the coverage value for the given exon for each sample
		// This means, if most samples split the original exon into three sub exons, and some samples
		// split the original exon in just two sub exons use the 3-part separation.
		//#####################################################################################################
		TreeMap<String, TreeSet<Exon>> vcAmbigiousExonSetsPerCondition = GetMostFrequentSetOfAmbigiousExons(exon);
		
		//#####################################################################################################
		// Now we have list of ambigious exons per sample, but we lost track of the target exon of the group
		// Thus, we need to get the target exon for each ambigious exon set.
		// (where target means 'key' to the indistuingishible exon set)
		//#####################################################################################################
		TreeMap<String, Exon> mapTargetExonsPerSample = new TreeMap<String, Exon>();
		
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			// get the set of ambigious exons for the current condition
			TreeSet<Exon> vcTarget = vcAmbigiousExonSetsPerCondition.get(strCondition);
			
			// now find the corresponding target exon
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				if(m_vcSelectedSamples.contains(strSample))
				{
					// check all target exons for all sets
					for(Exon tarEx : m_mapIndistinguishableExonPerSample.keySet())
					{
						if(m_mapIndistinguishableExonPerSample.get(tarEx).containsKey(strSample)) // not all samples might have this set of ambigious exons
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
							
							// alright, we found the correct target exon
							if(bIsSame)
								mapTargetExonsPerSample.put(strSample, tarEx);
						}
					}
				}
			}
		}

		// Now we know which exons should be returned for the coverage value of the query exon
		double pCoverages[] = null;
		pCoverages = new double[nSamples];

		int nCurSample = 0;
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				if(m_vcSelectedSamples.contains(strSample))
				{
					if(!mapTargetExonsPerSample.containsKey(strSample))
						continue;
					
					Exon targetExon = mapTargetExonsPerSample.get(strSample);
					
					double fCov = m_mapExonCoveragePerSample.get(targetExon).get(strSample);
					pCoverages[nCurSample] = fCov;
					nCurSample++;
					
					if(mapExonCoverageToSample != null)
						mapExonCoverageToSample.put(strSample, fCov);
				}
			}
		}

		return pCoverages;
	}
	
	//#############################################################################################
	//    Hiding irrelevant isoforms is done in two steps:
	//
	//    First, isoforms containing exons that are not supported by junction reads are removed.
	//
	//    Second, the coverage for the exons of the remaining isoforms are calculate and
	//    isoforms with too low coverage per base or coverage length are discarded.
	//#############################################################################################
	public boolean HideIrrelevantIsoforms(boolean bSkipFirstAndLastExon, boolean bDebug) throws IOException
	{
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		int nConditions = mapSamplesToConditions.size();
		
		//##################################
		//    remove irrelevant isoforms
		//##################################
		CalculateRawExonCoverage();

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
				
			// get exons for valid isoform
			TreeMap<String, int[]> mapExonsToValidIsoforms = new TreeMap<String, int[]>(); 
			for(String strIsoform : m_vcValidIsoforms)
			{
				for(String strIncludedIsoform : m_mapExonsToIsoforms.keySet())
				{
					if(strIsoform.equals(strIncludedIsoform) || strIsoform.equals(strIncludedIsoform.split("\\.")[0]))
					{
						mapExonsToValidIsoforms.put(strIsoform, m_mapExonsToIsoforms.get(strIncludedIsoform));
						break;
					}
				}
			}
	
			// isoforms have been removed, thus recalculate the exon groups
			m_pExonGroups = RecalculateExonGroups(mapExonsToValidIsoforms);
		
			// calculate exon expression. Remove unexpressed exons
			if(!CalculateExonCoverage())
				return false;
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
		
		// get exons for valid isoform
		TreeMap<String, int[]> mapExonsToValidIsoforms = new TreeMap<String, int[]>(); 
		for(String strIsoform : m_vcValidIsoforms)
		{
			for(String strIncludedIsoform : m_mapExonsToIsoforms.keySet())
			{
				if(strIsoform.equals(strIncludedIsoform) || strIsoform.equals(strIncludedIsoform.split("\\.")[0]))
				{
					mapExonsToValidIsoforms.put(strIsoform, m_mapExonsToIsoforms.get(strIncludedIsoform));
					break;
				}
			}
		}

		// remove isoforms that include exons with too low coverage
		for(String strIsoform : mapExonsToValidIsoforms.keySet())
		{
			int[] pExonIDs = mapExonsToValidIsoforms.get(strIsoform);
			
			for(int nExID : pExonIDs)
			{				
				int nExonID = Math.abs(nExID - (m_pExons.length-1));
				Exon exon = m_pExons[nExonID];
				
				boolean bIsFirstExon = false;
				boolean bIsLastExon	 = false;
				
				if(nExID == pExonIDs[0])
					bIsFirstExon = true;

				if(nExID == pExonIDs[pExonIDs.length-1])
					bIsLastExon = true;
				
				if(bSkipFirstAndLastExon && (bIsFirstExon || bIsLastExon))
					continue;

				int nInvalidInConditions = 0;
				
				TreeMap<String, Double> mapExonCoverageToSample = new TreeMap<String, Double>();
				if(bIsFirstExon || bIsLastExon)
				{					
					ExonGroup[] grps = m_Gene.computeOverlappingExonGroups();
					for(ExonGroup grp : grps)
					{
						if(grp.groupContainsExon(exon.getExonID()))
						{
							for(Exon curExon : grp.getExons())
							{
								// Only merge first with other first exons and last with other last exons
								Iterator<String> it = m_Gene.getGeneProductNames();
								boolean bOkay = false;
								while(it.hasNext())
								{
									String strCurrentIsoform = it.next();
									if(!vcIrrelevantIsoforms.contains(strCurrentIsoform.split("\\.")[0]))
									{
										Exon pIsoformExons[] = m_Gene.getSortedExonsForGeneProduct(strCurrentIsoform);
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
									TreeMap<String, Double> mapExonCoverageToSampleTmp = new TreeMap<String, Double>();
									GetFinalCoverageForExon(curExon, mapExonCoverageToSampleTmp);
									
//										System.out.println("adding coverage");									
//										System.out.println(Arrays.toString(pCoverageCurrent));
									
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
				}
				else
				{
					GetFinalCoverageForExon(exon, mapExonCoverageToSample);
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
						if(exon.getCodingStart() == 36014863 && exon.getCodingStop() == 36014880)
							System.out.println(strCondition + ": " + fMedianCov + " < " + m_nMinCovPerBase + " " + Arrays.toString(pCoverage));
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

		return true;
	}

	public void RunAnalysis(String strFileGTF, String strFileProject, String strConditionType, String strFirstGene, boolean bSkipFirstAndLastExon, boolean bDebug, int nThreads) throws IOException
	{
		boolean bOkay 	= true;
		if(strFirstGene != null)
			bOkay = false;
		
		System.out.println("GTF file: " + strFileGTF);
		System.out.println("project: " + strFileProject);
		System.out.println("Condition type: " + strConditionType);
		System.out.println("skipping first and last exons: " + bSkipFirstAndLastExon);
		System.out.println("debug mode: " + bDebug);
		System.out.println("starting at: " + strFirstGene);
		
		m_nThreads = nThreads;
		
		System.out.println("max number of threads: " + m_nThreads);
		
		// init project
		try
		{
			m_projectModel.Init(strFileProject, -1, -1, -1.0, -1, true);
		}
		catch(Exception e)
		{
			System.out.println("failed to open project file: " + strFileProject);
			return;
		}
		
		for(String strSample : m_projectModel.GetSamples())
			m_vcSelectedSamples.add(strSample);
		
		String strFileName = m_projectModel.GetProjectName() + "_" + m_nMinJunctionReads +"_" + m_nMinCovPerBase + "_" + m_fMinCoveredBases + "_" + m_fVariableExonThreshold + "_" + bSkipFirstAndLastExon + "_analysis_output.txt";
		
		// open output file
		PrintWriter pWriter = null;
		if(!bOkay)
		{
			pWriter = new PrintWriter(new FileWriter(strFileName, true));
		}
		else
		{
			pWriter = new PrintWriter(new FileWriter(strFileName, false));
			
			String strTimeStamp = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(Calendar.getInstance().getTime());
			pWriter.println("## " + strTimeStamp);
			pWriter.println("## minimum number of junction reads: " 			+ m_nMinJunctionReads);
			pWriter.println("## minimum coverage per base per exon: " 			+ m_nMinCovPerBase);
			pWriter.println("## minimum fraction of covered bases per exon: " 	+ m_fMinCoveredBases);
			pWriter.println("## coverage ratio threshold: " 					+ m_fVariableExonThreshold + "\n");
			pWriter.println("#ID\ttype\tENSEMBL_ID\tgene_symbol\tgrp_position\tconditionA_ratio\tconditionB_ratio\talt_exon_cov_per_base_A\talt_exon_cov_per_base_B\talt_exon_cov_per_base_diff\tabs_ratio_changes\trel_ratio_changes\tp.values_ratio\tex_start\tex_end\tconditionA_PSI\tconditionB_PSI\tincl_jun_start\tincl_jun_end\texcl_jun_start\texcl_jun_end\teffect_PSI\tp.value_PSI\tuses_novel_junction");
		}
		
		m_strSelectedConditionType	= strConditionType;
		
		if(strFileGTF.toLowerCase().contains("gff") || strFileGTF.toLowerCase().contains("gtf"))
		{
			int nID = 0;
			
			// analyze all genes in the GTF file
			GTFParser parser = new GTFParser(strFileGTF);
			while(true)
			{				
				m_Gene				= null;
				m_pExons 			= null;
				m_pExonGroups		= null;
				
				m_mapBigWigCoverageToSamples = null;
				m_vcValidIsoforms.clear();
				m_mapExonsToIsoforms.clear();
				
				m_mapExonCoveragePerSample					= null;
				m_mapExonCoveragePerSampleRawFromBigWig		= null;		
				m_mapIndistinguishableExonPerSample 		= null;
				m_mapExonCoverageLengthPerSampleRawFromBigWig = null;
				
				GTFGene gtf_gene = parser.nextGene();
				if(gtf_gene == null)
					break;
				
				m_Gene = gtf_gene.createGene();

				if(bDebug || !bOkay)
				{					
					if(m_Gene.getGeneName() != null)
					{
						if(!m_Gene.getGeneName().equals(strFirstGene))				
							continue;
					}
					else if(!m_Gene.getGeneID().equals(strFirstGene))
					{
						continue;
					}
					
					bOkay = true;
				}
				
				if(!bOkay)
					continue;
				
				if(m_Gene.getGeneName() == null)
				{
					GeneIdentifier gid = GetGeneIdentifierForGene(m_Gene.getGeneID());
					if(gid != null)
					{
						m_Gene.setGeneName(gid.m_strApprovedGeneSymbol);
					}
					else
						m_Gene.setGeneName("?");
				}
				
				System.out.println("current gene: " + m_Gene.getGeneName() + " (" + m_Gene.getGeneID() + ")");
				
				// get all exons
				Exon pExonsTmp[] = m_Gene.getSortedExons();
				
				// reverse exon order
				m_pExons = null;
				m_pExons = new Exon[pExonsTmp.length];
				for(int i=0; i<pExonsTmp.length; i++)
				{
					m_pExons[i] = pExonsTmp[pExonsTmp.length-i-1];
				}
				
				// skip 1-exon-transcripts
				if(m_pExons.length == 1)
					continue;

				// get exon groups				
				m_pExonGroups = m_Gene.computeOverlappingExonGroupsNormalOrder();
			
				//#####################################
				//    get all isoforms for the gene
				//#####################################
				if(bDebug) System.out.println("getting valid isoforms");

				Iterator<String> it = m_Gene.getGeneProductNames();
				while(it.hasNext())
				{
					String strGeneProductID = it.next();
					int[] pExonIDs = m_Gene.getExonsForGeneProduct(strGeneProductID);
					
					m_mapExonsToIsoforms.put(strGeneProductID, pExonIDs);
					
					// add all isoforms on gene change
					m_vcValidIsoforms.add(strGeneProductID.split("\\.")[0]);
				}
/*
				// recalculate expression values
				if(!ProcessCountData())
					continue;
*/			
				if(bDebug) System.out.println("isoforms: " + m_vcValidIsoforms);

				// get bigwig coverage
				if(bDebug) System.out.println("getting bigwig coverage");
				GetBigWigCoverageForGene();
				
				// remove irrelevant isoforms
				if(bDebug) System.out.println("removing irrelevant isoforms");
				if(!HideIrrelevantIsoforms(bSkipFirstAndLastExon, bDebug))
				{
					System.out.println("failed to hide irrelevant isoforms");
					continue;
				}
				
				if(bDebug) System.out.println("valid isoforms after: " + m_vcValidIsoforms);

				// recalculate exon groups based on valid isoforms
				TreeMap<String, int[]> mapExonsToValidIsoforms = new TreeMap<String, int[]>();
				
				// get exons for valid isoform 
				for(String strIsoform : m_vcValidIsoforms)
				{
					for(String strIncludedIsoform : m_mapExonsToIsoforms.keySet())
					{
						if(strIsoform.equals(strIncludedIsoform) || strIsoform.equals(strIncludedIsoform.split("\\.")[0]))
						{
							// get exons
							int[] pExons = m_mapExonsToIsoforms.get(strIncludedIsoform);
							
							if(pExons.length > 2)
							{
								mapExonsToValidIsoforms.put(strIsoform, pExons);
							}
														
							break;
						}
					}
				}
				
				if(bDebug) System.out.println("exon to valid isoforms: " + mapExonsToValidIsoforms);				

				// isoforms have been removed, thus recalculate the exon groups
				if(bDebug) System.out.println("recalculating exon groups");
				m_pExonGroups = null;
				m_pExonGroups = RecalculateExonGroups(mapExonsToValidIsoforms);
				
				if(bDebug) System.out.println(mapExonsToValidIsoforms);
				
				// calculate variable exons
				if(bDebug) System.out.println("calculating variable exons");
				AlternativeSplicingHit res = CalculateVariableExons();
				
				if(bDebug)
				{
					if(res != null)
					{
						System.out.println("########### variable exon coverage ###########");
						System.out.println(res);
					}
				}
				
				// calculate and add number of PSI score results and the best p-value				
//				PSI_score_container psi_scores = CalculatePSIScores();
				TreeSet<SimpleSpliceScore> psi_scores = CalculateNovelJunctionPSIScores(bDebug);
				
				// calculate and add number of PSI score results and the best p-value				
//				PSI_score_container psi_scores_novel = CalculateNovelJunctionPSIScores();
				
				if(bDebug)
				{
					System.out.println("########### psi_scores ###########");
					System.out.println(psi_scores);
				}

				// add to output
				if(bDebug) System.out.println("writing output");			
				TreeSet<SimpleSpliceScore> vcUsedSpliceEvent = new TreeSet<SimpleSpliceScore>(); 
				
				if(res != null)
				{
					int nResults = res.GetNumHits();
					
					for(int k=0; k<nResults; k++)
					{
						AlternativeSplicingExon hit = res.GetResult(k);
						int nStart = hit.GetStart();
						int nEnd   = hit.GetEnd();
					
						boolean bOverlapWithSpliceScore = false;
						for(SimpleSpliceScore r : psi_scores)
						{
							// check for overlap
							if(r.m_Exon.getCodingStart() == nStart && r.m_Exon.getCodingStop() == nEnd &&
									( (r.m_strConditionA.equals(hit.GetCondition(true)) && r.m_strConditionB.equals(hit.GetCondition(false)) ||
									  (r.m_strConditionB.equals(hit.GetCondition(true)) && r.m_strConditionA.equals(hit.GetCondition(false)) ))))
							{
								String strOut = "ID_" + nID + "\tcombined";
								strOut += "\t" + res.GetDataForExonGroup(hit.GetID());
								strOut += String.format(Locale.ENGLISH, "\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%.2f\t%.3e\t%b%n", r.m_Exon.getCodingStart(), r.m_Exon.getCodingStop(), r.m_strConditionA, r.m_strConditionB, r.m_JunctionInclusion.m_nStart, r.m_JunctionInclusion.m_nEnd, r.m_JunctionExclusion.m_nStart, r.m_JunctionExclusion.m_nEnd, r.m_fInclusionChange, r.m_fPValue, r.m_bIsNovel);
								pWriter.print(strOut);
								bOverlapWithSpliceScore = true;
								vcUsedSpliceEvent.add(r);
								nID++;
							}
						}
						
						if(!bOverlapWithSpliceScore)
						{
							String strOut = "ID_" + nID + "\tratio_only";
							strOut += "\t" + res.GetDataForExonGroup(hit.GetID());
							strOut += String.format(Locale.ENGLISH, "\t?\t?\t?\t?\t?\t?\t?\t?\t?\t?\t?%n");
							pWriter.print(strOut);
							nID++;
						}
					}
				}
				
				for(SimpleSpliceScore r : psi_scores)
				{
					if(!vcUsedSpliceEvent.contains(r))
					{
						String strOut = "ID_" + nID + "\tsplit_read_only\t"; 
						strOut += String.format(Locale.ENGLISH, "%s\t%s\t?\t?\t?\t?\t?\t?\t?\t?\t?\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%.2f\t%.3e\t%b%n", m_Gene.getGeneID(), m_Gene.getGeneName(), r.m_Exon.getCodingStart(), r.m_Exon.getCodingStop(), r.m_strConditionA, r.m_strConditionB, r.m_JunctionInclusion.m_nStart, r.m_JunctionInclusion.m_nEnd, r.m_JunctionExclusion.m_nStart, r.m_JunctionExclusion.m_nEnd, r.m_fInclusionChange, r.m_fPValue, r.m_bIsNovel);
						pWriter.print(strOut);
						nID++;
					}
				}

				pWriter.flush();
				
				if(bDebug)
					break;
			}
		}
		else
		{
			// NOT SUPPORTED
		}
		
		String strTimeStamp = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(Calendar.getInstance().getTime());
		pWriter.println("finished: " + strTimeStamp);
		
		pWriter.close();
	}
	
	public void SetParameters(int nMinJunctionReads, int nMinCovPerBase, double fMinFracCoveredBases, double fMinRatio, int nThreads, String strFileEntrezGene, String strFileGTF)
	{
		m_nMinJunctionReads					= nMinJunctionReads;
		m_nMinCovPerBase					= nMinCovPerBase;
		m_fMinCoveredBases					= fMinFracCoveredBases;
		m_fVariableExonThreshold			= fMinRatio;
		m_nThreads							= nThreads;
		
		m_strFileGTF						= strFileGTF;
		m_nReferenceType 					= GTF_REFERENCE_FILE;
		
		m_vcGeneIdentifier = new Vector<GeneIdentifier>();
		
		try
		{
			ReadGeneIdentifier(strFileEntrezGene);
			for(GeneIdentifier gid : m_vcGeneIdentifier)
			{
				gid.m_bIsValid = true;
			}
		}
		catch (FileNotFoundException e)
		{
			System.out.println(e.getMessage());
		}
	}

	public void LoadGeneAnnotation(String strFile) throws Exception
	{
		if(strFile.toLowerCase().endsWith(".gtf") || strFile.toLowerCase().endsWith(".gff"))
		{
			m_nReferenceType = GTF_REFERENCE_FILE;
			System.out.println("GTF reference detected");

			// validate gene identifier
			String strFileIdx = m_strFileGTF + ".idx";
			File pFile = new File(strFileIdx);
			if(!pFile.exists())
			{
				// missing index file, thus create it
				RandomAccessGFFReader reader = new RandomAccessGFFReader(new File(m_strFileGTF), pFile);
				reader.close();
			}
			
			Scanner pIn = new Scanner(pFile);
			TreeSet<String> vcGeneIDs = new TreeSet<String>();
			while(pIn.hasNextLine())
			{
				String strLine = pIn.nextLine();
				vcGeneIDs.add(strLine.split("\\s+")[0]);
			}
			pIn.close();
			
			for(GeneIdentifier gid : m_vcGeneIdentifier)
			{
				if(vcGeneIDs.contains(gid.m_strEnsemblGeneID) || vcGeneIDs.contains(gid.m_strEntrezGeneID))
					gid.m_bIsValid = true;
				else
					gid.m_bIsValid = false;
			}
		}
		else if(strFile.toLowerCase().endsWith(".refflat"))
		{
			m_nReferenceType = REFFLAT_REFERENCE_FILE;
			System.out.println("refflat reference detected");
			
			// validate gene identifier
			Scanner pIn = new Scanner(new File(strFile));
			TreeSet<String> vcGeneIDs = new TreeSet<String>();
			while(pIn.hasNextLine())
			{
				String strLine = pIn.nextLine();
				vcGeneIDs.add(strLine.split("\\s+")[1]);
			}
			pIn.close();

			for(GeneIdentifier gid : m_vcGeneIdentifier)
			{
				if(vcGeneIDs.contains(gid.m_strRefGeneID))
					gid.m_bIsValid = true;
				else
					gid.m_bIsValid = false;
			}
		}
	}

	public Comparator<Row> getAscRatingComparator()
	{
		return new Comparator<Row>()
		{
			@Override
			public int compare(Row arg0, Row arg1)
			{
				int nRating1 = Integer.parseInt(arg0.getId().split("_")[2]);
				int nRating2 = Integer.parseInt(arg1.getId().split("_")[2]);
				
				if(nRating1 < nRating2)
					return -1;
				else if(nRating1 > nRating2)
					return 1;
				else
					return 0;
			}
		};
	}
	
	public Comparator<Row> getDescRatingComparator()
	{
		return new Comparator<Row>()
		{
			@Override
			public int compare(Row arg0, Row arg1)
			{
				int nRating1 = Integer.parseInt(arg0.getId().split("_")[2]);
				int nRating2 = Integer.parseInt(arg1.getId().split("_")[2]);
				
				if(nRating1 > nRating2)
					return -1;
				else if(nRating1 < nRating2)
					return 1;
				else
					return 0;
			}
		};
	}

	public int DrawOverlappingTranscripts(int nX, int nY, int nIntronLength, double fShrinkageFactor, int nVerticalSpaceBetweenIsoforms, Graphics2D graph, FontMetrics fontMetrics) throws FileNotFoundException
	{
		TreeSet<Gene> vcOverlappingTranscripts = GetOverlappingTranscripts();
		
		if(vcOverlappingTranscripts == null)
			return nY;
		
		// opposite strand transcripts
		TreeMap<ExonGroup, TreeSet<CountElement>> mapExonGroupToOverlappingRegionsOrange = new TreeMap<ExonGroup, TreeSet<CountElement>>();
		TreeMap<ExonGroup, TreeSet<CountElement>> mapExonGroupToOverlappingRegionsRed    = new TreeMap<ExonGroup, TreeSet<CountElement>>();
		
		// same strand transcripts
		TreeMap<ExonGroup, TreeSet<CountElement>> mapExonGroupToOverlappingRegionsBlue    	  = new TreeMap<ExonGroup, TreeSet<CountElement>>();
		TreeMap<ExonGroup, TreeSet<CountElement>> mapExonGroupToOverlappingRegionsDarkBlue    = new TreeMap<ExonGroup, TreeSet<CountElement>>();
		
		for(Gene g : vcOverlappingTranscripts)
		{
			TreeMap<ExonGroup, TreeSet<CountElement>> map1 = mapExonGroupToOverlappingRegionsOrange;
			TreeMap<ExonGroup, TreeSet<CountElement>> map2 = mapExonGroupToOverlappingRegionsRed;
			
			if(g.isPlusStrand() == m_Gene.isPlusStrand())
			{
				map1 = mapExonGroupToOverlappingRegionsBlue;
				map2 = mapExonGroupToOverlappingRegionsDarkBlue;
			}
			
			// see which exon groups are overlapping
			for(ExonGroup grp : m_pExonGroups)
			{
				// check if the gene overlapps the exon group (warning level orange)
				if(g.getStart() <= grp.getGenomicStopOfGroup() && g.getStop() >= grp.getGenomicStartOfGroup())
				{
					if(map1.containsKey(grp))
					{
						TreeSet<CountElement> vcWarningOrange = map1.get(grp);
						
						CountElement r = new CountElement();
						r.m_nStart = Math.max(g.getStart(), grp.getGenomicStartOfGroup());
						r.m_nEnd = Math.min(g.getStop(), grp.getGenomicStopOfGroup());
						vcWarningOrange.add(r);
					}
					else
					{
						TreeSet<CountElement> vcWarningOrange = new TreeSet<CountElement>();

						CountElement r = new CountElement();
						r.m_nStart = Math.max(g.getStart(), grp.getGenomicStartOfGroup());
						r.m_nEnd = Math.min(g.getStop(), grp.getGenomicStopOfGroup());
						vcWarningOrange.add(r);
						
						map1.put(grp, vcWarningOrange);
					}
					
					// check if an exon overlaps the exon group (warning level red)
					for(Exon ex : g.getSortedExons())
					{
						if(ex.getCodingStart() <= grp.getGenomicStopOfGroup() && ex.getCodingStop() >= grp.getGenomicStartOfGroup())
						{
							if(map2.containsKey(grp))
							{
								TreeSet<CountElement> vcWarningRed = map2.get(grp);
								
								CountElement r = new CountElement();
								r.m_nStart = Math.max(ex.getCodingStart(), grp.getGenomicStartOfGroup());
								r.m_nEnd = Math.min(ex.getCodingStop(), grp.getGenomicStopOfGroup());
								vcWarningRed.add(r);
							}
							else
							{
								TreeSet<CountElement> vcWarningRed = new TreeSet<CountElement>();

								CountElement r = new CountElement();
								r.m_nStart = Math.max(ex.getCodingStart(), grp.getGenomicStartOfGroup());
								r.m_nEnd = Math.min(ex.getCodingStop(), grp.getGenomicStopOfGroup());
								vcWarningRed.add(r);
								
								map2.put(grp, vcWarningRed);
							}
						}
					}
				}
			}
		}
		
		// merge overlapping regions
		for(int i=0; i<4; i++)
		{
			TreeMap<ExonGroup, TreeSet<CountElement>> map = null;
			
			switch(i)
			{
				case 0: map = mapExonGroupToOverlappingRegionsOrange; break;
				case 1: map = mapExonGroupToOverlappingRegionsRed; break;
				case 2: map = mapExonGroupToOverlappingRegionsBlue; break;
				case 3: map = mapExonGroupToOverlappingRegionsDarkBlue; break;
			}
			
			for(ExonGroup grp : map.keySet())
			{
				boolean bMerged = true;
				
				while(bMerged)
				{
					bMerged = false;
					
					TreeSet<CountElement> vcRegions = map.get(grp);
					
					outer_loop:
					for(CountElement r1: vcRegions)
					{
						for(CountElement r2: vcRegions)
						{
							if(!r1.equals(r2))
							{
								if(r1.m_nStart <= r2.m_nEnd && r1.m_nEnd >= r2.m_nStart)
								{
									r1.m_nStart = Math.min(r1.m_nStart, r2.m_nStart);
									r1.m_nEnd   = Math.max(r1.m_nEnd, r2.m_nEnd);
									
									bMerged = true;
									
									vcRegions.remove(r2);
									break outer_loop;
								}
							}
						}
					}
				}
			}
		}
		
		for(int i=0; i<4; i++)
		{
		
			TreeMap<ExonGroup, TreeSet<CountElement>> map = null;
			Color pClrs[] 		= null;
			float pFractions[] 	= null;
			
			switch(i)
			{
				case 0:
				{
					map = mapExonGroupToOverlappingRegionsOrange;
					
					//  add description
					graph.setColor(Color.BLACK);
					String strText = "AS transcipts";
					int nTextWidth = fontMetrics.stringWidth(strText);
					int nTextHeight = fontMetrics.getHeight();
					graph.drawString(strText, nX - nTextWidth - 5, nY + 10 + (int)(nTextHeight*0.25));
					
					// set colors for gradient
					pClrs = new Color[3];
					pClrs[0] = Color.ORANGE;
					pClrs[1] = Color.RED;
					pClrs[2] = Color.ORANGE;
					
					pFractions = new float[3];
					pFractions[0] = 0.0f;
					pFractions[1] = 0.5f;
					pFractions[2] = 1.0f;
					break;
				}
					
				case 1:
				{
					map = mapExonGroupToOverlappingRegionsRed;
					
					pClrs = new Color[3];				
					pClrs[0] = new Color(150, 0, 0);
					pClrs[1] = Color.RED;
					pClrs[2] = new Color(150, 0, 0);
					
					pFractions = new float[3];
					pFractions[0] = 0.2f;
					pFractions[1] = 0.5f;
					pFractions[2] = 0.8f;
					break;
				}
					
				case 2:
				{
					map = mapExonGroupToOverlappingRegionsBlue;
					
					// change height
					nY += nVerticalSpaceBetweenIsoforms;
					
					//  add description
					graph.setColor(Color.BLACK);
					String strText = "S transcipts";
					int nTextWidth = fontMetrics.stringWidth(strText);
					int nTextHeight = fontMetrics.getHeight();
					graph.drawString(strText, nX - nTextWidth - 5, nY + 10 + (int)(nTextHeight*0.25));
					
					// set colors
					pClrs = new Color[3];
					pClrs[0] = new Color(180, 180, 255);
					pClrs[1] = new Color(100, 100, 255);
					pClrs[2] = new Color(180, 180, 255);
					
					pFractions = new float[3];
					pFractions[0] = 0.0f;
					pFractions[1] = 0.5f;
					pFractions[2] = 1.0f;
					break;
				}
					
				case 3:
				{
					map = mapExonGroupToOverlappingRegionsDarkBlue;
					
					pClrs = new Color[3];
					pClrs[0] = new Color(0, 0, 150);
					pClrs[1] = new Color(0, 0, 230);
					pClrs[2] = new Color(0, 0, 150);
					
					pFractions = new float[3];
					pFractions[0] = 0.0f;
					pFractions[1] = 0.5f;
					pFractions[2] = 1.0f;
					break;
				}
			}
		
			// draw regions
			int nXTmp = nX;
			for(ExonGroup grp : m_pExonGroups)
			{
				int nGrpWidth = (int)(grp.getExonGroupLengthInBp() * fShrinkageFactor);
				
				if(!map.containsKey(grp))
				{
					nXTmp += nGrpWidth + nIntronLength;
					continue;
				}
				
				for(CountElement r : map.get(grp))
				{
					int nOffset = (int) ((r.m_nStart - grp.getGenomicStartOfGroup()) * fShrinkageFactor);
					int nWidth  = (int) ((r.m_nEnd - r.m_nStart +1) * fShrinkageFactor);
					
					LinearGradientPaint clrGradient = new LinearGradientPaint(nXTmp+nOffset, nY, nXTmp+nOffset, nY+20, pFractions, pClrs);
					
					// color region
					graph.setPaint(clrGradient);
					graph.fillRect(nXTmp+nOffset, nY, nWidth, 20);
					
					graph.setPaint(Color.BLACK);
					graph.drawRect(nXTmp+nOffset, nY, nWidth, 20);
					
					// add orientation
					if(nWidth >= 20) // only for exons with a minimum size
					{
						String strText = "";
						
						int nTextLength = (int)Math.floor((nWidth-12)/7);
						
						for(int k=0; k<nTextLength; k++)
						{
							if((m_Gene.isPlusStrand() && (i == 2 || i == 3)) || !m_Gene.isPlusStrand() && (i == 0 || i == 1))
								strText += ">";
							else
								strText += "<";
						}
						
						int nTextWidth = fontMetrics.stringWidth(strText);
						int nTextHeight = fontMetrics.getHeight();
						
						int nXOffset = (int)((nWidth - nTextWidth) * 0.5);
						
						if(i == 3) graph.setPaint(Color.WHITE);
						graph.drawString(strText, nXTmp + nOffset + nXOffset, (int)(nY+6+nTextHeight*0.5));
					}
				}
				
				nXTmp += nGrpWidth + nIntronLength;
			}
		}
		
		nY += nVerticalSpaceBetweenIsoforms;
		
		return nY;
	}
		
	public TreeSet<Gene> GetOverlappingTranscripts() throws FileNotFoundException
	{
		if(m_Gene == null)
		{
			System.out.println("GetOverlappingTranscripts() -> no gene selected");
			return null;
		}
		
		TreeSet<Gene> vcGenes = null;

		// check if reference file exists
		File pFile = new File(m_strFileGTF);
		if(!pFile.exists())
		{
			System.out.println("Failed to open file: " + m_strFileGTF);
			return null;
		}

		// check whether it is a GTF/GFF file or refFlat file
		if(m_nReferenceType == GTF_REFERENCE_FILE)
		{
			// check if the index file exists
			File pFileIndex = new File(m_strFileGTF + ".idx");

			try
			{
				RandomAccessGFFReader gffReader = new RandomAccessGFFReader(new File(m_strFileGTF), pFileIndex);
				vcGenes = gffReader.GetGenesForRange(m_Gene.getChromosome(), m_Gene.getStart()-2000, m_Gene.getStop()+2000);
			}
			catch(Exception e)
			{
				System.out.println("failed to open GTF file: " + m_strFileGTF);
				System.out.println(e.getMessage());
			}
		}
		//TODO
		else if(m_nReferenceType == REFFLAT_REFERENCE_FILE)
		{
			Scanner pIn = new Scanner(pFile);

			while(pIn.hasNextLine())
			{
				String strLine = pIn.nextLine();
				
				if(strLine.startsWith("#"))
					continue;
				
				/*
				String pSplit[] = strLine.split("\\s+");
				
				String strGeneID = pSplit[1];
				String strRef = pSplit[2];
				int nStart  = Integer.parseInt(pSplit[4])+1; // must be made 1-based
				int nEnd	= Integer.parseInt(pSplit[5]);  
				String strGeneSymbol = pSplit[12];
				
				if(strGeneID.equals(strGene) || (m_Gene != null && m_Gene.getGeneName().equals(strGeneSymbol) && m_Gene.getStart() <= nEnd && m_Gene.getStop() >= nStart && m_Gene.getChromosome().equals(strRef)))
				{
					boolean bFirstStrand = false;
					if(pSplit[3].equals("+"))
						bFirstStrand = true;
					
					Gene trans = new Gene(strGeneID, strGeneSymbol, "?", nStart, nEnd, strRef, bFirstStrand);
					trans.ParseFromRefFlat(pSplit);

					if(m_Gene == null)
						m_Gene = trans;
					
					// add new isoform to current gene
					m_Gene.addGeneProduct(strGeneID,  "", trans.getSortedExons());
				}
				*/
			}
			
			pIn.close();
			
			return vcGenes;
		}
		else
		{
			Messagebox.show("Invalid reference file detected");
			return null;
		}
		
		// remove current gene
		if(vcGenes.contains(m_Gene))
			vcGenes.remove(m_Gene);

		return vcGenes;
	}

/*
	public PSI_score_container CalculatePSIScores() throws IOException
	{
		// result container
		PSI_score_container results = new PSI_score_container();
		
		// get samples per Condition
		TreeMap<String, TreeSet<String>> mapSamplesToConditions =  m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		TreeSet<Integer> vcJunctionsStartPositions = new TreeSet<Integer>();
		TreeSet<Integer> vcJunctionsEndPositions   = new TreeSet<Integer>();
		
		for(int i=0; i<m_pExons.length; i++)
		{
			// first exon has junction start position only
			if(i==0)
			{
				vcJunctionsStartPositions.add(m_pExons[i].getGenomicStop());
			}
			// last exon has exon end position only
			else if(i == m_pExons.length-1)
			{
				vcJunctionsEndPositions.add(m_pExons[i].getGenomicStart());
			}
			else
			{
				vcJunctionsStartPositions.add(m_pExons[i].getGenomicStop());
				vcJunctionsEndPositions.add(m_pExons[i].getGenomicStart());
			}
		}
		
		// check all exon permutations and calculate all junction paths between these exons
		for(int nPathStart : vcJunctionsStartPositions)
		{
			for(int nPathEnd : vcJunctionsEndPositions)
			{
				Vector<TreeSet<CountElement>> vcJunctionPaths = new Vector<TreeSet<CountElement>>();

				// skip those events where the exon end is before the exon start
				if(nPathEnd < nPathStart)
					continue;
				
				Iterator<String> it = m_Gene.getGeneProductNames();
				
				// search for isoforms containing both exons
				while(it.hasNext())
				{
					String strIsoform = it.next();		
					
//					Exon[] pExons = m_Gene.getSortedExonsForGeneProduct(strIsoform);
					Exon[] pExons = m_Gene.getUnsortedExonsForGeneProduct(strIsoform);
					
					int nValid = 0;
					for(Exon ex : pExons)
					{
						if(ex.getCodingStop() == nPathStart || ex.getCodingStart() == nPathEnd)
							nValid += 1;
						
						if(nValid == 2)
							break;
					}

					// isoform must contain both exons
					if(nValid == 2)
					{
						// get all junctions for the current gene product
						TreeSet<Junction> vcJunctions = m_Gene.getSpliceJunctionInformationForGeneProductAsJunctions(strIsoform);
						
						// get junction path from exon 1 to exon 2
						TreeSet<CountElement> vcIncludedJunctions = new TreeSet<CountElement>();
						
						for(Junction jun : vcJunctions)
						{
							int nStart = jun.getCodingStart();
							int nEnd   = jun.getCodingStop();
							
							if(nStart >= nPathStart && nEnd <= nPathEnd)
							{
								CountElement e = new CountElement();
								e.m_nStart = nStart;
								e.m_nEnd = nEnd;
								vcIncludedJunctions.add(e);
							}
						}
						
						boolean bNewPath = true;
						for(int i=0; i<vcJunctionPaths.size(); i++)
						{
							// skip paths with unequal lengths
							if(vcJunctionPaths.get(i).size() != vcIncludedJunctions.size())
								continue;
							
							int nOverlap = 0;
							for(CountElement e : vcIncludedJunctions)
							{
								if(vcJunctionPaths.get(i).contains(e))
								{
									nOverlap++;
								}
							}
							
							if(nOverlap == vcIncludedJunctions.size())
								bNewPath = false;
						}
						
						if(bNewPath)
							vcJunctionPaths.add(vcIncludedJunctions);
					}
				}

				if(vcJunctionPaths.size() <= 1)
					continue;
			
				TreeMap<String, Integer> mapSamplesToIndex = new TreeMap<String, Integer>(); 
				
				// get junction counts
				TreeMap<String, TreeMap<String, Integer>> mapJunctionCounts = m_projectModel.GetJunctionCountsForGene(m_Gene);
				
				// add counts to junction paths
				for(TreeSet<CountElement> vcJunctions : vcJunctionPaths)
				{
					for(CountElement jun : vcJunctions)
					{
						for(String strJunction : mapJunctionCounts.keySet())
						{
							String pSplit[] = strJunction.split("_");
							
							int nStart = Integer.parseInt(pSplit[0]);
							int nEnd   = Integer.parseInt(pSplit[1]);
							
							if(nStart == jun.m_nStart && nEnd == jun.m_nEnd)
							{
								if(mapSamplesToIndex.isEmpty())
								{
									int nIdx = 0;
									for(String strSample : mapJunctionCounts.get(strJunction).keySet())
									{
										mapSamplesToIndex.put(strSample, nIdx);
										nIdx++;
									}
								}
								
								for(String strSample : mapJunctionCounts.get(strJunction).keySet())
								{
									int nValue = mapJunctionCounts.get(strJunction).get(strSample);
									jun.m_vcCounts.add(nValue);
								}
							}
						}
					}
				}
				
				if(mapSamplesToIndex.size() == 0)
				{
					// no counts and thus no samples for this junction
					continue;
				}
				
				int nPaths = vcJunctionPaths.size();
				
				// discard paths with insufficient coverage (<5 median number of reads in all conditions)
				boolean bModified = true;
				while(bModified)
				{
					bModified = false;
					
					loop_paths:
					for(int i=0; i<nPaths; i++)
					{
						for(CountElement e : vcJunctionPaths.get(i))
						{
							if(e.m_vcCounts.size() == 0)
								continue;
							
							int nInvalid = 0;
							for(String strCondition : mapSamplesToConditions.keySet())
							{
								int nSamples = mapSamplesToConditions.get(strCondition).size();
								double pValues[] = new double[nSamples];
	
								int nSampleIdx = 0;
								for(String strSample : mapSamplesToConditions.get(strCondition))
								{
									int nIdx = mapSamplesToIndex.get(strSample);
										
									pValues[nSampleIdx] = e.m_vcCounts.get(nIdx);
									nSampleIdx++;
								}
								
								double fMedian = StatUtils.percentile(pValues,  50.0);
								
								if(fMedian < 5)
									nInvalid++;
							}
							
							// remove path
							if(nInvalid == mapSamplesToConditions.keySet().size())
							{
								vcJunctionPaths.remove(i);
								nPaths--;
								bModified = true;
								break loop_paths;
							}
						}
					}
				}
				
				// calculate psi scores
				for(int i=0; i<nPaths; i++)
				{
					loop_current_path:
					for(int j=i+1; j<nPaths; j++)
					{
						TreeSet<CountElement> vcPathA = new TreeSet<CountElement>();
						TreeSet<CountElement> vcPathB = new TreeSet<CountElement>();
						TreeSet<CountElement> vcPathShared = new TreeSet<CountElement>();
						
						for(CountElement e : vcJunctionPaths.get(i))
							vcPathA.add(e);
						
						for(CountElement e : vcJunctionPaths.get(j))
							vcPathB.add(e);
						
						// remove junctions that are present in both paths
						bModified = true;
						while(bModified)
						{
							bModified = false;
							
							outer_loop:
							for(CountElement e1 : vcPathA)
							{
								for(CountElement e2 : vcPathB)
								{
									if(e1.m_nStart == e2.m_nStart && e1.m_nEnd == e2.m_nEnd)
									{
										vcPathShared.add(e1);
										vcPathA.remove(e1);
										vcPathB.remove(e2);
										bModified = true;
										break outer_loop;
									}
								}
							}
						}
						
						if(vcPathA.size() < 1 || vcPathB.size() < 1)
							continue loop_current_path;
						
						Vector<String> vcConditionsInOrder 		= new Vector<String>();
						Vector<double[]> vcScoresPerConditions 	= new Vector<double[]>();
						
						for(String strCondition : mapSamplesToConditions.keySet())
						{
							int nSamples = mapSamplesToConditions.get(strCondition).size();
							double pScores[] = new double[nSamples];
							
							int nCurrentSampleIdx = 0;
							for(String strSample : mapSamplesToConditions.get(strCondition))
							{
								int nIdx = mapSamplesToIndex.get(strSample);
								
								double nCountsA = Integer.MAX_VALUE;
								double nCountsB = Integer.MAX_VALUE;
								
								// get counts for path A
								for(CountElement e : vcPathA)
								{
									nCountsA = Math.min(nCountsA, e.m_vcCounts.get(nIdx));
								}
								
								// get counts for path B
								for(CountElement e : vcPathB)
								{
									nCountsB = Math.min(nCountsB, e.m_vcCounts.get(nIdx));
								}
								
								// calculate ratio
								double fScore = 0.5;
								
								if(nCountsA + nCountsB != 0)
									fScore = nCountsA / (nCountsA + nCountsB);
								
								pScores[nCurrentSampleIdx] = fScore;
								nCurrentSampleIdx++;
							}
							
							vcScoresPerConditions.add(pScores);
							vcConditionsInOrder.add(strCondition);
						}
						
						// use one-way anova to identify differences between groups
						OneWayAnova anova = new OneWayAnova();
						double fPValue = anova.anovaPValue(vcScoresPerConditions);

						if(fPValue < 0.05)
						{
							results.AddResult(vcPathA, vcPathB, vcPathShared, fPValue, m_Gene, vcConditionsInOrder, vcScoresPerConditions);
						}
					}
				}
			}
		}
		
		return results;
	}
*/
	public TreeSet<SimpleSpliceScore> CalculateNovelJunctionPSIScores(boolean bDebug) throws IOException
	{	
		// result container
		TreeSet<SimpleSpliceScore> results = new TreeSet<SimpleSpliceScore>();
		
		// get samples per Condition
		TreeMap<String, TreeSet<String>> mapSamplesToConditions =  m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		// get all splice junctions from file
		TreeMap<CountElement, TreeMap<String, Integer>> mapJunctionCounts = m_projectModel.GetJunctionsForRangeAsCountElements(m_Gene.getChromosome(), m_Gene.getStart(), m_Gene.getStop());
		
		// get size factors
		TreeMap<String, Double> mapSizeFactorsToSamples			= m_projectModel.GetSizeFactors();
		
		// get number of conditions
		int nConditions = mapSamplesToConditions.keySet().size();
		
		// remove junctions with insufficient coverage
		TreeSet<CountElement> vcInvalid = new TreeSet<CountElement>();
		
		for(CountElement jun : mapJunctionCounts.keySet())
		{
			int nInvalid = 0;
			
			TreeMap<String, double[]> mapCoverageToConditions = new TreeMap<String, double[]>();
			for(String strCondition : mapSamplesToConditions.keySet())
			{
				int nSamples = mapSamplesToConditions.get(strCondition).size();
				mapCoverageToConditions.put(strCondition, new double[nSamples]);
			}
			
			for(String strCondition : mapSamplesToConditions.keySet())
			{							
				int nSampleIdx = 0;
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					if(mapJunctionCounts.get(jun).containsKey(strSample))
					{
						mapCoverageToConditions.get(strCondition)[nSampleIdx] = mapJunctionCounts.get(jun).get(strSample);									
					}
					else
					{
						mapCoverageToConditions.get(strCondition)[nSampleIdx] = 0;
					}
					nSampleIdx++;
				}
				double fMean = StatUtils.mean(mapCoverageToConditions.get(strCondition));
				
				if(fMean < m_nMinJunctionReads)
					nInvalid++;
			}
			
			if(nInvalid == mapSamplesToConditions.keySet().size())
			{
				vcInvalid.add(jun);
			}
		}
		
		if(bDebug) System.out.println("discarding " + vcInvalid.size() + " junctions due to insufficient coverage");
		for(CountElement jun : vcInvalid)
		{
			mapJunctionCounts.remove(jun);
		}
		
		if(bDebug)
		{
			System.out.println("remaining junctions:");
			for(CountElement e : mapJunctionCounts.keySet())
			{
				System.out.println(e);
			}
		}
		
		// test each exon in the specified path
		for(Exon ex : m_pExons)
		{
			if(bDebug) System.out.println("checking exon: " + ex);
			
			//#####################################################################
			//    Get number of reads for all inclusion and exclusion junctions
			//#####################################################################
			TreeMap<CountElement, double[]> mapInclReads = new TreeMap<CountElement, double[]>();
			TreeMap<CountElement, double[]> mapExclReads = new TreeMap<CountElement, double[]>();
			
			for(CountElement jun : mapJunctionCounts.keySet())
			{
				// inclusion junction
				if(jun.m_nStart == ex.getCodingStop() || jun.m_nEnd == ex.getCodingStart())
				{
					double pInclReads[] = new double[nConditions];

					int nConditionIdx = 0;
					for(String strCondition : mapSamplesToConditions.keySet())
					{
						double fReads = 0.0;
						for(String strSample : mapSamplesToConditions.get(strCondition))
						{
							if(mapJunctionCounts.get(jun).containsKey(strSample))
							{
								fReads += mapJunctionCounts.get(jun).get(strSample) * mapSizeFactorsToSamples.get(strSample);
							}
						}
			
						pInclReads[nConditionIdx] = fReads;								
						nConditionIdx++;
					}
					
					mapInclReads.put(jun, pInclReads);
				}
				// exclusion junction
				else if(jun.m_nStart < ex.getCodingStart() && jun.m_nEnd > ex.getCodingStop())
				{
					double pExclReads[] = new double[nConditions];

					int nConditionIdx = 0;
					for(String strCondition : mapSamplesToConditions.keySet())
					{
						double fReads = 0.0;
						for(String strSample : mapSamplesToConditions.get(strCondition))
						{
							if(mapJunctionCounts.get(jun).containsKey(strSample))
							{
								fReads += mapJunctionCounts.get(jun).get(strSample) * mapSizeFactorsToSamples.get(strSample);
							}
						}
			
						pExclReads[nConditionIdx] = fReads;								
						nConditionIdx++;
					}
					
					mapExclReads.put(jun, pExclReads);
				}
			}
			
			if(bDebug)
			{
				System.out.println("inclusion junctions: " + mapInclReads);
				System.out.println("exclusion junctions: " + mapExclReads);
			}
			
			if(mapInclReads.size() < 1 || mapExclReads.size() < 1)
				continue;
			
			/*
			if(ex.getCodingStart() == 63353397 && ex.getCodingStop() == 63353472)
			{
				System.out.println("inclusion");
				for(Map.Entry<CountElement, double[]> e : mapInclReads.entrySet())
					System.out.println(e.getKey() + " " + Arrays.toString(e.getValue()));
				
				System.out.println("exclusion");
				for(Map.Entry<CountElement, double[]> e : mapExclReads.entrySet())
					System.out.println(e.getKey() + " " + Arrays.toString(e.getValue()));
			}
			*/
			
			//#######################################################################
			//    Find the junction supporting the in/exclusion of the exon that has
			//    the largest count difference between two given conditions.
			//#######################################################################
			String pConditions[] = new String[nConditions];
			mapSamplesToConditions.keySet().toArray(pConditions);
			
			for(int i=0; i<nConditions; i++)
			{
				for(int j=i+1; j<nConditions; j++)
				{
					String strConditionA = pConditions[i];
					String strConditionB = pConditions[j];
					
					// find inclusion/exclusion combination with the largest difference between conditions
					TreeSet<SimpleSpliceScore> vcBestSpliceScores = new TreeSet<SimpleSpliceScore>();
					
					// do this twice, once for known junctions only, and if no significant hit was identified, also consider novel junctions
					for(int n=0; n<2; n++)
					{
						for(CountElement inclJunction : mapInclReads.keySet())
						{
							if(n==0 && !inclJunction.m_bKnown)
								continue;
								
							for(CountElement exclJunction : mapExclReads.keySet())
							{
								if(n==0 && !exclJunction.m_bKnown)
									continue;
								
//								int nSamplesA = mapSamplesToConditions.get(strConditionA).size();
//								int nSamplesB = mapSamplesToConditions.get(strConditionB).size();
								
								// test the inclusion junction reads against the exclusion junction reads
								Vector<Double> vcScoresA = new Vector<Double>();
								Vector<Double> vcScoresB = new Vector<Double>();
								
								TreeMap<String, Integer> mapInclReadsToSamples = mapJunctionCounts.get(inclJunction);
								TreeMap<String, Integer> mapExclReadsToSamples = mapJunctionCounts.get(exclJunction);
								
//								int nSampleIdx = 0;
								for(String strSample : mapSamplesToConditions.get(strConditionA))
								{
									if(!mapInclReadsToSamples.containsKey(strSample) || !mapExclReadsToSamples.containsKey(strSample))
										continue;
									
									double fSum = mapInclReadsToSamples.get(strSample) + mapExclReadsToSamples.get(strSample);
									
									if(fSum == 0.0)
									{
										// skip it
									}
									else
									{
										vcScoresA.add(mapInclReadsToSamples.get(strSample) / fSum);
									}
//									nSampleIdx++;
								}
								
//								nSampleIdx = 0;
								for(String strSample : mapSamplesToConditions.get(strConditionB))
								{
									if(!mapInclReadsToSamples.containsKey(strSample) || !mapExclReadsToSamples.containsKey(strSample))
										continue;
									
									double fSum = mapInclReadsToSamples.get(strSample) + mapExclReadsToSamples.get(strSample);
				
									if(fSum == 0.0)
									{
										// skip it
									}
									else
									{
										vcScoresB.add(mapInclReadsToSamples.get(strSample) / fSum);
									}
								}
								
								// there must be at least two values for each group
								if(vcScoresA.size() > 1 && vcScoresB.size() > 1)
								{
									double pScoresA[] = new double[vcScoresA.size()];
									double pScoresB[] = new double[vcScoresB.size()];
									
									for(int nIdx=0; nIdx<vcScoresA.size(); nIdx++)
										pScoresA[nIdx] = vcScoresA.get(nIdx);
									
									for(int nIdx=0; nIdx<vcScoresB.size(); nIdx++)
										pScoresB[nIdx] = vcScoresB.get(nIdx);
									
									TTest test = new TTest();
									double fPValue = test.tTest(pScoresA, pScoresB);
									double fIncLevel = Math.abs(StatUtils.mean(pScoresA) - StatUtils.mean(pScoresB));
		/*								
										if(bDebug)
										{
											
			//								if(ex.getCodingStart() == 35223218 && ex.getCodingStop() <= 35223334)
											{
												System.out.println(inclJunction);
												System.out.println(exclJunction);
												System.out.println("inclusion reads: " + mapInclReadsToSamples);
												System.out.println("exclusion reads: " + mapExclReadsToSamples);
												System.out.println("p-value: " + fPValue + " " + Arrays.toString(pScoresA) + " " + Arrays.toString(pScoresB));
											}
										}
		*/							
									if(fPValue <= 0.05)
									{
										SimpleSpliceScore res = new SimpleSpliceScore(m_Gene.getGeneID(), inclJunction, exclJunction, ex, fPValue, fIncLevel, strConditionA, strConditionB);
										res.GetValidIsoforms(m_Gene);
										vcBestSpliceScores.add(res);
									}
								}
							}
						}
					}
					
					if(vcBestSpliceScores.size() > 0 )
					{
						// get best known splice junction hit and best novel splice junction hit
						double fPValueKnown = Double.MAX_VALUE;
						double fPValueNovel = Double.MAX_VALUE; 
						SimpleSpliceScore bestKnown = null;
						SimpleSpliceScore bestNovel = null;
						
						/*
						if(bDebug)
						{
							if(ex.getCodingStart() == 35223218 && ex.getCodingStop() <= 35223334)
							{
								System.out.println("BEST scores:");
								System.out.println(vcBestSpliceScores);
							}
						}
						*/
						
						for(SimpleSpliceScore score : vcBestSpliceScores)
						{
							if(!score.m_bIsNovel)
							{
								if(score.m_fPValue < fPValueKnown)
								{
									fPValueKnown = score.m_fPValue;
									bestKnown = score;
								}
							}
							else
							{
								if(score.m_fPValue < fPValueNovel)
								{
									fPValueNovel = score.m_fPValue;
									bestNovel = score;
								}
							}
						}
						
						if(bestKnown != null)
							results.add(bestKnown);
						
						if(bestNovel != null)
						{
							if(bestKnown == null || (!bestKnown.equals(bestNovel) && bestKnown.m_fPValue > bestNovel.m_fPValue))
								results.add(bestNovel);
						}
						
						/*
						if(bDebug)
						{
							if(ex.getCodingStart() == 35223218 && ex.getCodingStop() <= 35223334)
							{
								System.out.println("selected scores");
								System.out.println(bestKnown);
								System.out.println(bestNovel);
							}
						}
						*/
					}
				}
			}
		}
		return results;
	}
	
	Color getHeatMapColor(float value)
	{
		 // http://www.andrewnoske.com/wiki/Code_-_heatmaps_and_color_gradients
		int NUM_COLORS = 4;
		float color[][] = new float[NUM_COLORS][3]; // { {0,0,1}, {0,1,0}, {1,1,0}, {1,0,0} };
		color[0][0] = 0;
		color[0][1] = 0;
		color[0][2] = 1;
		
		color[1][0] = 0;
		color[1][1] = 1;
		color[1][2] = 0;
		
		color[2][0] = 1;
		color[2][1] = 1;
		color[2][2] = 0;
		
		color[3][0] = 1;
		color[3][1] = 0;
		color[3][2] = 0;
	    //A static array of 4 colors:  (blue,   green,  yellow,  red) using {r,g,b} for each.
	 
	    int idx1;        // |-- Our desired color will be between these two indexes in "color".
	    int idx2;        // |
	    float fractBetween = 0.0f;  // Fraction between "idx1" and "idx2" where our value is.
	 
	    if(value <= 0)
	    {
	    	idx1 = idx2 = 0;    // accounts for an input <=0
	    }
	    else if(value >= 1)
	    {
	    	idx1 = idx2 = NUM_COLORS-1;     // accounts for an input >=0
	    }
	    else
	    {
	    	value = value * (NUM_COLORS-1);        // Will multiply value by 3.
	    	idx1  = (int)Math.floor(value);                  // Our desired color will be after this index.
	    	idx2  = idx1+1;                        // ... and before this index (inclusive).
	    	fractBetween = value - (float)idx1;    // Distance between the two indexes (0-1).
	    }
	 
	    float red   = (color[idx2][0] - color[idx1][0])*fractBetween + color[idx1][0];
	    float green = (color[idx2][1] - color[idx1][1])*fractBetween + color[idx1][1];
	    float blue  = (color[idx2][2] - color[idx1][2])*fractBetween + color[idx1][2];
	    
	    Color clrRes = new Color(red, green, blue);
	    return clrRes;
	}
	
	public void CreateCoverageMatrix() throws IOException
	{
		PrintWriter pOut = new PrintWriter(new File("C:/Data/tmp/" + m_Gene.getGeneName() + "_all_exon_groups.txt"));
		ExonGroup pGroups[] = m_Gene.computeOverlappingExonGroups();

		TreeMap<String, TreeSet<String>> mapSamplesToConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		// write header
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				pOut.print("\t" + strSample);
			}			
		}
		pOut.print("\n");

		for(ExonGroup grp : pGroups)
		{
			pOut.print(m_Gene.getChromosome() + ":" + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup());
			
			for(String strCondition : mapSamplesToConditions.keySet())
			{
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					double pValues[] = GetCoverageForExonGroup(grp, strSample, true);			
					double fMeanCoverage = StatUtils.mean(pValues);

					pOut.print("\t" + fMeanCoverage);	
				}
			}
			pOut.print("\n");
		}
		
		pOut.close();
		
		// output matrix for selected exon groups
		pOut = new PrintWriter(new File("C:/Data/tmp/" + m_Gene.getGeneName() + "_selected_exon_groups.txt"));

		// write header
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				pOut.print("\t" + strSample);
			}			
		}
		pOut.print("\n");

		for(ExonGroup grp : m_pExonGroups)
		{
			pOut.print(m_Gene.getChromosome() + ":" + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup());
			
			for(String strCondition : mapSamplesToConditions.keySet())
			{
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					double pValues[] = GetCoverageForExonGroup(grp, strSample, true);			
					double fMeanCoverage = StatUtils.mean(pValues);

					pOut.print("\t" + fMeanCoverage);	
				}
			}
			pOut.print("\n");
		}
		
		pOut.close();
		
		// output matrix for all exons
		Exon pExons[] = m_Gene.getSortedExons();

		pOut = new PrintWriter(new File("C:/Data/tmp/" + m_Gene.getGeneName() + "_all_exons.txt"));

		// write header
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				pOut.print("\t" + strSample);
			}			
		}
		pOut.print("\n");

		for(Exon ex : pExons)
		{
			pOut.print(m_Gene.getChromosome() + ":" + ex.getCodingStart() + "-" + ex.getCodingStop());
			
			for(String strCondition : mapSamplesToConditions.keySet())
			{
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					double pValues[] = GetCoverageForExon(ex, strSample, true);			
					double fMeanCoverage = StatUtils.mean(pValues);

					pOut.print("\t" + fMeanCoverage);	
				}
			}
			pOut.print("\n");
		}
		
		pOut.close();
	}
	
	public TreeMap<ExonGroup, Double> CalculateVariableExonsOld() throws IOException
	{
		TreeMap<ExonGroup, Double> mapResult = new TreeMap<ExonGroup, Double>();
		
		TreeMap<String, String> mapBigWigFiles 	= m_projectModel.GetBigWigFilesForSamples();
		
		if(mapBigWigFiles == null)
		{
			Messagebox.show("ERROR: no valid bigwig or bam files detected");
			return null;
		}
		
		// get reference name
		String strRef = m_Gene.getChromosome();
		
		// get samples per condition
		final TreeMap<String, TreeSet<String>> vcSamplesAndConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		//########################################
		//      prepare coverage containers
		//########################################
		// key 1 = exon_name, key 2 = sample, value = coverage array
		TreeMap<String, TreeMap<String, double[]>> mapCoverageToSamplesAndExons = new TreeMap<String, TreeMap<String, double[]>>();
		for(ExonGroup grp : m_pExonGroups)
		{			
			int nStart	= grp.getGenomicStartOfGroup();
			int nEnd	= grp.getGenomicStopOfGroup();
			
			String strName = strRef + ":" + nStart + "-" + nEnd;
			
			mapCoverageToSamplesAndExons.put(strName, new TreeMap<String, double[]>());
		}
		
		//#####################################
		//    obtain coverage for all exons 
		//#####################################
		for(String strCondition : vcSamplesAndConditions.keySet())
		{
			for(String strSample : vcSamplesAndConditions.get(strCondition))
			{
				if(!m_vcSelectedSamples.contains(strSample))
					continue;
	
				for(ExonGroup grp : m_pExonGroups)
				{
					int nStart	= grp.getGenomicStartOfGroup();
					int nEnd	= grp.getGenomicStopOfGroup();
					String strGrpName = strRef + ":" + nStart + "-" + nEnd;
					
					double[] pCoverage = GetCoverageForExonGroup(grp, strSample, true);
					
					if(m_bCoveragePlotUseLog2)
					{
						for(int i=0; i<pCoverage.length; i++)
						{
//							pCoverage[i] = Math.Log2(pCoverage[i]+1);
							pCoverage[i] = Math.log(pCoverage[i]+1) / Math.log(2);
						}
					}
					
					mapCoverageToSamplesAndExons.get(strGrpName).put(strSample, pCoverage);
				}
			}
		}
		
		// prepare fraction map
		TreeMap<String, TreeMap<ExonGroup, double[]>> mapCovFractionsPerExonGroupAndCondition = new TreeMap<String, TreeMap<ExonGroup, double[]>>();
		
		// get correct group coverage data
		for(ExonGroup grp : m_pExonGroups)
		{
			String strName = strRef + ":" + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup();
			TreeMap<String, double[]> mapData = mapCoverageToSamplesAndExons.get(strName);
			
			// key = condition, value = coverage per position for each sample
			TreeMap<String, double[][]> mapCoveragePerSamplePerCondition = new TreeMap<String, double[][]>();
	
			for(String strCondition : vcSamplesAndConditions.keySet())
			{
				int nValidSamples = 0;
				for(String strSample : vcSamplesAndConditions.get(strCondition))
				{
					if(m_vcSelectedSamples.contains(strSample))
						nValidSamples += 1;
				}

				if(nValidSamples > 0)
				{
					double[][] pValues =  new double[grp.getExonGroupLengthInBp()+1][nValidSamples];
					
					int nCurrentSample = 0;
					for(String strSample : vcSamplesAndConditions.get(strCondition))
					{
						if(!m_vcSelectedSamples.contains(strSample))
							continue;
						
						double pData[] =  mapData.get(strSample);						
						for(int x=0; x<pData.length; x++)
						{
							double fValue = pData[x];
							pValues[x][nCurrentSample] = fValue;
						}
						
						nCurrentSample += 1;
					}
					
					if(nCurrentSample != 0)
						mapCoveragePerSamplePerCondition.put(strCondition, pValues);
				}
			}
			
			TreeMap<String, double[]> mapFractionsToCondition = new TreeMap<String, double[]>();
			for(int x=0; x<grp.getExonGroupLengthInBp(); x++)
			{
				// get median coverage for all conditions
				TreeMap<String, Double> mapMedianCoverageToCondition = new TreeMap<String, Double>();
				
				// get coverage sum
				double fTotalCoverage = 0;
	
				for(String strCondition : mapCoveragePerSamplePerCondition.keySet())
				{
					double pValues[][] = mapCoveragePerSamplePerCondition.get(strCondition);					
					double pCov[] = new double[pValues[0].length];
					
					for(int nSample = 0; nSample<pValues[0].length; nSample++)
					{
						pCov[nSample] = pValues[x][nSample]; 
					}
	
					double fMedian = StatUtils.percentile(pCov, 50.0);
					mapMedianCoverageToCondition.put(strCondition, fMedian);
					fTotalCoverage += fMedian;
				}
				
				if(fTotalCoverage == 0)
					continue;
	
				// calculate contribution of each condition to the coverage
				for(String strCondition : mapMedianCoverageToCondition.keySet())
				{
					double fFraction = mapMedianCoverageToCondition.get(strCondition) / fTotalCoverage;
					
					if(mapFractionsToCondition.containsKey(strCondition))
					{
						double[] pFractions = mapFractionsToCondition.get(strCondition);
						pFractions[x] = fFraction;
					}
					else
					{
						double[] pFractions = new double[grp.getExonGroupLengthInBp()];
						pFractions[x] = fFraction;
						mapFractionsToCondition.put(strCondition, pFractions);
					}
				}
			}
			
			for(String strCondition : mapFractionsToCondition.keySet())
			{
				if(mapCovFractionsPerExonGroupAndCondition.containsKey(strCondition))
				{
					TreeMap<ExonGroup, double[]> mapFractionsToExonGroup = mapCovFractionsPerExonGroupAndCondition.get(strCondition);
					mapFractionsToExonGroup.put(grp, mapFractionsToCondition.get(strCondition));
					
//					System.out.println(grp + " -> " +  Arrays.toString(mapFractionsToCondition.get(strCondition)));
				}
				else
				{
					TreeMap<ExonGroup, double[]> mapFractionsToExonGroup = new TreeMap<ExonGroup, double[]>();
					mapFractionsToExonGroup.put(grp, mapFractionsToCondition.get(strCondition));
					
					mapCovFractionsPerExonGroupAndCondition.put(strCondition, mapFractionsToExonGroup);
					
//					System.out.println(grp + " -> " +  Arrays.toString(mapFractionsToCondition.get(strCondition)));
				}
			}
		}
		
		// perform test for each exon group to identify alternatively spliced exons
		for(String strCondition : mapCovFractionsPerExonGroupAndCondition.keySet())
		{
			for(ExonGroup grp1 : mapCovFractionsPerExonGroupAndCondition.get(strCondition).keySet())
			{
				double[] pFractionsGroup1 = mapCovFractionsPerExonGroupAndCondition.get(strCondition).get(grp1);
				double fMedianGroup1 = StatUtils.percentile(pFractionsGroup1, 50.0);
				
//				System.out.println(grp1 + " 1-> " +  Arrays.toString(mapCovFractionsPerExonGroupAndCondition.get(strCondition).get(grp1)));

				// get average coverage ratio for all other exons
				Vector<Double> vcAllGroups = new Vector<Double>();
				for(ExonGroup grp2 : mapCovFractionsPerExonGroupAndCondition.get(strCondition).keySet())
				{					
					if(grp1.getGroupID() != grp2.getGroupID())
					{
//						System.out.println(grp2 + " 2-> " +  Arrays.toString(mapCovFractionsPerExonGroupAndCondition.get(strCondition).get(grp2)));
						
						double fFractionsGroup2[] = mapCovFractionsPerExonGroupAndCondition.get(strCondition).get(grp2);
						vcAllGroups.add(StatUtils.mean(fFractionsGroup2));
					}
				}
				double pAllGroups[] = new double[vcAllGroups.size()];
				for(int i=0; i<vcAllGroups.size(); i++)
					pAllGroups[i] = vcAllGroups.get(i);
				
//				double fMeanAll  = StatUtils.mean(pAllGroups);		
				double fMedianAll  = StatUtils.percentile(pAllGroups, 50);
				
				double fRes = Math.abs(fMedianGroup1 - fMedianAll);
				if(fRes > m_fVariableExonThreshold)
				{				
					if(mapResult.containsKey(grp1))
					{
						if(fRes > mapResult.get(grp1))
							mapResult.put(grp1, fRes);
					}
					else
					{
						mapResult.put(grp1, fRes);
					}
//					System.out.println(strCondition + " " + grp1.getGenomicStartOfGroup() + "-" + grp1.getGenomicStopOfGroup() + " is significantly different" );
				}
				
//				System.out.println("grp: " + fMeanGrp1);
//				System.out.println("all: " + fMeanAll);
//				System.out.println(strCondition + " " + grp1.getGenomicStartOfGroup() + "-" + grp1.getGenomicStopOfGroup() + " " + ((fMeanGrp1-fMeanAll)/fMeanAll));
			}
		}
		
		return mapResult;
	}
	
	public void ConvertHitList(String strFileProject) throws IOException
	{
		// init project
		try
		{
			if(!m_projectModel.Init(strFileProject, -1, -1, -1.0, -1, true))
				return;
		}
		catch(Exception e)
		{
			System.out.println("failed to process project file: " + strFileProject);
			System.out.println(ExceptionUtils.getStackTrace(e));
			return;
		}
		
		TreeSet<AlternativeSplicingHit> vcHits = m_projectModel.LoadHitListOld(m_strPathHitLists);
		m_projectModel.SaveHitListToFile(m_strPathHitLists, vcHits);
	}

	public void ReduceSampleSize()
	{
		int nMaxTotalSamples = 50;
		
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		int nTotalSamples = m_projectModel.GetSamples().size();
		
		m_vcSelectedSamples.clear();
		
		if(nTotalSamples > nMaxTotalSamples)
		{
			int nMaxSamplesPerCondition = nMaxTotalSamples / mapSamplesToConditions.keySet().size();

			for(String strCondition : mapSamplesToConditions.keySet())
			{
				TreeSet<String> vcSamples = mapSamplesToConditions.get(strCondition);
				if(vcSamples.size() < nMaxSamplesPerCondition)
				{
					for(String strSample : vcSamples)
						m_vcSelectedSamples.add(strSample);
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
						m_vcSelectedSamples.add(pSamples[idx]);
					}
				}
			}
		}
		else
		{
			for(String strSample : m_projectModel.GetSamples())
				m_vcSelectedSamples.add(strSample);
		}

		UpdateComboboxForSampleSelection(true);
	}

}
