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
import org.zkoss.zul.Radio;
import org.zkoss.zul.Radiogroup;
import org.zkoss.zul.Separator;
import org.zkoss.zul.South;
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
import org.zkoss.zul.West;
import org.zkoss.zul.Window;

import BioKit.Exon;
import BioKit.ExonGroup;
import BioKit.GTFGene;
import BioKit.GTFParser;
import BioKit.Gene;
import BioKit.RandomAccessGFFReader;
import BioKit.Utils;

public class SplicingWebApp extends Window
{		
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
	
	static final int		GTF_REFERENCE_FILE 		= 1;
	static final int		REFFLAT_REFERENCE_FILE 	= 2;
	
	static final int		AS_TYPE_EXON_SKIPPING			= 1;
	static final int		AS_TYPE_ALT_START_UNIQUE_JUN	= 2;		// alternative start exon with unique junction
	static final int		AS_TYPE_ALT_END_UNIQUE_JUN		= 3;		// alternative end exon with unique junction
	static final int		AS_TYPE_ALT_START_SHARED_JUN	= 4;		// alternative start exon that shares the junction with a 'middle' exon
	static final int		AS_TYPE_ALT_END_SHARED_JUN		= 5;		// alternative start exon that shares the junction with a 'middle' exon
	static final int		AS_TYPE_ALT_START_UNIQUE_JUN_DOUBLE	= 6;	// this type is only set if there are at least two start exons that have coverage ratio changes
	static final int		AS_TYPE_ALT_END_UNIQUE_JUN_DOUBLE	= 7;	// this type is only set if there are at least two end exons that have coverage ratio changes
	static final int		AS_TYPE_ALT_START_SHARED_JUN_DOUBLE	= 8;	// this type is only set if there are at least two start exons that have coverage ratio changes
	static final int		AS_TYPE_ALT_END_SHARED_JUN_DOUBLE	= 9;	// this type is only set if there are at least two end exons that have coverage ratio changes
	static final int		AS_TYPE_RETAINED_INTRON				= 10;
	static final int		AS_TYPE_ALT_5_PRIME_EXON_END		= 11;
	static final int		AS_TYPE_ALT_3_PRIME_EXON_END		= 12;
	
	private GeneIdentifierHandler 	m_geneIdentifierHandler;
	
	// this layout defines regions for buttons/options (north), the splice graph (west) and main plots (center)
	private Borderlayout 			m_layout;
	private North					m_layoutNorth;
	private West					m_layoutWest;
	private Center					m_layoutCenter;

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
	
	// handle to application and popup window
	SplicingWebApp 					m_hWindow;
	Window							m_windowPopup;
	
	// plot factories
	PlotFactory						m_plotFactory;
	PlotFactoryGTEX					m_plotFactoryGTEX;
	
	ResultListHandler				m_resultListHandler;
	
	// stores information when clicking an element in the isoform plot
	private ClickEvent 			m_ClickEvent;
	Parameters					m_parameters;
	
	// grants access to the SQL junction and exon count data
	private ProjectModel 		m_projectModel;
	
	// available options
	private boolean				m_bSkipFirstAndLastExon;
	private boolean				m_bDetectIntronRetentionEvents;
	
	private Textbox				m_txtboxURL;
	private Bandbox				m_bandboxProjects;
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
	private Checkbox			m_checkboxShowQuartiles;
	private Checkbox			m_checkboxColorExonsAndJunctionsByCoverage;
	private Checkbox			m_checkboxShowSecondCoveragePlot;
	private	Checkbox			m_checkboxHighlightUniqueExons;
	private Checkbox			m_checkboxCoverageGrid;

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
	
	private TreeMap<String, TreeMap<Integer, TreeSet<String>>> 	m_mapAnalysedJunctionPaths;
	
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
	
	// The relative coverage for each path and sample is stored in this map: <sample, <path_id, ratio>>
	private TreeMap<String, TreeMap<Integer, Double>> 			m_mapCoveragePerJunctionPathAndSample;
	
	// values important for the popup windows	
	private TreeSet<String> 				m_vcPreviouslySelectedSamples;	
	private TreeSet<String>					m_vcPreviouslySelectedIsoforms;
	private AnalysisResultHandler			m_vcASResults;
	private AnalysisResult					m_selectedResult;	
	
	private String 												m_strPathToExonData;
	private String												m_strFileDEXSeqIDs;
	private TreeMap<String, TreeMap<String, Vector<Double>>> 	m_mapCountsPerExonAndTissue;
	private TreeMap<String, Double> 							m_mapEntropyToExonicPart;
	
	private class Parameters
	{
		ProjectModel 	m_project;
		long			m_nSelectedIsoforms;
		int				m_nProjectID;
		GeneIdentifier 	m_gid;
		boolean			m_bShowScreenshot;	// show just a screenshot instead of the web interface (or create it if necessary)
		TreeSet<Integer> m_vcSelectedSamples;
		String			m_strSelectedSamples;
		
		public Parameters()
		{
			m_project 			= null;
			m_nSelectedIsoforms = 0;
			m_nProjectID		= 0;
			m_gid 				= null;
			m_bShowScreenshot	= false;
			m_vcSelectedSamples	= new TreeSet<Integer>();
			m_strSelectedSamples = null;
		}
	};
	
	public SplicingWebApp(boolean bEmpty) throws IOException
	{
		Configure4JLogger();
//		ReadAppPaths();
		
		m_hWindow 			= this;
		m_gffReader			= null;
		
		m_dataSupplier 		= new DataSupplier();

		m_vcValidIsoforms	= new TreeSet<String>();
		
		m_mapCoveragePerJunctionPathAndSample = new TreeMap<String, TreeMap<Integer, Double>>();
		
		m_ClickEvent		= new ClickEvent();
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
		
		m_layoutNorth		= new North();
		m_layoutWest		= new West();
		m_layoutCenter		= new Center();
				
		m_layoutMainTop	 	= new North();
		m_layoutMainSouth	= new South();	
		m_layoutMainPlots 	= new Borderlayout();	
		m_layoutMainBottom 	= new Vlayout();
		
		m_vcASResults 		= new AnalysisResultHandler();
		
		m_bWindowSizeChanged = false;
	}

	public SplicingWebApp() throws Exception
	{
		Configure4JLogger();
		
		m_hWindow			= this;
		m_plotFactory	 	= new PlotFactory(m_hWindow);
		m_plotFactoryGTEX 	= new PlotFactoryGTEX(m_hWindow);
		m_ClickEvent		= new ClickEvent();
		m_gffReader			= null;

		ReadAppPaths();	
		
		m_dataSupplier 		= new DataSupplier();

		m_imgMapIsoforms	= new Imagemap();
		m_imgMapCoverage	= new Imagemap();
		m_imgMapColors		= null;
		
		m_layoutNorth		= new North();
		m_layoutWest		= new West();
		m_layoutCenter		= new Center();
				
		m_layoutMainTop	 	= new North();
		m_layoutMainSouth	= new South();		
		
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

		m_mapAnalysedJunctionPaths	 			= new TreeMap<String, TreeMap<Integer, TreeSet<String>>>();
		m_mapCoveragePerJunctionPathAndSample 	= new TreeMap<String, TreeMap<Integer, Double>>();
		
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
					strFileScreenshot = m_strScreenshotPath + "/" + m_parameters.m_nProjectID + "_" + m_parameters.m_gid.m_strEnsemblGeneID;
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
					Messagebox.show("ERROR: coult not open project file: " + m_parameters.m_project.GetFullPathOfProjectFile());
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
				
				m_bandboxSelectedGene.setText(m_parameters.m_gid.m_strEnsemblGeneID);
				m_bandboxSelectedGene.setValue(m_parameters.m_gid.m_strEnsemblGeneID);
				try
				{
					OnGeneChange(false, false);
				}
				catch(Exception e)
				{
					Messagebox.show("ERROR: could not open gene: " + m_parameters.m_gid.m_strApprovedGeneSymbol + " (" + m_parameters.m_gid.m_strEnsemblGeneID + ")");
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
						Messagebox.show("ERROR: could not select isoforms");
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
		//           add gui elements (top)
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
		//     add main plot layouts (center=coverage track, bottom = isoforms)
		//#################################################################################################
		m_layoutCenter.setVflex("min");
		m_layoutCenter.setParent(m_layout);

		m_layoutMainPlots = new Borderlayout();
		m_layoutMainPlots.setVflex("min");
		m_layoutMainPlots.setParent(m_layoutCenter);		

		m_layoutMainTop.setTitle("Coverage Track");
		m_layoutMainTop.setHeight("280px");
		m_layoutMainTop.setAutoscroll(true);
		m_layoutMainTop.setId("CoveragePlot");
		m_layoutMainTop.setAttribute("org.zkoss.zul.nativebar", "true");
		m_layoutMainTop.setParent(m_layoutMainPlots);

		m_layoutMainSouth.setTitle("Isoform View");
		m_layoutMainSouth.setHeight("600px");
		m_layoutMainSouth.setParent(m_layoutMainPlots);
		
		m_layoutMainBottom 	= new Vlayout();
		m_layoutMainBottom.setId("IsoformPlot");
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
			Messagebox.show("Invalid project identifier: " + nID);
		}
				
		// check if the gene is valid
		strParameter = Executions.getCurrent().getParameter("gene");
		if(strParameter == null)
			return;
		
		GeneIdentifier gid = m_geneIdentifierHandler.GetGeneIdentifierForGene(strParameter, null);
		if(gid == null)
		{
			Messagebox.show("unknown gene identifier: " + strParameter);			
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
	
	public RandomAccessGFFReader GetGFFReader()
	{
		return m_gffReader;
	}
	
	public North GetCoveragePlotRegion()
	{
		return m_layoutMainTop;
	}
	
	public West GetPlotRegion()
	{
		return m_layoutWest;
	}
	
	public Borderlayout GetBorderLayout()
	{
		return m_layout;
	}
	
	public Vlayout GetIsoformPlotRegion()
	{
		return m_layoutMainBottom;
	}
	
	public Imagemap GetImageMapForCoveragePlot()
	{
		return m_imgMapCoverage;
	}
	
	public Imagemap GetImageMapForIsoforms()
	{
		return m_imgMapIsoforms;
	}
	
	public ClickEvent GetClickEvent()
	{
		return m_ClickEvent;
	}
	
	public boolean DetectIntronRetentionEvents()
	{
		return m_bDetectIntronRetentionEvents;
	}
	
	public boolean ProjectHasPairedData()
	{
		return m_projectModel.ProjectHasPairedData();
	}
	
	public TreeSet<String> GetValidIsoforms()
	{
		return m_vcValidIsoforms;
	}
	
	public int GetMaximumThreads()
	{
		return m_nThreads;
	}
	
	public String GetSelectedConditionType()
	{
		return m_strSelectedConditionType;
	}
	
	public String GetSelectedCondition()
	{
		return m_strSelectedCondition;
	}
	
	public String GetPathToMMSeqExecutable()
	{
		return m_strMMSeqExecutable;
	}
	
	public String GetPathToMMSeqResults()
	{
		return m_strPathMMSeqResults;
	}
	
	public String GetPathToTemporaryDirectory()
	{
		return m_strTmpFolder;
	}
	
	public String GetScreenshotDirectory()
	{
		return m_strScreenshotPath;
	}
	
	public TreeSet<String> GetSelectedSamples()
	{
		return m_vcSelectedSamples;
	}
	
	public int GetMinimumCoveragePerBase()
	{
		return m_nMinCovPerBase;
	}
	
	public int GetMinimumJunctionReads()
	{
		return m_nMinJunctionReads;
	}
	
	public double GetMinimumCoveredBases()
	{
		return m_fMinCoveredBases;
	}
	
	public double GetVariableExonThreshold()
	{
		return m_fVariableExonThreshold;
	}
	
	public AnalysisResult GetSelectedResult()
	{
		return m_selectedResult;
	}
	
	public String GetGTFFile()
	{
		return m_strFileGTF;
	}
	
	public int GetReferenceType()
	{
		return m_nReferenceType;
	}
	
	public ProjectModel GetProjectModel()
	{
		return m_projectModel;
	}
	
	public DataSupplier GetDataSupplier()
	{
		return m_dataSupplier;
	}
	
	public AnalysisResultHandler GetResultHandler()
	{
		return m_vcASResults;
	}
	
	public GeneIdentifierHandler GetGeneIdentifierHandler()
	{
		return m_geneIdentifierHandler;
	}

	public void SetNumberOfThreads(int nThreads)
	{
		m_nThreads = nThreads;
		System.out.println("set max. number of threads to " + nThreads);
	}
		
	public void SetConditionType(String strConditionType)
	{
		m_strSelectedConditionType = strConditionType;
	}
		
	public void ClearCurrentGeneData()
	{
		m_dataSupplier.Clear();
		m_vcValidIsoforms.clear();
	}

	public boolean IsEntropyEnabled()
	{
		return m_bShowEntropyData;
	}
	
	public TreeMap<String, Double> GetEntropyPerExonicPart()
	{
		return m_mapEntropyToExonicPart;
	}
	
	public TreeMap<String, Vector<Double>> GetCountsForGTEXExonicPart(String strExonicPartID)
	{
		return m_mapCountsPerExonAndTissue.get(strExonicPartID);
	}
	
	public TreeMap<String, TreeMap<Integer, TreeSet<String>>> GetAnalysedJunctionPaths()
	{
		return m_mapAnalysedJunctionPaths;
	}
	
	public TreeMap<String, TreeMap<Integer, Double>> GetCoveragePerJunctionPathAndSample()
	{
		return m_mapCoveragePerJunctionPathAndSample;
	}
	
	public void ShowGTEXDataForExon(String strExon)
	{
		m_plotFactoryGTEX.ShowGTEXDataForExon(strExon);
	}
	
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

	private static final long serialVersionUID = 1L;

	// just to avoid the warning message...
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
		
		/*
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
				
				if(strGeneID.equals(strGene) || (gene != null && gene.getGeneName().equals(strGeneSymbol) && gene.getStart() <= nEnd && gene.getStop() >= nStart && gene.getChromosome().equals(strRef)))
				{
					boolean bFirstStrand = false;
					if(pSplit[3].equals("+"))
						bFirstStrand = true;
					
					Gene trans = new Gene(strGeneID, strGeneSymbol, "?", nStart, nEnd, strRef, bFirstStrand);
					trans.ParseFromRefFlat(pSplit);

					if(gene == null)
						gene = trans;
					
					// add new isoform to current gene
					gene.addGeneProduct(strGeneID,  "", trans.getSortedExons());
				}
			}
			
			pIn.close();
		}
		else
		{
			Messagebox.show("Invalid reference file detected");
		}
		*/

		if(gene == null)
		{
			Messagebox.show("Could not find " + gid.m_strApprovedGeneSymbol + " (" + gid.m_strEnsemblGeneID  + ") in GTF file.");
			return null;
		}
		
		return gene;
	}
	
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
	
	// the gene parameter is optional and may be null
	public void ProcessGeneInformation(GeneIdentifier gid, String strFileGTF, Gene gene) throws Exception
	{
		if(gene == null)
			gene = ReadGeneFromFile(gid, strFileGTF);
		
		if(gene == null)
			return;
		
		InitValidIsoforms(gene);
		
//		DetermineValidIsoforms();
		m_bIsoformSelectionChanged = true;
	
		// get coverage from bigwig files and junction counts, also detects invalid samples and conditions
		PrepareDataSupplier(gene);
		
		OnIsoformChange(true);
	}
	
	// get coverage from bigwig files and junction counts, also detects invalid samples and conditions
	public void PrepareDataSupplier(Gene gene)
	{
		// prepare exon and exon group arrays
		m_dataSupplier = new DataSupplier(m_hWindow, gene);
		
		m_dataSupplier.RetrieveCoverageData(m_hWindow);
		m_dataSupplier.RetrieveJunctionCounts(m_hWindow);
		m_dataSupplier.IdentifyInvalidJunctions(m_hWindow, false);
	}

	// calculates the exon coverage just using their posititon and the coverage in the bigwig files

	//################################################################################################
	//    Calculates the expression of an exon relative to all other exons in the same exon group.
	//    Returns an map with the median relative coverage (value) per condition (key).
	//################################################################################################
	public double[] GetCoverageForRegion(String strRef, int nStart, int nEnd, String strSample, boolean bAdjustForSizeFactor) throws IOException
	{		
		// check whether bigwig and bam files are available for the samples
		TreeMap<String, Double> mapSizeFactors = m_projectModel.GetSizeFactors();
		
		int nLength = nEnd - nStart + 1;
		double pCoverage[] = new double[nLength];
		for(int i=0; i<pCoverage.length; i++)
			pCoverage[i] = 0.0;
		
		// get size factor
		double fSizeFactor = mapSizeFactors.get(strSample);

		double pBigWigCoverage[] = null;
		
		pBigWigCoverage = m_dataSupplier.GetGeneCoverageArrayForSample(strSample);
		int nArrayStart = m_dataSupplier.GetGeneStart();
		
		if(pBigWigCoverage == null || nStart < m_dataSupplier.GetGeneStart() || nEnd > m_dataSupplier.GetGeneEnd())
		{
			nArrayStart = nStart;
			
			if(m_dataSupplier.GetGene() == null)
			{
				Messagebox.show("GetBigWigCoverageForGene - No gene selected");
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
				pCoverage[nIdx] *= fSizeFactor;
		}

		return pCoverage;
	}
	
	public double[] GetCoverageForExonGroup(ExonGroup grp, String strSample, boolean bAdjustForSizeFactor) throws IOException
	{
		//###################################
		//    get coverage for exon group
		//###################################
		String strRef 	= m_dataSupplier.GetReferenceName();
		int nGrpStart	= grp.getGenomicStartOfGroup();
		int nGrpEnd		= grp.getGenomicStopOfGroup();
		
		return(GetCoverageForRegion(strRef, nGrpStart, nGrpEnd, strSample, bAdjustForSizeFactor));
	}
	
	public double[] GetCoverageForExon(Exon ex, String strSample, boolean bAdjustForSizeFactor) throws IOException
	{
		//###################################
		//    get coverage for single exon
		//###################################
		String strRef 	= m_dataSupplier.GetReferenceName();
		int nStart		= ex.getCodingStart();
		int nEnd		= ex.getCodingStop();
		
		return(GetCoverageForRegion(strRef, nStart, nEnd, strSample, bAdjustForSizeFactor));
	}

	//#################################################################################
	//    returns isoforms that were choosen by the user by selecting certain exons
	//#################################################################################
	public TreeMap<String, int[]> GetSelectedIsoforms()
	{
		// <isoform_name, exon_ids>
		TreeMap<String, int[]> mapValidIsoforms = new TreeMap<String, int[]>();
		
		// get first and second selected exon position
		ExonGroup grpA = null;
		ExonGroup grpB = null;

		if(m_dataSupplier.GetExonGroups() != null)
		{
			for(ExonGroup grp : m_dataSupplier.GetExonGroups())
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
			String[] pIsoforms = m_dataSupplier.GetIsoformNames();
			
			for(String strIsoform : pIsoforms)
			{		
				Exon pExons[] = m_dataSupplier.GetExonsForIsoform(strIsoform);
										
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
					int[] pExonIDs = new int[pExons.length];
					int nIdx = 0;
					for(Exon ex : pExons)
					{
						pExonIDs[nIdx] = ex.getExonID();
						nIdx++;
					}
					mapValidIsoforms.put(strIsoform, pExonIDs);
				}
			}
		}
		return mapValidIsoforms;
	}

	public void OnProjectChange(String strProjectFile, boolean bHideMessage) throws Exception
	{
		m_strSelectedCondition 			= null;
		m_strSelectedConditionType		= null;
		
		m_projectModel.clear();
		m_projectModel.Init(strProjectFile, -1, -1, -1.0, -1, true);
		
		if(m_projectModel.GetSamples().size() > 50 && !bHideMessage)
		{
			Messagebox.show("The project has more than 50 samples. To speed up things, detection of retained introns was disabled. You can enable retained intron detection in the advanced options tab again.");
			m_bDetectIntronRetentionEvents = false;
			m_checkboxSkipIntronRetentionDetection.setChecked(m_bDetectIntronRetentionEvents);
		}

		m_bandboxSelectedGene.setText("Select Gene");
		
		m_vcASResults.Clear();
		m_dataSupplier.Clear();		
		
		UpdateComboboxesForCondition();
		
		// select first available condition type
		if(m_comboboxSelectedConditionType.getItemCount() > 0)
		{
			m_comboboxSelectedConditionType.setSelectedIndex(0);
			m_strSelectedConditionType = m_comboboxSelectedConditionType.getSelectedItem().getLabel();
			OnConditionTypeChange();
			
			// select first available condition
			m_comboboxSelectedCondition.setSelectedIndex(0);
			m_strSelectedCondition = m_comboboxSelectedCondition.getSelectedItem().getLabel();
		}
		
		//###########################################################################
		//    load previously calculated alternative splicing hits, if available
		//###########################################################################
		String strResultFile = m_projectModel.GetProjectName() + "_splicing_results.dat";
		File pIn = new File(m_strPathHitLists + "/" + strResultFile);
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
		
		m_mapCoveragePerJunctionPathAndSample 	= new TreeMap<String, TreeMap<Integer, Double>>();
		
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
		
		m_ClickEvent.clear();

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

		layoutV = new Vlayout();
		layoutV.setStyle("margin-left: 40px");
		layoutV.setParent(layoutH);
		
		// add label for project selection
		lab = new Label("Selected Project:");
		lab.setParent(layoutV);
		
		// add bandbox for project selection
		m_bandboxProjects = new Bandbox("Select Project");
		m_bandboxProjects.setParent(layoutV);
		m_bandboxProjects.setWidth("180px");
		
		// add bandpopup
		Bandpopup bandpopup = new Bandpopup();
		bandpopup.setParent(m_bandboxProjects);
		bandpopup.setWidth("75%");
		bandpopup.setHeight("75%");
		
		// create tree set of project information
		File pFolder = new File(m_strPathInput);
		for(File pFile : pFolder.listFiles())
		{
			String strFile = pFile.getName();
			
			if(strFile.endsWith("project"))// || strFile.endsWith("project2"))
			{
				ProjectModel project = new ProjectModel();
				project.Init(m_strPathInput + "/" + strFile, 0, 0, 0.0, 0, false);
				m_vcProjectInfos.add(project);
			}
		}
		
		Listbox listboxProjects = new Listbox();
		listboxProjects.setParent(bandpopup);		
		listboxProjects.setHflex("min");
		listboxProjects.setStyle("autoWidth:true");
		listboxProjects.setMold("paging");
		listboxProjects.setPageSize(10);
		
		Listhead listheadProject = new Listhead();
		listheadProject.setSizable(true);
		listheadProject.setParent(listboxProjects);
		
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
		
		listboxProjects.addEventListener(Events.ON_SELECT, new EventListener<Event>()
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
		
		// now add the sorted list of gene identifiers to the selection
		for(ProjectModel project : m_vcProjectInfos)
		{
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
			
			listboxProjects.appendChild(item);
		}
		
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
					TreeSet<GeneIdentifier> vcGeneIdentifier = new TreeSet<GeneIdentifier>();

					// get list of valid gene identifiers
					for(GeneIdentifier gid : m_geneIdentifierHandler.GetAllGeneIdentifiersStartingWith(strInput))
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

	public void AddOptionsMenu(Vlayout parentLayout)
	{
		Groupbox grpBox = new Groupbox();
		grpBox.setTitle("Options");
		grpBox.setMold("3d");
		grpBox.setParent(parentLayout);
		grpBox.setWidth("500px");
		grpBox.setStyle("margin-left: 10px;");
		grpBox.setVflex("true");
		
		Hlayout layoutH = new Hlayout();
		layoutH.setParent(grpBox);
		
		Vlayout layoutV = new Vlayout();
		layoutV.setHflex("true");
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
		
		m_checkboxUseLog2 = new Checkbox("Show log2 tranformed");
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
		
		Hbox layoutH2 = new Hbox();
		layoutH2.setParent(layoutV);
		layoutH2.setStyle("width: 100%; margin-top: 40px;");
		
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
		layoutV.setHeight("252px");
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
					Messagebox.show("No gene specified");
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
	
	public void AddAdvancedOptionsMenu(Vlayout parentLayout)
	{
		Groupbox grpBox = new Groupbox();
		grpBox.setTitle("Advanced Options");
		grpBox.setMold("3d");
		grpBox.setParent(parentLayout);
		grpBox.setWidth("440px");
		grpBox.setHeight("470px");

		Hlayout layoutH = new Hlayout();
		layoutH.setParent(grpBox);
		
		Vlayout layoutV = new Vlayout();
		layoutV.setWidth("220px");
		layoutV.setParent(layoutH);
		
		//#############################################################
		layoutH = new Hlayout();
		layoutH.setStyle("margin-top: 20px");
		layoutH.setParent(grpBox);
		
		layoutV = new Vlayout();
		layoutV.setWidth("220px");
		layoutV.setParent(layoutH);

		m_checkboxShowSecondCoveragePlot = new Checkbox("Attach coverage plot to isoforms");
		m_checkboxShowSecondCoveragePlot.setParent(layoutV);
		m_checkboxShowSecondCoveragePlot.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				m_plotFactory.ShowSecondCoveragePlot(m_checkboxShowSecondCoveragePlot.isChecked());
				m_plotFactory.RequestCoverageRedraw();
				m_plotFactory.DrawPlots();
			}
		});
		
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
		btnGetURL.setStyle("margin-top: 40px; margin-bottom: 10px;");
		
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
					strURL = strServer + ":" + nServerPort + strContentsPath + strRequestPath + "?ID=" + nProjectID + "&gene=" + strSelectedGene + "&isoforms=" + nSelectedIsoforms + "&samples=" + strSelectedSamples;
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
		
		if(m_selectedResult.GetGeneSymbol().isEmpty())
			m_bandboxSelectedGene.setValue(m_selectedResult.GetGeneID().split("\\.")[0]);
		else
			m_bandboxSelectedGene.setValue(m_selectedResult.GetGeneSymbol());
		
		if(!m_strFileGTF.equals(m_selectedResult.GetGTFFile()))
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
		
		boolean bSameGene = m_selectedResult.GetGeneID().equals(strOldGeneID);
		if(!bSameGene)
		{
			OnGeneChange(false, false);
			m_selectedResult = res;
		}
		
		m_vcValidIsoforms = m_selectedResult.GetValidIsoforms(m_dataSupplier);

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

		// check whether the selected gene is a new gene
		OnIsoformChange(!bSameGene);
		
		// check if any of the junctions is not annotated
		if(m_selectedResult.HasPSIScore())
		{
			if(m_dataSupplier.IsNovelJunction(m_selectedResult.GetInclusionJunction()))
			{
				m_plotFactory.AddUnknownJunction(m_selectedResult.GetInclusionJunction());
			}
			
			if(m_dataSupplier.IsNovelJunction(m_selectedResult.GetExclusionJunction()))
			{
				m_plotFactory.AddUnknownJunction(m_selectedResult.GetExclusionJunction());
			}
		}
	}
	
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
			m_treeSelectedSamples.setHeight("470px");
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
						
						if(m_dataSupplier.GetGene() != null)
							m_plotFactory.DrawPlots();
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

	//########################################################################
	//    Calculates the isoform expression value per sample
	//    Returns a map of the form <isoform_id, <sample, coverage_value>>
	//########################################################################
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
//		System.out.println(mapIsoformsToPositions);
		
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
	public void DetermineValidIsoforms() throws IOException, ClassNotFoundException, InstantiationException, IllegalAccessException
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
				
				int nClickedEdgeStart = Integer.parseInt(pSplitA[1])-1;
				int nClickedEdgeEnd   = Integer.parseInt(pSplitA[2])-1;
				
				ExonGroup pExonGroups[] = m_dataSupplier.GetExonGroups();
				String vcIsoforms[] 	= m_dataSupplier.GetIsoformNames();
				
				for(String strIsoform : vcIsoforms)
				{						
					Exon pExons[] = m_dataSupplier.GetExonsForIsoform(strIsoform);
					
					Vector<Integer> vcExonGroupIDs = new Vector<Integer>();

					// get exon group for each exon
					for(Exon ex : pExons)
					{
						for(ExonGroup grp : pExonGroups)
						{
							if(grp.groupContainsExon(ex.getExonID()))
							{
								vcExonGroupIDs.add(grp.getGroupID());
								break;
							}
						}
					}

					for(int i=0; i<vcExonGroupIDs.size(); i++)
					{
						if(vcExonGroupIDs.get(i) == nClickedEdgeStart && vcExonGroupIDs.get(i+1) == nClickedEdgeEnd)
						{
							m_vcValidIsoforms.add(strIsoform);
							break;
						}
					}

				}
			}
		}
		else
		{
			//##############################################
			//       select all isoforms for the gene
			//##############################################
			m_vcValidIsoforms.clear();
			for(Treeitem item : m_treeSelectedIsoforms.getSelectedItems())
			{
				m_vcValidIsoforms.add(item.getValue().toString());
			}
		}
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

	// uses Benjamini-Hochberg multiple-testing correction to adjust p-values of PSI scores
	protected TreeSet<SimpleSpliceScore> AdjustPValuesPsiScores(TreeSet<SimpleSpliceScore> vcScores, Vector<Double> vcPValues)
	{
		int nTests = vcPValues.size();
		
		// sort p-values asdencing
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

	public void LoadGeneAnnotation(String strFile)
	{
		if(strFile.toLowerCase().endsWith(".gtf") || strFile.toLowerCase().endsWith(".gff"))
		{
			m_strFileGTF = strFile;
			String strFileIdx = m_strFileGTF + ".idx";
			File pFile = new File(strFileIdx);
			
			m_nReferenceType = GTF_REFERENCE_FILE;
			System.out.println("GTF reference detected");

			if(m_gffReader == null || !m_gffReader.GetFileName().equals(m_strFileGTF))
			{
				// validate gene identifier
				try
				{
					m_gffReader = new RandomAccessGFFReader(new File(m_strFileGTF), pFile);
				}
				catch(Exception e)
				{
					System.out.println("failed to open GFF file: " + m_strFileGTF);
					e.printStackTrace();
					return;
				}
			}

			Scanner pIn = null;
			
			try
			{
				pIn = new Scanner(pFile);
			}
			catch(Exception e)
			{
				System.out.println("failed to open file: " + strFileIdx);
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
			
			for(GeneIdentifier gid : m_geneIdentifierHandler.GetAllGeneIdentifiers())
			{
				if(vcGeneIDs.contains(gid.m_strEnsemblGeneID) || vcGeneIDs.contains(gid.m_strEntrezGeneID))
					gid.m_bIsValid = true;
				else
					gid.m_bIsValid = false;
			}
		}
		else if(strFile.toLowerCase().endsWith(".refflat"))
		{
			/*
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

			for(GeneIdentifier gid : m_geneIdentifierHandler.GetAllGeneIdentifiers())
			{
				if(vcGeneIDs.contains(gid.m_strRefGeneID))
					gid.m_bIsValid = true;
				else
					gid.m_bIsValid = false;
			}
			*/
		}
	}

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
			e.printStackTrace();
			System.out.println("failed to topen file: " + m_strFileDEXSeqIDs);
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
		File pFolder = new File(m_strPathToExonData);
		if(!pFolder.exists())
		{
			Messagebox.show("ERROR: Invalid GTEX folder specified: " + m_strPathToExonData);
			return null;
		}
		
		for(String strPos : mapIDsToPositions.keySet())
		{
			String strFile = m_strPathToExonData + "//GTEX_" + strGTEXID + "_" + mapIDsToPositions.get(strPos) + ".tsv";
			File pFile = new File(strFile);
			
			if(!pFile.exists())
			{
				Messagebox.show("ERROR: No GTEX data available for exonic part: " + strGTEXID + " " + strPos + " " + mapIDsToPositions.get(strPos) + " -> " + strFile);
				continue;
			}
			
			try
			{
				pScanner = new Scanner(pFile);
			}
			catch(Exception e)
			{
				System.out.println("failed to open file: " + strFile);
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
		String strRef = m_dataSupplier.GetReferenceName();
		
		// get samples per condition
		final TreeMap<String, TreeSet<String>> vcSamplesAndConditions = m_projectModel.GetSamplesPerCondition(m_strSelectedConditionType);
		
		//########################################
		//      prepare coverage containers
		//########################################
		// key 1 = exon_name, key 2 = sample, value = coverage array
		TreeMap<String, TreeMap<String, double[]>> mapCoverageToSamplesAndExons = new TreeMap<String, TreeMap<String, double[]>>();
		for(ExonGroup grp : m_dataSupplier.GetExonGroups())
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
	
				for(ExonGroup grp : m_dataSupplier.GetExonGroups())
				{
					int nStart	= grp.getGenomicStartOfGroup();
					int nEnd	= grp.getGenomicStopOfGroup();
					String strGrpName = strRef + ":" + nStart + "-" + nEnd;
					
					double[] pCoverage = GetCoverageForExonGroup(grp, strSample, true);
					
					if(m_plotFactory.IsLog2Enabled())
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
		for(ExonGroup grp : m_dataSupplier.GetExonGroups())
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

	public void InitProjectModel(String strFileProject)
	{
		try
		{
			m_projectModel.Init(strFileProject, -1, -1, -1.0, -1, true);
		}
		catch(Exception e)
		{
			System.out.println("failed to open project file: " + strFileProject);
			return;
		}
	}

	public void SaveSplicingResults()
	{
		String strResultFile = m_projectModel.GetProjectName() + "_splicing_results.dat";
		m_vcASResults.SaveToFile(m_strPathHitLists, strResultFile);
	}
	
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

}
