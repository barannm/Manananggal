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
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.GradientPaint;
import java.awt.Graphics2D;
import java.awt.LinearGradientPaint;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collections;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import javax.imageio.ImageIO;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.OneWayAnova;
import org.zkoss.zk.ui.Executions;
import org.zkoss.zk.ui.event.Event;
import org.zkoss.zk.ui.event.EventListener;
import org.zkoss.zk.ui.event.Events;
import org.zkoss.zk.ui.event.MouseEvent;
import org.zkoss.zul.Area;
import org.zkoss.zul.Column;
import org.zkoss.zul.Columns;
import org.zkoss.zul.Grid;
import org.zkoss.zul.Hlayout;
import org.zkoss.zul.Imagemap;
import org.zkoss.zul.Label;
import org.zkoss.zul.Messagebox;
import org.zkoss.zul.Popup;
import org.zkoss.zul.Row;
import org.zkoss.zul.Rows;
import org.zkoss.zul.Vlayout;
import org.zkoss.zul.Window;

import BioKit.Exon;
import BioKit.ExonGroup;
import BioKit.Gene;
import BioKit.RandomAccessGFFReader;

/**
 *    The plot factory handles most of the drawing concerning coverage plots and isoform plots.
 *    The only exception are the plots associated with GTEX data that have their own PlotFactory.
 */
public class PlotFactory
{
	SplicingWebApp					m_app;                          // app object, used to get data on the current gene
	private int 					m_nMaxWidth;                    // maximum plot width (determined by client width by default)
	private TreeMap<String, Color> 	m_mapColorsToConditions;        // defines the colors per condition
	private int						m_nClientWindowWidth;           // store client window width
	private int						m_nClientWindowHeight;          // store client window height
	
	private boolean					m_bCoveragePlotUseLog2;         // show log2 coverage
	private boolean					m_bCoveragePlotShowMean;        //  show mean coverage 
	private boolean					m_bCoveragePlotShowGeometricMean;   // show geometric mean coverage
	private boolean					m_bCoveragePlotShowMedian;      // show median coverage
	private boolean					m_bCoveragePlotShowRelativeCoverage; // show relative coverage
	private boolean					m_bCoveragePlotShowQuartiles;   // draw the upper and lower quartiles in the coverage plot
	private boolean					m_bCoverageRequiresRedraw;      // indicates that the coverage plot must be redrawn
	private boolean					m_bShowUniqueFeatures;          // colors exons and junctions that are unique to a single isoform
	private boolean					m_bShowAmbigiousUniqueExons;    // color exons for which coverage could not be determine unambiguously
	private boolean					m_bShowCoverageGrid;            // show the light grey grid in the coverage window
	private boolean					m_bSwitchStrands; 	            // flip the image (second strand transcripts will be drawn from left to right)

	private boolean					m_bColorExonsAndJunctionsByCoverage; // color exons and junctions based on their coverage values

	TreeSet<CountElement>			m_vcUnknownJunctions;	// junctions that do not match any known isoform

	/** Helper class that stores drawing offsets */
	private class DrawingOffsets
	{
		int m_nBorderOffset;
		int m_nIsoformNameOffset;
		int m_nMaxExonLength;
		int m_nMMSeqTextLength;
		int m_nIntronLength;
		int m_nTotalWidth;
		
		public String toString()
		{
			String strOut = "border offset: " + m_nBorderOffset + "\n" +
					        "isoform offset: " + m_nIsoformNameOffset + "\n" +
					        "max exon length: " + m_nMaxExonLength + "\n" +
					        "MMSEQ text length: " + m_nMMSeqTextLength + "\n" +
					        "intron length: " + m_nIntronLength + "\n" +
					        "total width: " + m_nTotalWidth + "\n";
			return strOut;
		}
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
	    new Color(  0, 255, 153, 127), 	// Teal-Teal-Cyan							21
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

	public PlotFactory(SplicingWebApp app)
	{
		m_app				  					= app;

		m_bCoveragePlotUseLog2					= false;
		m_bCoveragePlotShowRelativeCoverage 	= false;
		m_bCoveragePlotShowGeometricMean		= false;
		m_bCoveragePlotShowMedian				= true;
		m_bCoveragePlotShowMean					= false;
		m_bCoveragePlotShowQuartiles			= true;		
		
		m_bShowCoverageGrid						= true;
		m_bColorExonsAndJunctionsByCoverage		= true;
		
		m_bCoverageRequiresRedraw				= true;
		m_bSwitchStrands						= false;

		m_vcUnknownJunctions					= new TreeSet<CountElement>();
	}
	
	/** Returns the client window width and height */
	public int[] GetClientSize()
	{
		int pResult[] = new int[2];
		pResult[0] = m_nClientWindowWidth;
		pResult[1] = m_nClientWindowHeight;
		
		return pResult;
	}

	/**
	 *    Specifies which junctions are unannotated. These junctions
	 *    will be colored in a different color .
	 */
	public void SetUnknownJunctions(TreeSet<CountElement> vcJunctions)
	{
		m_vcUnknownJunctions.clear();
		m_vcUnknownJunctions.addAll(vcJunctions);
	}
	
	/** Adds an unannotated junction */
	public void AddUnknownJunction(CountElement jun)
	{
		m_vcUnknownJunctions.add(jun);
	}
	
	/** Clears the list of unannotated junctions */
	public void ClearUnknownJunctions()
	{
		m_vcUnknownJunctions.clear();
	}
	
	/**
	 *    Sets the client window width and height.
	 *    This function is invoked once, when the GUI is initialized.
	 */
	public void SetClientSize(int nClientWidth, int nClientHeight)
	{
		m_nMaxWidth			  = nClientWidth;
		m_nClientWindowWidth  = nClientWidth;
		m_nClientWindowHeight = nClientHeight;
	}
	
	/**
	 *   Sets the maximum plot width. By default this value equals the client window width, 
	 *   but it may be larger (or smaller).
	 */
	public void SetMaxWidth(int nWidth)
	{
		m_nMaxWidth = nWidth;
	}
	
	/** Show the coverage plot as log2 transformed values */
	public void ShowAsLog2(boolean bState)
	{
		m_bCoveragePlotUseLog2 = bState;
	}
	
	/** Show the mean coverage value for each base */
	public void ShowMeanCoverage(boolean bState)
	{
		m_bCoveragePlotShowMean = bState;
		
		if(bState == true)
		{
			m_bCoveragePlotShowMedian			= false;
			m_bCoveragePlotShowGeometricMean	= false;
		}
	}
	
	/** Show the geometric mean coverage value for each base */
	public void ShowGeometricMeanCoverage(boolean bState)
	{
		m_bCoveragePlotShowGeometricMean = bState;
		
		if(bState == true)
		{
			m_bCoveragePlotShowMean				= false;
			m_bCoveragePlotShowMedian			= false;
		}
	}
	
	/** Show the median coverage value for each base */
	public void ShowMedianCoverage(boolean bState)
	{
		m_bCoveragePlotShowMedian = bState;
		
		if(bState == true)
		{
			m_bCoveragePlotShowMean				= false;
			m_bCoveragePlotShowGeometricMean	= false;
		}
	}
	
	/** Show the relative coverage plot */
	public void ShowRelativeCoverage(boolean bState)
	{
		m_bCoveragePlotShowRelativeCoverage = bState;
	}
	
	/** Show the upper and lower quartiles. Not available if showing relative coverage values */
	public void ShowQuartiles(boolean bState)
	{
		m_bCoveragePlotShowQuartiles = bState;
	}
	
	/** Highlight unique isoform elements (e.g. exons only present in a single isoform) */
	public void ShowUniqueFeatures(boolean bState)
	{
		m_bShowUniqueFeatures		= bState;
		m_bShowAmbigiousUniqueExons = bState;
		
		m_bColorExonsAndJunctionsByCoverage = !bState;
	}
	
	/** Add expression dependent colors to exons and junctions */
	public void ShowColoredExonsAndJunctions(boolean bState)
	{
		m_bColorExonsAndJunctionsByCoverage = bState;
	}
	
	/** Show light gray grid in the coverage plots */
	public void ShowCoverageGrid(boolean bState)
	{
		m_bShowCoverageGrid = bState;
	}
	
	/** Set the color for a specific condition */
	public void UpdateColorSelection(String strCondition, Color clr)
	{
		m_mapColorsToConditions.put(strCondition, clr);
	}
	
	/** Request a redraw of the coverage on the next call of DrawPlots()*/
	public void RequestCoverageRedraw()
	{
		m_bCoverageRequiresRedraw = true;
	}
	
	/** Always show the gene in left (5') to right (3') direction. */
	public void ToggleSwitchStrand(boolean bState)
	{
		m_bSwitchStrands = bState;
	}
	
	/** Returns whether the relative coverage plot should be drawn */
	public boolean IsRelativeCoverageEnabled()
	{
		return m_bCoveragePlotShowRelativeCoverage;
	}
	
	/** Returns whether log2 coverage values should be plotted */
	public boolean IsLog2Enabled()
	{
		return m_bCoveragePlotUseLog2;
	}
	
	/** Returns the color for a specific condition */
	public Color GetColorForCondition(String strCondition)
	{
		if(m_mapColorsToConditions.containsKey(strCondition))
			return m_mapColorsToConditions.get(strCondition);
		
		return null;
	}
	
	/** Returns whether the median coverage should be plotted */
	public boolean IsMedianEnabled()
	{
		return m_bCoveragePlotShowMedian; 
	}
	
	/** Returns whether the mean coverage should be plotted */
	public boolean IsMeanEnabled()
	{
		return m_bCoveragePlotShowMean;
	}
	
	/** Returns whether the geometric mean coverage should be plotted */
	public boolean IsGeometricMeanEnabled()
	{
		return m_bCoveragePlotShowGeometricMean;
	}
	
	/** Returns whether the coverage quartiles should be plotted */
	public boolean IsQuartilesEnabled()
	{
		return m_bCoveragePlotShowQuartiles;
	}
	
	/** Retrives the color selection from the SplicingWebApp and stores it in the member variable m_mapColorsToConditions */
	public void GetColorsPerCondition()
	{
		ProjectModel model 					= m_app.GetProjectModel();
		String strSelectedConditionType 	= m_app.GetSelectedConditionType();
		
		if(m_mapColorsToConditions == null)
			m_mapColorsToConditions = new TreeMap<String, Color>();
		else
			m_mapColorsToConditions.clear();

		TreeMap<String, TreeSet<String>> vcSamplesAndConditions = model.GetSamplesPerCondition(strSelectedConditionType);
		
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
		
		// overwrite with any colors that might have been saved previously
		String strIn = model.GetFullPathOfProjectFile() + ".color_scheme";
		Scanner pIn = null;
		
		try
		{
			pIn = new Scanner(new File(strIn));
			
			while(pIn.hasNextLine())
			{
				String strLine = pIn.nextLine();
				String pSplit[] = strLine.split("\t");
				
				Color clr = new Color(Integer.parseInt(pSplit[1]), Integer.parseInt(pSplit[2]), Integer.parseInt(pSplit[3]));
				m_mapColorsToConditions.put(pSplit[0], clr);
			}
			pIn.close();
		}
		catch(FileNotFoundException e)
		{
			// don't throw an error
			return;
		}
	}
	
	/** Computes the drawing offsets based on the specified maximum width of the plots */
	private DrawingOffsets ComputeDrawingOffsets()
	{
		ProjectModel model 					= m_app.GetProjectModel();
		String strSelectedConditionType 	= m_app.GetSelectedConditionType();
		TreeSet<String> vcSelectedSamples 	= m_app.GetSelectedSamples();
		DataSupplier dataSupplier 			= m_app.GetDataSupplier();
		
		DrawingOffsets res = new DrawingOffsets();
		
		TreeMap<String, TreeSet<String>> vcSamplesAndConditions = model.GetSamplesPerCondition(strSelectedConditionType);
		
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
		Font fontBoldText	= new Font("Lucida Sans", Font.BOLD, 13);
		graph.setFont(fontBoldText);
		
		// define selected conditions
		TreeSet<String> vcSelectedConditions = new TreeSet<String>();
		for(String strCondition : vcSamplesAndConditions.keySet())
		{
			for(String strSample : vcSelectedSamples)
			{
				if(vcSamplesAndConditions.get(strCondition).contains(strSample))
					vcSelectedConditions.add(strCondition);
			}
		}
		
		// get longest condition name and adjust isoform offset if necessary
		int nIsoformNameOffset 	= 120; 	// save 100 px for the isoform and condition names (or more if condition names are longer)
		int nMMSeqTextLength 	= 0; 	// this size is added to the graph for MMseq results
		for(String strCondition : vcSelectedConditions)
		{
			FontMetrics fontMetrics = graph.getFontMetrics(fontBoldText);
			int nStringWidth = fontMetrics.stringWidth(strCondition);
			
			nIsoformNameOffset = Math.max(nIsoformNameOffset, nStringWidth+20);
			nMMSeqTextLength += nStringWidth + 20;
		}
		
		// calculate exon group positions (maximum screen size - transcript lengths - 2xborder - isoform name length)
		int nExonGroups = dataSupplier.GetExonGroups().length;
		int nMaxExonLength = nMaxWidth - (nIntronLength*nExonGroups) - (nOffsetX*2) - nIsoformNameOffset;
		
		res.m_nBorderOffset 		= nOffsetX;
		res.m_nIntronLength			= nIntronLength;
		res.m_nIsoformNameOffset	= nIsoformNameOffset;
		res.m_nMaxExonLength		= nMaxExonLength;
		res.m_nMMSeqTextLength		= nMMSeqTextLength;
		res.m_nTotalWidth			= nMaxWidth + nMMSeqTextLength;
		
		return res;
	}
	
	/**
	 *    Adds a y-axis ranging from 0 to fMaximum value to a graph at position fPosX/fPosY. The height of the axis is determined
	 *    by fHeight. The axis title is specified by strTitle.
	 */
	private double AddYAxis(double fPosX, double fPosY, double fHeight, double fMaximumValue, String strTitle, Graphics2D graph)
	{
		// draw axis
		Line2D.Double line = new Line2D.Double(fPosX, fPosY, fPosX, fPosY+fHeight);
		graph.draw(line);
		
		// set font
		Font fontNormalText		 = new Font("Lucida Sans", Font.PLAIN, 13);
		FontMetrics fontMetrics  = graph.getFontMetrics(fontNormalText);
		graph.setFont(fontNormalText);
		graph.setColor(Color.BLACK);
	
		// add y-axis label
		// set font
		int nStringWidth  = fontMetrics.stringWidth(strTitle);
		int nStringHeight  = fontMetrics.getHeight();
		AffineTransform trans = graph.getTransform();
		graph.translate((int)(fPosX-nStringHeight-30), (int)(fPosY + fHeight*0.5 + nStringWidth*0.5));
		graph.rotate(-Math.PI/2);
		graph.drawString(strTitle, 0,0);
		graph.setTransform(trans);

		// add ticks to y-axis
		double fStep = 0;
		
		if(fMaximumValue > 10)
		{
			fStep = Math.ceil(fMaximumValue / 10.0);	// results in ~10 ticks
			if(fStep > 9999) 		fStep = fStep - (fStep%1000);
			else if(fStep > 999) 	fStep = fStep - (fStep%100);
			else if(fStep > 10) 	fStep = fStep - (fStep%10);
			else					fStep = 10.0;
		}
		else
		{
			fStep = fMaximumValue / 10.0;
		}
		
		fontNormalText	= new Font("Lucida Sans", Font.PLAIN, 11);
		fontMetrics  	= graph.getFontMetrics(fontNormalText);
		graph.setFont(fontNormalText);
		
		fMaximumValue = Math.ceil(fMaximumValue / fStep)*fStep;

		// add ticks
		double fVal = 0.0;
		while(fVal <= fMaximumValue)
		{
			String strValue = "";
			if(fVal % 1 == 0 || fVal > 10)
			{
				strValue = String.format(Locale.ENGLISH, "%.0f", fVal);
			}
			else
			{
				strValue = String.format(Locale.ENGLISH, "%.2f", fVal);
			}
			
			double fYPos = fPosY + (fHeight - (fVal / fMaximumValue * fHeight));
			
			line = new Line2D.Double(fPosX-2, fYPos, fPosX, fYPos);
			graph.draw(line);
							
			nStringWidth  = fontMetrics.stringWidth(strValue);
			nStringHeight = fontMetrics.getHeight();
			
			graph.drawString(strValue, (int)(fPosX-3-nStringWidth), (int)(fYPos + nStringHeight*0.25));
			
			fVal += fStep;
		}
		
		return fMaximumValue;
	}
	
	/**
	 *    Adds an x-axis ranging from 0 to fMaximum value to a graph at position fPosX/fPosY. The width of the axis is determined
	 *    by fWidth. The axis title is specified by strTitle.
	 */
	private double AddXAxis(double fPosX, double fPosY, double fWidth, double fMaximumValue, String strTitle, Graphics2D graph)
	{
		// draw axis
		Line2D.Double line = new Line2D.Double(fPosX, fPosY, fPosX+fWidth, fPosY);
		graph.draw(line);
		
		// set font
		Font fontNormalText		 = new Font("Lucida Sans", Font.PLAIN, 13);
		FontMetrics fontMetrics  = graph.getFontMetrics(fontNormalText);
				
		// add ticks to y-axis
		double fStep = 0;
		
		if(fMaximumValue > 10)
		{
			fStep = Math.ceil(fMaximumValue / 10.0);	// results in ~10 ticks
			if(fStep > 9999) 		fStep = fStep - (fStep%1000);
			else if(fStep > 999) 	fStep = fStep - (fStep%100);
			else if(fStep > 10) 	fStep = fStep - (fStep%10);
			else					fStep = 10.0;
		}
		else
		{
			fStep = fMaximumValue / 10.0;
		}
		
		// add x-axis label
		int nStringWidth  = fontMetrics.stringWidth(strTitle);
		int nStringHeight = fontMetrics.getHeight();
		if(fMaximumValue > 99)
			graph.drawString(strTitle, (int)(fPosX + fWidth*0.5 - nStringWidth*0.5), (int)(fPosY+nStringHeight*2+20));
		else
			graph.drawString(strTitle, (int)(fPosX + fWidth*0.5 - nStringWidth*0.5), (int)(fPosY+nStringHeight+20));

		fontNormalText	= new Font("Lucida Sans", Font.PLAIN, 11);
		fontMetrics  	= graph.getFontMetrics(fontNormalText);
		nStringHeight = fontMetrics.getHeight();
		graph.setFont(fontNormalText);
		graph.setColor(Color.BLACK);
		
		// add ticks
		double fVal = 0.0;
		int nEvenOdd = 0;
		while(fVal <= fMaximumValue)
		{
			String strValue = "";
			if(fVal % 1 == 0)
			{
				strValue = String.format(Locale.ENGLISH, "%.0f", fVal);				
			}
			else
			{
				strValue = String.format(Locale.ENGLISH, "%.2f", fVal);
			}
			
			double fOffsetY = 0.0;
			if(fVal > 99 && nEvenOdd%2 == 1)
			{
				fOffsetY = nStringHeight;
			}

			double fXPos = fPosX + (int) Math.ceil(fVal / fMaximumValue * fWidth);
			
			line = new Line2D.Double(fXPos, fPosY, fXPos, fPosY+3);
			graph.draw(line);
							
			nStringWidth = fontMetrics.stringWidth(strValue);			
			graph.drawString(strValue, (int)(fXPos - nStringWidth*0.5), (int)(fPosY+nStringHeight+3+fOffsetY));
			
			fVal += fStep;
			nEvenOdd++;
		}
		
		return fMaximumValue;
	}
	
	/**
	 *    Draws the coverage plot. Different visualizations are possible depending on the options selected (e.g. mean or median
	 *    coverage plots).
	 */
	public void DrawExtendedCoveragePlot(boolean bAsPopupWindow) throws IOException
	{
		ProjectModel model 						= m_app.GetProjectModel();
		String strSelectedConditionType 		= m_app.GetSelectedConditionType();
		DataSupplier dataSupplier 				= m_app.GetDataSupplier();
		TreeMap<String, String> mapBigWigFiles 	= model.GetBigWigFilesForSamples();
		
		// switch strand for second strand genes
		boolean bShowSecondStrand = false;
		if(m_bSwitchStrands && !dataSupplier.GetGene().isPlusStrand())
			bShowSecondStrand = true;
		
		if(mapBigWigFiles == null)
		{
			System.out.println("ERROR: no valid bigwig or bam files detected");
			return;
		}
		
		//#########################################
		//    prepare popup window if requested
		//#########################################
		Window windowPopup = null;
		if(bAsPopupWindow)
		{
			windowPopup = new Window();
			windowPopup.setParent(m_app);			
			windowPopup.setTitle("Junction Heatmap");
			windowPopup.setSizable(true);
			windowPopup.setClosable(true);
			windowPopup.setMaximizable(true);
			windowPopup.setBorder(true);
			windowPopup.setLeft("10px");
			windowPopup.setTop("10px");
			windowPopup.setVisible(true);
			windowPopup.doPopup();
			windowPopup.setTopmost();
		}
		
		// get reference name
		String strRef = dataSupplier.GetReferenceName();
		
		// get samples per condition
		final TreeMap<String, TreeSet<String>> mapSamplesToConditions = model.GetSelectedSamplesPerCondition(strSelectedConditionType, m_app.GetSelectedSamples());
		
		ExonGroup pExonGroups[] = dataSupplier.GetExonGroups();
		if(m_vcUnknownJunctions.size() > 0)
		{
			TreeSet<Exon> vcExons = new TreeSet<Exon>();
			for(ExonGroup grp : pExonGroups)
			{
				for(Exon ex : grp.getExons())
					vcExons.add(ex);
			}
			
			// add potential novel exons
			int nExonID = 10000;
			for(CountElement jun : m_vcUnknownJunctions)
			{
				Exon exon = new Exon(jun.m_nStart-100, jun.m_nStart);
				exon.setID(nExonID);
				vcExons.add(exon);
				nExonID++;
				
				exon = new Exon(jun.m_nEnd, jun.m_nEnd+100);
				exon.setID(nExonID);
				vcExons.add(exon);
				
				nExonID++;
			}
	
			pExonGroups = dataSupplier.RecalculateExonGroups(vcExons);
		}

		//########################################
		//      prepare coverage containers
		//########################################
		// key 1 = exon_name, key 2 = sample, value = coverage array
		TreeMap<String, TreeMap<String, double[]>> mapCoverageToSamplesAndExons = new TreeMap<String, TreeMap<String, double[]>>();
		for(ExonGroup grp : pExonGroups)
		{			
			int nStart	= grp.getGenomicStartOfGroup();
			int nEnd	= grp.getGenomicStopOfGroup();
			int nLength = nEnd-nStart+1;
			
			String strName = strRef + ":" + nStart + "-" + nEnd;
			
			mapCoverageToSamplesAndExons.put(strName, new TreeMap<String, double[]>());

			for(String strCondition : mapSamplesToConditions.keySet())
			{
				for(String strSample: mapSamplesToConditions.get(strCondition))
				{
					TreeMap<String, double[]> tmp = mapCoverageToSamplesAndExons.get(strName);
					tmp.put(strSample, new double[nLength]);
				}
			}
		}
		
		//######################################################
		//                 prepare graph
		//######################################################
		DrawingOffsets offsets = ComputeDrawingOffsets();

		int nPlotHeight = 200;
		int nMargin 	= 10;
		int nSpacer		= 20;
		int nMaxHeight  = nPlotHeight+nMargin*2;	// 200 for the track, 20 as border
		
		if(bAsPopupWindow)
		{
			nMaxHeight = nMargin*2 + (nPlotHeight+nSpacer)*mapSamplesToConditions.size(); // 10 for the border, 210 for the plot and 20 as spacer between plots
			windowPopup.setWidth((offsets.m_nTotalWidth-10) +"px");
			windowPopup.setHeight((nMaxHeight+40) +"px");
		}
		
		// prepare graph
		BufferedImage img = new BufferedImage(offsets.m_nTotalWidth, nMaxHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D graph = img.createGraphics();
		graph.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		// set font size
		Font fontBoldText	= new Font("Lucida Sans", Font.BOLD, 13);
		Font fontNormalText = new Font("Lucida Sans", Font.PLAIN, 11);
		
		graph.setFont(fontBoldText);
		
		// fill background
		graph.setColor(Color.white);
		graph.fillRect(0, 0, offsets.m_nTotalWidth, nMaxHeight);
		
		// reserve some space for the isoform name
		int nOffsetX = offsets.m_nBorderOffset + offsets.m_nIsoformNameOffset;
		
		int nCombinedExonLength = 0;
		for(ExonGroup grp : pExonGroups)
		{
			nCombinedExonLength += grp.getExonGroupLengthInBp();
		}
		
		// draw condition name in associated color		
		int nOffset = 13;
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			graph.setColor(m_mapColorsToConditions.get(strCondition));
			graph.drawString(strCondition, 5, nMargin + nOffset);			
			
			if(bAsPopupWindow)
				nOffset += nSpacer + nPlotHeight;
			else
				nOffset += 13;
		}

		// calculate shrinkage factor for window width
		double fXShrinkageFactor = (double)offsets.m_nMaxExonLength / (double)nCombinedExonLength;
				
		//####################################################################
		//     obtain coverage for all exons and estimate maximum coverage
		//####################################################################
		double fMaxCoverage = 1.0;
		
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				for(ExonGroup grp : pExonGroups)
				{
					int nStart	= grp.getGenomicStartOfGroup();
					int nEnd	= grp.getGenomicStopOfGroup();
					String strName = strRef + ":" + nStart + "-" + nEnd;
					
					// size factor adjusted values
					double pCoverage[] = dataSupplier.GetCoverageArrayForExonGroupAndSample(m_app, grp, strSample);
					
					// if there is no specific coverage array for the exon group (in case the junction extends past the
					// gene boundaries) it is necessary to create one
					if(pCoverage == null)
					{
						pCoverage = m_app.GetCoverageForRegion(dataSupplier.GetReferenceName(), grp.getGenomicStartOfGroup(), grp.getGenomicStopOfGroup(), strSample, true);
					}
					
					if(m_bCoveragePlotUseLog2)
					{
						for(int i=0; i<pCoverage.length; i++)
						{
							pCoverage[i] = Math.log(pCoverage[i]+1) / Math.log(2);
						}
					}
					
					mapCoverageToSamplesAndExons.get(strName).put(strSample, pCoverage);
					
					if(!m_bCoveragePlotShowMedian && !m_bCoveragePlotShowGeometricMean && !m_bCoveragePlotShowMean)
						fMaxCoverage = Math.max(fMaxCoverage, StatUtils.max(pCoverage)); 
				}
			}
		}
				
		if(m_bCoveragePlotShowMedian || m_bCoveragePlotShowGeometricMean || m_bCoveragePlotShowMean)
		{
			for(ExonGroup grp : pExonGroups)
			{
				// get correct group coverage data
				String strName = strRef + ":" + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup();
				TreeMap<String, double[]> mapData = mapCoverageToSamplesAndExons.get(strName);
	
				// combine coverage per condition
				for(String strCondition : mapSamplesToConditions.keySet())
				{
					graph.setColor(this.m_mapColorsToConditions.get(strCondition));
					
					int nSamples = mapSamplesToConditions.get(strCondition).size();

					double[][] pValues =  new double[grp.getExonGroupLengthInBp()][nSamples];
					
					int nCurrentSample = 0;
					for(String strSample : mapSamplesToConditions.get(strCondition))
					{
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
		
		double fYShrinkageFactor = nPlotHeight / fMaxCoverage;

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
			graph.drawLine(nYAxisPos,  10, nYAxisPos, nPlotHeight);

			for(int i=0; i<11; i++)
			{
				int nYPos = (int) Math.ceil((nPlotHeight+nMargin) - (i*0.1 * nPlotHeight));
				
				graph.drawLine(nYAxisPos-2, nYPos, nYAxisPos, nYPos);
				
				String strValue = String.format(Locale.ENGLISH, "%1f", 0.1*i);
				FontMetrics fontMetrics = graph.getFontMetrics(fontNormalText);
				int nStringWidth = fontMetrics.stringWidth(strValue);
				int nStringHeight = fontMetrics.getHeight();
				
				graph.drawString(strValue, (int)(nYAxisPos-3-nStringWidth), (int)(nYPos + nStringHeight*0.25));
			}
			
			int nFirstExonGroup = 0;
			int nLastExonGroup  = pExonGroups.length;
			int nChange = +1;
			
			// flip the plot of we're showing the second strand
			if(bShowSecondStrand)
			{
				nFirstExonGroup = pExonGroups.length-1;
				nLastExonGroup  = -1;
				nChange = -1;
			}
			
			for(int nGrp=nFirstExonGroup; nGrp!=nLastExonGroup; nGrp+=nChange)
			{
				ExonGroup grp = pExonGroups[nGrp];
				
				int nGrpStart = nOffsetX;
				
				// get correct group coverage data
				String strName = strRef + ":" + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup();
				TreeMap<String, double[]> mapData = mapCoverageToSamplesAndExons.get(strName);
				
				int nWidth = (int)(grp.getExonGroupLengthInBp() * fXShrinkageFactor);
				
				// key = condition, value = coverage per position for each sample
				TreeMap<String, double[][]> mapCoveragePerSamplePerCondition = new TreeMap<String, double[][]>();

				for(String strCondition : mapSamplesToConditions.keySet())
				{
					int nSamples = mapSamplesToConditions.get(strCondition).size();

					if(nSamples > 0)
					{
						double[][] pValues =  new double[grp.getExonGroupLengthInBp()][nSamples];
						
						int nCurrentSample = 0;
						for(String strSample : mapSamplesToConditions.get(strCondition))
						{
							double pData[] = mapData.get(strSample);						
							for(int x=0; x<pData.length; x++)
							{
								double fValue = pData[x];
								pValues[x][nCurrentSample] = fValue;
							}
							
							nCurrentSample += 1;
						}
						
						// reverse array if second strand is to be plotted
						if(bShowSecondStrand)
						{
							ArrayUtils.reverse(pValues);
						}
						
						if(nCurrentSample != 0)
							mapCoveragePerSamplePerCondition.put(strCondition, pValues);
					}
				}
				
				for(int x=1; x<grp.getExonGroupLengthInBp(); x++)
				{			
					// get mean coverage for all conditions
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
						double fMean = StatUtils.mean(pCov);
						mapMedianCoverageToCondition.put(strCondition, fMean);
						fTotalCoverage += fMean;
					}
					
					if(fTotalCoverage == 0)
						continue;

					// calculate contribution of each condition to the coverage
					double fPlotBottom = nPlotHeight+nMargin;
					for(String strCondition : mapMedianCoverageToCondition.keySet())
					{
						graph.setColor(this.m_mapColorsToConditions.get(strCondition));						
						double fValue = (mapMedianCoverageToCondition.get(strCondition) / fTotalCoverage)*nPlotHeight;
						
						Shape rect = new Rectangle2D.Double(nGrpStart+(x*fXShrinkageFactor), fPlotBottom-fValue, 1, fValue);
						graph.draw(rect);
						fPlotBottom -= fValue;
					}
				}				

				nOffsetX += nWidth + offsets.m_nIntronLength;
			}
		}
		else
		{
			//######################################################
			//                   draw y-axis
			//######################################################
			graph.setFont(fontNormalText);
			graph.setColor(Color.BLACK);
			int nPlotBottom = nPlotHeight+nMargin;
			
			int nTimes = 1;
			if(bAsPopupWindow)
				nTimes = mapSamplesToConditions.size();
			
			Color clrLightGrey = new Color(245,245,245);
			
			int nOffsetY = 0;
			for(int i=0; i<nTimes; i++)
			{
				// y-axis
				int nYAxisPos = nOffsetX-10;
				graph.drawLine(nYAxisPos,  nMargin+nOffsetY, nYAxisPos, nPlotBottom+nOffsetY);
				int nStep = (int)Math.ceil(fMaxCoverage / 10.0);	// results in ~10 ticks
				
				if(nStep > 9999) nStep = nStep - (nStep%1000);
				else if(nStep > 999) nStep = nStep - (nStep%100);
				else if(nStep > 10) nStep = nStep - (nStep%10);
	
				int nVal = 0;		
				while(nVal < fMaxCoverage)
				{
					int nYPos = (int) Math.ceil(nPlotBottom - (nVal / fMaxCoverage * nPlotHeight));
					graph.drawLine(nYAxisPos-2, nYPos+nOffsetY, nYAxisPos, nYPos+nOffsetY);
					
					if(m_bShowCoverageGrid)
					{
						graph.setColor(clrLightGrey);
						graph.drawLine(nYAxisPos+1, nYPos+nOffsetY, nYAxisPos + offsets.m_nTotalWidth, nYPos+nOffsetY);
						graph.setColor(Color.black);
					}
					
					String strValue = ""+nVal;
					FontMetrics fontMetrics = graph.getFontMetrics(fontNormalText);
					int nStringWidth = fontMetrics.stringWidth(strValue);
					int nStringHeight = fontMetrics.getHeight();
					
					graph.drawString(strValue, (int)(nYAxisPos-3-nStringWidth), (int)(nYPos + nOffsetY + nStringHeight*0.25));
					
					nVal += nStep;
				}
				
				nOffsetY += nSpacer + nPlotHeight;
			}
			
			int nFirstExonGroup = 0;
			int nLastExonGroup  = pExonGroups.length;
			int nChange = +1;
			
			// flip the plot of we're showing the second strand
			if(bShowSecondStrand)
			{
				nFirstExonGroup = pExonGroups.length-1;
				nLastExonGroup  = -1;
				nChange = -1;
			}
			
			for(int nGrp=nFirstExonGroup; nGrp!=nLastExonGroup; nGrp+=nChange)
			{
				ExonGroup grp = pExonGroups[nGrp];
				int nGrpStart = nOffsetX;
				
				// get correct group coverage data
				String strName = strRef + ":" + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup();
				TreeMap<String, double[]> mapData = mapCoverageToSamplesAndExons.get(strName);
				
				int nWidth = (int)(grp.getExonGroupLengthInBp() * fXShrinkageFactor);
				
				TreeMap<String, double[][]> mapCoverageValuesToConditions = new TreeMap<String, double[][]>();
				
				//##########################################
				//##   combine coverage per condition     ##
				//##   and draw coverage for all samples  ##
				//##   if required.                       ##
				//##########################################
				for(String strCondition : mapSamplesToConditions.keySet())
				{
					graph.setColor(m_mapColorsToConditions.get(strCondition));
					
					int nSamples = mapSamplesToConditions.get(strCondition).size();
					
					if(nSamples > 0)
					{
						double[][] pValues =  new double[grp.getExonGroupLengthInBp()][nSamples];
						
						int nCurrentSample = 0;
	
						for(String strSample : mapSamplesToConditions.get(strCondition))
						{							
							double pData[] =  mapData.get(strSample);
							
							// reverse the array if drawing the second strand
							if(bShowSecondStrand)
							{
								ArrayUtils.reverse(pData);
							}
												
							GeneralPath polyline = new GeneralPath(GeneralPath.WIND_EVEN_ODD, pValues.length);

							for(int x=0; x<pData.length; x++)
							{
								double fValue = pData[x];
								
								if(m_bCoveragePlotShowMedian || m_bCoveragePlotShowGeometricMean || m_bCoveragePlotShowMean || bAsPopupWindow)
								{
									pValues[x][nCurrentSample] = fValue;
								}
								else
								{
									double fX = nGrpStart+(x*fXShrinkageFactor);
									double fNextX = nGrpStart+(x+1)*fXShrinkageFactor;
									
									if(x == 0)
									{
										polyline.moveTo(fX, nPlotBottom);
									}
									
									polyline.lineTo(fX, 	nPlotBottom-fValue*fYShrinkageFactor);
									polyline.lineTo(fNextX, nPlotBottom-fValue*fYShrinkageFactor);
								}
								
							}
							
							if(!m_bCoveragePlotShowMedian && !m_bCoveragePlotShowGeometricMean && !m_bCoveragePlotShowMean)
							{							
								// move back to 0
								polyline.lineTo(nGrpStart+(pData.length*fXShrinkageFactor), nPlotBottom);
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
						nOffsetY = 0;
						for(String strCondition : mapSamplesToConditions.keySet())
						{
							if(mapCoverageValuesToConditions.containsKey(strCondition))
							{
								double[][] pValues = mapCoverageValuesToConditions.get(strCondition);
							
								Color clrAlpha = new Color(m_mapColorsToConditions.get(strCondition).getRed(), m_mapColorsToConditions.get(strCondition).getGreen(), m_mapColorsToConditions.get(strCondition).getBlue(), 40);
								graph.setColor(clrAlpha);
								
								GeneralPath polyline = new GeneralPath(GeneralPath.WIND_EVEN_ODD, pValues.length);

								// upper quartile
								for(int x=0; x<pValues.length; x++)
								{									
									// get median coverage
									double fQ75 = 	StatUtils.percentile(pValues[x], 75);
									double fX = 	nGrpStart+(x*fXShrinkageFactor);
									double fNextX = nGrpStart+(x+1)*fXShrinkageFactor;
									
									if(x == 0)
									{
										polyline.moveTo(fX, nPlotBottom+nOffsetY-fQ75*fYShrinkageFactor);
									}
									else
									{
										polyline.lineTo(fX, nPlotBottom+nOffsetY-fQ75*fYShrinkageFactor);
									}
									
									polyline.lineTo(fNextX, nPlotBottom+nOffsetY-fQ75*fYShrinkageFactor);
								}

								// lower quartile
								for(int x=(pValues.length-1); x>=0; x--)
								{									
									// get median coverage
									double fQ25 = StatUtils.percentile(pValues[x], 25);
									double fX = nGrpStart+(x+1)*fXShrinkageFactor;
									double fNextX = nGrpStart+x*fXShrinkageFactor;
									
									polyline.lineTo(fX, 	nPlotBottom+nOffsetY-fQ25*fYShrinkageFactor);
									polyline.lineTo(fNextX, nPlotBottom+nOffsetY-fQ25*fYShrinkageFactor);
								}
								
								polyline.closePath();
								graph.fill(polyline);
							}
							
							if(bAsPopupWindow)
								nOffsetY += nSpacer+nPlotHeight;
						}
					}

					//####################################
					//##           draw lines           ##
					//####################################
					nOffsetY = 0;
					for(String strCondition : mapSamplesToConditions.keySet())
					{
						if(mapCoverageValuesToConditions.containsKey(strCondition))
						{
							double[][] pValues = mapCoverageValuesToConditions.get(strCondition);
			
							GeneralPath polyline = new GeneralPath(GeneralPath.WIND_EVEN_ODD, pValues.length);
							
							for(int x=0; x<pValues.length; x++)
							{
								// get median coverage
								double fQ50 	= 0.0;
								double fX 		= nGrpStart+(x*fXShrinkageFactor);
								double fNextX 	= nGrpStart+(x+1)*fXShrinkageFactor;
								
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
									polyline.moveTo(fX, nPlotBottom+nOffsetY);
								}
								
								polyline.lineTo(fX, 	nPlotBottom+nOffsetY-fQ50*fYShrinkageFactor);
								polyline.lineTo(fNextX, nPlotBottom+nOffsetY-fQ50*fYShrinkageFactor);
							}
							
							graph.setColor(this.m_mapColorsToConditions.get(strCondition));
							polyline.lineTo(nGrpStart+(pValues.length*fXShrinkageFactor), nPlotBottom+nOffsetY);
							graph.draw(polyline);
						}
						
						if(bAsPopupWindow)
							nOffsetY += nSpacer+nPlotHeight;
					}
				}

				nOffsetX += nWidth + offsets.m_nIntronLength;
			}
		}

		if(bAsPopupWindow)
		{
			Imagemap imgMapPopup = new Imagemap();			
			imgMapPopup.setWidth(offsets.m_nTotalWidth+"px");
			imgMapPopup.setHeight(nMaxHeight+"px");
			imgMapPopup.setContent(img);
			imgMapPopup.setParent(windowPopup);
		}
		else
		{
			Imagemap imgMapCoverage = m_app.GetImageMapForCoveragePlot();
			imgMapCoverage.getChildren().clear();
			imgMapCoverage.setWidth(offsets.m_nTotalWidth+"px");
			imgMapCoverage.setHeight(nMaxHeight+"px");
			imgMapCoverage.setContent(img);
			imgMapCoverage.setParent(m_app.GetCoveragePlotRegion());
			
			imgMapCoverage.addEventListener(Events.ON_CLICK, new EventListener<Event>()
			{
				@Override
				public void onEvent(Event event) throws Exception
				{
					MouseEvent evnt = (MouseEvent) event;
					evnt.stopPropagation();

					DrawExtendedCoveragePlot(true);
				}
			});
		}
		
		m_bCoverageRequiresRedraw = false;
	}

	/**
	 *    Draws the isoform plot. Different information levels are added depending on the options selected.
	 *    For example, tissue specific expression (entropy based) or exon skipping events  may be indicated
	 *    on the meta exon track. The function became quite large, I might consider to split it into multiple
	 *    smaller functions...
	 */
	public void DrawIsoforms() throws IOException
	{
		ProjectModel projectModel				= m_app.GetProjectModel();
		AnalysisResultHandler resultHandler 	= m_app.GetResultHandler();
		DataSupplier dataSupplier 				= m_app.GetDataSupplier();
		String strSelectedConditionType 		= m_app.GetSelectedConditionType();
		String strSelectedCondition				= m_app.GetSelectedCondition();
		TreeSet<String> vcSelectedSamples 		= m_app.GetSelectedSamples();
		TreeSet<String> vcValidIsoforms			= m_app.GetValidIsoforms();
		AnalysisResult  selectedResult			= m_app.GetSelectedResult();
				
		TreeMap<String, Double> mapEntropyToExonicParts = m_app.GetEntropyPerExonicPart();
		
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = projectModel.GetSamplesPerCondition(strSelectedConditionType);
		
		// switch strand for second strand genes
		boolean bShowSecondStrand = false;
		if(m_bSwitchStrands && !dataSupplier.GetGene().isPlusStrand())
			bShowSecondStrand = true;
		
		TreeMap<String, Exon[]> mapExonsToUnknownIsoforms = new TreeMap<String, Exon[]>();
		ExonGroup pExonGroups[] = dataSupplier.GetExonGroups();
		if(m_vcUnknownJunctions.size() > 0)
		{
			TreeSet<Exon> vcExons = new TreeSet<Exon>();
			for(ExonGroup grp : pExonGroups)
			{
				for(Exon ex : grp.getExons())
					vcExons.add(ex);
			}
			
			// add potential novel exons
			int nExonID = 10000;
			int nIsoformID = 1;
			for(CountElement jun : m_vcUnknownJunctions)
			{
				Exon pNewExons[] = new Exon[2];
				Exon exon = new Exon(jun.m_nStart-100, jun.m_nStart);
				exon.setID(nExonID);
				vcExons.add(exon);
				pNewExons[0] = exon;
				nExonID++;
				
				exon = new Exon(jun.m_nEnd, jun.m_nEnd+100);
				exon.setID(nExonID);
				vcExons.add(exon);
				pNewExons[1] = exon;
				
				mapExonsToUnknownIsoforms.put("UNKNOWN_" + nIsoformID, pNewExons);
				
				nExonID++;
				nIsoformID++;
			}
	
			pExonGroups = dataSupplier.RecalculateExonGroups(vcExons);
		}

		double fMaxJunCoverage = 0.0;
		if(m_bColorExonsAndJunctionsByCoverage)
		{
			// get maximum junction coverage
			for(CountElement jun : dataSupplier.GetJunctions())
			{
				// get junction counts
				TreeMap<String, Integer> mapCountsToSamples = dataSupplier.GetCountsForJunction(jun);
				
				// get median for the current condition
				int nSamples = 0;
				
				// determine number of valid samples
				for(String strSample : mapSamplesToConditions.get(strSelectedCondition))
				{
					if(vcSelectedSamples.contains(strSample))
						nSamples++;
				}
				
				double pCoverages[] = new double[nSamples];			
				double fCoverage = 0.0;
				
				int nIdx = 0;
				for(String strSample : mapSamplesToConditions.get(strSelectedCondition))
				{
					if(vcSelectedSamples.contains(strSample))
					{
						if(mapCountsToSamples.containsKey(strSample))
							pCoverages[nIdx] = mapCountsToSamples.get(strSample);
						else
							pCoverages[nIdx] = 0.0;
						nIdx++;
					}
				}
				
				fCoverage = StatUtils.percentile(pCoverages, 50.0);
				fMaxJunCoverage = Math.max(fMaxJunCoverage, fCoverage);
			}
		}
		
		// get mmseq results (if available)
		TreeMap<String, TreeMap<String, Double>> mapIsoformExpressionValues = projectModel.GetIsoformExpressions(dataSupplier.GetGeneID(), m_app.GetPathToMMSeqResults());

		int nVerticalSpaceBetweenIsoforms = 30;
		
		// if numbers are added to the graph, we need more space between the isoforms
		if(m_bColorExonsAndJunctionsByCoverage)
			nVerticalSpaceBetweenIsoforms += 20;
		
		// define selected conditions
		TreeSet<String> vcSelectedConditions = new TreeSet<String>();
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			for(String strSample : vcSelectedSamples)
			{
				if(mapSamplesToConditions.get(strCondition).contains(strSample))
					vcSelectedConditions.add(strCondition);
			}
		}
		
		DrawingOffsets offsets = ComputeDrawingOffsets();
		
		// an additional 30 pixel as margin, and three more isoforms for the meta transcript and overlapping sense + antisense transcripts
		int nMaxHeight  = (vcValidIsoforms.size()+mapExonsToUnknownIsoforms.size()+3) * nVerticalSpaceBetweenIsoforms + 20 + 30;
		
		// prepare graph
		BufferedImage img	= new BufferedImage(offsets.m_nTotalWidth, nMaxHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D graph	= img.createGraphics();
		graph.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
		graph.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		Imagemap imgMapIsoforms = m_app.GetImageMapForIsoforms();
		imgMapIsoforms.getChildren().clear();
		
		// fill background
		graph.setColor(Color.white);
		graph.fillRect(0, 0, offsets.m_nTotalWidth, nMaxHeight);
		
		// reserve some space for the isoform name
		int nOffsetX = offsets.m_nIsoformNameOffset + offsets.m_nBorderOffset;

		int nCombinedExonLength = 0;
		for(ExonGroup grp : pExonGroups)
		{
			nCombinedExonLength += grp.getExonGroupLengthInBp();
		}

		// calculate shrinkage factor
		double fShrinkageFactor = (double)offsets.m_nMaxExonLength / (double)nCombinedExonLength;
		
		// calculate start position (=value) of each exon (=key), this is required when drawing first strand transcripts
		TreeMap<Integer, Integer> mapExonStarts = new TreeMap<Integer, Integer>(); 	// key = exon_id, value = exon_start_pos
			
		int nX = nOffsetX;
		for(int nGrp=0; nGrp<pExonGroups.length; nGrp++)
		{
			ExonGroup grp = pExonGroups[nGrp];
			if(bShowSecondStrand)
			{
				grp = pExonGroups[pExonGroups.length - nGrp -1];
			}
			int nGrpStart = nX;
			
			for(Exon ex : grp.getExons())
			{
				int nExonID = ex.getExonID();

				if(bShowSecondStrand)
				{
					// all exons must end at the same position
					
					double fVal1 = (grp.getGenomicStopOfGroup() - grp.getGenomicStartOfGroup()) * fShrinkageFactor;
					double fVal2 = (ex.getGenomicStop() - grp.getGenomicStartOfGroup()) * fShrinkageFactor;
					
					int nOffset = (int) (fVal1 - fVal2);

					mapExonStarts.put(nExonID, nGrpStart + nOffset);
				}
				else
				{
					// all exons must start at the same position
					int nOffset = (int) ((ex.getGenomicStart() - grp.getGenomicStartOfGroup()) * fShrinkageFactor);
					mapExonStarts.put(nExonID, nGrpStart + nOffset);
				}
			}

			nX += (grp.getExonGroupLengthInBp() * fShrinkageFactor) + offsets.m_nIntronLength;
		}

		Color clrGradient3 = new Color(217,217,255);	// BLUE
		Color clrGradient4 = new Color(120,120,230);	// BLUE
		Color clrGradient5 = new Color(255,255,255);	// WHITE
		Color clrGradient6 = new Color(150,150,150);	// GRAY
		
		Color clrGradient7 = new Color(255,240,200);	// LIGHT ORANGE
		Color clrGradient8 = new Color(230,180,0);		// ORANGE
		Color clrGradient9 = new Color(200, 80, 0);		// DAKR ORANGE

		//###################################
		//              prepare font
		//###################################
		Font fontNumbers = new Font("Lucida Sans", Font.BOLD, 11);
		graph.setFont(fontNumbers);
		FontMetrics fontMetrics = graph.getFontMetrics(fontNumbers);
		
		// draw condition name in associated color
		graph.setColor(m_mapColorsToConditions.get(strSelectedCondition));
		graph.drawString(strSelectedCondition, 5, 12);
		
		//###################################
		//    load incomplete exon images
		//###################################
		String strPath = Executions.getCurrent().getDesktop().getWebApp().getRealPath("/");
		BufferedImage imgMark = null;
		try
		{
			imgMark = ImageIO.read(new File(strPath + "/img/noncontinuous_left.png"));
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.out.println("could not open image: " + strPath + "/img/noncontinuous_left.png");
			return;
		}
		BufferedImage imgMarkLeft = new BufferedImage(42, 99, BufferedImage.TYPE_INT_ARGB);
		
		Graphics2D g = imgMarkLeft.createGraphics();
		g.drawImage(imgMark, 0, 0, 20, 30, null);
		g.dispose();
		
		try
		{
			imgMark = ImageIO.read(new File(strPath + "/img/noncontinuous_right.png"));
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.out.println("could not open image: " + strPath + "/img/noncontinuous_right.png");
			return;
		}
		BufferedImage imgMarkRight = new BufferedImage(42, 99, BufferedImage.TYPE_INT_ARGB);
		
		g = imgMarkRight.createGraphics();
		g.drawImage(imgMark, 0, 0, 20, 30, null);
		g.dispose();

		//#######################################################
		//    get list of unique features of the current gene
		//#######################################################
		TreeSet<Exon>  vcUniqueExons		= dataSupplier.GetUniqueExons(null);
		TreeSet<String> vcUniqueJunctions	= dataSupplier.GetUniqueJunctions(null);

		//###################################
		//          specify offsets
		//###################################
		int nY = 20; 	// top margin
		nX = nOffsetX;	// left margin
		
		//##############################################
		//         add overlapping transcripts
		//##############################################
		nY = DrawOverlappingTranscripts(nX, nY, offsets.m_nIntronLength, fShrinkageFactor, nVerticalSpaceBetweenIsoforms, graph);
		
		//#######################################################
		//           add description to meta transcript
		//#######################################################
		nX = nOffsetX;

		graph.setColor(Color.BLACK);
		String strText = "meta exons";
		int nTextWidth = fontMetrics.stringWidth(strText);
		int nTextHeight = fontMetrics.getHeight();
		graph.drawString(strText, nX - nTextWidth - 5, nY + 10 + (int)(nTextHeight*0.25));
		
		//##############################################
		//            draw meta transcript
		//##############################################
		nX = nOffsetX;
		int nFirstGroup  = 0;
		int nLastGroup   = pExonGroups.length;
		int nGroupChange = +1;
		
		if(bShowSecondStrand)
		{
			nFirstGroup  = pExonGroups.length-1;
			nLastGroup   = -1;
			nGroupChange = -1;
		}
		
		for(int nGrp=nFirstGroup; nGrp!=nLastGroup; nGrp+=nGroupChange)
		{
			ExonGroup grp = pExonGroups[nGrp];
			
			int nExonGrpID = grp.getGroupID()+1;

			int nWidth = (int)(grp.getExonGroupLengthInBp() * fShrinkageFactor);
			
			Color pClrs[] = new Color[3];
			pClrs[0] = new Color(120, 120, 120);
			pClrs[1] = new Color(200, 200, 200);
			pClrs[2] = new Color(120, 120, 125);
			
			float pFractions[] = new float[3];
			pFractions[0] = 0.0f;
			pFractions[1] = 0.5f;
			pFractions[2] = 1.0f;

			// show alternatively spliced exon groups			
			AnalysisResult res = null;
			for(Exon ex : grp.getExons())
			{
				AnalysisResult res2 = resultHandler.GetASResultForExon(dataSupplier.GetReferenceName(), ex, SplicingWebApp.AS_TYPE_EXON_SKIPPING);

				if(res2 != null && res2.GetResultTypeAsString().equals("combined"))
				{
					if( (res2 != null && res == null) || (res != null && res2.GetPValueA() < res.GetPValueA()))
					{
						res = res2;
					}
				}
			}
			
			if(res != null)
			{
				pClrs[0] = new Color(255, 180, 180);
				pClrs[1] = new Color(255,  60,  60);
				pClrs[2] = new Color(255, 180, 180);

				pFractions[0] = 0.0f;
				pFractions[1] = 0.5f;
				pFractions[2] = 1.0f;
			}
			
			LinearGradientPaint clrGradient = new LinearGradientPaint(nX, nY, nX, nY+20, pFractions, pClrs);
				
			// draw exon
			graph.setPaint(clrGradient);
			graph.fillRect(nX, nY, nWidth, 20);
			
			// add entropy levels if available
			if(m_app.IsEntropyEnabled() && mapEntropyToExonicParts != null)
			{
				// clear gradient
				graph.setPaint(Color.WHITE);
				graph.fillRect(nX, nY, nWidth, 20);
				
				for(Map.Entry<String, Double> e : mapEntropyToExonicParts.entrySet())
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
						area.setParent(imgMapIsoforms);
					}
				}
			}
			
			// draw exon border
			graph.setColor(Color.BLACK);
			graph.drawRect(nX, nY, nWidth, 20);
			
			fontNumbers = new Font("Lucida Sans", Font.BOLD, 10);
			graph.setFont(fontNumbers);
			fontMetrics = graph.getFontMetrics(fontNumbers);
			
			// add orientation
			strText          = "";
			nTextWidth 		 = 0;
			int nLetterWidth = fontMetrics.stringWidth(">");
			
			while(nTextWidth < nWidth-20)
			{
				char chStrand = '+';
				if(bShowSecondStrand)
					chStrand = '-';
				
				if(dataSupplier.GetStrand() == chStrand)
				{
					strText += ">";
					nTextWidth += nLetterWidth;
				}
				else
				{
					strText += "<";
					nTextWidth += nLetterWidth;
				}
			}

			nTextHeight = fontMetrics.getHeight();
			
			int nXOffset = (int)((nWidth - nTextWidth) * 0.5);
			
			graph.drawString(strText, nX + nXOffset, (int)(nY+6+nTextHeight*0.5));
			
			// add clickable area
			Area area = new Area();
			if(res != null)
				area.setId("highlighted_exon_group " + nExonGrpID + " " + grp.getGenomicStartOfGroup() + " " + grp.getGenomicStopOfGroup());
			else
				area.setId("exon_group " + nExonGrpID + " " + grp.getGenomicStartOfGroup() + "-" + grp.getGenomicStopOfGroup());
			area.setShape("rectangle");
			area.setCoords(nX + ", " + nY + ", " + (nX+nWidth) + ", " + (nY+20));
			
			String strToolTip = String.format(Locale.ENGLISH, "exon group %d\n%s: %,d - %,d\n", nExonGrpID, dataSupplier.GetReferenceName(), grp.getGenomicStartOfGroup(), grp.getGenomicStopOfGroup());
			
			if(res != null)
			{
				strToolTip += String.format(Locale.ENGLISH, "ratio change: %.4f", res.GetAbsoluteChangeA());
			}
			
			area.setTooltiptext(strToolTip);
			area.setParent(imgMapIsoforms);
			
			nX += nWidth + offsets.m_nIntronLength;
		}
		
		nY += nVerticalSpaceBetweenIsoforms;
		//##############################################
		//               draw isoforms
		//##############################################
		TreeSet<String> vcIsoforms = new TreeSet<String>();
		vcIsoforms.addAll(vcValidIsoforms);
		vcIsoforms.addAll(mapExonsToUnknownIsoforms.keySet());

		for(String strIsoform : vcIsoforms)
		{
			int nPreviousExonEnd   = -1;
			int nPreviousExonGrpID = -1;
			Exon previousExon      = null;
			boolean bIsFirstExon   = true;
			
			//#######################################################
			//            add icon to unselect the isoform
			//#######################################################
			String strImageString = Executions.getCurrent().getDesktop().getWebApp().getRealPath("/img/red_cross.png");
			BufferedImage img2 = ImageIO.read(new File(strImageString));
			graph.drawImage(img2, 2, nY -4 + (int)(nTextHeight*0.5), 16, 16, null);
			
			Area area = new Area();
			area.setId("remove_" + strIsoform);
			area.setShape("rectangle");
			area.setCoords(2 + ", " + (nY -4 + (int)(nTextHeight*0.5)) + ", " + (2+16) + ", " + (nY-4 + (int)(nTextHeight*0.5)+16));
			area.setTooltiptext("remove isoform: " + strIsoform);
			area.setParent(imgMapIsoforms);
			
			//#######################################################
			//                add transcript name
			//#######################################################
			graph.setColor(Color.BLACK);
			nTextWidth = fontMetrics.stringWidth(strIsoform);
			nTextHeight = fontMetrics.getHeight();
			graph.drawString(strIsoform, nOffsetX - nTextWidth - 5, nY + 8 + (int)(nTextHeight*0.5));
			
			Exon[] pIsoformExons = null;
			boolean bIsNovelIsoform = false;
			
			if(mapExonsToUnknownIsoforms.containsKey(strIsoform))
			{
				pIsoformExons = mapExonsToUnknownIsoforms.get(strIsoform);
				bIsNovelIsoform = true;
			}
			else
			{
				pIsoformExons = dataSupplier.GetExonsForIsoform(strIsoform);
			}
			
			// get maximum exon value
			double fMaxAbsoluteCoverage = 0.0;
			if(m_bColorExonsAndJunctionsByCoverage)
			{
				for(Exon exon : pIsoformExons)
				{					
					TreeMap<String, Double> mapAbsoluteCoveragesPerSample = dataSupplier.GetAbsoluteCoveragePerSample(m_app, exon);
					
					int nSamples = 0;
					for(String strSample : mapSamplesToConditions.get(strSelectedCondition))
					{
						if(vcSelectedSamples.contains(strSample))
							nSamples++;
					}
					double pCoverages[] = new double[nSamples];
					
					int nIdx = 0;
					for(String strSample : mapSamplesToConditions.get(strSelectedCondition))
					{
						if(vcSelectedSamples.contains(strSample))
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
			
			int nFirstExon = 0;
			int nLastExon  = pIsoformExons.length;
			int nChange    = +1;
			if(bShowSecondStrand)
			{
				nFirstExon = pIsoformExons.length-1;
				nLastExon  = -1;
				nChange    = -1;
			}
			
			for(int nExonIdx=nFirstExon; nExonIdx!=nLastExon; nExonIdx+=nChange)
			{
				Exon exon = pIsoformExons[nExonIdx];

				// get ambigious exons within an exon group per condition
				TreeMap<String, TreeSet<Exon>> mapAmbigiousExonsPerCondition = dataSupplier.GetMostFrequentSetOfAmbigiousExons(m_app, exon);
				
				// get exon group ID
				int nExonGrpID = -1;
				for(ExonGroup grp : pExonGroups)
				{
					if(grp.groupContainsExon(exon.getExonID()))
					{
						nExonGrpID = grp.getGroupID()+1;
					}
				}
				
				boolean bIsHighlighted = false;
				if(selectedResult != null)
				{
					if(selectedResult.HasAltExonA())
					{
						if(exon.getGenomicStart() == selectedResult.GetStartA() && exon.getGenomicStop() == selectedResult.GetEndA())
							bIsHighlighted = true;
					}
					
					if(selectedResult.HasAltExonB())
					{
						if(exon.getGenomicStart() == selectedResult.GetStartB() && exon.getGenomicStop() == selectedResult.GetEndB())
							bIsHighlighted = true;
					}
				}

				nX = mapExonStarts.get(exon.getExonID());
				int nWidth = (int)(exon.getLength() * fShrinkageFactor);

				// exon fill color depends on strand orientation
				GradientPaint clrGradient = null;
				clrGradient = new GradientPaint((int)(nX+nWidth*0.5), nY, clrGradient3, (int)(nX+nWidth*0.5), nY+20, clrGradient4);
				
				if(m_bShowUniqueFeatures)
				{
					if(vcUniqueExons.contains(exon))
					{
						if(mapAmbigiousExonsPerCondition != null && m_bShowAmbigiousUniqueExons)
						{
							TreeSet<Exon> vcAmbigiousExons = mapAmbigiousExonsPerCondition.get(strSelectedCondition);
						
							if(vcAmbigiousExons.size() == 1)
							{
								clrGradient = new GradientPaint((int)(nX+nWidth*0.5), nY, clrGradient7, (int)(nX+nWidth*0.5), nY+20, clrGradient8);
							}
							else
							{
								for(Exon ex : vcAmbigiousExons)
								{
									if(ex.equals(exon))
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
				if(m_bColorExonsAndJunctionsByCoverage)
				{
					graph.setColor(Color.BLACK);
					
					TreeMap<String, Double> mapAbsoluteCoveragesPerSample = dataSupplier.GetAbsoluteCoveragePerSample(m_app, exon);
					
					//##################################
					//            add text
					//##################################
					strText = "";
					
					int nSamples = 0;
					for(String strSample : mapSamplesToConditions.get(strSelectedCondition))
					{
						if(vcSelectedSamples.contains(strSample))
							nSamples++;
					}
					
					double pCoverages[] = new double[nSamples];
					
					int nIdx = 0;
					for(String strSample : mapSamplesToConditions.get(strSelectedCondition))
					{
						if(vcSelectedSamples.contains(strSample))
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
				
				// do not color exons that are not highlighted, if there are any highlighted exons
				if(!bIsHighlighted && selectedResult != null)
				{
					graph.setPaint(Color.WHITE);
				}
				
				//###################################
				//             draw exon
				//###################################
				
				// determine whether the exon is a UTR exon
				int nCodingStart = dataSupplier.GetCodingStartForIsoform(strIsoform);
				int nCodingEnd	 = dataSupplier.GetCodingEndForIsoform(strIsoform);

				// full UTR exon
				if( (dataSupplier.GetStrand() == '+' && (exon.getGenomicStart() > nCodingEnd   || exon.getGenomicStop() < nCodingStart)) ||
					(dataSupplier.GetStrand() != '+' && (exon.getGenomicStart() > nCodingStart || exon.getGenomicStop() < nCodingEnd))   ||
					(nCodingStart == Integer.MIN_VALUE && nCodingEnd == Integer.MAX_VALUE)     || (nCodingStart == Integer.MAX_VALUE && nCodingEnd == Integer.MIN_VALUE))
				{
					graph.fillRect(nX, nY+4, nWidth, 12);

					graph.setColor(Color.BLACK);
					graph.drawRect(nX, nY+4, nWidth, 12);
				}
				// full CDS exon
				else if(((dataSupplier.GetStrand() == '+' && exon.getGenomicStart() >= nCodingStart && exon.getGenomicStop() <= nCodingEnd) ||
						 (dataSupplier.GetStrand() != '+' && exon.getGenomicStart() >= nCodingEnd   && exon.getGenomicStop() <= nCodingStart)))
				{
					graph.fillRect(nX, nY, nWidth, 20);

					graph.setColor(Color.BLACK);
					graph.drawRect(nX, nY, nWidth, 20);
				}
				// partial CDS exon
				else
				{
					Polygon polyExon = new Polygon();
					char chStrand = '+';
					
					int nMidPoint = 0;
					boolean bIsCoding = false;
					
					if(nCodingStart >= exon.getGenomicStart() && nCodingStart <= exon.getGenomicStop())
					{
						nMidPoint = nCodingStart - exon.getGenomicStart();
						
						if(bShowSecondStrand)
						{
							nMidPoint = (exon.getGenomicStop() -exon.getGenomicStart()+1) - (nCodingStart - exon.getGenomicStart());
							chStrand = '-';
						}
						
						if(dataSupplier.GetStrand() == chStrand)
						{
						
							if(exon.getGenomicStart() >= nCodingStart)
								bIsCoding = false;
							else
								bIsCoding = true;
						}
						else
						{
							if(exon.getGenomicStart() >= nCodingStart)
								bIsCoding = true;
							else
								bIsCoding = false;
						}
					}
					else
					{
						nMidPoint = nCodingEnd - exon.getGenomicStart();
						
						if(bShowSecondStrand)
						{
							nMidPoint = (exon.getGenomicStop() -exon.getGenomicStart()+1) - (nCodingEnd - exon.getGenomicStart());
							chStrand = '-';
						}
						
						if(dataSupplier.GetStrand() == chStrand)
						{
							if(exon.getGenomicStart() >= nCodingEnd)
								bIsCoding = true;
							else
								bIsCoding = false;
						}
						else
						{
							if(exon.getGenomicStart() >= nCodingEnd)
								bIsCoding = false;
							else
								bIsCoding = true;
						}
					}
					
					int nMidPointCoord = nX + (int)(nMidPoint * fShrinkageFactor);

					// draw top|left to top|middle
					if(!bIsCoding)
					{
						polyExon.addPoint(nX, nY);
						polyExon.addPoint(nMidPointCoord, nY);
						bIsCoding = false;
					}
					else
					{
						polyExon.addPoint(nX, nY+4);
						polyExon.addPoint(nMidPointCoord, nY+4);
						bIsCoding = true;
					}
					
					// draw top|middle to top|right to bottom|right
					if(bIsCoding)
					{
						polyExon.addPoint(nMidPointCoord, nY);
						polyExon.addPoint(nX + nWidth, nY);
						polyExon.addPoint(nX + nWidth, nY+20);
					}
					else
					{
						polyExon.addPoint(nMidPointCoord, nY+4);
						polyExon.addPoint(nX + nWidth, nY+4);
						polyExon.addPoint(nX + nWidth, nY+16);
					}

					// draw bottom|right to bottom|middle
					if(bIsCoding)
					{
						polyExon.addPoint(nMidPointCoord, nY+20);
						polyExon.addPoint(nMidPointCoord, nY+16);
						bIsCoding = false;
					}
					else
					{
						polyExon.addPoint(nMidPointCoord, nY+16);
						polyExon.addPoint(nMidPointCoord, nY+20);
						bIsCoding = true;
					}
					
					// draw bottom|middle to bottom|left
					if(bIsCoding)
					{
						polyExon.addPoint(nX, nY+20);
						
						// close the polygon
						polyExon.addPoint(nX, nY);
					}
					else
					{
						polyExon.addPoint(nX, nY+16);
						
						// close the polygon
						polyExon.addPoint(nX, nY+4);
					}
					
					graph.fillPolygon(polyExon);
					graph.setColor(Color.BLACK);
					graph.drawPolygon(polyExon);
				}

				area = new Area();
				area.setId("exon_" + nExonGrpID + " " + strIsoform + " " + exon.getExonID());
				area.setShape("rectangle");
				area.setCoords(nX + ", " + nY + ", " + (nX+nWidth) + ", " + (nY+20));				
//				area.setTooltiptext(strIsoform + " exon_" + nExonGrpID + ", " + exon.getGenomicStart() + "-" + exon.getGenomicStop() + ", (" + exon.getExonID() + ")");
				area.setTooltiptext(String.format(Locale.ENGLISH, "%s\nexon %d\n%s: %,d - %,d (%d)", strIsoform, nExonGrpID, dataSupplier.GetReferenceName(), exon.getGenomicStart(), exon.getGenomicStop(), exon.getExonID()));
				area.setParent(imgMapIsoforms);
				
				// add marker for incomplete exons
				if(bIsNovelIsoform)
				{
					if(nExonIdx == 0)
					{
						int nPos = nX - 14;
						graph.drawImage(imgMarkLeft, nPos, nY-5, null);
					}
					else
					{
						int nPos = (nX+nWidth) - 6;
						graph.drawImage(imgMarkRight, nPos, nY-5, null);
					}
				}

				//#############################################
				//               draw introns 
				//#############################################
				if(bIsFirstExon)
					bIsFirstExon = false;
				else
				{
					int nIntronStart  = nPreviousExonEnd;
					int nIntronLength = nX-nPreviousExonEnd;	// = length - 1
					
					CountElement junction = new CountElement();
					
					if(bShowSecondStrand)
					{
						junction.m_nStart 	= exon.getGenomicStop();
						junction.m_nEnd		= previousExon.getGenomicStart();
					}
					else
					{
						junction.m_nStart 	= previousExon.getGenomicStop();
						junction.m_nEnd		= exon.getGenomicStart();
					}
					
					if(strIsoform.contains("UNKNOWN"))
						junction.m_chStrand = '?';
					else
						junction.m_chStrand = dataSupplier.GetStrand();
					String strJunctionName = (junction.m_nStart+1) + "-" + (junction.m_nEnd-1);
					
					clrGradient = new GradientPaint((int)(nIntronStart+(nIntronLength)*0.5), nY+7, clrGradient5, (int)(nIntronStart+(nIntronLength)*0.5), nY+13, clrGradient6);
					
					if(m_bShowUniqueFeatures)
					{
						if(vcUniqueJunctions.contains(strJunctionName))
							clrGradient = new GradientPaint((int)(nIntronStart+(nIntronLength)*0.5), nY+7, clrGradient7, (int)(nIntronStart+(nIntronLength)*0.5), nY+13, clrGradient8);
					}
					
					bIsHighlighted = false;
					if(selectedResult != null && selectedResult.HasPSIScore())
					{
						SimpleSpliceScore score = selectedResult.GetPSIScore();
						String pSplit[] = strJunctionName.split("-");
						int nJunctionStart = Integer.parseInt(pSplit[0].trim())-1;
						int nJunctionEnd = Integer.parseInt(pSplit[1].trim())+1;
						
						if(score.m_JunctionInclusion.m_nStart == nJunctionStart && score.m_JunctionInclusion.m_nEnd == nJunctionEnd)
							bIsHighlighted = true;
						
						if(score.m_JunctionExclusion.m_nStart == nJunctionStart && score.m_JunctionExclusion.m_nEnd == nJunctionEnd)
							bIsHighlighted = true;
					}

					if(m_bColorExonsAndJunctionsByCoverage)
					{
						//##################################
						//             add text
						//##################################
						double fCoverage = 0.0;
						
						//TODO adjust for size factors?
						TreeMap<String, Integer> mapJunCountsToSamples = dataSupplier.GetCountsForJunction(junction);
						
						// get median for the current condition
						int nSamples = 0;
						for(String strSample : mapSamplesToConditions.get(strSelectedCondition))
						{
							if(vcSelectedSamples.contains(strSample))
								nSamples++;
						}

						double pCoverages[] = new double[nSamples];
						
						int nIdx = 0;
						for(String strSample : mapSamplesToConditions.get(strSelectedCondition))
						{
							if(vcSelectedSamples.contains(strSample))
							{
								if(mapJunCountsToSamples == null)
									pCoverages[nIdx] = 0.0;
								else
								{
									if(mapJunCountsToSamples.containsKey(strSample))
										pCoverages[nIdx] = mapJunCountsToSamples.get(strSample);
									else
										pCoverages[nIdx] = 0.0;
								}
								nIdx++;
							}
						}
						
						fCoverage = StatUtils.percentile(pCoverages, 50.0);
						
						graph.setPaint(Color.BLUE);
						strText = String.format(Locale.ENGLISH, "%.0f", fCoverage);

						nTextWidth = fontMetrics.stringWidth(strText);
						nTextHeight = fontMetrics.getHeight();
						graph.drawString(strText, nIntronStart + (int)(nIntronLength*0.5) - (int)(nTextWidth*0.5), nY+18+nTextHeight);
						graph.setPaint(Color.BLACK);
						
						//##################################
						//             add color
						//##################################
						double fRelCoverage = fCoverage/fMaxJunCoverage;
						
						if(fRelCoverage == 0 && !bIsHighlighted)
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
					if(!bIsHighlighted && selectedResult != null)
					{
						graph.setPaint(Color.WHITE);
					}
					
					// fill the junction
					graph.fillRect(nIntronStart, nY+7, nIntronLength, 6);
					
					// intron border
					graph.setColor(Color.BLACK);
					graph.drawRect(nIntronStart, nY+7, nIntronLength, 6);
					
					area = new Area();
//					area.setId("intron_" + nPreviousExonGrpID + "_" + nExonGrpID + " " + strIsoform + " " + previousExon.getExonID() + " " + exon.getExonID());
					area.setId("intron_" + dataSupplier.GetReferenceName() + "_" + previousExon.getGenomicStop() + "_" + exon.getGenomicStart());
					area.setShape("rectangle");
					area.setCoords(nIntronStart + ", " + (nY+7) + ", " + (nIntronStart+nIntronLength) + ", " + (nY+13));
										
					String strToolTip = String.format(Locale.ENGLISH, "%s\nintron ex%d - ex%d\n%s: %,d - %,d", strIsoform, nPreviousExonGrpID, nExonGrpID, dataSupplier.GetReferenceName(), previousExon.getGenomicStop()+1, exon.getGenomicStart()-1);
					area.setTooltiptext(strToolTip);
					area.setParent(imgMapIsoforms);
				}
				
				nPreviousExonGrpID 	= nExonGrpID;
				nPreviousExonEnd 	= nX + nWidth;
				previousExon		= exon;
			}

			//#######################################
			//    add mmseq estimate if available
			//#######################################
			if(mapIsoformExpressionValues != null)
			{				
				TreeMap<String, Double> mapIsoformExpressionToSamples = null;
				if(mapIsoformExpressionValues.containsKey(strIsoform))
				{
					mapIsoformExpressionToSamples = mapIsoformExpressionValues.get(strIsoform);
				}
				else
				{
					System.out.println("failed to get mmseq results for isoform: " + strIsoform);
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
						fontNumbers = new Font("Lucida Sans", Font.PLAIN, 11);
						graph.setFont(fontNumbers);
					}
					
					nTextWidth = fontMetrics.stringWidth(strText);
					nTextHeight = fontMetrics.getHeight();
					graph.drawString(strText, nPreviousExonEnd+10, nY + 9 + (int)(nTextHeight*0.5));
					
					if(!bHasLargerPercentage)
					{
						fontNumbers = new Font("Lucida Sans", Font.BOLD, 11);  
						graph.setFont(fontNumbers);
					}
				}
			}
			
			nY += nVerticalSpaceBetweenIsoforms;
		}
				
		imgMapIsoforms.setWidth(offsets.m_nTotalWidth+"px");
		imgMapIsoforms.setHeight(nMaxHeight+"px");
		imgMapIsoforms.setContent(img);

		Vlayout layoutIsoformPlot = m_app.GetIsoformPlotRegion();
//		layoutIsoformPlot.setHeight("100%");
		imgMapIsoforms.setParent(layoutIsoformPlot);
		
		Area areaBackground = new Area();
		areaBackground.setCoords("0, 0, " + offsets.m_nTotalWidth + ", " + nMaxHeight);
		areaBackground.setId("background_isoforms");
		areaBackground.setParent(imgMapIsoforms);
		
		imgMapIsoforms.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				MouseEvent evnt = (MouseEvent) event;
				evnt.stopPropagation();
				
				DataSupplier dataSupplier	= m_app.GetDataSupplier();

				if(evnt.getArea().startsWith("highlighted_exon_group"))
				{
					String pSplit[] = evnt.getArea().split("\\s");
					
					int nGroupStart  = Integer.parseInt(pSplit[2]);
					int nGroupEnd	 = Integer.parseInt(pSplit[3]);					
					
					Popup popup = new Popup();
					popup.setParent(m_app);
					popup.open(m_app, "after_pointer");
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
					
					AnalysisResultHandler resultHandler = m_app.GetResultHandler();
					
					TreeSet<AnalysisResult> vcResults = resultHandler.GetASResultsForRegion(dataSupplier.GetReferenceName(), nGroupStart, nGroupEnd, SplicingWebApp.AS_TYPE_EXON_SKIPPING);
					
					for(AnalysisResult res : vcResults)
					{
						Row row = new Row();
						row.setParent(rows);
						
						Label label = new Label(res.GetConditionA());
						label.setParent(row);
						
						label = new Label(res.GetConditionB());
						label.setParent(row);
						
						label = new Label(String.format(Locale.ENGLISH, "%.3f", res.GetAbsoluteChangeA()));
						label.setParent(row);
						
						label = new Label(String.format(Locale.ENGLISH, "%.3f", res.GetRelativeChangeA()));
						label.setParent(row);
						
						label = new Label(String.format(Locale.ENGLISH, "%.3e", res.GetPValueA()));
						label.setParent(row);
					}
				}
				else if(evnt.getArea().startsWith("remove"))
				{
					String strIsoform = evnt.getArea().split("_")[1];	
					m_app.RemoveIsoformFromSelection(strIsoform);					
				}
				else if(evnt.getArea().startsWith("exonic_part"))
				{
					m_app.m_plotFactoryGTEX.ShowGTEXDataForExon(evnt.getArea());
					return;
				}
				else if(evnt.getArea().startsWith("intron"))
				{
					String pSplit[] = evnt.getArea().split("_");
					int nJunStart = Integer.parseInt(pSplit[2]);
					int nJunEnd = Integer.parseInt(pSplit[3]);
					
					DrawJunctionBoxplot(nJunStart, nJunEnd);
				}
			}
		});
	}
	
	/**
	 *    Draws a heatmap of the junction coverage.
	 *    Actually, two heatmaps are drawn. The first shows all junctions that pass the minimum junction coverage criteria and
	 *    colors them relative to the highest covered junction.
	 *    The second heatmap shows the same junctions but colors them relative to the highest value in each row, and uses an anova
	 *    to identify junctions that show a high variance. Only p-values that pass BH-multiple testing correction are shown.
	 */
	public void DrawJunctionHeatmap() throws IOException
	{
		ProjectModel projectModel				= m_app.GetProjectModel();
		DataSupplier dataSupplier 				= m_app.GetDataSupplier();
		String strSelectedConditionType 		= m_app.GetSelectedConditionType();
		int nMinJunctionReads					= m_app.GetMinimumJunctionReads();
		
		if(dataSupplier.GetGene() == null)
		{
			Messagebox.show("No gene selected!");
			return;
		}
		
		//###################################
		//       prepare popup window
		//###################################
		Window windowPopup = new Window();
		windowPopup.setParent(m_app);
		windowPopup.setWidth((m_nClientWindowWidth*0.5) +"px");
		windowPopup.setHeight((m_nClientWindowHeight*0.5) +"px");
		windowPopup.setTitle("Junction Heatmap");
		windowPopup.setSizable(true);
		windowPopup.setClosable(true);
		windowPopup.setMaximizable(true);
		windowPopup.setBorder(true);
		windowPopup.setPosition("center,center");
		windowPopup.setVisible(true);
		windowPopup.doPopup();
		windowPopup.setTopmost();

		//####################################################
		//           add options to options box
		//####################################################	
		Hlayout layout = new Hlayout();
		
		layout.setHflex("100%");
		layout.setHeight("100%");
		layout.setStyle("overflow:auto;");

		TreeMap<String, Double> mapSizeFactors	= projectModel.GetSizeFactors();
		TreeMap<String, String> mapBigWigFiles 	= projectModel.GetBigWigFilesForSamples();
		
		if(mapBigWigFiles == null)
		{
			System.out.println("ERROR: no valid bigwig or bam files detected");
			return;
		}
		
		// get samples per Condition
		TreeMap<String, TreeSet<String>> mapSamplesToConditions =  projectModel.GetSamplesPerCondition(strSelectedConditionType);

		// get reference name
		String strRef = dataSupplier.GetReferenceName();
		
		//##########################################################
		//                 get junction counts
		//##########################################################
		// get maximum junction count value
		int nSamples = projectModel.GetSamples().size();
		double fMaxCount = 0.0;
		for(CountElement jun : dataSupplier.GetJunctions())
		{
			TreeMap<String, Integer> mapCountsToSamples = dataSupplier.GetCountsForJunction(jun);
			
			for(String strSample : mapCountsToSamples.keySet())
			{
				double fValue = mapCountsToSamples.get(strSample) / mapSizeFactors.get(strSample);
				fMaxCount = Math.max(fMaxCount, fValue);
			}
		}
		
		//##########################################################
		//          get junctions with sufficient coverage
		//##########################################################
		Vector<CountElement> vcValidJunctions = new Vector<CountElement>();

		for(CountElement jun : dataSupplier.GetJunctions())
		{
			double fRowMaxCount = 0;
			int nValidSamples 	= 0;
						
			TreeMap<String, Integer> mapCountsToSamples = dataSupplier.GetCountsForJunction(jun);

			for(String strSample : mapCountsToSamples.keySet())
			{
				double fValue = mapCountsToSamples.get(strSample) / mapSizeFactors.get(strSample);
				fRowMaxCount = Math.max(fRowMaxCount, fValue);
				
				if(fValue >= nMinJunctionReads)
					nValidSamples++;
			}
			
			if(fRowMaxCount > 0 && nValidSamples >= 3)
			{
				vcValidJunctions.add(jun);
			}
		}
		
		//######################################################
		//                 prepare graph
		//######################################################
		int nNameOffset = 180; // subtract 100 px for the isoform names
		int nMargin = 20;
		int nSymbolWidth = 20+60;
		int nMaxWidth	= nMargin+30 + nSamples*20 + nNameOffset + nSymbolWidth;
		int nMaxHeight  = (vcValidJunctions.size()*40) + 50;
		
		// prepare graph
		BufferedImage img = new BufferedImage(nMaxWidth, nMaxHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D graph = img.createGraphics();
		
		if(m_nClientWindowWidth*0.5 > nMaxWidth)
			windowPopup.setWidth((nMaxWidth+30)+"px");
		
		if(m_nClientWindowHeight*0.5 > nMaxHeight)
			windowPopup.setHeight((nMaxHeight+60)+"px");
		
		// set font size
		Font fontNormalText = new Font("Lucida Sans", Font.PLAIN, 11);
		graph.setFont(fontNormalText);
		
		FontMetrics fontMetrics = graph.getFontMetrics(fontNormalText);
		
		// fill background
		graph.setColor(Color.white);
		graph.fillRect(0, 0, nMaxWidth, nMaxHeight);
		
		// reserve some space for the isoform name
		int nOffsetX = nMargin + nNameOffset;
		
		//###################################################
		//                   draw heatmap
		//###################################################
		int nY = 20;
		
		for(CountElement jun : vcValidJunctions)
		{
			TreeMap<String, Integer> mapCountsToSamples = dataSupplier.GetCountsForJunction(jun);
			
			double fRowMaxCount = 0;
			for(String strSample : mapCountsToSamples.keySet())
			{
				double fValue = mapCountsToSamples.get(strSample) / mapSizeFactors.get(strSample);
				fRowMaxCount = Math.max(fRowMaxCount, fValue);
			}

			int nStart = jun.GetStart()+1;
			int nEnd   = jun.GetEnd()-1;
			
			String strText = strRef + ":" + nStart + "-" + nEnd;
			graph.setColor(Color.BLACK);
 			graph.drawString(strText, nMargin, nY+fontMetrics.getHeight());
			
			int nX = nOffsetX;
			
			for(String strCondition : mapSamplesToConditions.keySet())
			{
				graph.setColor(m_mapColorsToConditions.get(strCondition));
	 			graph.drawString(strCondition, nX, 10);
				
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					// determine color
					double fValue = 0.0;
					
					if(mapCountsToSamples.containsKey(strSample))
					{
						fValue = mapCountsToSamples.get(strSample) / mapSizeFactors.get(strSample);
					}
					
					Shape rect = new Rectangle2D.Double(nX, nY, 20, 20);
					graph.setColor(new Color((int)(fValue/fMaxCount*255.0), 0, 0));
					graph.fill(rect);
					
					nX += 20;
				}
				
				// gap for each condition
				nX += 2;
			}
			
			nY += 20;
		}
		
		
		int nTmpY = nY+20;
		int nTmpX = 0;
		Vector<Double> vcPValues = new Vector<Double>(); // store anova p-values
		
		nY += 20;
		for(CountElement jun : vcValidJunctions)
		{
			TreeMap<String, Integer> mapCountsToSamples = dataSupplier.GetCountsForJunction(jun);
			
			double fRowMaxCount = 0;
			for(String strSample : mapCountsToSamples.keySet())
			{
				double fValue = mapCountsToSamples.get(strSample) / mapSizeFactors.get(strSample);
				fRowMaxCount = Math.max(fRowMaxCount, fValue);
			}
			
			int nStart = jun.GetStart()+1;
			int nEnd   = jun.GetEnd()-1;
			
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
					double fValue = 0.0;
					
					if(mapCountsToSamples.containsKey(strSample))
					{
						fValue = mapCountsToSamples.get(strSample) / mapSizeFactors.get(strSample);
					}
					
					pCounts[nIdx] = fValue;
					float fRes = (float)(fValue / fRowMaxCount);
					
					// determine color
					Shape rect = new Rectangle2D.Double(nX, nY, 20, 20);
					graph.setColor(new Color((int)(fRes*255.0), 0, 0));
					
					graph.fill(rect);
					
					nX += 20;
					
					nIdx++;
				}
				
				// gap for each condition
				nX += 2;
				
				vcCounts.add(pCounts);
			}
			
			nTmpX = nX;
			
			// calculate anova	
			OneWayAnova anova = new OneWayAnova();
			double fPValue = anova.anovaPValue(vcCounts);
			
			// store p-value
			vcPValues.add(fPValue);

			nY += 20;
		}
			
		//###################################
		//    multiple testing correction
		//###################################
		int nTests = vcPValues.size();
		
		// create backup of p-vlaues to keep the order
		Vector<Double> vcOrgPValues = new Vector<Double>();
		vcOrgPValues.setSize(vcPValues.size());
		Collections.copy(vcOrgPValues, vcPValues);
		
		// sort p-values ascendingly
		Collections.sort(vcPValues);
		
		// find the last significant hit with highest p-value
		double fThreshold = 0.0f;
		for(int i=0; i<nTests; i++)
		{
			double fCriticalValue = ((double)i / (double)nTests) * 0.05;
			double fPValue = vcPValues.get(i);
	
			if(fPValue < fCriticalValue)
				fThreshold = fPValue;
		}

		// indicate significant p-values
		nY = nTmpY;
		for(double fPValue : vcOrgPValues)
		{
			if(fPValue > fThreshold)
			{
				nY += 20;
				continue;
			}
			
			String strPath = Executions.getCurrent().getDesktop().getWebApp().getRealPath("/");
			BufferedImage imgExclamation = ImageIO.read(new File(strPath + "/img/exclamation_sign.png"));
			BufferedImage resizedImg = new BufferedImage(20, 18, BufferedImage.TYPE_INT_ARGB);
			
			Graphics2D g = resizedImg.createGraphics();
			g.drawImage(imgExclamation, 0, 0, 20, 18, null);
			g.dispose();				
			
			graph.drawImage(resizedImg, nTmpX+5, nY, null);
			
			graph.setColor(Color.black);
			graph.drawString(String.format(Locale.ENGLISH, "p.value %.3f", fPValue), nTmpX+30, nY+13);
			
			nY += 20;
		}
		
		layout.setParent(windowPopup);
		
		Imagemap imgMap = new Imagemap();
		imgMap.setWidth(nMaxWidth+"px");
		imgMap.setHeight((nMaxHeight+60)+"px");
		imgMap.setContent(img);
		imgMap.setParent(layout);
	}

	/**
	 *    Draws boxplots for a single splice junction
	**/
	public void DrawJunctionBoxplot(int nJunStart, int nJunEnd) throws IOException
	{
		ProjectModel projectModel				= m_app.GetProjectModel();
		DataSupplier dataSupplier 				= m_app.GetDataSupplier();
		String strSelectedConditionType 		= m_app.GetSelectedConditionType();
		
		// get size factors
		TreeMap<String, Double> mapSizeFactorsToSamples = projectModel.GetSizeFactors();
		
		// get samples per Condition
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = projectModel.GetSamplesPerCondition(strSelectedConditionType);
		
		if(dataSupplier.GetGene() == null)
		{
			Messagebox.show("No gene selected!");
			return;
		}
		
		//##########################################################
		//                 get junction counts
		//##########################################################
		
		// get number of samples
		int nSamples = projectModel.GetSamples().size();
		
		// get adjusted junction counts for each sample, grouped by condition
		TreeMap<String, double[]> mapCountsToConditions = new TreeMap<String, double[]>();
		Vector<String> vcDataLabels = new Vector<String>();
		
		double fMaxCount = 0.0;
		for(CountElement jun : dataSupplier.GetJunctions())
		{
			if(jun.m_nStart == nJunStart && jun.m_nEnd == nJunEnd)
			{
				TreeMap<String, Integer> mapCountsToSamples = dataSupplier.GetCountsForJunction(jun);
				
				for(String strCondition : mapSamplesToConditions.keySet())
				{
					Vector<Double> vcCounts = new Vector<Double>();
					
					for(String strSample : mapSamplesToConditions.get(strCondition))
					{
						double fCount = 0.0;
						if(mapCountsToSamples.containsKey(strSample))
							fCount = mapCountsToSamples.get(strSample) / mapSizeFactorsToSamples.get(strSample);
						vcCounts.add(fCount);
						fMaxCount = Math.max(fMaxCount, fCount);
						
						vcDataLabels.add(strSample);
					}
					
					double pCounts[] = new double[vcCounts.size()];
					for(int i=0; i<vcCounts.size(); i++)
						pCounts[i] = vcCounts.get(i);
					
					mapCountsToConditions.put(strCondition, pCounts);
				}
				break;
			}
		}

		//###################################
		//       prepare popup window
		//###################################
		Window windowPopup = new Window();
		windowPopup.setParent(m_app);
		windowPopup.setWidth("600px");
		windowPopup.setHeight("605px");
		windowPopup.setTitle("Junction Boxplots");
		windowPopup.setSizable(true);
		windowPopup.setClosable(true);
		windowPopup.setMaximizable(true);
		windowPopup.setBorder(true);
		windowPopup.setPosition("center,center");
		windowPopup.setVisible(true);
		windowPopup.doPopup();
		windowPopup.setTopmost();
		
		//###################################################
		//                   draw heatmap
		//###################################################
		
		BufferedImage img = BoxPlot(mapCountsToConditions, nJunStart+"-"+nJunEnd, null, 65, 5, 405, "adj. count");
		
		windowPopup.setWidth( Math.min(m_nClientWindowWidth *0.75, (img.getWidth()+10))+"px");
		windowPopup.setHeight(Math.min(m_nClientWindowHeight*0.75, (img.getHeight()+55))+"px");
					
		Imagemap imgMap = new Imagemap();
		imgMap.setWidth(img.getWidth()+"px");
		imgMap.setHeight(img.getHeight()+"px");
		imgMap.setContent(img);
		imgMap.setParent(windowPopup);
	}
	
	/**
	 *    Generates a box plot for the given count arrays per condition. The plot will be located at nX/nX with height nHeight on the graph.
	 *    Returns the width [in pixel] of the plot.
	 */
	public BufferedImage BoxPlot(TreeMap<String, double[]> mapCountsToConditions, String strTitle, Graphics2D graph, int nX, int nY, int nHeight, String strAxisTitleY)
	{
		// settings for the box plot
		int nBarWidth 	= 80;
		int nSpacer   	= 20;
		int nTotalBarDistance = nBarWidth+nSpacer;
		
		int nWidth = mapCountsToConditions.keySet().size() * nTotalBarDistance;

		// prepare a new graph if required
		BufferedImage img = null;
		if(graph == null)
		{
			img = new BufferedImage(nWidth+nX, nHeight+nY+45, BufferedImage.TYPE_INT_RGB); // +30 for the x-axis text
			graph = img.createGraphics();
			
			// fill background
			graph.setColor(Color.white);
			graph.fillRect(0, 0, nWidth+nX, nHeight+nY+45);
		}
		
		// process data
		double fMaxValue = 1.0;
		for(String strCondition : mapCountsToConditions.keySet())
		{
			for(double fVal : mapCountsToConditions.get(strCondition))
			{
				fMaxValue = Math.max(fVal, fMaxValue);
			}
		}
		
		// set font
		Font fontNormalText = new Font("Lucida Sans", Font.BOLD, 14);
		FontMetrics fontMetrics = graph.getFontMetrics(fontNormalText);
		
		graph.setFont(fontNormalText);
		graph.setColor(Color.BLACK);
		
		// add title
		int nStringHeight 	= fontMetrics.getHeight();
		int nTitleYPos  	= nY-10;
		
		int nStringWidth  = fontMetrics.stringWidth(strTitle);
		graph.drawString(strTitle, nX+(int)(nWidth*0.5 - nStringWidth*0.5), nTitleYPos);

		//################
		//    add axes
		//################
		// add X Axis
		graph.drawLine(nX, nY+nHeight, nX+nWidth,  nY+nHeight);
		
		// set font
		fontNormalText = new Font("Lucida Sans", Font.PLAIN, 11);
		fontMetrics = graph.getFontMetrics(fontNormalText);
		
		graph.setFont(fontNormalText);
		graph.setColor(Color.BLACK);
		
		// add ticks to x-axis
		int j = 0;
		boolean bPreviousTextHadOffset = false;
		for(String strCondition : mapCountsToConditions.keySet())
		{
			int nXPos = nX + j*nTotalBarDistance + (int)(nTotalBarDistance*0.5);
			graph.drawLine(nXPos, nY+nHeight, nXPos, nY+nHeight+3);
			
			nStringWidth  = fontMetrics.stringWidth(strCondition);
			
			if(nStringWidth > (nBarWidth+nSpacer*0.5))
			{
				if(bPreviousTextHadOffset)
					bPreviousTextHadOffset = false;
				else
					graph.drawString(strCondition, (int)(nXPos-nStringWidth*0.5), nY+nHeight+nStringHeight*2+3);				
			}
			else
				graph.drawString(strCondition, (int)(nXPos-nStringWidth*0.5), nY+nHeight+nStringHeight+3);
			j++;
		}
		
//		AddXAxis();
		fMaxValue = AddYAxis(nX, nY, nHeight, fMaxValue, strAxisTitleY, graph);
		
		// add boxplots
		j = 0;
		for(String strCondition : mapCountsToConditions.keySet())
		{
			double pCounts[] = null;
			
			pCounts = mapCountsToConditions.get(strCondition);
			
			double fPercentile90 = StatUtils.percentile(pCounts, 90.0);
			double fPercentile75 = StatUtils.percentile(pCounts, 75.0);
			double fPercentile50 = StatUtils.percentile(pCounts, 50.0);
			double fPercentile25 = StatUtils.percentile(pCounts, 25.0);
			double fPercentile10 = StatUtils.percentile(pCounts, 10.0);

			int nXPos = nX + j*nTotalBarDistance + (int)(nTotalBarDistance*0.5 - nBarWidth*0.5);

			Color clr1 = null;
			if(m_mapColorsToConditions.containsKey(strCondition))
			{
				clr1 = m_mapColorsToConditions.get(strCondition);
			}
			else
				clr1 = m_pColors[j];
			
			Color clr2 = clr1.darker();
			Point ptStart = new Point(nXPos, 0);
			Point ptEnd   = new Point(nXPos + nBarWidth, 0);
			float pDist[] = {0.0f,0.2f, 0.8f, 1.0f};
			Color pClrs[]  = {clr2, clr1, clr1, clr2};
			LinearGradientPaint clr = new LinearGradientPaint(ptStart, ptEnd, pDist, pClrs);
			
			// draw upper box
			int nTop 	= nY + (int) Math.ceil(nHeight - (fPercentile75 / fMaxValue * nHeight));
			int nBottom = nY + (int) Math.ceil(nHeight - (fPercentile50 / fMaxValue * nHeight));
			graph.setPaint(clr);
			graph.fillRect(nXPos, nTop, nBarWidth, nBottom-nTop);
			graph.setColor(Color.black);
			graph.drawRect(nXPos, nTop, nBarWidth, nBottom-nTop);
			
			// draw lower box
			nTop 	= nY + (int) Math.ceil(nHeight - (fPercentile50 / fMaxValue * nHeight));
			nBottom = nY + (int) Math.ceil(nHeight - (fPercentile25 / fMaxValue * nHeight));
			graph.setPaint(clr);
			graph.fillRect(nXPos, nTop, nBarWidth, nBottom-nTop);
			graph.setColor(Color.black);
			graph.drawRect(nXPos, nTop, nBarWidth, nBottom-nTop);
			
			// plot upper whiskers
			nTop 	= nY + (int) Math.ceil(nHeight - (fPercentile90 / fMaxValue * nHeight));
			nBottom = nY + (int) Math.ceil(nHeight - (fPercentile75 / fMaxValue * nHeight));
			graph.drawLine(nXPos + (int)(nBarWidth*0.5), nTop, nXPos + (int)(nBarWidth*0.5), nBottom);
			graph.drawLine(nXPos + (int)(nBarWidth*0.25), nTop, nXPos + (int)(nBarWidth*0.75), nTop);
			
			// plot lower whiskers
			nTop 	= nY + (int) Math.ceil(nHeight - (fPercentile25 / fMaxValue * nHeight));
			nBottom = nY + (int) Math.ceil(nHeight - (fPercentile10 / fMaxValue * nHeight));
			graph.drawLine(nXPos + (int)(nBarWidth*0.5), nTop, nXPos + (int)(nBarWidth*0.5), nBottom);
			graph.drawLine(nXPos + (int)(nBarWidth*0.25), nBottom, nXPos + (int)(nBarWidth*0.75), nBottom);
			
			j++;
		}
		
		return img;
	}
	
	/**
	 *    Generates a 2D (XY) scatter plot using the counts specified for two junctions given as count arrays per condition. 
	 */
	public Vector<Area> Plot2D(TreeMap<String, double[]> mapCountsToConditionsX, TreeMap<String, double[]> mapCountsToConditionsY, Vector<String> vcDataLabels, String strAxisLabelX, String strAxisLabelY, String strTitle, Graphics2D graph, int nX, int nY, int nWidth, int nHeight, boolean bAxesInPercent)
	{
		Vector<Area> vcAreas = new Vector<Area>();
		
		// get maximum value for both axes
		double fMaxValueX = 1.0;
		for(String strCondition : mapCountsToConditionsX.keySet())
		{
			for(double fVal : mapCountsToConditionsX.get(strCondition))
			{
				fMaxValueX = Math.max(fVal, fMaxValueX);
			}
		}
		
		double fMaxValueY = 1.0;
		for(String strCondition : mapCountsToConditionsY.keySet())
		{
			for(double fVal : mapCountsToConditionsY.get(strCondition))
			{
				fMaxValueY = Math.max(fVal, fMaxValueY);
			}
		}
		
		//########################
		//        add title
		//########################
		// set font
		Font fontNormalText = new Font("Lucida Sans", Font.PLAIN, 14);
		FontMetrics fontMetrics = graph.getFontMetrics(fontNormalText);
		
		graph.setFont(fontNormalText);
		graph.setColor(Color.BLACK);
		int nTitleYPos  	= nY-10;
		
		int nStringWidth  = fontMetrics.stringWidth(strTitle);
		graph.drawString(strTitle, nX+(int)(nWidth*0.5 - nStringWidth*0.5), nTitleYPos);
				
		//########################
		//        add axes
		//########################		
		fMaxValueX = AddXAxis(nX, nY+nHeight, nWidth, fMaxValueX, strAxisLabelX, graph);
		fMaxValueY = AddYAxis(nX, nY, 		 nHeight, fMaxValueY, strAxisLabelY, graph);

		//########################
		//   add data points
		//########################
		for(String strCondition : mapCountsToConditionsX.keySet())
		{
			double pValuesX[] = mapCountsToConditionsX.get(strCondition);
			double pValuesY[] = mapCountsToConditionsY.get(strCondition);
			
			graph.setColor(m_mapColorsToConditions.get(strCondition));
			
			for(int i=0; i<pValuesX.length; i++)
			{
				double fValX = pValuesX[i];
				double fValY = pValuesY[i];

				double fX = nX + (fValX / fMaxValueX) * nWidth;
				double fY = nY + nHeight - (fValY / fMaxValueY) * nHeight;
				
				Ellipse2D.Double dataPoint = new Ellipse2D.Double();
				dataPoint.x = fX - 4.0;
				dataPoint.y = fY - 4.0;
				dataPoint.width = 8.0;
				dataPoint.height = 8.0;
				
				// add clickable area
				Area area = new Area();
				area.setShape("circle");
				area.setCoords(dataPoint.x + ", " + dataPoint.y + ", " + 8);				
				area.setTooltiptext("Sample: " + vcDataLabels.get(i));
				vcAreas.add(area);
									
				graph.fill(dataPoint);
			}
		}
		
		return vcAreas;
	}
	
	/**
	 *    Opens a popup window that zooms in on an alternative splicing event. The surrounding exons and the alternative exon are shown
	 *    accompanied by a coverage plot above the exons (and some adjacent intron bases).
	 *    ATTENTION: because the viewer was never intended to show zoomed in views of splicing events, this function uses some workarounds,
	 *    is thus rather complicated and requires some simplifications.
	 */
	public void DrawSpliceJunctionCloseup(AnalysisResult result)
	{
		ProjectModel projectModel				= m_app.GetProjectModel();
		DataSupplier dataSupplier 				= m_app.GetDataSupplier();
		String strSelectedConditionType 		= m_app.GetSelectedConditionType();
		TreeSet<String> vcValidIsoforms 		= m_app.GetValidIsoforms();
		TreeSet<String> vcSelectedSamples		= m_app.GetSelectedSamples();
//		int nMinJunctionReads					= m_app.GetMinimumJunctionReads();
		Exon pExons[]							= dataSupplier.GetExons();		
		Gene gene 								= dataSupplier.GetGene();
		TreeSet<Exon> vcMiddleExons				= dataSupplier.GetMiddleExons();		
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = projectModel.GetSamplesPerCondition(strSelectedConditionType);
		
		class Region implements Comparable<Region>
		{
			double m_fStart;
			double m_fEnd;
			double m_fLength;
			double m_fOriginalStart;
			double m_fOriginalLength;
			double m_fOffsetY;
			boolean m_bIsLeftTerminalExon;
			boolean m_bIsRightTerminalExon;
			boolean m_bIsIncomplete;
			boolean m_bIsNovel;
			boolean m_bIsTarget;
			
			public Region(double fStart, double fEnd, double fOffsetY, boolean bIsLeftTerminalExon, boolean bIsRightTerminalExon, double fOriginalStart, double fOriginalLength)
			{
				m_fStart 		= fStart;
				m_fEnd	 		= fEnd;
				m_fLength		= fEnd - fStart + 1;
				m_fOffsetY 		= fOffsetY;
				m_bIsLeftTerminalExon  = bIsLeftTerminalExon;
				m_bIsRightTerminalExon = bIsRightTerminalExon;
				m_fOriginalLength = fOriginalLength;
				m_fOriginalStart  = fOriginalStart;
				m_bIsNovel		= false;
				m_bIsTarget		= false;
				
				if(m_fLength < fOriginalLength || fStart > fOriginalStart)
					m_bIsIncomplete = true;
				else
					m_bIsIncomplete = false;
			}
			
			public Region(double fStart, double fEnd, double fOffsetY, boolean bIsLeftTerminalExon, boolean bIsRightTerminalExon, boolean bIsIncomplete)
			{
				m_fStart 		= fStart;
				m_fEnd	 		= fEnd;
				m_fLength		= fEnd - fStart + 1;
				m_fOffsetY 		= fOffsetY;
				m_bIsLeftTerminalExon  = bIsLeftTerminalExon;
				m_bIsRightTerminalExon = bIsRightTerminalExon;
				m_bIsIncomplete	= bIsIncomplete;
				m_bIsNovel		= false;
				m_bIsTarget		= false;
			}
			
			public void SetIsNovel(boolean bIsNovel)
			{
				m_bIsNovel = bIsNovel;
			}
			
			public void SetIsTarget(boolean bIsTarget)
			{
				m_bIsTarget = bIsTarget;
			}

			@Override
			public int compareTo(Region other)
			{
				if(m_fStart < other.m_fStart)
					return -1;
				
				if(m_fStart > other.m_fStart)
					return 1;
				
				if(m_fEnd < other.m_fEnd)
					return -1;
				
				if(m_fEnd > other.m_fEnd)
					return 1;
				
				if(m_fOffsetY < other.m_fOffsetY)
					return -1;
				
				if(m_fOffsetY > other.m_fOffsetY)
					return 1;
				
				return 0;
			}
			
			public boolean equals(Region other)
			{
				if(m_fStart == other.m_fStart && m_fEnd == other.m_fEnd && m_fOffsetY == other.m_fOffsetY)
					return true;
				
				return false;
			}
			
			public boolean overlaps(Region other)
			{
				if(m_fStart <= other.m_fEnd && m_fEnd >= other.m_fStart)
					return true;
				
				return false;
			}
			
			public void extendByExons(TreeSet<Region> vcModifiedExons)
			{
				for(Region ex : vcModifiedExons)
				{
					if(ex.overlaps(this))
					{
						// extend the start of 'this' exons
						if(m_fStart > ex.m_fStart)
							m_fStart = Math.max(ex.m_fStart, m_fOriginalStart);
						
						// extend the end of 'this' exons
						if(m_fEnd < ex.m_fEnd)
							m_fEnd = Math.min(ex.m_fEnd, m_fOriginalStart + m_fOriginalLength-1);
						
						// extend the start of 'the other' exons
						if(ex.m_fStart > m_fStart)
							ex.m_fStart = Math.max(m_fStart, ex.m_fOriginalStart);
						
						// extend the end of 'the other' exons
						if(ex.m_fEnd < m_fEnd)
							ex.m_fEnd = Math.min(m_fEnd, ex.m_fOriginalStart + ex.m_fOriginalLength-1);
					}
				}
			}
			
			public String toString()
			{
				return m_fStart + " - " + m_fEnd + "\n - " + m_fOffsetY + "\n - leftTermExon: " + m_bIsLeftTerminalExon + "\n - rightTermExon: " + m_bIsRightTerminalExon + "\n - incompl.:" + m_bIsIncomplete + "\n - novel: " + m_bIsNovel + "\n - target: " + m_bIsTarget + "\n - org. start: " + m_fOriginalStart  + "\n - org. length: " + m_fOriginalLength + "\n";
			}
		};
		
		//###################################
		//       prepare popup window
		//###################################
		int nMaxWidth  = (int)(m_nClientWindowWidth*0.7);		
		int nWindowHeight = (int)(m_nClientWindowHeight*0.7);
		
		Window windowPopup = new Window();
		windowPopup.setParent(m_app);
		windowPopup.setWidth(nMaxWidth +"px");
		windowPopup.setHeight((nWindowHeight+20) +"px");
		windowPopup.setTitle("Detailed View");
		windowPopup.setSizable(true);
		windowPopup.setClosable(true);
		windowPopup.setMaximizable(true);
		windowPopup.setBorder(true);
		windowPopup.setPosition("center,center");
		windowPopup.setVisible(true);
		windowPopup.doPopup();
		windowPopup.setTopmost();
		
		Hlayout layout = new Hlayout();
		layout.setWidth("100%");
		layout.setHeight("100%");
		layout.setStyle("overflow:auto;");
		layout.setParent(windowPopup);
		
		int nNameOffset = 180; // subtract 100 px for the isoform names
		int nMargin = 10;
		
		// settings for the plots
		int nInterPlotDistance = 120;
		int nPlotHeight = 300;
		int nBarWidth 	= 80;
		int nSpacer   	= 20;
		int nTotalBarDistance = nBarWidth+nSpacer;
		
		TreeSet<String> vcValidConditions = new TreeSet<String>();
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			// determine number of valid samples
			int nValidSamples = 0;
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				if(vcSelectedSamples.contains(strSample))
					nValidSamples++;
			}
			
			if(nValidSamples == 0)
				continue;
			
			vcValidConditions.add(strCondition);
		}
		
		int nConditions	  = vcValidConditions.size();		
		int nBoxPlotWidth = nConditions * nTotalBarDistance;
		int n2DPlotWidth  = 300;
		int nCoveragePlotHeight = 200;
		
		int nMaxHeight = nPlotHeight + nCoveragePlotHeight + 260;

		// prepare graph
		int nMaxPlotWidth = Math.max(nMaxWidth, (nBoxPlotWidth+nInterPlotDistance)*3 + (n2DPlotWidth+nInterPlotDistance) + nMargin + nNameOffset);
		
		if(!result.HasPSIScore())
			return;
		
		// get location of alternatively spliced exon groups
		int nStartA = result.GetInclusionJunctionStart();
		int nEndA	= result.GetInclusionJunctionEnd();
		int nStartB = result.GetExclusionJunctionStart();
		int nEndB	= result.GetExclusionJunctionEnd();
		
		// get all exon groups involved in this splicing event
		Exon exLeftA  = null;
		Exon exLeftB  = null;
		Exon exRightA = null;
		Exon exRightB = null;
		
		// get exons linked to the inclusion and exclusion junctions
		// keep the largest middle exon, or if none is available, keep the largest terminal exon
		for(int i=0; i<2; i++)
		{
			if(exLeftA != null && exRightA != null && exLeftB != null && exRightB != null)
				break;
			
			for(Exon ex : pExons)
			{			
				if(ex.getGenomicStop() == nStartA && (exLeftA == null || (vcMiddleExons.contains(ex) && i==0)))
				{
					if(exLeftA == null || (exLeftA != null && exLeftA.getGenomicLength() < ex.getGenomicLength()))
						exLeftA = ex;
				}
							
				if(ex.getGenomicStart() == nEndA && (exRightA == null || (vcMiddleExons.contains(ex) && i==0)))
				{
					if(exRightA == null || (exRightA != null && exRightA.getGenomicLength() < ex.getGenomicLength()))
						exRightA = ex;
				}
				
				if(ex.getGenomicStop() == nStartB && (exLeftB == null || (vcMiddleExons.contains(ex) && i==0)))
				{
					if(exLeftB == null || (exLeftB != null && exLeftB.getGenomicLength() < ex.getGenomicLength()))
						exLeftB = ex;
				}
				
				if(ex.getGenomicStart() == nEndB && (exRightB == null || (vcMiddleExons.contains(ex) && i==0)))
				{
					if(exRightB == null || (exRightB != null && exRightB.getGenomicLength() < ex.getGenomicLength()))
						exRightB = ex;
				}
			}
		}
		
		if(exLeftA == null)
		{
			exLeftA = new Exon(nStartA-200, nStartA);
		}
		
		if(exLeftB == null)
		{
			exLeftB = new Exon(nStartB-200, nStartB);
		}
		
		if(exRightA == null)
		{
			exRightA = new Exon(nEndA, nEndA+200);
		}
		
		if(exRightB == null)
		{
			exRightB = new Exon(nEndB, nEndB+200);
		}
		
		//######################################################
		//                 prepare graph
		//######################################################
		BufferedImage img = new BufferedImage(nMaxPlotWidth, nMaxHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D graph = img.createGraphics();
		
		Font fontNumbers = new Font("Lucida Sans", Font.BOLD, 11);
		graph.setFont(fontNumbers);
		
		graph.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		// load noncontinous-mark images
		String strPath = Executions.getCurrent().getDesktop().getWebApp().getRealPath("/");
		BufferedImage imgMark = null;
		try
		{
			imgMark = ImageIO.read(new File(strPath + "/img/noncontinuous_left.png"));
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.out.println("could not open image: " + strPath + "/img/noncontinuous_left.png");
			return;
		}
		BufferedImage imgMarkLeft = new BufferedImage(42, 99, BufferedImage.TYPE_INT_ARGB);
		
		Graphics2D g = imgMarkLeft.createGraphics();
		g.drawImage(imgMark, 0, 0, 20, 30, null);
		g.dispose();
		
		try
		{
			imgMark = ImageIO.read(new File(strPath + "/img/noncontinuous_right.png"));
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.out.println("could not open image: " + strPath + "/img/noncontinuous_right.png");
			return;
		}
		BufferedImage imgMarkRight = new BufferedImage(42, 99, BufferedImage.TYPE_INT_ARGB);
		
		g = imgMarkRight.createGraphics();
		g.drawImage(imgMark, 0, 0, 20, 30, null);
		g.dispose();

		// fill background
		graph.setColor(Color.white);
		graph.fillRect(0, 0, nMaxPlotWidth, nMaxHeight);
		
		// subtract the margin from the canvas width and reserve some space for the exon names, also subtract something for the window border and scrollbars, etc.
		nMaxWidth -= (nMargin*2 + nNameOffset + 100);
		
		// reserve some space for the isoform name
		int nOffsetX = nMargin + nNameOffset;
		
		// calculate exon extension length
		int nExtensionLength 	= Integer.MAX_VALUE;
		int nMinDrawLength		= 100;
		
		if(result.GetType() == SplicingWebApp.AS_TYPE_ALT_3_PRIME_EXON_END || result.GetType() == SplicingWebApp.AS_TYPE_ALT_5_PRIME_EXON_END )
		{
			nExtensionLength = 0;
			
			if(!exLeftA.equals(exLeftB) && exLeftA.intersects(exLeftB))
			{
				nExtensionLength = Math.max(exLeftA.getGenomicStop(), exLeftB.getGenomicStop()) - Math.min(exLeftA.getGenomicStop(), exLeftB.getGenomicStop());
			}
			if(!exRightA.equals(exRightB) && exRightA.intersects(exRightB))
			{
				nExtensionLength = Math.max(exRightA.getGenomicStart(), exRightB.getGenomicStart()) - Math.min(exRightA.getGenomicStart(), exRightB.getGenomicStart());
			}
			
			// try to draw 2x the extension length from each exon, but try to draw at least 50 bp
			nMinDrawLength = Math.max(nExtensionLength * 2, nMinDrawLength);	
		}
		
		int nOffsetY = nMargin + nCoveragePlotHeight + 20;
		
		// Colors
		Color clrTargetExonTop 				= new Color(217,217,255, 255);	// BLUE
		Color clrTargetExonBottom 			= new Color(120,120,230, 255);	// BLUE
		Color clrTargetNovelExonTop 		= new Color(255,217,217, 255);	// RED
		Color clrTargetNovelExonBottom 		= new Color(230,120,120, 255);	// RED
		
		Color clrNonTargetExonTop 			= new Color(217,217,255, 100);	// BLUE
		Color clrNonTargetExonBottom 		= new Color(120,120,230, 100);	// BLUE
		Color clrNonTargetNovelExonTop 		= new Color(255,217,217, 100);	// RED
		Color clrNonTargetNovelExonBottom 	= new Color(230,120,120, 100);	// RED
		
		// draw legend
		int nStartX = nMargin;
		int nStartY = nMargin;
		int nWidth	= 10;
		GradientPaint clr = new GradientPaint(nStartX, nStartY, clrTargetExonTop, nStartX+20, nStartY+20, clrTargetExonBottom);
		graph.setPaint(clr);
		
		graph.fillRect(nStartX, nStartY, nWidth, nWidth);
		graph.setColor(Color.black);
		graph.drawRect(nStartX, nStartY, nWidth, nWidth);
		
		graph.drawString("known target", nStartX + 25, nStartY+10);

		//---
		nStartY += nWidth + 5;
		clr = new GradientPaint(nStartX, nStartY, clrTargetNovelExonTop, nStartX+nWidth, nStartY+nWidth, clrTargetNovelExonBottom);
		graph.setPaint(clr);
		
		graph.fillRect(nStartX, nStartY, nWidth, nWidth);
		graph.setColor(Color.black);
		graph.drawRect(nStartX, nStartY, nWidth, nWidth);
		
		graph.drawString("novel target", nStartX + 25, nStartY+10);
		
		//---
		nStartY += nWidth + 5;
		clr = new GradientPaint(nStartX, nStartY, clrNonTargetExonTop, nStartX+nWidth, nStartY+nWidth, clrNonTargetExonBottom);
		graph.setPaint(clr);
		
		graph.fillRect(nStartX, nStartY, nWidth, nWidth);
		graph.setColor(Color.black);
		graph.drawRect(nStartX, nStartY, nWidth, nWidth);
		
		graph.drawString("known non-target", nStartX + 25, nStartY+10);
		
		//---
		nStartY += nWidth + 5;
		clr = new GradientPaint(nStartX, nStartY, clrNonTargetNovelExonTop, nStartX+nWidth, nStartY+nWidth, clrNonTargetNovelExonBottom);
		graph.setPaint(clr);
		
		graph.fillRect(nStartX, nStartY, nWidth, nWidth);
		graph.setColor(Color.black);
		graph.drawRect(nStartX, nStartY, nWidth, nWidth);
		
		graph.drawString("novel non-target", nStartX + 25, nStartY+10);
		
		// draw isoform labels
		graph.setColor(Color.black);
		graph.drawString("inclusion path", nMargin, nOffsetY+16);
		graph.drawString("exclusion path", nMargin, nOffsetY+56);
		
		//############################################
		//    specify regions that should be drawn
		//############################################
		double fOffsetY = nOffsetY;		
		int nRegionStart = Math.min(result.GetInclusionJunction().m_nStart, result.GetExclusionJunction().m_nStart);
		int nRegionEnd   = Math.max(result.GetInclusionJunction().m_nEnd,   result.GetExclusionJunction().m_nEnd);

		TreeSet<Region> vcModifiedExons = new TreeSet<Region>();

		// get the best matching isoforms for the inclusion and exclusion junction
		for(int i=0; i<2; i++)
		{
			TreeSet<Exon> vcTargetExons = new TreeSet<Exon>();
			CountElement junTarget = result.GetInclusionJunction();
			if(i == 1)
				junTarget = result.GetExclusionJunction();

			// calculate average coverage for all junctions
			TreeSet<CountElement> vcJunctions = dataSupplier.GetJunctions();
			TreeMap<CountElement, Double> mapAvgCoverageToJunctions = new TreeMap<CountElement, Double>();
			
			for(CountElement jun : vcJunctions)
			{
				TreeMap<String, Integer> mapCountsToSamples = dataSupplier.GetCountsForJunction(jun);
				
				int nCount   = 0;
				for(String strSample : vcSelectedSamples)
				{
					if(mapCountsToSamples.containsKey(strSample))
						nCount += mapCountsToSamples.get(strSample);
				}
				
				double fAverageCoverage = nCount / vcSelectedSamples.size();
				mapAvgCoverageToJunctions.put(jun, fAverageCoverage);
			}
			
			// find all isoforms that contain the alternative spliced exon
			TreeSet<String> vcCandidateIsoforms = new TreeSet<String>();
			for(String strIsoform : vcValidIsoforms)
			{
				// find all isoforms that include the target junction
				TreeSet<CountElement> pIsoformJunctions = dataSupplier.GetJunctionsForIsoform(strIsoform);
				for(CountElement jun : pIsoformJunctions)
				{
					if(jun.GetStart() == junTarget.m_nStart && jun.GetEnd() == junTarget.m_nEnd)
					{
						vcCandidateIsoforms.add(strIsoform);
						break;
					}
				}
			}
/*			
			// now find the isoform where the coverage of the junctions fit best to the coverage of the
			// inclusion/exclusion junction
			double fAvgTargetJunctionCoverage = mapAvgCoverageToJunctions.get(junTarget);
			
			String strTargetIsoform = "";
			double fMinDiff = Double.MAX_VALUE;
			for(String strIsoform : vcCandidateIsoforms)
			{
				TreeSet<CountElement> vcIsoformJunctions = dataSupplier.GetJunctionsForIsoform(strIsoform);
				
				double fDifferenceSum = 0.0;
				double fJunCount = 0.0;
				
				for(CountElement jun : vcIsoformJunctions)
				{
					// restrict to junctions in the target region
					if(jun.m_nStart >= nRegionStart && jun.m_nEnd <= nRegionEnd && !jun.equals(junTarget))
					{
						if(mapAvgCoverageToJunctions.containsKey(jun))
							fDifferenceSum += Math.abs(mapAvgCoverageToJunctions.get(jun) - fAvgTargetJunctionCoverage);
						else
							fDifferenceSum += fAvgTargetJunctionCoverage;
						fJunCount++;
					}
				}

				double fDiff = Double.MAX_VALUE;
				if(i == 1 && fJunCount == 0.0) // allow single junction isoforms for exclusion junctions
					fDiff = 0;
				else if(fJunCount != 0.0)
					fDiff = fDifferenceSum / fJunCount;
				
				if(fDiff < fMinDiff)
				{
					strTargetIsoform = strIsoform;
					fMinDiff = fDiff;
				}
			}
*/
			// use isoform with most exons
			String strTargetIsoform = "";
			double fMaxExons = -1;
			for(String strIsoform : vcCandidateIsoforms)
			{
				TreeSet<CountElement> vcIsoformJunctions = dataSupplier.GetJunctionsForIsoform(strIsoform);
				
				double fJunCount = 0.0;
				
				for(CountElement jun : vcIsoformJunctions)
				{
					// restrict to junctions in the target region
					if(jun.m_nStart >= nRegionStart && jun.m_nEnd <= nRegionEnd && !jun.equals(junTarget))
					{
						fJunCount++;
					}
				}
				
				if(i == 1 && fJunCount == 0.0) // allow single junction isoforms for exclusion junctions
					fJunCount = 1;

				if(fJunCount > fMaxExons)
				{
					fMaxExons = fJunCount;
					strTargetIsoform = strIsoform;
				}
			}
			
			if(!strTargetIsoform.isEmpty())
			{
				for(Exon ex : gene.getSortedExonsForGeneProduct(strTargetIsoform))
				{
					// add all exons within the target region
					if(ex.getCodingStart() >= nRegionStart && ex.getCodingStop() <= nRegionEnd)
						vcTargetExons.add(ex);
					
					// and add exons that end or start at the border
					if(ex.getCodingStop() == nRegionStart || ex.getCodingStart() == nRegionEnd)
						vcTargetExons.add(ex);
				}
			}
			else
			{
				// add novel exon to the left
				int nExonLength 	= nMinDrawLength;
				int nExEnd	 		= junTarget.m_nStart;
				int nExStart 		= nExEnd - nExonLength+1;
				Region modifiedExon = null; 
				if(i == 0)
					modifiedExon = new Region(nExStart, nExEnd, fOffsetY, true, false, nExStart, nExonLength);
				else
					modifiedExon = new Region(nExStart, nExEnd, fOffsetY+40, true, false, nExStart, nExonLength);
				modifiedExon.SetIsNovel(true);
				modifiedExon.m_bIsIncomplete = true;
				
				// extend exon length if other exons are overlapping and longer
				modifiedExon.extendByExons(vcModifiedExons);				
				vcModifiedExons.add(modifiedExon);

				// add novel exon to the right
				nExStart 		= junTarget.m_nEnd;
				nExEnd	 		= nExStart + nExonLength;				

				// extend exon length if other exons are overlapping and longer
				if(i == 0)		
					modifiedExon = new Region(nExStart, nExEnd, fOffsetY, false, true, nExStart, nExonLength);
				else
					modifiedExon = new Region(nExStart, nExEnd, fOffsetY+40, false, true, nExStart, nExonLength);					
				
				modifiedExon.SetIsNovel(true);
				modifiedExon.m_bIsIncomplete = true;
				
				// extend exon length if other exons are overlapping and longer
				modifiedExon.extendByExons(vcModifiedExons);
				vcModifiedExons.add(modifiedExon);
			}
			
			// generate regions
			int nExonLength 	= -1;
			int nExStart 		= -1;			
			int nExEnd	 		= -1;
			Region modifiedExon	= null;
					
			int nCurEx = 0;
			for(Exon ex : vcTargetExons)
			{
				nExonLength 	= ex.getGenomicLength();
				nExStart 		= ex.getCodingStart();
				nExEnd	 		= ex.getCodingStop();
		
				boolean bIsLeftTerminalExon  = false;
				boolean bIsRightTerminalExon = false;
				if(nCurEx == 0)
				{
					nExonLength 	= Math.min(nMinDrawLength, ex.getGenomicLength());
					nExStart 		= ex.getGenomicStop() - nExonLength+1;
					nExEnd	 		= ex.getGenomicStop();
					bIsLeftTerminalExon = true;
				}
				if(nCurEx == vcTargetExons.size()-1)
				{
					nExonLength 	= Math.min(nMinDrawLength, ex.getGenomicLength());
					nExStart 		= ex.getGenomicStart();		
					nExEnd	 		= nExStart + nMinDrawLength;
					bIsRightTerminalExon = true;
				}
				
				if(i == 0)
					modifiedExon	= new Region(nExStart, nExEnd, fOffsetY, bIsLeftTerminalExon, bIsRightTerminalExon, ex.getGenomicStart(), ex.getGenomicLength());
				else
					modifiedExon	= new Region(nExStart, nExEnd, fOffsetY+40, bIsLeftTerminalExon, bIsRightTerminalExon, ex.getGenomicStart(), ex.getGenomicLength());

				// extend exon length if other exons are overlapping and longer
				modifiedExon.extendByExons(vcModifiedExons);
				vcModifiedExons.add(modifiedExon);

				nCurEx++;
			}
		}

		//############################################################################
		//                            specify regions
		//############################################################################
		// a region is defined by all overlapping exons
		// regions are separated by an intron of undefined length (40% plot width)
		//############################################################################

		TreeSet<Region> vcRegions	= new TreeSet<Region>();
		TreeSet<Region> vcUsedExons = new TreeSet<Region>();
		
		double fMinBasePosition = Double.MAX_VALUE;
		double fMaxBasePosition = 0.0;
		
		// repeat until all exons were used
		for(Region ex : vcModifiedExons)
		{
			// skip already used exons
			if(vcUsedExons.contains(ex))
				continue;
			
			double fRegionStart = 0.0;
			double fRegionEnd	= 0.0;
			
			fRegionStart= ex.m_fStart;
			fRegionEnd	= ex.m_fEnd;
			
			Region region = new Region(fRegionStart, fRegionEnd, 0.0, false, false, 0.0, 0.0);
			
			boolean bModified = true;
			
			TreeSet<Region> vcTmpUsedExons = new TreeSet<Region>();
			
			while(bModified)
			{
				bModified = false;
				
				for(Region ex2 : vcModifiedExons)
				{
					if(ex.equals(ex2) || vcTmpUsedExons.contains(ex2))
						continue;
					
					fRegionStart= ex2.m_fStart;
					fRegionEnd	= ex2.m_fEnd;
	
					// cut large terminal exons
					if(fRegionStart <= region.m_fEnd && fRegionEnd >= region.m_fStart)
					{
						region.m_fStart = Math.min(fRegionStart, 	region.m_fStart);
						region.m_fEnd	= Math.max(fRegionEnd, 		region.m_fEnd);
						
						vcTmpUsedExons.add(ex2);
						bModified = true;
						break;
					}
				}
			}
			
			for(Region reg : vcTmpUsedExons)
				vcUsedExons.add(reg);
			
			fMinBasePosition = Math.min(fMinBasePosition, region.m_fStart);
			fMaxBasePosition = Math.max(fMaxBasePosition, region.m_fEnd);

			vcRegions.add(region);
		}
		
		//######################################################################
		//    calculate how many pixels belong to each intron and region (exon)
		//######################################################################
		// 40% will be used by introns
		double fIntronPixelWidth = (nMaxWidth * 0.4) / (vcRegions.size()-1);
		
		// 60% will be used by regions
		double fRegionPixelWidth = (nMaxWidth * 0.6);

		//####################################################
		//    define pixel positions for regions and exons
		//####################################################
		TreeSet<Region>	vcExonsScreenPositions 	 = new TreeSet<Region>();
		TreeMap<Region, Region> mapExonBasePositionToExonScreenPosition = new TreeMap<Region, Region>(); 
		TreeMap<Region, Region>	mapScreenPositionRegionsToBasePositionRegions = new TreeMap<Region, Region>();
		
		double fScreenStart = nOffsetX;
		
		// get a rough estimate of the pixels per base value
		double fTotalRegionLengthInBp = 0;
		for(Region region : vcRegions)
		{
			fTotalRegionLengthInBp += region.m_fLength;
		}
		double fPixelsPerBase = fRegionPixelWidth / fTotalRegionLengthInBp;
		
		// define region screen positions
		// right now, regions are just exonic, but we want to add some coverage from the introns as well
		for(Region region : vcRegions)
		{
			// extend regions by intronic coverage, keep 20 pixel as spacer
			// add at most 100 bp!
			double fIntronLength = Math.min(100.0, (fIntronPixelWidth*0.4) / fPixelsPerBase);

			// the first region has only downstream intron coverage
			if(region.m_fStart == fMinBasePosition)
			{
				region.m_fEnd += fIntronLength;
			}
			// the last region has only upstream intron coverage
			else if(region.m_fEnd == fMaxBasePosition)
			{
				region.m_fStart -= fIntronLength;
			}
			// middle exons have upstream and downstream intron coverage
			else
			{
				region.m_fStart -= fIntronLength;
				region.m_fEnd   += fIntronLength;
			}
		}
/*		
		System.out.println("regions before merging:");
		for(Region reg : vcRegions)
			System.out.println((int)reg.m_fStart + " - " + (int)reg.m_fEnd);
*/	
		// merge overlapping regions (basespace)
		boolean bMerged = true;
		while(bMerged)
		{
			bMerged = false;
			
			outer_loop: for(Region region : vcRegions)
			{
				for(Region region2 : vcRegions)
				{
					// skip if it's the same region
					if(region.equals(region2))
						continue;
					
					// merge overlapping region
					if(region.overlaps(region2))
					{
						region.m_fStart = Math.min(region.m_fStart, region2.m_fStart);
						region.m_fEnd   = Math.max(region.m_fEnd, region2.m_fEnd);
						region.m_bIsLeftTerminalExon  = region.m_bIsLeftTerminalExon  | region2.m_bIsLeftTerminalExon;
						region.m_bIsRightTerminalExon = region.m_bIsRightTerminalExon | region2.m_bIsRightTerminalExon;
						region.m_fLength = region.m_fEnd - region.m_fStart +1;
						vcRegions.remove(region2);
						bMerged = true;
						break outer_loop;
					}
				}
			}
		}
/*	
		System.out.println("regions after merging:");
		for(Region reg : vcRegions)
			System.out.println((int)reg.m_fStart + " - " + (int)reg.m_fEnd);
*/
		// combine length of regions and calculate the number of bases per pixel
		fTotalRegionLengthInBp = 0;
		for(Region region : vcRegions)
		{
			fTotalRegionLengthInBp += region.m_fLength;
		}
		fPixelsPerBase = fRegionPixelWidth / fTotalRegionLengthInBp;
		
		// calculate position of regions (screen space)
		double fScreenStartTmp = fScreenStart;
		for(Region region : vcRegions)
		{
			double fRegionScreenSize = region.m_fLength * fPixelsPerBase;
			double fScreenEnd = fScreenStartTmp + fRegionScreenSize;

			Region screenPos = new Region(fScreenStartTmp, fScreenEnd, 0.0, false, false, 0.0, 0.0);
			mapScreenPositionRegionsToBasePositionRegions.put(region, screenPos);
			
			fScreenStartTmp += fRegionScreenSize + fIntronPixelWidth;
		}

		for(Region region : vcRegions)
		{
			double fRegionScreenSize = region.m_fLength * fPixelsPerBase;
			
			// define start position for each exon in the region
			for(Region ex : vcModifiedExons)
			{
				if(region.overlaps(ex))
				{
					double fExonStart = Math.min(ex.m_fStart, region.m_fStart);
					fExonStart = Math.max(fExonStart, ex.m_fOriginalStart);					
					
					double fExonEnd = Math.max(ex.m_fEnd, region.m_fEnd);
					fExonEnd		= Math.min(ex.m_fOriginalStart+ex.m_fOriginalLength-1, fExonEnd);
					double fExonWidth = fExonEnd - fExonStart +1;

					// transform to pixel coordinates
					fExonStart = fScreenStart + (fExonStart - region.m_fStart) * fPixelsPerBase;
					fExonEnd   = fExonStart + fExonWidth*fPixelsPerBase;
					
					Region screenPos = new Region(fExonStart, fExonEnd, ex.m_fOffsetY, ex.m_bIsLeftTerminalExon, ex.m_bIsRightTerminalExon, ex.m_bIsIncomplete);

					screenPos.SetIsNovel(ex.m_bIsNovel);
					
					vcExonsScreenPositions.add(screenPos);
					mapExonBasePositionToExonScreenPosition.put(screenPos, ex);
				}
			}
			
			fScreenStart += fRegionScreenSize + fIntronPixelWidth;
		}
		
		//####################################################
		//        define pixel positions for introns
		//####################################################
		TreeMap<Double, TreeSet<Region>> mapExonsToYOffset = new TreeMap<Double, TreeSet<Region>>();
		for(Region ex : vcExonsScreenPositions)
		{
			double fCurOffset = ex.m_fOffsetY;
			if(mapExonsToYOffset.containsKey(fCurOffset))
			{
				mapExonsToYOffset.get(fCurOffset).add(ex);
			}
			else
			{
				TreeSet<Region> vcTmp = new TreeSet<Region>();
				vcTmp.add(ex);
				mapExonsToYOffset.put(fCurOffset, vcTmp);
			}
		}

		TreeSet<Region> vcIntronScreenPositions = new TreeSet<Region>();
		TreeMap<Region, Region> mapIntronBasePositionToScreenPosition = new TreeMap<Region, Region>(); 
		for(double fCurOffset : mapExonsToYOffset.keySet())
		{
			TreeSet<Region> vcExons = mapExonsToYOffset.get(fCurOffset);
			
			double fIntronStart = -1.0;	// screen position
			double fIntronEnd   = -1.0;	// screen position
			int nIntronStart 	= -1;	// base position
			int nIntronEnd		= -1;	// base position	
			
			for(Region ex : vcExons)
			{
				// get base position
				Region regionBasePosition = mapExonBasePositionToExonScreenPosition.get(ex);
				
				if(fIntronStart == -1.0)
				{
					fIntronStart = ex.m_fEnd+1;
					nIntronStart = (int)regionBasePosition.m_fEnd;
				}
				else
				{
					fIntronEnd = ex.m_fStart-1;
					nIntronEnd = (int)regionBasePosition.m_fStart;
					
					Region intron = new Region(fIntronStart, fIntronEnd, fCurOffset, false, false, false);
					vcIntronScreenPositions.add(intron);
					
					if(	(nIntronStart == result.GetInclusionJunctionStart() && nIntronEnd == result.GetInclusionJunctionEnd()) ||
						(nIntronStart == result.GetExclusionJunctionStart() && nIntronEnd == result.GetExclusionJunctionEnd()) )
							intron.SetIsTarget(true);
					
					// check whether it is a novel junction
					boolean bIsNovel = true;
					isoform_loop : for(String strIsoform : dataSupplier.GetIsoformNames())
					{
						for(CountElement jun : dataSupplier.GetJunctionsForIsoform(strIsoform))
						{
							if(jun.m_nStart == nIntronStart && jun.m_nEnd == nIntronEnd)
							{
								bIsNovel = false;
								break isoform_loop;
							}
						}
					}
					
					Region intronBasePosition = new Region(nIntronStart, nIntronEnd, fCurOffset, false, false, false);
					mapIntronBasePositionToScreenPosition.put(intron, intronBasePosition);
					
					intron.SetIsNovel(bIsNovel);
					
					fIntronStart = ex.m_fEnd+1;
					nIntronStart = (int)regionBasePosition.m_fEnd;
				}
			}
		}
		
		/*
		System.out.println("BASE Regions: " + mapScreenPositionRegionsToBasePositionRegions.size());
		for(Region reg : mapScreenPositionRegionsToBasePositionRegions.keySet())
			System.out.println((int)reg.m_fStart + " - " + (int)reg.m_fEnd + " -> " + (int)mapScreenPositionRegionsToBasePositionRegions.get(reg).m_fStart + " " + (int)mapScreenPositionRegionsToBasePositionRegions.get(reg).m_fEnd);
		System.out.println("BASE position exons:");
		for(Region reg : mapExonBasePositionToExonScreenPosition.keySet())
			System.out.println((int)reg.m_fStart + " - " + (int)reg.m_fEnd + " -> " + (int)mapExonBasePositionToExonScreenPosition.get(reg).m_fStart + " " + (int)mapExonBasePositionToExonScreenPosition.get(reg).m_fEnd);
		
		System.out.println("SCREEN position exons:");
		for(Region reg : vcExonsScreenPositions)
			System.out.println((int)reg.m_fStart + " - " + (int)reg.m_fEnd);
		*/

		//######################################################################
		//                          draw coverage
		//######################################################################
		nOffsetY = nMargin;
		
		//#################################
		//    prepare coverage container
		//#################################		
		TreeMap<Region, TreeMap<String, double[][]>> mapCoverageToConditionsToRegionStart = new TreeMap<Region, TreeMap<String, double[][]>>();
		
		double fYAxisPos = Double.MAX_VALUE;

		for(Region regionBasePosition : mapScreenPositionRegionsToBasePositionRegions.keySet())
		{
			int nBasePositionStart 		= (int)regionBasePosition.m_fStart;
			int nBasePositionEnd   		= (int)regionBasePosition.m_fEnd;
			int nBaseLength		  	  	= nBasePositionEnd - nBasePositionStart + 1;
			
			double fScreenPositionStart = mapScreenPositionRegionsToBasePositionRegions.get(regionBasePosition).m_fStart;
						
			fYAxisPos = Math.min(fYAxisPos, fScreenPositionStart);
			
			TreeMap<String, double[][]> mapCoverageToConditions = new TreeMap<String, double[][]>();
		
			for(String strCondition : mapSamplesToConditions.keySet())
			{
				// determine number of valid samples
				int nValidSamples = 0;
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					if(vcSelectedSamples.contains(strSample))
						nValidSamples++;
				}
				
				if(nValidSamples == 0)
					continue;

				// obtain coverage for the region
				double pCoverage[][] = new double[nBaseLength][nValidSamples];

				int nSampleIdx = 0;
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					if(!vcSelectedSamples.contains(strSample))
						continue;
					
					// get coverage for the whole gene
					double pSampleCoverage[] = dataSupplier.GetGeneCoverageArrayForSample(strSample);
					int nArrayStart = dataSupplier.GetGeneStart();
					
					// try to read coverage from file if necessary
					if(pSampleCoverage == null)
					{
						nArrayStart = -1;
					}
										
					// if there is no specific coverage array for the exon group (in case the junction extends past the
					// gene boundaries) it is necessary to create one
					if(pSampleCoverage == null || nBasePositionStart < nArrayStart || nBasePositionEnd > dataSupplier.GetGeneEnd())
					{
						try
						{
							pSampleCoverage = m_app.GetCoverageForRegion(dataSupplier.GetReferenceName(), nBasePositionStart, nBasePositionEnd, strSample, true);
							nArrayStart = nBasePositionStart;
						}
						catch(Exception e)
						{
							e.printStackTrace();
							System.out.println("ERROR: could not retrieve coverage for: " + dataSupplier.GetReferenceName() + ":" + nBasePositionStart + "-" + nBasePositionEnd);	
							Messagebox.show("ERROR: could not retrieve coverage for: " + strSample + "" + e.getLocalizedMessage());
						}
					}

					// get coverage for the specified region
					for(int nPos=nBasePositionStart; nPos<=nBasePositionEnd; nPos++)
					{
						double fValue = pSampleCoverage[nPos - nArrayStart];
						if(m_bCoveragePlotUseLog2)
						{
							fValue = Math.log(fValue+1) / Math.log(2);
						}

						int nArrayPos = nPos - nBasePositionStart;	
						pCoverage[nArrayPos][nSampleIdx] = fValue;
					}
					
					nSampleIdx++;
				}
				
				mapCoverageToConditions.put(strCondition, pCoverage);
			}
			
			mapCoverageToConditionsToRegionStart.put(regionBasePosition, mapCoverageToConditions);
		}
		
		// find maximum value
		double fMaxCoverage = 0.0;
		for(Region regionBasePosition : mapScreenPositionRegionsToBasePositionRegions.keySet())
		{
			int nBasePositionStart 		= (int)regionBasePosition.m_fStart;
			int nBasePositionEnd   		= (int)regionBasePosition.m_fEnd;
			int nBaseLength		  	  	= nBasePositionEnd - nBasePositionStart + 1;
			TreeMap<String, double[][]> mapCoverageToConditions = mapCoverageToConditionsToRegionStart.get(regionBasePosition);
			
			for(String strCondition : mapSamplesToConditions.keySet())
			{
				if(!vcValidConditions.contains(strCondition))
					continue;
				
				double[][] pCoverage = mapCoverageToConditions.get(strCondition);
				
				for(int i=0; i<nBaseLength; i++)
				{
					double fVal = StatUtils.percentile(pCoverage[i], 75.0);
					fMaxCoverage = Math.max(fMaxCoverage, fVal);
				}
			}
		}
		
		// plot y-axis
		fMaxCoverage = AddYAxis(fYAxisPos, nOffsetY, nCoveragePlotHeight, fMaxCoverage, "Coverage", graph);
			
		for(Region regionBasePosition : mapScreenPositionRegionsToBasePositionRegions.keySet())
		{
			int nBasePositionStart 		= (int)regionBasePosition.m_fStart;
			int nBasePositionEnd   		= (int)regionBasePosition.m_fEnd;
			int nBaseLength		  	  	= nBasePositionEnd - nBasePositionStart + 1;
			double fScreenPositionStart	= mapScreenPositionRegionsToBasePositionRegions.get(regionBasePosition).m_fStart;
			
			TreeMap<String, double[][]> mapCoverageToConditions = mapCoverageToConditionsToRegionStart.get(regionBasePosition);
			
			// plot coverage
			for(String strCondition : mapSamplesToConditions.keySet())
			{
				if(!vcValidConditions.contains(strCondition))
					continue;
				
				double[][] pCoverage = mapCoverageToConditions.get(strCondition);
				
				double p25Coverage[] = new double[nBaseLength];
				double p50Coverage[] = new double[nBaseLength];
				double p75Coverage[] = new double[nBaseLength];
				
				for(int i=0; i<nBaseLength; i++)
				{
					p25Coverage[i] = StatUtils.percentile(pCoverage[i], 25.0);
					p50Coverage[i] = StatUtils.percentile(pCoverage[i], 50.0);
					p75Coverage[i] = StatUtils.percentile(pCoverage[i], 75.0);
				}
				
				// plot quartiles
				Color clrAlpha = new Color(m_mapColorsToConditions.get(strCondition).getRed(), m_mapColorsToConditions.get(strCondition).getGreen(), m_mapColorsToConditions.get(strCondition).getBlue(), 40);
				graph.setColor(clrAlpha);
				
				GeneralPath polyline = new GeneralPath(GeneralPath.WIND_EVEN_ODD, nBaseLength);

				for(int x=0; x<nBaseLength; x++)
				{
					double fY = nMargin + ((fMaxCoverage - p75Coverage[x]) / fMaxCoverage) * nCoveragePlotHeight;
					
					if(x==0)
						polyline.moveTo(fScreenPositionStart+x*fPixelsPerBase, fY);
					else
						polyline.lineTo(fScreenPositionStart+x*fPixelsPerBase, fY);

					polyline.lineTo(fScreenPositionStart+(x+1)*fPixelsPerBase, fY);
				}

				for(int x=nBaseLength-1; x>=0; x--)
				{									
					double fY = nMargin + ((fMaxCoverage - p25Coverage[x]) / fMaxCoverage) * nCoveragePlotHeight;
					
					polyline.lineTo(fScreenPositionStart+(x+1)*fPixelsPerBase, fY);
					polyline.lineTo(fScreenPositionStart+x*fPixelsPerBase, fY);
				}
				
				polyline.closePath();
				graph.fill(polyline);
				
				// plot median
				graph.setColor(m_mapColorsToConditions.get(strCondition));
				polyline = new GeneralPath(GeneralPath.WIND_EVEN_ODD, p50Coverage.length);

				for(int x=0; x<nBaseLength; x++)
				{
					double fY = nMargin + ((fMaxCoverage - p50Coverage[x]) / fMaxCoverage) * nCoveragePlotHeight;
					
					if(x==0)
						polyline.moveTo(fScreenPositionStart+x*fPixelsPerBase, fY);
					else
						polyline.lineTo(fScreenPositionStart+x*fPixelsPerBase, fY);

					polyline.lineTo(fScreenPositionStart+(x+1)*fPixelsPerBase, fY);
				}
				
				graph.draw(polyline);
			}
		}
		
		//######################################################################
		//                           draw exons
		//######################################################################	
		for(Region ex : vcExonsScreenPositions)
		{			
			if(ex.m_bIsNovel)
			{
				clr = new GradientPaint(0, (int)ex.m_fOffsetY, clrTargetNovelExonTop, 0, (int)ex.m_fOffsetY+20, clrTargetNovelExonBottom);
				graph.setPaint(clr);
			}
			else
			{
				clr = new GradientPaint(0, (int)ex.m_fOffsetY, clrTargetExonTop, 0, (int)ex.m_fOffsetY+20, clrTargetExonBottom);
				graph.setPaint(clr);
			}
			
			Rectangle2D.Double rect = new Rectangle2D.Double(ex.m_fStart, ex.m_fOffsetY, ex.m_fLength, 20); 
			graph.fill(rect);
			graph.setColor(Color.BLACK);
			graph.draw(rect);
			
			// add marker for shortened exons
			if(ex.m_bIsIncomplete)
			{
				if(ex.m_bIsLeftTerminalExon)
				{
					int nPos = (int)ex.m_fStart - 14;
					graph.drawImage(imgMarkLeft, nPos, (int)ex.m_fOffsetY-5, null);
				}
				else if(ex.m_bIsRightTerminalExon)
				{
					int nPos = (int)ex.m_fEnd - 6;
					graph.drawImage(imgMarkRight, nPos, (int)ex.m_fOffsetY-5, null);
				}
			}
		}
		
		//######################################################################
		//                           draw introns
		//######################################################################
		for(Region intron : vcIntronScreenPositions)
		{
			if(intron.m_bIsNovel)
			{
				if(intron.m_bIsTarget)
				{
					clr = new GradientPaint(0, (int)intron.m_fOffsetY+1, clrTargetNovelExonTop, 0, (int)intron.m_fOffsetY+13, clrTargetNovelExonBottom);
					graph.setPaint(clr);
				}
				else
				{
					clr = new GradientPaint(0, (int)intron.m_fOffsetY+1, clrNonTargetNovelExonTop, 0, (int)intron.m_fOffsetY+13, clrNonTargetNovelExonBottom);
					graph.setPaint(clr);
				}
			}
			else
			{
				if(intron.m_bIsTarget)
				{
					clr = new GradientPaint(0, (int)intron.m_fOffsetY+1, clrTargetExonTop, 0, (int)intron.m_fOffsetY+13, clrTargetExonBottom);
					graph.setPaint(clr);
				}
				else
				{
					clr = new GradientPaint(0, (int)intron.m_fOffsetY+1, clrNonTargetExonTop, 0, (int)intron.m_fOffsetY+13, clrNonTargetExonBottom);
					graph.setPaint(clr);
				}
			}
			
			Rectangle2D.Double rect = new Rectangle2D.Double(intron.m_fStart, intron.m_fOffsetY+7, intron.m_fLength, 6);
			graph.fill(rect);
			graph.setColor(Color.BLACK);
			graph.draw(rect);
			
			// add label
			String strJunctionName = dataSupplier.GetReferenceName() + ":" + (int)(mapIntronBasePositionToScreenPosition.get(intron).m_fStart+1) + "-" + (int)(mapIntronBasePositionToScreenPosition.get(intron).m_fEnd-1);
			if(intron.m_bIsTarget)
			{
				Font font = new Font("Lucida Sans", Font.BOLD, 11);
				graph.setFont(font);
				FontMetrics metrics = graph.getFontMetrics(font);
				
				graph.setColor(Color.black);
				
				int nStringWidth = metrics.stringWidth(strJunctionName);
				graph.drawString(strJunctionName, (int)(intron.m_fStart + intron.m_fLength*0.5 - nStringWidth *0.5), (int)intron.m_fOffsetY+30);
			}
			else
			{
				// can't show all junction names
				/*
				Font font = new Font("Lucida Sans", Font.PLAIN, 11);
				graph.setFont(font);
				FontMetrics metrics = graph.getFontMetrics(font);
				
				graph.setColor(Color.DARK_GRAY);
				
				int nStringWidth = metrics.stringWidth(strJunctionName);
				graph.drawString(strJunctionName, (int)(intron.m_fStart + intron.m_fLength*0.5 - nStringWidth *0.5), (int)intron.m_fOffsetY+30);
				*/
			}
			
			
		}
		
		//####################################################
		//             calculate junction coverage
		//####################################################
		// get size factors
		TreeMap<String, Double> mapSizeFactorsToSamples = projectModel.GetSizeFactors();
		
		// get junction coverages
		TreeMap<String, Integer> vcCountsA = dataSupplier.GetCountsForJunction(result.GetInclusionJunction());
		TreeMap<String, Integer> vcCountsB = dataSupplier.GetCountsForJunction(result.GetExclusionJunction());
		
		double fMaxCountA = 0;
		double fMaxCountB = 0;
		
		TreeMap<String, double[]> mapCountsToConditionsA = new TreeMap<String, double[]>();
		TreeMap<String, double[]> mapCountsToConditionsB = new TreeMap<String, double[]>();
		TreeMap<String, double[]> mapRatiosToConditions  = new TreeMap<String, double[]>();
		
		Vector<String> vcDataLabels = new Vector<String>();

		for(String strCondition : mapSamplesToConditions.keySet())
		{
			if(!vcValidConditions.contains(strCondition))
				continue;
			
			Vector<Double> vcValuesA = new Vector<Double>();
			Vector<Double> vcValuesB = new Vector<Double>();
			
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				if(vcSelectedSamples.contains(strSample))
				{
					double fCount = 0.0;
					if(vcCountsA.containsKey(strSample))
						fCount = vcCountsA.get(strSample) / mapSizeFactorsToSamples.get(strSample);
					vcValuesA.add(fCount);
					
					if(fCount > fMaxCountA)
						fMaxCountA = fCount;
					
					fCount = 0.0;
					if(vcCountsB.containsKey(strSample))
						fCount = vcCountsB.get(strSample) / mapSizeFactorsToSamples.get(strSample);
					vcValuesB.add(fCount);
					
					if(fCount > fMaxCountB)
						fMaxCountB = fCount;
					
					vcDataLabels.add(strSample);
				}
			}
			
			double pValuesA[] = new double[vcValuesA.size()];
			for(int i=0; i<vcValuesA.size(); i++)
			{
				pValuesA[i] = vcValuesA.get(i);
			}
			
			double pValuesB[] = new double[vcValuesB.size()];
			for(int i=0; i<vcValuesB.size(); i++)
			{
				pValuesB[i] = vcValuesB.get(i);
			}
			
			double pCountsRatio[] = new double[pValuesA.length];
			for(int k=0; k<pValuesA.length; k++)
			{
				pCountsRatio[k] = pValuesA[k] / (pValuesA[k] + pValuesB[k]); 
			}
			
			mapCountsToConditionsA.put(strCondition, pValuesA);
			mapCountsToConditionsB.put(strCondition, pValuesB);
			mapRatiosToConditions.put(strCondition, pCountsRatio);
		}

		//######################################################
		//                  draw box plots
		//######################################################
		nOffsetX = nMargin + 50;
		nOffsetY = nMargin + 360;
		
		String strJunctionNameA = dataSupplier.GetReferenceName() + ":" + (result.GetInclusionJunctionStart()+1) + "-" + (result.GetInclusionJunctionEnd()-1);
		String strJunctionNameB = dataSupplier.GetReferenceName() + ":" + (result.GetExclusionJunctionStart()+1) + "-" + (result.GetExclusionJunctionEnd()-1);
		
		// draw a plot for each junction and one plot for the ratios
		BoxPlot(mapCountsToConditionsA, strJunctionNameA, graph, nOffsetX, nOffsetY, nPlotHeight, "Coverage");
		nOffsetX += nBoxPlotWidth + nInterPlotDistance;
		
		BoxPlot(mapCountsToConditionsB, strJunctionNameB, graph, nOffsetX, nOffsetY, nPlotHeight, "Coverage");
		nOffsetX += nBoxPlotWidth + nInterPlotDistance;
	
		BoxPlot(mapRatiosToConditions, "Ratio", graph, nOffsetX, nOffsetY, nPlotHeight, "Ratio");
		nOffsetX += nBoxPlotWidth + nInterPlotDistance;
		
		Vector<Area> vcAreas = Plot2D(mapCountsToConditionsA, mapCountsToConditionsB, vcDataLabels, strJunctionNameA, strJunctionNameB, "2D plot", graph, nOffsetX, nOffsetY, n2DPlotWidth, nPlotHeight, false);

		Imagemap imgMap = new Imagemap();
		imgMap.setWidth(nMaxPlotWidth+"px");
		imgMap.setHeight(nMaxHeight+"px");
		imgMap.setContent(img);
		imgMap.setParent(layout);
				
		for(Area area : vcAreas)
		{
			area.setParent(imgMap);
		}
	}
	
	/**
	 *    Similar to above, but instead of showing alternative splicing events that involve additional/missing exons, this plot
	 *    focuses on retained introns (which also makes the plot function less complicated)
	 */
	public void DrawRetainedIntronCloseup(AnalysisResult result)
	{
		ProjectModel projectModel				= m_app.GetProjectModel();
		DataSupplier dataSupplier 				= m_app.GetDataSupplier();
		String strSelectedConditionType 		= m_app.GetSelectedConditionType();
		TreeSet<String> vcSelectedSamples		= m_app.GetSelectedSamples();
		TreeMap<String, TreeSet<String>> mapSamplesToConditions = projectModel.GetSamplesPerCondition(strSelectedConditionType);
		
		//###################################
		//       prepare popup window
		//###################################
		int nMaxWidth  		= (int)(m_nClientWindowWidth*0.5);
		int nWindowHeight 	= (int)(m_nClientWindowHeight*0.5);
		
		Window windowPopup = new Window();
		windowPopup.setParent(m_app);
		windowPopup.setWidth((nMaxWidth+20) +"px");
		windowPopup.setHeight(nWindowHeight +"px");
		windowPopup.setTitle("Junction Heatmap");
		windowPopup.setSizable(true);
		windowPopup.setClosable(true);
		windowPopup.setMaximizable(true);
		windowPopup.setBorder(true);
		windowPopup.setPosition("center,center");
		windowPopup.setVisible(true);
		windowPopup.doPopup();
		windowPopup.setTopmost();
		
		Hlayout layout = new Hlayout();
		layout.setWidth("100%");
		layout.setHeight("100%");
		layout.setStyle("overflow:auto;");
		layout.setParent(windowPopup);
		
		int nNameOffset = 180; // subtract 100 px for the isoform names
		int nMargin = 10;
		
		// settings for the plots
		int nPlotHeight = 300;

		int nCoveragePlotHeight = 200;
		
		int nMaxHeight = nPlotHeight + nCoveragePlotHeight + 260;
		
		// Colors
		Color clrExonTop 				= new Color(217,217,255,255);	// BLUE
		Color clrExonBottom 			= new Color(120,120,230,255);	// BLUE
		
		//######################################################
		//                 prepare graph
		//######################################################
		BufferedImage img = new BufferedImage(nMaxWidth, nMaxHeight, BufferedImage.TYPE_INT_RGB);
		Graphics2D graph = img.createGraphics();
		
		Font fontNumbers = new Font("Lucida Sans", Font.BOLD, 11);
		graph.setFont(fontNumbers);
		
		graph.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		// load noncontinous-mark images
		String strPath = Executions.getCurrent().getDesktop().getWebApp().getRealPath("/");
		BufferedImage imgMark = null;
		try
		{
			imgMark = ImageIO.read(new File(strPath + "/img/noncontinuous_left.png"));
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.out.println("could not open image: " + strPath + "/img/noncontinuous_left.png");
			return;
		}
		BufferedImage imgMarkLeft = new BufferedImage(42, 99, BufferedImage.TYPE_INT_ARGB);
		
		Graphics2D g = imgMarkLeft.createGraphics();
		g.drawImage(imgMark, 0, 0, 20, 30, null);
		g.dispose();
		
		try
		{
			imgMark = ImageIO.read(new File(strPath + "/img/noncontinuous_right.png"));
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.out.println("could not open image: " + strPath + "/img/noncontinuous_right.png");
			return;
		}
		BufferedImage imgMarkRight = new BufferedImage(42, 99, BufferedImage.TYPE_INT_ARGB);
		
		g = imgMarkRight.createGraphics();
		g.drawImage(imgMark, 0, 0, 20, 30, null);
		g.dispose();

		// fill background
		graph.setColor(Color.white);
		graph.fillRect(0, 0, nMaxWidth, nMaxHeight);
		
		// subtract the margin from the canvas width and reserve some space for the exon names
		
		// reserve some space for the isoform name
		int nOffsetX = nMargin + nNameOffset;
		
		//######################################################################
		//                          draw coverage
		//######################################################################
		int    nOffsetY 		= nMargin;
		double fMaxCoverage 	= 0.0;
		
		// get location of alternatively spliced exon groups
		int nBasePositionStart 		= result.GetStartA() - 100;
		int nBasePositionEnd   		= result.GetEndA() + 100;
		int nBaseLength		  	  	= nBasePositionEnd - nBasePositionStart + 1;
		
		double fPixelsPerBase = (double)(nMaxWidth - (nMargin*2 + nNameOffset)) / (double)nBaseLength;
		
		//#################################
		//    prepare coverage container
		//#################################		
		TreeMap<String, double[][]> mapCoverageToConditions = new TreeMap<String, double[][]>();
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			// determine number of valid samples
			int nValidSamples = 0;
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				if(vcSelectedSamples.contains(strSample))
					nValidSamples++;
			}

			// obtain coverage for the region
			double pCoverage[][] = new double[nBaseLength][nValidSamples];

			int nSampleIdx = 0;
			for(String strSample : mapSamplesToConditions.get(strCondition))
			{
				if(!vcSelectedSamples.contains(strSample))
					continue;
				
				// get coverage for the whole gene
				double pSampleCoverage[] = dataSupplier.GetGeneCoverageArrayForSample(strSample);
				
				// get coverage for the specified region
				for(int nPos=nBasePositionStart; nPos<=nBasePositionEnd; nPos++)
				{
					double fValue = pSampleCoverage[nPos - dataSupplier.GetGeneStart()];

					if(m_bCoveragePlotUseLog2)
					{
						fValue = Math.log(fValue+1) / Math.log(2);
					}

					int nArrayPos = nPos - nBasePositionStart;	
					pCoverage[nArrayPos][nSampleIdx] = fValue;

					fMaxCoverage = Math.max(fMaxCoverage, fValue);
				}
				
				nSampleIdx++;
			}
			
			mapCoverageToConditions.put(strCondition, pCoverage);
		}
		
		// plot y-axis
		fMaxCoverage = AddYAxis(nOffsetX, nOffsetY, nCoveragePlotHeight, fMaxCoverage, "Coverage", graph);
		
		// plot coverage
		for(String strCondition : mapSamplesToConditions.keySet())
		{
			double[][] pCoverage = mapCoverageToConditions.get(strCondition);
			
			double p25Coverage[] = new double[nBaseLength];
			double p50Coverage[] = new double[nBaseLength];
			double p75Coverage[] = new double[nBaseLength];
			
			for(int i=0; i<nBaseLength; i++)
			{
				p25Coverage[i] = StatUtils.percentile(pCoverage[i], 25.0);
				p50Coverage[i] = StatUtils.percentile(pCoverage[i], 50.0);
				p75Coverage[i] = StatUtils.percentile(pCoverage[i], 75.0);
			}
			
			// plot quartiles
			Color clrAlpha = new Color(m_mapColorsToConditions.get(strCondition).getRed(), m_mapColorsToConditions.get(strCondition).getGreen(), m_mapColorsToConditions.get(strCondition).getBlue(), 40);
			graph.setColor(clrAlpha);
			
			GeneralPath polyline = new GeneralPath(GeneralPath.WIND_EVEN_ODD, nBaseLength);

			for(int x=0; x<nBaseLength; x++)
			{
				double fY = nMargin + ((fMaxCoverage - p75Coverage[x]) / fMaxCoverage) * nCoveragePlotHeight;
				
				if(x==0)
					polyline.moveTo(nOffsetX+x*fPixelsPerBase, fY);
				else
					polyline.lineTo(nOffsetX+x*fPixelsPerBase, fY);

				polyline.lineTo(nOffsetX+(x+1)*fPixelsPerBase, fY);
			}

			for(int x=nBaseLength-1; x>=0; x--)
			{									
				double fY = nMargin + ((fMaxCoverage - p25Coverage[x]) / fMaxCoverage) * nCoveragePlotHeight;
				
				polyline.lineTo(nOffsetX+(x+1)*fPixelsPerBase, fY);
				polyline.lineTo(nOffsetX+x*fPixelsPerBase, fY);
			}
			
			polyline.closePath();
			graph.fill(polyline);
			
			// plot median
			graph.setColor(m_mapColorsToConditions.get(strCondition));
			polyline = new GeneralPath(GeneralPath.WIND_EVEN_ODD, p50Coverage.length);

			for(int x=0; x<nBaseLength; x++)
			{
				double fY = nMargin + ((fMaxCoverage - p50Coverage[x]) / fMaxCoverage) * nCoveragePlotHeight;
				
				if(x==0)
					polyline.moveTo(nOffsetX+x*fPixelsPerBase, fY);
				else
					polyline.lineTo(nOffsetX+x*fPixelsPerBase, fY);

				polyline.lineTo(nOffsetX+(x+1)*fPixelsPerBase, fY);
			}
			
			graph.draw(polyline);
		}

		//######################################################################
		//                           draw exons
		//######################################################################
		nOffsetY = nMargin + nCoveragePlotHeight + 20;

		double fIntronStart = 0.0;
		double fIntronEnd	= 0.0;
		
		double fExonStart 	= nOffsetX;
		double fExonEnd 	= nOffsetX + 100 * fPixelsPerBase;
		double fExonLength	= fExonEnd - fExonStart+1;
		
		fIntronStart = fExonEnd+1;
		
		GradientPaint clr = new GradientPaint(0, nOffsetY, clrExonTop, 0, nOffsetY+20, clrExonBottom);
		graph.setPaint(clr);
			
		Rectangle2D.Double rect = new Rectangle2D.Double(fExonStart, nOffsetY, fExonLength, 20); 
		graph.fill(rect);
		graph.setColor(Color.BLACK);
		graph.draw(rect);
			
		// add marker for shortened exons
		int nPos = (int)fExonStart - 14;
		graph.drawImage(imgMarkLeft, nPos, nOffsetY-5, null);
		
		nOffsetX += fExonLength;
		
		fExonStart 	= nOffsetX + (result.GetEndA()-result.GetStartA()+1) * fPixelsPerBase;
		fExonEnd 	= fExonStart + 100 * fPixelsPerBase;
		fExonLength	= fExonEnd - fExonStart+1;
		
		fIntronEnd = fExonStart-1;

		graph.setPaint(clr);
			
		rect = new Rectangle2D.Double(fExonStart, nOffsetY, fExonLength, 20); 
		graph.fill(rect);
		graph.setColor(Color.BLACK);
		graph.draw(rect);

		nPos = (int)fExonEnd - 6;
		graph.drawImage(imgMarkRight, nPos, nOffsetY-5, null);
		
		//######################################################################
		//                           draw introns
		//######################################################################
		double fIntronLength = fIntronEnd - fIntronStart +1;
		
		graph.setPaint(clr);
		rect = new Rectangle2D.Double(fIntronStart, nOffsetY+7, fIntronLength, 6);
		graph.fill(rect);
		graph.setColor(Color.BLACK);
		graph.draw(rect);
		
		// add label
		String strJunctionName = dataSupplier.GetReferenceName() + ":" + result.GetStartA() + "-" + result.GetEndA();

		Font font = new Font("Lucida Sans", Font.BOLD, 11);
		graph.setFont(font);
		FontMetrics metrics = graph.getFontMetrics(font);

		int nStringWidth = metrics.stringWidth(strJunctionName);
		graph.drawString(strJunctionName, (int)(fIntronStart + fIntronLength*0.5 - nStringWidth *0.5), nOffsetY+30);

//		Plot2D(mapCountsToConditionsA, mapCountsToConditionsB, strJunctionNameA, strJunctionNameB, "2D plot", graph, nOffsetX, nOffsetY, n2DPlotWidth, nPlotHeight, false);

		Imagemap imgMap = new Imagemap();
		imgMap.setWidth(nMaxWidth+"px");
		imgMap.setHeight(nMaxHeight+"px");
		imgMap.setContent(img);
		imgMap.setParent(layout);
	}
	
	/** Adds overlapping sense and antisense transcripts to the isoform plot */
	public int DrawOverlappingTranscripts(int nX, int nY, int nIntronLength, double fShrinkageFactor, int nVerticalSpaceBetweenIsoforms, Graphics2D graph) throws FileNotFoundException
	{
		DataSupplier dataSupplier 				= m_app.GetDataSupplier();
		
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
			
			if(g.getStrandStringID() == dataSupplier.GetStrand())
			{
				map1 = mapExonGroupToOverlappingRegionsBlue;
				map2 = mapExonGroupToOverlappingRegionsDarkBlue;
			}
			
			// see which exon groups are overlapping
			for(ExonGroup grp : dataSupplier.GetExonGroups())
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
						if(ex.getGenomicStart() <= grp.getGenomicStopOfGroup() && ex.getGenomicStop() >= grp.getGenomicStartOfGroup())
						{
							if(map2.containsKey(grp))
							{
								TreeSet<CountElement> vcWarningRed = map2.get(grp);
								
								CountElement r = new CountElement();
								r.m_nStart = Math.max(ex.getGenomicStart(), grp.getGenomicStartOfGroup());
								r.m_nEnd = Math.min(ex.getGenomicStop(), grp.getGenomicStopOfGroup());
								vcWarningRed.add(r);
							}
							else
							{
								TreeSet<CountElement> vcWarningRed = new TreeSet<CountElement>();

								CountElement r = new CountElement();
								r.m_nStart = Math.max(ex.getGenomicStart(), grp.getGenomicStartOfGroup());
								r.m_nEnd = Math.min(ex.getGenomicStop(), grp.getGenomicStopOfGroup());
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
			Font fontNumbers = new Font("Lucida Sans", Font.BOLD, 11);
			graph.setFont(fontNumbers);
			FontMetrics fontMetrics = graph.getFontMetrics(fontNumbers);
			
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
					String strText = "AS transcripts";
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
					String strText = "S transcripts";
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
			
			fontNumbers = new Font("Lucida Sans", Font.BOLD, 10);
			graph.setFont(fontNumbers);
			fontMetrics = graph.getFontMetrics(fontNumbers);

			// draw regions
			int nXTmp = nX;
			for(ExonGroup grp : dataSupplier.GetExonGroups())
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
					
					String strText = "";
					int nTextWidth = 0;
					int nLetterWidth = fontMetrics.stringWidth(">");
					
					// add orientation
					while(nTextWidth < nWidth-20)
					{
						if(((dataSupplier.GetStrand() == '+' && (i == 2 || i == 3))) || (dataSupplier.GetStrand() == '-' && (i == 0 || i == 1)))
						{
							strText += ">";
							nTextWidth += nLetterWidth;
						}
						else
						{
							strText += "<";
							nTextWidth += nLetterWidth;
						}
					}

					int nTextHeight = fontMetrics.getHeight();
					
					int nXOffset = (int)((nWidth - nTextWidth) * 0.5);
					
					if(i == 3) graph.setPaint(Color.WHITE);
					graph.drawString(strText, nXTmp + nOffset + nXOffset, (int)(nY+6+nTextHeight*0.5));
				}
				
				nXTmp += nGrpWidth + nIntronLength;
			}
		}
		
		nY += nVerticalSpaceBetweenIsoforms;
		
		return nY;
	}
	
	/** Returns a list of overlapping transcripts */
	public TreeSet<Gene> GetOverlappingTranscripts() throws FileNotFoundException
	{
		DataSupplier dataSupplier 				= m_app.GetDataSupplier();
		String strFileGTF						= m_app.GetGTFFile();
		
		if(dataSupplier.GetGene() == null)
		{
			System.out.println("GetOverlappingTranscripts() -> no gene selected");
			return null;
		}
		
		TreeSet<Gene> vcGenes = null;

		// check if reference file exists
		File pFile = new File(strFileGTF);
		if(!pFile.exists())
		{
			System.out.println("Failed to open file: " + strFileGTF);
			return null;
		}
		
		RandomAccessGFFReader gffReader = m_app.GetGFFReader();	
		if(gffReader == null || !gffReader.GetFileName().equals(strFileGTF))
		{
			m_app.LoadGeneAnnotation(strFileGTF);
		}

		vcGenes = gffReader.GetGenesForRange(dataSupplier.GetReferenceName(), dataSupplier.GetGeneStart()-2000, dataSupplier.GetGeneEnd()+2000);
		
		// remove current gene
		if(vcGenes.contains(dataSupplier.GetGene()))
			vcGenes.remove(dataSupplier.GetGene());

		return vcGenes;
	}
	
	/** Invokes drawing of the coverage (if requested) and isoform plot */
	public boolean DrawPlots() throws IOException, SQLException
	{
		DataSupplier dataSupplier = m_app.GetDataSupplier();
		
		if(dataSupplier.GetGene() == null)
		{
			Messagebox.show("No gene specified");
			return false;
		}
		
		// clear old plots
		m_app.GetIsoformPlotRegion().getChildren().clear();		
		
		if(m_bCoverageRequiresRedraw)
		{
			m_app.GetCoveragePlotRegion().getChildren().clear();

			// prepares the coverage image map
			DrawExtendedCoveragePlot(false);
		}

		// prepare the isoform image map
		DrawIsoforms();
		
		// put image maps into layout
		m_app.GetImageMapForIsoforms();
		
		return true;
	}

}
