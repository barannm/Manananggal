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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.Map;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import org.apache.commons.math3.stat.StatUtils;
import org.zkoss.zk.ui.Components;
import org.zkoss.zk.ui.event.Event;
import org.zkoss.zk.ui.event.EventListener;
import org.zkoss.zk.ui.event.Events;
import org.zkoss.zul.Borderlayout;
import org.zkoss.zul.Button;
import org.zkoss.zul.East;
import org.zkoss.zul.Imagemap;
import org.zkoss.zul.Messagebox;
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

public class PlotFactoryGTEX
{
	SplicingWebApp 												m_app;
	private String 												m_strPathToGeneData;
	private TreeSet<String> 									m_vcPreviouslySelectedGTEXTissues;
	private Window 												m_windowPopup;
	
	
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
		
	PlotFactoryGTEX(SplicingWebApp app)
	{
		m_app = app;
		
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
	
	public void SetPaths(String strPathToGeneData)
	{
		m_strPathToGeneData = strPathToGeneData;
	}

	public TreeMap<String, Vector<Double>> GetGTEXDataForGene()
	{
		DataSupplier dataSupplier = m_app.GetDataSupplier();
		
		GeneIdentifier id = dataSupplier.GetGeneIdentifier();
		if(id == null)
			return null;
		
		String strEntrezID = id.m_strEntrezGeneID;

		File pFolder = new File(m_strPathToGeneData);
		if(!pFolder.exists())
		{
			Messagebox.show("ERROR: Invalid GTEX folder specified");
			return null;
		}

		String strFile = m_strPathToGeneData + "/GTEX_" + strEntrezID + ".tsv";
		File pFile = new File(strFile);
		
		if(!pFile.exists())
		{
			Messagebox.show("ERROR: No GTEX data available for gene: " + id);
			return null;
		}
		
		// store counts per condition <tissue, counts> 
		TreeMap<String, Vector<Double>> mapCountsToTissues = new TreeMap<String, Vector<Double>>(); 
		Scanner pScanner = null;
		
		try
		{
			pScanner = new Scanner(pFile);
		}
		catch(Exception e)
		{
			System.out.println("could not open file: " + strFile);
			e.printStackTrace();
		}
		
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
		
		int nWidthPerTissue = 20;
		if(nTissues > 0)
		{
			nWidthPerTissue = ((nMaxWidth-nLeftMargin) - 5*nTissues) / nTissues;	// use 20 pixel as spacer between tissues and 40 pixels for the border/axis of the graph
		}

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
		m_windowPopup.setParent(m_app);
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
		m_windowPopup.setHeight((nMaxHeight + 40) + "px");
		m_windowPopup.setWidth ((nMaxWidth  + 220) + "px");
		m_windowPopup.setSizable(false);

		regionPlot.setParent(layout);
		regionOptions.setParent(layout);
		
		Imagemap imgMap = new Imagemap();
		imgMap.setWidth(nMaxWidth+"px");
		imgMap.setHeight(nMaxHeight+"px");
		imgMap.setContent(img);
		imgMap.setParent(regionPlot);
	}
	
	public void ShowGTEXDataForGene(String strGene)
	{
		TreeMap<String, Vector<Double>> mapCountsPerTissue = GetGTEXDataForGene();
		if(mapCountsPerTissue != null)
		{
			DrawGTEXBoxPlot(mapCountsPerTissue);
		}
	}
	
	public void ShowGTEXDataForExon(String strExonicPart)
	{
		String pSplit[] = strExonicPart.split(" ");
		TreeMap<String, Vector<Double>> mapCountsPerTissue = m_app.GetCountsForGTEXExonicPart(pSplit[1]);
		if(mapCountsPerTissue != null)
		{
			DrawGTEXBoxPlot(mapCountsPerTissue);
		}
	}
}
