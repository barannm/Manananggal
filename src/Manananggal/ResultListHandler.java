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
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import org.zkoss.zk.ui.Component;
import org.zkoss.zk.ui.event.Event;
import org.zkoss.zk.ui.event.EventListener;
import org.zkoss.zk.ui.event.Events;
import org.zkoss.zk.ui.event.OpenEvent;
import org.zkoss.zul.Area;
import org.zkoss.zul.Bandbox;
import org.zkoss.zul.Bandpopup;
import org.zkoss.zul.Button;
import org.zkoss.zul.Cell;
import org.zkoss.zul.Checkbox;
import org.zkoss.zul.Column;
import org.zkoss.zul.Columns;
import org.zkoss.zul.Grid;
import org.zkoss.zul.Groupbox;
import org.zkoss.zul.Hlayout;
import org.zkoss.zul.Image;
import org.zkoss.zul.Imagemap;
import org.zkoss.zul.Label;
import org.zkoss.zul.Listbox;
import org.zkoss.zul.Listitem;
import org.zkoss.zul.Messagebox;
import org.zkoss.zul.Row;
import org.zkoss.zul.Rows;
import org.zkoss.zul.Tab;
import org.zkoss.zul.Tabpanel;
import org.zkoss.zul.Tabpanels;
import org.zkoss.zul.Tabs;
import org.zkoss.zul.Textbox;
import org.zkoss.zul.Vlayout;
import org.zkoss.zul.Window;

/**
 *    This class handles the permanent and temporary result lists displayed in the GUI.
 *    It manages the filters and allows for transfer of temporary results to the permanent
 *    result list.
*/
public class ResultListHandler
{
	static final int GRID_TYPE_PERMANENT	= 0;
	static final int GRID_TYPE_TEMPORARY	= 1;
	static final int GRID_TYPE_IMPORT		= 2;
	
	private Tabs 				m_Tabs;
	private Tabpanels 			m_Panels;
	
	private Grid				m_gridHitList;
	private Grid				m_gridTmpHitList;
	private Grid				m_gridImportedResults;
	private ResultFilterRule	m_resultFilterRulePermanentList;
	private ResultFilterRule	m_resultFilterRuleTemporaryList;
	private ResultFilterRule	m_resultFilterRuleImportedList;
	private MyFileUploadDlg		m_FileUploadDlg;
	
	private Vector<AnalysisResult> m_vcImportedResults;
	private AnalysisResultHandler m_resultHandler;
	
	SplicingWebApp				m_App;
	
	public ResultListHandler(Tabs tabs, Tabpanels panels, SplicingWebApp app)
	{
		m_Tabs			= tabs;
		m_Panels		= panels;
		m_App			= app;
		
		m_resultFilterRulePermanentList = new ResultFilterRule();
		m_resultFilterRuleTemporaryList = new ResultFilterRule();
		m_resultFilterRuleImportedList  = new ResultFilterRule();

		AddPermanentResultList();
		AddTemporaryResultList();
		AddImportTab();
		
		m_vcImportedResults = new Vector<AnalysisResult>();
		m_FileUploadDlg = new MyFileUploadDlg(m_App);
	}

	/** Clears all data in the result table */
	public void Clear()
	{
		if(m_gridHitList != null)
		{
			m_gridHitList.setMold("default");
			
			// clear old rows
			m_gridHitList.getChildren().clear();
			m_gridHitList.setMold("paging");
			m_gridHitList.setPageSize(10);
		}

		m_vcImportedResults.clear();
	}
	
	/**
	 *    Adds a tab for the permanent (stored) AS results and
	 *    a button to save the results to a file.
	 */
	private void AddPermanentResultList()
	{
		Tab tab = new Tab("Curated AS Results");
		tab.setParent(m_Tabs);
		
		Tabpanel panel = new Tabpanel();
		panel.setVflex("min");
		panel.setParent(m_Panels);
		
		Vlayout layoutV = new Vlayout();		
		layoutV.setParent(panel);
		
		m_gridHitList = new Grid();
		m_gridHitList.setId("hitlist");
		m_gridHitList.setParent(layoutV);
		m_gridHitList.setVflex("min");
		m_gridHitList.setHflex("min");
		m_gridHitList.setMold("paging");
		m_gridHitList.setPageSize(10);
		
		m_gridHitList.addEventListener("onPaging", new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				m_Tabs.invalidate(); // if we don't force a redraw of the tabs, the grid will be truncated
			}
			
		});
		
		m_gridHitList.addEventListener(Events.ON_SORT, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				m_Tabs.invalidate(); // if we don't force a redraw of the tabs, the grid will be truncated
			}
			
		});
		layoutV.setVflex("min");
		
		Button btnSaveHitList = new Button("Save changes");
		btnSaveHitList.setParent(layoutV);
		btnSaveHitList.setSclass("button orange");
		btnSaveHitList.setStyle("margin-left: 4px;");
		btnSaveHitList.setWidth("200px");
		btnSaveHitList.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				m_App.SaveSplicingResults();
				Messagebox.show("Changes saved.");
			}
		});
	}
	
	/**
	 *    Adds a tab for the temporary AS results and adds a
	 *    button to transfer temporary AS results to the permanent
	 *    AS result list.
	 */
	private void AddTemporaryResultList()
	{
		Tab tab = new Tab("Temporary AS Results");
		tab.setParent(m_Tabs);
		
		Tabpanel panel = new Tabpanel();
		panel.setVflex("min");
		panel.setParent(m_Panels);
		
		Vlayout layoutV = new Vlayout();
		layoutV.setParent(panel);
		
		m_gridTmpHitList = new Grid();
		m_gridTmpHitList.setId("tmpHitList");
		m_gridTmpHitList.setParent(layoutV);
		m_gridTmpHitList.setVflex("min");
		m_gridTmpHitList.setHflex("min");
		m_gridTmpHitList.setMold("paging");
		m_gridTmpHitList.setPageSize(10);
		
		m_gridTmpHitList.addEventListener("onPaging", new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				m_Tabs.invalidate(); // if we don't force a redraw of the tabs, the grid will be truncated
			}
			
		});

		layoutV.setVflex("min");
		
		Button btnSaveHitList = new Button("Add selection to curated results");
		btnSaveHitList.setParent(layoutV);
		btnSaveHitList.setSclass("button orange");
		btnSaveHitList.setStyle("margin-left: 4px;");
		btnSaveHitList.setWidth("200px");
		btnSaveHitList.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				TreeSet<AnalysisResult> vcTmpResults = m_resultHandler.GetAllTemporaryResults();
				TreeSet<AnalysisResult> vcResults 	 = new TreeSet<AnalysisResult>();
				
				// get event IDs
				List<Row> rows = m_gridTmpHitList.getRows().getChildren();
				
				for(Row row : rows)
				{
					Checkbox checkbox = (Checkbox)row.getFirstChild().getFirstChild();
					
					if(checkbox.isChecked())
					{
						int nResultID = (int)checkbox.getValue();
						
						for(AnalysisResult res : vcTmpResults)
						{
							if(res.GetID() == nResultID)
							{
								vcResults.add(res);
							}
						}
					}					
				}
				
				m_resultHandler.AddCuratedResults(vcResults);
				UpdatePermanentResultList(m_resultHandler);
			}
		});
	}

	/**
	 *    Adds a tab for the import of AS events generated
	 *    outside the web application (e.g. by Manananggal,
	 *    DEXSeq, rMATS, JSplice or Cuffdiff)
	 */
	private void AddImportTab()
	{
		Tab tab = new Tab("Imported Results");
		tab.setParent(m_Tabs);
		
		Tabpanel panel = new Tabpanel();
		panel.setVflex("min");
		panel.setHflex("min");
		panel.setParent(m_Panels);
		
		Vlayout layoutV = new Vlayout();
		layoutV.setParent(panel);
		
		m_gridImportedResults = new Grid();
		m_gridImportedResults.setId("importedHitList");
		m_gridImportedResults.setParent(layoutV);
		m_gridImportedResults.setVflex("min");
		m_gridImportedResults.setHflex("min");
		m_gridImportedResults.setMold("paging");
		m_gridImportedResults.setPageSize(10);
		
		m_gridImportedResults.addEventListener("onPaging", new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				m_Tabs.invalidate(); // if we don't force a redraw of the tabs, the grid will be truncated
			}
			
		});

		layoutV.setVflex("min");
		layoutV.setHflex("min");
		
		Button btnUpload = new Button("Import AS results");
		btnUpload.setParent(layoutV);
		btnUpload.setSclass("button orange");
		btnUpload.setStyle("margin-left: 4px;");
		btnUpload.setWidth("200px");
		
		btnUpload.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				if(!m_App.GetProjectModel().IsReady())
				{
					Messagebox.show("You need to select a project first");
					return;
				}
				
				m_FileUploadDlg.Show();
			}
			
		});
	}
	
	/**
	 *    Adds the menu for the splicing type filter to the corresponding
	 *    column in the result table.
	 */
	private void AddSplicingTypeFilter(Column parent, int nGridType, int nDataType, int nColIdx)
	{		
		Bandbox bandBox = new Bandbox();
		bandBox.setId("bbSplicingTypeFilter_" + nGridType + "_" + nDataType + "_" + nColIdx);
		bandBox.setStyle("position: absolute; right: 0px; margin-right: 10px");
		bandBox.setWidth("30px");
		bandBox.setParent(parent);
		
		Bandpopup bandPopup = new Bandpopup();
		bandPopup.setHflex("min");
		bandPopup.setParent(bandBox);
		
		Hlayout layout = new Hlayout();
		layout.setParent(bandPopup);
		
		Button btnSelectAll = new Button("Select All");
		btnSelectAll.setId("btnSelectAllASTypes_" + nGridType + "_" + nDataType + "_" + nColIdx);
		btnSelectAll.setParent(layout);
		btnSelectAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				for(int i=0; i<13; i++)
				{
					Button btn = (Button)event.getTarget();
					int nGridType = Integer.parseInt(btn.getId().split("_")[1]);
					int nDataType = Integer.parseInt(btn.getId().split("_")[2]);
					
					switch(nGridType)
					{
						case GRID_TYPE_PERMANENT:
						{
							m_resultFilterRulePermanentList.ShowASType(i);
							UpdatePermanentResultList(null);
							break;
						}
						
						case GRID_TYPE_TEMPORARY:
						{
							m_resultFilterRuleTemporaryList.ShowASType(i);
							UpdateTemporaryResultList(null);
							break;
						}
						
						case GRID_TYPE_IMPORT:
						{
							m_resultFilterRuleImportedList.ShowASType(i);
							UpdateImportedResultList(null, nDataType, false);
							break;
						}
					}
				}
			}
		});
		
		Button btnHideAll = new Button("Unselect All");
		btnHideAll.setId("btnUnselectAllASTypes_" + nGridType + "_" + nDataType + "_" + nColIdx);
		btnHideAll.setParent(layout);
		btnHideAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				for(int i=0; i<13; i++)
				{
					Button btn = (Button)event.getTarget();
					int nGridType = Integer.parseInt(btn.getId().split("_")[1]);
					int nDataType = Integer.parseInt(btn.getId().split("_")[2]);
					
					switch(nGridType)
					{
						case GRID_TYPE_PERMANENT:
						{
							m_resultFilterRulePermanentList.HideASType(i);
							UpdatePermanentResultList(null);
							break;
						}
						
						case GRID_TYPE_TEMPORARY:
						{
							m_resultFilterRuleTemporaryList.HideASType(i);
							UpdateTemporaryResultList(null);
							break;
						}
						
						case GRID_TYPE_IMPORT:
						{
							m_resultFilterRuleImportedList.HideASType(i);
							UpdateImportedResultList(null, nDataType, false);
							break;
						}
					}
				}
			}
		});
		
		Listbox listBox = new Listbox();
		listBox.setWidth("300px");
		listBox.setMultiple(true);
		listBox.setCheckmark(true);
		listBox.setParent(bandPopup);
		
		bandBox.addEventListener(Events.ON_OPEN, new EventListener<OpenEvent>()
		{
			public void onEvent(OpenEvent event) throws Exception
			{
				if(event.isOpen())
					return;
				
				Bandbox bandbox = (Bandbox)event.getTarget();
				Bandpopup popup = (Bandpopup)bandbox.getFirstChild();
				
				
				int nGridType = Integer.parseInt(bandbox.getId().split("_")[1]);
				int nDataType = Integer.parseInt(bandbox.getId().split("_")[2]);
				
				ResultFilterRule resultFilterRule = null;
				switch(nGridType)
				{
					case GRID_TYPE_PERMANENT: resultFilterRule = m_resultFilterRulePermanentList; break;
					case GRID_TYPE_TEMPORARY: resultFilterRule = m_resultFilterRuleTemporaryList; break;
					case GRID_TYPE_IMPORT:    resultFilterRule = m_resultFilterRuleImportedList;  break;
				}
				
				for(Component c : popup.getChildren())
				{
					if(c instanceof Listbox)
					{
						Listbox listbox = (Listbox)c;
						
						for(Listitem item : listbox.getItems())
						{					
							if(item.isSelected())
								resultFilterRule.ShowASType((int)item.getValue());
							else
								resultFilterRule.HideASType((int)item.getValue());
						}
						
						switch(nGridType)
						{
							case GRID_TYPE_PERMANENT: UpdatePermanentResultList(null); break;
							case GRID_TYPE_TEMPORARY: UpdateTemporaryResultList(null); break;
							case GRID_TYPE_IMPORT:    UpdateImportedResultList(null, nDataType, false);  break;
						}

						return;
					}
				}
			}
		});
		
		ResultFilterRule resultFilterRule = null;
		switch(nGridType)
		{
			case GRID_TYPE_PERMANENT: resultFilterRule = m_resultFilterRulePermanentList; break;
			case GRID_TYPE_TEMPORARY: resultFilterRule = m_resultFilterRuleTemporaryList; break;
			case GRID_TYPE_IMPORT:    resultFilterRule = m_resultFilterRuleImportedList;  break;
		}
		
		Listitem item = new Listitem("Exon skipping");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_EXON_SKIPPING);
		if(resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_EXON_SKIPPING))
			item.setSelected(true);
		
		item = new Listitem("retained introns");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_RETAINED_INTRON);
		if(resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_RETAINED_INTRON))
			item.setSelected(true);
		
		item = new Listitem("alt. 5' exon end");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_5_PRIME_EXON_END);
		if(resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_5_PRIME_EXON_END))
			item.setSelected(true);
				
		item = new Listitem("alt. 3' exon end");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_3_PRIME_EXON_END);
		if(resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_3_PRIME_EXON_END))
			item.setSelected(true);
		
		item = new Listitem("alt. start (unique junction; double)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN_DOUBLE);
		if(resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN_DOUBLE))
			item.setSelected(true);
		
		item = new Listitem("alt. end (unique junction; double)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN_DOUBLE);
		if(resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN_DOUBLE))
			item.setSelected(true);
		
		item = new Listitem("alt. start (shared junction; double)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN_DOUBLE);
		if(resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN_DOUBLE))
			item.setSelected(true);
		
		item = new Listitem("alt. end (shared junction; double)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN_DOUBLE);
		if(resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN_DOUBLE))
			item.setSelected(true);
		
		item = new Listitem("alt. start (unique junction)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN);
		if(resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN))
			item.setSelected(true);
		
		item = new Listitem("alt. end (unique junction)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN);
		if(resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN))
			item.setSelected(true);
		
		item = new Listitem("alt. start (shared junction)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN);
		if(resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN))
			item.setSelected(true);
		
		item = new Listitem("alt. end (shared junction)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN);
		if(resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN))
			item.setSelected(true);
	}

	/**
	 *    Adds a menu for the result type filter to the corresponding
	 *    column in the result table.
	 */
	private void AddResultTypeFilter(Column parent, int nGridType, int nDataType, int nColIdx)
	{		
		Bandbox bandBox = new Bandbox();
		bandBox.setId("bbResultTypeFilter_" + nGridType + "_" + nDataType + "_" + nColIdx);
		bandBox.setStyle("position: absolute; right: 0px; margin-right: 10px");
		bandBox.setWidth("30px");
		bandBox.setParent(parent);
		
		Bandpopup bandPopup = new Bandpopup();
		bandPopup.setHflex("min");
		bandPopup.setParent(bandBox);
		
		Hlayout layout = new Hlayout();
		layout.setParent(bandPopup);
		
		Button btnSelectAll = new Button("Select All");
		btnSelectAll.setId("btnSelectAllResultTypes_" + nGridType + "_" + nDataType + "_" + nColIdx);
		btnSelectAll.setParent(layout);
		btnSelectAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{				
				for(int i=0; i<4; i++)
				{
					Button btn = (Button)event.getTarget();
					int nGridType = Integer.parseInt(btn.getId().split("_")[1]);
					int nDataType = Integer.parseInt(btn.getId().split("_")[2]);
					
					switch(nGridType)
					{
						case GRID_TYPE_PERMANENT:
						{
							m_resultFilterRulePermanentList.ShowResultType(i);
							UpdatePermanentResultList(null);
							break;
						}
						
						case GRID_TYPE_TEMPORARY:
						{
							m_resultFilterRuleTemporaryList.ShowResultType(i);
							UpdateTemporaryResultList(null);
							break;
						}
						
						case GRID_TYPE_IMPORT:
						{
							m_resultFilterRuleImportedList.ShowResultType(i);
							UpdateImportedResultList(null, nDataType, false);
							break;
						}
					}
				}
			}
		});
		
		Button btnHideAll = new Button("Unselect All");
		btnHideAll.setParent(layout);
		btnHideAll.setId("btnUnselectAllResultTypes_" + nGridType + "_" + nDataType + "_" + nColIdx);
		btnHideAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				for(int i=0; i<4; i++)
				{
					Button btn = (Button)event.getTarget();
					int nGridType = Integer.parseInt(btn.getId().split("_")[1]);
					int nDataType = Integer.parseInt(btn.getId().split("_")[2]);
					
					switch(nGridType)
					{
						case GRID_TYPE_PERMANENT:
						{
							m_resultFilterRulePermanentList.HideResultType(i);
							UpdatePermanentResultList(null);
							break;
						}
						
						case GRID_TYPE_TEMPORARY:
						{
							m_resultFilterRuleTemporaryList.HideResultType(i);
							UpdateTemporaryResultList(null);
							break;
						}
						
						case GRID_TYPE_IMPORT:
						{
							m_resultFilterRuleImportedList.HideResultType(i);
							UpdateImportedResultList(null, nDataType, false);
							break;
						}
					}
				}
			}
		});
		
		Listbox listBox = new Listbox();
		listBox.setWidth("300px");
		listBox.setMultiple(true);
		listBox.setCheckmark(true);
		listBox.setParent(bandPopup);
		
		bandBox.addEventListener(Events.ON_OPEN, new EventListener<OpenEvent>()
		{
			public void onEvent(OpenEvent event) throws Exception
			{
				if(event.isOpen())
					return;
				
				Bandbox bandbox = (Bandbox)event.getTarget();
				Bandpopup popup = (Bandpopup)bandbox.getFirstChild();
				
				int nGridType = Integer.parseInt(bandbox.getId().split("_")[1]);
				int nDataType = Integer.parseInt(bandbox.getId().split("_")[2]);
				
				ResultFilterRule resultFilterRule = null;
				switch(nGridType)
				{
					case GRID_TYPE_PERMANENT: resultFilterRule = m_resultFilterRulePermanentList; break;
					case GRID_TYPE_TEMPORARY: resultFilterRule = m_resultFilterRuleTemporaryList; break;
					case GRID_TYPE_IMPORT:    resultFilterRule = m_resultFilterRuleImportedList;  break;
				}
				
				for(Component c : popup.getChildren())
				{
					if(c instanceof Listbox)
					{
						Listbox listbox = (Listbox)c;
						
						for(Listitem item : listbox.getItems())
						{					
							if(item.isSelected())
								resultFilterRule.ShowResultType((int)item.getValue());
							else
								resultFilterRule.HideResultType((int)item.getValue());
						}
						
						switch(nGridType)
						{
							case GRID_TYPE_PERMANENT: UpdatePermanentResultList(null); break;
							case GRID_TYPE_TEMPORARY: UpdateTemporaryResultList(null); break;
							case GRID_TYPE_IMPORT:    UpdateImportedResultList(null, nDataType, false);  break;
						}

						return;
					}
				}
			}
		});
		
		ResultFilterRule resultFilterRule = null;
		switch(nGridType)
		{
			case GRID_TYPE_PERMANENT: resultFilterRule = m_resultFilterRulePermanentList; break;
			case GRID_TYPE_TEMPORARY: resultFilterRule = m_resultFilterRuleTemporaryList; break;
			case GRID_TYPE_IMPORT:    resultFilterRule = m_resultFilterRuleImportedList;  break;
		}
		
		Listitem item = new Listitem("combined");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(AnalysisResult.RESULT_EVIDENCE_TYPE_COMBINED);
		if(resultFilterRule.IsShowingResultType(AnalysisResult.RESULT_EVIDENCE_TYPE_COMBINED))
			item.setSelected(true);
		
		item = new Listitem("ratio only");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(AnalysisResult.RESULT_EVIDENCE_TYPE_RATIO_ONLY);
		if(resultFilterRule.IsShowingResultType(AnalysisResult.RESULT_EVIDENCE_TYPE_RATIO_ONLY))
			item.setSelected(true);
		
		item = new Listitem("split read only");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(AnalysisResult.RESULT_EVIDENCE_TYPE_SPLIT_READ_ONLY);
		if(resultFilterRule.IsShowingResultType(AnalysisResult.RESULT_EVIDENCE_TYPE_SPLIT_READ_ONLY))
			item.setSelected(true);
	}

	/**
	 *    Adds a menu to filter for strings in a given column of the
	 *    result table. Because multiple string filters exist, each
	 *    requires a unique ID that is specified by the filter ID.
	 *    See ResultFilterRule for valid filter IDs.
	 */
	private void AddStringFilter(Column parent, int nFilterID, int nGridType, int nDataType, int nColIdx)
	{
		ResultFilterRule resultFilterRule = null;
		switch(nGridType)
		{
			case GRID_TYPE_PERMANENT: resultFilterRule = m_resultFilterRulePermanentList; break;
			case GRID_TYPE_TEMPORARY: resultFilterRule = m_resultFilterRuleTemporaryList; break;
			case GRID_TYPE_IMPORT:    resultFilterRule = m_resultFilterRuleImportedList;  break;
		}
		
		TreeSet<String> vcStrings = resultFilterRule.GetTextFilter(nFilterID).GetStrings();
		
		Bandbox bandBox = new Bandbox();
		bandBox.setStyle("position: absolute; right: 0px; margin-right: 10px");
		bandBox.setId("bbStringFilter_" + nGridType + "_" + nDataType + "_" + nColIdx);
		bandBox.setWidth("30px");
		bandBox.setParent(parent);
		bandBox.setValue(""+nFilterID);
		
		Bandpopup bandPopup = new Bandpopup();
		bandPopup.setHflex("min");
		bandPopup.setParent(bandBox);
		
		Hlayout layout = new Hlayout();
		layout.setParent(bandPopup);
		
		Button btnSelectAll = new Button("Select All");
		btnSelectAll.setId("btnStringFilter_" + nGridType + "_" + nDataType + "_" + nColIdx);
		btnSelectAll.setParent(layout);
		btnSelectAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Bandbox bandbox = (Bandbox)event.getTarget().getParent().getParent().getParent();
				int nFilterID = Integer.parseInt(bandbox.getValue());
				
				Button btn = (Button)event.getTarget();
				int nGridType = Integer.parseInt(btn.getId().split("_")[1]);
				int nDataType = Integer.parseInt(btn.getId().split("_")[2]);
				
				switch(nGridType)
				{
					case GRID_TYPE_PERMANENT:
					{
						m_resultFilterRulePermanentList.GetTextFilter(nFilterID).SelectAll();
						UpdatePermanentResultList(null);
						break;
					}
					
					case GRID_TYPE_TEMPORARY:
					{
						m_resultFilterRuleTemporaryList.GetTextFilter(nFilterID).SelectAll();
						UpdateTemporaryResultList(null);
						break;
					}
					
					case GRID_TYPE_IMPORT:
					{
						m_resultFilterRuleImportedList.GetTextFilter(nFilterID).SelectAll();
						UpdateImportedResultList(null, nDataType, false);
						break;
					}
				}
			}
		});
		
		Button btnHideAll = new Button("Unselect All");
		btnHideAll.setParent(layout);
		btnHideAll.setId("btnStringFilterUnselectAll_" + nGridType + "_" + nDataType + "_" + nColIdx);
		btnHideAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Bandbox bandbox = (Bandbox)event.getTarget().getParent().getParent().getParent();
				int nFilterID = Integer.parseInt(bandbox.getValue());
				
				Button btn = (Button)event.getTarget();
				int nGridType = Integer.parseInt(btn.getId().split("_")[1]);
				int nDataType = Integer.parseInt(btn.getId().split("_")[2]);
				
				switch(nGridType)
				{
					case GRID_TYPE_PERMANENT:
					{
						m_resultFilterRulePermanentList.GetTextFilter(nFilterID).UnselectAll();
						UpdatePermanentResultList(null);
						break;
					}
					
					case GRID_TYPE_TEMPORARY:
					{
						m_resultFilterRuleTemporaryList.GetTextFilter(nFilterID).UnselectAll();
						UpdateTemporaryResultList(null);
						break;
					}
					
					case GRID_TYPE_IMPORT:
					{
						m_resultFilterRuleImportedList.GetTextFilter(nFilterID).UnselectAll();
						UpdateImportedResultList(null, nDataType, false);
						break;
					}
				}
			}
		});
		
		Listbox listBox = new Listbox();
		listBox.setWidth("300px");
		listBox.setMultiple(true);
		listBox.setCheckmark(true);
		listBox.setParent(bandPopup);
		
		bandBox.addEventListener(Events.ON_OPEN, new EventListener<OpenEvent>()
		{
			public void onEvent(OpenEvent event) throws Exception
			{
				if(event.isOpen())
					return;
				
				Bandbox bandbox = (Bandbox)event.getTarget();
				Bandpopup popup = (Bandpopup)bandbox.getFirstChild();
				int nFilterID = Integer.parseInt(bandbox.getValue());
				
				ResultFilterRule resultFilterRule = null;
				int nGridType = Integer.parseInt(bandbox.getId().split("_")[1]);
				int nDataType = Integer.parseInt(bandbox.getId().split("_")[2]);
				switch(nGridType)
				{
					case GRID_TYPE_PERMANENT: resultFilterRule = m_resultFilterRulePermanentList; break;
					case GRID_TYPE_TEMPORARY: resultFilterRule = m_resultFilterRuleTemporaryList; break;
					case GRID_TYPE_IMPORT:    resultFilterRule = m_resultFilterRuleImportedList;  break;
				}
				
				resultFilterRule.GetTextFilter(nFilterID).SelectAll();
				
				for(Component c : popup.getChildren())
				{
					if(c instanceof Listbox)
					{
						Listbox listbox = (Listbox)c;
						
						for(Listitem item : listbox.getItems())
						{
							if(item.isSelected())
								resultFilterRule.GetTextFilter(nFilterID).SelectString((String)item.getValue());
							else
								resultFilterRule.GetTextFilter(nFilterID).UnselectString((String)item.getValue());
						}
						
						switch(nGridType)
						{
							case GRID_TYPE_PERMANENT: UpdatePermanentResultList(null); break;
							case GRID_TYPE_TEMPORARY: UpdateTemporaryResultList(null); break;
							case GRID_TYPE_IMPORT:    UpdateImportedResultList(null, nDataType, false); break;
						}

						return;
					}
				}
			}
		});
		
		for(String strString : vcStrings)
		{
			Listitem item = new Listitem(strString);
			item.setParent(listBox);
			item.setCheckable(true);
			item.setValue(strString);
			
			if(resultFilterRule.GetTextFilter(nFilterID).IsStringSelected(strString))
				item.setSelected(true);
		}
	}

	/**
	 *    Adds a menu to filter for values in a given column of the
	 *    result table. Because multiple number filters exist, each
	 *    requires a unique ID that is specified by the filter ID.
	 *    See ResultFilterRule for valid filter IDs.
	 */
	private void AddNumberFilter(Component parent, int nFilterID, int nGridType, int nDataType, int nColIdx)
	{
		ResultFilterRule resultFilterRule = null;
		switch(nGridType)
		{
			case GRID_TYPE_PERMANENT: resultFilterRule = m_resultFilterRulePermanentList; break;
			case GRID_TYPE_TEMPORARY: resultFilterRule = m_resultFilterRuleTemporaryList; break;
			case GRID_TYPE_IMPORT:    resultFilterRule = m_resultFilterRuleImportedList;  break;
		}
		
		Bandbox bandBox = new Bandbox();
		bandBox.setStyle("position: absolute; right: 0px; margin-right: 10px");
		bandBox.setWidth("30px");
		bandBox.setParent(parent);
		bandBox.setValue(""+nFilterID);
		
		bandBox.addEventListener(Events.ON_BLUR, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Bandbox bandBox = (Bandbox)event.getTarget();
				bandBox.close();
			}
		});
		
		Bandpopup bandPopup = new Bandpopup();
		bandPopup.setHflex("min");
		bandPopup.setParent(bandBox);
		
		Vlayout layoutV = new Vlayout();
		layoutV.setParent(bandPopup);
		
		Hlayout layout = new Hlayout();
		layout.setParent(layoutV);
		
		Label labMinValue = new Label("min value: ");
		labMinValue.setParent(layout);
		
		Textbox boxMinValue = new Textbox("" + resultFilterRule.GetNumberFilter(nFilterID).GetMinValue());
		boxMinValue.setConstraint("/[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?/: Only numbers are allowed, including decimal numbers and scientific writing, e.g. 5.0E-3");
		boxMinValue.setId("tbNumberFilterMin_" + nGridType + "_" + nDataType + "_" + nColIdx);
		boxMinValue.setParent(layout);
		boxMinValue.addEventListener(Events.ON_CHANGE, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Textbox textbox = (Textbox)event.getTarget();
				Bandbox bandbox = (Bandbox)textbox.getParent().getParent().getParent().getParent();
				int nFilterID = Integer.parseInt(bandbox.getValue());
				
				float fValue = Float.parseFloat(textbox.getText());
				
				int nGridType = Integer.parseInt(textbox.getId().split("_")[1]);
				int nDataType = Integer.parseInt(textbox.getId().split("_")[2]);
				switch(nGridType)
				{
					case GRID_TYPE_PERMANENT:
					{
						m_resultFilterRulePermanentList.GetNumberFilter(nFilterID).SetMinValue(fValue);
						UpdatePermanentResultList(null);
						break;
					}
					case GRID_TYPE_TEMPORARY:
					{
						m_resultFilterRuleTemporaryList.GetNumberFilter(nFilterID).SetMinValue(fValue);
						UpdateTemporaryResultList(null);
						break;
					}
					case GRID_TYPE_IMPORT:
					{
						m_resultFilterRuleImportedList.GetNumberFilter(nFilterID).SetMinValue(fValue);
						UpdateImportedResultList(null, nDataType, false);
						break;
					}
				}
			}
		});
		
		layout = new Hlayout();
		layout.setParent(layoutV);
		
		Label labMaxValue = new Label("max value: ");
		labMaxValue.setParent(layout);
		
		Textbox boxMaxValue = new Textbox("" + resultFilterRule.GetNumberFilter(nFilterID).GetMaxValue());
		boxMaxValue.setConstraint("/[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?/: Only numbers are allowed, including decimal numbers and scientific writing, e.g. 5.0E-3");
		boxMaxValue.setId("tbNumberFilterMax_" + nGridType + "_" + nDataType + "_" + nColIdx);
		boxMaxValue.setParent(layout);
		boxMaxValue.addEventListener(Events.ON_CHANGE, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Textbox textbox = (Textbox)event.getTarget();
				Bandbox bandbox = (Bandbox)textbox.getParent().getParent().getParent().getParent();
				int nFilterID = Integer.parseInt(bandbox.getValue());
				
				float fValue = Float.parseFloat(textbox.getText());
				
				int nGridType = Integer.parseInt(textbox.getId().split("_")[1]);
				int nDataType = Integer.parseInt(textbox.getId().split("_")[2]);
				switch(nGridType)
				{
					case GRID_TYPE_PERMANENT:
					{
						m_resultFilterRulePermanentList.GetNumberFilter(nFilterID).SetMaxValue(fValue);
						UpdatePermanentResultList(null);
						break;
					}
					case GRID_TYPE_TEMPORARY:
					{
						m_resultFilterRuleTemporaryList.GetNumberFilter(nFilterID).SetMaxValue(fValue);
						UpdateTemporaryResultList(null);
						break;
					}
					case GRID_TYPE_IMPORT:
					{
						m_resultFilterRuleImportedList.GetNumberFilter(nFilterID).SetMaxValue(fValue);
						UpdateImportedResultList(null, nDataType, false);
						break;
					}
				}
			}
		});
		
		Button btnDisable = new Button("Disable filter");
		btnDisable.setParent(layoutV);
		btnDisable.setId("btnDisableFilter_" + nGridType + "_" + nDataType + "_" + nColIdx);
		btnDisable.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Bandbox bandbox = (Bandbox)event.getTarget().getParent().getParent().getParent();
				int nFilterID = Integer.parseInt(bandbox.getValue());
				
				Button btnTarget = (Button)event.getTarget();
				int nGridType = Integer.parseInt(btnTarget.getId().split("_")[1]);
				int nDataType = Integer.parseInt(btnTarget.getId().split("_")[2]);
				switch(nGridType)
				{
					case GRID_TYPE_PERMANENT:
					{
						m_resultFilterRulePermanentList.GetNumberFilter(nFilterID).Disable();
						UpdatePermanentResultList(null);
						break;
					}
					case GRID_TYPE_TEMPORARY:
					{
						m_resultFilterRuleTemporaryList.GetNumberFilter(nFilterID).Disable();
						UpdateTemporaryResultList(null);
						break;
					}
					case GRID_TYPE_IMPORT:
					{
						m_resultFilterRuleImportedList.GetNumberFilter(nFilterID).Disable();
						UpdateImportedResultList(null, nDataType, false);
						break;
					}
				}
			}
		});
	}
	
	/**
	 *    Updates the contents of the permanent result list.
	 *    This function adds the table to the permanent result list tab
	 *    and fills it with data.
	 *    It is invoked when a new project is opened or if the visible 
	 *    contents of the result table changes (e.g. by changing the
	 *    filter rules).
	 */
	public void UpdatePermanentResultList(AnalysisResultHandler resultHandler) throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		// update results
		boolean bInit = false;
		if(resultHandler != null)
		{
			m_resultHandler = resultHandler;
			bInit = true;
		}
		
		NumberFormat nf = NumberFormat.getNumberInstance(Locale.ENGLISH);
		DecimalFormat format = (DecimalFormat)nf;
		format.applyPattern("0.####E0");
		
		m_gridHitList.setMold("default");
		
		// clear old rows
		m_gridHitList.getChildren().clear();
		m_gridHitList.setMold("paging");
		m_gridHitList.setPageSize(10);
		
		Column colPSI = AddColumnHeadersToGridForManananggalResults(m_gridHitList, GRID_TYPE_PERMANENT, 0, bInit);
		
		Rows rows = new Rows();
		rows.setParent(m_gridHitList);
		
		//##################################
		//             add rows
		//##################################
		TreeSet<AnalysisResult> vcResults = m_resultHandler.GetAllResults();
		
		for(AnalysisResult res : vcResults)
		{
			if(!m_resultFilterRulePermanentList.bIsValidResult(res))
				continue;
			
			Row row = new Row();
			row.setId("ID_" + res.GetID());
			row.setValue(res.GetRating());
			row.setParent(rows);
			
			AddCellToRow(row, "", "remove_", "30px", "float: center", res, 1);
			AddCellToRow(row, "", "view_", "40px", "float: center", res, 2);
			AddCellToRow(row, "", "detailView_", "40px", "float: center", res, 3);
			AddCellToRow(row, "", "details_", "40px", "float: center", res, 4);
			AddCellToRow(row, "", "rating_", "110px", "float: center", res, 5);
			
			AddCellToRow(row, res.GetGeneID(),             "", "140px", "float: center", res, 0);
			AddCellToRow(row, res.GetGeneSymbol(),         "", "140px", "float: center", res, 0);
			AddCellToRow(row, res.GetConditionA(),         "", "160px", "float: center", res, 0);
			AddCellToRow(row, res.GetConditionB(),         "", "160px", "float: center", res, 0);
			AddCellToRow(row, res.GetASTypeAsString(),       "", "170px", "float: center", res, 0);
			AddCellToRow(row, res.GetResultTypeAsString(), "", "120px", "float: center", res, 0);
			
			if(res.HasAltExonA())
			{
				String strText = String.format(Locale.ENGLISH, "%.2f%%", res.GetAbsoluteChangeA()*100);
				if(Math.abs(res.GetAbsoluteChangeA()) < 0.05)
					AddCellToRow(row, strText, "", "210px", "float: right; color: red;", res, 0);
				else
					AddCellToRow(row, strText, "", "210px", "float: right", res, 0);
			}
			else
			{
				AddCellToRow(row, "NA", "", "210px", "float: right", res, 0);
			}
			
			if(res.HasAltExonB())
			{
				String strText = String.format(Locale.ENGLISH, "%.2f%%", res.GetAbsoluteChangeB()*100);
				if(Math.abs(res.GetAbsoluteChangeB()) < 0.05)
					AddCellToRow(row, strText, "", "210px", "float: right; color: red;", res, 0);
				else
					AddCellToRow(row, strText, "", "210px", "float: right", res, 0);
			}
			else
			{
				AddCellToRow(row, "NA", "", "210px", "float: right", res, 0);
			}
			
			if(res.HasPSIScore())
			{
				String strText = String.format(Locale.ENGLISH, "%.2f%%", res.GetInclusionChange()*100);
				if(Math.abs(res.GetInclusionChange()) < 0.05)
					AddCellToRow(row, strText, "", "150px", "float: right; color: red;", res, 0);
				else
					AddCellToRow(row, strText, "", "150px", "float: right", res, 0);
			}
			else
			{
				AddCellToRow(row, "NA", "", "150px", "float: right", res, 0);
			}
			
			if(res.HasAltExonA())
			{
				AddCellToRow(row, format.format(res.GetPValueA()), "", "210px", "float: right", res, 0);	
			}
			else
			{
				AddCellToRow(row, "NA", "", "210px", "float: right", res, 0);
			}
	
			if(res.HasAltExonB())
			{
				AddCellToRow(row, format.format(res.GetPValueB()), "", "210px", "float: right", res, 0);	
			}
			else
			{
				AddCellToRow(row, "NA", "", "210px", "float: right", res, 0);
			}
	
			if(res.HasPSIScore())
			{
				AddCellToRow(row, format.format(res.GetPValuePSI()), "", "150px", "float: right", res, 0);
				AddCellToRow(row, String.format(Locale.ENGLISH, "%s", res.HasSignificantPSIScore()), "", "150px", "float: right", res, 0);
				AddCellToRow(row, String.format(Locale.ENGLISH, "%s", res.UsesNovelJunction()), "", "170px", "float: right", res, 0);
			}
			else
			{
				AddCellToRow(row, "NA", "", "150px", "float: right", res, 0);
				AddCellToRow(row, "NA", "", "150px", "float: center", res, 0);			
				AddCellToRow(row, "NA", "", "170px", "float: center", res, 0);
			}
			
			int nEventStart = Integer.MAX_VALUE;
			int nEventEnd	= Integer.MIN_VALUE;
			if(res.HasAltExonA())
			{
				nEventStart = Math.min(res.GetStartA(), nEventStart);
				nEventEnd 	= Math.max(res.GetEndA(),   nEventEnd);
			}
			if(res.HasAltExonB())
			{
				nEventStart = Math.min(res.GetStartB(), nEventStart);
				nEventEnd	= Math.max(res.GetEndB(),   nEventEnd);
			}
			if(res.HasPSIScore())
			{
				nEventStart = Math.min(res.GetInclusionJunctionStart(), nEventStart);
				nEventEnd 	= Math.max(res.GetInclusionJunctionEnd(),   nEventEnd);
				
				nEventStart = Math.min(res.GetExclusionJunctionStart(), nEventStart);
				nEventEnd 	= Math.max(res.GetExclusionJunctionEnd(),   nEventEnd);
			}

			AddCellToRow(row, String.format(Locale.ENGLISH, "%s:%,d - %,d", res.GetReferenceName(), nEventStart, nEventEnd), "", "200px", "float: left", res, 0);
			AddCellToRow(row, res.GetComment(), "commentField_", "400px", "float: left", res, 6);
		}

		colPSI.sort(true);
		
		m_Tabs.invalidate();
	}
	
	/**
	 *    This function is invoked by clicking the 'detail icon'
	 *    located in front of each result. It opens a popup window
	 *    that shows some detailed information on the event.
	 */
	public void ShowDetailPopup(AnalysisResult res)
	{		
		NumberFormat nf = NumberFormat.getNumberInstance(Locale.ENGLISH);
		DecimalFormat format = (DecimalFormat)nf;
		format.applyPattern("0.####E0");
		
		int nBoxWidth = 505;
		
		Window windowPopup = new Window();
		windowPopup.setWidth("544px");
		windowPopup.setHeight("600px");
		windowPopup.setParent(m_App);
		windowPopup.setTitle("Details");
		windowPopup.setSizable(true);
		windowPopup.setClosable(true);
		windowPopup.setMaximizable(true);
		windowPopup.setBorder(true);
		windowPopup.setPosition("center,center");
		windowPopup.setVisible(true);
		windowPopup.doPopup();
		windowPopup.setTopmost();
		windowPopup.setContentStyle("overflow:auto");
		
		Vlayout layout = new Vlayout();
		layout.setParent(windowPopup);
		
		Groupbox grpBox = new Groupbox();
		grpBox.setTitle("General");
		grpBox.setMold("3d");
		grpBox.setParent(layout);
		grpBox.setWidth((nBoxWidth + 10) + "px");
		grpBox.setOpen(false);
		
		Grid grid = new Grid();
		grid.setParent(grpBox);
		grid.setWidth(nBoxWidth + "px");
		
		//##################################
		//           add header
		//##################################
		Columns cols = new Columns();
		cols.setParent(grid);
		//	cols.setMenupopup("auto");
		
		Column col = new Column("");
		col.setWidth("300px");
		col.setParent(cols);
		
		col = new Column("value");
		col.setWidth("200px");
		col.setParent(cols);
		
		//##################################
		//             add rows
		//##################################
		Rows rows = new Rows();
		rows.setParent(grid);

		Row row = new Row();
		row.setParent(rows);
		Label label = new Label("Condition A:");
		label.setParent(row);
		label = new Label(String.format(Locale.ENGLISH, "%s", res.GetConditionA()));
		label.setParent(row);
		
		row = new Row();
		row.setParent(rows);
		label = new Label("Condition B:");
		label.setParent(row);
		label = new Label(String.format(Locale.ENGLISH, "%s", res.GetConditionB()));
		label.setParent(row);

		for(int i=0; i<2; i++)
		{
			AlternativeSplicingExon ex = null;
			
			if(i == 0 && res.HasAltExonA())
				ex = res.GetAltExonA();
			else if(res.HasAltExonB())
				ex = res.GetAltExonB();
			
			if(ex != null)
			{
				grpBox = new Groupbox();
				if(i == 0)
					grpBox.setTitle("1st alternative exon");
				else
					grpBox.setTitle("2nd alternative exon");
				
				grpBox.setMold("3d");
				grpBox.setParent(layout);
				grpBox.setWidth((nBoxWidth + 10) + "px");
				grpBox.setOpen(false);
				
				grid = new Grid();
				grid.setParent(grpBox);
				grid.setWidth(nBoxWidth + "px");
				
				//##################################
				//           add header
				//##################################
				cols = new Columns();
				cols.setParent(grid);
				//	cols.setMenupopup("auto");
				
				col = new Column("");
				col.setWidth("300px");
				col.setParent(cols);
				
				col = new Column("value");
				col.setWidth("200px");
				col.setParent(cols);
				
				//##################################
				//             add rows
				//##################################
				rows = new Rows();
				rows.setParent(grid);
				
				row = new Row();
				row.setParent(rows);
				label = new Label("alternative exon:");
				label.setParent(row);
				label = new Label(String.format(Locale.ENGLISH, "%d - %d", ex.m_nStart, ex.m_nEnd));
				label.setParent(row);
				
				row = new Row();
				row.setParent(rows);
				label = new Label("absolute ratio change:");
				label.setParent(row);
				label = new Label(String.format(Locale.ENGLISH, "%.2f", ex.m_fFractionChangeAbsolute));
				label.setParent(row);
				
				row = new Row();
				row.setParent(rows);
				label = new Label("relative ratio change:");
				label.setParent(row);
				label = new Label(String.format(Locale.ENGLISH, "%.2f", ex.m_fFractionChangeRelative));
				label.setParent(row);
				
				row = new Row();
				row.setParent(rows);
				label = new Label("p-value: ");
				label.setParent(row);
				label = new Label(format.format(ex.m_fPValue));
				label.setParent(row);
				
				row = new Row();
				row.setParent(rows);
				label = new Label("cov. per base in " + ex.m_strConditionA + ":");
				label.setParent(row);
				label = new Label(String.format(Locale.ENGLISH, "%.4f", ex.m_fAltExonCovPerBaseA));
				label.setParent(row);
				
				row = new Row();
				row.setParent(rows);
				label = new Label("cov. per base in " + ex.m_strConditionB + ":");
				label.setParent(row);
				label = new Label(String.format(Locale.ENGLISH, "%.4f", ex.m_fAltExonCovPerBaseB));
				label.setParent(row);
				
				// add boxplot			
				TreeMap<String, double[]> mapRatios = new TreeMap<String, double[]>();
				mapRatios.put("tested_exon", ex.m_pFractionTestedExon);
				mapRatios.put("other_exons", ex.m_pFractionOtherExons);
				
				// settings for the box plot
				int nBarWidth 	= 80;
				int nSpacer   	= 20;
				int nTotalBarDistance = nBarWidth+nSpacer;
				int nPlotWidth = nTotalBarDistance*2+40;
				
				BufferedImage img = new BufferedImage(nPlotWidth+40, 360, BufferedImage.TYPE_INT_RGB);
				Graphics2D graph = img.createGraphics();
				graph.setColor(Color.WHITE);
				graph.fillRect(0, 0, nPlotWidth+40, 360);
				m_App.m_plotFactory.BoxPlot(mapRatios, "Ratios " + ex.m_nStart+"-"+ex.m_nEnd, graph, 30, 30, 300, "Ratio");
				
				Imagemap imgMap = new Imagemap();
				imgMap.setWidth(nPlotWidth+40 + "px");
				imgMap.setHeight("360px");
				imgMap.setContent(img);
				imgMap.setParent(grpBox);
			}
		}
		
		if(res.HasPSIScore())
		{
			SimpleSpliceScore score = res.GetPSIScore();
			
			grpBox = new Groupbox();
			grpBox.setTitle("Junctions");
			
			grpBox.setMold("3d");
			grpBox.setParent(layout);
			grpBox.setWidth((nBoxWidth + 10) + "px");
			grpBox.setOpen(false);
			
			grid = new Grid();
			grid.setParent(grpBox);
			grid.setWidth(nBoxWidth + "px");
			
			//##################################
			//           add header
			//##################################
			cols = new Columns();
			cols.setParent(grid);
			//	cols.setMenupopup("auto");
			
			col = new Column("");
			col.setWidth("300px");
			col.setParent(cols);
			
			col = new Column("value");
			col.setWidth("200px");
			col.setParent(cols);
			
			//##################################
			//             add rows
			//##################################
			rows = new Rows();
			rows.setParent(grid);
			
			row = new Row();
			row.setParent(rows);
			label = new Label("inclusion junction:");
			label.setParent(row);
			label = new Label(String.format(Locale.ENGLISH, "%d - %d", score.m_JunctionInclusion.m_nStart, score.m_JunctionInclusion.m_nEnd));
			label.setParent(row);
			
			row = new Row();
			row.setParent(rows);
			label = new Label("exclusion junction:");
			label.setParent(row);
			label = new Label(String.format(Locale.ENGLISH, "%d - %d", score.m_JunctionExclusion.m_nStart, score.m_JunctionExclusion.m_nEnd));
			label.setParent(row);
			
			row = new Row();
			row.setParent(rows);
			label = new Label("includes novel junction");
			label.setParent(row);
			if(score.m_bIsNovel)
				label = new Label("yes");
			else
				label = new Label("no");
			label.setParent(row);
			
			row = new Row();
			row.setParent(rows);
			label = new Label("inclusion change");
			label.setParent(row);
			label = new Label(String.format(Locale.ENGLISH, "%.2f", score.m_fInclusionChange));
			label.setParent(row);
			
			row = new Row();
			row.setParent(rows);
			label = new Label("p-value");
			label.setParent(row);
			label = new Label(format.format(score.m_fPValue));
			label.setParent(row);
			
			row = new Row();
			row.setParent(rows);
			label = new Label("significant");
			label.setParent(row);
			if(score.m_bSignificant)
				label = new Label("yes");
			else
				label = new Label("no");
			label.setParent(row);
			
			// draw box plots
			DataSupplier dataSupplier = m_App.GetDataSupplier();
			ProjectModel projectModel = m_App.GetProjectModel();
			TreeMap<String, TreeSet<String>> mapSamplesToConditions = projectModel.GetSamplesPerCondition(m_App.GetSelectedConditionType());
			TreeSet<String> vcSelectedSamples = m_App.GetSelectedSamples();
			
			//####################################################
			//             calculate junction coverage
			//####################################################
			
			// get size factors
			TreeMap<String, Double> mapSizeFactorsToSamples = projectModel.GetSizeFactors();
			
			// get junction coverages
			TreeMap<String, Integer> vcCountsA = dataSupplier.GetCountsForJunction(score.m_JunctionInclusion);
			TreeMap<String, Integer> vcCountsB = dataSupplier.GetCountsForJunction(score.m_JunctionExclusion);
			
			TreeMap<String, double[]> mapCountsToConditionsA = new TreeMap<String, double[]>();
			TreeMap<String, double[]> mapCountsToConditionsB = new TreeMap<String, double[]>();
			TreeMap<String, double[]> mapRatiosToConditions  = new TreeMap<String, double[]>();
			Vector<String> vcDataLabels = new Vector<String>();

			for(String strCondition : mapSamplesToConditions.keySet())
			{
				Vector<Double> vcValuesA = new Vector<Double>();
				Vector<Double> vcValuesB = new Vector<Double>();
				
				for(String strSample : mapSamplesToConditions.get(strCondition))
				{
					if(vcSelectedSamples.contains(strSample))
					{
						vcDataLabels.add(strSample);
						
						double fCount = 0.0;
						if(vcCountsA.containsKey(strSample))
							fCount = vcCountsA.get(strSample) / mapSizeFactorsToSamples.get(strSample);
						vcValuesA.add(fCount);
						
						fCount = 0.0;
						if(vcCountsB.containsKey(strSample))
							fCount = vcCountsB.get(strSample) / mapSizeFactorsToSamples.get(strSample);
						vcValuesB.add(fCount);
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
			
			// settings for the box plot
			int nBarWidth 	= 80;
			int nSpacer   	= 20;
			int nTotalBarDistance = nBarWidth+nSpacer;
			int nPlotWidth = nTotalBarDistance*2+40;
			int nXOffset = 60;
			
			String strJunctionNameA = dataSupplier.GetReferenceName() + ":" + (score.m_JunctionInclusion.m_nStart+1) + "-" + (score.m_JunctionInclusion.m_nEnd-1);
			String strJunctionNameB = dataSupplier.GetReferenceName() + ":" + (score.m_JunctionExclusion.m_nStart+1) + "-" + (score.m_JunctionExclusion.m_nEnd-1);
			
			BufferedImage img = new BufferedImage(nPlotWidth+nXOffset, 360, BufferedImage.TYPE_INT_RGB);
			Graphics2D graph = img.createGraphics();
			graph.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			graph.setColor(Color.WHITE);
			graph.fillRect(0, 0, nPlotWidth+nXOffset, 360);
			m_App.m_plotFactory.BoxPlot(mapCountsToConditionsA, "junction: " + strJunctionNameA, graph, nXOffset, 30, 300, "Coverage");
			
			Imagemap imgMap = new Imagemap();
			imgMap.setWidth(nPlotWidth+nXOffset + "px");
			imgMap.setHeight("360px");
			imgMap.setContent(img);
			imgMap.setParent(grpBox);
			
			img = new BufferedImage(nPlotWidth+nXOffset, 360, BufferedImage.TYPE_INT_RGB);
			graph = img.createGraphics();
			graph.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			graph.setColor(Color.WHITE);
			graph.fillRect(0, 0, nPlotWidth+nXOffset, 360);
			m_App.m_plotFactory.BoxPlot(mapCountsToConditionsB, "junction: " + strJunctionNameA, graph, nXOffset, 30, 300, "Coverage");
			
			imgMap = new Imagemap();
			imgMap.setWidth(nPlotWidth+nXOffset + "px");
			imgMap.setHeight("360px");
			imgMap.setContent(img);
			imgMap.setParent(grpBox);
			
			img = new BufferedImage(nPlotWidth+nXOffset, 360, BufferedImage.TYPE_INT_RGB);
			graph = img.createGraphics();
			graph.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			graph.setColor(Color.WHITE);
			graph.fillRect(0, 0, nPlotWidth+nXOffset, 360);
			m_App.m_plotFactory.BoxPlot(mapRatiosToConditions, "Ratios", graph, nXOffset, 30, 300, "Ratio");
			
			imgMap = new Imagemap();
			imgMap.setWidth(nPlotWidth+nXOffset + "px");
			imgMap.setHeight("360px");
			imgMap.setContent(img);
			imgMap.setParent(grpBox);
			
			nPlotWidth = 300+100;
			
			img = new BufferedImage(nPlotWidth, 400, BufferedImage.TYPE_INT_RGB);
			graph = img.createGraphics();
			graph.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			graph.setColor(Color.WHITE);
			graph.fillRect(0, 0, nPlotWidth, 400);
			Vector<Area> vcAreas = m_App.m_plotFactory.Plot2D(mapCountsToConditionsA, mapCountsToConditionsB, vcDataLabels, strJunctionNameA, strJunctionNameB, "2D plot", graph, nXOffset, 30, 300, 300, false);
			
			imgMap = new Imagemap();
			imgMap.setWidth(nPlotWidth + "px");
			imgMap.setHeight("400px");
			imgMap.setContent(img);
			imgMap.setParent(grpBox);
			
			for(Area area : vcAreas)
			{
				area.setParent(imgMap);
			}
		}
	}
	
	/**
	 *    Updates the contents of the temporary result list.
	 *    This function adds the table to the temporary result list tab
	 *    and fills it with data.
	 *    It is invoked each time a new gene is selected, the "reanalyze"
	 *    button is clicked or if the filter rules change.
	 */
	public void UpdateTemporaryResultList(AnalysisResultHandler resultHandler) throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		// update results
		boolean bInit = false;
		if(resultHandler != null)
		{
			m_resultHandler = resultHandler;
			bInit = true;
		}
		
		NumberFormat nf = NumberFormat.getNumberInstance(Locale.ENGLISH);
		DecimalFormat format = (DecimalFormat)nf;
		format.applyPattern("0.####E0");
				
		m_resultHandler.ClearTemporaryResults();
		m_resultHandler.UpdateTemporaryResults();
		
		m_gridTmpHitList.setMold("default");
		
		// clear old rows
		m_gridTmpHitList.getChildren().clear();
		m_gridTmpHitList.setMold("paging");
		m_gridTmpHitList.setPageSize(10);
		
		Column colPSI = AddColumnHeadersToGridForManananggalResults(m_gridTmpHitList, GRID_TYPE_TEMPORARY, 0, bInit);
		
		Rows rows = new Rows();
		rows.setParent(m_gridTmpHitList);
		
		//##################################
		//             add rows
		//##################################
		TreeSet<AnalysisResult> vcResults = m_resultHandler.GetAllTemporaryResults();
		
		for(AnalysisResult res : vcResults)
		{
			if(!m_resultFilterRuleTemporaryList.bIsValidResult(res))
				continue;
			
			Row row = new Row();
			row.setId("TMPID_" + res.GetID());
			row.setValue(res.GetRating());
			row.setParent(rows);
			
			AddCellToRow(row, "", "", "40px", "float: center", res, 7);
			AddCellToRow(row, "", "TMPview_", "40px", "float: center", res, 2);
			AddCellToRow(row, "", "TMPdetailView_", "40px", "float: center", res, 3);
			AddCellToRow(row, "", "TMPdetails_", "40px", "float: center", res, 4);
			AddCellToRow(row, "", "TMPrating_", "110px", "float: center", res, 5);
			
			AddCellToRow(row, res.GetGeneID(),             "", "140px", "float: center", res, 0);
			AddCellToRow(row, res.GetGeneSymbol(),         "", "140px", "float: center", res, 0);
			AddCellToRow(row, res.GetConditionA(),         "", "160px", "float: center", res, 0);
			AddCellToRow(row, res.GetConditionB(),         "", "160px", "float: center", res, 0);
			AddCellToRow(row, res.GetASTypeAsString(),       "", "170px", "float: center", res, 0);
			AddCellToRow(row, res.GetResultTypeAsString(), "", "120px", "float: center", res, 0);
			
			if(res.HasAltExonA())
			{
				String strText = String.format(Locale.ENGLISH, "%.2f%%", res.GetAbsoluteChangeA()*100);
				if(Math.abs(res.GetAbsoluteChangeA()) < 0.05)
					AddCellToRow(row, strText, "", "210px", "float: right; color: red;", res, 0);
				else
					AddCellToRow(row, strText, "", "210px", "float: right", res, 0);
			}
			else
			{
				AddCellToRow(row, "NA", "", "210px", "float: right", res, 0);
			}
			
			if(res.HasAltExonB())
			{
				String strText = String.format(Locale.ENGLISH, "%.2f%%", res.GetAbsoluteChangeB()*100);
				if(Math.abs(res.GetAbsoluteChangeB()) < 0.05)
					AddCellToRow(row, strText, "", "210px", "float: right; color: red;", res, 0);
				else
					AddCellToRow(row, strText, "", "210px", "float: right", res, 0);
			}
			else
			{
				AddCellToRow(row, "NA", "", "210px", "float: right", res, 0);
			}
			
			if(res.HasPSIScore())
			{
				String strText = String.format(Locale.ENGLISH, "%.2f%%", res.GetInclusionChange()*100);
				if(Math.abs(res.GetInclusionChange()) < 0.05)
					AddCellToRow(row, strText, "", "150px", "float: right; color: red;", res, 0);
				else
					AddCellToRow(row, strText, "", "150px", "float: right", res, 0);
			}
			else
			{
				AddCellToRow(row, "NA", "", "150px", "float: right", res, 0);
			}
			
			if(res.HasAltExonA())
			{
				AddCellToRow(row, format.format(res.GetPValueA()), "", "210px", "float: right", res, 0);	
			}
			else
			{
				AddCellToRow(row, "NA", "", "210px", "float: right", res, 0);
			}
	
			if(res.HasAltExonB())
			{
				AddCellToRow(row, format.format(res.GetPValueB()), "", "210px", "float: right", res, 0);	
			}
			else
			{
				AddCellToRow(row, "NA", "", "210px", "float: right", res, 0);
			}
	
			if(res.HasPSIScore())
			{
				AddCellToRow(row, format.format(res.GetPValuePSI()), "", "150px", "float: right", res, 0);
				AddCellToRow(row, String.format(Locale.ENGLISH, "%s", res.HasSignificantPSIScore()), "", "150px", "float: right", res, 0);
				AddCellToRow(row, String.format(Locale.ENGLISH, "%s", res.UsesNovelJunction()), "", "170px", "float: right", res, 0);
			}
			else
			{
				AddCellToRow(row, "NA", "", "150px", "float: right", res, 0);
				AddCellToRow(row, "NA", "", "150px", "float: center", res, 0);			
				AddCellToRow(row, "NA", "", "170px", "float: center", res, 0);
			}
			
			int nEventStart = Integer.MAX_VALUE;
			int nEventEnd	= Integer.MIN_VALUE;
			if(res.HasAltExonA())
			{
				nEventStart = Math.min(res.GetStartA(), nEventStart);
				nEventEnd 	= Math.max(res.GetEndA(),   nEventEnd);
			}
			if(res.HasAltExonB())
			{
				nEventStart = Math.min(res.GetStartB(), nEventStart);
				nEventEnd	= Math.max(res.GetEndB(),   nEventEnd);
			}
			if(res.HasPSIScore())
			{
				nEventStart = Math.min(res.GetInclusionJunctionStart(), nEventStart);
				nEventEnd 	= Math.max(res.GetInclusionJunctionEnd(),   nEventEnd);
				
				nEventStart = Math.min(res.GetExclusionJunctionStart(), nEventStart);
				nEventEnd 	= Math.max(res.GetExclusionJunctionEnd(),   nEventEnd);
			}

			AddCellToRow(row, String.format(Locale.ENGLISH, "%s:%,d - %,d", res.GetReferenceName(), nEventStart, nEventEnd), "", "200px", "float: left", res, 0);
			AddCellToRow(row, res.GetComment(), "TMPCommentField_", "400px", "float: left", res, 6);
		}

		colPSI.sort(true);
		
		m_Tabs.invalidate();
	}

	/**
	 *    Updates the contents of the imported result list.
	 */
	public void UpdateImportedResultList(Vector<AnalysisResult> vcResults, int nDataType, boolean bInit) throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		if(bInit)
		{
			m_vcImportedResults.clear();
			m_vcImportedResults = vcResults;
		}
		
		// clear old rows
		m_gridImportedResults.setMold("default");
		m_gridImportedResults.getChildren().clear();
		m_gridImportedResults.setMold("paging");
		m_gridImportedResults.setPageSize(10);
		
		// reset filter rules on data initialization
		if(bInit) ModifyFilterValuesForImportedData();
		
		// prepare new rows
		Rows rows = new Rows();
		rows.setParent(m_gridImportedResults);
		
		Column colPSI = null;
		switch(nDataType)
		{
			case SplicingWebApp.IMPORTED_DATA_TYPE_MANA:
			{
				// add columns
				colPSI = AddColumnHeadersToGridForImportedManananggalResults(m_gridImportedResults, GRID_TYPE_IMPORT, nDataType);
				// add rows
				AddImportedManananggalDataToGrid(rows);
				break;
			}
				
			case SplicingWebApp.IMPORTED_DATA_TYPE_DEXSEQ:
			{
				// add columns
				colPSI = AddColumnHeadersToGridForImportedDEXSeqResults(m_gridImportedResults, GRID_TYPE_IMPORT, nDataType);
				// add rows
				AddImportedDEXSeqDataToGrid(rows);
				break;
			}
			
			case SplicingWebApp.IMPORTED_DATA_TYPE_RMATS:
			{
				// add columns
				colPSI = AddColumnHeadersToGridForImportedRMATSResults(m_gridImportedResults, GRID_TYPE_IMPORT, nDataType);
				// add rows
				AddImportedRMATSDataToGrid(rows);
				break;
			}
			
			case SplicingWebApp.IMPORTED_DATA_TYPE_JSPLICE:
			{
				// add columns
				colPSI = AddColumnHeadersToGridForImportedJSpliceResults(m_gridImportedResults, GRID_TYPE_IMPORT, nDataType);
				// add rows
				AddImportedJSpliceDataToGrid(rows);
				break;
			}
			
			case SplicingWebApp.IMPORTED_DATA_TYPE_CUFFDIFF:
			{
				// add columns
				colPSI = AddColumnHeadersToGridForImportedCuffdiffResults(m_gridImportedResults, GRID_TYPE_IMPORT, nDataType);
				// add rows
				AddImportedCuffdiffDataToGrid(rows);
				break;
			}
			
			default:
			{
				ErrorMessage.ShowError("Unknown file type of imported file");
				return;
			}
		}

		colPSI.sort(true);
		m_Tabs.invalidate();
	}
	
	/** Helper function used to sort rows by ascending user ratings. */
	private Comparator<Row> getAscRatingComparator()
	{
		return new Comparator<Row>()
		{
			@Override
			public int compare(Row arg0, Row arg1)
			{
				int nRating1 = (int)arg0.getValue();
				int nRating2 = (int)arg1.getValue();
				
				if(nRating1 < nRating2)
					return -1;
				else if(nRating1 > nRating2)
					return 1;
				else
					return 0;
			}
		};
	}
	
	/** Helper function used to sort rows by descending user ratings */
	private Comparator<Row> getDescRatingComparator()
	{
		return new Comparator<Row>()
		{
			@Override
			public int compare(Row arg0, Row arg1)
			{
				int nRating1 = (int)arg0.getValue();
				int nRating2 = (int)arg1.getValue();
				
				if(nRating1 > nRating2)
					return -1;
				else if(nRating1 < nRating2)
					return 1;
				else
					return 0;
			}
		};
	}

	/**
	 *    Helper function used to sort number containing columns in descending order
	 *    by the value specified in the column nIdx.
	 */
	private Comparator<Row> getDescNumberComparator(int nColIdx)
	{	
		final int m_nIdx = nColIdx;
		
		return new Comparator<Row>()
		{
			@Override
			public int compare(Row arg0, Row arg1)
			{
				Label lab1 = (Label)arg0.getChildren().get(m_nIdx).getChildren().get(0);
				Label lab2 = (Label)arg1.getChildren().get(m_nIdx).getChildren().get(0);
				
				String strValue1 = lab1.getValue();
				String strValue2 = lab2.getValue();
				
				if(!strValue1.equals("NA") && !strValue2.equals("NA"))
				{
					double fValue1 = Double.parseDouble(strValue1.replace("%", ""));//new BigDecimal(strValue1).longValueExact();
					double fValue2 = Double.parseDouble(strValue2.replace("%", ""));//new BigDecimal(strValue2).longValueExact();
					
					if(fValue1 < fValue2)
						return +1;
					else if(fValue1 > fValue2)
						return -1;
					else
						return 0;
				}
				else if(!strValue1.equals("NA") && strValue2.equals("NA"))
					return +1;
				else if(strValue1.equals("NA") && !strValue2.equals("NA"))
					return -1;
				else
					return 0;
			}
		};
	}
	
	/**
	 *    Helper function used to sort number containing columns in ascending order
	 *    by the value specified in the column nIdx.
	 */
	private Comparator<Row> getAscNumberComparator(int nColIdx)
	{
		final int m_nIdx = nColIdx;
			
		return new Comparator<Row>()
		{
			@Override
			public int compare(Row arg0, Row arg1)
			{				
				Label lab1 = (Label)arg0.getChildren().get(m_nIdx).getChildren().get(0);
				Label lab2 = (Label)arg1.getChildren().get(m_nIdx).getChildren().get(0);
				
				String strValue1 = lab1.getValue();
				String strValue2 = lab2.getValue();
				
				if(!strValue1.equals("NA") && !strValue2.equals("NA"))
				{
					double fValue1 = Double.parseDouble(strValue1.replace("%", ""));
					double fValue2 = Double.parseDouble(strValue2.replace("%", ""));
					
					if(fValue1 > fValue2)
						return +1;
					else if(fValue1 < fValue2)
						return -1;
					else
						return 0;
				}
				else if(!strValue1.equals("NA") && strValue2.equals("NA"))
					return -1;
				else if(strValue1.equals("NA") && !strValue2.equals("NA"))
					return +1;
				else
					return 0;
			}
		};
	}

	/**
	 *    Helper function used to sort number containing columns in descending order
	 *    by the value specified in the column nIdx.
	 */
	private Comparator<Row> getDescTextComparator(int nColIdx)
	{	
		final int m_nIdx = nColIdx;
		
		return new Comparator<Row>()
		{
			@Override
			public int compare(Row arg0, Row arg1)
			{
				Label lab1 = (Label)arg0.getChildren().get(m_nIdx).getChildren().get(0);
				Label lab2 = (Label)arg1.getChildren().get(m_nIdx).getChildren().get(0);
				
				String strValue1 = lab1.getValue();
				String strValue2 = lab2.getValue();
				
				return strValue1.compareTo(strValue2) * -1;
			}
		};
	}
	
	/**
	 *    Helper function used to sort number containing columns in ascending order
	 *    by the value specified in the column nIdx.
	 */
	private Comparator<Row> getAscTextComparator(int nColIdx)
	{
		final int m_nIdx = nColIdx;
			
		return new Comparator<Row>()
		{
			@Override
			public int compare(Row arg0, Row arg1)
			{				
				Label lab1 = (Label)arg0.getChildren().get(m_nIdx).getChildren().get(0);
				Label lab2 = (Label)arg1.getChildren().get(m_nIdx).getChildren().get(0);
				
				String strValue1 = lab1.getValue();
				String strValue2 = lab2.getValue();
				
				return strValue1.compareTo(strValue2);
			}
		};
	}

	/**
	 * Adds column headers for manananggal results to the grid.
	 * Also returns the PSI column for auto-sorting.
	 * GridType 0 = permanent results
	 * GridType 1 = temporary results
	 * GridType 2 = imported results
	 */
	private Column AddColumnHeadersToGridForManananggalResults(Grid grid, int nGridType, int nDataType, boolean bInit) throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		ResultFilterRule resultFilterRule = null;
		switch(nGridType)
		{
			case GRID_TYPE_PERMANENT: resultFilterRule = m_resultFilterRulePermanentList; break;
			case GRID_TYPE_TEMPORARY: resultFilterRule = m_resultFilterRuleTemporaryList; break;
			case GRID_TYPE_IMPORT:    resultFilterRule = m_resultFilterRuleImportedList;  break;
		}
		
		// prepare filter values on first open
		if(bInit)
		{
			TreeSet<String> vcConditionsA = new TreeSet<String>();
			TreeSet<String> vcConditionsB = new TreeSet<String>();

			TreeSet<AnalysisResult> vcResults = null;
			if(nGridType == GRID_TYPE_TEMPORARY)
				vcResults = m_resultHandler.GetAllTemporaryResults();
			else
				vcResults = m_resultHandler.GetAllResults();

			for(AnalysisResult res : vcResults)
			{
				vcConditionsA.add(res.GetConditionA());
				vcConditionsB.add(res.GetConditionB());
			}
			
			TextFilter filter;
			filter = new TextFilter(vcConditionsA);
			resultFilterRule.SetTextFilter(ResultFilterRule.FILTER_TYPE_CONDITION_A, filter);
			filter = new TextFilter(vcConditionsB);
			resultFilterRule.SetTextFilter(ResultFilterRule.FILTER_TYPE_CONDITION_B, filter);

			TreeSet<String> vcStrings = new TreeSet<String>();
			vcStrings.add("true");
			vcStrings.add("false");
			vcStrings.add("NA");
			filter = new TextFilter(vcStrings);
			resultFilterRule.SetTextFilter(ResultFilterRule.FILTER_TYPE_SIGNIFICANT_PSI, filter);
			filter = new TextFilter(vcStrings);
			resultFilterRule.SetTextFilter(ResultFilterRule.FILTER_TYPE_NOVEL_JUNCTION, filter);
		}
		
		//##################################
		//           add header
		//##################################
		Columns cols = new Columns();
		cols.setSizable(false);
		cols.setParent(grid);
			
		int nColIdx = 0;		
		AddColumnToGrid(cols, "", "30px", "center", false, nColIdx, -1, -1, nGridType, nDataType); nColIdx++; // space for result removal			
		AddColumnToGrid(cols, "", "40px", "center", false, nColIdx, -1, -1, nGridType, nDataType); nColIdx++; // space for 'view result'
		AddColumnToGrid(cols, "", "40px", "center", false, nColIdx, -1, -1, nGridType, nDataType); nColIdx++; // space for 'detailed view'
		AddColumnToGrid(cols, "", "40px", "center", false, nColIdx, -1, -1, nGridType, nDataType); nColIdx++; // space for 'detailed description'
		AddColumnToGrid(cols, "Rating",      "110px", "center", true, nColIdx, ResultFilterRule.OTHER_FILTER, ResultFilterRule.FILTER_RATING, nGridType, nDataType); nColIdx++; // space for 'detailed description'
		AddColumnToGrid(cols, "Gene ID",     "140px", "center", true, nColIdx, -1, -1, nGridType, nDataType); nColIdx++;
		AddColumnToGrid(cols, "Gene Symbol", "140px", "center", true, nColIdx, -1, -1, nGridType, nDataType); nColIdx++;
		AddColumnToGrid(cols, "Condition A", "160px", "center", true, nColIdx,  ResultFilterRule.TEXT_FILTER, ResultFilterRule.FILTER_TYPE_CONDITION_A, nGridType, nDataType); nColIdx++;
		AddColumnToGrid(cols, "Condition B", "160px", "center", true, nColIdx,  ResultFilterRule.TEXT_FILTER, ResultFilterRule.FILTER_TYPE_CONDITION_B, nGridType, nDataType); nColIdx++;
		AddColumnToGrid(cols, "AS Type",     "170px", "center", true, nColIdx,  ResultFilterRule.OTHER_FILTER, ResultFilterRule.FILTER_TYPE_AS_TYPE, nGridType, nDataType); nColIdx++;
		AddColumnToGrid(cols, "Result Type", "120px", "center", true, nColIdx, ResultFilterRule.OTHER_FILTER, ResultFilterRule.FILTER_TYPE_RESULT_TYPE, nGridType, nDataType); nColIdx++;
		AddColumnToGrid(cols, "Ratio Change (Exon A)",  "210px", "right", true, nColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_RATIO_CHANGE_A, nGridType, nDataType); nColIdx++;
		AddColumnToGrid(cols, "Ratio Change (Exon B)",  "210px", "right", true, nColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_RATIO_CHANGE_B, nGridType, nDataType); nColIdx++;		
		AddColumnToGrid(cols, "PSI Change",             "150px", "right", true, nColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_PSI_CHANGE, nGridType, nDataType); nColIdx++;
		AddColumnToGrid(cols, "p-value Ratio (Exon A)", "210px", "right", true, nColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_P_VALUE_RATIO_A, nGridType, nDataType); nColIdx++;
		AddColumnToGrid(cols, "p-value Ratio (Exon B)", "210px", "right", true, nColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_P_VALUE_RATIO_B, nGridType, nDataType); nColIdx++;
		Column colPSI = AddColumnToGrid(cols, "p-value PSI", "150px", "right", true, nColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_P_VALUE_PSI, nGridType, nDataType); nColIdx++;
		AddColumnToGrid(cols, "Significant",    "150px", "center", true, nColIdx, ResultFilterRule.TEXT_FILTER, ResultFilterRule.FILTER_TYPE_SIGNIFICANT_PSI, nGridType, nDataType); nColIdx++;
		AddColumnToGrid(cols, "Novel Junction", "170px", "center", true, nColIdx, ResultFilterRule.TEXT_FILTER, ResultFilterRule.FILTER_TYPE_NOVEL_JUNCTION, nGridType, nDataType); nColIdx++;
		AddColumnToGrid(cols, "Event Location", "200px","center",  false, nColIdx, -1, -1, nGridType, nDataType); nColIdx++;
		AddColumnToGrid(cols, "Comment",        "410px", "left",   false, nColIdx, -1, -1, nGridType, nDataType);
		
		return colPSI;
	}
	
	/**
	 * Adds column headers for imported manananggal results to the grid
	 */
	private Column AddColumnHeadersToGridForImportedManananggalResults(Grid grid, int nGridType, int nDataType) throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		//##################################
		//           add header
		//##################################
		Columns cols = new Columns();
		cols.setParent(grid);
		
		int ColIdx = 0;
		AddColumnToGrid(cols, "", "40px", "center", false, ColIdx, -1, -1, nGridType, nDataType); ColIdx++; // space for 'view result'
		AddColumnToGrid(cols, "Gene ID",     "140px", "center", true, ColIdx, -1, -1, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Gene Symbol", "140px", "center", true, ColIdx, -1, -1, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Condition A", "160px", "center", true, ColIdx, ResultFilterRule.TEXT_FILTER, ResultFilterRule.FILTER_TYPE_CONDITION_A, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Condition B", "160px", "center", true, ColIdx, ResultFilterRule.TEXT_FILTER, ResultFilterRule.FILTER_TYPE_CONDITION_B, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "AS Type",     "170px", "center", true, ColIdx, ResultFilterRule.OTHER_FILTER, ResultFilterRule.FILTER_TYPE_AS_TYPE, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Result Type", "120px", "center", true, ColIdx, ResultFilterRule.OTHER_FILTER, ResultFilterRule.FILTER_TYPE_RESULT_TYPE, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Ratio Change",  "210px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_RATIO_CHANGE_A, nGridType, nDataType); ColIdx++;	
		AddColumnToGrid(cols, "PSI Change",    "150px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_PSI_CHANGE, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "p-value Ratio", "210px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_P_VALUE_RATIO_A, nGridType, nDataType); ColIdx++;
		Column colPSI = AddColumnToGrid(cols, "p-value PSI", "150px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_P_VALUE_PSI, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Significant",    "150px", "center", true, ColIdx, ResultFilterRule.TEXT_FILTER, ResultFilterRule.FILTER_TYPE_SIGNIFICANT_PSI, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Novel Junction", "170px", "center", true, ColIdx, ResultFilterRule.TEXT_FILTER, ResultFilterRule.FILTER_TYPE_NOVEL_JUNCTION, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Event Location", "200px", "center", false, ColIdx, -1, -1, nGridType, nDataType);

		return colPSI;
	}
	
	/**
	 * Adds column headers for imported DEXSeq results to the grid
	 */
	private Column AddColumnHeadersToGridForImportedDEXSeqResults(Grid grid, int nGridType, int nDataType)
	{
		//##################################
		//           add header
		//##################################
		Columns cols = new Columns();
		cols.setParent(grid);
		
		int ColIdx = 0;
		AddColumnToGrid(cols, "", "40px", "center", false, ColIdx, -1, -1, nGridType, nDataType); ColIdx++; // space for 'view result'
		AddColumnToGrid(cols, "Gene ID",     "140px", "center", true, ColIdx, -1, -1, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Gene Symbol", "140px", "center", true, ColIdx, -1, -1, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Condition A", "160px", "center", true, ColIdx, ResultFilterRule.TEXT_FILTER, ResultFilterRule.FILTER_TYPE_CONDITION_A, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Condition B", "160px", "center", true, ColIdx, ResultFilterRule.TEXT_FILTER, ResultFilterRule.FILTER_TYPE_CONDITION_B, nGridType, nDataType); ColIdx++;
		Column colPSI = AddColumnToGrid(cols, "p-value", "140px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_P_VALUE_RATIO_A, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "padj.", "140px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_P_VALUE_RATIO_A, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "log2FC",  "140px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_RATIO_CHANGE_A, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Event Location", "200px", "center", false, ColIdx, -1, -1, nGridType, nDataType);

		return colPSI;
	}

	/**
	 * Adds column headers for imported rMATS results to the grid
	 */
	private Column AddColumnHeadersToGridForImportedRMATSResults(Grid grid, int nGridType, int nDataType)
	{
		//##################################
		//           add header
		//##################################
		Columns cols = new Columns();
		cols.setParent(grid);
		
		int ColIdx = 0;
		AddColumnToGrid(cols, "", "40px", "center", false, ColIdx, -1, -1, nGridType, nDataType); ColIdx++; // space for 'view result'
		AddColumnToGrid(cols, "Gene ID",     "140px", "center", true, ColIdx, -1, -1, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Gene Symbol", "140px", "center", true, ColIdx, -1, -1, nGridType, nDataType); ColIdx++;
		Column colPSI = AddColumnToGrid(cols, "p-value", "140px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_P_VALUE_RATIO_A, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "FDR", "140px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_P_VALUE_RATIO_A, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "inclusion change",  "140px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_RATIO_CHANGE_A, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Event Location", "200px", "center", false, ColIdx, -1, -1, nGridType, nDataType);

		return colPSI;
	}
	
	/**
	 * Adds column headers for imported JSplice results to the grid
	 */
	private Column AddColumnHeadersToGridForImportedJSpliceResults(Grid grid, int nGridType, int nDataType)
	{
		//##################################
		//           add header
		//##################################
		Columns cols = new Columns();
		cols.setParent(grid);
		
		int ColIdx = 0;
		AddColumnToGrid(cols, "", "40px", "center", false, ColIdx, -1, -1, nGridType, nDataType); ColIdx++; // space for 'view result'
		AddColumnToGrid(cols, "Gene ID",     "140px", "center", true, ColIdx, -1, -1, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Gene Symbol", "140px", "center", true, ColIdx, -1, -1, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "AS Type",     "170px", "center", true, ColIdx, ResultFilterRule.OTHER_FILTER, ResultFilterRule.FILTER_TYPE_AS_TYPE, nGridType, nDataType); ColIdx++;
		Column colPSI = AddColumnToGrid(cols, "adj. p-value", "140px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_P_VALUE_RATIO_A, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "log2FC",  "140px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_RATIO_CHANGE_A, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "Event Location", "200px", "center", false, ColIdx, -1, -1, nGridType, nDataType);

		return colPSI;
	}
	
	/**
	 * Adds column headers for imported Cuffdiff results to the grid
	 */
	private Column AddColumnHeadersToGridForImportedCuffdiffResults(Grid grid, int nGridType, int nDataType)
	{
		//##################################
		//           add header
		//##################################
		Columns cols = new Columns();
		cols.setParent(grid);
		
		int ColIdx = 0;
		AddColumnToGrid(cols, "", "40px", "center", false, ColIdx, -1, -1, nGridType, nDataType); ColIdx++; // space for 'view result'
		AddColumnToGrid(cols, "Gene Symbol", "140px", "center", true, ColIdx, -1, -1, nGridType, nDataType); ColIdx++;
		Column colPSI = AddColumnToGrid(cols, "p-value", "140px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_P_VALUE_RATIO_A, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "q-value", "140px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_P_VALUE_RATIO_A, nGridType, nDataType); ColIdx++;
		AddColumnToGrid(cols, "sqrt(JS)",  "140px", "right", true, ColIdx, ResultFilterRule.NUMERIC_FILTER, ResultFilterRule.FILTER_TYPE_RATIO_CHANGE_A, nGridType, nDataType);

		return colPSI;
	}

	/**
	 *   Adds a single column to the grid. May include a filter rule.
	 *   - Cols is the column parent object
	 *   - strLabel is the text shown in the column
	 *   - strColWidth specifies the width of the column, e.g. "100px"
	 *   - strAlign specifies the alignment of the text, e.g. "center"
	 *   - bEnableSorting adds sorting functionality to the column
	 *   - nComparatorColIdx specifies that the column (specified by this index) uses a default numeric comparator
	 *   - nFiterType defines the type of filter rule (-1 = none, 0 = numeric filter, 1 = text filter, 2 = other filter)
	 *   - nFilterSubtype specifies the filter sub type, e.g. FILTER_TYPE_CONDITION_A
	 */
	private Column AddColumnToGrid(Columns cols, String strLabel, String strColWidth, String strAlign, boolean bEnableSorting, int nComparatorColIdx, int nFilterType, int nFilterSubType, int nGridType, int nDataType)
	{
		Column col = new Column();		
		col.setAlign(strAlign);
		col.setParent(cols);
		col.setWidth(strColWidth);
		
		if(nFilterType == -1)
		{
			if(bEnableSorting)
			{
				col.setSortAscending(getAscTextComparator(nComparatorColIdx));
				col.setSortDescending(getDescTextComparator(nComparatorColIdx));
			}
			col.setLabel(strLabel);
		}
		else
		{
			Label lab = new Label(strLabel);
			lab.setParent(col);
			lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
			
			switch(nFilterType)
			{
				case ResultFilterRule.NUMERIC_FILTER:
				{
					col.setSortAscending(getAscNumberComparator(nComparatorColIdx));
					col.setSortDescending(getDescNumberComparator(nComparatorColIdx));
					AddNumberFilter(col, nFilterSubType, nGridType, nDataType, nComparatorColIdx);
					break;
				}
				
				case ResultFilterRule.TEXT_FILTER:
				{
					col.setSortAscending(getAscTextComparator(nComparatorColIdx));
					col.setSortDescending(getDescTextComparator(nComparatorColIdx));
					AddStringFilter(col, nFilterSubType, nGridType, nDataType, nComparatorColIdx);
					break;
				}
				
				case ResultFilterRule.OTHER_FILTER:
				{
					if(nFilterSubType == ResultFilterRule.FILTER_TYPE_AS_TYPE)
					{
						col.setSortAscending(getAscTextComparator(nComparatorColIdx));
						col.setSortDescending(getDescTextComparator(nComparatorColIdx));
						AddSplicingTypeFilter(col, nGridType, nDataType, nComparatorColIdx);						
					}
					else if(nFilterSubType == ResultFilterRule.FILTER_TYPE_RESULT_TYPE)
					{
						col.setSortAscending(getAscTextComparator(nComparatorColIdx));
						col.setSortDescending(getDescTextComparator(nComparatorColIdx));
						AddResultTypeFilter(col, nGridType, nDataType, nComparatorColIdx);
					}
					else if(nFilterSubType == ResultFilterRule.FILTER_RATING)
					{
						col.setSortAscending(getAscRatingComparator());
						col.setSortDescending(getDescRatingComparator());
					}
					break;
				}
			}
		}
		
		if(bEnableSorting)
		{
			try
			{
				col.setSort("auto");
			}
			catch(Exception e)
			{
				System.out.println("ERROR: Failed to set filter type for column: " + strLabel);
				e.printStackTrace();
				return null;
			}
			
			col.addEventListener(Events.ON_SORT, new EventListener<Event>()
			{
				@Override
				public void onEvent(Event event) throws Exception
				{
					m_Tabs.invalidate(); // if we don't force a redraw of the tabs, the grid will be truncated after sorting
				}
				
			});
		}
		
		return col;
	}

	/**
	 * Adds a cell to a row in the grid.
	 * nType: 0 = default, 1=remove button, 2=view button, 3=detail button, 4 = description button
	 * nGridType: GRID_TYPE_PERMANENT, GRID_TYPE_TEMPORARY, GRID_TYPE_IMPORT
	 */
	private void AddCellToRow(Row parent, String strText, String strPrefixID, String strWidth, String strStyle, AnalysisResult res, int nType)
	{
		Cell cell = new Cell();
		cell.setParent(parent);
		cell.setAlign("center");
		cell.setWidth(strWidth);

		switch(nType)
		{
			case 0:
			{
				Label label = new Label(strText);
				label.setStyle(strStyle);
				label.setParent(cell);
				break;
			}
		
			case 1:
			{
				String strImageString = "/img/red_cross.png";
				Image img = new Image(strImageString);
				img.setHeight("18px");
				img.setParent(cell);
				img.setId(strPrefixID + res.GetID());
				img.setStyle("display: block; margin-left: auto; margin-right: auto; margin-top: 4px; margin-bottom: auto; cursor:pointer;");
				
				img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
				{
					public void onEvent(Event event) throws Exception
					{
						Image img = (Image)event.getTarget();
		
						int nID = Integer.parseInt(img.getId().split("_")[1]);
						
						m_resultHandler.RemoveResult(nID);
						UpdatePermanentResultList(null);
					}
				});
				break;
			}
		
			case 2:
			{
				if(res.GetType() == SplicingWebApp.AS_TYPE_RETAINED_INTRON)
				{
					Label label = new Label("");			
					label.setParent(cell);
				}
				else
				{
					String strImageString = "/img/view.png";
					
					if(res.UsesNovelJunction())
						strImageString = "/img/view_incomplete.png";
					
					Image img = new Image(strImageString);
					img.setHeight("12px");
					img.setParent(cell);
					img.setId(strPrefixID + res.GetID());
					img.setStyle("display: block; margin-left: auto; margin-right: auto; margin-top: 4px; margin-bottom: auto; cursor:pointer;");
					
					img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
					{
						public void onEvent(Event event) throws Exception
						{
							Image img = (Image)event.getTarget();
		
							String pSplit[] = img.getId().split("_");
							String strGridType = pSplit[0];
							int nGridType = GRID_TYPE_PERMANENT;
							if(strGridType.startsWith("TMP"))
							{
								nGridType = GRID_TYPE_TEMPORARY;
							}
							else if(strGridType.startsWith("IMP"))
							{
								nGridType = GRID_TYPE_IMPORT;
							}
							int nID = Integer.parseInt(pSplit[1]);
							
							AnalysisResult res = null;
							if(nGridType == GRID_TYPE_PERMANENT)
							{
								res = m_resultHandler.GetResult(nID);
							}
							else if(nGridType == GRID_TYPE_TEMPORARY)
							{
								res = m_resultHandler.GetTemporaryResult(nID);
							}
							else if(nGridType == GRID_TYPE_IMPORT)
							{
								res = GetImportedResult(nID);
							}
							m_App.PrepareHitForVisualization(res);
							m_App.m_plotFactory.RequestCoverageRedraw();
							m_App.m_plotFactory.DrawPlots();
						}
					});
				}
				break;
			}
		
			case 3:
			{
				if(res.HasPSIScore() || res.GetType() == SplicingWebApp.AS_TYPE_RETAINED_INTRON)
				{
					String strImageString = "/img/magnifier_blue.png";
					
					if(res.UsesNovelJunction())
						strImageString = "/img/magnifier_red.png";
					
					Image img = new Image(strImageString);
					img.setHeight("20px");
					img.setParent(cell);
					img.setId(strPrefixID + res.GetID());
					img.setStyle("display: block; margin-left: auto; margin-right: auto; margin-top: 4px; margin-bottom: auto; cursor:pointer;");
					
					img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
					{
						public void onEvent(Event event) throws Exception
						{
							Image img = (Image)event.getTarget();

							String pSplit[] = img.getId().split("_");
							String strGridType = pSplit[0];
							int nGridType = GRID_TYPE_PERMANENT;
							if(strGridType.startsWith("TMP"))
							{
								nGridType = GRID_TYPE_TEMPORARY;
							}
							else if(strGridType.startsWith("IMP"))
							{
								nGridType = GRID_TYPE_IMPORT;
							}
							int nID = Integer.parseInt(pSplit[1]);
							
							AnalysisResult res = null;
							if(nGridType == GRID_TYPE_PERMANENT)
							{
								res = m_resultHandler.GetResult(nID);
							}
							else if(nGridType == GRID_TYPE_TEMPORARY)
							{
								res = m_resultHandler.GetTemporaryResult(nID);
							}
							else if(nGridType == GRID_TYPE_IMPORT)
							{
								res = GetImportedResult(nID);
							}
							
							m_App.PrepareHitForVisualization(res);
							
							if(res.GetType() == SplicingWebApp.AS_TYPE_RETAINED_INTRON)
							{
								m_App.m_plotFactory.DrawRetainedIntronCloseup(res);
							}
							else
							{
								m_App.m_plotFactory.DrawSpliceJunctionCloseup(res);
							}
						}
					});
				}
				else
				{
					Label label = new Label("");			
					label.setParent(cell);
				}
				break;
			}
		
			case 4:
			{
				String strImageString = "/img/detailed_desc.png";
				Image img = new Image(strImageString);
				img.setWidth("16px");
				img.setHeight("20px");
				img.setParent(cell);
				img.setId(strPrefixID + res.GetID());
				img.setStyle("display: block; margin-left: auto; margin-right: auto; margin-top: 4px; margin-bottom: auto; cursor:hand; cursor:pointer;");
				
				img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
				{
					public void onEvent(Event event) throws Exception
					{
						Image img = (Image)event.getTarget();
						
						String pSplit[] = img.getId().split("_");
						String strGridType = pSplit[0];
						int nGridType = GRID_TYPE_PERMANENT;
						if(strGridType.startsWith("TMP"))
						{
							nGridType = GRID_TYPE_TEMPORARY;
						}
						else if(strGridType.startsWith("IMP"))
						{
							nGridType = GRID_TYPE_IMPORT;
						}
						int nID = Integer.parseInt(pSplit[1]);
						
						AnalysisResult res = null;
						if(nGridType == GRID_TYPE_PERMANENT)
						{
							res = m_resultHandler.GetResult(nID);
						}
						else if(nGridType == GRID_TYPE_TEMPORARY)
						{
							res = m_resultHandler.GetTemporaryResult(nID);
						}
						else if(nGridType == GRID_TYPE_IMPORT)
						{
							res = GetImportedResult(nID);
						}

						m_App.PrepareHitForVisualization(res);
						ShowDetailPopup(res);
					}
				});
				break;
			}
			
			case 5:
			{
				Hlayout layoutH = new Hlayout();
				layoutH.setParent(cell);
				for(int i=0; i<5; i++)
				{
					Image img = null;
					
					if(i<res.GetRating())
					{
						img = new Image("/img/good_rating.png");
					}
					else
					{
						img = new Image("/img/neutral_rating.png");
					}
					
					img.setHeight("15px");
					img.setParent(layoutH);
					img.setId(strPrefixID + (i+1) + "_" + res.GetID());
					img.setStyle("cursor:hand; cursor:pointer;");
					
					img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
					{
						public void onEvent(Event event) throws Exception
						{
							Image img = (Image)event.getTarget();
							
							String pSplit[] = img.getId().split("_");
							String strGridType = pSplit[0];
							int nGridType = GRID_TYPE_PERMANENT;
							if(strGridType.startsWith("TMP"))
							{
								nGridType = GRID_TYPE_TEMPORARY;
							}
							else if(strGridType.startsWith("IMP"))
							{
								nGridType = GRID_TYPE_IMPORT;
							}
							int nNewRating = Integer.parseInt(pSplit[1]);
							int nID = Integer.parseInt(pSplit[2]);
							
							AnalysisResult res = null;
							if(nGridType == GRID_TYPE_PERMANENT)
							{
								res = m_resultHandler.GetResult(nID);
							}
							else if(nGridType == GRID_TYPE_TEMPORARY)
							{
								res = m_resultHandler.GetTemporaryResult(nID);
							}
							else if(nGridType == GRID_TYPE_IMPORT)
							{
								res = GetImportedResult(nID);
							}
							
							int nOldRating = res.GetRating();
							
							// remove rating if the same star was clicked again.
							if(nNewRating == nOldRating)
								nNewRating = 0;
							
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
		
							res.SetRating(nNewRating);
							Row row = (Row)img.getParent().getParent().getParent();
							row.setValue(nNewRating);
						}
					});
				}
				break;
			}
			
			case 6:
			{
				Textbox boxComment = new Textbox();
				boxComment.setParent(cell);
				boxComment.setInplace(true);
				boxComment.setWidth("99%");
				boxComment.setStyle("margin-top: -5px; margin-bottom: -5px;");
				boxComment.setHflex("1");
				boxComment.setId(strPrefixID + res.GetID());
				boxComment.setText(res.GetComment());
				boxComment.setTooltiptext(res.GetComment());
				
				boxComment.addEventListener(Events.ON_CHANGE, new EventListener<Event>()
				{
					public void onEvent(Event event) throws Exception
					{
						Textbox box = (Textbox)event.getTarget();
						
						String strID = box.getId();
						
						String pSplit[] = strID.split("_");
						int nID = Integer.parseInt(pSplit[1]);
						
						if(strID.startsWith("TMP"))
						{
							AnalysisResult res = m_resultHandler.GetTemporaryResult(nID);
							res.SetComment(box.getText());
						}
						else
						{
							AnalysisResult res = m_resultHandler.GetResult(nID);
							res.SetComment(box.getText());
						}
					}
				});
				break;
			}
			
			case 7:
			{
				Checkbox checkbox = new Checkbox();
				checkbox.setValue(res.GetID());
				checkbox.setHeight("20px");
				checkbox.setWidth("20px");
				checkbox.setParent(cell);
				break;
			}
		}
	}

	private void ModifyFilterValuesForImportedData()
	{
		// reset filters for the imported data
		m_resultFilterRuleImportedList = new ResultFilterRule();
		
		// collect condition names for the condition filter
		TreeSet<String> vcConditionsA = new TreeSet<String>();
		TreeSet<String> vcConditionsB = new TreeSet<String>();
		
		// process imported results
		for(AnalysisResult res : m_vcImportedResults)
		{
			vcConditionsA.add(res.GetConditionA());
			vcConditionsB.add(res.GetConditionB());
		}
		
		// prepare the filters
		// modify filters
		TextFilter filter;
		filter = new TextFilter(vcConditionsA);
		m_resultFilterRuleImportedList.SetTextFilter(ResultFilterRule.FILTER_TYPE_CONDITION_A, filter);
		filter = new TextFilter(vcConditionsB);
		m_resultFilterRuleImportedList.SetTextFilter(ResultFilterRule.FILTER_TYPE_CONDITION_B, filter);
			
		// these filters don't need to be changed afterwards, they are fine as is
		TreeSet<String> vcStrings = new TreeSet<String>();
		vcStrings.add("true");
		vcStrings.add("false");
		vcStrings.add("NA");
		filter = new TextFilter(vcStrings);
		m_resultFilterRuleImportedList.SetTextFilter(ResultFilterRule.FILTER_TYPE_SIGNIFICANT_PSI, filter);
		filter = new TextFilter(vcStrings);
		m_resultFilterRuleImportedList.SetTextFilter(ResultFilterRule.FILTER_TYPE_NOVEL_JUNCTION, filter);
	}
	
	private void AddImportedManananggalDataToGrid(Rows parent)
	{
		// prepare number format
		NumberFormat nf = NumberFormat.getNumberInstance(Locale.ENGLISH);
		DecimalFormat format = (DecimalFormat)nf;
		format.applyPattern("0.####E0");
		
		// process imported results
		for(AnalysisResult res : m_vcImportedResults)
		{			
			if(!m_resultFilterRuleImportedList.bIsValidResult(res))
				continue;
			
			Row row = new Row();
			row.setId("IMPID_" + res.GetID());
			row.setValue(res.GetRating());
			row.setParent(parent);

			AddCellToRow(row, "", "IMPview_", "40px", "float: center", res, 2);
			
			AddCellToRow(row, res.GetGeneID(),             "", "140px", "float: center", res, 0);
			AddCellToRow(row, res.GetGeneSymbol(),         "", "140px", "float: center", res, 0);
			AddCellToRow(row, res.GetConditionA(),         "", "160px", "float: center", res, 0);
			AddCellToRow(row, res.GetConditionB(),         "", "160px", "float: center", res, 0);
			AddCellToRow(row, res.GetASTypeAsString(),     "", "170px", "float: center", res, 0);
			AddCellToRow(row, res.GetResultTypeAsString(), "", "120px", "float: center", res, 0);
			
			if(res.HasAltExonA())
			{
				String strText = String.format(Locale.ENGLISH, "%.2f%%", res.GetAbsoluteChangeA()*100);
				if(Math.abs(res.GetAbsoluteChangeA()) < 0.05)
					AddCellToRow(row, strText, "", "210px", "float: right; color: red;", res, 0);
				else
					AddCellToRow(row, strText, "", "210px", "float: right", res, 0);
			}
			else
			{
				AddCellToRow(row, "NA", "", "210px", "float: right", res, 0);
			}
			
			if(res.HasPSIScore())
			{
				String strText = String.format(Locale.ENGLISH, "%.2f%%", res.GetInclusionChange()*100);
				if(Math.abs(res.GetInclusionChange()) < 0.05)
					AddCellToRow(row, strText, "", "150px", "float: right; color: red;", res, 0);
				else
					AddCellToRow(row, strText, "", "150px", "float: right", res, 0);
			}
			else
			{
				AddCellToRow(row, "NA", "", "150px", "float: right", res, 0);
			}
			
			if(res.HasAltExonA())
			{
				AddCellToRow(row, format.format(res.GetPValueA()), "", "210px", "float: right", res, 0);	
			}
			else
			{
				AddCellToRow(row, "NA", "", "210px", "float: right", res, 0);
			}
	
			if(res.HasPSIScore())
			{
				AddCellToRow(row, format.format(res.GetPValuePSI()), "", "150px", "float: right", res, 0);
				AddCellToRow(row, String.format(Locale.ENGLISH, "%s", res.HasSignificantPSIScore()), "", "150px", "float: right", res, 0);
				AddCellToRow(row, String.format(Locale.ENGLISH, "%s", res.UsesNovelJunction()), "", "170px", "float: right", res, 0);
			}
			else
			{
				AddCellToRow(row, "NA", "", "150px", "float: right", res, 0);
				AddCellToRow(row, "NA", "", "150px", "float: center", res, 0);			
				AddCellToRow(row, "NA", "", "170px", "float: center", res, 0);
			}
			
			int nEventStart = Integer.MAX_VALUE;
			int nEventEnd	= Integer.MIN_VALUE;
			if(res.HasAltExonA())
			{
				nEventStart = Math.min(res.GetStartA(), nEventStart);
				nEventEnd 	= Math.max(res.GetEndA(),   nEventEnd);
			}

			if(res.HasPSIScore())
			{
				nEventStart = Math.min(res.GetInclusionJunctionStart(), nEventStart);
				nEventEnd 	= Math.max(res.GetInclusionJunctionEnd(),   nEventEnd);
				
				nEventStart = Math.min(res.GetExclusionJunctionStart(), nEventStart);
				nEventEnd 	= Math.max(res.GetExclusionJunctionEnd(),   nEventEnd);
			}

			AddCellToRow(row, String.format(Locale.ENGLISH, "%s:%,d - %,d", res.GetReferenceName(), nEventStart, nEventEnd), "", "200px", "float: left", res, 0);
		}
	}
	
	private void AddImportedDEXSeqDataToGrid(Rows parent)
	{
		// prepare number format
		NumberFormat nf = NumberFormat.getNumberInstance(Locale.ENGLISH);
		DecimalFormat format = (DecimalFormat)nf;
		format.applyPattern("0.####E0");
		
		// process imported results
		for(AnalysisResult res : m_vcImportedResults)
		{			
			if(!m_resultFilterRuleImportedList.bIsValidResult(res))
				continue;
			
			Row row = new Row();
			row.setId("IMPID_" + res.GetID());
			row.setValue(res.GetRating());
			row.setParent(parent);

			AddCellToRow(row, "", "IMPview_", "40px", "float: center", res, 2);
			
			AddCellToRow(row, res.GetGeneID(),             "", "140px", "float: center", res, 0);
			AddCellToRow(row, res.GetGeneSymbol(),         "", "140px", "float: center", res, 0);
			AddCellToRow(row, res.GetConditionA(),         "", "160px", "float: center", res, 0);
			AddCellToRow(row, res.GetConditionB(),         "", "160px", "float: center", res, 0);
			
			String strText = format.format(res.GetImportedDataPValue());
			if(Math.abs(res.GetImportedDataPValue()) > 0.05)
				AddCellToRow(row, strText, "", "140px", "float: right; color: red;", res, 0);
			else
				AddCellToRow(row, strText, "", "140px", "float: right", res, 0);

			strText = format.format(res.GetImportedDataAdjustedPValue());
			if(Math.abs(res.GetImportedDataAdjustedPValue()) > 0.1)
				AddCellToRow(row, strText, "", "140px", "float: right; color: red;", res, 0);
			else
				AddCellToRow(row, strText, "", "140px", "float: right", res, 0);
			
			strText = String.format(Locale.ENGLISH, "%.4f", res.GetImportedDataLog2FC());
			AddCellToRow(row, strText, "", "140px", "float: right", res, 0);
			
			int nEventStart = res.GetStartA();
			int nEventEnd	= res.GetEndA();

			AddCellToRow(row, String.format(Locale.ENGLISH, "%s:%,d - %,d", res.GetReferenceName(), nEventStart, nEventEnd), "", "200px", "float: left", res, 0);
		}
	}
	
	private void AddImportedRMATSDataToGrid(Rows parent)
	{
		// prepare number format
		NumberFormat nf = NumberFormat.getNumberInstance(Locale.ENGLISH);
		DecimalFormat format = (DecimalFormat)nf;
		format.applyPattern("0.####E0");
		
		// process imported results
		for(AnalysisResult res : m_vcImportedResults)
		{			
			if(!m_resultFilterRuleImportedList.bIsValidResult(res))
				continue;
			
			Row row = new Row();
			row.setId("IMPID_" + res.GetID());
			row.setValue(res.GetRating());
			row.setParent(parent);

			AddCellToRow(row, "", "IMPview_", "40px", "float: center", res, 2);
			
			AddCellToRow(row, res.GetGeneID(),             "", "140px", "float: center", res, 0);
			AddCellToRow(row, res.GetGeneSymbol(),         "", "140px", "float: center", res, 0);
			
			String strText = format.format(res.GetImportedDataPValue());
			if(Math.abs(res.GetImportedDataPValue()) > 0.05)
				AddCellToRow(row, strText, "", "140px", "float: right; color: red;", res, 0);
			else
				AddCellToRow(row, strText, "", "140px", "float: right", res, 0);

			strText = format.format(res.GetImportedDataFDR());
			if(Math.abs(res.GetImportedDataFDR()) > 0.1)
				AddCellToRow(row, strText, "", "140px", "float: right; color: red;", res, 0);
			else
				AddCellToRow(row, strText, "", "140px", "float: right", res, 0);
			
			strText = String.format(Locale.ENGLISH, "%.4f", res.GetAbsoluteChangeA());
			if(Math.abs(res.GetAbsoluteChangeA()) < 0.05)
				AddCellToRow(row, strText, "", "140px", "float: right; color: red;", res, 0);
			else
				AddCellToRow(row, strText, "", "140px", "float: right", res, 0);
			
			int nEventStart = res.GetStartA();
			int nEventEnd	= res.GetEndA();

			AddCellToRow(row, String.format(Locale.ENGLISH, "%s:%,d - %,d", res.GetReferenceName(), nEventStart, nEventEnd), "", "200px", "float: left", res, 0);
		}
	}

	private void AddImportedJSpliceDataToGrid(Rows parent)
	{
		// prepare number format
		NumberFormat nf = NumberFormat.getNumberInstance(Locale.ENGLISH);
		DecimalFormat format = (DecimalFormat)nf;
		format.applyPattern("0.####E0");
		
		// process imported results
		for(AnalysisResult res : m_vcImportedResults)
		{			
			if(!m_resultFilterRuleImportedList.bIsValidResult(res))
				continue;
			
			Row row = new Row();
			row.setId("IMPID_" + res.GetID());
			row.setValue(res.GetRating());
			row.setParent(parent);

			AddCellToRow(row, "", "IMPview_", "40px", "float: center", res, 2);
			
			AddCellToRow(row, res.GetGeneID(),             "", "140px", "float: center", res, 0);
			AddCellToRow(row, res.GetGeneSymbol(),         "", "140px", "float: center", res, 0);
			AddCellToRow(row, res.GetASTypeAsString(),     "", "170px", "float: center", res, 0);
			
			String strText = format.format(res.GetImportedDataPValue());
			if(Math.abs(res.GetImportedDataPValue()) > 0.05)
				AddCellToRow(row, strText, "", "140px", "float: right; color: red;", res, 0);
			else
				AddCellToRow(row, strText, "", "140px", "float: right", res, 0);
			
			strText = String.format(Locale.ENGLISH, "%.4f", res.GetImportedDataLog2FC());
			if(Math.abs(res.GetImportedDataLog2FC()) < 2)
				AddCellToRow(row, strText, "", "140px", "float: right; color: red;", res, 0);
			else
				AddCellToRow(row, strText, "", "140px", "float: right", res, 0);
			
			int nEventStart = res.GetStartA();
			int nEventEnd	= res.GetEndA();

			AddCellToRow(row, String.format(Locale.ENGLISH, "%s:%,d - %,d", res.GetReferenceName(), nEventStart, nEventEnd), "", "200px", "float: left", res, 0);
		}
	}
	
	private void AddImportedCuffdiffDataToGrid(Rows parent)
	{
		// prepare number format
		NumberFormat nf = NumberFormat.getNumberInstance(Locale.ENGLISH);
		DecimalFormat format = (DecimalFormat)nf;
		format.applyPattern("0.####E0");
		
		// process imported results
		for(AnalysisResult res : m_vcImportedResults)
		{			
			if(!m_resultFilterRuleImportedList.bIsValidResult(res))
				continue;
			
			Row row = new Row();
			row.setId("IMPID_" + res.GetID());
			row.setValue(res.GetRating());
			row.setParent(parent);

			AddCellToRow(row, "", "IMPview_", "40px", "float: center", res, 2);

			AddCellToRow(row, res.GetGeneSymbol(),         "", "140px", "float: center", res, 0);
			
			String strText = format.format(res.GetImportedDataPValue());
			if(Math.abs(res.GetImportedDataPValue()) > 0.05)
				AddCellToRow(row, strText, "", "140px", "float: right; color: red;", res, 0);
			else
				AddCellToRow(row, strText, "", "140px", "float: right", res, 0);
			
			strText = format.format(res.GetImportedDataAdjustedPValue());
			if(Math.abs(res.GetImportedDataAdjustedPValue()) > 0.1)
				AddCellToRow(row, strText, "", "140px", "float: right; color: red;", res, 0);
			else
				AddCellToRow(row, strText, "", "140px", "float: right", res, 0);
			
			strText = String.format(Locale.ENGLISH, "%.4f", res.GetImportedDataLog2FC());
			AddCellToRow(row, strText, "", "140px", "float: right", res, 0);
		}
	}
	
	/**
	 * Retrieves an imported splicing result
	 */
	private AnalysisResult GetImportedResult(int nID)
	{
		for(AnalysisResult res : m_vcImportedResults)
			if(res.GetID() == nID)
				return res;
		
		return null;
	}

}
