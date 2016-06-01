package Manananggal;

import java.awt.Color;
import java.awt.Graphics2D;
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

//######################################################
//    This class handles the permanent and temporary
//    result lists displayed in the GUI
//######################################################
public class ResultListHandler
{
	private Tabs 				m_Tabs;
	private Tabpanels 			m_Panels;
	
	private Grid				m_gridHitList;
	private Grid				m_gridTmpHitList;
	private ResultFilterRule	m_resultFilterRule;
	
	private AnalysisResultHandler m_resultHandler;
	
	SplicingWebApp				m_App;
	
	public ResultListHandler(Tabs tabs, Tabpanels panels, SplicingWebApp app)
	{
		m_Tabs			= tabs;
		m_Panels		= panels;
		m_App			= app;
		
		m_resultFilterRule = new ResultFilterRule();

		AddPermanentResultList();
		AddTemporaryResultList();
	}

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
	}
	
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
//		m_gridHitList.setWidth("830px");
		m_gridHitList.setHeight("380px");
		m_gridHitList.setVflex("min");
		m_gridHitList.setHflex("min");
		m_gridHitList.setMold("paging");
		m_gridHitList.setPageSize(10);
		
//		String strHeight = Integer.parseInt(m_gridHitList.getHeight().replace("px", "")) + 60 + "px";
//		layoutV.setHeight(strHeight);
		layoutV.setHeight("440px");
		
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
//		m_gridTmpHitList.setWidth("830px");
		m_gridTmpHitList.setHeight("380px");
		m_gridTmpHitList.setVflex("min");
		m_gridTmpHitList.setHflex("min");
		m_gridTmpHitList.setMold("paging");
		m_gridTmpHitList.setPageSize(10);
		
//		String strHeight = Integer.parseInt(m_gridTmpHitList.getHeight().replace("px", "")) + 60 + "px";
//		layoutV.setHeight(strHeight);
		layoutV.setHeight("440px");
		
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
					Checkbox checkbox = (Checkbox)row.getFirstChild();
					
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
				
//				m_vcASResults.AddTemporaryResults(vcResults);
//				Messagebox.show("Results saved.");
			}
		});
	}
	
	private void AddSplicingTypeFilter(Column parent)
	{		
		Bandbox bandBox = new Bandbox();
		bandBox.setStyle("position: absolute; right: 0px; margin-right: 10px");
		bandBox.setWidth("30px");
		bandBox.setParent(parent);
		
		Bandpopup bandPopup = new Bandpopup();
		bandPopup.setHflex("min");
		bandPopup.setParent(bandBox);
		
		Hlayout layout = new Hlayout();
		layout.setParent(bandPopup);
		
		Button btnSelectAll = new Button("Select All");
		btnSelectAll.setParent(layout);
		btnSelectAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				for(int i=0; i<13; i++)
					m_resultFilterRule.ShowASType(i);
				
				UpdatePermanentResultList(null);
				UpdateTemporaryResultList(null);
			}
		});
		
		Button btnHideAll = new Button("Unselect All");
		btnHideAll.setParent(layout);
		btnHideAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				for(int i=0; i<13; i++)
					m_resultFilterRule.HideASType(i);
				
				UpdatePermanentResultList(null);
				UpdateTemporaryResultList(null);
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
				
				for(Component c : popup.getChildren())
				{
					if(c instanceof Listbox)
					{
						Listbox listbox = (Listbox)c;
						
						for(Listitem item : listbox.getItems())
						{					
							if(item.isSelected())
								m_resultFilterRule.ShowASType((int)item.getValue());
							else
								m_resultFilterRule.HideASType((int)item.getValue());
						}
						
						UpdatePermanentResultList(null);
						UpdateTemporaryResultList(null);
						return;
					}
				}
			}
		});
		
		Listitem item = new Listitem("Exon skipping");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_EXON_SKIPPING);
		if(m_resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_EXON_SKIPPING))
			item.setSelected(true);
		
		item = new Listitem("retained introns");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_RETAINED_INTRON);
		if(m_resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_RETAINED_INTRON))
			item.setSelected(true);
		
		item = new Listitem("alt. 5' exon end");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_5_PRIME_EXON_END);
		if(m_resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_5_PRIME_EXON_END))
			item.setSelected(true);
				
		item = new Listitem("alt. 3' exon end");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_3_PRIME_EXON_END);
		if(m_resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_3_PRIME_EXON_END))
			item.setSelected(true);
		
		item = new Listitem("alt. start (unique junction; double)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN_DOUBLE);
		if(m_resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN_DOUBLE))
			item.setSelected(true);
		
		item = new Listitem("alt. end (unique junction; double)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN_DOUBLE);
		if(m_resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN_DOUBLE))
			item.setSelected(true);
		
		item = new Listitem("alt. start (shared junction; double)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN_DOUBLE);
		if(m_resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN_DOUBLE))
			item.setSelected(true);
		
		item = new Listitem("alt. end (shared junction; double)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN_DOUBLE);
		if(m_resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN_DOUBLE))
			item.setSelected(true);
		
		item = new Listitem("alt. start (unique junction)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN);
		if(m_resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_START_UNIQUE_JUN))
			item.setSelected(true);
		
		item = new Listitem("alt. end (unique junction)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN);
		if(m_resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_END_UNIQUE_JUN))
			item.setSelected(true);
		
		item = new Listitem("alt. start (shared junction)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN);
		if(m_resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_START_SHARED_JUN))
			item.setSelected(true);
		
		item = new Listitem("alt. end (shared junction)");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN);
		if(m_resultFilterRule.IsShowingASType(SplicingWebApp.AS_TYPE_ALT_END_SHARED_JUN))
			item.setSelected(true);
	}

	private void AddResultTypeFilter(Column parent)
	{		
		Bandbox bandBox = new Bandbox();
		bandBox.setStyle("position: absolute; right: 0px; margin-right: 10px");
		bandBox.setWidth("30px");
		bandBox.setParent(parent);
		
		Bandpopup bandPopup = new Bandpopup();
		bandPopup.setHflex("min");
		bandPopup.setParent(bandBox);
		
		Hlayout layout = new Hlayout();
		layout.setParent(bandPopup);
		
		Button btnSelectAll = new Button("Select All");
		btnSelectAll.setParent(layout);
		btnSelectAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				for(int i=0; i<4; i++)
					m_resultFilterRule.ShowResultType(i);
				
				UpdatePermanentResultList(null);
				UpdateTemporaryResultList(null);
			}
		});
		
		Button btnHideAll = new Button("Unselect All");
		btnHideAll.setParent(layout);
		btnHideAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				for(int i=0; i<4; i++)
					m_resultFilterRule.HideResultType(i);
				
				UpdatePermanentResultList(null);
				UpdateTemporaryResultList(null);
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
				
				for(Component c : popup.getChildren())
				{
					if(c instanceof Listbox)
					{
						Listbox listbox = (Listbox)c;
						
						for(Listitem item : listbox.getItems())
						{					
							if(item.isSelected())
								m_resultFilterRule.ShowResultType((int)item.getValue());
							else
								m_resultFilterRule.HideResultType((int)item.getValue());
						}
						
						UpdatePermanentResultList(null);
						UpdateTemporaryResultList(null);
						return;
					}
				}
			}
		});
		
		Listitem item = new Listitem("combined");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(AnalysisResult.RESULT_EVIDENCE_TYPE_COMBINED);
		if(m_resultFilterRule.IsShowingResultType(AnalysisResult.RESULT_EVIDENCE_TYPE_COMBINED))
			item.setSelected(true);
		
		item = new Listitem("ratio only");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(AnalysisResult.RESULT_EVIDENCE_TYPE_RATIO_ONLY);
		if(m_resultFilterRule.IsShowingResultType(AnalysisResult.RESULT_EVIDENCE_TYPE_RATIO_ONLY))
			item.setSelected(true);
		
		item = new Listitem("split read only");
		item.setParent(listBox);
		item.setCheckable(true);
		item.setValue(AnalysisResult.RESULT_EVIDENCE_TYPE_SPLIT_READ_ONLY);
		if(m_resultFilterRule.IsShowingResultType(AnalysisResult.RESULT_EVIDENCE_TYPE_SPLIT_READ_ONLY))
			item.setSelected(true);
	}

	private void AddStringFilter(Column parent, int nFilterID)
	{
		TreeSet<String> vcStrings = m_resultFilterRule.GetTextFilter(nFilterID).GetStrings();
		
		Bandbox bandBox = new Bandbox();
		bandBox.setStyle("position: absolute; right: 0px; margin-right: 10px");
		bandBox.setWidth("30px");
		bandBox.setParent(parent);
		bandBox.setValue(""+nFilterID);
		
		Bandpopup bandPopup = new Bandpopup();
		bandPopup.setHflex("min");
		bandPopup.setParent(bandBox);
		
		Hlayout layout = new Hlayout();
		layout.setParent(bandPopup);
		
		Button btnSelectAll = new Button("Select All");
		btnSelectAll.setParent(layout);
		btnSelectAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Bandbox bandbox = (Bandbox)event.getTarget().getParent().getParent().getParent();
				int nFilterID = Integer.parseInt(bandbox.getValue());
				
				m_resultFilterRule.GetTextFilter(nFilterID).SelectAll();
				
				UpdatePermanentResultList(null);
				UpdateTemporaryResultList(null);
			}
		});
		
		Button btnHideAll = new Button("Unselect All");
		btnHideAll.setParent(layout);
		btnHideAll.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Bandbox bandbox = (Bandbox)event.getTarget().getParent().getParent().getParent();
				int nFilterID = Integer.parseInt(bandbox.getValue());
				
				m_resultFilterRule.GetTextFilter(nFilterID).UnselectAll();
				
				UpdatePermanentResultList(null);
				UpdateTemporaryResultList(null);
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
				
				m_resultFilterRule.GetTextFilter(nFilterID).SelectAll();
				
				for(Component c : popup.getChildren())
				{
					if(c instanceof Listbox)
					{
						Listbox listbox = (Listbox)c;
						
						for(Listitem item : listbox.getItems())
						{
							if(item.isSelected())
								m_resultFilterRule.GetTextFilter(nFilterID).SelectString((String)item.getValue());
							else
								m_resultFilterRule.GetTextFilter(nFilterID).UnselectString((String)item.getValue());
						}
						
						UpdatePermanentResultList(null);
						UpdateTemporaryResultList(null);
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
			
			if(m_resultFilterRule.GetTextFilter(nFilterID).IsStringSelected(strString))
				item.setSelected(true);
		}
	}

	private void AddNumberFilter(Component parent, int nFilterID)
	{
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
		
		Textbox boxMinValue = new Textbox("" + m_resultFilterRule.GetNumberFilter(nFilterID).GetMinValue());
		boxMinValue.setConstraint("/[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?/: Only numbers are allowed, including decimal numbers and scientific writing, e.g. 5.0E-3");
		
		boxMinValue.setParent(layout);
		boxMinValue.addEventListener(Events.ON_CHANGE, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Textbox textbox = (Textbox)event.getTarget();
				Bandbox bandbox = (Bandbox)textbox.getParent().getParent().getParent().getParent();
				int nFilterID = Integer.parseInt(bandbox.getValue());
				
				float fValue = Float.parseFloat(textbox.getText());
				m_resultFilterRule.GetNumberFilter(nFilterID).SetMinValue(fValue);
				
				UpdatePermanentResultList(null);
				UpdateTemporaryResultList(null);
			}
		});
		
		layout = new Hlayout();
		layout.setParent(layoutV);
		
		labMinValue = new Label("max value: ");
		labMinValue.setParent(layout);
		
		Textbox boxMaxValue = new Textbox("" + m_resultFilterRule.GetNumberFilter(nFilterID).GetMaxValue());
		boxMinValue.setConstraint("/[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?/: Only numbers are allowed, including decimal numbers and scientific writing, e.g. 5.0E-3");
		boxMaxValue.setParent(layout);
		boxMaxValue.addEventListener(Events.ON_CHANGE, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Textbox textbox = (Textbox)event.getTarget();
				Bandbox bandbox = (Bandbox)textbox.getParent().getParent().getParent().getParent();
				int nFilterID = Integer.parseInt(bandbox.getValue());
				
				float fValue = Float.parseFloat(textbox.getText());
				m_resultFilterRule.GetNumberFilter(nFilterID).SetMaxValue(fValue);
				
				UpdatePermanentResultList(null);
				UpdateTemporaryResultList(null);
			}
		});
		
		Button btnDisable = new Button("Disable filter");
		btnDisable.setParent(layoutV);
		btnDisable.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			public void onEvent(Event event) throws Exception
			{
				Bandbox bandbox = (Bandbox)event.getTarget().getParent().getParent().getParent();
				int nFilterID = Integer.parseInt(bandbox.getValue());
				
				m_resultFilterRule.GetNumberFilter(nFilterID).Disable();
				
				UpdateTemporaryResultList(null);
			}
		});
	}
	
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
		
		if(bInit)
		{
			TreeSet<String> vcConditionsA = new TreeSet<String>();
			TreeSet<String> vcConditionsB = new TreeSet<String>();

			TreeSet<AnalysisResult> vcResults = m_resultHandler.GetAllResults();

			for(AnalysisResult res : vcResults)
			{
				vcConditionsA.add(res.GetConditionA());
				vcConditionsB.add(res.GetConditionB());
			}
			
			TextFilter filter;
			filter = new TextFilter(vcConditionsA);
			m_resultFilterRule.SetTextFilter(ResultFilterRule.FILTER_TYPE_CONDITION_A, filter);
			filter = new TextFilter(vcConditionsB);
			m_resultFilterRule.SetTextFilter(ResultFilterRule.FILTER_TYPE_CONDITION_B, filter);
			
			TreeSet<String> vcStrings = new TreeSet<String>();
			vcStrings.add("true");
			vcStrings.add("false");
			vcStrings.add("NA");
			filter = new TextFilter(vcStrings);
			m_resultFilterRule.SetTextFilter(ResultFilterRule.FILTER_TYPE_SIGNIFICANT_PSI, filter);
			filter = new TextFilter(vcStrings);
			m_resultFilterRule.SetTextFilter(ResultFilterRule.FILTER_TYPE_NOVEL_JUNCTION, filter);
			
		}
		
		m_gridHitList.setMold("default");
		
		// clear old rows
		m_gridHitList.getChildren().clear();
		m_gridHitList.setMold("paging");
		m_gridHitList.setPageSize(10);
	
		//##################################
		//           add header
		//##################################
		Columns cols = new Columns();
	//	cols.setMenupopup("auto");
		
		// image for result removal
		Column col = new Column("");
		col.setWidth("30px");
		col.setParent(cols);
	
		// view result
		col = new Column("");
		col.setWidth("40px");
		col.setParent(cols);
		
		// detailed view
		col = new Column("");
		col.setWidth("40px");
		col.setParent(cols);
		
		// detailed description
		col = new Column("");
		col.setWidth("40px");
		col.setParent(cols);
		
		col = new Column("Rating");
		col.setAlign("center");
		col.setSortAscending(getAscRatingComparator());
		col.setSortDescending(getDescRatingComparator());
		col.setWidth("110px");
		col.setParent(cols);
	
		col = new Column("Gene ID");
		col.setSort("auto");
		col.setAlign("center");
		col.setWidth("140px");
		col.setParent(cols);
		
		col = new Column("Gene Symbol");
		col.setSort("auto");
		col.setAlign("center");
		col.setWidth("100px");
		col.setParent(cols);
		
		col = new Column("Condition A");
		col.setAlign("center");
		col.setWidth("160px");
		col.setParent(cols);
		AddStringFilter(col, ResultFilterRule.FILTER_TYPE_CONDITION_A);
		
		col = new Column("Condition B");
		col.setAlign("center");
		col.setWidth("160px");
		col.setParent(cols);
		AddStringFilter(col, ResultFilterRule.FILTER_TYPE_CONDITION_B);
				
		col = new Column("");
		col.setSort("auto");
		col.setAlign("center");
		col.setWidth("170px");
		col.setParent(cols);
		Label lab = new Label("AS Type");
		lab.setParent(col);
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		AddSplicingTypeFilter(col);
		
		col = new Column("");
		col.setSort("auto");
		col.setWidth("120px");
		col.setAlign("center");
		col.setParent(cols);
		lab = new Label("Result Type");
		lab.setParent(col);
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		AddResultTypeFilter(col);
		
		col = new Column("");
		col.setSort("auto");
		col.setAlign("right");
		col.setWidth("210px");
		col.setParent(cols);
		lab = new Label("ratio change (exon A)");
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		lab.setParent(col);
		AddNumberFilter(col, ResultFilterRule.FILTER_TYPE_RATIO_CHANGE_A);
		
		col = new Column("");
		col.setSort("auto");
		col.setAlign("right");
		col.setWidth("210px");
		col.setParent(cols);
		lab = new Label("ratio change (exon B)");
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		lab.setParent(col);
		AddNumberFilter(col, ResultFilterRule.FILTER_TYPE_RATIO_CHANGE_B);
		
		col = new Column("");
		col.setSort("auto");
		col.setAlign("right");
		col.setWidth("150px");
		col.setParent(cols);
		lab = new Label("PSI change");
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		lab.setParent(col);
		AddNumberFilter(col, ResultFilterRule.FILTER_TYPE_PSI_CHANGE);
		
		col = new Column("");
		col.setSortAscending(getAscNumberComparator(14));
		col.setSortDescending(getDescNumberComparator(14));
		col.setSort("auto");
		col.setAlign("right");
		col.setWidth("210px");
		col.setParent(cols);
		lab = new Label("p-value ratio (exon A)");
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		lab.setParent(col);
		AddNumberFilter(col, ResultFilterRule.FILTER_TYPE_P_VALUE_RATIO_A);
		
		col = new Column("");
		col.setSortAscending(getAscNumberComparator(15));
		col.setSortDescending(getDescNumberComparator(15));
		col.setAlign("right");
		col.setWidth("210px");
		col.setParent(cols);
		lab = new Label("p-value ratio (exon B)");
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		lab.setParent(col);
		AddNumberFilter(col, ResultFilterRule.FILTER_TYPE_P_VALUE_RATIO_B);
		
		Column colPSI = new Column("");
		colPSI.setSortAscending(getAscNumberComparator(16));
		colPSI.setSortDescending(getDescNumberComparator(16));
		colPSI.setAlign("right");
		colPSI.setWidth("150px");
		colPSI.setParent(cols);
		lab = new Label("p-value PSI");
		lab.setParent(colPSI);
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		AddNumberFilter(colPSI, ResultFilterRule.FILTER_TYPE_P_VALUE_PSI);
		
		col = new Column("significant");
		col.setSort("auto");
		col.setAlign("center");
		col.setWidth("150px");
		col.setParent(cols);
		AddStringFilter(col, ResultFilterRule.FILTER_TYPE_SIGNIFICANT_PSI);
		
		col = new Column("novel junction");
		col.setSort("auto");
		col.setAlign("center");
		col.setWidth("170px");		
		col.setParent(cols);
		AddStringFilter(col, ResultFilterRule.FILTER_TYPE_NOVEL_JUNCTION);
		
		col = new Column("event Location");
		col.setSort("auto");
		col.setAlign("left");
		col.setWidth("200px");
		col.setParent(cols);
		
		col = new Column("Comment");
		col.setWidth("410px");
		col.setAlign("left");
		col.setParent(cols);
	
		cols.setParent(m_gridHitList);
		
		Rows rows = new Rows();
		rows.setParent(m_gridHitList);
		
		//##################################
		//             add rows
		//##################################
		TreeSet<AnalysisResult> vcResults = m_resultHandler.GetAllResults();
		
		for(AnalysisResult res : vcResults)
		{
			if(!m_resultFilterRule.bIsValidResult(res))
				continue;
			
			Row row = new Row();
			row.setId("ID_" + res.GetID());
			row.setValue(res.GetRating());
			row.setParent(rows);
			
			String strImageString = "/img/red_cross.png";
			Image img = new Image(strImageString);
			img.setHeight("18px");
			img.setParent(row);
			img.setId("remove_" + res.GetID());
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
			
			if(res.GetType() == SplicingWebApp.AS_TYPE_RETAINED_INTRON)
			{
				Label label = new Label("");			
				label.setParent(row);
			}
			else
			{
				strImageString = "/img/view.png";
				
				if(res.UsesNovelJunction())
					strImageString = "/img/view_incomplete.png";
				
				img = new Image(strImageString);
				img.setHeight("12px");
				img.setParent(row);
				img.setId("view_" + res.GetID());
				img.setStyle("display: block; margin-left: auto; margin-right: auto; margin-top: 4px; margin-bottom: auto; cursor:pointer;");
				
				img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
				{
					public void onEvent(Event event) throws Exception
					{
						Image img = (Image)event.getTarget();
	
						int nID = Integer.parseInt(img.getId().split("_")[1]);
						
						AnalysisResult res = m_resultHandler.GetResult(nID);
						m_App.PrepareHitForVisualization(res);
						m_App.m_plotFactory.RequestCoverageRedraw();
						m_App.m_plotFactory.DrawPlots();
					}
				});
			}
			
			if(res.HasPSIScore() || res.GetType() == SplicingWebApp.AS_TYPE_RETAINED_INTRON)
			{
				strImageString = "/img/magnifier_blue.png";
				
				if(res.UsesNovelJunction())
					strImageString = "/img/magnifier_red.png";
				
				img = new Image(strImageString);
				img.setHeight("20px");
				img.setParent(row);
				img.setId("detailView_" + res.GetID());
				img.setStyle("display: block; margin-left: auto; margin-right: auto; margin-top: 4px; margin-bottom: auto; cursor:pointer;");
				
				img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
				{
					public void onEvent(Event event) throws Exception
					{
						Image img = (Image)event.getTarget();
	
						int nID = Integer.parseInt(img.getId().split("_")[1]);
						
						AnalysisResult res = m_resultHandler.GetResult(nID);
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
				label.setParent(row);
			}
			
			strImageString = "/img/detailed_desc.png";
			img = new Image(strImageString);
			img.setWidth("16px");
			img.setHeight("20px");
			img.setParent(row);
			img.setId("details_" + res.GetID());
			img.setStyle("display: block; margin-left: auto; margin-right: auto; margin-top: 4px; margin-bottom: auto; cursor:hand; cursor:pointer;");
			
			img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
			{
				public void onEvent(Event event) throws Exception
				{
					Image img = (Image)event.getTarget();

					int nID = Integer.parseInt(img.getId().split("_")[1]);
					AnalysisResult res = m_resultHandler.GetResult(nID);
					m_App.PrepareHitForVisualization(res);
					ShowDetailPopup(res);
				}
			});
			
			Hlayout layoutH = new Hlayout();
			layoutH.setParent(row);
			for(int i=0; i<5; i++)
			{
				img = null;
				
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
				img.setId("rating_" + (i+1) + "_" + res.GetID());
				img.setStyle("cursor:hand; cursor:pointer;");
				
				img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
				{
					public void onEvent(Event event) throws Exception
					{
						Image img = (Image)event.getTarget();
						
						String strID = img.getId();
						
						String pSplit[] = strID.split("_");
						int nNewRating = Integer.parseInt(pSplit[1]);
						int nID = Integer.parseInt(pSplit[2]);
						
						AnalysisResult res = m_resultHandler.GetResult(nID);
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
						Row row = (Row)img.getParent().getParent();
						row.setValue(nNewRating);
					}
				});
			}
			
			Label label = new Label(res.GetGeneID());
			label.setStyle("float: center");
			label.setParent(row);
			
			label = new Label(res.GetGeneSymbol());
			label.setStyle("float: center");
			label.setParent(row);
			
			label = new Label(res.GetConditionA());
			label.setStyle("float: center");
			label.setParent(row);
			
			label = new Label(res.GetConditionB());
			label.setStyle("float: center");
			label.setParent(row);
			
			label = new Label(res.GetTypeAsString());
			label.setStyle("float: center");
			label.setParent(row);
			
			label = new Label(res.GetResultTypeAsString());
			label.setStyle("float: center");
			label.setParent(row);
			
			if(res.HasAltExonA())
			{
				label = new Label(String.format(Locale.ENGLISH, "%.2f%%", res.GetAbsoluteChangeA()*100));
				if(Math.abs(res.GetAbsoluteChangeA()) < 0.05)
					label.setStyle("color: red");
				label.setStyle("float: right");
				label.setParent(row);
			}
			else
			{
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setStyle("float: right");
				label.setParent(row);
			}
			
			if(res.HasAltExonB())
			{
				label = new Label(String.format(Locale.ENGLISH, "%.2f%%", res.GetAbsoluteChangeB()*100));
				if(Math.abs(res.GetAbsoluteChangeB()) < 0.05)
					label.setStyle("color: red");
				label.setStyle("float: right");
				label.setParent(row);
			}
			else
			{
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setStyle("float: right");
				label.setParent(row);
			}
			
			if(res.HasPSIScore())
			{
				label = new Label(String.format(Locale.ENGLISH, "%.2f%%", res.GetInclusionChange()*100));
				if(res.GetInclusionChange() < 0.05)
					label.setStyle("color: red");
				label.setStyle("float: right");
				label.setParent(row);
			}
			else
			{
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setStyle("float: right");
				label.setParent(row);
			}
			
			if(res.HasAltExonA())
			{
				label = new Label(format.format(res.GetPValueA()));
				label.setStyle("float: right");
				label.setParent(row);				
			}
			else
			{
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setStyle("float: right");
				label.setParent(row);
			}
	
			if(res.HasAltExonB())
			{
				label = new Label(format.format(res.GetPValueB()));
				label.setStyle("float: right");
				label.setParent(row);				
			}
			else
			{
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setStyle("float: right");
				label.setParent(row);				
			}
	
			if(res.HasPSIScore())
			{
				label = new Label(format.format(res.GetPValuePSI()));
				label.setStyle("float: right");
				label.setParent(row);
				
				label = new Label(String.format(Locale.ENGLISH, "%s", res.HasSignificantPSIScore()));
				label.setStyle("align: center");
				label.setParent(row);
				
				label = new Label(String.format(Locale.ENGLISH, "%s", res.UsesNovelJunction()));
				label.setStyle("align: center");
				label.setParent(row);
			}
			else
			{
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setParent(row);
				
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setParent(row);
				
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setParent(row);
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
			
			label = new Label(String.format(Locale.ENGLISH, "%s:%,d - %,d", res.GetReferenceName(), nEventStart, nEventEnd));
			label.setStyle("float: left");
			label.setParent(row);
			
			label = new Label(res.GetComment());
			label.setStyle("float: left");
			label.setParent(row);
		}
		
		colPSI.sort(true);
	}
	
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
				m_App.m_plotFactory.BoxPlot(mapRatios, "Ratios " + ex.m_nStart+"-"+ex.m_nEnd, graph, 30, 30, 300);
				
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
							fCount = vcCountsA.get(strSample) * mapSizeFactorsToSamples.get(strSample);
						vcValuesA.add(fCount);
						
						fCount = 0.0;
						if(vcCountsB.containsKey(strSample))
							fCount = vcCountsB.get(strSample) * mapSizeFactorsToSamples.get(strSample);
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
			
			String strJunctionNameA = dataSupplier.GetReferenceName() + ":" + (score.m_JunctionInclusion.m_nStart+1) + "-" + (score.m_JunctionInclusion.m_nEnd-1);
			String strJunctionNameB = dataSupplier.GetReferenceName() + ":" + (score.m_JunctionExclusion.m_nStart+1) + "-" + (score.m_JunctionExclusion.m_nEnd-1);
			
			BufferedImage img = new BufferedImage(nPlotWidth+50, 360, BufferedImage.TYPE_INT_RGB);
			Graphics2D graph = img.createGraphics();
			graph.setColor(Color.WHITE);
			graph.fillRect(0, 0, nPlotWidth+50, 360);
			m_App.m_plotFactory.BoxPlot(mapCountsToConditionsA, "junction: " + strJunctionNameA, graph, 40, 30, 300);
			
			Imagemap imgMap = new Imagemap();
			imgMap.setWidth(nPlotWidth+50 + "px");
			imgMap.setHeight("360px");
			imgMap.setContent(img);
			imgMap.setParent(grpBox);
			
			img = new BufferedImage(nPlotWidth+50, 360, BufferedImage.TYPE_INT_RGB);
			graph = img.createGraphics();
			graph.setColor(Color.WHITE);
			graph.fillRect(0, 0, nPlotWidth+50, 360);
			m_App.m_plotFactory.BoxPlot(mapCountsToConditionsB, "junction: " + strJunctionNameA, graph, 40, 30, 300);
			
			imgMap = new Imagemap();
			imgMap.setWidth(nPlotWidth+50 + "px");
			imgMap.setHeight("360px");
			imgMap.setContent(img);
			imgMap.setParent(grpBox);
			
			img = new BufferedImage(nPlotWidth+50, 360, BufferedImage.TYPE_INT_RGB);
			graph = img.createGraphics();
			graph.setColor(Color.WHITE);
			graph.fillRect(0, 0, nPlotWidth+50, 360);
			m_App.m_plotFactory.BoxPlot(mapRatiosToConditions, "Ratios", graph, 30, 30, 300);
			
			imgMap = new Imagemap();
			imgMap.setWidth(nPlotWidth+50 + "px");
			imgMap.setHeight("360px");
			imgMap.setContent(img);
			imgMap.setParent(grpBox);
			
			nPlotWidth = 400;
			
			img = new BufferedImage(nPlotWidth+50, 480, BufferedImage.TYPE_INT_RGB);
			graph = img.createGraphics();
			graph.setColor(Color.WHITE);
			graph.fillRect(0, 0, nPlotWidth+50, 480);
			Vector<Area> vcAreas = m_App.m_plotFactory.Plot2D(mapCountsToConditionsA, mapCountsToConditionsB, vcDataLabels, strJunctionNameA, strJunctionNameB, "2D plot", graph, 40, 30, nPlotWidth, 400, false);
			
			imgMap = new Imagemap();
			imgMap.setWidth(nPlotWidth+50 + "px");
			imgMap.setHeight("480px");
			imgMap.setContent(img);
			imgMap.setParent(grpBox);
			
			for(Area area : vcAreas)
			{
				area.setParent(imgMap);
			}
		}
	}
	
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
		
		if(bInit)
		{
			TreeSet<String> vcConditionsA = new TreeSet<String>();
			TreeSet<String> vcConditionsB = new TreeSet<String>();

			TreeSet<AnalysisResult> vcResults = m_resultHandler.GetAllTemporaryResults();

			for(AnalysisResult res : vcResults)
			{
				vcConditionsA.add(res.GetConditionA());
				vcConditionsB.add(res.GetConditionB());
			}
			
			TextFilter filter;
			filter = new TextFilter(vcConditionsA);
			m_resultFilterRule.SetTextFilter(ResultFilterRule.FILTER_TYPE_CONDITION_A, filter);
			filter = new TextFilter(vcConditionsB);
			m_resultFilterRule.SetTextFilter(ResultFilterRule.FILTER_TYPE_CONDITION_B, filter);
			
			TreeSet<String> vcStrings = new TreeSet<String>();
			vcStrings.add("true");
			vcStrings.add("false");
			vcStrings.add("NA");
			filter = new TextFilter(vcStrings);
			m_resultFilterRule.SetTextFilter(ResultFilterRule.FILTER_TYPE_SIGNIFICANT_PSI, filter);
			filter = new TextFilter(vcStrings);
			m_resultFilterRule.SetTextFilter(ResultFilterRule.FILTER_TYPE_NOVEL_JUNCTION, filter);
			
		}
		
		m_gridTmpHitList.setMold("default");
		
		// clear old rows
		m_gridTmpHitList.getChildren().clear();
		m_gridTmpHitList.setMold("paging");
		m_gridTmpHitList.setPageSize(10);
	
		//##################################
		//           add header
		//##################################
		Columns cols = new Columns();
	//	cols.setMenupopup("auto");
		
		// combobox for saving
		Column col = new Column("");
		col.setWidth("30px");
		col.setParent(cols);
	
		// view result
		col = new Column("");
		col.setWidth("40px");
		col.setParent(cols);
		
		// detailed view
		col = new Column("");
		col.setWidth("40px");
		col.setParent(cols);
		
		// detailed description
		col = new Column("");
		col.setWidth("40px");
		col.setParent(cols);
		
		col = new Column("Rating");
		col.setAlign("center");
		col.setSortAscending(getAscRatingComparator());
		col.setSortDescending(getDescRatingComparator());
		col.setWidth("110px");
		col.setParent(cols);
	
		col = new Column("Gene ID");
		col.setAlign("center");
		col.setWidth("140px");
		col.setParent(cols);
		
		col = new Column("Gene Symbol");
		col.setAlign("center");
		col.setWidth("100px");
		col.setParent(cols);
		
		col = new Column("Condition A");
		col.setAlign("center");
		col.setWidth("160px");
		col.setParent(cols);
		AddStringFilter(col, ResultFilterRule.FILTER_TYPE_CONDITION_A);
		
		col = new Column("Condition B");
		col.setAlign("center");
		col.setWidth("160px");
		col.setParent(cols);
		AddStringFilter(col, ResultFilterRule.FILTER_TYPE_CONDITION_B);
				
		col = new Column("");
		col.setSort("auto");
		col.setAlign("center");
		col.setWidth("170px");
		col.setParent(cols);
		Label lab = new Label("AS Type");
		lab.setParent(col);
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		AddSplicingTypeFilter(col);
		
		col = new Column("");
		col.setSort("auto");
		col.setWidth("120px");
		col.setAlign("center");
		col.setParent(cols);
		lab = new Label("Result Type");
		lab.setParent(col);
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		AddResultTypeFilter(col);
		
		col = new Column("");
		col.setSort("auto");
		col.setAlign("right");
		col.setWidth("210px");
		col.setParent(cols);
		lab = new Label("ratio change (exon A)");
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		lab.setParent(col);
		AddNumberFilter(col, ResultFilterRule.FILTER_TYPE_RATIO_CHANGE_A);
		
		col = new Column("");
		col.setSort("auto");
		col.setAlign("right");
		col.setWidth("210px");
		col.setParent(cols);
		lab = new Label("ratio change (exon B)");
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		lab.setParent(col);
		AddNumberFilter(col, ResultFilterRule.FILTER_TYPE_RATIO_CHANGE_B);
		
		col = new Column("");
		col.setSort("auto");
		col.setAlign("right");
		col.setWidth("150px");
		col.setParent(cols);
		lab = new Label("PSI change");
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		lab.setParent(col);
		AddNumberFilter(col, ResultFilterRule.FILTER_TYPE_PSI_CHANGE);
		
		col = new Column("");
		col.setSortAscending(getAscNumberComparator(14));
		col.setSortDescending(getDescNumberComparator(14));
		col.setSort("auto");
		col.setAlign("right");
		col.setWidth("210px");
		col.setParent(cols);
		lab = new Label("p-value ratio (exon A)");
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		lab.setParent(col);
		AddNumberFilter(col, ResultFilterRule.FILTER_TYPE_P_VALUE_RATIO_A);
		
		col = new Column("");
		col.setSortAscending(getAscNumberComparator(15));
		col.setSortDescending(getDescNumberComparator(15));
		col.setAlign("right");
		col.setWidth("210px");
		col.setParent(cols);
		lab = new Label("p-value ratio (exon B)");
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		lab.setParent(col);
		AddNumberFilter(col, ResultFilterRule.FILTER_TYPE_P_VALUE_RATIO_B);
		
		Column colPSI = new Column("");
		colPSI.setSortAscending(getAscNumberComparator(16));
		colPSI.setSortDescending(getDescNumberComparator(16));
		colPSI.setAlign("right");
		colPSI.setWidth("150px");
		colPSI.setParent(cols);
		lab = new Label("p-value PSI");
		lab.setParent(colPSI);
		lab.setStyle("align: left; margin-right: 50px; font-weight: bold;");
		AddNumberFilter(colPSI, ResultFilterRule.FILTER_TYPE_P_VALUE_PSI);
		
		col = new Column("significant");
		col.setSort("auto");
		col.setAlign("center");
		col.setWidth("150px");
		col.setParent(cols);
		AddStringFilter(col, ResultFilterRule.FILTER_TYPE_SIGNIFICANT_PSI);
		
		col = new Column("novel junction");
		col.setSort("auto");
		col.setAlign("center");
		col.setWidth("170px");		
		col.setParent(cols);
		AddStringFilter(col, ResultFilterRule.FILTER_TYPE_NOVEL_JUNCTION);
		
		col = new Column("event Location");
		col.setSort("auto");
		col.setAlign("left");
		col.setWidth("200px");
		col.setParent(cols);
		
		col = new Column("Comment");
		col.setWidth("410px");
		col.setAlign("left");
		col.setParent(cols);
	
		cols.setParent(m_gridTmpHitList);
		
		Rows rows = new Rows();
		rows.setParent(m_gridTmpHitList);
		
		//##################################
		//             add rows
		//##################################
		TreeSet<AnalysisResult> vcResults = m_resultHandler.GetAllTemporaryResults();
		
		for(AnalysisResult res : vcResults)
		{
			if(!m_resultFilterRule.bIsValidResult(res))
				continue;
			
			Row row = new Row();
			row.setId("TMPID_" + res.GetID());
			row.setValue(res.GetRating());
			row.setParent(rows);
			
			Checkbox checkbox = new Checkbox();
			checkbox.setValue(res.GetID());
			checkbox.setHeight("20px");
			checkbox.setWidth("20px");
			checkbox.setParent(row);
			
			if(res.GetType() == SplicingWebApp.AS_TYPE_RETAINED_INTRON)
			{
				Label label = new Label("");			
				label.setParent(row);
			}
			else
			{
				String strImageString = "/img/view.png";
				
				if(res.UsesNovelJunction())
					strImageString = "/img/view_incomplete.png";
				
				Image img = new Image(strImageString);
				img.setHeight("12px");
				img.setParent(row);
				img.setId("TMPview_" + res.GetID());
				img.setStyle("display: block; margin-left: auto; margin-right: auto; margin-top: 4px; margin-bottom: auto; cursor:hand; cursor:pointer;");
				
				img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
				{
					public void onEvent(Event event) throws Exception
					{
						Image img = (Image)event.getTarget();
	
						int nID = Integer.parseInt(img.getId().split("_")[1]);
						AnalysisResult res = m_resultHandler.GetTemporaryResult(nID);
						m_App.PrepareHitForVisualization(res);
						m_App.m_plotFactory.RequestCoverageRedraw();
						m_App.m_plotFactory.DrawPlots();					
					}
				});
			}
			
			if(res.HasPSIScore() || res.GetType() == SplicingWebApp.AS_TYPE_RETAINED_INTRON)
			{
				String strImageString = "/img/magnifier_blue.png";
				
				if(res.UsesNovelJunction())
					strImageString = "/img/magnifier_red.png";
				
				Image img = new Image(strImageString);
				img.setHeight("20px");
				img.setParent(row);
				img.setId("TMPdetailView_" + res.GetID());
				img.setStyle("display: block; margin-left: auto; margin-right: auto; margin-top: 4px; margin-bottom: auto; cursor:hand; cursor:pointer;");
				
				img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
				{
					public void onEvent(Event event) throws Exception
					{
						Image img = (Image)event.getTarget();
	
						int nID = Integer.parseInt(img.getId().split("_")[1]);
						
						AnalysisResult res = null;
						res = m_resultHandler.GetTemporaryResult(nID);
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
				label.setParent(row);
			}
			
			String strImageString = "/img/detailed_desc.png";

			Image img = new Image(strImageString);
			img.setWidth("16px");
			img.setHeight("20px");
			img.setParent(row);
			img.setId("TMPdetails_" + res.GetID());
			img.setStyle("display: block; margin-left: auto; margin-right: auto; margin-top: 4px; margin-bottom: auto; cursor:hand; cursor:pointer;");
			
			img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
			{
				public void onEvent(Event event) throws Exception
				{
					Image img = (Image)event.getTarget();

					int nID = Integer.parseInt(img.getId().split("_")[1]);
					AnalysisResult res = m_resultHandler.GetTemporaryResult(nID);
					m_App.PrepareHitForVisualization(res);
					ShowDetailPopup(res);
				}
			});
			
			Hlayout layoutH = new Hlayout();
			layoutH.setParent(row);
			for(int i=0; i<5; i++)
			{
				img = null;
				
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
				img.setId("TMPrating_" + (i+1) + "_" + res.GetID());
				img.setStyle("cursor:hand; cursor:pointer;");
				
				img.addEventListener(Events.ON_CLICK, new EventListener<Event>()
				{
					public void onEvent(Event event) throws Exception
					{
						Image img = (Image)event.getTarget();
						
						String strID = img.getId();
						
						String pSplit[] = strID.split("_");
						int nNewRating = Integer.parseInt(pSplit[1]);
						int nID = Integer.parseInt(pSplit[2]);
						
						AnalysisResult res = m_resultHandler.GetTemporaryResult(nID);
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
						Row row = (Row)img.getParent().getParent();
						row.setValue(nNewRating);
					}
				});
			}
			
			Label label = new Label(res.GetGeneID());
			label.setStyle("float: center");
			label.setParent(row);
			
			label = new Label(res.GetGeneSymbol());
			label.setStyle("float: center");
			label.setParent(row);
			
			label = new Label(res.GetConditionA());
			label.setStyle("float: center");
			label.setParent(row);
			
			label = new Label(res.GetConditionB());
			label.setStyle("float: center");
			label.setParent(row);
			
			label = new Label(res.GetTypeAsString());
			label.setStyle("float: center");
			label.setParent(row);
			
			label = new Label(res.GetResultTypeAsString());
			label.setStyle("float: center");
			label.setParent(row);
			
			if(res.HasAltExonA())
			{
				label = new Label(String.format(Locale.ENGLISH, "%.2f%%", res.GetAbsoluteChangeA()*100));
				if(Math.abs(res.GetAbsoluteChangeA()) < 0.05)
					label.setStyle("color: red");
				label.setStyle("float: right");
				label.setParent(row);
			}
			else
			{
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setStyle("float: right");
				label.setParent(row);
			}
			
			if(res.HasAltExonB())
			{
				label = new Label(String.format(Locale.ENGLISH, "%.2f%%", res.GetAbsoluteChangeB()*100));
				if(Math.abs(res.GetAbsoluteChangeB()) < 0.05)
					label.setStyle("color: red");
				label.setStyle("float: right");
				label.setParent(row);
			}
			else
			{
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setStyle("float: right");
				label.setParent(row);
			}
			
			if(res.HasPSIScore())
			{
				label = new Label(String.format(Locale.ENGLISH, "%.2f%%", res.GetInclusionChange()*100));
				if(res.GetInclusionChange() < 0.05)
					label.setStyle("color: red");
				label.setStyle("float: right");
				label.setParent(row);
			}
			else
			{
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setStyle("float: right");
				label.setParent(row);
			}
			
			if(res.HasAltExonA())
			{
				label = new Label(format.format(res.GetPValueA()));
				label.setStyle("float: right");
				label.setParent(row);				
			}
			else
			{
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setStyle("float: right");
				label.setParent(row);
			}
	
			if(res.HasAltExonB())
			{
				label = new Label(format.format(res.GetPValueB()));
				label.setStyle("float: right");
				label.setParent(row);				
			}
			else
			{
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setStyle("float: right");
				label.setParent(row);				
			}
	
			if(res.HasPSIScore())
			{
				label = new Label(format.format(res.GetPValuePSI()));
				label.setStyle("float: right");
				label.setParent(row);
				
				label = new Label(String.format(Locale.ENGLISH, "%s", res.HasSignificantPSIScore()));
				label.setStyle("align: center");
				label.setParent(row);
				
				label = new Label(String.format(Locale.ENGLISH, "%s", res.UsesNovelJunction()));
				label.setStyle("align: center");
				label.setParent(row);
			}
			else
			{
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setParent(row);
				
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setParent(row);
				
				label = new Label(String.format(Locale.ENGLISH, "NA"));
				label.setParent(row);
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
			
			label = new Label(String.format(Locale.ENGLISH, "%s:%,d - %,d", m_App.m_dataSupplier.GetReferenceName(), nEventStart, nEventEnd));
			label.setStyle("float: left");
			label.setParent(row);
			
			Textbox boxComment = new Textbox();
			boxComment.setParent(row);
			boxComment.setInplace(true);
			boxComment.setWidth("99%");
			boxComment.setStyle("margin-top: -5px; margin-bottom: -5px;");
			boxComment.setHflex("1");
			boxComment.setId("tmpCommentField_" + res.GetID());
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
					
					AnalysisResult res = m_resultHandler.GetTemporaryResult(nID);
					res.SetComment(box.getText());
				}
			});
		}
		
		colPSI.sort(true);
	}

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

	private Comparator<Row> getDescNumberComparator(int nIdx)
	{	
		final int m_nIdx = nIdx;
		
		return new Comparator<Row>()
		{
			@Override
			public int compare(Row arg0, Row arg1)
			{
				Label lab1 = (Label)arg0.getChildren().get(m_nIdx);
				Label lab2 = (Label)arg1.getChildren().get(m_nIdx);
				
				String strValue1 = lab1.getValue();
				String strValue2 = lab2.getValue();
				
				if(!strValue1.equals("NA") && !strValue2.equals("NA"))
				{
					double fValue1 = Double.parseDouble(strValue1);//new BigDecimal(strValue1).longValueExact();
					double fValue2 = Double.parseDouble(strValue2);//new BigDecimal(strValue2).longValueExact();
					
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
	
	private Comparator<Row> getAscNumberComparator(int nIdx)
	{
		final int m_nIdx = nIdx;
			
		return new Comparator<Row>()
		{
			@Override
			public int compare(Row arg0, Row arg1)
			{
				Label lab1 = (Label)arg0.getChildren().get(m_nIdx);
				Label lab2 = (Label)arg1.getChildren().get(m_nIdx);
				
				String strValue1 = lab1.getValue();
				String strValue2 = lab2.getValue();
				
				if(!strValue1.equals("NA") && !strValue2.equals("NA"))
				{
					double fValue1 = Double.parseDouble(strValue1);//new BigDecimal(strValue1).longValueExact();
					double fValue2 = Double.parseDouble(strValue2);//new BigDecimal(strValue2).longValueExact();
					
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
}
