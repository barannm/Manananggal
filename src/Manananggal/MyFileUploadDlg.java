package Manananggal;

import java.util.Collections;
import java.util.Scanner;
import java.util.Vector;
import java.util.zip.GZIPInputStream;

import org.zkoss.util.media.Media;
import org.zkoss.zk.ui.event.Event;
import org.zkoss.zk.ui.event.EventListener;
import org.zkoss.zk.ui.event.Events;
import org.zkoss.zk.ui.event.UploadEvent;
import org.zkoss.zul.Button;
import org.zkoss.zul.Hlayout;
import org.zkoss.zul.Label;
import org.zkoss.zul.Messagebox;
import org.zkoss.zul.Radio;
import org.zkoss.zul.Radiogroup;
import org.zkoss.zul.Vlayout;
import org.zkoss.zul.Window;

public class MyFileUploadDlg
{	
	SplicingWebApp m_App;
	Window	m_windowPopup;
	Label	m_labSelectedFile;
	Label	m_labError;
	Media	m_media;
	int		m_nSelectedFileType;
	
	public MyFileUploadDlg(SplicingWebApp app)
	{
		m_App 	= app;
		m_media = null;
		m_nSelectedFileType = 0;
		
		m_windowPopup = new Window();
		m_windowPopup.setParent(m_App);
		m_windowPopup.setId("upldDlg");
		
		m_windowPopup.doPopup();
		m_windowPopup.setTitle("coverage");
		m_windowPopup.setSizable(false);
		m_windowPopup.setClosable(false);
		m_windowPopup.setMinimizable(false);
		m_windowPopup.setMaximizable(false);
		m_windowPopup.setBorder(true);
		
		m_windowPopup.setVflex("0");
		m_windowPopup.setHflex("0");

		m_windowPopup.setPosition("center,center");
		m_windowPopup.setVisible(true);
		
		m_windowPopup.setWidth (600 + "px");
		m_windowPopup.setHeight(400 + "px");
		
		Vlayout layoutV = new Vlayout();
		layoutV.setParent(m_windowPopup);
		
		Hlayout layoutH = new Hlayout();
		layoutH.setParent(layoutV);
		
		Button btnSelectFile = new Button("Select File");
		btnSelectFile.setId("btnUpload");
		btnSelectFile.setUpload("true");
		btnSelectFile.setParent(layoutH);
		
		m_labSelectedFile = new Label();
		m_labSelectedFile.setParent(layoutH);
		m_labSelectedFile.setStyle("margin-left: 10px; margin-top: 5px; display: inline-block;");
				
		Radiogroup grp = new Radiogroup();
		grp.setParent(layoutV);
		
		Radio radioBtn = new Radio("auto");		
		radioBtn.setParent(grp);
		radioBtn.setValue(0);
		radioBtn.setChecked(true);
		
		radioBtn = new Radio("Manananggal");
		radioBtn.setParent(grp);
		radioBtn.setValue(1);
		radioBtn.setStyle("margin-left: 10px");
		
		radioBtn = new Radio("rMATS");
		radioBtn.setParent(grp);
		radioBtn.setValue(2);
		radioBtn.setStyle("margin-left: 10px");
		
		radioBtn = new Radio("DEXSeq");
		radioBtn.setParent(grp);
		radioBtn.setValue(3);
		radioBtn.setStyle("margin-left: 10px");
		
		radioBtn = new Radio("JSplice");
		radioBtn.setParent(grp);
		radioBtn.setValue(4);
		radioBtn.setStyle("margin-left: 10px");
		
		radioBtn = new Radio("Cuffdiff");
		radioBtn.setParent(grp);
		radioBtn.setValue(5);
		radioBtn.setStyle("margin-left: 10px");
		
		Button btnImport = new Button("Import");
		btnImport.setParent(layoutV);
		
		m_labError = new Label();
		m_labError.setId("labError");
		m_labError.setParent(layoutV);
		m_labError.setStyle("color: red;");
		
		grp.addEventListener(Events.ON_CHECK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				Radio radioBtn = (Radio)event.getTarget();
				m_nSelectedFileType = radioBtn.getValue();
			}
		});
		
		btnSelectFile.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				// clear error state
				m_labError.setValue("");
				m_labError.invalidate();
				
				// clear selected file label
				m_labSelectedFile.setValue("");				
				m_labSelectedFile.invalidate();
				
				// clear previous data
				m_media = null;
			}
		});
		
		btnSelectFile.addEventListener(Events.ON_UPLOAD, new EventListener<UploadEvent>()
		{
			@Override
			public void onEvent(UploadEvent event) throws Exception
			{
				m_media = event.getMedia();				
				m_labSelectedFile.setValue(m_media.getName());
			}
		});
		
		btnImport.addEventListener(Events.ON_CLICK, new EventListener<Event>()
		{
			@Override
			public void onEvent(Event event) throws Exception
			{
				if(m_media == null)
				{
					m_labError.setValue("No file selected!");
					return;
				}
				
				Scanner pIn = null;
				if(m_media.getName().endsWith(".gz"))
				{
					GZIPInputStream pGZin = new GZIPInputStream(m_media.getStreamData());
					pIn = new Scanner(pGZin);	
				}
				else
				{
					pIn = new Scanner(m_media.getReaderData());
				}
				
				//#############################
				//      detect file type
				//#############################
				String strLine = "";
				
				// skip empty lines and comments
				while(strLine.isEmpty() || strLine.startsWith("#"))
					strLine = pIn.nextLine();
				
				int nFileType = -1;
				
				// create a map that stores the column names
				Vector<String> vcColumnNames = new Vector<String>();
				
				// process the header but also check the file type
				// please note that the manananggal header won't be processed here
				nFileType = DetectFileType(strLine, vcColumnNames);
				
				String[] pColumnNames = new String[vcColumnNames.size()];
				vcColumnNames.toArray(pColumnNames);

				// overwrite the file type if something specific was selected
				if(m_nSelectedFileType != 0)
					nFileType = m_nSelectedFileType;
				
				if(nFileType == -1)
				{
					m_labError.setValue("Could not determine file type!");
					pIn.close();
					return;
				}

				// pass the input file and the first line to the processing function
				ProcessData(pIn, strLine, nFileType, pColumnNames);
				pIn.close();

				m_windowPopup.setVisible(false);
			}
		});

	}

	public void Show()
	{		
		m_windowPopup.doPopup();
		m_windowPopup.setTopmost();
		m_windowPopup.setVisible(true);
		
		// clear error state
		m_labError.setValue("");
		m_labError.invalidate();
		
		// clear selected file label
		m_labSelectedFile.setValue("");				
		m_labSelectedFile.invalidate();
		
		// clear previous data
		m_media = null;
	}
	
	public int DetectFileType(String strLine, Vector<String> vcColumnNames)
	{
		String[] pSplit = strLine.split("\t");
		
		if(pSplit.length < 4)
			return -1;
		
		// store column names and column indices
		for(int i=0; i<pSplit.length; i++)
		{
			// for JSplice input, check whether the bugged column headers have been corrected in the meantime (Largest_relFCs and adjPvalues were merged into a single column header)
			if(pSplit[i].startsWith("Largest_relFCs") && !pSplit[i].equals("Largest_relFCs"))
			{
				vcColumnNames.add("Largest_relFCs");
				vcColumnNames.add("adjPvalues");
			}
			else
				vcColumnNames.add(pSplit[i]);
		}
		
		// Manananggal
		if(pSplit[0].startsWith("ID_") && (pSplit[3].equals("ratio_only") || pSplit[3].equals("ratio_only") || pSplit[3].equals("split_read_only") || pSplit[3].equals("combined") ))
			return SplicingWebApp.IMPORTED_DATA_TYPE_MANA;
		
		// rMATS
		if(pSplit[0].equals("ID") && pSplit[1].equals("GeneID") && pSplit[2].equals("geneSymbol"))
			return SplicingWebApp.IMPORTED_DATA_TYPE_RMATS;
		
		// DEXSeq - Check twice, because maybe the file includes row names
		if((pSplit[0].equals("groupID") && pSplit[1].equals("featureID")) || (pSplit[1].equals("groupID") && pSplit[2].equals("featureID")))
			return SplicingWebApp.IMPORTED_DATA_TYPE_DEXSEQ;
		
		// JSplice file
		if(pSplit[0].equals("Gene_name|ID") && pSplit[1].equals("ASM_type"))
			return SplicingWebApp.IMPORTED_DATA_TYPE_JSPLICE;
		
		// Cuffdiff
		if(pSplit[0].equals("test_id") && pSplit[1].equals("gene_id"))
			return SplicingWebApp.IMPORTED_DATA_TYPE_CUFFDIFF;
		
		return -1;
	}

	public void ProcessData(Scanner pIn, String strLine, int nFileType, String[] pColumnNames) throws ClassNotFoundException, InstantiationException, IllegalAccessException
	{
		Vector<AnalysisResult> vcResults = new Vector<AnalysisResult>();
		
		// also process the first line if it's manananggal output
		// because the header included '#' and was skipped
		if(nFileType == SplicingWebApp.IMPORTED_DATA_TYPE_MANA)
		{
			AnalysisResult res = new AnalysisResult();
			res.ParseManananggalOutput(strLine);
			vcResults.add(res);
		}
		
		int nID = 0;
		while_loop : while(pIn.hasNextLine())
		{
			strLine = pIn.nextLine();
			
			//## skip empty lines
			if(strLine.trim().isEmpty())
				continue;
			
			AnalysisResult res = new AnalysisResult();
			
			switch(nFileType)
			{
				case SplicingWebApp.IMPORTED_DATA_TYPE_MANA:
				{
					if(strLine.startsWith("finished"))
						break while_loop;
					res.ParseManananggalOutput(strLine);
					break;
				}
				
				case SplicingWebApp.IMPORTED_DATA_TYPE_DEXSEQ:
				{
					res.ParseDEXSeqOutput(strLine, nID, pColumnNames);
					nID++;
					break;
				}
				
				case SplicingWebApp.IMPORTED_DATA_TYPE_RMATS:
				{
					res.ParseMatsOutput(strLine, pColumnNames);
					break;
				}
				
				case SplicingWebApp.IMPORTED_DATA_TYPE_JSPLICE:
				{
					res.ParseJSpliceOutput(strLine, nID, pColumnNames);
					nID++;
					break;
				}
				
				case SplicingWebApp.IMPORTED_DATA_TYPE_CUFFDIFF:
				{
					res.ParseCuffdiffOutput(strLine, pColumnNames);
					break;
				}
			}

			if(res.GetID() != -1)
				vcResults.add(res);
		}
		
		if(vcResults.size() > 2000)
		{
			Messagebox.show("Detected more than 2000 results. Therefore, the result list was reduced to 2000 results based on the p-value.");
			
			Collections.sort(vcResults, AnalysisResult.ImportedPValueComparator);
			Vector<AnalysisResult> vcNewList = new Vector<AnalysisResult>();
			vcNewList.addAll(vcResults.subList(0, 1999));
			
			m_App.GetResultListHandler().UpdateImportedResultList(vcNewList, nFileType, true); 
			
		}
		else
		{
			m_App.GetResultListHandler().UpdateImportedResultList(vcResults, nFileType, true);
		}
	}
}