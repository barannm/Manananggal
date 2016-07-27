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

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.nio.ByteBuffer;
import java.util.TreeSet;

import BioKit.Exon;

//################################################################################
//    The AnalysisResultHandler is used to fill the result list of the web
//    interface and handles request, such as adding comments or changing the
//    rating of an alternative splicing event.
//    It also handles the transfer of temporary results (= results in RAM)
//    to the permanent (= stored as physical file) result list
//################################################################################
public class AnalysisResultHandler
{
	private TreeSet<AnalysisResult> 	m_vcResults;			// stores results that were previously saved to disk
	private TreeSet<AnalysisResult> 	m_vcTemporaryResults;	// stores current results that are not yet saved to disk
	
	// stores results for the currently selected isoforms, these results won't be visible unless
	// a redraw of the temporary hit list is forced by the user (Button: "update temporary result list")
	private TreeSet<AnalysisResult>		m_vcTemporaryForCurrentSelectedIsoforms;
	
	private int							m_nMaxTmpID;				// the largest result identifier
	private int							m_nMaxPermanentID;
	
	AnalysisResultHandler()
	{
		m_vcResults 			= new TreeSet<AnalysisResult>();
		m_vcTemporaryResults 	= new TreeSet<AnalysisResult>();
		m_nMaxTmpID				= 0;
		m_nMaxPermanentID		= 0;
	}
	
	//################################################################################
	//    This function retrieves all alternative splicing events of type "nType"
	//    that involved an exon, given by reference name and exon object.
	//
	//    Currently, this function is only invoked during visualization. It is used
	//    to identify exons that refer to exon skipping events that will be shown
	//    by meta exons colored in red.
	//################################################################################
	public AnalysisResult GetASResultForExon(String strRef, Exon ex, int nType)
	{
		// search temporary results first
		if(m_vcTemporaryResults.size() > 0)
		{
			for(AnalysisResult res : m_vcTemporaryResults)
			{
				if(res.IncludesASExon(strRef, ex, nType))
					return res;
			}
		}
		else
		{
			for(AnalysisResult res : m_vcResults)
			{
				if(res.IncludesASExon(strRef, ex, nType))
					return res;
			}
		}
		
		return null;
	}
	
	//################################################################################
	//    This function retrieves all alternative splicing events of type "nType"
	//    that overlap a region specified by reference name, start and end position.
	//
	//    Currently, this function is only invoked when an highlighted meta exon
	//    is clicked to show a popup window with information on the AS event
	//################################################################################
	public TreeSet<AnalysisResult> GetASResultsForRegion(String strRef, int nStart, int nEnd, int nType)
	{
		TreeSet<AnalysisResult> vcRes = new TreeSet<AnalysisResult>();
		
		// search temporary results first
		if(m_vcTemporaryResults.size() > 0)
		{
			for(AnalysisResult res : m_vcTemporaryResults)
			{
				if(res.OverlapsRegion(strRef, nStart, nEnd, nType))
					vcRes.add(res);
			}
		}
		else
		{
			for(AnalysisResult res : m_vcResults)
			{
				if(res.OverlapsRegion(strRef, nStart, nEnd, nType))
					vcRes.add(res);
			}
		}
		
		return vcRes;
	}
	
	//################################################################################
	//    Sets the list of alternative splicing events detected for the current
	//    isoform selection. This function is called when a new gene is analyzed
	//    or when the isoform selection changed and the "Reanalyze"-Button was clicked.
	//################################################################################
	public void AddTemporaryResultsForCurrentlySelectedIsoforms(TreeSet<AnalysisResult> vcResults)
	{
		m_vcTemporaryForCurrentSelectedIsoforms = vcResults;
	}
	
	//################################################################################
	//    Transfers a selection (used defined) of splicing events from the temporary
	//    result list to the permanent result list. The function is invoked when the
	//    "Add selection to curated results" Button is clicked.
	//
	//    To make changes permanent, an additional click on the "Save Changes" Button
	//    is required.
	//################################################################################
	public void AddCuratedResults(TreeSet<AnalysisResult> vcResults)
	{
		for(AnalysisResult res : vcResults)
		{
			AnalysisResult newRes = new AnalysisResult(res);
			newRes.SetID(m_nMaxPermanentID);
			m_vcResults.add(newRes);
			m_nMaxPermanentID++;
		}
	}
	
	//################################################################################
	//    Transfers temporary results from the 'invisible'
    //    m_vcTemporaryForCurrentSelectedIsoforms list to the 'visible' 
	//    m_vcTemporaryResults list. During this process, all results get a new ID
	//    starting with the maximum ID of permanently stored results to avoid any
	//    duplication.
	//################################################################################
	public void UpdateTemporaryResults()
	{
		if(m_vcTemporaryForCurrentSelectedIsoforms == null)
			return;
		
		int nID = m_nMaxTmpID;
		for(AnalysisResult res : m_vcTemporaryForCurrentSelectedIsoforms)
		{
			nID++;
			res.SetID(nID);
			m_vcTemporaryResults.add(res);
		}
	}
	
	public void SaveToFile(String strPath, String strFileName)
	{
		try
		{
			FileOutputStream pOut = new FileOutputStream(strPath + "/" + strFileName);
			
			ByteBuffer bb = ByteBuffer.allocate(Integer.BYTES);
			bb.putInt(m_vcResults.size());
			pOut.write(bb.array());
			
			for(AnalysisResult res : m_vcResults)
				res.WriteToFile(pOut);
			
			pOut.close();
		}
		catch(Exception e)
		{
			System.out.println("failed to write alternative splicing results to file: " + strPath + "/" + strFileName);
			e.printStackTrace();
		}
	}
	
	// loads AS results from disk
	public boolean LoadHitList(String strPath, String strFileName)
	{
		try
		{
			FileInputStream pIn = new FileInputStream(strPath + "/" + strFileName);
			
			byte pBytes[] = new byte[Integer.BYTES];
			pIn.read(pBytes);
			ByteBuffer bb = ByteBuffer.wrap(pBytes);
			int nResults = bb.getInt();
			
			System.out.println("reading " + nResults + " saved results.");
			m_nMaxPermanentID = 0;
			
			for(int i=0; i<nResults; i++)
			{
				AnalysisResult res = new AnalysisResult();
				res.ReadFromFile(pIn);
				m_vcResults.add(res);
				
				m_nMaxPermanentID = Math.max(m_nMaxPermanentID, res.GetID());
			}
			m_nMaxPermanentID++;

			pIn.close();
		}
		catch(Exception e)
		{
			System.out.println("failed to write alternative splicing results to file: " + strPath + "/" + strFileName);
			e.printStackTrace();
		}
		return false;
	}
	
	public TreeSet<AnalysisResult> GetAllResults()
	{
		return m_vcResults;
	}
	
	public TreeSet<AnalysisResult> GetAllTemporaryResults()
	{
		return m_vcTemporaryResults;
	}
	
	public void RemoveResult(int nID)
	{
		for(AnalysisResult res : m_vcResults)
		{
			if(res.GetID() == nID)
			{
				m_vcResults.remove(res);
				return;
			}
		}
	}
	
	public AnalysisResult GetResult(int nID)
	{
		for(AnalysisResult res : m_vcResults)
		{
			if(res.GetID() == nID)
			{
				return res;
			}
		}
		return null;
	}
	
	public AnalysisResult GetTemporaryResult(int nID)
	{
		for(AnalysisResult res : m_vcTemporaryResults)
		{
			if(res.GetID() == nID)
			{
				return res;
			}
		}
		return null;
	}

	public void Clear()
	{
		m_vcTemporaryResults.clear();
		m_nMaxTmpID = 0;
		
		m_vcResults.clear();
		m_nMaxPermanentID = 0;
	}
	
	public void ClearTemporaryResults()
	{
		m_vcTemporaryResults.clear();
	}
}
