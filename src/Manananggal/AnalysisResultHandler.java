package Manananggal;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.nio.ByteBuffer;
import java.util.TreeSet;

import BioKit.Exon;

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
	
	public void AddTemporaryResultsForCurrentlySelectedIsoforms(TreeSet<AnalysisResult> vcResults)
	{
		m_vcTemporaryForCurrentSelectedIsoforms = vcResults;
	}
	
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
