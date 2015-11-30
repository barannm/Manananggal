package Manananggal;

import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import BioKit.Exon;

// This class stores are tree where each node corresponds to either a junction or exon
// and points to x other exon/junctions

public class JunctionTree
{
	class JunctionNode
	{
		CountElement			m_CurElement;
		TreeSet<CountElement> 	m_vcConnections;

		JunctionNode()
		{
			m_CurElement = new CountElement();
			m_vcConnections = new TreeSet<CountElement>();
		}
		
		JunctionNode(CountElement curElement, TreeSet<CountElement> vcConnections)
		{
			m_CurElement	= curElement;
			m_vcConnections = vcConnections;
		}
		
		Vector<TreeSet<CountElement>> GetPath(int nPathEnd, TreeMap<CountElement, JunctionNode> mapNodes, TreeSet<CountElement> vcDeadEnds)
		{
			/*
			if(nPathEnd == 35236399)
				System.out.println("current: " + m_CurElement);
			*/
			
			Vector<TreeSet<CountElement>> vcRes = new Vector<TreeSet<CountElement>>();
			
			// this node has already been proven to be a dead end
			if(vcDeadEnds.contains(m_CurElement))
				return vcRes;
			
			if(m_vcConnections.size() == 0 || (m_CurElement.m_bIsExon && m_CurElement.m_nStart == nPathEnd))
			{
				if(m_CurElement.m_bIsExon && m_CurElement.m_nStart == nPathEnd)
				{
					TreeSet<CountElement> vcPaths = new TreeSet<CountElement>();
					vcPaths.add(m_CurElement);
					vcRes.add(vcPaths);
					
					/*
					if(nPathEnd == 35236399)
						System.out.println("end detected: " + m_CurElement);
					*/
					
					return vcRes;
				}
				else
				{
					if(nPathEnd == 35236399)
						System.out.println("dead end: " + m_CurElement);
					vcDeadEnds.add(m_CurElement);
					return vcRes;
				}
			}
			else
			{
				boolean bHasValidEnds = false;
				for(CountElement e : m_vcConnections)
				{
					if(vcDeadEnds.contains(e))
						continue;
						
					Vector<TreeSet<CountElement>> vcTmp = mapNodes.get(e).GetPath(nPathEnd, mapNodes, vcDeadEnds);
					
					if(vcTmp.size() == 0)
					{
						/*
						if(nPathEnd == 35236399)
							System.out.println("invalid node: " + m_CurElement);
						*/
						vcDeadEnds.add(e);
						continue;
					}
					
					bHasValidEnds = true;
					
					// process each subsequent path
					for(TreeSet<CountElement> vcPath : vcTmp)
					{
						TreeSet<CountElement> vcNew = new TreeSet<CountElement>();
						
/*
						if(nPathEnd == 35236399)
							System.out.println("adding: " + m_CurElement + " linked to " + m_vcConnections);
	*/					
						vcNew.add(m_CurElement);
						
						for(CountElement ce : vcPath)
							vcNew.add(ce);
						
						vcRes.add(vcNew);
					}
				}
				
				if(!bHasValidEnds)
					vcRes.clear();
				
				return vcRes;
			}
		}
	}

	TreeMap<CountElement, JunctionNode> 				m_mapConnections;		// maps each element to a node
	TreeMap<CountElement, TreeMap<String, Integer>> 	m_mapJunctionCounts;	// stores the junction counts
	
	JunctionTree()
	{
		m_mapConnections 	= new TreeMap<CountElement, JunctionNode>();
		m_mapJunctionCounts	= null;
	}
	
	JunctionTree(TreeMap<CountElement, TreeMap<String, Integer>> mapJunctionCounts, Exon[] pExons)
	{
		m_mapConnections 	= new TreeMap<CountElement, JunctionNode>();
		m_mapJunctionCounts = mapJunctionCounts;
		
		// add new junctions
		for(CountElement jun : mapJunctionCounts.keySet())
		{
			jun.m_bIsExon = false;
			
			// get all exons connected to the new junctions
			TreeSet<CountElement> vcDownStreamExons = new TreeSet<CountElement>();
			TreeSet<CountElement> vcUpStreamExons 	= new TreeSet<CountElement>();
			
			for(Exon ex : pExons)
			{
				CountElement e 	= new CountElement();
				e.m_nStart		= ex.getCodingStart();
				e.m_nEnd		= ex.getCodingStop();
				e.m_bIsExon		= true;
				
				if(e.m_nStart == jun.m_nEnd)
					vcDownStreamExons.add(e);
				
				if(e.m_nEnd == jun.m_nStart)
					vcUpStreamExons.add(e);
			}
			
			// add junction
			JunctionNode node = new JunctionNode(jun, vcDownStreamExons);
			node.m_CurElement = jun;
			m_mapConnections.put(jun, node);
			
			// add all downstream exons
			for(CountElement ex : vcDownStreamExons)
			{
				if(!m_mapConnections.containsKey(ex))
				{
					node = new JunctionNode(ex, new TreeSet<CountElement>());
					node.m_CurElement = ex;
					m_mapConnections.put(ex, node);
				}
			}
			
			// add all upstream exons
			for(CountElement ex : vcUpStreamExons)
			{
				if(m_mapConnections.containsKey(ex))
				{
					m_mapConnections.get(ex).m_vcConnections.add(jun);
				}
				else
				{
					TreeSet<CountElement> vcJunctions = new TreeSet<CountElement>();
					vcJunctions.add(jun);
					node = new JunctionNode(ex, vcJunctions);
					node.m_CurElement = ex;

					m_mapConnections.put(ex, node);
				}
			}
		}
	}

	public void GetPaths(int nPathStart, int nPathEnd, TreeMap<Integer, TreeSet<CountElement>> pathFragmentsA, TreeMap<Integer, TreeSet<CountElement>> pathFragmentsB)
	{
		Vector<TreeSet<CountElement>> res = new Vector<TreeSet<CountElement>>();
		
		// get start exon
		CountElement startExon = null; 
		for(CountElement e : m_mapConnections.keySet())
		{
			if(e.m_bIsExon && e.m_nEnd == nPathStart)
				startExon = e;
		}
		
		if(nPathStart == 35211612 && nPathEnd == 35236399)
			System.out.println("----- starting at: " + startExon);
		else
			return;
		
		if(startExon == null)
		{
			System.out.println("failed to get path for: " + nPathStart + " " + nPathEnd);
			return;
		}

		res = m_mapConnections.get(startExon).GetPath(nPathEnd, m_mapConnections, new TreeSet<CountElement>());

/*
		System.out.println("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
		System.out.println(res.get(0));
		System.out.println("-----------------------------------------------------------");
		System.out.println(res.get(1));
		System.out.println(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
*/	
		// post process paths if there are more than 50
		// split paths by common junctions
		// i.e. A-B-C-D-E-F vs. A---C-D---F becomes
		// A-B-C vs. A---C and D-E-F vs. D---F
//		if(res.size() > 50)
		{
			int nPathA = 0;
			
			for(TreeSet<CountElement> pathA : res)
			{
				int nPathB = 0;
				for(TreeSet<CountElement> pathB : res)
				{
					// skip equal paths
					if(nPathA != nPathB)
					{
						CountElement pElementsA[] = new CountElement[pathA.size()];
						pathA.toArray(pElementsA);
						
						CountElement pElementsB[] = new CountElement[pathB.size()];
						pathB.toArray(pElementsB);
						
						TreeMap<Integer, Integer> mapSharedElementIDs = new TreeMap<Integer, Integer>();
						int nJunctionsA = 0;
						for(int i=0; i<pElementsA.length; i++)
						{
							if(pElementsA[i].m_bIsExon)
								continue;
							
							nJunctionsA++;
							
							for(int j=0; j<pElementsB.length; j++)
							{
								if(pElementsA[i] == pElementsB[j])
									mapSharedElementIDs.put(i, j);
							}
						}
						// skip equal paths
						if(mapSharedElementIDs.size() == nJunctionsA)
							continue;
						
						System.out.println("junctions: " + nJunctionsA);
						System.out.println(nPathA + " -> " + Arrays.toString(pElementsA));
						System.out.println(nPathB + " -> " + Arrays.toString(pElementsB));
						System.out.println(mapSharedElementIDs);
						try {
							Thread.sleep(1000);
						} catch (InterruptedException e1) {
							// TODO Auto-generated catch block
							e1.printStackTrace();
						}
						
						int nLastA = 0;
						int nLastB = 0;
						for(Map.Entry<Integer, Integer> e : mapSharedElementIDs.entrySet())
						{
							int nID_A = e.getKey();
							int nID_B = e.getValue();
							
							TreeSet<CountElement> vcTmpA = new TreeSet<CountElement>();
							TreeSet<CountElement> vcTmpB = new TreeSet<CountElement>();
							
							String strID_A = "";
							String strID_B = "";
							
							for(int i=nLastA; i<nID_A; i++)
							{
								vcTmpA.add(pElementsA[i]);
								strID_A += pElementsA[i].m_nStart + "_" + pElementsA[i].m_nEnd; 
							}
							
							for(int i=nLastB; i<nID_B; i++)
							{
								vcTmpB.add(pElementsB[i]);
								strID_B += pElementsB[i].m_nStart + "_" + pElementsB[i].m_nEnd;
							}

							if(vcTmpA.size() > 1 && vcTmpB.size() > 1)
							{
								pathFragmentsA.put(strID_A.hashCode(), vcTmpA);
								pathFragmentsB.put(strID_B.hashCode(), vcTmpB);
							}
							else
							{
								/*
								System.out.println("unknown error");
								System.out.println("A: " + vcTmpA);
								System.out.println("B: " + vcTmpB);
								*/
							}
							
							nLastA = nID_A;
							nLastB = nID_B;
						}
					}
					nPathB++;
				}
				nPathA ++;
			}
		}
		/*
		System.out.println("fragments A: " + pathFragmentsA);
		System.out.println("fragments B: " + pathFragmentsB);
		*/
	}
	
	@Override
	public String toString()
	{
		String strOut = "";
		for(CountElement e : m_mapConnections.keySet())
		{
			strOut += e.toString() + "\n";
			strOut += "connected to:\n";
			for(CountElement ce : m_mapConnections.get(e).m_vcConnections)
			{
					strOut += "\t" + ce + " " + m_mapJunctionCounts.get(ce) + "\n";
			}
		}
		
		return(strOut);
	}
	
	public TreeMap<String, Integer> GetCountsForElement(CountElement e)
	{
		return(m_mapJunctionCounts.get(e));
	}
}