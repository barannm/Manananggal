package Manananggal;
import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.IntBuffer;
import java.nio.charset.Charset;
import java.util.Arrays;
import java.util.TreeMap;
import java.util.TreeSet;

import org.zkoss.zul.Messagebox;

import BioKit.Gene;

public class RandomJunctionReader
{
	private class Index implements Comparable<Index>
	{
		String	m_strKey;
		String 	m_strRef;
		int		m_nStart;
		int		m_nEnd;
		long	m_nOffset;
		
		@SuppressWarnings("unused")
		Index()
		{
			m_strKey 	= "?";
			m_nStart	= -1;
			m_nEnd		= -1;
			m_nOffset	= -1;
		}
		
		Index(String strRef, int nStart, int nEnd, long nOffset)
		{
			m_strRef	= strRef;
			m_nStart	= nStart;
			m_nEnd		= nEnd;
			m_nOffset	= nOffset;
		}
		
		@SuppressWarnings("unused")
		Index(String strKey, long nOffset)
		{
			m_strKey = strKey;
			m_nOffset = nOffset;
		}
		
		@Override
		public String toString()
		{
			return m_strKey + "_" + m_strRef + "_" + m_nStart + "_" + m_nEnd + " " + m_nOffset;
		}

//		@Override
		public int compareTo(Index other)
		{
			if(m_strRef.equals(other.m_strRef))
			{
				if(m_nStart == other.m_nStart)
				{
					return m_nEnd - other.m_nEnd;
				}
				
				return m_nStart - other.m_nStart;
			}

			return m_strRef.compareTo(other.m_strRef);
		}
	}
	
	@SuppressWarnings("unused")
	private class BinaryTree
	{
		String 				m_strKey;		// the name of the indexed element (either ref_start_end, or gene ID)
		TreeSet<Long> 		m_vcOffsets;	// there might be multiple offsets for a given key (i.e. genes will have multiple junctions)
		
		BinaryTree m_pLeft;
		BinaryTree m_pRight;
		
		BinaryTree()
		{
			m_strKey	= "?";
			m_vcOffsets = null;
			m_pLeft		= null;
			m_pRight	= null;
		}
		
		// the list of Indices MUST be sorted
		public void Init(Index pIndices[])
		{
			/*
			if(pIndices.length < 100)
				System.out.println(Arrays.toString(pIndices));
			System.out.println("------------------------------------");
			*/
			
			// check if leaf
			if(pIndices[0].m_strKey.equals(pIndices[pIndices.length-1].m_strKey))
			{
				m_strKey = pIndices[0].m_strKey;
								
				m_vcOffsets = new TreeSet<Long>();
				for(Index idx : pIndices)
				{
					m_vcOffsets.add(idx.m_nOffset);
				}
				
//				System.out.println("new leaf: " + m_strKey);
				//System.out.println(m_vcOffsets);
				
			}
			// just another node
			else
			{
				// get middle element
				int nMiddle = (int)Math.floor(pIndices.length *0.5);
				m_strKey = pIndices[nMiddle].m_strKey;
				
				// make sure to get the last entry with that key
				int i=0;
				while((nMiddle+i) < pIndices.length && pIndices[nMiddle+i].m_strKey.equals(m_strKey))
					i++;
				
				// new middle
				nMiddle += (i-1);

				// all elements less or equal to the middle element will be put in the left branch, the rest into the right
				Index pIndicesLeft[]  = Arrays.copyOfRange(pIndices, 0, nMiddle+1);
				Index pIndicesRight[] = Arrays.copyOfRange(pIndices, nMiddle+1, pIndices.length);

				// same list as before? then this is a node before a leaf!
				if(pIndicesLeft.length == pIndices.length || pIndicesRight.length == pIndices.length)
				{
					// sort all items by 'middle' key 
					TreeSet<String> vcKeys = new TreeSet<String>();
					
					TreeSet<Index> vcIndicesLeft  = new TreeSet<Index>();
					TreeSet<Index> vcIndicesRight = new TreeSet<Index>();
					
					for(Index idx : pIndices)
					{
						vcKeys.add(idx.m_strKey);
					}
					
					nMiddle = (int)Math.floor(vcKeys.size() * 0.5)-1;
					
					int j=0;
					for(Index idx : pIndices)
					{
						if(j == nMiddle)
							m_strKey = idx.m_strKey;
						j++;
					}
					
					for(Index idx : pIndices)
					{
						if(idx.m_strKey.compareTo(m_strKey) <= 0)
							vcIndicesLeft.add(idx);
						else
							vcIndicesRight.add(idx);
					}
					
					pIndicesLeft  = new Index[vcIndicesLeft.size() ];
					int nIdx = 0;
					for(Index idx : vcIndicesLeft)
					{
						pIndicesLeft[nIdx] = idx;
						nIdx++;
					}
					
					pIndicesRight = new Index[vcIndicesRight.size()];
					nIdx = 0;
					for(Index idx : vcIndicesRight)
					{
						pIndicesRight[nIdx] = idx;
						nIdx++;
					}
				}
				
				// add left branch
				if(pIndicesLeft.length > 0)
				{
					m_pLeft   = new BinaryTree();
					m_pLeft.Init(pIndicesLeft);
				}
				else
					m_pLeft = null;
				
				// add right branch
				if(pIndicesRight.length > 0)
				{
					m_pRight  = new BinaryTree();
					m_pRight.Init(pIndicesRight);
				}
				else
					m_pRight = null;
			}
		}
	
		public TreeSet<Long> Find(String strKey)
		{
			int nRes = strKey.compareTo(m_strKey);

			// check if this could be a leaf
			if(nRes == 0)
			{
				// leaves don't have any child
				if(this.m_pLeft == null && this.m_pRight == null)
				{
//					System.out.println("leaf reached");
					return m_vcOffsets;
				}
				else // otherwise this was not the real leaf
					return m_pLeft.Find(strKey);
			}
			else
			{
				if(nRes < 0)
				{
					if(m_pLeft == null)
					{
//						System.out.println("could not find key: " + strKey);
						return null;
					}
					
					return m_pLeft.Find(strKey);
				}
				else
				{
					if(m_pRight == null)
					{
//						System.out.println("could not find key: " + strKey);
						return null;
					}
					return m_pRight.Find(strKey);
				}
			}
		}
	}
	
	String m_strFileName;
//	BinaryTree m_treeIndex;
	TreeSet<Index> m_vcIndices;
	
	RandomJunctionReader()
	{
//		m_treeIndex = new BinaryTree();
		m_vcIndices = new TreeSet<Index>();
		m_strFileName = "?";
	}
	
	public void Init(String strFile) throws IOException
	{
		m_strFileName = strFile;
		BufferedInputStream pIn = new BufferedInputStream(new FileInputStream(strFile + ".idx"));

		while(true)
		{
			// read reference name
			byte pBytes[] = new byte[4];
			if(pIn.read(pBytes) == -1) break;
			ByteBuffer bb = ByteBuffer.wrap(pBytes);
			int nLength = bb.getInt();
			
			pBytes = new byte[nLength];
			if(pIn.read(pBytes) == -1) break;
			String strRef = new String(pBytes, Charset.forName("UTF-8"));
			
			// read start position
			pBytes = new byte[4];
			if(pIn.read(pBytes) == -1) break;
			bb = ByteBuffer.wrap(pBytes);
			int nStart = bb.getInt();
			
			// read end positions
			pBytes = new byte[4];
			if(pIn.read(pBytes) == -1) break;
			bb = ByteBuffer.wrap(pBytes);
			int nEnd = bb.getInt();
			
			// read file offset
			pBytes = new byte[8];			
			if(pIn.read(pBytes) == -1) break;
			bb = ByteBuffer.wrap(pBytes);
			long nOffset = bb.getLong();
			
//			Index idx = new Index(strRange, nOffset);
			Index idx = new Index(strRef, nStart, nEnd, nOffset);
			m_vcIndices.add(idx);	
		}

		pIn.close();
		
		System.out.println("read " + m_vcIndices.size() + " indices");
/*
		Index pIndices[] = new Index[m_vcIndices.size()];
		
		int nIdx = 0;
		for(Index idx : m_vcIndices)
		{
			pIndices[nIdx] = idx;
			nIdx +=1;
		}
		
		System.out.println("initializing index tree");
		m_treeIndex.Init(pIndices);
		System.out.println("..done");
*/
	}
	
	// strand is either -1 (-), 0 (?) or +1 (+)
	public TreeMap<String, TreeMap<String, Integer>> GetJunctionsForRange(String strRef, int nStart, int nEnd, int nStrand) throws IOException
	{	
		// for each junction (1st key = string) store the read counts (values) per sample (2nd key)
		TreeMap<String, TreeMap<String, Integer>> vcResult = new TreeMap<String, TreeMap<String, Integer>>();
		
		RandomAccessFile pIn = null;
		
		// open input file
		try
		{
			pIn = new RandomAccessFile(m_strFileName, "r");
		}
		catch(IOException ex)
		{
			Messagebox.show("failed to open file: " + m_strFileName);
			System.out.println("failed to open file: " + m_strFileName);
			System.out.println(ex.getMessage());
			return null;
		}
		
		if(pIn.length() == 0)
		{
			Messagebox.show("ERROR: Empty merged junction count file!");
			System.out.println("ERROR: Empty merged junction count file!");
			pIn.close();
			return null;
		}

		long nTotalFileLength = pIn.length();
		
		TreeMap<String, Integer> mapIdxToSamples = new TreeMap<String, Integer>();
		int nSamples = pIn.readInt();
		for(int i=0; i<nSamples; i++)
		{
			mapIdxToSamples.put(pIn.readUTF(), i);
		}

		idx_loop:
		for(Index idx : m_vcIndices)
		{
			if(!idx.m_strRef.equals(strRef))
				continue;
			
			if(idx.m_nEnd < nStart)
				continue;
			
			if(idx.m_nStart > nEnd)
				break;

			pIn.seek(idx.m_nOffset);
			
			while(pIn.getFilePointer() < nTotalFileLength)
			{
				String strCurRef = pIn.readUTF();
				int nCurStart = pIn.readInt();
				int nCurEnd   = pIn.readInt();
				
				pIn.readBoolean(); 	// skip known/unknown
				pIn.readUTF();		// skip geneID
				char chStrand = pIn.readChar();
				
				if(nCurStart > nEnd || !strCurRef.equals(strRef))
					break idx_loop;
				
				if(nCurStart <= nEnd && nCurEnd >= nStart)
				{
					TreeMap<String, Integer> mapCounts = new TreeMap<String, Integer>();

					byte pBuffer[] = new byte[mapIdxToSamples.keySet().size()*4];
					pIn.read(pBuffer);
					
					IntBuffer pBufferInt = ByteBuffer.wrap(pBuffer).asIntBuffer();
					int pCounts[] = new int[mapIdxToSamples.keySet().size()];
					pBufferInt.get(pCounts);
					
//					int pCounts[] = ByteBuffer.wrap(pBuffer).asIntBuffer().array();
					

					int i=0;
					for(String strSample : mapIdxToSamples.keySet())
					{
//						mapCounts.put(strSample, pIn.readInt());
						mapCounts.put(strSample, pCounts[i]);
						i++;
					}
				
					if(nStrand != 0)
					{
						if(chStrand == '?' || (nStrand == +1 && chStrand == '+') || (nStrand == -1 && chStrand == '-'))
							vcResult.put(nCurStart + "_" + nCurEnd, mapCounts);
					}
					else
					{
						vcResult.put(nCurStart + "_" + nCurEnd, mapCounts);
					}
				}
				else
				{
					// skip junction
					pIn.seek(pIn.getFilePointer() + mapIdxToSamples.keySet().size()*4);
				}
			}
		}

		pIn.close();
		
		return vcResult;
	}

	public TreeMap<CountElement, TreeMap<String, Integer>> GetJunctionsForRangeAsCountElements(String strRef, int nStart, int nEnd, int nStrand) throws IOException
	{	
		// for each junction (1st key = string) store the read counts (values) per sample (2nd key)
		TreeMap<CountElement, TreeMap<String, Integer>> vcResult = new TreeMap<CountElement, TreeMap<String, Integer>>();
		
		RandomAccessFile pIn = null;
		
		// open input file
		try
		{
			pIn = new RandomAccessFile(m_strFileName, "r");
		}
		catch(IOException ex)
		{
			System.out.println("failed to open file: " + m_strFileName);
			System.out.println(ex.getMessage());
			return null;
		}

		long nTotalFileLength = pIn.length();
		
		TreeMap<String, Integer> mapIdxToSamples = new TreeMap<String, Integer>();
		int nSamples = pIn.readInt();
		for(int i=0; i<nSamples; i++)
		{
			mapIdxToSamples.put(pIn.readUTF(), i);
		}

		idx_loop:
		for(Index idx : m_vcIndices)
		{
			if(!idx.m_strRef.equals(strRef))
				continue;
			
			if(idx.m_nEnd < nStart)
				continue;

			pIn.seek(idx.m_nOffset);
			
			while(pIn.getFilePointer() < nTotalFileLength)
			{
				String strCurRef = pIn.readUTF();
				int nCurStart = pIn.readInt();
				int nCurEnd   = pIn.readInt();
				pIn.readBoolean(); 	// skip known/unknown
				pIn.readUTF();		// skip geneID
				char chStrand = pIn.readChar();
				
				if(nCurStart > nEnd || !strCurRef.equals(strRef))
					break idx_loop;
				
				if(nCurStart <= nEnd && nCurEnd >= nStart)
				{
					TreeMap<String, Integer> mapCounts = new TreeMap<String, Integer>();

					byte pBuffer[] = new byte[mapIdxToSamples.keySet().size()*4];
					pIn.read(pBuffer);
					
					IntBuffer pBufferInt = ByteBuffer.wrap(pBuffer).asIntBuffer();
					int pCounts[] = new int[mapIdxToSamples.keySet().size()];
					pBufferInt.get(pCounts);

					int i=0;
					for(String strSample : mapIdxToSamples.keySet())
					{
						mapCounts.put(strSample, pCounts[i]);
						i++;
					}
					
					if(nStrand != 0)
					{
						if(chStrand == '?' || (nStrand == +1 && chStrand == '+') || (nStrand == -1 && chStrand == '-'))
						{
							CountElement e 	= new CountElement();
							e.m_nStart 		= nCurStart;
							e.m_nEnd 		= nCurEnd;
							e.m_chStrand 	= chStrand;
							vcResult.put(e,  mapCounts);
						}
					}
					else
					{
						CountElement e 	= new CountElement();
						e.m_nStart 		= nCurStart;
						e.m_nEnd 		= nCurEnd;
						e.m_chStrand 	= chStrand;
						vcResult.put(e,  mapCounts);
					}
				}
				else
				{
					// skip junction
					pIn.seek(pIn.getFilePointer() + mapIdxToSamples.keySet().size()*4);
				}
			}
		}

		pIn.close();
		
		return vcResult;
	}
	
	public TreeMap<String, TreeMap<String, Integer>> GetJunctionsForGene(Gene gene) throws IOException
	{
		TreeMap<String, TreeMap<String, Integer>> vcRes = null;
		
		if(gene.isPlusStrand())
			vcRes = GetJunctionsForRange(gene.getChromosome(), gene.getStart(), gene.getStop(), +1);
		else
			vcRes = GetJunctionsForRange(gene.getChromosome(), gene.getStart(), gene.getStop(), -1);
		
		return vcRes;
	}

}
