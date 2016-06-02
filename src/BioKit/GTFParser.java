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
 * @author Fabian Birzele
 */
package BioKit;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.util.HashMap;
import java.util.TreeMap;

/**
 * GTF parser that allows to read a GTF file gene by gene. Results are returned wrapped in GTF Gene instances. Those can be converted
 * to standard Biokit genes using an inbuilt method. Please note that the parser assumes that all transcripts belonging to a gene are 
 * found in subsequent lines!
 * 
 * @author Fabian Birzele
 *
 */
public class GTFParser 
{
	private RandomAccessFile reader;
	private GTFGene newGene;
	
	// bStrict will only read lines of type "exon" as exons (and thus skips CDS entries)
	public GTFParser(String file) throws IOException
	{
		reader = new RandomAccessFile(file, "r");
	}
	
	public GTFParser(RandomAccessFile newReader, long nOffset) throws IOException
	{
		reader = newReader;
		reader.seek(nOffset);		
	}
	
	public GTFGene nextGene() throws IOException
	{
		GTFGene gene = newGene;
		newGene = null;
		
		String line = null;
		String[] split = null;
		String[] featuresSplit = null;
		String chromosome = null;
		String type = null;
		String geneID = null;
		String transcriptID = null;
		String synonyms = null;
		String geneName = null;
		String geneType = null;
		int start = -1;
		int stop = -1;
		boolean isPlusStrand = true;
		
		TreeMap<String, Integer> mapCDSstartToIsoform	= new TreeMap<String, Integer>();
		TreeMap<String, Integer> mapCDSendToIsoform		= new TreeMap<String, Integer>();

		while((line = reader.readLine()) != null)
		{
			//System.out.println("read line: " + line);
			split = line.split("\\t");
			
			try
			{
				// collect start codons
				if(split[2].equals("start_codon"))
				{
					// get position of the start codon
					int nPos = 0;
					if(split[6].equals("+"))
						nPos = Integer.parseInt(split[3]);
					else
						nPos = Integer.parseInt(split[4]);
					
					// get transcript ID
					featuresSplit = split[8].split(";");
					
					for(String f : featuresSplit)
					{
						if(f.trim().startsWith("transcript_id"))
						{
							transcriptID = f.replaceAll("transcript_id", "").replaceAll("\"", "").trim();
							mapCDSstartToIsoform.put(transcriptID, nPos);
							break;
						}
					}
					continue;
				}
				// collect stop codons
				else if(split[2].equals("stop_codon"))
				{				
					// get position of the stop codon
					int nPos = 0;
					if(split[6].equals("+"))
						nPos = Integer.parseInt(split[4]);
					else
						nPos = Integer.parseInt(split[3]);
					
					// get transcript ID
					featuresSplit = split[8].split(";");
					
					for(String f : featuresSplit)
					{
						if(f.trim().startsWith("transcript_id"))
						{
							transcriptID = f.replaceAll("transcript_id", "").replaceAll("\"", "").trim();
							mapCDSendToIsoform.put(transcriptID, nPos);
							break;
						}
					}
					continue;
				}
				// skip everything else but exons
				else if(!split[2].equals("exon"))
					continue;
				
				chromosome = split[0];
				type = split[1];
				start = Integer.parseInt(split[3]);
				stop = Integer.parseInt(split[4]);
				
				if(split[6].equals("-"))
					isPlusStrand = false;
				else if(split[6].equals("+"))
					isPlusStrand = true;
				
				featuresSplit = split[8].split(";");
				
				for(String f : featuresSplit)
				{
					if(f.trim().startsWith("gene_id"))
					{
						geneID = f.replaceAll("gene_id", "").replaceAll("\"", "").trim();
					}
					else if(f.trim().startsWith("transcript_id"))
					{
						transcriptID = f.replaceAll("transcript_id", "").replaceAll("\"", "").trim();
					}
					else if(f.trim().startsWith("synonyms"))
					{
						synonyms = f.replaceAll("synonyms", "").replaceAll("\"", "").trim();
					}

					// add gene symbol if available
					String[] split2 = f.trim().split("\\s");
					if(split2.length == 2 && split2[0].equals("gene_name"))
					{
						geneName = split2[1].replaceAll("\"", "").trim();
					}
					
					// add gene type if available
					if(split2.length == 2 && split2[0].equals("gene_type"))
					{
						geneType = split2[1].replaceAll("\"", "").trim();
					}
				}
			
				// first line
				if(gene == null)
				{
					if(geneID.equals(""))
						geneID = transcriptID;

					gene = new GTFGene(chromosome, isPlusStrand, geneID, geneName, geneType, type);
					gene.setSynonymInformation(synonyms);
					gene.addExon(transcriptID, start, stop);
				}
				else if(geneID.equals(gene.getId()) && gene.getChromosome().equals(chromosome))
				{
					gene.addExon(transcriptID, start, stop);
				}
				else
				{
					// I believe everything but the break is not required here, because you set newGene to null each time you call nextGene() anyway <comment: Matthias>
					if(geneID.equals(""))
						geneID = transcriptID;
					
					newGene = new GTFGene(chromosome, isPlusStrand, geneID, geneName, geneType, type);
					newGene.setSynonymInformation(synonyms);
					newGene.addExon(transcriptID, start, stop);
					
					// it is a different gene, thus break
					break;
				}
			}
			catch(Exception ex)
			{
				System.err.println("GTFParser could not parse line: " + line);
			}
		}

		// add cds start and end positions to gene information
		gene.SetCDSPositions(mapCDSstartToIsoform, mapCDSendToIsoform);
		return gene;
	}
	
	public void closeReader()
	{
		try 
		{
			this.reader.close();
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
			System.exit(0);
		}
	}
	
	public static void main(String[] args) throws Exception
	{
		if(args.length == 0)
		{
			System.out.println("-prepareUCLCGTFFile transcriptToGeneMapping inputGTF outputGTFFile");
			System.exit(0);
		}
		
		if(args[0].equals("-prepareUCLCGTFFile"))
		{
			RandomAccessFile reader = new RandomAccessFile(args[1], "r");
			String line = null;
			
			HashMap<String, String> transcriptToGeneMap = new HashMap<String, String>();
			
			while((line = reader.readLine()) != null)
			{
				String[] split = line.split("\t");
				
				transcriptToGeneMap.put(split[0], split[1]);
			}
			
			reader.close();
			
			System.out.println("Transcript to Gene Map size: " + transcriptToGeneMap.size());
			
			reader = new RandomAccessFile(args[2], "r");
			
			PrintWriter writer = new PrintWriter(new FileWriter(args[3]));
			
			line = null;
			String[] split = null;
			String[] featuresSplit = null;
			String chromosome = null;
			String type = null;
			String transcriptID = null;
			int start = -1;
			int stop = -1;
			
			while((line = reader.readLine()) != null)
			{
				//System.out.println("read line: " + line);
				split = line.split("\\t");
				
				// only read exon lines
				if(!split[2].equals("exon"))
					continue;
				
				chromosome = split[0];
				type = split[1];
				start = Integer.parseInt(split[3]);
				stop = Integer.parseInt(split[4]);
				
				/*
				if(split[6].equals("-"))
					isPlusStrand = false;
				else if(split[6].equals("+"))
					isPlusStrand = true;
				*/
				featuresSplit = split[8].split(";");
				
				for(String f : featuresSplit)
				{
					if(f.trim().startsWith("gene_id"))
					{
//						geneID = f.replaceAll("gene_id", "").replaceAll("\"", "").trim();
					}
					else if(f.trim().startsWith("transcript_id"))
					{
						transcriptID = f.replaceAll("transcript_id", "").replaceAll("\"", "").trim();
					}
					else if(f.trim().startsWith("synonyms"))
					{
//						synonyms = f.replaceAll("synonyms", "").replaceAll("\"", "").trim();
					}
					
				}
				
				if(transcriptToGeneMap.containsKey(transcriptID))
					writer.println(chromosome + "\t" + type + "\texon\t"  + start + "\t" + stop + "\t"+split[5]+"\t" + split[6] 
					    + "\t" + split[7] + "\tgene_id\"" + transcriptToGeneMap.get(transcriptID) + "\"; transcript_id \"" 
					    + transcriptID + "\"");
			}
			reader.close();
			
			writer.flush();
			writer.close();
			
		} 
		
	}

	
}
