package BioKit;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashMap;

public class MMSeqControl 
{
	private String mmseqPath;
	private String tempDirectory;

	public MMSeqControl(String mmseqPath, String tempDirectory)
	{
		this.mmseqPath = mmseqPath;
		this.tempDirectory = tempDirectory;
	}
	
	public File generateSAMFileName(String geneID)
	{
		return new File(tempDirectory + File.separatorChar + geneID + "_" + (int)(Math.random()*100000) + ".sam");
	}
	
	public HashMap<String, Double> runMMSeq(Gene gene, File samFile, int readLength)
	{
		try
		{
			HashMap<String, Double> isoformLevels = new HashMap<String, Double>();
			
			System.out.println("MMSeq: write gtf file...");
			
			// write gtf file
			File gtfFile = new File(tempDirectory + File.separatorChar + gene.getGeneID() + ".gtf");
			
			for(int i=0; i<10; i++)
			{
				PrintWriter writer = new PrintWriter(new FileWriter(gtfFile));
				
				for(String transcript : gene.getArrayOfGeneProductNames())
				{
					Exon[] exons = (Exon[])gene.getGeneProductForName(transcript).getProperty(Gene.PROPERTY_EXONS);
					
					for(Exon e : exons)
					{
						writer.println(gene.getChromosome() + "\tProSAS\texon\t" + e.getGenomicStart() + "\t" + e.getGenomicStop() + "\t.\t" 
							+ gene.getStrandStringID() + "\t.\t" + "gene_id \"" + gene.getGeneID() + "\"; transcript_id \"" + transcript + "\";");
					}
				}
				
				writer.flush();
				writer.close();
				
				if(gtfFile.length() > 0)
				{
					System.out.println("GTF File contains data..." + gtfFile.length());
					break;
				}
				else
				{
					System.out.println("GTF File is empty... make another attempt " + (i+2) + " out of 10");
				}
			}
			
			
			System.out.println("MMseq: create hits file...");
			
			// generate hits file
			File hitsFile = new File(tempDirectory + File.separatorChar + gene.getGeneID() + "_" + (int)(Math.random()*100000) + ".hits");
			
			for(int i=0; i<10; i++)
			{
				MMSeqUtils.generateHitsFile(gtfFile.getAbsolutePath(), samFile.getAbsolutePath(), hitsFile.getAbsolutePath(), readLength, 1000, 0);
				
				if(hitsFile.length() > 0)
				{
					System.out.println("Hits File contains data..." + hitsFile.length());
					break;
				}
				else
				{
					System.out.println("Hits File is empty... make another attempt " + (i+2) + " out of 10");
				}
			}
			
			// run mmseq
			System.out.println("Start MMseq: " + mmseqPath + " " + hitsFile.getAbsolutePath() + " " + hitsFile.getAbsolutePath());
			Utils.runProcess(mmseqPath + " " + hitsFile.getAbsolutePath() + " " + hitsFile.getAbsolutePath(), false, false, false);
			
			File isoformLevelFile = new File(hitsFile.getAbsolutePath() + ".mmseq");
			
			try
			{
				BufferedReader reader = new BufferedReader(new FileReader(isoformLevelFile));
				String line = null;
//				double expressionSum = 0;
				
				while((line = reader.readLine()) != null)
				{
					try
					{
						String[] split = line.split("\\s+");
						
						System.out.println("Isoform level prediction: " + line);
						
						isoformLevels.put(split[0], Math.exp(Double.parseDouble(split[2])));
						
//						expressionSum += Math.exp(Double.parseDouble(split[2]));
					}
					catch(Exception ex)
					{
						// do nothing, ignore line
					}
				}
				reader.close();
			}
			catch(Exception ex)
			{
				System.out.println("MMSeq isoform prediction was not successful.");
				ex.printStackTrace();
			}
			
			for(String transcript : gene.getArrayOfGeneProductNames())
			{
				if(!isoformLevels.containsKey(transcript))
					isoformLevels.put(transcript, 0.0);
			}
			
			Utils.removeFile(gtfFile.getAbsolutePath());
			Utils.removeFile(isoformLevelFile.getAbsolutePath());
			Utils.removeFile(hitsFile.getAbsolutePath());
			Utils.removeFile(hitsFile.getAbsolutePath() + ".amalgamated.mmseq");
			Utils.removeFile(hitsFile.getAbsolutePath() + ".gene.mmseq");
			
			return isoformLevels;
		}
		catch(Exception ex)
		{
			ex.printStackTrace();
			return null;
		}
	}

}
