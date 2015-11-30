package BioKit;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Vector;

/**
 * 
 * Wrapper for GTF gene information. Memory efficient version to store gene information. 
 * 
 * Contains inbuilt method to convert to standard Biokit Gene (less memory efficient).
 * 
 * @author Fabian Birzele
 *
 */
public class GTFGene 
{
	private String id;
	private boolean isPlusStrand;
	private String chromosome;
	private String type;
	private String synonyms;
	private String geneName;
	private String geneType;	// review http://www.gencodegenes.org/gencode_biotypes.html for gene types
	
	private HashMap<String, Vector<Interval>> exonIntervals;
	
	public GTFGene(String chromosome, boolean isPlusStrand, String id, String geneName, String geneType, String type)
	{
		this.chromosome = chromosome;
		this.isPlusStrand = isPlusStrand;
		this.id = id;
		this.type = type;
		this.geneName = geneName;
		this.geneType = geneType;
		
		exonIntervals = new HashMap<String, Vector<Interval>>();
	}
	
	public void setSynonymInformation(String synonyms)
	{
		this.synonyms = synonyms;
	}
	
	public void addExon(String transcript, int start, int stop)
	{
		if(!exonIntervals.containsKey(transcript))
		{
			exonIntervals.put(transcript, new Vector<Interval>());
		}
		
		exonIntervals.get(transcript).add(new Interval(start, stop));
	}
	
	public String getId() 
	{
		return id;
	}

	public boolean isPlusStrand() 
	{
		return isPlusStrand;
	}

	public String getChromosome() 
	{
		return chromosome;
	}

	public String getType()
	{
		return type;
	}
	
	public String getGeneName()
	{
		return geneName;
	}
	
	public String getGeneType()
	{
		return geneType;
	}

	public HashMap<String, Vector<Interval>> getExonIntervals() 
	{
		return exonIntervals;
	}
	
	public String getSynonyms()
	{
		return synonyms;
	}
	
	public Vector<Exon> getExons()
	{
		HashMap<String, Exon> allExons = new HashMap<String, Exon>();
		Vector<Exon> exons = new Vector<Exon>();
		
		for(String t : exonIntervals.keySet())
		{
			Vector<Interval> intervals = exonIntervals.get(t);
			
			for(Interval i : intervals)
			{
				if(allExons.containsKey(i.getStart() + "-" + i.getStop()))
				{
					continue;
				}
				else
				{
					Exon exon = new Exon(i.getStart(), i.getStop());
					exon.setID((getChromosome() + "-" + i.getStart() + "-" + i.getStop()).hashCode());
					allExons.put(i.getStart() + "-" + i.getStop(), exon);
					exons.add(exon);
				}
			}
		}
		
		return exons;
	}
	
	public Gene createGene()
	{
		Gene gene = new Gene(id, geneName, geneType, Integer.MAX_VALUE, Integer.MIN_VALUE, chromosome, isPlusStrand);
		
		HashMap<String, Exon> allExons = new HashMap<String, Exon>();
		
		for(String t : exonIntervals.keySet())
		{
			Vector<Interval> intervals = exonIntervals.get(t);
			
			for(Interval i : intervals)
			{
				if(allExons.containsKey(i.getStart() + "-" + i.getStop()))
				{
					continue;
				}
				else
				{
					Exon exon = new Exon(i.getStart(), i.getStop());
					exon.setID((getChromosome() + "-" + i.getStart() + "-" + i.getStop()).hashCode());
					allExons.put(i.getStart() + "-" + i.getStop(), exon);
				}
			}
		}
		
		for(String t : exonIntervals.keySet())
		{
			Vector<Interval> intervals = exonIntervals.get(t);
			
			Exon[] exons = new Exon[intervals.size()];
			int counter = 0;
			
			for(Interval i : intervals)
			{
				Exon exon = allExons.get(i.getStart() + "-" + i.getStop());
				gene.addExon(exon);
				exons[counter] = exon;
				counter++;
			}
			
			gene.addGeneProduct(t, "", exons);
		}
		
		//System.out.println(gene.toString());
		
		return gene;
	}
	
	/**
	 * This method creates a gene from a GTF lightweight wrapper instance. With a specified genome file and samtools path, the method
	 * will also try to parse the exon sequence information for the gene.
	 * 
	 * @param genomeFile
	 * @param samtoolsPath
	 * @return
	 * @throws Exception
	 */
	public Gene createGeneIncludingSequenceInformation(String genomeFile, String samtoolsPath) throws Exception
	{
		Gene gene = createGene();
		
		File sequenceFile = new File("temp" + (int)(Math.random()*100000) + ".fasta");
		
		Utils.runProcessInShell(samtoolsPath +" faidx "+genomeFile+" " +gene.getChromosome() + ":"+ +gene.getStart() + "-" + +gene.getStop() + " > " + sequenceFile.getAbsolutePath(), false, false, null);	
		System.out.println(samtoolsPath +" faidx "+genomeFile+" " +gene.getChromosome() + ":"+ +gene.getStart() + "-" + +gene.getStop() + " > " + sequenceFile.getAbsolutePath());
		Thread.sleep(100);
				
		BufferedReader reader = new BufferedReader(new FileReader(sequenceFile.getAbsolutePath()));
		String line = null;
				
		StringBuffer sequence = new StringBuffer();
				
		while((line = reader.readLine()) != null)
		{
			if(line.startsWith(">"))
				continue;
					
			sequence.append(line.trim());
		}
		
		reader.close();
				
		int geneLength = gene.getStop()-gene.getStart()+1;
				
		if(geneLength != sequence.length())
		{
			Utils.runProcessInShell(samtoolsPath +" faidx "+genomeFile+" " +gene.getChromosome() + ":"+ gene.getStart() + "-" + gene.getStop() + " > " + sequenceFile.getAbsolutePath(), false, false, null);
					
			Thread.sleep(5000);
					
			reader = new BufferedReader(new FileReader(sequenceFile.getAbsolutePath()));
			line = null;
				
			sequence = new StringBuffer();
				
			while((line = reader.readLine()) != null)
			{
				if(line.startsWith(">"))
					continue;
					
				sequence.append(line.trim());
			}
			
			reader.close();
				
			geneLength = gene.getStop()-gene.getStart()+1;
				
			if(geneLength != sequence.length())
			{
				throw new Exception("Can not extract sequence for gene " + gene.getGeneID() + " " + chromosome + " " + gene.getStart() + " " + gene.getStop() + " Expected length: " + geneLength + ", Got length: " + sequence.length());
			}
		}
		
		deleteFile(sequenceFile.getAbsolutePath());
				
		Exon[] exons = gene.getSortedExons();
				
		int offset = gene.getStart();
				
		for(Exon e : exons)
		{
			String exonSequence = sequence.substring(e.getGenomicStart()-offset, e.getGenomicStop()-offset+1);
					
			if(!gene.isPlusStrand())
				exonSequence = (new DNA("exon", exonSequence)).getReverseCompDNASequence().getSequence();
					
			if(e.getGenomicLength() == exonSequence.length())
			{
				e.setGenomicNucleotideSequence(exonSequence);
			}
			else
			{
				System.err.println("\tExtraction of exon sequence failed..." + gene.getGeneID() + " " + e.getGenomicStart() + "-" 
					+ e.getGenomicStop() + " " + e.getGenomicLength() + " vs. " + exonSequence.length() + ", exon sequence is not available...");
				System.err.println("Call for samtools was: " +samtoolsPath +" faidx "+genomeFile+" " +gene.getChromosome() + ":"+ gene.getStart() + "-" + gene.getStop() + " > " + sequenceFile.getAbsolutePath());
			}
		}
		
		return gene;
	}
	
	public Gene createGeneIncludingSequenceInformation(RandomAccessFastaReader reader) throws Exception
	{
		Gene gene = createGene();
		
		HashMap<String, String> exonSequences = new HashMap<String, String>();
		
		//System.out.println("Process gene: " + gene.getGeneID());
		
		for(String value : reader.getSequencesForIdHashcode(gene.getGeneID().hashCode()))
		{
			String[] split = value.split("\t");
			exonSequences.put(split[0], split[1]);
			
			//System.out.println("Added sequence: " + value);
		}
				
		Exon[] exons = gene.getSortedExons();
				
		for(Exon e : exons)
		{
			String sequence = exonSequences.get(e.getGenomicStart() + "-" + e.getGenomicStop());
			
			if(sequence != null && e.getGenomicLength() == sequence.length())
			{
				e.setGenomicNucleotideSequence(sequence);
			}
			else
			{
				System.err.println("\tExtraction of exon sequence failed..." + gene.getGeneID() + " " + e.getGenomicStart() + "-" 
					+ e.getGenomicStop() + ", exon sequence is not available...");
			}
		}
		
		return gene;
	}
	
	private void deleteFile(String file)
	{
		File f = new File(file);
		
		while(f.exists())
		{
			f.delete();
		}
		
		//System.out.println("File "+file+" has been deleted...");
	}
	
	public String getStrandStringId()
	{
		if(isPlusStrand)
			return "+";
		else
			return "-";
	}
	
	public String toString()
	{
		StringBuffer b = new StringBuffer();
		b.append(id + "\t" + chromosome + "\t" + getStrandStringId() + "\t" + type);
		
		for(String t : exonIntervals.keySet())
		{
			b.append("\t" + t + ":" + exonIntervals.get(t).size());
		}
		
		return b.toString();
	}

	private class Interval
	{
		private int start;
		private int stop;
		
		public Interval(int start, int stop)
		{
			this.start = start;
			this.stop = stop;
		}

		public int getStart() 
		{
			return start;
		}

		public int getStop() 
		{
			return stop;
		}
	}

}
