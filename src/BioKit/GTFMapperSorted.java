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

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;
import java.util.zip.GZIPInputStream;

public class GTFMapperSorted 
{
	private BufferedReader reader;
	private String curLine;	
	private HashMap<Integer, Integer> geneToNumLines;
	private Vector<SlimExonicRegionGene> curGenes;
	private int maxGenomicPosCovered;
	private String currentChromosome;
	private int geneCounter;
	private boolean setCorrespondingExonsForExonicRegionForGenes;
	private boolean setCorrespondingTranscriptsForExonicRegionForGenes;
	/* Consecutive exonic regions overlap by one position.
	 * Therefore it can happen that there belong two
	 * exonic regions to one gene position.
	 */
	private HashMap<Integer, HashMap<Integer, Vector<Integer>>> geneToPosToExonicRegionID;
	private HashMap<Integer, HashMap<Integer, Integer>> geneToPosToStartPosNextExonicRegion;
	private HashMap<Integer, HashMap<Integer, GeneticRegion>> geneToExonicRegionIDsToExonicRegion;
	private HashMap<Integer, String> geneIDToGene;
	private HashMap<String, Integer> geneToGeneID;
	
	public GTFMapperSorted(String sortedGTFFile)
	{
		this.geneToNumLines = new HashMap<Integer, Integer>();
		init(sortedGTFFile);		
		this.curGenes = new Vector<SlimExonicRegionGene>();
		this.reader = getGTFFileReader(sortedGTFFile);
		this.curLine = null;
		this.maxGenomicPosCovered = -1;
		this.currentChromosome = null;
		this.geneCounter = 0;
		this.geneToPosToExonicRegionID = new HashMap<Integer, HashMap<Integer,Vector<Integer>>>();
		this.geneToPosToStartPosNextExonicRegion = new HashMap<Integer, HashMap<Integer,Integer>>();
		this.geneToExonicRegionIDsToExonicRegion = new HashMap<Integer, HashMap<Integer,GeneticRegion>>();
		this.geneIDToGene = new HashMap<Integer, String>();
		this.geneToGeneID = new HashMap<String, Integer>();
		this.setCorrespondingExonsForExonicRegionForGenes = true;
		this.setCorrespondingTranscriptsForExonicRegionForGenes = true;
	}
	
	private void init(String gtfFile)
	{
		try
		{
			System.out.println("Initialize GTFMapperSorted");
			int numGenesFound = 0;
			this.reader = getGTFFileReader(gtfFile);
			String curLine = "";
			String[] curInfo = null;
			
			String[] featuresSplit = null;
			String geneName = null;
			HashMap<String, Integer> geneToID = new HashMap<String, Integer>();
			while((curLine = this.reader.readLine()) != null)
			{
				curInfo = curLine.split("\\t");
								
				// only read exon lines
				if(!curInfo[2].equals("exon"))
					continue;
				
				featuresSplit = curInfo[8].split(";");
				
				for(String f : featuresSplit)
				{
					if(f.trim().startsWith("gene_id"))
					{
						geneName = f.replaceAll("gene_id", "").replaceAll("\"", "").trim();
						break;
					}					
				}
				if (!geneToID.containsKey(geneName))
				{
					numGenesFound++;
					geneToID.put(geneName, numGenesFound);
					this.geneToNumLines.put(numGenesFound, 0);
				}
				this.geneToNumLines.put(geneToID.get(geneName), this.geneToNumLines.get(geneToID.get(geneName)) + 1);
			}	
			this.reader.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.exit(0);
		}
	}
	
	private void prepareGene(SlimExonicRegionGene gene)
	{
//		this.curGenes.add(gene);
		GeneticRegion[] exonicRegions = gene.getExonicRegions();
		GeneticRegion curExonicRegion = null;
		GeneticRegion nextExonicRegion = null;
		HashMap<Integer, Vector<Integer>> curGenePosToExonicRegionID = new HashMap<Integer, Vector<Integer>>();
		HashMap<Integer, Integer> curGenePosToStartPosNextExonicRegion = new HashMap<Integer, Integer>();
		HashMap<Integer, GeneticRegion> exonicRegionIDToExonicRegion = new HashMap<Integer, GeneticRegion>();
		//DEBUG
		boolean debug = false;
		int interestingPosition = -1;
//		if (gene.getName().equals("ENSG00000237613_1")) debug = true;
		if (debug)
		{
			System.out.println("PREPARE: Adding new gene: " + gene.getGeneID() + " --> " + gene.getStart() + "-" + gene.getStop());
			System.out.println("PREPARE: Exonic regions:");
			for (GeneticRegion exonicRegion : exonicRegions)
			{
				System.out.println(exonicRegion.getCodingStart() + "-" + exonicRegion.getCodingStop());
			}			
		}
		//END DEBUG
		
		for (int i = 0; i < exonicRegions.length; i++)
		{
			curExonicRegion = exonicRegions[i];
			//DEBUG
			if (debug) System.out.println("PREPARE: Current exonic region: " 
					+ curExonicRegion.getCodingStart() + "-" + curExonicRegion.getCodingStop());
			//END DEBUG
			if ((i+1) < exonicRegions.length)
			{
				nextExonicRegion = exonicRegions[i+1];
				//DEBUG
				if (debug) System.out.println("PREPARE: Next exonic region: " 
						+ nextExonicRegion.getCodingStart() + "-" + nextExonicRegion.getCodingStop());
				//END DEBUG
			}
			else
			{
				//DEBUG
				if (debug) System.out.println("PREPARE: No more exonic regions set next exonic region to null");
				//END DEBUG
				/* null means that there exist no further exonic regions
				 * for the current gene.
				 */
				nextExonicRegion = null;
			}
			
			exonicRegionIDToExonicRegion.put(curExonicRegion.getIntegerID(), curExonicRegion);
			if (nextExonicRegion != null)
			{
				exonicRegionIDToExonicRegion.put(nextExonicRegion.getIntegerID(), nextExonicRegion);
			}
						
			int startCurExonicRegion = curExonicRegion.getGenomicStart();
			int stopCurExonicRegion = curExonicRegion.getGenomicStop();
			//DEBUG
			if (debug)
			{
				System.out.println("PREPARE: Setting current and next regions for pos of region " + startCurExonicRegion + "-" + stopCurExonicRegion);
			}
			//END DEBUG
			for (int pos = startCurExonicRegion; pos < (stopCurExonicRegion + 1); pos++)
			{
				//DEBUG
				if (debug && pos == interestingPosition)
				{
					System.out.println("PREPARE: Parameter settings:");
					System.out.println("PREPARE: Current region for " + pos + ": " + curExonicRegion.getGenomicStart() + "-" + curExonicRegion.getGenomicStop());
					if (nextExonicRegion != null)
						System.out.println("PREPARE: Next region for " + pos + ": " + nextExonicRegion.getGenomicStart() + "-" + nextExonicRegion.getGenomicStop());
					else System.out.println("PREPARE: Next region is null.");
				}
				
				if (debug && pos == interestingPosition) System.out.println("\nPREPARE: setting current region");
				//END DEBUG
				
				if (!curGenePosToExonicRegionID.containsKey(pos))
				{
					//DEBUG
					if (debug && pos == interestingPosition)
					{
						GeneticRegion tmp = exonicRegionIDToExonicRegion.get(curExonicRegion.getIntegerID());
						System.out.println("PREPARE: For " + pos + " set current region to " + curExonicRegion.getCodingStart() + "-" + curExonicRegion.getCodingStop() + "/" + tmp.getCodingStart() + "-" + tmp.getCodingStop());
					}
					//END DEBUG
					curGenePosToExonicRegionID.put(pos, new Vector<Integer>());
					curGenePosToExonicRegionID.get(pos).add(curExonicRegion.getIntegerID());
				}
				else
				{
					//DEBUG
					if (debug && pos == interestingPosition)
					{
						System.out.println("Current region is already set for position " + pos);
					}
					//END DEBUG
				}
				
				if ((nextExonicRegion != null) && pos == nextExonicRegion.getCodingStart())
				{
					//DEBUG
					if (debug && pos == interestingPosition)
					{
						System.out.println("PREPARE: " + pos + " is equals start position of next region.");
						GeneticRegion tmp = exonicRegionIDToExonicRegion.get(nextExonicRegion.getIntegerID());
						System.out.println("PREPARE: For " + pos + " add to current regions the region to " + nextExonicRegion.getCodingStart() + "-" + nextExonicRegion.getCodingStop() + "/" + tmp.getCodingStart() + "-" + tmp.getCodingStop());
					}
					//END DEBUG
					
					curGenePosToExonicRegionID.get(pos).add(nextExonicRegion.getIntegerID());
				}
				
				/* example of consecutive exonic regions:
				 * 1-5, 5-10, 10-15
				 * at position 5 the possible current
				 */
				//DEBUG
				if (debug && pos == interestingPosition) System.out.println("\n setting next region");
				//END DEBUG
				if (pos == startCurExonicRegion)
				{
					//DEBUG
					if (debug && pos == interestingPosition)
					{
						System.out.println("PREPARE: " + pos + " is equals start position.");
						GeneticRegion tmp = exonicRegionIDToExonicRegion.get(curExonicRegion.getIntegerID());
						System.out.println("PREPARE: For " + pos + " set next region to " + curExonicRegion.getCodingStart() + "-" + curExonicRegion.getCodingStop() + "/" + tmp.getCodingStart() + "-" + tmp.getCodingStop());
					}
					//END DEBUG
					curGenePosToStartPosNextExonicRegion.put(pos, curExonicRegion.getIntegerID());
				}
				else if (nextExonicRegion == null)
				{
					//DEBUG
					if (debug && pos == interestingPosition)
					{
						System.out.println("PREPARE: For " + pos + " set next region is null");
					}
					//END DEBUG
					curGenePosToStartPosNextExonicRegion.put(pos, null);
				}
				else 
				{
					//DEBUG
					if (debug && pos == interestingPosition)
					{
						System.out.println("PREPARE: " + pos + " is inbetween position.");
						GeneticRegion tmp = exonicRegionIDToExonicRegion.get(nextExonicRegion.getIntegerID());
						System.out.println("PREPARE: For " + pos + " set next region to " + nextExonicRegion.getCodingStart() + "-" + nextExonicRegion.getCodingStop() + "/" + tmp.getCodingStart() + "-" + tmp.getCodingStop());
					}
					//END DEBUG
					curGenePosToStartPosNextExonicRegion.put(pos, nextExonicRegion.getIntegerID());
				}
				//DEBUG
				if (debug && pos == interestingPosition)
				{
					System.out.println("PREPARE: DATA for pos " + pos + ":");
					System.out.print("PREPARE: Current region:");
					for (Integer exonicRegionId : curGenePosToExonicRegionID.get(pos))
					{
						System.out.print("\t" + exonicRegionIDToExonicRegion.get(exonicRegionId).getCodingStart() + "-" + exonicRegionIDToExonicRegion.get(exonicRegionId).getCodingStop());
					}
					System.out.println();
					System.out.println("PREPARE: Next region:\t" + exonicRegionIDToExonicRegion.get(curGenePosToStartPosNextExonicRegion.get(pos)).getCodingStart() + "-" + exonicRegionIDToExonicRegion.get(curGenePosToStartPosNextExonicRegion.get(pos)).getCodingStop());
				}
				//END DEBUG
			}
		}
		this.geneToPosToExonicRegionID.put(gene.getGeneID(), curGenePosToExonicRegionID);
		this.geneToPosToStartPosNextExonicRegion.put(gene.getGeneID(), curGenePosToStartPosNextExonicRegion);
		this.geneToExonicRegionIDsToExonicRegion.put(gene.getGeneID(), exonicRegionIDToExonicRegion);
	}
	
	private void readNextGenes(int regionStart, int regionStop, String chromosomeOfCurRegion)
	{
//		System.out.println("Read next gene(s)...");
		//DEBUG
//		System.out.println("Current region: " + regionStart + "-" + regionStop);
		//END DEBUG
		/* It is expected that not only the GTF file is sorted,
		 * but that also the range are asked in a sorted manner.
		 * Therefore genes lying before the region start are
		 * not needed anymore and can be removed.
		 */
		removeUnecessaryGenes(regionStart, chromosomeOfCurRegion);
		try
		{
			//DEBUG
			boolean debug = false;
//			if (regionStart == 879599 && regionStop == 879649 && chromosomeOfCurRegion.equals("chr1")) debug = true;
//			if (regionStart == 879800 && regionStop == 879849 && chromosomeOfCurRegion.equals("chr1")) debug = true;
//			if (regionStart == 93245 && regionStop == 93294 && chromosomeOfCurRegion.equals("chr10")) debug = true;
//			if (regionStart == 127926 && regionStop == 127975 && chromosomeOfCurRegion.equals("chr11")) debug = true;
			//END DEBUG
			String[] curInfo = null;
			
			String[] featuresSplit = null;
			String geneName = null;
			String transcriptName = null;
			
			String chromosome = null;
			Integer geneID = null;
			int curExonStart = -1;
			int curExonStop = -1;
			String exonID = null;
			boolean isPlusStrand = true;
			boolean startedToReadGenes = false;
			HashMap<Integer, Integer> geneIDToIndex = new HashMap<Integer, Integer>();
			HashMap<Integer, Integer> genesToNumLinesRead = new HashMap<Integer, Integer>();
			/* Since often exons are shared among transcripts they are mentioned several
			 * times gene, but for the different transcripts. In this case the transcript
			 * information is ignored and the exon is needed only once. To make sure each
			 * exon is only assigned once to its gene the exons already regarded are stored
			 * in this hase.
			 */
			HashMap<Integer, HashSet<String>> geneIDToExonsAlreadyRead = new HashMap<Integer, HashSet<String>>();
			
			HashMap<Integer, SlimExonicRegionGene> geneIDToGene = new HashMap<Integer, SlimExonicRegionGene>();
			boolean includeNextGene = false;
			/* The reading is continued as long as the end of the file is not reached
			 * and until either all genes that have been started to read have been
			 * finished or until the region which needs to be mapped is covered.
			 * Do not switch conditions in while loop otherwise some lines won't
			 * be regarded!!
			 */
			while((!startedToReadGenes || !genesToNumLinesRead.isEmpty() || this.maxGenomicPosCovered < regionStop 
					|| !this.currentChromosome.equals(chromosomeOfCurRegion) || includeNextGene) &&
					((this.curLine = this.reader.readLine()) != null))
			{
				//DEBUG
				if (includeNextGene && debug)
				{
					System.out.println("READ GENE: First line of next gene to include:");
					System.out.println("READ GENE: \t" + curLine);
				}
				//END DEBUG
				includeNextGene = false;
				curInfo = curLine.split("\\t");
							
				// only read exon lines
				if(!curInfo[2].equals("exon"))
					continue;
				
				curExonStart = Integer.parseInt(curInfo[3]);
				curExonStop = Integer.parseInt(curInfo[4]);
				exonID = curExonStart + "-" + curExonStop;
				featuresSplit = curInfo[8].split(";");
				
				for(String f : featuresSplit)
				{
					if(f.trim().startsWith("gene_id"))
					{
						geneName = f.replaceAll("gene_id", "").replaceAll("\"", "").trim();
					}					
					else if(f.trim().startsWith("transcript_id"))
					{
						transcriptName = f.replaceAll("transcript_id", "").replaceAll("\"", "").trim();
					}	
				}
				
				if (!this.geneToGeneID.containsKey(geneName))
				{
					//DEBUG
					if (debug) System.out.println("READ GENE: Started reading gene " + geneName);
					//END DEBUG
					
					chromosome = curInfo[0];
					if (curInfo[6].equals("-")) isPlusStrand = false;
					
					geneID = ++this.geneCounter;
					//DEBUG
					if (debug) System.out.println("READ GENE: geneName " + geneName + " - geneID " + geneID);
					//END DEBUG
					this.geneToGeneID.put(geneName, geneID);
					this.geneIDToGene.put(geneID, geneName);
					geneIDToIndex.put(geneID, this.curGenes.size());
					genesToNumLinesRead.put(geneID, 1);
					geneIDToExonsAlreadyRead.put(geneID, new HashSet<String>());
					geneIDToGene.put(geneID, 
							new SlimExonicRegionGene(geneName, geneID, chromosome, isPlusStrand));
					
					geneIDToGene.get(geneID).setSetCorrespondingExonsForExonicRegion(setCorrespondingExonsForExonicRegionForGenes);
					geneIDToGene.get(geneID).setSetCorrespondingTranscriptsForExonicRegion(setCorrespondingTranscriptsForExonicRegionForGenes);
					this.curGenes.add(geneIDToGene.get(geneID));
				}	
				else
				{
					geneID = this.geneToGeneID.get(geneName);
					genesToNumLinesRead.put(geneID, genesToNumLinesRead.get(geneID) + 1);
				}
								
				/* This exon was already encountered for this gene
				 * and therefore will be ignored.
				 */
				if (!geneIDToExonsAlreadyRead.get(geneID).contains(exonID))
				{
					Exon exon = new Exon(curExonStart, curExonStop, curExonStart, curExonStop);
					exon.setReference(chromosome);
					if (this.setCorrespondingTranscriptsForExonicRegionForGenes)
					{
						exon.addMappedTranscript(transcriptName);
					}
					geneIDToGene.get(geneID).addExon(exon);
					geneIDToExonsAlreadyRead.get(geneID).add(exonID);
					//DEBUG
					if (debug) 
					{
						System.out.println("\nREAD GENE: lines read = " + genesToNumLinesRead.get(geneID) + " out of " + this.geneToNumLines.get(geneID));
						System.out.println("READ GENE: Add exon " + exonID + " to gene " + geneName);
					}
					//END DEBUG
				}
				else
				{
					if (this.setCorrespondingTranscriptsForExonicRegionForGenes)
					{
						geneIDToGene.get(geneID).getExons().lastElement().addMappedTranscript(transcriptName);
						//DEBUG
						if (debug) 
						{
							System.out.println("DEBUG\tHave already seen exon " + exonID + " for gene " + geneName);
							System.out.println("DEBUG\tAdd transcript " + transcriptName + " to exon " + exonID);
							System.out.println("DEBUG\tThis exons needs to be the last one added to gene. Order of exons:");
							for (Exon exon : geneIDToGene.get(geneID).getExons())
								System.out.println(exon.getDescription());
							System.out.println("DEBUG\tTranscripts of last exon " + geneIDToGene.get(geneID).getExons().lastElement().getDescription() + ": ");
							for (String transcript : geneIDToGene.get(geneID).getExons().lastElement().getMappedTranscripts())
							{
								System.out.println(transcript);
							}
						}
						//END DEBUG
					}
				}
				
				if (genesToNumLinesRead.get(geneID).equals(this.geneToNumLines.get(geneID)))
				{
					//DEBUG
					if (debug) System.out.println("READ GENE: Done reading gene " + geneName);
					//END DEBUG
					
					if (geneIDToGene.get(geneID).getStop() < regionStart || 
							!chromosomeOfCurRegion.equals(geneIDToGene.get(geneID).getChromosome()))
					{
						//DEBUG
						if (debug) System.out.println("READ GENE: Ignore gene " + geneName + " (" + geneIDToGene.get(geneID).getStart() + "-" + geneIDToGene.get(geneID).getStop() + ") since it does not overlap with region " + regionStart + "-" + regionStop + "!");
						//END DEBUG
						genesToNumLinesRead.remove(geneID);
						this.curGenes.remove(geneIDToGene.get(geneID));
						geneIDToGene.remove(geneID);
						geneIDToIndex.remove(geneID);
						geneIDToExonsAlreadyRead.remove(geneID);
						continue;
					}
					
					//DEBUG
					if (debug) System.out.println("READ GENE: Gene " + geneName  + " (" + geneIDToGene.get(geneID).getStart() + "-" + geneIDToGene.get(geneID).getStop() + ") overlaps with region " + regionStart + "-" + regionStop + "!");
					//END DEBUG
					
					prepareGene(geneIDToGene.get(geneID));
					if (this.maxGenomicPosCovered < geneIDToGene.get(geneID).getStop())
					{
						this.maxGenomicPosCovered = geneIDToGene.get(geneID).getStop();
					}
					/* In a sorted GTF file there exists a switch between the chromosomes
					 * which means in term of positions the counting starts from the
					 * beginning. Therefore the max genomic position covered has to be 
					 * newly set when chromosomes change. 
					 */
					else if (!this.currentChromosome.equals(geneIDToGene.get(geneID).getChromosome()))
					{
						this.maxGenomicPosCovered = geneIDToGene.get(geneID).getStop();
					}
					this.currentChromosome = geneIDToGene.get(geneID).getChromosome();
					genesToNumLinesRead.remove(geneID);
					geneIDToGene.remove(geneID);
					geneIDToIndex.remove(geneID);
					geneIDToExonsAlreadyRead.remove(geneID);
				}
				startedToReadGenes = true;
				
				/* If the subsequent conditions are fullfilled, all genes that were started to be read
				 * are finished and the current region is covered. However it can occur that the next gene
				 * also overlaps with the current region, but the first exon starts shortly after the last
				 * exon of the last gene read. This gene also needs to be read for appropriate mapping of
				 * the current region to the transcriptome. 
				 * Example:
				 * region: 
				 * Current Mapping
				 * ERR030880.4230505       chr1    879681-879730
				 * 50M
				 * 
				 * last exon first gene:    chr1    exon    879288  879955  gene_id "ENSG00000187634_1";
				 * first exon second gene:  chr1    exon    879584  880180  gene_id "ENSG00000188976_1";
				 * 
				 * The first gene already completely covers the specified region. The genes do not overlap 
				 * with each other, but both overlap with region. Therefore, when all genes started to be 
				 * read are read and the region is covered the next gene needs to tested if still overlaps 
				 * with the region. If not one has to go back to its start and read next time needed.
				 */
				if (startedToReadGenes && genesToNumLinesRead.isEmpty() && this.maxGenomicPosCovered >= regionStop)
				{
					// only read exon lines					
					this.reader.mark(1000);
					this.curLine = this.reader.readLine();
					curInfo = curLine.split("\\t");
					
					while (!curInfo[2].equals("exon"))
					{
						this.curLine = this.reader.readLine();
						curInfo = curLine.split("\\t");
					}
				
					chromosome = curInfo[0];
					curExonStart = Integer.parseInt(curInfo[3]);
					curExonStop = Integer.parseInt(curInfo[4]);
					
					// Check whether next gene overlaps with 
					// region read so far from GTF file
					if (curExonStart < this.maxGenomicPosCovered && chromosome.equals(chromosomeOfCurRegion))
					{
						//DEBUG
						if (debug)
						{
							System.out.println("READ GENE: Current line: " + this.curLine);
							System.out.println("READ GENE: Exon start probed gene: " + curExonStart);
							System.out.println("READ GENE: current max genomic pos covered" + this.maxGenomicPosCovered);
							System.out.println("READ GENE: INCLUDE PROBED GENE");
						}
						//END DEBUG
						includeNextGene = true;
					}
					else
					{
						//DEBUG
						if (debug)
						{
							System.out.println("READ GENE: Current line: " + this.curLine);
							System.out.println("READ GENE: Exon start probed gene: " + curExonStart);
							System.out.println("READ GENE: current max genomic pos covered" + this.maxGenomicPosCovered);
							System.out.println("READ GENE: DO NOT INCLUDE PROBED GENE!!!");
						}
						//END DEBUG
					}
					reader.reset();
				}
			}
			// Close stream when file is read completely.
			// In this case and only this case curLine is null.
			if (curLine == null) this.reader.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.exit(0);
		}
	}
	
	public HashMap<String, Vector<MappableElement>> mapPositionToGeneticElements(String regionChromosome, 
			int regionStart, int regionStop, boolean regionIsPlusStrand, boolean regardStrand)
	{
		//DEBUG
		boolean debug = false; 
//		if (regionStart == 879599 && regionStop == 879648 && regionChromosome.equals("chr1")) debug = true;
//		if (regionStart == 879800 && regionStop == 879849 && regionChromosome.equals("chr10")) debug = true;
//		if (regionStart == 93245 && regionStop == 93294 && regionChromosome.equals("chr10")) debug = true;
//		if (regionStart == 127926 && regionStop == 127975 && regionChromosome.equals("chr11")) debug = true;
		//END DEBUG
		if (debug) System.out.println("DEBUG\tSearching overlap for region: " + regionStart + "-" + regionStop);
		HashMap<String, Vector<MappableElement>> geneToMappableElements = 
				new HashMap<String, Vector<MappableElement>>(); 
		
		/* Make sure the region needed to be mapped is covered by
		 * the current genes.
		 */
		//DEBUG
		if (debug)
		{
			System.out.println("DEBUG: Max genomic position covered: " + this.maxGenomicPosCovered);
			System.out.println("DEBUG: Current chromosome: " + this.currentChromosome);
			System.out.println("DEBUG: Current genes");
			for (SlimExonicRegionGene gene : this.curGenes)
			{
				System.out.println("DEBUG: " + gene.getName() + ":" + gene.getChromosome() + ":" + gene.getStart() + "-" + gene.getStop());
			}
		}
		if (this.maxGenomicPosCovered < regionStop && debug)
		{
			System.out.println("DEBUG: Max genomic position covered: " + this.maxGenomicPosCovered);
			System.out.println("DEBUG: Next gene(s) has to be read!");
		}
		if ((this.maxGenomicPosCovered < regionStop || !this.currentChromosome.equals(regionChromosome))
				&& debug) 
		{
			System.out.println("DEBUG: Max genomic position covered: " + this.maxGenomicPosCovered);
			System.out.println("DEBUG: Current chromosome: " + this.currentChromosome);
			System.out.println("Next gene has to be read!");
		}
		//END DEBUG
		if (this.maxGenomicPosCovered < regionStop || !this.currentChromosome.equals(regionChromosome)) readNextGenes(regionStart, regionStop, regionChromosome);
		//DEBUG
		
//		System.out.println("Genes to search: " + this.curGenes.size());
//		System.out.println("Searching overlap for region: " + regionStart + "-" + regionStop);

		//END DEBUG
		
		for (int i = 0; i < this.curGenes.size(); i++)			
		{
			SlimExonicRegionGene gene = this.curGenes.get(i);

			//DEBUG
			if (debug) System.out.println("DEBUG: Range current gene " + gene.getName() + ": " + gene.getStart() + "-" + gene.getStop());
			//END DEBUG
			
			if (!gene.getChromosome().equals(regionChromosome))
			{
				//DEBUG
				if (debug) System.out.println("DEBUG: Gene ignored. Gene has different reference than region: gene -> " 
						+ gene.getChromosome() + ", region -> " + regionChromosome);
				//END DEBUG
				continue;
			}
			if (regardStrand && ((gene.isPlusStrand() && !regionIsPlusStrand)
					|| (!gene.isPlusStrand() && regionIsPlusStrand)))
			{
				//DEBUG
				if (debug) System.out.println("DEBUG: Gene ignored. Gene has different strand than region: gene is + strand -> " 
						+ gene.isPlusStrand() + ", region is + strand -> " + regionIsPlusStrand);
				//END DEBUG
				continue;
			}
			
			/* It can happen that the last genes read lied before
			 * the region. Then the gene (the current gene was read), but
			 * this lies after the mapped region. Thus, the region does not
			 * map to any gene
			 */
			if (regionStop < gene.getStart())
			{
				//DEBUG
				if (debug) System.out.println("DEBUG: Gene and genes afterwards ignored since region stops before gene start: gene -> " 
						+ gene.getStart() + ", region -> " + regionStart);
				//END DEBUG
				break;
			}
			
			//DEBUG
			if (debug) System.out.println("DEBUG: Region to be mapped: " + regionStart + "-" + regionStop);
			//END DEBUG
			
			HashMap<Integer, Vector<Integer>> posToExonicRegionID = this.geneToPosToExonicRegionID.get(gene.getGeneID());
			HashMap<Integer, Integer> posToNextExonicRegionID = this.geneToPosToStartPosNextExonicRegion.get(gene.getGeneID());
			HashMap<Integer, GeneticRegion> exonicRegionIDsToexonicRegions  = this.geneToExonicRegionIDsToExonicRegion.get(gene.getGeneID());
			int curStart = -1;
			int curStop = -1;
			GeneticRegion curGeneticRegion = null;
			boolean ignoreGene = false;
		
			Vector<Integer> mappedExonicRegions = new Vector<Integer>();
			Vector<Integer> allMappedExonicRegions = new Vector<Integer>();
			/* If there exist no exonic regions for the start position
			 * then there exists no exon in this gene the read could
			 * have originated from and it is ignored.
			 */
			if (!posToExonicRegionID.containsKey(regionStart))
			{
				//DEBUG
				if (debug)
				{					
					System.out.println("DEBUG " + gene.getName() + ": Gene is ignored since it does not fully overlap with region!");
					System.out.println("DEBUG " + gene.getName() + ": " + regionStart + " has no exonic region ID information");
				}
				//END DEBUG
				continue;
			}
			
			
			//DEBUG
			if (debug)
			{					
				System.out.println("DEBUG " + gene.getName() + ": Gene overlaps with region!");
				System.out.println("DEBUG " + gene.getName() + ": " + regionStart + " has exonic region ID information");
			}
			//END DEBUG
			
			mappedExonicRegions.addAll(posToExonicRegionID.get(regionStart));
			
			/* if regionStart == regionStop holds the exonic regions before or 
			 * after a junction are seeked. In this case both exonic regions
			 * found at that position are returned since it not known whether
			 * the junction start or the junction stop is given.
			 */
			if (regionStart == regionStop)
			{
				for (Integer exonicRegionID : mappedExonicRegions)
				{
					curGeneticRegion = exonicRegionIDsToexonicRegions.get(exonicRegionID);
					allMappedExonicRegions.add(curGeneticRegion.getElementID());
					addMappableElementToGene(geneToMappableElements, exonicRegionIDsToexonicRegions, allMappedExonicRegions, 
							gene, curGeneticRegion.getCodingStart(), 
							curGeneticRegion.getCodingStop());
				}
				continue;
			}
			
			curGeneticRegion = exonicRegionIDsToexonicRegions.get(mappedExonicRegions.firstElement());
			
			/* Consecutive exonic regions overlap by one position. In case that the start of the
			 * requested region exactly lies on the border of two consecutive exonic regions
			 * the exonic region lying before the region start is ignored. The aim is to find
			 * the exonic regions that cover the requested region as precisely as possible.
			 */
			if (curGeneticRegion.getCodingStop() == regionStart && mappedExonicRegions.size() > 1)
			{
				mappedExonicRegions.remove(0);
				curGeneticRegion = exonicRegionIDsToexonicRegions.get(mappedExonicRegions.firstElement());
			}
														
			while(curGeneticRegion != null)
			{	
				//DEBUG
				if (debug) 
				{
					System.out.println("DEBUG: Regions in vector");
					for (Integer regionID : mappedExonicRegions)
					{
						System.out.println("DEBUG: \t" + exonicRegionIDsToexonicRegions.get(regionID).getCodingStart() + "-" + exonicRegionIDsToexonicRegions.get(regionID).getCodingStop());
					}
					
					System.out.println("DEBUG: Checking overlap with exonic region: " + 
							curGeneticRegion.getCodingStart() + "-" + curGeneticRegion.getCodingStop());
				}
				//END DEBUG
				if (curStart == -1)
				{
					curStart = curGeneticRegion.getCodingStart();
					curStop = curGeneticRegion.getCodingStop(); 
					allMappedExonicRegions.add(curGeneticRegion.getElementID());
					//DEBUG
					if (debug) 
					{
						System.out.println("DEBUG\tAdded exonic region to all mapped exonic regions");
						System.out.println("DEBUG\tCurrent:" + curGeneticRegion.getElementID());
						System.out.println("DEBUG\tTotal:" + allMappedExonicRegions);
					}
					//END DEBUG
				}
				else
				{
					/* Check whether the current exonic region and the one stored are
					 * consecutive, namely overlap by one position. If this is the case
					 * they belong to the same exon and are summarized as one exonic region.
					 */
					if (curStop == curGeneticRegion.getCodingStart())
					{
						curStop = curGeneticRegion.getCodingStop();
						allMappedExonicRegions.add(curGeneticRegion.getElementID());
						//DEBUG
						if (debug) 
						{
							System.out.println("DEBUG\tAdded exonic region to all mapped exonic regions");
							System.out.println("DEBUG\tCurrent:" + curGeneticRegion.getElementID());
							System.out.println("DEBUG\tTotal:" + allMappedExonicRegions);
						}
						//END DEBUG
					}
					/* If the current exonic region does not overlap with the previous then
					 * they do not lie on the same exon. Hence the gene has a junction within
					 * the specified region. However, since the specified region is expected 
					 * to originate from one exon the current gene is ignored.
					 */
					else
					{
						//DEBUG
						if (debug)
						{
							System.out.println("DEBUG: Gene " + gene.getName() + " is ignored due to junction within region " + regionStart + "-" + regionStop);
							System.out.println("DEBUG: Gene junction: " + curStop + ":" + curGeneticRegion.getCodingStart());
						}
						//END DEBUG
						ignoreGene = true;
						break;
					}
				}
				
				/* Check whether subsequent exonic regions are already considered.
				 * If not continue with subsequent exonic region if there exist 
				 * further exonic regions for this genes and in they still lie
				 * in the specified region. If the start position of the next
				 * exon is null the current exonic region is the last.
				 */
				mappedExonicRegions.remove(0);
				if (mappedExonicRegions.isEmpty() && curStop < gene.getStop() && curStop < regionStop)
				{
					//DEBUG
//					System.out.printf("No more exonic regions to process looking for next (curStop=%d, geneStop=%d, regionStop=%d)\n",
//							curStop, gene.getStop(), regionStop);
					
					//END DEBUG
					if (posToNextExonicRegionID.get(curStop) != null)
					{
						mappedExonicRegions.add(posToNextExonicRegionID.get(curStop));
						//DEBUG
//						GeneticRegion tmp = exonicRegionIDsToExonicRegion.get(posToNextExonicRegionID.get(curStop));
//						System.out.println("next exonic region for position " + curStop + ": " + tmp.getCodingStart() + "-" + tmp.getCodingStop());
						//END DEBUG
					}
					//DEBUG
//					System.out.println("Found " + mappedExonicRegions.size() + " exonic regions for curStop pos: " + curStop);
					//END DEBUG
				}
				
				if (!mappedExonicRegions.isEmpty())
				{
					curGeneticRegion = exonicRegionIDsToexonicRegions.get(mappedExonicRegions.firstElement());
					//DEBUG
//					System.out.println("next current exonic region for position " + curStop + ": " + curExonicRegion.getCodingStart() + "-" + curExonicRegion.getCodingStop());
					//END DEBUG
				}
				else curGeneticRegion = null;
			}
			
			if (ignoreGene) continue;
			
			/* Only if the current region is clearly covered by consecutive 
			 * exonic regions then this exonic region and gene are collected.
			 */
			if (curStart <= regionStart && curStop >= regionStop)
			{
				addMappableElementToGene(geneToMappableElements, exonicRegionIDsToexonicRegions, allMappedExonicRegions, gene, curStart, curStop);
			}
		}
		//DEBUG
		if (geneToMappableElements.isEmpty())
		{
//			System.out.println(debugString);
//			System.out.println("DEBUG Unmapped region: " + regionChromosome + ":" + regionStart + "-" + regionStop);
//			System.out.println("DEBUG Unmapped region Current genes: ");
//			for (SlimExonicRegionGene gene : this.curGenes)
//			{
//				System.out.println("DEBUG Unmapped region \t" + gene.getName() + ":" + gene.getStart() + "-" + gene.getStop());
//				for (GeneticRegion exonicRegion : gene.getExonicRegions())
//				{
//					System.out.println("DEBUG Unmapped region \t\t" + exonicRegion.getCodingStart() + "-" + exonicRegion.getCodingStop());					
//				}
//			}
//			System.out.println("DEBUG Unmapped region");
			
//			System.exit(0);
		}
		if (debug && geneToMappableElements.isEmpty())
		{
//			System.out.println(debugString);
//			System.out.println("DEBUG Unmapped region: " + regionChromosome + ":" + regionStart + "-" + regionStop);
//			System.out.println("DEBUG Unmapped region Current genes: ");
//			for (SlimExonicRegionGene gene : this.curGenes)
//			{
//				System.out.println("DEBUG Unmapped region \t" + gene.getName() + ":" + gene.getStart() + "-" + gene.getStop());
//				for (GeneticRegion exonicRegion : gene.getExonicRegions())
//				{
//					System.out.println("DEBUG Unmapped region \t\t" + exonicRegion.getCodingStart() + "-" + exonicRegion.getCodingStop());					
//				}
//			}
//			System.out.println("DEBUG Unmapped region");
//			System.exit(0);
		}
 		//END DEBUG
		return geneToMappableElements;
	}
		
	private BufferedReader getGTFFileReader(String pathToGTFFile)
	{
		BufferedReader reader = null;
		try
		{
			if (pathToGTFFile.endsWith(".gz")) 
			{
				GZIPInputStream gzip = new GZIPInputStream(new FileInputStream(pathToGTFFile));
				reader = new BufferedReader(new InputStreamReader(gzip));
				gzip.close();
			}
			else reader = new BufferedReader(new FileReader(pathToGTFFile));
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.exit(0);
		}
		return reader;
	}
	
	private void removeUnecessaryGenes(int regionStart, String chromosome)
	{
		int numOfFirstGenesToBeRemoved = 0;
		for (int i = 0; i < this.curGenes.size(); i++)
		{
			if (this.curGenes.get(i).getStop() < regionStart 
					|| !chromosome.equals(this.curGenes.get(i).getChromosome()))
			{
				numOfFirstGenesToBeRemoved++;
			}
			else break;
		}
		
		removeNumFirstGenes(numOfFirstGenesToBeRemoved);
	}
	
	private void removeNumFirstGenes(int numOfGenesToRemove)
	{
		for (int i = 0; i < numOfGenesToRemove; i++)
		{
			//DEBUG
//			System.out.println("Removing gene: " + this.curGenes.get(0).getName());
			//END DEBUG
			int geneID = curGenes.get(0).getGeneID();
			this.geneToNumLines.remove(geneID);
			this.curGenes.remove(0);
			this.geneToExonicRegionIDsToExonicRegion.remove(geneID);
			this.geneToPosToExonicRegionID.remove(geneID);
			this.geneToPosToStartPosNextExonicRegion.remove(geneID);
			this.geneToGeneID.remove(geneIDToGene.get(geneID));
			this.geneIDToGene.remove(geneID);
			//DEBUG
//			System.out.printf("Removed gene. Num of genes: \n\tcurGenes=%d, \n\tgeneToExonicRegionIDsToExonicRegion=%d, " +
//					"\n\tgeneToPosToExonicRegionID=%d, \n\tgeneToPosToStartPosNextExonicRegion=%d, \n\tgeneToNumLines=%d\n", this.curGenes.size(),
//					this.geneToExonicRegionIDsToExonicRegion.size(), this.geneToPosToExonicRegionID.size(),
//					this.geneToPosToStartPosNextExonicRegion.size(), this.geneToNumLines.keySet().size());
			//END DEBUG
		}
	}	
	
	/* Adds exonic region covered. Consecutive exonic regions are merged!
	 * It is assumed that the region specified by curStart and curStop in 
	 * the parameters is continous. Example: 
	 * exonic region 1: 5-10
	 * exonic region 2: 10-30
	 * --> 5-30
	 */	
	private void addMappableElementToGene(HashMap<String, Vector<MappableElement>> geneToMappableFeatures,
			HashMap<Integer, GeneticRegion> exonicRegionIDsToExonicRegions, 
			Vector<Integer> allMappedExonicRegions, SlimExonicRegionGene gene, int curStart, int curStop)
	{
		GeneticRegion exonicRegion = new GeneticRegion(curStart, curStop, curStart, curStop);
		exonicRegion.setReference(gene.getChromosome());
		
		
		//DEBUG
		if (allMappedExonicRegions.isEmpty())
		{
			System.out.println("DEBUG\tallMappedExonicRegions is empty!!");
			System.out.println("DEBUG\tCurrent exonic region: " + exonicRegion.getDescription());
			System.exit(0);
		}
		//END DEBUG
		
		if (this.setCorrespondingExonsForExonicRegionForGenes 
				|| this.setCorrespondingTranscriptsForExonicRegionForGenes)
		{
		  HashSet<Exon> exons = new HashSet<Exon>();
		  for (Integer exonicRegionID : allMappedExonicRegions)
		  {
			  exons.addAll(exonicRegionIDsToExonicRegions.get(exonicRegionID).getMappedExons());
		  }
		  
		  /* Iterate over all exons of all mapped exonic regions and check which exons
		   * fully overlap the newly created exonic region. Annotate these exons for the
		   * exonic regions and also the corresponding transcripts.
		   */
		  for (Exon exon : exons)
		  {
			  if (!(exon.getCodingStart() <= exonicRegion.getCodingStart() && exon.getCodingStop() >= exonicRegion.getCodingStop()))
			  {
				  continue;
			  }
			  exonicRegion.addMappedExon(exon);
			  if (this.setCorrespondingTranscriptsForExonicRegionForGenes)
			  {
				  //DEBUG
//					  System.out.println("DEBUG\tAdding transcripts of exon" + exon.getDescription());
				  //END DEBUG
				  for (String transcript : exon.getMappedTranscripts()) 
				  {
					  //DEBUG
//						  System.out.println("DEBUG\t\t" + transcript);
					  //END DEBUG
					  exonicRegion.addMappedTranscript(transcript);
				  }
			  }
		  }
		}		
		
		if (!geneToMappableFeatures.containsKey(gene.getName()))
		{
			geneToMappableFeatures.put(gene.getName(), new Vector<MappableElement>());
		}
		geneToMappableFeatures.get(gene.getName()).add(exonicRegion);
		//DEBUG
//		System.out.println("Added exonic region " + curStart + "-" + curStop + " to current gene " + geneName);
		//END DEBUG
		
	}
	
	public Vector<SlimExonicRegionGene> getCurGenes()
	{
		return this.curGenes;
	}

	public void setSetCorrespondingExonsForExonicRegionForGenes(
			boolean setCorrespondingExonsForExonicRegionForGenes) 
	{
		this.setCorrespondingExonsForExonicRegionForGenes = setCorrespondingExonsForExonicRegionForGenes;
	}

	public void setSetCorrespondingTranscriptsForExonicRegionForGenes(
			boolean setCorrespondingTranscriptsForExonicRegionForGenes) 
	{
		this.setCorrespondingTranscriptsForExonicRegionForGenes = setCorrespondingTranscriptsForExonicRegionForGenes;
	}	
	
	public Vector<GeneticRegion> getExonicRegionsForGene(String gene)
	{
		Vector<GeneticRegion> exonicRegionsForGene = new Vector<GeneticRegion>();
		if (geneToExonicRegionIDsToExonicRegion.containsKey(this.geneToGeneID.get(gene)))
		{
			exonicRegionsForGene.addAll(geneToExonicRegionIDsToExonicRegion.get(this.geneToGeneID.get(gene)).values());
		}
		return exonicRegionsForGene;
	}
}
