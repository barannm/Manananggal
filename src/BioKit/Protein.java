package BioKit;

import java.io.File;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

/**
 * A protein instance stores all information associated with a protein in the biokit toolkit. As a standard, the instance
 * holds an array containing the AminoAcids of the protein and the sequence as well as proteinID and pdbFileName. 
 * Additionally one can add properties to the protein identified by a name. The property can be any object, so that
 * it could be an multiple alignment, an additional 1D property like Solvent Accessibility...
 * 
 * @author Fabian Birzele
 */
public class Protein implements Sequence, Serializable
{
	/**
	 * Serial Version UID created with the following instance variables:
	 * 
	 * private AminoAcid[] aminoAcids;
	 * private HashMap properties;
	 */
	private static final long serialVersionUID = -6978959387880198678L;
	
	// standard property ids
	public static String PROPERTY_PDB_ID = "PDBID";
	public static String PROPERTY_PDB_FILE = "PDBFILE";
	public static String PROPERTY_PROTEIN_ID = "PROTEINID";
	public static String PROPERTY_ENSEMBL_GENE_ID = "ENSEMBL_GENE_ID";
	public static String PROPERTY_SEQUENCE = "SEQUENCE";
	public static String PROPERTY_SCOP_CLASS = "SCOPCLASS";
	public static String PROPERTY_SCOP_INFORMATION ="SCOPINFORMATION";
	public static String PROPERTY_SOLVENT_ACCESSIBILITY_AREA = "SAAREA";
	public static String PROPERTY_SOLVENT_ACCESSIBILITY_TWO_STATE = "SATWOSTATE";
	public static String PROPERTY_SOLVENT_ACCESSIBILITY_PREDICTION = "SAPREDICTION";
	public static String PROPERTY_PSSM = "PSSM";
	public static String PROPERTY_VORONOI_TESSELLATION = "VORONOI_TESSELLATION";
	public static String PROPERTY_VORONOI_ALL_ATOM_TESSELLATION = "VORONOI_ALL_ATOM_TESSELLATION";
	public static String PROPERTY_PROSITE_PATTERN_MATCHES = "PROSITE_PATTERN_MATCH_LIST";
	public static String PROPERTY_POTENTIAL_COMPONENTS = "POTENTIAL_COMPONENTS";
	public static String PROPERTY_3D_PATTERNS = "3D_PATTERNS";
	public static String PROPERTY_DSSP = "DSSP";
	
	// instance variables
	private AminoAcid[] aminoAcids;
	private HashMap<String, Object> properties;

	/**
	 * Initializes the protein with its id and its pdb id. Those ids can be the same but do not necessarily be the same...
	 * 
	 * @param proteinID
	 * @param pdbID
	 */
	public Protein(String proteinID, String pdbID)
	{
		addProperty(Protein.PROPERTY_PDB_ID, pdbID);
		addProperty(Protein.PROPERTY_PROTEIN_ID, proteinID);
	}
	
	/**
	 * Method initializes a protein with its id and its pdb id. Those ids can be the same but do not necessarily be the same...
	 * Also the protein is initialized with the specified sequence.
	 * 
	 * @param proteinID
	 * @param pdbID
	 * @param sequence
	 */
	public Protein(String proteinID, String pdbID, String sequence)
	{
		addProperty(Protein.PROPERTY_PDB_ID, pdbID);
		addProperty(Protein.PROPERTY_PROTEIN_ID, proteinID);
		
		// replace * if exists...
		if(sequence.endsWith("*"))
			sequence = sequence.replaceAll("\\*", "");
		
		char[] seq = sequence.toCharArray();
		aminoAcids = new AminoAcid[seq.length];
		
		for(int i=0; i<seq.length; i++)
		{
			aminoAcids[i] = new AminoAcid(seq[i], 'C', i);
			aminoAcids[i].setResiduePositionInSequence(i);
		}
		
		addProperty(Protein.PROPERTY_SEQUENCE, sequence);
	}
	
	/**
	 * Method returns the index of the amino acid that is assigned to the specified residue number in the 
	 * pdb file. If such a residue does not exist, null is returned...
	 * 
	 * @param residueNumber
	 * @return
	 */
	public int getAminoAcidIndexCorrespondingToResidueNumberInPDB(int residueNumber)
	{
		for(int i=0; i<aminoAcids.length; i++)
		{
			if(aminoAcids[i].getResidueNumber() == residueNumber)
				return i;
		}
		
		return -1;
	}
	
	/**
	 * Method returns the indes of the amino acid that is assigned to the specified residue number in the pdb file.
	 * If such a residue does not exist, null is returned
	 * @param residueNumber
	 * @param chain
	 * @return
	 */
	public int getAminoAcidIndexCorrespondingToResidueNumberOnChainInPDB(int residueNumber, char chain)
	{
		for(int i=0; i<aminoAcids.length; i++)
		{
			if(aminoAcids[i].getResidueNumber() == residueNumber && aminoAcids[i].getChainID() == chain)
				return i;
		}
		
		//System.out.println("Residue " + residueNumber + " not found on chain " + chain);
		
		return -1;
	}
	
	/**
	 * Method sets the chain id for all amino acids of the protein
	 * @param chainID
	 */
	public void setChainID(char chainID)
	{
		for(int i=0; i<aminoAcids.length; i++)
		{
			aminoAcids[i].setChainID(chainID);
		}
	}
	
	/**
	 * Method returns a vector containing all chain ids found in the protein
	 * 
	 * @return
	 */
	public Vector<Character> getChainsOfProtein()
	{
		Vector<Character> chains = new Vector<Character>();
		
		for(int i=0; i<aminoAcids.length; i++)
		{
			if(!chains.contains(new Character(aminoAcids[i].getChainID())))
			{
				chains.add(new Character(aminoAcids[i].getChainID()));
			}
		}
		
		return chains;
	}
	
	/**
	 * Method returns the start position of the specified chain in the protein amino acid array.
	 * If the chain does not exist, -1 is returned.
	 * 
	 * @param chain
	 * @return
	 */
	public int getStartOfChainInProteinSequence(char chain)
	{
		for(int i=0; i<aminoAcids.length; i++)
		{
			if(aminoAcids[i].getChainID() == chain)
			{
				return i;
			}
		}
		
		return -1;
	}
	
	/**
	 * Method returns the end position of the specified chain in the protein amino acid array.
	 * If the chain does not exist, -1 is returned.
	 * 
	 * @param chain
	 * @return
	 */
	public int getEndOfChainInProteinSequence(char chain)
	{
		boolean chainFound = false;
		
		// if there is only one chain => end is the last amino acid
		if(getChainsOfProtein().size() == 1)
			return aminoAcids.length-1;
		
		int lastResidueOfChain = -1;
		
		for(int i=0; i<aminoAcids.length; i++)
		{
			if(aminoAcids[i].getChainID() == chain)
			{
				chainFound = true;
				lastResidueOfChain = i;
			}
			// new chain after wanted chain
			else if(chainFound && aminoAcids[i].getChainID() != chain)
				return i-1;
		}
		
		if(chainFound)
			return lastResidueOfChain;
		
		return -1;
	}
	
	/**
	 * Method returns the part of the secondary structure sequence of the protein that belongs to the specified
	 * chain id. If the chain does not exist in the protein, null is returned!
	 * 
	 * @param chain
	 * @return
	 */
	public String getChainSecondaryStructureSequence(char chain)
	{
		StringBuffer buffer = new StringBuffer();
		boolean chainExists = false;
		
		for(int i=0; i<aminoAcids.length; i++)
		{
			if(aminoAcids[i].getChainID() == chain)
			{
				buffer.append(aminoAcids[i].getSecondaryStructureElement());
				
				chainExists = true;
			}
		}
		
		if(chainExists)
			return buffer.toString();
		else
			return null;
	}
	
	/**
	 * Method returns the part of the sequence of the protein that belongs to the specified
	 * chain id. If the chain does not exist in the protein, null is returned!
	 * 
	 * @param chain
	 * @return
	 */
	public String getChainSequence(char chain)
	{
		StringBuffer buffer = new StringBuffer();
		boolean chainExists = false;
		
		for(int i=0; i<aminoAcids.length; i++)
		{
			if(aminoAcids[i].getChainID() == chain)
			{
				buffer.append(aminoAcids[i].getOneLetterCode());
				chainExists = true;
			}
		}
		
		if(chainExists)
			return buffer.toString();
		else
			return null;
	}
	
	/**
	 * Method returns the amino acids of the protein that belong to the corresponding chain
	 * 
	 * @param chain
	 * @return
	 */
	public AminoAcid[] getChainAminoAcids(char chain)
	{
		Vector<AminoAcid> aminoAcidsVector = new Vector<AminoAcid>();
		
		for(int i=0; i<aminoAcids.length; i++)
		{
			if(aminoAcids[i].getChainID() == chain)
			{
				aminoAcidsVector.add(aminoAcids[i]);
			}
		}
		
		return (AminoAcid[])aminoAcidsVector.toArray(new AminoAcid[aminoAcidsVector.size()]);
	}

	/**
	 * @return Returns the pdbID
	 */
	public String getPdbID()
	{
		return (String)properties.get(Protein.PROPERTY_PDB_ID);
	}
	
	/**
	 * Method returns the pdb file
	 * 
	 * @return
	 */
	public File getPdbFile()
	{
		return (File)properties.get(Protein.PROPERTY_PDB_FILE);
	}
	
	/**
	 * Method sets the pdb file of the protein
	 * 
	 * @param pdbFilePath
	 */
	public void setPdbFile(File pdbFile)
	{
		addProperty(Protein.PROPERTY_PDB_FILE, pdbFile);
	}
	
	/**
	 * @return Returns the proteinID.
	 */
	public String getProteinID()
	{
		return (String)properties.get(Protein.PROPERTY_PROTEIN_ID);
	}
	/**
	 * @return Returns the sequence.
	 */
	public String getSequence()
	{
		return (String)properties.get(Protein.PROPERTY_SEQUENCE);
	}
	
	/**
	 * Method is required by the Sequence interface. The one letter code of the amino acid
	 * at the specified position is returned
	 * 
	 * @param position
	 * @return
	 */
	public String getSequenceElement(int position)
	{
		return "" + aminoAcids[position].getOneLetterCode();
	}
	
	/**
	 * Method returns the amino acid stored at the specified position
	 * 
	 * @param position
	 * @return
	 */
	public Object getSequenceObjectAt(int position)
	{
		return aminoAcids[position];
	}

	/**
	 * Method returns the length of the protein, i.e. the number of amino acids
	 * 
	 * @return
	 */
	public int getSequenceLength()
	{
		return aminoAcids.length;
	}

	public String getSequenceWithModifiedAminoAcids() {
		StringBuffer sequence = new StringBuffer("");
		for ( AminoAcid aminoAcid : aminoAcids ) {
			if ( aminoAcid.isModfiedAminoAcid() ) {
				sequence.append( aminoAcid.getBaseOneLetterCode() );
			}
			else {
				sequence.append( aminoAcid.getOneLetterCode() );
			}
		}
		return sequence.toString();
	}
	
	/**
	 * Method returns the ID of the protein
	 */
	public String getID()
	{
		return getProteinID();
	}
	
	/**
	 * Method returns the secondary structure string of the protein
	 * 
	 * @return
	 */
	public String getSecondaryStructure()
	{
		StringBuffer buffer = new StringBuffer();
		
		for(int i=0; i<aminoAcids.length; i++)
		{
			buffer.append(aminoAcids[i].getSecondaryStructureElement());
		}
		
		return buffer.toString();
	}

	/**
	 * Method sets the amino acid array and the amino acid sequence
	 * 
	 * @param aminoAcids
	 * @param sequence
	 */
	public void setAminoAcids(AminoAcid[] aminoAcids, String sequence)
	{
		this.aminoAcids = aminoAcids;
		
		for(int i=0; i<this.aminoAcids.length; i++)
		{
			this.aminoAcids[i].setResiduePositionInSequence(i);
		}
		
		addProperty(Protein.PROPERTY_SEQUENCE, sequence);
	}
	
	/**
	 * This method removes the residue at the specified position from the protein.
	 * Please note that this method should be used with care. Alignment instances
	 * or other instances that reference only positions in the amino acid array 
	 * might cause problems afterwards.
	 * 
	 * @param position
	 */
	public void removeResidue(int position)
	{
		AminoAcid[] newSequence = new AminoAcid[aminoAcids.length - 1];

		System.arraycopy(aminoAcids, 0, newSequence, 0, position);
		System.arraycopy(aminoAcids, position+1, newSequence, position, newSequence.length-position);
		
		aminoAcids = newSequence;
		
		StringBuffer newSeqString = new StringBuffer();
		
		for(int i=0; i<aminoAcids.length; i++)
		{
			newSeqString.append(aminoAcids[i].getOneLetterCode());
		}
		
		addProperty(Protein.PROPERTY_SEQUENCE, newSeqString.toString());
	}
	
	/**
	 * This method removes a set of residues from the protein. It should be used if more
	 * than one residue is removed to make sure that no index errors occur due to the
	 * changing sequence lenght in case that more residues are removed in a row.
	 * 
	 * @param positions
	 */
	public void removeResidues(Vector<Integer> positions)
	{
		int removed = 0;
		for(Integer removal : positions)
		{
			aminoAcids[removal] = null;
			removed++;
		}
		
		AminoAcid[] newSequence = new AminoAcid[aminoAcids.length - removed];
		int newPosition = 0;
		
		StringBuffer newSeqString = new StringBuffer();
		
		for(int i=0; i<aminoAcids.length; i++)
		{
			if(aminoAcids[i] != null)
			{
				newSequence[newPosition] = aminoAcids[i];
				newSeqString.append(aminoAcids[i].getOneLetterCode());
				aminoAcids[i].setResiduePositionInSequence(newPosition);
				newPosition++;
			}
		}
		
		aminoAcids = newSequence;
		
		addProperty(Protein.PROPERTY_SEQUENCE, newSeqString.toString());
	}

	/**
	 * Method returns the amino acid array
	 * 
	 * @return
	 */
	public AminoAcid[] getAminoAcids()
	{
		return aminoAcids;
	}
	
	/**
	 * Method adds the given property to the property map of the protein instance
	 * using the given name as key. If there is already a property with that name, 
	 * the old property is overwritten.
	 * 
	 * @param name
	 * @param property
	 */
	public void addProperty(String name, Object property)
	{
		if(properties == null)
			properties = new HashMap<String, Object>();
		
		properties.put(name, property);
	}
	
	/**
	 * Method returns the property object associated with the specified name. If no property
	 * with that name exists, null is returned.
	 * 
	 * @param name
	 * @return
	 */
	public Object getProperty(String name)
	{
		return properties.get(name);
	}
	
	/**
	 * This method renumbers the amino acids of the protein such that their residue numbers
	 * correspond to their positions in the amino acid array of the protein instance.
	 */
	public void renumberResidues()
	{
		for(int i=0;i<aminoAcids.length; i++)
		{
			aminoAcids[i].setResidueNumber(i);
		}
	}
	
	/**
	 * This method clones the protein structure information of the protein, i.e. its amino acid sequence.
	 * It does not clone all properties assigned to the protein. Those are only transfered by their 
	 * references
	 * 
	 * @return
	 */
	public Protein cloneStructure()
	{
		AminoAcid[] clonedAcids = new AminoAcid[this.aminoAcids.length];
		
		for(int i=0; i<aminoAcids.length; i++)
		{
			clonedAcids[i] = aminoAcids[i].clone();
		}
		
		Protein clone = new Protein(this.getProteinID(), this.getPdbID());
		clone.setAminoAcids(clonedAcids, this.getSequence());
		
		for(String p : properties.keySet())
		{
			clone.addProperty(p, properties.get(p));
		}
		
		return clone;
	}

	/**
	 * Returns the Ca-coordinates, a set of x-y-z coordinates per line.
	 * 
	 * @return String the Ca-coordinates, a set of x-y-z coordinates per line
	 */
	public String toCaCoordinates()
	{
		StringBuffer caCoordinates = new StringBuffer();
		for (int i = 0; i < aminoAcids.length; i++)
		{
			caCoordinates.append(aminoAcids[i].getXCACoordinate());
			caCoordinates.append(' ');
			caCoordinates.append(aminoAcids[i].getYCACoordinate());
			caCoordinates.append(' ');
			caCoordinates.append(aminoAcids[i].getZCACoordinate());
			caCoordinates.append('\n');
			
		}
		return caCoordinates.toString();
	}

	/**
	 * Method returns a string representation of the protein which looks as follows:
	 * 
	 * > proteinID pdbFileName
	 * AminoAcid SecondaryStructureElement XCA YCA ZCA
	 * NextAminoAcid...
	 */
	public String toString()
	{
		StringBuffer finalRepresentation = new StringBuffer();

		finalRepresentation.append("> ");
		finalRepresentation.append((String)properties.get(Protein.PROPERTY_PROTEIN_ID));
		finalRepresentation.append("\t");
		finalRepresentation.append((String)properties.get(Protein.PROPERTY_PDB_ID));
		finalRepresentation.append("\n");

		for (int i = 0; i < aminoAcids.length; i++)
		{	
			finalRepresentation.append(i+1);
			finalRepresentation.append("\t");
			finalRepresentation.append(aminoAcids[i].getOneLetterCode());
			finalRepresentation.append("\t");
			finalRepresentation.append(aminoAcids[i].getSecondaryStructureElement());
			finalRepresentation.append("\t");
			finalRepresentation.append(aminoAcids[i].getChainID());
			finalRepresentation.append("\t");
			finalRepresentation.append(aminoAcids[i].getXCACoordinate());
			finalRepresentation.append("\t");
			finalRepresentation.append(aminoAcids[i].getYCACoordinate());
			finalRepresentation.append("\t");
			finalRepresentation.append(aminoAcids[i].getZCACoordinate());
			finalRepresentation.append("\n");
		}

		return finalRepresentation.toString();
	}
	
	/**
	 * Method returns an iterator object on the protein amino
	 * acid sequence. null references in the alignment are ignored.
	 * 
	 * @return Iterator
	 */
	public Iterator<?> getSequenceIterator()
	{
		return new SequenceIterator(aminoAcids, this);
	}


	/**
         * This Class implements an iterator for a given sequence of amino acids. 
	 * null - references in the alignment are ignored so always the next amino
	 * acid is returned.
	 * 
	 * @author Fabian Birzele
	 *
	 */
	private class SequenceIterator implements Iterator<Object>
	{
		private AminoAcid[] sequence;
		private Protein protein;
		private int currentPosition;

		/**
		 * Initializes the alignment iterator
		 * 
		 * @param sequence
		 * @param alignment
		 */
		public SequenceIterator(AminoAcid[] sequence, Protein protein)
		{
			this.sequence = sequence;
			this.protein = protein;

			currentPosition = -1;
		}

		/**
		 * Tests if there is a next nun null position in the alignment
		 * 
		 * @see java.util.Iterator#hasNext()
		 */
		public boolean hasNext()
		{
			if (currentPosition < sequence.length)
			{
				for (int i = currentPosition + 1; i < sequence.length; i++)
				{
					if (sequence[i] != null)
					{
						return true;
					}
				}
			}

			return false;
		}

		/**
		 * returns the next non null amino acid in the alignment
		 * 
		 * @see java.util.Iterator#next()
		 */
		public Object next()
		{
			for (int i = currentPosition + 1; i < sequence.length; i++)
			{
				if (sequence[i] != null)
				{
					currentPosition = i;

					return sequence[i];
				}
			}

			return null;
		}

		/**
		 * Method removed the current amino acid form the alignment (in query and template)
		 * 
		 * @see java.util.Iterator#remove()
		 */
		public void remove()
		{
			// remove amino acid from alignment
			protein.removeResidue(currentPosition);

			// remove amino acid from iterator sequence array also
			AminoAcid[] tempSequence = new AminoAcid[sequence.length - 1];

			System.arraycopy(sequence, 0, tempSequence, 0, currentPosition);
			System.arraycopy(sequence, currentPosition+1, tempSequence, currentPosition, tempSequence.length-currentPosition);
		
			sequence = tempSequence;

			currentPosition--;
		}
	}
}
