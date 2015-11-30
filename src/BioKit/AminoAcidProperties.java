package BioKit;


public class AminoAcidProperties
{
	/**
	 * Method returns the hydrophobicity of an amino acid with the specified one letter code based
	 * on the Kyte and Doolittle table.
	 * 
	 * @param oneLetterCode
	 * @return
	 */
	public static float getKyteDoolittleHydrophobicity(char oneLetterCode)
	{
		// make sure that one letter code comes in upper case form
		char olc = Character.toUpperCase(oneLetterCode);

		if(olc =='A')
		{
			return 1.8f;
		}
		else if(olc == 'I')
		{
			return 4.5f;
		}
		else if(olc == 'L')
		{
			return 3.8f;
		}
		else if(olc == 'M')
		{
			return 1.9f;
		}
		else if(olc == 'F')
		{
			return 2.8f;
		}
		else if(olc == 'P')
		{
			return -1.6f;
		}
		else if(olc == 'V')
		{
			return  4.2f;
		}
		else if(olc == 'R')
		{
			return -4.5f;
		}
		else if(olc == 'D')
		{
			return -3.5f;
		}
		else if(olc == 'E')
		{
			return -3.5f;
		}
		else if(olc == 'K')
		{
			return -3.9f;
		}
		else if(olc == 'N')
		{
			return -3.5f;
		}
		else if(olc == 'C')
		{
			return 2.5f;
		}
		else if(olc == 'Q')
		{
			return -3.5f;
		}
		else if(olc == 'H')
		{
			return -3.2f;
		}
		else if(olc == 'S')
		{
			return -0.8f;
		}
		else if(olc == 'T')
		{
			return -0.7f;
		}
		else if(olc == 'W')
		{
			return  -0.9f;
		}
		else if(olc == 'Y')
		{
			return -1.3f;
		}
		else if(olc == 'G')
		{
			return -0.4f;
		}

		return 0;
	}
	
	/**
	 * Method returns the the maximal solvent accessiblity of the amino acid in an extended conformation Ala - X - Ala
	 * in A^2 as described in 
	 * 
	 *  ASAView: Database and tool for solvent accessibility representation in proteins by Shandar Ahmad, BMC Bioinformatics 2004, 5, 51
	 *  
	 *  110.2 (Ala), 144.1 (Asp),140.4 (Cys), 174.7 (Glu), 200.7 (Phe), 78.7 (Gly), 181.9 (His), 185.0 (Ile), 205.7 (Lys), 183.1 (Leu), 200.1 (Met),
	 *  146.4 (Asn), 141.9 (Pro), 178.6 (Gln), 229.0 (Arg), 117.2 (Ser), 138.7 (Thr), 153.7 (Val), 240.5 (Trp), and 213.7 (Tyr)
	 * 
	 * @param oneLetterCode
	 * @return
	 */
	public static float getMaximalSolventAccessibility(char oneLetterCode)
	{
		// make sure that one letter code comes in upper case form
		char olc = Character.toUpperCase(oneLetterCode);

		if(olc =='A')
		{
			return 110.2f;
		}
		else if(olc == 'I')
		{
			return 185.0f;
		}
		else if(olc == 'L')
		{
			return 183.1f;
		}
		else if(olc == 'M')
		{
			return 200.1f;
		}
		else if(olc == 'F')
		{
			return 200.7f;
		}
		else if(olc == 'P')
		{
			return 141.9f;
		}
		else if(olc == 'V')
		{
			return  153.7f;
		}
		else if(olc == 'R')
		{
			return 229.0f;
		}
		else if(olc == 'D')
		{
			return 144.1f;
		}
		else if(olc == 'E')
		{
			return 174.7f;
		}
		else if(olc == 'K')
		{
			return 205.7f;
		}
		else if(olc == 'N')
		{
			return 146.4f;
		}
		else if(olc == 'C')
		{
			return 140.4f;
		}
		else if(olc == 'Q')
		{
			return 178.6f;
		}
		else if(olc == 'H')
		{
			return 181.9f;
		}
		else if(olc == 'S')
		{
			return 117.2f;
		}
		else if(olc == 'T')
		{
			return 138.7f;
		}
		else if(olc == 'W')
		{
			return 240.5f;
		}
		else if(olc == 'Y')
		{
			return 213.7f;
		}
		else if(olc == 'G')
		{
			return 78.7f;
		}

		return 0;
	}
	
	/**
	 * This method translates the standard three letter amino acid codes into
	 * standard one letter amino acid codes.
	 * 
	 * @param threeLetterCode
	 * @return char
	 */
	public static char threeLetterToOneLetterCode(String threeLetterCode)
	{
		if(threeLetterCode.equalsIgnoreCase("ALA"))
		{
			return 'A';
		}
		else if(threeLetterCode.equalsIgnoreCase("VAL"))
		{
			return 'V';
		}
		else if(threeLetterCode.equalsIgnoreCase("LEU"))
		{
			return 'L';
		}
		else if(threeLetterCode.equalsIgnoreCase("ILE"))
		{
			return 'I';
		}
		else if(threeLetterCode.equalsIgnoreCase("PRO"))
		{
			return 'P';
		}
		else if(threeLetterCode.equalsIgnoreCase("TRP"))
		{
			return 'W';
		}
		else if(threeLetterCode.equalsIgnoreCase("PHE"))
		{
			return 'F';
		}
		else if(threeLetterCode.equalsIgnoreCase("MET"))
		{
			return 'M';
		}
		else if(threeLetterCode.equalsIgnoreCase("GLY"))
		{
			return 'G';
		}
		else if(threeLetterCode.equalsIgnoreCase("SER"))
		{
			return 'S';
		}
		else if(threeLetterCode.equalsIgnoreCase("THR"))
		{
			return 'T';
		}
		else if(threeLetterCode.equalsIgnoreCase("TYR"))
		{
			return 'Y';
		}
		else if(threeLetterCode.equalsIgnoreCase("CYS"))
		{
			return 'C';
		}
		else if(threeLetterCode.equalsIgnoreCase("ASN"))
		{
			return 'N';
		}
		else if(threeLetterCode.equalsIgnoreCase("GLN"))
		{
			return 'Q';
		}
		else if(threeLetterCode.equalsIgnoreCase("LYS"))
		{
			return 'K';
		}
		else if(threeLetterCode.equalsIgnoreCase("ARG"))
		{
			return 'R';
		}
		else if(threeLetterCode.equalsIgnoreCase("HIS"))
		{
			return 'H';
		}
		else if(threeLetterCode.equalsIgnoreCase("ASP"))
		{
			return 'D';
		}
		else if(threeLetterCode.equalsIgnoreCase("GLU"))
		{
			return 'E';
		}
		else
		{
			return 'X';
		}
	}
	
	/**
	 * This method takes the standard three letter code of an amino acid and returns an
	 * array of ordered atom names that belong to the amino acid (including side chain 
	 * atoms). Names refer to names used in the pdb...
	 * 
	 * @param threeLetterCode
	 * @return char
	 */
	public static String[] threeLetterToAtomNameArray(String threeLetterCode)
	{
		if(threeLetterCode.equalsIgnoreCase("ALA"))
		{
			return new String[]{"N", "CA", "C", "O", "CB"};
		}
		else if(threeLetterCode.equalsIgnoreCase("VAL"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG1", "CG2"};
		}
		else if(threeLetterCode.equalsIgnoreCase("LEU"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"};
		}
		else if(threeLetterCode.equalsIgnoreCase("ILE"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"};
		}
		else if(threeLetterCode.equalsIgnoreCase("PRO"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG", "CD"};
		}
		else if(threeLetterCode.equalsIgnoreCase("TRP"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"};
		}
		else if(threeLetterCode.equalsIgnoreCase("PHE"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"};
		}
		else if(threeLetterCode.equalsIgnoreCase("MET"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG", "SD", "CE"};
		}
		else if(threeLetterCode.equalsIgnoreCase("GLY"))
		{
			return new String[]{"N", "CA", "C", "O"};
		}
		else if(threeLetterCode.equalsIgnoreCase("SER"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "OG"};
		}
		else if(threeLetterCode.equalsIgnoreCase("THR"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "OG1", "CG2"};
		}
		else if(threeLetterCode.equalsIgnoreCase("TYR"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"};
		}
		else if(threeLetterCode.equalsIgnoreCase("CYS"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "SG"};
		}
		else if(threeLetterCode.equalsIgnoreCase("ASN"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"};
		}
		else if(threeLetterCode.equalsIgnoreCase("GLN"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"};
		}
		else if(threeLetterCode.equalsIgnoreCase("LYS"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"};
		}
		else if(threeLetterCode.equalsIgnoreCase("ARG"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"};
		}
		else if(threeLetterCode.equalsIgnoreCase("HIS"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"};
		}
		else if(threeLetterCode.equalsIgnoreCase("ASP"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"};
		}
		else if(threeLetterCode.equalsIgnoreCase("GLU"))
		{
			return new String[]{"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"};
		}
		else
		{
			return new String[]{};
		}
	}
	
	/**
	 * Method returns true if a specified atom is part of the standard atom set of the amino acid.
	 * Otherwise false is returned.
	 * 
	 * @param aminoAcidThreeLetterCode
	 * @param atomName
	 * @return
	 */
	public static boolean isStandardAminoAcidAtom(String aminoAcidThreeLetterCode, String atomName)
	{
		String[] atoms = threeLetterToAtomNameArray(aminoAcidThreeLetterCode);
		
		for(int i=0;i<atoms.length;i++)
		{
			if(atomName.equals(atoms[i]))
				return true;
		}
		
		return false;
	}
	
	/**
	 * Method takes the standard amino acid one letter code and translates it to the standard
	 * upper case three letter code.
	 * 
	 * @param oneLetterCode
	 * @return String
	 */
	public static String oneLetterToThreeLetterCode(char oneLetterCode)
	{
		// make sure that one letter code comes in upper case form
		char olc = Character.toUpperCase(oneLetterCode);

		if(olc =='A')
		{
			return "ALA";
		}
		else if(olc == 'I')
		{
			return "ILE";
		}
		else if(olc == 'L')
		{
			return "LEU";
		}
		else if(olc == 'M')
		{
			return "MET";
		}
		else if(olc == 'F')
		{
			return "PHE";
		}
		else if(olc == 'P')
		{
			return "PRO";
		}
		else if(olc == 'V')
		{
			return "VAL";
		}
		else if(olc == 'R')
		{
			return "ARG";
		}
		else if(olc == 'D')
		{
			return "ASP";
		}
		else if(olc == 'E')
		{
			return "GLU";
		}
		else if(olc == 'K')
		{
			return "LYS";
		}
		else if(olc == 'N')
		{
			return "ASN";
		}
		else if(olc == 'C')
		{
			return "CYS";
		}
		else if(olc == 'Q')
		{
			return "GLN";
		}
		else if(olc == 'H')
		{
			return "HIS";
		}
		else if(olc == 'S')
		{
			return "SER";
		}
		else if(olc == 'T')
		{
			return "THR";
		}
		else if(olc == 'W')
		{
			return "TRP";
		}
		else if(olc == 'Y')
		{
			return "TYR";
		}
		else if(olc == 'G')
		{
			return "GLY";
		}

		return "XXX";
	}
	
	/**
	 * Method returns the property groups associated with an amino acid (specified by the given
	 * one letter code). The property groups follow the Taylor classification scheme and contain:
	 * Small, Hydrophobic, Aliphatic, Positive, Polar, Negative and Aromatic. An amino acid can
	 * belong to more than one group.
	 * 
	 * @param oneLetterCode
	 * @return
	 */
	public static String[] getPropertyGroupsOfAminoAcid(char oneLetterCode)
	{
		char olc = Character.toUpperCase(oneLetterCode);

		if (olc == 'A')
		{
			return new String[]{"Small[STDNGA]", "Hydrophobic[VILMFA]"};
		}
		else if (olc == 'I')
		{
			return new String[]{"Aliphatic[VIL]", "Hydrophobic[VILMFA]"};
		}
		else if (olc == 'L')
		{
			return new String[]{"Aliphatic[VIL]", "Hydrophobic[VILMFA]"};
		}
		else if (olc == 'M')
		{
			return new String[]{"Hydrophobic[VILMFA]"};
		}
		else if (olc == 'F')
		{
			return new String[]{"Aromatic[FYW]", "Hydrophobic[VILMFA]"};
		}
		else if (olc == 'P')
		{
			return new String[]{};
		}
		else if (olc == 'V')
		{
			return new String[]{"Aliphatic[VIL]", "Hydrophobic[VILMFA]"};
		}
		else if (olc == 'R')
		{
			return new String[]{"Positive[HKR]", "Polar[DENQRS]"};
		}
		else if (olc == 'D')
		{
			return new String[]{"Negative[DE]", "Small[STDNGA]",  "Polar[DENQRS]"};
		}
		else if (olc == 'E')
		{
			return new String[]{"Negative[DE]",  "Polar[DENQRS]"};
		}
		else if (olc == 'K')
		{
			return new String[]{"Positive[HKR]"};
		}
		else if (olc == 'N')
		{
			return new String[]{"Small[STDNGA]",  "Polar[DENQRS]"};
		}
		else if (olc == 'C')
		{
			return new String[]{};
		}
		else if (olc == 'Q')
		{
			return new String[]{ "Polar[DENQRS]"};
		}
		else if (olc == 'H')
		{
			return new String[]{"Positive[HKR]"};
		}
		else if (olc == 'S')
		{
			return new String[]{"Small[STDNGA]",  "Polar[DENQRS]"};
		}
		else if (olc == 'T')
		{
			return new String[]{"Small[STDNGA]"};
		}
		else if (olc == 'W')
		{
			return new String[]{"Aromatic[FYW]"};
		}
		else if (olc == 'Y')
		{
			return new String[]{"Aromatic[FYW]"};
		}
		else if (olc == 'G')
		{
			return new String[]{"Small[STDNGA]"};
		}
		
		return new String[]{};
	}
}
