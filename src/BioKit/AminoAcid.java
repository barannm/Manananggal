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

import java.io.Serializable;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;


/**
 * This class represents one amino acid. The class must be able to store
 * 1D (sequence), 2D (secondary structure information) as well as 3D
 * (x, y and z coordinates) information.
 * 
 * The standard x, y and z coordinates are the CA coordinates of the
 * amino acid stored in the pdb file. The structure supports adding 
 * more atoms if wanted.
 * 
 * @author Fabian Birzele
 */
public class AminoAcid implements Comparable<AminoAcid>, Serializable
{
	/**
	 * Serial Version UID created on 21.04.2006. The class contained the following instance variables.
	 * 
	 * private float xCACoordinate;
	 * private float yCACoordinate;
	 * private float zCACoordinate;
	 * private int residueNumber;
	 * private char oneLetterCode;
	 * private char secondaryStructureElement;
	 * private HashMap<String, Atom> additionalAtoms;
	 */
	private static final long serialVersionUID = 2924713600229221989L;

	// stores the residue number of the amino acid. This is NOT necessarily equivalent to the position of the amino acid in the 
	// internal protein representation, i.e. the amino acid array of the protein.
	private int residueNumber;
	
	// this variable stores the position of the residue in the amino acid array of the protein which holds this amino acid.
	private int residuePositionInArray;

	// one letter code of the amino acid
	private char oneLetterCode;
	
	// chain id of the amino acid
	private char chainID = ' ';

	// Secondary structure assignment. Please notice that in the case of
	// the query this assignment is only a prediction! 
	private char secondaryStructureElement;

    //  coordinates of the amino acid in the 3D structure. Please notice
	// that in the case of the template those coordinates are real pdb
	// coordinates, in the case of the query those coordinates are derived
	// from the corresponding aligned amino acid in the template
    private Atom CAAtom = null;

	// Hash map holding all atoms
	private HashMap<String, Atom> atoms = new HashMap<String, Atom>();

	// information about a pattern the amino acid is described in
	// e.g. a PROSITE pattern
	// using the AC to describe the pattern
	// & the active site residue within the pattern (indexing from 0)
	// e.g. [PS00135,0]
	private String[] patternOccurrenceInfo;
	
	// solvent accessibility of amino acid
	private float solventAccessibility;
	
        private boolean is_heteroatom = false;
        private char baseOneLetterCode = 'X';
        private String residueName = null; // name of the residue if it is heteroatom or modified residue

        public void setHeteroAtomName( String residueName, String backmappedOneLetterCode ) {
            this.is_heteroatom = true;
            this.residueName = residueName;
            if ( backmappedOneLetterCode != null ) {
            	baseOneLetterCode = backmappedOneLetterCode.charAt(0);
            }
        }

        public String getResidueName()
        {
            return ( residueName != null ) ? residueName : getThreeLetterCode();
        }

        public boolean isHeteroAtom()
        {
            return is_heteroatom;
        }

        public boolean isModfiedAminoAcid()
        {
            return oneLetterCode == 'X' && baseOneLetterCode != 'X';
        }

        public char getBaseOneLetterCode()
        {
            return baseOneLetterCode;
        }
        
        public char getOneLetterCode( boolean base ) {
        	if ( oneLetterCode == 'X' ) {
        		if ( base ) {
        			return baseOneLetterCode;
        		}
        	}
       		return oneLetterCode;
       	}

        /**
	 * Constructor, initializes the variables of an amino acid when all the information
	 * is available
	 * 
	 * @param xCoordinate					pdb x CA coordinate
	 * @param yCoordinate					pdb y CA coordinate
	 * @param zCoordinate					pdb z CA coordinate
	 * @param oneLetterCode				amino acid one letter code
	 * @param secondaryStructureElement	amino acid secondary structure element assignment / prediction
	 */
        public AminoAcid(float xCoordinate, float yCoordinate, float zCoordinate, char oneLetterCode, char secondaryStructureElement, int residueNumber)
        {
            this(xCoordinate, yCoordinate, zCoordinate, "CA", oneLetterCode, secondaryStructureElement, residueNumber);
        }
	
        /**
	 * Constructor, initializes the variables of an amino acid when all the information
	 * is available
	 * 
	 * @param xCoordinate					pdb x CA coordinate
	 * @param yCoordinate					pdb y CA coordinate
	 * @param zCoordinate					pdb z CA coordinate
         * @param atomtype                                      pdb atomtype it defaults to CA but can be different (for example
         *                                                      when parsing heteroatom
	 * @param oneLetterCode				amino acid one letter code
	 * @param secondaryStructureElement	amino acid secondary structure element assignment / prediction
	 */
        public AminoAcid(float xCoordinate, float yCoordinate, float zCoordinate, String atomtype, char oneLetterCode, char secondaryStructureElement, int residueNumber)
	{
                addAtom(atomtype, xCoordinate, yCoordinate, zCoordinate);
		this.residueNumber = residueNumber;
		this.oneLetterCode = Character.toUpperCase(oneLetterCode);
		this.secondaryStructureElement = secondaryStructureElement;
	}

	/**
	 * Constructor initializes the variables of an amino acid when only 1D and 2D information is
	 * available. The 3D information might be set afterwards. The default coordinate value is 0.0
	 * for all amino acids
	 * 
	 * @param oneLetterCode
	 * @param secondaryStructureElement
	 */
	public AminoAcid(char oneLetterCode, char secondaryStructureElement, int residueNumber)
	{
                //do we need this???
                addAtom("CA", Float.MIN_VALUE, Float.MIN_VALUE, Float.MIN_VALUE);
		
		this.residueNumber = residueNumber;
		this.oneLetterCode = Character.toUpperCase(oneLetterCode);
		this.secondaryStructureElement = secondaryStructureElement;
	}
	
	public String toString()
	{
        Atom dummy = CAAtom;
        if (dummy == null)
        {
            Iterator<String> it = getAtomNames();
            if (it.hasNext())
            {
                //take the first atom
                dummy = getAtom(it.next());
            }
            
        }
        return oneLetterCode + ", "+residueNumber + " " + chainID +" ("+dummy+")" ;
                
	}

	/**
	 * Method returns the one letter code value
	 * 
	 * @return char
	 */
	public char getOneLetterCode()
	{
		return oneLetterCode;
	}

	/**
	 * Method returns the residue number. This number is NOT necessarily equivalent to the position of the amino acid in the 
	 * internal protein representation, i.e. the amino acid array of the protein.
	 * 
	 * @return int
	 */
	public int getResidueNumber()
	{
		return residueNumber;
	}
	
	/**
	 * Method sets the residue number of the amino acid
	 * 
	 * @param residueNumber
	 */
	public void setResidueNumber(int residueNumber)
	{
		this.residueNumber = residueNumber;
	}
	
	/**
	 * Method returns the residue position in the sequence
	 * 
	 * @return
	 */
	public int getResiduePositionInSequence()
	{
		return residuePositionInArray;
	}
	
	/**
	 * Method sets the residue position in the sequence
	 * 
	 * @param residuePositionInSequence
	 */
	public void setResiduePositionInSequence(int residuePositionInSequence)
	{
		this.residuePositionInArray = residuePositionInSequence;
	}
	
	/**
	 * Method sets the solvent accessibility of the amino acid. The Solvent accessibility in Angstroem should
	 * be normalized by the maximally possible solvent accessibility surface of the amino acid in Angstroem.
	 * This can be accessed in the class AminoAcidProperties and the method getMaximalSolventAccessibility.
	 * 
	 * @param accessibility
	 */
	public void setSolventAccessibility(float accessibility)
	{
		solventAccessibility = accessibility;
	}
	
	/**
	 * Method returns the solvent accessibility of the amino acid
	 * 
	 * @return
	 */
	public float getSolventAccessibility()
	{
		return solventAccessibility;
	}
	
	/**
	 * Method returns the chain id of the amino acid
	 * 
	 * @return
	 */
	public char getChainID()
	{
		return chainID;
	}

	/**
	 * Method sets the chain id of the amino acid
	 * 
	 * @param chainID
	 */
	public void setChainID(char chainID)
	{
                this.chainID = chainID;
	}

	/**
	 * Method compares two amino acids due to their residue numbers.
	 * 
	 * @param otherAminoAcid
	 * @return
	 */
	public int compareTo(AminoAcid otherAminoAcid)
	{
		int otherID = ((AminoAcid)otherAminoAcid).getResidueNumber();
		
		if(residueNumber < otherID)
			return -1;
		else if(residueNumber == otherID && ((AminoAcid)otherAminoAcid).chainID == chainID)
			return 0;
		
		return 1;
	}

	/**
	 * Method returns the three letter code to the one letter value of the 
	 * amino acid
	 * 
	 * @return String three letter code of the amino acid
	 */
	public String getThreeLetterCode()
	{
		return AminoAcidProperties.oneLetterToThreeLetterCode(oneLetterCode);
	}

	/**
	 * Method returns the secondary structure assignment of the amino acid
	 * 
	 * @return char
	 */
	public char getSecondaryStructureElement()
	{
		return secondaryStructureElement;
	}
	
	/**
	 * Method sets the secondary structure element to a new value
	 * @param element
	 */
	public void setSecondaryStructureElement(char element)
	{
		secondaryStructureElement = element;
	}

	/**
	 * Method returns the x coordinate of the amino acid
	 * 
	 * @return float
	 */
	public float getXCACoordinate()
	{       
                return (CAAtom!=null)?CAAtom.getXCoordinate():Float.MIN_VALUE;
	}

	/**
	 * Method returns the y coordinate of the amino acid
	 * 
	 * @return float
	 */
	public float getYCACoordinate()
	{
	    return (CAAtom!=null)?CAAtom.getYCoordinate():Float.MIN_VALUE;    
            
	}

	/**
	 * Method returns the z coordinate of the amino acid
	 * 
	 * @return float
	 */
	public float getZCACoordinate()
	{
	    return (CAAtom!=null)?CAAtom.getZCoordinate():Float.MIN_VALUE;
	}

        
        private void initCA()
        {
            if (CAAtom!=null) return;
            addAtom("CA", Float.MIN_VALUE, Float.MIN_VALUE, Float.MIN_VALUE);
            
        }

	/**
	 * Sets the xCoordinate.
	 * 
	 * @param xCoordinate The xCoordinate to set
	 */
	public void setXCACoordinate(float xCoordinate)
	{
            initCA();
            CAAtom.setXCoordinate(xCoordinate);
            
	}

	/**
	 * Sets the yCoordinate.
	 * 
	 * @param yCoordinate The yCoordinate to set
	 */
	public void setYCACoordinate(float yCoordinate)
	{
            initCA();
            CAAtom.setYCoordinate(yCoordinate);
            
	}

	/**
	 * Sets the zCoordinate.
	 * 
	 * @param zCoordinate The zCoordinate to set
	 */
	public void setZCACoordinate(float zCoordinate)
	{
            initCA();
            CAAtom.setZCoordinate(zCoordinate);
	    
	}


        public Collection<Atom> getAtoms()
        {
            return atoms.values();
        }
	/**
	 * Method returns an iterator instance of the atom names stored
	 * in this amino acid instance
	 * 
	 * @return Iterator<String>
	 */
	public Iterator<String> getAtomNames()
	{
            return atoms.keySet().iterator();
		
	}

	/**
	 * This method returns the corresponding Atom instance held by this amino 
	 * acid instance for the specified atom name. If the amino acid does not 
	 * contain this atom name, null is returned!
	 * 
	 * @param name
	 * @return Atom
	 */
	public Atom getAtom(String name)
	{
            return atoms.get(name);
	}

	/**
	 * Method adds an Atom (different than the CA atom) to the amino acid.
	 * 
	 * @param name
	 * @param xCoordinate
	 * @param yCoordinate
	 * @param zCoordinate
	 */
	public void addAtom(String name, float xCoordinate, float yCoordinate, float zCoordinate)
	{
		// initialize hash map when method is called first
                Atom a = new Atom(name, xCoordinate, yCoordinate, zCoordinate);
                atoms.put(name, a);
                if (name.equals("CA"))
                {
                    CAAtom = a;
                }
		
	}
	
	/**
	 * Method adds an Atom to the amino acid. The specified atom instance is copied, no reference is stored!
	 * 
	 * @param atom
	 */
	public void addAtom(Atom atom)
	{
		// initialize hash map when method is called first
		atoms.put(atom.getAtomName(), new Atom(atom.getAtomName(), atom.getXCoordinate(), atom.getYCoordinate(), atom.getZCoordinate()));
		
		if (atom.getAtomName().equals("CA"))
                {
                    CAAtom = atoms.get(atom.getAtomName());
                }
	}

	/**
	 * This method transfers the coordinates of the other amino acid to this amino 
	 * acid. If the amino acids are of the same type, all atoms are transfered.
	 * If they have different types, only the CA and CB information (if available)
	 * are transfered. If this amino acid is a Glycine, no CB information is transfered.
	 * 
	 * @param otherAcid
	 */
	public void transferAtoms(AminoAcid otherAcid)
	{
                 
                Atom otherAtom = otherAcid.getAtom("CA");
                if (otherAtom!=null)
                {
                    addAtom(otherAtom.getAtomName(), otherAtom.getXCoordinate(), otherAtom.getYCoordinate(), otherAtom.getZCoordinate());
                }
                
	
                // same Amino acids, transfer all available coordinates
                if (oneLetterCode == otherAcid.getOneLetterCode())
                {
                    Iterator<String> otherAtomsIterator = otherAcid.getAtomNames();

                    while (otherAtomsIterator.hasNext())
                    {
                        otherAtom = otherAcid.getAtom(otherAtomsIterator.next());
                        addAtom(otherAtom.getAtomName(), otherAtom.getXCoordinate(), otherAtom.getYCoordinate(), otherAtom.getZCoordinate());
                    }
                    
    
                }
                else
                // different amino acids, transfer only CB information, if available
		{
                        Atom otherCBAtom = otherAcid.getAtom("CB");

                        if (otherCBAtom != null)
                        {
                                addAtom(otherCBAtom.getAtomName(), otherCBAtom.getXCoordinate(), otherCBAtom.getYCoordinate(), otherCBAtom.getZCoordinate());
                        }
                }
		
	}

        public float calculateAtomDistance(AminoAcid otherAcid, String ownatomtype, String otheratomtype)
        {
            Atom thisAtom = this.getAtom(ownatomtype);
	    Atom otherAtom = otherAcid.getAtom(otheratomtype);

            if(thisAtom==null || otherAtom == null) 
                return Float.MAX_VALUE; //is this a good idea? show we not throw an exception??

            return thisAtom.calculateDistance(otherAtom);
            
        }
	/**
	 * Method calculates the distance of the CA atoms of this amino acid with the
	 * other specified amino acid
	 * 
	 * @param otherAcid
	 * @return float
	 */
	public float calculateCADistance(AminoAcid otherAcid)
	{
	    return calculateAtomDistance(otherAcid,"CA", "CA");
            
	}
	
	/**
	 * Method returns the CB distance of the two amino acids. For Glycin the CA atom is used for 
	 * distance calculation
	 * 
	 * @param otherAcid
	 * @return
	 */
	public float calculateCBDistance(AminoAcid otherAcid)
	{
                String own_type = "CA";
                String other_type = "CA";
		if(oneLetterCode != 'G')
		    own_type="CB";
		
		if(otherAcid.getOneLetterCode() != 'G')
                    other_type="CB";
                
                return calculateAtomDistance(otherAcid, own_type, other_type);
	}
	
	/**
	 * Method artificially places the CB atom of the amino acid and returns the additional atom to the calling class
	 * 
	 * @return
	 */
	public Atom computeArtificialCBPosition()
	{                                
		float cbX = 0;
		float cbY = 0;
		float cbZ = 0;
		
		float[] v1 = new float[]{0,0,0};
		float[] v2 = new float[]{0,0,0};
		float[] v3 = new float[]{0,0,0};
		
		v1[0] = getXCACoordinate() - getAtom("N").getXCoordinate();
		v1[1] = getYCACoordinate() - getAtom("N").getYCoordinate();
		v1[2] = getZCACoordinate() - getAtom("N").getZCoordinate();

		v2[0] = getXCACoordinate() - getAtom("C").getXCoordinate();
		v2[1] = getYCACoordinate() - getAtom("C").getYCoordinate();
		v2[2] = getZCACoordinate() - getAtom("C").getZCoordinate();

		v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
		v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
		v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

		float gamma = 0.56f;

		cbX = getXCACoordinate() + gamma * (v1[0] +v2[0] + v3[0]);
		cbY = getYCACoordinate() + gamma * (v1[1] +v2[1] + v3[1]);
		cbZ = getZCACoordinate() + gamma * (v1[2] +v2[2] + v3[2]);
		
		return new Atom("CB", cbX, cbY, cbZ);
	}

        public class ClosestContact
        {
            public String atomtype1;
            public String atomtype2;
            public float distance;
        }
        
        public float calculateMinimalDistance(AminoAcid otherAcid)
        {
            ClosestContact csc = calculateClosestContact(otherAcid);
            return csc.distance;
        }

	/**
	 * This method calculates the minimal distance between this amino acid and
	 * the given amino acid. This means that the distances of all atoms of this
	 * amino acid with the other amino acids atoms are calculated and the minimal
	 * distance is returned.
	 * 
	 * @param otherAcid
	 * @return float
	 */
	public ClosestContact calculateClosestContact(AminoAcid otherAcid)
	{
		// calculate the CA distance first
                ClosestContact csc = new ClosestContact();
		csc.distance = Float.MAX_VALUE;
                

		// calculate all other distances between the atoms of this amino acid and 
		// all the other atoms of the other amino acid
                
		Iterator<String> otherAtomsIterator = otherAcid.getAtomNames();
                
                // iterate over all other atoms
                while (otherAtomsIterator.hasNext())
		{
                    Atom otherAtom = otherAcid.getAtom(otherAtomsIterator.next());
                    
                    Iterator<String> thisAtomsIterator = getAtomNames();
					
                    while (thisAtomsIterator.hasNext())
                    {
                            Atom thisAtom = getAtom(thisAtomsIterator.next());

                            float newDistance = thisAtom.calculateDistance(otherAtom);

                            if (newDistance < csc.distance)
                            {
                                    csc.distance = newDistance;
                                    csc.atomtype1=thisAtom.getAtomName();
                                    csc.atomtype2=otherAtom.getAtomName();
                            }
                    }
                }
		return csc;
	}
	
	/**
	 * Method creates a clone of this instance
	 */
	public AminoAcid clone()
	{
		// initialize clone
		AminoAcid clone = new AminoAcid(this.oneLetterCode, this.secondaryStructureElement, this.residueNumber);
		clone.setResiduePositionInSequence(this.residuePositionInArray);
		
		clone.setChainID(this.chainID);
		
		for(String a : atoms.keySet())
		{
			clone.addAtom(atoms.get(a).clone());
		}
		
		clone.setPatternOccurrenceInfo(this.patternOccurrenceInfo);
		
		if(is_heteroatom)
			setHeteroAtomName(this.residueName, this.baseOneLetterCode+"");
		
		return clone;
	}

	/**
	 * Method returns pattern reference information of the amino acid
	 * 
	 * @return
	 */
	public String[] getPatternOccurrenceInfo()
	{
		return patternOccurrenceInfo;
	}

	/**
	 * Method sets pattern reference information of the amino acid
	 * 
	 * @param patternReferenceInfo
	 */
	public void setPatternOccurrenceInfo(String[] patternReferenceInfo)
	{
		this.patternOccurrenceInfo = patternReferenceInfo;
	}	
}