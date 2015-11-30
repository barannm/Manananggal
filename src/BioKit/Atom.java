
package BioKit;

import java.io.Serializable;


/**
 * This class represents a simple atom instance that holds a name,
 * as well as three coordinates (x, y and z). Atom instances can
 * be used to hold information about amino acid atoms.
 * 
 * @author Fabian Birzele
 */
public class Atom implements Serializable
{
	/**
	 * Serial Version UID created on 21.04.2006 with the following instance variables:
	 * 
	 * private String atomName;
	 * private float xCoordinate;
	 * private float yCoordinate;
	 * private float zCoordinate;
	 */
	private static final long serialVersionUID = 8398495386396709951L;
	
	private String atomName;
	private float xCoordinate;
	private float yCoordinate;
	private float zCoordinate;
	
	public Atom(String atomName, float xCoordinate, float yCoordinate, float zCoordinate)
	{
		this.atomName = atomName;
		this.xCoordinate = xCoordinate;
		this.yCoordinate = yCoordinate;
		this.zCoordinate = zCoordinate;
	}
	
	
	/**
	 * Returns the name of the atom
	 * 
	 * @return String
	 */
	public String getAtomName()
	{
		return atomName;
	}

	/**
	 * Returns the x coordinate of the atom
	 * 
	 * @return float
	 */
	public float getXCoordinate()
	{
		return xCoordinate;
	}

	/**
	 * Returns the y coordinate of the atom
	 * 
	 * @return float
	 */
	public float getYCoordinate()
	{
		return yCoordinate;
	}

	/**
	 * Returns the z coordinate of the atom
	 * 
	 * @return float
	 */
	public float getZCoordinate()
	{
		return zCoordinate;
	}
	
	/**
	 * Returns the coordinates wrapped in a Vector3D instance
	 */
	public Vector3D getCoordinates()
	{
		return new Vector3D(0,xCoordinate,yCoordinate,zCoordinate);
	}

	/**
	 * Sets the atomName.
	 * 
	 * @param atomName The atomName to set
	 */
	public void setAtomName(String atomName)
	{
		this.atomName = atomName;
	}

	/**
	 * Sets the xCoordinate.
	 * 
	 * @param xCoordinate The xCoordinate to set
	 */
	public void setXCoordinate(float xCoordinate)
	{
		this.xCoordinate = xCoordinate;
	}

	/**
	 * Sets the yCoordinate.
	 * 
	 * @param yCoordinate The yCoordinate to set
	 */
	public void setYCoordinate(float yCoordinate)
	{
		this.yCoordinate = yCoordinate;
	}

	/**
	 * Sets the zCoordinate.
	 * 
	 * @param zCoordinate The zCoordinate to set
	 */
	public void setZCoordinate(float zCoordinate)
	{
		this.zCoordinate = zCoordinate;
	}
	
	/**
	 * Method calculates the distance of this atom to the other atom specified
	 * 
	 * @param otherAtom
	 * @return float
	 */
   	public float calculateDistance(Atom otherAtom)
   	{
	   	float distanceX = xCoordinate - otherAtom.getXCoordinate();
	   	float distanceY = yCoordinate - otherAtom.getYCoordinate();
	   	float distanceZ = zCoordinate - otherAtom.getZCoordinate();
		
	   	return (float)Math.sqrt(distanceX*distanceX + distanceY*distanceY + distanceZ*distanceZ);
   	}
   	
   	public String toString()
   	{
   		return atomName + "\t" + xCoordinate + "\t" + yCoordinate + "\t" + zCoordinate;
   	}
   	
   	/**
   	 * Method returns a clone of this instance
   	 */
   	public Atom clone()
   	{
   		Atom clone = new Atom(this.atomName, this.xCoordinate, this.yCoordinate, this.zCoordinate);
   		
   		return clone;
   	}
}
