package BioKit;

import java.io.Serializable;

/**
 * This class represents a vector in 3D space. The class provides a set of vector operations.
 * Vectors are compared (compareTo - method) and identified based on their id set in the 
 * constructor. 
 * 
 * @author Fabian Birzele
 */
public class Vector3D implements Comparable<Vector3D>, Serializable, Cloneable
{
	/**
	 * Serial Version UID, created on 21.04.2006 with the following instance variables:
	 * 
	 * private float xCoordinate;
	 * private float yCoordinate;
	 * private float zCoordinate;
	 * private int id;
	 */
	private static final long serialVersionUID = -7892177393515345337L;

	public static final double RADIANS_TO_DEGREES_CONSTANT = 180.0 / Math.PI;
	
	private float xCoordinate;
	private float yCoordinate;
	private float zCoordinate;
	
	private int id;
	
	/**
	 * Constructor, initializes the vertex
	 * 
	 * @param id
	 * @param xCoordinate
	 * @param yCoordinate
	 * @param zCoordinate
	 */
	public Vector3D(int id, float xCoordinate, float yCoordinate, float zCoordinate)
	{
		this.id = id;
		
		this.xCoordinate = xCoordinate;
		this.yCoordinate = yCoordinate;
		this.zCoordinate = zCoordinate;
	}
	
	/**
	 * Method compares two vertices based on their IDs
	 * 
	 * @param otherVertex
	 * @return
	 */
	public int compareTo(Vector3D otherVertex)
	{
		int otherID = ((Vector3D)otherVertex).getId();
		
		if(id < otherID)
			return -1;
		else if(id == otherID)
			return 0;
		
		return 1;
	}

	/**
	 * Method returns the ID of the vertex
	 * 
	 * @return
	 */
	public int getId()
	{
		return id;
	}
	
	/**
	 * Method sets the x coordinate of the vertex
	 */
	public void setXCoordinate(float xCoordinate)
	{
		this.xCoordinate = xCoordinate;
	}

	/**
	 * Method sets the y coordinate of the vertex
	 * 
	 * @return
	 */
	public void setYCoordinate(float yCoordinate)
	{
		this.yCoordinate = yCoordinate;
	}

	/**
	 * Method sets the z coordinate of the vertex
	 * 
	 * @return
	 */
	public void setZCoordinate(float zCoordinate)
	{
		this.zCoordinate = zCoordinate;
	}

	/**
	 * Method returns the x coordinate of the vertex
	 * 
	 * @return
	 */
	public float getXCoordinate()
	{
		return xCoordinate;
	}

	/**
	 * Method returns the y coordinate of the vertex
	 * 
	 * @return
	 */
	public float getYCoordinate()
	{
		return yCoordinate;
	}

	/**
	 * Method returns the z coordinate of the vertex
	 * 
	 * @return
	 */
	public float getZCoordinate()
	{
		return zCoordinate;
	}
	
	public String toString()
	{
		return "ID: " + id + ", coordinates: (" + xCoordinate + ", " + yCoordinate + ", " + zCoordinate + "), length: " + length(); 
	}
	
	/**
	 * Method calculates the length of the vertex
	 * 
	 * @return
	 */
	public float length()
	{
		return (float)Math.sqrt(xCoordinate*xCoordinate + yCoordinate*yCoordinate + zCoordinate*zCoordinate);
	}
	
	public void normalize()
	{
		float length = length();
		
		xCoordinate /= length;
		yCoordinate /= length;
		zCoordinate /= length;
	}
	
	/**
	 * Method returns a cloned instance of this vector
	 */
	public Object clone()
	{
		return new Vector3D(id, xCoordinate, yCoordinate, zCoordinate);
	}

	/**
	 * Subtract the specified vector from this vector. This vector is changed with
	 * that operation!
	 * 
	 * @param second
	 */
	public void substractFrom(Vector3D second)
	{
		xCoordinate -= second.getXCoordinate();
		yCoordinate -= second.getYCoordinate();
		zCoordinate -= second.getZCoordinate();
	}
	
	/**
	 * Multiplies Vector with a scalar.
	 */
	public void multiply(double scalar)
	{
		xCoordinate *= scalar;
		yCoordinate *= scalar;
		zCoordinate *= scalar;
	}
	
	/**
	 * Method rotates the vector using the specified rotation matrix
	 * 
	 * @param rotationMatrix
	 */
	public void rotateVector(double [][] rotationMatrix)
	{
		float xCoordinateTemp = (float)(rotationMatrix[0][0] * xCoordinate + rotationMatrix[0][1] * yCoordinate + rotationMatrix[0][2] * zCoordinate);
		float yCoordinateTemp = (float)(rotationMatrix[1][0] * xCoordinate + rotationMatrix[1][1] * yCoordinate + rotationMatrix[1][2] * zCoordinate);
		float zCoordinateTemp = (float)(rotationMatrix[2][0] * xCoordinate + rotationMatrix[2][1] * yCoordinate + rotationMatrix[2][2] * zCoordinate);
		
		xCoordinate = xCoordinateTemp;
		yCoordinate = yCoordinateTemp;
		zCoordinate = zCoordinateTemp;
	}
	
	/**
	 * Substracts vector from this vector.
	 * This vector is changed during translation!
	 *
	 * @param traslationVector
	 */
	public void substract(Vector3D traslationVector){
		xCoordinate -= traslationVector.getXCoordinate();
		yCoordinate -= traslationVector.getYCoordinate();
		zCoordinate -= traslationVector.getZCoordinate();
	}

	
	/**
	 * Method returns a new Vector3D moved along a second vector
	 * This vector is changed during translation!
	 * 
	 * @param traslationVector
	 * 
	 * @return new Vector3D - moved along the translation vector
	 */
	public void add(Vector3D traslationVector){
		xCoordinate += traslationVector.getXCoordinate();
		yCoordinate += traslationVector.getYCoordinate();
		zCoordinate += traslationVector.getZCoordinate();
	}
	
	/**
	 * Method rotates the vector using the specified rotation matrix
	 * 
	 * @return Array with x,y,z-coordinate
	 */
	public double[] getCoordinates(){
		double[] returnArray = {getXCoordinate(), getYCoordinate(), getZCoordinate()};
		return returnArray;
	}
	
	/**
	 * Unary minus. This Vector is CHANGED during the operation!
	 * 
	 */
	public void unMinus(){
		xCoordinate = -xCoordinate;
		yCoordinate = -yCoordinate;
		zCoordinate = -zCoordinate;
	}
	
	/**
	 * Print Method
	 * 
	 */
	public void print(){
		System.out.println("x: " + xCoordinate + " y: " + yCoordinate + " z: " + zCoordinate);
	}
	
	/**
	 * Method subtracts the second argument from the first argument and returns a new 
	 * vertex.
	 * 
	 * @param one
	 * @param two
	 * @return
	 */
	public static Vector3D substract(Vector3D one, Vector3D two)
	{
		return new Vector3D(0, one.getXCoordinate() - two.getXCoordinate(), one.getYCoordinate() - two.getYCoordinate(), one.getZCoordinate() - two.getZCoordinate());
	}
	
	/**
	 * Method adds the second argument to the first argument and returns a new 
	 * vertex.
	 * 
	 * @param one
	 * @param two
	 * @return
	 */
	public static Vector3D add(Vector3D one, Vector3D two)
	{
		return new Vector3D(0, one.getXCoordinate() + two.getXCoordinate(), one.getYCoordinate() + two.getYCoordinate(), one.getZCoordinate() + two.getZCoordinate());
	}
	
	/**
	 * Method calculates the cross product of the two vertices (vectors)
	 * 
	 * @param one
	 * @param two
	 * @return
	 */
	public static Vector3D crossProduct(Vector3D one, Vector3D two)
	{
		float x = one.getYCoordinate()*two.getZCoordinate() - one.getZCoordinate()*two.getYCoordinate();
		float y = one.getZCoordinate()*two.getXCoordinate() - one.getXCoordinate()*two.getZCoordinate();
		float z = one.getXCoordinate()*two.getYCoordinate() - one.getYCoordinate()*two.getXCoordinate();
		
		return new Vector3D(0, x, y, z);
	}
	
	/**
	 * Method calculates the distance of the two vectors and returns it
	 * 
	 * @param one
	 * @param two
	 * @return
	 */
	public static double distance(Vector3D one, Vector3D two)
	{
		float x = one.getXCoordinate() - two.getXCoordinate();
		float y = one.getYCoordinate() - two.getYCoordinate();
		float z = one.getZCoordinate() - two.getZCoordinate();
		
		return Math.sqrt(x*x + y*y + z*z);
	}
	
	/**
	 * Calculates the dot product of the two vertices
	 * 
	 * @param one
	 * @param two
	 * @return
	 */
	public static double dotProduct(Vector3D one, Vector3D two)
	{
		return one.getXCoordinate()*two.getXCoordinate() + one.getYCoordinate()*two.getYCoordinate() + one.getZCoordinate()*two.getZCoordinate();
	}
	
	/**
	 * Calculates the vector resulting from a scalar multiplication
	 * 
	 * @param scale
	 * @param vertex
	 * @return
	 */
	public static Vector3D scalarMultiplication(double scalar, Vector3D one){
		return new Vector3D(0,one.getXCoordinate()*(float)scalar,one.getYCoordinate()*(float)scalar,one.getZCoordinate()*(float)scalar);
	}
	
	/**
	 * Calculates the angle between two vectors
	 * 
	 * @param one First Vector
	 * @param two Second Vector
	 * @return Returns the angle in radians.
	 */
	public static double angle(Vector3D one, Vector3D two){
		return Math.acos(dotProduct(one, two)/one.length()/two.length());
	}
	
	/**
	 * Method rotates a vector object using the given rotation matrix. A new Vector3D
	 * instance containing the rotated coordinates is returned
	 * 
	 * @param rotationMatrix	rotation matrix to rotate the point
	 * @param point			point to be rotated
	 * @return Vector3D			rotatedPoint
	 */
	public static Vector3D generateNewRotatedVector(double [][] rotationMatrix, Vector3D point)
	{
		double x = rotationMatrix[0][0] * point.getXCoordinate() + rotationMatrix[0][1] * point.getYCoordinate() + rotationMatrix[0][2] * point.getZCoordinate();
		double y = rotationMatrix[1][0] * point.getXCoordinate() + rotationMatrix[1][1] * point.getYCoordinate() + rotationMatrix[1][2] * point.getZCoordinate();
		double z = rotationMatrix[2][0] * point.getXCoordinate() + rotationMatrix[2][1] * point.getYCoordinate() + rotationMatrix[2][2] * point.getZCoordinate();
		
		return new Vector3D(0, (float)x, (float)y, (float)z); 
	}
	
	/**
	 * Method calculates unary minus of the vector and returns a new vector instance containing this
	 * information
	 * 
	 * @param one
	 * @return
	 */
	public static Vector3D unMinus(Vector3D one){
		return new Vector3D(0,-one.xCoordinate, -one.yCoordinate, -one.zCoordinate);
	}
	
	/**
	 * Static method to normalize a vector. A new, normalized vector instance is returned
	 * 
	 * @param one
	 * @return
	 */
	public static Vector3D normalize(Vector3D one){
		return new Vector3D(0, one.xCoordinate/one.length(), one.yCoordinate/one.length(), one.zCoordinate/one.length());
	}
	
	/**
	 * Method converts radians to degrees
	 * 
	 * @param radians
	 * @return
	 */
	public static double radiansToDegrees(double radians)
	{
		return radians * RADIANS_TO_DEGREES_CONSTANT;
	}
	
}
