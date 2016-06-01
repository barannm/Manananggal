package BioKit;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * This class stores compressed sequences in the main memory. Additionally to encoding sequences
 * as byte arrays, the class encodes byte arrays using GZIP writer and reader classes.
 * 
 * @author Florian Oefinger
 */
public abstract class CompressedSequence
{
	// storeage array to store compressed sequence as byte array
	protected byte[] storage;

	/**
	 * Method must be provided by implementing classes to allow compression
	 * of a specified sequence to a byte array.
	 * 
	 * @param sequence
	 * @return
	 */
	protected abstract byte[] encode(String sequence);

	/**
	 * Method must be provided by implementing classes to allow
	 * decompression of a specified byte array to a sequence representation
	 * 
	 * @param b
	 * @return
	 */
	protected abstract String decode(byte[] b);

	/**
	 * Constructor, used to zip a sequence using GZIP writer
	 * 
	 * @param sequence
	 * @throws IOException
	 */
	public CompressedSequence(String sequence) throws IOException
	{
		ByteArrayOutputStream stream = new ByteArrayOutputStream();
		GZIPOutputStream out = new GZIPOutputStream(stream);
		
		try
		{
			out.write(encode(sequence));
		}
		catch (NullPointerException e)
		{
			System.out.println(sequence);
		}
		out.close();
		
		storage = stream.toByteArray();
	}

	/**
	 * Method returns a string representation of the sequence zipped in this instance.
	 */
	public String toString()
	{
		try
		{
			GZIPInputStream in = new GZIPInputStream(new ByteArrayInputStream(storage));
			byte[] b = new byte[1];
			StringBuffer buff = new StringBuffer();
			while (in.read(b) != -1)
			{
				buff.append(decode(b));
			}
			return buff.toString();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		return null;
	}
}