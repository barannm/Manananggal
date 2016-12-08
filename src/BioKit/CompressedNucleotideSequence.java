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

import java.io.IOException;
import java.util.HashMap;

/**
 * Class is used to compress nucleotide sequences
 * @author Florian Oefinger
 *
 */
public class CompressedNucleotideSequence extends CompressedSequence
{
	private static HashMap<String, Byte> encoding;
	private static HashMap<Byte, String> decoding;
	private static String[] alphabet;
	
	static
	{
		alphabet = new String[] { "A", "C", "G", "T", "N" };
		encoding = new HashMap<String, Byte>();
		decoding = new HashMap<Byte, String>();
		
		int index = 0;
		
		for (int a = 0; a < alphabet.length; a++)
		{
			for (int b = 0; b < alphabet.length; b++)
			{
				for (int c = 0; c < alphabet.length; c++)
				{
					encoding.put(alphabet[a] + alphabet[b] + alphabet[c], new Byte((byte) index));
					decoding.put(new Byte((byte) index), alphabet[a] + alphabet[b] + alphabet[c]);
					index++;
				}
			}
		}
		
		for (int a = 0; a < alphabet.length; a++)
		{
			for (int b = 0; b < alphabet.length; b++)
			{
				encoding.put(alphabet[a] + alphabet[b], new Byte((byte) index));
				decoding.put(new Byte((byte) index), alphabet[a] + alphabet[b]);
				index++;
			}
		}
		
		for (int a = 0; a < alphabet.length; a++)
		{
			encoding.put(alphabet[a], new Byte((byte) index));
			decoding.put(new Byte((byte) index), alphabet[a]);
			index++;
		}

	}

	/**
	 * Constructor, initializes a compressed sequence with a specified string
	 * 
	 * @param sequence
	 * @throws IOException
	 */
	public CompressedNucleotideSequence(String sequence) throws IOException
	{
		super(sequence.toUpperCase().replaceAll("[^ACGTN]", "N"));
	}

	/**
	 * Method encodes a nucleotide sequence and returns the encoded byte array
	 */
	protected byte[] encode(String sequence)
	{
		byte[] res;
		if (sequence.length() % 3 == 0)
		{
			res = new byte[sequence.length() / 3];
		}
		else
		{
			res = new byte[sequence.length() / 3 + 1];
		}
		for (int i = 0; i < sequence.length() / 3; i++)
		{
			res[i] = encoding.get(sequence.substring(i * 3, i * 3 + 3).toUpperCase());
		}
		if (sequence.length() % 3 == 2)
		{
			res[res.length - 1] = encoding.get(sequence.substring(sequence.length() - 2).toUpperCase());
		}
		else if (sequence.length() % 3 == 1)
		{
			res[res.length - 1] = encoding.get(sequence.substring(sequence.length() - 1).toUpperCase());
		}
		return res;
	}

	/**
	 * Method decodes the nucleotide sequence represented in the byte array and returns the string
	 * representation to the calling method.
	 * 
	 */
	protected String decode(byte[] enc)
	{
		StringBuffer buff = new StringBuffer();
		for (byte b : enc)
		{
			buff.append(decoding.get(b));
		}
		return buff.toString();
	}
}
