package Manananggal;

public class NumberFilter
{
	private float m_fMinValue;
	private float m_fMaxValue;
	
	public NumberFilter()
	{
		m_fMinValue = Float.NEGATIVE_INFINITY;
		m_fMaxValue	= Float.POSITIVE_INFINITY;
	}
	
	public void SetMinValue(float fValue)
	{
		m_fMinValue = fValue;
	}
	
	public void SetMaxValue(float fValue)
	{
		m_fMaxValue = fValue;
	}
	
	public double GetMinValue()
	{
		return m_fMinValue;
	}
	
	public double GetMaxValue()
	{
		return m_fMaxValue;
	}
	
	public void Disable()
	{
		m_fMinValue = Float.NEGATIVE_INFINITY;
		m_fMaxValue	= Float.POSITIVE_INFINITY;
	}
	
	public boolean IsInRange(double fVal)
	{		
		if(fVal >= m_fMinValue && fVal <= m_fMaxValue)
			return true;

		return false;
	}
	
	public boolean IsInRange(int nVal)
	{
		if(nVal >= m_fMinValue && nVal <= m_fMaxValue)
			return true;
		
		return false;
	}
}
