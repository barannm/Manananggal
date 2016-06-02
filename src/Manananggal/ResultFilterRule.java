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
 * @author Matthias Barann
 */
package Manananggal;

public class ResultFilterRule
{	
	static final int FILTER_TYPE_CONDITION_A 		= 0;
	static final int FILTER_TYPE_CONDITION_B 		= 1;
	static final int FILTER_TYPE_SIGNIFICANT_PSI	= 2;
	static final int FILTER_TYPE_NOVEL_JUNCTION		= 3;
	private static final int MAX_TEXT_FILTERS		= 4;
	
	static final int FILTER_TYPE_RATIO_CHANGE_A		= 0;
	static final int FILTER_TYPE_RATIO_CHANGE_B		= 1;
	static final int FILTER_TYPE_PSI_CHANGE			= 2;
	static final int FILTER_TYPE_P_VALUE_RATIO_A	= 3;
	static final int FILTER_TYPE_P_VALUE_RATIO_B	= 4;
	static final int FILTER_TYPE_P_VALUE_PSI		= 5;
	private static final int MAX_NUMBER_FILTERS		= 6;

	private boolean 		m_bShowASTypes[];
	private boolean 		m_bShowResultType[];
	private TextFilter[] 	m_pTextFilters;
	private NumberFilter[] 	m_pNumberFilters;
	
	public ResultFilterRule()
	{
		m_bShowASTypes = new boolean[13]; // 0 is not used
		for(int i=0; i<13; i++)
			m_bShowASTypes[i] = true;
		
		m_bShowResultType = new boolean[4]; // 0 is not used
		for(int i=0; i<4; i++)
			m_bShowResultType[i] = true;
		
		m_pTextFilters = new TextFilter[MAX_TEXT_FILTERS];
		m_pNumberFilters = new NumberFilter[MAX_NUMBER_FILTERS];
		
		for(int i=0; i<MAX_NUMBER_FILTERS; i++)
		{
			NumberFilter filter = new NumberFilter(); 
			m_pNumberFilters[i] = filter;
		}
	}
	
	public void SetTextFilter(int nFilterID, TextFilter filter)
	{
		m_pTextFilters[nFilterID] = filter;
	}
	
	public TextFilter GetTextFilter(int nFilterID)
	{
		return m_pTextFilters[nFilterID];
	}
	
	public void SetNumberFilter(int nFilterID, NumberFilter filter)
	{
		m_pNumberFilters[nFilterID] = filter;
	}
	
	public NumberFilter GetNumberFilter(int nFilterID)
	{
		return m_pNumberFilters[nFilterID];
	}
	
	//#############################################
	public void ShowASType(int nSplicingType)
	{
		m_bShowASTypes[nSplicingType] = true;
	}
	
	public void HideASType(int nSplicingType)
	{
		m_bShowASTypes[nSplicingType] = false;
	}
	
	public boolean IsShowingASType(int nSplicingType)
	{
		return m_bShowASTypes[nSplicingType];
	}
	
	//#############################################
	public void ShowResultType(int nResultType)
	{
		m_bShowResultType[nResultType] = true;
	}
	
	public void HideResultType(int nResultType)
	{
		m_bShowResultType[nResultType] = false;
	}
	
	public boolean IsShowingResultType(int nResultType)
	{
		return m_bShowResultType[nResultType];
	}
	
	//#############################################
	public boolean bIsValidResult(AnalysisResult result)
	{
		if(!IsShowingASType(result.GetType()))
			return false;
		
		if(!IsShowingResultType(result.GetResultType()))
			return false;
		
		for(int nFilterID = 0; nFilterID < MAX_TEXT_FILTERS; nFilterID++)
		{
			switch(nFilterID)
			{
				case FILTER_TYPE_CONDITION_A:
				{
					if(!m_pTextFilters[nFilterID].IsStringSelected(result.GetConditionA()))
						return false;
					break;
				}
				
				case FILTER_TYPE_CONDITION_B:
				{
					if(!m_pTextFilters[nFilterID].IsStringSelected(result.GetConditionB()))
						return false;
					break;
				}
				
				case FILTER_TYPE_SIGNIFICANT_PSI:
				{
					String strQuery = "NA";
					if(result.HasSignificantPSIScore())
						strQuery = "true";
					else
					{
						if(result.HasPSIScore())
							strQuery = "false";
					}
						 
					if(!m_pTextFilters[nFilterID].IsStringSelected(strQuery))
						return false;
					break;
				}
				
				case FILTER_TYPE_NOVEL_JUNCTION:
				{
					String strQuery = "NA";
					if(result.UsesNovelJunction())
						strQuery = "true";
					else
					{
						if(result.HasPSIScore())
							strQuery = "false";
					}
						 
					if(!m_pTextFilters[nFilterID].IsStringSelected(strQuery))
						return false;
					
					break;
				}
			}
		}
		
		for(int nFilterID = 0; nFilterID < MAX_NUMBER_FILTERS; nFilterID++)
		{
			double fQueryValue = 0.0;

			switch(nFilterID)
			{
				case FILTER_TYPE_RATIO_CHANGE_A:
				{
					if(result.HasAltExonA())
						fQueryValue = result.GetAbsoluteChangeA();
					break;
				}
				case FILTER_TYPE_RATIO_CHANGE_B:
				{
					if(result.HasAltExonB())
						fQueryValue = result.GetAbsoluteChangeB();
					break;
				}
				case FILTER_TYPE_PSI_CHANGE:
				{
					if(result.HasPSIScore())
						fQueryValue = result.GetInclusionChange();
					break;
				}
				case FILTER_TYPE_P_VALUE_RATIO_A:
				{
					if(result.HasAltExonA())
						fQueryValue = result.GetPValueA();
					break;
				}
				case FILTER_TYPE_P_VALUE_RATIO_B:
				{
					if(result.HasAltExonB())
						fQueryValue = result.GetPValueB();
					break;
				}
				case FILTER_TYPE_P_VALUE_PSI:
				{
					if(result.HasPSIScore())
						fQueryValue = result.GetPValuePSI();
					break;
				}
			}
			
			if(!m_pNumberFilters[nFilterID].IsInRange(fQueryValue))
				return false;
		}
		
		return true;
	}
};