#include "CCPiLabelQuantificationResult.h"

CCPiLabelQuantificationResult::CCPiLabelQuantificationResult()
{
	initialise();
}

CCPiLabelQuantificationResult::~CCPiLabelQuantificationResult()
{
}

void CCPiLabelQuantificationResult::initialise()
{
}

std::vector<std::string> CCPiLabelQuantificationResult::GetQuantityNames()
{
	std::vector<std::string> result;
	for(std::map<std::string,std::map<int,double>>::iterator itr = QuantityMap.begin(); itr!=QuantityMap.end(); itr++)
	{
		result.push_back(itr->first);
	}
	return result;
}

void CCPiLabelQuantificationResult::SetValue(std::string quantityName, int labelIndex, double value)
{
	(QuantityMap[quantityName])[labelIndex] = value;
	rowIndexList.push_back(labelIndex);
	rowIndexList.unique();
}

double CCPiLabelQuantificationResult::GetValue(std::string quantityName, int labelIndex)
{
	return (QuantityMap[quantityName])[labelIndex];
}

std::list<int> CCPiLabelQuantificationResult::GetLabelIndexes()
{
	return rowIndexList;
}

