/**
 * This class holds the results from the label quantification algorithm
 * Author: Mr. Srikanth Nagella
 * Date: 18.06.2014
 */

#ifndef CCPILABELQUANTIFICATIONRESULT_H
#define CCPILABELQUANTIFICATIONRESULT_H

#include <string>
#include <vector>
#include <map>
#include <list>

class CCPiLabelQuantificationResult
{
public:
	CCPiLabelQuantificationResult();
	~CCPiLabelQuantificationResult();
	std::vector<std::string> GetQuantityNames();
	void SetValue(std::string quantityName, int labelIndex, double value);
	double GetValue(std::string quantityName, int labelIndex);
	std::list<int> GetLabelIndexes();
private:
	void initialise();
	std::map< std::string,std::map<int, double> > QuantityMap;
	std::list<int> rowIndexList;
};

#endif
