/**
 * This header file is a data class for the nexus reader. This allows reading of data from command line and via the GUI interface
 * @author: Mr. Srikanth Nagella
 */

#ifndef CCPIAVIZONEXUSCMDREADER_H
#define CCPIAVIZONEXUSCMDREADER_H

#include <hxcore/HxData.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortFilename.h>
#include <hxcore/HxPortMultiMenu.h>
#include <hxcore/HxPortIntTextN.h>
#include "api.h"


class CCPI_API CCPiAvizoNexusCmdReader : public HxData
{
	HX_HEADER(CCPiAvizoNexusCmdReader);

public:
//	CCPiAvizoNexusCmdReader();
//	~CCPiAvizoNexusCmdReader();

	/** Port providing a button to click to run the module */
	HxPortDoIt portAction;
	/** Path to the nexus file*/
	HxPortFilename portFilename;
	/** Path for the dataset*/
	HxPortMultiMenu portDatasetPath;
	/** Dimension one*/
	HxPortIntTextN portDimensionOne;
	/** Dimension two*/
	HxPortIntTextN portDimensionTwo;
	/** Dimension three*/
	HxPortIntTextN portDimensionThree;
	/** Perform the calculation in this module. Called by Avizo. */
	virtual void compute();

private:
	void resetDimensionWidgets();
	std::vector<long> getStartValues();
	std::vector<long> getStrideValues();
	std::vector<long> getCountValues();
};
#endif
