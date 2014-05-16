/**
 * This is an interface between CCPi core algorithms and External application specific calls.
 * such as Avizo, etc.
 * Author: Mr. Srikanth Nagella
 * Date  : 14.05.2014
 */

#ifndef CCPIUSERAPPLICATIONINTERFACE_H
#define CCPIUSERAPPLICATIONINTERFACE_H

#include <string>

class CCPiUserApplicationInterface
{
public:
	virtual void LogMessage(std::string message)=0;
	virtual void SetStatusMessage(std::string message)=0;
	virtual void SetProgressValue(float value)=0;
	virtual bool isCancel()=0;
};

#endif
