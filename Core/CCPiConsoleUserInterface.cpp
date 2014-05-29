#include "CCPiConsoleUserInterface.h"
#include <iostream>

void CCPiConsoleUserInterface::LogMessage(std::string message)
{
	std::cout << message <<std::endl;
}

void CCPiConsoleUserInterface::SetStatusMessage(std::string message)
{
	std::cout<< "Status: "<<message<<std::endl;
}

void CCPiConsoleUserInterface::SetProgressValue(float value)
{
	std::cout<< "Progress: "<<value<<std::endl;
}

bool CCPiConsoleUserInterface::isCancel()
{
	return false;
}
