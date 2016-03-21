#ifndef CCPINEXUSSUBSETWIDGET_H
#define CCPINEXUSSUBSETWIDGET_H

#include "CCPiDefines.h"
#include <QDialog>
#include <QGridLayout>
#include <QLineEdit>
#include <QPushButton>
#include <vector>

class CCPI_EXPORT CCPiNexusSubsetWidget : public QDialog
{
	Q_OBJECT

public:
	CCPiNexusSubsetWidget(QWidget *parent, std::vector<std::string> names, std::vector<long> dimensionValues);
	std::vector<long> getSelectedStartValues(); //returns the selected start values
	std::vector<long> getSelectedCountValues(); // returns the selected count values
	std::vector<long> getSelectedStrideValues(); //returns the selected stride values

private slots:
	void okButtonClicked(); //Validate the given inputs when ok button is pressed. if any errors display dialog
	void cancelButtonClicked(); //cancel the dialog and set the values of the start,end and stride to defaults.(entire range)

private:
	std::vector<long> selectedStartValues; // selected start values
	std::vector<long> selectedCountValues; // selected count values
	std::vector<long> selectedStrideValues; // selected Stride values
	std::vector<long> rangeValues; // Original dimensions
	std::vector<std::string> dimensionNames;
	QPushButton *okButton;
	QPushButton *cancelButton;
	QGridLayout *layout;
	std::vector<QLineEdit*> *startWidgets;
	std::vector<QLineEdit*> *countWidgets;
	std::vector<QLineEdit*> *strideWidgets;
};
#endif