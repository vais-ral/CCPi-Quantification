#include "CCPiNexusSubsetWidget.h"
#include <QGridLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <iostream>

CCPiNexusSubsetWidget::CCPiNexusSubsetWidget(QWidget *parent, std::vector<std::string> names, std::vector<long> rangeValues)
						: QDialog(parent), dimensionNames(names), rangeValues(rangeValues)
{
	layout = new QGridLayout();
	if(dimensionNames.size()==0){
		for(int i=0;i<rangeValues.size();i++)
		{
			QString dimname = QString("dim")+QString::number(i);
			QByteArray ba = dimname.toUtf8();
			std::string strdimname(ba.data());
			dimensionNames.push_back(strdimname);
		}
	}
	//Initalize startWidgets
	startWidgets = new std::vector<QLineEdit *>();
	countWidgets = new std::vector<QLineEdit *>();
	strideWidgets = new std::vector<QLineEdit *>();

	//Add the header
	layout->addWidget(new QLabel("Dimension Name"),0,0);
	layout->addWidget(new QLabel("Start"),0,1);
	layout->addWidget(new QLabel("Count"),0,2);
	layout->addWidget(new QLabel("Stride"),0,3);

	for(int i=0;i<dimensionNames.size();i++)
	{
		layout->addWidget(new QLabel(QString(dimensionNames[i].c_str())), i+1, 0);
		QLineEdit *startWidget = new QLineEdit();
		startWidgets->push_back(startWidget);
		QLineEdit *endWidget = new QLineEdit();
		countWidgets->push_back(endWidget);
		QLineEdit *strideWidget = new QLineEdit();
		strideWidgets->push_back(strideWidget);

		//set values
		startWidget->setText("0");
		endWidget->setText(QString::number(rangeValues[i]));
		strideWidget->setText("1");

		//add to layout
		layout->addWidget(startWidget, i+1, 1);
		layout->addWidget(endWidget, i+1, 2);
		layout->addWidget(strideWidget, i+1, 3);
	}
	//Add buttons
	okButton = new QPushButton("OK");
	cancelButton = new QPushButton("Cancel");
	layout->addWidget(okButton, dimensionNames.size()+2, 2);
	layout->addWidget(cancelButton, dimensionNames.size()+2, 3);

	//Connect the signals and slots
	QObject::connect(okButton,SIGNAL(clicked()),  this, SLOT(okButtonClicked()));
	QObject::connect(cancelButton,SIGNAL(clicked()),  this, SLOT(cancelButtonClicked()));

	//add the layout to the dialog
	this->setLayout(layout);

	//set the default values
	for(int i=0;i<dimensionNames.size();i++)
	{
		selectedStartValues.push_back(0L);
		selectedCountValues.push_back(rangeValues[i]);
		selectedStrideValues.push_back(1L);
	}
}

std::vector<long> CCPiNexusSubsetWidget::getSelectedStartValues()
{
	return selectedStartValues;
}

std::vector<long> CCPiNexusSubsetWidget::getSelectedCountValues()
{
	return selectedCountValues;
}

std::vector<long> CCPiNexusSubsetWidget::getSelectedStrideValues()
{
	return selectedStrideValues;
}

void CCPiNexusSubsetWidget::okButtonClicked()
{
	bool success = true;
	//Check if the start value is > end value
	for(int i=0;i<selectedStartValues.size();i++)
	{
		if((*startWidgets)[i]->text().toLong()>rangeValues[i]) success = false; //start values should always be less than total count
		long totalElements = (*startWidgets)[i]->text().toLong()+((*countWidgets)[i]->text().toLong()*(*strideWidgets)[i]->text().toLong());
		if(totalElements>rangeValues[i]) success = false;
	}
	if(!success)
	{ //Display message box
		QMessageBox::critical(this,"Subset Selection","Input range is incorrect. Please check the values");
		return;
	}
	//everything is ok. now copy the values
	for(int i=0;i<selectedStartValues.size();i++)
	{
		selectedStartValues[i]=(*startWidgets)[i]->text().toLong();
		selectedCountValues[i]=(*countWidgets)[i]->text().toLong();
		selectedStrideValues[i]=(*strideWidgets)[i]->text().toLong();
	}
	close();
}

void CCPiNexusSubsetWidget::cancelButtonClicked()
{
	close();
}