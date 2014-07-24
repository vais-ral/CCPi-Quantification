/**
 * This is the nexus Dialog box
 */

#ifndef CCPINEXUSWIDGETDIALOG_H
#define CCPINEXUSWIDGETDIALOG_H

#include "CCPiDefines.h"
#include "ui_CCPiNexusWidget.h"
#include "CCPiNexusTreeModel.h"
#include <vector>

class CCPI_EXPORT CCPiNexusWidgetDialog: public QDialog
{
	Q_OBJECT
public:
	explicit CCPiNexusWidgetDialog(std::string filename, QWidget *parent=0);
	~CCPiNexusWidgetDialog();
	void SetModel(QAbstractItemModel *);
	std::vector<std::string> GetSelectedDataSetList(){return SelectedDatasetList;}
private:
	void accept();
	std::vector<std::string> SelectedDatasetList;
	Ui::CCPiNexusWidgetDialog *ui;
    CCPiNexusTreeModel *model;
	std::string filename;
};

#endif
