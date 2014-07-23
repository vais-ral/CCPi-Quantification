#include "CCPiNexusWidgetDialog.h"
#include "CCPiNexusTreeItem.h"
#include <iostream>

CCPiNexusWidgetDialog::CCPiNexusWidgetDialog(std::string filename, QWidget *parent) : QDialog(parent), ui(new Ui::CCPiNexusWidgetDialog)
{
	ui->setupUi(this);
	this->filename = filename;
	model = new CCPiNexusTreeModel(filename);
	ui->treeView->setModel(model);
}

CCPiNexusWidgetDialog::~CCPiNexusWidgetDialog()
{
	delete ui;
	delete model;
}

void CCPiNexusWidgetDialog::SetModel(QAbstractItemModel *model)
{
	ui->treeView->setModel(model);
}

void CCPiNexusWidgetDialog::accept()
{
	QModelIndexList selectedRows = ui->treeView->selectionModel()->selectedRows();
	foreach(QModelIndex index, selectedRows)
	{
		CCPiNexusTreeItem *item = static_cast<CCPiNexusTreeItem*>(index.internalPointer());
		//std::cout<<item->data(0).toString().toUtf8().constData()<<std::endl;
		SelectedDatasetList.push_back(item->data(0).toString().toUtf8().constData());
	}
	QDialog::accept();
}