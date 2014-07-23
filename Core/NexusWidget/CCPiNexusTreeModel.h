/**
 * Tree model for Nexus data
 *
 */
#ifndef CCPINEXUSTREEMODEL_H
#define CCPINEXUSTREEMODEL_H

#include <QAbstractItemModel>
#include <QModelIndex>
#include <QVariant>

#include "hdf5.h"
#include <string>
class CCPiNexusTreeItem;
class CCPiNexusTreeModel;


class CCPiNexusTreeModel : public QAbstractItemModel
{
    Q_OBJECT

public:
    CCPiNexusTreeModel( std::string filename, QObject *parent = 0);
    ~CCPiNexusTreeModel();

    QVariant data(const QModelIndex &index, int role) const;
    Qt::ItemFlags flags(const QModelIndex &index) const;
    QVariant headerData(int section, Qt::Orientation orientation,
                        int role = Qt::DisplayRole) const;
    QModelIndex index(int row, int column,
                      const QModelIndex &parent = QModelIndex()) const;
    QModelIndex parent(const QModelIndex &index) const;
    int rowCount(const QModelIndex &parent = QModelIndex()) const;
    int columnCount(const QModelIndex &parent = QModelIndex()) const;

private:
    void setupModelData(std::string filename,CCPiNexusTreeItem *parent);

    CCPiNexusTreeItem *rootItem;
};


#endif
