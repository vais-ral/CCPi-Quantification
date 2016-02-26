/**
 * This class is model for the Nexus data set
 *
 *
 */

#include <QtGui>
#include <qstring.h>

#include "CCPiNexusTreeItem.h"
#include "CCPiNexusTreeModel.h"

#include "hdf5.h"
#include <qmutex.h>
#include <string>
#include <sstream> 
#include <map>

QMutex HDFDataMutex;
std::map<std::string,std::string> HDFTreeDataMap;

CCPiNexusTreeModel::CCPiNexusTreeModel(std::string filename, QObject *parent)
    : QAbstractItemModel(parent)
{
    QList<QVariant> rootData;
    rootData << "Data" << "Dimensions";
    rootItem = new CCPiNexusTreeItem(rootData);
    setupModelData(filename, rootItem);
}

CCPiNexusTreeModel::~CCPiNexusTreeModel()
{
    delete rootItem;
}

int CCPiNexusTreeModel::columnCount(const QModelIndex &parent) const
{
    if (parent.isValid())
        return static_cast<CCPiNexusTreeItem*>(parent.internalPointer())->columnCount();
    else
        return rootItem->columnCount();
}

QVariant CCPiNexusTreeModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();

    if (role != Qt::DisplayRole)
        return QVariant();

    CCPiNexusTreeItem *item = static_cast<CCPiNexusTreeItem*>(index.internalPointer());

    return item->data(index.column());
}

Qt::ItemFlags CCPiNexusTreeModel::flags(const QModelIndex &index) const
{
    if (!index.isValid())
        return 0;

    return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}

QVariant CCPiNexusTreeModel::headerData(int section, Qt::Orientation orientation,
                               int role) const
{
    if (orientation == Qt::Horizontal && role == Qt::DisplayRole)
        return rootItem->data(section);

    return QVariant();
}

QModelIndex CCPiNexusTreeModel::index(int row, int column, const QModelIndex &parent)
            const
{
    if (!hasIndex(row, column, parent))
        return QModelIndex();

    CCPiNexusTreeItem *parentItem;

    if (!parent.isValid())
        parentItem = rootItem;
    else
        parentItem = static_cast<CCPiNexusTreeItem*>(parent.internalPointer());

    CCPiNexusTreeItem *childItem = parentItem->child(row);
    if (childItem)
        return createIndex(row, column, childItem);
    else
        return QModelIndex();
}


QModelIndex CCPiNexusTreeModel::parent(const QModelIndex &index) const
{
    if (!index.isValid())
        return QModelIndex();

    CCPiNexusTreeItem *childItem = static_cast<CCPiNexusTreeItem*>(index.internalPointer());
    CCPiNexusTreeItem *parentItem = childItem->parent();

    if (parentItem == rootItem)
        return QModelIndex();

    return createIndex(parentItem->row(), 0, parentItem);
}

int CCPiNexusTreeModel::rowCount(const QModelIndex &parent) const
{
    CCPiNexusTreeItem *parentItem;
    if (parent.column() > 0)
        return 0;

    if (!parent.isValid())
        parentItem = rootItem;
    else
        parentItem = static_cast<CCPiNexusTreeItem*>(parent.internalPointer());

    return parentItem->childCount();
}

herr_t op_func_o(hid_t loc_id, const char *name, const H5O_info_t *info, void *operator_data)
{
    if (name[0] == '.')         /* Root group, do not print '.' */
        printf ("  (Group)\n");
    else
        switch (info->type) {
            case H5O_TYPE_GROUP:
//				HDFTreeDataMap.insert(std::pair<std::string,std::string>(name,name));
                printf ("%s  (Group)\n", name);
                break;
            case H5O_TYPE_DATASET:
				printf("%s (Dataset)\n", name);
				{
				hid_t dataset_id = H5Dopen(loc_id,name,H5P_DEFAULT);
				hid_t dataspace = H5Dget_space(dataset_id);
				int rank = H5Sget_simple_extent_ndims(dataspace);
				hsize_t *dims_out = new hsize_t[rank];
				H5Sget_simple_extent_dims(dataspace,dims_out,NULL);
				std::ostringstream dimsString;
				for(int i=0;i<rank;i++)
				{
					if(i !=rank-1)
						dimsString <<dims_out[i]<<"x";
					else
						dimsString <<dims_out[i];
				}
				HDFTreeDataMap.insert(std::pair<std::string,std::string>(name,dimsString.str()));
				}
                break;
            case H5O_TYPE_NAMED_DATATYPE:
//				HDFTreeDataMap.insert(std::pair<std::string,std::string>(name,name));
                printf ("%s  (Datatype)\n", name);
                break;
            default:
//				HDFTreeDataMap.insert(std::pair<std::string,std::string>(name,name));
                printf ("%s  (Unknown)\n", name);
        }
	return 0;
}

herr_t op_func_l(hid_t loc_id, const char *name, const H5L_info_t *info, void *operator_data)
{
	herr_t status;
	H5O_info_t infobuf;
	switch(info->type){
	case H5L_TYPE_HARD:
		printf("Hard link");
		break;
	case H5L_TYPE_SOFT:
		printf("Soft link");
		break;
	case H5L_TYPE_EXTERNAL:
		printf("External link");
		{
#define HDF5_MAXSTRLEN 1024
			char namebuf[HDF5_MAXSTRLEN];
			std::string fullname;
			// Open group by name group_name
			hid_t group_id = H5Gopen2(loc_id, name, H5P_DEFAULT);

			// Count how many datasets in the group
			hsize_t num_datasets;
			status = H5Gget_num_objs(group_id, &num_datasets);

			// Iterate through group collecting all dataset names
			for(hsize_t i = 0; i < num_datasets; i++)
			{
				H5Gget_objname_by_idx(group_id, i, namebuf, HDF5_MAXSTRLEN);
				fullname = name;
				fullname +="/";
				fullname +=namebuf;
				status = H5Oget_info_by_name(loc_id, fullname.c_str(), &infobuf, H5P_DEFAULT);
				op_func_o(loc_id, fullname.c_str(), &infobuf, operator_data);
			}

			// Close group
			status = H5Gclose(group_id);
		}
		return 0;
	case H5L_TYPE_ERROR:
		printf("Error");
		break;
	}
	status = H5Oget_info_by_name(loc_id, name, &infobuf, H5P_DEFAULT);
	return op_func_o(loc_id, name, &infobuf, operator_data);
}


void CCPiNexusTreeModel::setupModelData(std::string filename, CCPiNexusTreeItem *parent)
{

	hid_t file;
	herr_t	status;
	std::map<std::string,std::string> resultTreeMap;

	file = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);

	HDFDataMutex.lock();
	HDFTreeDataMap.clear();
//	status = H5Ovisit (file, H5_INDEX_NAME, H5_ITER_NATIVE, op_func_o , NULL);
//	for(std::map<std::string,std::string>::iterator itr=HDFTreeDataMap.begin();itr != HDFTreeDataMap.end();itr++)
//		resultTreeMap.insert(std::pair<std::string,std::string>(itr->first,itr->second));
//	HDFTreeDataMap.clear();
	hid_t lapl_id = H5Pcreate(H5P_LINK_ACCESS); 
	status = H5Pset_elink_prefix( lapl_id, ".");
	status = H5Lvisit_by_name (file, ".", H5_INDEX_NAME, H5_ITER_NATIVE, op_func_l , NULL, lapl_id);
	for(std::map<std::string,std::string>::iterator itr=HDFTreeDataMap.begin();itr != HDFTreeDataMap.end();itr++)
		resultTreeMap.insert(std::pair<std::string,std::string>(itr->first,itr->second));
	HDFDataMutex.unlock();

    QList<CCPiNexusTreeItem*> parents;
    parents << parent;


	for(std::map<std::string,std::string>::iterator itr=resultTreeMap.begin();itr!=resultTreeMap.end();itr++)
	{
		QList<QVariant> columnData;
		columnData << itr->first.c_str();
		columnData << itr->second.c_str();

		parents.last()->appendChild(new CCPiNexusTreeItem(columnData, parents.last()));
	}
}


