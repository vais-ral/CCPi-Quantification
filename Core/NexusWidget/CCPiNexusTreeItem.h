/**
 * This is the nexus data item in the tree. currently only stores the 3D datasets.
 *
 */

#ifndef CCPINEXUSTREEITEM_H
#define CCPINEXUSTREEITEM_H

#include <QList>
#include <QVariant>

class CCPiNexusTreeItem
{
public:
    CCPiNexusTreeItem(const QList<QVariant> &data, CCPiNexusTreeItem *parent = 0);
    ~CCPiNexusTreeItem();

    void appendChild(CCPiNexusTreeItem *child);

    CCPiNexusTreeItem *child(int row);
    int childCount() const;
    int columnCount() const;
    QVariant data(int column) const;
    int row() const;
    CCPiNexusTreeItem *parent();

private:
    QList<CCPiNexusTreeItem*> childItems;
    QList<QVariant> itemData;
    CCPiNexusTreeItem *parentItem;
};


#endif
