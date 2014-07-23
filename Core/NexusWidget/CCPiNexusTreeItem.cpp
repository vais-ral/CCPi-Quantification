/**
 * Nexus Tree item to store nexus 3D data information.
 *
 */

#include <QStringList>

#include "CCPiNexusTreeItem.h"

CCPiNexusTreeItem::CCPiNexusTreeItem(const QList<QVariant> &data, CCPiNexusTreeItem *parent)
{
    parentItem = parent;
    itemData = data;
}

CCPiNexusTreeItem::~CCPiNexusTreeItem()
{
    qDeleteAll(childItems);
}

void CCPiNexusTreeItem::appendChild(CCPiNexusTreeItem *item)
{
    childItems.append(item);
}


CCPiNexusTreeItem *CCPiNexusTreeItem::child(int row)
{
    return childItems.value(row);
}


int CCPiNexusTreeItem::childCount() const
{
    return childItems.count();
}

int CCPiNexusTreeItem::columnCount() const
{
    return itemData.count();
}

QVariant CCPiNexusTreeItem::data(int column) const
{
    return itemData.value(column);
}

CCPiNexusTreeItem *CCPiNexusTreeItem::parent()
{
    return parentItem;
}

int CCPiNexusTreeItem::row() const
{
    if (parentItem)
        return parentItem->childItems.indexOf(const_cast<CCPiNexusTreeItem*>(this));

    return 0;
}
