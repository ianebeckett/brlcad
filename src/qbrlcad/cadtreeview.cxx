#include "cadtreeview.h"

#include <QPainter>
#include <qmath.h>

void GObjectDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    /*
        // rect with width proportional to value
        QRect rect(option.rect.adjusted(4,4,-4,-4));

        // draw the value bar
        painter->fillRect(rect, QBrush(QColor("steelblue")));
*/
        // calculate text min width
        QString text = index.data().toString();
        // draw value text centered
        painter->drawText(option.rect, text, QTextOption(Qt::AlignLeft));
}

QSize GObjectDelegate::sizeHint(const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    return option.fontMetrics.size(Qt::TextSingleLine, index.data().toString());
}


CADTreeView::CADTreeView(QWidget *pparent) : QTreeView(pparent)
{
}

void CADTreeView::header_state()
{
    header()->setStretchLastSection(true);
    int header_length = header()->length();
    header()->setStretchLastSection(false);
    header()->setSectionResizeMode(0, QHeaderView::Fixed);
    resizeColumnToContents(0);
    int column_width = columnWidth(0);
    header()->setSectionResizeMode(0, QHeaderView::ResizeToContents);
    if (column_width > header_length) {
	header()->setStretchLastSection(false);
    } else {
	header()->setStretchLastSection(true);
    }
}

void CADTreeView::resizeEvent(QResizeEvent *)
{
    header_state();
}

void CADTreeView::tree_column_size(const QModelIndex &)
{
    header_state();
}


/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */

