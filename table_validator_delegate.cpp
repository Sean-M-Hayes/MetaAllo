/**************************************************************************
**    Copyright 2017 Sean M. Hayes
**    This file is part of MetaAllo.
**
**    MetaAllo is free software: you can redistribute it and/or modify
**    it under the terms of the GNU General Public License as published by
**    the Free Software Foundation, either version 3 of the License, or
**    (at your option) any later version.

**    MetaAllo is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**    GNU General Public License for more details.

**    You should have received a copy of the GNU General Public License
**    along with MetaAllo.  If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/
#include <table_validator_delegate.h>
#include <QLineEdit>

//Int

TableValidatorDelegate_Int::TableValidatorDelegate_Int(QObject *parent) : QStyledItemDelegate(parent)
{
}

QWidget* TableValidatorDelegate_Int::createEditor(QWidget* parent, const QStyleOptionViewItem &option,const QModelIndex &index) const
{
    QLineEdit* editor = new QLineEdit(parent);
    QRegExp single_int_regex("^\\d+$");
    QRegExpValidator *single_int_validator = new QRegExpValidator(single_int_regex);

    editor->setValidator(single_int_validator);
    return editor;
}

void TableValidatorDelegate_Int::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    int value = index.model()->data(index,Qt::EditRole).toInt();
    QLineEdit* line = static_cast<QLineEdit*>(editor);
    line->setText(QString().setNum(value));
}

void TableValidatorDelegate_Int::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    QLineEdit* line = static_cast<QLineEdit*>(editor);
    QString value = line->text();
    model->setData(index,value);
}

void TableValidatorDelegate_Int::updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    editor->setGeometry(option.rect);
}

//Double

TableValidatorDelegate_Double::TableValidatorDelegate_Double(QObject *parent) : QStyledItemDelegate(parent)
{
}

QWidget* TableValidatorDelegate_Double::createEditor(QWidget* parent, const QStyleOptionViewItem &option,const QModelIndex &index) const
{
    QLineEdit* editor = new QLineEdit(parent);

    QRegExp single_double_regex("^\\d*\\.?\\d+$");
    QRegExpValidator *single_double_validator = new QRegExpValidator(single_double_regex);

    editor->setValidator(single_double_validator);
    return editor;
}

void TableValidatorDelegate_Double::setEditorData(QWidget *editor, const QModelIndex &index) const
{
    double value = index.model()->data(index,Qt::EditRole).toDouble();
    QLineEdit* line = static_cast<QLineEdit*>(editor);
    line->setText(QString().setNum(value));
}

void TableValidatorDelegate_Double::setModelData(QWidget *editor, QAbstractItemModel *model, const QModelIndex &index) const
{
    QLineEdit* line = static_cast<QLineEdit*>(editor);
    QString value = line->text();
    model->setData(index,value);
}

void TableValidatorDelegate_Double::updateEditorGeometry(QWidget *editor, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    editor->setGeometry(option.rect);
}
