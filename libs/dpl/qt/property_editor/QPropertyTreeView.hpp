/*
 * This file is part of Dmytro Petrovskyy Library (DPL).
 *
 * Copyright (c) 2023
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * DPL is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * DPL is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with RRM. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <dpl/qt/property_editor/PropertyEntities.hpp>

#include <QTreeView>
#include <QHeaderView>
#include <QKeyEvent>
#include <QStyledItemDelegate>
#include <QPainter>

namespace dpl::qt::property_editor
{
  class PropertyDelegate : public QStyledItemDelegate
  {
  Q_OBJECT   
    
  public:
    QTreeView* tree_view;

    PropertyDelegate(QObject* parent = nullptr) : QStyledItemDelegate(parent) {}

    QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem&, const QModelIndex& index) const override {
      auto* editor = try_ptr<PropertyItem>(index)->CreateEditor(parent, index);
      // editor->installEventFilter(const_cast<PropertyDelegate*>(this));
      return editor;
    }


    // void setEditorData(QWidget *widget, const QModelIndex &index) const override {
    //   auto item = get_ptr<PropertyItem>(index);
    //   item->SetEditorData(widget, index);
    // }
  
    void setModelData(QWidget* widget, QAbstractItemModel*, const QModelIndex& index) const override {
      try_ptr<PropertyItem>(index)->SetModelData(widget, index);
    }

    // bool eventFilter(QObject* object, QEvent* event) {
      // auto* editor = qobject_cast<QWidget*>(object);
      // if (editor && event->type() == QEvent::KeyPress) {
      //   auto* key_event = static_cast<QKeyEvent*>(event);
      //   if (key_event->key() == Qt::Key_Return) {
      //     std::cout << "HERE";
      //     emit commitData(editor); //save changes
      //     auto* line_edit = dynamic_cast<QLineEdit*>(editor);
      //     
      //     std::cout << editor->metaObject()->className() << " " << editor->children().count();
      //     std::cout << " " << editor->children()[0]->metaObject()->className();
      //     std::cout << " " << editor->children()[1]->metaObject()->className();
      //     if (line_edit) {
      //       std::cout << "THIS IS LINE EDIT";
      //       line_edit->selectAll();
      //     }
      //     return true;
      //   }
      // }
      // return false;
    // }


protected:
  // void initStyleOption(QStyleOptionViewItem* option, const QModelIndex& index) const override {
  //   QStyledItemDelegate::initStyleOption(option, index);
  //
  //   if (true/*index.column() == 1*/) {
  //     auto* item = dynamic_cast<PropertyItem*>(get_ptr<PropertyBase>(index));
  //     if (item && item->IsReadOnly())
  //       option->backgroundBrush = QColor{245, 245, 245};
  //   }
  // }

public:
    void paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const override {
      if (index.column() == 0) {
        if (auto* item = try_ptr<PropertyItem>(index)) {
          if (item->IsReadOnly()) {
            painter->save();
            
            painter->setBrush(/*QColor{255, 0, 0}*/QColor{34, 160, 223});
            painter->setPen(::Qt::NoPen);
            painter->setRenderHint(QPainter::Antialiasing, true);

            auto indent = static_cast<const QTreeView*>(option.widget)->indentation() / 2;
            auto size = 3;
            painter->drawEllipse(QPointF(option.rect.left() - indent, option.rect.top() + ceil(option.rect.height()/2.)), size, size);
            

            painter->restore();
          }
        }
      }
            
      if (try_ptr<PropertyCategory>(index)) {                
        if (index.column() == 1) {
          auto* tree_model = index.model();
          auto secondColRect = static_cast<const QTreeView*>(option.widget)->visualRect(
            tree_model->index(index.row(), 0, tree_model->parent(index)));
          
          auto option_row = option;
          option_row.rect.setLeft(secondColRect.left());                            
                    
          QStyledItemDelegate::paint(painter, option_row, index);
          // QStyledItemDelegate::paint(painter, option, index);
        }              
      }
      else
        QStyledItemDelegate::paint(painter, option, index);
    }
  
    void updateEditorGeometry(QWidget* editor, const QStyleOptionViewItem& option, const QModelIndex&) const override {
      editor->setGeometry(option.rect);
    }
  };

  
  class QPropertyTreeView : public QTreeView
  {
    Q_OBJECT


    PropertyDelegate delegate_;
    PropertyModel model_;

    void setModel(PropertyModel*) {}
    void setModel(QAbstractItemModel*) override {}
    
  public:
    

    explicit QPropertyTreeView(QWidget* parent = nullptr) : QTreeView{parent} {
      delegate_.tree_view = this;
      
      setEditTriggers(NoEditTriggers);
      setFocusPolicy(Qt::NoFocus);

      // connect(this, &QTreeView::clicked, this, [this](const QModelIndex& index) {
      //   using namespace property_editor;
      //
      //   if (index.column() == 1 && !dynamic_cast<PropertyCategory*>(get_ptr<PropertyBase>(index))
      //     /*index.parent() != QModelIndex()*/)
      //     edit(index);
      // });
      
      header()->setDefaultAlignment(Qt::AlignCenter);
      // setColumnWidth(0, rendering_settings_.font_size*12);
      setColumnWidth(1, 0);
      header()->setStretchLastSection(true);

      QTreeView::setModel(&model_);
      setItemDelegate(&delegate_);

      setFrameShape(NoFrame);
    }

    PropertyModel* model() {
      return &model_;
    }
    
  protected:
    void mouseReleaseEvent(QMouseEvent* event) override {
      QTreeView::mouseReleaseEvent(event);

      if (event->button() == Qt::LeftButton)
        if (auto index = indexAt(event->pos());
          index.column() == 1 && !try_ptr<PropertyCategory>(index))
          edit(index);
    }

    void mouseDoubleClickEvent(QMouseEvent* event) override {
      QTreeView::mouseDoubleClickEvent(event);

      if (event->button() == Qt::LeftButton)
        if (auto index = indexAt(event->pos()); index.column() == 0)
          if (auto* item = try_ptr<PropertyItem>(index))
            item->SetActive();
    }
  };
}
