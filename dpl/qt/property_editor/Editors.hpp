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

#include <QApplication>
#include <QComboBox>
#include <QKeyEvent>
#include <QLineEdit>
#include <QMenu>
// #include <QCore>


namespace dpl::qt::property_editor
{    
  class QComboBoxCustom : public QComboBox
  {
    Q_OBJECT

  public:
    QComboBoxCustom(QWidget* parent = nullptr) : QComboBox{parent} {
      connect(this, QOverload<int>::of(&QComboBox::activated), this, [this](){ 
        QApplication::postEvent(this, new QKeyEvent(QKeyEvent::KeyPress, ::Qt::Key_Enter, ::Qt::NoModifier));                
      });
    }

  protected:
    void showEvent(QShowEvent* e) override {
      QComboBox::showEvent(e);
      
      showPopup();
    }

    // void keyPressEvent(QKeyEvent *e) override {
    //   if (e->key() == Qt::Key_Enter || e->key() == Qt::Key_Return)
    //     accept_ = true;
    //   QLineEdit::keyPressEvent(e);
    // }

    // void focusOutEvent(QFocusEvent *e) override {        
    //   QApplication::postEvent(this, new QKeyEvent(QKeyEvent::KeyPress, Qt::Key_Escape, Qt::NoModifier));                
    // }
  };



  class QLineEditCustom : public QLineEdit
  {
    Q_OBJECT

  public:
    QLineEditCustom(QWidget* parent = nullptr) : QLineEdit{parent} {
    }

    bool accept_ = false;
    bool reset_ = false;
    
    std::unique_ptr<QMenu> menu_;    
    
  protected:      
    void focusInEvent(QFocusEvent* e) override {      
      QLineEdit::focusInEvent(e);

      if (e->reason() != ::Qt::FocusReason::PopupFocusReason)
        setCursorPosition(cursorPositionAt(mapFromGlobal(QCursor::pos())));
    }

    void focusOutEvent(QFocusEvent *e) override {
      QLineEdit::focusOutEvent(e); 
      
      if (e->reason() != ::Qt::FocusReason::PopupFocusReason)
        QApplication::postEvent(this, new QKeyEvent(QKeyEvent::KeyPress, ::Qt::Key_Escape, ::Qt::NoModifier));  
    }


    void contextMenuEvent(QContextMenuEvent* e) override {      
      menu_.reset(createStandardContextMenu());
      menu_->addSeparator();
      
      connect(  
        menu_->addAction("Reset"), &QAction::triggered, this, [this]() {
          reset_ = true;
          QApplication::postEvent(this, new QKeyEvent(QKeyEvent::KeyPress, ::Qt::Key_Enter, ::Qt::NoModifier));
        });

      menu_->popup(e->globalPos());
    }

    void keyPressEvent(QKeyEvent *e) override {
      QLineEdit::keyPressEvent(e);
      
      if (!reset_ && (e->key() == ::Qt::Key_Enter || e->key() == ::Qt::Key_Return)) {
        accept_ = true;
      }
      
    }
  };


  
}




