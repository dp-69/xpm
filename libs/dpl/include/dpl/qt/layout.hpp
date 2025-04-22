/*
 * This file is part of Dmytro Petrovskyy Library (dpl).
 *
 * Copyright (c) 2024
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * dpl is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * dpl is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with dpl. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <QBoxLayout>
#include <QSplitter>

namespace dpl::qt::layout
{
  struct spacing
  {
    int value;
    operator int() const { return value; }
  };

  struct stretch
  {
    int value;
    operator int() const { return value; }
  };

  enum class dir
  {
    horizontal,
    vertical
  };
  
  struct box
  {
    QBoxLayout* qbox;

    box(dir d, int spacing = 0) {
      qbox = new QBoxLayout{d == dir::horizontal ? QBoxLayout::LeftToRight : QBoxLayout::TopToBottom};
      qbox->setSpacing(spacing);
      qbox->setContentsMargins(0, 0, 0, 0);
    }

    operator QBoxLayout*() const { return qbox; }
  };

  struct splitter
  {
    QSplitter* qsplitter;

    splitter(dir d) {
      qsplitter = new QSplitter{d == dir::horizontal ? Qt::Horizontal : Qt::Vertical};
    }

    operator QSplitter*() const { return qsplitter; }
  };

  struct widget
  {
    QWidget* qwidget;

    widget(QLayout* l) {
      qwidget = new QWidget{};
      qwidget->setContentsMargins(0, 0, 0, 0);
      qwidget->setLayout(l);
    }

    operator QWidget*() const { return qwidget; }
  };

  namespace resolve
  {
    inline void add(QBoxLayout* l, const std::tuple<QWidget*, stretch>& arg) {
      l->addWidget(std::get<QWidget*>(arg), std::get<stretch>(arg));
    }
    
    inline void add(QBoxLayout* l, QWidget* arg) { l->addWidget(arg); }
    inline void add(QBoxLayout* l, QLayout* arg) { l->addLayout(arg); }
    inline void add(QBoxLayout* l, QWidget& arg) { l->addWidget(&arg); }
    inline void add(QBoxLayout* l, QLayout& arg) { l->addLayout(&arg); }
    inline void add(QBoxLayout* l, stretch arg) { l->addStretch(arg); }
    inline void add(QBoxLayout* l, spacing arg) { l->addSpacing(arg); }

    inline void add(QSplitter* s, QWidget* w) { s->addWidget(w); }
    inline void add(QSplitter* s, stretch str) { s->setStretchFactor(s->count() - 1, str); }
  }

  template <typename Item>
  box operator<<(box b, Item&& arg) {
    resolve::add(b, std::forward<Item>(arg));
    return b;
  }

  template <typename Item>
  box operator|(box b, Item&& arg) {
    resolve::add(b, std::forward<Item>(arg));
    return b;
  }

  template <typename Item>
  splitter operator<<(splitter s, Item&& arg) {
    resolve::add(s, std::forward<Item>(arg));
    return s;
  }
}
