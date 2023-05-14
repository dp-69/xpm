/*
 * This file is part of Rapid Reservoir Modelling Software.
 *   | https://rapidreservoir.org/
 *   | https://bitbucket.org/rapidreservoirmodelling/rrm/
 *
 * Copyright (c) 2022
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * RRM is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * RRM is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with RRM. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#pragma once

#include <dpl/static_vector.hpp>
#include <dpl/qt/property_editor/Editors.hpp>

#include <utility>
#include <vector>
#include <string>
#include <regex>

#include <QAbstractItemModel>
#include <QComboBox>
#include <QPlainTextEdit>
#include <QLineEdit>
#include <QWidget>
#include <QVariant>
#include <QTimer>
#include <QHBoxLayout>



namespace dpl::qt::property_editor
{
  const inline std::string NumberRegEx = "([-+.eE\\d]+)";
  const inline std::string SeparatorRegEx = "[,;\\s]+";

  class PropertyCategory;

  class PropertyBase
  {
    friend class PropertyCategory;
    friend class PropertyModel;

  protected:
    auto* GetVisibleChildren(int i) {      
      for (auto* c : children)
        if (c->IsVisible() && i-- == 0)
          return c;
    }

    int VisibleChildrenCount() {
      int i = 0;
      for (auto* c : children)
        if (c->IsVisible())
          ++i;
      
      return i;      
    }
    
  public:
    PropertyCategory* parent;
    std::vector<PropertyBase*> children;
    
    int VisibleIndex();
    
    virtual ~PropertyBase() = default;

    virtual bool IsVisible() = 0;    
    virtual std::string Name() = 0;
    virtual std::string Tooltip() { return {}; }
  };

  inline auto* get_ptr(const QModelIndex& index) {
    return static_cast<PropertyBase*>(index.internalPointer());
  }  

  template<typename T>
  auto* try_ptr(const QModelIndex& index) {
    return dynamic_cast<T*>(get_ptr(index));
  }  

 
  
  


  class PropertyItem : public PropertyBase
  {        
  public:
    virtual void SetModelData(QWidget *widget, const QModelIndex& index) = 0;    
    virtual QWidget* CreateEditor(QWidget *parent, const QModelIndex& index) = 0;

    virtual bool IsReadOnly() = 0;
    virtual bool IsActive() = 0;
    virtual void SetActive() {}

    virtual QVariant GetDisplayValueQVariant(Qt::ItemDataRole role) = 0;

    QColor color;
  };
  

  class PropertyCategory : public PropertyBase
  {
    friend class PropertyModel;
    
    std::string name_;
    std::string tooltip_;
    
  public:
    explicit PropertyCategory(
      std::string_view name = {},
      std::string_view tootip = {}
    ) : name_{name}, tooltip_{tootip} {}

    std::string Name() override { return name_; }

    std::string Tooltip() override { return tooltip_; }  

    bool IsVisible() override { return true; }
  };


  inline int PropertyBase::VisibleIndex() {    
    int i = 0;
    for (auto* ptr : parent->children) {
      if (ptr == this)
        return i;
      
      if (ptr->IsVisible())
        ++i;
    }
  }



  static QWidget* CreateLineEditor(QWidget *parent, const QString& text) {
    auto* widget = new QWidget{parent};
    auto* layout = new QHBoxLayout;
    layout->setContentsMargins(1, 0, 0, 0);
    widget->setLayout(layout);
    widget->setAutoFillBackground(true);
    auto* editor = new QLineEditCustom{parent};
    editor->setFrame(false);
    editor->setText(text);
    layout->addWidget(editor);
    widget->setFocusProxy(editor);       
    return widget;
  }
  
  template <typename Scalar> struct ScalarConvert {};

  template<> struct ScalarConvert<bool> {
    static bool FromString(const QString& str, bool& value) {
      value = QVariant{str}.toBool();
      return true;
    }

    static auto ToQString(bool value) { return QString::number(value); }
    static auto ToString(bool value) { return std::to_string(value); } 
  };
  
  template<> struct ScalarConvert<double> {
    static bool FromString(const QString& str, double& value) {
      bool valid;
      auto str_copy = str;
      value = str_copy.replace(",", "").toDouble(&valid);
      return valid;
    }

    static auto ToQString(double value) { return QString::number(value); }
    static auto ToString(double value) { return std::to_string(value); } 
  };

  template<> struct ScalarConvert<int> {
    static bool FromString(const QString& str, int& value) {
      bool valid;
      value = str.toInt(&valid);
      return valid;
    }
    
    static auto ToQString(int value) { return QString::number(value); }
    static auto ToString(int value) { return std::to_string(value); } 
  };

  template<> struct ScalarConvert<QString> {
    static auto FromString(const QString& str, QString& value) {
      value = str;
      return true;
    }
    
    static auto ToQString(const QString& value) { return value; }
    static auto ToString(const QString& value) { return value.toStdString(); }
  };

  template<> struct ScalarConvert<std::string> {
    static auto FromString(const QString& str, std::string& value) {
      value = str.toStdString();
      return true;
    }

    static auto ToQString(const std::string& value) { return QString::fromStdString(value); }
    static auto ToString(const std::string& value) { return value; }
  };


  
  
  

  template<typename Scalar>
  struct PropertyItemTextEdit
  {
    template<typename Inter>
    static QWidget* CreateEditor(QWidget *parent, const QModelIndex& index, Inter& inter) {
      if constexpr (std::is_same_v<Scalar, double>) {
        auto val = inter.Get();
        
        return CreateLineEditor(parent,
          std::abs(val) == std::numeric_limits<double>::max() ? "" : QString::fromStdString(inter.Format(val)));        
      }
      else
        return CreateLineEditor(parent,
          QString::fromStdString(inter.Format(inter.Get()))
        );
    }

    template<typename Inter>
    static void SetModelData(QWidget *widget, const QModelIndex &index, Inter& inter) {
      if (auto* editor = static_cast<QLineEditCustom*>(widget->focusProxy()); editor->reset_)
        inter.Reset();
      else if (editor->accept_) {
        if (auto text = editor->text(); text.isEmpty())
          inter.Reset();
        else if (Scalar value; ScalarConvert<Scalar>::FromString(text, value))
          inter.Set(value);
      }
    }


    template<typename Inter>
    static QVariant GetModelValue(Inter& inter, Qt::ItemDataRole role) {
      if (role != Qt::ItemDataRole::DisplayRole)
        return {};

      auto value = inter.Get();

      if constexpr (std::is_same_v<Scalar, double>)
        if (value == std::numeric_limits<double>::max())
          return {};

      auto str = inter.Format(value);
      if (inter.Calculated())
        str.append(" : auto");

      return QString::fromStdString(str);
    }
  };



  
  template<typename ValueType> struct PropertyItemTraits {};

  template<> struct PropertyItemTraits<std::string> : PropertyItemTextEdit<std::string> {};
  template<> struct PropertyItemTraits<int> : PropertyItemTextEdit<int> {};
  template<> struct PropertyItemTraits<double> : PropertyItemTextEdit<double> {};

  template<> struct PropertyItemTraits<bool>
  {
    template<typename Inter>
    static QWidget* CreateEditor(QWidget *parent, const QModelIndex& index, Inter& inter) {         
      inter.Set(!inter.Get());          
      emit const_cast<QAbstractItemModel*>(index.model())->dataChanged(index, index);
      return nullptr;
    }

    template<typename Inter>
    static void SetModelData(QWidget *widget, const QModelIndex &index, Inter& inter) { }

    template<typename Inter>
    static QVariant GetModelValue(Inter& inter, Qt::ItemDataRole role) {
      if (role != Qt::ItemDataRole::CheckStateRole)
        return {};

      return inter.Get() ? Qt::Checked : Qt::Unchecked;
    }
  };



  
  template<typename Scalar, int Count>
  struct PropertyItemTraits<dpl::vector_n<Scalar, Count>>
  {
    template<typename Inter>
    static QWidget* CreateEditor(QWidget *parent, const QModelIndex&, Inter& inter) {
      return CreateLineEditor(parent, QString::fromStdString(GetTextStream(inter).str()));
    }

    template<typename Inter>
    static void SetModelData(QWidget *widget, const QModelIndex&, Inter& inter) {
      if (auto* editor = static_cast<QLineEditCustom*>(widget->focusProxy()); editor->reset_)
        inter.Reset();
      else if (editor->accept_) {
        dpl::vector_n<Scalar, Count> arr;

        if (auto text = editor->text(); text.isEmpty())
          inter.Reset();
        else {
          if (Scalar value; ScalarConvert<Scalar>::FromString(text, value)) {
            arr = value;
            inter.Set(arr);
          }
          else {
            std::stringstream ss;
            ss << NumberRegEx;
            dpl::sfor<Count - 1>([&]() { ss << SeparatorRegEx << NumberRegEx; });

            if (auto m = QRegularExpression(QString::fromStdString(ss.str())).match(text); m.hasMatch()) {
              dpl::sfor<Count>([&](auto i) { ScalarConvert<Scalar>::FromString(m.captured(i + 1), arr[i]); });
              inter.Set(arr);        
            }
            else if constexpr (Count == 3)
              if (auto m32 = QRegularExpression(
                QString::fromStdString(NumberRegEx + SeparatorRegEx + NumberRegEx)).match(text); m32.hasMatch()) {
                ScalarConvert<Scalar>::FromString(m32.captured(1), arr[0]);
                ScalarConvert<Scalar>::FromString(m32.captured(1), arr[1]);
                ScalarConvert<Scalar>::FromString(m32.captured(2), arr[2]);
                inter.Set(arr);
              }
          }     
        }                        
      }
    }

    template<typename Inter>
    static auto GetTextStream(Inter& inter) {
      auto arr = inter.Get();

      std::stringstream ss;
      
      dpl::sfor<Count>([&](auto i) {
        if constexpr (i == 0)
          ss << '[' << inter.template Format<0>(arr[i]);
        else
          ss << ", " << inter.template Format<i>(arr[i]);
      });
      ss << ']';

      return ss;
    }

    template<typename Inter>
    static QVariant GetModelValue(Inter& inter, Qt::ItemDataRole role) {
      if (role != Qt::ItemDataRole::DisplayRole)
        return {};

      auto str = GetTextStream(inter);
      if (inter.Calculated())
        str << " : auto";        

      return QString::fromStdString(str.str());
    }
  };






  







  








  template<typename Functor>
  class FunctionGenericItem : public PropertyItem
  {
    using Type = typename Functor::Type;
    using Traits = PropertyItemTraits<Type>;

    Functor functor_;
    
  public:
    explicit FunctionGenericItem(Functor functor) : functor_{functor} {}   

    std::string Name() override {
      return functor_.Name();
    }

    std::string Tooltip() override {
      return functor_.Tooltip();
    }  

    bool IsVisible() override {
      return functor_.IsVisible();
    }

    bool IsReadOnly() override {
      return functor_.IsReadOnly();
    }

    bool IsActive() override {
      return functor_.IsActive();
    }

    void SetActive() override {
      functor_.SetActive();
    }
    

    QVariant GetDisplayValueQVariant(Qt::ItemDataRole role) override {
      return Traits::GetModelValue(functor_, role);
    }
    
    QWidget* CreateEditor(QWidget *parent, const QModelIndex& index) override {
      return Traits::CreateEditor(parent, index, functor_);
    }

    void SetModelData(QWidget* widget, const QModelIndex& index) override {
      Traits::SetModelData(widget, index, functor_);
    }
  };
  

  
  



  template<typename Functor>
  class ComboBoxItem : public PropertyItem
  {
    Functor functor_;
    using ComboItem = typename Functor::ComboItem;
    
  public:
    std::string Name() override {
      return functor_.Name();
    }

    explicit ComboBoxItem(const Functor& inter) : functor_(inter) {}

    QVariant GetDisplayValueQVariant(Qt::ItemDataRole role) override {
      if (role == Qt::DisplayRole)
        return QString::fromStdString(Functor::Caption(functor_.Get()));

      return {};
    }
    
    bool IsVisible() override { return true; }

    bool IsReadOnly() override { return false; }

    bool IsActive() override { return false; }

    QWidget* CreateEditor(QWidget *parent, const QModelIndex& index) override {
      auto* editor = new QComboBoxCustom(parent);

      if constexpr (std::is_same_v<ComboItem, int>)
        for (int i = 0, count = functor_.Count(); i < count; ++i)
          editor->addItem(QString::fromStdString(Functor::Caption(i)));
      else
        for (auto* item : functor_.Items())
          editor->addItem(QString::fromStdString(Functor::Caption(item)), ItemToQVariant<ComboItem>(item));

      editor->setCurrentText(QString::fromStdString(Functor::Caption(functor_.Get())));
      return editor;            
    }
  
  
    void SetModelData(QWidget* widget, const QModelIndex& index) override {            
      auto combo = static_cast<QComboBox*>(widget);

      if constexpr (std::is_same_v<ComboItem, int>)
        functor_.Set(combo->currentIndex());
      else
        functor_.Set(ItemFromQVariant<ComboItem>(combo->itemData(combo->currentIndex())));
    }
  };


  class PropertyModel : public QAbstractItemModel
  {
  Q_OBJECT

    std::vector<std::unique_ptr<PropertyBase>> all_items_;
    PropertyCategory* root_;

    template <typename Functor, typename = std::void_t<>>
    struct parse_property_item {
      using type = FunctionGenericItem<Functor>;
    };

    template <typename Functor>
    struct parse_property_item<Functor, std::void_t<typename Functor::ComboItem>> {
      using type = ComboBoxItem<Functor>;
    };

    template <typename Functor>
    using property_item_t = typename parse_property_item<Functor>::type;
  public:      

    explicit PropertyModel(QObject* parent = nullptr) : QAbstractItemModel(parent) {
      root_ = new PropertyCategory;
      all_items_.emplace_back(root_);
    }


    // TODO 
    // template <typename Getter, typename Setter>
    // PropertyItem* AddBool(PropertyCategory* category, Functor functor) {
    //
    //
    //   // FunctorPropertyItem<>
    //
    // //   auto add_bool = [this, &tree_model](
    // //   dpl::qt::property_editor::PropertyCategory* cat, std::string_view name, const auto& get, const auto& set) {
    // //   tree_model.AddItem(cat, FunctorPropertyItem<
    // //     decltype(get), bool, decltype(set), FinaliserRender>{fd_widget_, name, get, set});
    // // };
    // //
    // // auto add_bool_visual = [add_bool, visual_category](
    // //   std::string_view name, const auto& get, const auto& set) {
    // //   add_bool(visual_category, name, get, set);
    // // };
    //   
    //   
    //
    //   
    //
    //   PropertyItem* ptr = new property_item_t<Functor>{functor};
    //   ptr->parent = category;
    //   category->children.push_back(ptr);
    //   all_items_.emplace_back(ptr);
    //   return ptr;
    // }
    
    template <typename Functor>
    PropertyItem* AddItem(PropertyCategory* category, Functor functor) {
      PropertyItem* ptr = new property_item_t<Functor>{functor};
      ptr->parent = category;
      category->children.push_back(ptr);
      all_items_.emplace_back(ptr);
      return ptr;
    }

    template <typename Functor>
    auto* AddItem(Functor functor) {
      return AddItem(root_, functor);
    }

    auto* AddCategory(PropertyCategory* category, std::string_view name, std::string_view tooltip = {}) {      
      auto* ptr = new PropertyCategory{name, tooltip};
      ptr->parent = category;
      category->children.push_back(ptr);
      all_items_.emplace_back(ptr);
      return ptr;
    }

    auto* AddCategory(std::string_view name) {      
      return AddCategory(root_, name);
    }

    QModelIndex index(int row, int column, const QModelIndex& parent) const override {
      return createIndex(row, column,
        (parent == QModelIndex{} ? root_ : get_ptr(parent))->GetVisibleChildren(row));
    }

    QModelIndex parent(const QModelIndex& index) const override {
      auto* base = get_ptr(index);
      return !base || base->parent == root_ ? QModelIndex{} : createIndex(base->VisibleIndex(), 0, base->parent);
      // [base] is sometimes [nullptr] in Julio's test cases.
    }

    int rowCount(const QModelIndex& parent) const override {
      return (parent == QModelIndex{} ? root_ : get_ptr(parent))->VisibleChildrenCount();            
    }

    int columnCount(const QModelIndex& parent) const override {
      return 2;
    }

    QModelIndex index(PropertyBase* data) const {
      return createIndex(data->VisibleIndex(), 0, data);
    }


    QVariant data(const QModelIndex& index, int /*Qt::ItemDataRole*/ role) const override {
      auto* base = get_ptr(index);

      if (role == Qt::ToolTipRole)
        return QString::fromStdString(base->Tooltip());

      if (index.column() == 0) {
        if (auto* item = dynamic_cast<PropertyItem*>(base); item) {
          if (role == Qt::FontRole) {
            QFont font;
            font.setUnderline(item->IsActive());
            return font;
          }

          if (role == Qt::DisplayRole)
            return QString::fromStdString(base->Name());
          
          if (item->color.isValid()) {
            if (role == Qt::BackgroundRole)
              return item->color;
            if (role == Qt::ForegroundRole)
              return item->color.lightnessF() < 0.5 ? QColor(Qt::white) : QColor(Qt::black);
          }
        }
      }
      else if (index.column() == 1) {
        if (auto* item = dynamic_cast<PropertyItem*>(base); item)
          return item->GetDisplayValueQVariant(static_cast<::Qt::ItemDataRole>(role));

        /* Category caption is assigned to the 2nd column to disregard
           its size when autofitting the 1st column width. */
        if (role == ::Qt::DisplayRole)
          return QString::fromStdString(base->Name());
      }

      return {};
    }


    ::Qt::ItemFlags flags(const QModelIndex& index) const override {            
      if (index.column() == 1 && try_ptr<PropertyItem>(index))
        return ::Qt::ItemIsEnabled | ::Qt::ItemIsEditable;

      return ::Qt::ItemIsEnabled;
    }

    QVariant headerData(int section, ::Qt::Orientation orientation, int role) const override {
      if (role == ::Qt::DisplayRole)
        return section == 0 ? "Property" : "Value";

      return {};
    }
  };
}