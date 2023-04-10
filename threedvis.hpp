#pragma once

#include <QWidget>


namespace xpm
{
  class XPMWidget : public QWidget
  {
    
  };
}

// // Global module fragment where #includes can happen
// module;
// #include <iostream>
//
// // first thing after the Global module fragment must be a module command
// export module fooMODULE;
//
// // char const* world() { return "hello world\n"; }
//
// export class foo {
// public:
//   foo();
//   ~foo();
//   void helloworld();
// };
//
// foo::foo() = default;
// foo::~foo() = default;
// void foo::helloworld() { std::cout << "hello world\n"; }
//
// // export module threedvisMODULE;
// //
// // #include <QWidget>;
// //
// //
// //
// //
// // namespace xpm
// // {
// //   export class FDWidget : public QWidget
// //   {
// //     
// //   };
// // }
