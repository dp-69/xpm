#pragma once

#include <QColor>
#include <nlohmann/json.hpp>

namespace dpl::qt
{
  inline void try_parse(const nlohmann::json& j, QColor& color) {
    if (auto c = j.find("rgb"), end = j.end(); c != end)
      color = QColor::fromRgb((*c)[0], (*c)[1], (*c)[2]);
    else if (c = j.find("rgb_f"); c != end)
      color = QColor::fromRgbF((*c)[0], (*c)[1], (*c)[2]);
    else if (c = j.find("hsl"); c != end)
      color = QColor::fromHsl((*c)[0], (*c)[1], (*c)[2]);
    else if (c = j.find("hsl_f"); c != end)
      color = QColor::fromHslF((*c)[0], (*c)[1], (*c)[2]);
  }
}
