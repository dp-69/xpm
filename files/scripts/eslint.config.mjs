/* eslint-disable quotes */

import globals from "globals";
import pluginJs from "@eslint/js";

export default [
  {languageOptions: { globals: globals.browser }},
  pluginJs.configs.recommended,
  {
    rules: {
      // "no-multiple-empty-lines" : ["warn", { "max" : 4}],
      // "no-console" : "warn",
      "arrow-spacing": ["error", { "before": true, "after": true }],
      "dot-notation": "error",
      "indent": ["warn", 2, { SwitchCase: 1 }],
      "quotes": ["error", "single"],
      "no-path-concat": "error",
      "no-undef": "error",
      "no-redeclare": "error",
      "no-unused-vars": "warn",
      "no-unreachable": "warn",
      "no-shadow": "warn",
      "prefer-const": "warn",
      "comma-dangle": "warn",
      "no-var": "error",
      // "no-eq-null" : "error",
      "eqeqeq": "warn",
      "no-use-before-define": "error",
      "semi": ["error", "always"]
    }
  }
];