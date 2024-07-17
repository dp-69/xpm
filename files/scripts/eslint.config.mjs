import globals from "globals";
import pluginJs from "@eslint/js";


export default [
  {languageOptions: { globals: globals.browser }},
  pluginJs.configs.recommended,
  {
    rules: {
      "no-unused-vars": "warn",
      "no-unreachable": "warn",      
      "no-shadow" : "warn",
      "comma-dangle" : "warn",
      // "no-eq-null" : "error",
      "eqeqeq" : "warn",
      "no-use-before-define": "error",
      "semi": ["error", "always"]
    }
  }
];