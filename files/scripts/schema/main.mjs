/*
 * This file is part of Extensive Pore Modelling (xpm).
 *   | https://github.com/dp-69/xpm
 *
 * Copyright (c) 2024
 *   | Dmytro Petrovskyy, PhD
 *   | dmytro.petrovsky@gmail.com
 *   | https://www.linkedin.com/in/dmytro-petrovskyy/
 *
 * xpm is free software: you can redistribute it and/or modify              
 * it under the terms of the GNU General Public License as published by     
 * the Free Software Foundation, either version 3 of the License, or        
 * (at your option) any later version.                                      
 *                                                                         
 * xpm is distributed in the hope that it will be useful,                   
 * but WITHOUT ANY WARRANTY; without even the implied warranty of           
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            
 * GNU General Public License for more details.                             
 *                                                                         
 * You should have received a copy of the GNU General Public License        
 * along with xpm. If not, see <http://www.gnu.org/licenses/>.
 *
 */

import fs from 'fs';
import stripJsonComments from 'strip-json-comments';

// function optional({required : r}, name) {  // eslint-disable-line no-unused-vars
//   return r !== undefined && r.includes(name) ? '' : '*opt* ';
// }

const valueMap = {
  'tolerance' : x => x.toExponential(2)
};

function recordStr(key, value, level) {
  function typeStr(x) {
    switch (x.type) {
      case 'string':  return 'string';
      case 'boolean': return 'bool';
      case 'integer': return 'int';
      case 'number':  return 'real';
      case 'object':  return 'object';
      case 'array': 
        return x.maxItems !== undefined
          ? typeStr(x.items) + `[${x.maxItems}]`
          : 'array'; /* && x.minItems == x.maxItems */
  
      default:
        return `${x.type}`;
    }
  }
  
  const indentStr = '<br/>\n&ensp;&ensp;&ensp;';

  function descStr({ description: d }/* , level */) {
    return d !== undefined                                      
      ? (level === 0 ? '\n' : indentStr) + d.replaceAll('\n', indentStr)
      : '';
  }

  function valueStr({ enum: e, default: d }) {
    return (
      e !== undefined  // ReSharper disable once QualifiedExpressionMaybeNull
        ? indentStr + `*values* : ${e.map(x => `\`"${x}"\``).join(' | ')}`
        : d !== undefined
          ? indentStr + `*default* : \`${                                   /* eslint-disable indent */
              typeof d === 'string'                                                       
                ? `"${d}"`                                                
                : (valueMap[key] !== undefined ? valueMap[key](d) : d)}\``  /* eslint-enable indent */
          : '');    
  }


  return (
    '\n'
    + ' '.repeat(level < 2 ? 0 : 2*(level - 1))
    + (level === 0 ? '## ' : '- ')
    + `\`${key}\` : ` + typeStr(value)
    + descStr(value, level) + valueStr(value) + '\n'); // <br/>
}

let output =
`# JSON specification for [\`xpm\`](https://github.com/dp-69/xpm) input

Paragraph indent resembles the nesting levels in the JSON input file.

This file has been automatically generated from a machine-readable [JSON schema](https://raw.githubusercontent.com/dp-69/xpm/main/files/xpm.schema.json).
`;

const details = {
  'darcy':
    'Acceptable data table format is `[[x1, y1], [x2, y2], ...]` or `[x1, y1, x2, y2, ...]` where `x(i)` < `x(i+1)`.'
};

function proc(x, level = 0) {
  function keyValues({ type: t, properties: p, items: i }) {
    return (  
      t === 'array'
        ? (i !== undefined ? keyValues(i) : [])
        : (p !== undefined ? Object.entries(p) : []));
  }

  for (const [key, value] of keyValues(x)) {    
    output += recordStr(key, value, level);  
    proc(value, level + 1);

    if (details[key] !== undefined)
      output += `\n> ${details[key]}\n`;
  } 
}

proc(
  JSON.parse(fs.readFileSync('../../xpm.schema.json', 'utf8'))
);

output +=
`
# Examples

## Two-phase image
\`\`\`json
${
  stripJsonComments(
    fs.readFileSync('example1.json', 'utf8')
  ).replace(/^\s*$\n/gm, '')
}
\`\`\`
`;

fs.writeFileSync('xpm.schema.md', output);

console.log('done.');


// console.log(md_file);