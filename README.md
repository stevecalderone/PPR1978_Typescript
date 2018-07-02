# PPR1978_Typescript
'*********************************************************************************************
MIT License

Copyright (c) 2018 Steve G. Calderone

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'*********************************************************************************************
This code is implemented and shared in a read only Google Sheet here:
https://docs.google.com/spreadsheets/d/1HAIj2jcMhk3-gon62P7JgCQ3EdnSpwoLPlI1hyzopug/edit?usp=sharing

Some PPR1978 examples are included.

To use this sheet make a copy (File -> Make a copy...). The sheet can then be edited.

To use this VSCode developed code in Google Sheets follow these steps:

This code (eosCode.js) was developed in VSCode and needs to be modified for use in Google Sheets. To use this code in Google Sheets preform a find and replace in eosCode.js to delete the following lines or sections of code:

1 - 
In eosCode.js delete all instances of 'exports.* = *;'
for example:
'exports.phaseCp = phaseCp;'

2
Delete the following code located at the top of eosCode.js
'Object.defineProperty(exports, "__esModule", { value: true });'
'var findCubicRoots_js_1 = require("./findCubicRoots.js");'
'var Math = require("mathjs");'

3
In eosCode.js delete two instances of 'findCubicRoots_js_1.'

Change
Z = findCubicRoots_js_1.getCubicRoot(1, (b - 1), a - 3 * Math.pow(b, 2) - 2 * b, (-a * b + Math.pow(b, 2) + Math.pow(b, 3)), phase);
to
Z = getCubicRoot(1, (b - 1), a - 3 * Math.pow(b, 2) - 2 * b, (-a * b + Math.pow(b, 2) + Math.pow(b, 3)), phase);

and 
Z = findCubicRoots_js_1.getCubicRoot(1, (b - 1), a - 3 * Math.pow(b, 2) - 2 * b, (-a * b + Math.pow(b, 2) + Math.pow(b, 3)), phase);
to
Z = getCubicRoot(1, (b - 1), a - 3 * Math.pow(b, 2) - 2 * b, (-a * b + Math.pow(b, 2) + Math.pow(b, 3)), phase);

4 - In findCubicRoots.js delete:
Object.defineProperty(exports, "__esModule", { value: true });
//import { errorSub } from "./errorCode";

Within Google Sheets script copy then paste this revised code (eosCode.js & findCubicRoots.js) into the code.gs tab or a new tab or tabs.

The gsheetsCode.js was developed in Google sheets and does not require any revisions. Paste the gsheetsCode.js code into code.gs or a new tab.

Save code.gs or the new tab(s) and begin calculations.

References

1
Peng Robinson Equation of State
A New Two-Constant Equation of State
Ding-Yu Peng, Donald B. Robinson
Ind. Eng. Chem. Fundamen., 1976, 15 (1), pp 59–64
DOI:  10 0.1021 / i160057a011
Publication Date: February 1976

2
Spreadsheet for Thermodynamics Instruction
Phillip Savage - University of Michigan
'ChE classroom'
Fall 1995
http://ufdcimages.uflib.ufl.edu/AA/00/00/03/83/00128/AA00000383_00128_00262.pdf

3
Flash Routine Reference:
CHEMENG 120 taught by Professor Musgrave during the Spring '04 term at Stanford.
http://documentslide.com/documents/lecture-5-isothermal-flash-calculations.html

4
Implementation of Departure Functions on worksheets Departure_Pt1, 2, 3 & 4
Dr. Phillip Savage - Penn State University
http://www.che.psu.edu/department/directory-detail.aspx?LandOn=Gen&q=pes15

5
Bubble and Dew Point Routines adapted from:
It’s not as easy as it looks - Revisiting Peng–Robinson equation of state convergence issues for dew point, bubble point and flash alculations
Vamshi Krishna Kandula,(a) John C. Telotte (b) and F. Carl Knopf (corresponding author) (a)
(a) Chemical Engineering Department, Louisiana State University, USA
e -mail: cknopf@ southalabama.edu
(b) Chemical Engineering Department, Florida A&M University – Florida State University, USA

International Journal of Mechanical Engineering Education, Volume 41, Number 3 (July 2013), © Manchester University Press
 http://journals.sagepub.com/doi/pdf/10.7227/IJMEE.41.3.2

6
Cubic Equation VBA Code
Dr.Tomas b.Co
Michigan Tech University
https://www.mtu.edu/chemical/department/faculty/co/

7
Bicubic Interpolation Code (for LeeKeslerZ() function)
https://mathformeremortals.wordpress.com/

8
modArraySupport module and array handling techniques
Chip Pearson, chip@cpearson.com, www.cpearson.com

9
Derivation of the enthalpy departure function
https://shareok.org/bitstream/handle/11244/12606/Thesis-1996-R231e.pdf?sequence=1
by Abhishek Rastogi
APPENDIX a
A DETAILED DERIVATION OF PENG-ROBINSON EQUATION OF STATE ENTHALPY DEPARTURE FUNCTION

10
PData worksheet physical properties adapted from:
Properties Databank 1.0
Pedro Fajardo
ppfk@ yahoo.com
http://www.cheresources.com/invision/files/file/125-physical-properties-ms-excel-add-in/

11
Thermodynamic Properties Involving Derivatives
Using the Peng-Robinson Equation of State
R.M. Pratt Ph. D.
The National University of Malaysia
(now Associate Professor of Mathematics and Sciences, Fresno Pacific University)
http://ufdcimages.uflib.ufl.edu/AA/00/00/03/83/00150/AA00000383_00150_00112.pdf

12
Chemical Equilibrium by Gibbs Energy Minimization on Spreadsheets
Int. J. Engng Ed. Vol. 16, No. 4, pp. 335±339, 2000
Y.LWIN
Department of Chemical Engineering, Rangoon Institute of Technology, Insein P. O., Rangoon, Burma.
e -mail: ylwin@ yahoo.com
http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.476.5931&rep=rep1&type=pdf

13
Getting a Handle on Advanced Cubic Equations of State
www.cepmagazine.org, November 2002, CEP
Chorng H. Twu, Wayne D. Sim and Vince Tassone, Aspen Technology, Inc.
http://people.clarkson.edu/~wwilcox/Design/adv-ceos.pdf

14
VLE predictions with the Peng.Robinson equation of state and
temperature dependent kij calculated through a group contribution method
Jean-Noel Jaubert., Fabrice Mutelet
Laboratoire de Thermodynamique des Milieux Polyphases, Institut National Polytechnique de Lorraine, Ecole Nationale Superieure des Industries
Chimiques, 1rue Grandville, 54000 Nancy, France
Received 24 January 2004; accepted 25 June 2004

15
Prediction of Thermodynamic Properties of Alkyne-Containing
Mixtures with the E.PPR78 Model
Xiaochun Xu, Jean-Noe.l Jaubert, Romain Privat, and Philippe Arpentinierö
Ecole Nationale Supe.rieure des Industries Chimiques, Laboratoire Re.actions et Ge.nie des Proce.de.s (UMR CNRS 7274),
Universite. de Lorraine, 1 rue Grandville, 54000 Nancy, France
Centre de Recherche Paris Saclay, Air Liquide, 1 chemin de la porte des loges, BP 126, 78354 Jouy-en-Josas, France

16
Efficient flash calculations for chemical process
design — extension of the Boston–Britt
‘‘Inside–out’’ flash algorithm to extreme
conditions and new flash types
Vipul S. Parekh and Paul M. Mathias
Computers Chem. Engng Vol. 22, No. 10, pp. 1371—1380, 1998
( 1998 Published by Elsevier Science Ltd.
All rights reserved. Printed in Great Britain
