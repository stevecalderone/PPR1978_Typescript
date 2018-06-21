/*********************************************************************************************
MIT License

Copyright (c) [2018] [Steve G. Calderone]

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
'*********************************************************************************************/

"use strict";

/**
* Validates that the required pure component data named ranges are available.
*
* @customfunction
*/
function validateWorkbook(){
  
  var PDataExists  = false
  var PDataPropertiesExist = false
  var PDataPropertyNamesExist = false

  var myWorkbook = SpreadsheetApp.getActiveSpreadsheet()
  var myRanges = myWorkbook.getNamedRanges()
  var numberOfSheets = SpreadsheetApp.getActiveSpreadsheet().getNumSheets()
  var mySheets = SpreadsheetApp.getActiveSpreadsheet().getSheets()
  
  for(var i = 0; i < myRanges.length; i++){
    if(myRanges[i].getName() === "PData_Properties"){
       PDataPropertiesExist = true
    }
    if(myRanges[i].getName() === "PData_PropertyNames"){
       PDataPropertyNamesExist = true
    }
  }
    if(PDataPropertiesExist && PDataPropertyNamesExist){
      return true
    }else{
      return false
    }
}
/*********************************************************************/
/**
* Calculates constant pressure heat capacity (j/mol/K) given T (C), P (pbara), mole amounts (kg-moles), phase (vapor or liquid) and a dataSet.
*
* @customfunction
*/
function createDatasetHeader(){
  try{

    var LiquidPhaseExists = false
    var int_Species = 0
    var int_MW = 0
    var int_TC = 0
    var int_PC = 0
    var int_TMN1 = 0
    var int_Omega = 0
    var int_ZC = 0
    var int_TB = 0
    var int_HVap = 0
    var int_Hf298 = 0
    var int_Gf298 = 0
    var int_S298 = 0
    var int_HHV = 0
    var int_LHV = 0
    var int_CpDataType = 0
    var NumberOfDatasetRows = 0
    var pDataArray = []
    var myErrorMsg = ""
    var LiquidSpeciesFound = false
    var LiquidIndex  = 0
    var fcnName = "CreateDataset"
    var Phase = "vapor"
    var SfNotFound = false
    var GfNotFound = false
    var columnHeaders = []

    pDataArray = SpreadsheetApp.getActiveSpreadsheet().getRangeByName("PData_Properties").getValues()
    
    fcnName = "CreateDataset"
    
    if(validateWorkbook === false){
      myErrorMsg = "Pdata worksheet, or Cell A1:A1 in PData worksheet does not contain 'PData' or Pdata named ranges 'PData_Properties' or 'PData_PropertyNames' do not exist"
      return myErrorMsg
    }
    
    columnHeaders.push("(MW)", "(TC, K)", "(PC, bara)", "(OMEGA)", "(ZC)", "(TB, K)", "(Hvap, kJ/kg-mole)", "(Hf298, kJ/g-mole)", "(S298, J/g-mole/K)",
                        "(Gf298, kJ/g-mole)", "(CPData)", "(NIST-TMN1)", "(NIST-TMX1)", "(NIST-A1)", "(NIST-B1)", "(NIST-C1)", "(NIST-D1)",
                        "(NIST-E1)", "(NIST-F1)", "(NIST-G1)", "(NIST-H1)", "(NIST-TMN2)", "(NIST-TMX2)", "(NIST-A2)", "(NIST-B2)", "(NIST-C2)", "(NIST-D2)",
                        "(NIST-E2)", "(NIST-F2)", "(NIST-G2)", "(NIST-H2)", "(NIST-TMN3)", "(NIST-TMX3)", "(NIST-A3)", "(NIST-B3)", "(NIST-C3)", "(NIST-D3)", "(NIST-E3)",
                        "(NIST-F3)", "(NIST-G3)", "(NIST-H3)", "(NIST-TMN4)", "(NIST-TMX4)", "(NIST-A4)", "(NIST-B4)", "(NIST-C4)", "(NIST-D4)", "(NIST-E4)", "(NIST-F4)",
                        "(NIST-G4)", "(NIST-H4)", "(NIST-TMN5)", "(NIST-TMX5)", "(NIST-A5)", "(NIST-B5)", "(NIST-C5)", "(NIST-D5)", "(NIST-E5)", "(NIST-F5)", "(NIST-G5)",
                        "(NIST-H5)", "(NIST-TMN6)", "(NIST-TMX6)", "(NIST-A6)", "(NIST-B6)", "(NIST-C6)", "(NIST-D6)", "(NIST-E6)", "(NIST-F6)", "(NIST-G6)", "(NIST-H6)");


    for(var i = 0; i < columnHeaders.length; i++){
      for(var j = 0; j < pDataArray.length; j++){
        if(columnHeaders[i] === pDataArray[0][j]){
          columnHeaders[i] = "found"
        }
      }
    }
    
    var missingHeaders = []
    for(var i = 0; i < columnHeaders.length; i++){
      if(columnHeaders[i] !== "found"){
       missingHeaders.push(columnHeaders[i]) 
      }
    }
    
    if(missingHeaders.length !== 0){
      return ["The Following column headers are missing from the PData worksheet: ", missingHeaders]
    }

    
    columnHeaders = []
    
    columnHeaders.push("(MW)", "(TC, K)", "(PC, bara)", "(OMEGA)", "(ZC)", "(TB, K)", "(Hvap, kJ/kg-mole)", "(Hf298, kJ/g-mole)", "(S298, J/g-mole/K)",
                        "(Gf298, kJ/g-mole)", "(CPData)", "(NIST-TMN1)", "(NIST-TMX1)", "(NIST-A1)", "(NIST-B1)", "(NIST-C1)", "(NIST-D1)",
                        "(NIST-E1)", "(NIST-F1)", "(NIST-G1)", "(NIST-H1)", "(NIST-TMN2)", "(NIST-TMX2)", "(NIST-A2)", "(NIST-B2)", "(NIST-C2)", "(NIST-D2)",
                        "(NIST-E2)", "(NIST-F2)", "(NIST-G2)", "(NIST-H2)", "(NIST-TMN3)", "(NIST-TMX3)", "(NIST-A3)", "(NIST-B3)", "(NIST-C3)", "(NIST-D3)", "(NIST-E3)",
                        "(NIST-F3)", "(NIST-G3)", "(NIST-H3)", "(NIST-TMN4)", "(NIST-TMX4)", "(NIST-A4)", "(NIST-B4)", "(NIST-C4)", "(NIST-D4)", "(NIST-E4)", "(NIST-F4)",
                        "(NIST-G4)", "(NIST-H4)", "(NIST-TMN5)", "(NIST-TMX5)", "(NIST-A5)", "(NIST-B5)", "(NIST-C5)", "(NIST-D5)", "(NIST-E5)", "(NIST-F5)", "(NIST-G5)",
                        "(NIST-H5)", "(NIST-TMN6)", "(NIST-TMX6)", "(NIST-A6)", "(NIST-B6)", "(NIST-C6)", "(NIST-D6)", "(NIST-E6)", "(NIST-F6)", "(NIST-G6)", "(NIST-H6)");
    
    var columnHeaderIndexes = []
        for(var i = 0; i < columnHeaders.length; i++){
      for(var j = 0; j < pDataArray[0].length; j++){
        if(columnHeaders[i].toString() === pDataArray[0][j].toString()){
         columnHeaderIndexes.push(j)
        }
      }
    }
    
    if(columnHeaderIndexes.length !== columnHeaders.length){
     return "There is something wrong with the column headers" 
    }
    
    var omegaIndex = columnHeaderIndexes[columnHeaders.indexOf("(OMEGA)")]
    var outputArray = []
    columnHeaderIndexes.splice(columnHeaders.indexOf("(ZC)") + 1,0,"k(i)","b(i)")
    columnHeaders.splice(columnHeaders.indexOf("(ZC)") + 1,0,"k(i)","b(i)")
    outputArray.push(columnHeaderIndexes)
    outputArray.push(columnHeaders)
    return outputArray
    
      }catch(myErrorHandler){
    return myErrorHandler.message
  }

}
/***********************************************************************************/
/**
* Returns a 2 line header for the specie decomposition.
*
* @customfunction
*/
function createDecompHeader(){
  try{
    var numOfPredictiveGroups = 31
    var columnHeaders = []
    var pDataArray = SpreadsheetApp.getActiveSpreadsheet().getRangeByName("PData_Properties").getValues();
    var group1Index = pDataArray[0].indexOf("(GROUP 1)");

    
    var columnHeaderIndexes = ["Indexes"]
    for(var j = group1Index; j < numOfPredictiveGroups + group1Index; j++){
      columnHeaderIndexes.push(j)
    }
    columnHeaders.push(columnHeaderIndexes)
    
    columnHeaders.push(["Groups","CH3","CH2","CH","C","CH4","C2H6","CHaro","Caro","C, aro-fused rings","CH2, cyclic","CH/C, cyclic","CO2","N2","H2S","SH","H2O","C2H4","CH2/CH, alkenic","C, alkenic","CH/C, cycloalkenic","H2","CO","He","Ar","SO2","O2","NO","COS","NH3","NO2/N2O4","N2O"],
                       ["Groups Names","(GROUP 1)","(GROUP 2)","(GROUP 3)","(GROUP 4)","(GROUP 5)","(GROUP 6)","(GROUP 7)","(GROUP 8)","(GROUP 9)","(GROUP 10)","(GROUP 11)","(GROUP 12)","(GROUP 13)","(GROUP 14)","(GROUP 15)","(GROUP 16)","(GROUP 17)","(GROUP 18)","(GROUP 19)","(GROUP 20)","(GROUP 21)","(GROUP 22)","(GROUP 23)","(GROUP 24)","(GROUP 25)","(GROUP 26)","(GROUP 27)","(GROUP 28)","(GROUP 29)","(GROUP 30)","(GROUP 31)"])
    
    return columnHeaders
    
  }catch(myErrorMessage){
  return myErrorMessage.message
  }
}
/***************************************************************************/
/**
* Adds a species to the decompoosition.
*
* @customfunction
*/
function addSpeciesToDecomp(Species, decompHeader){
  try{
    var outputArray = []  ;
    var numOfPredictiveGroups = 31;
    var sum = 0;
    var group1Index = decompHeader[0][1];
    
    var pDataArray = SpreadsheetApp.getActiveSpreadsheet().getRangeByName("PData_Properties").getValues();
    
    for(var i = group1Index; i < numOfPredictiveGroups + group1Index; i++){
      if(typeof(Species[0][0].valueOf()) === "number"){
        if(pDataArray[Species[0][0]][i].toString === "" || typeof(pDataArray[Species[0][0]][i]) === null){        
          outputArray.push(0)
        }else{
          if(Number(pDataArray[Species[0][0]][i]) === NaN){
            myErrorMsg = "Some decomposition values are not numeric.";
            throw new Error("Function name: " + fcnName + ", " + myErrorMsg);
          }else{ 
            outputArray.push(Number(pDataArray[Species[0][0]][i]));
            sum = sum + Number(pDataArray[Species[0][0]][i])
          }
        }
        
      }else{
        outputArray.push(0)
      }
    }
    
    if(sum !== 0){
      for(var i = 0; i < numOfPredictiveGroups; i++){            
        outputArray[i] = outputArray[i]/sum;
      }
    };
    
    return outputArray;
  }catch(myErrorHandler){
    return myErrorHandler.message;
  }
}

/***************************************************************************/
/**
* Validates that the species are found in the PData_Properties named range and returns an array of the row indexes for the species.
*
* @customfunction
*/
function validateSpecies(Species){
  try{
    var liquidsExist = false
    var liquidsIndex = 0
    var numberOfSpecies = 0
    var speciesArray = []
    var pDataArray = []
    var speciesIndexes = []
    var myErrorMsg = ""
    var fcnName = "validateSpecies"
    var sizeArray = []
    var curSpecies = ""
    var speciesFound = false
    var arLength = 0
    
    if(!validateWorkbook()){
      myErrorMsg = "There is something wrong with the PData workbook"
      throw new Error("Function name: " + fcnName + ", " + myErrorMsg);
    }
    
    pDataArray = SpreadsheetApp.getActiveSpreadsheet().getRangeByName("PData_Properties").getValues()
    speciesArray = Species.slice(0)
    sizeArray = size(speciesArray)
    
    if (typeof (Species) === "string" ) {            
      arLength = 1
    }else{
      arLength = speciesArray.length
    }
      
    
    for(var i = 0; i < arLength; i++){
      speciesFound = false
      if (sizeArray[i][0 /* arrayTest */] === true && sizeArray[i][1 /* valueORlength1 */] === 1) {
        
        curSpecies = speciesArray[i].toString();
      }else{
        if(sizeArray[i][1 /* arrayTest */] !== false){
          if (sizeArray[i][2 /* arrayTest */] === false) {            
            curSpecies = speciesArray.toString();
          }
        }else{
          return "The species parameter must be a single dimensional range." 
          throw new Error("Function name: " + fcnName + ", " + myErrorMsg);
        }
        
      }
      for(var j = 0; j < pDataArray.length; j++){
        if(curSpecies.toString() === pDataArray[j][0].toString()){
          speciesIndexes.push(j)
          speciesFound = true
          j = pDataArray.length - 1
        }
      }
      if(speciesFound === false){
        speciesIndexes.push("not found")
      }
    }
    return speciesIndexes
  }catch(myErrorHandler){
    return myErrorHandler.message
  }
}
/********************************************************************************/
/**
* Validates that the required pure component data named ranges are available.
*
* @customfunction
*/
function addSpeciesToDataset(Species, dataSetHeader, Phase, errorMsgOn){
  if (Phase === void 0) { Phase = "vapor"; }
  if (errorMsgOn === void 0) { errorMsgOn = "false"; }
  try{
    var outputArray = []  ;
    var omegaIndex = dataSetHeader[1].indexOf("(OMEGA)");
    var kiIndex = dataSetHeader[1].indexOf("k(i)");
    var biIndex = dataSetHeader[1].indexOf("b(i)");
    var tcIndex = dataSetHeader[1].indexOf("(TC, K)");
    var pcIndex = dataSetHeader[1].indexOf("(PC, bara)");
    var phase = Phase.toString().toLowerCase() ;
    
    var pDataArray = SpreadsheetApp.getActiveSpreadsheet().getRangeByName("PData_Properties").getValues();
    
    var omega = pDataArray[Species[0][0]][dataSetHeader[0][omegaIndex]];
    var tc = pDataArray[Species[0][0]][dataSetHeader[0][tcIndex]];
    var pc = pDataArray[Species[0][0]][dataSetHeader[0][pcIndex]];
    
    for(var i = 0; i < dataSetHeader[0].length; i++){            
      if(i === 1){
        outputArray.push(pDataArray[Species[0][0]][dataSetHeader[0][i]]);
      }else if(i === kiIndex && phase === "vapor"){
        if(omega < 0.491){
          outputArray.push(0.37464 + 1.54226 * omega - 0.26992 * Math.pow(omega, 2));
        }else{
          outputArray.push(0.379642 + 1.487503 * omega - 0.164423 * Math.pow(omega, 2) + 0.016666 * Math.pow(omega, 3)); //'<= m or k for PR1978
        }
      }else if(i === biIndex && phase === "vapor"){
        outputArray.push(0.0778 * gasLawR * tc / pc);      
      }else if(i > 0 && i < 12 /* cpData */ && i !== biIndex && i !== kiIndex && phase === "vapor"){
        outputArray.push(pDataArray[Species[0][0]][dataSetHeader[0][i]]);
      }else if(i > 0 && i < 12 /* cpData */ && i !== biIndex && i !== kiIndex && phase === "liquid"){
        outputArray.push(""); 
      }else{
        outputArray.push(pDataArray[Species[0][0]][dataSetHeader[0][i]]);
      }
      
    };
  
    return outputArray;
    
  }catch(myErrorHandler){
    return myErrorHandler.message;
  }
  
  
}

  


