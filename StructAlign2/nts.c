/* #include <stdio.h>
#include <string.h>
#include <stdlib.h> */

unsigned int convertNT(char *res, char *nt)
{
  char result;
  result = '?';
  if ( strcmp(res, "02I") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "08Q") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "08T") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "0AD") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "0C") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "0DC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "0DG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "0DT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "0G") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "0KL") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "0KX") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "0KZ") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "0L3") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "0L4") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "0L5") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "0L6") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "0L7") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "0OH") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "0OJ") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "0R5") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "0R6") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "0R8") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "0U") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "0UH") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "10C") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "125") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "126") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "127") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "12A") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "18M") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "18Q") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "1AP") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "1CC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "1FC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "1FZ") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "1GC") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "1MA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "1MG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "1RN") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "1RT") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "1SC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "1TL") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "1TW") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "23G") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "23T") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "2AD") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "2AR") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "2AT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "2AU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "2BA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "2BD") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "2BP") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "2BT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "2BU") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "2DA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "2DT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "2EG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "2FE") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "2FI") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "2GT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "2IA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "2L8") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "2LA") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "2LF") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "2MA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "2MG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "2MU") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "2NT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "2OT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "2PR") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "2ST") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "365") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "3AD") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "3AT") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "3AU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "3AY") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "3D1") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "3DA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "3ME") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "3TD") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "40A") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "40C") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "40G") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "40T") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "47C") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "4BD") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "4DG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "4DU") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "4OC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "4PC") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "4SC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "4SU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "574") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "5AA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "5AT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "5BU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "5CF") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "5CG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "5CM") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "5FC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "5FU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "5GP") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "5HC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "5HT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "5HU") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "5IC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "5IU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "5MC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "5MU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "5NC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "5OC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "5PC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "5PY") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "5SE") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "5SI") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "5UA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "63G") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "63H") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "64P") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "64T") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "6AP") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "6CF") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "6GO") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "6GU") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "6HA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "6HC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "6HG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "6HT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "6IA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "6MA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "6MP") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "6MZ") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "6OG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "6PO") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "70U") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "7DA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "7GU") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "7MG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "84T") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "8AG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "8AN") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "8BA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "8DG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "8FG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "8MG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "8OG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "9DG") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "A1P") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "A23") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "A2F") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "A2L") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "A2M") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "A38") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "A3P") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "A40") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "A43") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "A44") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "A47") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "A5L") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "A5M") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "A5O") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "A66") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "A6A") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "A6C") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "A6G") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "A6U") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "A9Z") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "ABR") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "ABS") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "ACP") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "AD2") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "ADE") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "ADI") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "ADK") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "ADP") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "AET") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "AF2") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "AG9") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "AGD") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "AGS") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "AMO") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "AMP") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "ANG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "ANP") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "ANZ") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "AP7") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "APC") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "APN") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "AS") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "AT7") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "ATD") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "ATL") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "ATM") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "ATP") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "AVC") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "AZT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "B7C") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "BA2") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "BGM") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "BGR") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "BLS") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "BOE") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "BRG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "BRU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "C2E") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "C2L") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C2S") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C31") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C34") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C36") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C37") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C38") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C42") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C43") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C45") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C46") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C49") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C4S") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C5P") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C66") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "C6G") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "CAR") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CBR") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CBV") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CCC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CFL") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CFZ") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CG1") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "CH") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CH1") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CM0") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "CMP") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "CMR") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CP1") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CPN") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CSG") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CSL") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CSM") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "CTG") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "CTP") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CUD") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CX2") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "CYT") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "D3N") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "D3T") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "D5M") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "DAD") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "DCM") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "DCP") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "DCT") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "DCZ") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "DDG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "DDS") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "DDY") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "DG3") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "DG8") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "DGP") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "DGT") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "DI") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "DJF") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "DNR") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "DOC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "DTP") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "DU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "DUP") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "DUR") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "DUT") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "DUZ") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "DX4") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "DZ4") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "DZM") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "E") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "E1X") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "EDA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "EDC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "EEM") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "EFG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "EHG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "EIT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "EPE") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "F2A") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "F3A") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "F3H") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "F4H") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "F5H") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "F6H") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "FA2") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "FAG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "FAX") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "FDG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "FHA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "FHG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "FHU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "FMG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "FMU") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "FOX") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "FYA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "G1C") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "G1M") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "G2C") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "G2L") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "G2M") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "G2P") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "G2S") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "G31") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "G36") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "G38") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "G46") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "G47") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "G48") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "G49") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "G7M") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GAO") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GBR") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GCK") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "GCP") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GDO") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GDP") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GF2") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GFC") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GFF") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GFH") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GFL") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GFM") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GGH") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GH3") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GMP") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GMS") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GMU") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "GN7") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GNE") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GNG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GNP") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GOM") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "GPN") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GRB") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GRC") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GS") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GSR") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GSS") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GSU") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "GTF") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "GTP") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GUN") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "GX1") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "H2U") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "HEU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "HGL") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "HMU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "HN0") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "HN1") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "HPA") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "HXB") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "HXZ") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "I") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "IC") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "IG") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "IGU") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "IMC") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "IPN") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "IU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "JDT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "KAG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "KIR") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "L8P") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "LCA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "LCC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "LCG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "LGP") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "LHC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "LKC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "LLP") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "LMS") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "M1G") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "M2G") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "M5M") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "MA6") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "MA7") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "MAD") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "MAU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "MCY") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "MDJ") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "MDK") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "MDQ") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "MDU") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "MDV") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "ME6") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "MFT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "MG1") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "MGT") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "MIA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "MMT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "MNU") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "MRG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "MSP") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "MTU") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "N2G") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "N5C") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "N5M") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "N6G") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "N6M") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "N79") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "NEA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "NMS") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "NMT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "NYM") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "O2C") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "O2G") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "OGX") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "OHU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "OMC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "OMG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "OMU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "ONE") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "OXG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "P") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "P2T") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "P2U") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "P5P") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "PBT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "PDU") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "PGN") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "PGP") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "PLR") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "PPU") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "PPW") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "PPZ") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "PQ0") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "PQ1") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "PRF") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "PRN") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "PSD") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "PST") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "PSU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "PU") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "PUY") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "PVX") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "PYO") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "QBT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "QSI") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "QUO") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "R") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "RCE") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "RIA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "RMP") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "RPC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "RSP") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "RSQ") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "RUS") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "S2M") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "S4A") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "S4C") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "S4G") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "S6G") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "SAH") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "SAM") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "SC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "SDG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "SFG") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "SMP") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "SMT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "SOS") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "SPT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "SRA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "SSA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "SSJ") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "SSU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "SUR") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "T23") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "T2S") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "T2T") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "T32") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "T38") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "T39") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "T3P") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "T41") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "T48") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "T49") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "T4S") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "T5O") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "T5S") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "T66") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "T6A") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "TAF") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TCP") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TCY") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "TDR") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "TDY") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TED") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TEP") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "TFE") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TFF") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TFO") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "TFT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TGP") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "THM") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "THP") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "THX") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "TLB") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TLC") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TLN") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TM2") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TMP") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TNV") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "TP1") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TPC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "TPN") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TSB") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "TSP") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TT") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TTD") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TTE") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TTI") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "TTM") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TTP") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "TX2") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "U2L") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "U31") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "U33") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "U34") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "U36") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "U37") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "U3H") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "U5P") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "U8U") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UAR") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UBB") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UBD") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UBI") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UCL") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UD5") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UDP") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UF2") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UFP") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UFR") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "UFT") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UMP") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UMS") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UMX") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UPC") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UPE") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UPG") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UPS") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UPV") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UR3") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "URA") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "URI") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "URX") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "US1") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "US2") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "US3") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "US4") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "US5") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "USM") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UTP") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UVX") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "UZR") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "VAA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "X") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "XAD") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "XAE") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "XAL") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "XAN") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "XAR") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "XCL") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "XCR") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "XCT") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "XCY") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "XG4") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "XGA") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "XGL") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "XGR") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "XGU") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "XJS") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "XTF") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "XTH") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "XTL") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "XTR") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "XUA") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "XUG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "Y") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "YCO") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "YG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "YMP") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "YYG") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "Z") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "ZAD") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "ZAN") == 0 ) {
  	result = 'A';}
  if ( strcmp(res, "ZBC") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "ZBU") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "ZCY") == 0 ) {
  	result = 'C';}
  if ( strcmp(res, "ZDU") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "ZGU") == 0 ) {
  	result = 'G';}
  if ( strcmp(res, "ZHP") == 0 ) {
  	result = 'U';}
  if ( strcmp(res, "ZP4") == 0 ) {
  	result = 'T';}
  if ( strcmp(res, "ZTH") == 0 ) {
  	result = 'T';}

(*nt) = result;
if (result != '?') {
	return 0;}
else {
	return 1;}
}

/* void main()
{
char a;
char res[4];
unsigned int result;
result = convertNT("08Q", &a);
printf("%u\n%c\n", result, a);

a = 'g';
strcpy(res, "ZAN");
if ( strcmp(res, "ZAN")==0 )
{
	a = 'T';
}
printf("%c", a);
} */
