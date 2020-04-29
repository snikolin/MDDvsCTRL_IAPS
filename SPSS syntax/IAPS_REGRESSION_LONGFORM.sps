* Encoding: UTF-8.
*** Longform dataset ***

USE ALL.
COMPUTE filter_$=(group=0 & outlier = 1 & valence = 0).
VARIABLE LABELS filter_$ 'group=0 & outlier = 1 & valence = 0 (FILTER)'.
VALUE LABELS filter_$ 0 'Not Selected' 1 'Selected'.
FORMATS filter_$ (f1.0).
FILTER BY filter_$.
EXECUTE.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LPP
  /METHOD=ENTER strategy_2x2.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LPP
  /METHOD=ENTER DASS_anx.

*** Using a significant time window ***

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LPP_sig
  /METHOD=ENTER strategy_2x2 .

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LPP_sig
  /METHOD=ENTER  DASS_anx.

*** Analysis on complete dataset ***

USE ALL. 

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LPP
  /METHOD=ENTER group valence strategy_2x2.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LPP
  /METHOD=ENTER group valence DASS_anx.

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LPP_sig
  /METHOD=ENTER group valence strategy_2x2 .

REGRESSION
  /MISSING LISTWISE
  /STATISTICS COEFF OUTS R ANOVA
  /CRITERIA=PIN(.05) POUT(.10)
  /NOORIGIN 
  /DEPENDENT LPP_sig
  /METHOD=ENTER group valence DASS_anx.

*** MRMM method testing ***

MIXED LPP BY group valence strategy_2x2
  /FIXED = group valence group*valence strategy_2x2 | SSTYPE(3)
  /METHOD = ML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= valence | SUBJECT(PID) COVTYPE(AD1).

MIXED LPP BY group valence DASS_anx
  /FIXED = group valence group*valence DASS_anx | SSTYPE(3)
  /METHOD = ML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= valence | SUBJECT(PID) COVTYPE(AD1).


MIXED LPP_sig BY group valence strategy_2x2
  /FIXED = group valence group*valence strategy_2x2 | SSTYPE(3)
  /METHOD = ML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= valence | SUBJECT(PID) COVTYPE(AD1).

MIXED LPP_sig BY group valence DASS_anx
  /FIXED = group valence group*valence DASS_anx | SSTYPE(3)
  /METHOD = ML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= valence | SUBJECT(PID) COVTYPE(AD1).

