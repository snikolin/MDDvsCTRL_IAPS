* Encoding: UTF-8.

FILTER BY outlier.

*** Build MRMMs ***

MIXED LPP BY group valence
  /FIXED = group valence group*valence | SSTYPE(3)
  /METHOD = ML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= valence | SUBJECT(PID) COVTYPE(AD1).

MIXED LPP_sig BY group valence
  /FIXED = group valence group*valence | SSTYPE(3)
  /METHOD = ML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= valence | SUBJECT(PID) COVTYPE(AD1).

* UN or AR1

*** Run MRMMs ***

MIXED LPP BY group valence
  /FIXED = group valence group*valence | SSTYPE(3)
  /METHOD = REML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= valence | SUBJECT(PID) COVTYPE(UN) 
  /EMMEANS = TABLES(group*valence) COMPARE(group) ADJ(LSD).
*  /SAVE RESID (RESID_LPP).

MIXED LPP_sig BY group valence
  /FIXED = group valence group*valence | SSTYPE(3)
  /METHOD = REML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= valence | SUBJECT(PID) COVTYPE(UN) 
  /EMMEANS = TABLES(group*valence) COMPARE(group) ADJ(LSD).
*  /SAVE RESID (RESID_LPPsig).

*** Add Strategy to MRMM ***

MIXED LPP BY group valence strategy_2x2
  /FIXED = group valence group*valence strategy_2x2 | SSTYPE(3)
  /METHOD = REML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= valence | SUBJECT(PID) COVTYPE(UN) 
  /EMMEANS = TABLES(group*valence) COMPARE(group) ADJ(LSD)
  /SAVE RESID (RESID_LPP_strat).

MIXED LPP_sig BY group valence strategy_2x2
  /FIXED = group valence group*valence strategy_2x2 | SSTYPE(3)
  /METHOD = REML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= valence | SUBJECT(PID) COVTYPE(UN) 
  /EMMEANS = TABLES(group*valence) COMPARE(group) ADJ(LSD)
  /SAVE RESID (RESID_LPPsig_strat).

*** Add Anxiety to MRMM ***

MIXED LPP BY group valence DASS_anx
  /FIXED = group valence group*valence DASS_anx | SSTYPE(3)
  /METHOD = REML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= valence | SUBJECT(PID) COVTYPE(UN) 
  /EMMEANS = TABLES(group*valence) COMPARE(group) ADJ(LSD)
  /SAVE RESID (RESID_LPP_anx).

MIXED LPP_sig BY group valence DASS_anx
  /FIXED = group valence group*valence DASS_anx | SSTYPE(3)
  /METHOD = REML
  /PRINT = G R SOLUTION TESTCOV
  /REPEATED= valence | SUBJECT(PID) COVTYPE(UN) 
  /EMMEANS = TABLES(group*valence) COMPARE(group) ADJ(LSD)
  /SAVE RESID (RESID_LPPsig_anx).


*** Examine residuals ***

EXAMINE VARIABLES = RESID_LPP BY group
  /ID=PID
  /PLOT BOXPLOT STEMLEAF HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES = RESID_LPPsig BY group
  /ID=PID
  /PLOT BOXPLOT STEMLEAF HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES = RESID_LPP_strat BY group
  /ID=PID
  /PLOT BOXPLOT STEMLEAF HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES = RESID_LPPsig_strat BY group
  /ID=PID
  /PLOT BOXPLOT STEMLEAF HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES = RESID_LPP_anx BY group
  /ID=PID
  /PLOT BOXPLOT STEMLEAF HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.

EXAMINE VARIABLES = RESID_LPPsig_anx BY group
  /ID=PID
  /PLOT BOXPLOT STEMLEAF HISTOGRAM
  /COMPARE GROUPS
  /STATISTICS DESCRIPTIVES
  /CINTERVAL 95
  /MISSING LISTWISE
  /NOTOTAL.
