# Digital Predistorter with Real-Valued Feedback Employing Forward Model Estimation

Matlab scripts which perform DPD with real-valued feedback
(only real part or imaginary part of IQ demodulator output).
Firstly, the forward PA model is estimated and subsequently it is used
to train DPD coefficients using the indirect learning architecture.

The scripts supplement a research paper Digital Predistorter with Real-Valued
Feedback Employing Forward Model Estimation by Jan KRAL, Tomas GOTTHANS,
Roman MARSALEK, and Michal HARVANEK.

[RUN_ANALYSIS_01.m](RUN_ANALYSIS_01.m) - main script, runs all the analysis presented
  in the paper. During the analysis, the main charts are plot. 

[Run_Diff_DPD_Arch.m](Run_Diff_DPD_Arch.m) - script which executes different DPD 
  architectures including the proposed technique using the forward model extraction
  and ILA

[PLOT_01.m](PLOT_01.m) - plots results in the form presented in the paper

# Citation
If used for research, plese cite our paper
