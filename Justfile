alias a: all

render:
  quarto render case_control_E1_pairwise.Rmd
  quarto render case_control_E2_pairwise.Rmd
  quarto render case_control_E3_pairwise.Rmd
  quarto render case_control_E4_pairwise.Rmd

uka:
  Rscript uka_analysis.R

creeden:
  Rscript creedenzymatic_analysis.R
  Rscript generate_quartile_plots.R

all: render uka creeden