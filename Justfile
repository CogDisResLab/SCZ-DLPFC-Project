alias a: all

render:
  quarto render cell-level-run1.Rmd
  quarto render cell-level-run2.Rmd
  quarto render cell-level-run3.Rmd
  quarto render cell-level-run4.Rmd

uka:
  Rscript uka_analysis_single_sample.R

creeden:
  Rscript creedenzymatic_analysis.R
  Rscript generate_quartile_plots.R

all: render uka creeden
