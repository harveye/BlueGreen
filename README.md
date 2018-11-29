# Blue-Green project 

## Code and R data repository for the Harvey et al., manuscript

### Repository structure:
- `data`: Raw data 
- `scripts`: R scripts (**see below for a description**)
- `R`: Custom R functions called by the main code from the `scripts` folder
- `doc`: Documents that are not generated by the R code 
- `figs`: Figures generated by the code
- `output`: Documents generated by the code other than figures (tables, clean data files)

### Scripts order: 
- 1- BG_DatMan.R 
   - input: Raw data
   - output: Clean data for analyses ('cl.data.RDS' in `output` folder)
- 2- BG_data_processing_SVM_final.R (here for information purpose only - cannot be run)
   - input: Data from video counts (not present in this repository)
   - output: Counts by species (the script runs the SVM model to ID species)

### Markdown files: 
- Overview.Rmd : A notebook documenting all exploratory analyses (please note that LRR calculations in this notebook are wrong and were not fixed - see Final_analysis.Rmd for final LRR calculation and figures)
- Final_analysis.Rmd : Document only the analytical pipeline from the article with final analyses (RDA, LRR, mixed-effect models)

