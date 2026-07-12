@echo off
cd /d d:\repos\FCSimple
C:\Users\joshu\miniconda3\envs\fcview\Lib\R\bin\Rscript.exe -e "rmarkdown::render('inst/analysis/fcview_sccomp_guide.Rmd', output_format='pdf_document', output_dir='inst/analysis')" > inst\analysis\build_log.txt 2>&1
echo Exit code: %ERRORLEVEL% >> inst\analysis\build_log.txt