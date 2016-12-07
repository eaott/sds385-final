# sds385-final
This is the final project for @jgscott's SDS 385 course for myself and @raghavshroff.
The main write-up can be found in `final-write-up.pdf`.

## Data
There are two sets of data for this project, `master.csv` and `master_countedByKegg.csv`.
These are the "master" and "grouped" data sets in the final paper.

## Code
The code we used for analysis can be found in `final_master.Rmd` and `final_kegg.Rmd`,
which are RMarkdown documents containing both R code and markdown text. The code
is effectively the same, save some constants, and is separated more for our
own convenience rather than long-term use. Each of these is associated
with a `.nb.html` and `.pdf` file that contain the formatted output.

We also used @jgscott's SGD logistic regression code as it better handled the intercept
compared to @eaott's code from earlier in the course which shrunk the intercept toward zero.
