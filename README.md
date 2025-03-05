# TesCA-hormones
Clinical and biological implications of differential expression of sex hormone-related genes in testicular cancer

## Results

Differential expression and co-regulation of genes related to sex hormone metabolism in testicular cancers as well as classification of testicular tumors by differing expression was done in two publicly available cohorts: TCGA [^1] and GSE99420 [^2]. 
The resulting hormonal clusters of cancers were thoroughly characterized in terms of clinics, tumor biology, and predicted therapy response.

![result_summary](https://github.com/PiotrTymoszuk/TesCA-hormones/assets/80723424/55cedc63-5b12-4576-a563-3068fc1bf314)


Manuscript parts (Result draft, Figures, Tables, and Method) are available as [a Word file](https://github.com/PiotrTymoszuk/TesCA-hormones/blob/main/paper/figures_tables_methods.docx).
Supplementary Material is provided as [a Word file](https://github.com/PiotrTymoszuk/TesCA-hormones/blob/main/paper/supplementary_material.docx) and [an Excel file with Supplementary Tables](https://github.com/PiotrTymoszuk/TesCA-hormones/blob/main/paper/supplementary_tables.xlsx).

You may also access [Figures](https://github.com/PiotrTymoszuk/TesCA-hormones/tree/main/paper/figures) and [Supplementary Figures](https://github.com/PiotrTymoszuk/TesCA-hormones/tree/main/paper/supplementary%20figures) as high-resolution PDF files.

## Usage

The analysis pipeline requires some development packages, the easiest way to install them is to use `devtools`:

```r

devtools::install_github('PiotrTymoszuk/trafo')
devtools::install_github('PiotrTymoszuk/clustTools')
devtools::install_github('PiotrTymoszuk/soucer')
devtools::install_github('PiotrTymoszuk/figur')
devtools::install_github('PiotrTymoszuk/microViz')
devtools::install_github('PiotrTymoszuk/coxExtensions')
devtools::install_github('PiotrTymoszuk/htGLMNET')
devtools::install_github('PiotrTymoszuk/ExDA')

```
To launch the entire pipeline, source the `exec.R` file:

```r

source('exec.R')

```

## Terms of use

To reference and use analysis results, please cite our GitHub repository and [our publication (DOI: 10.1186/s12610-025-00254-5)](https://pubmed.ncbi.nlm.nih.gov/40011822/). 

## Contact

Data and analysis requests should be addressed to [Dr. Renate Pichler](mailto:renate.pichler@i-med.ac.at). [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com) is the maintainer of the repository.

## References

[^1]: Liu J, Lichtenberg T, Hoadley KA, Poisson LM, Lazar AJ, Cherniack AD, Kovatich AJ, Benz CC, Levine DA, Lee A V., et al. An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. Cell (2018) 173:400-416.e11. doi:10.1016/J.CELL.2018.02.052

[^2]: Lewin J, Soltan Ghoraie L, Bedard PL, Hamilton RJ, Chung P, Moore M, Jewett MAS, Anson-Cartwright L, Virtanen C, Winegarden N, et al. Gene expression signatures prognostic for relapse in stage I testicular germ cell tumours. BJU Int (2018) 122:814–822. doi:10.1111/BJU.14372
