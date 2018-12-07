\name{VIDTAXA}
\alias{VIDTAXA}
\encoding{latin1}
\title{IDENTIFICATION OF TAXA BASED ON MORPHOLOGICAL VARIABILITY}
\description{Identification of the different taxa based on the morphological variability observed in a Principal Components Analysis or a Correspondence Analysis.} 
\usage{
VIDTAXA(data, var, labels, cat=NULL, analysis="PCA", por=80, k=NULL,
pthreshold=0.05, ellipse=FALSE, convex=FALSE, dim=c(1,2), size=c(1,5),
showCluster=TRUE, VIF=FALSE, threshold=10, method="overlap", minimum=TRUE,
ResetPAR=TRUE, PAR=NULL, PCA=NULL, SCATTERPLOT=NULL, HCLUST=NULL,
CLUSTER=NULL, BOXPLOT=NULL, mfrowBOXPLOT=NULL, LabelCat=NULL, COLOR=NULL,
COLORC=NULL, COLORB=NULL, PCH=NULL, XLIM=NULL, YLIM=NULL, XLAB=NULL, YLAB=NULL,
ylabBOXPLOT=NULL, LEGEND=NULL, MTEXT= NULL, TEXTvar=NULL, TEXTlabels=NULL,
arrows=TRUE, larrow=0.7, colArrows="black", quadratic=FALSE, file1="Output.txt",
file2="Cat loadings.csv", file3="Descriptive statistics of clusters.csv",
file4="Original data and cluster number.csv", file5="Var loadings-Linear.csv",
file6="Cat loadings-Linear.csv", file7="Table cross-validation-Linear.csv",
file8="Cases cross-validation-Linear.csv", file9="Table cross-validation-Quadratic.csv",
file10="Cases cross-validation-Quadratic.csv", file11="Plots VARSEDIG.pdf",
file12="U Mann-Whitney test.csv", na="NA", dec=",", row.names=TRUE)
}
\arguments{
  \item{data}{Data file.}
  \item{var}{Variables that are included in the analysis.}
  \item{labels}{Variable that allows to display a label for each case.}
  \item{cat}{Optionally, it is possible to specify a variable to show a grouping in the plot of the Principal Components or Correspondence analyses.}
  \item{analysis}{If it is "PCA" a Principal Components analysis is carried out, whereas a Correspondence analysis is performed if the selection is "CA".}
  \item{por}{Cut-off threshold specifying the cumulative variance percentage, to determine how many axes are selected from the Principal Components or Correspondence analyses. By default it is 80\%, which means that the axes are selected until reaching an accumulated variance percentage of 80\%.}
  \item{k}{Number of clusters in which the Dendrogram is divided. If it is NULL, the algorithm select automatically the maximum number of clusters in which the Dendrogram can be divided, which are those groups that are statistically different in at least one variable according to the U Mann-Whitney test.}
  \item{pthreshold}{Threshold probability of the U Mann-Whitney test.}
\item{ellipse}{If it is TRUE, the ellipses with the levels of significance to the 0.5 (inner ellipse) and 0.95 (outer ellipse) of each category of the variable \emph{cat} are depicted. These levels of significance can be modified by entering the function \link[car]{scatterplot} using the argument \emph{SCATTERPLOT} and modifying the argument \emph{levels=c(0.5,0.95)}. If it is TRUE, the ellipses of the clusters in the Discriminant analysis and in the polar coordinate plot of the VARSEDIG algorithm are also calculated.}
  \item{convex}{If it is TRUE, the convex hull is calculated for each category in the plot of the Principal Components or Correspondence analyses, but only if some variable has been selected in the argument \emph{cat}. If TRUE, the convex hull of the clusters is also calculated in the Discriminant analysis and in the polar coordinate plot of the VARSEDIG algorithm.}
  \item{dim}{Vector with two values indicating the axes that are shown in the plot of the Principal Components or Correspondence analyses.}
  \item{size}{Size range of bubbles. Two values: minimum and maximum size.}
  \item{showCluster}{If it is TRUE, the number of each cluster is shown in the Dendrogram.}
  \item{VIF}{If it is TRUE, the inflation factor of the variance (VIF) is used to select the highly correlated variables and, therefore, not correlated variables are excluded from the Principal Components analysis.}
  \item{threshold}{Cut-off value for the VIF.}
  \item{method}{Three different methods for prioritizing the variables according to their capacity for discrimination can be used in the VARSEDIG algorithm. If the method is "overlap", a density curve is obtained for each variable and the overlap of the area under the curve between the two groups of the variable \emph{group} is estimated for all variables. Those variables with lower overlap should have better discrimination capacities and, hence, all variables are ordered from lowest to highest overlap; in other words, from the highest to lowest discrimination capacity. If the method is "Monte-Carlo", a Monte-Carlo test is performed comparing all values of group 1 with group 2, and all values of group 2 with 1. The variables are prioritized from the variable with the lowest mean of all p-values (highest discrimination capacity) to the variable with the highest mean of all p-values (lowest discrimination capacity). If the method is "logistic regression", then a binomial logistic regression is calculated and only significant variables are selected for further analyses with the regression performed by steps using the Akaike Information Criterion (AIC).}
  \item{minimum}{If it is TRUE, the algorithm is designed to find a significant discrimination between both groups with the minimum possible number of significant variables. Therefore, only the variables with higher discrimination capacity are selected. If it is FALSE, the algorithm selects all significant variables, and not only those with higher discrimination capacity. This argument is only valid with the methods "Monte-Carlo" and "overlap" and it is useful in those cases that discrimination between the groups is difficult and requires to include as many as variables as possible.}
   \item{ResetPAR}{If it is FALSE, the default condition of the function PAR are not placed and those defined by the user on previous graphics are maintained.}
  \item{PAR}{It accesses the PAR function that allows to modify many different aspects of the graphs.}
  \item{PCA}{It accesses the \link[stats]{prcomp} function of the stats package.}
  \item{SCATTERPLOT}{It accesses the function \link[car]{scatterplot} of the car package.}
  \item{HCLUST}{You may access the function \link[stats]{hclust} of the stats package.}
  \item{CLUSTER}{Access to the function that allows to modify the graphic representation of the Dendrogram.}
  \item{BOXPLOT}{Allows to specify the characteristics of the boxplot.}
  \item{mfrowBOXPLOT}{It allows to specify the boxplot panel. It is a vector with two numbers, for example c(2,5) which means that the boxplots are put in 2 rows and 5 columns.}
  \item{LabelCat}{It allows to specify a vector with the names of the clusters in the boxplots. They must be as many as clusters.}
\item{COLOR}{It allows to modify the colours of the graphic in the in the plot of the Principal Components or Correspondence analyses, but they must be as many as different groups have the variable \emph{cat}.}
\item{COLORC}{It allows to modify the colours of the clusters in the Dendrogram, but they must be as many as clusters.}
\item{COLORB}{It allows to modify the colours of the clusters in the boxplots, but they must be as many as clusters.}
  \item{PCH}{Vector with the symbols in the plot of the Principal Components or Correspondence analyses, which must be as many as different groups have the variable \emph{cat}. If it is NULL they are calculated automatically starting with the symbol 15.}
  \item{XLIM, YLIM}{Vectors with the axes limits \emph{X} and \emph{Y} in the plot of the Principal Components or Correspondence analyses.}
  \item{XLAB, YLAB}{Legends of the axes \emph{X} and \emph{Y} in the plot of the Principal Components or Correspondence analyses.}
  \item{ylabBOXPLOT}{You can specify a vector with the legends of the axes \emph{Y} of the boxplots. They should be as many as the number of variables.}
  \item{LEGEND}{It allows to include or to modify a legend in the plot of the Principal Components or Correspondence analyses.}
  \item{MTEXT}{It allows to add text in the margins in the plot of the Principal Components or Correspondence analyses.}
  \item{TEXTvar}{It allows to modify the labels of the variables in the plot of the Principal Components or Correspondence analyses.}
  \item{TEXTlabels}{It allows to modify the labels of the cases in the plot of the Principal Components or Correspondence analyses plot.}
  \item{arrows}{If it is TRUE the arrows are shown in the scatterplot in the plot of the Principal Components or Correspondence analyses.}
  \item{larrow}{It modifies the length of the arrows in the plot of the Principal Components or Correspondence analyses.}
  \item{colArrows}{Colours of the arrows in the plot of the Principal Components or Correspondence analyses.}
  \item{quadratic}{If TRUE, a Quadratic Discriminant Analysis is performed, in addition to the Linear Discriminant Analysis.}
  \item{file1}{TXT FILE. Name of the output file with the results.}
  \item{file2}{CSV FILE. Name of the output file with the coordinates of the cases in the plot of the Principal Components or Correspondence analyses.}
  \item{file3}{CSV FILE. Name of the output file with the descriptive statistics of each variable for each of the clusters obtained in the Dendrogram.}
  \item{file4}{CSV FILE. Name of the output file with the original data of the variables and the cluster to which each case belongs.}
  \item{file5}{CSV FILE. Name of the output file with the coordinates of the variables in the Linear Discriminant Analysis plot.}
  \item{file6}{CSV FILE. Name of the output file with the coordinates of the categories in the Linear Discriminant Analysis plot.}
  \item{file7}{CSV FILE. Name of the output file with the prediction table using the cross-validation of the Linear Discriminant Analysis.}
  \item{file8}{CSV FILE. Name of the output file with the group to which each case belongs and the prediction of the Discriminant Analysis using the cross-validation of the Linear Discriminant Analysis.}
  \item{file9}{CSV FILE. Name of the output file with the predictions table using the cross-validation of the Quadratic Discriminant Analysis.}
  \item{file10}{CSV FILE. Name of the output file with the group to which each case belongs and the prediction of the Discriminant Analysis using the cross-validation of the Quadratic Discriminant Analysis.}
  \item{file11}{PDF File. Name of the output file with the graphics obtained from the VARSEDIG algorithm.}
  \item{file12}{CSV FILE. Name of the output file with the obtained probabilities of comparing all the variables among all the clusters with the U Mann-Whitney test.}
  \item{na}{CSV FILES. Text that is used in the cells without data.}
  \item{dec}{CSV FILES. It defines if a comma "," or a dot "." is used as decimal separator.}
  \item{row.names}{CSV FILES. Logical value that defines if identifiers are put in rows or a vector with a text for each of the rows.}
}
\value{
It is obtained:

1. A TXT file with the VIF (if the argument \emph{VIF=TRUE}), the correlations between variables, the Kaiser-Meyer-Olkin (KMO) test, the Bartlett sphericity test and the results of the Principal Components or Correspondence analyses The file is called by default "Output.TXT".

2. A CSV FILE with the coordinates for each case of the Principal Components or Correspondence analyses. The file is called by default "Cat loadings.CSV".

3. A CSV FILE with the descriptive statistics of each variable for each of the clusters obtained in the Dendrogram. The file is called by default "Descriptive statistics of clusters.CSV".

4. A CSV FILE with the original data of the variables and the cluster to which each case belongs. The file is called by default "Original data and cluster number.CSV".

5. A CSV FILE with the coordinates of the variables in the Linear Discriminant Analysis plot. The file is called by default "Var loadings-Linear.csv"

6. A CSV FILE with the coordinates of the categories in the Linear Discriminant Analysis plot. The file is called by default "Cat loadings-Linear.csv".

7. A CSV FILE with the predictions table using the cross-validation of Linear Discriminant Analysis. The file is called by default "Table cross-validation-Linear.csv".

8. A CSV FILE with the group to which each case belongs and the prediction of the Discriminant Analysis using the cross-validation of the Linear Discriminant Analysis. The file is called by default "Cases cross-validation-Linear.csv".

9. A CSV file with the predictions table using the cross-validation of the Quadratic Discriminant Analysis. The file is called by default "Table cross-validation-Quadratic.csv".

10. A CSV file with the group to which each case belongs and the prediction of the Discriminant Analysis using the cross-validation of the Quadratic Discriminant Analysis. The file is called by default "Cases cross-validation-Quadratic.csv".

11. A CSV file with the obtained probabilities of comparing all the variables among all the clusters with the U Mann-Whitney test.

12. A PDF file with the graphics obtained from the VARSEDIG algorithm.

13. A scatterplot of the Principal Components or Correspondence analyses.

14. A Dendrogram grouping by clusters according to the scores of the Principal Components or Correspondence analyses.

15. A graphic panel with a boxplot for each variable comparing the values of these variables between each of the clusters obtained in the Dendrogram.

16. A Graph of the Discriminant Analysis showing the influence of the variables on the discriminant axis I, differentiating the different clusters.

17. A graph of the Discriminant Analysis showing the scores of the discriminant axes I and II, differentiating the different clusters.

18. A bubble chart with the number of variables that are statistically different between clusters.
}
\details{
The aim of this analysis is to determine what statistically different groups are formed by applying a Principal Components or Correspondence analyses.

The first axis in a Principal Components analysis or Correspondence analysis is the linear combination of the original variables that has maximum variance. The second component is the linear combination of the original variables with maximum variance with the added condition that it is independent of the first (orthogonal), and so on, all the main components can be obtained, which, being independent of each other, contain different information. The independence or absence of correlation means that the new variables or components do not share common information. Each main component, therefore, explains the maximum possible residual variability (which has not already been explained above). Therefore, in a Principal Components or Correspondence analyses the cases are differentiated according to the variables that have greater variability. The idea of the analysis is to determine if statistically different groups are formed associated to the variability observed in the variables.

This analysis can be useful to find different groups when you really do not know what they are. For example, find different species using morphometric variables, without really knowing how many potential species there are and to what species each individual belongs. However, it is important to note that only different groups will be detected if the variables that have more variability give rise to different groups. It is possible that a variable does not present a great variability, but it is important for discriminating groups. This type of differentiation based on variables that do not have high variance, would not be detected in this analysis.

To detect the potential groups being formed, a Dendrogram is applied to the scores obtained from the axes that absorb a greater variance. By default, the axes that absorb 80\% of the variability are chosen, but this value can be modified by the user.

Subsequently, a Discriminant Analysis is carried out to determine if the clusters that have been generated are well discriminated, that is, to determine the number of correctly identified cases in each cluster.

Next, a U Mann-Whitney test is performed to determine if there are significant differences in the variables between the clusters.

Finally, the algorithm of the VARSEDIG function is applied (see for more details (Guisande, 2018). With this algorithm it is possible to determine if all the cases of each cluster are statistically different from the other clusters.

The idea of this function is to find the largest possible number of clusters with the highest discrimination percentage. To do this the user should perform tests, modifying the cut-off threshold by specifying the cumulative variance percentage to determine how many axes are selected from the Main Components (by default \emph{by=80}) and the variables to be included, eliminating those that are not correlated and are not useful in the Principal Components or Correspondence analyses, as well as those that have little discrimination power in the Discriminant Analysis.

\strong{FUNCTIONS}

The Correspondence analysis was performed with the \link[ca]{ca} function of the package ca (Greenacre & Pardo, 2006; Greenacre, 2007; Nenadic & Greenacre, 2007; Greenacre, 2013). The Principal Components Analysis was performed with the \link[stats]{prcomp} function of the stats package. The \link[usdm]{vif} function of the usdm package was used for the calculation of VIF (Naimi et al., 2014; Naimi, 2017). To perform the \emph{biplot} graph the \link[car]{scatterplot} function of the car package was used (Fox et al., 2018). The arrows are depicted with the function \link[IDPmisc]{Arrows} of the package IDPmisc (Locher & Ruckstuhl, 2014). The convex hull is estimated with the function \link[grDevices]{chull} of the package grDevices. KMO test was performed with the function \link[psych]{KMO} of the package psych (Revelle, 2018). Bartlett's test sphericity was performed with the function \link[REdaS]{bart_spher} of the package REdaS (Maier, 2015). The U Mann-Whitney test is performed with the \emph{wilcox.test} function of the base stats package. The comparison between clusters with the VARSEDIG algorithm is done with the VARSEDIM function of the VARSEDIG package (Guisande et al., 2016: Guisande, 2018). The Linear Discriminant Analysis was performed with the functions \link[candisc]{candisc} of the candisc package (Friendly, 2007; Friendly & Fox, 2017) and \link[MASS]{lda} of the MASS package (Venables & Ripley, 2002; Ripley et al., 2018). The Quadratic Discriminant Analysis was performed with the function \link[MASS]{qda} of the MASS package (Venables & Ripley, 2002; Ripley et al., 2018). The graph with one dimension in the Discriminant analysis was performed with the function \link[candisc]{plot.cancor} of the candisc package (Friendly, 2007; Friendly & Fox, 2017).

\strong{EXAMPLE}

The example consisted of analysing the morphometric variability of several species of scorpaeniformes. The aim is to find how many groups are statistically different based on the morphometric variability observed in the Principal Components analysis. For purposes only of graphic presentation in the Principal Components, the genus is used as a category \emph{cat="Genus"}. It is important to highlight that the category is not used for any statistical analysis and it is simply used to group the cases with ellipses or with the convex hull in the Principal Components graphic.

The analysis is performed by eliminating the variables that are not correlated, for which it is specified \emph{VIF=TRUE}. Therefore, the first result obtained is the VIF values of the variables. Those variables with a VIF lower than the threshold are no included in the Principal Components analysis.
\tabular{c}{
\if{html}{\figure{R1.jpg}{options: width=2.5in}}
\if{latex}{\figure{R1.jpg}{options: width=2.5in}}
}
The second statistic obtained is the KMO test, which tells us if the variables are adequate for the Principal Components. The value must be greater than 0.5. Therefore, all variables that do not have a value greater than 0.5, could be eliminated from the analysis. In the case that the value is exactly 0.5, it means that it is not possible to estimate the KMO.
\tabular{c}{
\if{html}{\figure{R2.jpg}{options: width=5.0in}}
\if{latex}{\figure{R2.jpg}{options: width=5.0in}}
}
The next statistic that appears is Bartlett's test of sphericity, which tests whether the correlation matrix is an identity matrix, which would indicate that the factor model is inappropriate. A value \emph{p} of the contrast smaller than the level of significance allows rejecting the hypothesis and concluding that there is correlation. Therefore, for the Principal Components analysis to be valid, the probability must be less than 0.05, as it is in this case.
\tabular{c}{
\if{html}{\figure{R3.jpg}{options: width=2.8in}}
\if{latex}{\figure{R3.jpg}{options: width=2.8in}}
}
Figure VIDTAXA.1 shows that the variability observed in the Principal Components analysis allows to clearly differentiate among the genera.
\tabular{c}{
\strong{Figure VIDTAXA.1.} Principal Components analysis showing the \cr
variability observed in the genera.\cr
\if{html}{\figure{F2.jpg}{options: width=5.5in}}
\if{latex}{\figure{F2.jpg}{options: width=5.5in}}
}
The first axis accounts for 54\%, the second for 25.3\% and the third for 8.5\% of the variance observed. The first three axes explain 87.8\% of the variance. Since the default value of \emph{by=80} was selected, these three Principal Component axes are selected.

Figure VIDTAXA.2 shows the Dendrogram where 6 clusters are grouped, which are the six  genera used in this example.
\tabular{c}{
\strong{Figure VIDTAXA.2.} Dendrogram with the scores of the axes selected\cr
from the Principal Components analysis.\cr
\if{html}{\figure{F3.jpg}{options: width=5.0in}}
\if{latex}{\figure{F3.jpg}{options: width=5.0in}}
}
Figure VIDTAXA.3 shows the differences between clusters for each of the variables. It is clear, for instance, the difference in M21 for cluster 1, in M6 for cluster 5, etc.
\tabular{c}{
\strong{Figure VIDTAXA.3.} Boxplot obtained for each of the variables\cr
with the averaged values for each cluster.\cr
\if{html}{\figure{F4.jpg}{options: width=5.2in}}
\if{latex}{\figure{F4.jpg}{options: width=5.2in}}
}
The Discriminant Analysis shows that it is possible to correctly discriminate 100\% of cases by cross-validation with the Linear method. The first discriminant axis explains most of the variability and discriminates well between the 6 clusters (Figure VIDTAXA.4). Many variables seem to be important for the discrimination since the arrows are not small. Figure VIDTAXA.5 shows the first two discriminant axes and shows the differences between the 6 clusters.
\tabular{c}{
\strong{Figure VIDTAXA.4.} Axis I of the Discriminant analysis\cr
\if{html}{\figure{F5.jpg}{options: width=3.7in}}
\if{latex}{\figure{F5.jpg}{options: width=3.7in}}
}
\tabular{c}{
\strong{Figure VIDTAXA.5.} Axes I and II of the Discriminant analysis\cr
\if{html}{\figure{F6.jpg}{options: width=3.7in}}
\if{latex}{\figure{F6.jpg}{options: width=3.7in}}
}
The next test to determine if the clusters are statistically different was the comparison of the variables between the clusters. The results of the U Mann-Whitney test are shown in Figure VIDTAXA.6. For clusters to be different, there must be at least one statistically different variable when comparing each cluster with all the others. In the graph it is noted that in the comparison between all the clusters there is always a point, that is, there is always at least one variable that is different. In fact, between cluster 2 and cluster 4, the smaller number of statistically different variables was observed, a total of 14 variables. Therefore, from the comparison of the variables between clusters with the U Mann-Whitney test, it is concluded that the clusters are statistically different from each other.
\tabular{c}{
\strong{Figure VIDTAXA.6.} Plot where the bubbles show the number of variables,\cr
that are statistically different (p <= 0.05) between clusters.\cr
\if{html}{\figure{F7.jpg}{options: width=4.3in}}
\if{latex}{\figure{F7.jpg}{options: width=4.3in}}
}
Finally, in a pdf, the plots obtained from applying the VARSEDIG algorithm are saved, whose objective is to compare all the clusters with each other.

Figure VIDTAXA.7 shows the example of the comparison of cluster 1 with 2. It is observed that the variable that discriminate significantly between both clusters is M22 (upper right panel). The Monte-Carlo test showed that the individuals that most resembles cluster 2 in cluster 1 (lower left panel) does not have significant differences in the polar coordinate axes X and Y (p = 0.1).

The individual that most resembles cluster 1 to cluster 2 (bottom right panel), it is very close to the significance threshold on both the polar coordinate axes X and Y (p = 0.077). Therefore, it cannot be concluded that cluster 1 and 2 are different. The same process would be done to compare the rest of the clusters.
\tabular{c}{
\strong{Figure VIDTAXA.8.} Plots obtained from the algorithm VARSEDIG.\cr
It is shown the comparison between the cluster 1 and 2.\cr
\if{html}{\figure{F8.jpg}{options: width=5.5in}}
\if{latex}{\figure{F8.jpg}{options: width=5.5in}}
}
Therefore, according to the Discriminant Analysis and the tests performed with the U Mann-Whitney test, the clusters are statistically different from each other, but the VARSEDIG algorithm showed that not all clusters are statistically different. However, it is very important to emphasize that the VARSEDIG algorithm considers two statistically different groups if the case that most resembles each group is statistically different using the Monte-Carlo test. The Monte-Carlo test needs a large number of cases in each group for detecting significant differences. That is, it is possible that, as it was shown in the comparison of cluster 1 with cluster 2, the cases of both groups that resemble each other are not within the point cloud of the other group, but due to the low number of cases in each group, it is not possible to determine that the difference is not due to chance.
}
\references{
Fox, J., Weisberg, S., Adler, D., Bates, D., Baud-Bovy, G., Ellison, S., Firth, D., Friendly, M., Gorjanc, G., Graves, S., Heiberger, R., Laboissiere, R., Monette, G., Murdoch, D., Nilsson, H., Ogle, D., Ripley, B., Venables, W. & Zeileis, A. (2018) Companion to Applied Regression. R package version 3.0-0. Available at: \url{https://CRAN.R-project.org/package=car}.

Friendly, M. & Fox, J. (2017) Visualizing Generalized Canonical Discriminant and Canonical Correlation Analysis. R package version 0.8-0. Available at: \url{https://CRAN.R-project.org/package=candisc}.

Friendly, M. (2007). HE plots for Multivariate General Linear Models. \emph{Journal of Computational and Graphical Statistics}, 16: 421-444.

Greenacre, M. (2007) \emph{Correspondence Analysis in Practice}. Second Edition. London: Chapman & Hall / CRC.

Greenacre, M. (2013).  Simple, Multiple and Joint Correspondence Analysis. R package version 0.53. Available at: \url{https://CRAN.R-project.org/package=ca}.

Greenacre, M.J. & Pardo, R. (2006) Subset correspondence analysis: visualizing relationships among a selected set of response categories from a questionnaire survey. \emph{Sociological Methods and Research}, 35: 193-218.

Guisande, C. (2018) An Algorithm for Morphometric Characters Selection and Statistical Validation in Morphological Taxonomy. R package version 1.8. Available at: \url{https://CRAN.R-project.org/package=VARSEDIG}.

Locher, R. & Ruckstuhl, A. (2014) Utilities of Institute of Data Analyses and Process Design. R package version 1.1.17. Available at: \url{https://CRAN.R-project.org/package=IDPmisc}.

Maier, M.J. (2015) Companion Package to the Book 'R: Einfuehrung durch angewandte Statistik. R package version 0.9.3. Available at: \url{https://CRAN.R-project.org/package=REdaS}.

Naimi, B. (2017) Uncertainty analysis for species distribution models. R package version 1.1-18. Available at: \url{https://CRAN.R-project.org/package=usdm}.

Naimi, B., Hamm, N.A.S., Groen, T.A., Skidmore, A.K., & Toxopeus, A.G. (2014) Where is positional uncertainty a problem for species distribution modelling? \emph{Ecography}, 37: 191-203.

Nenadic, O. & Greenacre, M. (2007) Correspondence analysis in R, with two- and three-dimensional graphics: The ca package. \emph{Journal of Statistical Software}, 20: 1-13.

Revelle,W. (2018) Procedures for Psychological, Psychometric, and Personality Research. R package version 1.8.4. Available at: \url{https://CRAN.R-project.org/package=psych}.

Ripley, B., Venables, B., Bates, D.M., Hornik, K., Gebhardt, A. & Firth, D. (2018) Support Functions and Datasets for Venables and Ripley's MASS. R package version 7.3-50. Available at: \url{https://CRAN.R-project.org/package=MASS}.

Rizopoulos, D. (2006) ltm: An R package for latent variable modelling and item response theory analyses. \emph{Journal of Statistical Software}, 17: 1-25.

Rizopoulos, D. (2018) Latent Trait Models under IRT. R package version 1.1-1. Available at: \url{https://CRAN.R-project.org/package=ltm}.

Venables, W.N. & Ripley, B.D. (2002) \emph{Modern Applied Statistics with S}. Springer, fourth edition, New York. \url{https://www.stats.ox.ac.uk/pub/MASS4}.
}
\examples{
\dontrun{

data(scorpaeniformes)

VIDTAXA(data=scorpaeniformes, var=c("M2","M3","M4","M5","M6","M7",
"M8","M9","M10","M11","M12","M13","M14","M15","M16","M19","M20",
"M21","M22","M23","M24","M25","M26","M27"), labels="Genus",
cat="Genus", VIF=TRUE, convex=TRUE)
}
}
\keyword{VIDTAXA}
