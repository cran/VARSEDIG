\author{Cástor Guisande González, Universidad de Vigo, Spain.}
\name{VARSEDIM}
\alias{VARSEDIM}
\encoding{latin1}
\title{Variable selection to discriminate many taxonomic groups}
\description{This function performs an algorithm for morphometric characters selection and statistical validation in morphological taxonomy among many taxonomic groups.}
\usage{
VARSEDIM(data, variables, group, method="overlap", stepwise=TRUE,
VARSEDIG=TRUE, minimum=TRUE, kernel="gaussian", cor=TRUE, ellipse=FALSE,
convex=TRUE, file="Plots VARSEDIG.pdf", na="NA", dec=",", row.names=FALSE)
}
\arguments{
  \item{data}{Data file.}
  \item{variables}{Variables to be selected.}
  \item{group}{Variable with the groups to be discriminated.}
  \item{method}{Three different methods for prioritizing the variables according to their capacity for discrimination can be used. If the method is "overlap", a density curve is obtained for each variable and the overlap of the area under the curve between the two groups of the variable \emph{group} is estimated for all variables. Those variables with lower overlap should have better discrimination capacities and, hence, all variables are ordered from lowest to highest overlap; in other words from the highest to lowest discrimination capacity. If the method is "Monte-Carlo", a Monte-Carlo test is performed comparing all values of group 1 with group 2, and all values of group 2 with 1. The variables are prioritized from the variable with the lowest mean of all p-values (highest discrimination capacity) to the variable with the highest mean of all p-values (lowest discrimination capacity). If the method is "logistic regression", then a binomial logistic regression is calculated and if the argument stepwise=TRUE (default option), then only significant variables are selected for further analyses with the regression performed by steps using the Akaike Information Criterion (AIC).}
  \item{stepwise}{If TRUE, the logistic regression is applied by steps, in order to eliminate those variables that are not significant. The Akaike information criterion (\emph{AIC}) is used to define what are the variables that are excluded.}
  \item{VARSEDIG}{If it is TRUE, the variables are added for the estimation of polar coordinates in the priority order according to the method "overlap", "Monte-Carlo", or "logistic regression" and the variable is selected if it significantly contributes to discriminate between both groups. See details section for further information.}
  \item{minimum}{If it is TRUE, the algorithm is designed to find a significant discrimination between both groups with the minimum possible number of significant variables. Therefore, only the variables with higher discrimination capacity are selected. If it is FALSE, the algorithm selects all significant variables,and not only those with higher discrimination capacity. This argument is only valid with the methods "Monte-Carlo" and "overlap" and it is useful in those cases that discrimination between the groups is difficult and requires to include as many as variables as possible.}
  \item{kernel}{A character string giving the smoothing kernel to be used for estimating the overlap of the area under the curve between groups. This must be one of "gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine" or "optcosine". For further details about the estimation of the density curve see the details section of the function \link[stats]{density} of base stats package.}
  \item{cor}{If it is TRUE the variables are ordered according to the correlation between them when estimating the polar coordinates. Therefore, the next variable to another variable is the one that has a greater positive correlation.}
\item{ellipse}{If it is TRUE the ellipses with the levels of significance to the 0.5 (inner ellipse) and 0.95 (outer ellipse) of each category of the variable \emph{group} is depicted. These levels of significance can be modified by entering the function \link[car]{scatterplot} using the argument \emph{SCATTERPLOT} and modifying the argument \emph{levels=c(0.5,0.95)}.}
  \item{convex}{If it is TRUE the convex hull is depicted for each category.}
  \item{file}{PDF FILE. Filename with the plots of the function \link[VARSEDIG]{VARSEDIG}.}
  \item{na}{CSV FILE. Text that is used in the cells without data.}
  \item{dec}{CSV FILE. It defines if the comma "," is used as decimal separator or the dot ".".}
  \item{row.names}{CSV FILE. Logical value that defines if identifiers are put in rows or a vector with a text for each of the rows.}
}
\value{It is obtained a PDF file with the plots of the function \link[VARSEDIG]{VARSEDIG}.}
\details{
The difference with the function \link[VARSEDIG]{VARSEDIG} is that all the different taxa of the variable group are compared with each other, instead of just comparing two taxa. It uses the same algorithm described in the function \link[VARSEDIG]{VARSEDIG}.
}
\examples{
\dontrun{
data(characiformes)
VARSEDIM(data=characiformes, variables= c("M2", "M3", "M4",  "M5", "M6",
"M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16", "M17",
"M18", "M19", "M20", "M21", "M22", "M23", "M24", "M25", "M26", "M27", "M28"),
group="Genus")
}
}
\keyword{VARSEDIM}
