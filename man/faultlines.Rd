\name{faultlines}
\alias{faultlines}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Diversity Faultline Measures
}
\description{
Wrapper function for calculating various diversity faultline measures as proposed by Thatcher,  Jehn,  and 
  Zanutto (2003),  by Gibson and Vermeulen (2003),  by Shaw (2004),  by Trezzini (2008),  by 
  Bezrukova,  Jehn,  and Zanutto (2009),  by van Knippenberg,  Dawson,  West,  and Homan (2011) 
  and by Lawrence and Zyphur (2011).  The function can also calculate the cluster-based ASW measure
  proposed by Meyer and Glenz (2012).
}
\usage{
faultlines(data, 
           group.par = NA, 
           attr.type = NA, 
           attr.weight = NA, 
           rescale = NA, 
           method = "asw", 
           metric = "euclid", 
           maxgroups = 6, 
           i.level = FALSE,
           cores = 1,
           usesghomo = FALSE,
           by.attr = FALSE,
           ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data frame containing the diversity attributes of individual team members from one or more teams. All attributes in this data frame are used for calculating diversity faultline measure(s). Individuals are listed in rows and diversity attributes (dimensions) are listed in named columns. If the data frame contains more than one team, it must contain a variable denoting team membership (i.e., a team number for each individual) in one column. The name of ths column has to be passed to the \code{group.par} parameter.
}
  \item{group.par}{Specifies the column (by name) of \code{data} that contains a numeric team-identifier in the cases where the dataset contains data from multiple teams.
}
  \item{attr.type}{Vector of strings specifying the attribute type (\code{"numeric"} or \code{"nominal"}) for each column of \code{data}, except for the column denoted by \code{group.par} if the data set contains more than one team. Note that the methods \code{"trezzini"} and \code{"shaw"} require all attributes to be \code{"nominal"}.
}
  \item{attr.weight}{Vector of numeric weights, each usually in the range of 0 and 1, specifying each attribute's relative weight for each column of \code{data}, except for the column denoted by \code{group.par} if the data set contains more than one team. These weights are used with the methods \code{"thatcher"}, \code{"bezrukova"}, and \code{"asw"} only, and ignored for all other methods. See 'Details'.
}
  \item{rescale}{Either one of \code{NA} (the default) or \code{"std"}.  If \code{rescale = "std"} is selected,  all attributes of type \code{"numeric"} (specified by the \code{attr.type}-parameter),  are divided by the sample's standard-deviation.
}
  \item{method}{Chracter string specifying the diversity faultline measure to calculate.  This must be one of \code{"asw"},   \code{"thatcher"},   \code{"shaw"},   \code{"bezrukova"},  \code{"trezzini"},  
  	 \code{"knippenberg"},   \code{"lcca"} or \code{"gibson"}.
}
  \item{metric}{String denoting the employed metric (Euclidean or Mahalanobis) for calculating distances between individuals.  Must be either \code{"euclid"} or \code{"mahal"}.  See 'Details'. 
}
  \item{maxgroups}{Maximum number of subgroups that ASW can detect for each team. The default is 6. Only employed if \code{method = "asw"}, see 'Details'.
}
  \item{i.level}{Logical, indicates if faultlines should be computed on the individual level instead of on team-level. This is implemented for the "asw"-method only.
}
  \item{cores}{calculate faultlines for multiple groups in parallel, using the \code{mclapply}-functino from the \code{parallels} package. Specify the number of cores to be used.
}
  \item{usesghomo}{Experimental, logical. Specifies whether outgroup homogeneity should be takten into account.
}
  \item{by.attr}{Logical. If set to \code{TRUE}, faultline calculation is based on the weighted average of normed single-attribute metacontrast ratios.
}
  \item{...}{Other parameters.
}
}
\details{
	The \code{faultlines()} function is the main function of the package and is intended for calculating the faultlines; most other functions in the package are workhorse and sub functions that are called by \code{faultlines()}.
	
	The function is run over a data set (a data frame passed to the function with the parameter \code{data}) containing the members of one or more teams as rows and their diversity attributes that are used for calculating a given faultline measure as columns. Note that all columns will be used for calculating the faultline. Thus, this will most likely be a subset of the 'full' data frame that was gathered in a given research. If the data set contains more than one team, it must contain a column that specifies a team number for each team member, thus indicating team membership. The name of this column must be passed to the argument \code{group.par}.
	
	For each diversity attribute contained in the data frame, the user must specify its scale (either numeric or nominal) in the same order of the variables in the data frame as a character vector and pass it to the function with the parameter \code{attr.type}. For example, if the data set contains the variables age (numeric in years), ethnicity (character factor), and gender (character factor), this parameter must be specified as \code{attr.type = c("numeric", "nominal", "nominal")}. Note that the faultline measures proposed by Shaw (2004) and Trezzini (2008) require that all attributes are nominal. Thus, prior to calculating diversity faultline strengths with these two methods, the user must recode numeric attributes such as age to factors with levels such as 'young', 'middle-aged', and 'old' and specify the attribute type of this variable as nominal.
	
	For faultline measures \code{"thatcher"},  \code{"bezrukova"},  or \code{"asw"} that are capable of dealing with numeric attributes such as age or tenure, an optional weightign factor for each diversity attribute can bre specified using the \code{attr.weight} parameter.  These weights indicate how strong a difference of 1 (in case of numeric attributes) or a different category (in case of nominal attributes) is fractured into the faultline.  In the example case of age (in years),  gender,  and ethnicity,  specifying this parameter as \code{attr.weight = c(0.1, 1, 1)} means that an age difference of ten years is equally weighted as a difference in gender,  which is equally weighted as a difference in ethnicity.  Note that these are the default values for Thatcher's et al. (2003) Fau that are probably used in most papers,  but these appear to be arbitrary.  More research is required with regard to the choice of these weights in a given context.
	
	The \code{metric} parameter lets the user specify whether Euclidean or Mahalanobis distances should be employed in determining how different team members are from each other.  This metric is only employed by the methods \code{"thatcher"},  \code{"bezrukova"},  and \code{"asw"}.  Note that the former two methods (Bezrukova et al., 2009; Thatcher et al., 2003) were introduced based on Euclidean distances,  which assume that diversity attributes are uncorrelated.  (Meyer et al., 2012) showed that correlations between diversity attributes (e.g.,  between age and tenure) can have a significant influence on diversity faultline measures.  They thus suggested to employ Mahalanobis distances to control for such correlations.  They explicitly included this option in the calculation of the ASW measure,  but invoking it for Thatcher's et al. Fau or for Bezrukova's Faultline Distance measure is purely experimental.  Employing Mahalanobis distances for the latter two measures will thus deliver a measure that has not been described in the literature.  Furthermore, calculating Mahalanobis distances requires an inversion of the variance-covariance-martix of attributes.  Using Mahalanobis metrics is therefore restricted to data sets with invertible variance-/covariance matrices,  i.e., to numeric attributes only. 
	
	For calculating faultlines on the individual level, specify \code{i.level = TRUE}.  This will cause the program to find faultlines that maximize each team-member's individual silhouette-width,  instead of maximizing the team's average silhouette-width.  Thus,  in a team with n members, n individual faultlines will be calculated.  In addition to an n-sized vector of faultline strengths,  the result includes an n x n adjacency matrix,  indicating in each cell the portion of faultlines that assign the corresponding pair of team-members to the same subgroup.
}

\value{
A data frame with one row for each team with the following columns:
\item{team}{The number of the given team}

\item{teamsize}{The size of the team, i.e., the number of team members}

\item{fl.value}{The diversity faultline strength (or in the case of \code{method  =  "bezruvkova"},  faultline strength * distance) for the group as identified by the employed measure.}

\item{mbr_to_subgroups}{The member-to-subgroup-association.  A character vector of numbers seperated by commas.  These are the team members' subgroup numbers,  in the order of the original data files.}

\item{number_of_subgroups}{The number of homogeneous subgroups inside the team as identified by the faultline measure.  The measures \code{"thatcher"} and \code{"bezrukova"} alsways assume a two-subgroup solution.  The measures by \code{"gibson"},  \code{"shaw"},  \code{"trezzini"} and \code{"knippenberg"} do not identify subgroups and hence return \code{NA}.  Only \code{"asw"} will return a value larger than two.}

\item{adjm}{...}

\item{fl.ind}{...}
}
\references{
Bezrukova, K., Jehn, K. A., Zanutto, E. L., & Thatcher, S. M. B. (2009).  Do workgroup faultlines help or hurt? A moderated model of faultlines, team identification, and group performance.  Organization Science,  20,  35-50. 
  doi:10.1287/orsc.1080.0379

Gibson, C., & Vermeulen, F. (2003). A healthy divide: Subgroups as a stimulus for team learning behavior. Administrative Science Quarterly,  48,  202-239. 
  doi:10.2307/3556657

Lawrence, B., & Zyphur, M. (2011). Identifying organizational faultlines with latent class cluster analysis. Organizational Research Methods,  14,  32-57. 
  doi:10.1177/1094428110376838

Meyer, B., & Glenz, A.  (in press).  Diversity faultline measures:  A computational comparison and a new approach to multiple subgroups.  Organizational Research Methods.

Shaw, J. (2004). The development and analysis of a measure of group faultlines. Organizational Research Methods, 7, 66-100.  
  doi:10.1177/1094428103259562

Thatcher, S., Jehn, K. A., & Zanutto, E. (2003). Cracks in diversity research: the effects of diversity faultlines on conflict and performance. Group Decision and Negotiation,  12,  217-241. 
  doi:10.1023/A:1023325406946

Trezzini, B. (2008). Probing the group faultline concept: An evaluation of measures of patterned multi-dimensional group diversity. Quality and Quantity,  42,  339-368. 
  doi:10.1007/s11135-006-9049-z

van Knippenberg, D., Dawson, J., West, M., & Homan, A. (2011). Diversity faultlines, shared objectives, and top management team performance. Human Relations,  64,  307-336. 
  doi:10.1177/0018726710378384
}
\author{
Bertolt Meyer,  b.meyer@psychologie.uzh.ch
Andreas Glenz,  a.glenz@psychologie.uzh.ch

}
\note{
%%  ~~further notes~~
}

\section{Warning }{This is a beta version of the package. Please report any bugs to the authors.}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
library(parallel)
# create a data set consisting of two teams with six members each
# with information on age, gender, and ethnicity for each team member
# and a variable teamid that identifies team membership in the first column
mydata <- data.frame(teamid = c(rep(1,6),rep(2,6)),  
                     age = c(44,18,40,33,33,50,22,23,39,42,57,51),  
                     gender = c("f","m","f","f","m","f","f","f","m","m","m","m"),  
                     ethnicity = c("A","B","A","D","C","B","A","A","B","B","C","C"))
 
# specify the scale for each attribute                     
my_types <- c("numeric", "nominal", "nominal")

# specify the weight for each attribute
my_weights <- c(0.1, 1, 1)

# calculate ASW and save the result in a data frame my_asw
my_ASW <- faultlines(data = mydata,  
                     group.par = "teamid",  
                     attr.type = my_types,   
                     attr.weight = my_weights,   
                     method = "asw",
                     cores=1)
my_ASW

# calculate ASW on the individual level and save the result in a data frame my_asw_i
my_ASW_i <- faultlines(data = mydata,  
                       group.par = "teamid",  
                       attr.type = my_types,  
                       attr.weight = my_weights,  
                       method = "asw",  
                       i.level = TRUE,
                       cores=1)
my_ASW_i

# NOT RUN:
# my_ASW
#   team teamsize  fl.value mbr_to_subgroups number_of_subgroups
# 1    1        6 0.3317999 1, 2, 1, 2, 2, 1                   2
# 2    2        6 0.8054895 1, 1, 2, 2, 3, 3                   3

# calculate Thatcher et al.'s Fau and save the result in a data frame my_faus
my_faus <- faultlines(data = mydata,  
                      group.par = "teamid",  
                      attr.type = my_types,  
                      attr.weight = my_weights,  
                      method = "thatcher",
                      cores=1)

# NOT RUN:
# my_faus
#   team teamsize  fl.value mbr_to_subgroups number_of_subgroups
# 1    1        6 0.5513439 1, 2, 1, 2, 2, 1                   2
# 2    2        6 0.7747787 1, 1, 2, 2, 2, 2                   2


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. Call RShowDoc("KEYWORDS") from within R
% to see stabdard keywords
\keyword{cluster}
\keyword{multivariate}% __ONLY ONE__ keyword per line
