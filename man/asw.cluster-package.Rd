\name{asw.cluster-package}
\alias{asw.cluster-package}
\alias{asw.cluster}
\docType{package}
\title{Diversity Faultlines}
\description{
Package for calculating diversity faultline measures as suggested by Thatcher, Jehn, and Zanutto (2003), by Gibson and Vermeulen (2003), by Shaw (2004), by Trezzini (2008), by Bezrukova, Jehn, and Zanutto (2009), by van Knippenberg, Dawson, West, and Homan (2011) and by Lawrence and Zyphur (2011).  The package can also calculate a group's Faultline with the maximal average silhouette width (ASW) by using agglomerative and partitioning clustering algorithms suggested by Meyer and Glenz (2013).
}
\details{
\tabular{ll}{
Package: \tab asw.cluster\cr
Type: \tab Package\cr
Version: \tab 2.10.2\cr
Date: \tab 2023-02-23\cr
License: \tab For Non-commercial Use Only\cr
}
  The \code{\link{faultlines}} function is the main function of the package and is intended for calculating the faultlines; most other functions in the package are workhorse and sub functions that are called by \code{faultlines()}.
}
\author{
Andreas Glenz, University of Zurich, Switzerland
Bertolt Meyer, TU-Chemnitz, Germany

Maintainer: Andreas Glenz <a.glenz@uzh.ch>

}
\references{

Bezrukova, K., Jehn, K. A., Zanutto, E. L., & Thatcher, S. M. B. (2009). Do workgroup faultlines help or hurt? A moderated model of faultlines, team identification, and group performance. Organization Science, 20, 35-50. doi:10.1287/orsc.1080.0379

Gibson, C., & Vermeulen, F. (2003). A healthy divide: Subgroups as a stimulus for team learning behavior. Administrative Science Quarterly, 48, 202-239. doi:10.2307/3556657

Lawrence, B., & Zyphur, M. (2011). Identifying organizational faultlines with latent class cluster analysis. Organizational Research Methods, 14, 32-57. doi:10.1177/1094428110376838

Meyer, B., & Glenz, A.  (2013). Team faultline measures: A computational comparison and a new approach to multiple subgroups. Organizational Research Methods. Advance online publication. doi:10.1177/1094428113484970

Shaw, J. (2004). The development and analysis of a measure of group faultlines. Organizational Research Methods, 7, 66-100. doi:10.1177/1094428103259562

Thatcher, S., Jehn, K. A., & Zanutto, E. (2003). Cracks in diversity research: the effects of diversity faultlines on conflict and performance. Group Decision and Negotiation, 12, 217-241. doi:10.1023/A:1023325406946

Trezzini, B. (2008). Probing the group faultline concept: An evaluation of measures of patterned multi-dimensional group diversity. Quality and Quantity, 42, 339-368. doi:10.1007/s11135-006-9049-z

van Knippenberg, D., Dawson, J., West, M., & Homan, A. (2011). Diversity faultlines, shared objectives, and top management team performance. Human Relations, 64, 307-336. doi:10.1177/0018726710378384
}

\keyword{ package }
\seealso{
\code{\link{faultlines}}
}

