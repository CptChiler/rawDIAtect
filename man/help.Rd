\name{rawDIAtect}
\alias{rawDIAtect}
\title{rawDIAtect(path_in = "path", pep_filter = 3)}
\usage{
rawDIAtect(path_in = "path")
}
\description{
Data analysis tool for MS DIA-NN outputs using the rawDIAtect workflow.
}
\examples{
Minimal input:
rawDIAtect(path_in = "/Quants")

Expert input:
rawDIAtect(path_in = "/Quants", pep_filter = 3 , path_out = "/AMRs" , whitelist = 0 , Exp_name = "Induced")
}
