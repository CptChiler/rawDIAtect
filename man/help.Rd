\name{rawDIAtect}
\alias{rawDIAtect}
\title{rawDIAtect}
\usage{
rawDIAtect(path_in = "path/to/DIA-NN/outputs")
}
\description{
Data analysis tool for MS DIA-NN outputs using the rawDIAtect workflow.
}
\examples{
Minimal input:
rawDIAtect(path_in = "/Quants")

iso_diff = 2 (0 = Off)
pep_filter = 3 (<1 = Off)
path_out = (Default = path_in)

Expert input:
rawDIAtect(path_in = "/Quants", pep_filter = 3 ,iso_diff = 2, path_out = "/AMRs", Exp_name = "Induced")
}
