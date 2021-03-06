% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pso.R
\name{pso}
\alias{pso}
\title{Non-Linear Pseudo Second-Order Adsorption Kinetic Model}
\usage{
pso(t, qt, qe)
}
\arguments{
\item{t}{the numerical value for contact time}

\item{qt}{the numerical value for the amount adsorbed at time t}

\item{qe}{the numerical value for the amount adsorbed at equilibrium. If this parameter is not defined, it will be estimated.}
}
\value{
the non linear regression and the parameters for the pseudo-second order non-linear model analysis
}
\description{
The Pseudo-Second-Order Adsorption Kinetic Model is an empirical rate equation known to be the simplified second-order expression of the Pseudo-First Order Adsorption Kinetic Model. It is widely applied to adsorption systems, from biomass to nanomaterials as adsorbent and from heavy metals to pharmaceuticals as adsorbate or contaminant (Revellame, Fortela, Sharp, Hernandez, and Zappi, 2020).
}
\examples{
\donttest{
t <- c(0,15,30,45,60,75,90,105,120)
qt <- c(0.000,3.718,3.888,4.102,4.274,4.402,4.444,4.488,4.616)
qe <- 4.68
pso(t,qt,qe)}
}
\references{
Ho, Y. S., &; Mckay, G. (1999) <doi:10.1016/S0032-9592(98)00112-5> Pseudo-second order model for sorption processes. In Process Biochemistry (Vol. 34).

Revellame, E. D., Fortela, D. L., Sharp, W., Hernandez, R., &; Zappi, M. E. (2020) <doi:10.1016/j.clet.2020.100032> Adsorption kinetic modeling using pseudo-first order and pseudo-second order rate laws: A review. In Cleaner Engineering and Technology (Vol. 1). Elsevier Ltd.
}
\author{
Jeff Ryan S. Magalong

Joshua Z. DelaCruz

Jeann M. Bumatay

Chester C. Deocaris
}
