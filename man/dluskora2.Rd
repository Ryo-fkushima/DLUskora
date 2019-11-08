\name{dluskora2}
\alias{dluskora2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{dluskora2:A simulator of diffusion-limited REE uptake by garnets (a little different from dluskora1)
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
"dluskora2" resembles dluskora1 but this can be  more useful if you'd like to simulate REE profile with complicated growth of garnets. You can set radius-time path as four points and three straight lines.
}
\usage{
dluskora2(fac1, Q, syR, c_ave, D_0, R, T_1, T_2, K_d, Mr, Mt, ft, garsize, Tm_1, Tm_2, gmratio)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fac1}{
This is a "diffusion length" and equal to D_0/A [cm]. A is the garnet growth coefficient (r = A t^(fg)). The meanings of D_0 and fg are shown below.
}
\item{Q}{
Activation energy of REE diffusion around the garnet [J/mol]
}
\item{syR}{
The size of the system concluding only a big garnet [cm]
}
\item{c_ave}{
The initial concentration of REE on the system [ppm]
}
\item{D_0}{
A pre-exponential factor of REE diffusion around the garnet [cm^2/year]
}
\item{R}{
The universal gas constant [m^2 kg s^(-2) K^(-1) mol^(-1)]
}
\item{T_1, T_2}{
This function can simulate a situation where the temperature around the system increased and you can set the initial and final temperature. T_1 is the initial temperature [degree Celsius] and T_2 is the final temperature [degree Celsius].
}
\item{K_d}{
partition coefficient > 1 (REE concentration of garnet over REE concentration of the matrix)
}
\item{Mr, Mt}{
They are the numbers of spatial(radial) and time meshes respectively.
}

\item{ft}{
This is an exponent of the temperature-increase law. If you set ft as 1, the temperature will increase with the linier increase rate law: T = B t.
}
\item{garsize}{
The final size of the garnet [cm]
}
\item{Tm_1, Tm_2}{
}
\item{gmratio}{
The final size of the garnet [cm]
}


}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{
%% ~put references to the literature/web site here ~
Crank (1975) The mathematics of diffusion. Oxford University Press, p414.

Skora et al. (2006) Diffusion-limited REE uptake by eclogite garnets and its consequences for Lu-Hf and Sm-Nd geochronology. Contrib. Mineral Petrol., 152:703-720.
}
\author{
%%  ~~who you are~~
Ryo Fukushima
}
\note{
%%  ~~further notes~~
>K_d and A do not depend on temperature in this function. Temperature increase affects only the REE diffusion coefficient around the garnet.

>Diffusion in the garnet is ignored.

>This program is based on a Crank-Nicholson scheme. If you set too large Mt, the result would not be appropriate.

>This function returns "mass_gain_percent" with the REE profile plot. Please check whether its absolute value is very small(~e-4). Your result is not appropriate when the value is too big.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
dluskora2(fac1=(1.78 *10^(21)), Q=300000, syR=0.60, c_ave=30,
D_0=(4 * 10^(13)), R=8.3, T_1=450, T_2=600, K_d=15, Mr=100,
Mt=45, ft=1, garsize=0.27, Tm_1=500, Tm_2=550, gmratio=0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{REE}% use one of  RShowDoc("KEYWORDS")
\keyword{diffusion-limited}% __ONLY ONE__ keyword per line
\keyword{garnet}
\keyword{eclogite}
