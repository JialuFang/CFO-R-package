% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CFO.next.R
\name{CFO.next}
\alias{CFO.next}
\title{Determination of the dose level for next cohort in the calibration-free odds (CFO) design for phase I trials}
\usage{
CFO.next(target, cys, cns, currdose, 
       prior.para = list(alp.prior = target, bet.prior = 1 - target),
       cutoff.eli = 0.95, early.stop = 0.95)
}
\arguments{
\item{target}{the target DLT rate.}

\item{cys}{the cumulative numbers of DLTs observed at the left, current, and right dose levels.}

\item{cns}{the cumulative numbers of patients treated at the left, current, and right dose levels.}

\item{currdose}{the current dose level.}

\item{prior.para}{the prior parameters for a beta distribution, where set as \code{list(alp.prior = target, bet.prior = 1 - target)}
by default, \code{alp.prior} and \code{bet.prior} represent the parameters of the prior distribution for
the true DLT rate at any dose level. This prior distribution is specified as Beta(\code{alpha.prior}, \code{beta.prior}).}

\item{cutoff.eli}{the cutoff to eliminate overly toxic doses for safety. We recommend
the default value of \code{cutoff.eli = 0.95} for general use.}

\item{early.stop}{the threshold value for early stopping. The default value \code{early.stop = 0.95}
generally works well.}
}
\value{
The \code{CFO.next()} function returns a list object comprising the following elements:
\itemize{
\item target: the target DLT rate.
\item cys: the cumulative counts of DLTs observed at the left, current, and right dose levels.
\item cns: the cumulative counts of patients treated at the left, current, and right dose levels.
\item decision: the decision in the CFO design, where \code{left}, \code{stay}, and \code{right} represent the
movement directions, and \code{stop} indicates stopping the experiment.
\item currdose: the current dose level.
\item nextdose: the recommended dose level for the next cohort. \code{nextdose = 99} indicates that the trial is
terminated due to early stopping.
\item overtox: the situation regarding which positions experience over-toxicity. The dose level indicated
by \code{overtox} and all the dose levels above experience over-toxicity. \code{overtox = NA} signifies that
the occurrence of over-toxicity did not happen.
\item toxprob: the expected toxicity probability, \eqn{Pr(p_k > \phi | x_k, m_k)}, at the left, current, and
right dose levels, where \eqn{p_k}, \eqn{x_k}, and \eqn{m_k} is the dose-limiting toxicity (DLT) rate, the
numbers of observed DLTs, and the numbers of patients at dose level \eqn{k}.
}
}
\description{
In the CFO design for phase I trials, the function is used to determine the dose movement based on the toxicity outcomes of the enrolled cohorts.
}
\details{
The CFO design determines the dose level for the next cohort by assessing evidence from the current
dose level and its adjacent levels. This evaluation is based on odds ratios denoted as \eqn{O_k}, where
\eqn{k = L, C, R} represents left, current (central), and right dose levels. Additionally, we define \eqn{\overline{O}_k = 1/O_k}.
The ratio \eqn{O_C / \overline{O}_{L}} indicates the inclination for de-escalation, while \eqn{\overline{O}_C / O_R}
quantifies the tendency for escalation. Threshold values \eqn{\gamma_L} and \eqn{\gamma_R} are chosen to
minimize the probability of making incorrect decisions. The decision process is summarized in Table 1
of Jin and Yin (2022).
The early stopping and dose elimination rules are implemented to ensure patient safety. If the data suggest excessive
toxicity at the current dose level, we exclude that dose level and those higher levels. If the lowest dose level is overly toxic,
the trial will be terminated according to the early stopping rule.
}
\note{
When the current dose level is the lowest or highest (i.e., at the boundary), the parts in \code{cys},
\code{cns}, and \code{toxprob} where there is no data are filled with \code{NA}. \cr
The dose level indicated by \code{overtox} and all the dose levels above experience over-toxicity, and these dose levels will be eliminated.
}
\examples{
## determine the dose level for the next cohort of new patients
cys <- c(0, 1, 0); cns <- c(3, 6, 0)
decision <- CFO.next(target=0.2, cys=cys, cns=cns, currdose=3)
summary(decision)

cys <- c(NA, 3, 0); cns <- c(NA, 3, 0)
decision <- CFO.next(target=0.2, cys=cys, cns=cns, currdose=1)
summary(decision)

cys <- c(0, 3, NA); cns <- c(3, 3, NA)
decision <- CFO.next(target=0.2, cys=cys, cns=cns, currdose=7)
summary(decision)

}
\references{
Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials.
\emph{Statistical Methods in Medical Research}, 31(6), 1051-1066.
}
\author{
Jialu Fang, Wenliang Wang, Ninghao Zhang, and Guosheng Yin
}
