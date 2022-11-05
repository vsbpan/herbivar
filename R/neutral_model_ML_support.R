#'@title Neutral 'Bite Size' Distribution Based On Allometric Scaling Laws
#'@description Density, distribution function, quantile function, and random generation for the neutral 'bite size' distribution based on allometric scaling laws.
#' @param x,q a vector of quantities
#' @param p a vector of probabilities
#' @param n number of observations to generate.
#' @param min.phi the minimum bite size in terms of proportion leaf herbivory. Defaults to 0.005 (0.5%).
#' @param max.phi the maximum bite size in terms of proportion leaf herbivory Defaults to 1 (100%).
#' @param a the combined allometric scaling coefficient. Defaults to 14/9.
#' @param log,log.p logical; if \code{TRUE}, probabilities p are given as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P(X \leq x)} otherwise, \eqn{P(X > x)}
#' @details
#' The neutral 'bite size' distribution has density function
#' \deqn{P(\phi) = \frac{1-\alpha}{\phi_M^{1-\alpha} - \phi_m^{1-\alpha}} \phi^{-\alpha}} when \eqn{\alpha \neq 1} and
#' \deqn{P(\phi) = \frac{1}{\ln{\phi_M} - \ln{\phi_m}} \phi^{-1}} when \eqn{\alpha = 1}
#' where \eqn{\phi_m} is the minimum 'bite size' and \eqn{\phi_M} is the maximum 'bite size' in terms of proportion leaf herbivory, and \eqn{\alpha} is the combined allometric scaling coefficient defined as \deqn{\alpha = - \frac{\alpha_N + \alpha_S + 1 - \alpha_I}{\alpha_I}.} \eqn{\alpha_N} is the allometric scaling exponent between population density \eqn{N} and body mass \eqn{W} such that \eqn{N \propto W^{\alpha_N}}. A priori value is \eqn{-\frac{3}{4}} according to Damuth’s rule (Damuth 1981). \eqn{\alpha_S} is the allometric scaling exponent of species richness \eqn{S_j} among a body mass class \eqn{W_j} such that \eqn{S_j \propto W_{j}^{\alpha_S}}. A priori value is \eqn{-\frac{2}{3}} according to Hutchinson & MacArthur (1959) and May (1978).\eqn{\alpha_I} is the allometric scaling exponent of whole body metabolic rate \eqn{I_j} among a body mass class \eqn{W_j} such that \eqn{I_j \propto W_{j}^{\alpha_I}}. A priori value is \eqn{\frac{3}{4}} according to Kleiber's law (1932).
#'
#' The cumulative density function is
#' \deqn{P(\phi \leq q) = \frac{q^{1-\alpha} - \phi_m^{1-\alpha}}{\phi_M^{1-\alpha} - \phi_m^{1-\alpha}}} when \eqn{\alpha \neq 1} and
#' \deqn{P(\phi \leq q) = \frac{\ln{q} - \ln{\phi_m}}{\ln{\phi_M} - \ln{\phi_m}}} when \eqn{\alpha = 1}.
#'
#' @return a vector of numeric values
#' @references
#' Damuth, J. 1981. Population density and body size in mammals. Nature 290:699–700.
#'
#' Hutchinson, G. E., and R. H. MacArthur. 1959. A Theoretical Ecological Model of Size Distributions Among Species of Animals. The American Naturalist 93:117–125.
#'
#' Kleiber, M. 1932. Body size and metabolism. Hilgardia 6:315–353.
#'
#' May, R. M. 1978. The dynamics and diversity of insect faunas. Page in L. A. Mound and N. Waloff, editors. Diversity of insect faunas. Published for the Royal Entomological Society by Blackwell Scientific Publications; Distributed by Halsted Press, Oxford England; New York.
#'
#' @rdname allo
#' @export
dallo<-function(x,min.phi=0.005,max.phi=1,a=14/9,log = FALSE){
  # Check if phi is within bound. If so, apply pdf formula, otherwise 0
  if(a == 1){
    prob<-ifelse((x<=max.phi)*(x>=min.phi),1/(log(max.phi)-log(min.phi))/x,0)
  } else {
    prob<-ifelse((x<=max.phi)*(x>=min.phi),(1-a)/(max.phi^(1-a)-min.phi^(1-a))/x^a,0)
  }
  if(log){
    prob<-log(prob)
  }
  return(prob)
}

#' @rdname allo
#' @export
rallo<-function(n,min.phi=0.005,max.phi=1,a=14/9){
  cdf.prob<-runif(n,min = 0,max = 1)
  # inverse transform sampling
  phi <- qallo(cdf.prob, min.phi = min.phi, max.phi = max.phi, a = a)
  return(phi)
}

#' @rdname allo
#' @export
pallo <- function(q, min.phi = 0.005, max.phi = 1, a = 14/9,
                 lower.tail = TRUE, log.p = FALSE){
  if(a == 1){
    p <- ifelse(q > 1,
                1,
                ifelse(q < min.phi,
                       0,
                       (log(q) - log(min.phi)) / (log(max.phi) - log(min.phi))
                )
    )
  } else {
    p <- ifelse(q > 1,
                1,
                ifelse(q < min.phi,
                       0,
                       (q^(1-a) - min.phi^(1-a)) / (max.phi^(1-a) - min.phi^(1-a))
                )
    )
  }

  if(!lower.tail){
   p <- 1-p
  }

  if(log.p){
    p <- log(p)
  }
  return(p)
}

#' @rdname allo
#' @export
qallo <- function(p, min.phi = 0.005, max.phi = 1, a = 14/9,
                lower.tail = TRUE, log.p = FALSE){
  if(log.p){
    p <- exp(p)
  }
  if(!lower.tail){
    p <- 1 - p
  }
  if(a == 1){
    q <- ifelse(p > 1 | p < 0,
                NaN,
                exp(p*(log(max.phi)-log(min.phi))+log(min.phi))
                )
  } else {
    q <- ifelse(p > 1 | p < 0,
                NaN,
                (p*(max.phi^(1-a)-min.phi^(1-a))+min.phi^(1-a))^(1/(1-a))
                )
  }

  return(q)
}

#' @title Generate Neutral Herbivory Distribution Based On Allometric Scaling Laws Using Simulation
#' @description This is the workhorse of \code{ralloT()} that generate draws from the neutral herbivory distribution using simulation. The CDF of the neutral herbivory model is too computationally intensive to calculate, so inverse transform sampling is not implemented.
#' @details
#' \eqn{\lambda} of a Poisson distribution is calculated from parameters \eqn{\phi_m}, \eqn{\phi_M}, \eqn{a}, and \eqn{\overline{\phi_{T}'}} and used to draw random number of feeding events on a leaf \eqn{k_i}. \deqn{\lambda = \frac{\overline{\phi_{T}'}}{\overline{\phi}}}
#' \eqn{\phi_{Ti}} is then obtained via \deqn{\phi_{Ti} = \sum^{k_i}_{j=1} \phi_j} if \eqn{\phi_{Ti} \leq 1}, otherwise \deqn{\phi_{Ti} = 1,} where
#' \deqn{P(\phi) = \frac{1-\alpha}{\phi_M^{1-\alpha} - \phi_m^{1-\alpha}} \phi^{-\alpha}} when \eqn{\alpha \neq 1} and
#' \deqn{P(\phi) = \frac{1}{\ln{\phi_M} - \ln{\phi_m}} \phi^{-1}} when \eqn{\alpha = 1}
#' \eqn{\overline{\phi_{T}'}} is the mean herbivory of the distribution when herbivores are not plant limited. That is, if \eqn{\overline{\phi_{T}}}, which has support \eqn{[x, \infty]}, is not truncated to \eqn{[0, 1]}. \eqn{\overline{\phi_{T}}} and \eqn{\overline{\phi_{T}'}} are often close to each other. For more details on the 'bite size' distribution see \code{?rallo()}.
#' @param mean.phi.T The mean herbivory of the distribution when herbivores are not plant limited. See details.
#' @param min.phi the minimum bite size in terms of proportion leaf herbivory. Defaults to 0.005 (0.5%).
#' @param max.phi the maximum bite size in terms of proportion leaf herbivory Defaults to 1 (100%).
#' @param a the combined allometric scaling coefficient. Defaults to 14/9.
#' @param n.sim the number of random numbers to draw.
#' @param truncate if \code{TRUE} (default), truncate generated values of \eqn{\phi_T} to \eqn{[0, 1]}.
#' @return a vector of numeric values
#' @export
.allometry.herb.quasi.sim<-function(mean.phi.T,min.phi=0.005,max.phi=1,a=14/9,
                                   n.sim=1000,truncate=T){
  #Forward approximate simulation (not exact!)
  lambda<-get_lambda(mean.phi.T = mean.phi.T,min.phi = min.phi, max.phi = max.phi, a = a)
  k<-rpois(n.sim,lambda)
  phi_T<-vapply(X = k,
                FUN = function(x) sum(
                  rallo(n = x,
                        min.phi = min.phi,
                        max.phi = max.phi,
                        a = a)),
                FUN.VALUE = numeric(1))
  if(truncate){
    phi_T[phi_T>1]<-1 # Cut off phi_T if goes over 1. May underestimate mean.phi.T
  }
  return(phi_T)
}


#' @title Neutral Herbivory Distribution Based On Allometric Scaling Laws
#' @description Density, distribution function, quantile function, and random generation for the neutral herbivory distribution based on allometric scaling laws.
#' @details
#' The neutral herbivory distribution is a type of compound Poisson distribution taking the form:
#' \deqn{\phi_{Ti} = \sum^{k_i}_{j=1} \phi_j,} if \eqn{\phi_{Ti} \leq 1}, otherwise \deqn{\phi_{Ti} = 1,} where \deqn{P(\phi) = \frac{1-\alpha}{\phi_M^{1-\alpha} - \phi_m^{1-\alpha}} \phi^{-\alpha},} when \eqn{\alpha \neq 1} and
#' \deqn{P(\phi) = \frac{1}{\ln{\phi_M} - \ln{\phi_m}} \phi^{-1},} when \eqn{\alpha = 1} and \deqn{k \sim Pois(\lambda = \frac{\overline{\phi_{T}'}}{\overline{\phi}})}
#' \eqn{\overline{\phi_{T}'}} is the mean herbivory of the distribution when herbivores are not plant limited. That is, if \eqn{\overline{\phi_{T}}}, which has support \eqn{[x, \infty]}, is not truncated to \eqn{[0, 1]}. \eqn{\overline{\phi_{T}}} and \eqn{\overline{\phi_{T}'}} are often close to each other. For more details on the 'bite size' distribution see \code{?rallo}.
#'
#' The PDF of \eqn{P(\phi_T)} does not have a closed form solution and is therefore numerically approximated. We first marginalize out \eqn{k} via
#' \deqn{P(\phi_T)=\sum_{k=0}^\infty P(\phi_T|k) P(k|\lambda ).}
#' As \eqn{k} increases, \eqn{P(k|\lambda)} becomes vanishingly small, so the infinite sum can be cut off at a rather low \code{k.max} value without loosing much accuracy. \code{k.max.tolerance} sets the tolerance threshold for the cut off in the calculation. More precisely, the function checks to see whether \eqn{P(k.max | \lambda)} is less than the threshold of triviality. If exceeded, the function would throw a warning.
#'
#' The calculation of the conditional probability distribution \eqn{P(\phi_T | k)} is given by
#' \deqn{P(\phi_T | k) = P(\phi_T^k)}
#' when \eqn{\phi_T < 1}. The truncated probabilities \eqn{P(\phi_T > 1 | k)} is added back to \eqn{P(\phi_T = 1 | k)} such that \eqn{P(\phi_T \leq 1 | k)} sums up to 1. Here, \eqn{P(\phi_T^k)} is the kth convolution of \eqn{P(\phi)}. \eqn{P(\phi_T^k)} is computationally expensive. Thankfully, we can use the Fast Fourier Transformation (FFT) to simplify the calculation significantly. Convolution of functions is simply the product of those functions in the frequency domain that is then back transformed into the time domain.
#' \deqn{P(\phi_T^k) = \overbrace{P(\phi)*P(\phi)*...*P(\phi)}^{\text{k times}} = \mathcal{F}^{-1} [\mathcal{F}[P(\phi)]^k]}
#' Argument \code{by} sets the grid resolution of the discrete Fourier transform. Usually, a value below 0.001 is required to achieve reasonable accuracy. Eights times as many zeros are added to the end of the PDF to improve the accuracy and efficiency of the FFT (i.e. zero padding). To avoid overflow and improve computational efficiency, convolutions with \eqn{k > 100} are approximated with a normal distribution, given the central limit theorem (see \code{?get_phi_bar()} and \code{?get_phi_var()}):
#' \deqn{P(\phi_T^k) \approx \mathcal{N}(\mu = k\overline\phi,\sigma=\sqrt{k\mathbb{Var}[\phi]})}
#'
#' The CDF of the neutral herbivory distribution is calculated numerically by adding up the density. Because the neutral herbivory distribution is a mix of discrete and continuous, the CDF is the sum of the discrete portion \deqn{P(\phi_T = 0, \overline{\phi_T'}, \phi_m, \phi_M, \alpha) = e^{-\lambda}}
#'  and the integral of the continuous portion \deqn{\int_0^q P(\phi_T = q, \overline{\phi_T'}, \phi_m, \phi_M, \alpha) d \phi_T,} and when \eqn{q = 1}, \deqn{P(\phi_T \leq q) = 1.}
#'
#' Because the CDF of the neutral herbivory distribution is too computationally intensive, \code{ralloT()} generates values using random draws from the Poisson and 'bite size' distribution. It is essentially a wrapper for the function \code{.allometry.herb.quasi.sim()}. For the same reason, \code{qalloT()} estimates the quantile function by generating \code{n.sim} draws from the neutral herbivory distribution then finding the empirical cumulative density function.
#'
#'
#' @param x,q a vector of values of proportion herbivory
#' @param p a vector of probabilities
#' @param n the number of observations to generate
#' @param mean.phi.T The mean herbivory of the distribution when herbivores are not plant limited. If  \code{NULL}, supplied \code{lambda} value is converted to \code{mean.phi.T}. See details.
#' @param min.phi the minimum bite size in terms of proportion leaf herbivory. Defaults to 0.005 (0.5%).
#' @param max.phi the maximum bite size in terms of proportion leaf herbivory Defaults to 1 (100%).
#' @param a the combined allometric scaling coefficient. Defaults to 14/9.
#' @param lambda An alternative parameterization of mean.phi.T. see \code{?get_lambda()}
#' @param k.max The maximum number of convolutions of the neutral 'bite size' distribution in numerical approximation of the PDF of the neutral herbivory distribution. Default is 50.
#' @param by The grid resolution used in the FFT convolutions. Default is 0.001.
#' @param k.max.tolerance the tolerance threshold of maximum convolution cut off (ideally probabilities above \code{k.max} convolutions is vanishingly small).
#' @param k.fft.limit the maximum number of convolutions performed by FFT. For k convolutions above this limit, a Gaussian approximation is used. See details for more information.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are \eqn{P(X \leq x)} otherwise, \eqn{P(X > x)}
#' @param log,log.p if \code{TRUE} (default is \code{FASLE}), return probabilities on the log scale.
#' @param parallel if \code{TRUE}, use parallel computing supported via the package \code{foreach} to speed up the computation. The number of parallel processes depends on the argument \code{cores}. Ignored if \code{cores} is set  above 1.
#' @param cores the number of parallel processes in parallel computing.
#' @param n.sim the number of simulations used to numerically estimate the quantile function. Defaults to 1000.
#' @param ... addition arguments passed to \code{dalloT()} in \code{palloT()}.
#' @return A vector of numeric values
#' @rdname alloT
#' @export
ralloT<-function(n, mean.phi.T = NULL, min.phi = 0.005, max.phi = 1, a = 14/9, lambda = NULL){
  if(is.null(mean.phi.T)){
    if(is.null(lambda)){
      stop("Either 'lambda' or 'mean.phi.T must be specified'.")
    }
    mean.phi.T <- get_mean_phi_T(lambda,
                                 min.phi = min.phi,
                                 max.phi = max.phi,
                                 a = a)
  }

  phi.T<-.allometry.herb.quasi.sim(mean.phi.T = mean.phi.T,
                                  n.sim = n,
                                  min.phi = min.phi,
                                  max.phi = max.phi,
                                  a = a,
                                  truncate = T)
  # while(any(phi.T>1)){
  #   phi.T[phi.T>1] <- .allometry.herb.quasi.sim(mean.phi.T = mean.phi.T,
  #                                              n.sim = sum(phi.T>1),
  #                                              min.phi = min.phi,
  #                                              max.phi = max.phi,
  #                                              a = a,
  #                                              truncate = F)
  # }
  # This method has less numeric error than sampling using dalloT() directly.
  # If using sample(), the sampling grid resolution needs to be at least <= 0.0001
  # which is computationally intensive. If the resolution is too low, there would be
  # a downward bias on the variance and mean approximation
  return(phi.T)
}

#' @rdname alloT
#' @export
dalloT<-function(x, mean.phi.T = NULL, min.phi = 0.005, max.phi = 1, a = 14/9, lambda = NULL,
                 k.max = 50, by = 0.001, k.max.tolerance = 1e-5,
                 k.fft.limit = 100, log = FALSE, parallel = FALSE, cores = 1){
  if(is.null(mean.phi.T)){
    if(is.null(lambda)){
      stop("Either 'lambda' or 'mean.phi.T must be specified'.")
    }
    mean.phi.T <- get_mean_phi_T(lambda,
                                 min.phi = min.phi,
                                 max.phi = max.phi,
                                 a = a)
  }

  if(cores>1){
    parallel <- TRUE
    cluster <- makeCluster(cores)
    doParallel::registerDoParallel(cluster)

  } else {
    parallel <- FALSE
  }

  if(max.phi<=min.phi){
    prob<-rep(0,length(x))
  } else {
    # if(by > 0.01){
    #   warning("Approximation resolution too low; results are crappy")
    # }
    # mean.phi.T is the mean before cutting off at 1 (!!).
    # mean.phi.T is an overestimation of the true mean, especially at higher values
    k.vec <- seq(0,k.max,by = 1) # Set the vector of k that is calculated.
    # 50 usually does a good approximation

    phi <- seq(from=0,to = 1,by = by) # Generate some x
    p.phi <- dallo(phi,
                   min.phi = min.phi,
                   max.phi = max.phi,
                   a = a,
                   log = FALSE) # Calculate the pdf
    #Add zero padding to get the length to a power of 2 (for computational efficiency)
    p.phi.length<-2^ceiling(log2(length(p.phi))+2)
    p.phi<-c(p.phi,rep(0,p.phi.length-length(p.phi)))

    #Rescale p.phi to avoid over or under flow
    p.phi<-p.phi/sum(p.phi)
    Ft.p.phi <- fft(p.phi,inverse = FALSE) # Apply FFT

    #Pre-calculate indices
    phi.T.index<-findInterval(x = x, vec = phi)

    #Find lambda outside of .dk.cond.lambda() to save computation
    lambda <- get_lambda(mean.phi.T,
                         min.phi = min.phi,
                         max.phi = max.phi,
                         a = a)

    #Numeric summation of joint prob from k=0 to k=max_k for each observed x (phi.T)
    prob<-as.vector(tcrossprod(.dk.cond.lambda(k = k.vec,
                                              k.max.tolerance = k.max.tolerance,
                                              lambda = lambda,
                                              log = FALSE),
                               .dalloT.cond.k.fft.conv(phi.T = x,
                                                      k = k.vec,
                                                      fft.vec = Ft.p.phi,
                                                      phi = phi,
                                                      log = FALSE,
                                                      min.phi = min.phi,
                                                      max.phi = max.phi,
                                                      a = a,
                                                      k.fft.limit = k.fft.limit,
                                                      phi.T.index = phi.T.index,
                                                      parallel = parallel)
    ))
  } # Parallel computing only possible with registered cluster

  #Fix scaling issue
  prob <- ifelse(x == 0 | x == 1,
                 prob,
                 prob/by)

  if(log){
    prob<-log(prob)
  }
  return(prob)

  on.exit(
    try({
      if(!is.null(cluster)){
        doParallel::stopImplicitCluster()
        stopCluster(cluster)
      }
    })
  ) # Stop all connections on exit
}

#' @rdname alloT
#' @export
palloT<-function(q, mean.phi.T = NULL, min.phi = 0.005, max.phi = 1, a = 14/9, lambda = NULL,
                 by = 0.001, k.max = 50, k.max.tolerance = 1e-5,
                 k.fft.limit = 100, parallel = FALSE, cores = 1,
                 lower.tail = TRUE, log.p = FALSE, ...){
  if(is.null(mean.phi.T)){
    if(is.null(lambda)){
      stop("Either 'lambda' or 'mean.phi.T must be specified'.")
    }
    mean.phi.T <- get_mean_phi_T(lambda,
                                 min.phi = min.phi,
                                 max.phi = max.phi,
                                 a = a)
  }
  lambda <- get_lambda(mean.phi.T,
                           min.phi = min.phi,
                           max.phi = max.phi,
                           a = a)

  if(any((q > 0) & ((by > min.phi & q < min.phi) | (by > q)))){
    warning("Numeric resolution not high enough, result is unreliable. Lower the 'by' setting.")
  }

  x <- seq(0,1, by = by) # Set up a grid to numerically integrate at 'by' resolution

  x.prob <- dalloT(x = x, # Generate probability
         lambda = lambda,
         min.phi = min.phi,
         max.phi = max.phi,
         a = a,
         k.max = k.max,
         by = by,
         k.max.tolerance = k.max.tolerance,
         k.fft.limit = k.fft.limit,
         log = FALSE,
         parallel = parallel,
         cores = cores)
  x.prob[-c(1,length(x.prob))] <- x.prob[-c(1,length(x.prob))] * by # Fix scaling issue for the continuous interval
  x.prob <- c(0, x.prob) # Append prob = 0 to any q less than 0
  x.cdf <- cumsum(x.prob) # Find cdf grid
  indices <- findInterval(q, x)+1

  p <- x.cdf[indices] # Index q from the cdf grid
  p <- ifelse(indices == length(x.prob), 1, p) # Fix precision issue

  if(!lower.tail){
   p <-  1 - p
  }

  if(log.p){
    p <- log(p)
  }

  return(p)
}

#' @rdname alloT
#' @export
qalloT<-function(p, mean.phi.T = NULL, min.phi = 0.005, max.phi = 1, a = 14/9, lambda = NULL,
                 n.sim = 1000, lower.tail = TRUE, log.p = FALSE){
  if(log.p){
    p <- exp(p)
  }

  if(!lower.tail){
    p <- 1 - p
  }

  if(is.null(mean.phi.T)){
    if(is.null(lambda)){
      stop("Either 'lambda' or 'mean.phi.T must be specified'.")
    }
    mean.phi.T <- get_mean_phi_T(lambda,
                                 min.phi = min.phi,
                                 max.phi = max.phi,
                                 a = a)
  }

  phi.T<-.allometry.herb.quasi.sim(mean.phi.T = mean.phi.T,
                                  n.sim = n.sim,
                                  min.phi = min.phi,
                                  max.phi = max.phi,
                                  a = a,
                                  truncate = TRUE)
  q <- ifelse(
    p > 1 | p < 0,
    NaN,
    quantile(phi.T, probs = p, names = FALSE)
  )
  return(q)
}



#' @title Calculate the kth Convolution Probability Matrix of P(phi) Via Gaussian Approximation
#' @description Internal function of \code{.dalloT.cond.k.fft.conv()} used to calculate \eqn{P(\phi^k)} via Gaussian approximation for very large values of \eqn{k}.
#' @details
#' Owing to CLT,
#' \deqn{P(\phi_T^k) \approx \mathcal{N}(\mu = k\overline\phi,\sigma=\sqrt{k\mathbb{Var}[\phi]})}
#' Also see \code{?get_phi_bar()} and \code{?get_phi_var()}
#' @return a length(phi.T) X length(k) matrix of probabilities.
#' @param phi.T a vector of cumulative proportion herbivory
#' @param k a vector of integer values indicating the number of convolutions to perform
#' @param min.phi the minimum bite size in terms of proportion leaf herbivory. Defaults to 0.005 (0.5%).
#' @param max.phi the maximum bite size in terms of proportion leaf herbivory Defaults to 1 (100%).
#' @param a the combined allometric scaling coefficient. Defaults to 14/9.
#' @param log if \code{TRUE} (default is \code{FASLE}), return probabilities on the log scale.
#' @export
.dalloT.cond.k.gauss.approx<-function(phi.T,k,min.phi=0.005,max.phi=1,a=14/9,log=FALSE){
  phi.bar<-get_phi_bar(min.phi = min.phi, max.phi = max.phi, a = a)
  phi.var<-get_phi_var(min.phi = min.phi, max.phi = max.phi, a = a)

  prob.mat<-(outer(phi.T, k, Vectorize(FUN = function(phi.T, k){
    dnorm(x = phi.T,mean = k*phi.bar,sd = sqrt(k*phi.var),log = log)
  }))) # Returns a matrix. Row is phi.T and column is k.
  return(prob.mat)
}


#' @title Perform k Convolutions of Fourier Transformed Vector
#' @description Internal function of \code{.dalloT.cond.k.fft.conv()} used to perform \eqn{k} convolutions of a supplied Fourier transformed vector. Parallel computing is supported and requires an already registered cluster.
#' @details
#' \deqn{P(\phi_T^k) = \overbrace{P(\phi)*P(\phi)*...*P(\phi)}^{\text{k times}} = \mathcal{F}^{-1} [\mathcal{F}[P(\phi)]^k]}
#' @param k a vector of integers indicating the number of convolutions to apply to a Fourier transformed vector
#' @param fft.vec a Fourier transformed vector on the frequency domain
#' @param parallel if \code{TRUE} (default is \code{FALSE}), convolve \code{fft.vec} in parallel. Automatically turned off if \code{length(ff.vec)} is less than 50000, as no computational efficiency would be gained.
#' @return A \code{length(fft.vec)} x \code{length(k)} matrix of conditional probabilities
#' @export
.convolve.dist <- function(k, fft.vec, parallel = FALSE){
  fft.vec.length<-length(fft.vec)
  if(parallel &&
     (50000 < fft.vec.length)){ # Only long fft.vec benefits from parallel
    cond.prob.mat <- foreach::foreach(i=k, .combine = cbind) %dopar%{
      (Re(fft(fft.vec^i,inverse = T))/fft.vec.length)
    }
  } else {
    cond.prob.mat<-sapply(k,function(x){
      (Re(fft(fft.vec^x,inverse = T))/fft.vec.length)
    })
  }
  return(cond.prob.mat)
}

#' @title Numerically Estimate Conditional Probability Matrix of P(phi.T | k)
#' @description Internal function of \code{dalloT()} used to estimate the conditional probability matrix \eqn{P(\phi_T|k)} via FFT or Gaussian approximation.
#' @param phi.T a vector of cumulative proportion herbivory
#' @param k a vector of integers indicating the number of convolutions to apply to a Fourier transformed vector
#' @param fft.vec a Fourier transformed vector on the frequency domain, specifically that of the probability density of phi.
#' @param phi a numeric vector of the proportion herbivory grid used to generate the \code{fft.vec} with \code{dallo()}. Supplied to find closest matched indices for \code{phi.T} for approximate conditional probabilities of \eqn{P(\phi_T | k)}.
#' @param log if \code{TRUE}, returns probabilities on the log scale.
#' @param min.phi the minimum bite size in terms of proportion leaf herbivory. Defaults to 0.005 (0.5%).
#' @param max.phi the maximum bite size in terms of proportion leaf herbivory Defaults to 1 (100%).
#' @param a the combined allometric scaling coefficient. Defaults to 14/9.
#' @param k.fft.limit the maximum number of convolutions performed by FFT. For k convolutions above this limit, a Gaussian approximation is used.
#' @param parallel if \code{TRUE} (default is \code{FALSE}), convolve \code{fft.vec} in parallel. Automatically turned off if \code{length(ff.vec)} is less than 50000, as no computational efficiency would be gained.
#' @param phi.T.index a vector of integers indicating the location where phi.T most closely match the grid of probabilities. If \code{NULL} (default), the indices are calculated using \code{phi.T}, \code{min.phi}, \code{max.phi}, and \code{a}. Supplying a pre-calculated index vector saves a lot of computation when the function is used repeatedly.
#' @return a \code{length(phi)} X \code{length(k)} matrix of conditional probabilities.
#' @export
.dalloT.cond.k.fft.conv<-function(phi.T, k, fft.vec, phi, log = FALSE,
                                 min.phi = 0.005, max.phi = 1,
                                 a = 14/9, k.fft.limit = 100,
                                 parallel = FALSE,phi.T.index = NULL){
  if(max(k)>k.fft.limit){
    message("Underflow in FFT; switching to gaussian approximation for k > ",k.fft.limit)
    k1 <- seq(min(k),k.fft.limit, by = 1)
    k2 <- seq((k.fft.limit+1),max(k),by = 1)
    # The repeated convolutions approaches a gaussian distribution when k is high
    # Automatically switch to gaussian approximation above some set threshold.
    # 100 seems to be a good limit
    # Might be able to save some computation by feeding in observed data directly to
    # get the prob, rather than searching in a big matrix
    # The approximation is not needed for biologically reasonable parameter space (!!)
    prob.max.gauss.approx<-.dalloT.cond.k.gauss.approx(phi.T = phi,
                                                      k = k2,
                                                      min.phi = min.phi,
                                                      max.phi = max.phi,
                                                      a = a,
                                                      log = FALSE)


    cond.prob.mat<-.convolve.dist(k = k1,
                                 fft.vec =  fft.vec,
                                 parallel = parallel)[seq_along(phi),]
    # Values over the upper boundary are cut off with the indexing

    cond.prob.mat<-cbind(cond.prob.mat,prob.max.gauss.approx)
  } else {
    cond.prob.mat<-.convolve.dist(k = k,
                                 fft.vec =  fft.vec,
                                 parallel = parallel)[seq_along(phi),]
    # Returns a matrix. Row is phi.T and column is k.
    # Values over the upper boundary are cut off with the indexing
  }
  cond.prob.mat[1,] <- 0 # P(phi_T = 0 | k > 0) = 0 Discrete probability
  cond.prob.mat[1,1] <- 1 # P(phi_T = 0 | k = 0) = 1 Discrete probability


  # Cut off probabilities added back to phi=1
  cond.prob.mat[nrow(cond.prob.mat),]<-cond.prob.mat[nrow(cond.prob.mat),]+1-colSums(cond.prob.mat)
  #cond.prob.mat <- cond.prob.mat/ colSums(cond.prob.mat)

  if(is.null(phi.T.index)){
    #Check if there is pre-calculated indices
    phi.T.index<-findInterval(x = phi.T, vec = phi)
  }

  # For each pairwise combination of phi.T and k, index the corresponding probability and returns a matrix
  prob.mat<-outer(phi.T.index, k,
                  Vectorize(
                    FUN = function(phi.T.index, k){
                      cond.prob.mat[
                        phi.T.index,
                        (k+1)]
                    })) # slow
  prob.mat[(phi.T>1)|(phi.T<0)|(phi.T > 0 & phi.T < min.phi),]<-0 # Set values outside of bounds to zero

  prob.mat[prob.mat<0]<-0 # Fix some numeric errors
  if(log){
    prob.mat<-log(prob.mat)
  }
  return(prob.mat) # Out puts a matrix of P(phi.T | k)
  # phi.T is the rows and k is the columns
}


#' @title Numerically Estimate Conditional Probability Matrix of P(k | lambda)
#' @description Internal function of \code{dalloT()} used to estimate the conditional probability matrix  \eqn{P(k|\lambda)}. Basically a wrapper function for \code{dpois()}.
#' @param mean.phi.T The mean herbivory of the distribution when herbivores are not plant limited.
#' @param k a vector of integers indicating the number of convolutions to apply to a Fourier transformed vector
#' @param min.phi the minimum bite size in terms of proportion leaf herbivory. Defaults to 0.005 (0.5%).
#' @param max.phi the maximum bite size in terms of proportion leaf herbivory Defaults to 1 (100%).
#' @param a the combined allometric scaling coefficient. Defaults to 14/9.
#' @param log if \code{TRUE}, returns probabilities on the log scale.
#' @param lambda if \code{NULL} (default), lambda would be estimated from the parameters given. Supplying \code{lambda} is useful to avoid repeated computation when the function is used repeatedly.
#' @param k.max.tolerance the tolerance threshold of maximum convolution cut off (ideally probabilities above \code{k.max} convolutions is vanishingly small). See details of \code{?dalloT()}
#' @return a \code{(length(phi.T))} X \code{length(k)} matrix of conditional probabilities.
#' @export
.dk.cond.lambda<-function(k,mean.phi.T=NULL,min.phi=NULL,max.phi=NULL,
                         a=NULL,k.max.tolerance=1e-5,lambda=NULL,log=FALSE){
  if(is.null(lambda)){
    lambda <- get_lambda(mean.phi.T,
                                 min.phi = min.phi,
                                 max.phi = max.phi,
                                 a = a)
  }
  prob<-dpois(x = k,lambda = lambda)
  if((prob[length(prob)]>k.max.tolerance)||
     (prob[length(prob)]>prob[length(prob)-1])){
    warning("Approximation may be unreliable; increase k.max")
  }
  if(log){
    prob<-log(prob)
  }

  return(prob) # Out puts a vector of P(k | lambda)
}

#' @title Fit Neutral Herbivory Distribution Using Maximum Likelihood Estimation
#' @description Estimate unknown parameter(s) in the neutral herbivory model from a vector of observed herbivory data using Maximum Likelihood Estimation (MLE). The likelihood function is of the form
#' \deqn{\mathcal{L}(\phi_{T1},\phi_{T2},...,\phi_{Tn}|\overline{\phi_T},\phi_M,\phi_m,\alpha)=\prod_i^n P(\phi_{Ti}|\overline{\phi_T},\phi_M,\phi_m,\alpha).}
#' @details
#' The order of \code{optim.var}, \code{upper}, \code{lower}, and \code{init} must match exactly. The initial value and bounds must be on the scale of the transformed variable.
#' @param data.vec A vector of numeric data. \code{NA}s are ignored. If the whole vector does not contain any non-zero values or contains values outside of \eqn{[0,1]}, the function would throw an error.
#' @param optim.vars A vector of character string indicating the name of the variables to be estimated. Acceptable values are "mean.phi.T" (default), "min.phi", "max.phi", and "a".
#' @param init A numeric vector of initial values in the order of optim.vars for the optimizer to start searching. Defaults to 0.
#' @param method The method of maximum likelihood estimation algorithm. Acceptable methods are "Nelder-Mead" (Default), "BFGS", "L-BFGS-B", "Brent", "nlminb", and "nlm". Defaults to "Brent" for one dimensional optimization. For details see \code{?stats::optim()}, \code{?stats::nlminb()}, and \code{?stats::nlm()}
#' @param k.max The maximum number of convolutions of the neutral 'bite size' distribution in numerical approximation of the PDF of the neutral herbivory distribution. Default is 50.
#' @param by The grid resolution used in the FFT convolutions. Default is 0.001.
#' @param k.max.tolerance the tolerance threshold of maximum convolution cut off (ideally probabilities above \code{k.max} convolutions is vanishingly small).
#' @param k.fft.limit The maximum number of convolutions performed by FFT. For k convolutions above this limit, a Gaussian approximation is used. See details for more information.
#' @param param.vals A named vector of the default parameter values of the neutral herbivory distribution. Set to \code{NA} if the variable is declared in the \code{optim.vars} argument (i.e. is estimated). "mean.phi.T" defaults to 0.1, "min.phi" to 0.005, "max.phi" to 1, and "a" to 14/9.
#' @param param.val.trans A vector of functions that transform each estimated parameter value in the optimization process (the optimizer find the values on the transformed scale). Set to \code{NA} if the variable is not declared in the \code{optim.vars} argument (i.e. is not estimated). Default transformation is x/100 for "mean.phi.T", "min.phi", and "max.phi". Default transformation for "a" is identity.
#' @param lower,upper A numeric vector indicating the bounds of each estimated variable for the "L-BFGS-B" and "Brent" method. Defaults to -Inf and Inf. If the length of supplied vector is shorter than the number of estimated variables, the first element of the vector is set as the bound for all variables. If the method is "Brent" and the \code{optim.vars = "mean.phi.T"}, the bounds default to 0 and 100.
#' @param cores The number of parallel processes in parallel computing.
#' @param id A character string supplied for book keeping purposes. Default is \code{NULL}. Useful for when storing a large number of "fit_allo_herb" objects in a list.
#' @param ... Additional arguments that are passed to the optimizer functions \code{optim()}, \code{nlminb()}, or \code{nlm()}.
#' @return
#' An 'allo_herb_fit' object with the following slots:
#'
#' \code{theta.names}: A vector of character string naming the fitted parameters
#'
#' \code{par}: the value of fitted parameters on the transformed scale in the order of theta.names
#'
#' \code{se}: standard error of fitted parameters on the transformed scale. Quadratic approximation of the standard error will be implemented in the future.
#'
#' \code{loglik}: log likelihood of the data given the parameter combinations
#'
#' \code{hessian}: The Hessian matrix
#'
#' \code{iters}: The number of times the likelihood function or its gradient is evaluated.
#'
#' \code{message}: message passed from the optimizer function
#'
#' \code{convergence}: 0 indicates successful convergence. Error codes are passed from the optimizer function. 1 indicate that the iteration limit has been reached and convergence has failed. 2 indicate that the Hessian matrix is non-positive definite.
#'
#' \code{method}: Method of optimization
#'
#' \code{init}: A vector of the used initial values for optimization
#'
#' \code{data}: The data to which the model is fitted to
#'
#' \code{param.val.trans}: A vector of functions that transform each estimated parameter value in the optimization process (the optimizer find the values on the transformed scale).
#'
#' \code{param.vals}: the default parameter values assumed for the model
#'
#' \code{id}: name of model supplied via the id argument
#'
#' \code{n}: the sample size
#'
#' \code{df}: the degrees of freedom of the model
#'
#' \code{herbivar.version}: version of herbivar used to fit the model
#'
#' \code{num.approx.param}: a named vector of values indicating the parameters used in the numerical approximation of the density function of the neutral herbivory distribution.
#' @export
fit_allo_herb<-function(data.vec,
                             optim.vars=c("mean.phi.T"),
                             init = 0,
                             method=c("Nelder-Mead","BFGS",
                                      "L-BFGS-B","Brent","nlminb","nlm"),
                             k.max = 50,
                             by = 0.001,
                             k.max.tolerance = 1e-5,
                             k.fft.limit = 100,
                             param.vals=c("mean.phi.T" = 0.1,
                                          "min.phi" = 0.005,
                                          "max.phi" = 1,
                                          "a" = 14/9),
                             param.val.trans = c(
                               "mean.phi.T" = function(x) {x/100},
                               "min.phi" = function(x) {x/100},
                               "max.phi" = function(x) {x/100},
                               "a" = function(x) {x}),
                             lower = -Inf,
                             upper = Inf,
                             cores = 1,
                             id = NULL, # A name for book keeping purposes
                             ...){
  if(!all(optim.vars %in% names(param.vals))){
     stop("'",
          paste0(optim.vars[!optim.vars %in% names(param.vals)], collapse = " and "),
          "' not valid model parameter")
  }

  if(length(lower)!=length(optim.vars)){
    lower<-rep(lower[1],length(optim.vars))
  }
  if(length(upper)!=length(optim.vars)){
    upper<-rep(upper[1],length(optim.vars))
  }

  method <- match.arg(method)

  if((length(optim.vars)<2) & (method!= "Brent") & (method != "nlminb")){
    method<-"Brent"
    message("Overwriting optimizer to Brent for one dimensional optimization")
  }

  if(length(init)!=length(optim.vars)){
    stop("Initial values length not equal to number of fitted parameters")
  }

  param.val.trans <- lapply(param.val.trans, function(x){
    match.fun(x)
  })

  param.vals[(names(param.vals)%in%optim.vars)]<-NA # Fitted values set to NA

  if(length(optim.vars)!=sum(is.na(param.vals))){
    stop("Number of fitted parameters does not equal to the number of unknown parameters")
  }

  if(method == "Brent" &&
     (optim.vars == "mean.phi.T")) {
    if(is.infinite(lower)){
      lower <- 0
    }
    if(is.infinite(upper)){
      upper <- 100
    }
  }

  min.phi <- param.vals["min.phi"]
  max.phi <- param.vals["max.phi"]
  lower.min.phi <- param.val.trans[["min.phi"]](lower[optim.vars == "min.phi"])
  lower.max.phi <- param.val.trans[["max.phi"]](lower[optim.vars == "max.phi"])
  upper.min.phi <- param.val.trans[["min.phi"]](upper[optim.vars == "min.phi"])

  data.vec <- .herb_data_check(data.vec,
                               min.phi = ifelse(is.na(min.phi),
                                                lower.min.phi,
                                                min.phi),
                               allow.zero = TRUE)

  if(method == "L-BFGS-B" || method == "Brent"){
    if(("max.phi"%in%optim.vars) &&
       !is.na(min.phi) &&
       (lower.max.phi < min.phi)){
      stop("Lower bound of max.phi needs to be higher than min.phi")
    }
    if(("min.phi"%in%optim.vars) &&
       !is.na(max.phi) &&
       (upper.min.phi > max.phi)){
      stop("Upper bound of min.phi needs to be higher than max.phi")
    }
  }


  if(by > 0.01){
    warning("Approximation resolution too low; results are crappy")
  } else if((!is.na(param.vals["max.phi"]) && !is.na(param.vals["min.phi"])) &&
            (param.vals["max.phi"]-param.vals["min.phi"] < by)){
    stop("Approximation resolution too low; set 'by' to a lower number")
  }

  if(cores > 1){
    parallel <- TRUE
    cluster <- makeCluster(cores)
    doParallel::registerDoParallel(cluster)

  } else {
    parallel <- FALSE
  }

  nll <- function(theta) {
    nll.out<-sum(dalloT(x = data.vec,
                        mean.phi.T = ifelse(
                          is.na(param.vals["mean.phi.T"]),
                          param.val.trans[["mean.phi.T"]](eval(parse(
                            text =
                              paste0("theta[",
                                     which(optim.vars=="mean.phi.T"),"]")))),
                          param.vals["mean.phi.T"]),
                        min.phi = ifelse(
                          is.na(param.vals["min.phi"]),
                          param.val.trans[["min.phi"]](eval(parse(
                            text =
                              paste0("theta[",
                                     which(optim.vars=="min.phi"),"]")))),
                          param.vals["min.phi"]),
                        max.phi = ifelse(
                          is.na(param.vals["max.phi"]),
                          param.val.trans[["max.phi"]](eval(parse(
                            text =
                              paste0("theta[",
                                     which(optim.vars=="max.phi"),"]")))),
                          param.vals["max.phi"]),
                        a = ifelse(
                          is.na(param.vals["a"]),
                          param.val.trans[["a"]](eval(parse(
                            text =
                              paste0("theta[",
                                     which(optim.vars=="a"),"]")))),
                          param.vals["a"]),
                        k.max = k.max,
                        by = by,
                        k.max.tolerance = k.max.tolerance,
                        k.fft.limit = k.fft.limit,
                        log = TRUE,
                        parallel = parallel))*-1
    return(nll.out)
  }

    ML.fit<-optim2(init = init,
                  fn = nll,
                  hessian = TRUE,
                  method = method,
                  lower = lower,
                  upper = upper,
                  ...)
    hessian <- ML.fit$hessian
    loglik <- -ML.fit$value
    iters <- ML.fit$counts

  allo.herb.fit.out<-list("theta.names" = optim.vars,
                          "par" = ML.fit$par,
                          "se" = tryCatch(sqrt(diag(solve(hessian))),
                                          error = function(e) {
                                            rep(NA,length(diag(hessian)))
                                          }),
                          "loglik" = loglik,
                          "hessian" = hessian,
                          "message" = ML.fit$message,
                          "iters" = iters,
                          "convergence" = ML.fit$convergence, # > 0 indicates issues
                          "method" = method,
                          "init" = init,
                          "data" = data.vec,
                          "param.vals" = param.vals,
                          "param.val.trans" = param.val.trans,
                          "id" = id,
                          "n" = length(data.vec),
                          "df" = length(optim.vars),
                          "herbivar.version" = herbivar.version(silent = TRUE),
                          "num.approx.param" = c("k.max" = k.max,
                                                 "resolution" = by,
                                                 "k.max.tolerance" = k.max.tolerance,
                                                 "k.fft.limit" = k.fft.limit)

  )
  if(any(is.na(allo.herb.fit.out$se))){
    warning("Non-positive definite hessian matrix")
    #Usually increasing 'by' resolution fixes this issue. But SLOW!
    allo.herb.fit.out$convergence <- 2
  } else if(allo.herb.fit.out$convergence!=0){
    warning("Model did not converge","  ",allo.herb.fit.out$convergence)
  }

  allo.herb.fit.out<-structure(allo.herb.fit.out,
                               class=c("allo_herb_fit","list"))
  return(allo.herb.fit.out)

  on.exit(
    try({
      if(!is.null(cluster)){
        doParallel::stopImplicitCluster()
        stopCluster(cluster)
      }
    })
  ) # Stop all connections on exit
  # Use closeAllConnections() to kill zombie workers
}

#' @title Print Values
#' @description Prints fitted object and returns an invisible matrix array of the fitted model coefficients
#' @param x An object of class 'allo_herb_fit'
#' @param ... additional arguments
#' @param trans A string or function of transformation applied on the model coefficients
#' @param digits A numeric value indicating the number of digits to be displayed in the model coefficients. Default is 3.
#' @return A matrix array
#' @export
print.allo_herb_fit <- function(x, ..., trans = NULL, digits = 3){
  if(!inherits(x,"allo_herb_fit")){
    stop("Needs to by of object type 'allo_herb_fit'")
  }
  if(is.null(trans)){
    FUN <- function(x) x
  } else {
    FUN <- match.fun(trans)
  }

  cat("Fitted result on", x$id,"\n")
  cat("method =",x$method,"\n")
  cat("n =",length(x$data),"\n")
  cat("\n")
  print(c("loglik" = x$loglik, "AIC" = AIC(x), "AICc" = AICc(logLik(x)), "BIC" = BIC(x)))
  cat("\n")
  print(c(x$iters,"convergence"=x$convergence))
  cat("\n")
  print(noquote(vapply(
    x$param.val.trans[x$theta.names],
    function(z){
      gsub(".*\\{| |\\}","",paste0(deparse(z),collapse = ""))
    },
    character(1)
  )))
  cat("\n")
  x$param.vals<-round(x$param.vals,digits = 3)
  x$param.vals[is.na(x$param.vals)]<-"fitted"
  print(noquote(x$param.vals))
  cat("\n")
  out<-round(
    apply(as.matrix(
      data.frame(row.names = x$theta.names,
                 "Estimate" = x$par,
                 "Std." = x$se,
                 "lower" = x$par-x$se*1.96,
                 "upper" = x$par+x$se*1.96)
    ), 2, FUN = FUN),
    digits = digits)
  if(is.vector(out)){
    out<-t(as.matrix(out))
    rownames(out)<-x$theta.names
  }
  print(out)
}

#' @title Extract Log-Likelihood
#' @description extract log-likelihood
#' @param object a fitted model object
#' @param ... additional arguments
#' @return a 'logLik' object
#' @rdname logLik
#' @export
logLik.allo_herb_fit <- function(object, ...){
  out <- object$loglik
  out <- structure(out,
                   class= "logLik",
                   nall = length(object$data),
                   nobs= length(object$data),
                   df = length(object$theta.names))
  return(out)
}

#' @title Akaike's An Information Criterion
#' @description calculate AIC or BIC for fitted model objects. see \code{stats::AIC()} for more details.
#' @param object a fitted model object for which there exist a \code{logLik} method to extract the log-likelihood, or an object of class 'logLik'.
#' @param ... additional arguments passed to \code{stats::AIC()} or \code{stats::BIC()} for additional fitted model objects.
#' @param k an atomic numeric value for the penalty. Defaults to 2.
#' @return if just one object is provided, a numeric value of AIC or BIC is returned. If multiple objects are provided, a data.frame with rows corresponding to the objects and columns representing the number of parameters in the model and the AIC or BIC is returned.
#' @rdname AIC
#' @export
AIC.allo_herb_fit <- function(object, ..., k = 2){
  dots <- as.list(match.call())[-1]
  dots <- dots[!names(dots) == "k"]

  if(length(dots) == 1){
    return(AIC(logLik(object),..., k = k))
  } else {
    out <- do.call("rbind",
                   lapply(unname(dots), function(x){
                     loglik <- logLik(eval(x))
                     df <- attributes(loglik)$df
                     aic <- AIC(loglik)
                     return(data.frame("model"= deparse(x), "AIC" = aic, "df"= df))
                   }))

    out$dAIC <- out$AIC - min(out$AIC)
    out <- out[order(out$AIC),]
    return(out)
  }
}

#' @title Corrected Akaike's An Information Criterion (AICc)
#' @description Calculate AICc taking into account small sample size bias
#' @param object The object from which log likelihood can be extracted
#' @param ... additional fitted model objects.
#' @param k the penalty on model complexity. Default is 2.
#' @details
#' The formulat for AICc is given by:
#' \deqn{AIC_c = -2\ln{\mathcal{L}(x)} + k  d \frac{n}{n - d - 1}}
#' \eqn{d} is the number of parameters, \eqn{k} is the penalty, \eqn{n} is the sample size, and \eqn{\mathcal{L}(x)} is the likelihood.
#'
#' @return if just one object is provided, a numeric value of AICc is returned. If multiple objects are provided, a data.frame with rows corresponding to the objects and columns representing the number of parameters in the model and the AICc is returned.
#' @export
AICc <- function(object, ..., k = 2){
  dots <- as.list(match.call())[-1]
  dots <- dots[!names(dots) == "k"]

  if(length(dots) == 1){
    logLik <- logLik(object)
    n <- attributes(logLik)$nobs
    df <- attributes(logLik)$df
    aicc<-as.numeric(logLik) * -2 + k*df*(n/(n-df-1))
    if(n <= (df+1)){
      warning("Sample size too small.")
    }
    return(aicc)
  } else {
    out <- do.call("rbind",
                   lapply(unname(dots), function(x){
                     loglik <- logLik(eval(x))
                     n <- attributes(loglik)$nobs
                     df <- attributes(loglik)$df
                     aicc<-as.numeric(loglik) * -2 + k*df*(n/(n-df-1))
                     if(n <= (df+1)){
                       warning("Sample size too small.")
                     }
                     return(data.frame("model"= deparse(x), "AICc" = aicc, "df"= df))
                   }))

    out$dAICc <- out$AICc - min(out$AICc)
    out <- out[order(out$AICc),]
    return(out)
  }
}


#' @title Bayesian information criterion
#' @rdname AIC
#' @export
BIC.allo_herb_fit <- function(object,...){
  dots <- as.list(match.call())[-1]

  if(length(dots) == 1){
    return(BIC(logLik(object),...))
  } else {
    out <- do.call("rbind",
                   lapply(unname(dots), function(x){
                     loglik <- logLik(eval(x))
                     df <- attributes(loglik)$df
                     bic <- BIC(loglik)
                     return(data.frame("model"= deparse(x), "BIC" = bic, "df"= df))
                   }))

    out$dBIC <- out$BIC - min(out$BIC)
    out <- out[order(out$BIC),]
    return(out)
  }
}

#' @title Extract Model Coefficients
#' @description extract model coefficients from allo_herb_fit objects
#' @param object An 'allo_herb_fit' object
#' @param ... additional arguments
#' @param backtransform A logic value indicating whether to back transform coefficients from the scale the coefficient was estimated at to the scale the coefficient is parameterized as for the neutral model. Default is TRUE.
#' @return A named vector
coef.allo_herb_fit<-function(object, ..., backtransform = TRUE){
  if(backtransform){
    out<-vapply(seq_along(object$par), function(i){
      vapply(object$par[i],object$param.val.trans[[object$theta.names[i]]],numeric(1))
    }, numeric(1))
  } else {
    out <- object$par
  }
  names(out) <- object$theta.name
  return(out)
}




#' @title Calculate Attack Rate From Unlimited Mean Cumulative Proportion Herbivory
#' @description Convenience function that reparameterizes unlimited mean cumulative proportion herbivory \eqn{\overline{\phi_T'}} as attack rate \eqn{\lambda}.
#' @details
#' \deqn{\lambda = \frac{\overline{\phi_{T}'}}{\overline{\phi}}}
#' @param mean.phi.T The mean herbivory of the distribution when herbivores are not plant limited.
#' @param min.phi the minimum bite size in terms of proportion leaf herbivory. Defaults to 0.005 (0.5%).
#' @param max.phi the maximum bite size in terms of proportion leaf herbivory Defaults to 1 (100%).
#' @param a the combined allometric scaling coefficient. Defaults to 14/9.
#' @return a vector of numeric values
get_lambda <- function(mean.phi.T, min.phi = 0.005, max.phi = 1, a = 14/9){
  mean.phi.T / get_phi_bar(min.phi = min.phi, max.phi = max.phi, a = a)
}

#' @title Calculate Unlimited Mean Cumulative Proportion Herbivory From Attack Rate
#' @description Convenience function that reparameterizes attack rate \eqn{\lambda} as unlimited mean cumulative proportion herbivory \eqn{\overline{\phi_T'}}.
#' @details
#' \deqn{\overline{\phi_T'} = \overline{\phi}\lambda}
#' @param lambda the attack rate on the plant or leaf
#' @param min.phi the minimum bite size in terms of proportion leaf herbivory. Defaults to 0.005 (0.5%).
#' @param max.phi the maximum bite size in terms of proportion leaf herbivory Defaults to 1 (100%).
#' @param a the combined allometric scaling coefficient. Defaults to 14/9.
#' @return a vector of numeric values
get_mean_phi_T <- function(lambda, min.phi = 0.005, max.phi = 1, a = 14/9){
  lambda * get_phi_bar(min.phi = min.phi, max.phi = max.phi, a = a)
}

#' @title Calculate Mean Proportion Bite Size
#' @description Calculate the theoretical mean of the proportion bite size distribution
#' @param min.phi the minimum bite size in terms of proportion leaf herbivory. Defaults to 0.005 (0.5%).
#' @param max.phi the maximum bite size in terms of proportion leaf herbivory Defaults to 1 (100%).
#' @param a the combined allometric scaling coefficient. Defaults to 14/9.
#' @details
#'
#' The mean of the bite size distribution can be found via:
#'
#' When \eqn{\phi \neq 1,2}
#' \deqn{\overline{\phi} = \frac{1-\alpha}{2-\alpha} \frac{\phi_M^{2-\alpha} - \phi_m^{2-\alpha}}{\phi_M^{1-\alpha} - \phi_m^{1-\alpha}}}
#'
#' when \eqn{\phi = 1}
#' \deqn{\overline{\phi} = \frac{1}{2-\alpha} \frac{\phi_M^{2-\alpha} - \phi_m^{2-\alpha}}{\ln{\phi_M} - \ln{\phi_m}}}
#'
#' when \eqn{\phi = 2}
#' \deqn{\overline{\phi} = (1-\alpha) \frac{\ln{\phi_M} - \ln{\phi_m}}{\phi_M^{1-\alpha} - \phi_m^{1-\alpha}}}
#'
#' @return a numeric value of the mean
#'
get_phi_bar <- function(min.phi = 0.005, max.phi = 1, a = 14/9){
  if(a == 1){
    phi.bar <- 1 / (2-a) * (max.phi^(2-a) - min.phi^(2-a))/(log(max.phi)-log(min.phi))
  } else if(a == 2){
    phi.bar <- (1-a) * (log(max.phi) - log(min.phi))/(max.phi^(1-a) - min.phi^(1-a))
  } else {
    phi.bar <- (1-a)/(2-a) * (max.phi^(2-a) - min.phi^(2-a))/(max.phi^(1-a) - min.phi^(1-a))
  }
  return(phi.bar)
}

#' @title Calculate Variation of Proportion Bite Size
#' @description Calculate the theoretical variance of the proportion bite size distribution
#' @param min.phi the minimum bite size in terms of proportion leaf herbivory. Defaults to 0.005 (0.5%).
#' @param max.phi the maximum bite size in terms of proportion leaf herbivory Defaults to 1 (100%).
#' @param a the combined allometric scaling coefficient. Defaults to 14/9.
#' @details
#'
#' The variance of the bite size distribution can be found via:
#'
#' When \eqn{\phi \neq 1,3}
#' \deqn{\mathbb{V}ar[\phi] = \frac{1-\alpha}{3-\alpha} \frac{\phi_M^{3-\alpha} - \phi_m^{3-\alpha}}{\phi_M^{1-\alpha} - \phi_m^{1-\alpha}} - \overline{\phi}^2}
#'
#' when \eqn{\phi = 1}
#' \deqn{\mathbb{V}ar[\phi] = \frac{1}{3-\alpha} \frac{\phi_M^{3-\alpha} - \phi_m^{3-\alpha}}{\ln{\phi_M} - \ln{\phi_m}} - \overline{\phi}^2}
#'
#' when \eqn{\phi = 3}
#' \deqn{\mathbb{V}ar[\phi] = (1-\alpha) \frac{\ln{\phi_M} - \ln{\phi_m}}{\phi_M^{1-\alpha} - \phi_m^{1-\alpha}} - \overline{\phi}^2}
#'
#' @return a numeric value of the variance
#'
get_phi_var <- function(min.phi = 0.005, max.phi = 1, a = 14/9){
  phi.bar <- get_phi_bar(min.phi = min.phi, max.phi = max.phi, a = a)
  if(a == 1){
    phi.var <- 1 / (3-a) * (max.phi^(3-a) - min.phi^(3-a))/(log(max.phi)-log(min.phi)) - phi.bar^2
  } else if(a == 3){
    phi.var <- (1-a) * (log(max.phi) - log(min.phi))/(max.phi^(1-a) - min.phi^(1-a)) - phi.bar^2
  } else {
    phi.var<-(1-a)/(3-a)*(max.phi^(3-a)-min.phi^(3-a))/(max.phi^(1-a)-min.phi^(1-a))-phi.bar^2
  }
  return(phi.var)
}























