kl_rate.R - function that gives an evolutionary rate from given inputs
(kirkpatrick and land 1989 equation 7)



# function into which you put specific parameters

kl_rate <- function(
	m,		# maternal effect on offspring trait
	Goo, 	# additive genetic variance (of offspring trait?)
	Eoo,	#1 - m - Goo # environmental variance of offspring trait
	beta.o,	# selection coefficient on offspring trait

	p,		# maternal effect on maternal effector
	Gmm,	# additive genetic variance (of mother trait?)
	Emm,	# environmental variance of offspring trait
	beta.m,	# selection coefficient on mother trait

	Gmo, # genetic covariance (between maternal and offspring trait)
	Emo){
	
	beta <- matrix(c(beta.m,beta.o))

	M <- matrix(c(p,m,0,0),ncol=2)

	G <- matrix(c(Gmm,Gmo,Gmo,Goo),ncol=2) #Gom = Gmo

	E <- matrix(c(Emm,Emo,Emo,Eoo),ncol=2) #Eom = Emo

	I <- diag(1,ncol=2, nrow=2)

	Caz <- G %*% (solve(I - 0.5* t(M))) #solve creates the inverse of a matrix


	kl.R <- solve(I-M) %*% Caz %*% beta
	return(kl.R)
	}








# function into which you put matrices

kl_rate_M <- function(M, G, E, beta){

## check if inputs are correct (i.e. matrix for M, G and E, vector for beta, both dimensions the same, all the same size)
## probably neater way of doing this
## also need some meaningful errors

	dimM <- dim(M)
	dimG <- dim(G)
	dimE <- dim(E)
	dimBeta <- length(beta)

stopifnot(length(dimM)==2, length(dimG)==2, length(dimE)==2)

stopifnot(dimM[1]==dimM[2], dimG[1]==dimG[2], dimE[1]==dimE[2])

stopifnot(dimM[1]==dimG[1], dimG[1]==dimE[1], dimE[1]==dimBeta)


	I <- diag(1,ncol=dimBeta, nrow=dimBeta)

	Caz <- G %*% (solve(I - 0.5* t(M))) #solve creates the inverse of a matrix


	kl.R <- solve(I-M) %*% Caz %*% beta
	return(kl.R)
	}





## test

M <- matrix(c(0.5,0.5,0,0), nrow=2, ncol=2)
G <- matrix(c(0.5,0,0,0.5), nrow=2, ncol=2)
E <- matrix(c(0.25,0,0,0.25), nrow=2, ncol=2) 
beta <- c(0,0.5)

kl_rate_M(M,G,E,beta)