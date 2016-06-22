
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
