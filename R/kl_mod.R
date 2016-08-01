# We need to modify 


# Parameters
p =  2
m = 1.5
Goo = 0.78 # Can't be neg
Eoo = 5  # Can't be neg
beta.o = 0.8
beta.m = 0.2 
Gmo = 0
Emo = 0
Gmm = 
KLmod <- function(
m,		# maternal effect on offspring trait
Goo, 	# additive genetic variance (of offspring trait?)
Eoo,	#1 - m - Goo # environmental variance of offspring trait
beta.o,	# selection coefficient on offspring trait

p,		# maternal effect on maternal effector
Gmm,	# additive genetic variance (of mother trait?)
#Emm,	# environmental variance of offspring trait
beta.m,	# selection coefficient on mother trait

Gmo, # genetic covariance (between maternal and offspring trait)
Emo,

gens, # number of generations
gens.sel) # generation at which selection stops
{


Emm <-  1 - Gmm - p^2 - 2*p*Gmm/(2-p) #1 - p^2 - Gmm,		#1 - p - Gmm # environmental variance of offspring trait

	
	M <- matrix(c(p,m,0,0),ncol=2)
#	M <- matrix(c(p,0,m,0),ncol=2)


	G <- matrix(c(Gmm,Gmo,Gmo,Goo),ncol=2) #Gom = Gmo

	E <- matrix(c(Emm,Emo,Emo,Eoo),ncol=2) #Eom = Emo

	I <- diag(1,ncol=2, nrow=2)

	Caz <- G %*% (solve(I - 0.5* t(M))) #solve creates the inverse of a matrix


 Pmo <- Pom <- ( ((m * (1 + 2*p)) / ((1 - p^2) * (2 - p))) * Gmm ) + (2/(2-p))*Gmo + Emo + ( (m*p)/(1-p^2) )*Emm

 Pmm <- ((2+p)*Gmm + (2-p)*Emm) / ((1-p^2)*(2-p))

 Poo <- Goo + Eoo + ((2*m)/(2-p))*Gmo + (((m^2)*(2+p))/((2-p)*(1-p^2)))*Gmm + ((m^2)/(1-p^2))*Emm

 P <- matrix(c(Pmm,Pmo,Pom,Poo),ncol=2)
print(P)


#	P <- 	G + E + M %*% P %*% t(M) + 0.5* M %*% t(Caz) + 0.5* Caz %*% t(M)


	beta <- matrix(c(c(0,0),rep(c(beta.m,beta.o), gens.sel-1),rep(c(0,0), gens-gens.sel)),ncol=gens)

	z.hat <- matrix(NA, nrow=2,ncol=gens)
	z.hat[,1] <- c(0,0)
 
 
	# MODEL
	for(i in 2:gens) z.hat[,i] <- (Caz + M %*% P) %*% beta[,i] + M %*% z.hat[,i-1] - M %*% P %*% beta[,i-1] 

	out <- apply(z.hat,1,cumsum)
	colnames(out) <- c("m","o")
	return(out)
}

