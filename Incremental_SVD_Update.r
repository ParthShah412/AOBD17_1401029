cat("\014")					# clears the command console

RankOneUpdate <- function(U,S,V,A,B)
{
	V <- rbind(V,matrix(0,1,ncol(V)));

	# eq 6
	M <- t(U) %*% A;
	P <- A - (U%*%M);
	Ra = norm(P);
	P <- (1/Ra) * P;

	# eq 7
	N <- t(V) %*% B;
	Q <- B - (V%*%N);
	Rb <- norm(Q);
	Q <- (1/Rb) * Q;

	K <- S
	K <- rbind(K,0);
	K <- cbind(K,0);

	M1 <- rbind(M,Ra);
	M2 <- rbind(N,Rb);
	K <- K + M1 %*% t(M2);

	KSVD <- svd(K);
	K_U <- KSVD$u;
	K_S <- diag(KSVD$d);
	K_V <- KSVD$v;

	U_new <- cbind(U,P);
	V_new <- cbind(V,Q);

	U_new <- U_new %*% K_U;
	V_new <- V_new %*% K_V;
	S_new <- K_S;
	result <- list(U=U_new,S=S_new,V=V_new);
	rm(KSVD,M1,M2,K,N,Q,Rb,M,P,Ra);
	return(result);
}

{
	v_h = 108;
	v_w = 192;
	p = v_h * v_w;
	FPS = 25;
	seconds = 1;
	q = FPS * seconds;
	max_rank = q/2;

	Amp = 255						# Amp defines the maximum value a matrix element can take

	sizeX = p*q;

	X <- matrix(sample.int(Amp, sizeX, replace = TRUE),nrow = p,ncol = q,byrow = TRUE)
	SVD <- svd(X);
	U <- SVD$u;
	S <- diag(SVD$d);
	V <- SVD$v;

	rm(SVD,sizeX);

	B <- matrix(0,q+1,1)
	B[q+1] = 1;

	ERR <- matrix(0,1,max_rank);

	for (I in (1:max_rank))
	{
		A <- matrix(sample.int(Amp, p, replace = TRUE),nrow = p,ncol = 1,byrow = TRUE)
		
		result <- RankOneUpdate(U,S,V,A,B);
		Un <- result$U;
		Sn <- result$S;
		Vn <- result$V;
		rm(result,U,S,V);
		U <- Un; S <- Sn; V <- Vn;
		
		X <- cbind(X,matrix(0,nrow(X),1));
		X <- X + (A%*%t(B));
		
		B <- rbind(0,B);
		
		Xn <- (Un %*% Sn)%*%t(Vn);
		ERR[I] <- norm(X-Xn);
	}
	rm(A);
	rm(Amp);
	rm(Un,Sn,Vn);	 # remove this line in order to compare obtained and desired eigens
	rm(U,S,V);	 	 # remove this line in order to compare obtained and desired eigens
	rm(X,Xn);		 # remove this line in order to compare X and Xn
	par(mar = rep(2, 4));
	plot(1:max_rank, ERR);
}
rm(B,I,max_rank,p,q);
rm(RankOneUpdate)