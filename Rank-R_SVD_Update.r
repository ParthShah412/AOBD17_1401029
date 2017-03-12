##########################################################################################
# Phase 0 : Setting up the environment and parameters.
##########################################################################################
cat("\014")					# clears the command console

v_w = 12
v_h = 3
p = v_w * v_h				# dimensions of input data (rows = h * w)

FPS = 30;                # for the sake of my laptops sanity
seconds_elapsed = 30/10;
q = FPS*seconds_elapsed;	# number of input samples to X.. let q > max_c
max_c = 5;

sizeX = p * q				# defining the sizes of X , A and B

ITERS = 5;

for (J in (1:ITERS))
{
  for (c in(1:max_c))
  {
# if c is more than q, the error will be much more than expected
# dont ever let c > q

sizeA = p * c
sizeB = q * c

Amp = 8						# Amp defines the maximum value a matrix element can take

##########################################################################################
# Phase 1: Defining the initial input X
##########################################################################################
X <- matrix(sample.int(Amp, sizeX, replace = TRUE),nrow = p,ncol = q,byrow = TRUE)
SVD <- svd(X)
U_init <- SVD$u
d_init <- SVD$d; S_init <- diag(d_init)		# save D and diagonal matrix of D
V_init <- SVD$v
rm(SVD, d_init)

##########################################################################################
# Phase 2: Defining the new samples A and B with sizes defined as in paper
##########################################################################################
A <- matrix(sample.int(Amp, sizeA, replace = TRUE),nrow = p,ncol = c,byrow = TRUE)
B <- matrix(sample.int(Amp, sizeB, replace = TRUE),nrow = q,ncol = c,byrow = TRUE)

rm(sizeX, sizeA, sizeB, Amp)		# remove variables not required now



##########################################################################################
# Phase 3: Find orthogonal bases P and Q, and find RA and RB .... equation 2 of paper
##########################################################################################

I <- diag(p) - (U_init %*% t(U_init));
IA <- I %*% A;
QR.qr <- qr(IA);
P <- qr.Q(QR.qr); 				# have successfully found P (orthonormal "Q" matrix of (I-UUT)A=QR
RA <- (t(P) %*% I) %*% A;		# have succesfuly found RA and P
rm(IA);							# remove variables not required now

I <- diag(q) - (V_init %*% t(V_init));
IB <- I %*% B
QR.qr <- qr(IB)
Q <- qr.Q(QR.qr) 				# have successfully found Q (orthonormal "Q" matrix of (I-VVT)B=QR
RB <- (t(Q) %*% I) %*% B 		# have successfully found RB and Q
rm(I, IB, QR.qr)				# remove variables not required now


##########################################################################################
# Phase 4: Find "K" and get to equation (4) of the paper
##########################################################################################

K <- S_init

for (dum in (1:c - 1))			  # pad K with zeros, to add second apart as given in paper
{
  K <- rbind(cbind(K, 0), 0);
}
rm(dum);

temp <- t(U_init) %*% A
M1 <- rbind(temp, RA)

temp <- t(V_init) %*% B
M2 <- rbind(temp, RB)

K <- K + (M1 %*% t(M2))				# is implemenataiton of equation 4
rm(M1, M2)							# remove unneeded variables, keep "temp" for next phase


##########################################################################################
# Phase 5: Calculating the OBTAINED values of U,S,V by our formula
##########################################################################################

orth_K <- svd(K)
Udash <- orth_K$u
Vdash <- orth_K$v
Sdash <- orth_K$d

S_obtained <- diag(Sdash)						# OBTAINED S MATRIX
temp <- cbind(U_init, P)
U_obtained <- temp %*% Udash					# OBTAINED U MATRIX

temp <- cbind(V_init, Q)
temp <- temp %*% Vdash
V_obtained <- (temp);
#V_obtained <- temp;								# OBTAINED V MATRIX

rm(p, q, c)
rm(K, orth_K, temp)
rm(Udash, Sdash, Vdash)
rm(P, RA, Q, RB)
rm(U_init, S_init, V_init)

X_Dash <- (U_obtained %*% S_obtained) %*% t(V_obtained);

##########################################################################################
# Phase 6: Verification by comparing desired and obtained values of U,S,V
##########################################################################################

X <- X + (A %*% t(B))
ErrorMat <- X - X_Dash;
Norm_Error <- norm(ErrorMat,'f');

desired_SVD <- svd(X)
S_desired <- desired_SVD$d

S_desired <- diag(S_desired)						# DESIRED S MATRIX
U_desired <- desired_SVD$u							# DESIRED U MATRIX
V_desired <- desired_SVD$v							# DESIRED V MATRIX
rm(desired_SVD)
}
}
c=5;
X = matrix(c(1:c));
plot(x,err)