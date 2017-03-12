v_w = 4;
v_h = 3;

FPS = 25;
seconds = 2/25;


# Define the samples of the input, consider 2 inputs for "training data"
samples = FPS*seconds;

# Define the dimensions of the input, consider 2-D (x,y) type input
dims = v_w * v_h;

#This product will be used to allocate matrix by length
size = samples * dims

# P is the data/matrix used to find the 
max_iter=10
# "50" max_iter defines the max number of points to evaluate over, after which code will stop
P <- matrix(sample.int(max_iter, size, replace = TRUE),nrow = samples,ncol = dims,byrow = TRUE)

X <- cbind(1,P[1,])			# Form the X matrix of rows [1 x1; 1 x2...]
Y <- cbind(P[2,])			# Form the Y matrix of rows [y1; y2; ...]

# Now, we perform calculations as expected
B <- t(X)%*%Y				# B = Xt Y
A <- t(X)%*%X				# A = Xt X

temp <- solve(A)%*%B		# temp = A(-1) B

# Add points from training dataset to plot
plot(X[,2],Y)					# add the individual points (part of online data) to the output graph

# Plot the line represented by y=ax+b
abline(a=temp[1],b=temp[2])		# For y=mx+c, temp(0) is the intercept (c) and temp(1) is the slope (m)

for (dum in (1:max_iter-1)) 
{
	# Create and store new points
	N <- matrix(sample.int(max_iter, size, replace = TRUE),nrow = samples,ncol = 1,byrow = TRUE)
	P <- cbind(P,N)
	X <- cbind(1,P[1,])
	Y <- cbind(P[2,])
	
	# Perform necessary calculations
	Bdash <- t(X)%*%Y
	Adash <- t(X)%*%X
	temp <- solve(A+Adash)%*%(B+Bdash)
	A <- Adash
	B <- Bdash

	# Add the new points obtained online, to the plot
	Sys.sleep(0.5)
	plot(X[,2],Y)
	abline(a=temp[1],b=temp[2])
}