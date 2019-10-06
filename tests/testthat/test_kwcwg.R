context("Kw-CWG distribution")

require(microbenchmark)
require(elfDistr)

kwcwg.pdf = function(x, alpha, beta, gamma, a, b){
	#cat("Called with", "(", x, ")", alpha, beta, gamma, a, b, "\n")
	
	# If parameters are not within valid range, return 0
	if(  sum(x < 0) > 0
	  || alpha < 0 || alpha > 1
	  || beta < 0
	  || gamma < 0
	  || a < 0
	  || b < 0
	){
		return(rep(0, length(x)));
	}
	
	# Original function
	# return(
	#	alpha**a * beta * gamma * a * b * (gamma * x)**(beta-1) * exp(-(gamma*x)**beta) *
	#	(
	#		(1 - exp(-(gamma*x)**beta))**(a-1) /
	#		(alpha + (1 - alpha)*exp(-(gamma*x)**beta))**(a+1)
	#	) *
	#	(
	#		1 -
	#		(alpha**a*(1 - exp(-(gamma*x)**beta))**a) /
	#		(alpha + (1-alpha)*exp(-(gamma*x)**beta))**a
	#	)**(b-1)
	#)
	
	# A common term in the equation
	aux1 = exp(-(gamma*x)**beta)
	
	# Here we will factor f(x) as being A * (B/C) * (1 - D/E)^(b-1)
	A = alpha**a * beta * gamma * a * b * (gamma * x)**(beta-1) * aux1
	B = 1 - aux1
	C = alpha + (1 - alpha)*aux1
	D = (alpha**a * (1 - aux1)**a)
	E = (alpha + (1-alpha)*aux1)**a
	
	# cat("A", A, "\n")
	# cat("B", B, "\n")
	# cat("C", C, "\n")
	# cat("D", D, "\n")
	# cat("E", E, "\n")
	
	result = A * (B**(a-1)/C**(a+1)) * (1 - D/E)**(b-1)
	
	return(result)
}

func_generator = function(alpha, beta, gamma, a, b){
	f = function(x){ dkwcwg(x, alpha, beta, gamma, a, b) }
	return(f)
}

test_that("Integrates to 1", {
	f = func_generator(0.1, 1, 1, 1, 1)
	expect_equal(integrate(f, 1e-10, 30)$value, 1, tolerance=0.001)

	f = func_generator(0.8, 0.2, 1, 1, 0.4)
	expect_equal(integrate(f, 1e-6, Inf)$value, 1, tolerance=0.03) # Could not get better than this, numerically

	f = func_generator(0.5, 0.5, 0.2, 0.1, 1)
	expect_equal(integrate(f, 1e-10, 150)$value, 1, tolerance=0.001)

	f = func_generator(0.5, 5, 2, 3, 1)
	expect_equal(integrate(f, 0, 1)$value, 1, tolerance=0.001)

	f = func_generator(0.5, 1, 2, 3, 10)
	expect_equal(integrate(f, 0, 2)$value, 1, tolerance=0.001)
})

test_that("Benchmark", {
	savedOptions = options()
	setOption("width", 150)

	x = seq(0.1, 2, length=150000)

	result = microbenchmark(
		kwcwg.pdf(1, 0.5, 1, 1, 1, 1),
		dkwcwg(1, 0.5, 1, 1, 1, 1),
		unit="ms",
		control=list(warmup=1000)
	)
	cat("\n\n============================\n")
	result = summary(result)
	print(result)
	expect_gt(result[1, "mean"], result[2, "mean"])

	result = microbenchmark(
		kwcwg.pdf(x, 0.5, 1, 1, 1, 1),
		dkwcwg(x, 0.5, 1, 1, 1, 1),
		unit="s",
		control=list(warmup=10)
	)
	cat("\n\n============================\n")
	result = summary(result)
	print(result)
	expect_gt(result[1, "mean"], result[2, "mean"])

	options(savedOptions)
})


