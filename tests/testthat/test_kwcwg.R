context("KW-CWG distribution")
require(elfDistr)

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


