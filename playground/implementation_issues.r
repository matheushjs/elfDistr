require(Rcpp)

dggamma1 = function(x, a, b, k, log=F){
	result = log(b) - lgamma(k) + (b*k - 1)*log(x) - (b*k)*log(a) - (x/a)**b;
	if(!log) result = exp(result);
	return(result);
}

dggamma2 = function(x, a, b, k, log=F){
	maxLength = max(length(x), length(a), length(b), length(k));
	x = rep_len(x, maxLength)
	a = rep_len(a, maxLength)
	b = rep_len(b, maxLength)
	k = rep_len(k, maxLength)

	result = log(b) - lgamma(k) + (b*k - 1)*log(x) - (b*k)*log(a) - (x/a)**b;
	if(!log) result = exp(result);
	return(result);
}

dggamma3 = function(x, a, b, k, log=F){
	maxLength = max(length(x), length(a), length(b), length(k));
	if(length(x) != 1 && length(x) != maxLength) x = rep_len(x, maxLength)
	if(length(a) != 1 && length(a) != maxLength) a = rep_len(a, maxLength)
	if(length(b) != 1 && length(b) != maxLength) b = rep_len(b, maxLength)
	if(length(k) != 1 && length(k) != maxLength) k = rep_len(k, maxLength)

	result = log(b) - lgamma(k) + (b*k - 1)*log(x) - (b*k)*log(a) - (x/a)**b;
	if(!log) result = exp(result);
	return(result);
}

cppFunction("
	NumericVector dggamma4(
		NumericVector &vx,
		NumericVector &va,
		NumericVector &vb,
		NumericVector &vk,
		const bool &log_prob = false
	){
		int maxLength = std::max({
			vx.length(),
			va.length(),
			vb.length(),
			vk.length()
		});
		NumericVector p(maxLength);
		for(int i = 0; i < maxLength; i++){
			const double x = vx[i % vx.length()];
			const double a = va[i % va.length()];
			const double b = vb[i % vb.length()];
			const double k = vk[i % vk.length()];
			p[i] = std::log(b) - R::lgammafn(k) + (b*k - 1)*std::log(x) - (b*k)*std::log(a) - std::pow(x/a, b);
		}
		if(!log_prob)
			p = Rcpp::exp(p);
		return p;
	}
");

cppFunction("
	NumericVector dggamma5(
		NumericVector &vx,
		NumericVector &va,
		NumericVector &vb,
		NumericVector &vk,
		const bool &log_prob = false
	){
		int maxLength = std::max({
			vx.length(),
			va.length(),
			vb.length(),
			vk.length()
		});
		NumericVector p(maxLength);
		for(int i = 0; i < maxLength; i++){
			const double b = vb[i % vb.length()];
			p[i] = std::log(b);
		}
		for(int i = 0; i < maxLength; i++){
			const double x = vx[i % vx.length()];
			const double a = va[i % va.length()];
			const double b = vb[i % vb.length()];
			p[i] = p[i] - std::pow(x/a, b);
		}
		for(int i = 0; i < maxLength; i++){
			const double x = vx[i % vx.length()];
			const double a = va[i % va.length()];
			const double b = vb[i % vb.length()];
			const double k = vk[i % vk.length()];
			p[i] = p[i] - R::lgammafn(k) + (b*k - 1)*std::log(x) - (b*k)*std::log(a);
		}
		if(!log_prob)
			p = Rcpp::exp(p);
		return p;
	}
");

x = seq(0, 10, length=10000);
result = microbenchmark(
	dggamma1(x, 1, 1, 1),
	dggamma2(x, 1, 1, 1),
	dggamma3(x, 1, 1, 1),
	dggamma4(x, 1, 1, 1),
	dggamma5(x, 1, 1, 1),
	unit="ms",
	times=1000,
	control=list(warmup=100)
)
print(summary(result))
