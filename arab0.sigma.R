
nf = estimate.norm.factors(arab0)
nb.da <- prepare.nb.data(arab0, norm.factors=nf)

x = model.matrix(~ as.factor(group0))
dispersion.0 = estimate.dispersion(nb.da, x=x, model="NBQ")

row.mean <- apply(nb.da$rel.frequencies, 1, mean)
y.0 <- dispersion.0$estimates[,1]
x.0 = log(row.mean)
fit <- glm( y.0~ poly(x.0, degree=2), family=Gamma(link="log"))
phi.hat <- fitted(fit)
phi <- expandAsMatrix(phi.hat, c(dim(arab0)[1], dim(arab0)[2]))
length(phi.hat)

# this is because dispersion is not estimated from NBQ

system.time(sigma.0 <- estimate.dispersion.var.edgeR(nb.da, phi, x = x))

sigma.0 # 2.93


plot(x.0, log(y.0), pch=20, main="arab0")
range(dispersion.0$estimates)
range(phi.hat)
points(x.0, log(phi.hat), pch=8, col="red")


plot(xx0, log(yy0), pch=20, main="dipsersion est. by edgeR")
points(xx0, log(ll0$phi[, 2]), pch=8, col="red")


plot(xx1, log(yy0), pch=20, main="dipsersion est. by edgeR")
points(xx1, log(ll0$phi[, 2]), pch=8, col="red")

