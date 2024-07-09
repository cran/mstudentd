# Dimension p = 1

Sigma1 <- 0.5
Sigma2 <- 1

kl1_12 <- kldstudent(nu1 = 1, Sigma1, nu2 = 1, Sigma2, eps = 1e-16)
kl1_21 <- kldstudent(nu1 = 1, Sigma2, nu2 = 1, Sigma1, eps = 1e-16)

lambda <- 0.5

test_that("kl works (dim 1)", {
  expect_equal(
    round(as.numeric(kl1_21), 15),
    log( (1 + sqrt(lambda))^2 / (4*sqrt(lambda)) )
  )
  expect_equal(
    round(as.numeric(kl1_21), 15),
    log( (1 + sqrt(lambda))^2 / (4*sqrt(lambda)) )
  )
})

# Dimension p = 2

Sigma1 <- diag(0.5, nrow = 2)
Sigma2 <- diag(1, nrow = 2)

kl2_12 <- kldstudent(nu1 = 1, Sigma1, nu2 = 1, Sigma2, eps = 1e-16)
kl2_21 <- kldstudent(nu1 = 1, Sigma2, nu2 = 1, Sigma1, eps = 1e-16)

lambda <- as.complex(0.5)

test_that("kl works (dim 2)", {
  expect_equal(
    round(as.numeric(kl2_12), 15),
    Re(-log(lambda) + 3/sqrt(1-1/lambda) * log(sqrt(lambda) + sqrt(lambda-1)) - 3)
  )
  expect_equal(
    round(as.numeric(kl2_21), 15),
    Re(log(lambda) + 3/sqrt(1-lambda) * log(sqrt(1/lambda) + sqrt(1/lambda-1)) - 3)
  )
})

# Dimension p = 2; 2nd example

Sigma1 <- matrix(c(0.5, 0, 0, 1), nrow = 2)
Sigma2 <- diag(nrow = 2)

lambda <- 0.5

kl2 <- kldstudent(nu1 = 1, Sigma1, nu2 = 1, Sigma2, eps = 1e-16)

test_that("kl works (dim 2, one of the eigenvalues = 1)", {
  expect_equal(
    round(as.numeric(kl2), 15),
    log(lambda) - 3/2 * 1/sqrt(1-lambda) * log((1 - sqrt(1-lambda))/(1 + sqrt(1-lambda))) - 3
  )
})

nu1 <- 2; nu2 <- 4

#Dimension p = 3

Sigma1 <- 4*rbind(c(1, 0.6, 0.2), c(0.6, 1, 0.3), c(0.2, 0.3, 1))
Sigma2 <- rbind(c(1, 0.3, 0.1), c(0.3, 1, 0.4), c(0.1, 0.4, 1))

lambda <- 0.5

kl3_12 <- kldstudent(nu1, Sigma1, nu2, Sigma2, eps = 1e-8)
kl3_21 <- kldstudent(nu1, Sigma2, nu2, Sigma1, eps = 5e-5)
test_that("kl works (dim 3)", {
  expect_equal(
    attr(kl3_12, "epsilon"), 1e-8
  )
  expect_equal(
    round(as.numeric(kl3_12), 16), 0.9297752865860369
  )

  expect_equal(
    attr(kl3_21, "epsilon"), 5e-5
  )
  expect_equal(
    round(as.numeric(kl3_21), 16), 0.4074954441658625
  )
})

# Dimension p = 3, 2nd example

Sigma1 <- 2*rbind(c(1, 0.6, 0.2), c(0.6, 1, 0.3), c(0.2, 0.3, 1))
Sigma2 <- rbind(c(1, 0.3, 0.1), c(0.3, 1, 0.4), c(0.1, 0.4, 1))

kl3 <- kldstudent(nu1, Sigma1, nu2, Sigma2, eps = 1e-16)

test_that("kl12 works (dim 3, 2nd)", {
  expect_equal(
    attr(kl3, "epsilon"), 1e-16
  )
  expect_equal(
    round(as.numeric(kl3), 16), 0.3979439491689158
  )
})

# Dimension p = 4

Sigma1 <- 4*rbind(c(1, 0.6, 0.2, 0),
                  c(0.6, 1, 0.3, 0),
                  c(0.2, 0.3, 1, 0),
                  c(0, 0, 0, 1))
Sigma2 <- rbind(c(1, 0.3, 0.1, 0),
                c(0.3, 1, 0.4, 0),
                c(0.1, 0.4, 1, 0),
                c(0, 0, 0, 1))

kl4_12 <- kldstudent(nu1, Sigma1, nu2, Sigma2, eps = 1e-6)
kl4_21 <- kldstudent(nu1, Sigma2, nu2, Sigma1, eps = 1e-6)

test_that("kl12 works (dim 4)", {
  expect_equal(
    attr(kl4_12, "epsilon"), 1e-06
  )
  expect_equal(
    round(as.numeric(kl4_12), 16), 1.039925196101446
  )

  expect_equal(
    attr(kl4_21, "epsilon"), 1e-06
  )
  expect_equal(
    round(as.numeric(kl4_21), 16), 0.5359743613606762
  )
})
