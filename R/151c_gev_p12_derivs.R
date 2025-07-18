######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_p12_fd=function (x, t1, t2, v1, v2, v3, v4, v5) 
{
    .e2 <- exp(-(t2 * v4 + v3))
    .e5 <- x - (t1 * v2 + v1)
    .e6 <- 1/v5
    .e8 <- v5 * .e2 * .e5
    .e9 <- 1 + .e8
    .e10 <- 1 + .e6
    .e12 <- .e9^.e10
    .e13 <- .e9^(.e6 + 2)
    .e14 <- exp(-.e9^-.e6)
    .e19 <- v5 * .e10/.e13 - 1/.e9^(2 * .e10)
    .e23 <- .e2 * .e19 * .e5 - 1/.e12
    .e24 <- .e2^2
    .e25 <- log1p(.e8)
    c(v1 = .e14 * .e24 * .e19, v2 = t1 * .e14 * .e24 * .e19, 
        v3 = .e14 * .e2 * .e23, v4 = t2 * .e14 * .e2 * .e23, 
        v5 = ((.e25/(v5 * .e12) - (.e25/(v5 * .e9^.e6) - .e2 * 
            .e5/.e12)/.e12)/v5 - .e10 * .e2 * .e5/.e13) * .e14 * 
            .e2)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p12_fdd=function (x, t1, t2, v1, v2, v3, v4, v5) 
{
    .e2 <- exp(-(t2 * v4 + v3))
    .e5 <- x - (t1 * v2 + v1)
    .e6 <- 1/v5
    .e8 <- v5 * .e2 * .e5
    .e9 <- 1 + .e8
    .e10 <- 1 + .e6
    .e11 <- .e6 + 2
    .e12 <- .e9^.e10
    .e13 <- 2 * .e10
    .e14 <- .e9^.e11
    .e15 <- log1p(.e8)
    .e16 <- v5 * .e10
    .e17 <- .e9^.e13
    .e18 <- .e9^.e6
    .e19 <- 1/.e17
    .e20 <- .e2 * .e5
    .e21 <- v5 * .e18
    .e23 <- 2 * .e11
    .e24 <- exp(-.e9^-.e6)
    .e25 <- .e20/.e12
    .e28 <- .e16/.e14 - .e19
    .e30 <- v5 * .e9^(.e10 - .e23) * .e11
    .e31 <- 1/.e12
    .e32 <- 1/.e14
    .e34 <- .e15/.e21 - .e25
    .e36 <- .e9^(.e6 - .e13)
    .e38 <- .e30 - 2/.e9^(1 + .e13)
    .e39 <- v5^2
    .e40 <- v5 * .e12
    .e41 <- .e2^2
    .e44 <- .e2 * .e28 * .e5 - .e31
    .e46 <- .e34/.e12
    .e47 <- .e19 + .e16 * (.e2 * .e38 * .e5 - (.e36 + .e32))
    .e48 <- .e15/.e40
    .e49 <- (.e48 - .e46)/v5
    .e53 <- .e10 * .e18 * .e2 * .e5 - .e12 * .e15/.e39
    .e56 <- .e10 * .e2 * .e5/.e14
    .e60 <- .e12 * .e11 * .e2 * .e5 - .e14 * .e15/.e39
    .e61 <- .e9^(.e6 - 1)
    .e62 <- .e9^.e23
    .e63 <- .e40^2
    .e64 <- .e21^2
    .e65 <- .e53/.e17
    .e66 <- (.e49 - .e56)/.e12
    .e85 <- .e47 * .e2 * .e5 - (1 + .e25) * .e44
    .e87 <- (2 * (.e10 * .e9^(.e13 - 1) * .e2 * .e5) - 2 * (.e17 * 
        .e15/.e39))/.e9^(4 * .e10) + .e32
    .e90 <- (.e16 * (.e21 * .e15/.e63 - .e36 * .e34) - ((.e31 + 
        v5 * (.e61 * .e15/.e64 - .e10 * .e36 * .e2 * .e5) - .e31)/.e12 + 
        .e32))/v5
    .e92 <- .e47 - .e44/.e12
    .e93 <- .e2^3
    .e96 <- v5 * .e60 * .e10/.e62
    .e98 <- .e16 * .e38 - .e28/.e12
    .e102 <- .e16 * .e2 * .e38 * .e5 - (2 + .e25) * .e28
    .e104 <- .e30 * .e2 * .e5
    .e105 <- t1 * .e24
    .e107 <- .e65 + (.e87 - .e96) * .e2 * .e5 - .e44 * .e34/v5
    .e110 <- (.e10 * (2/.e14 - .e104) + .e90 - .e66) * .e2 * 
        .e5 - .e49
    .e112 <- .e87 - (.e34 * .e28/v5 + .e96)
    .e113 <- .e90 - (.e66 + .e10 * (.e104 - .e32))
    .e115 <- .e105 * .e93 * .e98
    .e116 <- t1 * t2
    .e119 <- t2 * .e85 * .e24 * .e2
    c(v1 = c(v1 = .e24 * .e93 * .e98, v2 = .e115, v3 = .e92 * 
        .e24 * .e41, v4 = t2 * .e92 * .e24 * .e41, v5 = .e113 * 
        .e24 * .e41), v2 = c(v1 = .e115, v2 = t1^2 * .e24 * .e93 * 
        .e98, v3 = t1 * .e92 * .e24 * .e41, v4 = .e116 * .e92 * 
        .e24 * .e41, v5 = t1 * .e113 * .e24 * .e41), v3 = c(v1 = .e24 * 
        .e41 * .e102, v2 = .e105 * .e41 * .e102, v3 = .e85 * 
        .e24 * .e2, v4 = .e119, v5 = .e110 * .e24 * .e2), v4 = c(v1 = t2 * 
        .e24 * .e41 * .e102, v2 = .e116 * .e24 * .e41 * .e102, 
        v3 = .e119, v4 = t2^2 * .e85 * .e24 * .e2, v5 = t2 * 
            .e110 * .e24 * .e2), v5 = c(v1 = .e112 * .e24 * .e41, 
        v2 = t1 * .e112 * .e24 * .e41, v3 = .e107 * .e24 * .e2, 
        v4 = t2 * .e107 * .e24 * .e2, v5 = (((.e65 + .e56 - .e49) * 
            .e34 + (.e46 + .e20/.e14 - .e48)/v5 - (((.e65 + 1/.e40) * 
            .e2 * .e5 - (.e61 * .e2 * .e5 + .e18 - .e18 * .e15/v5) * 
            .e15/.e64)/.e12 + (.e12 + v5 * .e53) * .e15/.e63))/v5 + 
            (.e60 * .e10/.e62 + 1/(.e39 * .e14)) * .e2 * .e5) * 
            .e24 * .e2))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_p12_pd=function (x, t1, t2, v1, v2, v3, v4, v5) 
{
    .e2 <- exp(-(t2 * v4 + v3))
    .e5 <- x - (t1 * v2 + v1)
    .e7 <- v5 * .e2 * .e5
    .e8 <- 1 + .e7
    .e9 <- 1/v5
    .e11 <- .e8^(1 + .e9)
    .e12 <- exp(-.e8^-.e9)
    .e13 <- .e12 * .e2
    c(v1 = -(.e13/.e11), v2 = -(t1 * .e12 * .e2/.e11), v3 = -(.e13 * 
        .e5/.e11), v4 = -(t2 * .e12 * .e2 * .e5/.e11), v5 = -(.e12 * 
        (log1p(.e7)/(v5 * .e8^.e9) - .e2 * .e5/.e11)/v5))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p12_pdd=function (x, t1, t2, v1, v2, v3, v4, v5) 
{
    .e2 <- exp(-(t2 * v4 + v3))
    .e5 <- x - (t1 * v2 + v1)
    .e6 <- 1/v5
    .e8 <- v5 * .e2 * .e5
    .e9 <- 1 + .e8
    .e10 <- 1 + .e6
    .e11 <- .e9^.e10
    .e12 <- 2 * .e10
    .e14 <- exp(-.e9^-.e6)
    .e15 <- .e9^.e6
    .e17 <- .e2 * .e5/.e11
    .e18 <- log1p(.e8)
    .e19 <- .e9^(.e6 - .e12)
    .e21 <- v5 * .e10 * .e19
    .e22 <- v5 * .e15
    .e26 <- .e21 * .e2 * .e5 - (1 + .e17)/.e11
    .e27 <- .e18/.e22
    .e28 <- .e9^.e12
    .e29 <- .e27 - .e17
    .e30 <- 1/.e11
    .e31 <- (.e10 * .e15 * .e2 * .e5 - .e11 * .e18/v5^2)/.e28
    .e32 <- .e9^(.e6 - 1)
    .e33 <- .e22^2
    .e34 <- v5 * .e11
    .e35 <- .e31 + .e29/.e34
    .e45 <- .e30 + v5 * (.e32 * .e18/.e33 - .e10 * .e19 * .e2 * 
        .e5) - (.e29/.e11 + .e30)
    .e47 <- .e2^2
    .e48 <- t1 * .e14
    .e51 <- t2 * .e14 * .e2 * .e26
    .e52 <- .e21 - 1/.e28
    .e54 <- .e14 * .e2 * .e26
    .e55 <- -.e54
    .e56 <- -(.e48 * .e2 * .e26)
    .e57 <- -(.e48 * .e47 * .e52)
    .e58 <- -(t1 * t2 * .e14 * .e2 * .e26)
    .e59 <- -(.e51 * .e5)
    .e60 <- -.e51
    .e62 <- .e35 * .e14 * .e2
    .e64 <- .e45 * .e14 * .e2
    c(v1 = c(v1 = -(.e14 * .e47 * .e52), v2 = .e57, v3 = .e55, 
        v4 = .e60, v5 = -(.e64/v5)), v2 = c(v1 = .e57, v2 = -(t1^2 * 
        .e14 * .e47 * .e52), v3 = .e56, v4 = .e58, v5 = -(t1 * 
        .e45 * .e14 * .e2/v5)), v3 = c(v1 = .e55, v2 = .e56, 
        v3 = -(.e54 * .e5), v4 = .e59, v5 = -(.e64 * .e5/v5)), 
        v4 = c(v1 = .e60, v2 = .e58, v3 = .e59, v4 = -(t2^2 * 
            .e14 * .e2 * .e26 * .e5), v5 = -(t2 * .e45 * .e14 * 
            .e2 * .e5/v5)), v5 = c(v1 = .e62, v2 = t1 * .e35 * 
            .e14 * .e2, v3 = .e62 * .e5, v4 = t2 * .e35 * .e14 * 
            .e2 * .e5, v5 = -(((.e31 + 1/.e34) * .e2 * .e5 - 
            ((.e32 * .e2 * .e5 + .e15 - .e15 * .e18/v5) * .e18/.e33 + 
                (1 + .e27 - .e17) * .e29/v5)) * .e14/v5)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p12_logfdd=function (x, t1, t2, v1, v2, v3, v4, v5) 
{
    .e2 <- exp(-(t2 * v4 + v3))
    .e5 <- x - (t1 * v2 + v1)
    .e7 <- v5 * .e2 * .e5
    .e8 <- 1 + .e7
    .e9 <- 1/v5
    .e10 <- .e8^.e9
    .e11 <- 1 + .e9
    .e12 <- 1/.e10
    .e13 <- v5 * .e11
    .e14 <- .e13 - .e12
    .e15 <- .e8^(.e9 - 1)
    .e17 <- v5 * .e14 - .e12
    .e19 <- (1/.e15 + .e2 * .e17 * .e5)/.e8 - .e13
    .e20 <- log1p(.e7)
    .e23 <- .e15 * .e2 * .e5 - .e10 * .e20/v5
    .e25 <- .e12 - 1
    .e26 <- 2 * .e11
    .e27 <- v5 * .e8^(2/v5 - 1)
    .e29 <- (.e23/.e27 - .e2 * .e14 * .e5)/.e8 + 1
    .e31 <- (((.e20/.e10 - v5 * .e25)/v5 + .e12)/.e8 - .e13 * 
        .e8^(.e9 - .e26) * .e2 * .e5)/v5 + .e11 * (.e7/.e8 - 
        1)/.e8
    .e32 <- .e8^2
    .e33 <- .e2^2
    .e35 <- t2 * .e19 * .e2
    .e36 <- .e19 * .e2
    .e37 <- .e29 * .e2
    .e38 <- .e31 * .e2
    .e39 <- .e36/.e8
    .e42 <- t1 * .e19 * .e2/.e8
    .e45 <- t1 * .e33 * .e17/.e32
    .e49 <- t1 * t2 * .e19 * .e2/.e8
    .e51 <- .e35 * .e5/.e8
    .e52 <- .e35/.e8
    .e53 <- v5^2
    c(v1 = c(v1 = .e33 * .e17/.e32, v2 = .e45, v3 = .e39, v4 = .e52, 
        v5 = -.e38), v2 = c(v1 = .e45, v2 = t1^2 * .e33 * .e17/.e32, 
        v3 = .e42, v4 = .e49, v5 = -(t1 * .e31 * .e2)), v3 = c(v1 = .e39, 
        v2 = .e42, v3 = .e36 * .e5/.e8, v4 = .e51, v5 = -(.e38 * 
            .e5)), v4 = c(v1 = .e52, v2 = .e49, v3 = .e51, v4 = t2^2 * 
        .e19 * .e2 * .e5/.e8, v5 = -(t2 * .e31 * .e2 * .e5)), 
        v5 = c(v1 = .e37/.e8, v2 = t1 * .e29 * .e2/.e8, v3 = .e37 * 
            .e5/.e8, v4 = t2 * .e29 * .e2 * .e5/.e8, v5 = -(((((2/.e10 - 
            1) * .e2 * .e5 - .e23 * .e20/.e27)/.e8 - 2 * (.e25 * 
            .e20/v5))/v5 + (.e11 * .e10 * .e2 * .e5 - .e8^.e11 * 
            .e20/.e53) * .e2 * .e5/.e8^.e26)/v5 - (.e11 * .e2 * 
            .e5/.e8 + 1/.e53) * .e2 * .e5/.e8)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
gev_p12_logfddd=function (x, t1, t2, v1, v2, v3, v4, v5) 
{
    .e2 <- exp(-(t2 * v4 + v3))
    .e5 <- x - (t1 * v2 + v1)
    .e7 <- v5 * .e2 * .e5
    .e8 <- 1 + .e7
    .e9 <- 1/v5
    .e10 <- .e8^.e9
    .e11 <- 1 + .e9
    .e12 <- .e9 - 1
    .e13 <- 1/.e10
    .e14 <- .e8^.e12
    .e15 <- log1p(.e7)
    .e16 <- v5 * .e11
    .e17 <- 2/v5
    .e18 <- .e16 - .e13
    .e19 <- .e8^.e11
    .e20 <- .e17 - 1
    .e22 <- .e14 * .e2 * .e5
    .e23 <- 1/.e14
    .e25 <- .e10 * .e15/v5
    .e26 <- .e22 - .e25
    .e28 <- v5 * .e18 - .e13
    .e29 <- .e8^.e20
    .e31 <- .e2 * .e28 * .e5
    .e32 <- v5 * .e29
    .e33 <- 2 * .e11
    .e34 <- 2/.e14
    .e35 <- v5/.e19
    .e36 <- .e26/.e32
    .e37 <- .e9 - .e33
    .e38 <- 2 * .e12
    .e39 <- .e8^2
    .e40 <- v5^2
    .e42 <- .e8^(.e9 - (2 + .e38)) * .e12
    .e44 <- (1/.e19 + .e35) * .e2 * .e5
    .e45 <- 2 * .e16
    .e46 <- .e13 - 1
    .e47 <- .e8^(.e17 - 2)
    .e48 <- (.e15/.e10 - v5 * .e46)/v5
    .e49 <- .e8^.e37
    .e50 <- .e48 + .e13
    .e51 <- .e7/.e8
    .e52 <- 2/.e10
    .e53 <- .e14 * .e15
    .e54 <- .e8^(.e9 - 2)
    .e56 <- .e2 * .e18 * .e5
    .e57 <- .e2^2
    .e58 <- .e36 - .e56
    .e65 <- ((.e52 + v5 * (.e42 + (.e23 + .e34 + .e31)/.e8 - 
        .e45) - .e44) * .e2 * .e5 - .e23)/.e8 + v5 * (((.e23 + 
        .e31)/.e8 - .e16) * .e2 * .e5/.e8 + 1 + .e9)
    .e66 <- .e23 + 2 * .e31
    .e67 <- .e51 - 1
    .e68 <- .e8^.e17
    .e69 <- .e66 + .e34
    .e70 <- .e32^2
    .e71 <- .e2 * .e5
    .e72 <- 2 * .e28
    .e73 <- v5 * .e47
    .e75 <- (.e53 + v5 * .e14)/v5 - (.e14 + v5 * .e54 * .e12 * 
        .e2 * .e5)
    .e76 <- .e8^(.e9 - (1 + .e33))
    .e77 <- .e71/.e10
    .e78 <- v5 * .e26
    .e79 <- t1 * t2
    .e84 <- .e54 * .e12 * .e2 * .e5 - .e53/.e40
    .e85 <- .e58/.e8
    .e88 <- .e75/.e73 + .e77 + v5 * .e58
    .e90 <- (((.e13 + v5 * .e50)/.e8 + (.e15/.e19 - 2 * .e35)/v5)/.e8 + 
        v5 * (.e49 + v5 * .e76 * .e37 * .e2 * .e5) * .e11)/v5 + 
        2 * (.e16 * .e67/.e39)
    .e91 <- (.e50/.e8 - .e16 * .e49 * .e2 * .e5)/v5
    .e92 <- .e26/.e68
    .e96 <- .e13 + v5 * (.e42 + .e69/.e8 - .e45) - .e44
    .e100 <- .e78 * .e47 * .e20/.e70
    .e102 <- v5 * (.e72 - .e13) - .e13
    .e104 <- (((.e88 - .e34)/.e8 + v5 * (.e85 + .e17 + 3 + .e100)) * 
        .e2 * .e5 - .e36)/.e8 - 1
    .e106 <- (((.e36 - .e69)/.e8 + v5 * ((.e92 + 2)/v5 + 3)) * 
        .e2 * .e5 - .e84/.e8^.e38)/.e8 - 1
    .e109 <- .e90 * .e2 * .e5 - (.e91 + .e11 * .e67/.e8)
    .e110 <- t1^2
    .e111 <- t2^2
    .e112 <- .e11 * .e10
    .e113 <- .e19 * .e15
    .e114 <- .e113/.e40
    .e116 <- t2 * .e65 * .e2
    .e117 <- .e8^3
    .e122 <- .e2 * .e102 * .e5/.e8 - .e72
    .e123 <- .e2^3
    .e125 <- .e26 * .e15/.e32
    .e127 <- .e112 * .e2 * .e5
    .e128 <- .e127 - .e114
    .e129 <- .e8^.e33
    .e131 <- (.e52 - 1) * .e2 * .e5
    .e134 <- .e79 * .e65 * .e2/.e8
    .e135 <- .e116/.e8
    .e137 <- .e111 * .e65 * .e2
    .e141 <- .e84 * .e2 * .e5 - ((.e22 - (.e10 + .e25)) * .e15/v5 + 
        .e22)/v5
    .e145 <- .e11 * .e2 * .e5/.e8 + 1/.e40
    .e146 <- .e8^(.e33 - 1)
    .e149 <- .e29 + v5 * (.e47 * .e20 * .e2 * .e5 - 2 * (.e29 * 
        .e15/.e40))
    .e150 <- .e131 - .e125
    .e151 <- v5 * .e68
    .e153 <- ((((.e112 * .e15 + .e10)/v5 - (.e22 + .e10 + .e10) * 
        .e11) * .e2 * .e5 + .e114)/.e146 + ((2 * .e77 + v5 * 
        .e150 - ((.e75 * .e15 - .e78/.e8)/.e73 + .e34))/.e8 + 
        2 - (2 * .e48 + .e40 * .e26 * .e47 * .e20 * .e15/.e70))/v5 + 
        2 * (v5 * .e128 * .e11 * .e2 * .e5/.e129))/v5 - (.e11 * 
        (.e51 - 2) + v5 * .e145) * .e2 * .e5/.e8
    .e155 <- (((.e26 * (1/.e68 - (1/.e68 + .e15/.e151)) + 1 + 
        .e71/.e19 - .e50)/v5 - .e50 * .e2 * .e5/.e8)/.e8 - (.e91 + 
        (.e49 + v5 * (.e76 * .e37 * .e2 * .e5 + .e49 * .e15/.e40) * 
            .e11) * .e2 * .e5))/v5 + (.e11 * (2 - 2 * .e51) * 
        .e2 * .e5/.e8 - .e67/.e40)/.e8
    .e160 <- (.e88 - .e23)/.e8 + v5 * (.e85 + .e9 + 2 + .e100)
    .e164 <- (.e36 - .e66)/.e8 + v5 * ((.e92 + 1)/v5 + 2)
    .e166 <- .e141/.e32 - (((.e26 * (1/.e29 + 2/.e29)/v5 - 2 * 
        .e56)/.e8 + 2) * .e2 * .e5 + .e26 * .e149/.e70)
    .e167 <- .e65 * .e2
    .e170 <- t1 * .e96 * .e57/.e39
    .e173 <- .e79 * .e96 * .e57/.e39
    .e175 <- t2 * .e104 * .e2
    .e177 <- t2 * .e106 * .e2
    .e179 <- t2 * .e109 * .e2
    .e180 <- .e104 * .e2
    .e181 <- .e106 * .e2
    .e182 <- .e109 * .e2
    .e183 <- .e167/.e8
    .e186 <- t1 * .e65 * .e2/.e8
    .e189 <- t1 * .e123 * .e102/.e117
    .e193 <- t1 * .e111 * .e65 * .e2/.e8
    .e196 <- .e110 * .e123 * .e102/.e117
    .e197 <- .e110 * t2
    .e199 <- .e116 * .e5/.e8
    .e201 <- .e137 * .e5/.e8
    .e202 <- .e137/.e8
    .e203 <- -.e182
    .e204 <- -(t1 * .e109 * .e2)
    .e206 <- -(.e79 * .e109 * .e2)
    .e208 <- -.e179
    .e209 <- .e180/.e8
    .e210 <- .e153 * .e2
    .e211 <- .e155 * .e2
    .e212 <- .e181/.e8
    .e213 <- .e166 * .e2
    .e215 <- .e46 * .e15/v5
    .e217 <- .e96 * .e57/.e39
    .e218 <- c(v1 = .e189, v2 = .e196, v3 = .e170, v4 = .e173, 
        v5 = -(t1 * .e90 * .e57))
    .e219 <- c(v1 = .e135, v2 = .e134, v3 = .e199, v4 = .e201, 
        v5 = -(.e179 * .e5))
    .e222 <- t1 * .e104 * .e2/.e8
    .e225 <- t1 * .e106 * .e2/.e8
    .e228 <- t1 * .e160 * .e57/.e39
    .e231 <- t1 * .e164 * .e57/.e39
    .e234 <- t1 * .e57 * .e122/.e39
    .e237 <- .e79 * .e104 * .e2/.e8
    .e240 <- .e79 * .e106 * .e2/.e8
    .e243 <- .e79 * .e57 * .e122/.e39
    .e246 <- .e110 * .e96 * .e57/.e39
    .e249 <- .e197 * .e96 * .e57/.e39
    .e251 <- .e175 * .e5/.e8
    .e252 <- .e175/.e8
    .e254 <- .e177 * .e5/.e8
    .e255 <- .e177/.e8
    .e258 <- t2 * .e96 * .e57/.e39
    c(v1 = c(v1 = c(v1 = .e123 * .e102/.e117, v2 = .e189, v3 = .e217, 
        v4 = .e258, v5 = -(.e90 * .e57)), v2 = .e218, v3 = c(v1 = .e217, 
        v2 = .e170, v3 = .e183, v4 = .e135, v5 = .e203), v4 = c(v1 = .e258, 
        v2 = .e173, v3 = .e135, v4 = .e202, v5 = .e208), v5 = c(v1 = .e160 * 
        .e57/.e39, v2 = .e228, v3 = .e209, v4 = .e252, v5 = -(.e210/.e8))), 
        v2 = c(v1 = .e218, v2 = c(v1 = .e196, v2 = t1^3 * .e123 * 
            .e102/.e117, v3 = .e246, v4 = .e249, v5 = -(.e110 * 
            .e90 * .e57)), v3 = c(v1 = .e170, v2 = .e246, v3 = .e186, 
            v4 = .e134, v5 = .e204), v4 = c(v1 = .e173, v2 = .e249, 
            v3 = .e134, v4 = .e193, v5 = .e206), v5 = c(v1 = .e228, 
            v2 = .e110 * .e160 * .e57/.e39, v3 = .e222, v4 = .e237, 
            v5 = -(t1 * .e153 * .e2/.e8))), v3 = c(v1 = c(v1 = .e57 * 
            .e122/.e39, v2 = .e234, v3 = .e183, v4 = .e135, v5 = .e203), 
            v2 = c(v1 = .e234, v2 = .e110 * .e57 * .e122/.e39, 
                v3 = .e186, v4 = .e134, v5 = .e204), v3 = c(v1 = .e183, 
                v2 = .e186, v3 = .e167 * .e5/.e8, v4 = .e199, 
                v5 = -(.e182 * .e5)), v4 = .e219, v5 = c(v1 = .e209, 
                v2 = .e222, v3 = .e180 * .e5/.e8, v4 = .e251, 
                v5 = -(.e210 * .e5/.e8))), v4 = c(v1 = c(v1 = t2 * 
            .e57 * .e122/.e39, v2 = .e243, v3 = .e135, v4 = .e202, 
            v5 = .e208), v2 = c(v1 = .e243, v2 = .e197 * .e57 * 
            .e122/.e39, v3 = .e134, v4 = .e193, v5 = .e206), 
            v3 = .e219, v4 = c(v1 = .e202, v2 = .e193, v3 = .e201, 
                v4 = t2^3 * .e65 * .e2 * .e5/.e8, v5 = -(.e111 * 
                  .e109 * .e2 * .e5)), v5 = c(v1 = .e252, v2 = .e237, 
                v3 = .e251, v4 = .e111 * .e104 * .e2 * .e5/.e8, 
                v5 = -(t2 * .e153 * .e2 * .e5/.e8))), v5 = c(v1 = c(v1 = .e164 * 
            .e57/.e39, v2 = .e231, v3 = .e212, v4 = .e255, v5 = -.e211), 
            v2 = c(v1 = .e231, v2 = .e110 * .e164 * .e57/.e39, 
                v3 = .e225, v4 = .e240, v5 = -(t1 * .e155 * .e2)), 
            v3 = c(v1 = .e212, v2 = .e225, v3 = .e181 * .e5/.e8, 
                v4 = .e254, v5 = -(.e211 * .e5)), v4 = c(v1 = .e255, 
                v2 = .e240, v3 = .e254, v4 = .e111 * .e106 * 
                  .e2 * .e5/.e8, v5 = -(t2 * .e155 * .e2 * .e5)), 
            v5 = c(v1 = .e213/.e39, v2 = t1 * .e166 * .e2/.e39, 
                v3 = .e213 * .e5/.e39, v4 = t2 * .e166 * .e2 * 
                  .e5/.e39, v5 = -(((((.e26 * .e11 - .e10/v5) * 
                  .e2 * .e5 - ((.e127 - (.e113/v5 + 2 * .e19)/v5) * 
                  .e15 + .e10 * .e2 * .e5)/v5)/(v5 * .e129) - 
                  .e128 * (2 * (.e11 * .e146 * .e2 * .e5) - 2 * 
                    (.e129 * .e15/.e40))/.e8^(4 * .e11)) * .e2 * 
                  .e5 + (((.e149 * .e15/.e70 - 2 * (.e71/.e151)) * 
                  .e26 - ((.e141 * .e15 + .e26 * .e2 * .e5/.e8)/.e73 + 
                  (.e128/.e8^(.e33 - 2) + .e131 - .e125) * .e2 * 
                    .e5)/.e8)/.e8 - (2 * ((.e46 * .e2 * .e5 - 
                  .e125)/.e8 - .e215) + 2 * (.e150/.e8) - 4 * 
                  .e215)/v5)/v5)/v5 + (2 * (.e145 * .e2 * .e5/.e8) + 
                  2/v5^3) * .e2 * .e5/.e8))))
}
############################################################
#' The first derivative of the density for DMGS
#' @returns Vector
#' @inheritParams manf
gev_p12_f1fa=function(x,t01,t02,v1,v2,v3,v4,v5){

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p12_fd,"x")
	f1=vf(x,t01,t02,v1,v2,v3,v4,v5)
	return(f1)
}
############################################################
#' The first derivative of the density for WAIC
#' @returns Vector
#' @inheritParams manf
gev_p12_f1fw=function(x,t1,t2,v1,v2,v3,v4,v5){

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p12_fd,c("x","t1","t2"))
	f1=vf(x,t1,t2,v1,v2,v3,v4,v5)
	return(f1)
}
############################################################
#' The second derivative of the density for DMGS
#' @returns Matrix
#' @inheritParams manf
gev_p12_f2fa=function(x,t01,t02,v1,v2,v3,v4,v5){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p12_fdd,"x")
	temp1=vf(x,t01,t02,v1,v2,v3,v4,v5)
	f2=deriv_copyfdd(temp1,nx,dim=5)
	return(f2)
}
############################################################
#' The second derivative of the density for WAIC
#' @returns Matrix
#' @inheritParams manf
gev_p12_f2fw=function(x,t1,t2,v1,v2,v3,v4,v5){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p12_fdd,c("x","t1","t2"))
	temp1=vf(x,t1,t2,v1,v2,v3,v4,v5)
	f2=deriv_copyfdd(temp1,nx,dim=5)
	return(f2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
gev_p12_mu1fa=function(alpha,t01,t02,v1,v2,v3,v4,v5){
	x=qgev((1-alpha),mu=v1+v2*t01,sigma=exp(v3+v4*t02),xi=v5)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p12_pd,"x")
	mu1=-vf(x,t01,t02,v1,v2,v3,v4,v5)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
gev_p12_mu2fa=function(alpha,t01,t02,v1,v2,v3,v4,v5){
	x=qgev((1-alpha),mu=v1+v2*t01,sigma=exp(v3+v4*t02),xi=v5)
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p12_pdd,"x")
	temp1=vf(x,t01,t02,v1,v2,v3,v4,v5)
	mu2=-deriv_copyfdd(temp1,nx,dim=5)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
gev_p12_ldda=function(x,t1,t2,v1,v2,v3,v4,v5){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p12_logfdd,c("x","t1","t2"))
	temp1=vf(x,t1,t2,v1,v2,v3,v4,v5)
	ldd=deriv_copyldd(temp1,nx,dim=5)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
gev_p12_lddda=function(x,t1,t2,v1,v2,v3,v4,v5){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p12_logfddd,c("x","t1","t2"))
	temp1=vf(x,t1,t2,v1,v2,v3,v4,v5)
	lddd=deriv_copylddd(temp1,nx,dim=5)
	return(lddd)
}
