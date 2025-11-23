######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_p1c_fd=function (x, ta, tb, tc, v1, v2, v3, v4, v5, v6) 
{
    .e1 <- 1/v6
    .e6 <- x - (ta * v2 + tb * v3 + tc * v4 + v1)
    .e8 <- v6 * .e6/v5
    .e9 <- 1 + .e8
    .e10 <- 1 + .e1
    .e12 <- .e9^(.e1 + 2)
    .e13 <- exp(-.e9^-.e1)
    .e16 <- v5^2
    .e19 <- v6 * .e10/.e12 - 1/.e9^(2 * .e10)
    .e20 <- .e9^.e10
    .e21 <- log1p(.e8)
    c(v1 = .e13 * .e19/.e16, v2 = ta * .e13 * .e19/.e16, v3 = tb * 
        .e13 * .e19/.e16, v4 = tc * .e13 * .e19/.e16, v5 = (.e19 * 
        .e6/v5 - 1/.e20) * .e13/.e16, v6 = ((.e21/(v6 * .e20) - 
        (.e21/(v6 * .e9^.e1) - .e6/(v5 * .e20))/.e20)/v6 - .e10 * 
        .e6/(v5 * .e12)) * .e13/v5)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p1c_fdd=function (x, ta, tb, tc, v1, v2, v3, v4, v5, v6) 
{
    .e1 <- 1/v6
    .e6 <- x - (ta * v2 + tb * v3 + tc * v4 + v1)
    .e8 <- v6 * .e6/v5
    .e9 <- 1 + .e8
    .e10 <- 1 + .e1
    .e11 <- .e1 + 2
    .e12 <- .e9^.e10
    .e13 <- 2 * .e10
    .e14 <- .e9^.e11
    .e15 <- v6 * .e10
    .e16 <- log1p(.e8)
    .e17 <- .e9^.e13
    .e18 <- .e9^.e1
    .e19 <- 1/.e17
    .e21 <- exp(-.e9^-.e1)
    .e23 <- .e15/.e14 - .e19
    .e24 <- 2 * .e11
    .e25 <- v6 * .e18
    .e26 <- v5 * .e12
    .e30 <- v5^3
    .e33 <- v6 * .e9^(.e10 - .e24) * .e11 - 2/.e9^(1 + .e13)
    .e36 <- .e15 * .e33 - .e23/.e12
    .e39 <- .e16/.e25 - .e6/.e26
    .e40 <- v6 * .e12
    .e41 <- v5 * .e14
    .e42 <- 1/.e12
    .e43 <- v6^2
    .e44 <- 1/.e14
    .e45 <- .e9^(.e1 - .e13)
    .e46 <- .e39/.e12
    .e47 <- .e16/.e40
    .e48 <- v5^2
    .e50 <- .e10 * .e6/.e41
    .e51 <- (.e47 - .e46)/v6
    .e54 <- .e23 * .e6/v5 - .e42
    .e58 <- .e12 * .e11 * .e6/v5 - .e14 * .e16/.e43
    .e59 <- .e9^(.e1 - 1)
    .e60 <- .e51 - .e50
    .e61 <- .e26^2
    .e62 <- .e41^2
    .e63 <- .e40^2
    .e64 <- .e25^2
    .e76 <- (2 * (.e10 * .e9^(.e13 - 1) * .e6/v5) - 2 * (.e17 * 
        .e16/.e43))/.e9^(4 * .e10) + .e44
    .e81 <- .e19 + v6 * (.e33 * .e6/v5 - (.e45 + .e44)) * .e10 - 
        .e54/.e12
    .e84 <- v6 * .e58 * .e10/.e9^.e24
    .e86 <- .e15 * .e18 * .e6
    .e88 <- .e15 * (.e25 * .e16/.e63 - .e45 * .e39) - .e44
    .e90 <- .e40 * .e11 * .e6
    .e93 <- v6 * .e59 * .e16/.e64
    .e98 <- (.e88/v5 - ((.e42 + .e93 - .e42)/v5 - .e86/.e61)/.e12)/v6 - 
        (.e60/.e26 + .e10 * (.e90/.e62 - 1/.e41))
    .e102 <- .e10 * .e18 * .e6/v5 - .e12 * .e16/.e43
    .e103 <- .e76 - (.e39 * .e23/v6 + .e84)
    .e106 <- .e36 * .e6/v5 - 2 * .e23
    .e107 <- .e102/.e17
    .e110 <- ta * .e21 * .e36/.e30
    .e114 <- ta * tb * .e21 * .e36/.e30
    .e118 <- ta * tc * .e21 * .e36/.e30
    .e121 <- tb * .e21 * .e36/.e30
    .e125 <- tb * tc * .e21 * .e36/.e30
    .e128 <- tc * .e21 * .e36/.e30
    c(v1 = c(v1 = .e21 * .e36/.e30, v2 = .e110, v3 = .e121, v4 = .e128, 
        v5 = .e81 * .e21/.e30, v6 = .e98 * .e21/v5), v2 = c(v1 = .e110, 
        v2 = ta^2 * .e21 * .e36/.e30, v3 = .e114, v4 = .e118, 
        v5 = ta * .e81 * .e21/.e30, v6 = ta * .e98 * .e21/v5), 
        v3 = c(v1 = .e121, v2 = .e114, v3 = tb^2 * .e21 * .e36/.e30, 
            v4 = .e125, v5 = tb * .e81 * .e21/.e30, v6 = tb * 
                .e98 * .e21/v5), v4 = c(v1 = .e128, v2 = .e118, 
            v3 = .e125, v4 = tc^2 * .e21 * .e36/.e30, v5 = tc * 
                .e81 * .e21/.e30, v6 = tc * .e98 * .e21/v5), 
        v5 = c(v1 = .e106 * .e21/.e30, v2 = ta * .e106 * .e21/.e30, 
            v3 = tb * .e106 * .e21/.e30, v4 = tc * .e106 * .e21/.e30, 
            v5 = (.e81 * .e6/v5 - 2 * .e54) * .e21/.e30, v6 = (((.e14 - 
                .e90/v5) * .e10/.e62 + (.e88/.e48 - ((.e12 - 
                .e86/v5)/.e61 + (.e93 - .e42)/.e48)/.e12)/v6 - 
                .e60/(.e48 * .e12)) * .e6 - .e60/v5) * .e21/v5), 
        v6 = c(v1 = .e103 * .e21/.e48, v2 = ta * .e103 * .e21/.e48, 
            v3 = tb * .e103 * .e21/.e48, v4 = tc * .e103 * .e21/.e48, 
            v5 = (.e107 + (.e76 - .e84) * .e6/v5 - .e54 * .e39/v6) * 
                .e21/.e48, v6 = (((.e107 + .e50 - .e51) * .e39 + 
                (.e46 + .e6/.e41 - .e47)/v6 - ((.e12 + v6 * .e102) * 
                .e16/.e63 + ((1/(v5 * v6 * .e12) + v5 * .e102/.e61) * 
                .e6 - (.e59 * .e6/v5 + .e18 - .e18 * .e16/v6) * 
                .e16/.e64)/.e12))/v6 + (1/(v5 * .e43 * .e14) + 
                v5 * .e58 * .e10/.e62) * .e6) * .e21/v5))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_p1c_pd=function (x, ta, tb, tc, v1, v2, v3, v4, v5, v6) 
{
    .e5 <- x - (ta * v2 + tb * v3 + tc * v4 + v1)
    .e7 <- v6 * .e5/v5
    .e8 <- 1 + .e7
    .e9 <- 1/v6
    .e11 <- .e8^(1 + .e9)
    .e12 <- exp(-.e8^-.e9)
    .e13 <- v5 * .e11
    c(v1 = -(.e12/.e13), v2 = -(ta * .e12/.e13), v3 = -(tb * 
        .e12/.e13), v4 = -(tc * .e12/.e13), v5 = -(.e12 * .e5/(v5^2 * 
        .e11)), v6 = -(.e12 * (log1p(.e7)/(v6 * .e8^.e9) - .e5/.e13)/v6))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p1c_pdd=function (x, ta, tb, tc, v1, v2, v3, v4, v5, v6) 
{
    .e5 <- x - (ta * v2 + tb * v3 + tc * v4 + v1)
    .e6 <- 1/v6
    .e8 <- v6 * .e5/v5
    .e9 <- 1 + .e8
    .e10 <- 1 + .e6
    .e11 <- .e9^.e10
    .e12 <- .e9^.e6
    .e13 <- v5 * .e11
    .e15 <- exp(-.e9^-.e6)
    .e16 <- .e13^2
    .e17 <- v5^2
    .e19 <- v6 * .e10 * .e12
    .e20 <- log1p(.e8)
    .e21 <- .e9^(2 * .e10)
    .e22 <- v6 * .e12
    .e23 <- .e5/.e13
    .e26 <- .e19/.e16 - 1/(.e17 * .e21)
    .e27 <- .e20/.e22
    .e28 <- .e27 - .e23
    .e29 <- .e17 * .e11
    .e30 <- .e19 * .e5
    .e31 <- 1/.e11
    .e32 <- v5 * v6
    .e36 <- .e10 * .e12 * .e5/v5 - .e11 * .e20/v6^2
    .e37 <- .e9^(.e6 - 1)
    .e38 <- .e29^2
    .e39 <- .e22^2
    .e40 <- (.e11 - .e30/v5)/.e16
    .e42 <- .e28/.e11 + .e31
    .e44 <- v5 * .e36/.e16
    .e45 <- .e32 * .e11
    .e48 <- v6 * .e37 * .e20/.e39
    .e49 <- .e40 + .e5/(v5^3 * .e21)
    .e52 <- (.e31 + .e48 - .e42)/v5 - .e30/.e16
    .e54 <- .e28/.e45 + .e44
    .e59 <- .e32 * .e10 * .e12 * .e5/.e38 - (.e23 + 1)/.e29
    .e60 <- ta * .e15
    .e61 <- tb * .e15
    .e62 <- tc * .e15
    .e63 <- -(.e60 * .e26)
    .e64 <- -(ta * tb * .e15 * .e26)
    .e65 <- -(ta * tc * .e15 * .e26)
    .e66 <- -(.e61 * .e26)
    .e67 <- -(tb * tc * .e15 * .e26)
    .e68 <- -(.e62 * .e26)
    c(v1 = c(v1 = -(.e15 * .e26), v2 = .e63, v3 = .e66, v4 = .e68, 
        v5 = -(.e15 * .e59), v6 = -(.e52 * .e15/v6)), v2 = c(v1 = .e63, 
        v2 = -(ta^2 * .e15 * .e26), v3 = .e64, v4 = .e65, v5 = -(.e60 * 
            .e59), v6 = -(ta * .e52 * .e15/v6)), v3 = c(v1 = .e66, 
        v2 = .e64, v3 = -(tb^2 * .e15 * .e26), v4 = .e67, v5 = -(.e61 * 
            .e59), v6 = -(tb * .e52 * .e15/v6)), v4 = c(v1 = .e68, 
        v2 = .e65, v3 = .e67, v4 = -(tc^2 * .e15 * .e26), v5 = -(.e62 * 
            .e59), v6 = -(tc * .e52 * .e15/v6)), v5 = c(v1 = .e49 * 
        .e15, v2 = ta * .e49 * .e15, v3 = tb * .e49 * .e15, v4 = tc * 
        .e49 * .e15, v5 = ((2 * .e13 - .e30)/.e38 + .e5/(v5^4 * 
        .e21)) * .e15 * .e5, v6 = -((.e40 + (.e48 - .e42)/.e17) * 
        .e15 * .e5/v6)), v6 = c(v1 = .e54 * .e15, v2 = ta * .e54 * 
        .e15, v3 = tb * .e54 * .e15, v4 = tc * .e54 * .e15, v5 = (.e28/(.e17 * 
        v6 * .e11) + .e17 * .e36/.e38) * .e15 * .e5, v6 = -(((1/.e45 + 
        .e44) * .e5 - ((.e37 * .e5/v5 + .e12 - .e12 * .e20/v6) * 
        .e20/.e39 + (1 + .e27 - .e23) * .e28/v6)) * .e15/v6)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p1c_logfdd=function (x, ta, tb, tc, v1, v2, v3, v4, v5, v6) 
{
    .e5 <- x - (ta * v2 + tb * v3 + tc * v4 + v1)
    .e6 <- v6 * .e5
    .e7 <- .e6/v5
    .e8 <- 1 + .e7
    .e9 <- 1/v6
    .e10 <- .e8^.e9
    .e11 <- 1 + .e9
    .e12 <- v5 * .e8
    .e13 <- 1/.e10
    .e14 <- v6 * .e11
    .e15 <- .e12^2
    .e16 <- .e14 - .e13
    .e17 <- .e8^(.e9 + 2)
    .e18 <- v6 * .e16
    .e19 <- v5^2
    .e22 <- .e18/.e15 - 1/(.e19 * .e17)
    .e23 <- log1p(.e7)
    .e24 <- .e8^.e11
    .e25 <- v5 * .e24
    .e29 <- .e8^(.e9 - 1) * .e5/v5 - .e10 * .e23/v6
    .e30 <- .e25^2
    .e31 <- .e16 * .e5
    .e32 <- .e13 - 1
    .e33 <- 2/v6
    .e35 <- (.e29/(v6 * .e8^.e33) + 1)/.e12 - .e31/.e15
    .e37 <- .e16/.e15 + .e5/(v5^3 * .e17)
    .e39 <- .e23/.e10 - v6 * .e32
    .e41 <- .e14 * .e10 * .e5
    .e43 <- ((.e39/v6 + .e13)/.e12 - .e41/.e30)/v6 + .e11 * (.e6/.e15 - 
        1/.e12)
    .e47 <- .e18 * .e5/.e15 - (.e5/.e25 + .e14 - .e13)/.e12
    .e48 <- ta * .e22
    .e50 <- ta * tb * .e22
    .e52 <- ta * tc * .e22
    .e53 <- tb * .e22
    .e55 <- tb * tc * .e22
    .e56 <- tc * .e22
    .e57 <- v6^2
    c(v1 = c(v1 = .e22, v2 = .e48, v3 = .e53, v4 = .e56, v5 = .e47/v5, 
        v6 = -.e43), v2 = c(v1 = .e48, v2 = ta^2 * .e22, v3 = .e50, 
        v4 = .e52, v5 = ta * .e47/v5, v6 = -(ta * .e43)), v3 = c(v1 = .e53, 
        v2 = .e50, v3 = tb^2 * .e22, v4 = .e55, v5 = tb * .e47/v5, 
        v6 = -(tb * .e43)), v4 = c(v1 = .e56, v2 = .e52, v3 = .e55, 
        v4 = tc^2 * .e22, v5 = tc * .e47/v5, v6 = -(tc * .e43)), 
        v5 = c(v1 = -.e37, v2 = -(ta * .e37), v3 = -(tb * .e37), 
            v4 = -(tc * .e37), v5 = -(((.e31/.e12 - 1)/v5 + .e37 * 
                .e5)/v5), v6 = -((((.e24 - .e41/v5)/.e30 + .e39/(.e19 * 
                v6 * .e8))/v6 - .e11/.e15) * .e5)), v6 = c(v1 = .e35, 
            v2 = ta * .e35, v3 = tb * .e35, v4 = tc * .e35, v5 = .e35 * 
                .e5/v5, v6 = -(((((2/.e10 - 1) * .e5/v5 - .e29 * 
                .e23/(v6 * .e8^(.e33 - 1)))/.e8 - 2 * (.e32 * 
                .e23/v6))/v6 + v5 * (.e11 * .e10 * .e5/v5 - .e24 * 
                .e23/.e57) * .e5/.e30)/v6 - (.e11 * .e5/.e15 + 
                1/(v5 * .e57 * .e8)) * .e5)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
gev_p1c_logfddd=function (x, ta, tb, tc, v1, v2, v3, v4, v5, v6) 
{
    .e5 <- x - (ta * v2 + tb * v3 + tc * v4 + v1)
    .e6 <- v6 * .e5
    .e7 <- .e6/v5
    .e8 <- 1 + .e7
    .e9 <- 1/v6
    .e10 <- 1 + .e9
    .e11 <- .e8^.e9
    .e12 <- v5 * .e8
    .e13 <- .e8^.e10
    .e14 <- .e12^2
    .e15 <- v5 * .e13
    .e16 <- log1p(.e7)
    .e17 <- .e9 + 2
    .e18 <- 1/.e11
    .e19 <- .e9 - 1
    .e20 <- .e8^.e19
    .e21 <- v6 * .e10
    .e22 <- 2/v6
    .e23 <- v5^2
    .e24 <- .e21 - .e18
    .e25 <- .e8^.e17
    .e26 <- .e20 * .e5
    .e27 <- .e8^.e22
    .e28 <- .e26/v5
    .e29 <- v5 * v6
    .e31 <- .e11 * .e16/v6
    .e32 <- .e28 - .e31
    .e33 <- .e29 * .e8
    .e34 <- .e15^2
    .e35 <- (.e23 * .e25)^2
    .e36 <- .e12 * .e24
    .e37 <- v6 * .e27
    .e38 <- (2 * (.e33 * .e24/.e14) - 1/.e15)/.e14
    .e39 <- v5^3
    .e40 <- 2 * (.e36 * .e5/.e14)
    .e41 <- .e38 - .e15 * .e17/.e35
    .e42 <- v6^2
    .e43 <- .e5/.e15
    .e44 <- .e22 - 1
    .e45 <- .e39 * .e25
    .e46 <- .e18 - 1
    .e47 <- .e8^.e44
    .e49 <- .e16/.e11 - v6 * .e46
    .e50 <- .e32/.e27
    .e51 <- .e23 * .e8
    .e52 <- .e20 * .e16
    .e53 <- .e8^(.e9 - 2)
    .e55 <- .e21 * .e11 * .e5
    .e56 <- 2/.e11
    .e57 <- .e23 * v6
    .e58 <- .e49/v6
    .e59 <- .e23 * .e13
    .e63 <- v6 * ((.e50 + 1)/v6 + 2 - .e40)
    .e65 <- (.e52 + v6 * .e20)/v6 - (.e20 + v6 * .e53 * .e19 * 
        .e5/v5)
    .e66 <- .e45^2
    .e67 <- .e37^2
    .e69 <- .e5/.e59 + 2 * (.e36/.e14)
    .e70 <- v6 * .e32
    .e71 <- 2 * (.e33 * .e5/.e14)
    .e72 <- .e10 * .e11
    .e74 <- .e65/.e37 + 2 * (.e70 * .e47/.e67)
    .e75 <- .e32/.e37
    .e76 <- .e58 + .e18
    .e79 <- .e16/.e13 - 2 * (v6/.e13)
    .e80 <- .e76/.e14
    .e81 <- ta * tb
    .e82 <- .e29 * .e10
    .e86 <- .e13 * .e17 * .e5/v5 - .e25 * .e16/.e42
    .e87 <- v5 * .e25
    .e90 <- v6 * .e13 * .e17 * .e5
    .e91 <- .e74/.e51
    .e92 <- .e75 + 1
    .e93 <- .e13 - .e55/v5
    .e94 <- .e38 + .e57 * .e13 * .e17 * .e5/.e66
    .e96 <- .e79/v6 + 1/.e13
    .e97 <- .e71 - 2
    .e99 <- v6 * .e69/.e14
    .e100 <- .e13 * .e16
    .e101 <- .e100/.e42
    .e102 <- ta * tc
    .e103 <- tb * tc
    .e105 <- .e72 * .e5/v5
    .e106 <- .e91 + (.e43 + .e63 - .e18)/.e14
    .e108 <- (.e96/.e51 + v6 * (.e80 - ((2 * (.e82 * .e8^(1 + 
        3/v6)/.e34) - .e20/v5) * .e5 - .e11) * .e10/.e34))/v6 + 
        .e21 * .e97/.e14
    .e109 <- .e105 - .e101
    .e110 <- .e8^(1 + .e22)
    .e112 <- (2 * .e87 - .e90)/.e35 - .e99
    .e113 <- .e94 - 1/.e45
    .e115 <- (.e63 - .e18)/.e14 + .e23 * .e86/.e35
    .e123 <- v6 * (.e56 + v6 * (.e40 - (2 + .e22)) - 2 * .e43)/.e14 - 
        (.e55/.e34 - 2/.e15)/.e12
    .e124 <- ta^2
    .e125 <- tb^2
    .e126 <- tc^2
    .e127 <- v6 * .e47
    .e128 <- .e8^(.e22 - 2)
    .e130 <- .e43 + .e21 - .e18
    .e131 <- .e57 * .e8
    .e132 <- .e69/.e14
    .e133 <- .e28 + .e11
    .e134 <- .e93/.e34
    .e135 <- 1/.e27
    .e136 <- 2 * .e58
    .e137 <- v5 * .e109
    .e138 <- .e39 * .e8
    .e143 <- .e32 * .e16/.e127
    .e145 <- (v6 * .e24 * .e5/.e14 - .e130/.e12)/v5
    .e147 <- .e43 + v6 * ((.e50 + 2)/v6 + 3 - .e40) - .e56
    .e148 <- 2 * (v5 * .e93 * .e110/.e34)
    .e149 <- .e11 + .e31
    .e152 <- (.e56 - 1) * .e5/v5 - .e143
    .e154 <- .e92/.e12
    .e157 <- (.e53 * .e19 * .e5/v5 - .e52/.e42) * .e5/v5 - ((.e28 - 
        .e149) * .e16/v6 + .e28)/v6
    .e158 <- .e12^4
    .e159 <- (v5 * .e42 * .e8)^2
    .e160 <- .e131^2
    .e161 <- .e127^2
    .e164 <- .e81 * tc * v6 * .e41
    .e166 <- .e81 * v6 * .e41
    .e168 <- .e102 * v6 * .e41
    .e170 <- .e103 * v6 * .e41
    .e171 <- v6 * .e128
    .e172 <- (.e74/.e138 + .e132) * .e5
    .e176 <- (.e157/.e37 - .e32 * (.e27 + 2 * (.e47 * .e5/v5) - 
        2 * (.e27 * .e16/v6))/.e67)/.e12 - (2 + 2 * .e75 - .e40) * 
        .e5/.e14
    .e177 <- (.e92 - .e40)/.e14
    .e178 <- .e92/.e14
    .e179 <- (.e72 * .e16 + .e11)/v6
    .e180 <- .e132 + v5 * (3 * .e87 - .e90) * .e5/.e66
    .e181 <- .e133 - .e11
    .e182 <- (2 * (.e5/(v5 * .e11)) + v6 * .e152 - ((.e65 * .e16 - 
        .e70/.e8)/.e171 + 2/.e20))/.e8
    .e183 <- .e24/.e14
    .e184 <- .e136 + .e42 * .e32 * .e128 * .e44 * .e16/.e161
    .e185 <- .e16/.e37
    .e187 <- .e137 * .e5/.e34
    .e190 <- .e39 * .e86 * .e5/.e66
    .e192 <- (((.e179 - .e133 * .e10) * .e5 + .e137 * (2 * (.e82 * 
        .e110 * .e5/.e34) - 1))/.e34 + (.e182 + 2 - .e184)/.e33)/v6 - 
        (.e10 * .e97/.e14 + v6^3/.e159) * .e5
    .e194 <- (((.e32 * (.e135 - (.e135 + .e185)) + .e43 + 2 - 
        (.e136 + .e56))/.e12 + v6 * (.e72/.e34 - 1/.e14) * .e5)/v6 - 
        (((.e28 - (.e31 + 2 * (.e57 * .e109 * .e110/.e34))) * 
            .e10 + .e11)/.e34 + .e80) * .e5)/v6 + .e10 * (2 - 
        .e71) * .e5/.e14
    .e196 <- (.e91 + .e147/.e14) * .e5 - .e154
    .e197 <- .e172 - .e178
    .e200 <- ((.e134 + 1/.e59)/.e12 - .e99) * .e5 + .e130/.e14 - 
        .e145
    .e202 <- ((.e96/.e138 + v6 * (.e26/.e23 + .e148) * .e10/.e34) * 
        .e5 - .e80)/v6 + (1 - .e71) * .e10/.e14
    .e204 <- .e177 - .e190
    .e211 <- ((.e79/(.e39 * v6 * .e8) + v6 * ((.e181/v5 + .e148) * 
        .e10/.e34 + .e29 * .e49/.e160))/v6 - 2 * (.e82 * .e8/.e158)) * 
        .e5 + .e10/.e14 - (.e134 + .e49/.e131)/v6
    .e214 <- (.e94 - 2/.e45) * .e5 + .e145 - .e183
    .e217 <- .e147 * .e5/.e14 - (.e92 - .e187)/.e12
    .e219 <- v5 * .e10 * .e8
    .e220 <- .e24 * .e5
    .e223 <- ta * .e125 * v6 * .e41
    .e226 <- ta * .e126 * v6 * .e41
    .e228 <- ta * v6 * .e41
    .e231 <- .e124 * tb * v6 * .e41
    .e234 <- .e124 * tc * v6 * .e41
    .e236 <- .e124 * v6 * .e41
    .e239 <- tb * .e126 * v6 * .e41
    .e241 <- tb * v6 * .e41
    .e244 <- .e125 * tc * v6 * .e41
    .e246 <- .e125 * v6 * .e41
    .e248 <- tc * v6 * .e41
    .e250 <- .e126 * v6 * .e41
    .e251 <- .e219 * .e5
    .e253 <- -(ta * .e113)
    .e255 <- -(.e81 * .e113)
    .e257 <- -(.e102 * .e113)
    .e259 <- -(tb * .e113)
    .e261 <- -(.e103 * .e113)
    .e263 <- -(tc * .e113)
    .e264 <- (.e154 - .e220/.e14)/v5
    .e265 <- .e101 + 2 * (.e23 * .e109 * .e93 * .e13/.e34)
    .e267 <- .e46 * .e16/v6
    .e269 <- 1/.e42 + 2 * (.e251/.e14)
    .e270 <- 2 * .e12
    .e271 <- c(v1 = .e166, v2 = .e231, v3 = .e223, v4 = .e164, 
        v5 = .e81 * .e123/v5, v6 = -(.e81 * .e108))
    .e272 <- c(v1 = .e168, v2 = .e234, v3 = .e164, v4 = .e226, 
        v5 = .e102 * .e123/v5, v6 = -(.e102 * .e108))
    .e273 <- c(v1 = .e228, v2 = .e236, v3 = .e166, v4 = .e168, 
        v5 = ta * .e123/v5, v6 = -(ta * .e108))
    .e274 <- c(v1 = .e170, v2 = .e164, v3 = .e244, v4 = .e239, 
        v5 = .e103 * .e123/v5, v6 = -(.e103 * .e108))
    .e275 <- c(v1 = .e241, v2 = .e166, v3 = .e246, v4 = .e170, 
        v5 = tb * .e123/v5, v6 = -(tb * .e108))
    .e276 <- c(v1 = .e248, v2 = .e168, v3 = .e170, v4 = .e250, 
        v5 = tc * .e123/v5, v6 = -(tc * .e108))
    .e277 <- ta * .e106
    .e278 <- ta * .e112
    .e279 <- ta * .e115
    .e280 <- .e81 * .e106
    .e281 <- .e81 * .e112
    .e282 <- .e81 * .e115
    .e283 <- .e102 * .e106
    .e284 <- .e102 * .e112
    .e285 <- .e102 * .e115
    .e286 <- tb * .e106
    .e287 <- tb * .e112
    .e288 <- tb * .e115
    .e289 <- .e103 * .e106
    .e290 <- .e103 * .e112
    .e291 <- .e103 * .e115
    .e292 <- tc * .e106
    .e293 <- tc * .e112
    .e294 <- tc * .e115
    c(v1 = c(v1 = c(v1 = v6 * .e41, v2 = .e228, v3 = .e241, v4 = .e248, 
        v5 = .e123/v5, v6 = -.e108), v2 = .e273, v3 = .e275, 
        v4 = .e276, v5 = c(v1 = -.e113, v2 = .e253, v3 = .e259, 
            v4 = .e263, v5 = -(.e214/v5), v6 = -.e211), v6 = c(v1 = .e106, 
            v2 = .e277, v3 = .e286, v4 = .e292, v5 = .e196/v5, 
            v6 = -.e192)), v2 = c(v1 = .e273, v2 = c(v1 = .e236, 
        v2 = ta^3 * v6 * .e41, v3 = .e231, v4 = .e234, v5 = .e124 * 
            .e123/v5, v6 = -(.e124 * .e108)), v3 = .e271, v4 = .e272, 
        v5 = c(v1 = .e253, v2 = -(.e124 * .e113), v3 = .e255, 
            v4 = .e257, v5 = -(ta * .e214/v5), v6 = -(ta * .e211)), 
        v6 = c(v1 = .e277, v2 = .e124 * .e106, v3 = .e280, v4 = .e283, 
            v5 = ta * .e196/v5, v6 = -(ta * .e192))), v3 = c(v1 = .e275, 
        v2 = .e271, v3 = c(v1 = .e246, v2 = .e223, v3 = tb^3 * 
            v6 * .e41, v4 = .e244, v5 = .e125 * .e123/v5, v6 = -(.e125 * 
            .e108)), v4 = .e274, v5 = c(v1 = .e259, v2 = .e255, 
            v3 = -(.e125 * .e113), v4 = .e261, v5 = -(tb * .e214/v5), 
            v6 = -(tb * .e211)), v6 = c(v1 = .e286, v2 = .e280, 
            v3 = .e125 * .e106, v4 = .e289, v5 = tb * .e196/v5, 
            v6 = -(tb * .e192))), v4 = c(v1 = .e276, v2 = .e272, 
        v3 = .e274, v4 = c(v1 = .e250, v2 = .e226, v3 = .e239, 
            v4 = tc^3 * v6 * .e41, v5 = .e126 * .e123/v5, v6 = -(.e126 * 
                .e108)), v5 = c(v1 = .e263, v2 = .e257, v3 = .e261, 
            v4 = -(.e126 * .e113), v5 = -(tc * .e214/v5), v6 = -(tc * 
                .e211)), v6 = c(v1 = .e292, v2 = .e283, v3 = .e289, 
            v4 = .e126 * .e106, v5 = tc * .e196/v5, v6 = -(tc * 
                .e192))), v5 = c(v1 = c(v1 = .e112, v2 = .e278, 
        v3 = .e287, v4 = .e293, v5 = .e200/v5, v6 = -.e202), 
        v2 = c(v1 = .e278, v2 = .e124 * .e112, v3 = .e281, v4 = .e284, 
            v5 = ta * .e200/v5, v6 = -(ta * .e202)), v3 = c(v1 = .e287, 
            v2 = .e281, v3 = .e125 * .e112, v4 = .e290, v5 = tb * 
                .e200/v5, v6 = -(tb * .e202)), v4 = c(v1 = .e293, 
            v2 = .e284, v3 = .e290, v4 = .e126 * .e112, v5 = tc * 
                .e200/v5, v6 = -(tc * .e202)), v5 = c(v1 = .e180, 
            v2 = ta * .e180, v3 = tb * .e180, v4 = tc * .e180, 
            v5 = (.e180 * .e5 + 2 * (((.e220/.e12 - 1)/v5 + (.e183 + 
                .e5/.e45) * .e5)/v5))/v5, v6 = -(((.e79 * .e5/(v5^4 * 
                v6 * .e8) + (v6 * .e181 * .e10 * .e5/.e23 - 2 * 
                (v5 * .e93^2 * .e13/.e34))/.e34 - v6 * (.e270 - 
                .e6) * .e49/.e160)/v6 + 2 * (.e219/.e158)) * 
                .e5)), v6 = c(v1 = .e197, v2 = ta * .e197, v3 = tb * 
            .e197, v4 = tc * .e197, v5 = (.e172 - (.e264 + .e178)) * 
            .e5/v5, v6 = -(((((.e179 + (.e11 - .e133) * .e10) * 
            .e5/v5 - .e265)/.e34 + (.e182 + 1 - .e184)/.e131)/v6 + 
            2 * (.e251/.e158) + .e42/.e159) * .e5))), v6 = c(v1 = c(v1 = .e115, 
        v2 = .e279, v3 = .e288, v4 = .e294, v5 = .e217/v5, v6 = -.e194), 
        v2 = c(v1 = .e279, v2 = .e124 * .e115, v3 = .e282, v4 = .e285, 
            v5 = ta * .e217/v5, v6 = -(ta * .e194)), v3 = c(v1 = .e288, 
            v2 = .e282, v3 = .e125 * .e115, v4 = .e291, v5 = tb * 
                .e217/v5, v6 = -(tb * .e194)), v4 = c(v1 = .e294, 
            v2 = .e285, v3 = .e291, v4 = .e126 * .e115, v5 = tc * 
                .e217/v5, v6 = -(tc * .e194)), v5 = c(v1 = -.e204, 
            v2 = -(ta * .e204), v3 = -(tb * .e204), v4 = -(tc * 
                .e204), v5 = -((.e264 + .e177 - .e190) * .e5/v5), 
            v6 = -(((((.e32 * (.e135 - .e185) + .e43 + 1 - .e76)/.e51 - 
                .e134)/v6 + (((.e149 - .e28) * .e10 - .e11) * 
                .e5/v5 - .e265)/.e34 - v5 * .e49 * (.e12 + .e6)/.e160)/v6 + 
                .e269/.e14) * .e5)), v6 = c(v1 = .e176, v2 = ta * 
            .e176, v3 = tb * .e176, v4 = tc * .e176, v5 = .e176 * 
            .e5/v5, v6 = -((((((.e47 + v6 * (.e128 * .e44 * .e5/v5 - 
            2 * (.e47 * .e16/.e42))) * .e16/.e161 - 2 * (.e5/(.e29 * 
            .e27))) * .e32 - ((.e157 * .e16 + .e32 * .e5/.e12)/.e171 + 
            .e152 * .e5/v5)/.e8)/.e8 - ((2 * ((.e46 * .e5/v5 - 
            .e143)/.e8 - .e267) + 2 * (.e152/.e8) - 4 * .e267)/v6 + 
            .e187))/v6 + v5 * (((.e32 * .e10 - .e11/v6) * .e5/v5 - 
            ((.e105 - (.e100/v6 + 2 * .e13)/v6) * .e16 + .e11 * 
                .e5/v5)/v6)/v6 - 2 * (.e23 * .e109^2 * .e13/.e34)) * 
            .e5/.e34)/v6 + (.e269 * .e5/.e14 + v6 * (.e270 + 
            .e6)/.e159) * .e5))))
}
############################################################
#' The first derivative of the density for DMGS
#' @returns Vector
#' @inheritParams manf
gev_p1c_f1fa=function(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6){

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_fd,"x")
	f1=vf(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6)
	return(f1)
}
############################################################
#' The first derivative of the density for WAIC
#' @returns Vector
#' @inheritParams manf
gev_p1c_f1fw=function(x,ta,tb,tc,v1,v2,v3,v4,v5,v6){

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_fd,c("x","ta","tb","tc"))
	f1=vf(x,ta,tb,tc,v1,v2,v3,v4,v5,v6)
	return(f1)
}
############################################################
#' The second derivative of the density for DMGS
#' @returns Matrix
#' @inheritParams manf
gev_p1c_f2fa=function(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6){
	nx=length(x)

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_fdd,"x")
	temp1=vf(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6)
	f2=deriv_copyfdd(temp1,nx,dim=6)
	return(f2)
}
############################################################
#' The second derivative of the density for WAIC
#' @returns Matrix
#' @inheritParams manf
gev_p1c_f2fw=function(x,ta,tb,tc,v1,v2,v3,v4,v5,v6){
	nx=length(x)

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_fdd,c("x","ta","tb","tc"))
	temp1=vf(x,ta,tb,tc,v1,v2,v3,v4,v5,v6)
	f2=deriv_copyfdd(temp1,nx,dim=6)
	return(f2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
gev_p1c_mu1fa=function(alpha,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6){
	x=qgev((1-alpha),mu=v1+v2*t0a+v3*t0b+v4*t0c,sigma=v5,xi=v6)

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_pd,"x")
	mu1=-vf(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
gev_p1c_mu2fa=function(alpha,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6){
	x=qgev((1-alpha),mu=v1+v2*t0a+v3*t0b+v4*t0c,sigma=v5,xi=v6)
	nx=length(x)

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_pdd,"x")
	temp1=vf(x,t0a,t0b,t0c,v1,v2,v3,v4,v5,v6)
	mu2=-deriv_copyfdd(temp1,nx,dim=6)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
gev_p1c_ldda=function(x,ta,tb,tc,v1,v2,v3,v4,v5,v6){
	nx=length(x)

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_logfdd,c("x","ta","tb","tc"))
	temp1=vf(x,ta,tb,tc,v1,v2,v3,v4,v5,v6)
	ldd=deriv_copyldd(temp1,nx,dim=6)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
gev_p1c_lddda=function(x,ta,tb,tc,v1,v2,v3,v4,v5,v6){
	nx=length(x)

	v6=movexiawayfromzero(v6)

	vf=Vectorize(gev_p1c_logfddd,c("x","ta","tb","tc"))
	temp1=vf(x,ta,tb,tc,v1,v2,v3,v4,v5,v6)
	lddd=deriv_copylddd(temp1,nx,dim=6)
	return(lddd)
}
