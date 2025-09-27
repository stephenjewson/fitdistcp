######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_p1b_fd=function (x, ta, tb, v1, v2, v3, v4, v5) 
{
    .e1 <- 1/v5
    .e5 <- x - (ta * v2 + tb * v3 + v1)
    .e7 <- v5 * .e5/v4
    .e8 <- 1 + .e7
    .e9 <- 1 + .e1
    .e11 <- .e8^(.e1 + 2)
    .e12 <- exp(-.e8^-.e1)
    .e13 <- .e8^.e9
    .e16 <- v4^2
    .e19 <- v5 * .e9/.e11 - 1/.e8^(2 * .e9)
    .e20 <- log1p(.e7)
    c(v1 = .e12 * .e19/.e16, v2 = ta * .e12 * .e19/.e16, v3 = tb * 
        .e12 * .e19/.e16, v4 = (.e19 * .e5/v4 - 1/.e13) * .e12/.e16, 
        v5 = ((.e20/(v5 * .e13) - (.e20/(v5 * .e8^.e1) - .e5/(v4 * 
            .e13))/.e13)/v5 - .e9 * .e5/(v4 * .e11)) * .e12/v4)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p1b_fdd=function (x, ta, tb, v1, v2, v3, v4, v5) 
{
    .e1 <- 1/v5
    .e5 <- x - (ta * v2 + tb * v3 + v1)
    .e7 <- v5 * .e5/v4
    .e8 <- 1 + .e7
    .e9 <- 1 + .e1
    .e10 <- .e1 + 2
    .e11 <- .e8^.e9
    .e12 <- 2 * .e9
    .e13 <- .e8^.e10
    .e14 <- log1p(.e7)
    .e15 <- v5 * .e9
    .e16 <- .e8^.e1
    .e17 <- .e8^.e12
    .e18 <- 1/.e17
    .e20 <- exp(-.e8^-.e1)
    .e21 <- v5 * .e16
    .e22 <- v4 * .e11
    .e24 <- .e15/.e13 - .e18
    .e25 <- 2 * .e10
    .e31 <- .e14/.e21 - .e5/.e22
    .e32 <- v4^3
    .e35 <- v5 * .e8^(.e9 - .e25) * .e10 - 2/.e8^(1 + .e12)
    .e36 <- v5 * .e11
    .e37 <- v4 * .e13
    .e38 <- v5^2
    .e39 <- 1/.e11
    .e41 <- 1/.e13
    .e43 <- .e15 * .e35 - .e24/.e11
    .e44 <- .e8^(.e1 - .e12)
    .e45 <- .e31/.e11
    .e46 <- .e14/.e36
    .e47 <- v4^2
    .e49 <- .e9 * .e5/.e37
    .e50 <- (.e46 - .e45)/v5
    .e53 <- .e24 * .e5/v4 - .e39
    .e57 <- .e11 * .e10 * .e5/v4 - .e13 * .e14/.e38
    .e58 <- .e8^(.e1 - 1)
    .e59 <- .e50 - .e49
    .e60 <- .e22^2
    .e61 <- .e37^2
    .e62 <- .e36^2
    .e63 <- .e21^2
    .e68 <- .e9 * .e16 * .e5/v4 - .e11 * .e14/.e38
    .e79 <- (2 * (.e9 * .e8^(.e12 - 1) * .e5/v4) - 2 * (.e17 * 
        .e14/.e38))/.e8^(4 * .e9) + .e41
    .e84 <- .e18 + v5 * (.e35 * .e5/v4 - (.e44 + .e41)) * .e9 - 
        .e53/.e11
    .e87 <- v5 * .e57 * .e9/.e8^.e25
    .e89 <- .e15 * .e16 * .e5
    .e91 <- .e15 * (.e21 * .e14/.e62 - .e44 * .e31) - .e41
    .e93 <- .e36 * .e10 * .e5
    .e96 <- v5 * .e58 * .e14/.e63
    .e101 <- (.e91/v4 - ((.e39 + .e96 - .e39)/v4 - .e89/.e60)/.e11)/v5 - 
        (.e59/.e22 + .e9 * (.e93/.e61 - 1/.e37))
    .e102 <- .e79 - (.e31 * .e24/v5 + .e87)
    .e105 <- .e43 * .e5/v4 - 2 * .e24
    .e106 <- .e68/.e17
    .e109 <- ta * .e20 * .e43/.e32
    .e113 <- ta * tb * .e20 * .e43/.e32
    .e116 <- tb * .e20 * .e43/.e32
    c(v1 = c(v1 = .e20 * .e43/.e32, v2 = .e109, v3 = .e116, v4 = .e84 * 
        .e20/.e32, v5 = .e101 * .e20/v4), v2 = c(v1 = .e109, 
        v2 = ta^2 * .e20 * .e43/.e32, v3 = .e113, v4 = ta * .e84 * 
            .e20/.e32, v5 = ta * .e101 * .e20/v4), v3 = c(v1 = .e116, 
        v2 = .e113, v3 = tb^2 * .e20 * .e43/.e32, v4 = tb * .e84 * 
            .e20/.e32, v5 = tb * .e101 * .e20/v4), v4 = c(v1 = .e105 * 
        .e20/.e32, v2 = ta * .e105 * .e20/.e32, v3 = tb * .e105 * 
        .e20/.e32, v4 = (.e84 * .e5/v4 - 2 * .e53) * .e20/.e32, 
        v5 = (((.e13 - .e93/v4) * .e9/.e61 + (.e91/.e47 - ((.e11 - 
            .e89/v4)/.e60 + (.e96 - .e39)/.e47)/.e11)/v5 - .e59/(.e47 * 
            .e11)) * .e5 - .e59/v4) * .e20/v4), v5 = c(v1 = .e102 * 
        .e20/.e47, v2 = ta * .e102 * .e20/.e47, v3 = tb * .e102 * 
        .e20/.e47, v4 = (.e106 + (.e79 - .e87) * .e5/v4 - .e53 * 
        .e31/v5) * .e20/.e47, v5 = (((.e106 + .e49 - .e50) * 
        .e31 + (.e45 + .e5/.e37 - .e46)/v5 - ((.e11 + v5 * .e68) * 
        .e14/.e62 + ((1/(v4 * v5 * .e11) + v4 * .e68/.e60) * 
        .e5 - (.e58 * .e5/v4 + .e16 - .e16 * .e14/v5) * .e14/.e63)/.e11))/v5 + 
        (1/(v4 * .e38 * .e13) + v4 * .e57 * .e9/.e61) * .e5) * 
        .e20/v4))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_p1b_pd=function (x, ta, tb, v1, v2, v3, v4, v5) 
{
    .e4 <- x - (ta * v2 + tb * v3 + v1)
    .e6 <- v5 * .e4/v4
    .e7 <- 1 + .e6
    .e8 <- 1/v5
    .e10 <- .e7^(1 + .e8)
    .e11 <- exp(-.e7^-.e8)
    .e12 <- v4 * .e10
    c(v1 = -(.e11/.e12), v2 = -(ta * .e11/.e12), v3 = -(tb * 
        .e11/.e12), v4 = -(.e11 * .e4/(v4^2 * .e10)), v5 = -(.e11 * 
        (log1p(.e6)/(v5 * .e7^.e8) - .e4/.e12)/v5))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p1b_pdd=function (x, ta, tb, v1, v2, v3, v4, v5) 
{
    .e4 <- x - (ta * v2 + tb * v3 + v1)
    .e5 <- 1/v5
    .e7 <- v5 * .e4/v4
    .e8 <- 1 + .e7
    .e9 <- 1 + .e5
    .e10 <- .e8^.e9
    .e11 <- .e8^.e5
    .e12 <- v4 * .e10
    .e14 <- exp(-.e8^-.e5)
    .e15 <- log1p(.e7)
    .e16 <- .e12^2
    .e17 <- v4^2
    .e19 <- v5 * .e9 * .e11
    .e20 <- v5 * .e11
    .e21 <- .e8^(2 * .e9)
    .e22 <- .e4/.e12
    .e23 <- .e15/.e20
    .e25 <- .e23 - .e22
    .e27 <- .e19/.e16 - 1/(.e17 * .e21)
    .e28 <- .e17 * .e10
    .e29 <- .e19 * .e4
    .e30 <- 1/.e10
    .e31 <- v4 * v5
    .e35 <- .e9 * .e11 * .e4/v4 - .e10 * .e15/v5^2
    .e36 <- .e8^(.e5 - 1)
    .e37 <- .e28^2
    .e38 <- .e20^2
    .e39 <- (.e10 - .e29/v4)/.e16
    .e41 <- .e25/.e10 + .e30
    .e43 <- v4 * .e35/.e16
    .e44 <- .e31 * .e10
    .e47 <- v5 * .e36 * .e15/.e38
    .e48 <- .e39 + .e4/(v4^3 * .e21)
    .e51 <- (.e30 + .e47 - .e41)/v4 - .e29/.e16
    .e53 <- .e25/.e44 + .e43
    .e54 <- ta * .e14
    .e55 <- tb * .e14
    .e60 <- .e31 * .e9 * .e11 * .e4/.e37 - (.e22 + 1)/.e28
    .e61 <- -(.e54 * .e27)
    .e62 <- -(ta * tb * .e14 * .e27)
    .e63 <- -(.e55 * .e27)
    c(v1 = c(v1 = -(.e14 * .e27), v2 = .e61, v3 = .e63, v4 = -(.e14 * 
        .e60), v5 = -(.e51 * .e14/v5)), v2 = c(v1 = .e61, v2 = -(ta^2 * 
        .e14 * .e27), v3 = .e62, v4 = -(.e54 * .e60), v5 = -(ta * 
        .e51 * .e14/v5)), v3 = c(v1 = .e63, v2 = .e62, v3 = -(tb^2 * 
        .e14 * .e27), v4 = -(.e55 * .e60), v5 = -(tb * .e51 * 
        .e14/v5)), v4 = c(v1 = .e48 * .e14, v2 = ta * .e48 * 
        .e14, v3 = tb * .e48 * .e14, v4 = ((2 * .e12 - .e29)/.e37 + 
        .e4/(v4^4 * .e21)) * .e14 * .e4, v5 = -((.e39 + (.e47 - 
        .e41)/.e17) * .e14 * .e4/v5)), v5 = c(v1 = .e53 * .e14, 
        v2 = ta * .e53 * .e14, v3 = tb * .e53 * .e14, v4 = (.e25/(.e17 * 
            v5 * .e10) + .e17 * .e35/.e37) * .e14 * .e4, v5 = -(((1/.e44 + 
            .e43) * .e4 - ((.e36 * .e4/v4 + .e11 - .e11 * .e15/v5) * 
            .e15/.e38 + (1 + .e23 - .e22) * .e25/v5)) * .e14/v5)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p1b_logfdd=function (x, ta, tb, v1, v2, v3, v4, v5) 
{
    .e4 <- x - (ta * v2 + tb * v3 + v1)
    .e5 <- v5 * .e4
    .e6 <- .e5/v4
    .e7 <- 1 + .e6
    .e8 <- 1/v5
    .e9 <- .e7^.e8
    .e10 <- 1 + .e8
    .e11 <- v4 * .e7
    .e12 <- 1/.e9
    .e13 <- v5 * .e10
    .e14 <- .e11^2
    .e15 <- .e13 - .e12
    .e16 <- .e7^(.e8 + 2)
    .e17 <- log1p(.e6)
    .e18 <- v5 * .e15
    .e19 <- .e7^.e10
    .e20 <- v4^2
    .e23 <- .e18/.e14 - 1/(.e20 * .e16)
    .e24 <- v4 * .e19
    .e28 <- .e7^(.e8 - 1) * .e4/v4 - .e9 * .e17/v5
    .e29 <- .e24^2
    .e30 <- .e15 * .e4
    .e31 <- .e12 - 1
    .e32 <- 2/v5
    .e34 <- (.e28/(v5 * .e7^.e32) + 1)/.e11 - .e30/.e14
    .e36 <- .e15/.e14 + .e4/(v4^3 * .e16)
    .e38 <- .e17/.e9 - v5 * .e31
    .e40 <- .e13 * .e9 * .e4
    .e42 <- ((.e38/v5 + .e12)/.e11 - .e40/.e29)/v5 + .e10 * (.e5/.e14 - 
        1/.e11)
    .e46 <- .e18 * .e4/.e14 - (.e4/.e24 + .e13 - .e12)/.e11
    .e47 <- ta * .e23
    .e49 <- ta * tb * .e23
    .e50 <- tb * .e23
    .e51 <- v5^2
    c(v1 = c(v1 = .e23, v2 = .e47, v3 = .e50, v4 = .e46/v4, v5 = -.e42), 
        v2 = c(v1 = .e47, v2 = ta^2 * .e23, v3 = .e49, v4 = ta * 
            .e46/v4, v5 = -(ta * .e42)), v3 = c(v1 = .e50, v2 = .e49, 
            v3 = tb^2 * .e23, v4 = tb * .e46/v4, v5 = -(tb * 
                .e42)), v4 = c(v1 = -.e36, v2 = -(ta * .e36), 
            v3 = -(tb * .e36), v4 = -(((.e30/.e11 - 1)/v4 + .e36 * 
                .e4)/v4), v5 = -((((.e19 - .e40/v4)/.e29 + .e38/(.e20 * 
                v5 * .e7))/v5 - .e10/.e14) * .e4)), v5 = c(v1 = .e34, 
            v2 = ta * .e34, v3 = tb * .e34, v4 = .e34 * .e4/v4, 
            v5 = -(((((2/.e9 - 1) * .e4/v4 - .e28 * .e17/(v5 * 
                .e7^(.e32 - 1)))/.e7 - 2 * (.e31 * .e17/v5))/v5 + 
                v4 * (.e10 * .e9 * .e4/v4 - .e19 * .e17/.e51) * 
                  .e4/.e29)/v5 - (.e10 * .e4/.e14 + 1/(v4 * .e51 * 
                .e7)) * .e4)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
gev_p1b_logfddd=function (x, ta, tb, v1, v2, v3, v4, v5) 
{
    .e4 <- x - (ta * v2 + tb * v3 + v1)
    .e5 <- v5 * .e4
    .e6 <- .e5/v4
    .e7 <- 1 + .e6
    .e8 <- 1/v5
    .e9 <- 1 + .e8
    .e10 <- .e7^.e8
    .e11 <- v4 * .e7
    .e12 <- .e7^.e9
    .e13 <- .e11^2
    .e14 <- log1p(.e6)
    .e15 <- .e8 - 1
    .e16 <- 1/.e10
    .e17 <- .e7^.e15
    .e18 <- v4 * .e12
    .e19 <- 2/v5
    .e20 <- .e8 + 2
    .e21 <- v5 * .e9
    .e22 <- v4^2
    .e23 <- .e21 - .e16
    .e24 <- .e17 * .e4
    .e25 <- .e24/v4
    .e26 <- .e7^.e20
    .e27 <- .e7^.e19
    .e29 <- .e10 * .e14/v5
    .e30 <- .e25 - .e29
    .e31 <- v4 * v5
    .e32 <- .e18^2
    .e33 <- v5 * .e27
    .e34 <- .e11 * .e23
    .e35 <- .e31 * .e7
    .e36 <- v4^3
    .e37 <- v5^2
    .e38 <- (.e22 * .e26)^2
    .e39 <- 2 * (.e34 * .e4/.e13)
    .e40 <- (2 * (.e35 * .e23/.e13) - 1/.e18)/.e13
    .e41 <- .e19 - 1
    .e42 <- .e4/.e18
    .e43 <- .e7^.e41
    .e44 <- .e16 - 1
    .e45 <- .e36 * .e26
    .e47 <- .e14/.e10 - v5 * .e44
    .e48 <- .e40 - .e18 * .e20/.e38
    .e49 <- .e17 * .e14
    .e50 <- .e7^(.e8 - 2)
    .e52 <- .e21 * .e10 * .e4
    .e53 <- .e30/.e27
    .e54 <- 2/.e10
    .e55 <- .e22 * v5
    .e56 <- .e47/v5
    .e57 <- .e22 * .e12
    .e58 <- .e22 * .e7
    .e59 <- .e9 * .e10
    .e60 <- .e30/.e33
    .e62 <- (.e49 + v5 * .e17)/v5 - (.e17 + v5 * .e50 * .e15 * 
        .e4/v4)
    .e63 <- .e45^2
    .e64 <- .e33^2
    .e66 <- .e4/.e57 + 2 * (.e34/.e13)
    .e67 <- v5 * .e30
    .e71 <- 2 * (.e35 * .e4/.e13)
    .e72 <- ta * tb
    .e73 <- v5 * ((.e53 + 1)/v5 + 2 - .e39)
    .e75 <- .e62/.e33 + 2 * (.e67 * .e43/.e64)
    .e76 <- .e60 + 1
    .e77 <- .e12 - .e52/v4
    .e78 <- .e12 * .e14
    .e79 <- .e56 + .e16
    .e82 <- .e14/.e12 - 2 * (v5/.e12)
    .e83 <- .e79/.e13
    .e84 <- .e78/.e37
    .e85 <- .e31 * .e9
    .e87 <- .e59 * .e4/v4
    .e88 <- .e87 - .e84
    .e92 <- .e12 * .e20 * .e4/v4 - .e26 * .e14/.e37
    .e93 <- v4 * .e26
    .e96 <- v5 * .e12 * .e20 * .e4
    .e97 <- .e75/.e58
    .e98 <- .e7^(1 + .e19)
    .e99 <- .e40 + .e55 * .e12 * .e20 * .e4/.e63
    .e101 <- .e82/v5 + 1/.e12
    .e102 <- .e71 - 2
    .e103 <- ta^2
    .e104 <- tb^2
    .e106 <- v5 * .e66/.e13
    .e107 <- v5 * .e43
    .e108 <- .e7^(.e19 - 2)
    .e109 <- .e97 + (.e42 + .e73 - .e16)/.e13
    .e111 <- (.e101/.e58 + v5 * (.e83 - ((2 * (.e85 * .e7^(1 + 
        3/v5)/.e32) - .e17/v4) * .e4 - .e10) * .e9/.e32))/v5 + 
        .e21 * .e102/.e13
    .e113 <- (2 * .e93 - .e96)/.e38 - .e106
    .e114 <- .e99 - 1/.e45
    .e116 <- (.e73 - .e16)/.e13 + .e22 * .e92/.e38
    .e119 <- .e42 + .e21 - .e16
    .e124 <- .e55 * .e7
    .e127 <- v5 * (.e54 + v5 * (.e39 - (2 + .e19)) - 2 * .e42)/.e13 - 
        (.e52/.e32 - 2/.e18)/.e11
    .e128 <- .e66/.e13
    .e129 <- .e25 + .e10
    .e130 <- .e77/.e32
    .e132 <- .e30 * .e14/.e107
    .e133 <- 1/.e27
    .e134 <- 2 * .e56
    .e135 <- v4 * .e88
    .e136 <- .e36 * .e7
    .e141 <- .e10 + .e29
    .e144 <- (.e54 - 1) * .e4/v4 - .e132
    .e145 <- (v5 * .e23 * .e4/.e13 - .e119/.e11)/v4
    .e147 <- .e42 + v5 * ((.e53 + 2)/v5 + 3 - .e39) - .e54
    .e148 <- 2 * (v4 * .e77 * .e98/.e32)
    .e150 <- .e72 * v5 * .e48
    .e152 <- .e76/.e11
    .e155 <- (.e50 * .e15 * .e4/v4 - .e49/.e37) * .e4/v4 - ((.e25 - 
        .e141) * .e14/v5 + .e25)/v5
    .e156 <- .e11^4
    .e157 <- (v4 * .e37 * .e7)^2
    .e158 <- .e124^2
    .e159 <- .e107^2
    .e160 <- v5 * .e108
    .e161 <- (.e75/.e136 + .e128) * .e4
    .e165 <- (.e155/.e33 - .e30 * (.e27 + 2 * (.e43 * .e4/v4) - 
        2 * (.e27 * .e14/v5))/.e64)/.e11 - (2 + 2 * .e60 - .e39) * 
        .e4/.e13
    .e166 <- (.e76 - .e39)/.e13
    .e167 <- .e76/.e13
    .e168 <- (.e59 * .e14 + .e10)/v5
    .e169 <- .e128 + v4 * (3 * .e93 - .e96) * .e4/.e63
    .e170 <- .e129 - .e10
    .e171 <- (2 * (.e4/(v4 * .e10)) + v5 * .e144 - ((.e62 * .e14 - 
        .e67/.e7)/.e160 + 2/.e17))/.e7
    .e172 <- .e23/.e13
    .e173 <- .e134 + .e37 * .e30 * .e108 * .e41 * .e14/.e159
    .e174 <- .e14/.e33
    .e176 <- .e135 * .e4/.e32
    .e178 <- v4 * .e9 * .e7
    .e181 <- .e36 * .e92 * .e4/.e63
    .e183 <- (((.e168 - .e129 * .e9) * .e4 + .e135 * (2 * (.e85 * 
        .e98 * .e4/.e32) - 1))/.e32 + (.e171 + 2 - .e173)/.e35)/v5 - 
        (.e9 * .e102/.e13 + v5^3/.e157) * .e4
    .e185 <- (((.e30 * (.e133 - (.e133 + .e174)) + .e42 + 2 - 
        (.e134 + .e54))/.e11 + v5 * (.e59/.e32 - 1/.e13) * .e4)/v5 - 
        (((.e25 - (.e29 + 2 * (.e55 * .e88 * .e98/.e32))) * .e9 + 
            .e10)/.e32 + .e83) * .e4)/v5 + .e9 * (2 - .e71) * 
        .e4/.e13
    .e187 <- (.e97 + .e147/.e13) * .e4 - .e152
    .e188 <- .e161 - .e167
    .e191 <- ((.e130 + 1/.e57)/.e11 - .e106) * .e4 + .e119/.e13 - 
        .e145
    .e193 <- ((.e101/.e136 + v5 * (.e24/.e22 + .e148) * .e9/.e32) * 
        .e4 - .e83)/v5 + (1 - .e71) * .e9/.e13
    .e195 <- .e166 - .e181
    .e202 <- ((.e82/(.e36 * v5 * .e7) + v5 * ((.e170/v4 + .e148) * 
        .e9/.e32 + .e31 * .e47/.e158))/v5 - 2 * (.e85 * .e7/.e156)) * 
        .e4 + .e9/.e13 - (.e130 + .e47/.e124)/v5
    .e205 <- (.e99 - 2/.e45) * .e4 + .e145 - .e172
    .e208 <- .e147 * .e4/.e13 - (.e76 - .e176)/.e11
    .e209 <- .e23 * .e4
    .e212 <- ta * .e104 * v5 * .e48
    .e214 <- ta * v5 * .e48
    .e217 <- .e103 * tb * v5 * .e48
    .e219 <- .e103 * v5 * .e48
    .e221 <- tb * v5 * .e48
    .e223 <- .e104 * v5 * .e48
    .e224 <- .e178 * .e4
    .e226 <- -(ta * .e114)
    .e228 <- -(.e72 * .e114)
    .e230 <- -(tb * .e114)
    .e231 <- (.e152 - .e209/.e13)/v4
    .e232 <- .e84 + 2 * (.e22 * .e88 * .e77 * .e12/.e32)
    .e234 <- .e44 * .e14/v5
    .e236 <- 1/.e37 + 2 * (.e224/.e13)
    .e237 <- 2 * .e11
    .e238 <- c(v1 = .e150, v2 = .e217, v3 = .e212, v4 = .e72 * 
        .e127/v4, v5 = -(.e72 * .e111))
    .e239 <- c(v1 = .e214, v2 = .e219, v3 = .e150, v4 = ta * 
        .e127/v4, v5 = -(ta * .e111))
    .e240 <- c(v1 = .e221, v2 = .e150, v3 = .e223, v4 = tb * 
        .e127/v4, v5 = -(tb * .e111))
    .e241 <- ta * .e109
    .e242 <- ta * .e113
    .e243 <- ta * .e116
    .e244 <- .e72 * .e109
    .e245 <- .e72 * .e113
    .e246 <- .e72 * .e116
    .e247 <- tb * .e109
    .e248 <- tb * .e113
    .e249 <- tb * .e116
    c(v1 = c(v1 = c(v1 = v5 * .e48, v2 = .e214, v3 = .e221, v4 = .e127/v4, 
        v5 = -.e111), v2 = .e239, v3 = .e240, v4 = c(v1 = -.e114, 
        v2 = .e226, v3 = .e230, v4 = -(.e205/v4), v5 = -.e202), 
        v5 = c(v1 = .e109, v2 = .e241, v3 = .e247, v4 = .e187/v4, 
            v5 = -.e183)), v2 = c(v1 = .e239, v2 = c(v1 = .e219, 
        v2 = ta^3 * v5 * .e48, v3 = .e217, v4 = .e103 * .e127/v4, 
        v5 = -(.e103 * .e111)), v3 = .e238, v4 = c(v1 = .e226, 
        v2 = -(.e103 * .e114), v3 = .e228, v4 = -(ta * .e205/v4), 
        v5 = -(ta * .e202)), v5 = c(v1 = .e241, v2 = .e103 * 
        .e109, v3 = .e244, v4 = ta * .e187/v4, v5 = -(ta * .e183))), 
        v3 = c(v1 = .e240, v2 = .e238, v3 = c(v1 = .e223, v2 = .e212, 
            v3 = tb^3 * v5 * .e48, v4 = .e104 * .e127/v4, v5 = -(.e104 * 
                .e111)), v4 = c(v1 = .e230, v2 = .e228, v3 = -(.e104 * 
            .e114), v4 = -(tb * .e205/v4), v5 = -(tb * .e202)), 
            v5 = c(v1 = .e247, v2 = .e244, v3 = .e104 * .e109, 
                v4 = tb * .e187/v4, v5 = -(tb * .e183))), v4 = c(v1 = c(v1 = .e113, 
            v2 = .e242, v3 = .e248, v4 = .e191/v4, v5 = -.e193), 
            v2 = c(v1 = .e242, v2 = .e103 * .e113, v3 = .e245, 
                v4 = ta * .e191/v4, v5 = -(ta * .e193)), v3 = c(v1 = .e248, 
                v2 = .e245, v3 = .e104 * .e113, v4 = tb * .e191/v4, 
                v5 = -(tb * .e193)), v4 = c(v1 = .e169, v2 = ta * 
                .e169, v3 = tb * .e169, v4 = (.e169 * .e4 + 2 * 
                (((.e209/.e11 - 1)/v4 + (.e172 + .e4/.e45) * 
                  .e4)/v4))/v4, v5 = -(((.e82 * .e4/(v4^4 * v5 * 
                .e7) + (v5 * .e170 * .e9 * .e4/.e22 - 2 * (v4 * 
                .e77^2 * .e12/.e32))/.e32 - v5 * (.e237 - .e5) * 
                .e47/.e158)/v5 + 2 * (.e178/.e156)) * .e4)), 
            v5 = c(v1 = .e188, v2 = ta * .e188, v3 = tb * .e188, 
                v4 = (.e161 - (.e231 + .e167)) * .e4/v4, v5 = -(((((.e168 + 
                  (.e10 - .e129) * .e9) * .e4/v4 - .e232)/.e32 + 
                  (.e171 + 1 - .e173)/.e124)/v5 + 2 * (.e224/.e156) + 
                  .e37/.e157) * .e4))), v5 = c(v1 = c(v1 = .e116, 
            v2 = .e243, v3 = .e249, v4 = .e208/v4, v5 = -.e185), 
            v2 = c(v1 = .e243, v2 = .e103 * .e116, v3 = .e246, 
                v4 = ta * .e208/v4, v5 = -(ta * .e185)), v3 = c(v1 = .e249, 
                v2 = .e246, v3 = .e104 * .e116, v4 = tb * .e208/v4, 
                v5 = -(tb * .e185)), v4 = c(v1 = -.e195, v2 = -(ta * 
                .e195), v3 = -(tb * .e195), v4 = -((.e231 + .e166 - 
                .e181) * .e4/v4), v5 = -(((((.e30 * (.e133 - 
                .e174) + .e42 + 1 - .e79)/.e58 - .e130)/v5 + 
                (((.e141 - .e25) * .e9 - .e10) * .e4/v4 - .e232)/.e32 - 
                v4 * .e47 * (.e11 + .e5)/.e158)/v5 + .e236/.e13) * 
                .e4)), v5 = c(v1 = .e165, v2 = ta * .e165, v3 = tb * 
                .e165, v4 = .e165 * .e4/v4, v5 = -((((((.e43 + 
                v5 * (.e108 * .e41 * .e4/v4 - 2 * (.e43 * .e14/.e37))) * 
                .e14/.e159 - 2 * (.e4/(.e31 * .e27))) * .e30 - 
                ((.e155 * .e14 + .e30 * .e4/.e11)/.e160 + .e144 * 
                  .e4/v4)/.e7)/.e7 - ((2 * ((.e44 * .e4/v4 - 
                .e132)/.e7 - .e234) + 2 * (.e144/.e7) - 4 * .e234)/v5 + 
                .e176))/v5 + v4 * (((.e30 * .e9 - .e10/v5) * 
                .e4/v4 - ((.e87 - (.e78/v5 + 2 * .e12)/v5) * 
                .e14 + .e10 * .e4/v4)/v5)/v5 - 2 * (.e22 * .e88^2 * 
                .e12/.e32)) * .e4/.e32)/v5 + (.e236 * .e4/.e13 + 
                v5 * (.e237 + .e5)/.e157) * .e4))))
}
############################################################
#' The first derivative of the density for DMGS
#' @returns Vector
#' @inheritParams manf
gev_p1b_f1fa=function(x,t0a,t0b,v1,v2,v3,v4,v5){

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_fd,"x")
	f1=vf(x,t0a,t0b,v1,v2,v3,v4,v5)
	return(f1)
}
############################################################
#' The first derivative of the density for WAIC
#' @returns Vector
#' @inheritParams manf
gev_p1b_f1fw=function(x,ta,tb,v1,v2,v3,v4,v5){

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_fd,c("x","ta","tb"))
	f1=vf(x,ta,tb,v1,v2,v3,v4,v5)
	return(f1)
}
############################################################
#' The second derivative of the density for DMGS
#' @returns Matrix
#' @inheritParams manf
gev_p1b_f2fa=function(x,t0a,t0b,v1,v2,v3,v4,v5){
	nx=length(x)

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_fdd,"x")
	temp1=vf(x,t0a,t0b,v1,v2,v3,v4,v5)
	f2=deriv_copyfdd(temp1,nx,dim=5)
	return(f2)
}
############################################################
#' The second derivative of the density for WAIC
#' @returns Matrix
#' @inheritParams manf
gev_p1b_f2fw=function(x,ta,tb,v1,v2,v3,v4,v5){
	nx=length(x)

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_fdd,c("x","ta","tb"))
	temp1=vf(x,ta,tb,v1,v2,v3,v4,v5)
	f2=deriv_copyfdd(temp1,nx,dim=5)
	return(f2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
gev_p1b_mu1fa=function(alpha,t0a,t0b,v1,v2,v3,v4,v5){
	x=qgev((1-alpha),mu=v1+v2*t0a+v3*t0b,sigma=v4,xi=v5)

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_pd,"x")
	mu1=-vf(x,t0a,t0b,v1,v2,v3,v4,v5)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
gev_p1b_mu2fa=function(alpha,t0a,t0b,v1,v2,v3,v4,v5){
	x=qgev((1-alpha),mu=v1+v2*t0a+v3*t0b,sigma=v4,xi=v5)
	nx=length(x)

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_pdd,"x")
	temp1=vf(x,t0a,t0b,v1,v2,v3,v4,v5)
	mu2=-deriv_copyfdd(temp1,nx,dim=5)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
gev_p1b_ldda=function(x,ta,tb,v1,v2,v3,v4,v5){
	nx=length(x)

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_logfdd,c("x","ta","tb"))
	temp1=vf(x,ta,tb,v1,v2,v3,v4,v5)
	ldd=deriv_copyldd(temp1,nx,dim=5)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
gev_p1b_lddda=function(x,ta,tb,v1,v2,v3,v4,v5){
	nx=length(x)

	v5=movexiawayfromzero(v5)

	vf=Vectorize(gev_p1b_logfddd,c("x","ta","tb"))
	temp1=vf(x,ta,tb,v1,v2,v3,v4,v5)
	lddd=deriv_copylddd(temp1,nx,dim=5)
	return(lddd)
}
