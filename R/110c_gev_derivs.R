######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_fd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- 1/v3
    .e4 <- v3 * .e1/v2
    .e5 <- 1 + .e4
    .e6 <- 1 + .e2
    .e7 <- .e5^.e6
    .e9 <- .e5^(.e2 + 2)
    .e10 <- exp(-.e5^-.e2)
    .e13 <- log1p(.e4)
    .e14 <- v2^2
    .e17 <- v3 * .e6/.e9 - 1/.e5^(2 * .e6)
    c(v1 = .e10 * .e17/.e14, v2 = (.e17 * .e1/v2 - 1/.e7) * .e10/.e14, 
        v3 = ((.e13/(v3 * .e7) - (.e13/(v3 * .e5^.e2) - .e1/(v2 * 
            .e7))/.e7)/v3 - .e6 * .e1/(v2 * .e9)) * .e10/v2)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_fdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- 1/v3
    .e4 <- v3 * .e1/v2
    .e5 <- 1 + .e4
    .e6 <- 1 + .e2
    .e7 <- .e5^.e6
    .e8 <- .e2 + 2
    .e9 <- log1p(.e4)
    .e10 <- .e5^.e8
    .e11 <- 2 * .e6
    .e12 <- .e5^.e2
    .e13 <- v3 * .e12
    .e14 <- .e5^.e11
    .e15 <- v2 * .e7
    .e16 <- v3 * .e6
    .e18 <- 1/.e14
    .e20 <- .e9/.e13 - .e1/.e15
    .e21 <- v3 * .e7
    .e22 <- v3^2
    .e24 <- exp(-.e5^-.e2)
    .e25 <- v2 * .e10
    .e27 <- .e16/.e10 - .e18
    .e28 <- 1/.e7
    .e29 <- 1/.e10
    .e30 <- 2 * .e8
    .e31 <- .e20/.e7
    .e32 <- .e9/.e21
    .e33 <- v2^2
    .e37 <- .e6 * .e12 * .e1/v2 - .e7 * .e9/.e22
    .e39 <- .e6 * .e1/.e25
    .e42 <- .e5^(.e2 - .e11)
    .e43 <- (.e32 - .e31)/v3
    .e46 <- .e27 * .e1/v2 - .e28
    .e48 <- v2^3
    .e51 <- v3 * .e5^(.e6 - .e30) * .e8 - 2/.e5^(1 + .e11)
    .e55 <- .e7 * .e8 * .e1/v2 - .e10 * .e9/.e22
    .e56 <- .e5^(.e2 - 1)
    .e57 <- .e43 - .e39
    .e58 <- .e15^2
    .e59 <- .e25^2
    .e60 <- .e21^2
    .e61 <- .e13^2
    .e62 <- .e37/.e14
    .e74 <- (2 * (.e6 * .e5^(.e11 - 1) * .e1/v2) - 2 * (.e14 * 
        .e9/.e22))/.e5^(4 * .e6) + .e29
    .e80 <- .e18 + v3 * (.e51 * .e1/v2 - (.e42 + .e29)) * .e6 - 
        .e46/.e7
    .e83 <- v3 * .e55 * .e6/.e5^.e30
    .e85 <- .e16 * .e12 * .e1
    .e87 <- .e16 * .e51 - .e27/.e7
    .e89 <- .e16 * (.e13 * .e9/.e60 - .e42 * .e20) - .e29
    .e91 <- .e21 * .e8 * .e1
    .e94 <- v3 * .e56 * .e9/.e61
    c(v1 = c(v1 = .e24 * .e87/.e48, v2 = .e80 * .e24/.e48, v3 = ((.e89/v2 - 
        ((.e28 + .e94 - .e28)/v2 - .e85/.e58)/.e7)/v3 - (.e57/.e15 + 
        .e6 * (.e91/.e59 - 1/.e25))) * .e24/v2), v2 = c(v1 = (.e87 * 
        .e1/v2 - 2 * .e27) * .e24/.e48, v2 = (.e80 * .e1/v2 - 
        2 * .e46) * .e24/.e48, v3 = (((.e10 - .e91/v2) * .e6/.e59 + 
        (.e89/.e33 - ((.e7 - .e85/v2)/.e58 + (.e94 - .e28)/.e33)/.e7)/v3 - 
        .e57/(.e33 * .e7)) * .e1 - .e57/v2) * .e24/v2), v3 = c(v1 = (.e74 - 
        (.e20 * .e27/v3 + .e83)) * .e24/.e33, v2 = (.e62 + (.e74 - 
        .e83) * .e1/v2 - .e46 * .e20/v3) * .e24/.e33, v3 = (((.e62 + 
        .e39 - .e43) * .e20 + (.e31 + .e1/.e25 - .e32)/v3 - ((.e7 + 
        v3 * .e37) * .e9/.e60 + ((1/(v2 * v3 * .e7) + v2 * .e37/.e58) * 
        .e1 - (.e56 * .e1/v2 + .e12 - .e12 * .e9/v3) * .e9/.e61)/.e7))/v3 + 
        (1/(v2 * .e22 * .e10) + v2 * .e55 * .e6/.e59) * .e1) * 
        .e24/v2))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_pd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e3 <- v3 * .e1/v2
    .e4 <- 1 + .e3
    .e5 <- 1/v3
    .e7 <- .e4^(1 + .e5)
    .e8 <- exp(-.e4^-.e5)
    .e9 <- v2 * .e7
    c(v1 = -(.e8/.e9), v2 = -(.e8 * .e1/(v2^2 * .e7)), v3 = -(.e8 * 
        (log1p(.e3)/(v3 * .e4^.e5) - .e1/.e9)/v3))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_pdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e3 <- v3 * .e1/v2
    .e4 <- 1/v3
    .e5 <- 1 + .e3
    .e6 <- 1 + .e4
    .e7 <- .e5^.e6
    .e8 <- .e5^.e4
    .e9 <- v2 * .e7
    .e10 <- log1p(.e3)
    .e12 <- exp(-.e5^-.e4)
    .e13 <- v3 * .e8
    .e14 <- v2^2
    .e15 <- .e1/.e9
    .e16 <- .e9^2
    .e17 <- .e10/.e13
    .e18 <- .e17 - .e15
    .e20 <- v3 * .e6 * .e8
    .e21 <- .e14 * .e7
    .e22 <- .e20 * .e1
    .e26 <- .e6 * .e8 * .e1/v2 - .e7 * .e10/v3^2
    .e27 <- .e5^(.e4 - 1)
    .e28 <- .e5^(2 * .e6)
    .e29 <- .e21^2
    .e30 <- .e13^2
    .e31 <- 1/.e7
    .e32 <- v2 * v3
    .e33 <- (.e7 - .e22/v2)/.e16
    .e35 <- .e18/.e7 + .e31
    .e37 <- v2 * .e26/.e16
    .e38 <- .e32 * .e7
    .e41 <- v3 * .e27 * .e10/.e30
    c(v1 = c(v1 = -(.e12 * (.e20/.e16 - 1/(.e14 * .e28))), v2 = -(.e12 * 
        (.e32 * .e6 * .e8 * .e1/.e29 - (.e15 + 1)/.e21)), v3 = -(((.e31 + 
        .e41 - .e35)/v2 - .e22/.e16) * .e12/v3)), v2 = c(v1 = (.e33 + 
        .e1/(v2^3 * .e28)) * .e12, v2 = ((2 * .e9 - .e22)/.e29 + 
        .e1/(v2^4 * .e28)) * .e12 * .e1, v3 = -((.e33 + (.e41 - 
        .e35)/.e14) * .e12 * .e1/v3)), v3 = c(v1 = (.e18/.e38 + 
        .e37) * .e12, v2 = (.e18/(.e14 * v3 * .e7) + .e14 * .e26/.e29) * 
        .e12 * .e1, v3 = -(((1/.e38 + .e37) * .e1 - ((.e27 * 
        .e1/v2 + .e8 - .e8 * .e10/v3) * .e10/.e30 + (1 + .e17 - 
        .e15) * .e18/v3)) * .e12/v3)))
}
############################################################
#' First derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_logfd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e3 <- v3 * .e1/v2
    .e4 <- 1 + .e3
    .e5 <- 1/v3
    .e6 <- 1 + .e5
    .e8 <- 1/.e4^.e5
    .e9 <- v2 * .e4
    .e11 <- v3 * .e6 - .e8
    c(v1 = .e11/.e9, v2 = (.e11 * .e1/.e9 - 1)/v2, v3 = -(((.e8 - 
        1) * log1p(.e3)/v3 - .e1/(v2 * .e4^.e6))/v3 + .e6 * .e1/.e9))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_logfdd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- v3 * .e1
    .e3 <- .e2/v2
    .e4 <- 1 + .e3
    .e5 <- 1/v3
    .e6 <- .e4^.e5
    .e7 <- 1 + .e5
    .e8 <- v2 * .e4
    .e9 <- 1/.e6
    .e10 <- v3 * .e7
    .e11 <- .e8^2
    .e12 <- log1p(.e3)
    .e13 <- .e10 - .e9
    .e14 <- .e4^.e7
    .e15 <- v2 * .e14
    .e19 <- .e4^(.e5 - 1) * .e1/v2 - .e6 * .e12/v3
    .e20 <- .e4^(.e5 + 2)
    .e21 <- .e15^2
    .e22 <- .e13 * .e1
    .e23 <- .e9 - 1
    .e24 <- 2/v3
    .e26 <- (.e19/(v3 * .e4^.e24) + 1)/.e8 - .e22/.e11
    .e28 <- .e13/.e11 + .e1/(v2^3 * .e20)
    .e30 <- .e12/.e6 - v3 * .e23
    .e31 <- v2^2
    .e33 <- .e10 * .e6 * .e1
    .e34 <- v3 * .e13
    .e35 <- v3^2
    c(v1 = c(v1 = .e34/.e11 - 1/(.e31 * .e20), v2 = (.e34 * .e1/.e11 - 
        (.e1/.e15 + .e10 - .e9)/.e8)/v2, v3 = -(((.e30/v3 + .e9)/.e8 - 
        .e33/.e21)/v3 + .e7 * (.e2/.e11 - 1/.e8))), v2 = c(v1 = -.e28, 
        v2 = -(((.e22/.e8 - 1)/v2 + .e28 * .e1)/v2), v3 = -((((.e14 - 
            .e33/v2)/.e21 + .e30/(.e31 * v3 * .e4))/v3 - .e7/.e11) * 
            .e1)), v3 = c(v1 = .e26, v2 = .e26 * .e1/v2, v3 = -(((((2/.e6 - 
        1) * .e1/v2 - .e19 * .e12/(v3 * .e4^(.e24 - 1)))/.e4 - 
        2 * (.e23 * .e12/v3))/v3 + v2 * (.e7 * .e6 * .e1/v2 - 
        .e14 * .e12/.e35) * .e1/.e21)/v3 - (.e7 * .e1/.e11 + 
        1/(v2 * .e35 * .e4)) * .e1)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
gev_logfddd=function (x, v1, v2, v3) 
{
    .e1 <- x - v1
    .e2 <- v3 * .e1
    .e3 <- .e2/v2
    .e4 <- 1 + .e3
    .e5 <- 1/v3
    .e6 <- 1 + .e5
    .e7 <- .e4^.e5
    .e8 <- log1p(.e3)
    .e9 <- v2 * .e4
    .e10 <- .e4^.e6
    .e11 <- .e5 - 1
    .e12 <- .e4^.e11
    .e13 <- 2/v3
    .e14 <- .e9^2
    .e15 <- 1/.e7
    .e16 <- .e12 * .e1
    .e17 <- .e16/v2
    .e19 <- .e7 * .e8/v3
    .e20 <- v2 * .e10
    .e21 <- v3 * .e6
    .e22 <- .e4^.e13
    .e23 <- .e17 - .e19
    .e24 <- v2^2
    .e25 <- .e5 + 2
    .e26 <- .e20^2
    .e27 <- .e21 - .e15
    .e28 <- v3 * .e22
    .e29 <- v3^2
    .e30 <- .e13 - 1
    .e31 <- .e4^.e25
    .e32 <- .e4^.e30
    .e33 <- .e9 * .e27
    .e34 <- v2^3
    .e35 <- .e15 - 1
    .e36 <- v2 * v3
    .e37 <- .e6 * .e7
    .e39 <- .e8/.e7 - v3 * .e35
    .e40 <- .e23/.e28
    .e41 <- .e10 * .e8
    .e42 <- .e41/.e29
    .e43 <- .e12 * .e8
    .e44 <- .e4^(.e5 - 2)
    .e45 <- .e1/.e20
    .e46 <- 2 * (.e33 * .e1/.e14)
    .e47 <- .e34 * .e31
    .e49 <- .e21 * .e7 * .e1
    .e50 <- .e40 + 1
    .e52 <- .e37 * .e1/v2
    .e53 <- .e10 - .e49/v2
    .e54 <- 2/.e7
    .e55 <- .e36 * .e4
    .e56 <- .e24 * v3
    .e57 <- v3 * .e32
    .e58 <- .e52 - .e42
    .e59 <- .e39/v3
    .e60 <- .e24 * .e10
    .e62 <- (.e43 + v3 * .e12)/v3 - (.e12 + v3 * .e44 * .e11 * 
        .e1/v2)
    .e63 <- .e4^(.e13 - 2)
    .e64 <- .e47^2
    .e65 <- .e28^2
    .e67 <- .e1/.e60 + 2 * (.e33/.e14)
    .e68 <- v3 * .e23
    .e70 <- .e23 * .e8/.e57
    .e71 <- .e56 * .e4
    .e73 <- .e62/.e28 + 2 * (.e68 * .e32/.e65)
    .e74 <- .e23/.e22
    .e75 <- .e67/.e14
    .e76 <- .e4^(1 + .e13)
    .e77 <- .e17 + .e7
    .e78 <- .e7 + .e19
    .e81 <- (.e54 - 1) * .e1/v2 - .e70
    .e82 <- .e59 + .e15
    .e83 <- 2 * (.e55 * .e1/.e14)
    .e86 <- .e8/.e10 - 2 * (v3/.e10)
    .e88 <- v2 * .e6 * .e4
    .e89 <- .e24 * .e4
    .e91 <- .e50/.e9
    .e92 <- .e53/.e26
    .e95 <- (.e44 * .e11 * .e1/v2 - .e43/.e29) * .e1/v2 - ((.e17 - 
        .e78) * .e8/v3 + .e17)/v3
    .e96 <- .e82/.e14
    .e100 <- .e10 * .e25 * .e1/v2 - .e31 * .e8/.e29
    .e101 <- (2 * (.e55 * .e27/.e14) - 1/.e20)/.e14
    .e102 <- .e9^4
    .e103 <- (v2 * .e29 * .e4)^2
    .e104 <- (.e24 * .e31)^2
    .e105 <- .e71^2
    .e106 <- .e27 * .e1
    .e107 <- .e57^2
    .e109 <- .e45 + .e21 - .e15
    .e110 <- 1/.e22
    .e111 <- 2 * .e59
    .e112 <- v2 * .e58
    .e113 <- .e88 * .e1
    .e114 <- v2 * .e31
    .e115 <- .e36 * .e6
    .e116 <- .e34 * .e4
    .e119 <- v3 * .e10 * .e25 * .e1
    .e120 <- v3 * .e63
    .e121 <- (.e73/.e116 + .e75) * .e1
    .e122 <- (.e91 - .e106/.e14)/v2
    .e125 <- .e73/.e89
    .e127 <- (.e95/.e28 - .e23 * (.e22 + 2 * (.e32 * .e1/v2) - 
        2 * (.e22 * .e8/v3))/.e65)/.e9 - (2 + 2 * .e40 - .e46) * 
        .e1/.e14
    .e134 <- (.e50 - .e46)/.e14
    .e135 <- .e50/.e14
    .e136 <- (.e37 * .e8 + .e7)/v3
    .e138 <- .e75 + v2 * (3 * .e114 - .e119) * .e1/.e64
    .e139 <- .e42 + 2 * (.e24 * .e58 * .e53 * .e10/.e26)
    .e140 <- .e77 - .e7
    .e142 <- .e35 * .e8/v3
    .e143 <- (2 * (.e1/(v2 * .e7)) + v3 * .e81 - ((.e62 * .e8 - 
        .e68/.e4)/.e120 + 2/.e12))/.e4
    .e144 <- .e101 + .e56 * .e10 * .e25 * .e1/.e64
    .e146 <- .e86/v3 + 1/.e10
    .e147 <- .e27/.e14
    .e148 <- (v3 * .e27 * .e1/.e14 - .e109/.e9)/v2
    .e150 <- .e45 + v3 * ((.e74 + 2)/v3 + 3 - .e46) - .e54
    .e152 <- 1/.e29 + 2 * (.e113/.e14)
    .e153 <- .e111 + .e29 * .e23 * .e63 * .e30 * .e8/.e107
    .e154 <- 2 * (v2 * .e53 * .e76/.e26)
    .e155 <- 2 * .e9
    .e156 <- .e83 - 2
    .e157 <- .e8/.e28
    .e159 <- .e112 * .e1/.e26
    .e162 <- .e34 * .e100 * .e1/.e64
    .e163 <- v3 * ((.e74 + 1)/v3 + 2 - .e46)
    .e165 <- v3 * .e67/.e14
    c(v1 = c(v1 = c(v1 = v3 * (.e101 - .e20 * .e25/.e104), v2 = (v3 * 
        (.e54 + v3 * (.e46 - (2 + .e13)) - 2 * .e45)/.e14 - (.e49/.e26 - 
        2/.e20)/.e9)/v2, v3 = -((.e146/.e89 + v3 * (.e96 - ((2 * 
        (.e115 * .e4^(1 + 3/v3)/.e26) - .e12/v2) * .e1 - .e7) * 
        .e6/.e26))/v3 + .e21 * .e156/.e14)), v2 = c(v1 = -(.e144 - 
        1/.e47), v2 = -(((.e144 - 2/.e47) * .e1 + .e148 - .e147)/v2), 
        v3 = -(((.e86/(.e34 * v3 * .e4) + v3 * ((.e140/v2 + .e154) * 
            .e6/.e26 + .e36 * .e39/.e105))/v3 - 2 * (.e115 * 
            .e4/.e102)) * .e1 + .e6/.e14 - (.e92 + .e39/.e71)/v3)), 
        v3 = c(v1 = .e125 + (.e45 + .e163 - .e15)/.e14, v2 = ((.e125 + 
            .e150/.e14) * .e1 - .e91)/v2, v3 = -((((.e136 - .e77 * 
            .e6) * .e1 + .e112 * (2 * (.e115 * .e76 * .e1/.e26) - 
            1))/.e26 + (.e143 + 2 - .e153)/.e55)/v3 - (.e6 * 
            .e156/.e14 + v3^3/.e103) * .e1))), v2 = c(v1 = c(v1 = (2 * 
        .e114 - .e119)/.e104 - .e165, v2 = (((.e92 + 1/.e60)/.e9 - 
        .e165) * .e1 + .e109/.e14 - .e148)/v2, v3 = -(((.e146/.e116 + 
        v3 * (.e16/.e24 + .e154) * .e6/.e26) * .e1 - .e96)/v3 + 
        (1 - .e83) * .e6/.e14)), v2 = c(v1 = .e138, v2 = (.e138 * 
        .e1 + 2 * (((.e106/.e9 - 1)/v2 + (.e147 + .e1/.e47) * 
        .e1)/v2))/v2, v3 = -(((.e86 * .e1/(v2^4 * v3 * .e4) + 
        (v3 * .e140 * .e6 * .e1/.e24 - 2 * (v2 * .e53^2 * .e10/.e26))/.e26 - 
        v3 * (.e155 - .e2) * .e39/.e105)/v3 + 2 * (.e88/.e102)) * 
        .e1)), v3 = c(v1 = .e121 - .e135, v2 = (.e121 - (.e122 + 
        .e135)) * .e1/v2, v3 = -(((((.e136 + (.e7 - .e77) * .e6) * 
        .e1/v2 - .e139)/.e26 + (.e143 + 1 - .e153)/.e71)/v3 + 
        2 * (.e113/.e102) + .e29/.e103) * .e1))), v3 = c(v1 = c(v1 = (.e163 - 
        .e15)/.e14 + .e24 * .e100/.e104, v2 = (.e150 * .e1/.e14 - 
        (.e50 - .e159)/.e9)/v2, v3 = -((((.e23 * (.e110 - (.e110 + 
        .e157)) + .e45 + 2 - (.e111 + .e54))/.e9 + v3 * (.e37/.e26 - 
        1/.e14) * .e1)/v3 - (((.e17 - (.e19 + 2 * (.e56 * .e58 * 
        .e76/.e26))) * .e6 + .e7)/.e26 + .e96) * .e1)/v3 + .e6 * 
        (2 - .e83) * .e1/.e14)), v2 = c(v1 = -(.e134 - .e162), 
        v2 = -((.e122 + .e134 - .e162) * .e1/v2), v3 = -(((((.e23 * 
            (.e110 - .e157) + .e45 + 1 - .e82)/.e89 - .e92)/v3 + 
            (((.e78 - .e17) * .e6 - .e7) * .e1/v2 - .e139)/.e26 - 
            v2 * .e39 * (.e9 + .e2)/.e105)/v3 + .e152/.e14) * 
            .e1)), v3 = c(v1 = .e127, v2 = .e127 * .e1/v2, v3 = -((((((.e32 + 
        v3 * (.e63 * .e30 * .e1/v2 - 2 * (.e32 * .e8/.e29))) * 
        .e8/.e107 - 2 * (.e1/(.e36 * .e22))) * .e23 - ((.e95 * 
        .e8 + .e23 * .e1/.e9)/.e120 + .e81 * .e1/v2)/.e4)/.e4 - 
        ((2 * ((.e35 * .e1/v2 - .e70)/.e4 - .e142) + 2 * (.e81/.e4) - 
            4 * .e142)/v3 + .e159))/v3 + v2 * (((.e23 * .e6 - 
        .e7/v3) * .e1/v2 - ((.e52 - (.e41/v3 + 2 * .e10)/v3) * 
        .e8 + .e7 * .e1/v2)/v3)/v3 - 2 * (.e24 * .e58^2 * .e10/.e26)) * 
        .e1/.e26)/v3 + (.e152 * .e1/.e14 + v3 * (.e155 + .e2)/.e103) * 
        .e1))))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
gev_f1fa=function(x,v1,v2,v3){

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_fd,"x")
	f1=vf(x,v1,v2,v3)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
gev_f2fa=function(x,v1,v2,v3){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_fdd,"x")
	temp1=vf(x,v1,v2,v3)
	f2=deriv_copyfdd(temp1,nx,dim=3)
	return(f2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
gev_mu1fa=function(alpha,v1,v2,v3){
	x=qgev((1-alpha),mu=v1,sigma=v2,xi=v3)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_pd,"x")
	mu1=-vf(x,v1,v2,v3)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
gev_mu2fa=function(alpha,v1,v2,v3){
	x=qgev((1-alpha),mu=v1,sigma=v2,xi=v3)
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_pdd,"x")
	temp1=vf(x,v1,v2,v3)
	mu2=-deriv_copyfdd(temp1,nx,dim=3)
	return(mu2)
}
############################################################
#' The first derivative of the normalized log-likelihood
#' @returns Vector
#' @inheritParams manf
gev_lda=function(x,v1,v2,v3){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_logfd,"x")
	ld=vf(x,v1,v2,v3)
	return(ld)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
gev_ldda=function(x,v1,v2,v3){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_logfdd,"x")
	temp1=vf(x,v1,v2,v3)
	ldd=deriv_copyldd(temp1,nx,dim=3)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
gev_lddda=function(x,v1,v2,v3){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_logfddd,"x")
	temp1=vf(x,v1,v2,v3)
	lddd=deriv_copylddd(temp1,nx,dim=3)
	return(lddd)
}
############################################################
#' The combined derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
gev_ld12a=function(x,v1,v2,v3){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf1=Vectorize(gev_logfd,"x")
	vf2=Vectorize(gev_logfdd,"x")
	ld1=vf1(x,v1,v2,v3) 							#a matrix which is 3 by nx
	temp2=vf2(x,v1,v2,v3) 						#a matrix which is 9 by nx
	ld2=deriv_copyld2(temp2,nx,dim=3) #a matrix which is 3 by 3 by nx
	dim=3
	ld12=array(0,c(dim,dim,dim))
	for (j in 1:dim){
		for (s in 1:dim){
			for (r in 1:dim){
				for (i in 1:nx){
					ld12[j,s,r]=ld12[j,s,r]+(ld2[j,s,i]*ld1[r,i]/nx)
				}
			}
		}
	}
	return(ld12)
}
