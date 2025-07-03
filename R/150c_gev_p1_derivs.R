######################################################################
#' First derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_p1_fd=function (x, t, v1, v2, v3, v4) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- 1/v4
    .e6 <- v4 * .e3/v3
    .e7 <- 1 + .e6
    .e8 <- 1 + .e4
    .e10 <- .e7^.e8
    .e11 <- .e7^(.e4 + 2)
    .e12 <- exp(-.e7^-.e4)
    .e15 <- v3^2
    .e18 <- v4 * .e8/.e11 - 1/.e7^(2 * .e8)
    .e19 <- log1p(.e6)
    c(v1 = .e12 * .e18/.e15, v2 = t * .e12 * .e18/.e15, v3 = (.e18 * 
        .e3/v3 - 1/.e10) * .e12/.e15, v4 = ((.e19/(v4 * .e10) - 
        (.e19/(v4 * .e7^.e4) - .e3/(v3 * .e10))/.e10)/v4 - .e8 * 
        .e3/(v3 * .e11)) * .e12/v3)
}
######################################################################
#' Second derivative of the density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p1_fdd=function (x, t, v1, v2, v3, v4) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- 1/v4
    .e6 <- v4 * .e3/v3
    .e7 <- 1 + .e6
    .e8 <- 1 + .e4
    .e9 <- .e7^.e8
    .e10 <- .e4 + 2
    .e11 <- .e7^.e10
    .e12 <- 2 * .e8
    .e13 <- log1p(.e6)
    .e14 <- .e7^.e4
    .e15 <- v4 * .e8
    .e16 <- .e7^.e12
    .e17 <- v4 * .e14
    .e18 <- v3 * .e9
    .e19 <- 1/.e16
    .e21 <- exp(-.e7^-.e4)
    .e23 <- .e15/.e11 - .e19
    .e26 <- .e13/.e17 - .e3/.e18
    .e27 <- v4 * .e9
    .e28 <- 2 * .e10
    .e29 <- v3 * .e11
    .e30 <- v4^2
    .e31 <- 1/.e9
    .e34 <- 1/.e11
    .e36 <- v3^3
    .e39 <- v4 * .e7^(.e8 - .e28) * .e10 - 2/.e7^(1 + .e12)
    .e40 <- .e7^(.e4 - .e12)
    .e41 <- .e26/.e9
    .e43 <- .e13/.e27
    .e44 <- v3^2
    .e46 <- .e15 * .e39 - .e23/.e9
    .e48 <- .e8 * .e3/.e29
    .e49 <- (.e43 - .e41)/v4
    .e52 <- .e23 * .e3/v3 - .e31
    .e56 <- .e8 * .e14 * .e3/v3 - .e9 * .e13/.e30
    .e60 <- .e9 * .e10 * .e3/v3 - .e11 * .e13/.e30
    .e61 <- .e7^(.e4 - 1)
    .e62 <- .e49 - .e48
    .e63 <- .e18^2
    .e64 <- .e29^2
    .e65 <- .e27^2
    .e66 <- .e17^2
    .e78 <- (2 * (.e8 * .e7^(.e12 - 1) * .e3/v3) - 2 * (.e16 * 
        .e13/.e30))/.e7^(4 * .e8) + .e34
    .e83 <- .e19 + v4 * (.e39 * .e3/v3 - (.e40 + .e34)) * .e8 - 
        .e52/.e9
    .e86 <- v4 * .e60 * .e8/.e7^.e28
    .e88 <- .e15 * .e14 * .e3
    .e90 <- .e15 * (.e17 * .e13/.e65 - .e40 * .e26) - .e34
    .e92 <- .e27 * .e10 * .e3
    .e95 <- v4 * .e61 * .e13/.e66
    .e96 <- .e56/.e16
    .e101 <- (.e90/v3 - ((.e31 + .e95 - .e31)/v3 - .e88/.e63)/.e9)/v4 - 
        (.e62/.e18 + .e8 * (.e92/.e64 - 1/.e29))
    .e102 <- .e78 - (.e26 * .e23/v4 + .e86)
    .e105 <- .e46 * .e3/v3 - 2 * .e23
    .e108 <- t * .e21 * .e46/.e36
    c(v1 = c(v1 = .e21 * .e46/.e36, v2 = .e108, v3 = .e83 * .e21/.e36, 
        v4 = .e101 * .e21/v3), v2 = c(v1 = .e108, v2 = t^2 * 
        .e21 * .e46/.e36, v3 = t * .e83 * .e21/.e36, v4 = t * 
        .e101 * .e21/v3), v3 = c(v1 = .e105 * .e21/.e36, v2 = t * 
        .e105 * .e21/.e36, v3 = (.e83 * .e3/v3 - 2 * .e52) * 
        .e21/.e36, v4 = (((.e11 - .e92/v3) * .e8/.e64 + (.e90/.e44 - 
        ((.e9 - .e88/v3)/.e63 + (.e95 - .e31)/.e44)/.e9)/v4 - 
        .e62/(.e44 * .e9)) * .e3 - .e62/v3) * .e21/v3), v4 = c(v1 = .e102 * 
        .e21/.e44, v2 = t * .e102 * .e21/.e44, v3 = (.e96 + (.e78 - 
        .e86) * .e3/v3 - .e52 * .e26/v4) * .e21/.e44, v4 = (((.e96 + 
        .e48 - .e49) * .e26 + (.e41 + .e3/.e29 - .e43)/v4 - ((.e9 + 
        v4 * .e56) * .e13/.e65 + ((1/(v3 * v4 * .e9) + v3 * .e56/.e63) * 
        .e3 - (.e61 * .e3/v3 + .e14 - .e14 * .e13/v4) * .e13/.e66)/.e9))/v4 + 
        (1/(v3 * .e30 * .e11) + v3 * .e60 * .e8/.e64) * .e3) * 
        .e21/v3))
}
######################################################################
#' First derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Vector
#' @inheritParams manf
gev_p1_pd=function (x, t, v1, v2, v3, v4) 
{
    .e3 <- x - (t * v2 + v1)
    .e5 <- v4 * .e3/v3
    .e6 <- 1 + .e5
    .e7 <- 1/v4
    .e9 <- .e6^(1 + .e7)
    .e10 <- exp(-.e6^-.e7)
    .e11 <- v3 * .e9
    c(v1 = -(.e10/.e11), v2 = -(t * .e10/.e11), v3 = -(.e10 * 
        .e3/(v3^2 * .e9)), v4 = -(.e10 * (log1p(.e5)/(v4 * .e6^.e7) - 
        .e3/.e11)/v4))
}
######################################################################
#' Second derivative of the cdf
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p1_pdd=function (x, t, v1, v2, v3, v4) 
{
    .e3 <- x - (t * v2 + v1)
    .e5 <- v4 * .e3/v3
    .e6 <- 1/v4
    .e7 <- 1 + .e5
    .e8 <- 1 + .e6
    .e9 <- .e7^.e8
    .e10 <- .e7^.e6
    .e11 <- v3 * .e9
    .e12 <- log1p(.e5)
    .e14 <- exp(-.e7^-.e6)
    .e15 <- v3^2
    .e16 <- .e11^2
    .e17 <- v4 * .e10
    .e18 <- .e3/.e11
    .e20 <- v4 * .e8 * .e10
    .e21 <- .e12/.e17
    .e22 <- .e7^(2 * .e8)
    .e23 <- .e21 - .e18
    .e24 <- .e15 * .e9
    .e25 <- .e20 * .e3
    .e26 <- 1/.e9
    .e27 <- v3 * v4
    .e31 <- .e8 * .e10 * .e3/v3 - .e9 * .e12/v4^2
    .e32 <- .e7^(.e6 - 1)
    .e33 <- .e24^2
    .e34 <- .e17^2
    .e37 <- .e20/.e16 - 1/(.e15 * .e22)
    .e38 <- (.e9 - .e25/v3)/.e16
    .e40 <- .e23/.e9 + .e26
    .e41 <- t * .e14
    .e43 <- v3 * .e31/.e16
    .e44 <- .e27 * .e9
    .e47 <- v4 * .e32 * .e12/.e34
    .e48 <- -(.e41 * .e37)
    .e49 <- .e38 + .e3/(v3^3 * .e22)
    .e52 <- (.e26 + .e47 - .e40)/v3 - .e25/.e16
    .e54 <- .e23/.e44 + .e43
    .e59 <- .e27 * .e8 * .e10 * .e3/.e33 - (.e18 + 1)/.e24
    c(v1 = c(v1 = -(.e14 * .e37), v2 = .e48, v3 = -(.e14 * .e59), 
        v4 = -(.e52 * .e14/v4)), v2 = c(v1 = .e48, v2 = -(t^2 * 
        .e14 * .e37), v3 = -(.e41 * .e59), v4 = -(t * .e52 * 
        .e14/v4)), v3 = c(v1 = .e49 * .e14, v2 = t * .e49 * .e14, 
        v3 = ((2 * .e11 - .e25)/.e33 + .e3/(v3^4 * .e22)) * .e14 * 
            .e3, v4 = -((.e38 + (.e47 - .e40)/.e15) * .e14 * 
            .e3/v4)), v4 = c(v1 = .e54 * .e14, v2 = t * .e54 * 
        .e14, v3 = (.e23/(.e15 * v4 * .e9) + .e15 * .e31/.e33) * 
        .e14 * .e3, v4 = -(((1/.e44 + .e43) * .e3 - ((.e32 * 
        .e3/v3 + .e10 - .e10 * .e12/v4) * .e12/.e34 + (1 + .e21 - 
        .e18) * .e23/v4)) * .e14/v4)))
}
############################################################
#' Second derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns Matrix
#' @inheritParams manf
gev_p1_logfdd=function (x, t, v1, v2, v3, v4) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- v4 * .e3
    .e5 <- .e4/v3
    .e6 <- 1 + .e5
    .e7 <- 1/v4
    .e8 <- .e6^.e7
    .e9 <- 1 + .e7
    .e10 <- v3 * .e6
    .e11 <- 1/.e8
    .e12 <- v4 * .e9
    .e13 <- .e10^2
    .e14 <- .e12 - .e11
    .e15 <- log1p(.e5)
    .e16 <- .e6^.e9
    .e17 <- .e6^(.e7 + 2)
    .e18 <- v3 * .e16
    .e19 <- v4 * .e14
    .e20 <- v3^2
    .e24 <- .e6^(.e7 - 1) * .e3/v3 - .e8 * .e15/v4
    .e25 <- .e18^2
    .e26 <- .e14 * .e3
    .e27 <- .e11 - 1
    .e29 <- 2/v4
    .e31 <- .e19/.e13 - 1/(.e20 * .e17)
    .e33 <- (.e24/(v4 * .e6^.e29) + 1)/.e10 - .e26/.e13
    .e35 <- .e14/.e13 + .e3/(v3^3 * .e17)
    .e37 <- .e15/.e8 - v4 * .e27
    .e39 <- .e12 * .e8 * .e3
    .e41 <- ((.e37/v4 + .e11)/.e10 - .e39/.e25)/v4 + .e9 * (.e4/.e13 - 
        1/.e10)
    .e43 <- t * .e31
    .e46 <- .e19 * .e3/.e13 - (.e3/.e18 + .e12 - .e11)/.e10
    .e47 <- v4^2
    c(v1 = c(v1 = .e31, v2 = .e43, v3 = .e46/v3, v4 = -.e41), 
        v2 = c(v1 = .e43, v2 = t^2 * .e31, v3 = t * .e46/v3, 
            v4 = -(t * .e41)), v3 = c(v1 = -.e35, v2 = -(t * 
            .e35), v3 = -(((.e26/.e10 - 1)/v3 + .e35 * .e3)/v3), 
            v4 = -((((.e16 - .e39/v3)/.e25 + .e37/(.e20 * v4 * 
                .e6))/v4 - .e9/.e13) * .e3)), v4 = c(v1 = .e33, 
            v2 = t * .e33, v3 = .e33 * .e3/v3, v4 = -(((((2/.e8 - 
                1) * .e3/v3 - .e24 * .e15/(v4 * .e6^(.e29 - 1)))/.e6 - 
                2 * (.e27 * .e15/v4))/v4 + v3 * (.e9 * .e8 * 
                .e3/v3 - .e16 * .e15/.e47) * .e3/.e25)/v4 - (.e9 * 
                .e3/.e13 + 1/(v3 * .e47 * .e6)) * .e3)))
}
############################################################
#' Third derivative of the log density
#' Created by Stephen Jewson
#' using Deriv() by Andrew Clausen and Serguei Sokol
#' @returns 3d array
#' @inheritParams manf
gev_p1_logfddd=function (x, t, v1, v2, v3, v4) 
{
    .e3 <- x - (t * v2 + v1)
    .e4 <- v4 * .e3
    .e5 <- .e4/v3
    .e6 <- 1 + .e5
    .e7 <- 1/v4
    .e8 <- 1 + .e7
    .e9 <- .e6^.e7
    .e10 <- v3 * .e6
    .e11 <- .e6^.e8
    .e12 <- log1p(.e5)
    .e13 <- .e10^2
    .e14 <- .e7 - 1
    .e15 <- .e6^.e14
    .e16 <- 2/v4
    .e17 <- 1/.e9
    .e18 <- v3 * .e11
    .e19 <- v4 * .e8
    .e20 <- .e7 + 2
    .e21 <- .e15 * .e3
    .e22 <- .e21/v3
    .e23 <- .e6^.e16
    .e24 <- v3^2
    .e26 <- .e9 * .e12/v4
    .e27 <- .e19 - .e17
    .e28 <- .e22 - .e26
    .e29 <- .e6^.e20
    .e30 <- .e18^2
    .e31 <- v4 * .e23
    .e32 <- v3 * v4
    .e33 <- .e10 * .e27
    .e34 <- v4^2
    .e35 <- v3^3
    .e36 <- .e16 - 1
    .e37 <- .e32 * .e6
    .e38 <- .e6^.e36
    .e39 <- .e17 - 1
    .e40 <- 2 * (.e33 * .e3/.e13)
    .e41 <- .e3/.e18
    .e43 <- .e12/.e9 - v4 * .e39
    .e44 <- .e35 * .e29
    .e45 <- .e8 * .e9
    .e46 <- .e15 * .e12
    .e47 <- .e6^(.e7 - 2)
    .e48 <- (.e24 * .e29)^2
    .e50 <- .e19 * .e9 * .e3
    .e51 <- .e28/.e31
    .e52 <- 2/.e9
    .e53 <- .e24 * v4
    .e54 <- (2 * (.e37 * .e27/.e13) - 1/.e18)/.e13
    .e55 <- .e43/v4
    .e56 <- .e24 * .e11
    .e57 <- .e11 * .e12
    .e58 <- .e28/.e23
    .e59 <- .e51 + 1
    .e61 <- (.e46 + v4 * .e15)/v4 - (.e15 + v4 * .e47 * .e14 * 
        .e3/v3)
    .e62 <- .e11 - .e50/v3
    .e63 <- .e57/.e34
    .e64 <- .e44^2
    .e65 <- .e31^2
    .e67 <- .e3/.e56 + 2 * (.e33/.e13)
    .e68 <- v4 * .e28
    .e70 <- .e45 * .e3/v3
    .e71 <- .e24 * .e6
    .e72 <- .e70 - .e63
    .e73 <- 2 * (.e37 * .e3/.e13)
    .e74 <- v4 * .e38
    .e76 <- .e61/.e31 + 2 * (.e68 * .e38/.e65)
    .e77 <- .e55 + .e17
    .e80 <- .e12/.e11 - 2 * (v4/.e11)
    .e81 <- t^2
    .e85 <- .e77/.e13
    .e86 <- .e6^(1 + .e16)
    .e87 <- .e6^(.e16 - 2)
    .e88 <- .e54 - .e18 * .e20/.e48
    .e89 <- .e32 * .e8
    .e90 <- v4 * ((.e58 + 1)/v4 + 2 - .e40)
    .e94 <- .e11 * .e20 * .e3/v3 - .e29 * .e12/.e34
    .e95 <- v3 * .e29
    .e96 <- .e53 * .e6
    .e99 <- v4 * .e11 * .e20 * .e3
    .e100 <- .e76/.e71
    .e102 <- .e28 * .e12/.e74
    .e103 <- .e67/.e13
    .e104 <- .e22 + .e9
    .e105 <- .e54 + .e53 * .e11 * .e20 * .e3/.e64
    .e107 <- .e80/v4 + 1/.e11
    .e109 <- .e41 + .e19 - .e17
    .e110 <- .e73 - 2
    .e112 <- v4 * .e67/.e13
    .e113 <- .e62/.e30
    .e114 <- .e9 + .e26
    .e117 <- (.e52 - 1) * .e3/v3 - .e102
    .e118 <- 1/.e23
    .e119 <- 2 * .e55
    .e120 <- v3 * .e72
    .e121 <- .e35 * .e6
    .e122 <- .e100 + (.e41 + .e90 - .e17)/.e13
    .e127 <- .e59/.e10
    .e129 <- (.e107/.e71 + v4 * (.e85 - ((2 * (.e89 * .e6^(1 + 
        3/v4)/.e30) - .e15/v3) * .e3 - .e9) * .e8/.e30))/v4 + 
        .e19 * .e110/.e13
    .e132 <- (.e47 * .e14 * .e3/v3 - .e46/.e34) * .e3/v3 - ((.e22 - 
        .e114) * .e12/v4 + .e22)/v4
    .e135 <- (2 * .e95 - .e99)/.e48 - .e112
    .e136 <- .e105 - 1/.e44
    .e137 <- .e10^4
    .e138 <- (v3 * .e34 * .e6)^2
    .e139 <- .e96^2
    .e141 <- (.e90 - .e17)/.e13 + .e24 * .e94/.e48
    .e143 <- .e74^2
    .e144 <- (v4 * .e27 * .e3/.e13 - .e109/.e10)/v3
    .e146 <- .e41 + v4 * ((.e58 + 2)/v4 + 3 - .e40) - .e52
    .e148 <- 2 * (v3 * .e62 * .e86/.e30)
    .e153 <- v3 * .e8 * .e6
    .e154 <- v4 * .e87
    .e157 <- v4 * (.e52 + v4 * (.e40 - (2 + .e16)) - 2 * .e41)/.e13 - 
        (.e50/.e30 - 2/.e18)/.e10
    .e158 <- (.e76/.e121 + .e103) * .e3
    .e162 <- (.e132/.e31 - .e28 * (.e23 + 2 * (.e38 * .e3/v3) - 
        2 * (.e23 * .e12/v4))/.e65)/.e10 - (2 + 2 * .e51 - .e40) * 
        .e3/.e13
    .e163 <- (.e59 - .e40)/.e13
    .e164 <- .e59/.e13
    .e165 <- (.e45 * .e12 + .e9)/v4
    .e166 <- .e103 + v3 * (3 * .e95 - .e99) * .e3/.e64
    .e167 <- .e104 - .e9
    .e168 <- (2 * (.e3/(v3 * .e9)) + v4 * .e117 - ((.e61 * .e12 - 
        .e68/.e6)/.e154 + 2/.e15))/.e6
    .e169 <- .e27 * .e3
    .e170 <- .e27/.e13
    .e171 <- .e119 + .e34 * .e28 * .e87 * .e36 * .e12/.e143
    .e172 <- .e12/.e31
    .e174 <- t * v4 * .e88
    .e176 <- .e81 * v4 * .e88
    .e178 <- .e120 * .e3/.e30
    .e179 <- .e153 * .e3
    .e182 <- .e35 * .e94 * .e3/.e64
    .e184 <- -(t * .e136)
    .e186 <- (((.e165 - .e104 * .e8) * .e3 + .e120 * (2 * (.e89 * 
        .e86 * .e3/.e30) - 1))/.e30 + (.e168 + 2 - .e171)/.e37)/v4 - 
        (.e8 * .e110/.e13 + v4^3/.e138) * .e3
    .e188 <- (((.e28 * (.e118 - (.e118 + .e172)) + .e41 + 2 - 
        (.e119 + .e52))/.e10 + v4 * (.e45/.e30 - 1/.e13) * .e3)/v4 - 
        (((.e22 - (.e26 + 2 * (.e53 * .e72 * .e86/.e30))) * .e8 + 
            .e9)/.e30 + .e85) * .e3)/v4 + .e8 * (2 - .e73) * 
        .e3/.e13
    .e190 <- (.e100 + .e146/.e13) * .e3 - .e127
    .e191 <- .e158 - .e164
    .e194 <- ((.e113 + 1/.e56)/.e10 - .e112) * .e3 + .e109/.e13 - 
        .e144
    .e195 <- (.e127 - .e169/.e13)/v3
    .e197 <- ((.e107/.e121 + v4 * (.e21/.e24 + .e148) * .e8/.e30) * 
        .e3 - .e85)/v4 + (1 - .e73) * .e8/.e13
    .e199 <- .e163 - .e182
    .e206 <- ((.e80/(.e35 * v4 * .e6) + v4 * ((.e167/v3 + .e148) * 
        .e8/.e30 + .e32 * .e43/.e139))/v4 - 2 * (.e89 * .e6/.e137)) * 
        .e3 + .e8/.e13 - (.e113 + .e43/.e96)/v4
    .e209 <- (.e105 - 2/.e44) * .e3 + .e144 - .e170
    .e212 <- .e146 * .e3/.e13 - (.e59 - .e178)/.e10
    .e213 <- .e63 + 2 * (.e24 * .e72 * .e62 * .e11/.e30)
    .e215 <- .e39 * .e12/v4
    .e217 <- 1/.e34 + 2 * (.e179/.e13)
    .e218 <- 2 * .e10
    .e219 <- c(v1 = .e174, v2 = .e176, v3 = t * .e157/v3, v4 = -(t * 
        .e129))
    .e220 <- t * .e122
    .e221 <- t * .e135
    .e222 <- t * .e141
    c(v1 = c(v1 = c(v1 = v4 * .e88, v2 = .e174, v3 = .e157/v3, 
        v4 = -.e129), v2 = .e219, v3 = c(v1 = -.e136, v2 = .e184, 
        v3 = -(.e209/v3), v4 = -.e206), v4 = c(v1 = .e122, v2 = .e220, 
        v3 = .e190/v3, v4 = -.e186)), v2 = c(v1 = .e219, v2 = c(v1 = .e176, 
        v2 = t^3 * v4 * .e88, v3 = .e81 * .e157/v3, v4 = -(.e81 * 
            .e129)), v3 = c(v1 = .e184, v2 = -(.e81 * .e136), 
        v3 = -(t * .e209/v3), v4 = -(t * .e206)), v4 = c(v1 = .e220, 
        v2 = .e81 * .e122, v3 = t * .e190/v3, v4 = -(t * .e186))), 
        v3 = c(v1 = c(v1 = .e135, v2 = .e221, v3 = .e194/v3, 
            v4 = -.e197), v2 = c(v1 = .e221, v2 = .e81 * .e135, 
            v3 = t * .e194/v3, v4 = -(t * .e197)), v3 = c(v1 = .e166, 
            v2 = t * .e166, v3 = (.e166 * .e3 + 2 * (((.e169/.e10 - 
                1)/v3 + (.e170 + .e3/.e44) * .e3)/v3))/v3, v4 = -(((.e80 * 
                .e3/(v3^4 * v4 * .e6) + (v4 * .e167 * .e8 * .e3/.e24 - 
                2 * (v3 * .e62^2 * .e11/.e30))/.e30 - v4 * (.e218 - 
                .e4) * .e43/.e139)/v4 + 2 * (.e153/.e137)) * 
                .e3)), v4 = c(v1 = .e191, v2 = t * .e191, v3 = (.e158 - 
            (.e195 + .e164)) * .e3/v3, v4 = -(((((.e165 + (.e9 - 
            .e104) * .e8) * .e3/v3 - .e213)/.e30 + (.e168 + 1 - 
            .e171)/.e96)/v4 + 2 * (.e179/.e137) + .e34/.e138) * 
            .e3))), v4 = c(v1 = c(v1 = .e141, v2 = .e222, v3 = .e212/v3, 
            v4 = -.e188), v2 = c(v1 = .e222, v2 = .e81 * .e141, 
            v3 = t * .e212/v3, v4 = -(t * .e188)), v3 = c(v1 = -.e199, 
            v2 = -(t * .e199), v3 = -((.e195 + .e163 - .e182) * 
                .e3/v3), v4 = -(((((.e28 * (.e118 - .e172) + 
                .e41 + 1 - .e77)/.e71 - .e113)/v4 + (((.e114 - 
                .e22) * .e8 - .e9) * .e3/v3 - .e213)/.e30 - v3 * 
                .e43 * (.e10 + .e4)/.e139)/v4 + .e217/.e13) * 
                .e3)), v4 = c(v1 = .e162, v2 = t * .e162, v3 = .e162 * 
            .e3/v3, v4 = -((((((.e38 + v4 * (.e87 * .e36 * .e3/v3 - 
            2 * (.e38 * .e12/.e34))) * .e12/.e143 - 2 * (.e3/(.e32 * 
            .e23))) * .e28 - ((.e132 * .e12 + .e28 * .e3/.e10)/.e154 + 
            .e117 * .e3/v3)/.e6)/.e6 - ((2 * ((.e39 * .e3/v3 - 
            .e102)/.e6 - .e215) + 2 * (.e117/.e6) - 4 * .e215)/v4 + 
            .e178))/v4 + v3 * (((.e28 * .e8 - .e9/v4) * .e3/v3 - 
            ((.e70 - (.e57/v4 + 2 * .e11)/v4) * .e12 + .e9 * 
                .e3/v3)/v4)/v4 - 2 * (.e24 * .e72^2 * .e11/.e30)) * 
            .e3/.e30)/v4 + (.e217 * .e3/.e13 + v4 * (.e218 + 
            .e4)/.e138) * .e3))))
}
############################################################
#' The first derivative of the density
#' @returns Vector
#' @inheritParams manf
gev_p1_f1fa=function(x,t,v1,v2,v3,v4){

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p1_fd)
	f1=vf(x,t,v1,v2,v3,v4)
	return(f1)
}
############################################################
#' The second derivative of the density
#' @returns Matrix
#' @inheritParams manf
gev_p1_f2fa=function(x,t,v1,v2,v3,v4){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p1_fdd)
	temp1=vf(x,t,v1,v2,v3,v4)
	f2=deriv_copyfdd(temp1,nx,dim=4)
	return(f2)
}
############################################################
#' Minus the first derivative of the cdf, at alpha
#' @returns Vector
#' @inheritParams manf
gev_p1_mu1fa=function(alpha,t,v1,v2,v3,v4){
	x=qgev((1-alpha),mu=v1+v2*t,sigma=v3,xi=v4)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p1_pd)
	mu1=-vf(x,t,v1,v2,v3,v4)
	return(mu1)
}
############################################################
#' Minus the second derivative of the cdf, at alpha
#' @returns Matrix
#' @inheritParams manf
gev_p1_mu2fa=function(alpha,t,v1,v2,v3,v4){
	x=qgev((1-alpha),mu=v1+v2*t,sigma=v3,xi=v4)
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p1_pdd)
	temp1=vf(x,t,v1,v2,v3,v4)
	mu2=-deriv_copyfdd(temp1,nx,dim=4)
	return(mu2)
}
############################################################
#' The second derivative of the normalized log-likelihood
#' @returns Matrix
#' @inheritParams manf
gev_p1_ldda=function(x,t,v1,v2,v3,v4){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p1_logfdd)
	temp1=vf(x,t,v1,v2,v3,v4)
	ldd=deriv_copyldd(temp1,nx,dim=4)
	return(ldd)
}
############################################################
#' The third derivative of the normalized log-likelihood
#' @returns 3d array
#' @inheritParams manf
gev_p1_lddda=function(x,t,v1,v2,v3,v4){
	nx=length(x)

	v3=movexiawayfromzero(v3)

	vf=Vectorize(gev_p1_logfddd)
	temp1=vf(x,t,v1,v2,v3,v4)
	lddd=deriv_copylddd(temp1,nx,dim=4)
	return(lddd)
}
