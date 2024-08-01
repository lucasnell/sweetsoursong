
# install.packages("diffeqr")

library(diffeqr)
library(sweetsoursong)
library(tidyverse)
library(patchwork)

de <- diffeqr::diffeq_setup()



m = 0.1
R = 10
d_yp = 1.5
d_b0 = 0.3
d_bp = 0.4
g_yp = 0.005
g_b0 = 0.02
g_bp = 0.001
L_0 = 0.01
P_max = 1.5
q = 0.5
s_0 = 0.5
h = 3
f_0 = 0.5
F_tilde = 1e3
u = 2
dt = 0.1
max_t = 250.0
Y0 = 1.0
B0 = 1.0
N0 = 1.0


# rm(m, R, d_yp, d_b0, d_bp, g_yp, g_b0, g_bp, L_0, P_max, q, s_0, h, f_0, F_tilde, u, dt, max_t, Y0, B0, N0)


one_plant_ode_jl <- function(m = 0.1,
                             R = 10,
                             d_yp = 1.5,
                             d_b0 = 0.3,
                             d_bp = 0.4,
                             g_yp = 0.005,
                             g_b0 = 0.02,
                             g_bp = 0.001,
                             L_0 = 0.01,
                             P_max = 1.5,
                             q = 0.5,
                             s_0 = 0.5,
                             h = 3,
                             f_0 = 0.5,
                             F_tilde = 1e3,
                             u = 2,
                             dt = 0.1,
                             max_t = 250.0,
                             Y0 = 1.0,
                             B0 = 1.0,
                             N0 = 1.0) {

    pars <- c(m,
              R,
              d_yp,
              d_b0,
              d_bp,
              g_yp,
              g_b0,
              g_bp,
              L_0,
              P_max,
              q,
              s_0^h,
              h,
              f_0^u,
              F_tilde,
              u)

    # f <- function(u,p,t) {
    #     du1 = p[1]*(u[2]-u[1])
    #     du2 = u[1]*(p[2]-u[3]) - u[2]
    #     du3 = u[1]*u[2] - p[3]*u[3]
    #     return(c(du1,du2,du3))
    # }
    f <- function(ybn, pars, t) {

        Y = ybn[[1]]
        B = ybn[[2]]
        N = ybn[[3]]
        .F = sum(ybn)

        m = pars[[1]]
        R = pars[[2]]
        d_yp = pars[[3]]
        d_b0 = pars[[4]]
        d_bp = pars[[5]]
        g_yp = pars[[6]]
        g_b0 = pars[[7]]
        g_bp = pars[[8]]
        L_0 = pars[[9]]
        P_max = pars[[10]]
        q = pars[[11]]
        s_0_h = pars[[12]]
        h = pars[[13]]
        f_0_u = pars[[14]]
        F_tilde = pars[[15]]
        u = pars[[16]]

        FF_u = (.F / (.F + F_tilde))^u
        phi = FF_u / (f_0_u + FF_u)
        psi = s_0_h / (s_0_h + (B / .F)^h)
        P = P_max * (q * psi + (1-q) * phi)

        PF = P / .F
        Lambda = PF / (L_0 + PF)

        gamma_y = g_yp * Lambda
        gamma_b = g_b0 + g_bp * Lambda

        delta_y = d_yp * Lambda
        delta_b = d_b0 + d_bp * Lambda

        disp_y = delta_y * Y / .F + gamma_y
        disp_b = delta_b * B / .F + gamma_b

        # double& dYdt(dxdt[0]);
        # double& dBdt(dxdt[1]);
        # double& dNdt(dxdt[2]);

        dYdt = disp_y * N - m * Y
        dBdt = disp_b * N - m * B
        dNdt = R - N * (m + disp_y + disp_b)

        return(c(dYdt, dBdt, dNdt))

    }


    ybn0 <- c(Y0, B0, N0)
    tspan <- list(0.0, max_t)
    saveat <- seq(0, max_t, dt)
    prob <- de$ODEProblem(f, ybn0, tspan, pars)
    sol <- de$solve(prob, saveat=saveat)


    udf <- sapply(sol$u,identity) |>
        t() |>
        as.data.frame() |>
        set_names(c("Y", "B", "N")) |>
        as_tibble() |>
        mutate(t = sol$t) |>
        select(t, everything())

    s_0_h = s_0^h
    f_0_u = f_0^u
    udf$P <- sapply(1:nrow(udf), \(i) {
        .F = sum(udf[i,c("Y", "B", "N")])
        F_u = .F^u
        phi = F_u / (f_0_u + F_u)
        psi = s_0_h / (s_0_h + (udf[["B"]][[i]] / .F)^h)
        udf[i,"P"] = P_max * (q * psi + (1-q) * phi)
    })

    return(udf)
}




a <- one_plant_ode_jl(m = 0.1,
                      R = 10,
                      d_yp = 1.5,
                      d_b0 = 0.3,
                      d_bp = 0.4,
                      g_yp = 0.005,
                      g_b0 = 0.02,
                      g_bp = 0.001,
                      L_0 = 0.01,
                      P_max = 1.5,
                      q = 0.5,
                      s_0 = 0.5,
                      h = 3,
                      f_0 = 0.5,
                      F_tilde = 1e3,
                      u = 2,
                      dt = 0.1,
                      max_t = 250.0,
                      Y0 = 1.0,
                      B0 = 1.0,
                      N0 = 1.0)
b <- one_plant_ode(m = 0.1,
                   R = 10,
                   d_yp = 1.5,
                   d_b0 = 0.3,
                   d_bp = 0.4,
                   g_yp = 0.005,
                   g_b0 = 0.02,
                   g_bp = 0.001,
                   L_0 = 0.01,
                   P_max = 1.5,
                   q = 0.5,
                   s_0 = 0.5,
                   h = 3,
                   f_0 = 0.5,
                   F_tilde = 1e3,
                   u = 2,
                   dt = 0.1,
                   max_t = 250.0,
                   Y0 = 1.0,
                   B0 = 1.0,
                   N0 = 1.0) |>
    as_tibble()

a; b

pa <- a |>
    pivot_longer(Y:P) |>
    ggplot(aes(t, value, color = name)) +
    geom_line()
pb <- b |>
    pivot_longer(Y:P) |>
    ggplot(aes(t, value, color = name)) +
    geom_line()
pa + pb


microbenchmark::microbenchmark(jl = one_plant_ode_jl(m = 0.1,
                                                     R = 10,
                                                     d_yp = 1.5,
                                                     d_b0 = 0.3,
                                                     d_bp = 0.4,
                                                     g_yp = 0.005,
                                                     g_b0 = 0.02,
                                                     g_bp = 0.001,
                                                     L_0 = 0.01,
                                                     P_max = 1.5,
                                                     q = 0.5,
                                                     s_0 = 0.5,
                                                     h = 3,
                                                     f_0 = 0.5,
                                                     F_tilde = 1e3,
                                                     u = 2,
                                                     dt = 0.1,
                                                     max_t = 250.0,
                                                     Y0 = 1.0,
                                                     B0 = 1.0,
                                                     N0 = 1.0),
                               cpp = one_plant_ode(m = 0.1,
                                                   R = 10,
                                                   d_yp = 1.5,
                                                   d_b0 = 0.3,
                                                   d_bp = 0.4,
                                                   g_yp = 0.005,
                                                   g_b0 = 0.02,
                                                   g_bp = 0.001,
                                                   L_0 = 0.01,
                                                   P_max = 1.5,
                                                   q = 0.5,
                                                   s_0 = 0.5,
                                                   h = 3,
                                                   f_0 = 0.5,
                                                   F_tilde = 1e3,
                                                   u = 2,
                                                   dt = 0.1,
                                                   max_t = 250.0,
                                                   Y0 = 1.0,
                                                   B0 = 1.0,
                                                   N0 = 1.0),
                               times = 10)

