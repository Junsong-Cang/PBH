/* mdot module, theory: 2303.06616
Original script in python written by ZZH @ NAOC
C version written by JSC @ SNS
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double get_xmid(double x1, double x2)
{
    double lx1, lx2, lx, x;
    lx1 = log(x1);
    lx2 = log(x2);
    lx = lx1 + (lx2 - lx1) / 2.0;
    x = exp(lx);
    return x;
}

double rb_eff_solve_kernel(double rb, double rbh, double rh, double x)
{
    double p, r;

    p = 0.75;
    r = rb / x + rbh / x * pow(x / rh, p) - 1.0;
    return r;
}

double Solve_rb_eff(double rb, double rh, double rbh)
{
    // solve f(x) = 0 for x in [rb, rh], f is defined in rb_eff_solve_kernel
    // obviously since both rb and rh are positive, it would be most efficient to solve in log-space, see get_xmid

    double x1, x, x2, lx1, lx, lx2, f1, fx, f2, dif, rtol, small;
    int count, go;

    // solve precision
    rtol = 1.0e-5;

    // Initialise
    small = 1.0e-200;
    count = 0;
    go = 1;
    x1 = rb;
    x2 = rh;
    x = get_xmid(x1, x2);

    f1 = rb_eff_solve_kernel(rb, rbh, rh, x1);
    fx = rb_eff_solve_kernel(rb, rbh, rh, x);
    f2 = rb_eff_solve_kernel(rb, rbh, rh, x2);

    if (f1 * f2 > 0)
    {
        fprintf(stderr, "Error in Solve_rb_eff: solution not found in [rb, rh]\n");
        exit(1);
    }

    while (go)
    {
        // Update x1, x2 and f1, f2
        if (f1 * fx < 0)
        {
            // f(x1) and f(x) have different signs, solution must be in [x1, x], same logic applies for [x, x2]
            // f1 = f1;
            f2 = fx;
            // x1 = x1;
            x2 = x;
        }
        else
        {
            f1 = fx;
            // f2 = f2;
            x1 = x;
            // x2 = x2;
        }

        x = get_xmid(x1, x2);
        fx = rb_eff_solve_kernel(rb, rbh, rh, x);
        count += 1;
        dif = fabs(x2 / x1 - 1);
        if (dif < rtol)
        {
            go = 0;
        }
        if (count > 100)
        {
            fprintf(stderr, "Error in Solve_rb_eff: solution not found after 100 iterations, something must have gone wrong.\n");
            exit(1);
        }
    }
    // printf("Solve_rb_eff : x = %E\n", x);
    return x;
}

void r_B_eff_fun(double z, double mbh, double veff, double rb, double *results)
{
    double G, msun, pc, mh, rbh, rb_eff, rh;
    G = 6.67259e-8;  // cm^3/g/s^2
    msun = 1.989e33; // solar mass in g
    pc = 3.086e18;   // parsec in cm

    mh = 3000 * mbh / (1. + z);                      // DM halo mass
    rh = 0.339 * 58 / (1. + z) * pow(mh, 1.0 / 3.0); // DM halo radius
    rbh = G * msun * mh / pow(veff * 1e5, 2) / pc;   // DM halo's Bondi radius
    rb_eff = rbh + rb;
    if (rb_eff < rh)
    {
        rb_eff = Solve_rb_eff(rb, rh, rbh);
    }

    results[0] = rb_eff;
    results[1] = rbh;
    results[2] = rh;
}

double lambda_zzh(double beta0, double gamma0, double r_B_eff, double r_h, double z, double lglambda_maccmax, double dlnx0)
{
    double p, dlnx, lnx_cut, dlnx_cut, x_h, m_h_rBeff, lglambda_max, lglambda_min, lglambda_h, lambda_h, lnx_cut_change, lnx, x, rho;
    double T, lnT, lnv, v, m_h_r, m_PBH, M_r, A, A0, B0, I, dlnv, dlnT;

    p = 0.75;
    dlnx = dlnx0;
    lnx_cut = 2.;
    dlnx_cut = -1e-8;

    x_h = r_h / r_B_eff;
    if (x_h < 1.0)
    {
        m_h_rBeff = 1.;
    }
    else
    {
        m_h_rBeff = pow(1 / x_h, p);
    }

    lglambda_max = 1.;
    lglambda_min = -8.;
    lglambda_h = (lglambda_max + lglambda_min) / 2.;
    lambda_h = pow(10.0, lglambda_h);

    if (lglambda_maccmax < lglambda_max)
    {
        lglambda_max = lglambda_maccmax;
        lglambda_h = lglambda_maccmax;
        lambda_h = pow(10.0, lglambda_h);
    }

    while (1)
    {
        lnx_cut_change = 2.;
        lnx = 3.;
        x = exp(lnx);
        rho = 1.;
        T = 1.;
        lnT = log(T);
        v = lambda_h / rho / pow(x, 2.0);
        lnv = log(v);

        while (1)
        {
            if (x > x_h)
            {
                m_h_r = 1.;
            }
            else
            {
                m_h_r = pow(x / x_h, p);
            }
            m_PBH = (1. + z) / 3000.;
            M_r = (m_PBH + m_h_r) / (m_PBH + m_h_rBeff);

            A = pow(v, 2.0) / T;
            A0 = (M_r / x - beta0 * v * x) / T;
            B0 = gamma0 * x / v * (1. / T - 1);
            I = (A0 - B0 - 2 * A) / (5. / 3. - A);
            dlnv = -(2 - I) * dlnx;
            dlnT = -(B0 + 2. / 3. * I) * dlnx;
            lnx = lnx + dlnx;
            lnv = lnv + dlnv;
            lnT = lnT + dlnT;
            x = exp(lnx);
            v = exp(lnv);
            T = exp(lnT);

            if (lnx < lnx_cut_change)
            {
                if (dlnx > dlnx0)
                {
                    dlnx = dlnx * 10;
                }
                lnx_cut_change = lnx_cut_change - 1;
            }

            if (lnT < 2 * lnv - 0.510825623765991) // np.log(5./3.) = 0.510825623765991
            {
                if ((lnx > lnx_cut) && (dlnx < dlnx_cut))
                {
                    dlnx = dlnx * 1e-1;
                    break;
                }
                lglambda_max = lglambda_h;
                dlnx = dlnx0;
                break;
            }
            else if ((lnx < -8.) || ((10. / 3. * T / x - M_r / pow(x, 2.0) + gamma0 * (1 - T) / v + beta0 * v) < 0))
            {
                if ((lnx > lnx_cut) && (dlnx < dlnx_cut))
                {
                    dlnx = dlnx * 1e-1;
                    break;
                }
                lglambda_min = lglambda_h;
                dlnx = dlnx0;
                break;
            }
        }
        lglambda_h = (lglambda_max + lglambda_min) / 2;
        lambda_h = pow(10.0, lglambda_h);
        if (lglambda_max - lglambda_min < 0.01)
        {
            break;
        }
    }
    return pow(10.0, lglambda_max);
}

double mdot_kernel_c(double z, double M_PBH, double x_e, double T_k, double Omega_b, int Use_halo)
{
    double h, Omb_LCDM, k_B, m_H, G, M_SUN, pc, c, sigma_T, m_e, rhoc, Y_He, n_H_z0, n_He_z0, p, rho_crit, rho_b, L_Edd, T_k_z, x_e_z;
    double n_H, n_He, n_tot, mu, cs, v_L, r_B, v_eff, Pi, rho_CMB, t_B, beta0, gamma0, lambda_ad, lambda_iso, lambda0, lambda1, M_acc;
    double r1[3], r_B_eff, r_Bh, r_h, t_B_eff, beta0_h, gamma0_h, m_acc;
    double M_acc_max, lambda1_max;
    double m_acc_max = 100.0;
    // parameters
    h = 0.6766;
    Omb_LCDM = 0.04897468161;
    k_B = 1.3806542e-16;
    m_H = 1.66e-24;   // H mass in g
    G = 6.67259e-8;   // cm^3/g/s^2
    M_SUN = 1.989e33; // solar mass in g
    pc = 3.086e18;    // parsec in cm
    c = 2.99792458e10;
    sigma_T = 6.65e-25;
    Pi = 3.141592653589793;
    m_e = 9.1093897e-28;              // electron mass in g
    rhoc = (2.7752e11) * pow(h, 2.0); // M_SUN/Mpc^3
    Y_He = 0.245;
    n_H_z0 = rhoc * Omb_LCDM * M_SUN * (1 - Y_He) / m_H / pow(1e6 * pc, 3.0);
    n_He_z0 = rhoc * Omb_LCDM * M_SUN * Y_He / (4 * m_H) / pow(1e6 * pc, 3.0);

    p = 0.75;                                     // power-law DM halo density profile
    rho_crit = 1.879e-29 * pow(h, 2.0);           // g/cm^3
    rho_b = Omega_b * rho_crit * pow(1 + z, 3.0); // g/cm^3
    L_Edd = 1.26e38 * M_PBH;                      // erg/s
    T_k_z = T_k;
    x_e_z = x_e;

    n_H = n_H_z0 * pow(1 + z, 3);
    n_He = n_He_z0 * pow(1 + z, 3.0);
    n_tot = n_H * (1 + x_e_z) + n_He;
    mu = rhoc * Omega_b * M_SUN / pow(1e6 * pc, 3.0) * pow(1 + z, 3.0) / n_tot / m_H;

    cs = sqrt(5. / 3. * k_B * T_k_z / (mu * m_H)) / 1e5;    // km/s
    v_L = fmin(1., z / 1.0e3) * 30.0;                       // km/s
    v_eff = fmax(cs, sqrt(cs * v_L));                       // km/s
    r_B = G * M_SUN * M_PBH / pow(v_eff * 1.0e5, 2.0) / pc; // pc

    // lambda
    if (Use_halo)
    {
        // PBH with DM halo
        r_B_eff_fun(z, M_PBH, v_eff, r_B, r1);
        r_B_eff = r1[0];
        r_Bh = r1[1];
        r_h = r1[2];
        rho_CMB = 4.15e-5 * pow(h, -2.0) * rho_crit * pow(1.0 + z, 4.0) * pow(c, 2.0);
        t_B_eff = r_B_eff * pc / (v_eff * 1e5);
        beta0_h = 4 * x_e_z * sigma_T * rho_CMB * t_B_eff / (3 * m_H * c);
        gamma0_h = 8 * x_e_z * sigma_T * rho_CMB * t_B_eff / (3 * m_e * c * (1 + x_e_z));
        M_acc_max = m_acc_max / pow(c, 2.0) * L_Edd;
        lambda1_max = M_acc_max / 4.0 / Pi / rho_b / pow(r_B_eff * pc, 2.0) / (v_eff * 1e5);
        lambda1 = lambda_zzh(beta0_h, gamma0_h, r_B_eff, r_h, z,log10(lambda1_max), -1e-4);
        M_acc = 4 * Pi * lambda1 * rho_b * pow(r_B_eff * pc, 2.0) * (v_eff * 1e5); // #g/s
    }
    else
    {
        rho_CMB = 4.15e-5 * pow(h, -2) * rho_crit * pow(1 + z, 4) * pow(c, 2.0);
        t_B = r_B * pc / (v_eff * 1e5);
        beta0 = 4 * x_e_z * sigma_T * rho_CMB * t_B / (3 * m_H * c);
        gamma0 = 8 * x_e_z * sigma_T * rho_CMB * t_B / (3 * m_e * c * (1 + x_e_z));
        lambda_ad = 0.25 * pow(0.6, 1.5);
        lambda_iso = 0.25 * exp(1.5);
        lambda0 = (lambda_ad + (lambda_iso - lambda_ad) * pow((pow(gamma0, 2) / (88. + pow(gamma0, 2.))), 0.22) * exp(4.5 / (3. + pow(beta0, 0.75))) / pow(sqrt(1. + beta0) + 1., 2.0)) / lambda_iso;

        M_acc = 4 * Pi * lambda0 * rho_b * pow(r_B * pc, 2.0) * (v_eff * 1e5); // g/s
    }
    m_acc = M_acc * pow(c, 2.0) / L_Edd;
    return m_acc;
}

/*
int main()
{
    double z, M_PBH, x_e, T_k, Omega_b, r;
    int Use_halo, idx, n;
    Use_halo = 1;
    z = 10;
    M_PBH = 10000;
    x_e = 1e-4;
    T_k = 10;
    Omega_b = 0.04897468161;
    Use_halo = 1;
    n = 10000;
    for (idx = 0; idx < n; idx++)
    {
        r = mdot_kernel(z, M_PBH, x_e, T_k, Omega_b, Use_halo);
    }

    printf("%f\n", r);
    return 0;
}
*/
