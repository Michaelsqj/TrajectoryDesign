#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mex.h"

#define PI 3.14159265358979

double findgtwist(double gx, double gy, double kx, double ky, double phi)
{

    double adk, aik, anglediff, outsin, outcos, Gtwist;

    adk = atan2(gy, gx);
    aik = atan2(ky, kx) + PI / 2;
    anglediff = aik - adk;
    anglediff = (anglediff > 2 * PI) ? anglediff - 2 * PI : anglediff;
    anglediff = (anglediff > PI) ? anglediff - PI : anglediff;
    anglediff = (anglediff < PI) ? anglediff + 2 * PI : anglediff;
    anglediff = (anglediff < 0) ? anglediff + PI : anglediff;
    outsin = sin(anglediff);
    outcos = cos(anglediff);
    if (outsin == 0)
        Gtwist = -1;
    else
        Gtwist = outcos / (outsin * sqrt(1 + tan(phi) * tan(phi)));

    return Gtwist;
}

int rootof(double a, double b, double c, double *s1, double *s2)
{
    double descrim;
    if (a == 0)
    {
        if (b == 0)
        {
            *s1 = 0;
            *s2 = 0;
            return 0;
        }
        else
        {
            *s1 = -c / b;
            *s2 = -c / b;
            return 1;
        }
    }
    else
    {
        descrim = b * b - 4 * a * c;
        if (descrim < 0)
        {
            *s1 = 0;
            *s2 = 0;
            return 0;
        }
        else
        {
            *s1 = (-b + sqrt(descrim)) / 2 / a;
            *s2 = (-b - sqrt(descrim)) / 2 / a;
            return 2;
        }
    }
}

int circcirc(double d, double e, double k, double s, double type, double *x, double *y)
{

    double s1;
    double s2;
    double t;
    double ade;
    double tx1, tx2, ty1, ty2;
    double a, b, c;
    int numsols;
    if (e == 0)
        ade = 100;
    else
        ade = fabs(d / e);

    if ((d == 0) && (e == 0))
    {
        *x = 0;
        *y = 0;
        return 0;
    }
    else if (ade < 1)
    {
        t = e;
        e = d;
        d = t;
        a = 4 * d * d + 4 * e * e;
        b = (-4 * e * e * e - 4 * k * k * e + 4 * e * s * s - 4 * d * d * e);
        c = d * d * d * d + k * k * k * k - 2 * e * e * s * s + 2 * d * d * e * e - 2 * d * d * s * s + s * s * s * s - 2 * k * k * s * s - 2 * k * k * d * d + e * e * e * e + 2 * k * k * e * e;
        numsols = rootof(a, b, c, &s1, &s2);
        if (numsols > 0)
        {
            tx1 = s1;
            ty1 = -0.5 * (2 * e * s1 - k * k - d * d - e * e + s * s) / d;
            tx2 = s2;
            ty2 = -0.5 * (2 * e * s2 - k * k - d * d - e * e + s * s) / d;
        }
        else
        {
            *y = 0;
            *x = 0;
            return 0;
        }
        t = e;
        e = d;
        d = t;
    }
    else
    {
        a = 4 * d * d + 4 * e * e;
        b = (-4 * e * e * e - 4 * k * k * e + 4 * e * s * s - 4 * d * d * e);
        c = d * d * d * d + k * k * k * k - 2 * e * e * s * s + 2 * d * d * e * e - 2 * d * d * s * s + s * s * s * s - 2 * k * k * s * s - 2 * k * k * d * d + e * e * e * e + 2 * k * k * e * e;
        numsols = rootof(a, b, c, &s1, &s2);
        if (numsols > 0)
        {
            ty1 = s1;
            tx1 = -0.5 * (2 * e * s1 - k * k - d * d - e * e + s * s) / d;
            ty2 = s2;
            tx2 = -0.5 * (2 * e * s2 - k * k - d * d - e * e + s * s) / d;
        }
        else
        {
            *y = 0;
            *x = 0;
            return 0;
        }
    }

    if ((type * (-tx1 * e + ty1 * d)) < 0)
    {
        *x = tx2;
        *y = ty2;
    }
    else
    {
        *x = tx1;
        *y = ty1;
    }
    return numsols;
}

int intlinecirc(double a, double b, double d, double e, double r, double *s1, double *s2)
{

    return rootof((a * a + b * b), (-2 * e * b - 2 * d * a), (d * d + e * e - r * r), s1, s2);
}

int perpgrad(double a, double b, double c, double d, double e, double f,
             double kx, double ky, double kz, double phi, double smax,
             double gmax, double doplot, double *x, double *y, double *z)
{

    double vpx, vpy, kxn, kyn, kzn, kxy;
    double mx, my;
    int numsols;

    vpx = -e / (sqrt(d * d + e * e));
    vpy = d / (sqrt(d * d + e * e));

    *x = d + smax * vpx;
    *y = e + smax * vpy;

    kxn = kx + *x;
    kyn = ky + *y;
    kzn = tan(phi) * sqrt(kxn * kxn + kyn * kyn);
    *z = kzn - kz;

    if (((*z - f) < -smax) || ((*z - f) > 0))
    {
        if ((*z - f) > 0)
            *z = f;
        else
            *z = f - smax;

        kxy = (kz + *z) / tan(phi);

        numsols = circcirc((kx + d), (ky + e), kxy, smax, 1, &mx, &my);
        if (numsols == 0)
        {
            mexPrintf("Major Error! Circcirc failed\n");
            return 0;
        }
        else
        {
            *x = mx - kx;
            *y = my - ky;
            return numsols;
        }
    }
    return 2;
}

int linematch(double a, double b, double c, double d, double e, double f, double kx, double ky, double kz,
              double phi, double smax, double gmax, double doplot, double *x, double *y, double *z)
{

    int numsols, numsolg, numsolsmax, numsolsmin;
    double g1, gext2, smaxz, smaxz2, sminz, sminz2;
    double s1, s2;
    double sactual;
    double kxy;
    double kxn, kyn, kzn;
    int fail;
    numsols = intlinecirc(a, b, d, e, smax, &s1, &s2);
    if (s2 < 0)
        s2 = 0;

    numsolg = intlinecirc(a, b, 0, 0, gmax, &g1, &gext2);

    kxy = (kz + ((gmax > (f + smax)) ? (f + smax) : gmax)) / tan(phi);
    numsolsmax = intlinecirc(a, b, -kx, -ky, kxy, &smaxz, &smaxz2);

    kxy = (kz + (((f - smax) > 0) ? f - smax : 0)) / tan(phi);
    numsolsmin = intlinecirc(a, b, -kx, -ky, kxy, &sminz, &sminz2);

    sactual = (s1 > g1) ? g1 : s1;
    sactual = (sactual > smaxz) ? smaxz : sactual;
    if ((sactual < ((s2 > sminz) ? s2 : sminz)) || !numsols || !numsolg || !numsolsmax || !numsolsmin)
        fail = 1;
    else
        fail = 0;

    *x = a * sactual;
    *y = b * sactual;
    kxn = kx + *x;
    kyn = ky + *y;
    kzn = tan(phi) * sqrt(kxn * kxn + kyn * kyn);
    *z = kzn - kz;

    return fail;
}

int gtasub(double *pk,
           double *pg,
           int ai,
           int *perpstartai,
           int *perpendai,
           int numperp,
           double phi,
           double FOV1,
           double FOV2,
           double FOV3,
           double FOV4,
           double dcf,
           double kmax,
           double NINT,
           double GTWISTADJUST,
           int numpt,
           double Ts,
           double Smax,
           double Gmax,
           int *paiout,
           int *pnummatch)
{

    int fail;
    double krai;
    double krxy;

    double *pkx, *pky, *pkz, *pgx, *pgy, *pgz;
    double dscale, Gtwist, theta, xout, yout, zout;
    double out, rangle;
    double FOVa;
    double FOVb;
    int perptest;
    int bi;

    pkx = pk;
    pky = pk + numpt;
    pkz = pk + 2 * numpt;
    pgx = pg;
    pgy = pg + numpt;
    pgz = pg + 2 * numpt;

    fail = 0;
    *pnummatch = 0;
    krai = sqrt((pkx[ai]) * (pkx[ai]) + (pky[ai]) * (pky[ai]) + (pkz[ai]) * (pkz[ai]));
    krxy = sqrt((pkx[ai]) * (pkx[ai]) + (pky[ai]) * (pky[ai]));

    while ((krai < kmax) && (ai < numpt - 1) && (fail == 0))
    {
        /*FOVa = FOV1-(FOV1-FOV3)*krai/kmax;
    FOVb = FOV2-(FOV2-FOV4)*krai/kmax;*/
        FOVa = FOV1;
        FOVb = FOV2;
        dscale = ((1 - dcf) + dcf * krai / kmax);
        dscale = 1 - 4 * FOVa * FOVa / FOV1 / FOV1 * dscale * dscale;
        dscale = (dscale < 0) ? 0 : dscale;
        dscale = 1 / (1 + sqrt(dscale));
        /*dscale = ((1-dcf)+dcf*krai/kmax);*/
        /*Gtwist = pow(2*PI*krai*cos(phi)*dscale/NINT,2)-pow(1/FOVa,2);*/
        Gtwist = pow(2 * PI * krxy * dscale / NINT, 2) - pow(1 / FOVa, 2);
        Gtwist = (Gtwist < 0) ? 0 : Gtwist;
        Gtwist = sqrt(Gtwist) * FOVb;
        Gtwist = Gtwist * GTWISTADJUST;

        perptest = 0;
        for (bi = 0; bi < numperp; bi++)
            perptest = perptest || ((ai >= perpstartai[bi]) && (ai < perpendai[bi]));
        if (perptest)
        {

            *pnummatch = 0;
            perpgrad(-pgx[ai], pgy[ai], pgz[ai], pgx[ai], pgy[ai], pgz[ai], pkx[ai] / 4258.0 / Ts, pky[ai] / 4258.0 / Ts, pkz[ai] / 4258.0 / Ts, phi, Smax, Gmax, 0, &xout, &yout, &zout);
            pgx[ai + 1] = xout;
            pgy[ai + 1] = yout;
            pgz[ai + 1] = zout;
            out = findgtwist(pgx[ai + 1], pgy[ai + 1], pkx[ai], pky[ai], phi);
            if (out < 0)
                fail = 2;
            else if (out < Gtwist)
                fail = 1;
        }
        else
        {
            theta = (Gtwist == 0) ? PI / 2.0 : atan(1 / Gtwist / sqrt(1 + tan(phi) * tan(phi)));
            rangle = atan2(pky[ai], pkx[ai]) + (PI / 2 - theta);
            fail = linematch(cos(rangle), sin(rangle), 0, pgx[ai], pgy[ai], pgz[ai], pkx[ai] / 4258.0 / Ts, pky[ai] / 4258.0 / Ts, pkz[ai] / 4258.0 / Ts, phi, Smax, Gmax, 0, &xout, &yout, &zout);
            if (fail != 1)
            {
                pgx[ai + 1] = xout;
                pgy[ai + 1] = yout;
                pgz[ai + 1] = zout;
                *pnummatch = *pnummatch + 1;
            }
            else
            {
                out = findgtwist(pgx[ai], pgy[ai], pkx[ai], pky[ai], phi);
                if (out >= Gtwist)
                {
                    fail = 0;
                    pgx[ai + 1] = pgx[ai];
                    pgy[ai + 1] = pgy[ai];
                }
            }
        }

        if (fail != 1)
        {
            pkx[ai + 1] = pkx[ai] + pgx[ai + 1] * 4258 * Ts;
            pky[ai + 1] = pky[ai] + pgy[ai + 1] * 4258 * Ts;
            pkz[ai + 1] = tan(phi) * sqrt(pkx[ai + 1] * pkx[ai + 1] + pky[ai + 1] * pky[ai + 1]);
            pgz[ai + 1] = (pkz[ai + 1] - pkz[ai]) / 4258.0 / Ts;
            ai++;
            krai = sqrt((pkx[ai]) * (pkx[ai]) +
                        (pky[ai]) * (pky[ai]) +
                        (pkz[ai]) * (pkz[ai]));
            krxy = sqrt((pkx[ai]) * (pkx[ai]) + (pky[ai]) * (pky[ai]));
        }
    }
    *paiout = ai;
    return fail;
}

int whirlcone(double phi,
              double FOV1,
              double FOV2,
              double FOV3,
              double FOV4,
              double dcf,
              double kmax,
              double NINT,
              double GTWISTADJUST,
              double numpt,
              double Ts,
              double Smax,
              double Gmax,
              double *pg,
              double *pk,
              double *pai)
{

    int perpstartai[1000];
    int perpendai[1000];
    int numperp;
    int nummatch;
    int aiout;
    int result;
    int saimin;
    int saimax;
    int curai;
    int minai;

    if (phi > PI / 4)
    {
        pg[0] = Smax * tan(PI / 2 - phi);
        pg[(int)numpt] = 0;
        pg[2 * (int)numpt] = Smax;
    }
    else
    {
        pg[0] = Smax;
        pg[(int)numpt] = 0;
        pg[2 * (int)numpt] = Smax * tan(phi);
    }
    pk[0] = pg[0] * 4258 * Ts;
    pk[(int)numpt] = pg[(int)numpt] * 4258 * Ts;
    pk[2 * (int)numpt] = pg[2 * (int)numpt] * 4258 * Ts;
    numperp = 0;
    minai = 0;

    result = gtasub(pk, pg, minai, perpstartai, perpendai, numperp, phi, FOV1, FOV2, FOV3, FOV4, dcf, kmax, NINT, GTWISTADJUST, numpt, Ts, Smax, Gmax, &aiout, &nummatch);
    saimin = 1;
    saimax = aiout;
    curai = aiout;
    minai = 1;
    numperp = 1;
    perpstartai[numperp - 1] = 10000;
    perpendai[numperp - 1] = 10000;

    while (result != 0)
    {
        while ((saimax - saimin) > 1)
        {
            perpstartai[numperp - 1] = (saimax - saimin) / 2 + saimin;
            result = gtasub(pk, pg, minai, perpstartai, perpendai, numperp, phi, FOV1, FOV2, FOV3, FOV4, dcf, kmax, NINT, GTWISTADJUST, numpt, Ts, Smax, Gmax, &aiout, &nummatch);
            if (result == 1)
                saimax = perpstartai[numperp - 1];
            else
                saimin = perpstartai[numperp - 1];
        }
        perpstartai[numperp - 1] = saimin;
        result = gtasub(pk, pg, minai, perpstartai, perpendai, numperp, phi, FOV1, FOV2, FOV3, FOV4, dcf, kmax, NINT, GTWISTADJUST, numpt, Ts, Smax, Gmax, &aiout, &nummatch);
        result = 1;
        nummatch = 0;
        perpendai[numperp - 1] = curai;
        while ((result == 1) && (nummatch < 4))
        {
            perpendai[numperp - 1]++;
            result = gtasub(pk, pg, minai, perpstartai, perpendai, numperp, phi, FOV1, FOV2, FOV3, FOV4, dcf, kmax, NINT, GTWISTADJUST, numpt, Ts, Smax, Gmax, &aiout, &nummatch);
            if (perpendai[numperp - 1] > numpt)
                return -1;
        }
        saimin = curai;
        saimax = aiout;
        numperp++;
        perpendai[numperp - 1] = 10000;
        minai = curai;
        curai = aiout;
    }
    numperp--;
    *pai = aiout + 1;
    return 1;
}

/* mexFunction is the gateway routine for the MEX-file. */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    int mrows;
    int ncols;
    double *pr;
    double phi;
    double FOV1;
    double FOV2;
    double FOV3;
    double FOV4;
    double dcf;
    double kmax;
    double NINT;
    double GTWISTADJUST;
    double numpt;
    double Ts;
    double Smax;
    double Gmax;
    double *pg;
    double *pk;
    double *pai;

    pr = mxGetPr(prhs[0]);
    phi = *pr;

    pr = mxGetPr(prhs[1]);
    FOV1 = *pr++;
    if (mxGetNumberOfElements(prhs[1]) < 2)
        FOV2 = FOV1;
    else
        FOV2 = *pr++;

    if (mxGetNumberOfElements(prhs[1]) < 4)
    {
        FOV3 = FOV1;
        FOV4 = FOV2;
    }
    else
    {
        FOV3 = *pr++;
        FOV4 = *pr;
    }

    pr = mxGetPr(prhs[2]);
    dcf = *pr;

    pr = mxGetPr(prhs[3]);
    kmax = *pr;

    pr = mxGetPr(prhs[4]);
    NINT = *pr;

    pr = mxGetPr(prhs[5]);
    GTWISTADJUST = *pr;

    if (nrhs < 7)
        numpt = 2000;
    else
    {
        pr = mxGetPr(prhs[6]);
        numpt = *pr;
        if (numpt > 10000)
            numpt = 10000;
    }

    /*mexPrintf("Found numpt = %g\n",numpt);*/

    if (nrhs < 8)
        Ts = 0.000002;
    else
    {
        pr = mxGetPr(prhs[7]);
        Ts = *pr;
    }

    if (nrhs < 9)
        Smax = Ts * 15000;
    else
    {
        pr = mxGetPr(prhs[8]);
        Smax = *pr;
    }

    if (nrhs < 10)
        Gmax = 3.89;
    else
    {
        pr = mxGetPr(prhs[9]);
        Gmax = *pr;
    }

    mrows = 4;
    ncols = 3;
    /* Look at each input (right-hand-side) argument. */

    plhs[0] = mxCreateDoubleMatrix((int)numpt, 3, mxREAL);
    plhs[1] = mxCreateDoubleMatrix((int)numpt, 3, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

    pg = mxGetPr(plhs[0]);
    pk = mxGetPr(plhs[1]);
    pai = mxGetPr(plhs[2]);

    whirlcone(phi, FOV1, FOV2, FOV3, FOV4, dcf, kmax, NINT, GTWISTADJUST, numpt, Ts, Smax, Gmax, pg, pk, pai);
}