import math
import numpy as np
import pandas as pd
import brine_prop as bp


def third_type_boundary(
    nlayer,
    alfa,
    tbrine,
    tnext,
    tside0,
    dtime,
    n2e,
    twell,
    nc,
    z,
    mat,
    am,
    nodes_in_elements,
):

    if 1:
        nw = int(twell[nlayer, 0])
        r = nc[nw, 1]
        dr = nc[nw + len(z), 1] - nc[nw, 1]

        if (nodes_in_elements[nw, 2] >= 0) and (nodes_in_elements[nw, 3] >= 0):
            x1 = nodes_in_elements[nw, 2]
            x2 = nodes_in_elements[nw, 3]
            a = (
                (nc[nw + 1, 2] - nc[nw, 2]) * am[mat[x1], 1]
                + (nc[nw, 2] - nc[nw - 1, 2]) * am[mat[x2], 1]
            ) / (nc[nw + 1, 2] - nc[nw - 1, 2])
            lbd = (
                (nc[nw + 1, 2] - nc[nw, 2]) * am[mat[x1], 2]
                + (nc[nw, 2] - nc[nw - 1, 2]) * am[mat[x2], 2]
            ) / (nc[nw + 1, 2] - nc[nw - 1, 2])

        if (nodes_in_elements[nw, 2] >= 0) and (nodes_in_elements[nw, 3] < 0):
            x = nodes_in_elements[nw, 2]
            a = am[mat[x], 1]
            lbd = am[mat[x], 2]

        if (nodes_in_elements[nw, 2] < 0) and (nodes_in_elements[nw, 3] >= 0):
            x = nodes_in_elements[nw, 3]
            a = am[mat[x], 1]
            lbd = am[mat[x], 2]

        m = dr**2.0 / (a * dtime)
        c = (alfa * dr) / lbd
        tsn = (1 / m) * (
            2.0 * c * tbrine
            + (m - 2.0 * (c + 1.0 - dr / (2.0 * r))) * tside0
            + 2.0 * tnext * (1.0 - dr / (2.0 * r))
        )
    else:
        tsn = tbrine

    return tsn


def prod_well(
    m,
    t,
    p,
    s,
    tr,
    trock_matrix,
    dtime,
    fi,
    acc,
    n2e,
    nc,
    z,
    mat,
    am,
    nodes_in_elements,
):
    g = 9.80655
    out = np.zeros([len(tr), 7], dtype=float)
    trock = np.zeros([len(tr)], dtype=float)
    trock2 = np.zeros([len(tr)], dtype=float)
    tbrine1 = np.zeros([len(tr)], dtype=float)
    tbrine2 = np.zeros([len(tr)], dtype=float)
    pbrine1 = np.zeros([len(tr)], dtype=float)
    pbrine2 = np.zeros([len(tr)], dtype=float)

    for i in range(len(tr), 0, -1):
        n = i - 1
        trock[n] = tr[n, 3]

        if n == len(tr) - 1:
            tbrine1[n] = t
            pbrine1[n] = p

        else:
            tbrine1[n] = tbrine2[n + 1]
            pbrine1[n] = pbrine2[n + 1]

        w = m / (
            bp.density(tbrine1[n] + 273.15, pbrine1[n], s)
            * 3.141592654
            * tr[n, 1] ** 2.0
        )
        pDynamic = (
            0.5 * bp.density(tbrine1[n] + 273.15, pbrine1[n], s) * w**2.0
        )
        roB = bp.density(tbrine1[n] + 273.15, pbrine1[n] - pDynamic, s)
        cB = bp.c_heat(tbrine1[n] + 273.15, pbrine1[n] - pDynamic, s)
        if n > 0:
            hzone = tr[n, 2] - tr[n - 1, 2]
        else:
            hzone = 0.0
        dp_flow_res = bp.dp_flow(
            tbrine1[n] + 273.15, pbrine1[n], s, fi, w, hzone, tr[n, 1] * 2.0
        )
        alfaB = bp.heat_convection(
            tbrine1[n] + 273.15,
            pbrine1[n] - pDynamic,
            s,
            2.0 * tr[n, 1],
            w,
            hzone,
            tr[n, 3] + 273.15,
        )
        A = alfaB * (math.pi * 2.0 * tr[n, 1] * hzone) / ((m + 1e-3) * cB)
        if hzone > 0.0:
            tbrine2[n] = (
                2.0 * tbrine1[n] + A * (tr[n, 3] + tr[n - 1, 3] - tbrine1[n])
            ) / (A + 2.0)
        else:
            tbrine2[n] = tbrine1[n]
        pbrine2[n] = (
            pbrine1[n] - (g * roB * (tr[n, 2] - tr[n - 1, 2])) - dp_flow_res
        )

        if n == 0:
            pbrine2[n] = pbrine1[n]

        if n == 1:
            print(
                f"tbrine1[n]:{n}, {round(tbrine1[n],3)}, {round(pbrine1[n]/1e5,3)}, {dp_flow_res}"
            )

        trock2[n] = third_type_boundary(
            n,
            alfaB,
            tbrine1[n],
            trock_matrix[int(tr[n, 0] + len(z))],
            trock[n],
            dtime,
            n2e,
            tr,
            nc,
            z,
            mat,
            am,
            nodes_in_elements,
        )

        if i == 1:
            PR = [
                tbrine1[0],
                pbrine1[0],
                w[0],
                pDynamic[0],
                roB[0],
                cB[0],
                dp_flow_res[0],
                m,
                alfaB[0],
            ]

        out[n, 0] = trock[n]
        if n == len(tr) - 1:
            out[n, 1] = t
        else:
            out[n, 1] = trock2[n]

        out[n, 2] = tbrine1[n]
        out[n, 3] = tbrine2[n]
        out[n, 4] = pbrine1[n]
        out[n, 5] = pbrine2[n]
        out[n, 6] = alfaB

    df = pd.DataFrame(out)
    df.to_excel(
        "prod_well_output.xlsx",
        index=False,
        header=[
            "trock",
            "trock2",
            "tbrine1",
            "tbrine2",
            "pbrine1",
            "pbrine2",
            "alfaB",
        ],
    )

    return [out, PR]


def inj_well(
    m,
    t,
    p,
    s,
    tr,
    trock_matrix,
    dtime,
    fi,
    acc,
    n2e,
    nc,
    z,
    mat,
    am,
    nodes_top,
    nodes_in_elements,
):
    g = 9.80655
    out = np.zeros([len(tr), 7], dtype=float)
    trock = np.zeros([len(tr)], dtype=float)
    trock2 = np.zeros([len(tr)], dtype=float)
    tbrine1 = np.zeros([len(tr)], dtype=float)
    tbrine2 = np.zeros([len(tr)], dtype=float)
    pbrine1 = np.zeros([len(tr)], dtype=float)
    pbrine2 = np.zeros([len(tr)], dtype=float)

    for i in range(0, len(tr)):
        n = i

        trock[n] = tr[n, 3]

        if n == 0:
            tbrine1[n] = t
            pbrine1[n] = p
        else:
            tbrine1[n] = tbrine2[n - 1]
            pbrine1[n] = pbrine2[n - 1]

        w = m / (
            bp.density(tbrine1[n] + 273.15, pbrine1[n], s)
            * (math.pi * tr[n, 1] ** 2.0)
        )
        pDynamic = (
            0.5 * bp.density(tbrine1[n] + 273.15, pbrine1[n], s) * w**2.0
        )
        roB = bp.density(tbrine1[n] + 273.15, pbrine1[n] - pDynamic, s)
        cB = bp.c_heat(tbrine1[n] + 273.15, pbrine1[n] - pDynamic, s)
        if n > 0:
            hzone = tr[n, 2] - tr[n - 1, 2]
        else:
            hzone = 0.0
        dp_flow_res = bp.dp_flow(
            tbrine1[n] + 273.15, pbrine1[n], s, fi, w, hzone, tr[n, 1] * 2.0
        )
        alfaB = bp.heat_convection(
            tbrine1[n] + 273.15,
            pbrine1[n] - pDynamic,
            s,
            2.0 * tr[n, 1],
            w,
            hzone,
            tr[n, 3] + 273.15,
        )

        A = alfaB * (math.pi * 2.0 * tr[n, 1] * hzone) / (m * cB)

        if n < (len(tr) - 1):
            tbrine2[n] = (
                2.0 * tbrine1[n] + A * (tr[n, 3] + tr[n + 1, 3] - tbrine1[n])
            ) / (A + 2.0)
        else:
            tbrine2[n] = (
                2.0 * tbrine1[n] + A * (2.0 * tr[n, 3] - tbrine1[n])
            ) / (A + 2.0)

        if n == (len(tr) - 1):
            pbrine2[n] = pbrine1[n] + (g * roB * hzone) - dp_flow_res
        else:
            pbrine2[n] = pbrine1[n] + (g * roB * hzone) - dp_flow_res

        trock2[n] = third_type_boundary(
            n,
            alfaB,
            tbrine1[n],
            trock_matrix[int(tr[n, 0]) + len(z)],
            trock[n],
            dtime,
            n2e,
            tr,
            nc,
            z,
            mat,
            am,
            nodes_in_elements,
        )

        out[n, 0] = trock[n]

        if n == 0:
            out[n, 1] = t
        else:
            out[n, 1] = trock2[n]

        out[n, 2] = tbrine1[n]
        out[n, 3] = tbrine2[n]
        out[n, 4] = pbrine1[n]
        out[n, 5] = pbrine2[n]
        out[n, 6] = alfaB

    return out


def twell_convert_twell0(twell0, z):
    X = np.zeros([len(twell0), 4])
    for nw in range(0, len(twell0)):
        X[nw, 0] = twell0[nw, 0]
        X[nw, 1] = twell0[nw, 1]
        X[nw, 2] = twell0[nw, 2]
        X[nw, 3] = z[nw - len(twell0), 1]
    return X
