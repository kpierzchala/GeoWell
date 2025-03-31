import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import fdm as fdm
import wells as pw

# from openpyxl import load_workbook


def nodes_bottom(nc, z, r):
    nn = np.zeros([len(r)], dtype=int)
    for n in range(0, len(r)):
        nn[n] = len(z) * (n + 1) - 1
    return nn


def nodes_top(nc, z, r):
    nn = np.zeros([len(r)], dtype=int)
    for n in range(0, len(r)):
        nn[n] = len(z) * n
    return nn


def prt_mtx(mtx):
    formatted_mtx = [[f"{element:.4f}" for element in row] for row in mtx]
    for row in formatted_mtx:
        print(" ".join(row))
    return 1


def rows(a):
    out = np.shape(a)[0]
    return out


def cols(a):
    out = np.shape(a)[1]
    return out


def interp_1d(a, x):

    if x == a[len(a) - 1, 0]:
        y = a[len(a) - 1, 1]
    else:
        n = 0
        while a[n, 0] <= x:
            n = n + 1
        n = n - 1
        if a[n, 0] == x:
            y = a[n, 1]
        else:
            y = ((a[n + 1, 1] - a[n, 1]) / (a[n + 1, 0] - a[n, 0])) * (
                x - a[n, 0]
            ) + a[n, 1]
    return y


def seek_actv_nodes(nc, fix_nodes):
    actv_nodes = []

    for i in range(1, len(nc) - 1):
        if i not in fix_nodes:
            actv_nodes.append(i)
    return actv_nodes


def divide_z(x1, x2, nz):
    dz = (x2 - x1) / nz
    z2 = []
    z2.append(x1)
    wsp_x = x1

    for i in range(1, int(nz + 1)):
        z2.append(wsp_x + dz)
        wsp_x = wsp_x + dz

    return z2


def divide_from_to(x1, x2, dx1, ndx):
    xa = x1
    dxa = dx1
    wsp_x = []
    wsp_x.append(x1)

    while xa + dxa * ndx < x2:
        xa = xa + dxa * ndx
        dxa = dxa * ndx
        wsp_x.append(xa)
    last = (x2 + wsp_x[len(wsp_x) - 2]) / 2.0
    wsp_x[len(wsp_x) - 1] = last

    wsp_x.append(x2)
    return wsp_x


def well_geometry_in(file_name):
    data = pd.read_excel(
        file_name, sheet_name="geo", usecols="A:I", skiprows=0
    )
    out = data.values
    return out


def prof_temp_start(file_name, IN):

    if np.any(np.isnan(IN)):

        # if temperature distribution is empty (in the input file) then it will be calculated and filled by the Z coordinate
        PrmMat = pd.read_excel(
            file_name, sheet_name="mat", usecols="A:D", skiprows=0
        )
        prm = PrmMat.values

        # estimation of geothermal power denisity q [W/m²]
        # below R is the thermal reseistivity factor [(m² K)/W], not the radius
        R = 0.0
        for i in range(0, rows(IN)):
            R = (IN[i, 1] - IN[i, 0]) / prm[int(IN[i, 7]) - 1, 1] + R
        q = (1.0 / R) * (IN[np.shape(IN)[0] - 1, 3] - IN[0, 2])

        # we have q [W/m2], now it is the time to calculate temperatures counting from the Earth surface

        IN[0, 3] = (
            IN[0, 2] + q * (IN[0, 1] - IN[0, 0]) / prm[int(IN[0, 7]) - 1, 1]
        )

        for i in range(1, rows(IN)):
            IN[i, 2] = IN[i - 1, 3]
            IN[i, 3] = (
                IN[i, 2]
                + q * (IN[i, 1] - IN[i, 0]) / prm[int(IN[i, 7]) - 1, 1]
            )
    else:
        IN = IN

    # thermal profile curve
    th = np.zeros((rows(IN) + 1, 2))
    th[0, 0] = IN[0, 0]
    th[0, 1] = IN[0, 2]
    # due to in the input file we operate on diameters we recalculate it to radius, simply dividing by 2 (collumns 4, 6 and 8, counting from 0)
    for i in range(0, len(IN)):
        IN[i, 4] = IN[i, 4] * 0.5
        IN[i, 6] = IN[i, 6] * 0.5
        IN[i, 8] = IN[i, 8] * 0.5
    # the loop below allows you to see temperature profile at the begining
    for i in range(0, rows(IN)):
        th[i + 1, 0] = IN[i, 1]
        th[i + 1, 1] = IN[i, 3]
    if 0:  # if 1 - shows the temperature profile in the file defined a well
        plt.plot(th[:, 0], th[:, 1], marker="o", linestyle="-")
        plt.grid(True)
        plt.show()
    return IN

    if np.any(np.isnan(IN)):
        prm_mat = pd.read_excel(
            file_name, sheet_name="mat", usecols="A:D", skiprows=0
        )
        prm = prm_mat.values
        R = 0.0
        for i in range(0, rows(IN)):
            R = (IN[i, 1] - IN[i, 0]) / prm[int(IN[i, 7]) - 1, 1] + R
        q = (1.0 / R) * (IN[np.shape(IN)[0] - 1, 3] - IN[0, 2])

        IN[0, 3] = (
            IN[0, 2] + q * (IN[0, 1] - IN[0, 0]) / prm[int(IN[0, 7]) - 1, 1]
        )

        for i in range(1, rows(IN)):
            IN[i, 2] = IN[i - 1, 3]
            IN[i, 3] = (
                IN[i, 2]
                + q * (IN[i, 1] - IN[i, 0]) / prm[int(IN[i, 7]) - 1, 1]
            )

        th = np.zeros((rows(IN) + 1, 2))
        th[0, 0] = IN[0, 0]
        th[0, 1] = IN[0, 2]

        for i in range(0, len(IN)):
            IN[i, 4] = IN[i, 4] * 0.5
            IN[i, 6] = IN[i, 6] * 0.5
            IN[i, 8] = IN[i, 8] * 0.5

        for i in range(0, rows(IN)):
            th[i + 1, 0] = IN[i, 1]
            th[i + 1, 1] = IN[i, 3]
        if 0:
            plt.plot(th[:, 0], th[:, 1], marker="o", linestyle="-")
            plt.grid(True)
            plt.show()
        return IN


def skip_small_differences(a, difference, n_digit):
    b = np.array([])
    b = np.append(b, a[0])
    for i in range(1, len(a)):
        if abs(a[i] - a[i - 1]) > difference:
            b = np.append(b, round(a[i], n_digit))

    return b


def grid_gener(geom, max_z, nr):
    # --- R generation, R coordinates in r2 variable
    # IN - data origined form geometry definition, input file. Please remember thah in the input file we have diameter, not the radius.
    # maxZ - maximal allowed step on Z axis,
    # nR - increasement of R behind the defined zone (in the IN variable)
    r = np.zeros((2 * rows(geom)) + 1, dtype=float)
    r[0] = 0.0
    for i in range(0, rows(geom)):
        r[i + 1] = geom[i, 4]
    for i in range(rows(geom), 2 * rows(geom)):
        r[i + 1] = geom[i - rows(geom), 6]
    #
    r2 = np.array([])  # define the vector
    r2 = np.unique(r)  # delete duplicates
    r2 = np.sort(r2, axis=0)  # sort the data
    # let's divide the rest of the segment
    rNext = divide_from_to(
        r2[rows(r2) - 1],
        geom[0, cols(geom) - 1],
        (r2[rows(r2) - 1] - r2[rows(r2) - 2]),
        nr,
    )
    for i in range(1, rows(rNext)):
        r2 = np.append(r2, rNext[i])
    #
    if 0:  # if = 1, show the grid as x-y plot by R
        plt.plot(r2, r2, marker="o", linestyle="-")
        plt.grid(True)
        plt.show()
    # --- Z coordinates
    z = np.array([])
    z2 = np.array([])
    nk = np.array([])
    dz = np.array([])
    z = np.append(z, geom[0, 0])
    for i in range(1, rows(geom)):
        z = np.append(z, geom[i - 1, 1])
    # adding the last part
    z = np.append(z, geom[i, 1])

    for i in range(0, rows(geom)):
        if (z[i + 1] - z[i]) / max_z < 1.0:
            DZ = 1.0
            dz = np.append(dz, (z[i + 1] - z[i]))
            nk = np.append(nk, DZ)
        else:
            DZ = int((z[i + 1] - z[i]) / max_z)
            dz = np.append(dz, (z[i + 1] - z[i]) / DZ)
            nk = np.append(nk, int((z[i + 1] - z[i]) / max_z))

    for i in range(0, rows(nk)):
        z2 = np.append(z2, divide_z(z[i], z[i + 1], nk[i]))
    # delete repetitions
    # due to values are rounded values ​​with small differences may be repeated,
    z2hlp = skip_small_differences(z2, 1.0e-4, 4)
    # we ignore them here and round values to E-4 meter (4-th digit)
    z3 = np.array([])
    for element in z2hlp:
        if element not in z3:
            z3 = np.append(z3, element)
    z3 = np.sort(z3, axis=0)
    #
    if 0:  # if = 1, show the grid as x-y plot, by Z
        plt.plot(z3, z3, marker="o", linestyle="-")
        plt.plot(r2, r2, marker="x", linestyle="dashdot")
        plt.grid(True)
        plt.show()
    # R nodes coordinates versus R, Z nodes coordinates cersus Z
    R = r2
    Z = z3
    nNod = rows(R) * rows(Z)  # number of nodes in the grid
    out = []
    out.append(R)
    out.append(Z)
    # print(R)
    # print(Z)
    #
    return out

    r = np.zeros((2 * rows(geom)) + 1, dtype=float)

    r[0] = 0.0
    for i in range(0, rows(geom)):
        r[i + 1] = geom[i, 4]
    for i in range(rows(geom), 2 * rows(geom)):
        r[i + 1] = geom[i - rows(geom), 6]

    r2 = np.array([])
    r2 = np.unique(r)
    r2 = np.sort(r2, axis=0)
    rNext = divide_from_to(
        r2[rows(r2) - 1],
        geom[0, cols(geom) - 1],
        (r2[rows(r2) - 1] - r2[rows(r2) - 2]),
        nr,
    )
    for i in range(1, rows(rNext)):
        r2 = np.append(r2, rNext[i])

    if 0:
        plt.plot(r2, r2, marker="o", linestyle="-")
        plt.grid(True)
        plt.show()
    z = np.array([])
    z2 = np.array([])
    nk = np.array([])
    dz = np.array([])
    z = np.append(z, geom[0, 0])
    for i in range(1, rows(geom)):
        z = np.append(z, geom[i - 1, 1])
    z = np.append(z, geom[i, 1])
    for i in range(0, rows(geom)):
        if (z[i + 1] - z[i]) / max_z < 1.0:
            DZ = 1.0
            dz = np.append(dz, (z[i + 1] - z[i]))
            nk = np.append(nk, DZ)
        else:
            DZ = int((z[i + 1] - z[i]) / max_z)
            dz = np.append(dz, (z[i + 1] - z[i]) / DZ)
            nk = np.append(nk, int((z[i + 1] - z[i]) / max_z))
    for i in range(0, rows(nk)):
        z2 = np.append(z2, divide_z(z[i], z[i + 1], nk[i]))
    z2hlp = skip_small_differences(z2, 1.0e-4, 4)

    z3 = np.array([])
    for element in z2hlp:
        if element not in z3:
            z3 = np.append(z3, element)
    z3 = np.sort(z3, axis=0)

    if 0:
        plt.plot(z3, z3, marker="o", linestyle="-")
        plt.plot(r2, r2, marker="x", linestyle="dashdot")
        plt.grid(True)
        plt.show()

    R = r2
    Z = z3
    # nNod = rows(R) * rows(Z)
    out = []
    out.append(R)
    out.append(Z)

    return out


def nodes_to_elements(nc, r, z):
    n2e = np.zeros([(len(r) - 1) * (len(z) - 1), 4], dtype=int)
    for i in range(0, len(n2e)):
        for j in range(0, 4):
            n2e[i, j] = -1

    nCols = 1
    nel_in_col = 0
    for i in range(0, len(n2e)):
        nel_in_col = nel_in_col + 1
        if nel_in_col > (len(z) - 1):
            nel_in_col = 1
            nCols = nCols + 1
        n2e[i, 0] = ((nCols - 1) * len(z) + (nel_in_col - 1) + len(z) + 2) - 1
        n2e[i, 1] = ((nCols - 1) * len(z) + (nel_in_col - 1) + len(z) + 1) - 1
        n2e[i, 2] = ((nCols - 1) * len(z) + (nel_in_col - 1) + 1) - 1
        n2e[i, 3] = ((nCols - 1) * len(z) + (nel_in_col - 1) + 2) - 1

    return n2e


def find_node(node_number, n2e):
    indices = np.where(n2e == node_number)
    no_el = np.full(4, -1, dtype=int)
    for i, j in zip(*indices):
        if 0 <= j < 4:
            no_el[j] = i
    return no_el


def nodes_in_elements(n2e):
    n_in_el = np.zeros([n2e[-1][0] + 1, 4], dtype=int)
    for i in range(0, n2e[-1][0] + 1):
        temp = find_node(i, n2e)
        for j in range(0, 4):
            n_in_el[i, j] = temp[j]
    return n_in_el


def nodes_coordinates(r, z):
    nNod = rows(r) * rows(z)
    nc = np.zeros((nNod, 3))
    n = -1
    for i in range(0, rows(r)):
        for j in range(0, rows(z)):
            n = n + 1
            nc[n, 0] = n
            nc[n, 1] = r[i]
            nc[n, 2] = z[j]
    return nc


def initial_temp(nc, IN):
    th0 = np.zeros([rows(IN), 2], dtype=float)
    for i in range(0, rows(IN) - 1):
        th0[i, 0] = IN[i, 0]
        th0[i, 1] = IN[i, 2]
    th0[rows(IN) - 1, 0] = IN[rows(IN) - 1, 1]
    th0[rows(IN) - 1, 1] = IN[rows(IN) - 1, 3]

    t = np.zeros([rows(nc)], dtype=float)
    t[0] = IN[0, 2]
    t[len(nc) - 1] = IN[len(IN) - 1, 3]
    for i in range(1, rows(nc) - 1):
        t[i] = interp_1d(th0, round(nc[i, 2], 10))
    return t


def a_mat_fix(file_name):
    data = pd.read_excel(
        file_name, sheet_name="mat", usecols="A:D", skiprows=0
    )
    out = data.values

    a = np.zeros([len(out), 3], dtype=float)

    for i in range(0, len(out)):
        a[i, 0] = i
        a[i, 1] = out[i, 1] / (out[i, 2] * out[i, 3])
        a[i, 2] = out[i, 1]

    return a


def working_scheme(file_name):
    data = pd.read_excel(
        file_name, sheet_name="work", usecols="A:D", skiprows=0
    )
    out = data.values
    return out


def read_others(file_name):
    data = pd.read_excel(
        file_name, sheet_name="others", usecols="B:B", skiprows=0
    )
    out = data.values
    return out


def check_mat(IN, Grid, n2e, nc):
    nEl = (len(Grid[0]) - 1) * (len(Grid[1]) - 1)
    mat = np.zeros(nEl, dtype=int)
    for nz in range(0, len(IN)):
        for ne in range(0, nEl):
            if (IN[nz, 1] >= nc[n2e[ne, 0], 2]) and (
                IN[nz, 0] <= nc[n2e[ne, 1], 2]
            ):
                if (nc[n2e[ne, 2], 1] >= IN[nz, 4]) and (
                    nc[n2e[ne, 0], 1] <= IN[nz, 6]
                ):
                    mat[ne] = IN[nz, 5]
                elif (nc[n2e[ne, 2], 1] >= IN[nz, 6]) and (
                    nc[n2e[ne, 0], 1] <= IN[nz, 8]
                ):
                    mat[ne] = IN[nz, 7]
    return mat


def first_type_boundary(nc, n2e, mat):
    FT = np.empty([], dtype=int)
    FT2 = np.empty([], dtype=int)
    NB = np.zeros([4], dtype=int)
    for i in range(1, len(nc)):
        NB = find_node(i, n2e)
        if NB[0] < 0 or NB[1] < 0 or NB[2] < 0 or NB[3] < 0:
            FT = np.append(FT, i)
        if (
            mat[find_node(i, n2e)[0]] == 0
            or mat[find_node(i, n2e)[1]] == 0
            or mat[find_node(i, n2e)[2]] == 0
            or mat[find_node(i, n2e)[3]] == 0
        ):
            FT = np.append(FT, i)
    FT2 = list(set(FT))
    return FT2


def nodes_zones(Z, nc, n2e, mat):
    WF = np.empty([], dtype=int)
    NZ = []

    for i in range(0, len(nc)):
        if (
            mat[find_node(i, n2e)[0]] == 0
            or mat[find_node(i, n2e)[1]] == 0
            or mat[find_node(i, n2e)[2]] == 0
            or mat[find_node(i, n2e)[3]] == 0
        ):
            WF = np.append(WF, i)
    WF = list(set(WF))

    for i in range(0, len(Z)):
        X = []
        for j in range(0, len(nc) - 1):
            if nc[j, 2] == Z[i] and j in WF:
                X.append(j)
        NZ.append(X)
    if 0:
        nH2O = []
        j = 0
        for n in range(0, len(NZ)):
            for k in range(0, len(NZ[n])):
                nH2O.append(NZ[n][k])
        nH2O = sorted(nH2O)
        print(nH2O)

    return NZ


def eq_temp_in_zone(t, NinZ):
    for i in range(0, len(NinZ) - 1):
        for j in range(0, len(NinZ[i]) - 1):
            t[NinZ[i][j]] = t[NinZ[i][len(NinZ[i]) - 1]]
    return t


def temp_after_dtime(nw, n2e, nc, mat, am, t, dtime, nodes_in_elements):
    NinE = nodes_in_elements[nw]

    nw0 = n2e[NinE[1], 0]
    nw1 = n2e[NinE[2], 1]
    nw2 = n2e[NinE[0], 1]
    nw3 = n2e[NinE[0], 3]

    dx0 = nc[nw0, 2] - nc[nw, 2]
    dx1 = nc[nw1, 1] - nc[nw, 1]
    dx2 = nc[nw, 2] - nc[nw2, 2]
    dx3 = nc[nw, 1] - nc[nw3, 1]

    a_vals = [am[mat[NinE[i]], 1] for i in range(4)]

    dtspeed = fdm.dtdtime(
        [t[nw0], t[nw1], t[nw2], t[nw3], t[nw]],
        [dx0, dx1, dx2, dx3],
        a_vals,
        nc[nw, 2],
    )

    dt = dtspeed * dtime
    tend = t[nw] + dt
    return tend


def well_param(wws, time):
    A = np.zeros([3], dtype=float)
    n = 0
    while wws[n, 0] <= time:
        n = n + 1
    n = n - 1
    A[0] = time
    A[1] = wws[n, 1]
    A[2] = wws[n, 2]
    return A


def temp_in_well(node_in_zone, t, nc):
    tz = np.empty([len(node_in_zone), 4], dtype=float)
    for i in range(0, len(node_in_zone)):
        tz[i, 0] = node_in_zone[i][len(node_in_zone[i]) - 1]
        tz[i, 1] = nc[node_in_zone[i][len(node_in_zone[i]) - 1], 1]
        tz[i, 2] = nc[node_in_zone[i][len(node_in_zone[i]) - 1], 2]
        tz[i, 3] = t[node_in_zone[i][len(node_in_zone[i]) - 1]]
    return tz


def max_dtime(nc, n2e, mat, am, quality):
    if quality <= 4.0:
        quality = 4.0

    dtime = np.zeros([len(mat)], dtype=float)
    for n in range(0, len(mat)):
        if mat[n] == 0:
            dtime[n] = 1.0e10
        else:
            Dr = nc[n2e[n, 0], 1] - nc[n2e[n, 3], 1]
            a = am[mat[n], 1]
            dtime[n] = Dr**2.0 / (quality * a)
    out = min(dtime)
    return out


def nodes_filled_with_brine(node_in_well, n2e, mat):
    NFB = []
    for i in range(0, len(node_in_well)):
        hlp = []
        for j in range(0, len(node_in_well[i])):
            nodeNumber = node_in_well[i][j]
            noEls = find_node(nodeNumber, n2e)
            onlyBrine = True
            for k in range(0, 4):
                if (noEls[k] >= 0) and onlyBrine:
                    if mat[noEls[k]] > 0:
                        onlyBrine = False
            if onlyBrine:
                hlp.append(node_in_well[i][j])
        NFB.append(hlp)
    return NFB


def disribute_temp_changes(
    t0, twellNew, twell0, node_in_zone, node_with_brine
):
    tnew = np.zeros([len(t0)], dtype=float)

    for i in range(0, len(t0)):
        tnew[i] = t0[i]

    for k in range(0, len(twell0)):
        NodeNumber = int(twell0[k, 0])
        tnew[NodeNumber] = twellNew[k, 1]

    tnew = eq_temp_in_zone(tnew, node_in_zone)

    for i in range(0, len(twellNew)):
        for j in range(0, len(node_with_brine[i])):
            NodeNumber = node_with_brine[i][j]
            tnew[NodeNumber] = twellNew[i, 1]
    return tnew


def prod_calculations(
    wws,
    dtime,
    dtime_max,
    t0,
    twell0,
    n2e,
    nc,
    Z,
    R,
    mat,
    am,
    NinWell,
    NfbB,
    actvN,
    others,
    nodes_top,
    nodes_in_elements,
):
    tHeadTime = []
    show_graphic_results = 1
    PipeRoughness = others[1]
    S = others[0]
    accuracy = others[2]
    DtimeDynamic = dtime
    time = 0.0
    n = 0
    k = 0
    twell_1 = twell0
    counter = 0
    list_excel = []
    while (time + DtimeDynamic) <= wws[len(wws) - 1, 0]:
        if (time + dtime) > wws[k, 0]:
            DtimeDynamic = wws[k, 0] - time
        else:
            DtimeDynamic = dtime
        time = time + DtimeDynamic

        tHeadTime.append([time, twell_1[0, 3]])

        twellNe = pw.prod_well(
            wws[k, 1],
            wws[k, 2],
            wws[k, 3],
            S,
            twell0,
            t0,
            DtimeDynamic,
            PipeRoughness,
            accuracy,
            n2e,
            nc,
            Z,
            mat,
            am,
            nodes_in_elements,
        )
        twellNew = twellNe[0]
        PR = twellNe[1]
        PR.insert(0, time)
        list_excel.append(PR)

        twellNew[0, 0] = twellNew[0, 1]

        tnew = disribute_temp_changes(t0, twellNew, twell0, NinWell, NfbB)
        for nw in range(0, len(actvN)):
            tnew[actvN[nw]] = temp_after_dtime(
                actvN[nw],
                n2e,
                nc,
                mat,
                am,
                tnew,
                DtimeDynamic,
                nodes_in_elements,
            )
        t0 = tnew

        for nw in nodes_top:
            t0[nw] = t0[nw + 1]
            tnew[nw] = tnew[nw + 1]

        twell_1 = pw.twell_convert_twell0(twell0, twellNew)

        twell0 = twell_1

        if (time + DtimeDynamic) >= wws[len(wws) - 1, 0]:
            print(
                "- 100.00% of expected time (",
                round(time / 3600.0, 2),
                "hr) | time step length: ",
                round(DtimeDynamic, 2),
                " sec | bottom / wellhead temp. ",
                round(twell_1[len(twell_1) - 1, 3], 2),
                "/",
                round(twell_1[0, 3], 2),
                "°C",
            )
        else:
            print(
                "- ",
                round(100.0 * time / wws[len(wws) - 1, 0], 2),
                "% of time (",
                round(time / 3600.0, 2),
                "hr) | time step length: ",
                round(DtimeDynamic, 2),
                " sec | bottom / wellhead temp. ",
                round(twell_1[len(twell_1) - 1, 3], 2),
                "/",
                round(twell_1[0, 3], 2),
                "°C",
            )

        n = n + 1
        if time >= wws[k, 0]:
            k = k + 1
        counter = counter + 1
        if counter == 5:
            dtime = dtime + 10.0
            counter = 0
        if dtime > dtime_max:
            dtime = dtime_max
        result = pd.DataFrame(
            list_excel,
            columns=[
                "time",
                "temperature",
                "pressure",
                "velocity",
                "pDynamic",
                "density",
                "spec heat",
                "press_drop",
                "mass flow",
                "alfaB",
            ],
        )
        result.to_excel("wellhead_out_last_time.xlsx", index=False)

    if show_graphic_results:
        tHead_time = np.array(tHeadTime)
        fig1, ax = plt.subplots()
        ax.plot(tHead_time[:, 0], tHead_time[:, 1], marker="o", linestyle="-")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Wellhead temperature [°C]")
        ax.set_title("Wellhead temperature vs time")
        ax.grid(True)
    if show_graphic_results:
        tvsH = np.zeros([len(Z), 2])
        n = -1
        for nw in range(0, len(Z)):
            n = n + 1
            tvsH[n, 0] = nc[nw, 2] * (-1.0)
            tvsH[n, 1] = twell0[nw, 3]
        fig2, ax = plt.subplots()
        ax.plot(tvsH[:, 1], tvsH[:, 0], marker="o", linestyle="-")
        ax.set_xlabel("Temperature [°C]")
        ax.set_ylabel("Depth [m]")
        ax.set_title("Temperature vs depth at the end time")
        ax.grid(True)

    if show_graphic_results:
        nR = int(len(Z) / 3)
        if nR > len(Z) - 1:
            nR = len(Z) - 1
        tvsR = np.zeros([len(R), 2])
        n = -1
        for nw in range(0, len(R)):
            n = n + 1
            nwr = nR + nw * len(Z)
            tvsR[n, 0] = nc[nwr, 1]
            tvsR[n, 1] = tnew[nwr]
        fig3, ax = plt.subplots()
        ax.plot(tvsR[:, 0], tvsR[:, 1], marker="o", linestyle="-")
        ax.set_xlabel("Radius [m]")
        ax.set_ylabel("Temperature [°C]")
        ax.set_title(
            "Temperature of rock at depth of "
            + str(round(nc[nR, 2], 2))
            + " m at the end"
        )
        ax.grid(True)
    if show_graphic_results:
        x = np.zeros([len(R), len(Z)])
        y = np.zeros([len(R), len(Z)])
        z = np.zeros([len(R), len(Z)])
        n = -1
        for a in range(0, len(R)):
            for b in range(0, len(Z)):
                n = n + 1
                x[a, b] = nc[n, 1]
                y[a, b] = (nc[n, 2] - min(Z)) * -1.0
                z[a, b] = tnew[n]
        fig4, ax = plt.subplots()
        contour = ax.contourf(x, y, z, levels=25)
        ax.set_xlabel("Radius [m]")
        ax.set_ylabel("Depth [m]")
        ax.set_title("Temperature distribution vs R and H [°C]")
        ax.grid(True)
        cbar = fig4.colorbar(contour, ax=ax)
        cbar.set_label("Temperature [°C]")
        plt.show()
    return tnew


def inj_well_guess_press_head(
    m,
    t,
    pBottom,
    S,
    tr,
    tRockMatrix,
    dtime,
    fi,
    acc,
    n2e,
    nc,
    Z,
    mat,
    am,
    nodes_top,
    nodes_in_elements,
):
    pHeadMin = 101325.0
    pHeadMax = pBottom
    pHeadAver = (pHeadMin + pHeadMax) * 0.5

    pBottomMin = pw.inj_well(
        m,
        t,
        pHeadMin,
        S,
        tr,
        tRockMatrix,
        dtime,
        fi,
        acc,
        n2e,
        nc,
        Z,
        mat,
        am,
        nodes_top,
        nodes_in_elements,
    )[len(Z) - 1, 5]
    pBottomMax = pw.inj_well(
        m,
        t,
        pHeadMax,
        S,
        tr,
        tRockMatrix,
        dtime,
        fi,
        acc,
        n2e,
        nc,
        Z,
        mat,
        am,
        nodes_top,
        nodes_in_elements,
    )[len(Z) - 1, 5]
    pBottomAver = pw.inj_well(
        m,
        t,
        pHeadAver,
        S,
        tr,
        tRockMatrix,
        dtime,
        fi,
        acc,
        n2e,
        nc,
        Z,
        mat,
        am,
        nodes_top,
        nodes_in_elements,
    )[len(Z) - 1, 5]

    while (pBottomMax - pBottomMin) / pBottomAver > acc:
        if pBottomAver >= pBottom:
            pHeadMax = pHeadAver
        if pBottomAver < pBottom:
            pHeadMin = pHeadAver
        pHeadAver = (pHeadMin + pHeadMax) * 0.5
        pBottomMin = pw.inj_well(
            m,
            t,
            pHeadMin,
            S,
            tr,
            tRockMatrix,
            dtime,
            fi,
            acc,
            n2e,
            nc,
            Z,
            mat,
            am,
            nodes_top,
            nodes_in_elements,
        )[len(Z) - 1, 5]
        pBottomMax = pw.inj_well(
            m,
            t,
            pHeadMax,
            S,
            tr,
            tRockMatrix,
            dtime,
            fi,
            acc,
            n2e,
            nc,
            Z,
            mat,
            am,
            nodes_top,
            nodes_in_elements,
        )[len(Z) - 1, 5]
        pBottomAver = pw.inj_well(
            m,
            t,
            pHeadAver,
            S,
            tr,
            tRockMatrix,
            dtime,
            fi,
            acc,
            n2e,
            nc,
            Z,
            mat,
            am,
            nodes_top,
            nodes_in_elements,
        )[len(Z) - 1, 5]
    pHead = pHeadAver
    return pHead


def inj_calculations(
    wws,
    dtime,
    dtime_max,
    t0,
    twell0,
    n2e,
    nc,
    Z,
    R,
    mat,
    am,
    NinWell,
    NfbB,
    actvN,
    others,
    nodes_top,
    nodes_bottom,
    nodes_in_elements,
):
    tHeadTime = []
    pHeadTime = []
    show_graphic_results = 1

    for i in range(0, len(wws)):
        wws[i, 1] = wws[i, 1] * (-1.0)

    PipeRoughness = others[1]
    S = others[0]
    accuracy = others[2]

    DtimeDynamic = dtime
    time = 0.0
    n = 0
    k = 0
    twell_1 = twell0
    counter = 0

    while (time + DtimeDynamic) <= wws[len(wws) - 1, 0]:
        if (time + dtime) > wws[k, 0]:
            DtimeDynamic = wws[k, 0] - time
        else:
            DtimeDynamic = dtime

        time = time + DtimeDynamic

        if (time + DtimeDynamic) >= wws[len(wws) - 1, 0]:
            print(
                "- 100.00% of expected time (",
                round(time / 3600.0, 2),
                "hr) | time step length: ",
                round(DtimeDynamic, 2),
                " sec | wellhead / bottom temp. ",
                round(twell_1[0, 3], 2),
                "/",
                round(twell_1[len(twell_1) - 1, 3], 2),
                "°C",
            )
        else:
            print(
                "- ",
                round(100.0 * time / wws[len(wws) - 1, 0], 2),
                "% of time (",
                round(time / 3600.0, 2),
                "hr) | time step length: ",
                round(DtimeDynamic, 2),
                " sec | wellhead / bottom temp. ",
                round(twell_1[0, 3], 2),
                "/",
                round(twell_1[len(twell_1) - 1, 3], 2),
                "°C",
            )

        tHeadTime.append([time, twell_1[len(Z) - 1, 3]])

        pHead = inj_well_guess_press_head(
            wws[k, 1],
            wws[k, 2],
            wws[k, 3],
            S,
            twell0,
            t0,
            DtimeDynamic,
            PipeRoughness,
            accuracy,
            n2e,
            nc,
            Z,
            mat,
            am,
            nodes_top,
            nodes_in_elements,
        )
        pHeadTime.append([time, pHead])
        twellNew = pw.inj_well(
            wws[k, 1],
            wws[k, 2],
            pHead,
            S,
            twell0,
            t0,
            DtimeDynamic,
            PipeRoughness,
            accuracy,
            n2e,
            nc,
            Z,
            mat,
            am,
            nodes_top,
            nodes_in_elements,
        )
        tnew = disribute_temp_changes(t0, twellNew, twell0, NinWell, NfbB)
        for nw in range(0, len(actvN)):
            tnew[actvN[nw]] = temp_after_dtime(
                actvN[nw],
                n2e,
                nc,
                mat,
                am,
                tnew,
                DtimeDynamic,
                nodes_in_elements,
            )
        t0 = tnew

        for nw in nodes_top:
            t0[nw] = t0[nw + 1]
            tnew[nw] = tnew[nw + 1]
        for nw in nodes_bottom:
            t0[nw] = t0[nw - 1]
            tnew[nw] = tnew[nw - 1]

        twell_1 = pw.twell_convert_twell0(twell0, twellNew)
        twell0 = twell_1

        n = n + 1
        if time >= wws[k, 0]:
            k = k + 1

        counter = counter + 1
        if counter == 5:
            dtime = dtime + 5.0
            counter = 0
        if dtime > dtime_max:
            dtime = dtime_max

    if show_graphic_results:
        tHead_time = np.array(tHeadTime)
        fig1, ax = plt.subplots()
        ax.plot(tHead_time[:, 0], tHead_time[:, 1], marker="o", linestyle="-")
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Temperature at the liner depth [°C]")
        ax.set_title("Temperature at the liner depth vs time")
        ax.grid(True)
    if show_graphic_results:
        pHead_time = np.array(pHeadTime)
        fig2, ax = plt.subplots()
        ax.plot(
            pHead_time[:, 0],
            pHead_time[:, 1] / 1.0e6,
            marker="o",
            linestyle="-",
        )
        ax.set_xlabel("Time [s]")
        ax.set_ylabel("Pressure on the wellhead [MPa]")
        ax.set_title("Pressure on the wellhead vs time")
        ax.grid(True)
    if show_graphic_results:
        tvsH = np.zeros([len(Z), 2])
        n = -1
        for nw in range(0, len(Z)):
            n = n + 1
            tvsH[n, 0] = nc[nw, 2] * (-1.0)
            tvsH[n, 1] = tnew[nw]
        fig4, ax = plt.subplots()
        ax.plot(tvsH[:, 1], tvsH[:, 0], marker="o", linestyle="-")
        ax.set_xlabel("Temperature [°C]")
        ax.set_ylabel("Depth [m]")
        ax.set_title("Temperature vs depth at the end")
        ax.grid(True)

    if show_graphic_results:
        nR = int(len(Z) / 3)
        if nR > len(Z) - 1:
            nR = len(Z) - 1
        tvsR = np.zeros([len(R), 2])
        n = -1
        for nw in range(0, len(R)):
            n = n + 1
            nwr = nR + nw * len(Z)
            tvsR[n, 0] = nc[nwr, 1]
            tvsR[n, 1] = tnew[nwr]
        fig4, ax = plt.subplots()
        ax.plot(tvsR[:, 0], tvsR[:, 1], marker="o", linestyle="-")
        ax.set_xlabel("Radius [m]")
        ax.set_ylabel("Temperature [°C]")
        ax.set_title(
            "Temperature or rock at depth of "
            + str(round(nc[nR, 2], 2))
            + " m at the end"
        )
        ax.grid(True)
    if show_graphic_results:
        x = np.zeros([len(R), len(Z)])
        y = np.zeros([len(R), len(Z)])
        z = np.zeros([len(R), len(Z)])
        n = -1
        for a in range(0, len(R)):
            for b in range(0, len(Z)):
                n = n + 1
                x[a, b] = nc[n, 1]
                y[a, b] = (nc[n, 2] - min(Z)) * -1.0
                z[a, b] = tnew[n]
        fig5, ax = plt.subplots()
        contour = ax.contourf(x, y, z, levels=25)
        ax.set_xlabel("Radius [m]")
        ax.set_ylabel("Depth [m]")
        ax.set_title("Temperature distribution vs R and h [°C]")
        ax.grid(True)
        cbar = fig5.colorbar(contour, ax=ax)
        cbar.set_label("Temperature [°C]")
        plt.show()
    return tnew
